!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
*=====================================================================*
      SUBROUTINE CCREO2CON(IREO,NREO,IRDOTS,RCONS,
     &                     MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate contributions from second-order orbital
*              reorthogonalization and the coupling of reorthog.
*              and relaxation to response functions involving
*              perturbation-dependent basis sets.
*
*              IREO    --  array with the perturbation indeces
*              NREO    --  length of IREO
*              IRDOTS  --  matrix containing the second indeces
*              RCONS   --  matrix with the contributions
*              MXVEC   --  leading dimension of IRDOTS and RCONS
*
*     Christof Haettig 11-6-1999
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccr1rsp.h>
#include <ccfro.h>
#include <ccroper.h>
#include <iratdef.h>
#include <inftap.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYM0, LUFCK
      PARAMETER( ISYM0 = 1)
      CHARACTER LABEL0*(8)
      PARAMETER( LABEL0 = 'HAM0    ' )

      LOGICAL LORX, LFOCK0
      INTEGER LWORK, NREO, MXVEC

      INTEGER IREO(NREO), IRDOTS(MXVEC,NREO)

#if defined (SYS_CRAY)
      REAL TWO, FREQ, RCONS(MXVEC,NREO), WORK(LWORK)
#else
      DOUBLE PRECISION TWO, FREQ, RCONS(MXVEC,NREO), WORK(LWORK)
#endif
      PARAMETER (TWO = 2.0D0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA, LABELB
      LOGICAL LORXA, LORXB, LPDBSA, LPDBSB, NOKAPPA
      INTEGER NCMOT, NASHT, N2ASHX, LCINDX, IDXREO, IKAPPA
      INTEGER KCMO, KFOCK0, KOVERLP, KEND1, LWRK1, IFOCK, IADRF, ISYM
      INTEGER ISYMA, ISYMB, IOPERA, IOPERB, KR2EFF, KRMATA, KAPASQ
      INTEGER KAPPAA, KAPBSQ, KAPPAB, KSCR1, KEND2, LWRK2, KRMATB
      INTEGER KQMATPA, KQMATPB
      INTEGER IOPT, IVEC, IKAPPB, ISAMA, ISAMB

* external functions:
      INTEGER IEFFFOCK
      INTEGER ILSTSYM, ILSTSYMRLX
      INTEGER IROPER
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
*     check, if there is anything at all to do:
*---------------------------------------------------------------------*

      IF (NREO.LE.0) RETURN

*---------------------------------------------------------------------*
*     get some constants from sirius common block:
*---------------------------------------------------------------------*

      CALL CC_SIRINF(NCMOT,NASHT,N2ASHX,LCINDX) 

*---------------------------------------------------------------------*
*     allocate memory for perturbation-independent stuff needed:
*---------------------------------------------------------------------*
      KCMO    = 1
      KFOCK0  = KCMO    + NCMOT
      KOVERLP = KFOCK0  + N2BST(ISYM0)
      KEND1   = KOVERLP + N2BST(ISYM0)
      LWRK1   = LWORK   - KEND1

      IF (LWRK1.LT.0) THEN
         CALL QUIT('Insufficient work space in CCREO2CON.')
      END IF 

*---------------------------------------------------------------------*
*     read MO coefficients from file:
*---------------------------------------------------------------------*

      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NCMOT,WORK(KCMO))
      CALL GPCLOSE(LUSIFC,'KEEP')                  

*---------------------------------------------------------------------*
*     read overlap matrix from file:
*---------------------------------------------------------------------*

      IF (LWRK1.LT.NBAST) THEN
         CALL QUIT('Insufficient work space in CCREO2CON.')
      END IF 

      CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND1),NBAST) 
      CALL CCSD_SYMSQ(WORK(KEND1),ISYM0,WORK(KOVERLP))

*---------------------------------------------------------------------*
*     read zeroth-order effective Fock matrix if available:
*---------------------------------------------------------------------*
      IFOCK = IEFFFOCK(LABEL0,ISYM,1)

      IF (LEXPFCK(2,IFOCK)) THEN
         IADRF = IADRFCK(1,IFOCK)

         LUFCK = -1
         CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
         CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK0),IADRF,N2BST(ISYM0))
         CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')   

         CALL CC_EFFCKMO(WORK(KFOCK0),ISYM0,WORK(KCMO),WORK(KOVERLP),
     &                   WORK(KEND1),LWRK1)        

         LFOCK0 = .TRUE.
      ELSE
         LFOCK0 = .FALSE.
      END IF

*---------------------------------------------------------------------*
*     start loop over first index:
*---------------------------------------------------------------------*
      DO IDXREO = 1, NREO

         IKAPPA = IREO(IDXREO)
         ISYMA  = ILSTSYMRLX('R1',IKAPPA)
         LABELA = LRTHFLBL(IKAPPA)
         IOPERA = IROPER(LABELA,ISYMA)
         ISAMA  = ISYMAT(IOPERA)
         LORXA  = .TRUE.
         LPDBSA = LPDBSOP(IOPERA)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCREO2CON> IDXREO:',IDXREO
            WRITE (LUPRI,*) 'CCREO2CON> IKAPPA:',IKAPPA
            WRITE (LUPRI,*) 'CCREO2CON> LABELA:',LABELA
            WRITE (LUPRI,*) 'CCREO2CON> LORXA :',LORXA
            WRITE (LUPRI,*) 'CCREO2CON> LPDBSA:',LPDBSA
         END IF

         KR2EFF  = KEND1
         KQMATPA = KR2EFF  + N2BST(ISYM0)
         KRMATA  = KQMATPA + N2BST(ISYMA)
         KAPASQ  = KRMATA  + N2BST(ISYMA)
         KAPPAA  = KAPASQ  + N2BST(ISYMA)
         KQMATPB = KAPPAA  + 2*NALLAI(ISYMA)
         KRMATB  = KQMATPB + N2BST(ISYMA)
         KAPBSQ  = KRMATB  + N2BST(ISYMA)
         KAPPAB  = KAPBSQ  + N2BST(ISYMA)
         KSCR1   = KAPPAB  + NALLAI(ISYMA)
         KEND2   = KSCR1   + N2BST(ISYMA)

         LWRK2  = LWORK  - KEND2

         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CCREO2CON.')
         END IF
 
         IF (LORXA) THEN
           CALL CC_RDHFRSP('R1 ',IKAPPA,ISYMA,WORK(KAPPAA))
           CALL CCKAPPASQ(WORK(KAPASQ),WORK(KAPPAA),ISYMA,'N')
         END IF

         IF (LPDBSA) THEN
           CALL CC_GET_RMAT(WORK(KSCR1),IOPERA,1,ISYMA,
     &                      WORK(KEND2),LWRK2)

           NOKAPPA = .TRUE.
           CALL CC_QMAT(WORK(KQMATPA),WORK(KRMATA),WORK(KSCR1),DUMMY,
     &                ISAMA,ISYMA,NOKAPPA,WORK(KCMO),WORK(KEND2),LWRK2)
         END IF

*---------------------------------------------------------------------*
*        loop over second index:
*---------------------------------------------------------------------*
         IVEC = 1

         DO WHILE (IRDOTS(IVEC,IDXREO).NE.0 .AND. IVEC.LE.MXVEC)

C          WRITE (LUPRI,*) 'CCREO2CON> IVEC = ',IVEC
      
           IKAPPB = IRDOTS(IVEC,IDXREO)
           ISYMB  = ILSTSYMRLX('R1',IKAPPB)
           LABELB = LRTHFLBL(IKAPPB)
           IOPERB = IROPER(LABELB,ISYMB)
           ISAMB  = ISYMAT(IOPERB)
           LORXB  = .TRUE.
           LPDBSB = LPDBSOP(IOPERB)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'CCREO2CON> IVEC  :',IVEC
              WRITE (LUPRI,*) 'CCREO2CON> IKAPPB:',IKAPPB
              WRITE (LUPRI,*) 'CCREO2CON> LABELB:',LABELB
              WRITE (LUPRI,*) 'CCREO2CON> LORXB :',LORXB
              WRITE (LUPRI,*) 'CCREO2CON> LPDBSB:',LPDBSB
           END IF

           IF (ISYMB.NE.ISYMA) THEN
              WRITE (LUPRI,*) IDXREO, ISYMA, IVEC, ISYMB
              CALL QUIT('symmetry mismatch in CCREO2CON.')
           END IF

           IF (LORXB) THEN
             CALL CC_RDHFRSP('R1 ',IKAPPB,ISYMB,WORK(KAPPAB))
             CALL CCKAPPASQ(WORK(KAPBSQ),WORK(KAPPAB),ISYMB,'N')
           END IF

           IF (LPDBSB) THEN
             CALL CC_GET_RMAT(WORK(KSCR1),IOPERB,1,ISYMB,
     &                        WORK(KEND2),LWRK2)

             NOKAPPA = .TRUE.
             CALL CC_QMAT(WORK(KQMATPB),WORK(KRMATB),WORK(KSCR1),DUMMY,
     &                ISAMB,ISYMB,NOKAPPA,WORK(KCMO),WORK(KEND2),LWRK2)
           END IF

           CALL DZERO(WORK(KR2EFF),N2BST(ISYM0))
           
           IF (LORXB.AND.LPDBSA) THEN
              CALL CC_MMOMMO('N','N',+1.0D0,WORK(KAPBSQ),ISYMB,
     &                       WORK(KRMATA),ISYMA,1.0D0,WORK(KR2EFF),1)
              CALL CC_MMOMMO('N','N',-1.0D0,WORK(KRMATA),ISYMA,
     &                       WORK(KAPBSQ),ISYMB,1.0D0,WORK(KR2EFF),1)
           END IF

           IF (LORXA.AND.LPDBSB) THEN
              CALL CC_MMOMMO('N','N',+1.0D0,WORK(KAPASQ),ISYMA,
     &                       WORK(KRMATB),ISYMB,1.0D0,WORK(KR2EFF),1)
              CALL CC_MMOMMO('N','N',-1.0D0,WORK(KRMATB),ISYMB,
     &                       WORK(KAPASQ),ISYMA,1.0D0,WORK(KR2EFF),1)
           END IF

           IF (LPDBSA.AND.LPDBSB) THEN
              CALL QUIT('CCREO2CON NOT YET IMPLEMENTED '//
     &             'FOR 2. DERIVATIVES.')
           END IF

           RCONS(IVEC,IDXREO) = TWO*DDOT(N2BST(ISYM0),WORK(KR2EFF),1,
     &                                                WORK(KFOCK0),1)
C          WRITE (LUPRI,*) 'CCREO2CON>',RCONS(IVEC,IDXREO)

           IVEC = IVEC + 1
         END DO

      END DO

*---------------------------------------------------------------------*
*     print the results:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
     &        'CCREO2CON> results for X intermediate contribs.:'
         IF (MXVEC.NE.0) THEN
            DO IDXREO = 1, NREO
               WRITE (LUPRI,*) 'IDXREO = ',IDXREO
               IVEC = 1
               DO WHILE(IRDOTS(IVEC,IDXREO).NE.0 .AND. IVEC.LE.MXVEC)
                  WRITE(LUPRI,'(A,2I5,2X,E18.12)') 'CCREO2CON> ',
     &              IREO(IDXREO),IRDOTS(IVEC,IDXREO),RCONS(IVEC,IDXREO)
                  IVEC = IVEC + 1
               END DO
            END DO
         ELSE
            WRITE (LUPRI,*) 'MXVEC.EQ.0 --> nothing calculated.'
         END IF
      END IF

      RETURN
      END
*======================================================================*
