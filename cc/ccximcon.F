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
      SUBROUTINE CCXIMCON(IXIM,NXIM,LISTQ,IXDOTS,XCONS,
     &                    MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate contributions from the generalized relaxation
*              & reorthogonalization X intermediate to response 
*              functions involving orbital relaxation or perturbation-
*              dependent basis sets.
*
*              IXIM    --  array with the perturbation indeces
*              NXIM    --  length of IXIM
*              LISTQ   --  type of vectors the X intermediates are to
*                          be dotted on
*              IXDOTS  --  matrix of indeces for the vectors with which
*                          the dot products are to be calculated
*              XCONS   --  matrix with the dot product results
*              MXVEC   --  leading dimension of IXDOTS and XCONS
*
*     Christof Haettig 1-6-1999
*
*     N.B.: this routine is not yet adapted for QMATP diff. from QMATH
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
#include <inftap.h>
#include <iratdef.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYM0, LUFCK
      PARAMETER( ISYM0 = 1) 
      CHARACTER LABEL0*(8)
      PARAMETER( LABEL0 = 'HAM0    ' )

      LOGICAL LORX, LFOCK0, LORXQ, LPDBS
      CHARACTER*(3) LISTQ
      INTEGER ISYXIM, LWORK, NXIM, MXVEC

      INTEGER IXIM(NXIM), IXDOTS(MXVEC,NXIM)

#if defined (SYS_CRAY)
      REAL FREQ, XCONS(MXVEC,NXIM), WORK(LWORK)
      REAL HALF, ONE, TWO, ZERO, FTRACE
#else
      DOUBLE PRECISION FREQ, XCONS(MXVEC,NXIM), WORK(LWORK)
      DOUBLE PRECISION HALF, ONE, TWO, ZERO, FTRACE
#endif
      PARAMETER( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0, HALF = 0.5D0 )

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABEL, LABELQ
      LOGICAL NOKAPPA
      INTEGER NCMOT, NASHT, N2ASHX, LCINDX
      INTEGER IFOCK, IEXPV, IADRF, KOVERLP, KFOCK0, KFOCK1, KCMO, KAPPA
      INTEGER IKAPPA, IOPT, IOPER, ISYM, IDXXIM, KOFF1, IFCK1
      INTEGER KEND1,LWRK1,KEND2,LWRK2,KXIM,KSCR1,KRMAT,KQMATP,KQMATH
      INTEGER IVEC, IFILE, ISYMQ, IOPERQ, KEND3, LWRK3
      INTEGER ISYM1, ISYM2, KOFF2, KQTRP, IREAL

* external functions:
      INTEGER IEFFFOCK
      INTEGER IEXPECT
      INTEGER IROPER
      INTEGER I1DXFCK
      INTEGER ILSTSYMRLX
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
*     check, if there is anything at all to do:
*---------------------------------------------------------------------*

      IF (NXIM.LE.0) RETURN

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
         CALL QUIT('Insufficient work space in CCXIMCON.')
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
         CALL QUIT('Insufficient work space in CCXIMCON.')
      END IF 

      CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND1),NBAST) 
      CALL CCSD_SYMSQ(WORK(KEND1),ISYM0,WORK(KOVERLP))

*---------------------------------------------------------------------*
*     read zeroth-order effective Fock matrix if available:
*---------------------------------------------------------------------*
      call quit('hjaaj: code error, ISYM not defined, please fix!')
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

         IF (LOCDBG) THEN
            WRITE (LUPRI,*)
     &           'CCXIMCON> zeroth-order effective Fock in ao/AO:'
            CALL CC_PRONELAO(WORK(KFOCK0),ISYM0)
         END IF       

      ELSE
         LFOCK0 = .FALSE.
      END IF

*---------------------------------------------------------------------*
*     start loop over X intermediates:
*---------------------------------------------------------------------*
      DO IDXXIM = 1, NXIM

         IKAPPA = IXIM(IDXXIM)
         LABEL  = LRTHFLBL(IKAPPA)
         IOPER  = IROPER(LABEL,ISYXIM)
         LORX   = .TRUE.
         LPDBS  = LPDBSOP(IOPER)
         ISYXIM = ILSTSYMRLX('R1',IKAPPA)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCXIMCON> IDXXIM:',IDXXIM
            WRITE (LUPRI,*) 'CCXIMCON> IKAPPA:',IKAPPA
            WRITE (LUPRI,*) 'CCXIMCON> LABEL :',LABEL
            WRITE (LUPRI,*) 'CCXIMCON> LORX  :',LORX
         END IF

         KXIM   = KEND1
         KRMAT  = KXIM   + N2BST(ISYXIM)
         KQMATP = KRMAT  + N2BST(ISYXIM)
         KQMATH = KQMATP + N2BST(ISYXIM)
         KQTRP  = KQMATH + N2BST(ISYXIM)
         KAPPA  = KQTRP  + N2BST(ISYXIM)
         KEND2  = KAPPA  + 2*NALLAI(ISYXIM)

         LWRK2  = LWORK  - KEND2

         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CCXIMCON.')
         END IF
 
         ! read first-order effective Fock matrix from file
         call quit('hjaaj: code error, ISYM not defined, please fix!')
         IFOCK = IEFFFOCK(LABEL,ISYM,1)
         IEXPV = IEXPECT(LABEL,ISYM)
         IADRF = IADRFCK(1,IFOCK)  

         CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
         CALL GETWA2(LUFCK,FILFCKEFF,WORK(KXIM),IADRF,N2BST(ISYXIM))
         CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')      

         IF (LOCDBG) THEN
            FTRACE = ZERO
            IF (ISYXIM.EQ.1) THEN
               DO ISYM = 1, NSYM
                  KOFF1 = KXIM + IAODIS(ISYM,ISYM)
                  DO I = 1, NBAS(ISYM)
                    FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                  END DO
               END DO
            END IF
            WRITE (LUPRI,*) 'LABEL:',LABEL
            WRITE (LUPRI,*) 'ISYXIM,IFOCK,IEXPV:',ISYXIM,IFOCK,IEXPV
            WRITE (LUPRI,*) 'FTRACE of read matrix:',FTRACE
            WRITE (LUPRI,*) 'one-electron expect:',EXPVALUE(1,IEXPV)
            WRITE (LUPRI,*) 'two-electron expect:',EXPVALUE(2,IEXPV)
         END IF
                            
         ! transform effective Fock matrix to MO basis
         CALL CC_EFFCKMO(WORK(KXIM),ISYXIM,WORK(KCMO),WORK(KOVERLP),
     &                   WORK(KEND2),LWRK2)  

         IF (LOCDBG) THEN
            FTRACE = ZERO
            IF (ISYXIM.EQ.1) THEN
               DO ISYM = 1, NSYM
                  KOFF1 = KXIM + IAODIS(ISYM,ISYM)
                  DO I = 1, NBAS(ISYM)
                    FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                  END DO
               END DO
            END IF
            WRITE (LUPRI,*) 'LABEL:',LABEL
            WRITE (LUPRI,*) 'ISYXIM:',ISYXIM
            WRITE (LUPRI,*)
     &           'FTRACE of matrix generated in CC_EFCKMO:',FTRACE
         END IF

         ! scale with 2 to get contribution to X intermediate
         CALL DSCAL(N2BST(ISYXIM),TWO,WORK(KXIM),1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*)
     &           'CCXIMCON> direct contribution to X intermediate:'
            CALL CC_PRONELAO(WORK(KXIM),ISYXIM)
         END IF       

         IF (LORX .OR. LPDBS) THEN

           KSCR1  = KEND2
           KFOCK1 = KSCR1  + N2BST(ISYXIM)
           KEND3  = KFOCK1 + N2BST(ISYXIM)
           LWRK3  = LWORK  - KEND3

           IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insufficient memory in CCXIMCON.')
           END IF

           IF (LORX) THEN
             CALL CC_RDHFRSP('R1 ',IKAPPA,ISYXIM,WORK(KAPPA))
           ELSE
             CALL DZERO(WORK(KAPPA),2*NALLAI(ISYXIM))
           END IF

           IOPER = IROPER(LABEL,ISYXIM)
           CALL CC_GET_RMAT(WORK(KRMAT),IOPER,1,ISYXIM,
     &                      WORK(KEND3),LWRK3)

           NOKAPPA = .FALSE.
           IREAL   = ISYMAT(IOPER)
           CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),
     &                  WORK(KRMAT),WORK(KAPPA),
     &                  IREAL,ISYXIM,NOKAPPA,
     &                  WORK(KCMO),WORK(KEND3),LWRK3)

           DO ISYM1 = 1, NSYM
              ISYM2 = MULD2H(ISYM1,ISYXIM)
              KOFF1 = KQMATP + IAODIS(ISYM1,ISYM2)
              KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
              CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                    WORK(KOFF1),WORK(KOFF2))
           END DO                         
           CALL DSCAL(N2BST(ISYXIM),-ONE,WORK(KQTRP ),1) 

           CALL CC_MMOMMO('N','N',ONE,WORK(KFOCK0),ISYM0,
     &                    WORK(KQTRP ),ISYXIM,ZERO,WORK(KSCR1),ISYXIM)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &             'CCXIMCON> relax. contrib. 1 to X intermediate:'
              CALL CC_PRONELAO(WORK(KSCR1),ISYXIM)
              FTRACE = ZERO
              DO ISYM = 1, NSYM
                 KOFF1 = KSCR1 + IAODIS(ISYM,ISYM)
                 DO I = 1, NBAS(ISYM)
                   FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                 END DO
              END DO
              WRITE (LUPRI,*) 'trace:',FTRACE
           END IF         

           ! add to result matrix:
           CALL DAXPY(N2BST(ISYXIM),ONE,WORK(KSCR1),1,WORK(KXIM),1)     

           ! read contribution from 'one-index' transformed density
           ! from file
           FREQ  = FRQLRTHF(IKAPPA)
           IFCK1 = I1DXFCK('HAM0    ','R1 ',LABEL,FREQ,ISYXIM)
           IADRF = IADR1DXF(1,IFCK1)
           CALL WOPEN2(LUFCK,FIL1DXFCK,64,0)
           CALL GETWA2(LUFCK,FIL1DXFCK,WORK(KFOCK1),IADRF,N2BST(ISYXIM))
           CALL WCLOSE2(LUFCK,FIL1DXFCK,'KEEP')    

           CALL CC_EFFCKMO(WORK(KFOCK1),ISYXIM,WORK(KCMO),WORK(KOVERLP),
     &                     WORK(KEND3),LWRK3)      

           IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &             'CCXIMCON> relax. contrib. 2 to X intermediate:'
              CALL CC_PRONELAO(WORK(KFOCK1),ISYXIM)
              FTRACE = ZERO
              DO ISYM = 1, NSYM
                 KOFF1 = KFOCK1 + IAODIS(ISYM,ISYM)
                 DO I = 1, NBAS(ISYM)
                   FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                 END DO
              END DO
              WRITE (LUPRI,*) 'trace:',FTRACE
           END IF

           ! add to result matrix:
           CALL DAXPY(N2BST(ISYXIM),ONE,WORK(KFOCK1),1,WORK(KXIM),1)

         END IF

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCXIMCON> final result for X intermediate:'
            CALL CC_PRONELAO(WORK(KXIM),ISYXIM)
            FTRACE = ZERO
            DO ISYM = 1, NSYM
               KOFF1 = KXIM-1 + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I)
               END DO
            END DO
            WRITE (LUPRI,*) 'trace of X intermediate:',FTRACE
         END IF                            

*---------------------------------------------------------------------*
*        calculate all required dot products with this X intermediate:
*---------------------------------------------------------------------*
         IVEC = 1

         DO WHILE (IXDOTS(IVEC,IDXXIM).NE.0 .AND. IVEC.LE.MXVEC)

      
           IFILE = IXDOTS(IVEC,IDXXIM)
           ISYMQ = ILSTSYMRLX(LISTQ,IFILE)

           IF (ISYMQ.NE.ISYXIM) THEN
              WRITE (LUPRI,*) IDXXIM, ISYXIM, LISTQ, IFILE, ISYMQ
              CALL QUIT('symmetry mismatch in CCXIMCON.')
           END IF

           IF (LISTQ(1:3).EQ.'R1 ') THEN
              LORXQ  = .TRUE.
              LABELQ = LRTHFLBL(IFILE)
           ELSE
              CALL QUIT('Unknown LISTQ in CCXIMCON.')
           END IF

           IOPERQ = IROPER(LABELQ,ISYMQ)
           IREAL  = ISYMAT(IOPERQ)

           IF (LORXQ) THEN
             CALL CC_RDHFRSP(LISTQ,IFILE,ISYMQ,WORK(KAPPA))
           ELSE
             CALL DZERO(WORK(KAPPA),2*NALLAI(ISYMQ))
           END IF

           CALL CC_GET_RMAT(WORK(KRMAT),IOPERQ,1,ISYMQ,
     &                      WORK(KEND2),LWRK2)

           NOKAPPA = .FALSE.
           CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),
     &                  WORK(KRMAT),WORK(KAPPA),
     &                  IREAL,ISYMQ,NOKAPPA,
     &                  WORK(KCMO),WORK(KEND2),LWRK2)

           DO ISYM1 = 1, NSYM
              ISYM2 = MULD2H(ISYM1,ISYMQ)
              KOFF1 = KQMATH + IAODIS(ISYM1,ISYM2)
              KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
              CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                    WORK(KOFF1),WORK(KOFF2))
           END DO
           CALL DCOPY(N2BST(ISYMQ),WORK(KQTRP),1,WORK(KQMATH),1)
           CALL DSCAL(N2BST(ISYMQ),HALF,WORK(KQMATH),1)

           DO ISYM1 = 1, NSYM
              ISYM2 = MULD2H(ISYM1,ISYMQ)
              KOFF1 = KQMATP + IAODIS(ISYM1,ISYM2)
              KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
              CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                    WORK(KOFF1),WORK(KOFF2))
           END DO
           CALL DCOPY(N2BST(ISYMQ),WORK(KQTRP),1,WORK(KQMATP),1)
           CALL DSCAL(N2BST(ISYMQ),HALF,WORK(KQMATP),1)

           XCONS(IVEC,IDXXIM) = 
     &                  DDOT(N2BST(ISYMQ),WORK(KQMATH),1,WORK(KXIM),1)+
     &      DBLE(IREAL)*DDOT(N2BST(ISYMQ),WORK(KQMATP),1,WORK(KXIM),1) 


           IVEC = IVEC + 1
         END DO

      END DO

*---------------------------------------------------------------------*
*     print the results:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'CCXIMCON> results for '//
     &        'X intermediate contribs.:'
         IF (MXVEC.NE.0) THEN
            DO IDXXIM = 1, NXIM
               WRITE (LUPRI,*) 'IDXXIM = ',IDXXIM
               IVEC = 1
               DO WHILE(IXDOTS(IVEC,IDXXIM).NE.0 .AND. IVEC.LE.MXVEC)
                  WRITE(LUPRI,'(A,2I5,2X,E18.12)') 'CCXIMCON> ',
     &              IXIM(IDXXIM),IXDOTS(IVEC,IDXXIM),XCONS(IVEC,IDXXIM)
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
