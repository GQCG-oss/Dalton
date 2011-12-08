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
      SUBROUTINE  CCFBINT1(ITRAN,   LISTL,   IDLSTL, LISTR, IDLSTR,
     &                     XAINT,   YAINT,   MAINT, 
     &                     ZETA1,   ZETA2,   ISYCTR,
     &                     TA1AMP,  TA2AMP,  ISYMTA,   
     &                     XLAMDP,  XLAMDH,  ISYLAM,
     &                     XLAMDPQ, XLAMDHQ, ISYHOP,
     &                     BZDENA,  BZDENAB,
     &                     FILPQAMO,LUPQAMO,  IADRPQAMO,IADRPQA,
     &                     FILPQA0, LUPQA0,   IADRPQA0, IADRPQAI0, 
     &                     FILPQA1, LUPQA1,   IADRPQA1, IADRPQAI1,
     &                     LRELAX,  LTWOEL,  WORK,    LWORK)
*---------------------------------------------------------------------*
* Purpose:
*
*     Precalculate some intermediates for FbTa result vector depending
*     on ZETA, T^A and IOPER: 
*     -- the XA and YA intermediates (in memory, allocated outside)
*     -- the effective density BZDENA for the rho^BZA intermediate (disk)
*     -- the effective density BZDENAB for the rho^BZQA intermediate (memory)
*     -- the PA and QA intermediates, both with 4 MO indices and 
*                                     with 1 index in AO (0 and bar)
*
* Note: flags LRELAX & LTWOEL not yet carried through...?
*
*     Sonia Coriani, February 1999 (based on CCETAINT1)
*     Version: 08/10-1999
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LRELAX, LTWOEL
      CHARACTER*(*) LISTL, LISTR
c
      CHARACTER*(*) FILPQAMO, FILPQA0, FILPQA1
c
      INTEGER ITRAN, IDLSTL, ISYCTR, IDLSTR, ISYMTA
      INTEGER ISYLAM, ISYHOP, LWORK
c
      INTEGER LUPQAMO, IADRPQA,   IADRPQAMO(MXCORB,*)
      INTEGER LUPQA0,  IADRPQAI0, IADRPQA0(MXCORB,*)
      INTEGER LUPQA1,  IADRPQAI1, IADRPQA1(MXCORB,*)
      
#if defined (SYS_CRAY)
      REAL XAINT(*), YAINT(*), MAINT(*), ZETA1(*), ZETA2(*)
      REAL BZDENA(*), BZDENAB(*)
      REAL XLAMDP(*), XLAMDH(*), XLAMDPQ(*), XLAMDHQ(*)
      REAL TA1AMP(*), TA2AMP(*), WORK(*)
      REAL DUMMY, ONE, TWO
#else
      DOUBLE PRECISION XAINT(*), YAINT(*), MAINT(*), ZETA1(*), ZETA2(*)
      DOUBLE PRECISION BZDENA(*), BZDENAB(*)
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), XLAMDPQ(*), XLAMDHQ(*)
      DOUBLE PRECISION TA1AMP(*), TA2AMP(*), WORK(*)
      DOUBLE PRECISION DUMMY, ONE, TWO
#endif
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)

      CHARACTER MODEL*(10)
      INTEGER IOPT, ISYINT, ISYINTQ, KCHI, KCHIQ, KEND1, LWRK1, IDEL
      INTEGER IDUMMY, KZ1A, KCHIA
      INTEGER ILLSTOLD
      SAVE    ILLSTOLD
      INTEGER IRLSTOLD
      SAVE    IRLSTOLD

*---------------------------------------------------------------------*
* do some tests and set symmetries:
*---------------------------------------------------------------------*
* if ITRAN=1, make sure that everything will be initialzed:
*
      IF (ITRAN.EQ.1) THEN
         ILLSTOLD = IDLSTL - 1
         IRLSTOLD = IDLSTR - 1
      END IF

* set symmetries for intermediates:

      ISYINT  = MULD2H(ISYCTR,ISYMTA)     !ZxTA
      ISYINTQ = MULD2H(ISYINT,ISYHOP)     !ZxTAxOperB

*---------------------------------------------------------------------*
* if only left vector (IDLSTL) has changed, 
*    read new multipliers Z1+Z2 into memory and square Z2 up:
*---------------------------------------------------------------------*
      IF ((IDLSTL .NE. ILLSTOLD).AND.(IDLSTR.EQ.IRLSTOLD)) THEN

         IF (LWORK .LT. NT2AM(ISYCTR)) THEN
            CALL QUIT('Insufficient memory in CCFBINT1 (a1)')
         END IF 

         IOPT = 3
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,ZETA1,WORK) 

         CALL CC_T2SQ(WORK,ZETA2,ISYCTR)           

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCFBINT1> Only NEW Zeta vector' 
            WRITE (LUPRI,*) 'CCFBINT1> the zeta2 vector:'
            CALL CC_PRP(WORK,WORK,ISYCTR,0,1)
            WRITE (LUPRI,*) 'CCFBINT1> the tA2amp vector:'
            CALL CC_PRP(WORK,TA2AMP,ISYMTA,0,1)
         END IF
      END IF
*---------------------------------------------------------------------*
* if only right vector (IDLSTR) has changed, 
*    read new response amplitudes TA1+TA2 into memory:
*---------------------------------------------------------------------*

      IF ((IDLSTL .EQ. ILLSTOLD).AND.(IDLSTR.NE.IRLSTOLD)) THEN

         IF (LWORK .LT. NT2AM(ISYMTA)) THEN
            CALL QUIT('Insufficient memory in CCFBINT1 (a2)')
         END IF 

         IOPT = 3  
         CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,TA1AMP,TA2AMP)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCFBINT1> Only NEW T^A vector' 
            WRITE (LUPRI,*) 'CCFBINT1> the zeta2 vector:'
            CALL CC_PRP(WORK,WORK,ISYCTR,0,1)
            WRITE (LUPRI,*) 'CCFBINT1> the tA2amp vector:'
            CALL CC_PRP(WORK,TA2AMP,ISYMTA,0,1)
         END IF
*
* scale with 2 the tA_2 part 
*
         IF ( (IOPT.EQ.3) .AND. (.NOT. LISTR(1:2).EQ.'R0') ) THEN
            CALL CCLR_DIASCL(TA2AMP,TWO,ISYMTA)
         END IF
*
      END IF
*---------------------------------------------------------------------*
* if both vectors (IDLSTR,IDLSTL) have changed, 
*                 read new Zeta and response amplitudes into memory:
*---------------------------------------------------------------------*

      IF ((IDLSTL .NE. ILLSTOLD).AND.(IDLSTR.NE.IRLSTOLD)) THEN

c         IF (LWORK .LT. (NT2AM(ISYCTR)+NT2AM(ISYMTA))) THEN
         IF (LWORK .LT. NT2AM(ISYCTR)) THEN
            CALL QUIT('Insufficient memory in CCFBINT1 (a3)')
         END IF 

         IOPT = 3  !(both singles and doubles)
         CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,TA1AMP,TA2AMP)

c         KSCR  = 1
c         KEND1 = KSCR + NT2AM(ISYCTR)
         
         IOPT = 3
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,ZETA1,WORK)

         CALL CC_T2SQ(WORK,ZETA2,ISYCTR)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCFBINT1> Both NEW T^A and Zeta vectors' 
            WRITE (LUPRI,*) 'CCFBINT1> the zeta2 vector:(packed)'
            CALL CC_PRP(WORK,WORK,ISYCTR,0,1)
            WRITE (LUPRI,*) 'CCFBINT1> the tA2amp vector:'
            CALL CC_PRP(WORK(1),TA2AMP,ISYMTA,0,1)
         END IF
*
* scale with 2 the tA_2 part
*
         IF ( (IOPT.EQ.3) .AND. (.NOT. LISTR(1:2).EQ.'R0') ) THEN
            CALL CCLR_DIASCL(TA2AMP,TWO,ISYMTA)
         END IF
*
       END IF
*---------------------------------------------------------------------*
* if either left or right vector changed, calculate 
* new XA, YA and MA intermediates:
*---------------------------------------------------------------------*
      IF (( IDLSTL .NE. ILLSTOLD ).OR.( IDLSTR .NE. IRLSTOLD )) THEN

          CALL CC_XI(XAINT,ZETA2,ISYCTR,TA2AMP,ISYMTA,WORK,LWORK)  

          CALL CC_YI(YAINT,ZETA2,ISYCTR,TA2AMP,ISYMTA,WORK,LWORK)  

          IF (CCSD .AND. LRELAX) THEN
             CALL CC_MI(MAINT,ZETA2,ISYCTR,TA2AMP,ISYMTA,WORK,LWORK)  
          END IF

          IF (LOCDBG) THEN
             WRITE (LUPRI,'(//A)') 'CCFBINT1> XA-intermediate:'
             WRITE (LUPRI,'(5G15.6)') (XAINT(I),I=1,NMATIJ(ISYINT))
             WRITE (LUPRI,'(//A)') 'CCFBINT1> YA-intermediate:'
             WRITE (LUPRI,'(5G15.6)') (YAINT(I),I=1,NMATAB(ISYINT))
             IF (CCSD .AND. LRELAX) THEN
                WRITE (LUPRI,'(//A)') 'CCFBINT1> MA-intermediate:'
                WRITE (LUPRI,'(5G15.6)') (MAINT(I),I=1,N3ORHF(ISYINT))
             END IF
          END IF

      END IF

*---------------------------------------------------------------------*
* calculate Chi^A matrices (calN)=(breve-Zeta1 - XA) intermediate
* transform Zeta^ci to ZetaA^ki with T^A_ck
* (cheap, do always!)
*---------------------------------------------------------------------*
      IF (CCSD) THEN
        KCHIA = 1
        KEND1 = KCHIA + NMATIJ(ISYINT) 
        LWRK1 = LWORK - KEND1
        
        IF (LWRK1 .LT. NT2AM(ISYCTR)) THEN
          CALL QUIT('Insufficient memory in CCFBINT1 (a)')
        END IF 
        CALL CCLT_Z1A(ZETA1,ISYCTR,TA1AMP,ISYMTA, ISYINT,WORK(KCHIA) )
        IF (LOCDBG) THEN
          WRITE (LUPRI,'(//A)') 'CCFBINT1> ChiA-intermediate 1:'
          WRITE (LUPRI,'(5G15.6)') (WORK(KCHIA+I-1),I=1,NMATIJ(ISYINT))
        END IF
        CALL DAXPY(NMATIJ(ISYINT),-ONE,XAINT,1,WORK(KCHIA),1)
        IF (LOCDBG) THEN
          WRITE (LUPRI,'(//A)') 'CCFBINT1> ChiA-intermediate:'
          WRITE (LUPRI,'(5G15.6)') (WORK(KCHIA+I-1),I=1,NMATIJ(ISYINT))
        END IF
      END IF

*---------------------------------------------------------------------*
* calculate the effective BZA density (only for new TA or ZETA)
* D_{k\a}^{ij} (backtr. with Lambda^p) (does not depend on Lambda-bar)
* D_{k\a}^{ij} -> BZDENA (in memory)
*---------------------------------------------------------------------*
      IF (CCSD .AND. LRELAX) THEN
        IF ((IDLSTL .NE. ILLSTOLD).OR.(IDLSTR .NE. IRLSTOLD)) THEN

         IF (LOCDBG) THEN 
           WRITE (LUPRI,*) 'CCFBINT1: Print LambdaP for symmetry 1'
           CALL CC_PRLAM(XLAMDP,XLAMDH,ISYLAM)
           WRITE (LUPRI,*) 'FBINT1: Print ZETA2 for symmetry 1'
           CALL CC_PRSQ(ZETA1,ZETA2,ISYCTR,0,1)
           CALL DZERO(BZDENA,NT2AOIJ(ISYINT))
           CALL DZERO(WORK(KEND1),LWRK1)
         END IF

         CALL CC_BFDENF( ZETA2,  ISYCTR, MAINT, ISYINT,
     *                   XLAMDP, ISYLAM, WORK(KCHIA),  ISYINT,
     *                   TA1AMP, ISYMTA, BZDENA, WORK(KEND1), LWRK1)
         IF (LOCDBG) THEN 
           WRITE (LUPRI,*) 'FBINT1: Print BZDENA for symmetry 1'
           CALL OUTPUT(BZDENA,1,NT1AO(1),1,NMATIJ(1),
     *                          NT1AO(1),NMATIJ(1),1,LUPRI)
           WRITE (LUPRI,*) 'returned 1 from CC_BFDENF... '
           CALL FLSHFO(LUPRI)
         END IF
        ENDIF
      END IF

*---------------------------------------------------------------------*
* calculate the effective BZAB density which does depend on the
* external perturbation for _each_ transformation (ITRAN, as it
* depends on IOPER): (or pass IOPER?)
* 2) Dbar_{k\a}^{ij} -> BZDENAB (send in LambdapQ zero)
*---------------------------------------------------------------------*
      IF (CCSD .AND. LRELAX) THEN
         CALL CC_BFDENF( ZETA2,  ISYCTR, MAINT, ISYINT,
     *                   XLAMDPQ, ISYHOP, WORK(KCHIA),  ISYINT,
     *                   TA1AMP,  ISYMTA, BZDENAB, WORK(KEND1), LWRK1)
         IF (LOCDBG) THEN 
           WRITE (LUPRI,*) 'CCFBINT1: Print BZDENQA for symmetry 1'
           CALL OUTPUT(BZDENAB,1,NT1AO(1),1,NMATIJ(1),
     *               NT1AO(1),NMATIJ(1),1,LUPRI)
           WRITE (LUPRI,*) 'returned 2 from CC_BFDENF... '
           CALL FLSHFO(LUPRI)
         END IF
      ENDIF
*---------------------------------------------------------------------*
* calculate one-index backtransformed PA and QA intermediates used 
* in the F and G terms contribution to the result vector:
*---------------------------------------------------------------------*
      IF (CCSD .AND. LRELAX) THEN

        IF (( IDLSTL .NE. ILLSTOLD ).OR.( IDLSTR .NE. IRLSTOLD )) THEN

          IOPT = 1
          CALL CC_PQIMO(ZETA2,ISYCTR,TA2AMP,ISYMTA,YAINT,IOPT,
     *                  FILPQAMO,LUPQAMO,IADRPQAMO(1,ITRAN),IADRPQA,
     *                  ITRAN,WORK(KEND1),LWRK1)

          CALL CC_PQIAO(FILPQAMO,LUPQAMO,IADRPQAMO(1,ITRAN),ISYINT,
     *                  FILPQA0, LUPQA0, IADRPQA0(1,ITRAN), IADRPQAI0,
     *                  ITRAN, XLAMDH, ISYLAM,WORK(KEND1), LWRK1)


         ELSE
          DO IDEL = 1, NBAST
            !dimensioned bigger that Nvir
            IADRPQAMO(IDEL,ITRAN) = IADRPQAMO(IDEL,ITRAN-1)
            IADRPQA0(IDEL,ITRAN)  = IADRPQA0(IDEL,ITRAN-1)
          END DO
        END IF

C Backtransformation with Lambda-bar always performed
c this is not necessary. I could pass LNEWOP!!!!!!!!!!!!!

        CALL CC_PQIAO(FILPQAMO,LUPQAMO,IADRPQAMO(1,ITRAN),ISYINT,
     *                FILPQA1, LUPQA1, IADRPQA1(1,ITRAN), IADRPQAI1,
     *                ITRAN,XLAMDHQ,ISYHOP,WORK(KEND1), LWRK1)

      END IF
*---------------------------------------------------------------------*
* save present IDLSTL in IDLSTOLD and return:
*---------------------------------------------------------------------*
      ILLSTOLD = IDLSTL
      IRLSTOLD = IDLSTR
  
      RETURN

      END 
*=====================================================================*
