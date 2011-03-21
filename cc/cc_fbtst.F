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
*---------------------------------------------------------------------*
c/* Deck CC_FBTST */
*=====================================================================*
       SUBROUTINE CC_FBTST(WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose: provide some tests for the F{B}T{A} transformation 
*          calculated from a (one- and two-electron) derivative operator
*          B refers to the (two-electron) perturbation operator in F{B}, 
*          A to the other perturbation in T{A}
* Sonia Coriani / Christof Haettig, 1998/99
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <cco1rsp.h>
#include <ccx1rsp.h>
#include <ccfro.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MXFBTRAN
      PARAMETER (MXFBTRAN = 20)
      INTEGER LWORK

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FREQB, FREQA, WRKDLM, DUMMY
      REAL DDOT
      REAL ZERO, ONE, TWO, FOUR, FIVE
      REAL RDUM
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FREQB, FREQA, WRKDLM, DUMMY
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, FIVE
      DOUBLE PRECISION RDUM
#endif 
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0, FIVE = 5.0D0)

      INTEGER IFBTRAN(5,MXFBTRAN), NFBTRAN
      INTEGER ISYM0, ISYHOP, ISYOPA, IOPERB
      INTEGER ITEST, IREAL, IRELAX, IOPT, ILEFT, IRIGHT, IDXR1
      INTEGER IDETA1, IDRHSR1, IDXR1HF
      INTEGER IDUM, IDTA1

      CHARACTER*(3) LISTL, LISTR
      CHARACTER*(8) LABELB, LABELA, FILFBTA
      CHARACTER*(10) MODEL

      LOGICAL LORXB, LORXA, LRSP, LTWO, FD_ON_FMAT, FD_ON_ETA

      INTEGER KT0, KCMOPQ, KCMOHQ
      INTEGER KEND0, KEND1A, LWRK1A
      INTEGER KEND1, KDUM, KEND2, LWRK0, LWRK1, LWRK2, KSCR2
      INTEGER KKAPPA, KRMAT, KRHO1, KRHO2, KOMEGA1, KOMEGA2 
      INTEGER KRHS1, KRHS2, KT2AMP0, KT1AMP0, IORDER

      INTEGER IOPTRES, N2VEC

* external function:
      INTEGER ILSTSYM
      INTEGER IROPER
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IL1ZETA
      INTEGER IRHSR1
      INTEGER IETA1

*---------------------------------------------------------------------*
* set up information for rhs vectors for the different tests:
*---------------------------------------------------------------------*
*  number of simultaneous transformations:
      NFBTRAN = 1
*---------------------------------------------------------------------*
*  test 1: use zeroth-order Hamiltonian as perturbation. 
*          This gives a transformed vector which is 5 * {F T^A}
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c      ITEST = 1
c      LISTL = 'L0'
c      FILFBTA = 'CC_FBMAT'      !output file (debug)
       !The B operator is set uqual to zero-order Hamilton operator
c      LABELB = 'HAM0    '
c      LORXB  = .TRUE.
c      LTWO   = .TRUE.
c      FREQB  = ZERO
c      ISYHOP = 1
       !The A operator and T{A}
c      IOPTRES = 1
c      LISTR = 'R1'
c      ISYOPA = 1
c      LABELA = 'ZDIPLEN '
c      LORXA  = .TRUE.
c      FREQA  = ZERO
c
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 2: test the F^B matrix transformation for a non-relaxed 
*          one-electron perturbation against the old implementation.
*          (F^A) ---> working
*          (does only require that the integrals for the operator
*           are available on file)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c      ITEST = 2
c      LISTL = 'L0'
c      IOPTRES = 1
c      FILFBTA = 'CC_FBMAT'
       !The B operator is a non relaxed one-electron operator
c      LABELB = 'ZDIPLEN '
c      LORXB  = .FALSE.
c      LTWO   = .FALSE.
c      FREQB  = ZERO
c      ISYHOP = 1
       !the A operator and  T{A}
c      ISYOPA = 1
c      LABELA = 'ZDIPLEN '
c      LORXA  = .FALSE.
c      FREQA  = ZERO
c      LISTR = 'R1'
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 3: set differentiated integrals to zero and test only the
*          orbital relaxation & reorthogonalization contributions
*          to the F^A matrix transformed vector.
*          (requires that the CPHF equations for the `reference' 
*          operator, specified here, have been solved and the Kappa
*          vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
      ITEST = 3
      LISTL = 'L0'
      IOPTRES = 1
      FILFBTA = 'CC_FBMAT'      !output file (debug)
      !The B operator is DUMMYOP. It will be defined later (vide infra)
      LABELB = 'ZDIPLEN '
      LORXB  = .TRUE.
      LTWO   = .FALSE.
      FREQB  = ZERO
      ISYHOP = 1
      !the A operator and T{A}
      ISYOPA = 1
      LISTR = 'R1'
      LABELA = 'ZDIPLEN '
      LORXA  = .TRUE.
      FREQA  = ZERO
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 4: test the xi and eta vectors for a `orbital-relaxed'
*          one-electron perturbation.
*          (requires that the CPHF equations for the operator specified 
*           have been solved and the Kappa vector is available on disc)
*---------------------------------------------------------------------*
c     ITEST = 4
c     LABEL  = 'ZDIPLEN '
c     LORX   = .TRUE.
c     LTWO   = .FALSE.
c     FREQ   = ZERO
c     ISYHOP = 1
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 5: test use of London orbitals in FBTA (LAO)
*---------------------------------------------------------------------*
c      ITEST = 5
c      LISTL = 'L0'
c      IOPTRES = 1
c      FILFBTA = 'CC_FBMAT'     !output file (debug)
c
c       !The B operator is a component of magnetic field (x)
c      LABELB = 'dh/dBX  '      
c      LORXB  = .TRUE.
c      LTWO   = .TRUE.
c      FREQB  = ZERO
c      FREQB  = 0.09d+00
c      ISYHOP = 1
       !The A operator and T{A} amplitude
c      ISYOPA = 1
c      LISTR = 'R1'             
c      LABELA = 'YDIPLEN '     
c      LORXA  = .TRUE.
c      FREQA  = ZERO
c      FREQA  = 0.09d+00
*---------------------------------------------------------------------*
* Special test case ITEST = 5
* Force calculation of those vectors we need to be able to do FD on eta
*---------------------------------------------------------------------*
c      ! allow extensions in the vector lists for the next few lines
c      ! to calculate the (X1)^B (eta) vectors ourselves for ITEST 5
c
c      IF (ITEST.EQ.5) THEN
c        LOPROPN = .TRUE.
c        LO1OPN  = .TRUE.
c        LX1OPN  = .TRUE.
c        LR1OPN  = .TRUE.
c        IRELAX = 0
c        IF (LORXB) THEN
c          IRELAX = IR1KAPPA(LABELB,FREQB,ISYHOP)
c        END IF
c
c        LPDBSOP(IROPER(LABELB,ISYHOP)) = LTWO
c
c      ! set the IXETRAN array for one (XI,ETA) pair
c        IXETRAN(1,1) = IROPER(LABELB,ISYHOP)
c        IXETRAN(2,1) = 0
c        IXETRAN(3,1) = IRHSR1(LABELB,LORXB,FREQB,ISYHOP)
c        IXETRAN(4,1) = IETA1(LABELB,LORXB,FREQB,ISYHOP)
c        IXETRAN(5,1) = IRELAX
c      ! disallow again extension in the vector lists
c        LOPROPN = .FALSE.
c        LO1OPN  = .FALSE.
c        LX1OPN  = .FALSE.
c
c
*---------------------------------------------------------------------*
* call CC_XIETA to calculate the Xi and Eta vectors:
*---------------------------------------------------------------------*
c        LISTL  = 'L0 '
c        FILXI  = 'O1 '
c        FILETA = 'X1 '
c        IOPTRES = 3
c
c        KEND0 = 1
c        KEND1 = KEND0 + 1
c        LWRK1 = LWORK - KEND1
c        IF (LWRK1 .LT. 0) THEN
c           CALL QUIT('Insufficient work space in CC_XETST. (1)')
c        END IF
c
c        IDUM = 0
c        WORK(KEND1-1) = WORK(IDUM)
c        WRITE (LUPRI,*) 'WORK(0) before entering CC_XIETA:',WORK(KEND1-1)
c
c        CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL,
c     &                FILXI,  IXDOTS, XCONS,
c     &                FILETA, IEDOTS, ECONS,
c     &                .FALSE.,0,  WORK(KEND1), LWRK1 )
c
c        WRITE (LUPRI,*) 'returned from CC_XIETA... WORK(0) :',WORK(KEND1)
c        CALL FLSHFO(LUPRI)
c
c      END IF
*--------------------------------------------------------------------*
* End of special test case
*--------------------------------------------------------------------*

c     get the index of the kappa vector (kappa{B}) which is saved on the 
c     same file on which the corresponding T amplitude responses (T{B})
c     would be saved a value of 0 indicates `unrelaxed'

      IRELAX = 0
      IF (LORXB) THEN
        IRELAX = IR1KAPPA(LABELB,FREQB,ISYHOP)      !it could be frequency dep
      END IF                                             

c     set logical for perturbation-dependent basis sets (for B); this will
c     overwrite the default already saved on the common block

      LPDBSOP(IROPER(LABELB,ISYHOP)) = LTWO               

c     set the IFBTRAN array for one F{B} T{A} result vector

      IFBTRAN(1,1) = IROPER(LABELB,ISYHOP)               ! index B operator
      IFBTRAN(2,1) = 0                                   ! index left vector       (L0)
      IFBTRAN(3,1) = IR1TAMP(LABELA,LORXA,FREQA,ISYOPA)  ! index right vector T{A} (R1)
      IFBTRAN(4,1) = 0                                   ! index result vector (return in memory)
      IFBTRAN(5,1) = IRELAX                              ! index kappa vector

c     choose finite difference path:
      FD_ON_FMAT = .FALSE.   ! f.d. on usual F matrix
      FD_ON_ETA  = .TRUE.    ! f.d. on Eta vector
 
      IF (ITEST.EQ.3) THEN
 
         ! for test 3, replace now the operator labels on the 
         ! common blocks by 'DUMMYOP ', which should be interpreted 
         ! inside of CC_XIETA/CC_FBTA as a zero one-electr. operator,
         ! but associated with a orb.-relax. (kappa) vector.

         ! index for a R1 type vect for LABELB operator (B)
         IDXR1   = IR1TAMP(LABELB,LORXB,FREQB,ISYHOP)
         ! index for correponding CPHF vector
         IDXR1HF = IR1KAPPA(LABELB,FREQB,ISYHOP)
         KT0   = 1
         KEND1 = KT0   + 2*NALLAI(ISYHOP)
         LWRK1 = LWORK - KEND1
         CALL CC_RDHFRSP('R1 ',IDXR1HF,ISYHOP,WORK(KT0))

         IDRHSR1 = IRHSR1(LABELB,LORXB,FREQB,ISYHOP)          !index for O1, RHS of R1 equations (TB), dvs csi^B (O1)
         IDETA1  = IETA1(LABELB,LORXB,FREQB,ISYHOP)           !index for X1, the Eta^B vector (X1)
         LABELB  = 'DUMMYOP '
         LBLOPR(IFBTRAN(1,1)) = LABELB                        ! operator common block
         LRTLBL(IDXR1)        = LABELB                        ! common block for R1{B}/kappa{B}
         LBLO1(IDRHSR1) = LABELB                              !operator label of IDRHSR1 (O1)
         LBLX1(IDETA1)  = LABELB                              !operator label of IDETA1 (X1)
*
* substitute DUMMY to ZDIPLEN on file  for Kappa (use IOPT = 4)
*
         IOPT = 4
         CALL CC_WRRSP('R1 ',IDXR1HF,ISYHOP,IOPT,MODEL,WORK(KT0),
     &                 DUMMY,DUMMY,WORK(KEND1),LWRK1)
 
         WRITE (LUPRI,*) 'CC_FBTST: orbital relaxation vector:'
         CALL OUTPUT(WORK(KT0),1,2*NALLAI(ISYHOP),1,1,
     &                           2*NALLAI(ISYHOP),1,1,LUPRI)

      END IF
*---------------------------------------------------------------------*
* Print general infos 
*---------------------------------------------------------------------*
      WRITE (LUPRI,*) 'Test case nr  ITEST  =',ITEST
      WRITE (LUPRI,*) 'B operator   LABELB  =',LABELB
      WRITE (LUPRI,*) '2-el operator? LTWO  =',LTWO
      WRITE (LUPRI,*) 'relaxed oper? LORXB  =',LORXB
      WRITE (LUPRI,*) 'of frequency? FREQB  =',FREQB
      WRITE (LUPRI,*) 'IROPER =',IFBTRAN(1,1)
      WRITE (LUPRI,*) 'ILEFT  =',IFBTRAN(2,1)
      WRITE (LUPRI,*) 'IRIGHT =',IFBTRAN(3,1)
      WRITE (LUPRI,*) 
     &    'IRELAX =',IFBTRAN(5,1),' LRTLBL(RELAX):', LRTLBL(IRELAX)
      WRITE (LUPRI,*) 'A op. in T{A} LABELA =',LABELA, 
     &                             ' LRTLBL(R1):', LRTLBL(IFBTRAN(3,1))
      WRITE (LUPRI,*) 'relaxed?      LORXA  =',LORXA
      WRITE (LUPRI,*) 'of frequency? FREQA  =',FREQA
 
      N2VEC = 1
      IF (CCS) N2VEC = 0

*---------------------------------------------------------------------*
* call CC_FBTA  to calculate F{B}T{A} transformed vector:
*---------------------------------------------------------------------*
      
      KEND0 = 1
      LWRK0 = LWORK
      CALL CC_FBTA(IFBTRAN, NFBTRAN, IOPTRES, LISTL, LISTR,
     &                      FILFBTA, IDUM, RDUM,0,
     &                      WORK(KEND0), LWRK0)

c      CALL CC_FBTA(IFBTRAN, NFBTRAN, IOPTRES, LISTL, LISTR,
c     &                      FILFBTA, IFBDOTS, FBCON, MXVEC,
c     &                      WORK(KEND0), LWRK0)

*---------------------------------------------------------------------*
*     calculate the reference vector and compare 
*     with the results from the CC_FBTA routine:
*---------------------------------------------------------------------*

* ------------------------------------------------------------------- *
*  ITEST = 1     Test against F matrix
* ------------------------------------------------------------------- *
      IF (ITEST.EQ.1) THEN
        ILEFT  = IFBTRAN(2,1)
        IRIGHT = IFBTRAN(3,1)
        CALL AROUND('ITEST =1, 5 times the Fmatrix transformation:')
        CALL CCFTRANSC(LISTL, ILEFT, LISTR, IRIGHT, WORK(KEND0),
c        CALL CC_FTRAN(LISTL, ILEFT, LISTR, IRIGHT, WORK(KEND0),
     &                                                   LWRK0)
      END IF

      IF (ITEST.EQ.1) RETURN

* -------------------------------------------------------------------- *
* ITEST > 1      Test against others
* -------------------------------------------------------------------- *
      KKAPPA  = 1
      KRMAT   = KKAPPA  + 2*NALLAI(ISYHOP)
      KCMOPQ  = KRMAT   + N2BST(ISYHOP)
      KCMOHQ  = KCMOPQ  + NGLMDT(ISYHOP)
      KT1AMP0 = KCMOHQ  + NGLMDT(ISYHOP)
      KOMEGA1 = KT1AMP0 + NT1AMX
      KOMEGA2 = KOMEGA1 + NT1AM(ISYHOP)
      KRHS1   = KOMEGA2 + NT2AM(ISYHOP)
      KRHS2   = KRHS1   + NT1AM(ISYHOP)
      KEND1   = KRHS2   + NT2AM(ISYHOP)
      LWRK1   = LWORK   - KEND1

      KT2AMP0 = KOMEGA2 + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))
      KSCR2   = KT2AMP0 + NT2AMX
      KEND1A  = KSCR2   + NT2AMX + NT1AMX
      LWRK1A  = LWORK   - KEND1A

      IF (LWRK1.LT.0 .OR. LWRK1A.LT.0) THEN
         CALL QUIT('Insufficient memory in CC_FBTST.')
      END IF
 
c
      IF (LORXB) THEN
        IF (LABELB.EQ.'HAM0    ') THEN
          CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYHOP))
        ELSE
      CALL FLSHFO(LUPRI)
          CALL CC_RDHFRSP('R1 ',IRELAX,ISYHOP,WORK(KKAPPA))
      CALL FLSHFO(LUPRI)
        END IF
      ELSE
        CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYHOP))
      END IF

*     ------------------------------
*     get orbital connection matrix:
*     ------------------------------
       IF (LTWO) THEN
        IOPERB = IROPER(LABELB,ISYHOP)
        IORDER= 1
        CALL CC_GET_RMAT(WORK(KRMAT),IOPERB,IORDER,
     &                   ISYHOP,WORK(KEND1),LWRK1)
       END IF 

*     ------------------------------------
*     construct the derivative CMO vector:
*     ------------------------------------
      IF (LORXB.OR.LTWO) THEN
         IREAL = ISYMAT(IROPER(LABELB,ISYHOP))
         IOPT = 0
         CALL CC_LAMBDAQ(DUMMY,DUMMY,WORK(KCMOPQ),WORK(KCMOHQ),ISYHOP,
     &                   DUMMY,WORK(KKAPPA),WORK(KRMAT),IREAL,IOPT,
     &                   WORK(KEND1),LWRK1)
      END IF

      IF (ITEST.EQ.2) THEN
*         ------------------------------------------------------------
*         call the usual FA matrix routine to calculate the reference
*         result for the FA transformed vector:
*         ------------------------------------------------------------
c          ILEFT  = IFBTRAN(2,1)
c          IRIGHT = IFBTRAN(3,1)
c          CALL CCLR_FA(LABELB, ISYHOP,  ! inp: label/symmetry A
c     &                 LISTR,  IRIGHT,  ! inp: B resp. amplit.
c     &                 LISTL,  ILEFT,   ! inp: C resp. zeta vec.
c     &                 WORK(KEND1),   LWRK1   )! work space
c
      END IF

      IF (LORXB.OR.LTWO) THEN
*
*        ------------------------------------------------------------
*        evaluate orbital relaxation and reorthogonalization contrib.
*        to the F matrix transformed vector by 
*           a) finite difference on the usual F matrix (mht C)
*           b) finite difference on the Eta{A} vector (CC_XIETA) (mht T)
*
*        (does only work for totally symmetric perturbations... 
*        actually only tested without symmetry...)
*        ------------------------------------------------------------
         IF (ISYHOP.NE.1) THEN
            WRITE (LUPRI,*) 'finite difference test of CC_FBTA not '
            WRITE (LUPRI,*) 
     &           'available for non-totally sym. perturbations.'
            CALL QUIT(
     &           'CC_FDFBMAT only implemented for tot. sym. perturb.')
         END IF

         IF (FD_ON_FMAT) THEN
            ILEFT  = IFBTRAN(2,1)
            IRIGHT = IFBTRAN(3,1)

            IF (IREAL.EQ.-1) THEN
               WRITE (LUPRI,*) 
     &               'CC_FDFBMAT not applicable for IREAL = -1.'
            END IF
            WRITE(LUPRI,*)  'NOW CALL CC_FDFBMAT'
            CALL CC_FDFBMAT(WORK(KCMOHQ),LISTL,ILEFT,LISTR,IRIGHT,
     &                      WORK(KOMEGA1),WORK(KEND1),LWRK1)
         END IF

         IF (FD_ON_ETA) THEN
           IRIGHT = IFBTRAN(3,1)
           CALL CC_FDFBMAT2(LISTR,IRIGHT,WORK(KOMEGA1),WORK(KOMEGA2),
     &                      LABELB,IRELAX,WORK(KEND1),LWRK1)
           CALL AROUND('fin. diff. orbital relaxation F^B T^A vector:')
           CALL CC_PRP(WORK(KOMEGA1),WORK(KOMEGA2),ISYHOP,1,N2VEC)
         END IF

         CALL CC_PRP(WORK(KOMEGA1),WORK(KOMEGA2),ISYHOP,1,N2VEC)

      ELSE
         WRITE (LUPRI,*) 'no orb. relax. or reorth. contribution '//
     &        'to F^B T^A.'
         CALL DZERO(WORK(KOMEGA1),NT1AM(ISYHOP))
         CALL DZERO(WORK(KOMEGA2),NT2AM(ISYHOP))
      END IF

C     ------------------------------------------------------
C     calculate the contribution from one-electron h^1:
C     (not needed, if we do finite difference on Eta vector)
C     ------------------------------------------------------
      IF (FD_ON_FMAT .AND.
     & (LABELB(1:8).NE.'DUMMYOP '.AND. LABELB(1:8).NE.'HAM0    ')) THEN

         CALL QUIT('one-electron contr. not yet available in CC_FBTA.')

         CALL AROUND('one-electron h^1 contribution to Xi vector:')
         CALL CC_PRP(WORK(KRHS1),WORK(KRHS2),ISYHOP,1,N2VEC)
         CALL DAXPY(NT1AM(ISYHOP),ONE,WORK(KRHS1),1,WORK(KOMEGA1),1)
         CALL DAXPY(NT2AM(ISYHOP),ONE,WORK(KRHS2),1,WORK(KOMEGA2),1)

         CALL AROUND('num. orb.-relax. + analyt. one-electr. F^A T^B:')
         CALL CC_PRP(WORK(KOMEGA1),WORK(KOMEGA2),ISYHOP,1,N2VEC)

      END IF

C     -------------------------------------------------
C     read the result vector from CC_FBTA and compare:
C     -------------------------------------------------
c      DO ITRAN = 1, NFBTRAN
c        IOPT = 3
c        IVEC = IFBTRAN(3,ITRAN)
c        PRINT *,'CC_FBTST: IVEC = ',IVEC
c        CALL CC_RDRSP('O1 ',IVEC,ISYHOP,IOPT,MODEL,
c     &                WORK(KRHS1),WORK(KRHS2))
c
c        CALL AROUND('analytical F^B T^A vector:')
c        CALL CC_PRP(WORK(KRHS1),WORK(KRHS2),ISYHOP,1,N2VEC)
c
c        CALL DAXPY(NT1AM(ISYHOP),-1.0d0,WORK(KOMEGA1),1,WORK(KRHS1),1)
c        CALL DAXPY(NT2AM(ISYHOP),-1.0d0,WORK(KOMEGA2),1,WORK(KRHS2),1)
c
c        WRITE(LUPRI,*) 'Norm of difference between analytical F^B T^A '
c     >             // 'vector and the numerical result:'
c        WRITE (LUPRI,*) 'singles excitation part:',
c     >     DSQRT(DDOT(NT1AM(ISYHOP),WORK(KRHS1),1,WORK(KRHS1),1))
c        IF (.NOT.CCS) WRITE (LUPRI,*) 'double excitation part: ',
c     >     DSQRT(DDOT(NT2AM(ISYHOP),WORK(KRHS2),1,WORK(KRHS2),1))
c
c        WRITE (LUPRI,*) 'difference vector:'
c        CALL CC_PRP(WORK(KRHS1),WORK(KRHS2),ISYHOP,1,N2VEC)
c      END DO

      RETURN
      END 
**=====================================================================*
