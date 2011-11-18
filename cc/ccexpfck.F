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
      SUBROUTINE CCEFFFOCK1(I1DXORD,AODENCC,AODENHF,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate the one-electron contributions to 
*                1) CC effective Fock matrices
*                2) CC expection values
*                3) CPHF RHS vectors
*                4) 1-index transformed effective CC Fock matrices
*              
*
*     Christof Haettig, Februar 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccfro.h>
#include <ccroper.h>
#include <ccropr2.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccr1rsp.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
       
      INTEGER LUFCK, ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER I1DXORD, LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL AODENCC(*), AODENHF(*)
      REAL FTRACE, XNORM
      REAL ZERO, HALF, ONE, TWO
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION AODENCC(*), AODENHF(*)
      DOUBLE PRECISION FTRACE, XNORM
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)

      CHARACTER LABEL*8, LSTRLX*3, FILFCK*8, MODEL*10
      LOGICAL LEXPEC, LFOCK
      INTEGER NFOCKLBL, IDLST, ISYMOPR, ISYFCK, ISYRLX, IDXRLX, ISYRMAT
      INTEGER KKAPPA, KRMAT, KAODSY, K1INT, KEND0, LWRK0, KEND1, LWRK1
      INTEGER KQAOS, KBAODSY, KFCKEF, KOVERLP
      INTEGER IFIELD, IOPT, IOPRMAT, IREAL, LEN, ISYM, KOFF1
      INTEGER ISTRTFCK, IOPER, IMATRIX, IERR

* external functions:
      INTEGER IROPER
      INTEGER ILSTSYM, ILSTSYMRLX
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

      CALL QENTER('EFFFOCK1')
*---------------------------------------------------------------------*
* if LOCDBG set, print some debug output:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'entered CCEFFFOCK1... I1DXORD = ',I1DXORD
         XNORM = DDOT(N2BST(ISYM0),AODENCC,1,AODENCC,1)
         WRITE (LUPRI,*) 'CCEFFFOCK1> the AO CC density matrix:',XNORM
         CALL CC_PRONELAO(AODENCC,1)
         WRITE (LUPRI,*) 'CCEFFFOCK1> the AO SCF density matrix:'
         CALL CC_PRONELAO(AODENHF,1)
      END IF

*---------------------------------------------------------------------*
* some initializations:
*---------------------------------------------------------------------*
      ! set file name for file with effective Fock matrices and
      ! get the length of the loop of Fock matrices / exp. values
      IF      (I1DXORD.EQ.0) THEN
        FILFCK   = FILFCKEFF
        NFOCKLBL = NEXPFCKLBL
      ELSE IF (I1DXORD.EQ.1) THEN
        FILFCK   = FIL1DXFCK
        NFOCKLBL = N1DXFLBL
      ELSE
        CALL QUIT('Illegal value of I1DXORD in CCEFFFOCK1.')
      END IF                           

      ! open file for effective Fock matrices
      LUFCK = -1
      CALL WOPEN2(LUFCK,FILFCK,64,0)        

      ! initialize start address for eff. Fock matrices on file
      ISTRTFCK = 1

*---------------------------------------------------------------------*
* symmetrize the CC density matrix and read overlap matrix:
*---------------------------------------------------------------------*
      KAODSY  = 1
      KEND0   = KAODSY  + N2BST(ISYM0)
      IF (I1DXORD.EQ.1) THEN
        KOVERLP = KEND0
        KEND0   = KOVERLP + N2BST(ISYM0)
      END IF
      LWRK0   = LWORK  - KEND0

      IF (LWRK0 .LT. 0) THEN 
        WRITE (LUPRI,*) 'Insufficient memory in CCEFFFOCK1. (0)' 
        CALL QUIT('Insufficient memory in CCEFFFOCK1. (0)' )
      END IF

      ! symmetrize the density matrix:  DSY(pq) = D(pq) + D(qp)
      DO ISYM = 1, NSYM
         KOFF1 = 1 + IAODIS(ISYM,ISYM)
         CALL TRSREC(NBAS(ISYM),NBAS(ISYM),
     &               AODENCC(KOFF1),WORK(KAODSY-1+KOFF1))
      END DO
      CALL DAXPY(N2BST(ISYM0),ONE,AODENCC,1,WORK(KAODSY),1)


      ! if we are going to calculate one-index transformed effective
      ! Fock matrices, we read here the overlap matrix: 
      IF (I1DXORD.EQ.1) THEN
         IF (LWRK0.LT.NBAST) 
     *      CALL QUIT('Insufficient memory in CCEFFFOCK1.')
         CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND0),NBAST)
         CALL CCSD_SYMSQ(WORK(KEND0),ISYM0,WORK(KOVERLP))
      END IF             

*---------------------------------------------------------------------*
* start loop over complete list of Fock matrices / expectation values
* and set labels, symmetries, indeces etc.
*---------------------------------------------------------------------*

      DO IDLST = 1, NFOCKLBL

        IF (I1DXORD.EQ.0) THEN
           LABEL   = LBLEXPFCK(IDLST)  ! operator label
           ISYMOPR = ISYEXPFCK(IDLST)  ! operator symmetry
           LEXPEC  = LEXPFCK(1,IDLST)  ! flag for expectation value
           LFOCK   = LEXPFCK(2,IDLST)  ! flag for effective Fock mat.
           ISYFCK  = ISYMOPR           ! symmetry of eff. Fock mat.

           ! initialize EXPVALUE & IADRFCK
           EXPVALUE(1,IDLST) = ZERO
           EXPVALUE(2,IDLST) = ZERO
           EXPVALUE(3,IDLST) = ZERO
           EXPVALUE(4,IDLST) = ZERO
           IADRFCK(1,IDLST)  = 0
           IADRFCK(2,IDLST)  = 0                   
        ELSE IF (I1DXORD.EQ.1) THEN
           LABEL  = LBL1DXFCK(IDLST)        ! operator label
           IOPER  = IROPER(LABEL,ISYMOPR)   ! operator index/symmetry
           LEXPEC = .FALSE.                 ! flag for expect. value
           LFOCK  = .TRUE.                  ! flag for eff. Fock mat.
           LSTRLX = LST1DXFCK(IDLST)        ! list type for relx. vec.
           IDXRLX = IRELAX1DX(IDLST)        ! list index for relx. vec.
           ISYRLX = ILSTSYMRLX(LSTRLX,IDXRLX)  ! symmetry of relx. vec.
           ISYFCK = MULD2H(ISYMOPR,ISYRLX)  ! symmetry of eff Fock mat.

           ! initialize IADR1DXF
           IADR1DXF(1,IDLST) = 0
           IADR1DXF(2,IDLST) = 0
           IF (LOCDBG) THEN
             WRITE(LUPRI,*)'CCEFFFOCK1> I1DXORD=1'
             WRITE(LUPRI,*)'LABEL,IOPER:',LABEL,IOPER
             WRITE(LUPRI,*)'LSTRLX,IDXRLX,ISYRLX:',LSTRLX,IDXRLX,ISYRLX
           END IF
        ELSE
           CALL QUIT('Illegal value of I1DXORD in CCEFFFOCK1.')
        END IF                                            

      IF (LEXPEC.OR.LFOCK) THEN

*---------------------------------------------------------------------*
* allocate memory for some intermediates:
*---------------------------------------------------------------------*
      K1INT  = KEND0
      KFCKEF = K1INT  + N2BST(ISYMOPR)
      KEND1  = KFCKEF + N2BST(ISYFCK)
 
      IF (I1DXORD.GE.1) THEN
        KKAPPA  = KEND1
        KRMAT   = KKAPPA  + 2*NALLAI(ISYRLX)
        KQAOS   = KRMAT   + N2BST(ISYRLX)
        KBAODSY = KQAOS   + N2BST(ISYRLX)
        KEND1   = KBAODSY + N2BST(ISYRLX)
      END IF

      LWRK1 = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN 
        WRITE (LUPRI,*) 'Insufficient memory in CCEFFFOCK1. (1)' 
        CALL QUIT('Insufficient memory in CCEFFFOCK1. (1)' )
      END IF

*---------------------------------------------------------------------*
* get the one-electron integral matrix:
*---------------------------------------------------------------------*
      IF (LABEL(1:8).EQ.'HAM0    ') THEN

         ! hamiltonian used in the energy calculation:
         ! zeroth-order hamiltonian + finite fields
         CALL CCRHS_ONEAO(WORK(K1INT),WORK(KEND1),LWRK1)
         DO IFIELD = 1, NFIELD
            CALL CC_ONEP(WORK(K1INT),WORK(KEND1),LWRK1,
     &                   EFIELD(IFIELD),ISYMOPR,LFIELD(IFIELD))
         END DO

      ELSE IF (LABEL(1:8).EQ.'DUMMYOP ') THEN

         ! dummy zero operator (for tests)
         CALL DZERO(WORK(K1INT),N2BST(ISYMOPR))
                                                        
      ELSE

         ! the standard case: retrieve integrals using CCPRPAO
         CALL CCPRPAO(LABEL,.TRUE.,WORK(K1INT),ISYM,IMATRIX,IERR,
     &                     WORK(KEND1),LWRK1)
         IF (IERR.LT.0) THEN
           CALL DZERO(WORK(K1INT),N2BST(ISYMOPR))
         END IF
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'norm^2(K1INT):',
     &       DDOT(N2BST(ISYMOPR),WORK(K1INT),1,WORK(K1INT),1)
c       WRITE (LUPRI,*)LABEL,' 1el integral matrix:'
c       CALL CC_PRONELAO(WORK(K1INT),ISYMOPR)
      END IF                                 

*----------------------------------------------------------------------*
* calculate the one-electron expectation values 
*----------------------------------------------------------------------*
      IF (LEXPEC .AND. ISYMOPR.EQ.1) THEN

         ! with correlated density
         EXPVALUE(1,IDLST) = HALF * 
     &                DDOT(N2BST(ISYM0),WORK(K1INT),1,WORK(KAODSY),1)

         ! with SCF density
         EXPVALUE(3,IDLST) = DDOT(N2BST(ISYM0),WORK(K1INT),1,AODENHF,1)

      END IF                                                

*----------------------------------------------------------------------*
* get the orbital relaxation vector, the connection matrix and the
* one-index transformed ao density matrix:
*----------------------------------------------------------------------*
      IF (I1DXORD.GE.1) THEN

         ! retrieve orbital relaxation vector
         IF (IDXRLX.GE.1) THEN
           CALL CC_RDHFRSP(LSTRLX,IDXRLX,ISYRLX,WORK(KKAPPA))
         ELSE
           CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYRLX))
         END IF

         ! retrieve orbital connection matrix
         IF (LSTRLX(1:2).EQ.'R1') THEN
            IOPRMAT = IROPER(LRTHFLBL(IDXRLX),ISYRMAT)
         ELSE
            WRITE (LUPRI,*) 'IOPRMAT problem in CCEFFFOCK1.'
            CALL QUIT('IOPRMAT problem in CCEFFFOCK1.')
         END IF
         CALL CC_GET_RMAT(WORK(KRMAT),IOPRMAT,1,ISYRLX,
     &                    WORK(KEND1),LWRK1)

         IREAL = ISYMAT(IOPRMAT)
         CALL CC_QAOS(WORK(KQAOS),WORK(KRMAT),WORK(KKAPPA),
     &                IREAL,ISYRLX,WORK(KOVERLP),WORK(KEND1),LWRK1)

         ! transform the second index of the density matrix
         ! to the barred (one-index transformed) ao basis
         CALL CC_MAOMAO('N','T',ONE, WORK(KAODSY), ISYM0,
     &                               WORK(KQAOS),  ISYRLX,
     &                         ZERO, WORK(KBAODSY),ISYRLX)      

          
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'IDXRLX, LSTRLX:',IDXRLX, LSTRLX
           WRITE (LUPRI,*) 'ISYRLX, LRTLBL:',ISYRLX, LRTLBL(IDXRLX)
           WRITE (LUPRI,*) 'norm^2(kappa):',
     &       DDOT(2*NALLAI(ISYRLX),WORK(KKAPPA),1,WORK(KKAPPA),1)
           WRITE (LUPRI,*) 'norm^2(krmat):',
     &       DDOT(N2BST(ISYRLX),WORK(KRMAT),1,WORK(KRMAT),1)
           WRITE (LUPRI,*) 'norm^2(kqaos):',
     &       DDOT(N2BST(ISYRLX),WORK(KQAOS),1,WORK(KQAOS),1)
           WRITE (LUPRI,*) 'norm^2(kaodsy):',
     &       DDOT(N2BST(ISYM0),WORK(KAODSY),1,WORK(KAODSY),1)
           WRITE (LUPRI,*) 'norm^2(kbaodsy):',
     &       DDOT(N2BST(ISYRLX),WORK(KBAODSY),1,WORK(KBAODSY),1)
           CALL FLSHFO(LUPRI)
         END IF
      END IF                        

*---------------------------------------------------------------------*
* calculate the one-electron contributions to the eff. Fock matrices,
* and write them to file:
*---------------------------------------------------------------------*
      IF (LFOCK) THEN

         IF      (I1DXORD.EQ.0) THEN

            ! effective Fock matrix calculated from the usual
            ! CC density matrix
            CALL CC_MAOMAO('N','N',HALF,WORK(KAODSY), ISYM0,
     &                                  WORK(K1INT),  ISYMOPR,
     &                             ZERO,WORK(KFCKEF), ISYFCK)

         ELSE IF (I1DXORD.EQ.1) THEN   

            ! effective Fock matrix calculated from the one-index
            ! transformed density matrix
            CALL CC_MAOMAO('N','N',HALF,WORK(KBAODSY),ISYRLX,
     &                                  WORK(K1INT),  ISYMOPR,
     &                             ZERO,WORK(KFCKEF), ISYFCK)

         ELSE
            CALL QUIT('Illegal value for I1DXORD in CCEFFFOCK1.')
         END IF                                                  


         LEN = N2BST(ISYFCK)
         CALL PUTWA2(LUFCK,FILFCK,WORK(KFCKEF),ISTRTFCK,LEN)
         IF (I1DXORD.EQ.0) IADRFCK(1,IDLST)  = ISTRTFCK
         IF (I1DXORD.EQ.1) IADR1DXF(1,IDLST) = ISTRTFCK
         ISTRTFCK = ISTRTFCK + LEN
 
         IF (I1DXORD.EQ.0) THEN
            CALL PUTWA2(LUFCK,FILFCK,WORK(K1INT),ISTRTFCK,LEN)
            IADRFCK(2,IDLST)  = ISTRTFCK
            ISTRTFCK = ISTRTFCK + LEN
         END IF                                                   


         IF (LOCDBG) THEN
           FTRACE = ZERO
           IF (ISYFCK.EQ.1) THEN
             DO ISYM = 1, NSYM
               KOFF1 = KFCKEF + IAODIS(ISYM,ISYM)
               DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
               END DO
             END DO
           END IF                                   

           WRITE (LUPRI,*) 
     *           'one-electron contribution to effective Fock:'
           WRITE (LUPRI,*) 'LABEL :',LABEL
           WRITE (LUPRI,*) 'ISYFCK:',ISYFCK
           WRITE (LUPRI,*) 'IDLST :',IDLST
           WRITE (LUPRI,*) 'TRACE :',FTRACE
           WRITE (LUPRI,*) 'NORM  :',
     &       DDOT(N2BST(ISYFCK),WORK(KFCKEF),1,WORK(KFCKEF),1)
           IF (I1DXORD.EQ.0) THEN
             WRITE (LUPRI,*) 'IADRFCK(1,IDLST) = ',IADRFCK(1,IDLST)
             WRITE (LUPRI,*) 'IADRFCK(2,IDLST) = ',IADRFCK(2,IDLST)
c            WRITE (LUPRI,*) 'HF Fock matrix (one-electron contrib.):'
c            CALL CC_PRONELAO(WORK(K1INT),ISYFCK)
           ELSE IF (I1DXORD.EQ.1) THEN
             WRITE (LUPRI,*) 'IADR1DXF(1,IDLST) = ',IADR1DXF(1,IDLST)
           END IF
           IF (LEXPEC .AND. ISYMOPR.EQ.1)
     &        WRITE (LUPRI,*) 
     *           'EXPVAL:',EXPVALUE(1,IDLST),EXPVALUE(3,IDLST)
         END IF

      END IF ! (LFOCK) 

*---------------------------------------------------------------------*
* end loop over list of Fock matrices / expectation values and return
*---------------------------------------------------------------------*

      END IF ! (LEXPEC.OR.LFOCK)
      END DO ! IDLST

      ! close file for effective Fock matrices
      CALL WCLOSE2(LUFCK,FILFCK,'KEEP')

      CALL QEXIT('EFFFOCK1')
      RETURN
      END
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCEFFFOCK2(I1DXORD,LABEL,ISCOOR,ICOOR,ICORSY,
     &                      MXCOMP,ISYRLX,
     &                      G,ISYMG,D,ISYMD,LUFCK,FILFCK,ISYFCK,IDLST,
     &                      LFOCK, LEXPEC, ISYDEN, ISYMOPR,
     &                      XLAMDP,XLAMDH,QAOS,D1HFAO,DNS1D,XKAPAO,
     &                      XKAPAOB1,D1HFAOB1,XKAPAOB2,D1HFAOB2,
     &                      D2GAI,  D2GAB,  D2GIJ,  D2GIA,
     &                      D2GAIB1,D2GABB1,D2GIJB1,D2GIAB1,
     &                      D2GAIB2,D2GABB2,D2GIJB2,D2GIAB2,
     &                      DIRINT,LDERINT,IRECNR,NRECS,MXBUF,NFILES,
     &                      WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate the two-electron contributions to 
*                1) CC effective Fock matrices
*                2) CC expection values
*                3) CPHF RHS vectors
*                4) 1-index transformed effective CC Fock matrices
*              
*
*     Christof Haettig, Februar 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <ccroper.h>
#include <ccropr2.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccr1rsp.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
       
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER LABEL*8, FILFCK*8
      LOGICAL LEXPEC, LFOCK, DIRINT, LDERINT(8,3)
      INTEGER I1DXORD, LWORK, LUFCK, ISYFCK, ISYMOPR, IDLST, NFILES
      INTEGER ISYMG, ISYMD, ISCOOR,ICOOR,ICORSY, MXCOMP, ISYDEN, ISYRLX
      INTEGER IRECNR(*), NRECS(*), MXBUF

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XLAMDP(*), XLAMDH(*)
      REAL D1HFAO(*), QAOS(*), DNS1D(*), XKAPAO(*)
      REAL XKAPAOB1(*), D1HFAOB1(*), XKAPAOB2(*), D1HFAOB2(*)
      REAL D2GAI(*),  D2GAB(*),  D2GIJ(*),  D2GIA(*)
      REAL D2GAIB1(*),D2GABB1(*),D2GIJB1(*),D2GIAB1(*)
      REAL D2GAIB2(*),D2GABB2(*),D2GIJB2(*),D2GIAB2(*)
      REAL FTRACE, XNORM
      REAL ZERO, HALF, ONE, TWO
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XLAMDP(*),XLAMDH(*)
      DOUBLE PRECISION D1HFAO(*), XKAPAO(*), QAOS(*), DNS1D(*)
      DOUBLE PRECISION XKAPAOB1(*),D1HFAOB1(*),XKAPAOB2(*),D1HFAOB2(*)
      DOUBLE PRECISION D2GAI(*),  D2GAB(*),  D2GIJ(*),  D2GIA(*)
      DOUBLE PRECISION D2GAIB1(*),D2GABB1(*),D2GIJB1(*),D2GIAB1(*)
      DOUBLE PRECISION D2GAIB2(*),D2GABB2(*),D2GIJB2(*),D2GIAB2(*)
      DOUBLE PRECISION FTRACE, XNORM
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)

      LOGICAL SQRINT
      INTEGER ITYPE, NUMD, NUMG, IGAM, IDEL, ISTRTFCK, ISYM
      INTEGER ISYDIS, ISYMAB, ISYMPQ, ISYPQ1, ICON, ISYMP, ISYMQ
      INTEGER KAODEN, KAODSY, KOFF1, KOFF2, KOFF1D
      INTEGER KB1AODEN, KB2AODEN, KBAODSY, KEND0, KEND1, LWRK1
      INTEGER KXINT, KXTRP, KINTAO, KFCKEF, KFCK1D, KD2HFAO
      INTEGER KOFFAD, KOFFGB, KOFFAB, NBASA, ISYMA, ISYMB

* external functions:
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
* if LOCDBG set, print some debug output:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'entered CCEFFFOCK2... I1DXORD = ',I1DXORD
         WRITE (LUPRI,*) 'LABEL,LEXPEC,LFOCK:',LABEL,LEXPEC,LFOCK
         WRITE (LUPRI,*) 'ICORSY,ISYRLX,ISYMOPR:',ICORSY,ISYRLX,ISYMOPR
         WRITE (LUPRI,*) 'ISCOOR,ICOOR:',ISCOOR,ICOOR
         WRITE (LUPRI,*) 'ISYMG,ISYMD:',ISYMG,ISYMD
         WRITE (LUPRI,*) 'ISYFCK,ISYDEN:',ISYFCK,ISYDEN
      END IF

*---------------------------------------------------------------------*
* some intializations and memory allocation:
*---------------------------------------------------------------------*

      IF      ( LABEL(1:5).EQ.'HAM0 ' ) THEN
        ITYPE    = 0
        SQRINT = .FALSE.
      ELSE IF ( LABEL(1:5).EQ.'1DHAM' ) THEN
        ITYPE    = 1
        SQRINT = .FALSE.
      ELSE IF ( LABEL(1:5).EQ.'dh/dB' ) THEN
        ITYPE = 5
        SQRINT = .TRUE.
      ELSE
        WRITE (LUPRI,*) 'LABEL:',LABEL
        WRITE (LUPRI,*) 'ITYPE unknown in CCEFFFOCK2.'
        CALL QUIT('ITYPE unknown in CCEFFFOCK2.')
      END IF                                         

      IGAM = G + IBAS(ISYMG)
      IDEL = D + IBAS(ISYMD)

      NUMD = 1
      NUMG = 1

      ISYDIS = MULD2H(ISYMD,ICORSY)  ! 3-index distrib. of der. ints.
      ISYMAB = MULD2H(ISYMG,ISYDIS)  ! 2-index distrib. of der. ints.
      ISYMPQ = MULD2H(ISYMG,ISYDEN)  ! distrib. of usual density
      IF (I1DXORD.EQ.1)
     &  ISYPQ1 = MULD2H(ISYMPQ,ISYRLX) ! distrib. of 1-index trans.den.

      KFCKEF = 1
      KFCK1D = KFCKEF + N2BST(ISYFCK)
      KEND0  = KFCK1D + N2BST(ISYFCK)

      KAODEN   = KEND0
      KAODSY   = KAODEN  + N2BST(ISYMPQ)
      KEND0    = KAODSY  + N2BST(ISYMPQ)

      KD2HFAO  = KEND0
      KEND0    = KD2HFAO + N2BST(ISYMPQ)

      IF (I1DXORD.EQ.1) THEN
        KB1AODEN  = KEND0
        KB2AODEN  = KB1AODEN + N2BST(ISYPQ1)
        KBAODSY   = KB2AODEN + N2BST(ISYPQ1)
        KEND0     = KBAODSY  + N2BST(ISYPQ1)
      END IF

      KXINT  = KEND0
      IF (SQRINT) THEN
        KEND1  = KXINT  + N2BST(ISYMAB)
      ELSE
        KEND1  = KXINT  + NNBST(ISYMAB)
      END IF
      KINTAO = KEND1
      KEND1  = KINTAO + N2BST(ISYMAB)

      LWRK1  = LWORK  - KEND1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK,'Needed:',KEND1
         CALL QUIT('Insufficient space in CCEFFFOCK2. (1)')
      ENDIF

*---------------------------------------------------------------------*
* read derivative two-electron integrals and square up:
*---------------------------------------------------------------------*
c      WRITE(LUPRI,*)'calling CC_RDAOD for derivative integrals'
      CALL CC_RDAOD(WORK(KXINT),ISCOOR,ICOOR,ICORSY,
     &              IGAM,IDEL,NUMG,NUMD,
     &              WORK(KEND1),LWRK1,IRECNR,NRECS,MXBUF,NFILES,
     &              DIRINT,ITYPE,MXCOMP,LDERINT)
c      WRITE(LUPRI,*)'returned from CC_RDAOD...'

      IF (SQRINT) THEN
         CALL DCOPY(N2BST(ISYMAB),WORK(KXINT),1,WORK(KINTAO),1)
      ELSE
         CALL CCSD_SYMSQ(WORK(KXINT),ISYMAB,WORK(KINTAO))
      END IF
      
*---------------------------------------------------------------------*
* backtransform last remaining 2 indeces 2e- density to AO basis:
*---------------------------------------------------------------------*
C     IF (NSYM.EQ.1) THEN
        CALL DZERO(WORK(KAODEN),N2BST(ISYMPQ))
        CALL CC_DENAO(WORK(KAODEN),ISYMPQ,
     *                D2GAI,D2GAB,D2GIJ,D2GIA,ISYMPQ,
     *                XLAMDP,1,XLAMDH,1,WORK(KEND1),LWRK1)
                                                                   
        IF (I1DXORD.EQ.1) THEN
         CALL DZERO(WORK(KB1AODEN),N2BST(ISYPQ1))
         CALL CC_DENAO(WORK(KB1AODEN),ISYPQ1,
     *                 D2GAIB1,D2GABB1,D2GIJB1,D2GIAB1,ISYPQ1,
     *                 XLAMDP,1,XLAMDH,1,WORK(KEND1),LWRK1)

         CALL DZERO(WORK(KB2AODEN),N2BST(ISYPQ1))
         CALL CC_DENAO(WORK(KB2AODEN),ISYPQ1,
     *                 D2GAIB2,D2GABB2,D2GIJB2,D2GIAB2,ISYPQ1,
     *                 XLAMDP,1,XLAMDH,1,WORK(KEND1),LWRK1)
        END IF
C     END IF

*---------------------------------------------------------------------*
* build SCF 2-electron density matrix:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KD2HFAO),N2BST(ISYMPQ))

      ICON = 1
      CALL CC_D2EFF(WORK(KD2HFAO),G,ISYMG,IDEL,ISYMD,
     &              D1HFAO,D1HFAO,ICON)

*---------------------------------------------------------------------*
* add relaxation terms to set up effective 2e- densities:
*---------------------------------------------------------------------*
       ICON = 1
       CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,XKAPAO,D1HFAO,ICON)

       ICON = 1
       CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,D1HFAO,XKAPAO,ICON)

       IF (I1DXORD.EQ.1) THEN

        ICON = 1

        CALL CC_1IDX_D2EFF(WORK(KB1AODEN),ICON,G,ISYMG,D,ISYMD,ISYRLX,
     &               XKAPAO,D1HFAO,XKAPAOB1,D1HFAOB1,XKAPAOB2,D1HFAOB2)

        CALL CC_1IDX_D2EFF(WORK(KB2AODEN),ICON,G,ISYMG,D,ISYMD,ISYRLX,
     &               XKAPAO,D1HFAO,XKAPAOB1,D1HFAOB1,XKAPAOB2,D1HFAOB2)

       END IF

*---------------------------------------------------------------------*
* calculate the two-electron contributions to the expectation values:
*---------------------------------------------------------------------*
      IF (LEXPEC .AND. ISYMOPR.EQ.1) THEN
 
        EXPVALUE(2,IDLST) = EXPVALUE(2,IDLST) + HALF * 
     &       DDOT(N2BST(ISYMPQ),WORK(KAODEN),1,WORK(KINTAO),1)

        EXPVALUE(4,IDLST) = EXPVALUE(4,IDLST) + 
     &       DDOT(N2BST(ISYMPQ),WORK(KD2HFAO),1,WORK(KINTAO),1)

      END IF

*---------------------------------------------------------------------*
* calculate the two-electron contributions to effective Fock matrices:
*---------------------------------------------------------------------*
      IF (LFOCK) THEN

        ! symmetrize the usual density
        ! DSY(pq) = D(pq) + D(qp)
        DO ISYMP = 1, NSYM
           ISYMQ = MULD2H(ISYMP,ISYMPQ)
           KOFF1 = IAODIS(ISYMP,ISYMQ)
           KOFF2 = IAODIS(ISYMQ,ISYMP)
           CALL TRSREC(NBAS(ISYMP),NBAS(ISYMQ),
     &                 WORK(KAODEN+KOFF1),WORK(KAODSY+KOFF2))
        END DO
        CALL DAXPY(N2BST(ISYMPQ),ONE,WORK(KAODEN),1,WORK(KAODSY),1)


        IF (LOCDBG) THEN
          XNORM = DDOT(N2BST(ISYMPQ),WORK(KAODEN),1,WORK(KAODEN),1)
          WRITE (LUPRI,*) 'norm^2(kaoden):',XNORM
          XNORM = DDOT(N2BST(ISYMPQ),WORK(KAODSY),1,WORK(KAODSY),1)
          WRITE (LUPRI,*) 'norm^2(kaodsy):',XNORM
        END IF

c       --------------------------------------------------------
c       read effective Fock matrix from file and test its trace:
c       --------------------------------------------------------
        IF (I1DXORD.EQ.0) ISTRTFCK = IADRFCK(1,IDLST)
        IF (I1DXORD.EQ.1) ISTRTFCK = IADR1DXF(1,IDLST)
        
        CALL GETWA2(LUFCK,FILFCK,WORK(KFCKEF),ISTRTFCK,N2BST(ISYFCK))
     
        IF (I1DXORD.EQ.0) THEN
          ISTRTFCK = IADRFCK(2,IDLST)
          CALL GETWA2(LUFCK,FILFCK,WORK(KFCK1D),ISTRTFCK,N2BST(ISYFCK))
        END IF

        IF (LOCDBG .AND. ISYFCK.EQ.1) THEN
           FTRACE = 0.0D0
           DO ISYM = 1, NSYM
              KOFF1 = KFCKEF + IAODIS(ISYM,ISYM)
              DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
              END DO
           END DO
           WRITE (LUPRI,*) 
     *           'I1DXORD, IDLST, ISYFCK:',I1DXORD,IDLST, ISYFCK
           WRITE (LUPRI,*) 'read CC eff. Fock, TRACE = ',FTRACE
           IF (I1DXORD.EQ.0) THEN
             FTRACE = 0.0D0
             DO ISYM = 1, NSYM
                KOFF1 = KFCK1D + IAODIS(ISYM,ISYM)
                DO I = 1, NBAS(ISYM)
                   FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                END DO
             END DO
             WRITE (LUPRI,*) 'read CPHF rhs, TRACE = ',FTRACE
           END IF
        END IF                                       
        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'ccefffock2>  FCKEF matrix before ',label
          WRITE(LUPRI,*) 'ccefffock2>  norm:',
     &       DDOT(N2BST(ISYFCK),WORK(KFCKEF),1,WORK(KFCKEF),1)
          CALL CC_PRONELAO(WORK(KFCKEF),ISYFCK)
          IF (I1DXORD.EQ.0) THEN
            WRITE(LUPRI,*) 'ccefffock2>  FCK1D matrix before ',label
            CALL CC_PRONELAO(WORK(KFCK1D),ISYFCK)
          END IF
        END IF

c       -------------------------------------------------------
c       usual CC effective Fock matrices and SCF Fock matrices:
c       -------------------------------------------------------
        IF (I1DXORD.EQ.0) THEN

          ! effective Fock matrix calculated from
          ! the usual density matrix
          IF (SQRINT) THEN

            ! usual integrals / density matrix
            CALL CC_MAOMAO('T','N',HALF,WORK(KAODEN),ISYMPQ,
     &                     WORK(KINTAO),ISYMAB,ONE,WORK(KFCKEF),ISYFCK)

CORIGINAL
C           ! transposed integrals / density matrix
C           CALL CC_MAOMAO('N','N',HALF,WORK(KAODEN),ISYMPQ,
C    &                     WORK(KXTRP),ISYMAB,ONE,WORK(KFCKEF),ISYFCK)
CORIGINAL
            ! transposed integrals / density matrix
            CALL CC_MAOMAO('N','T',-HALF,WORK(KAODEN),ISYMPQ,
     &                     WORK(KINTAO),ISYMAB,ONE,WORK(KFCKEF),ISYFCK)
          ELSE

            ! use simple symmetrized density
            CALL CC_MAOMAO('N','N',HALF,WORK(KAODSY),ISYMPQ,
     &                     WORK(KINTAO),ISYMAB,ONE,WORK(KFCKEF),ISYFCK)

          END IF
       
          ! Coulomb contribution to HF Fock matrix
          IF (ISYMAB.EQ.1) THEN
            KOFF1D = IAODIS(ISYMG,ISYMD) + (D-1)*NBAS(ISYMG) + (G-1)
            WORK(KFCK1D+KOFF1D) = WORK(KFCK1D+KOFF1D) +
     &         TWO * DDOT(N2BST(ISYMAB),WORK(KINTAO),1,D1HFAO,1)
          END IF

          ! Exchange contribution to HF Fock matrix
          ISYMA = ISYMD
          ISYMB = MULD2H(ISYMA,ISYMAB)
          KOFFAD = IAODIS(ISYMA,ISYMD) + (D-1)*NBAS(ISYMA) + 1
          KOFFGB = IAODIS(ISYMG,ISYMB) + (G-1)
          KOFFAB = IAODIS(ISYMA,ISYMB)
          NBASA  = MAX(NBAS(ISYMA),1)
          IF (MULD2H(ISYMG,ISYMB).NE.ISYFCK) 
     &      CALL QUIT('symmetry problem...')

          CALL DGEMV('T',NBAS(ISYMA),NBAS(ISYMB),
     &               -ONE,WORK(KINTAO+KOFFAB),NBASA,D1HFAO(KOFFAD),1,
     &               +ONE,WORK(KFCK1D+KOFFGB),NBAS(ISYMG))
          
C         DO A = 1, NBAS(ISYMA)
C          DO B = 1, NBAS(ISYMB)
C            WORK(KFCK1D+KOFFGB+NBAS(ISYMG)*(B-1)) = 
C    &          WORK(KFCK1D+KOFFGB+NBAS(ISYMG)*(B-1)) -
C    &           WORK(KINTAO+KOFFAB+NBAS(ISYMA)*(B-1)+A-1) * 
C    &             D1HFAO(KOFFAD+NBAS(ISYMA)*(D-1)+A-1) 
C            IF (B.EQ.3 .AND. G.EQ.1) THEN
C              WRITE(LUPRI,*) 'Integral:', 
C    &             WORK(KINTAO+KOFFAB+NBAS(ISYMA)*(B-1)+A-1)
C              WRITE(LUPRI,*) 'density:', 
C    &             D1HFAO(KOFFAD+NBAS(ISYMA)*(D-1)+A)
C              WRITE(LUPRI,*) 'fock matrix:', 
C    &              WORK(KFCK1D+KOFFGB+NBAS(ISYMG)*(B-1))
C            END IF 
C          END DO
C         END DO

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'ccefffock2>  ISYMAB,ISYMA:',ISYMAB,ISYMA
            WRITE(LUPRI,*) 'ccefffock2>  FCK1D matrix after 1',label
            CALL CC_PRONELAO(WORK(KFCK1D),ISYFCK)
            WRITE(LUPRI,*) 'ccefffock2>  integral matrix at 1'
            CALL CC_PRONELAO(WORK(KINTAO),ISYFCK)
            WRITE(LUPRI,*) 'ccefffock2>  density matrix at 1'
            CALL CC_PRONELAO(D1HFAO,1)
          END IF
c       -------------------------------------------------
c       1-index transformed effective Fock matrices:
c       -------------------------------------------------
        ELSE IF (I1DXORD.EQ.1) THEN

          ! transpose the lambda^q transformed d2
          DO ISYMP = 1, NSYM
             ISYMQ = MULD2H(ISYMP,ISYPQ1)
             KOFF1 = IAODIS(ISYMP,ISYMQ)
             KOFF2 = IAODIS(ISYMQ,ISYMP)
             CALL TRSREC(NBAS(ISYMP),NBAS(ISYMQ),
     &                   WORK(KB1AODEN+KOFF1),WORK(KBAODSY+KOFF2))
          END DO

          ! add the lambda^q2 transformed d2
          CALL DAXPY(N2BST(ISYPQ1),ONE,WORK(KB2AODEN),1,WORK(KBAODSY),1)

          IF (LOCDBG) THEN
           XNORM = DDOT(N2BST(ISYPQ1),WORK(KB1AODEN),1,WORK(KB1AODEN),1)
           WRITE (LUPRI,*) 'norm^2(kb1aoden):',XNORM
           XNORM = DDOT(N2BST(ISYPQ1),WORK(KB2AODEN),1,WORK(KB2AODEN),1)
           WRITE (LUPRI,*) 'norm^2(kb2aoden):',XNORM
           XNORM = DDOT(N2BST(ISYPQ1),WORK(KBAODSY),1,WORK(KBAODSY),1)
           WRITE (LUPRI,*) 'norm^2(kbaodsy) - (1):',XNORM
          END IF

          ! add contribution from usual AODEN with
          ! the second index transformed to the
          ! barred (one-index transformed) ao basis
          CALL CC_MAOMAO('N','T',ONE,WORK(KAODSY),ISYMPQ,
     &                   QAOS,ISYRLX,ONE,WORK(KBAODSY),ISYPQ1)

          IF (LOCDBG) THEN
           XNORM = DDOT(N2BST(ISYRLX),QAOS,1,QAOS,1)
           WRITE (LUPRI,*) 'norm^2(kqaos) :',XNORM
           XNORM = DDOT(N2BST(ISYMPQ),WORK(KAODSY),1,WORK(KAODSY),1)
           WRITE (LUPRI,*) 'norm^2(kaodsy):',XNORM
           XNORM = DDOT(N2BST(ISYPQ1),WORK(KBAODSY),1,WORK(KBAODSY),1)
           WRITE (LUPRI,*) 'norm^2(kbaodsy) - (2):',XNORM
           XNORM = DDOT(N2BST(ISYMAB),WORK(KINTAO),1,WORK(KINTAO),1)
           WRITE (LUPRI,*) 'norm^2(kintao):',XNORM
           WRITE (LUPRI,*) 'ISYPQ1,ISYMAB,ISYFCK:',ISYPQ1,ISYMAB,ISYFCK
          END IF

          ! effective Fock matrix calculated from the
          ! one-index transformed density matrix
          CALL CC_MAOMAO('N','N',HALF,WORK(KBAODSY),ISYPQ1,
     &                   WORK(KINTAO),ISYMAB,ONE,WORK(KFCKEF),ISYFCK)
                                                                        
        END IF

c       --------------------------------------------------------
c       write effective Fock matrix to file and test its trace:
c       --------------------------------------------------------
        ! write matrix back to file
        IF (I1DXORD.EQ.0) ISTRTFCK = IADRFCK(1,IDLST)
        IF (I1DXORD.EQ.1) ISTRTFCK = IADR1DXF(1,IDLST)
        CALL PUTWA2(LUFCK,FILFCK,WORK(KFCKEF),ISTRTFCK,N2BST(ISYFCK))

        IF (LOCDBG .AND. ISYFCK.EQ.1) THEN
           FTRACE = 0.0D0
           DO ISYM = 1, NSYM
              KOFF1 = KFCKEF + IAODIS(ISYM,ISYM)
              DO I = 1, NBAS(ISYM)
                 FTRACE = FTRACE + WORK(KOFF1 +
     &                    NBAS(ISYM)*(I-1)+I-1)
              END DO
           END DO
           WRITE (LUPRI,*) 'written: FTRACE = ',FTRACE
        END IF
        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'ccefffock2>  FCKEF matrix before ',label
          WRITE(LUPRI,*) 'ccefffock2>  norm:',
     &       DDOT(N2BST(ISYFCK),WORK(KFCKEF),1,WORK(KFCKEF),1)
          CALL CC_PRONELAO(WORK(KFCKEF),ISYFCK)
        END IF

c       --------------------------------------------------------
c       calculate the contribution of the zeroth-order integrals
c       to the derivative SCF Fock matrix:
c       --------------------------------------------------------
        IF (I1DXORD.EQ.0) THEN

          ISYMAB = MULD2H(ISYMG,ISYMD)

          KXINT  = KEND0
          KINTAO = KXINT  + NNBST(ISYMAB)
          KEND1  = KINTAO + N2BST(ISYMAB)
          LWRK1  = LWORK  - KEND1

          IF (LWRK1.LT.0) THEN
             CALL QUIT('Insufficient space in CCEFFFOCK2. (2)')
          END IF

c          WRITE(LUPRI,*)'calling CC_RDAOD for undiff. integrals'
          CALL CC_RDAOD(WORK(KXINT),0,1,1,IGAM,IDEL,NUMG,NUMD,
     &                  WORK(KEND1),LWRK1,IRECNR,NRECS,MXBUF,NFILES,
     &                  DIRINT,0,MXCOMP,LDERINT)
c          WRITE(LUPRI,*)'returned from CC_RDAOD...'

          CALL CCSD_SYMSQ(WORK(KXINT),ISYMAB,WORK(KINTAO))

          IF (ISYFCK.NE.ISYMOPR) THEN
            WRITE (LUPRI,*) 'ISYFCK,ISYMOPR:',ISYFCK,ISYMOPR
            CALL QUIT('Symmetry mismatch in CCEFFFOCK2.')
          END IF

          IF (ISYMAB.EQ.ISYMOPR) THEN
            KOFF1D = IAODIS(ISYMG,ISYMD) + (D-1)*NBAS(ISYMG) + (G-1)
            WORK(KFCK1D+KOFF1D) = WORK(KFCK1D+KOFF1D) +
     &         TWO * DDOT(N2BST(ISYMAB),WORK(KINTAO),1,DNS1D,1)
          END IF

          ISYMA  = MULD2H(ISYMD,ISYMOPR)
          ISYMB  = MULD2H(ISYMA,ISYMAB)
          KOFFAD = IAODIS(ISYMA,ISYMD) + (D-1)*NBAS(ISYMA) + 1
          KOFFGB = IAODIS(ISYMG,ISYMB) + (G-1)
          KOFFAB = IAODIS(ISYMA,ISYMB)
          NBASA  = MAX(NBAS(ISYMA),1)

          CALL DGEMV('T',NBAS(ISYMA),NBAS(ISYMB),
     &               -ONE,WORK(KINTAO+KOFFAB),NBASA,DNS1D(KOFFAD),1,
     &               +ONE,WORK(KFCK1D+KOFFGB),NBAS(ISYMG))

c         WRITE(LUPRI,*) 'ccefffock2>  DNS1D matrix for ',label
c         CALL CC_PRONELAO(DNS1D,ISYFCK)
c         XNORM = DDOT(N2BST(ISYFCK),DNS1D,1,DNS1D,1)
c         WRITE(LUPRI,*) 'ccefffock2> norm^2 of DNS1D:',XNORM
          IF (LOCDBG) THEN
             WRITE(LUPRI,*) 'ccefffock2>  FCK1D matrix after 2',label
             CALL CC_PRONELAO(WORK(KFCK1D),ISYFCK)
          END IF

          IF (LOCDBG) THEN
             FTRACE = 0.0D0
             IF (ISYFCK.EQ.1) THEN
              DO ISYM = 1, NSYM
                KOFF1 = KFCK1D + IAODIS(ISYM,ISYM)
                DO I = 1, NBAS(ISYM)
                   FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                END DO
              END DO
             END IF
             WRITE (LUPRI,*) 'write CPHF rhs, TRACE = ',FTRACE
          END IF

          ISTRTFCK = IADRFCK(2,IDLST)
          CALL PUTWA2(LUFCK,FILFCK,WORK(KFCK1D),ISTRTFCK,N2BST(ISYFCK))

        END IF

      END IF

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'leaving CCEFFFOCK2... I1DXORD = ',I1DXORD
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_expfck_init */
      SUBROUTINE CC_EXPFCK_INIT(MOLGRD) 
*---------------------------------------------------------------------*
*
*     Purpose: Initialize the list of expectation values and effective
*              fock matrices. Purpose of the initialization is to
*              enforce the particular order of the operator lables
*              needed for the two-electron derivative integrals.
*              Extension of the list with simple one-electron operators
*              and maintainance of the logical flags for expect. values
*              and eff. Fock matrices is done via the integer functions
*              IEXPECT and IEFFFOCK.
*
*     Christof Haettig, Januar 1999
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccexpfck.h>
#include <ccroper.h>
#include <ccorb.h>
#include <ccsections.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
 
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      LOGICAL L1DHAM, MOLGRD, LDHDB
      CHARACTER*8 LABEL
      INTEGER I, IDX, IATOM, ICOOR, ICORSY, ISCOOR
 
#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER( ZERO = 0.0D0 )
 

*---------------------------------------------------------------------*
*     first element: zeroth-order Hamiltonian:
*---------------------------------------------------------------------*
      IDX = 1
      LBLEXPFCK(IDX) = 'HAM0    '
      ISYEXPFCK(IDX) = 1
      ORDEXPFCK(IDX) = 0

*---------------------------------------------------------------------*
*     if needed, first derivative integrals in the appropriate order:
*---------------------------------------------------------------------*
      ! there are at the moment three ways to ask for things which
      ! require first derivatve integrals: 
      !    1) gradient calc. enforced by DALTON    --> MOLGRD=.TRUE.
      !    2) gradient calc. requested by CC input --> CCDERI=.TRUE.
      !    3) ask for properties which need '1DHAM***' operators
      L1DHAM = MOLGRD .OR. CCDERI
      DO I = 1, NRSOLBL
         IF ( LBLOPR(I)(1:5).EQ.'1DHAM' ) L1DHAM = .TRUE.
      END DO
      
      IF (L1DHAM) THEN
   
         DO IATOM  = 1, NUCIND
         DO ICOOR  = 1, 3
         DO ICORSY = 1, NSYM
     
            ISCOOR = IPTCNT(3*(IATOM-1)+ICOOR,ICORSY-1,1)

            IF (ISCOOR.GT.0) THEN
               IDX = IDX + 1
               LABEL = '1DHAM'//CHRNOS(ISCOOR/100)
     &                        //CHRNOS(MOD(ISCOOR,100)/10)
     &                        //CHRNOS(MOD((MOD(ISCOOR,100)),10))
               LBLEXPFCK(IDX)  = LABEL
               ISYEXPFCK(IDX)  = ICORSY
               ORDEXPFCK(IDX)  = 1
               IF (LOCDBG) WRITE(2,*) 'LABEL:',LABEL
            END IF

         END DO
         END DO
         END DO

      END IF

*---------------------------------------------------------------------*
*     if needed, first magnetic derivative (London) integrals 
*---------------------------------------------------------------------*
      LDHDB = .FALSE.
      DO I = 1, NRSOLBL
         IF ( LBLOPR(I)(1:5).EQ.'dh/dB' ) LDHDB = .TRUE.
      END DO
      
      IF (LDHDB) THEN
        IDX = IDX + 1
        LBLEXPFCK(IDX)  = 'dh/dBX  '
        ISYEXPFCK(IDX)  = ISYMAX(1,2)+1
        ORDEXPFCK(IDX)  = 1

        IDX = IDX + 1
        LBLEXPFCK(IDX)  = 'dh/dBY  '
        ISYEXPFCK(IDX)  = ISYMAX(2,2)+1
        ORDEXPFCK(IDX)  = 1

        IDX = IDX + 1
        LBLEXPFCK(IDX)  = 'dh/dBZ  '
        ISYEXPFCK(IDX)  = ISYMAX(3,2)+1
        ORDEXPFCK(IDX)  = 1
      END IF
   
*---------------------------------------------------------------------*
*     set the length of the list to the number of operator we have
*     put on it and set all expect. value & eff. Fock flags to .FALSE.
*---------------------------------------------------------------------*
      NEXPFCKLBL = IDX
    
      DO IDX = 1, NEXPFCKLBL
         LEXPFCK(1,IDX) = .FALSE.
         LEXPFCK(2,IDX) = .FALSE.
      END DO
      
*---------------------------------------------------------------------*
*     print debug output and return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A)') ' CC_EXPFCK_INIT> EXPFCK list:'
        DO IDX = 1, NEXPFCKLBL
          WRITE (LUPRI,'(I5,3X,A8,2I5,2L3)') IDX,LBLEXPFCK(IDX),
     &                   ISYEXPFCK(IDX),ORDEXPFCK(IDX),LEXPFCK(1,IDX), 
     *                   LEXPFCK(2,IDX)
        END DO
      END IF
 
      RETURN
      END
*=====================================================================*
*                      END OF SUBROUTINE CC_EXPFCK_INIT               *
*=====================================================================*
*=====================================================================*
C  /* Deck iefffock */
      INTEGER FUNCTION IEFFFOCK(NEWLBL,ISYM,IORD)
*---------------------------------------------------------------------*
*
* Purpose: maintain Fock matrix entries in the list of expectation
*          values and effective Fock matrices
*          (which at present depend only on an operator label.
*           for 3. and higher derivatives a dependency on a density
*           label or index might offer the simplest way for extension.
*           see function IEXPECT for expectation value entries.)
*
*   if Fock matrix is on the list return list index 
*   if Fock matrix is NOT on the list:
*
*        LEFCKOPN=.true.  --> extend list, and return index
*        LEFCKOPN=.false. --> return -1
*
*        NEWLBL -- operator label
*        ISYM   -- symmetry
*
* Christof Haettig, Januar 1999
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>

#include <ccexpfck.h>
 
      INTEGER ISYM, IORD
      CHARACTER NEWLBL*(8)

      INTEGER I

      DO I = 1, NEXPFCKLBL
         IF ( NEWLBL .EQ. LBLEXPFCK(I) ) THEN  ! found list entry

            IEFFFOCK = I

            IF ( LEXPFCK(2,I) ) THEN 
               ! fock matrix flag already set
               ! nothing more to do
               RETURN
            ELSE IF (LEFCKOPN) THEN
               ! fock matrix flag was not set, but list
               ! extension is allowed -> set flag and return
               LEXPFCK(2,I) = .TRUE.
               RETURN
            END IF

         END IF
      END DO  

      IF (LEFCKOPN) THEN
        NEXPFCKLBL = NEXPFCKLBL + 1

        IF (NEXPFCKLBL.GT.MAXEXPFCKLBL) THEN
          WRITE(LUPRI,'(2A,/A,I5,A,I5)') '@ NUMBER OF SPECIFIED',
     *    ' FOCK MATRICES EXCEEDS THE MAXIMUM ALLOWED',
     *    '@ MAXEXPFCKLBL =',MAXEXPFCKLBL,' NEXPFCKLBL= ',NEXPFCKLBL
          CALL QUIT(' IEFFFOCK: TOO MANY FOCK MATRICES SPECIFIED')
        END IF

        LBLEXPFCK(NEXPFCKLBL)  = NEWLBL
        ISYEXPFCK(NEXPFCKLBL)  = ISYM
        ORDEXPFCK(NEXPFCKLBL)  = IORD
        LEXPFCK(1,NEXPFCKLBL)  = .FALSE.
        LEXPFCK(2,NEXPFCKLBL)  = .TRUE.
        IEFFFOCK               = NEXPFCKLBL

      ELSE
        WRITE(LUPRI,'(/3A)')
     *   '@ WARNING: EFFECTIVE FOCK MATRIX FOR ',NEWLBL,
     *   ' IS NOT AVAILABLE.'
        WRITE(LUPRI,'(/A)') ' LIST OF EFFECTIVE FOCK MATRICES:'
        DO I = 1, NEXPFCKLBL
           IF (LEXPFCK(2,I)) WRITE(LUPRI,'(I5,3X,A8,I5)') 
     &         I, LBLEXPFCK(I), ISYEXPFCK(I)
        END DO
        IEFFFOCK = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck iexpect */
      INTEGER FUNCTION IEXPECT(NEWLBL,ISYM,IORD)
*---------------------------------------------------------------------*
*
* Purpose: maintain expect. value entries in the list of expectation
*          values and effective Fock matrices
*          (which at present depend only on an operator label.
*           for 3. and higher derivatives a dependency on a density
*           label or index might offer the simplest way for extension.
*           see function IEFFFOCK for Fock matrix entries.)
*
*   if expectation value is on the list return list index 
*   if expectation value is NOT on the list:
*
*        LEXPTOPN=.true.  --> extend list, and return index
*        LEXPTOPN=.false. --> return -1
*
*        NEWLBL -- operator label
*        ISYM   -- symmetry
*
* Christof Haettig, Januar 1999
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccexpfck.h>
 
      INTEGER ISYM, IORD
      CHARACTER NEWLBL*(8)

      INTEGER I

      DO I = 1, NEXPFCKLBL
         IF ( NEWLBL .EQ. LBLEXPFCK(I) ) THEN  ! found list entry

            IEXPECT = I

            IF ( LEXPFCK(1,I) ) THEN 
               ! expect. value flag already set
               ! nothing more to do
               RETURN
            ELSE IF (LEXPTOPN) THEN
               ! expect. value flag was not set, but list
               ! extension is allowed -> set flag and return
               LEXPFCK(1,I) = .TRUE.
               RETURN
            END IF

         END IF
      END DO  

      IF (LEXPTOPN) THEN
        NEXPFCKLBL = NEXPFCKLBL + 1

        IF (NEXPFCKLBL.GT.MAXEXPFCKLBL) THEN
          WRITE(LUPRI,'(2A,/A,I5,A,I5)') '@ NUMBER OF SPECIFIED',
     *    ' EXPECTATION VALUES EXCEEDS THE MAXIMUM ALLOWED',
     *    '@ MAXEXPFCKLBL =',MAXEXPFCKLBL,' NEXPFCKLBL= ',NEXPFCKLBL
          CALL QUIT(' IEXPECT: TOO MANY EXPECTATION VALUES SPECIFIED')
        END IF

        LBLEXPFCK(NEXPFCKLBL)  = NEWLBL
        ISYEXPFCK(NEXPFCKLBL)  = ISYM
        ORDEXPFCK(NEXPFCKLBL)  = IORD
        LEXPFCK(1,NEXPFCKLBL)  = .TRUE.
        LEXPFCK(2,NEXPFCKLBL)  = .FALSE.
        IEXPECT                = NEXPFCKLBL

      ELSE
        WRITE(LUPRI,'(/3A)')
     *   '@ WARNING: EXPECTATION VALUE FOR ',NEWLBL,
     *   ' IS NOT AVAILABLE.'
        WRITE(LUPRI,'(/A)') ' LIST OF EXPECTATION VALUES:'
        DO I = 1, NEXPFCKLBL
           IF (LEXPFCK(1,I)) WRITE(LUPRI,'(I5,3X,A8,I5)') 
     &         I, LBLEXPFCK(I), ISYEXPFCK(I)
        END DO
        IEXPECT = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck i1dxfck */
      INTEGER FUNCTION I1DXFCK(NEWLBL,LSTRLX,LABELX,FREQX,ISYMX)
*---------------------------------------------------------------------*
*
* Purpose: maintain list of effective Fock matrices calculated from
*          (single) one-index transformed integrals
*
*   if Fock matrix is on the list return list index
*   if Fock matrix is NOT on the list:
*
*        L1DXFOPN=.true.  --> extend list, and return index
*        L1DXFOPN=.false. --> return -1
*
*        NEWLBL -- operator label
*
* Christof Haettig, Februar 1999
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <cc1dxfck.h>
 
      INTEGER ISYMX
      CHARACTER NEWLBL*(8), LSTRLX*(3), LABELX*(8)

#if defined (SYS_CRAY)
      REAL FREQX, TOL
#else
      DOUBLE PRECISION FREQX, TOL
#endif
      PARAMETER(TOL=1.0D-12)            

      INTEGER I, INUM, ISYMO
      INTEGER IR1KAPPA, IROPER

      DO I = 1, N1DXFLBL
         IF ( NEWLBL .EQ. LBL1DXFCK(I) .AND.
     &        LSTRLX .EQ. LST1DXFCK(I) .AND.
     &        LABELX .EQ. LABRLX(I)    .AND.
     &        (ABS(FREQX-FRQRLX(I)).LT.TOL) .AND.
     &        ISYMX  .EQ. ISYMRLX(I)      ) THEN  ! found list entry

            I1DXFCK = I
            RETURN

         END IF
      END DO  

      IF (L1DXFOPN) THEN
        N1DXFLBL = N1DXFLBL + 1

        IF (N1DXFLBL.GT.MAX1DXFCKLBL) THEN
          WRITE(LUPRI,'(2A,/A,I5,A,I5)') '@ NUMBER OF SPECIFIED',
     *    ' FOCK MATRICES EXCEEDS THE MAXIMUM ALLOWED',
     *    '@ MAX1DXFCKLBL =',MAX1DXFCKLBL,' N1DXFLBL= ',N1DXFLBL
          CALL QUIT(' I1DXFCK: TOO MANY FOCK MATRICES SPECIFIED')
        END IF

        LBL1DXFCK(N1DXFLBL)  = NEWLBL
        LST1DXFCK(N1DXFLBL)  = LSTRLX
        LABRLX(N1DXFLBL)     = LABELX
        FRQRLX(N1DXFLBL)     = FREQX
        ISYMRLX(N1DXFLBL)    = ISYMX
        I1DXFCK              = N1DXFLBL

        IF (LST1DXFCK(N1DXFLBL)(1:2).EQ.'R1') THEN
           IRELAX1DX(N1DXFLBL) = IR1KAPPA(LABELX,FREQX,ISYMX)
        ELSE
           CALL QUIT('Unknown list in I1DXFCK.')
        END IF

        ISYMO = 0
        INUM = IROPER(LABELX,ISYMO)

      ELSE
        WRITE(LUPRI,'(/3A,A3,1X,A8,1P,D12.5,A)')
     *   '@ WARNING: 1-IDX-TRAN EFFECTIVE FOCK MATRIX FOR ',NEWLBL,
     *   'WITH ORB.-RELX. VECTOR ',LSTRLX,LABELX,FREQX,
     *   ' IS NOT AVAILABLE.'
        WRITE(LUPRI,'(/A)') ' LIST OF 1-IDX-TRAN EFFECTIVE '//
     &       'FOCK MATRICES:'
        DO I = 1, N1DXFLBL
           WRITE(LUPRI,'(I5,3X,A8,1X,A3,1X,A8,1P,D12.5,I5)') 
     &       I,LBL1DXFCK(I),LST1DXFCK(I),LABRLX(I),FRQRLX(I),ISYMRLX(I)
        END DO
        I1DXFCK = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ccefffockhlp1 */
      SUBROUTINE CCEFFFOCKHLP1(IATOM,I1DXORD,DOTHISATOM,LABELH)
*---------------------------------------------------------------------*
*
* test if expectation values or effective fock matrices are to be 
* calculate which need first derivative of the integrals with respect
* to coordinates of atom no. IATOM.
*
* Christof Haettig, Februar 2000
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <cc1dxfck.h>
#include <ccexpfck.h>
#include <ccroper.h>
#include <ccropr2.h>
 
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IATOM, I1DXORD
      LOGICAL DOTHISATOM
      CHARACTER*8 LABELH

      LOGICAL LEXPEC, LFOCK, LTWO
      CHARACTER*8 LABELA, LABELB, LABEL
      INTEGER IORDER, IDLST, IOPER, JATOM, ISIGN, ISYM, ISYMA, ISYMB
      INTEGER IOPERA, IOPERB

      INTEGER IROPER, IROPER2

      LABELH     = '??????? '
      DOTHISATOM = .FALSE.


      IF (LOCDBG) WRITE (LUPRI,*) 'entered CCEFFFOCKHLP1...'

c     -----------------------------------------------------------
c     expec. values or eff. fock matrices for ordinary operators:
c     -----------------------------------------------------------
      IF (I1DXORD.EQ.0) THEN
        DO IDLST = 1, NEXPFCKLBL
          LEXPEC = LEXPFCK(1,IDLST)
          LFOCK  = LEXPFCK(2,IDLST)
          IF (LEXPEC .OR. LFOCK) THEN
            LABEL  = LBLEXPFCK(IDLST)
            IORDER = ORDEXPFCK(IDLST)
            IF (IORDER.EQ.0 .OR. IORDER.EQ.1) THEN
              IOPER = IROPER(LABEL,ISYM)
              JATOM = IATOPR(IOPER)
              LTWO  = LPDBSOP(IOPER)
            ELSE IF (IORDER.EQ.2) THEN
              LABELA = '?'
              LABELB = '?'
              IOPER  = IROPER2(LABELA,LABELB,LABEL,ISIGN,ISYM)
              JATOM  = IATOPR2(IOPER)
              IOPERA = IROPER(LABELA,ISYMA)
              IOPERB = IROPER(LABELB,ISYMB)
              LTWO   = (LPDBSOP(IOPERA) .AND. LPDBSOP(IOPERB))
            ELSE
              CALL QUIT('CCEFFFOCKHLP1> illegal operator order.')
            END IF
            IF (LOCDBG) THEN
              WRITE(LUPRI,*) 'LABEL,IORDER,JATOM,LTWO:',
     &                        LABEL,IORDER,JATOM,LTWO
            END IF
            IF (JATOM.EQ.IATOM .AND. LTWO) THEN
               DOTHISATOM = .TRUE.
               IF (LABEL(1:4).EQ.'HAM0'  .OR.
     &             LABEL(1:5).EQ.'1DHAM' .OR.
     &             LABEL(1:5).EQ.'dh/dB') THEN
                 LABELH = LABEL
               END IF
            END IF
          END IF
        END DO

c     -----------------------------------------------------
c     eff. fock matrices for 1-index transformed operators:
c     -----------------------------------------------------
      ELSE IF (I1DXORD.EQ.1) THEN
        DO IDLST = 1, N1DXFLBL
          LABEL = LBL1DXFCK(IDLST)
          IOPER = IROPER(LABEL,ISYM)
          JATOM = IATOPR(IOPER)
          LTWO  = LPDBSOP(IOPER)
          IF (JATOM.EQ.IATOM .AND. LTWO) THEN
            DOTHISATOM = .TRUE.
            IF (LABEL(1:4).EQ.'HAM0'  .OR.
     &          LABEL(1:5).EQ.'1DHAM' .OR.
     &          LABEL(1:5).EQ.'dh/dB') THEN
              LABELH = LABEL
            END IF
          END IF
        END DO

      ELSE 
        CALL QUIT('CCEFFFOCKHLP1> Illegel value for I1DXORD.')
      END IF

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'leaving CCEFFFOCKHLP1 '
        WRITE(LUPRI,*) 'LABELH     : ',LABELH
        WRITE(LUPRI,*) 'DOTHISATOM : ',DOTHISATOM
        WRITE(LUPRI,*) 'I1DXORD    : ',I1DXORD
        WRITE(LUPRI,*) 'IATOM      : ',IATOM
      END IF

      RETURN
      END
*======================================================================*
*======================================================================*
C  /* Deck ccefffockhlp2 */
      SUBROUTINE CCEFFFOCKHLP2(IDLST,JATOM,IATOM,ISYMOPR,ISYRLX,ISYFCK,
     &                         ISCOOR,ICOOR,ICORSY,MXCOMP,
     &                         LEXPEC,LFOCK,LTWO,LABEL,IOPER,
     &                         KAPPA,RMAT,QAOS,T1AM,OVERLP,
     &                         XLAMDPQ,XLAMDHQ,CMOPQ,CMOHQ,
     &                         B1DHFAO,B2DHFAO,B1KABAO,B2KABAO,
     &                         XLAMDPQ2,XLAMDHQ2,CMO,KABAR,
     &                         WORK,LWORK)
*----------------------------------------------------------------------*
*
* do some tests and preparations for the calculation of 
* one-index transformed effective Fock matrices
*
* Christof Haettig, Februar 2000
*
*======================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <cc1dxfck.h>
#include <ccroper.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <ccr1rsp.h>
#include <ccorb.h>
#include <chrxyz.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
 
      LOGICAL LTWO, LEXPEC, LFOCK
      CHARACTER*8 LABEL
      INTEGER IDLST, JATOM, IATOM, LWORK, ISYRLX, IOPER, ISYMOPR
      INTEGER ISYFCK, ISCOOR, ICOOR, ICORSY, MXCOMP

#if defined (SYS_CRAY)
      REAL KAPPA(*), RMAT(*), QAOS(*), OVERLP(*), T1AM(*)
      REAL XLAMDPQ(*), XLAMDHQ(*), CMOPQ(*), CMOHQ(*)
      REAL B1DHFAO(*), B2DHFAO(*), B1KABAO(*), B2KABAO(*)
      REAL XLAMDPQ2(*), XLAMDHQ2(*), CMO(*), KABAR(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION KAPPA(*), RMAT(*), QAOS(*), OVERLP(*), T1AM(*)
      DOUBLE PRECISION XLAMDPQ(*), XLAMDHQ(*), CMOPQ(*), CMOHQ(*)
      DOUBLE PRECISION B1DHFAO(*), B2DHFAO(*), B1KABAO(*), B2KABAO(*)
      DOUBLE PRECISION XLAMDPQ2(*), XLAMDHQ2(*), CMO(*), KABAR(*)
      DOUBLE PRECISION WORK(LWORK)
#endif

      CHARACTER LSTRLX*3, MODEL*10, LAB1*5
      INTEGER IDXRLX, IOPT, IOPRMAT, ISYRMAT, JSCOOR, JCOOR
      INTEGER IREAL, KOFDIJ, KOFDAB, KOFDAI, KOFDIA, ISDEN
      INTEGER IROPER, ILSTSYM, ILSTSYMRLX


      CALL QENTER('FCKHLP2')

      LABEL  = LBL1DXFCK(IDLST) 
      IOPER  = IROPER(LABEL,ISYMOPR)
      JATOM  = IATOPR(IOPER)
      LTWO   = LPDBSOP(IOPER)
      LEXPEC = .FALSE.
      LFOCK  = .TRUE.
      LSTRLX = LST1DXFCK(IDLST)
      IDXRLX = IRELAX1DX(IDLST)
      ISYRLX = ILSTSYMRLX(LSTRLX,IDXRLX)
      ISYFCK = MULD2H(ISYRLX,ISYMOPR)
      
      ! if the atom doesn't fit, there is nothing more to do
      IF (JATOM.NE.IATOM) RETURN

      ! retrieve orbital relaxation vector
      IF (IDXRLX.GE.1) THEN
        CALL CC_RDHFRSP(LSTRLX,IDXRLX,ISYRLX,KAPPA)
      ELSE
        CALL DZERO(KAPPA,2*NALLAI(ISYRLX))
      END IF

      ! retrieve orbital connection matrix
      IF (LSTRLX(1:2).EQ.'R1') THEN
        IOPRMAT = IROPER(LRTHFLBL(IDXRLX),ISYRMAT)
      ELSE
        WRITE (LUPRI,*) 'IOPRMAT problem in CCEFFFOCKHLP2.'
        CALL QUIT('IOPRMAT problem in CCEFFFOCKHLP2.')
      END IF

      CALL CC_GET_RMAT(RMAT,IOPRMAT,1,ISYRLX,WORK,LWORK)

      ! build the Q^ao x S^AO matrix:
      IREAL = ISYMAT(IOPRMAT)
      CALL CC_QAOS(QAOS,RMAT,KAPPA,IREAL,ISYRLX,OVERLP,WORK,LWORK)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'RMAT for ',LABEL
        CALL CC_PRONELAO(RMAT,ISYRMAT)
        WRITE (LUPRI,*) 'Q^ao x S^AO matrix for ',LABEL
        CALL CC_PRONELAO(QAOS,ISYRMAT)
      END IF


      ! build the Lambda^Q matrices:
      IOPT  = 1
      IREAL = ISYMAT(IOPRMAT)
      CALL CC_LAMBDAQ(XLAMDPQ,XLAMDHQ,CMOPQ,CMOHQ,ISYRLX,T1AM,KAPPA,
     &                RMAT,IREAL,IOPT,WORK,LWORK)


      ! build the Lambda^Q matrices with left and right Q matrices
      ! (CMOQ matrices) interchanged
      CALL DCOPY(NGLMDT(ISYRLX),CMOPQ,1,XLAMDHQ2,1)
      CALL DCOPY(NGLMDT(ISYRLX),CMOHQ,1,XLAMDPQ2,1)
      CALL LAMDA2(XLAMDPQ2,XLAMDHQ2,ISYRLX,T1AM,ISYM0,
     &            CMOHQ,CMOPQ,ISYRLX)
      
      
      ! calculate HF density with alpha index transformed      
      CALL DZERO(B2DHFAO,N2BST(ISYRLX))
      CALL CC_DHFAO(B2DHFAO,ISYRLX,CMO,ISYM0,CMOHQ,ISYRLX,WORK,LWORK)

      ! calculate HF density with beta index transformed      
      CALL DZERO(B1DHFAO,N2BST(ISYRLX))
      CALL CC_DHFAO(B1DHFAO,ISYRLX,CMOPQ,ISYRLX,CMO,ISYM0,WORK,LWORK)


      ! set pointer for subblocks of kappa-bar-0
      KOFDIJ = 1
      KOFDAB = KOFDIJ + NMATIJ(1)
      KOFDAI = KOFDAB + NMATAB(1)
      KOFDIA = KOFDAI + NT1AMX


      ! calculate zeta-kappa-bar-0 density with alpha index
      ! transformed with CMOHQ
      CALL DZERO(B2KABAO,N2BST(ISYRLX))
      CALL CC_DENAO(B2KABAO,ISYRLX,KABAR(KOFDAI),KABAR(KOFDAB),
     &              KABAR(KOFDIJ),KABAR(KOFDIA),ISYM0,
     &              CMO,1,CMOHQ,ISYRLX,WORK,LWORK)

      ! calculate zeta-kappa-bar-0 density with beta index
      ! transformed with CMOHQ
      CALL DZERO(B1KABAO,N2BST(ISYRLX))
      CALL CC_DENAO(B1KABAO,ISYRLX,KABAR(KOFDAI),KABAR(KOFDAB),
     &              KABAR(KOFDIJ),KABAR(KOFDIA),ISYM0,
     &              CMOPQ,ISYRLX,CMO,1,WORK,LWORK)

      IF (JATOM.EQ.IATOM .AND. LTWO) THEN
        
        IF      (LABEL(1:8).EQ.'HAM0    ') THEN

          MXCOMP = 1
          ISCOOR = 0
          ICORSY = 1
          ICOOR  = 0

        ELSE IF (LABEL(1:5).EQ.'dh/dB') THEN

          MXCOMP = 3
          ICORSY = ISYMOPR
          DO JCOOR = 1, MXCOMP
             IF (CHRXYZ(JCOOR).EQ.LABEL(6:6)) THEN
                ICOOR = JCOOR
             END IF
          END DO
          ISCOOR = ICOOR
          IF (LOCDBG) WRITE (LUPRI,*)LABEL,ISCOOR,MXCOMP,ICOOR,ICORSY

        ELSE IF (LABEL(1:5).EQ.'1DHAM') THEN

          MXCOMP = 3
          ICORSY = ISYMOPR
          READ(LABEL,'(A5,I3)') LAB1,ISCOOR
          DO JCOOR  = 1, MXCOMP
             JSCOOR = IPTCNT(3*(IATOM-1)+JCOOR,ICORSY-1,1)
             IF (JSCOOR.EQ.ISCOOR) THEN
                ICOOR = JCOOR
             END IF
          END DO 
          IF (LOCDBG) WRITE (LUPRI,*)
     *                LABEL,ISCOOR,MXCOMP,ICOOR,ICORSY

        ELSE
          CALL QUIT('Illegal or unknown label in CCEFFFOCKHLP2.')
        END IF

      END IF
    
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'will enter two-electron part for:'
        WRITE (LUPRI,*) 'IDLST:',IDLST
        WRITE (LUPRI,*) 'LABEL,ISYMOPR,ISYRLX:',LABEL,ISYMOPR,ISYRLX
        WRITE (LUPRI,*) 'IOPER:',IOPER
        WRITE (LUPRI,*) 'IATOM,ISCOOR:',IATOM,ISCOOR
        WRITE (LUPRI,*) 'LTWO :',LTWO
        WRITE (LUPRI,*) 'ICOOR,ICORSY,MXCOMP:',ICOOR,ICORSY,MXCOMP
c       WRITE (LUPRI,*) 'B1DHFAO:'
c       CALL CC_PRONELAO(B1DHFAO,ISYRLX)
c       WRITE (LUPRI,*) 'B2DHFAO:'
c       CALL CC_PRONELAO(B2DHFAO,ISYRLX)
c       WRITE (LUPRI,*) 'B1KABAO:'
c       CALL CC_PRONELAO(B1KABAO,ISYRLX)
c       WRITE (LUPRI,*) 'B2KABAO:'
c       CALL CC_PRONELAO(B2KABAO,ISYRLX)
      END IF

      CALL QEXIT('FCKHLP2')
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ccefffockhlp3 */
      SUBROUTINE CCEFFFOCKHLP3(IDLST,IATOM,JATOM,ISYMOPR,ISYFCK,
     &                         LABEL,IOPER,IORDER,
     &                         MXCOMP,ISCOOR,ICORSY,ICOOR,
     &                         LEXPEC,LFOCK,LTWO,DNS1D,WORK,LWORK)
*---------------------------------------------------------------------*
*
* do some tests for the calculation of effective Fock matrices
* and expectation values
*
* Christof Haettig, March 2000
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccroper.h>
#include <ccropr2.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <ccr1rsp.h>
#include <chrxyz.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
 
      LOGICAL LTWO, LEXPEC, LFOCK
      INTEGER IDLST, JATOM, IATOM, MXCOMP, ISCOOR, ICORSY, ICOOR, LWORK
      INTEGER ISYMOPR, ISYFCK, IORDER, IOPER, ISIGN
      CHARACTER LABEL*8

#if defined (SYS_CRAY)
      REAL DNS1D(*), WORK(LWORK)
#else
      DOUBLE PRECISION DNS1D(*), WORK(LWORK)
#endif

      CHARACTER LABELA*8, LABELB*8, LAB1*5
      INTEGER IOPERA, IOPERB, JCOOR, JSCOOR, ISYMA, ISYMB, ISYM
      INTEGER IROPER, IROPER2

* test if anything to do:
      IF (LOCDBG) WRITE (LUPRI,*) 'entered CCEFFFOCKHLP3...'

      LEXPEC = LEXPFCK(1,IDLST)
      LFOCK  = LEXPFCK(2,IDLST)

      IF ( .NOT. (LEXPEC.OR.LFOCK) ) RETURN

* ok, enter routine:
      CALL QENTER('FCKHLP3')

      LABEL  = LBLEXPFCK(IDLST)
      IORDER = ORDEXPFCK(IDLST)

      IF (IORDER.EQ.0 .OR. IORDER.EQ.1) THEN
        IOPER = IROPER(LABEL,ISYM)
        JATOM = IATOPR(IOPER)
        LTWO  = LPDBSOP(IOPER)
      ELSE IF (IORDER.EQ.2) THEN
        LABELA = '?'
        LABELB = '?'
        IOPER  = IROPER2(LABELA,LABELB,LABEL,ISIGN,ISYM)
        JATOM  = IATOPR2(IOPER)
        IOPERA = IROPER(LABELA,ISYMA)
        IOPERB = IROPER(LABELB,ISYMB)
        LTWO   = ( LPDBSOP(IOPERA) .AND. LPDBSOP(IOPERB) )
      ELSE
        CALL QUIT('CCEFFFOCKHLP3> illegal operator order.')
      END IF

      ISYMOPR = ISYEXPFCK(IDLST)
      ISYFCK  = ISYMOPR

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'two-electron part for:'
        WRITE (LUPRI,*) 'IDLST:',IDLST
        WRITE (LUPRI,*) 'LABEL:',LABEL, ISYMOPR
        WRITE (LUPRI,*) 'IOPER:',IOPER, IORDER
        WRITE (LUPRI,*) 'IATOM:',IATOM, JATOM
        WRITE (LUPRI,*) 'LEXPEC, LFOCK:',LEXPEC,LFOCK
        WRITE (LUPRI,*) 'LTWO:',LTWO
        WRITE (LUPRI,*) 'LWORK:',LWORK
      END IF

      IF (LFOCK) THEN
        IF (LABEL.EQ.'HAM0    ') THEN
          CALL DZERO(DNS1D,N2BST(ISYMOPR))
        ELSE
          CALL CC_HFR1DEN(DNS1D,IOPER,IORDER,ISYMOPR,WORK,LWORK)
        END IF
      END IF

      IF (JATOM.EQ.IATOM .AND. LTWO) THEN
        
        IF      (LABEL(1:8).EQ.'HAM0    ') THEN

          MXCOMP = 1
          ISCOOR = 0
          ICORSY = 1
          ICOOR  = 0

        ELSE IF (LABEL(1:5).EQ.'dh/dB') THEN

          MXCOMP = 3
          ICORSY = ISYMOPR
          DO JCOOR = 1, MXCOMP
             IF (CHRXYZ(JCOOR).EQ.LABEL(6:6)) THEN
                ICOOR = JCOOR
             END IF
          END DO
          ISCOOR = ICOOR
          IF (LOCDBG) WRITE (LUPRI,*)LABEL,ISCOOR,MXCOMP,ICOOR,ICORSY

        ELSE IF (LABEL(1:5).EQ.'1DHAM') THEN

          MXCOMP = 3
          ICORSY = ISYMOPR
          READ(LABEL,'(A5,I3)') LAB1,ISCOOR
          DO JCOOR  = 1, MXCOMP
             JSCOOR = IPTCNT(3*(IATOM-1)+JCOOR,ICORSY-1,1)
             IF (JSCOOR.EQ.ISCOOR) THEN
                ICOOR = JCOOR
             END IF
          END DO 
          IF (LOCDBG) WRITE (LUPRI,*)LABEL,ISCOOR,MXCOMP,ICOOR,ICORSY

        ELSE
          CALL QUIT('Illegal or unknown label in CCEFFFOCKHLP3.')
        END IF

      END IF

      CALL QEXIT('FCKHLP3')
      RETURN
      END
*=====================================================================*
