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
c/* Deck cc_wrrsp */
      SUBROUTINE CC_WRRSP(LIST,IDXLST,ISYM,IOPT,MODFIL,
     &                    VEC0,VEC1,VEC2,WORK,LWORK)
C---------------------------------------------------------------------*
C
C   Purpose:  Write a response vector on file: 
C             - the LIST and the list index IDXLST are used to 
C               construct a file name and to access the information 
C               for the file header
C             - MODFILE is written to file header as a Label to
C               identify the CC model.
C             - IOPT is used as bit wise flags to determine 
C               which parts of the vector are to be written on file:
C                 1.bit  (1): singles part   VEC1
C                 2.bit  (2): doubles part   VEC2
C                 3.bit  (4): cphf part      VEC0
C                 4.bit  (8): CC3 eff. singles part VEC1
C                 5.bit (16): CC3 eff. doubles part VEC2
C                 6.bit (32): R12 doubles part VEC2
C
C  implemented LISTs:
C
C          'R0' : right zeroth-order amplitude vectors (IDXLST ignored)
C          'R1' : right first-order amplitude vectors
C          'R2' : right second-order amplitude vectors
C          'R3' : right third-order amplitude vectors 
C          'R4' : right fourth-order amplitude vectors 
C          'RC' : right first-order Cauchy vectors
C          'CR2': right second-order Cauchy vectors
C          'RE' : right eigenvectors
C          'ER1': first-order response of right eigenvectors
C          'ER2': second-order response of right eigenvectors
C
C          'L0' : left zeroth-order lamdba vectors (IDXLST ignored)
C          'L1' : left first-order lambda vectors
C          'L2' : left second-order lambda vectors 
C          'L3' : left third-order lambda vectors 
C          'L4' : left fourth-order lambda vectors 
C          'LC' : left first-order Cauchy vectors
C          'CL2': left second-order Cauchy vectors
C          'LE' : left eigenvectors
C          'EL1': first-order response of left eigenvectors
C          'EL2': second-order response of left eigenvectors
C          'PL1': projected L1 vectors 
C
C          'F1' : F-transformed right first-order response vectors.   
C          'F2' : F-transformed right second-order response vectors.   
C          'F3' : F-transformed right third-order response vectors.   
C          'F4' : F-transformed right fourth-order response vectors.   
C          'FC' : F-transformed right first-order Cauchy vectors.
C          'CF2': F-transformed right second-order Cauchy vectors.
C          'FR' : F-transformed right eigenvectors.
C
C          'XF1': F-transformed right first-order response vectors,
C                 F12*R2 contribution for Cholesky CC2.
C
C          'O1 ': rhs for first-order amplitude equations
C          'O1e': effective 'O1' vector for triples models
C          'O2' : rhs for second-order amplitude equations
C          'O3' : rhs for third-order amplitude equations
C          'O4' : rhs for fourth-order amplitude equations
C          'CO2': rhs for second-order right Cauchy equations
C          'EO1': rhs for first-order response of right eigenvectors
C          'EO2': rhs for second-order response of right eigenvectors
C
C          'eO1': effective rhs for first-order amplitude equations for
C                 Cholesky CC2
C
C          'X1 ': first-order chi/eta vectors
C          'X1e': effective 'X1' vector for triples models
C          'X2' : second-order chi vectors
C          'X3' : third-order chi vectors
C          'X4' : fourth-order chi vectors
C          'CX2': second-order Cauchy eta vectors
C          'EX1': chi vector for left eigenvector first-order response
C          'EX2': chi vector for left eigenvector second-order response
C
C          'M1': zeroth-order Lagrangian multipliers for ground state
C                -- excited state transition moments
C
C          'BE': rhs for zeroth-order excited state Lagrangian
C                multipliers. (LE*B*RE)
C          'E0': zeroth-order excited state Lagrangian multipliers.
C                (the same as the N2 vectors, but diagonal case only)
C
C          'BR': rhs for N2 vectors
C 
C          'N2': zeroth-order Lagrangian multiplers for
C                excited state -- excited state transition moments
C                and excited state properties
C
C          'D0': dummy vector (used f.x. in CCLR_FANA)
C
C  Christof Haettig, october 1996, restructered in november 1998
C  PL1 vectors introduced in March 2000 by Sonia Coriani
C  CC3 effective vectors introduced 2002 by Christof Haettig
C  R12 vectors introduced Jun 2003 by Christof Haettig
C---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfro.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccl1rsp.h>
#include <cco1rsp.h>
#include <ccx1rsp.h>
#include <ccr2rsp.h>
#include <ccl2rsp.h>
#include <cco2rsp.h>
#include <ccx2rsp.h>
#include <ccr3rsp.h>
#include <ccl3rsp.h>
#include <cco3rsp.h>
#include <ccx3rsp.h>
#include <ccr4rsp.h>
#include <ccl4rsp.h>
#include <cco4rsp.h>
#include <ccx4rsp.h>
#include <ccn2rsp.h>
#include <cclrmrsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccl2rsp.h>
#include <cccx2rsp.h>
#include <ccexci.h>
#include <ccpl1rsp.h>
#include <dummy.h>
Cholesky
#include <maxorb.h>
#include <cclr.h>
#include <ccdeco.h>
Cholesky
 
      CHARACTER SUBRNAME*(8)
      INTEGER LUSAVE
      PARAMETER (SUBRNAME = 'CC_WRRSP')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      CHARACTER LIST*(*), LISTI*(4)
      CHARACTER FILEX*(10)
      CHARACTER MODFIL*(10), MODEL*(10)
      CHARACTER LABELA*(8)

      LOGICAL LEXISTS
      INTEGER IOS, NVEC, LWORK, KORBITAL, KSINGLE, KDOUBLE, KR12DBLE,
     &        KSINGLEFF, KDOUBLEFF, NVECR, IBIT32, IDXVEC
      INTEGER IDXLST, ISYM, IMUL, IDXFIL, IOPT, KEND, LWRK, NWRITE

#if defined (SYS_CRAY)
      REAL VEC0(*), VEC1(*),VEC2(*), WORK(LWORK)
      REAL DDOT, XXNORM
#else
      DOUBLE PRECISION VEC0(*), VEC1(*),VEC2(*), WORK(LWORK)
      DOUBLE PRECISION DDOT, XXNORM
#endif

* external functions:
      INTEGER IDXSYM
      INTEGER ILSTSYM

      CALL QENTER('CC_WRRSP')
*---------------------------------------------------------------------*
*     get file name and modified LIST label and IDXVEC:
*---------------------------------------------------------------------*
      CALL CC_RWPRE(LIST,IDXLST,ISYM,LISTI,IDXVEC,FILEX)

*---------------------------------------------------------------------*
* check compatibility of IOPT with the MODFIL, if not reset IOPT:   
* set number of components that have to be written:
*---------------------------------------------------------------------*
      IBIT32 = 0
      IF (CCR12) IBIT32 = 32

      IF (MODFIL(1:3).EQ.'SCF') THEN
        IOPT   = IAND(IOPT,4) 
        NWRITE = 1
C
Cholesky 
C
      ELSE IF (MODFIL(1:3).EQ.'CC2' .AND. CHOINT) THEN
        IOPT   = IAND(IOPT,5)
        NWRITE = 2
C     
Cholesky
C 
      ELSE IF ((MODFIL(1:3).EQ.'CCS').AND.
     *    (.NOT.(MODFIL(1:4).EQ.'CCSD'))) THEN
        IOPT   = IAND(IOPT,1+4) 
        NWRITE = 2
      ELSE IF (MODFIL(1:3).EQ.'MP2'.OR.MODFIL(1:3).EQ.'CC2'.OR.
     *         MODFIL(1:3).EQ.'CCD'.OR.MODFIL(1:4).EQ.'CCSD'
     *        ) THEN
        IOPT   = IAND(IOPT,1+2+4+IBIT32)
        NWRITE = 4
      ELSE IF (MODFIL(1:3).EQ.'CC3'      .OR.
     *         MODFIL(1:8).EQ.'CCSDT-1a' .OR.
     *         MODFIL(1:8).EQ.'CCSDT-1B' .OR.
     *         MODFIL(1:7).EQ.'CCSD(T)'  .OR.
     *         MODFIL(1:5).EQ.'CC(3)'    .OR.
     *         MODFIL(1:8).EQ.'CCSDR(T)' .OR.
     *         MODFIL(1:8).EQ.'CCSDR(3)' .OR.
     *         MODFIL(1:9).EQ.'CCSDR(1A)'.OR.
     *         MODFIL(1:9).EQ.'CCSDR(1B)'   
     *        ) THEN
        IOPT   = IAND(IOPT,1+2+4+8+16+IBIT32)
        NWRITE = 6
      ELSE 
        WRITE (LUPRI,*) 'Unknown model in CC_WRRSP: ',MODFIL
        CALL QUIT('MODEL not yet implemented in CC_WRRSP.')
      END IF

*---------------------------------------------------------------------*
* find the multiplicity of the vector:
*---------------------------------------------------------------------*
      ! default is singlett
      IMUL = 1

      ! for excited states we get it from the common block:
      IF (LISTI(1:2).EQ.'RE'  .OR. LISTI(1:2).EQ.'LE') THEN
         IMUL = IMULTE(IDXVEC)
      END IF                              

      IF (LOCDBG) WRITE(LUPRI,*) 'CC_WRRSP> FILE:',FILEX 
*---------------------------------------------------------------------*
* open & rewind file:
*---------------------------------------------------------------------*
      INQUIRE(FILE=FILEX,EXIST=LEXISTS)

      LUSAVE = -1
      CALL GPOPEN(LUSAVE,FILEX,'UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)

      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

      KORBITAL  = 1
      KSINGLE   = 1
      KDOUBLE   = 1
      KSINGLEFF = 1
      KDOUBLEFF = 1
      KR12DBLE  = 1
      KEND      = 1

CCH   IF (LEXISTS .AND. LISTI(1:2).NE.'R0') THEN
      IF (LEXISTS) THEN
C
C        check model and how many vectors were stored on file
C
         READ(LUSAVE,END=888,ERR=888) NVEC, MODEL
C
C        for CPHF we want sometimes not change MODEL on file
C
         IF (IOPT.EQ.4 .AND. MODFIL(1:4).EQ.'SCF?') THEN
            MODFIL = MODEL
         END IF
C
C        don't keep higher excitations from previous runs
C
         NVEC = MIN(NVEC,NWRITE)
C
C        read CPHF vector... 
C
         IF ( NVEC.GE.1 .AND. IAND(IOPT,4).EQ.0 ) THEN
            KORBITAL = KEND
            KEND     = KORBITAL + 2*NALLAI(ISYM)
            IF (LWORK.LT.KEND-1) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP (cphf).')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KORBITAL-1+I),I=1,2*NALLAI(ISYM))
            IF (LOCDBG) WRITE(LUPRI,*) 'read orbital part'
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
            IF (LOCDBG) WRITE(LUPRI,*) 'skipped orbital part'
         END IF
         NVECR = 1
C
C        read CC singles vector... 
C
         IF (NVEC.GE.2 .AND. IAND(IOPT,1).EQ.0 ) THEN
            KSINGLE  = KEND
            KEND     = KSINGLE  + NT1AM(ISYM)
            IF (LWORK.LT.KEND-1) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP (singles).')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KSINGLE-1+I),I=1,NT1AM(ISYM))
            IF (LOCDBG) WRITE(LUPRI,*) 'read singles part'
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
            IF (LOCDBG) WRITE(LUPRI,*) 'skipped singles part'
         END IF
         NVECR = 2
C
C        read CC doubles vector... 
C
         IF (NVEC.GE.3 .AND. IAND(IOPT,2).EQ.0 ) THEN
            KDOUBLE  = KEND
            KEND     = KDOUBLE  + NT2AM(ISYM)
            IF (LWORK.LT.KEND-1) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP (doubles).')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KDOUBLE-1+I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              KEND     = KDOUBLE  + NT2AM(ISYM)+NT2AMA(ISYM)
              READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KDOUBLE+NT2AM(ISYM)-1+I),I=1,NT2AMA(ISYM))
            ELSE
              READ(LUSAVE,ERR=888,END=888) 
            END IF
            IF (LOCDBG) WRITE(LUPRI,*) 'read doubles part'
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
            READ(LUSAVE,ERR=888,END=888) 
            IF (LOCDBG) WRITE(LUPRI,*) 'skipped doubles part'
         END IF
         NVECR = 3
C
C        read CC R12 doubles vector... 
C
         IF (NVEC.GE.4 .AND. IAND(IOPT,32).EQ.0 .AND. CCR12) THEN
            KR12DBLE = KEND
            KEND     = KR12DBLE  + NTR12AM(ISYM)
            IF (LWORK.LT.KEND-1) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP (r12 dbls).')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KR12DBLE-1+I),I=1,NTR12AM(ISYM))
            IF (IMUL.EQ.3) THEN
              CALL QUIT('Triplet for R12 not implemented in CC_WRRSP.')
C             KEND     = KDOUBLE  + NT2AM(ISYM)+NT2AMA(ISYM)
C             READ(LUSAVE,ERR=888,END=888) 
C    &            (WORK(KDOUBLE+NT2AM(ISYM)-1+I),I=1,NT2AMA(ISYM))
            END IF
            IF (LOCDBG) WRITE(LUPRI,*) 'read r12 doubles part'
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
            READ(LUSAVE,ERR=888,END=888) 
            IF (LOCDBG) WRITE(LUPRI,*) 'skipped r12 doubles part'
         END IF
         NVECR = 4
C
C        read CC singles of effective vector... 
C
         IF (NVEC.GE.5 .AND. IAND(IOPT,8).EQ.0 ) THEN
            KSINGLEFF  = KEND
            KEND       = KSINGLEFF  + NT1AM(ISYM)
            IF (LWORK.LT.KEND-1) THEN
              CALL QUIT('Insufficient memory in CC_WRRSP (eff. sgls.).')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KSINGLEFF-1+I),I=1,NT1AM(ISYM))
            IF (LOCDBG) WRITE(LUPRI,*) 'read eff. singles part'
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
            IF (LOCDBG) WRITE(LUPRI,*) 'skipped eff. singles part'
         END IF
         NVECR = 5
C
C        read CC doubles of effective vector... 
C
         IF (NVEC.GE.6 .AND. IAND(IOPT,16).EQ.0 ) THEN
            KDOUBLEFF  = KEND
            KEND       = KDOUBLEFF  + NT2AM(ISYM)
            IF (LWORK.LT.KEND-1) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP (eff. dbls).')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KDOUBLEFF-1+I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              KEND     = KDOUBLEFF  + NT2AM(ISYM)+NT2AMA(ISYM)
              READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KDOUBLEFF+NT2AM(ISYM)-1+I),I=1,NT2AMA(ISYM))
            END IF
            IF (LOCDBG) WRITE(LUPRI,*) 'read eff. doubles part'
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
            READ(LUSAVE,ERR=888,END=888) 
            IF (LOCDBG) WRITE(LUPRI,*) 'skipped eff. doubles part'
         END IF
         NVECR = 6

888      CONTINUE
    
         ! for some vectors we don't have an 'effective' triples part
         ! and we won't write zero vectors instead:
         IF (NWRITE.EQ.6 .AND. IAND(IOPT,24).LT.8 .AND. NVEC.LT.5) THEN
           NWRITE = 4
         END IF

      ELSE
         ! for some vectors we don't have an 'effective' triples part
         ! and we won't write zero vectors instead:
         IF (NWRITE.EQ.6 .AND. IAND(IOPT,24).LT.8) NWRITE = 4
      END IF

      ! if there are not effective singles and doubles (vectors 5,6)
      ! and this is not an R12 run, there is no need to worry about
      ! the R12 part, which is stored on file in record 4
      IF (NWRITE.EQ.4 .AND. (.NOT.CCR12)) NWRITE = 3

      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

*=====================================================================*
* write a 'LISTI' and index specific header to the file:
*=====================================================================*

*.....................................................................*
* 'R0' - zeroth-order t amplitudes:
* 'L0' - zeroth-order lagrangian multipliers:
* 'D0' - dummy vector
*.....................................................................*
      IF ( LISTI(1:2) .EQ. 'R0' .OR. LISTI(1:2) .EQ. 'L0' .OR.
     &     LISTI(1:2) .EQ. 'D0'                               ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL

*.....................................................................*
* 'R1' - first-order response t amplitudes:
* 'F1' - F transformed first-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'R1' .OR. LISTI(1:2) .EQ. 'F1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), LRTLBL(IDXVEC), ISYLRT(IDXVEC), FRQLRT(IDXVEC),
     &    LORXLRT(IDXVEC)

*.....................................................................*
* 'XF1' - F transformed first-order response t amplitudes,
*         F12*R2 contribution for Cholesky CC2: 
* 'eO1' - eff. rhs for first-order response t amplitudes for
*         Cholesky CC2:
*.....................................................................*

      ELSE IF ( LIST(1:3) .EQ. 'XF1' .OR. LIST(1:3) .EQ. 'eO1') THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:3), LRTLBL(IDXLST), ISYLRT(IDXLST), FRQLRT(IDXLST),
     &    LORXLRT(IDXLST)

*.....................................................................*
* 'R2' - second-order response t amplitudes:
* 'F2' - F transformed second-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'R2' .OR. LISTI(1:2) .EQ. 'F2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLR2T(IDXVEC,1), ISYR2T(IDXVEC,1), FRQR2T(IDXVEC,1),
     &       LBLR2T(IDXVEC,2), ISYR2T(IDXVEC,2), FRQR2T(IDXVEC,2)

*.....................................................................*
* 'R3' - third-order response t amplitudes:
* 'F3' - F transformed third-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'R3' .OR. LISTI(1:2) .EQ. 'F3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLR3T(IDXVEC,1), ISYR3T(IDXVEC,1), FRQR3T(IDXVEC,1),
     &       LBLR3T(IDXVEC,2), ISYR3T(IDXVEC,2), FRQR3T(IDXVEC,2),
     &       LBLR3T(IDXVEC,3), ISYR3T(IDXVEC,3), FRQR3T(IDXVEC,3)

*.....................................................................*
* 'R4' - fourth-order response t amplitudes:
* 'F4' - F transformed fourth-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'R4' .OR. LISTI(1:2) .EQ. 'F4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLR4T(IDXVEC,1), ISYR4T(IDXVEC,1), FRQR4T(IDXVEC,1),
     &       LBLR4T(IDXVEC,2), ISYR4T(IDXVEC,2), FRQR4T(IDXVEC,2),
     &       LBLR4T(IDXVEC,3), ISYR4T(IDXVEC,3), FRQR4T(IDXVEC,3),
     &       LBLR4T(IDXVEC,4), ISYR4T(IDXVEC,4), FRQR4T(IDXVEC,4)

*.....................................................................*
* 'RE' - Right eigenvector.
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'RE' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), ISYEXC(IDXVEC), IMULTE(IDXVEC), EIGVAL(IDXVEC)

*.....................................................................*
* 'ER1' - Right eigenvector, first-order response
* 'EO1' - RHS for right eigenvector, first-order response
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'ER1' .OR. LISTI(1:3) .EQ. 'EO1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    ISTER1(IDXVEC), ISYSER1(IDXVEC), EIGER1(IDXVEC),
     &    LBLER1(IDXVEC), ISYOER1(IDXVEC), FRQER1(IDXVEC)

*.....................................................................*
* 'ER2' - Right eigenvector, second-order response
* 'EO2' - RHS for right eigenvector, second-order response
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'ER2' .OR. LISTI(1:3) .EQ. 'EO2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    ISTER2(IDXVEC),   ISYSER2(IDXVEC),   EIGER2(IDXVEC),
     &    LBLER2(IDXVEC,1), ISYOER2(IDXVEC,1), FRQER2(IDXVEC,1),
     &    LBLER2(IDXVEC,2), ISYOER2(IDXVEC,2), FRQER2(IDXVEC,2)

*.....................................................................*
* 'EL1' - Left eigenvector, first-order response
* 'EX1' - chi vector for left eigenvector, first-order response
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'EL1' .OR. LISTI(1:3) .EQ. 'EX1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    ISTEL1(IDXVEC), ISYSEL1(IDXVEC), EIGEL1(IDXVEC),
     &    LBLEL1(IDXVEC), ISYOEL1(IDXVEC), FRQEL1(IDXVEC)

*.....................................................................*
* 'EL2' - Left eigenvector, second-order response
* 'EX2' - chi vector for left eigenvector, second-order response
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'EL2' .OR. LISTI(1:3) .EQ. 'EX2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    ISTEL2(IDXVEC),   ISYSEL2(IDXVEC),   EIGEL2(IDXVEC),
     &    LBLEL2(IDXVEC,1), ISYOEL2(IDXVEC,1), FRQEL2(IDXVEC,1),
     &    LBLEL2(IDXVEC,2), ISYOEL2(IDXVEC,2), FRQEL2(IDXVEC,2)

*.....................................................................*
* 'RC' - Right Cauchy vector.
* 'FC' - F transformed right Cauchy vector.
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'RC' .OR. LISTI(1:2) .EQ. 'FC' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), LRCLBL(IDXVEC), ISYLRC(IDXVEC), ILRCAU(IDXVEC)

*.....................................................................*
* 'Rx' - x-th. order response t amplitudes:
*.....................................................................*
      ELSE IF ( LISTI(1:1) .EQ. 'R' ) THEN
         CALL QUIT ( LISTI(1:2)//'-th. order t amplit.'//
     &             ' not yet available' )

*.....................................................................*
* 'CR2' - second-order right Cauchy vectors
* 'CF2' - F-transformed second-order right Cauchy vectors
*.....................................................................*
      ELSE IF (LISTI(1:3).EQ.'CR2' .OR. LISTI(1:3).EQ.'CF2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    LBLCR2(IDXVEC,1), ISYCR2(IDXVEC,1), ICR2CAU(IDXVEC,1),
     &    LBLCR2(IDXVEC,2), ISYCR2(IDXVEC,2), ICR2CAU(IDXVEC,2)

*.....................................................................*
* 'CO2' - RHS for second-order right Cauchy vectors
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'CO2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    LBLCO2(IDXVEC,1), ISYCO2(IDXVEC,1), ICO2CAU(IDXVEC,1),
     &    LBLCO2(IDXVEC,2), ISYCO2(IDXVEC,2), ICO2CAU(IDXVEC,2)

*.....................................................................*
* 'CL2' - second-order left Cauchy vectors
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'CL2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    LBLCL2(IDXVEC,1), ISYCL2(IDXVEC,1), ICL2CAU(IDXVEC,1),
     &    LBLCL2(IDXVEC,2), ISYCL2(IDXVEC,2), ICL2CAU(IDXVEC,2)

*.....................................................................*
* 'CX2' - second-order Cauchy eta vectors
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'CX2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &    LBLCX2(IDXVEC,1), ISYCX2(IDXVEC,1), ICX2CAU(IDXVEC,1),
     &    LBLCX2(IDXVEC,2), ISYCX2(IDXVEC,2), ICX2CAU(IDXVEC,2)

*.....................................................................*
* 'E0' and 'BE' - zeroth order lagrangian multipliers and rhs.
*.....................................................................*
      ELSE IF ((LISTI(1:2).EQ.'E0').OR.(LISTI(1:2).EQ.'BE')) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL

*.....................................................................*
* 'N2' and 'BR' - zeroth order lagrangian multipliers and rhs.
*.....................................................................*
      ELSE IF ((LISTI(1:2).EQ.'N2').OR.(LISTI(1:2).EQ.'BR')) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3),
     &       IIN2(IDXVEC), ISYIN2(IDXVEC), FRQIN2(IDXVEC),
     &       IFN2(IDXVEC), ISYFN2(IDXVEC), FRQFN2(IDXVEC)

*.....................................................................*
* 'L1' - first-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'L1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), LRZLBL(IDXVEC), ISYLRZ(IDXVEC), FRQLRZ(IDXVEC),
     &    LORXLRZ(IDXVEC)

*.....................................................................*
* 'L2' - second-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'L2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLL2(IDXVEC,1), ISYL2(IDXVEC,1), FRQL2(IDXVEC,1),
     &       LBLL2(IDXVEC,2), ISYL2(IDXVEC,2), FRQL2(IDXVEC,2)

*.....................................................................*
* 'L3' - third-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'L3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLL3(IDXVEC,1), ISYL3(IDXVEC,1), FRQL3(IDXVEC,1),
     &       LBLL3(IDXVEC,2), ISYL3(IDXVEC,2), FRQL3(IDXVEC,2),
     &       LBLL3(IDXVEC,3), ISYL3(IDXVEC,3), FRQL3(IDXVEC,3)

*.....................................................................*
* 'L4' - third-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'L4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLL4(IDXVEC,1), ISYL4(IDXVEC,1), FRQL4(IDXVEC,1),
     &       LBLL4(IDXVEC,2), ISYL4(IDXVEC,2), FRQL4(IDXVEC,2),
     &       LBLL4(IDXVEC,3), ISYL4(IDXVEC,3), FRQL4(IDXVEC,3),
     &       LBLL4(IDXVEC,4), ISYL4(IDXVEC,4), FRQL4(IDXVEC,4)

*.....................................................................*
* 'LE' - Left eigenvector.
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'LE' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), ISYEXC(IDXVEC), IMULTE(IDXVEC), EIGVAL(IDXVEC)

*.....................................................................*
* 'LC' - first-order left Cauchy vector.
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'LC' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), LBLLC1(IDXVEC), ISYLC1(IDXVEC), ILC1CAU(IDXVEC)

*.....................................................................*
* 'Lx' - x-th. order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LISTI(1:1) .EQ. 'L' ) THEN
         CALL QUIT ( LISTI(1:2)//'-th. order lagrang. multip.'//
     &             ' not yet available' )

*.....................................................................*
* 'O1' - rhs vectors for first-order amplitude equations:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'O1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &     LBLO1(IDXVEC), ISYO1(IDXVEC), FRQO1(IDXVEC), LORXO1(IDXVEC)

*.....................................................................*
* 'O2' - rhs vectors for second-order amplitude equations:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'O2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLO2(IDXVEC,1), ISYO2(IDXVEC,1), FRQO2(IDXVEC,1),
     &       LBLO2(IDXVEC,2), ISYO2(IDXVEC,2), FRQO2(IDXVEC,2)

*.....................................................................*
* 'O3' - rhs vectors for third-order amplitude equations:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'O3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLO3(IDXVEC,1), ISYO3(IDXVEC,1), FRQO3(IDXVEC,1),
     &       LBLO3(IDXVEC,2), ISYO3(IDXVEC,2), FRQO3(IDXVEC,2),
     &       LBLO3(IDXVEC,3), ISYO3(IDXVEC,3), FRQO3(IDXVEC,3)

*.....................................................................*
* 'O4' - rhs vectors for fourth-order amplitude equations:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'O4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLO4(IDXVEC,1), ISYO4(IDXVEC,1), FRQO4(IDXVEC,1),
     &       LBLO4(IDXVEC,2), ISYO4(IDXVEC,2), FRQO4(IDXVEC,2),
     &       LBLO4(IDXVEC,3), ISYO4(IDXVEC,3), FRQO4(IDXVEC,3),
     &       LBLO4(IDXVEC,4), ISYO4(IDXVEC,4), FRQO4(IDXVEC,4)

*.....................................................................*
* 'X1' - first-order chi/eta vectors
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'X1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &     LBLX1(IDXVEC), ISYX1(IDXVEC), FRQX1(IDXVEC), LORXX1(IDXVEC)

*.....................................................................*
* 'X2' - second-order chi vectors:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'X2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLX2(IDXVEC,1), ISYX2(IDXVEC,1), FRQX2(IDXVEC,1),
     &       LBLX2(IDXVEC,2), ISYX2(IDXVEC,2), FRQX2(IDXVEC,2)

*.....................................................................*
* 'X3' - third-order chi vectors:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'X3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLX3(IDXVEC,1), ISYX3(IDXVEC,1), FRQX3(IDXVEC,1),
     &       LBLX3(IDXVEC,2), ISYX3(IDXVEC,2), FRQX3(IDXVEC,2),
     &       LBLX3(IDXVEC,3), ISYX3(IDXVEC,3), FRQX3(IDXVEC,3)

*.....................................................................*
* 'X4' - fourth-order chi vectors:
*.....................................................................*
      ELSE IF ( LISTI(1:2) .EQ. 'X4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LISTI(1:3), 
     &       LBLX4(IDXVEC,1), ISYX4(IDXVEC,1), FRQX4(IDXVEC,1),
     &       LBLX4(IDXVEC,2), ISYX4(IDXVEC,2), FRQX4(IDXVEC,2),
     &       LBLX4(IDXVEC,3), ISYX4(IDXVEC,3), FRQX4(IDXVEC,3),
     &       LBLX4(IDXVEC,4), ISYX4(IDXVEC,4), FRQX4(IDXVEC,4)

*.....................................................................*
* 'M1' and 'FR' - first order response t amplitudes:
*.....................................................................*
      ELSE IF ((LISTI(1:2).EQ.'M1').OR.(LISTI(1:2).EQ.'FR')) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), ILRM(IDXVEC), ISYLRM(IDXVEC), FRQLRM(IDXVEC)

*.....................................................................*
* 'PL1' - projected first-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LISTI(1:3) .EQ. 'PL1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LISTI(1:3), LBLPL1(IDXVEC), ISYPL1(IDXVEC), FRQPL1(IDXVEC),
     &    LORXPL1(IDXVEC)
*.....................................................................*
* unknown list:
*.....................................................................*
      ELSE
         CALL QUIT ('unknown list '//LISTI(1:2)//' in '//SUBRNAME)
      END IF

*---------------------------------------------------------------------*
* write the vectors, close file & return:
*---------------------------------------------------------------------*
      CALL CC_WVRSP(ISYM,IMUL,IOPT,NWRITE,LUSAVE,CCR12,VEC0,VEC1,VEC2,
     &              WORK(KORBITAL),WORK(KSINGLE),WORK(KDOUBLE),
     &              WORK(KR12DBLE),WORK(KSINGLEFF),WORK(KDOUBLEFF))

      CALL GPCLOSE(LUSAVE,'KEEP')

      CALL QEXIT('CC_WRRSP')
      RETURN

*---------------------------------------------------------------------*
* handle i/o errors:
*---------------------------------------------------------------------*
991   CONTINUE
      WRITE(LUPRI,'(2A)') ' an error occured while opening file ',FILEX
      GOTO 999

992   CONTINUE
      WRITE(LUPRI,'(2A)') ' i/o error while rewinding file ',FILEX
      GOTO 999

993   CONTINUE
      WRITE(LUPRI,'(A)') ' write error on unit LUSAVE in CC_WRRSP:',
     &     FILEX
      GOTO 999

995   CONTINUE
      WRITE(LUPRI,'(2A)') ' i/o error while closing file ',FILEX
      GOTO 999

999   CONTINUE
      WRITE(LUPRI,'(A,I5)') ' unit number    :',LUSAVE
      WRITE(LUPRI,'(A,I5)') ' returned IOSTAT:',IOS
      CALL QUIT ('fatal i/o error in '//SUBRNAME)


      END 
*---------------------------------------------------------------------*
*                    END OF SUBROUTINE CC_WRRSP
*---------------------------------------------------------------------*
c/* Deck cc_wvrsp */
      SUBROUTINE CC_WVRSP(ISYM,IMUL,IOPT,NWRITE,LUSAVE,CCR12,
     *                    VEC0,VEC1,VEC2,VS0,VS1,VS2,VS2P,VS3,VS4)
C---------------------------------------------------------------------*
C
C   Purpose:  Write a vector with symmetry ISYM to unit LUSAVE
C             for an explanation of IOPT see CC_WRRSP
C
C---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      INTEGER LUSAVE
C
      LOGICAL CCR12
      INTEGER ISYM, IOPT

      DIMENSION VEC0(*),VEC1(*),VEC2(*)
      DIMENSION VS0(*),VS1(*),VS2(*),VS2P(*),VS3(*),VS4(*)
C
      INTEGER IOS

* write vector components to file:
      ! CPHF part for orbital relaxed response 
      IF (NWRITE.GE.1) THEN
         IF ( IAND(IOPT,4).EQ.4 ) THEN
           WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC0(I),I=1,2*NALLAI(ISYM))
         ELSE
           WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS0(I),I=1,2*NALLAI(ISYM))
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF
  
      ! singles cluster amplitudes
      IF (NWRITE.GE.2) THEN
         IF ( IAND(IOPT,1).EQ.1 ) THEN
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC1(I),I=1,NT1AM(ISYM))
         ELSE 
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS1(I),I=1,NT1AM(ISYM))
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! doubles cluster amplitudes
      IF (NWRITE.GE.3) THEN
         IF ( IAND(IOPT,2).EQ.2 ) THEN
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC2(I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)
     &            (VEC2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))     
            ELSE
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)   
            END IF
         ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS2(I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)
     &            (VS2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))     
            ELSE
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)   
            END IF
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! R12 doubles cluster amplitudes
      IF (NWRITE.GE.4 .AND. CCR12) THEN
         IF ( IAND(IOPT,32).EQ.32 ) THEN
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC2(I),I=1,NTR12AM(ISYM))
            IF (IMUL.EQ.3) THEN
              CALL QUIT('Triplet not implemented for R12 in CC_WVRSP')
C             WRITE(LUSAVE,IOSTAT=IOS,ERR=993)
C    &            (VEC2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))     
            ELSE
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)   
            END IF
         ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS2P(I),I=1,NTR12AM(ISYM))
            IF (IMUL.EQ.3) THEN
              CALL QUIT('Triplet not implemented for R12 in CC_WVRSP')
C             WRITE(LUSAVE,IOSTAT=IOS,ERR=993)
C    &            (VS2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))     
            ELSE
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)   
            END IF
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! effective singles amplitudes for CC3 RHS vectors
      IF (NWRITE.GE.5) THEN
         IF ( IAND(IOPT,8).EQ.8 ) THEN
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC1(I),I=1,NT1AM(ISYM))
         ELSE 
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS3(I),I=1,NT1AM(ISYM))
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! effective doubles amplitudes for CC3 RHS vectors
      IF (NWRITE.GE.6) THEN
         IF ( IAND(IOPT,16).EQ.16 ) THEN
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC2(I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)
     &            (VEC2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))     
            ELSE
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)   
            END IF
         ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS4(I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)
     &            (VS4(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))     
            ELSE
              WRITE(LUSAVE,IOSTAT=IOS,ERR=993)   
            END IF
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      RETURN

*---------------------------------------------------------------------*
* handle i/o error:
*---------------------------------------------------------------------*
993   CONTINUE
      WRITE(LUPRI,'(A)') ' write error on unit LUSAVE in CC_WVRSP:'
      WRITE(LUPRI,'(A,I5)') ' returned IOSTAT:',ios
      CALL QUIT ('fatal i/o error in CC_WVRSP')

      END 
*---------------------------------------------------------------------*
*                    END OF SUBROUTINE CC_WVRSP
*---------------------------------------------------------------------*
c/* Deck cc_warsp */
      SUBROUTINE CC_WARSP(LIST,IDXLST,ISYM,IOPT,MODFIL,
     &                    VEC0,VEC1,VEC2,WORK,LWORK)
C---------------------------------------------------------------------*
C
C   Purpose:  Add a vector to a vector stored on file and
C             write the result vector back to the file
C
C             program is aborted if the file not yet exists to
C             prevent undefined results.
C
C             - IOPT is used as bit wise flag to determine 
C               which parts of the vector are to be written on file:
C                 1.bit  (1): singles part
C                 2.bit  (2): doubles part
C                 3.bit  (4): cphf part
C                 4.bit  (8): singles part of effective vector
C                 5.bit (16): doubles part of effective vector
C                 6.bit (32): R12 doubles part 
C
C               allowed values are 1, 2, 8, 16, 32 and 
C               the combinations 1+2, 8+16
C
C             for further input description see routines
C             CC_RDRSP and CC_WRRSP
C
C  Christof Haettig, April 1997
C
C---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <ccexci.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER IDXLST, ISYM, IOPT, IOPTRW, IMUL
      CHARACTER LIST*(*)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK), DDOT
      REAL VEC0(*), VEC1(*),VEC2(*)
#else
      DOUBLE PRECISION WORK(LWORK), DDOT
      DOUBLE PRECISION VEC0(*), VEC1(*),VEC2(*)
#endif
 
      CHARACTER MODFIL*(10)
      INTEGER KVEC0, KVEC1, KVEC2, KEND, LEND, NTAMP

*---------------------------------------------------------------------*
* allocate work space and read old vector from file:
*---------------------------------------------------------------------*
      IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
        WRITE (LUPRI,*) 'Symmetry out of range in CC_WARSP:'
        WRITE (LUPRI,*) 'ISYM:',ISYM
        WRITE (LUPRI,*) 'NSYM:',NSYM
        CALL QUIT('Symmetry out of range in CC_WARSP.')
      END IF

      ! default is singlett
      IMUL = 1

      ! for excited states we get it from the common block:
      IF (LIST(1:2).EQ.'RE'  .OR. LIST(1:2).EQ.'LE') THEN
         IMUL = IMULTE(IDXLST)
      END IF                              
      IF (IMUL .NE. 1) THEN
        CALL QUIT('CC_WARSP> triplet vector not yet implemented.')
      END IF

      IF (IOPT.NE.IAND(IOPT,1+2) .AND. IOPT.NE.IAND(IOPT,8+16) .AND.
     &    IOPT.NE.IAND(IOPT,32)  ) THEN
        WRITE(LUPRI,*) 'CC_WARSP> IOPT = ',IOPT
        CALL QUIT('IOPT not implemented in CC_WARSP.')
      END IF

   
      IF (IOPT.EQ.IAND(IOPT,1+2) .OR. IOPT.EQ.IAND(IOPT,8+16)) THEN
        KVEC0 = 1
        KVEC1 = KVEC0 + 2*NALLAI(ISYM)
        KVEC2 = KVEC1 + NT1AM(ISYM)
        KEND  = KVEC2 + NT2AM(ISYM)
      ELSE IF (IOPT.EQ.32) THEN
        KVEC1 = -999 999
        KVEC2 = 1
        KEND  = KVEC2 + NTR12AM(ISYM)
      ELSE
        CALL QUIT('Illegal case in CC_WARSP.')
      END IF
    
      LEND  = LWORK - KEND
      IF (LEND .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_WARSP.')
      END IF

      IF      (IOPT.EQ.IAND(IOPT,1+2)) THEN
         IOPTRW = 3
      ELSE IF (IOPT.EQ.IAND(IOPT,8+16)) THEN
         IOPTRW = 24
      ELSE IF (IOPT.EQ.32) THEN
         IOPTRW = 32
      ELSE
        CALL QUIT('Illegal case in CC_WARSP.')
      END IF

      CALL CC_RDRSP(LIST,IDXLST,ISYM,IOPTRW,MODFIL,
     &              WORK(KVEC1),WORK(KVEC2))
      
C     IF ( IAND(IOPT,4).EQ.4 ) THEN
C        CALL CC_RDHFRSP(LIST,IDXLST,ISYM,WORK(KVEC0))
C        IOPTRW = IOPTRW + 4
C     END IF

      IF (LOCDBG) THEN
        NTAMP = NT1AM(ISYM)
        IF (IOPT.NE.1 .AND. IOPT.NE.8) NTAMP = NTAMP + NT2AM(ISYM)
        WRITE(LUPRI,*) 'CC_WARSP> LIST,IDXLST,IOPT:',LIST,IDXLST,IOPTRW
        WRITE(LUPRI,*) 'CC_WARSP> norm^2(vec-old):',
     &    DDOT(NTAMP,WORK(KVEC1),1,WORK(KVEC1),1)
        WRITE(LUPRI,*) 'CC_WARSP> norm^2(vec-new):',
     &    DDOT(NTAMP,VEC1,1,VEC1,1)
      END IF

* check return value of read option:
      IF (IOPTRW .LT. IOPT) THEN
        WRITE (LUPRI,*) 'IOPT in CC_WARSP incompatible with vector '//
     &                  'on file.'
        WRITE (LUPRI,*) 'IOPTRW, IOPT:',IOPTRW,IOPT
        CALL QUIT('IOPT in CC_WARSP incompatible with vector on file.')
      END IF

*---------------------------------------------------------------------*
* add new vector to old vector and write result back to file:
*---------------------------------------------------------------------*
      IF ( IAND(IOPT,1+8).NE.0 ) THEN
        CALL DAXPY(NT1AM(ISYM),1.0d0,VEC1,1,WORK(KVEC1),1)
      END IF
      IF ( IAND(IOPT,2+16).NE.0 ) THEN
        CALL DAXPY(NT2AM(ISYM),1.0d0,VEC2,1,WORK(KVEC2),1)
      END IF
      IF ( IAND(IOPT,32).NE.0 ) THEN
        CALL DAXPY(NTR12AM(ISYM),1.0d0,VEC2,1,WORK(KVEC2),1)
      END IF
C     IF ( IAND(IOPT,4).NE.0 ) THEN
C       CALL DAXPY(2*NALLAI(ISYM),1.0d0,VEC0,1,WORK(KVEC0),1)
C     END IF

      CALL CC_WRRSP(LIST,IDXLST,ISYM,IOPTRW,MODFIL,
     &              WORK(KVEC0),WORK(KVEC1),WORK(KVEC2),
     &              WORK(KEND),LEND)
  
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CC_WARSP> norm^2(vec-written):',
     &    DDOT(NTAMP,WORK(KVEC1),1,WORK(KVEC1),1)
      END IF

      RETURN

      END 
*---------------------------------------------------------------------*
*                  END OF SUBROUTINE CC_WARSP
*---------------------------------------------------------------------*
