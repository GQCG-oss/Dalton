C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
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
C                 1.bit (1): singles part   VEC2
C                 2.bit (2): double part    VEC1
C                 3.bit (4): cphf part      VEC0
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
C
C          'O1' : rhs for first-order amplitude equations
C          'O2' : rhs for second-order amplitude equations
C          'O3' : rhs for third-order amplitude equations
C          'O4' : rhs for fourth-order amplitude equations
C          'CO2': rhs for second-order right Cauchy equations
C          'EO1': rhs for first-order response of right eigenvectors
C          'EO2': rhs for second-order response of right eigenvectors
C
C          'X1' : first-order chi/eta vectors
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
C
C---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit_dec.h>
#include <priunit.h>
#include <ccsdsym_dec.h>
#include <ccsdsym.h>
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
#include <dummy_dec.h>
#include <dummy.h>
 
      CHARACTER SUBRNAME*(8)
      INTEGER LUSAVE
      PARAMETER (SUBRNAME = 'CC_WRRSP')
 
      CHARACTER LIST*(*), LISTI*(4)
      CHARACTER FILEX*(10)
      CHARACTER MODFIL*(10), MODEL*(10)
      CHARACTER LABELA*(8)

      LOGICAL LEXISTS
      INTEGER IOS, NVEC, LWORK, KORBITAL, KSINGLE, KDOUBLE
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

C----------------------------------------------------------------------
C     Check symmetry and list index range:
C----------------------------------------------------------------------
      IF (ISYM .NE. ILSTSYM(LIST,IDXLST) ) THEN
        WRITE (LUPRI,*) 'symmetry mismatch in CC_WRRSP.'
        WRITE (LUPRI,*) 'LIST:',LIST
        WRITE (LUPRI,*) 'input symmetry:',ISYM
        WRITE (LUPRI,*) 'symmetry returned form ILSTSYM:',
     &       ILSTSYM(LIST,IDXLST)
        CALL QUIT('symmetry mismatch in CC_WRRSP.')
      END IF

C----------------------------------------------------------------------
C     Make symmetry-adapted filenames and 
C     find symmetry reduced list numbers.
C----------------------------------------------------------------------
      IF ( LIST(1:2).EQ.'R0' .OR. LIST(1:2).EQ.'L0' .OR.
     &     LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'L1' .OR.
     &     LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'L2' .OR.
     &     LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'L3' .OR.
     &     LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'L4' .OR.
     &     LIST(1:2).EQ.'O1' .OR. LIST(1:2).EQ.'X1' .OR.
     &     LIST(1:2).EQ.'O2' .OR. LIST(1:2).EQ.'X2' .OR.
     &     LIST(1:2).EQ.'O3' .OR. LIST(1:2).EQ.'X3' .OR.
     &     LIST(1:2).EQ.'O4' .OR. LIST(1:2).EQ.'X4' .OR.
     &     LIST(1:2).EQ.'F1' .OR. LIST(1:2).EQ.'F2' .OR.
     &     LIST(1:2).EQ.'F3' .OR. LIST(1:2).EQ.'F4' .OR.
     &     LIST(1:2).EQ.'RC' .OR. LIST(1:2).EQ.'LC' .OR.
     &     LIST(1:2).EQ.'RE' .OR. LIST(1:2).EQ.'LE' .OR.
     &     LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE' .OR.
     &     LIST(1:2).EQ.'N2' .OR. LIST(1:2).EQ.'BR' .OR.
     &     LIST(1:2).EQ.'FR' .OR. LIST(1:2).EQ.'M1' .OR.
     &     LIST(1:2).EQ.'FC' .OR. LIST(1:2).EQ.'D0'      ) THEN

        LISTI = LIST(1:2)//'_'
      ELSE
        LISTI = LIST(1:3)
      END IF

      IF (ISYM .EQ. 1) LISTI = LISTI(1:3)//'1'
      IF (ISYM .EQ. 2) LISTI = LISTI(1:3)//'2'
      IF (ISYM .EQ. 3) LISTI = LISTI(1:3)//'3'
      IF (ISYM .EQ. 4) LISTI = LISTI(1:3)//'4'
      IF (ISYM .EQ. 5) LISTI = LISTI(1:3)//'5'
      IF (ISYM .EQ. 6) LISTI = LISTI(1:3)//'6'
      IF (ISYM .EQ. 7) LISTI = LISTI(1:3)//'7'
      IF (ISYM .EQ. 8) LISTI = LISTI(1:3)//'8'
      IDXFIL = IDXSYM(LIST,ISYM,IDXLST)

*---------------------------------------------------------------------*
* build file name:
*---------------------------------------------------------------------*
      WRITE(FILEX,'(A2,A4,1X,I3)') 'CC', LISTI(1:4), IDXFIL
      DO I = 1, 10
        IF ( FILEX(I:I) .EQ. ' ' ) FILEX(I:I) = '_'
      END DO

* dummy vector:
      IF ( LIST(1:2) .EQ. 'D0' ) THEN
        FILEX = 'CCDUMMY'
      END IF

*---------------------------------------------------------------------*
* check compatibility of IOPT with the MODFIL, if not reset IOPT:   
* set number of components that have to be written:
*---------------------------------------------------------------------*
      IF (MODFIL(1:3).EQ.'SCF') THEN
        IOPT   = IAND(IOPT,4) 
        NWRITE = 1
      ELSE IF ((MODFIL(1:3).EQ.'CCS').AND.
     *    (.NOT.(MODFIL(1:4).EQ.'CCSD'))) THEN
        IOPT   = IAND(IOPT,5) 
        NWRITE = 2
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
        IOPT   = IAND(IOPT,7)
        NWRITE = 3
      ELSE IF (MODFIL(1:3).EQ.'MP2'.OR.MODFIL(1:3).EQ.'CC2'.OR.
     *         MODFIL(1:3).EQ.'CCD'.OR.MODFIL(1:4).EQ.'CCSD'
     *        ) THEN
        IOPT   = IAND(IOPT,7)
        NWRITE = 3
      ELSE 
        WRITE (LUPRI,*) 'unkown model in CC_WRRSP:',MODFIL
        CALL QUIT('triples not yet implemented in CC_WRRSP.')
      END IF

*---------------------------------------------------------------------*
* find the multiplicity of the vector:
*---------------------------------------------------------------------*
      ! default is singlett
      IMUL = 1

      ! for excited states we get it from the common block:
      IF (LIST(1:2).EQ.'RE'  .OR. LIST(1:2).EQ.'LE') THEN
         IMUL = IMULTE(IDXLST)
      END IF                              

*---------------------------------------------------------------------*
* open & rewind file:
*---------------------------------------------------------------------*
      INQUIRE(FILE=FILEX,EXIST=LEXISTS)

      LUSAVE = 0
      CALL GPOPEN(LUSAVE,FILEX,'UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)

      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

      KORBITAL = 1
      KSINGLE  = 1
      KDOUBLE  = 1
      KEND     = 1

      IF (LEXISTS .AND. LIST(1:2).NE.'R0') THEN
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
            IF (LWORK.LT.KEND) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP.')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KORBITAL-1+I),I=1,2*NALLAI(ISYM))
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
         END IF
C
C        read CC singles vector... 
C
         IF (NVEC.GE.2 .AND. IAND(IOPT,1).EQ.0 ) THEN
            KSINGLE  = KEND
            KEND     = KSINGLE  + NT1AM(ISYM)
            IF (LWORK.LT.KEND) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP.')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KSINGLE-1+I),I=1,NT1AM(ISYM))
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
         END IF
C
C        read CC doubles vector... 
C
         IF (NVEC.GE.3 .AND. IAND(IOPT,2).EQ.0 ) THEN
            KDOUBLE  = KEND
            KEND     = KDOUBLE  + NT2AM(ISYM)
            IF (LWORK.LT.KEND) THEN
               CALL QUIT('Insufficient memory in CC_WRRSP.')
            END IF
            READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KDOUBLE-1+I),I=1,NT2AM(ISYM))
            IF (IMUL.EQ.3) THEN
              KEND     = KDOUBLE  + NT2AM(ISYM)+NT2AMA(ISYM)
              READ(LUSAVE,ERR=888,END=888) 
     &            (WORK(KDOUBLE+NT2AM(ISYM)-1+I),I=1,NT2AMA(ISYM))
            END IF
         ELSE
            READ(LUSAVE,ERR=888,END=888) 
         END IF

888      CONTINUE

      END IF

      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

*=====================================================================*
* write a 'LIST' and index specific header to the file:
*=====================================================================*

*.....................................................................*
* 'R0' - zeroth-order t amplitudes:
* 'L0' - zeroth-order lagrangian multipliers:
* 'D0' - dummy vector
*.....................................................................*
      IF ( LIST(1:2) .EQ. 'R0' .OR. LIST(1:2) .EQ. 'L0' .OR.
     &     LIST(1:2) .EQ. 'D0'                               ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL

*.....................................................................*
* 'R1' - first-order response t amplitudes:
* 'F1' - F transformed first-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'R1' .OR. LIST(1:2) .EQ. 'F1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), LRTLBL(IDXLST), ISYLRT(IDXLST), FRQLRT(IDXLST),
     &    LORXLRT(IDXLST)

*.....................................................................*
* 'R2' - second-order response t amplitudes:
* 'F2' - F transformed second-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'R2' .OR. LIST(1:2) .EQ. 'F2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLR2T(IDXLST,1), ISYR2T(IDXLST,1), FRQR2T(IDXLST,1),
     &       LBLR2T(IDXLST,2), ISYR2T(IDXLST,2), FRQR2T(IDXLST,2)

*.....................................................................*
* 'R3' - third-order response t amplitudes:
* 'F3' - F transformed third-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'R3' .OR. LIST(1:2) .EQ. 'F3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLR3T(IDXLST,1), ISYR3T(IDXLST,1), FRQR3T(IDXLST,1),
     &       LBLR3T(IDXLST,2), ISYR3T(IDXLST,2), FRQR3T(IDXLST,2),
     &       LBLR3T(IDXLST,3), ISYR3T(IDXLST,3), FRQR3T(IDXLST,3)

*.....................................................................*
* 'R4' - fourth-order response t amplitudes:
* 'F4' - F transformed fourth-order response t amplitudes:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'R4' .OR. LIST(1:2) .EQ. 'F4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLR4T(IDXLST,1), ISYR4T(IDXLST,1), FRQR4T(IDXLST,1),
     &       LBLR4T(IDXLST,2), ISYR4T(IDXLST,2), FRQR4T(IDXLST,2),
     &       LBLR4T(IDXLST,3), ISYR4T(IDXLST,3), FRQR4T(IDXLST,3),
     &       LBLR4T(IDXLST,4), ISYR4T(IDXLST,4), FRQR4T(IDXLST,4)

*.....................................................................*
* 'RE' - Right eigenvector.
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'RE' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), ISYEXC(IDXLST), IMULTE(IDXLST), EIGVAL(IDXLST)

*.....................................................................*
* 'ER1' - Right eigenvector, first-order response
* 'EO1' - RHS for right eigenvector, first-order response
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'ER1' .OR. LIST(1:3) .EQ. 'EO1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    ISTER1(IDXLST), ISYSER1(IDXLST), EIGER1(IDXLST),
     &    LBLER1(IDXLST), ISYOER1(IDXLST), FRQER1(IDXLST)

*.....................................................................*
* 'ER2' - Right eigenvector, second-order response
* 'EO2' - RHS for right eigenvector, second-order response
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'ER2' .OR. LIST(1:3) .EQ. 'EO2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    ISTER2(IDXLST),   ISYSER2(IDXLST),   EIGER2(IDXLST),
     &    LBLER2(IDXLST,1), ISYOER2(IDXLST,1), FRQER2(IDXLST,1),
     &    LBLER2(IDXLST,2), ISYOER2(IDXLST,2), FRQER2(IDXLST,2)

*.....................................................................*
* 'EL1' - Left eigenvector, first-order response
* 'EX1' - chi vector for left eigenvector, first-order response
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'EL1' .OR. LIST(1:3) .EQ. 'EX1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    ISTEL1(IDXLST), ISYSEL1(IDXLST), EIGEL1(IDXLST),
     &    LBLEL1(IDXLST), ISYOEL1(IDXLST), FRQEL1(IDXLST)

*.....................................................................*
* 'EL2' - Left eigenvector, second-order response
* 'EX2' - chi vector for left eigenvector, second-order response
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'EL2' .OR. LIST(1:3) .EQ. 'EX2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    ISTEL2(IDXLST),   ISYSEL2(IDXLST),   EIGEL2(IDXLST),
     &    LBLEL2(IDXLST,1), ISYOEL2(IDXLST,1), FRQEL2(IDXLST,1),
     &    LBLEL2(IDXLST,2), ISYOEL2(IDXLST,2), FRQEL2(IDXLST,2)

*.....................................................................*
* 'RC' - Right Cauchy vector.
* 'FC' - F transformed right Cauchy vector.
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'RC' .OR. LIST(1:2) .EQ. 'FC' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), LRCLBL(IDXLST), ISYLRC(IDXLST), ILRCAU(IDXLST)

*.....................................................................*
* 'Rx' - x-th. order response t amplitudes:
*.....................................................................*
      ELSE IF ( LIST(1:1) .EQ. 'R' ) THEN
         CALL QUIT ( LIST(1:2)//'-th. order t amplit.'//
     &             ' not yet available' )

*.....................................................................*
* 'CR2' - second-order right Cauchy vectors
* 'CF2' - F-transformed second-order right Cauchy vectors
*.....................................................................*
      ELSE IF (LIST(1:3).EQ.'CR2' .OR. LIST(1:3).EQ.'CF2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    LBLCR2(IDXLST,1), ISYCR2(IDXLST,1), ICR2CAU(IDXLST,1),
     &    LBLCR2(IDXLST,2), ISYCR2(IDXLST,2), ICR2CAU(IDXLST,2)

*.....................................................................*
* 'CO2' - RHS for second-order right Cauchy vectors
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'CO2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    LBLCO2(IDXLST,1), ISYCO2(IDXLST,1), ICO2CAU(IDXLST,1),
     &    LBLCO2(IDXLST,2), ISYCO2(IDXLST,2), ICO2CAU(IDXLST,2)

*.....................................................................*
* 'CL2' - second-order left Cauchy vectors
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'CL2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    LBLCL2(IDXLST,1), ISYCL2(IDXLST,1), ICL2CAU(IDXLST,1),
     &    LBLCL2(IDXLST,2), ISYCL2(IDXLST,2), ICL2CAU(IDXLST,2)

*.....................................................................*
* 'CX2' - second-order Cauchy eta vectors
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'CX2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:3),
     &    LBLCX2(IDXLST,1), ISYCX2(IDXLST,1), ICX2CAU(IDXLST,1),
     &    LBLCX2(IDXLST,2), ISYCX2(IDXLST,2), ICX2CAU(IDXLST,2)

*.....................................................................*
* 'E0' and 'BE' - zeroth order lagrangian multipliers and rhs.
*.....................................................................*
      ELSE IF ((LIST(1:2).EQ.'E0').OR.(LIST(1:2).EQ.'BE')) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL

*.....................................................................*
* 'N2' and 'LBR' - zeroth order lagrangian multipliers and rhs.
*.....................................................................*
      ELSE IF ((LIST(1:2).EQ.'N2').OR.(LIST(1:2).EQ.'BR')) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2),
     &       IIN2(IDXLST), ISYIN2(IDXLST), FRQIN2(IDXLST),
     &       IFN2(IDXLST), ISYFN2(IDXLST), FRQFN2(IDXLST)

*.....................................................................*
* 'L1' - first-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'L1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), LRZLBL(IDXLST), ISYLRZ(IDXLST), FRQLRZ(IDXLST),
     &    LORXLRZ(IDXLST)

*.....................................................................*
* 'L2' - second-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'L2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLL2(IDXLST,1), ISYL2(IDXLST,1), FRQL2(IDXLST,1),
     &       LBLL2(IDXLST,2), ISYL2(IDXLST,2), FRQL2(IDXLST,2)

*.....................................................................*
* 'L3' - third-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'L3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLL3(IDXLST,1), ISYL3(IDXLST,1), FRQL3(IDXLST,1),
     &       LBLL3(IDXLST,2), ISYL3(IDXLST,2), FRQL3(IDXLST,2),
     &       LBLL3(IDXLST,3), ISYL3(IDXLST,3), FRQL3(IDXLST,3)

*.....................................................................*
* 'L4' - third-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'L4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLL4(IDXLST,1), ISYL4(IDXLST,1), FRQL4(IDXLST,1),
     &       LBLL4(IDXLST,2), ISYL4(IDXLST,2), FRQL4(IDXLST,2),
     &       LBLL4(IDXLST,3), ISYL4(IDXLST,3), FRQL4(IDXLST,3),
     &       LBLL4(IDXLST,4), ISYL4(IDXLST,4), FRQL4(IDXLST,4)

*.....................................................................*
* 'LE' - Left eigenvector.
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'LE' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), ISYEXC(IDXLST), IMULTE(IDXLST), EIGVAL(IDXLST)

*.....................................................................*
* 'LC' - first-order left Cauchy vector.
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'LC' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), LBLLC1(IDXLST), ISYLC1(IDXLST), ILC1CAU(IDXLST)

*.....................................................................*
* 'Lx' - x-th. order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LIST(1:1) .EQ. 'L' ) THEN
         CALL QUIT ( LIST(1:2)//'-th. order lagrang. multip.'//
     &             ' not yet available' )

*.....................................................................*
* 'O1' - rhs vectors for first-order amplitude equations:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'O1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &     LBLO1(IDXLST), ISYO1(IDXLST), FRQO1(IDXLST), LORXO1(IDXLST)

*.....................................................................*
* 'O2' - rhs vectors for second-order amplitude equations:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'O2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLO2(IDXLST,1), ISYO2(IDXLST,1), FRQO2(IDXLST,1),
     &       LBLO2(IDXLST,2), ISYO2(IDXLST,2), FRQO2(IDXLST,2)

*.....................................................................*
* 'O3' - rhs vectors for third-order amplitude equations:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'O3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLO3(IDXLST,1), ISYO3(IDXLST,1), FRQO3(IDXLST,1),
     &       LBLO3(IDXLST,2), ISYO3(IDXLST,2), FRQO3(IDXLST,2),
     &       LBLO3(IDXLST,3), ISYO3(IDXLST,3), FRQO3(IDXLST,3)

*.....................................................................*
* 'O4' - rhs vectors for fourth-order amplitude equations:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'O4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLO4(IDXLST,1), ISYO4(IDXLST,1), FRQO4(IDXLST,1),
     &       LBLO4(IDXLST,2), ISYO4(IDXLST,2), FRQO4(IDXLST,2),
     &       LBLO4(IDXLST,3), ISYO4(IDXLST,3), FRQO4(IDXLST,3),
     &       LBLO4(IDXLST,4), ISYO4(IDXLST,4), FRQO4(IDXLST,4)

*.....................................................................*
* 'X1' - rhs vectors for first-order amplitude equations:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'X1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &     LBLX1(IDXLST), ISYX1(IDXLST), FRQX1(IDXLST), LORXX1(IDXLST)

*.....................................................................*
* 'X2' - second-order chi vectors:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'X2' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLX2(IDXLST,1), ISYX2(IDXLST,1), FRQX2(IDXLST,1),
     &       LBLX2(IDXLST,2), ISYX2(IDXLST,2), FRQX2(IDXLST,2)

*.....................................................................*
* 'X3' - third-order chi vectors:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'X3' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLX3(IDXLST,1), ISYX3(IDXLST,1), FRQX3(IDXLST,1),
     &       LBLX3(IDXLST,2), ISYX3(IDXLST,2), FRQX3(IDXLST,2),
     &       LBLX3(IDXLST,3), ISYX3(IDXLST,3), FRQX3(IDXLST,3)

*.....................................................................*
* 'X4' - fourth-order chi vectors:
*.....................................................................*
      ELSE IF ( LIST(1:2) .EQ. 'X4' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, LIST(1:2), 
     &       LBLX4(IDXLST,1), ISYX4(IDXLST,1), FRQX4(IDXLST,1),
     &       LBLX4(IDXLST,2), ISYX4(IDXLST,2), FRQX4(IDXLST,2),
     &       LBLX4(IDXLST,3), ISYX4(IDXLST,3), FRQX4(IDXLST,3),
     &       LBLX4(IDXLST,4), ISYX4(IDXLST,4), FRQX4(IDXLST,4)

*.....................................................................*
* 'M1' and 'FR' - first order response t amplitudes:
*.....................................................................*
      ELSE IF ((LIST(1:2).EQ.'M1').OR.(LIST(1:2).EQ.'FR')) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:2), ILRM(IDXLST), ISYLRM(IDXLST), FRQLRM(IDXLST)

*.....................................................................*
* 'PL1' - projected first-order response lagrangian multipliers:
*.....................................................................*
      ELSE IF ( LIST(1:3) .EQ. 'PL1' ) THEN
        WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:3), LBLPL1(IDXLST), ISYPL1(IDXLST), FRQPL1(IDXLST),
     &    LORXPL1(IDXLST)
*.....................................................................*
* unknown list:
*.....................................................................*
      ELSE
         CALL QUIT ('unknown list '//LIST(1:2)//' in '//SUBRNAME)
      END IF

*---------------------------------------------------------------------*
* write the vectors, close file & return:
*---------------------------------------------------------------------*
      CALL CC_WVRSP( ISYM, IMUL, IOPT, NWRITE, LUSAVE, VEC0, VEC1, VEC2,
     &               WORK(KORBITAL),WORK(KSINGLE),WORK(KDOUBLE))

      CALL GPCLOSE(LUSAVE,'KEEP')

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
      SUBROUTINE CC_WVRSP(ISYM,IMUL,IOPT,NWRITE,LUSAVE,
     *                    VEC0,VEC1,VEC2,VS0,VS1,VS2)
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
      INTEGER ISYM, IOPT

      DIMENSION VEC0(*),VEC1(*),VEC2(*),VS0(*),VS1(*),VS2(*)
C
      INTEGER IOS


* write vector components to file:
      IF (NWRITE.GE.1) THEN
         IF ( IAND(IOPT,4).EQ.4 ) THEN
           WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC0(I),I=1,2*NALLAI(ISYM))
         ELSE
           WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS0(I),I=1,2*NALLAI(ISYM))
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF
  
      IF (NWRITE.GE.2) THEN
         IF ( IAND(IOPT,1).EQ.1 ) THEN
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VEC1(I),I=1,NT1AM(ISYM))
         ELSE 
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (VS1(I),I=1,NT1AM(ISYM))
         END IF
      ELSE
            WRITE(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

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
C                 2.bit  (2): double part
C                 3.bit  (4): cphf part
C
C             for further input description see routines
C             CC_RDRSP and CC_WRRSP
C
C  Christof Haettig, April 1997
C
C---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit_dec.h>
#include <priunit.h>
#include <ccorb_dec.h>
#include <ccorb.h>
#include <ccsdsym_dec.h>
#include <ccsdsym.h>
#include <ccfro.h>
 
      INTEGER IDXLST, ISYM, IOPT, IOPTRW
      CHARACTER LIST*(*)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL VEC0(*), VEC1(*),VEC2(*)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION VEC0(*), VEC1(*),VEC2(*)
#endif
 
      CHARACTER MODFIL*(10)
      INTEGER KVEC0, KVEC1, KVEC2, KEND, LEND

*---------------------------------------------------------------------*
* allocate work space and read old vector from file:
*---------------------------------------------------------------------*
      IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
        WRITE (LUPRI,*) 'Symmetry out of range in CC_WARSP:'
        WRITE (LUPRI,*) 'ISYM:',ISYM
        WRITE (LUPRI,*) 'NSYM:',NSYM
        CALL QUIT('Symmetry out of range in CC_WARSP.')
      END IF

      KVEC0 = 1
      KVEC1 = KVEC0 + 2*NALLAI(ISYM)
      KVEC2 = KVEC1 + NT1AM(ISYM)
      KEND  = KVEC2 + NT2AM(ISYM)
      LEND  = LWORK - KEND
    
      IF (LEND .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_WARSP.')
      END IF

      IF (IOPT.GT.7) THEN
        CALL QUIT('IOPT > 7 not implemented in CC_WARSP.')
      END IF

      IOPTRW = 3
      CALL CC_RDRSP(LIST,IDXLST,ISYM,IOPTRW,MODFIL,
     &              WORK(KVEC1),WORK(KVEC2))
      
      IF ( IAND(IOPT,4).EQ.4 ) THEN
         CALL CC_RDHFRSP(LIST,IDXLST,ISYM,MODFIL,WORK(KVEC0))
         IOPTRW = IOPTRW + 4
      END IF

* check return value of read option:
      IF (IOPTRW .LT. IOPT) THEN
        WRITE (LUPRI,*) 'IOPT in CC_WARSP incompatible with vector '//
     &                  'on file.'
        CALL QUIT('IOPT in CC_WARSP incompatible with vector on file.')
      END IF

*---------------------------------------------------------------------*
* add new vector to old vector and write result back to file:
*---------------------------------------------------------------------*
      IF ( IAND(IOPT,1).EQ.1 ) THEN
        CALL DAXPY(NT1AM(ISYM),1.0d0,VEC1,1,WORK(KVEC1),1)
      END IF
      IF ( IAND(IOPT,2).EQ.2 ) THEN
        CALL DAXPY(NT2AM(ISYM),1.0d0,VEC2,1,WORK(KVEC2),1)
      END IF

      CALL CC_WRRSP(LIST,IDXLST,ISYM,IOPTRW,MODFIL,
     &              WORK(KVEC0),WORK(KVEC1),WORK(KVEC2),
     &              WORK(KEND),LEND)
  
      RETURN

      END 
*---------------------------------------------------------------------*
*                  END OF SUBROUTINE CC_WARSP
*---------------------------------------------------------------------*
