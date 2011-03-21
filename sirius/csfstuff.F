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
#ifdef UNDEF
===========================================================================
961125: GTSXST: added dimension JCLSX(*) to identify as vector
                for GTSXS2 call (found from SGI compiler message)
        LRFEXD: added dimension SIGN(*) to identify as vector
                for LRFEX2 call (found from SGI compiler message)
C  /* Deck gtsxst */
/* Comdeck comments */
C CSFDIM --- CISIZE
C        --- IWEYLF
C
C CSDTMT --- SPNCOM
C        --- UPDET
C        --- CSFDET -- MSSTRIN
C
C
C CNFORD ---- CONFGN
C        ---- CNTOST -- CNDET
C                    -- DETSTR
C                    -- IABNUM
C
C  CSDIAG
C
C  CSDTVC
===========================================================================
#endif
C  /* Deck csdtmt */
      SUBROUTINE CSDTMT(IDFTP,ICFTP,DTOC,WORK,IPRINT)
C
C     PURPOSE:
C
C     CONSTRUCT LIST OF PROTO TYPE DETERMINANTS, IN IDFTP
C     CONSTRUCT LIST OF PROTO TYPE CSF'S, IN ICFTP
C     CONSTRUCT MATRIX EXPANDING PROTO TYPE CSF'S IN TERMS OF
C     PROTOTYPE DETERMINANTS IN DTOC
C
C     SUBROUTINE AND(OR) EXTERNAL FUNCTION CALLS:
C
C     SPNCOM,CSFDET
C
c Jeppe Olsen
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION IDFTP(*),ICFTP(*),DTOC(*)
      DIMENSION WORK(*)
#include <spinfo.h>
C
C .. SET UP DETERMINANTS AND UPPER DETERMINANTS
C
      DO 20 ITP = 1, NTYP
        IOPEN = MINOP+ITP-1
        IF( IPRINT.GT.20 ) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,'(A)')
     &    '  *************************************************'
          WRITE(LUPRI,'(A,I3,A)')
     &    '    INFORMATION FOR TYPE WITH ',IOPEN,' OPEN ORBITALS '
          WRITE(LUPRI,'(A)')
     &    '  *************************************************'
          WRITE(LUPRI,*)
        END IF
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NDTFTP(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NCSFTP(ITP-1)
        END IF
C
        IF( IOPEN .NE. 0 ) THEN
C. PROTO TYPE DETERMINANTS AND BRANCHING DIAGRAM FOR
C  FOR PROTO TYPE DETS
          IF( MS2+1 .EQ. MULTS ) THEN
            IFLAG = 2
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' SPNCOM IS ENTERED FROM CSDTMT'
            CALL SPNCOM(IOPEN,MS2,NNDET,IDFTP(IDTBS),WORK(1) ,
     &                  ICFTP(ICSBS),IFLAG,IPRINT)
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RETURN TO CSDTMT FROM SPNCOM'
          ELSE
            IFLAG = 1
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' SPNCOM IS ENTERED FROM CSDTMT'
            CALL SPNCOM(IOPEN,MS2,NNDET,IDFTP(IDTBS),WORK(1) ,
     &                  ICFTP(ICSBS),IFLAG,IPRINT)
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RETURN TO CSDTMT FROM SPNCOM'
            IFLAG = 3
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' SPNCOM IS ENTERED FROM CSDTMT'
            CALL SPNCOM(IOPEN,MULTS-1,NNDET,IDFTP(IDTBS),WORK(1) ,
     &                  ICFTP(ICSBS),IFLAG,IPRINT)
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RETURN TO CSDTMT FROM SPNCOM'
      END IF
         END IF
   20 CONTINUE
C. matrix expressing csf's in terms of determinants
      DO 30 ITP = 1, NTYP
        IOPEN = MINOP+ITP-1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
          ICDCBS =1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NDTFTP(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NCSFTP(ITP-1)
          ICDCBS = ICDCBS + NDTFTP(ITP-1)*NCSFTP(ITP-1)
        END IF
        IF(NDTFTP(ITP)*NCSFTP(ITP).EQ.0) GOTO 30
        IF( IPRINT.GT.20 ) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,'(A)')
     &    '  *************************************************'
          WRITE(LUPRI,'(A,I3,A)')
     &    '    INFORMATION FOR TYPE WITH ',IOPEN,' OPEN ORBITALS '
          WRITE(LUPRI,'(A)')
     &    '  *************************************************'
          WRITE(LUPRI,*)
        END IF
        IF(IOPEN .EQ. 0 ) THEN
          DTOC(ICDCBS) = 1.0D0
        ELSE
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' CSFDET IS ENTERED FROM CSDTMT'
          CALL CSFDET(IOPEN,IDFTP(IDTBS),NDTFTP(ITP),
     &               ICFTP(ICSBS),NCSFTP(ITP),DTOC(ICDCBS),
     &               WORK(1),IPRINT)
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RETURN TO CSDTMT FROM CSFDET'
        END IF
   30 CONTINUE
C
      RETURN
      END
C  /* Deck spncom */
      SUBROUTINE SPNCOM(NOPEN,MS2,NDET,IABDET,IWORK,
     &                  IABUPP,IFLAG,IPRINT)
C
C
C     PURPOSE:
C
C     COMPUTE ALL COMBINATIONS OF NOPEN UNPAIRED ELECTRONS.REQUIRED
C     SPIN PROJECTION MS/2.
C
C     SUBROUTINE AND(OR) EXTERNAL FUNCTION CALLS:
C
C     ISETVC,IWRTMA,ICOPVE
C
C    IFLAG = 1 : ONLY DETERMINANTS ( IN IABDET )
C    IFLAG = 2 : DETERMINANTS AND UPPER DETS
C    IFLAG = 3 : ONLY UPPER DETS
C
c Jeppe Olsen
#include <implicit.h>
#include <priunit.h>
      INTEGER   ADD
      DIMENSION IABDET(NOPEN,*),IABUPP(NOPEN,*)
      DIMENSION IWORK(*)
C
      NDET=0
      NUPPER = 0
C
C DETERMINANTS ARE CONSIDERED AS BINARY NUMBERS,1=ALPHA,0=BETA
C
      MAX=2 ** NOPEN
      CALL ISETVC(IWORK,0,NOPEN)
C
C LOOP OVER ALL POSSIBLE BINARY NUMBERS
      DO 200 I=1,MAX
C
C.. 1 : NEXT BINARY NUMBER
        ADD=1
        J=0
  190   CONTINUE
        J=J+1
        IF(IWORK(J).EQ.1) THEN
          IWORK(J)=0
        ELSE
          IWORK(J)=1
          ADD=0
        END IF
        IF( ADD .EQ. 1 ) GOTO 190
C
C.. 2 :  CORRECT SPIN PROJECTION ?
        NALPHA=0
        DO 180 J=1,NOPEN
          NALPHA=NALPHA+IWORK(J)
  180   CONTINUE
C
        IF(2*NALPHA-NOPEN.EQ.MS2) THEN
          IF (IFLAG .LT. 3 ) THEN
            NDET=NDET+1
            CALL ICOPVE(IWORK,IABDET(1,NDET),NOPEN)
          END IF
          IF (IFLAG .GT. 1 ) THEN
C UPPER DET ?
            MS2L = 0
            LUPPER = 1
C
            DO 10 IEL = 1,NOPEN
              IF (IWORK(IEL).EQ.1) THEN
                 MS2L = MS2L + 1
              ELSE
                 MS2L = MS2L - 1
              END IF
              IF( MS2L .LT. 0 ) LUPPER = 0
   10       CONTINUE
            IF( LUPPER .EQ. 1 ) THEN
              NUPPER = NUPPER + 1
              CALL ICOPVE(IWORK,IABUPP(1,NUPPER),NOPEN)
            END IF
          END IF
        END  IF
C
  200 CONTINUE
C
C
      IF( IPRINT.GT.30 .AND. IFLAG .NE.3 ) THEN
         XMSD2=MS2
         XMSD2=XMSD2/2
         WRITE(LUPRI,1010) NOPEN,NDET,XMSD2
 1010    FORMAT(I6,' UNPAIRED ELECTRONS GIVE',I6,/,
     +'           DETERMINANTS WITH SPIN PROJECTION',F5.1,//,
     +'  DETERMINANTS : '/
     +'  ============== ')
         DO 20 J=1,NDET
           WRITE(LUPRI,1020) J,(IABDET(K,J),K=1,NOPEN)
  20     CONTINUE
 1020    FORMAT(/I6,2X,30I2,/,(8X,30I2))
      END IF
C
      IF( IFLAG.GT.1 .AND. IPRINT.GT.20 ) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)') ' UPPER DETERMINANTS '
         WRITE(LUPRI,'(A)') ' ================== '
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' NUPPER ', NUPPER
         DO 22 J=1,NUPPER
           WRITE(LUPRI,1020) J,(IABUPP(K,J),K=1,NOPEN)
  22     CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck csfdet */
      SUBROUTINE CSFDET(NOPEN,IDET,NDET,ICSF,NCSF,CDC,WORK,IPRINT)
C
C
C     PURPOSE:
C
C     EXPAND CSF'S IN TERMS OF DETERMINANTS WITH
C     THE USE OF THE GRAEBENSTETTER METHOD ( I.J.Q.C.10,P142(1976) )
C
C     CALLING PARAMETERS:
C
C     ON INPUT:  NOPEN: NUMBER OF OPEN ORBITALS
C                IDET : OCCUPATION OF DETERMINANTS
C                NDET : NUMBER OF DETERMINANTS
C                IPRINT: FLAG TO INDICATE PRINTING OF INFO
C                ICSF : INTERMEDIATE SPIN COUPLINGS OF
C                        CSF'S IN BRANCHING DIAGRAM
C                WORK  : SCRATCH ARRAY
C     ON OUTPUT:  CDC :  NDET X NCSF MATRIX
C                        GIVING EXPANSION FROM SD'S TO CSF,S
C            CSF BASIS = SD BASIS * CDC
C
C     SUBROUTINE AND(OR) EXTERNAL FUNCTION CALLS:
C
C     WRTMA,MSSTRN
c
c Jeppe Olsen
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IDET(NOPEN,NDET),ICSF(NOPEN,NCSF)
      DIMENSION CDC(NDET,NCSF)
      DIMENSION WORK(*)
C
      KLFREE = 1
      KLMDET = KLFREE
      KLFREE = KLMDET + NDET * NOPEN
      KLSCSF =  KLFREE
      KLFREE = KLSCSF + NOPEN
C
C.. OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' MSSTRN IS ENTERED FROM CSFDET'
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' MSSTRN IS CALLED',NDET,'TIMES'
      DO 10 JDET = 1, NDET
        CALL MSSTRN
     &  (IDET(1,JDET),WORK(KLMDET+(JDET-1)*NOPEN),NOPEN,IPRINT)
   10 CONTINUE
C
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' MSSTRN IS ENTERED FROM CSFDET'
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' MSSTRN IS CALLED',NCSF,'TIMES'
      DO 1000 JCSF = 1, NCSF
       IF( IPRINT.GT.30 ) WRITE(LUPRI,*) ' ....OUTPUT FOR CSF ',JCSF
C
C OBTAIN INTERMEDIATE COUPLINGS FOR CSF
      CALL MSSTRN(ICSF(1,JCSF),WORK(KLSCSF),NOPEN,IPRINT)
C
      DO 900 JDET = 1, NDET
C EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
      COEF = 1.0D0
      SIGN = 1.0D0
      JDADD = (JDET-1)*NOPEN
#include <scalar.h>
      DO 700 IOPEN = 1, NOPEN
C
C + + CASE
        IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.1) THEN
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.0) THEN
C + - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.1) THEN
C - + CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
          SIGN  = - SIGN
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.0) THEN
C - - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
        END IF
  700 CONTINUE
       CDC(JDET,JCSF) = SIGN * SQRT(COEF)
  900 CONTINUE
 1000 CONTINUE
C
      IF( IPRINT.GT.30 ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,'(A,2I2)')
     &  '  THE CDC ARRAY FOR  NOPEN ',NOPEN
        WRITE(LUPRI,*)
        CALL WRTMAT(CDC,NDET,NCSF,NDET,NCSF,0)
       END IF
C
      RETURN
      END
C  /* Deck msstrn */
      SUBROUTINE MSSTRN(INSTRN,UTSTRN,NOPEN,IPRINT)
C
C
C     PURPOSE:
C
C      A STRING IS GIVEN IN FORM A SEQUENCE OF ZEROES
C      AND ONE ' S
C
C      REINTERPRET THIS AS :
C
C      1 : THE INPUT STRING IS A DETERMINANT AND THE
C         1'S INDICATE ALPHA ELECTRONS AND THE
C         0'S INDICATE BETA ELECTRONS .
C         UTSTRN IS THE MS-VALUES AT EACH VERTEX
C
C     2 : THE INPUT STRING IS A CSF GIVEN IN A
C         BRANCHING DIAGRAM, WHERE
C         1'S INDICATE UPWARDS SPIN COUPLEING
C         WHILE THE 0'S INDICATES DOWNWARDS SPIN COUPLING ,
C         REEXPRESS THIS AS S VALUES OF ALL COUPLINGS
C
C
C     Externals
C
C     WRTMA,IWRTMA
C
C Jeppe Olsen
c
C THE TWO PROCEDURES ARE IDENTICAL .

#include <implicit.h>
#include <priunit.h>
      DIMENSION INSTRN(NOPEN),UTSTRN(NOPEN)
C
      UTSTRN(1) = DFLOAT(INSTRN(1)) - 0.5D0
      DO 10 IOPEN = 2, NOPEN
        UTSTRN(IOPEN) = UTSTRN(IOPEN-1) +DFLOAT(INSTRN(IOPEN))-0.5D0
   10 CONTINUE
C
      IF( IPRINT.GT.30 ) THEN
         WRITE(LUPRI,*) ' ... OUTPUT FROM MSSTRN '
         WRITE(LUPRI,*) ' INSTRN AND UTSTRN'
         CALL IWRTMA(INSTRN,1,NOPEN,1,NOPEN)
         CALL WRTMAT(UTSTRN,1,NOPEN,1,NOPEN,0)
       END IF
C
      RETURN
      END
C  /* Deck cnford */
      SUBROUTINE CNFORD(ICTSDT,ICONF,IWORK,ORBSYM,
     &           SYMPRO,IREFSM,NORB,IDFTP,XNDXCI,NCNFTP,NEL,NTEST)
C
C Generate configurations in ICONF
C
C Generate determinants in configuration order and obtain
C sign array for switching between the two formats .
C
C Jeppe Olsen January 1989
C
#include <implicit.h>
      DIMENSION ICTSDT(*)
      INTEGER ORBSYM(*),SYMPRO(8,8)
      INTEGER IDFTP(*)
      DIMENSION XNDXCI(*)
      DIMENSION IWORK(*)
      DIMENSION NCNFTP(*)
      DIMENSION ICONF(*)
C NOTE : NCNFTP IS COLUMN FOR SYMMETRY GIVEN , NOT COMPLETE MATRIX.
C Dim of IWORK : MAX(2*NTYP+3*NORB,(MXDT+2)*NEL),
C where MXDT is the largest number of prototype determinants of
C a given type.
C
C.. 1 : Construct list of configurations,offset for each configuration
C       and type for each configuration
C     CALL CONFGN(NORB,IREFSM,ORBSYM,SYMPRO,ICONF,IWORK,
C    & NEL,NTEST)
C     CONFG2 (NORB,IREFSM,ORBSYM,SYMPRO,ICONF,
C    &        NEL,IIOC,IIOP,IICL,IPRINT)
      KL1 = 1
      KL2 = KL1 + NORB
      KL3 = KL2 + NORB
      CALL CONFG2(NORB,IREFSM,ORBSYM,SYMPRO,ICONF,NEL,
     &            IWORK(KL1),IWORK(KL2),IWORK(KL3),NTEST)
C
C       Obtain determinants for each configuration and determine
C.. 2 : the corresponding adress and phaseshift to reform into
C       string form and ordering.
C
c. Offset arrays for SD's in string order
      CALL CIOFF(IREFSM,1,XNDXCI,NTEST)
c
      CALL CNTOST(ICONF,ICTSDT,IWORK,IDFTP(1),IREFSM,
     &            XNDXCI,NORB,NEL,ORBSYM,NTEST)
C          CNTOST(ICONF,ICTSDT,IWORK,IPRODT,IREFSM,
C    &            XNDXCI,NORB,NEL,ORBSYM,NTEST)
      RETURN
      END
C  /* Deck csdiag */
      SUBROUTINE CSDIAG(CSFDIA,DETDIA,NCNFTP,NTYP,
     &                   ICTSDT,NDTFTP,NCSFTP,IFLAG,NTEST )
C
C.. obtain average CI diagonal elements and store in
C   CSFDIA as CSF diagonal
C
C if IFLAG .NE. 0 determinant diagonal is changed to
C average diagonal values .
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION CSFDIA(*),DETDIA(*)
      DIMENSION NCNFTP(NTYP),NDTFTP(NTYP),NCSFTP(NTYP)
      DIMENSION ICTSDT(*)
C
      ICSOFF = 1
      IDTOFF = 1
      JCNABS = 0
      DO 100 ITYP = 1, NTYP
        IDET = NDTFTP(ITYP)
        ICSF = NCSFTP(ITYP)
        ICNF = NCNFTP(ITYP)
#if defined (SYS_NEC)
	IF (ICNF .LE. 0) GOTO 100
#endif
        DO 80 JCNF = 1, ICNF
          JCNABS = JCNABS + 1
          EAVER = 0.0D0
          DO 70 JDET = 1, IDET
            EAVER = EAVER + DETDIA( ABS(ICTSDT(IDTOFF-1+JDET) ))
   70     CONTINUE
          IF( IDET .NE. 0 )EAVER = EAVER/IDET
          IF( IFLAG .NE. 0 ) THEN
C.          write  average values in determinant diagonal
            DO 75 JDET = 1, IDET
              DETDIA( ABS(ICTSDT(IDTOFF-1+JDET)) ) = EAVER
   75       CONTINUE
          END IF
          CALL SETVEC(CSFDIA(ICSOFF),EAVER,ICSF)
          ICSOFF = ICSOFF + ICSF
          IDTOFF = IDTOFF + IDET
   80   CONTINUE
  100 CONTINUE
C
      IF( NTEST .GT.20 ) THEN
        WRITE(LUPRI,*) 'CIDIAGONAL IN DET AND CSF BASIS '
        NCSTOT = ICSOFF-1
        NDTTOT = IDTOFF-1
        CALL WRTMAT(DETDIA,1,NDTTOT,1,NDTTOT,0)
        CALL WRTMAT(CSFDIA,1,NCSTOT,1,NCSTOT,0)
      END IF
      RETURN
      END
C  /* Deck csdtvc */
      SUBROUTINE CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,
     &                  IREFSM,ICOPY,NTEST)
C
C IWAY = 1 : CSF to DETERMINANT TRANSFORMATION
C IWAY = 2 : DETERMINANT TO CSF TRANSFORMATION
C
C ICOPY .NE. 0 : COPY OUTPUT INTO INPUT
C                SO INPUT BECOMES OUTPUT WHILE
C                OUTPUT REMAINS OUTPUT ( FOR THE MOMENT )
C
C MOTECC-90: This module, CSDTVC, and the algorithms used
C            are described in Chapter 8 Setion D.1 of MOTECC-90
C            "RAS-CI Expansions in a CSF basis"
C
#include <implicit.h>
#include <priunit.h>
#include <spinfo.h>
#include <ciinfo.h>
      DIMENSION CSFVEC(*),DETVEC(*)
      DIMENSION DTOCMT(*),ICTSDT(*)
C
      IF( NTEST .GE. 20 ) WRITE(LUPRI,*) ' >>> CSDTVC SPEAKING <<<'
      NDET = NDTASM(IREFSM)
      NCSF = NCSASM(IREFSM)
C
      IF(NTEST.GE.20) WRITE(LUPRI,*) ' NDET NCSF IREFSM',NDET,
     &                                 NCSF,IREFSM
      IF(IWAY .EQ. 1 ) THEN
C
c--renzo important for conversion cfg->det
         IF(NTEST .GE. 40) THEN
            WRITE(LUPRI,*) ' INPUT VECTOR IN CSF BASIS '
            CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF,0)
         END IF
C.. CSF to DET transformation
C
C. Multiply with  expansion matrix
        DO 100 ITYP = 1,NTYP
          IDET = NDTFTP(ITYP)
          ICSF = NCSFTP(ITYP)
          ICNF = NCNFTP(ITYP,IREFSM)
          IF(ITYP .EQ. 1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNFTP(ITYP-1,IREFSM)*NCSFTP(ITYP-1)
            IOFFDT = IOFFDT+NCNFTP(ITYP-1,IREFSM)*NDTFTP(ITYP-1)
            IOFFCD = IOFFCD + NDTFTP(ITYP-1)*NCSFTP(ITYP-1)
          END IF
          CALL MATML4(DETVEC(IOFFDT),DTOCMT(IOFFCD),CSFVEC(IOFFCS),
     &                IDET,ICNF,IDET,ICSF,ICSF,ICNF,0)
  100   CONTINUE
C. Change from csf to string ordering with sign changes
        CALL SCACSF(NDET,CSFVEC,DETVEC,ICTSDT)
        CALL COPVEC(CSFVEC,DETVEC,NDET)
C       CALL SCACSF(NDET,ADET,BCSF,IORD)
C       IF(ICOPY.NE.0) CALL COPVEC(DETVEC,CSFVEC,NDET)
Chj-890501: vector is currently already in both DETVEC and CSFVEC
c--renzo 
        IF(NTEST .GE. 40) THEN
           WRITE(LUPRI,*) ' OUTPUT VECTOR IN DET BASIS '
           CALL WRTMAT(DETVEC,1,NDET,1,NDET,0)
        END IF
      ELSE
         IF( NTEST .GE. 40 ) THEN
            WRITE(LUPRI,*) ' INPUT VECTOR ON DETERMINANT BASIS '
            CALL WRTMAT(DETVEC,1,NDET,1,NDET,0)
         END IF
C
C.. Determinant to csf transformation
C
C. To CSF ordering
        CALL GATCSF(NDET,DETVEC,CSFVEC,ICTSDT)
        CALL COPVEC(CSFVEC,DETVEC,NDET)
C       CALL GATCSF(NDET,ADET,BCSF,IORD)
C. Multiply with CIND expansion matrix
        DO 200 ITYP = 1,NTYP
          IDET = NDTFTP(ITYP)
          ICSF = NCSFTP(ITYP)
          ICNF = NCNFTP(ITYP,IREFSM)
          IF(ITYP .EQ. 1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNFTP(ITYP-1,IREFSM)*NCSFTP(ITYP-1)
            IOFFDT = IOFFDT+NCNFTP(ITYP-1,IREFSM)*NDTFTP(ITYP-1)
            IOFFCD = IOFFCD + NDTFTP(ITYP-1)*NCSFTP(ITYP-1)
          END IF
          CALL MATML4(CSFVEC(IOFFCS),DTOCMT(IOFFCD),DETVEC(IOFFDT),
     &                ICSF,ICNF,IDET,ICSF,IDET,ICNF,1)
  200   CONTINUE
        IF( ICOPY .NE. 0 ) CALL COPVEC(CSFVEC,DETVEC,NCSF)
        IF( NTEST .GE. 40 ) THEN
          WRITE(LUPRI,*) ' OUTPUT VECTOR IN CSF BASIS '
          CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF,0)
        END IF
       END IF
C
      RETURN
      END
C  /* Deck iabnum */
      FUNCTION IABNUM(IASTR,IBSTR,XNDXCI,ORBSYM)
C
C Encapsulation routine for IABNUS
C
#include <implicit.h>
#include <mxpdim.h>
#include <detbas.h>
#include <ciinfo.h>
#include <strnum.h>
#include <inforb.h>
      DIMENSION IASTR(*),IBSTR(*),XNDXCI(*)
      INTEGER   ORBSYM(*)
C
      IABNUM = IABNUS(IASTR,IBSTR,NAEL,NBEL,XNDXCI(KCOFF),MULD2H,
     &                XNDXCI(KIPNSA),XNDXCI(KIPNSB),XNDXCI(KZA),
     &                XNDXCI(KZB),ORBSYM,
     &                XNDXCI(KSTBAA),XNDXCI(KSTBAB),XNDXCI(KSTASA),
     &                XNDXCI(KTPFSA),XNDXCI(KTPFSB),
     &                XNDXCI(KNSSOA),XNDXCI(KNSSOB),XNDXCI(KISSOA),
     &                XNDXCI(KISSOB),XNDXCI(KIOCOC),XNDXCI(KICOOS),
     &                NOCTPA,NOCTPB,NASHT)
C          IABNUS(IASTR,IBSTR,NAEL,NBEL,IOFF,SYMPRO,IAORD,
C     &           IBORD,ZA,ZB,ORBSYM,ISTBAA,ISTBAB,NSTASA,
C     &           ITPFSA,ITPFSB,NSSOA,NSSOB,ISSOA,ISSOB,
C     &           IOCOC,ICOOS,NORB)
C
      RETURN
      END
C  /* Deck cntost */
      SUBROUTINE CNTOST(ICONF,ICTSDT,IWORK,IPRODT,IREFSM,
     &                  XNDXCI,NORB,NEL,ORBSYM,NTEST)
C
C Obtain pointer ICTSDT(I).
C ABS(ICTSDT(I) giving adress of determinant I in
C STRING ordering for determinant I in CSF ordering.
C Going between the two formats can involve a sign change . this is
C stored in sign of ICTSDT(I)
C The sign is thus to be multiplied with determinant
C vector ordered in CSF ordering.
C
#include <implicit.h>
#include <priunit.h>
#include <mxpdim.h>
#include <strnum.h>
#include <spinfo.h>
#include <ciinfo.h>
      DIMENSION ICONF(*),IPRODT(*), ICTSDT(*), IWORK(*),XNDXCI(*)
      INTEGER   ORBSYM(*)
C
C IWORK should at least be of length (MXDT+2)*NEL,
C where MXDT is the largest number of prototype determinants occuring
C in a single block.
C
C.. Local memory
C
      IF( NTEST .GE. 10 ) WRITE(LUPRI,*) ' >>> OUTPUT FROM CNTOST <<< '
      KLFREE = 1
C Largest number of dets for a given type
      MXDT = 0
      DO 10 ITYP = 1, NTYP
        MXDT   = MAX(MXDT,NDTFTP(ITYP) )
   10 CONTINUE
      IF(NTEST.GE.10) WRITE(LUPRI,*) ' MXDT IS... ', MXDT
      KLDTBL = KLFREE
      KLFREE = KLDTBL + MXDT*NEL
C Space for an alpha and an betastring
      KLIA = KLFREE
      KLFREE = KLIA + NAEL
      KLIB = KLFREE
      KLFREE = KLFREE + NBEL
C
C.. Loop over configurations and generate determinants in compact form
C
      ICNF = 0
      IDT  = 0
      JDTABS = 0
      DO 1000 ITYP = 1, NTYP
        IDET = NDTFTP(ITYP)
        IOPEN = ITYP + MINOP - 1
        ICL = (NEL - IOPEN) / 2
        IOCC = IOPEN + ICL
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
        ELSE
          ICNBS0 = ICNBS0 + NCNFTP(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
        END IF
        IF( NTEST .GE. 50 ) WRITE(LUPRI,*) ' ITYP ICNBS0 ',ITYP,ICNBS0
        IF( NTEST .GE. 50 ) WRITE(LUPRI,*) ' IDET IOPEN ICL ',IDET,
     &                                       IOPEN,ICL
C Base for prototype determinants
        IF( ITYP .EQ. 1 ) THEN
          IPBAS = 1
        ELSE
          IPBAS = IPBAS + NDTFTP(ITYP-1)*(IOPEN-1)
        END IF
C Determinants for this configuration
        DO 900  IC = 1, NCNFTP(ITYP,IREFSM)
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
          IF( NTEST .GE. 10 ) WRITE(LUPRI,*) ' IC ICNF ICNBS',IC,ICNF,
     &                                       ICNBS
C              CNDET(ICONF,IPDET, NDET,NEL,NORB,NOP,NCL,IDET)
          CALL CNDET(ICONF(ICNBS),IPRODT(IPBAS),IDET,
     &               NEL,IOCC,IOPEN,ICL,
     &               IWORK(KLDTBL),NTEST )
C Separate determinants into strings and determine string number .
          DO 800 JDET = 1,IDET
            JDTABS = JDTABS + 1
            CALL DETSTR(IWORK(KLDTBL+(JDET-1)*NEL),IWORK(KLIA),
     &           IWORK(KLIB),NEL,NAEL,NBEL,NORB,ISIGN,IWORK(KLFREE) )
C Find number of this determinant in string ordering
            NEWDET = IABNUM(IWORK(KLIA),IWORK(KLIB),XNDXCI,ORBSYM)
C Save sign for configuration to string order switch
            ICTSDT(JDTABS) = SIGN(NEWDET,ISIGN)
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
C
      IF( NTEST .GE. 10) THEN
        WRITE(LUPRI,*) ' ICTSDT ARRAY '
        CALL IWRTMA(ICTSDT,1,JDTABS,1,JDTABS)
      END IF
C
      RETURN
      END
C  /* Deck iabnus */
      FUNCTION IABNUS(IASTR,IBSTR,NAEL,NBEL,IOFF,SYMPRO,IAORD,
     &                IBORD,ZA,ZB,ORBSYM,ISTBAA,ISTBAB,NSTASA,
     &                ITPFSA,ITPFSB,NSSOA,NSSOB,ISSOA,ISSOB,
     &                IOCOC,ICOOS,NOCTPA,NOCTPB,NORB)
C
C
C A DETERMINANT IS GIVEN BY STRINGS IASTR,IBSTR .
C FIND NUMBER OF THIS DETERMINANT
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IOFF(*)
      DIMENSION IASTR(NAEL),IBSTR(NBEL)
      DIMENSION IAORD(*),IBORD(*)
      INTEGER SYMPRO(8,8),ORBSYM(*)
      INTEGER ISTBAA(*),ISTBAB(*)
      INTEGER ZA(*),ZB(*)
      INTEGER NSTASA(*)
      INTEGER ITPFSA(*),ITPFSB(*)
      INTEGER NSSOA(NOCTPA,*),NSSOB(NOCTPB,*)
      INTEGER ISSOA(NOCTPA,*),ISSOB(NOCTPB,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER ICOOS(NOCTPB,NOCTPA,*)
C
C JEPPE OLSEN NOVEMBER 1988
C
C.. NUMBER OF ALPHA AND BETA STRING
      NTEST = 0
      IF( NTEST .GT. 10 ) THEN
       WRITE(LUPRI,*) ' >>> IABNUS SPEAKING <<< '
       WRITE(LUPRI,*) ' NOCTPA,NOCTPB ', NOCTPA,NOCTPB
C      WRITE(LUPRI,*) ' ALPHA AND BETA STRING '
C      CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
C      CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
      END IF
C
C     FUNCTION IZNUM(IOCC,NEL,Z,MAXEL,MXCORB,NEWORD)
      IANUM = IZNUM(IASTR,NAEL,ZA,NAEL,NORB,IAORD,NTEST)
      IBNUM = IZNUM(IBSTR,NBEL,ZB,NBEL,NORB,IBORD,NTEST)
      IF( NTEST .GE. 10 ) WRITE(LUPRI,*) ' IANUM AND IBNUM ',IANUM,IBNUM
C
C.. SYMMETRY OF ALPHA AND BETA STRING
C
C     FUNCTION MSYMP2(ISTR,NELMNT,ISYM)
      IASYM = MSYMP2(IASTR,NAEL,ORBSYM,SYMPRO)
      IBSYM = MSYMP2(IBSTR,NBEL,ORBSYM,SYMPRO)
      IF( NTEST .GE.10) WRITE(LUPRI,*) ' IASYM IBSYM ',IASYM,IBSYM
C
        IATP = ITPFSA(IANUM)
        IBTP = ITPFSB(IBNUM)
        IF(NTEST.GE.10) WRITE(LUPRI,*) ' IATP,IBTP ', IATP,IBTP
        IAREL = IANUM - ISSOA(IATP,IASYM)+1
        IBREL = IBNUM - ISSOB(IBTP,IBSYM)+1
        IF(NTEST .GE.10) WRITE(LUPRI,*) ' IAREL IBREL ', IAREL,IBREL
        IABNUS = ICOOS(IBTP,IATP,IASYM) + (IBREL-1)*NSSOA(IATP,IASYM)
     &           + IAREL - 1
C
      IF ( NTEST .GE. 2 ) THEN
         WRITE(LUPRI,*) ' ALPHA AND BETA STRING '
         CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
         CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
         WRITE(LUPRI,*) ' CORRESPONDING DETERMINANT NUMBER ', IABNUS
      END IF
C
      RETURN
      END
C  /* Deck ordstr */
      SUBROUTINE ORDSTR(IINST,IOUTST,NELMNT,ISIGN )
C
C ORDER A STRING OF INTEGERS TO ASCENDING ORDER
C
C IINST : INPUT STRING IS IINST
C IOUTST : OUTPUT STRING IS IOUTST
C NELMNT : NUMBER OF INTEGERS IN STRING
C ISIGN :  SIGN OF PERMUTATION : + 1 : EVEN PERMUTATIONN
C                                - 1 : ODD  PERMUTATION
C
C THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
C ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IINST(NELMNT),IOUTST(NELMNT)
C
      CALL ICOPVE(IINST,IOUTST,NELMNT)
      ISIGN = 1
C
C       BEGIN TO ORDER
C
        JOE = 1
  10    I = JOE
  20    CONTINUE
        IF(I.EQ.NELMNT) GO TO 50
        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
        JOE = I + 1
  30    ISWAP = IOUTST(I)
        ISIGN = - ISIGN
        IOUTST(I) = IOUTST(I+1)
        IOUTST(I+1) = ISWAP
        IF(I.EQ.1) GO TO 10
        I = I - 1
        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
        GO TO 10
 40     I = I + 1
      GO TO 20
C
C     END ORDER
C
 50   CONTINUE
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*)  ' INPUT STRING ORDERED STRING ISIGN '
        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
        WRITE(LUPRI,*) ' ISIGN : ', ISIGN
      END IF
C
      RETURN
      END
C  /* Deck cndet */
      SUBROUTINE CNDET(ICONF,IPDET,NDET,NEL,NORB,NOP,NCL,
     &                 IDET,NTEST)
C
C A configuration ICONF in compressed form and a set of
C prototype determinants ,IPDET, are given .
C
C Construct the corresponding determinants in contracted  form .
C
C JEPPE OLSEN , NOVEMBER 1988
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION ICONF(NORB), IPDET(*       ), IDET(NEL,NDET)
C     DIMENSION ICONF(NORB), IPDET(NOP,NDET), IDET(NEL,NDET)
C     890417-hj: new IPDET definition to avoid error when NOP.eq.0
C
C POSITIVE NUMBER  : ALPHA ORBITAL
C NEGATIVE NUMBER  : BETA  ORBITAL
C
      IF( NTEST .GT. 20 ) THEN
        WRITE(LUPRI,*) ' DOUBLY OCCUPIED ORBITALS '
        CALL IWRTMA(ICONF,1,NCL,1,NCL)
        WRITE(LUPRI,*) ' OPEN ORBITALS '
        IF (NOP .GT. 0) THEN
          CALL IWRTMA(ICONF(1+NCL),1,NOP,1,NOP)
        ELSE
          WRITE (LUPRI,'(/A/)') ' none'
        END IF
        WRITE(LUPRI,*) ' NEL IN CNDET ', NEL
      END IF
C
C.. 1 DOUBLY OCCUPIED ORBITALS ARE PLACED FIRST
C
      DO 100 ICL = 1, NCL
        IBASE = 2 * (ICL-1)
        DO 90 JDET = 1, NDET
          IDET(IBASE+1,JDET) =  ICONF(ICL)
          IDET(IBASE+2,JDET) = -ICONF(ICL)
   90   CONTINUE
  100 CONTINUE
C
C..2  SINGLY OCCUPIED ORBITALS
C
      IADD = 2*NCL
      DO 200 JDET = 1, NDET
        DO 190 IOP = 1, NOP
C         IF( IPDET(IOP,JDET) .EQ. 1 )
          IF( IPDET(IOP+(JDET-1)*NOP) .EQ. 1 )
     &      IDET(IADD+IOP,JDET) =   ICONF(NCL +IOP)
C
C         IF( IPDET(IOP,JDET) .EQ. 0 )
          IF( IPDET(IOP+(JDET-1)*NOP) .EQ. 0 )
     &      IDET(IADD+IOP,JDET) = - ICONF(NCL +IOP)
C
  190   CONTINUE
  200 CONTINUE
C
      IF( NTEST .GE. 20) THEN
       WRITE(LUPRI,*) ' CONFIGURATION FROM DETCON '
       CALL IWRTMA(ICONF,1,NORB,1,NORB)
       IF (NOP.GT.0) THEN
         WRITE(LUPRI,* ) ' PROTO TYPE DETERMINANTS '
         CALL IWRTMA(IPDET,NOP,NDET,NOP,NDET)
       END IF
       WRITE(LUPRI,*) ' CORRESPONDING DETERMINANTS '
       CALL IWRTMA(IDET,NEL,NDET,NEL,NDET)
      END IF
C
      RETURN
      END
C  /* Deck detstr */
      SUBROUTINE DETSTR(IDET,IASTR,IBSTR,NEL,NAEL,NBEL,NORB,
     &ISIGN,IWORK)
C
C A DETERMINANT,IDET,IS GIVEN AS A SET OF OCCUPIED SPIN ORBITALS,
C POSITIVE NUMBER INDICATES ALPHA ORBITAL AND NEGATIVE NUMBER
C INDICATES BETA ORBITAL .
C
C FIND CORRESPONDING ALPHA STRING AND BETA STRING ,
C AND DETERMINE SIGN NEEDED TO CHANGE DETERMINANT
C INTO PRODUCT OF ORDERED ALPHA STRING AND
C BETA STRING
C
C JEPPE OLSEN NOVEMBER 1988
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION IDET(NEL)
      DIMENSION IASTR(NAEL),IBSTR(NBEL)
      DIMENSION IWORK(*)
C required length of IWORK : NEL
C
      NTEST = 0
C
C FIRST REORDER SPIN ORBITALS IN ASCENDING SEQUENCE
C THIS WILL AUTOMATICALLY SPLIT ALPHA AND BETASTRING
C
C     SUBROUTINE ORDSTR(IINST,IOUTST,NELMNT,ISIGN )
      CALL ORDSTR(IDET,IWORK,NEL,ISIGN)
C
C ALPHA STRING IS LAST NAEL ORBITALS
      CALL ICOPVE(IWORK(NBEL+1),IASTR,NAEL)
C
C BETA  STRING MUST BE COMPLETELY TURNED AROUND
      DO 10 IBEL = 1, NBEL
        IBSTR(IBEL) = -IWORK(NBEL+1-IBEL)
   10 CONTINUE
C SIGN CHANGE FOR SWITCH OF BETA ORBITALS
      NNBEL = (NBEL*(NBEL+1))/2
      ISIGN = ISIGN * (-1) ** NNBEL
C
      IF( NTEST .GE. 20) THEN
        WRITE(LUPRI,*) ' INPUT DETERMINANT '
        CALL IWRTMA(IDET,1,NEL,1,NEL)
        WRITE(LUPRI,*) ' CORRESPONDING ALPHA STRING '
        CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
        WRITE(LUPRI,*) ' CORRESPONDING BETA STRING '
        CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
        WRITE(LUPRI,*) ' ISIGN FOR SWITCH ', ISIGN
      END IF
C
      RETURN
      END
C  /* Deck confg2 */
      SUBROUTINE CONFG2 (NORB,IREFSM,ORBSYM,SYMPRO,ICONF,
     &                   NEL,IIOC,IIOP,IICL,IPRINT)
C                CONFG2 (NORB,IREFSM,ORBSYM,SYMPRO,ICONF,
C    &                   NEL,IWORK(KL1),IWORK(KL2),IWORK(KL3),IPRINT)
COLD  SUBROUTINE CONFGN (NORB,IREFSM,ORBSYM,SYMPRO,ICONF,
COLD &                          IWORK,NEL,NTEST)
C
C     PURPOSE:
C
C     TURBO CONFIGURATION GENERATOR
C     GENERATE ARRAY,ICONF,GIVING OCCUPATION OF EACH CONFIGURATION
C     FOR CI SPACE OF REFERENCE SYMMETRY IREFSM.
C     ICONF IS ORDERED SO ALL CONFIGURATIUNS OF THE SAME TYPE ARE
C     CONSECUTIVELY STORED .
C     ICONF IS WRITTEN SO CLOSED ORBITALS ARE GIVEN FIRST AND THEN SINGL
C     OCCUPIED ORBITALS
C
C     SUBROUTINE CALLS:
C
C     ISETVC,IWRTMA
C
C
c
cJeppe Olsen , 1989
c
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION ICONF(*),IIOC(*),IICL(*),IIOP(*)
      INTEGER SYMPRO(8,8)
      INTEGER ORBSYM(NORB)
C
#include <spinfo.h>
#include <ciinfo.h>
C
C.. Loop over types of configurations
C
      JCONF =  0
      ICFREE = 1
c Min number of doubly occupied orbitals in RAS 1
      NORB1 = IORB1L-IORB1F+1
C corrected def of MINCL1
      MINCL1 = MAX(0,NEL1MN-NORB1)
      IF( IPRINT.GT.30 ) WRITE(LUPRI,*)
     &' Min number of doubly occupied orbitals in RAS 1',MINCL1
      DO 5000 NOP = MINOP,MAXOP,2
        ITYPE = NOP-MINOP+1
        NCL = (NEL-NOP)/2
        IF( IPRINT.GT.30 )
     &  WRITE(LUPRI,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
        IF(IPRINT.GT.30)
     &  WRITE(LUPRI,*) ' CONFG2 : NEL1MN NEL1MX NEL3MN NEL3MX ',
     &                        NEL1MN,NEL1MX,NEL3MN,NEL3MX
C
C. first combination of double occupied orbitals
        CALL ISETVC(IIOC,0,NORB)
        DO 10 ICL = 1, NCL
          IICL(ICL) = ICL
          IIOC(ICL) = 2
   10   CONTINUE
        IFRSTC = 1
C.. Loop over double occupied orbital configurations
 2000   CONTINUE
C
C. next double occupied configuration
          IF (.NOT. ( IFRSTC .EQ. 0 .AND. NCL .NE. 0 ) ) GO TO 810
C
          DO 50 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 1 ) IIOC(IORB) = 0
   50     CONTINUE
C
          IPLACE = 0
  800     IPLACE = IPLACE + 1

          IPRORB = IICL(IPLACE)
          IIOC(IPRORB) = 0
          NEWORB = IPRORB+1
          IF((IPLACE .LT. NCL .AND. NEWORB .LT. IICL(IPLACE+1))
     &      .OR.
     &      IPLACE .EQ. NCL .AND.  NEWORB.LE. NORB ) THEN
C
            IICL(IPLACE) = NEWORB
            IIOC(NEWORB) = 2
          ELSE IF
     &    (.NOT.(IPLACE.EQ.NCL.AND.NEWORB.GE.NORB)) THEN
C
            IF(IPLACE .EQ. 1 ) THEN
              IICL(1) = 1
              IIOC(1) = 2
            ELSE
              IICL(IPLACE) = IICL(IPLACE-1) + 1
              IIOC(IICL(IPLACE)) = 2
            END IF
            GOTO 800
          ELSE
C. No more inactive configurations
             GOTO 2001
          END IF
#if defined (VAR_JEPPE)
... entering IF -- ENDIF block not permitted
    (go to 800)
          END IF
          IFRSTC = 0
#else
  810     IFRSTC = 0
#endif
C..       CHECK RAS1 and RAS 3
          IEL1C = 0
          IEL3C = 0
          ICL1  = 0
          DO  20 ICL = 1,NCL
            IORB = IICL(ICL)
            IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
               IEL1C = IEL1C + 2
               ICL1 = ICL1 + 1
            ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
               IEL3C = IEL3C + 2
            END IF
   20     CONTINUE
          IF(ICL1 .LT.MINCL1) GOTO 2000
          IF(IEL3C.GT.NEL3MX) GOTO 2000
            IIICHK = 1
            IF(ICL1 .LT.MINCL1.AND. IIICHK.EQ.1) THEN
c Next higher combination with a higher number of inactive orbitals
              DO 41 ICL = 1,ICL1+1
                IIOC(IICL(ICL)) = 0
                IICL(ICL) = ICL
                IIOC(ICL) = 2
   41         CONTINUE
              IPLACE=ICL1+1
              IF( IPLACE.GE.NCL) GOTO 2001
              GOTO 800
            END IF
            IF(IEL3C.GT.NEL3MX) GOTO 2000
C. Highest orbital not occupied
         MXMPTY = NORB
         IORB = NORB+1
C. begin while
   12      CONTINUE
           IORB = IORB - 1
           IF(IIOC(IORB) .EQ. 2 ) THEN
             MXMPTY = IORB-1
             IF( IORB .NE. 1 ) GOTO 12
           END IF
C. End while
          IF( IPRINT.GT.30 ) THEN
            WRITE(LUPRI,*) ' Next inactive configuration '
            CALL IWRTMA(IICL,1,NCL,1,NCL)
          END IF
C
C. first active configuration
          IORB = 0
          IOP = 0
          DO 30 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 0 ) THEN
              IOP = IOP + 1
              IF( IOP .GT. NOP ) GOTO 31
              IIOC(IORB) = 1
              IIOP(IOP) = IORB
            END IF
   30     CONTINUE
   31     CONTINUE
          IFRSTO = 1
C
C. Next open shell configuration
 1000     CONTINUE
#if defined (VAR_JEPPE)
... entering IF -- ENDIF block not permitted
    (go to 700)
            IF(IFRSTO .EQ. 0 .AND. NOP .NE. 0 ) THEN
#else
            IF (.NOT. (IFRSTO .EQ. 0 .AND. NOP .NE. 0 ) ) GO TO 710
#endif
            IPLACE = 0
  700       CONTINUE
              IPLACE = IPLACE + 1
              IPRORB = IIOP(IPLACE)
              NEWORB = IPRORB + 1
              IIOC(IPRORB) = 0
  690         CONTINUE
                IF(IIOC(NEWORB) .NE. 0 ) THEN
                   NEWORB = NEWORB + 1
                   IF(NEWORB .LE. MXMPTY ) GOTO 690
                END IF
C 691         End of loop
            IF(IPLACE .LT. NOP .AND. NEWORB .LT. IIOP(IPLACE+1)
     &      .OR.
     &      IPLACE .EQ. NOP .AND. NEWORB .LE. MXMPTY ) THEN
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
            ELSE IF (  IPLACE .NE. NOP ) THEN
              IF(IPLACE.EQ.1) THEN
                NEWORB = 1 - 1
              ELSE
                NEWORB = IIOP(IPLACE-1)
              END IF
 671          CONTINUE
                NEWORB = NEWORB + 1
              IF(IIOC(NEWORB) .NE. 0 .AND. NEWORB.LE.MXMPTY ) GOTO 671
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
              GOTO 700
            ELSE
C. No more active configurations , so
              IF( NCL .NE. 0 ) GOTO 2000
              IF( NCL .EQ. 0 ) GOTO 5001
            END IF
#if defined (VAR_JEPPE)
... entering IF -- ENDIF block not permitted
    (go to 700)
           END IF
           IFRSTO = 0
#else
  710      IFRSTO = 0
#endif
          IF( IPRINT.GT.30 ) THEN
            WRITE(LUPRI,*) ' Next active configuration '
            CALL IWRTMA(IIOP,1,NOP,1,NOP)
          END IF
C        RAS  CONSTRAINTS
         IOKAY=1
         IEL1 = IEL1C
         IEL3 = IEL3C
C..      CHECK RAS1 and RAS3
         DO  40 IOP = 1,NOP
           IORB = IIOP(IOP)
           IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
              IEL1 = IEL1 + 1
           ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
              IEL3 = IEL3 + 1
           END IF
   40    CONTINUE
           IR3CHK  = 1
           IF(IEL3.GT.NEL3MX.AND.IR3CHK.EQ.1) THEN
c. Direct generation of next string with a lower number of
c. electrons in RAS III
c. Procedure : Obtain substring that contains on electron
c  in RAS I, modify this to lowest possible string.
c. Add one to remaining string
c
c. Number of electrons in substring
             IFSTR3 = 0
             DO 5607 IOP = 1, NOP
               IF(IIOP(IOP).GE.IORB3F) THEN
                IFSTR3 = IOP
                GOTO 5608
               END IF
 5607        CONTINUE
 5608        CONTINUE
             IF(IFSTR3.NE.NOP) THEN

c. Lowest possible string with NOP electrons
               DO 5610 K = 1, IFSTR3
                 IIOC(IIOP(K)) = 0
 5610          CONTINUE
c
               KEL = 0
               KORB = 0
 5630          CONTINUE
                 KORB = KORB + 1
                 IF(IIOC(KORB).NE.2) THEN
                   KEL = KEL + 1
                   IIOC(KORB) = 1
                   IIOP(KEL) = KORB
                 END IF
               IF(KEL.NE.IFSTR3) GOTO  5630
               IPLACE = IFSTR3
C?           IF(MAXPRT.GT.0) THEN
C?             WRITE(LUPRI,*) ' Updated IIOC AND IIOP '
C?             CALL IWRTMA(IIOC,1,NORB,1,NORB)
C?             CALL IWRTMA(IIOP,1,NOP,1,NOP)
C?           END IF
               GOTO 700
             END IF
           END IF
         IF( IEL1 .LT. NEL1MN .OR. IEL3 .GT. NEL3MX
     &  .OR. IEL1 .GT. NEL1MX .OR. IEL3 .LT. NEL3MN ) GOTO  999
C. Spatial symmetry
         ISYM = 1
         DO 920 IOP = 1, NOP
           ISYM = SYMPRO(ISYM,ORBSYM(IIOP(IOP)))
  920    CONTINUE
         IF(ISYM.EQ.IREFSM) THEN
           IF( IPRINT.GT.30 )
     &     WRITE(LUPRI,1120) ( IIOC(I),I = 1,NORB )
 1120      FORMAT(1H0,'  configuration included ',15I3)
           JCONF=JCONF+1
C
           DO 60 ICL = 1, NCL
             ICONF(ICFREE-1+ICL) = IICL(ICL)
   60      CONTINUE
           DO 61 IOP = 1, NOP
             ICONF(ICFREE-1+NCL+IOP) = IIOP(IOP)
   61      CONTINUE
           ICFREE = ICFREE + NOP + NCL
         END IF

C
C** LOOP OVER active configurations end
C
  999   CONTINUE
        IF( NOP .EQ. 0 .AND. NCL .EQ. 0 ) GOTO 5001
        IF( NOP .EQ. 0 ) GOTO 2000
        GOTO 1000
 2001 CONTINUE
 5000 CONTINUE
 5001 CONTINUE

C
      IF( IPRINT.GT.30 ) THEN
        WRITE(LUPRI,'(/A,I3)') '  Configurations of symmetry ', IREFSM
        WRITE(LUPRI,*)        ' ================================='
        IBAS = 0
        DO 1200 IOPEN = MINOP,MAXOP
          ITYPE = IOPEN - MINOP + 1
          ICL = (NEL-IOPEN)/2
          IOC = IOPEN + ICL
          LICONF = NCNFTP(ITYPE,IREFSM)
          WRITE(LUPRI,'(/A,2I3)')
     &    '  Type with number of closed and open orbitals ',ICL,IOPEN
          WRITE(LUPRI,'(A,I7)')
     &    '  Number of configurations of this type',LICONF
          DO 1180 JCONF = 1,LICONF
            WRITE(LUPRI,'(3X,20I3)') (ICONF(IBAS+IORB),IORB=1,IOC)
            IBAS = IBAS + IOC
 1180     CONTINUE
 1200   CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck iweylf */
      FUNCTION IWEYLF(NOPEN,MULTS)
C
C NUMBER OF CSF:S WITH NOPEN ORBITALS AND TOTAL MULTIPLICITY
C MULTS ACCORDING TO WEYLS FORMULAE
C
C     (2S+1)/(NOPEN+1) * BION(NOPEN+1/0.5NOPEN-S)
C
#include <implicit.h>
#include <priunit.h>
C
      NTEST = 0
C
      IF(NOPEN.EQ.0 .AND. MULTS .EQ. 1 ) THEN
        NCSF = 1
      ELSEIF(MOD(MULTS-1,2) .NE. MOD(NOPEN,2) ) THEN
        NCSF = 0
      ELSEIF(MOD(MULTS-1,2) .EQ. MOD(NOPEN,2) ) THEN
        NCSF = MULTS*IBION(NOPEN+1,(NOPEN+1-MULTS)/2)/(NOPEN+1)
      END IF
C
      IWEYLF = NCSF
C
      IF(NTEST .NE. 0 ) THEN
        WRITE(LUPRI,'(A,4I4)')
     &  '  IWEYLF SAYS : NOPEN MULTS NCSF : ', NOPEN,MULTS,NCSF
      END IF
C
      RETURN
      END
C  /* Deck msymp2 */
      FUNCTION MSYMP2(ISTR,NELMNT,ORBSYM,SYMPRO)
C
C A STRING of objects given. find total symmetry
C
#include <implicit.h>
      INTEGER SYMPRO(8,8),ORBSYM(1)
      INTEGER ISTR(1)
C
      ITOTSM = 1
      DO 100 I = 1, NELMNT
       ITOTSM = SYMPRO(ITOTSM,ORBSYM(ISTR(I)))
  100 CONTINUE
C
      MSYMP2 = ITOTSM
C
      RETURN
      END
C  /* Deck nexci */
      SUBROUTINE NEXCI(OLSEN,ENER,NCVAR,D,XVEC,RES,
     &                 DIAG,IPRPHP,WRK,LWRK)
C
C THE PURPOSE OF THIS ROUTINE IS TO CREATE A NEW CSF TRIAL VECTOR
C FROM THE RESIDUAL
C IF OLSEN THE TRIAL VECTOR IS CREATED ACCORDING TO JEPPES FORMULA
C IF .NOT.OLSEN THE TRIAL VECTOR IS OBTAINED FROM BLOCK DAVIDSON
C
C MOTECC-90: The purpose of this module, NEXCI, and the algorithms used
C            are described in Chapter 8 Appendix 8C of MOTECC-90
C            "Iterative Solution of the CI Eigenvalue Problem"
C
C
C INPUT:
C ENER   ENERGY OF CURRENT APPROXIMATION
C NCVAR  NUMBER OF CONFIGURATION VARIABLES
C NPCVAR DIMENSION OF H(0) SPACE
C XVEC   SOLUTION VECTOR FOR CURRENT APPROXIMATION (ONLY NEEDED
C                                                WHEN OLSEN = .TRUE.
C RES    RESIDUAL ALWAYS NEEDED
C DIAG   DIAGONAL CI ELEMENTS, contain at the end h(0) information
C IPRPHP PRINT LEVEL
C
C OUTPUT:
C D  CONTAIN NEW CONFIGURATION TRIAL VECTOR (can use the same
C                                   locations as xvec,res,or diag)
C
#include <implicit.h>
#include <priunit.h>
C
#include <phpinf.h>
C
      LOGICAL OLSEN
C
      DIMENSION D(*),XVEC(*),RES(*),DIAG(*),WRK(*)
C
      IF (IPRPHP .GE. 50) THEN
         WRITE(LUPRI,*) 'NEWCI: OLSEN =',OLSEN
         WRITE(LUPRI,*) 'NEWCI: ENER  =',ENER
      END IF
      KFREE = 1
      LFREE = LWRK
      IF (OLSEN) THEN
         INVCOR = 1
      ELSE
         INVCOR = 0
      END IF
      CALL NEWDIR(D,XVEC,RES,DIAG,ENER,NCVAR,NPCVAR,
     &            DIAG(KPNTR),DIAG(KPEVAL),DIAG(KPEVEC),
     &            INVCOR,IPRPHP,WRK(KFREE),LFREE)
C     CALL NEWDIR(D,X,G,DIAG,E,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
C    &            INVCOR,NTEST,WORK,LFREE)
      RETURN
      END
C  /* Deck phpget */
      SUBROUTINE PHPGET(LSTSYM,NCVAR,XNDXCI,FCAC,H2AC,DIAG,
     &                  ECORE,IPWAY,IPRPHP,WRK,LWRK)
C
C L.R. 28-Oct-1990 hjaaj : reorder if SLREOR
C
C THE PURPOSE OF THIS ROUTINE IS TO CREATE ALL INFORMATION WHICH
C IS REQUIRED TO SET UP H(0) IN THE DIAGONAL BASIS TOGETHER
C WITH POINTERS THAT SHOW THE LOCATIONS OF THE DIAGONAL ELEMENTS
C THAT IS INCLUDED IN H(0) C FROM THE RESIDUAL
C
C INPUT:
C LSTSYM SYMMETRY OF THE CONSIDERED CONFIGURATION SPACE
C NCVAR  NUMBER OR CONFIGURATION VARIABLES
C DIAG   DIAGONAL CI ELEMENTS
C ECORE  core energy for CI ***OR*** -EACTIV for MC Hessian.
C IPWAY  method for selecting PHP
C  = 1   select up to MXPDIM lowest elements of DIAG
C  = 2   select up to MXPDIM first elements of DIAG (test option)
C  = 3   select up to MXPDIM elements of DIAG with highest abs.val.
C IPRPHP PRINT LEVEL
C
C OUTPUT:
C NPCVAR DIMENSION OF H(0)
C DIAG   HAS APPENDED AT THE END POINTERS TELLING WHICH DIAGONAL
C        ELEMENTS ARE CONTAINED IN H(0),EIGENVECTORS FOR H(0), AND
C        EIGENVALUES FOR H(0)
C
#include <implicit.h>
#include <priunit.h>
C
C Used from common blocks:
C  INFINP : FLAG(27),LSYM,?
C  PHPINF : MXPDIM,NPCVAR,...
C  CIINFO : ICOMBI,IPSIGN,NDTASM,NCSASM
C  DETBAS : KLTSOB,?
C  CBREOR : SLREOR
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <phpinf.h>
#include <mxpdim.h>
#include <ciinfo.h>
#include <detbas.h>
#include <csfbas.h>
#include <strnum.h>
#include <cbreor.h>
C
      DIMENSION XNDXCI(*),FCAC(*),H2AC(*),DIAG(*),WRK(*)
C
C
C     Initialize:
C
      NPCVAR = 0
      IDODGN = 1
C     ... always diagonalize PHP matrix in this version
C
C     Exit if zero dimension of PHP
C
      IF (MXPDIM .LE. 0) GO TO 9999
C
      KFREE = 1
      LFREE = LWRK
      KSAV  = KFREE
      LSAV  = LFREE
      CALL MEMGET('REAL',KPHP  ,MXPDIM*MXPDIM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KIPCNF,MXPDIM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEBD,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRJ   ,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRK   ,N2ASHX,WRK,KFREE,LFREE)
      CALL GETFIJ(WRK(KRJ),WRK(KRK),H2AC)
      CALL DSPTSI(NASHT,FCAC,WRK(KONEBD))
      IF (SLREOR) THEN
C        ... reorder from Sirius to Lunar order
         CALL REORMT(WRK(KONEBD),WRK(KFREE),NASHT,NASHT,
     &               XNDXCI(KLTSOB),XNDXCI(KLTSOB))
         CALL DCOPY(N2ASHX,WRK(KFREE),1,WRK(KONEBD),1)
         CALL REORMT(WRK(KRJ),WRK(KFREE),NASHT,NASHT,
     &               XNDXCI(KLTSOB),XNDXCI(KLTSOB))
         CALL DCOPY(N2ASHX,WRK(KFREE),1,WRK(KRJ),1)
         CALL REORMT(WRK(KRK),WRK(KFREE),NASHT,NASHT,
     &               XNDXCI(KLTSOB),XNDXCI(KLTSOB))
         CALL DCOPY(N2ASHX,WRK(KFREE),1,WRK(KRK),1)
      END IF
C
C     Store determinant CI off-sets in C arrays and HC arrays
C
      CALL CIOFF(LSTSYM,1,XNDXCI,IPRPHP)
      CALL CIOFF(LSTSYM,2,XNDXCI,IPRPHP)
C     CALL CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
C
C
C000127-HJAaJ:
C     Find out if determinants (ICCSF .eq. 0) or CSFs:
C
      IF ( FLAG(27) .OR. (LSTSYM .NE. LSYM) ) THEN
         ICCSF = 0
      ELSE
C
C        LSTSYM is thus the reference symmetry, LSYM.
C        ABACUS and RESPONSE use CSF for singlet operators
C        and det for triplet operators in reference symmetry,
C        we find what is the case for this call by comparing
C        to NCSASM and NDTASM:
C
         IF ( NCVAR .EQ. NCSASM( LSTSYM ) ) THEN
            ICCSF = 1
         ELSE IF ( NCVAR .EQ. NDTASM( LSTSYM ) ) THEN
            ICCSF = 0
         ELSE
            WRITE (LUPRI,'(/A,I3/A/4I12)')
     &         ' CSF ERROR in PHPGET, LSYM =',LSYM,
     &         ' LSTSYM, NCVAR, NCSASM(LSTSYM), NDTASM(LSTSYM):',
     &           LSTSYM, NCVAR, NCSASM(LSTSYM), NDTASM(LSTSYM)
            CALL QUIT('NCVAR disagrees with N*ASM(:) in PHPGET')
         END IF
      END IF
C
      IF ( ICCSF .EQ. 0 ) THEN
C     IF ( determinants ) THEN
C     980820-hjaaj: CSFs only for symm. LSYM in this version
         IHSYM = 1
C        ... PHP matrix is symmetric
         PSIGN = IPSIGN
         NPREV = 0
C        NOTE: if NPREV .gt. 0 then DIAG is destroyed
C        NPREV is number of solution vectors to analyze
         KFREE0 = KFREE
         CALL PHPDET(DIAG,IPWAY,MXPDIM,XNDXCI,WRK,KFREE0,LWRK,
     &               LSTSYM,IHSYM,WRK(KPHP),WRK(KONEBD),WRK(KRJ),
     &               WRK(KRK),NASHT,MULD2H,ECORE,ICOMBI,PSIGN,
     &               NPCVAR,DIAG(KPEVEC),DIAG(KPEVAL),IDODGN,
     &               NCVAR,DIAG(KPNTR),NPREV,H2AC,IPRPHP)
      ELSE
         IF (IPWAY .EQ. 3) THEN
            WRITE (LUPRI,*)
     *      '***ERROR, .RESPHP is not implemented for CSF basis'
            WRITE (LUPRI,*) 'Use .DETERMINANTS or remove .RESPHP'
            CALL QUIT('ERROR, .RESPHP is not implemented for CSF basis')
         END IF
         CALL PHPCSF(WRK(KPHP),DIAG(KPNTR),WRK(KIPCNF),MXPDIM,
     &               XNDXCI,XNDXCI(KDTOC),XNDXCI(KDFTP),
     &               XNDXCI(KICONF(1)),
     &               LSTSYM,WRK(KONEBD),ECORE,0,NASHT,
     &               WRK(KFREE),LFREE,NCVAR,NACTEL,NAEL,NBEL,
     &               IPWAY,NPCVAR,NPCNF,
     &               IDODGN,DIAG(KPEVEC),DIAG(KPEVAL),DIAG,IPRPHP,
     &               H2AC,WRK(KRJ),WRK(KRK),XNDXCI(KLTSOB))
C        CALL PHPCSF(PHP,IPCSF,IPCNF,MXPDIM,
C    &               XNDXCI,DTOC,IPRODT,ICONF,
C    &               IREFSM,ONEBOD,ECORE,NINOB,NACTOB,
C    &               SCR,LSCR,NCONF,NEL,NAEL,NBEL,IPWAY,NPCSF,NPCNF,
C    &               IDODIA,EIGVEC,EIGVAL,DIAG,NTEST,
C    &               FIJKL,RJ,RK,ILTSOB)
      END IF
      CALL MEMREL('PHPGET',WRK,KSAV,KSAV,KFREE,LFREE)
      IF (NPCVAR .NE. MXPDIM .AND. IPRPHP .GT. 5) THEN
         WRITE (LUPRI,*) 'Dimension of explicit CI matrix (PHP matrix):'
         WRITE (LUPRI,*) '   -- requested maximum',MXPDIM
         WRITE (LUPRI,*) '   -- actually used    ',NPCVAR
      END IF
 9999 CONTINUE
      RETURN
      END
C  /* Deck phpini */
      SUBROUTINE PHPINI(LPHPT,NCVAR,NOVAR,MAXPHP)
C
#include <implicit.h>
C
C SPACE ALLOCATION FOR NEXCI ROUTINES
C
#include <phpinf.h>
      MXPDIM = MAXPHP
      KDIACI = 1
      KDIAOR = KDIACI + NCVAR
      KPNTR  = KDIAOR + NOVAR
      KPEVAL = KPNTR  + MXPDIM
      KPEVEC = KPEVAL + MXPDIM
      KTOT   = KPEVEC + MXPDIM*MXPDIM
      LPHP   = KTOT   - 1
      LPHPT  = LPHP
      RETURN
      END
#ifdef UNDEF
/* Comdeck fra_jeppe */
C Nov. 1989
C Hej arhusianere , her kommer lidt preconditionerings rutiner.
C rutinen til at beregne CSF hamilton matricen hedder PHPCSF.
C Den skulle have veldokumenteret input.
C Jeg har sat ind kaldet til  SIRIUS versionen af DIHDJ
C i rutinen, men dette er ikke testet-
C CSF hamilton matricen bliver konstrueret som folger
C
C PHPCSF-WRTMAT,ISTVC2,FNDMNX,XTRDI (standard routines )
C       -GETCNF (occupation vector from configuration number )
C       -CNHCN  ( hamiltonian matrix for a pair of configurations )
C             - CNFSTR ( strings of a given configuration )
C             - DIHDJ  ( hamiltonian matrix between two sets of dets )
C             - DGMM2  ( sign changes for going between dets and csf
C             - MATML4 ( Transform from det to CSf basis )
C
C rutinen til at beregne den modificerede
C Davidson retning hedder NEWDIR ,ogs} denne er rimeligt kommenteret
#endif
C  /* Deck newdir */
      SUBROUTINE NEWDIR(D,X,G,DIAG,E,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &                  INVCOR,NTEST,WORK,LFREE)
c
c Calculate
c
c  D = (H0-E)** (-1) * G    (INVCOR = 0 )
c
c  D = (H0-E)** (-1) * G - ALPHA * (H0 - E)**(-1) * X
c
c       ALPHA = X(T)(H0-E)**(-1)*D / X(T)(H0-E)**(-1)*X  (INVCOR .NE.0)
c
c The latter correction corresponds to inverse iteration
c correction to Davidson
c
c Where H0 consists of a diagonal Diag
c and a block matrix of dimension NPRDIM.
c
c The block matrix is defined by
c ==============================
c
c  NPRDIM : Size of block matrix
c  IPNTR(I) : Scatter array, gives adress of subblock element
c             I in full matrix
c  PEIGVL   : Eigenvalues of subblock mateix
c  PEIGVC   : Eigenvectors of subblock matrix
c
c Input
c=======
c X : for eigenvalue problem X is current eigenvector
c     (for INVCOR = 0 X can be a dummy variable )
c G : for eigenvalue problem G = (H - E ) * X
c Diag : Diaginal of matrix
c E : Energy for shift
c NVAR : Dimension of full matrix
c NPRDIM,IPNTR,PEIGVL,PEIGVC : See above
c INVCOR : use(.NE.0) , do not use (.eq.0) inverse correction
c Modification
c Work : Scratch space , at least ??
c
c Output
c ======
c D as given above, code has been constructed so D
c can occupy the same place as either X,G,DIAG
c
c Scratch space
c===============
c Should at least be of length ???
c
c Externals   GPRCTV
c===========
c
#include <implicit.h>
#include <priunit.h>
c
      DIMENSION D(*), X(*),G(*),DIAG(*)
      DIMENSION IPNTR(*),PEIGVL(*),PEIGVC(*), WORK(*)
c
c
      IF( INVCOR .EQ. 0 ) THEN
c (H0 - E ) **(-1) * G , store in D
C       CALL GPRCTV(DIAG,VECIN,VECUT,NVAR,NPRDIM,IPNTR,
C    &              PEIGVL,PEIGVC,SHIFT,WORK,NTEST)
        CALL GPRCTV(DIAG,G,D,NVAR,NPRDIM,IPNTR,
     &              PEIGVL,PEIGVC,-E,WORK,XDUMMY,NTEST)
        IF(NTEST.GE.50) THEN
          WRITE(LUPRI,*) ' Information from NEWDIR'
          WRITE(LUPRI,*) ' ======================='
          WRITE(LUPRI,*) ' New direction = (H0 - E) **(-1) * RESIDUAL'
          WRITE(LUPRI,*) ' E = ',E
          CALL WRTMAT(D,1,NVAR,1,NVAR,0)
        END IF
      ELSE
        IF(NTEST.GT.5) THEN
          WRITE(LUPRI,*) ' Information from NEWDIR'
          WRITE(LUPRI,*) ' ======================='
          WRITE(LUPRI,*) ' E = ',E
        END IF
c (H0 - E ) **(-1) * G , store in G
        CALL GPRCTV(DIAG,G,G,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &              -E,WORK,XDUMMY,NTEST)
c X(T) (H0 - E) ** (-1) G
        XH0MEG = DDOT(NVAR,X,1,G,1)
        IF(NTEST.GE.50) THEN
          WRITE(LUPRI,*) '  (H0 - E) **(-1) * RESIDUAL '
          CALL WRTMAT(G,1,NVAR,1,NVAR,0)
        END IF
c (H0 - E ) **(-1) * X , store in X
        CALL GPRCTV(DIAG,X,X,NVAR,NPRDIM,IPNTR,PEIGVL,PEIGVC,
     &              -E,WORK,XH0MEX,NTEST)
        FACTOR = -XH0MEG/XH0MEX
        IF (NTEST.GT.5) THEN
           WRITE(LUPRI,*) ' X(T) (H0 - E) ** (-1) G = XH0MEG ', XH0MEG
           WRITE(LUPRI,*) ' X(T) (H0 - E) ** (-1) X = XH0MEX ', XH0MEX
           WRITE(LUPRI,*) ' FACTOR = -XH0MEG/XH0MEX ', FACTOR
        END IF
        IF(NTEST.GE.50) THEN
          WRITE(LUPRI,*) '  (H0 - E) **(-1) * X '
          CALL WRTMAT(X,1,NVAR,1,NVAR,0)
        END IF
c
        CALL VECSUM(D,G,X,1.0D0,FACTOR,NVAR)
        IF(NTEST.GE.50) THEN
          WRITE(LUPRI,*) ' New direction with inverse iteration '//
     &                   'correction'
          CALL WRTMAT(D,1,NVAR,1,NVAR,0)
        END IF
      END IF
c
      RETURN
      END
C  /* Deck gprctv */
      SUBROUTINE GPRCTV(DIAG,VECIN,VECUT,NVAR,NPRDIM,IPNTR,
     &                  PEIGVL,PEIGVC,SHIFT,WORK,XH0PSX,NTEST )
c
c Calculate inverted general preconditioner matrix times vector
c
c  Vecut=  (H0 + shift )-1 Vecin
c
c  and XH0PSX = X(T) (H0 + shift )** -1  X
c
c Where H0 consists of a diagonal Diag
c and a block matrix of dimension NPRDIM.
c
c Note : The diagonal elements in DIAG corresponding to
c        elements in the subspace are neglected,
c        i.e. their elements can have arbitrary value
c        without affecting the results
c
c The block matrix is defined by
c ==============================
c
c  NPRDIM : Size of block matrix
c  IPNTR(I) : Scatter array, gives adress of subblock element
c             I in full matrix
c  PEIGVL   : Eigenvalues of subblock mateix
c  PEIGVC   : Eigenvectors of subblock matrix
c
c Jeppe Olsen , Sept. 1989
c
c Input
c=======
c DIAG : Diagonal of matrix
c VECIN : Input vector
c NVAR : Dimension of full matrix
c NPRDIM,PEIGVL,PEIGVC : See above
c SHIFT : constant ADDED to diagonal
c WORK : Scratch array , at least 2*NPRDIM
c NTEST : Print level
c
c Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
c ==========
c
c Output
c========
c VECUT : Output vector (you guessed ?? ), can occupy same space
c         as VECIN or DIAG
c XH0PSX  = X(T)(H0+SHIFT)**(-1)X
c
#include <implicit.h>
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),PEIGVL(*),PEIGVC(*)
      DIMENSION WORK(*)
c
      IF(NPRDIM.NE.0) THEN
        CALL GATVEC(WORK(1),VECIN,IPNTR,NPRDIM)
c X(T)(DIAG+SHIFT)X in subspace, for later subtraction
        CALL GATVEC(WORK(1+NPRDIM),DIAG,IPNTR,NPRDIM)
        CALL DIAVC3(WORK(1+NPRDIM),WORK(1),
     &       WORK(1+NPRDIM),SHIFT,NPRDIM,X1)
       ELSE
         X1 = 0.0D0
       END IF
c
      CALL DIAVC3(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
c
      IF(NPRDIM .NE. 0 ) THEN
         CALL SCAVEC(VECUT,WORK(1),IPNTR,NPRDIM)
         CALL SBINTV(NPRDIM,PEIGVC,PEIGVL,SHIFT,
     &              IPNTR,VECUT,VECUT,WORK(1),WORK(1+NPRDIM),X3,NTEST)
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(LUPRI,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3
c
      RETURN
      END
C  /* Deck sbintv */
      SUBROUTINE SBINTV(NSBDIM,EIGVC,EIGVL,SHIFT,INDEX,VECI,VECO,X1,X2,
     &                  XHPSX,NTEST)
c
c INVERTED SHIFTED SUBSPACE MATRIX TIMES VECTOR
c
c Last revision, oct 1989
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION EIGVC(NSBDIM,NSBDIM),EIGVL(NSBDIM),INDEX(NSBDIM)
      DIMENSION X1(NSBDIM),X2(NSBDIM)
      DIMENSION VECI(*),VECO(*)
c
      CALL GATVEC(X1,VECI,INDEX,NSBDIM)
      CALL MATVCD(EIGVC,X1,X2,NSBDIM,NSBDIM,1)
      CALL DIAVC3(X1,X2,EIGVL,SHIFT,NSBDIM,XHPSX)
      CALL MATVCD(EIGVC,X1,X2,NSBDIM,NSBDIM,0)
      CALL SCAVEC(VECO,X2,INDEX,NSBDIM)
C
      IF( NTEST .GE. 50 ) THEN
        WRITE(LUPRI,*) ' OUTPUT FROM SBINTV, VECTOR IN GATHERED FORM '
        write(lupri,*) ' Shift =',SHIFT
        write(lupri,*) ' INDEX =',(INDEX(I),I=1,NSBDIM)
        CALL WRTMAT(X1,1,NSBDIM,1,NSBDIM,0)
      END IF
C
      RETURN
      END
C  /* Deck diavc3 */
      SUBROUTINE DIAVC3(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
C
C VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
C
C VDSV = SUM(I) VECIN(I) ** 2 /( DIAG(I) + SHIFT )
C
#include <implicit.h>
      DIMENSION VECOUT(*),VECIN(*),DIAG(*)
C
      THRES=1.0D-5
      VDSV = 0.0D0
      DO 100 I=1,NDIM
c
        DIVIDE=DIAG(I)+SHIFT
        IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
c
        VDSV = VDSV + VECIN(I) ** 2 /DIVIDE
        VECOUT(I)=VECIN(I)/DIVIDE
c
  100 CONTINUE
      RETURN
      END
C  /* Deck phpcsf */
      SUBROUTINE PHPCSF(PHP,IPCSF,IPCNF,MXPDIM,
     &                  XNDXCI,DTOC,IPRODT,ICONF,
     &                  IREFSM,ONEBOD,ECORE,NINOB,NACTOB,
     &                  SCR,LSCR,NCONF,NEL,NAEL,NBEL,IPWAY,NPCSF,NPCNF,
     &                  IDODIA,EIGVEC,EIGVAL,DIAG,NTEST,
     &                  FIJKL,RJ,RK,ILTSOB)
c
c Obtain primary subspace and obtain
c explicit representation of hamilton matrix in subspace
c
c ARGUMENTS :
c ===========
c   PHP    : hamilton matrix in subspace ( output)
c   IPCSF  : CSF's defining subspace (output)
c   IPCNF  : Configurations defining subspace (output )
c   MXPDIM : Largest allowed dimension of subspace (Input)
c   XNDXCI : List with Ci-information as obtained in DETINF(INPUT)
c   DTOC   : Transformation matrix between CSF's and DET's (input)
c   IPRODT : Prototype determinants (input)
c   ICONF  : List of configurations  (input)
c   IREFSM : symmetry of considered CI space (input)
c   Onebod : one body hamilton matrix in rectangular form (input)
c   ECORE  : Core energy (input)
c   NINOB  : Number of inactive orbitals(input)
c   NACTOB : Number of active orbitals (input)
c   SCR    : Scratch array of length LSCR
c   LSCR   : length of SCR
c   NCONF  : Number of configurations of symmetry IREFSM
c   IPWAY  : Defines way of choosing Primary space
c    = 2   : use the first configurations (until atmost
c            MXPDIM CSF's have been included )
c   NPCNF  : Number of primary configurations obtained (output)
c   NPCSF  : Number of primary CSF's obtained  (OUTPUT)
c   IDODIA : Diagonalize primary matrix
c   EIGVEC : Eigenvectors of primary subspace
c   EIGVAL : Eigenvalues of primary subspace
c   DIAG   : Hamilton diagonal over CSF's
c
c For SIRIUS interface :
c FIJKL,RJ,RK,ILTSOB
c
c Externals :
c    FNDMNX,WRTMAT,ISTVC2,XTRDI
c    GETCNF,CNHCN
c Jeppe Olsen , Summer of '89
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION PHP(*),IPCSF(*),IPCNF(*),XNDXCI(*),DIAG(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*),SCR(*)
      DIMENSION EIGVEC(*),EIGVAL(*)
      DIMENSION FIJKL(*),RJ(*),RK(*),ILTSOB(*)
#include <spinfo.h>
c
c* 1 : Obtain primary subspace
c
      IF(IPWAY .EQ. 2 ) THEN
        NPCNF = 0
        NPCSF = 0
        ICNF = 0
        ICSF = 0
        DO 100 ITYP = 1, NTYP
          NJCNF = NCNFTP(ITYP,IREFSM)
          NIRREP = NCSFTP(ITYP)
          DO 90 IICNF = 1, NJCNF
            ICNF = ICNF + 1
            IF(NPCSF+NIRREP .LE. MXPDIM) THEN
              NPCNF = NPCNF + 1
              IPCNF(NPCNF) = ICNF
              DO 80 IICSF = 1, NIRREP
                NPCSF = NPCSF + 1
                IPCSF(NPCSF) = NPCSF
   80         CONTINUE
            ELSE
              GOTO 101
            END IF
   90     CONTINUE
  100   CONTINUE
  101   CONTINUE
c
      ELSE IF( IPWAY .EQ. 1) THEN
c          Obtain lowest CSF's
c
c 1 :  local Diagonal elements over configurations
c
        IICSF = 1
        IICNF = 1
        KLDIA = 1
        KLFREE = KLDIA + NCONF
        DO 300 ITYP = 1, NTYP
          NJCNF = NCNFTP(ITYP,IREFSM)
          NIRREP = NCSFTP(ITYP)
          DO 290 ICNF = 1, NJCNF
            SCR(KLDIA-1+IICNF) = DIAG(IICSF)
            IICNF = IICNF + 1
            IICSF = IICSF + NIRREP
  290     CONTINUE
  300   CONTINUE
        IF( NTEST .GE. 15 ) THEN
          WRITE(LUPRI,*) ' CI diagonal over configurations '
          WRITE(LUPRI,*) ' ================================'
          CALL WRTMAT(SCR(KLDIA),1,NCONF,1,NCONF,0)
        END IF
c Largest element
        IMAX = IDMAX(NCONF,SCR(KLDIA),1)
        XMAX = SCR(KLDIA-1+IMAX)
C       XMAX = FNDMNX(SCR(KLDIA),NCONF,2)
c loop over lowest configurations
c
        NPCSF = 0
        NPCNF = 0
c
        IFINIT = 0
  400   CONTINUE
c
        XMIN = XMAX + 1.0D0
        IMIN = 0
c
        IICNF = 1
        ICSFOF = 1
        DO 500 ITYP = 1, NTYP
          NIRREP = NCSFTP(ITYP)
          DO 450 ICNF = 1,NCNFTP(ITYP,IREFSM)
            IF(SCR(KLDIA-1+IICNF).LT.XMIN) THEN
              XMIN = SCR(KLDIA-1+IICNF)
              IMIN = IICNF
              ICSFMN = ICSFOF
              NCSFMN = NIRREP
            END IF
c
            IICNF = IICNF + 1
            ICSFOF = ICSFOF + NIRREP
  450     CONTINUE
  500   CONTINUE
c
c* Next lowest element has been found
c
c 1 : add new configuration
        IF(NPCSF + NCSFMN .LE. MXPDIM ) THEN
          NPCNF = NPCNF + 1
          IPCNF(NPCNF) = IMIN
          CALL ISTVC2(IPCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
          NPCSF = NPCSF + NCSFMN
c Mask
          SCR(KLDIA-1+IMIN) = XMAX + 1
        ELSE
          IFINIT = 1
c 2 : No space for this configuration, remove previous
c     configurations with the same diagonal value
          IICNF = NPCNF+1
  600     CONTINUE
            IICNF = IICNF - 1
            DIAVAL = SCR(KLDIA-1+IPCNF(IICNF) )
            IF(ABS(DIAVAL-XMIN) .LE. 1.0D-10) THEN
              NPCNF = NPCNF -1
C             CALL GETCNF(SCR(KRCONF),IRTYP,IPCNF(ICNR),ICONF,IREFSM,NEL
              CALL GETCNF(SCR(KLFREE),ITYP,IPCNF(IICNF),
     &        ICONF,IREFSM,NEL)
              NPCSF = NPCSF - NCSFTP(ITYP)
              GOTO 600
            END IF
        END IF
        IF(IFINIT.EQ.0.AND.NPCSF.LT.NCONF) GOTO 400
      ELSE
        WRITE (LUPRI,*) ' ERROR in PHPCSF, illegal IPWAY =',IPWAY
        CALL QTRACE(6)
        CALL QUIT('Illegal IPWAY in PHPCSF')
      END IF
c
      IF(NTEST .GE. 2 ) THEN
        WRITE(LUPRI,*) ' Output from PHPCSF '
        WRITE(LUPRI,*) ' ================== '
        WRITE(LUPRI,*)
     &  ' Number of Configurations in primary subspace',NPCNF
        WRITE(LUPRI,*)
     &  ' Number of CSFs in primary subspace          ',NPCSF
        IF( NTEST .GE. 5 ) THEN
          WRITE(LUPRI,*) ' Configurations included : '
          CALL IWRTMA(IPCNF,1,NPCNF,1,NPCNF)
          WRITE(LUPRI,*) ' CSFs included : '
          CALL IWRTMA(IPCSF,1,NPCSF,1,NPCSF)
        END IF
      END IF
c
c* 2 : Construct Hamiltonian matrix in subspace
c
c* size of largest possible subblock
      MXCSFC = 0
      DO 10 ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP) )
   10 CONTINUE
      IF( NTEST .GE.2 ) WRITE(LUPRI,*)
     &' Largest number of CSFs for given configuration  ', MXCSFC
c
      KLFREE = 1
c
      KLCONF = KLFREE
      KLFREE = KLFREE + NEL
c
      KRCONF = KLFREE
      KLFREE = KLFREE + NEL
c
      KLPHPS = KLFREE
      KLFREE = KLFREE + MXCSFC ** 2
      LFREE  = LSCR - KLFREE + 1
c
      IILB = 1
      DO 200 ICNL = 1, NPCNF
        CALL GETCNF(SCR(KLCONF),ILTYP,IPCNF(ICNL),ICONF,IREFSM,NEL)
        NCSFL = NCSFTP(ILTYP)
        IIRB = 1
        DO 190 ICNR = 1, ICNL
          IF(NTEST .GE. 20 )
     &    write(LUPRI,*) ' >>> ICNL,ICNR <<< ',ICNL,ICNR
          CALL GETCNF(SCR(KRCONF),IRTYP,IPCNF(ICNR),ICONF,IREFSM,NEL)
          NCSFR = NCSFTP(IRTYP)
          CALL CNHCN(SCR(KLCONF),ILTYP,SCR(KRCONF),IRTYP,SCR(KLPHPS),
     &               SCR(KLFREE),KLFREE,XNDXCI,NEL,NAEL,NBEL,NINOB,
     &               ECORE,ONEBOD,IPRODT,DTOC,NACTOB,FIJKL,RJ,RK,ILTSOB)
c Copy to PHP matrix
          DO 160 IIL = 1, NCSFL
            IF(IILB.EQ.IIRB) THEN
              IIRMAX = IIL
            ELSE
              IIRMAX = NCSFR
            END IF
            DO 150 IIR = 1, IIRMAX
              IIRACT = IIRB - 1 + IIR
              IILACT = IILB - 1 + IIL
              ILRO = IILACT*(IILACT-1)/2 + IIRACT
              ILRI = (IIR-1)*NCSFL + IIL
              PHP(ILRO) = SCR(KLPHPS-1+ILRI)
  150       CONTINUE
  160     CONTINUE
          IIRB = IIRB + NCSFR
  190   CONTINUE
      IILB = IILB + NCSFL
 200  CONTINUE
c
      IF( NTEST .GE.20 ) THEN
        WRITE(LUPRI,*) ' The subspace hamiltonian '
        WRITE(LUPRI,*) ' ======================== '
        CALL PRSYM(PHP,NPCSF)
      END IF
c
c* 3: Diagonalize PHP
c
      IF( IDODIA .GT. 0 ) THEN
c
         CALL EIGEN(PHP,EIGVEC,NPCSF,0,1)
C.       EXTRACT EIGENVALUES
         CALL XTRCDI(PHP,EIGVAL,NPCSF,1)
c
         IF(NTEST.GT.0)
     &   WRITE(LUPRI,*) ' Range of eigenvalues in subspace ',
     &             EIGVAL(1),EIGVAL(NPCSF)
         IF(NTEST .GE. 2 ) THEN
           WRITE(LUPRI,*) ' Eigenvalues of subspace hamiltonian '
           CALL WRTMAT(EIGVAL,1,NPCSF,1,NPCSF,0)
           IF( NTEST .GE.20 ) THEN
             WRITE(LUPRI,*) ' Eigenvectors of subspace hamiltonian '
             CALL WRTMAT(EIGVEC,NPCSF,NPCSF,NPCSF,NPCSF,0)
           END IF
          END IF
      END IF
c
      RETURN
      END
C  /* Deck getcnf */
      SUBROUTINE GETCNF(KCNF,KTYP,K,ICONF,IREFSM,NEL)
c
c Obtain configuration number K .
c Occupation in KCNF
c Type in KTYP
c
c Jeppe Olsen , summer of '89
c
#include <implicit.h>
#include <priunit.h>
#include <spinfo.h>
c
      DIMENSION KCNF(*),ICONF(*)
c
      ICNFB1 = 1
      ICNFB2 = 1
      KTYP   = 0
      DO 100 JTYP = 1, NTYP
        JOP = JTYP - 1 + MINOP
        JCL = (NEL-JOP)/2
        JOCC = JOP + JCL
c
        NJCNF = NCNFTP(JTYP,IREFSM)
        IF(K.GE.ICNFB1.AND.K.LE.ICNFB1+NJCNF-1) THEN
          KREL = K - ICNFB1 + 1
          KADD = (KREL-1)*JOCC
          KTYP = JTYP
          CALL ICOPVE(ICONF(ICNFB2+KADD),KCNF,JOCC)
        END IF
c
        ICNFB1 = ICNFB1 + NJCNF
        ICNFB2 = ICNFB2 + NJCNF*JOCC
  100 CONTINUE
c
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' Output fron GETCNF '
        WRITE(LUPRI,*) ' ================== '
        WRITE(LUPRI,*) ' Input configuration number : ', K
        WRITE(LUPRI,*) ' Corresponding type : ', KTYP
        WRITE(LUPRI,*) ' Occupation : '
        NOCC = (NEL+KTYP-1+MINOP)/2
        CALL IWRTMA(KCNF,1,NOCC,1,NOCC)
      END IF
c
      RETURN
      END
C  /* Deck cnfstr */
      SUBROUTINE CNFSTR(ICONF,ITYP,IASTR,IBSTR,NORB,NAEL,NBEL,
     &                  IDET,IPRODT,XNDXCI,ISCR, SIGN)
c
c An orbital configuration ICONF is given,
c Obtain the corresponding alpha strings,IASTR
c        the corresponding beta  strings,IBSTR
c        the corresponding sign array   ,ISIGN
c
c Jeppe Olsen , Summer of '89
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION ICONF(*),XNDXCI(*),ISCR(*),IPRODT(*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*), SIGN(*)
#include <spinfo.h>
c
      NTEST =  0
      NEL = NAEL + NBEL
      IOPEN = ITYP-1+MINOP
      ICLOS = (NAEL + NBEL - IOPEN) / 2
      IOCC  = IOPEN + ICLOS
c
c* Spin orbital occupations of determinants of configuration
c
      KLFREE = 1
      KLDETS = KLFREE
      KLFREE = KLDETS + IDET*(NAEL+NBEL)
c Pointer for determinants of prototype ITYP
      IP = 1
      DO 10 JTYP = 1, ITYP - 1
        IP = IP + NDTFTP(JTYP)*(JTYP-1+MINOP)
   10 CONTINUE
      CALL CNDET(ICONF,IPRODT(IP),IDET,NAEL+NBEL,NORB,
     &           IOPEN,ICLOS,ISCR(KLDETS),NTEST)
c
c* Separate determinants into strings and determine sign change
c
      DO 100 JDET = 1, IDET
        CALL DETSTR(ISCR(KLDETS+(JDET-1)*NEL),
     &              IASTR(1,JDET),IBSTR(1,JDET),
     &              NEL,NAEL,NBEL,NORB,ISIGN      ,ISCR(KLFREE))
      SIGN(JDET) = DFLOAT(ISIGN)
  100 CONTINUE
c
      IF( NTEST .GE. 1 ) THEN
        WRITE(LUPRI,*) ' Output from CNFSTR '
        WRITE(LUPRI,*) ' ================== '
        WRITE(LUPRI,*) ' Input configuration '
        CALL IWRTMA(ICONF,1,IOCC,1,IOCC)
        WRITE(LUPRI,*) ' Corresponding alpha and beta strings'
        CALL IWRTMA(IASTR,NAEL,IDET,NAEL,IDET)
        CALL IWRTMA(IBSTR,NBEL,IDET,NBEL,IDET)
        WRITE(LUPRI,*) ' SIGN ARRAY '
        CALL WRTMAT( SIGN,1,IDET,1,IDET,0)
      END IF
c
      RETURN
      END
C  /* Deck cnhcn */
      SUBROUTINE CNHCN(ICNL,ITPL,ICNR,ITPR,CNHCNM,SCR,LSCR,XNDXCI,NEL,
     &                 NAEL,NBEL,NINOB,ECORE,ONEBOD,IPRODT,DTOC,
     &                 NORB, FIJKL,RJ,RK,ILTSOB)
c
c Obtain Hamilton matrix over CSF's of configurations ICNL,ICNR
c
c Jeppe Olsen , Summer of '89
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION ICNL(*),ICNR(*),XNDXCI(*),SCR(*)
      DIMENSION IPRODT(*),DTOC(*),ONEBOD(*), CNHCNM(*)
      DIMENSION FIJKL(*),RJ(*),RK(*),ILTSOB(*)
#include <spinfo.h>
c Length of SCR : LSCR
      NTEST = 0
c
c* 1 : Obtain determints corresponding to configurations
c
      IOPL = ITPL - 1 + MINOP
      IOPR = ITPR - 1 + MINOP
c
      ICLL = (NEL-IOPL)/2
      ICLR = (NEL-IOPR)/2
c
      NDETL = NDTFTP(ITPL)
      NDETR = NDTFTP(ITPR)
c
      NCSFL = NCSFTP(ITPL)
      NCSFR = NCSFTP(ITPR)
c
      KLFREE = 1
c
      KLDTLA = KLFREE
      KLFREE = KLFREE + NDETL*NAEL
c
      KLDTLB = KLFREE
      KLFREE = KLFREE + NDETL*NBEL
c
      KLISL  = KLFREE
      KLFREE = KLFREE + NDETL
c
      KLDTRA = KLFREE
      KLFREE = KLFREE + NDETR*NAEL
c
      KLDTRB = KLFREE
      KLFREE = KLFREE + NDETR*NBEL
c
      KLISR  = KLFREE
      KLFREE = KLFREE + NDETR
c
      CALL CNFSTR(ICNL,ITPL,SCR(KLDTLA),SCR(KLDTLB),
     &            NORB,NAEL,NBEL,NDETL,IPRODT,XNDXCI,SCR(KLFREE),
     &            SCR(KLISL) )
      CALL CNFSTR(ICNR,ITPR,SCR(KLDTRA),SCR(KLDTRB),
     &            NORB,NAEL,NBEL,NDETR,IPRODT,XNDXCI,SCR(KLFREE),
     &            SCR(KLISR) )
c
c* Hamiltonin matrix over determinants of ths configuration
c
      KLDHD  = KLFREE
      KLFREE = KLFREE + NDETL*NDETR
      LFREE  = LSCR - KLFREE + 1
c
C     CALL DIHDJ(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,
C    &   IWORK,LWORK,NORB,ONEBOD,RJ,RK,FIJKL,ILTSOB,
C    &   HAMIL,ISYM,ECORE,ICOMBI,PSIGN)
      CALL DIHDJ(SCR(KLDTLA),SCR(KLDTLB),NDETL,
     &           SCR(KLDTRA),SCR(KLDTRB),NDETR,
     &           NAEL,NBEL,SCR(KLFREE),LFREE,NORB,ONEBOD,
     &           RJ,RK,FIJKL,ILTSOB,
     &           SCR(KLDHD),0,ECORE,0,1.0D0)
c
c* Transform matrix to CSF basis
c
c* : sign changes
      CALL DGMM2 (SCR(KLDHD),SCR(KLDHD),SCR(KLISL),1,NDETL,NDETR)
      CALL DGMM2 (SCR(KLDHD),SCR(KLDHD),SCR(KLISR),2,NDETL,NDETR)
      IPL = 1
      DO 200 JTYP = 1, ITPL - 1
        JCSF = NCSFTP(JTYP)
        JDET = NDTFTP(JTYP)
        IPL = IPL + JCSF*JDET
  200 CONTINUE
      IF(ITPR.EQ. ITPL) THEN
        IPR = IPL
      ELSE
        IPR = 1
        DO 300 JTYP = 1, ITPR - 1
          JCSF = NCSFTP(JTYP)
          JDET = NDTFTP(JTYP)
          IPR = IPR + JCSF*JDET
  300   CONTINUE
      END IF
c
      KLCHD = KLFREE
      KLFREE = KLFREE + NCSFL*NDETR
      CALL MATML4(SCR(KLCHD),DTOC(IPL),SCR(KLDHD),
     &            NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,1)
      CALL MATML4(CNHCNM,SCR(KLCHD),DTOC(IPR),
     &            NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,0)
c
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*)
     &  ' CSF-Hamiltonian matrix between two configurations'
        CALL WRTMAT(CNHCNM,NCSFL,NCSFR,NCSFL,NCSFR,0)
      END IF
c
      RETURN
      END
C  /* Deck weight */
      SUBROUTINE WEIGHT(Z,IAB,ISCR,IPRINT)
c
c     Purpose:
c
c     Construct vertex weights
c
c     Externals
c
c     RSMXMN,GRAPW
c
#include <implicit.h>
#include <priunit.h>
      INTEGER Z(*), ISCR(*)
c
#include <mxpdim.h>
#include <strnum.h>
#include <ciinfo.h>
c
c IAB = 1 : alpha electrons
c IAB = 2 : beta  electrons
c
      NORB = NORB1 + NORB2 + NORB3
      IF( IAB .EQ. 1 ) THEN
        NELEC = NAEL
        NELECO = NBEL
      ELSE
        NELEC = NBEL
        NELECO = NAEL
      END IF
c
c.   default values for NEL1MX and NEL3MN  ( Includes ALPHA AND BETA)
COLD    NEL1MX = MIN(2*NORB1,NAEL+NBEL)
COLD    NEL3MN = MAX(0,NAEL+NBEL-2*(NORB1+NORB2) )
c. number of electrons for alpha or beta case
        MXEL1  = MIN(NEL1MX,NORB1,NELEC)
        MXEL1O  = MIN(NEL1MX,NORB1,NELECO)
        MNEL1 = MAX(0,NEL1MN-MXEL1O)
        IF( IPRINT.GT.10 ) WRITE(LUPRI,*)' MXEL1,MXEL1O,MNEL1',
     &                                 MXEL1,MXEL1O,MNEL1
c
        MXEL3 = MIN(NEL3MX,NORB3,NELEC)
        MXEL3O = MIN(NEL3MX,NORB3,NELECO)
        MNEL3 = MAX(0,NEL3MN-MXEL3O)
        IF( IPRINT.GT.10 ) WRITE(LUPRI,*)' MXEL3,MXEL3O,MNEL3',
     &                                 MXEL3,MXEL3O,MNEL3
c
        KLFREE = 1
        KLMAX = KLFREE
        KLFREE = KLFREE + NORB
c
        KLMIN = KLFREE
        KLFREE = KLFREE + NORB
c
        KW = KLFREE
        KLFREE = KW + (NELEC+1)*(NORB+1)
c
c.Max and min arrays for strings
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RSMXMN IS ENTERED FROM WEIGHT'
        CALL RSMXMN(ISCR(KLMAX),ISCR(KLMIN),NORB1,NORB2,NORB3,
     &              NELEC,MNEL1,MXEL1,MNEL3,MXEL3,IPRINT)
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RETURN TO WEIGHT FROM RSMXMN'
c. Arc weights
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' GRAPW IS ENTERED FROM WEIGHT'
        CALL GRAPW
     &  (ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NELEC,IPRINT)
      IF( IPRINT.GT.5 ) WRITE(LUPRI,*)' RETURN TO WEIGHT FROM GRAPW'
c
      RETURN
      END
C  /* Deck rsmxmn */
      SUBROUTINE RSMXMN(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     &                  MIN1,MAX1,MIN3,MAX3,IPRINT)
c
c     Author:        J. Olsen, Univ. of Lund, Sweden, April 1987
c
c     Externals:
c
c     IWRTMA
c
c
c     CONSTRUCT ACCUMULATED MAX AND MIN ARRAYS FOR A RAS SPACE
c     OF STRINGS
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION  MINEL(*),MAXEL(*)
c
      NORB = NORB1 + NORB2 + NORB3
      MIN1A = MIN1
      MAX1A = MAX1
c
      MIN2A = NEL - MAX3
      MAX2A = NEL - MIN3
c
      MIN3A = NEL
      MAX3A = NEL
c
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF(NORB1 .GT. 0 )
     &    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF(NORB1+NORB2 .GT. 0 )
     &    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
  100 CONTINUE
c
      IF( IPRINT.GT.20 ) THEN
         WRITE(LUPRI,*) ' MINEL : '
         CALL IWRTMA(MINEL,1,NORB,1,NORB)
         WRITE(LUPRI,*) ' MAXEL : '
         CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
C
      RETURN
      END
C  /* Deck grapw */
      SUBROUTINE GRAPW(W,Y,MINEL,MAXEL,NORB,NEL,IPRINT)
c
c       Author:        J. Olsen, Univ. of Lund, Sweden, April 1987
c
c       Externals:
c
c       ISETVC,IWRTMA
c
c      A GRAPH OF STRINGS HAS BEEN DEFINED FROM
c      MINEL(I) IS THE SMALLEST ALLOWED NUMBER OF ELECTRONS IN
c      ORBITALS 1 THROUGH I
c      MAXEL(I) IS THE LARGEST ALLOWED  NUMBER OF ELECTRONS IN
c      ORBITALS 1 THROUGH I
c
c      SET UP VERTEX WEIGHTS W
c      SET UP ARC WEIGHTS    Y
c
c      REVERSE LEXICAL ORDEREING IS USED WITH
c      WEIGHTS OF UNOOCUPIED ORBITALS SETR TO ZERO
c
#include <implicit.h>
#include <priunit.h>
      INTEGER W(NORB+1,NEL+1)
      INTEGER Y(NORB,NEL)
      INTEGER MAXEL(NORB),MINEL(NORB)
c
      CALL ISETVC(W,0,(NEL+1)*(NORB+1) )
      CALL ISETVC(Y,0,NEL*NORB)
c. VERTEX WEIGHTS
c
c WEIGHT FOR VETEX (IEL,IORB) IS STORED IN W(IORB+1,IEL+1)
c
c HEAD
      W(1,1) = 1
      DO 300 IEL = 0, NEL
        DO 200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            IF( IEL .GT. 0 ) THEN
              W(IORB+1,IEL+1) = W(IORB-1+1,IEL+1)
     &                        + W(IORB-1+1,IEL-1+1)
            ELSE
              W(IORB+1,1) = W(IORB-1+1,1)
            END IF
          END IF
  200   CONTINUE
  300 CONTINUE
c
      IF( IPRINT.GT.20 ) THEN
        WRITE(LUPRI,'(A)') ' MATRIX FOR VERTEX WEIGHTS '
        CALL IWRTMA(W,NORB+1,NEL+1,NORB+1,NEL+1)
       END IF
c
c. ARC  WEIGHTS
c
c WEIGHT FOR ARC CONNECTING VERTICES (IORB-1,IEL-1) AND (IORB,IEL)
c IS STORED IN Y(IORB,IEL)
c Y(IORB,IEL) = W(IORB-1,IEL)
      DO 1300 IEL = 1, NEL
        DO 1200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            Y(IORB,IEL) = W(IORB-1+1,IEL+1)
          END IF
 1200   CONTINUE
 1300 CONTINUE
c
      IF( IPRINT.GT.20 ) THEN
        WRITE(LUPRI,'(A)') ' Matrix for arc weights  '
        CALL IWRTMA(Y,NORB,NEL,NORB,NEL)
      END IF
c
      RETURN
      END
C  /* Deck detgn4 */
      SUBROUTINE DETGN4(NEL,NORB1,NORB2,NORB3,
     &                  NELMN1,NELMX1,NELMN3,NELMX3,Z,IREORD,
     &                  MAXSTR,ISTR,IOC,NORB,IPRINT)
c
c     Purpose:
c
c     Turbo routine for generating strings that fulfills RAS constraints
c
c     Externals
c
c     ISTVC2,NXTORD,IWRTMA,ICOPVE,ISCAVE
c
c     Jeppe Olsen , sometime, someplace  in 1989
c
#include <implicit.h>
#include <priunit.h>
c.Input
      INTEGER Z(*), IREORD(*), ISTR(NEL,MAXSTR)
c.Scratch
      INTEGER IOC(*)
c
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
c Loop over possible partitionings between RAS1,RAS2,RAS3
c910416-hjaaj: new IEL1MX and IEL3MX (was IELaMX = NELMXa)
      IEL1MX = MIN(NORB1,NELMX1+1,NEL)
      DO 1001 IEL1 = IEL1MX,NELMN1,-1
       IEL3MX = MIN(NORB3,NELMX3,NEL-IEL1)
      DO 1003 IEL3 = NELMN3,IEL3MX, 1
       IEL2 = NEL - IEL1-IEL3
       IF( IPRINT.GT.30 )WRITE(LUPRI,*) ' IEL1 IEL2 IEL3 ',
     &      IEL1,IEL2,IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
c Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1,IPRINT)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( IPRINT.GT.30 ) THEN
           WRITE(LUPRI,*) ' RAS 1 STRING '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         ENDIF
         IFRST2 = 1
         IFRST3 = 1
c Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,
     &         IORB2F,IORB2L,NONEW2,IPRINT)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
         IF( IPRINT.GT.30 ) THEN
           WRITE(LUPRI,*) ' RAS 1 2 STRING '
           CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
         ENDIF
           IFRST3 = 1
c Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3,IPRINT)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
         IF( IPRINT.GT.30 ) THEN
           WRITE(LUPRI,*) ' RAS 1 2 3 STRING '
           CALL IWRTMA(IOC,1,NEL,1,NEL)
         ENDIF
c Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
             IADRES = IZNUM( IOC ,NEL,Z,NEL,NORB,IREORD,IPRINT)
             IF(IADRES .GT. MAXSTR) THEN
               WRITE(LUPRI,*)
     &         ' TOO MANY STRINGS IN DETGN4, STRING NUMBER ENCOUNTERED',
     &         IADRES
               WRITE(LUPRI,*) ' MAX NUMBER OF STRINGS :',MAXSTR
               CALL QUIT('Too many strings in DETGN4')
             END IF
             IREORD(IADRES) = -NSTRIN
             IF (NEL .GT. 0) CALL ICOPVE(IOC(1),ISTR(1,NSTRIN),NEL)
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL2 .NE. 0 ) GOTO 902
           IF( IEL1 .NE. 0 ) GOTO 901
 1003 CONTINUE
 1001 CONTINUE
C
c* ZERO FORBIDDEN HANDY NUMBERS AND MAKE ALLOWED HANDY NUMBERS POSITIVE
      DO 200 ISTIN = 1,MAXSTR
       IF( IREORD(ISTIN) .GT. 0 ) IREORD(ISTIN) = 0
  200 CONTINUE
       CALL ISCAVE(IREORD(1),-1,MAXSTR)
C
      IF( IPRINT.GT.20 ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' *** Output from DETGN4 ***'
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' HANDY NUMBER TO PACKED NUMBER '
        CALL IWRTMA(IREORD,1,MAXSTR,1,MAXSTR)
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' NUMBER OF STRINGS GENERATED   ', NSTRIN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' OCCUPATION OF STRINGS '
        WRITE(LUPRI,*) '======================='
        MAXPR = 100
        NPR  = MIN(NSTRIN,MAXPR)
        IF( NSTRIN.GT.NPR ) WRITE(LUPRI,*)' ONLY FIRST',NPR,' STRINGS',
     &  ' ARE PRINTED'
        WRITE(LUPRI,*) ' '
        WRITE(LUPRI,*) ' STRING          OCCUPATION '
        DO 1234 ISTRIN = 1,NPR
          WRITE(LUPRI,'(I4,4X,(20I3))')
     &    ISTRIN,(ISTR(IEL,ISTRIN),IEL = 1,NEL)
1234    CONTINUE
      ENDIF
C
      RETURN
      END
C  /* Deck nxtord */
      SUBROUTINE NXTORD(INUM,NELMNT,MINVAL,MAXVAL,NONEW,IPRINT)
c
c
c     Purpose
c
c     AN ORDERED SET OF NUMBERS INUM(I),I=1,NELMNT IS
c        GIVEN IN STRICTLY ASCENDING ORDER.
c     VALUES OF INUM(*) IS RESTRICTED TO THE INTERVAL MINVAL,MAXVAL .
c     FIND NEXT HIGHER NUMBER.
c     NONEW = 1 ON RETURN INDICATES THAT NO ADDITIONAL NUMBERS
c        COULD BE OBTAINED.
c
c     Externals
c
c     IWRTMA
c
c Jeppe Olsen , 1990
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION INUM(*)
c
       IF( IPRINT.GT.30 ) THEN
         WRITE(LUPRI,*) ' Initial number in NXTORD '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
c
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF( IPLACE .LT. NELMNT .AND.
     &      INUM(IPLACE)+1 .LT. INUM(IPLACE+1)
     &  .OR.IPLACE.EQ. NELMNT .AND.
     &      INUM(IPLACE)+1.LE.MAXVAL) THEN
              INUM(IPLACE) = INUM(IPLACE) + 1
              NONEW = 0
              GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
              IF(IPLACE .EQ. 1 ) THEN
                INUM(IPLACE) = MINVAL
              ELSE
                INUM(IPLACE) = INUM(IPLACE-1) + 1
              END IF
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
              NONEW = 1
              GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
c
      IF( IPRINT.GT.30 ) THEN
        WRITE(LUPRI,*) ' New number in NXTORD:'
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
c
      RETURN
      END
C  /* Deck csfdim */
      SUBROUTINE CSFDIM(NACTOB,NACTEL,MULTP,MAXSYM,SYMPRO,ORBSYM,IWORK,
     &                  LLCSF,LCSFFO,KFREE,NCNSM,IPRINT)
c
c
c
c     Initializing routine for CSF-DET expansions
c     Set up common block /CSFDIM/ which gives information
c     on the number of dets,csf's for each symmetry
c
c     find local memory requirements for CSF routines
c     Largest local memory requirements in CNFORD,CSFDET is
c     returned in LLCSF
c
c DETERMINE BASE ADRESSES /CSFBAS/ (STARTING FROM KFREE ) FOR :
c             KDFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
c             KCFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
c             KDTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
c             KICONF(I) : SPACE FOR STORING  NCNSM
c                        CONFIGURATION EXPANSIONS
c Local memory requirements : IWORK(NACTOB)
c
c     Externals
c
c     CISIZ2,MEMADD,IWRTMA,IWEYLF
c
c Jeppe Olsen, 1989
c
c Note . In this version, MS2 is assumed to be a input parameter
c
#include <implicit.h>
#include <priunit.h>
c
      INTEGER SYMPRO(8,8)
      INTEGER ORBSYM(NACTOB)
      DIMENSION IWORK(*)
c
#include <ciinfo.h>
#include <spinfo.h>
#include <csfbas.h>
c
      KFREEO = KFREE
c
c.. Define parameters in SPINFO
c
      MULTS = MULTP
      NEL = NACTEL
c
c. ALLOWED NUMBER OF OPEN ORBITALS
      MINOP = ABS(MS2)
      MAXOP = MIN(NACTOB,2*NACTOB-NACTEL,NACTEL)
c For RAS must RAS constraints be checked to get
c the correct max number of open shells
      MAXOP = 0
      NORB1 = IORB1L-IORB1F+1
      NORB3 = IORB3L-IORB3F+1
      NORB2 = NACTOB-NORB1-NORB3
c
      DO 5 IEL1 = NEL1MN,2*NORB1
        DO 4 IEL3 = 0,NEL3MX
          IEL2 = NACTEL-IEL1-IEL3
          IF(IEL2 .LT. 0 ) GOTO 4
          IOP1 = MIN(NORB1,2*NORB1-IEL1,IEL1)
          IOP2 = MIN(NORB2,2*NORB2-IEL2,IEL2)
          IOP3 = MIN(NORB3,2*NORB3-IEL3,IEL3)
          IOP = IOP1 + IOP2 + IOP3
          MAXOP = MAX(MAXOP,IOP)
    4   CONTINUE
    5 CONTINUE
      IF( IPRINT.GT.10 ) THEN
         WRITE(LUPRI,*) ' MAXOP WITH RAS CONSTRAINTS :' ,MAXOP
      ENDIF
      NTYP = MAXOP-MINOP + 1
c
      IF( NTYP .GT. MTYP ) THEN
        WRITE(LUPRI,*) '  NUMBER OF CONFIGURATION TYPES TOO LARGE'
        WRITE(LUPRI,*) '  CHANGE PARAMETER MTYP TO AT LEAST ',NTYP
        WRITE(LUPRI,*) '  CURRENT VALUE OF MTYP ',MTYP
        CALL QUIT('MTYP IN LUSPIN TOO SMALL')
      END IF
c
c.. Number of sd's and csf's per configuration type
c
      DO 10 ITP = 1,NTYP
        IOPEN = MINOP+ITP - 1
        IAEL = (IOPEN + MS2 ) / 2
        IBEL = (IOPEN - MS2 ) / 2
        IF(IAEL+IBEL .EQ. IOPEN ) THEN
          NDTFTP(ITP) = IBION(IOPEN,IAEL)
          IF(IOPEN .GE. MULTS-1) THEN
            NCSFTP(ITP) = IWEYLF(IOPEN,MULTS)
          ELSE
            NCSFTP(ITP) = 0
          END IF
        ELSE
          NDTFTP(ITP) = 0
          NCSFTP(ITP) = 0
        END IF
   10 CONTINUE
      IF( IPRINT.GT.5 ) THEN
      WRITE(LUPRI,'(/A)') ' Information about prototype configurations '
      WRITE(LUPRI,'( A)') ' ========================================== '
      WRITE(LUPRI,'(/A)')
     &'  Open orbitals   Determinants    CSFs '
      DO 580 IOPEN = MINOP,MAXOP,2
        ITYPE = IOPEN - MINOP + 1
        WRITE(LUPRI,'(5X,I3,I16,I13)')
     &  IOPEN,NDTFTP(ITYPE),NCSFTP(ITYPE)
  580 CONTINUE
      ENDIF
c
c.. Number of determinants and CSF's for each reference symmetry
c
      KL1 = 1
      KL2 = KL1 + NACTOB
      KL3 = KL2 + NACTOB
      CALL CISIZ2(NACTOB,ORBSYM,SYMPRO,IWORK(1),MAXSYM,NEL,
     &            IWORK(KL1),IWORK(KL2),IWORK(KL3),IPRINT)
c
c.. Permanent and local memory for csf routines
c
c    memory for CSDTMT arrays.
c    Largest block of proto type determinants .
c    Largest number of prototype determinants
c    All configurations( of any specific symmetry )
c
      LIDT = 0
      LICS = 0
      LDTOC = 0
      MXPTBL = 0
      MXDT = 0
      LCONF = 0
      DO 11 ITP = 1,NTYP
        IOPEN = MINOP+ITP - 1
        LIDT = LIDT + NDTFTP(ITP) * IOPEN
        LICS = LICS + NCSFTP(ITP) * IOPEN
        LDTOC= LDTOC + NCSFTP(ITP)*NDTFTP(ITP)
        MXDT =   MAX(MXDT,NDTFTP(ITP) )
        MXPTBL = MAX(NDTFTP(ITP)*IOPEN,MXPTBL)
   11 CONTINUE
c. local memory for CSFDET
      LCSFDT = MXPTBL + MAXOP
c. local memory for CNFORD
      LCNFOR = MAX(2*NTYP+3*NACTOB,(MXDT+2)*NACTEL)
c. local memory for any routine used in construction of csf basis
      LLCSF = MAX(LCSFDT,LCNFOR)
c. Memory needed to store ICONF array
      LCONF = 0
      LDET = 0
      DO 30 ISYM = 1,MAXSYM
        LLCONF = 0
        LDET = MAX(LDET,NDTASM(ISYM))
        DO 25 ITYP = 1, NTYP
          IOPEN = ITYP+MINOP-1
          ICL = ( NEL-IOPEN)/2
          LLCONF = LLCONF + NCNFTP(ITYP,ISYM)*(IOPEN+ICL)
   25   CONTINUE
        LCONF = MAX(LCONF,LLCONF)
   30 CONTINUE
      IF( IPRINT.GT.10 ) THEN
       WRITE(LUPRI,'(/A,I12)')
     & '  Memory for holding largest list of configurations ',LCONF
       WRITE(LUPRI,'(/A,I12)')
     & '  Size of largest CI expansion (dets)',LDET
      ENDIF
c
c. permanent memory for csf proto type arrays
c
      CALL MEMADD(KDFTP,LIDT,KFREE,1)
      CALL MEMADD(KCFTP,LICS,KFREE,1)
      CALL MEMADD(KDTOC,LDTOC,KFREE,2)
      LPROTO = KFREE - KDFTP
      IF( IPRINT.GT.10 ) THEN
        WRITE(LUPRI,*)'  MEMORY ALLOCATION IN CSFDIM : '
        WRITE(LUPRI,*)'  POINTERS FOR PROTOTYPE INFORMATION '
        WRITE(LUPRI,*)'    KDFTP,  KCFTP,  KDTOC'
        WRITE(LUPRI,'(1X,3I12)')KDFTP,KCFTP,KDTOC
      END IF
c
c. PERMANENT ARRAYS FOR
c. HOLDING CONFIGURATION EXPANSIONS AND REORDER ARRAYS
c  ( sign and reorder arrays are assumed to be written in a single array)
      IF(NCNSM .GT. MXCNSM ) THEN
        WRITE(LUPRI,'(A,2I5)')
     &  '  TROUBLE IN CSFDIM NCNSM > MXCNSM : NCNSM,MXCNSM',
     &  NCNSM,MXCNSM
        CALL QUIT('CSFDIM : NCNSM  IS GREATER THAN MXCNSM')
      END IF
      DO 60 ICNSM = 1, NCNSM
        CALL MEMADD(KICONF(ICNSM),LCONF,KFREE,1)
        CALL MEMADD(KICTS(ICNSM),LDET ,KFREE,1)
   60 CONTINUE
C
      IF( IPRINT.GT.20 ) THEN
        WRITE(LUPRI,*)'   KICONF'
        CALL IWRTMA(KICONF,1,NCNSM,1,NCNSM)
        WRITE(LUPRI,*)'   KICTS '
        CALL IWRTMA(KICTS,1,NCNSM,1,NCNSM)
      END IF
C
      LCSFFO = KFREE - KFREEO
      IF( IPRINT.GT.10 ) THEN
        WRITE(LUPRI,'(/A,I12)')
     &  '  Memory needed for prototype information',LPROTO
        WRITE(LUPRI,'(/A,I12)')
     &  '  Memory needed for using CSFs  ',LCSFFO
      END IF
C
      RETURN
      END
C  /* Deck cisiz2 */
      SUBROUTINE CISIZ2 (NORB,ORBSYM,SYMPRO,IWORK,MAXSYM,NEL,
     &                   IICL,IIOP,IIOC,IPRINT)
c
c Obtain :
c
c     Number of CSF's and SD's
c     Information about all symmetries generated
c
c     Externals
c
c     ISETVC,IWRTMA
c
c     Jeppe Olsen , 1989
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION IWORK(*)
      INTEGER SYMPRO(8,8)
      INTEGER ORBSYM(NORB)
      DIMENSION IICL(*),IIOP(*),IIOC(*)
c
#include <spinfo.h>
#include <ciinfo.h>
c
      CALL ISETVC(NCNFTP,0,MAXSYM*MTYP)
      CALL ISETVC(NCSASM,0,MAXSYM)
      CALL ISETVC(NDTASM,0,MAXSYM)
      CALL ISETVC(NCNASM,0,MAXSYM)
c
      ILOOP = 0
c Min number of doubly occupied orbitals in RAS 1
      NORB1 = IORB1L-IORB1F+1
      MINCL1 = MAX(0,NEL1MN-NORB1)
      IF( IPRINT.GT.20 ) THEN
        WRITE(LUPRI,*)
     &  ' MIN NUMBER OF DOUBLY OCCUPIED ORBITALS IN RAS 1',MINCL1
      ENDIF
      DO 5000 NOP = MINOP,MAXOP,2
        ITYPE = NOP-MINOP+1
        NCL = (NEL-NOP)/2
        IF( IPRINT.GT.20 )
     &  WRITE(LUPRI,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
c
c. first combination of double occupied orbitals
        CALL ISETVC(IIOC,0,NORB)
        DO 10 ICL = 1, NCL
          IICL(ICL) = ICL
          IIOC(ICL) = 2
   10   CONTINUE
        IFRSTC = 1
c.  Loop over double occupied orbital configurations
 2000   CONTINUE
c
c. next double occupied configuration
          IF ( .NOT. ( IFRSTC .EQ. 0 .AND. NCL .NE. 0 ) ) GO TO 810
c
          DO 50 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 1 ) IIOC(IORB) = 0
   50     CONTINUE
c
          IPLACE = 0
  800     IPLACE = IPLACE + 1

          IPRORB = IICL(IPLACE)
          IIOC(IPRORB) = 0
          NEWORB = IPRORB+1
          IF((IPLACE .LT. NCL .AND. NEWORB .LT. IICL(IPLACE+1))
     &      .OR.
     &      IPLACE .EQ. NCL .AND.  NEWORB.LE. NORB ) THEN
            IICL(IPLACE) = NEWORB
            IIOC(NEWORB) = 2
          ELSE IF
     &    (.NOT.(IPLACE.EQ.NCL.AND.NEWORB.GE.NORB)) THEN
            IF(IPLACE .EQ. 1 ) THEN
              IICL(1) = 1
              IIOC(1) = 2
            ELSE
              IICL(IPLACE) = IICL(IPLACE-1) + 1
              IIOC(IICL(IPLACE)) = 2
            END IF
            GOTO 800
          ELSE
c. No more inactive configurations
             GOTO 2001
          END IF
  810     IFRSTC = 0
          IF( IPRINT.GT.30 ) THEN
            WRITE(LUPRI,'(''  next inac. conf.'',20I3,/,(18X,20I3))')
     &      (IICL(ICFUE),ICFUE=1,NCL)
          END IF
c.        CHECK RAS1 and RAS 3
          IEL1C = 0
          IEL3C = 0
          ICL1  = 0
          DO  20 ICL = 1,NCL
            IORB = IICL(ICL)
            IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
               IEL1C = IEL1C + 2
               ICL1 = ICL1 + 1
            ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
               IEL3C = IEL3C + 2
            END IF
   20     CONTINUE
          IIICHK = 1
          IF(ICL1 .LT.MINCL1.AND. IIICHK.EQ.1) THEN
c Next higher combination with a higher number of inactive orbitals
            DO 41 ICL = 1,ICL1+1
              IIOC(IICL(ICL)) = 0
              IICL(ICL) = ICL
              IIOC(ICL) = 2
   41       CONTINUE
            IPLACE=ICL1+1
            IF( IPLACE.GE.NCL) GOTO 2001
            GOTO 800
          END IF
          IF(IEL3C.GT.NEL3MX) GOTO 2000
c. Highest orbital not occupied
         MXMPTY = NORB
         IORB = NORB+1
c. begin while
   12      CONTINUE
           IORB = IORB - 1
           IF(IIOC(IORB) .EQ. 2 ) THEN
             MXMPTY = IORB-1
             IF( IORB .NE. 1 ) GOTO 12
           END IF
c. End while
C
c. first active configuration
          IORB = 0
          IOP = 0
          DO 30 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 0 ) THEN
              IOP = IOP + 1
              IF( IOP .GT. NOP ) GOTO 31
              IIOC(IORB) = 1
              IIOP(IOP) = IORB
            END IF
   30     CONTINUE
   31     CONTINUE
          IFRSTO = 1
c
c. Next open shell configuration
 1000     CONTINUE
          IF (.NOT. (IFRSTO .EQ. 0 .AND. NOP .NE. 0 ) ) GO TO 710
            IPLACE = 0
  700       CONTINUE
              IPLACE = IPLACE + 1
              IPRORB = IIOP(IPLACE)
              NEWORB = IPRORB + 1
              IIOC(IPRORB) = 0
  690         CONTINUE
                IF(IIOC(NEWORB) .NE. 0 ) THEN
                   NEWORB = NEWORB + 1
                   IF(NEWORB .LE. MXMPTY ) GOTO 690
                END IF
c 691         End of loop
            IF(IPLACE .LT. NOP .AND. NEWORB .LT. IIOP(IPLACE+1) .OR.
     &         IPLACE .EQ. NOP .AND. NEWORB .LE. MXMPTY ) THEN
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
            ELSE IF (  IPLACE .NE. NOP ) THEN
              IF(IPLACE.EQ.1) THEN
                NEWORB = 1 - 1
              ELSE
                NEWORB = IIOP(IPLACE-1)
              END IF
 671          CONTINUE
                NEWORB = NEWORB + 1
              IF(IIOC(NEWORB) .NE. 0 .AND. NEWORB.LE.MXMPTY ) GOTO 671
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
              GOTO 700
            ELSE
c. No more active configurations , so
              IF( NCL .NE. 0 ) GOTO 2000
              IF( NCL .EQ. 0 ) GOTO 5001
            END IF
  710     IFRSTO = 0
          IF( IPRINT.GT.30 ) THEN
            WRITE(LUPRI,'(''  next active conf.'',20I3,/,(19X,20I3))')
     &      (IICL(ICFUE),ICFUE=1,NCL)
          ENDIF
c        RAS  CONSTRAINTS
         IOKAY=1
         IEL1 = IEL1C
         IEL3 = IEL3C
c.       CHECK RAS1 and RAS3
         DO  40 IOP = 1,NOP
           IORB = IIOP(IOP)
           IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
              IEL1 = IEL1 + 1
           ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
              IEL3 = IEL3 + 1
           END IF
   40    CONTINUE
           IR3CHK  = 1
           IF(IEL3.GT.NEL3MX.AND.IR3CHK.EQ.1) THEN
c. Direct generation of next string with a lower number of
c. electrons in RAS III
c. Procedure : Obtain substring that contains on electron
c  in RAS I, modify this to lowest possible string.
c. Add one to remaining string
c
c. Number of electrons in substring
             IFSTR3 = 0
             DO 5607 IOP = 1, NOP
               IF(IIOP(IOP).GE.IORB3F) THEN
                IFSTR3 = IOP
                GOTO 5608
               END IF
 5607        CONTINUE
 5608        CONTINUE
             IF(IFSTR3.NE.NOP) THEN

c. Lowest possible string with NOP electrons
               DO 5610 K = 1, IFSTR3
                 IIOC(IIOP(K)) = 0
 5610          CONTINUE
c
               KEL = 0
               KORB = 0
 5630          CONTINUE
                 KORB = KORB + 1
                 IF(IIOC(KORB).NE.2) THEN
                   KEL = KEL + 1
                   IIOC(KORB) = 1
                   IIOP(KEL) = KORB
                 END IF
               IF(KEL.NE.IFSTR3) GOTO  5630
               IPLACE = IFSTR3
               GOTO 700
             END IF
           END IF
         IF( IEL1 .LT. NEL1MN .OR. IEL1 .GT. NEL1MX   .OR.
     &       IEL3 .LT. NEL3MN .OR. IEL3 .GT. NEL3MX ) GOTO  999
c. Spatial symmetry
         ISYM = 1
         DO 920 IOP = 1, NOP
           ISYM = SYMPRO(ISYM,ORBSYM(IIOP(IOP)))
  920    CONTINUE
C
         IF( ISYM.GT.MAXSYM .OR. ITYPE.GT.MTYP ) THEN
            WRITE(LUPRI,*)'CISIZ2 ERROR, ISYM.gt.MAXSYM .or. '//
     &                    'ITYPE.gt.MTYP'
            WRITE(LUPRI,*)'  ISYM,ITYPE',ISYM,ITYPE
            WRITE(LUPRI,*)'  MAXSYM,MTYP',MAXSYM,MTYP
            CALL QUIT('ENFORCE PROG. INTERRUPT IN SUBR. CISIZ2')
         ENDIF
CHJ      NCNF=NCNF+1
         NCNASM(ISYM) = NCNASM(ISYM)+1
         NCNFTP(ITYPE,ISYM)=NCNFTP(ITYPE,ISYM)+1
         IF( IPRINT.GT.30 ) THEN
CHJ         inserted next line 901031/hjaaj
            NCNF=NCNFTP(ITYPE,ISYM)
CHJ         ... 901031/hjaaj: NCNF was uninitialized, I guess
CHJ             it should be NCNFTP(ITYPE,ISYM)
            WRITE(LUPRI,*) '  ACTIVE CNF: NCNF,ISYM,TYPE',NCNF,ISYM,
     &                     ITYPE
            WRITE(LUPRI,1311) (IIOC(I),I=1,NORB)
         ENDIF
 1311    FORMAT('  occupation',20I2,/,(12X,20I2))
c
c** LOOP OVER CONFIGURATIONS, end
c
  999  CONTINUE
       ILOOP = ILOOP + 1
c
       IF( NOP .EQ. 0 .AND. NCL .EQ. 0 ) GOTO 5001
       IF( NOP .EQ. 0 ) GOTO 2000
       GOTO 1000
 2001 CONTINUE
 5000 CONTINUE
 5001 CONTINUE
c
c*. TOTAL NUMBER OF CSF'S AND DETERMINANTS
c
      DO 610 ISYM = 1, MAXSYM
        DO 600 ITYPE = 1, NTYP
c.Dets
          NDTASM(ISYM) = NDTASM(ISYM)
     &    + NDTFTP(ITYPE)*NCNFTP(ITYPE,ISYM)
c.CSF's
          NCSASM(ISYM) = NCSASM(ISYM)
     &    + NCSFTP(ITYPE)*NCNFTP(ITYPE,ISYM)
  600   CONTINUE
  610 CONTINUE
c
      IF( IPRINT.GT.5 ) THEN
        WRITE(LUPRI,*)
        ICSM = 0
        WRITE(LUPRI,'(/A)') ' INFORMATION ABOUT ACTUAL CONFIGURATIONS '
        WRITE(LUPRI,'( A)') ' ========================================'
        WRITE(LUPRI,'(/A)')
     &  '    SYMMETRY     CONFIGURATIONS     CSFS     DETERMINANTS  '
        WRITE(LUPRI,'(A)')
     &  '  =============  ============== ============ ============  '
        DO 570 ICSM = 1, MAXSYM
            WRITE(LUPRI,'(4X,I3,6X,I8,6X,I8,6X,I8)')
     &      ICSM,NCNASM(ICSM),NCSASM(ICSM),NDTASM(ICSM)
  570   CONTINUE
      ENDIF
c
      RETURN
      END
C  /* Deck cisig9 */
      SUBROUTINE CISIG9(TAIJ,TATO,NAEL,NASTR,NAEXCI,
     &                  TBIJ,TBTO,NBEL,NBSTR,NBEXCI,
     &                  NORB,C,HC,NCDET,NHCDET,NSTASA,NSTASB,
     &                  ICOFF,IHCOFF,MAXSYM,DISSYM,
     &                  ISTBAA,ISTBAB,TASYM,TBSYM,G,IJKL,
     &                  ICSM,IHCSM,PERMSM,
     &                  IOCOC,NOCTPA,NOCTPB,ICSO,IHCSO,ICOOS,
     &                  IHCOOS,NSSOA,ISSOA,NSSOB,ISSOB,
     &                  KLTP,ICREA,IANNI,
     &                  VEC2,FACTOR,C2,INDX1,INDX2,L,R,PSIGN,FAC2,
     &                  ITYPSA,ITYPSB,NCDTAS,NHDTAS,
     &                  KLCAN,ITPFOB,HCDUM,IDOH2,ECORE,TUVX,WIJKL,
     &                  IPRINT,ISTINA,ISTINB,MXSTAS,
     &                  ILTSOB,ISTLOB,ISPIN1,ISPIN2,IH8SM,INTFRM,
     &                  INDEX2,F3)
c
c    PURPOSE: CONSTRUCT SIGMA VECTOR
c
c     JEPPE OLSEN , OCT-89 VERSION OF RAS CI ROUTINE
c     SIRIUS version
c
c Version where lower IATP.GT.IBTP always is imposed
c before IASM .GT. IBSM
#include <implicit.h>
#include <cbgetdis.h>
#include <priunit.h>
      EXTERNAL BDJEP
#include <tstjep.h>
c
#include <mxpdim.h>
#include <maxash.h>
#include <mxsmob.h>
c
      LOGICAL PERMSM  , IH8SM
      INTEGER TAIJ,TATO,TASYM ,SYMPRO(8,8),DISSYM
      INTEGER TBIJ,TBTO,TBSYM
      REAL*8  IJKL
      DIMENSION TAIJ(*),TATO(*)
      DIMENSION TBIJ(*),TBTO(*)
      DIMENSION C(NCDET),HC(NHCDET)
      DIMENSION NSTASA(*),NSTASB(*),ICOFF(*),IHCOFF(*)
      DIMENSION G(NORB,NORB),IJKL(NORB*NORB,*),DISSYM(NORB*NORB)
      DIMENSION ISTBAA(* ),ISTBAB(* )
      DIMENSION TASYM(MAXSYM+1,NASTR) ,TBSYM(MAXSYM+1,NBSTR )
c
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ICSO(NOCTPA,MAXSYM),IHCSO(NOCTPA,MAXSYM)
      DIMENSION NSSOA(NOCTPA,MAXSYM),ISSOA(NOCTPA,MAXSYM)
      DIMENSION NSSOB(NOCTPB,MAXSYM),ISSOB(NOCTPB,MAXSYM)
      DIMENSION  ICOOS(NOCTPB,NOCTPA,MAXSYM)
      DIMENSION IHCOOS(NOCTPB,NOCTPA,MAXSYM)
      DIMENSION HCDUM(NHCDET)
      DIMENSION KLTP(*),ICREA(*),IANNI(*)
      DIMENSION VEC2(NBSTR),FACTOR(MXSTAS,*),C2(*)
      DIMENSION FAC2(MXSTAS,*)
      INTEGER   INDX1(*),INDX2(*),L(MXSTAS,*),R(MXSTAS,*)
      INTEGER   ITYPSA(NASTR),ITYPSB(NBSTR)
      DIMENSION NCDTAS(*),NHDTAS(*)
      DIMENSION KLCAN(*),ITPFOB(NORB)
      DIMENSION TUVX(*)
      INTEGER   INDEX2(MXPST,MXPTP),NJBFTP(MXPTP),NJAFTP(MXPTP)
      DIMENSION F3(MXPST,MXPTP)
      INTEGER   INDX4(MXPST), INDX5(MXPST),INDX6(MXPST)
      DIMENSION IFTP(MXPTP),NFTP(MXPTP),IC2OOS(MXPTP,MXPTP)
      DIMENSION WIJKL(*)
c.added
CHJ: The MAXASH below is max for NORB.
      DIMENSION IFTP2(MXPTP),IOST(2*MXPST)
      DIMENSION LORB(MAXASH),NKLSTA(MAXASH)
      DIMENSION NFTPA(MXPTP,MXSMOB),IFTPA(MXPTP,MXSMOB),
     &          INDX7(MXPTP,MXSMOB)
c
      DIMENSION ILTSOB(*),ISTLOB(*)
      DIMENSION ISTINA(NAEL,*),ISTINB(NBEL,*)
c
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
      IF(ISPIN1+ISPIN2.EQ.1) THEN
        CALL SCALVE(HC,-1.0D0,NHCDET)
      END IF
      SIGNLL = (-1.0D0)**(ISPIN1+ISPIN2)
c
c.Largest number of strings in a given SO class
      IMXSOA = IMXVEC(NSSOA,NOCTPA*MAXSYM)
      IMXSOB = IMXVEC(NSSOB,NOCTPB*MAXSYM)
      IMXSO = MAX(IMXSOA,IMXSOB)
c
      IF(ITSTLP.EQ.1) THEN
         WRITE(LUPRI,*)
     &'  CISIG9 : Notice only alpha-alpha+beta-beta loops'
      END IF
      IF(ITSTLP.EQ.2) THEN
         WRITE(LUPRI,*)
     &'  CISIG9 : Notice only alpha-beta loops'
      END IF
      NERR = 0
      IF(IMXSO.GT.MXPST) THEN
        WRITE(LUPRI,*)
     &  ' >> MXPST in CISIG9 insufficient,current value  <<',MXPST
        WRITE(LUPRI,*)
     &  ' >> Raise MXPST in CISIG9 to atleast << ',IMXSO
        NERR = NERR + 1
      END IF
      IF(NOCTPA.GT.MXPTP.OR.NOCTPB.GT.MXPTP) THEN
        WRITE(LUPRI,*)
     &  ' >> MXPTP in CISIG9 insufficient,current value  <<',MXPTP
        WRITE(LUPRI,*)
     &  ' >> Raise MXPTP in CISIG9 to atleast <<' ,MAX(NOCTPA,NOCTPB)
        NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) CALL QUIT('ERRORS in CISIG9')
c
      TEST = 1.0D-12
      IF( IDOH2 .NE. 0 ) THEN
c
c* 1 : DISTRIBUTE INTEGRALS TO G MATRIX
c
c                      G(I,J) =
c                     (G(J,I))
c     SUM(K,I.GE.K,IF I.EQ.K,K .GE. J) (IK|KJ)/(1+D(I,K)*D(J,K))
c
      DO 110 KK = 1, NORB
      DO 110 LL = 1, NORB
        CALL GETIN2(IJKL,LL,KK,TUVX,ILTSOB,WIJKL)
C            GETINT( RIJKL, K1, L1, FIJKL, ILTSOB, WIJKL )
        DO 100 I = KK+1 , NORB
          G(I,LL) = G(I,LL) - IJKL((I-1)*NORB+KK,1 )
  100   CONTINUE
        IF(KK.GT.LL) G(KK,LL) = G(KK,LL) - IJKL((KK-1)*NORB+KK,1)
        IF(KK.EQ.LL) G(KK,KK) = G(KK,KK) - IJKL((KK-1)*NORB+KK,1)/2.0D0
  110 CONTINUE
      END IF
c
COLD  CALL SETVEC(HC,0.0D0,NHCDET)
      IF (PERMSM) CALL SCALVE(HC,0.5D0,NHCDET)
c  SYMMETRY OF H ( AS EXCITATION OPERATOR )
      IHSM = SYMPRO(ICSM,IHCSM)
c
c* RAS CI NEW LOOP 2
c
c**************************************************************
c*                                                            *
c**            Term 1 to sigma vector                         *
c*                                                            *
c**************************************************************
c
c   SIGMA(IA,IB) = SUM((IJ).GE.(KL)) <JB|EB(IJ)EB(KL)|IB>
c                                   *    C(IA,JB)
c                                   * (J I | L K ) /(1+DELTA(IJ,KL))
c
c                + SUM(KL)           <JB|EB(KL)|IB>
c                                   *    C(IA,JB)
c                                   *    G(L,K)
c
C|    CALL QENTER('LOOP1',0)
      IF (ITSTLP.NE.2) THEN
      DO 2000 IASM = 1, MAXSYM
c
        IBSM = SYMPRO(IHCSM,IASM)
        JBSM = SYMPRO(ICSM,IASM)
c
        DO 1801 IBTP = 1, NOCTPB
c
          NIB   = NSSOB(IBTP,IBSM)
          IBSTR = ISSOB(IBTP,IBSM)
          IBSTP = IBSTR + NIB - 1
c.Number of dets in this block
          IF(PERMSM) THEN
            MNIATP = IBTP
          ELSE
            MNIATP = 1
          END IF
          NIABDT = 0
          DO 7000 IATP = MNIATP,NOCTPA
            IF(IOCOC(IATP,IBTP).EQ.1)
     &      NIABDT = NIABDT + NSSOA(IATP,IASM)*NIB
 7000     CONTINUE
          IF(NIABDT.EQ.0) GOTO 1801
c
          DO 1800 IB = IBSTR,IBSTP
            IBREL = IB - ISSOB(IBTP,IBSM)+1
c
c** Set up list of excitations
c         SUM( IJ .GE. KL )<JB|EB(IJ)EB(KL)|IB>*(JI|LK)/(1+D(IJ,KL) )
c       + SUM(KL)          <JB|EB(KL)|IB>      *G(L,K)
c
            KBSTR = ISTBAB(JBSM)
            CALL ISETVC(NJBFTP,0,NOCTPB)
            CALL ISETVC(INDX1,0,NSTASB(JBSM))
c
c* One-electron terms
c
            NEX = TBSYM(IHSM+1,IB)-TBSYM(IHSM,IB)
            IFRSTX = TBSYM(IHSM,IB)
            CALL SXFSTR(TBIJ(IFRSTX),TBTO(IFRSTX),
     &            NEX,NBSTR,KLTP,ITYPSB,KBSTR,MXPST,F3,
     &            INDX1,INDEX2,NJBFTP,G)
c
c*  Two-electron terms
c
            IF(IDOH2.NE.0)
     &      CALL DXFSTR(IB,IBSM,IHSM,
     &           TBIJ,TBTO,TBSYM,NBEXCI,ICREA,IANNI,NORB,
     &           NBSTR,KLTP,ITYPSB,ISTBAB,MAXSYM,
     &           MXPST,F3,INDX1,INDEX2,NJBFTP,IJKL,INTFRM,
     &           SYMPRO,TUVX,ILTSOB,WIJKL)
            JBSTR = ISTBAB(JBSM)
c
            DO 1750 IATP = MNIATP,NOCTPA
              NIA = NSSOA(IATP,IASM)
              IF(NIA.EQ.0.OR.IOCOC(IATP,IBTP).NE.1 ) GOTO 1750
              IHC = IHCOOS(IBTP,IATP,IASM) + (IBREL-1)*NIA
              IC  = ICOOS (1   ,IATP,IASM)
              CALL ISETVC(INDX4,0,NSTASB(JBSM))
c
              NJBEFF = 1
              NJBXCL= 0
c for given alpha type find the allowed beta types
              DO 1740 JBTP = 1, NOCTPB
                IF(IOCOC(IATP,JBTP) .NE. 1) THEN
                  NJBXCL = NJBXCL + NSSOB(JBTP,JBSM)
                ELSE
                 DO 1730 JB = 1, NJBFTP(JBTP)
                   FACTOR(JB+NJBEFF-1,1) = F3(JB,JBTP)
                   INDX4(JB+NJBEFF-1) = INDEX2(JB,JBTP)-NJBXCL
 1730            CONTINUE
                 NJBEFF = NJBEFF + NJBFTP(JBTP)
                END IF
 1740         CONTINUE
              NJBEFF = NJBEFF - 1
c Perform multiple SAXPY operation
              CALL MSAXPY(HC(IHC),FACTOR,C(IC),TEST,
     &             NIA,NJBEFF,INDX4,NJBEFF)
 1750       CONTINUE
 1800     CONTINUE
 1801   CONTINUE
c
 2000 CONTINUE
      END IF
c
      IF( IPRINT.GE. 5 ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' SIGMA VECTOR AFTER FIRST LOOP '
        WRITE(LUPRI,*)
        CALL WRTMAT(HC(1),1,NHCDET,1,NHCDET,0)
      END IF
c
C|    CALL QEXIT('LOOP1',0)
C|    CALL QENTER('LOOP2',0)
c**************************************************************
c*                                                            *
c**            Term 2 to sigma vector                         *
c*                                                            *
c**************************************************************
c    SIGMA(IA,IB) = SIGMA(IA,IB)  +                      )
c                   SUM(IJ.GE.KL) * (JA|EA(IJ) EA(KL) |IA)
c                                 * C(JA,IB)
c                                 * ( J I | L K )/(1+DELTA(IJ,KL))
c                 +
c                   SUM ( KL )     G(L,K)*(JA|EA(KL)|IA)
c
      IF (ITSTLP.NE.2) THEN
      DO 3001 IASM = 1, MAXSYM
c
c symmetry of strings in sigma(ia,ib) and c(ja,jb)
        IBSM = SYMPRO(IASM,IHCSM)
        JBSM = IBSM
        JASM = SYMPRO(ICSM,IBSM)
c Transpose symmetry block of C
        IBASE = 1
        DO 5607 JBTP = 1, NOCTPB
        DO 5607 JATP = 1, NOCTPA
          IC2OOS(JBTP,JATP) = IBASE
          IF(IOCOC(JATP,JBTP).EQ.1  ) THEN
            NRIN = NSSOA(JATP,JASM)
            NCIN = NSSOB(JBTP,JBSM)
            CALL TRPMAT(C(ICOOS(JBTP,JATP,JASM)),
     &           NRIN,NCIN,C2(IBASE) )
            IBASE = IBASE + NRIN * NCIN
          END IF
 5607   CONTINUE
c
        DO 3000 IATP = 1, NOCTPA
          NIA    = NSSOA(IATP,IASM)
c.Number of dets in this block
          IF(PERMSM) THEN
            MXIBTP = IATP-1
          ELSE
            MXIBTP = NOCTPB
          END IF
          NIABDT = 0
          DO 7020 IBTP = 1, MXIBTP
            IF(IOCOC(IATP,IBTP).EQ.1)
     &      NIABDT = NIABDT + NSSOB(IBTP,IBSM)*NIA
 7020     CONTINUE
          IF(NIABDT.EQ.0) GOTO 3000
c
          IASTR = ISSOA(IATP,IASM)
          IASTP = IASTR + NIA - 1
          DO 2800 IA = IASTR,IASTP
c
c* obtain all excitations out from string IA
c
            CALL ISETVC(NJAFTP,0,NOCTPA)
            CALL ISETVC(INDX1,0,NSTASA(JASM))
c
c* One electron excitations
c
            NEX = TASYM(IHSM+1,IA)-TASYM(IHSM,IA)
            IFRSTX = TASYM(IHSM,IA)
            KASTR = ISTBAA(JASM)
            CALL SXFSTR(TAIJ(IFRSTX),TATO(IFRSTX),
     &           NEX,NASTR,KLTP,ITYPSA,KASTR,MXPST,F3,
     &           INDX1,INDEX2,NJAFTP,G)
c
c*  Two electron excitations
c
            IF(IDOH2.NE.0)
     &      CALL DXFSTR(IA,IASM,IHSM,
     &           TAIJ,TATO,TASYM,NAEXCI,ICREA,IANNI,NORB,
     &           NASTR,KLTP,ITYPSA,ISTBAA,MAXSYM,
     &           MXPST,F3,INDX1,INDEX2,NJAFTP,IJKL,INTFRM,
     &           SYMPRO,TUVX,ILTSOB,WIJKL)
c
            DO 2750 IBTP = 1,MXIBTP
              NIB = NSSOB(IBTP,IBSM)
              IF(NIB.EQ.0.OR.IOCOC(IATP,IBTP).NE.1) GOTO 2750
              IHC = IHCOOS(IBTP,IATP,IASM) + IA-IASTR
              IC  = IC2OOS(IBTP,1)
              CALL ISETVC(INDX4,0,NSTASA(JASM))
c for the given beta type find the allowed alpha types
              NJAEFF = 0
              NJAXCL = 0
              DO 2740 JATP = 1, NOCTPA
                IF(IOCOC(JATP,IBTP) .NE. 1) THEN
                   NJAXCL = NJAXCL + NSSOA(JATP,JASM)
                ELSE
                  DO 2730 JA = 1, NJAFTP(JATP)
                    FACTOR(JA+NJAEFF,1) = F3(JA,JATP)
                    INDX4(JA+NJAEFF) = INDEX2(JA,JATP) - NJAXCL
 2730             CONTINUE
                  NJAEFF = NJAEFF + NJAFTP(JATP)
                END IF
 2740         CONTINUE
c Perform multiple SAXPY operation
              CALL MSAXTY(VEC2,FACTOR,C2(IC),TEST,
     &              NIB,NJAEFF,INDX4,NJAEFF)
c add to total sigma vector
              DO 2300 IB = 1, NIB
                IIIHC = IHC + ( IB - 1 ) * NIA
                HC(IIIHC) = HC(IIIHC) + SIGNLL*VEC2(IB)
 2300         CONTINUE
 2750       CONTINUE
c
 2800     CONTINUE
 3000   CONTINUE
 3001 CONTINUE
      IF(PERMSM) THEN
c.Only blocks of sigma vector at or below the diagonal were
c obtained.Obtain blocks above the diagonal
c.Only aplha-alpha part of type diagonal blocks were obtained.
c.Obtain beta-component by transposing
        DO 4856 IATP = 1, NOCTPA
          DO 4700 IBTP = 1, IATP
            IF(IOCOC(IATP,IBTP).EQ.1) THEN
              DO 5000 IASM = 1, MAXSYM
                IBSM = SYMPRO(IASM,IHCSM)
                IBASEI = IHCOOS(IBTP,IATP,IASM)
                IBASEO = IHCOOS(IATP,IBTP,IBSM)
                NA = NSSOA(IATP,IASM)
                NB = NSSOB(IBTP,IBSM)
                IF(IATP.NE.IBTP) THEN
                  CALL TRPMAT(HC(IBASEI),NA,NB,HC(IBASEO))
                  CALL SCALVE(HC(IBASEO),PSIGN,NA*NB)
                ELSE IF (IATP.EQ.IBTP.AND.IASM.LE.IBSM) THEN
                  CALL TRPMAT(HC(IBASEO),NB,NA,C2)
                  CALL VECSUM(HC(IBASEI),HC(IBASEI),C2,
     &                      1.0D0,PSIGN,NA*NB)
                  IF(IASM.NE.IBSM) THEN
                    CALL TRPMAT(HC(IBASEI),NA,NB,HC(IBASEO))
                    CALL SCALVE(HC(IBASEO),PSIGN,NA*NB)
                  END IF
                END IF
 5000         CONTINUE
            END IF
 4700     CONTINUE
 4856   CONTINUE
      END IF
      END IF
c
      IF( IPRINT.GE. 5 ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' SIGMA VECTOR AFTER SECONDI LOOP '
        WRITE(LUPRI,*)
        CALL WRTMAT(HC(1),1,NHCDET,1,NHCDET,0)
      END IF
C|    CALL QEXIT('LOOP2',0)
c**************************************************************
c*                                                            *
c**            Term 3 to sigma vector                         *
c*                                                            *
c**************************************************************
      IF(IDOH2.NE.0.AND.ITSTLP.NE.1) THEN
C|    CALL QENTER('LOOP3',0)
      NOO = NORB * NORB
c
      DO 4000 IASM = 1, MAXSYM
        IBSM = SYMPRO(IASM,IHCSM)
         IF(PERMSM.AND.IASM.LT.IBSM.AND.NOCTPA.EQ.1) GOTO 4000
c
        IRSTRT = ISTBAA(IASM)
        IF( NSTASA(IASM) .EQ. 0 ) GOTO 3999
c
c Number of determinants in this symmetry block
        NIAIB = 0
        DO 7645 IATP = 1,NOCTPA
        IF(PERMSM) THEN
          IF(IASM.GE.IBSM) THEN
            IBTPMX = IATP
          ELSE
C|          IBTPMX = IATP - 1
            IBTPMX = IATP
          END IF
        ELSE
          IBTPMX = NOCTPB
        END IF
        DO 7645 IBTP = 1,IBTPMX
          IF(IOCOC(IATP,IBTP) .EQ. 1 )
     &    NIAIB = NIAIB + NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
 7645   CONTINUE
        IF(NIAIB .EQ. 0 ) GOTO 3999
c
        DO 3800 JASM = 1, MAXSYM
          JBSM = SYMPRO(JASM,ICSM)
c
          IF(NSTASA(JASM) .EQ. 0 ) GOTO 3800
          KLSM = SYMPRO(IASM,JASM)
          LKSM = KLSM
          IJSM = SYMPRO(IBSM,JBSM)
c Number of dets in this symmetry block
          NJAJB = 0
          DO 7646 JATP = 1,NOCTPA
          DO 7646 JBTP = 1,NOCTPB
            IF(IOCOC(JATP,JBTP) .EQ. 1 )
     &      NJAJB = NJAJB + NSSOA(JATP,JASM)*NSSOB(JBTP,JBSM)
 7646     CONTINUE
          IF(NJAJB .EQ. 0 ) GOTO 3800
c
          CALL ISETVC(INDX1,0,NSTASB(JBSM) )
          DO 3710 K = 1, NORB
          IFRST = 1
          DO 3711 LTYPE = 1,3
c
c.Number of excitations, that has K as an creation operator
c.and is of symmetry KLSM
            LEX = 0
            KLEX0 = (K-1)*NORB
            DO 3709 KLEX = KLEX0+1, KLEX0+NORB
            IF(DISSYM(KLEX).EQ.KLSM.AND.ITPFOB(IANNI(KLEX)).EQ.LTYPE)
     &      THEN
              LEX = LEX + 1
              LORB(LEX) = IANNI(KLEX)
            END IF
 3709       CONTINUE
            IF(LEX.EQ.0) GOTO 3711
c. Obtain all strings of symmetry JSYM that contains orbital K
            IOFFF = ISTBAA(JASM)
            IF(IFRST.EQ.1) THEN
C|            CALL QENTER('OINST',0)
C             CALL OINSTR(ISTINA(1,IOFFF),NSTASA(JASM),NAEL,K,
C    &             LSTRIN,IOST,NOCTPA,ITYPSA(IOFFF),IFTP2)
              KTYP = ITPFOB(K)
              CALL OINST3(ISTINA(1,IOFFF),NSSOA(1,JASM),NAEL,
     &             K,KTYP,LSTRIN,IOST,NOCTPA,IFTP2,1)
              IF (LSTRIN .GT. 2*MXPST) THEN
                 WRITE(LUPRI,'(2(/A,I10))')
     &  ' >> MXPST in CISIG9 insufficient, current value  <<',MXPST,
     &  ' >> Raise (2*MXPST) in CISIG9 to at least        <<',LSTRIN
                 CALL QUIT('CISIG9 after OINST3: MXPST too small')
              END IF
              IFRST = 0
C|            CALL QEXIT('OINST',0)
            END IF
            IF (LSTRIN.EQ.0) GOTO 3710
c
c*  Gather C2(I,JB) = C(IOST(I),JB)
c
            ISIGN = 0
C|          CALL QENTER('RASCG',0)
            CALL RASCG2(C(ICOOS(1,1,JASM)),C2,NOCTPA,NSSOA(1,JASM),
     &           IFTP2(1),NOCTPB,NSSOB(1,JBSM),IC2OOS,IOST,FAC2,IOCOC,
     &           JBTPMX,ISIGN)
            JBTPMO = JBTPMX
C|          CALL QEXIT('RASCG',0)
c.Loop over batches of L orbitals
            NBATCH = LEX/MXSMOB
            IF(NBATCH*MXSMOB.LT.LEX) NBATCH = NBATCH + 1
            DO 3705 IBATCH = 1, NBATCH
               IF(IBATCH.EQ.1) THEN
                  LSTRT = 1
                ELSE
                  LSTRT = LSTRT + MXSMOB
                END IF
                LSTOP = MIN(LEX,LSTRT+MXSMOB-1)
                NL = LSTOP-LSTRT+1
c Set up |L(I)> = EA(KL)|R(I)> * SIGN(I) for L in batch
              DO 9010 LL = LSTRT,LSTOP
                LABS = LORB(LL)
                LREL = LL - LSTRT + 1
                KL = (K-1)*NORB + LABS
                LK = KLTP(KL)
                KL2 = 0
c
c
c The list is ordered according to L, so first all L strings of type
c 1 comes , then all L strings of type 2 etc .
c
c The L strings are numbered with base for given symmetry and
c type as base, while R is numbered with beginning of symmetry
c as base
c <L|E(KL)|R>  is obtained from lists <R|E(LK)|L>
c
                ISEL = 1
                ILSEL = 1
C|              CALL QENTER('LRFEX',0)
                CALL LRFEX5(LK,KL2,LKSM,R(1,LREL),L(1,LREL),
     &               FAC2(1,LREL),INDX7(1,LREL),NFTPA(1,LREL),ITYPSA,
     &               TAIJ,TATO,
     &               NAEXCI,NKLSTA(LREL),NSSOA(1,JASM),NOCTPA,NASTR,
     &               TASYM,MAXSYM,IRSTRT,IFTPA(1,LREL),
     &               ISEL,IFTP2,IOST,ILSEL,1,ISTBAA(JASM) )
C|              CALL QEXIT('LRFEX',0)
c. Integrals (**|KL)
C|              CALL QENTER('GETIN',0)
                IF (ISPIN1+ISPIN2.GE.1) DISTYP = DISTYP + 2
                CALL GETIN2(IJKL(1,LREL),LABS,IANNI(LK),TUVX,
     &               ILTSOB,WIJKL)
C|              CALL QEXIT('GETIN',0)
                IF (ISPIN1+ISPIN2.GE.1) DISTYP = DISTYP - 2
c
 9010         CONTINUE
c*  Loop over beta strings and calculate terms to SIGMA(R(I),IB)
c
c
c The sigma term can now be written
c sigma(R(I),IB) = <JB|EB(ij)|IB> * C2(L(I),JB)
c
              IF (ISPIN1+ISPIN2.EQ.1) THEN
                SIGN12 = -1.0D0
              ELSE
                SIGN12 = 1.0D0
              END IF
              JBSTR = ISTBAB(JBSM)
              DO 3501 IBTP = 1, NOCTPB
C|            IF(IBTP.EQ.1) THEN
C|              CALL QENTER('SEC1',0)
C|            ELSE IF(IBTP.EQ.2) THEN
C|              CALL QENTER('SEC2',0)
C|            ELSE IF(IBTP.EQ.3) THEN
C|              CALL QENTER('SEC3',0)
C|            END IF
c
                IF (PERMSM) THEN
                 IF(IASM.GE.IBSM) THEN
                   IATPMN = IBTP
                 ELSE
C|                 IATPMN = IBTP + 1
                   IATPMN = IBTP
                 END IF
                ELSE
                 IATPMN = 1
                END  IF
                JBTPMX = JBTPMO
c
                IATPMX = 0
                DO 7040 IATP = IATPMN,NOCTPA
                  IF(IOCOC(IATP,IBTP).EQ.1) IATPMX = IATP
 7040           CONTINUE
                IF(IATPMN.EQ.IATPMX.AND.IATPMN.EQ.IBTP.AND.PERMSM)
     &          THEN
                  IDIASM = 1
                  JBTPM2 = 0
                  DO 7050 JATP = 1,NOCTPA
                  DO 7050 JBTP = 1, JATP
                    IF(IOCOC(JATP,JBTP).EQ.1) JBTPM2= MAX(JBTPM2,JBTP)
 7050             CONTINUE
                  JBTPMX = MIN(JBTPMO,JBTPM2)
                ELSE
                  IDIASM = 0
                END IF
                MXNIA = 0
                DO 5011 LL = 1,NL
                  NIATOT = 0
                  DO 5010 IATP = IATPMN, NOCTPA
                    IF(IOCOC(IATP,IBTP) .EQ. 1 )
     &              NIATOT = NIATOT + INDX7(IATP,LL)
 5010             CONTINUE
                  MXNIA = MAX(MXNIA,NIATOT)
 5011           CONTINUE
                IF( MXNIA.EQ.0) GOTO 3502
c
                NIB = NSSOB(IBTP,IBSM)
                IBSTR = ISSOB(IBTP,IBSM)
                IBSTP = IBSTR + NIB - 1
c
                DO 3500 IB = IBSTR,IBSTP
                  IBREL = IB - ISSOB(IBTP,IBSM) + 1
c* all single excitations out from string IB
                  CALL ISETVC(NJBFTP,0,NOCTPB)
                  NJBEFF = 0
                  DO 3400 IJEX = TBSYM(IJSM,IB),TBSYM(IJSM+1,IB)-1
                    IJ = TBIJ(IJEX)
                    JB = TBTO(IJEX)
                    IF(JB. GT. NBSTR) THEN
                      JB = JB - NBSTR
                      SIGN = -SIGN12
                    ELSE
                      SIGN = SIGN12
                    END IF
                    JBREL = JB - JBSTR + 1
                    JBTP = ITYPSB (JB)
                    IF(JBTP .GT. JBTPMX) GOTO 3401
c
                    IF(INDX1(JBREL) .NE. 0 ) THEN
                      DO 9020 LL = 1,NL
                        F3(INDX1(JBREL),LL) = F3(INDX1(JBREL),LL)+
     &                                        SIGN*IJKL(IJ,LL)
 9020                 CONTINUE
                    ELSE
                      NJBEFF=NJBEFF+1
                      INDX1(JBREL) = NJBEFF
                      INDX6(NJBEFF) = JBREL
                      NJBFTP(JBTP) = NJBFTP(JBTP)+1
                      DO 9030 LL = 1, NL
                        F3(NJBEFF,LL) = SIGN*IJKL(IJ,LL)
 9030                 CONTINUE
                    END IF
c
 3400             CONTINUE
 3401             CONTINUE
                  NJBTOT = NJBEFF
c
                  DO 4800 IATP = IATPMN, NOCTPA
                    IF(IOCOC(IATP,IBTP) .NE. 1 ) GOTO 4800
                    NIA = NSSOA(IATP,IASM)
                    IF(NIA .EQ. 0 ) GOTO 4800
c type JATP so |IATP> =  E(KL) |JATP>
                    IADD = ISTBAA(IASM)-1
                    JATP = - 1
c To be modified
                    DO 4791 LL = 1, NL
                    DO 4790 KATP = 1, NOCTPA
                      IF(NFTPA(KATP,LL).EQ.0) GOTO 4790
                      IF(ITYPSA(R(IFTPA(KATP,LL),LL)+IADD).EQ.IATP)
     &                JATP= KATP
 4790               CONTINUE
 4791               CONTINUE
                    IF(JATP .EQ. -1) GOTO 4800
c
                    NIAGMX = 0
                    DO 4792 LL = 1, NL
                      NIAGMX = MAX(NIAGMX,NFTPA(JATP,LL) )
 4792               CONTINUE
                    IF(NIAGMX.EQ.0) GOTO 4800
                    IHC = IHCOOS(IBTP,IATP,IASM)
                    IASTR = ISSOA(IATP,IASM)
                    IC =IC2OOS(1,JATP)
c
                    NJBEFF = 1
                    NJBT = 1
                    NJBXCL = 0
                    DO 4785 JBTP = 1, JBTPMX
C?                    WRITE(LUPRI,*)' JATP,JBTP,IDIASM',JATP,JBTP,IDIASM
                      FCT = 1.0D0
                      IF(IDIASM.EQ.1.AND.JATP.LT.JBTP)GOTO 4785
                      IF(IDIASM.EQ.1.AND.JATP.NE.JBTP) FCT = 2.0D0
                      NJB = NJBFTP(JBTP)
                      IF(IOCOC(JATP,JBTP) .NE. 1 ) THEN
                        NJBXCL = NJBXCL + NSSOB(JBTP,JBSM)
                      ELSE
                        DO 9040 LL = 1, NL
                          CALL COPVEC(F3(NJBT,LL),FACTOR(NJBEFF,LL),NJB)
                          IF(IDIASM.EQ.1.AND.FCT.EQ.2.0D0)
     &                    CALL SCALVE(FACTOR(NJBEFF,LL),FCT,NJB)
 9040                   CONTINUE
                        CALL ICOPVE(INDX6(NJBT),INDX4(NJBEFF),NJB)
                        IF(NJBXCL .NE. 0 ) THEN
                          DO 4783 JB = NJBEFF, NJBEFF+NJB-1
                            INDX4 (JB) = INDX4(JB)-NJBXCL
 4783                     CONTINUE
                        END IF
                        NJBEFF = NJBEFF + NJB
                      END IF
                      NJBT = NJBT + NJB
 4785               CONTINUE
                    NJBEFF = NJBEFF - 1
c
                    IF( NJBEFF.EQ.0) GOTO 4800
                    NCGROW = IFTP2(JATP)
                    IIHC = IHC + ( IB-IBSTR)*NIA-ISSOA(IATP,IASM)
     &                   +                       ISSOA(1,IASM)-1
c
                    DO 9050 LL = 1, NL
                      NIAG =   NFTPA(JATP,LL)
                      IAGSTR = IFTPA(JATP,LL)
                      IAGSTP = IFTPA(JATP,LL) + NIAG - 1
                      CALL MSAXTY(VEC2(1),FACTOR(1,LL),
     &                     C2(IC),TEST,NCGROW,NJBEFF,INDX4,NJBEFF)
                      DO 3430 I = IAGSTR,IAGSTP
                        HC(IIHC+R(I,LL)) = HCDUM(IIHC+R(I,LL))
     &                              + FAC2(I,LL)*VEC2(L(I,LL))
 3430                 CONTINUE
 9050              CONTINUE
 4800           CONTINUE
c
                  DO 5035 JB = 1, NJBTOT
                    INDX1(INDX6(JB)) = 0
 5035             CONTINUE
c
 3500         CONTINUE
 3502         CONTINUE
C|            IF(IBTP.EQ.1) THEN
C|              CALL QEXIT('SEC1',0)
C|            ELSE IF(IBTP.EQ.2) THEN
C|              CALL QEXIT('SEC2',0)
C|            ELSE IF(IBTP.EQ.3) THEN
C|              CALL QEXIT('SEC3',0)
C|            END IF
 3501      CONTINUE
 3705  CONTINUE
 3711  CONTINUE
 3710  CONTINUE
 3800  CONTINUE
 3999 CONTINUE
c
c
 4000 CONTINUE
c
      IF( IPRINT.GE.10) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Sigma vector after 4000  '
        CALL WRTMAT(HC,1,NHCDET,1,NHCDET,0)
      END IF
      IF(PERMSM) THEN
        DO 6856 IATP = 1, NOCTPA
          DO 6700 IBTP = 1, IATP
            IF(IOCOC(IATP,IBTP).EQ.1) THEN
              DO 6001 IASM = 1, MAXSYM
                IBSM = SYMPRO(IASM,IHCSM)
                IBASEI = IHCOOS(IBTP,IATP,IASM)
                IBASEO = IHCOOS(IATP,IBTP,IBSM)
                NA = NSSOA(IATP,IASM)
                NB = NSSOB(IBTP,IBSM)
                IF(IATP.NE.IBTP) THEN
                  CALL TRPMAT(HC(IBASEI),NA,NB,HC(IBASEO))
                  CALL SCALVE(HC(IBASEO),PSIGN,NA*NB)
                ELSE IF (IATP.EQ.IBTP.AND.IASM.GE.IBSM) THEN
                  IF(NOCTPA.EQ.1.AND.IASM.NE.IBSM) THEN
                    CALL TRPMAT(HC(IBASEI),NA,NB,HC(IBASEO))
                    CALL SCALVE(HC(IBASEO),PSIGN,NA*NB)
                  ELSE IF(NOCTPA.GE.2) THEN
                    CALL TRPMAT(HC(IBASEO),NB,NA,C2)
                    CALL VECSUM(HC(IBASEI),HC(IBASEI),C2,
     &                          0.5D0,0.5D0*PSIGN,NA*NB)
                    IF(IASM.NE.IBSM) THEN
                      CALL TRPMAT(HC(IBASEI),NA,NB,HC(IBASEO))
                      CALL SCALVE(HC(IBASEO),PSIGN,NA*NB)
                    END IF
                  END IF
                END IF
 6001         CONTINUE
            END IF
 6700     CONTINUE
 6856   CONTINUE
      END IF
C|    CALL QEXIT('LOOP3',0)
      END IF
c ( end of IF ( IDOH2 .GT.0 ) THEN )
c
      IF(ISPIN1+ISPIN2.EQ.1) THEN
        CALL SCALVE(HC,-1.0D0,NHCDET)
      END IF
c
      IF( IPRINT.GE.  5) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Final sigma vector '
        CALL WRTMAT(HC,1,NHCDET,1,NHCDET,0)
      END IF
c
      RETURN
      END
C  /* Deck lrfex */
#if defined (VAR_SECSEC)
      SUBROUTINE LRFEX(IJ,JI,IJSM,L,R,SIGN,NLFTP,NRFTP,ITYPST,
     &           IIJ,ITO,NEXCI,NTERM,NSTFTP,NTYPE,NSTR,IEXPSM,
     &           MAXSYM,ILSTRT,IBRFTP)
c
c     PURPOSE: For two symmetry blocks of strings set up
c              - <L| E |R> for excitations E, that can be E(ij) or E(ji)
c                in the form of three arrays L,R,ISIGN
c              - |L(I)> = E |R(I)> * SIGN(I)
c                The list is ordered according to R, so first
c                all R strings of type 1 comes , then all R strings
c                of type 2 etc .
c
c     NOTE:    The R strings are numbered with base for given
c              symmetry and type as base, while L is numbered
c              with beginning of symmetry as base
c
c
c     Input
c     ======
c     IJ,JI : numbers of two single excitations ( can be zero )
c     IJSM  : symmetry of excitation IJ (JI)
c     IIJ   : single excitations for each string
c     ITO   : excited strings for each
c     NEXCI : (Max ) number of excitations per string
c     NTYPE : Number of types ( RAS types )
c     NSTR  : Total number of strings
c     IEXPSM : Base for excitations of given symmetry from given string
c     MAXSYM  : Number of symmetries of single excitations
c     ILSTRT : string number of first string in L-block
c
c     Output
c     ======
c     L : L array above
c     R : R array above
c     SIGN : Sign array above
c     NLFTP(*) : number of strings in L per type of L
c     NRFTP(*) : number of strings in R per type of R
c     NTERM :  number of strings in L and R
c     IBRFTP : First string in R of given type
c
c     Jeppe Olsen , Sept 1989
c
#include <implicit.h>
c
      INTEGER R
      DIMENSION ITO(NEXCI,*),IIJ(NEXCI,*),NSTFTP(*),ITYPST(*),
     &          IEXPSM(MAXSYM+1,*)
      DIMENSION L(*),R(*),SIGN(*),NLFTP(*),NRFTP(*),IBRFTP(*)
c
      NTERM = 0
      DO 400 IRTP = 1, NTYPE
        IBRFTP(IRTP) = NTERM + 1
c
        IF( IRTP .EQ. 1 ) THEN
          IRSTR = 1
        ELSE
          IRSTR = IRSTR + NSTFTP(IRTP-1)
        END IF
        IRSTP = IRSTR + NSTFTP(IRTP)-1
c
        NTRFTP = 0
        NTERMP = NTERM
        DO 300 IR = IRSTR,IRSTP
c Is E(IJ)|IR> or E(JI)|IR> nonvanishing
          DO 200 IEX = IEXPSM(IJSM,IR),IEXPSM(IJSM+1,IR)-1
            IF(IIJ(IEX,IR).EQ.IJ.OR.IIJ(IEX,IR).EQ.JI ) THEN
              NTRFTP = NTRFTP + 1
              NTERM  = NTERM + 1
              R(NTERM) = IR- IRSTR +1
              IF(ITO(IEX,IR) .GT. NSTR) THEN
                 L(NTERM) = ITO(IEX,IR)-NSTR-ILSTRT+1
                 SIGN(NTERM) = -1.0D0
              ELSE
                 L(NTERM) = ITO(IEX,IR)-ILSTRT + 1
                 SIGN(NTERM) =  1.0D0
              END IF
              GOTO 201
            END IF
 200      CONTINUE
 201      CONTINUE
 300    CONTINUE
c
        NRFTP(IRTP) = NTRFTP
 400  CONTINUE
c Number of L strings per type
      CALL ISETVC(NLFTP,0,NTYPE)
      DO 500 ILSTR = 1,NTERM
         NLFTP(ITYPST(L(ILSTR)+ILSTRT-1)) =
     &   NLFTP(ITYPST(L(ILSTR)+ILSTRT-1)) + 1
 500  CONTINUE
c
      RETURN
      END
#endif
C  /* Deck dxfstr */
      SUBROUTINE DXFSTR(ISTR,ISTRSM,IGSM,
     &           IEX,ITO,IEXPSM,NEXCI,ICREA,IANNI,NORB,
     &           NSTR,KLTP,ITYPST,ISTASM,MAXSYM,
     &           IDIM,FACTOR,INDEXG,INDEXS,NEXFTP,IJKL,INTFRM,SYMPRO,
     &           TUVX,ILTSOB,SCR)
c
c     purpose: Multiply a two electron operator G times string ISTR
c              G |ISTR> = sum(i,j,k,l,(ij).ge.(kl) )
c                         E(ij)E(kl) |ISTR> (ji|lk)/(1+delta(ij,kl) )
c     Input
c     =====
c     ISTR : String to be excited from
c     ISTRSM : Symmetry  of ISTR
c     IGSM : Symmetry of two-electron operator
c     IEX : List of excitations from all strings
c     ITO : List of strings arising from excitations in IEX
c     IEXPSM: Pointer to first excitation of a given symmetry
c     NEXCI: Max Number of excitations per string
c     ICREA : Creation operator of a given excitation
c     IANNI : Annihilation operator of a given excitation
c     NORB  : Total number of orbitals
c     NEX : Number of excitations from string
c     NSTR : Total number of strings
c     KLTP : KLTP(IJ) = JI
c     ITYPST : Array giving type of string
c     ISTASM : First string of given symmetry
c     MAXSYM : Number of symmetries
c     IDIM : Row dimension of FACTOR and INDEXS and INDEXG
c     IJKL : Real Double precision scratch space ,size NORB X NORB
c     INTFRM : Defines input mode for integrals : 1=> fetch with getint
c                                                 2=> fetch with geth2a
c                                                 3=> fetch direct from
c                                                     integral list
c     SYMPRO : D2H multiplication table
c     TUVX : List of two-electron integrals
c     ILTSOB : LUNAR to SIRIUS order of orbitals
c     SCR    : SCR array , length : NORB ** 2
c
c     Output
c     ======
c     INDEXS : INDEXS(I) is string number of excited string I
c     INDEXG : INDEXG(I) is place  in factor for string number I
c     ( INDEXS is scatter index and INDEXG is gather index
c       for going between full string indexing and local indexing
c       for exciting from specific string )
c     FACTOR : FACTOR(I,ITYP) Is constant multiplying excited string I
c              of type ITYP
c     NEXFTP : Number of excitations per string
c     G times string =
c     sum(type,i=1,nterms) factor(i,type) *string(index(i))
c
c     Remarks
c     =======
c     1 : output arrays are not zeroed before use
c     2 : No symmetry-checking is performed.This
c         should be be performed by proper arguments.
c     3 : Full string number is defined relative to first
c         string in symmetry block.
c
c     Jeppe Olsen , sept 1989
c
c Direct access of integrals version
c
c =================
c SIRIUS version ||
c =================
c
#include <implicit.h>
#include <priunit.h>
      INTEGER SYMPRO(8,8)
      REAL * 8 IJKL(*)
      DIMENSION IEX(*),ITO(*),IEXPSM(MAXSYM+1,NSTR),
     &          KLTP(*),ITYPST(*),ISTASM(*),ICREA(*),IANNI(*)
      DIMENSION INDEXS(IDIM,*),INDEXG(*),
     *          FACTOR(IDIM,*),NEXFTP(*)
      DIMENSION TUVX(*),ILTSOB(*),SCR(*)
CHJ: check INTFRM
      IF (INTFRM .NE. 1 .AND. INTFRM .NE. 3) THEN
         WRITE (LUPRI,*) 'Illegal INTFRM in DXFSTR, INTFRM =',INTFRM
         WRITE (LUPRI,*) 'Allowed values are 1 or 3'
         CALL QUIT('Illegal INTFRM in DXFSTR')
      END IF
c. Symmetry of excited block
      JSM = SYMPRO(ISTRSM,IGSM)
      JSTRT= ISTASM(JSM)
c Loop over symmetries of E(kl)
      DO 500 KLSM = 1, MAXSYM
c Symmetry of E(ij)
        IJSM = SYMPRO(KLSM,IGSM)
c Loop over excitations from ISTR of symmetry KLSM
        DO  400 KLEX  = IEXPSM(KLSM,ISTR),IEXPSM(KLSM+1,ISTR) - 1
          KL = IEX(KLEX)
          LK = KLTP(KL)
c
          KSTR = ITO(KLEX)
          SIGN = 1.0D0
          IF(KSTR .GT. NSTR) THEN
            KSTR = KSTR - NSTR
            SIGN = -1.0D0
          END IF
c Obtain integrals (J I | K L )
          IF(INTFRM .EQ. 1 ) THEN
            CALL GETIN2(IJKL,ICREA(LK),IANNI(LK),TUVX,ILTSOB,SCR)
          ELSE IF(INTFRM.EQ.3) THEN
            K = ILTSOB(ICREA(LK))
            L = ILTSOB(IANNI(LK))
            KK = MAX(L,K)
            LL = MIN(L,K)
            KLP = KK*(KK-1)/2 + LL
            KLOFF = (KLP-1)*NORB*(NORB+1)/2
          END IF
c Excitations out from string KSTR of symmetry IJSM
          DO 300 IJEX = IEXPSM(IJSM,KSTR),IEXPSM(IJSM+1,KSTR)-1
            IJ = IEX(IJEX)
            IF( KL .GT. IJ ) GOTO 300
C
            JSTR = ITO(IJEX)
            SIGN1 = SIGN
            IF( JSTR .GT. NSTR ) THEN
              JSTR = JSTR - NSTR
              SIGN1 = - SIGN
            END IF
            JSTRRL = JSTR - JSTRT + 1
            JTP = ITYPST(JSTR)
            IF(INTFRM.EQ.1) THEN
              XIJKL = IJKL(IJ)
            ELSE IF(INTFRM.EQ.3) THEN
              I = ILTSOB(ICREA(IJ))
              J = ILTSOB(IANNI(IJ))
              II = MAX(I,J)
              JJ = MIN(I,J)
              IJP = II*(II-1)/2 + JJ
              XIJKL = TUVX(IJP + KLOFF)
            END IF
            IF( IJ .EQ .KL ) XIJKL = 0.5D0 * XIJKL
c
            IF(INDEXG(JSTRRL) .NE. 0 ) THEN
              FACTOR(INDEXG(JSTRRL),JTP) =
     &        FACTOR(INDEXG(JSTRRL),JTP) + SIGN1 * XIJKL
            ELSE
              NEXFTP(JTP) = NEXFTP(JTP) + 1
              INDEXG(JSTRRL) =NEXFTP(JTP)
              INDEXS(NEXFTP(JTP),JTP) = JSTRRL
              FACTOR(INDEXG(JSTRRL),JTP) = SIGN1*XIJKL
            END IF
  300     CONTINUE
  400   CONTINUE
  500 CONTINUE
c
C1000 CONTINUE
c
      RETURN
      END
C  /* Deck sxfstr */
      SUBROUTINE SXFSTR(IEX,ITO,NEX,NSTR,KLTP,ITYPST,IFEXST,
     &           IDIM,FACTOR,INDEXG,INDEXS,NEXFTP,G)
c
c     PURPOSE: An operator G is given as
c              G = sum(i) H(T)(IEX(I)) E(IX(I))
c              Calculate G times a string
c
c     Input
c     =====
c     IEX : List of excitations from this string
c     ITO : List of excited strings
c     NEX : Number of excitations from string
c     NSTR : Total number of strings
c     KLTP : KLTP(IJ) = JI
c     ITYPST : Array giving type of string
c     IFEXST : Number of first string in excited symmetry block
c     IDIM : Row dimension of FACTOR and INDEXS and INDEXG
c
c     Output
c     ======
c     INDEXS : INDEXS(I) is string number of excited string I
c     INDEXG : INDEXG(I) is place  in factor for string number I
c     ( INDEXS is scatter index and INDEXG is gather index
c       for going between full string indexing and local indexing
c       for exciting from specific string )
c     FACTOR : FACTOR(I,ITYP) Is constant multiplying excited string I
c              of type ITYP
c     NEXFTP : Number of excitations per string
c     G times string =
c     sum(type,i=1,nterms) factor(i,type) *string(index(i))
c
c     Remarks
c     =======
c     1 : output arrays are not zeroed before use
c     2 : No symmetry-checking is performed.This
c         should be be performed by proper arguments.
c     3 : Full string number is defined relative to first
c         string in symmetry block.
c
c     Jeppe Olsen , sept '89
c
#include <implicit.h>
      DIMENSION IEX(*),ITO(*),KLTP(*),ITYPST(*)
      DIMENSION INDEXS(IDIM,*),INDEXG(*)
      DIMENSION FACTOR(IDIM,*),NEXFTP(*),G(*)
c
      DO 100 KLEX = 1,NEX
c
        LK = KLTP(IEX(KLEX))
c
        KSTR = ITO(KLEX)
        SIGN = 1.0D0
        IF(KSTR .GT. NSTR ) THEN
         SIGN = -1.0D0
         KSTR = KSTR - NSTR
        END IF
        KSTRRL = KSTR - IFEXST + 1
        KTP = ITYPST(KSTR)
c
        IF(INDEXG(KSTRRL) .NE. 0 ) THEN
          FACTOR(INDEXG(KSTRRL),KTP)= FACTOR(INDEXG(KSTRRL),KTP) +
     &    SIGN*G(LK)
        ELSE
          NEXFTP(KTP) = NEXFTP(KTP) + 1
          INDEXG(KSTRRL) = NEXFTP(KTP)
          INDEXS(NEXFTP(KTP),KTP ) = KSTRRL
          FACTOR(INDEXG(KSTRRL),KTP) = SIGN * G(LK)
        END IF
  100 CONTINUE
c
      RETURN
      END
C  /* Deck rascg2 */
      SUBROUTINE RASCG2(C,CG,NTYPR,NFRFTP,NGRFTP,NTYPC,NFCFTP,
     &                 IBOUT,IG,SIGN,IOCOC,MXGCTP,ISIGN)

c
c Gather a set of RAS blocks for a given symmetry
c
c ISIGN .NE. 0 CG(I,J) = C(IG(I),J)*SIGN(I)
c ISIGN .EQ. 0 CG(I,J) = C(IG(I),J)
c
c Input
c======
c C : fulls symmetry block of RAS vector
c NTYPR : Number of types of rows
c NFRFTP : Number of rows per type for full matrix
c NGRFTP : Number of rows per type for gathered matrix
c NTYPC : Number of types of columns
c NFCFTP : Number of columns per symetry for full matrix
c IG    : Gathering vector,in this version with start of
c         symmetry block as base
c SIGN  : sign shift array
c IOCOC : array telling whether a given block is included or not
c IBRFTP : First Row element of a given type
c
c Output
c ======
c CG : gathered matrix as above
c IBOUT : adress of first element in CG belonging to a
c         subblock defined by given row and column type
c MXGCTP : last column type with a nonvanishing number of rows
c
c Jeppe Olsen , October 1989
c
#include <implicit.h>
#include <mxpdim.h>
      DIMENSION C(*),CG(*),NFRFTP(*),NGRFTP(*),NFCFTP(*),
     &          IBOUT(MXPTP,MXPTP),IG(*),SIGN(*),
     &          IOCOC(NTYPR,NTYPC)
c
      MXGCTP = 0
      IBSGC = 1
      IBSFC = 1
      DO  200 IRTP = 1, NTYPR
c
       NFR = NFRFTP(IRTP)
       NGR = NGRFTP(IRTP)
       IF(IRTP .EQ. 1 ) THEN
          IBFR = 1
          IBGR = 1
       ELSE
          IBFR = IBFR + NFRFTP(IRTP-1)
          IBGR = IBGR + NGRFTP(IRTP-1)
       END IF
       DO  150 ICTP = 1, NTYPC
         IBOUT(ICTP,IRTP) = IBSGC
         IF ( IOCOC(IRTP,ICTP) .EQ. 1 ) THEN
c
            NFC = NFCFTP(ICTP)
            IF(NFC*NGR .NE. 0 ) MXGCTP = MAX(MXGCTP,ICTP)
            IF(NGR.EQ.0)GOTO  140
            IICI = IBSFC-1-NFR
            IICO = IBSGC-IBGR-NGR
c
            IF(ISIGN.NE. 0 ) THEN
              DO 20 IC = 1, NFC
                 IICO = IICO + NGR
                 IICI = IICI + NFR
                 IICIP = IICI - IBFR + 1
                 DO 10 IR = IBGR, IBGR+NGR-1
                      CG(IICO+IR ) = C(IICIP + IG(IR))*SIGN(IR)
   10            CONTINUE
   20         CONTINUE
            ELSE
              DO 21 IC = 1, NFC
                 IICO = IICO + NGR
                 IICI = IICI + NFR
                 IICIP = IICI - IBFR + 1
                 DO 11 IR = IBGR, IBGR+NGR-1
                      CG(IICO+IR ) = C(IICIP + IG(IR))
   11            CONTINUE
   21         CONTINUE
            END IF
c
  140       CONTINUE
            IBSFC = IBSFC + NFR*NFC
            IBSGC = IBSGC + NGR*NFC
          END IF
c
  150   CONTINUE
  200 CONTINUE
c
C     NTEST = 0
C     IF(NTEST.NE. 0 ) THEN
C       WRITE(LUPRI,*) ' Output matrix  form RASCG2 '
C            PRCLM4(A,NROWBL,NCOLBL,LROW,LCOL,IRC)
C       CALL PRCLM4(CG,NTYPR,MXGCTP,NGRFTP,NFCFTP,IOCOC)
C     END IF
c
      RETURN
      END
C  /* Deck lrfex5 */
      SUBROUTINE LRFEX5(IJ,JI,IJSM,L,R,SIGN,NLFTP,NRFTP,ITYPST,
     &           IIJ,ITO,NEXCI,NTERM,NSTFTP,NTYPE,NSTR,IEXPSM,
     &           NSMSX,ILSTRT,IBRFTP ,
     &           ISEL,ISELTP,ISELST,ILSEL,IAB,IRSTRT)
c
c For two symmetry blocks of strings set up
c
c      <L| E |R> for excitations E, that can be E(ij)
c in the form of three arrays L,R,ISIGN
c
c |L(I)> = E |R(I)> * SIGN(I)
c
c The list is ordered according to R, so first all R strings of type
c 1 comes , then all R strings of type 2 etc .
c
c The R strings are numbered with base for given symmetry and
c type as base, while L is numbered with beginning of symmetry
c as base
c
c
c If ISEL .ne. 0  only the strings in ISELST are searched
c (ISELTP gives number of stringsin ISELST per type )
c
c IF ILSEL .EQ. 0 , L cointains normal (full) string ordering
c If ILSEL .NE. 0 , L contains numbers relating to ISELST


c Input
c=======
c IJ : numbers of two single excitations ( can be zero )
c IJSM  : symmetry of excitation IJ (JI)
c IIJ   : single excitations for each string
c ITO   : excited strings for each
c NEXCI : (Max ) number of excitations per string
c NTYPE : Number of types ( RAS types )
c NSTR  : Total number of strings
c IEXPSM : Base for excitations of given symmetry from given string
c NSMSX  : Number of aymmetries of single excitations
c ILSTRT : string number of first string in L-block
c IRSTRT : string number of first string in R-block
c
c Output
c ======
c L : L array above
c R : R array above
c SIGN : Sign array above
c NLFTP(*) : number of strings in L per type of L
c NRFTP(*) : number of strings in R per type of R
c NTERM :  number of strings in L and R
c IBRFTP : First string in R of given type
c
c Jeppe Olsen , Sept 1989 , LRFEX5 version july 1990
c               Sept 1990 , ADDON added
c
#include <implicit.h>
      INTEGER R
      DIMENSION ITO(*),IIJ(*),NSTFTP(*),ITYPST(*),
     &          IEXPSM(NSMSX+1,*)
      DIMENSION L(*),R(*),SIGN(*),NLFTP(*),NRFTP(*),IBRFTP(*)
      DIMENSION ISELST(*),ISELTP(*)
c
#include <mxpdim.h>
#include <maxash.h>
      COMMON/ADDON/NL1FTP(MXOCTP,2),NL2FTP(MXOCTP,2),NL3FTP(MXOCTP,3),
     &             IJTYP(MAXASH**2),IJTST(9,MXPTP,2)
c
      CALL ISETVC(NLFTP,0,NTYPE)
      NTERM = 0
      IJTP = IJTYP(IJ)
      DO 400 IRTP = 1, NTYPE
        IBRFTP(IRTP) = NTERM + 1
        ILTP = IJTST(IJTP,IRTP,IAB)
c
        IF( ISEL .EQ. 0 ) THEN
          IF( IRTP .EQ. 1 ) THEN
            IRSTR = 1
          ELSE
            IRSTR = IRSTR + NSTFTP(IRTP-1)
          END IF
          IRSTP = IRSTR + NSTFTP(IRTP)-1
        ELSE IF ( ISEL .NE. 0 ) THEN
          IF( IRTP .EQ. 1 ) THEN
            IRSTR = 1
          ELSE
            IRSTR = IRSTR + ISELTP(IRTP-1)
          END IF
          IRSTP = IRSTR + ISELTP(IRTP)-1
        END IF
c
        NTRFTP = 0
        NTERMP = NTERM
        IF(ILTP.EQ.0) GOTO 301
        DO 300 IIR = IRSTR,IRSTP
          IF(ISEL.EQ.0) THEN
            IR = IIR
          ELSE
            IR = ISELST(IIR)
          END IF
c Is E(IJ)|IR> or E(JI)|IR> nonvanishing
          DO 200 IEX = IEXPSM(IJSM,IR+IRSTRT-1),
     &                 IEXPSM(IJSM+1,IR+IRSTRT-1)-1
            IF(IIJ(IEX).EQ.IJ) THEN
              NTRFTP = NTRFTP + 1
              NTERM  = NTERM + 1
              IF(ILSEL.EQ.0) THEN
                R(NTERM) = IR- IRSTR +1
              ELSE
                R(NTERM) = IIR- IRSTR +1
              END IF
c
              IF(ITO(IEX) .GT. NSTR) THEN
                 L(NTERM) = ITO(IEX)-NSTR-ILSTRT+1
                 SIGN(NTERM) = -1.0D0
              ELSE
                 L(NTERM) = ITO(IEX)-ILSTRT + 1
                 SIGN(NTERM) =  1.0D0
              END IF
              GOTO 201
            END IF
 200      CONTINUE
 201      CONTINUE
 300    CONTINUE
 301    CONTINUE
c
        NRFTP(IRTP) = NTRFTP
        IF(ILTP.NE.0) NLFTP(ILTP) = NTRFTP
 400  CONTINUE
c
      RETURN
      END
C  /* Deck oinstr */
#if defined (VAR_SECSEC)
      SUBROUTINE OINSTR(STRING,NSTRIN,NEL,IORB,LSTRIN,ISCAT,
     &                  NTYPE,ITPFST,NSCPTP)
c
c Obtain all strings containing orbital IORB
c
c =======
c  Input
c =======
c
c STRING : Array containing string occupations
c NSTRIN : Number of strings
c NEL    : Number of elecetrons in string
c IORB   : Orbital required to be occupied
c NTYPE  : Number of RAS types
c ITPFST : Type of each string
c
c ========
c. Output
c ========
c
c LSTRIN : NUMBER OF STRINGS WITH IORB OCCUPIED
c ISCAT(I) : EXPANDED NUMBER OF STRING I WITH IORB OCCUPIED
c NSCPTP   : Number of strings on ISCAT per type
c
c JEPPE OLSEN , SLC, DEC 1989, Revised July 1990
c
#include <implicit.h>
#include <priunit.h>
c.Input
      INTEGER STRING(NEL,NSTRIN),ITPFST(NSTRIN)
c.Output
      DIMENSION ISCAT(*),NSCPTP(*)
c
      LSTRIN = 0
      CALL ISETVC(NSCPTP,0,NTYPE)
      DO 200 ISTRIN = 1, NSTRIN
       DO 100 IEL = 1, NEL
         IF(STRING(IEL,ISTRIN).EQ.IORB) THEN
           LSTRIN = LSTRIN + 1
           ISCAT(LSTRIN) = ISTRIN
           NSCPTP(ITPFST(ISTRIN)) = NSCPTP(ITPFST(ISTRIN))+1
           GOTO 200
         END IF
  100  CONTINUE
  200 CONTINUE
c
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(LUPRI,*) ' STRINGS CONTAINING ORBITAL ' ,IORB
        CALL IWRTMA(ISCAT,1,LSTRIN,1,LSTRIN)
        WRITE(LUPRI,*) ' Number of strings per type : '
        CALL IWRTMA(NSCPTP,1,NTYPE,1,NTYPE)

      END IF
      RETURN
      END
#endif
C  /* Deck oinst3 */
      SUBROUTINE OINST3(STRING,NSSO,NEL,IORB,IORBTP,LSTRIN,ISCAT,
     &                  NTYPE,NSCPTP,IAB)
c
c Obtain all strings containing orbital IORB
c
c =======
c  Input
c =======
c
c STRING : Array containing string occupations
c NSSO   : Number of strings per type
c NEL    : Number of elecetrons in string
c IORB   : Orbital required to be occupied
c IORBTP : Type of orbital IORB ( 1,2,3 )
c NTYPE  : Number of RAS types
c IAB    : 1 => alpha string, 2 => beta string
c
c ========
c. Output
c ========
c
c LSTRIN : NUMBER OF STRINGS WITH IORB OCCUPIED
c ISCAT(I) : EXPANDED NUMBER OF STRING I WITH IORB OCCUPIED
c NSCPTP   : Number of strings on ISCAT per type
c
c JEPPE OLSEN , SLC, DEC 1989, Revised July 1990
c                                     Sept. 1990 ( Addon added
c
#include <implicit.h>
#include <priunit.h>
c.Input
      INTEGER STRING(NEL,*),NSSO(NTYPE)
c.Output
      DIMENSION ISCAT(*),NSCPTP(*)
c
#include <mxpdim.h>
#include <maxash.h>
      COMMON/ADDON/NL1FTP(MXOCTP,2),NL2FTP(MXOCTP,2),NL3FTP(MXOCTP,3),
     &             IJTYP(MAXASH**2),IJTST(9,MXPTP,2)

      CALL ISETVC(NSCPTP,0,NTYPE)
      LSTRIN = 0
      DO 300 ITYPE = 1, NTYPE
        NSTRIN = NSSO(ITYPE)
        IF(ITYPE .EQ.1 ) THEN
          IBASE = 1
        ELSE
          IBASE = IBASE + NSSO(ITYPE-1)
        END IF
c
        IF(IORBTP.EQ.1) THEN
          IOFF = 1
          NLFTP = NL1FTP(ITYPE,IAB)
        ELSE IF(IORBTP.EQ.2) THEN
          IOFF = NL1FTP(ITYPE,IAB) + 1
          NLFTP = NL2FTP(ITYPE,IAB)
        ELSE IF(IORBTP.EQ.3) THEN
          IOFF = NL1FTP(ITYPE,IAB)+NL2FTP(ITYPE,IAB) + 1
          NLFTP = NL3FTP(ITYPE,IAB)
        END IF
c
        IF(NLFTP.EQ.1) THEN
          DO 201 ISTRIN = IBASE, IBASE+NSTRIN-1
            IF(STRING(IOFF,ISTRIN).EQ.IORB) THEN
              LSTRIN = LSTRIN + 1
              ISCAT(LSTRIN) = ISTRIN
              NSCPTP(ITYPE) = NSCPTP(ITYPE)+1
            END IF
  201     CONTINUE
        ELSE IF(NLFTP.EQ.2) THEN
          DO 202 ISTRIN = IBASE, IBASE+NSTRIN-1
            IF(STRING(IOFF,ISTRIN)  .EQ.IORB  .OR.
     &         STRING(IOFF+1,ISTRIN).EQ.IORB) THEN
              LSTRIN = LSTRIN + 1
              ISCAT(LSTRIN) = ISTRIN
              NSCPTP(ITYPE) = NSCPTP(ITYPE)+1
            END IF
  202     CONTINUE
        ELSE IF (NLFTP.GE.3) THEN
          DO 203 ISTRIN = IBASE,IBASE+NSTRIN-1
            DO 100 IEL = 1, NLFTP
              IF(STRING(IOFF-1+IEL,ISTRIN)  .EQ.IORB) THEN
                LSTRIN = LSTRIN + 1
                ISCAT(LSTRIN) = ISTRIN
                NSCPTP(ITYPE) = NSCPTP(ITYPE)+1
              END IF
  100       CONTINUE
  203     CONTINUE
        END IF
  300 CONTINUE
c
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(LUPRI,*) ' STRINGS CONTAINING ORBITAL ' ,IORB
        CALL IWRTMA(ISCAT,1,LSTRIN,1,LSTRIN)
        WRITE(LUPRI,*) ' Number of strings per type : '
        CALL IWRTMA(NSCPTP,1,NTYPE,1,NTYPE)

      END IF
      RETURN
      END
C  /* Deck zijtyp */
      SUBROUTINE ZIJTYP(IJTYP,ITPFOB,NORB)
c
c Obtain type of each excitation
c
c.Input
#include <implicit.h>
#include <priunit.h>
      DIMENSION ITPFOB(NORB)
c.Output
      DIMENSION IJTYP(NORB * NORB )
c
      DO 100 IORB = 1, NORB
        ITYP = ITPFOB(IORB)
        DO 50 JORB = 1, NORB
          JTYP =  ITPFOB(JORB)
          IJ = (IORB-1)*NORB + JORB
          IJTYP(IJ) = (ITYP-1)*3 + JTYP
   50   CONTINUE
  100 CONTINUE
c
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' IJTYP array '
        CALL IWRTMA(IJTYP,NORB,NORB,NORB,NORB)
      END IF
c
      RETURN
      END
C  /* Deck sxtst */
      SUBROUTINE SXTST(NTYPE,NL1FTP,NL2FTP,NL3FTP,IJTST)
c
c Obtain type of excitation times string .
c The type of an excitation is defined  as
c (ITP-1)*3 + JTP , where ITP,JTP are the types of orbitals
c (1 to 3 ) in excitation IJ
c
c.Input
#include <implicit.h>
#include <priunit.h>
      DIMENSION NL1FTP(NTYPE),NL2FTP(NTYPE),NL3FTP(NTYPE)
c.Output
      DIMENSION IJTST(9,*)
c
      DO 300 ITP = 1, 3
       IF(ITP.EQ.1) THEN
        IDEL1 = 1
        IDEL2 = 0
        IDEL3 = 0
       ELSE IF (ITP.EQ.2) THEN
        IDEL1 = 0
        IDEL2 = 1
        IDEL3 = 0
       ELSE IF (ITP.EQ.3) THEN
        IDEL1 = 0
        IDEL2 = 0
        IDEL3 = 1
       END IF
       DO 200 JTP = 1, 3
        IF(JTP.EQ.1) THEN
         JDEL1 =-1
         JDEL2 = 0
         JDEL3 = 0
        ELSE IF (JTP.EQ.2) THEN
         JDEL1 = 0
         JDEL2 =-1
         JDEL3 = 0
        ELSE IF (JTP.EQ.3) THEN
         JDEL1 = 0
         JDEL2 = 0
         JDEL3 =-1
        END IF
        IJDEL1 = IDEL1 + JDEL1
        IJDEL2 = IDEL2 + JDEL2
        IJDEL3 = IDEL3 + JDEL3
        IJTP = (ITP-1)*3 + JTP
        DO 100 ISTTP = 1, NTYPE
cOccupation of IJ * string
          JSEL1 = JDEL1 + NL1FTP(ISTTP)
          JSEL2 = JDEL2 + NL2FTP(ISTTP)
          JSEL3 = JDEL3 + NL3FTP(ISTTP)
          IF(JSEL1.LT.0.OR.JSEL2.LT.0.OR.JSEL3.LT.0) GOTO 100
          IJSEL1 = IJDEL1 + NL1FTP(ISTTP)
          IJSEL2 = IJDEL2 + NL2FTP(ISTTP)
          IJSEL3 = IJDEL3 + NL3FTP(ISTTP)
          IJTYPE = 0
          DO 50 JTYPE = 1, NTYPE
           IF((NL1FTP(JTYPE).EQ.IJSEL1).AND.
     &        (NL2FTP(JTYPE).EQ.IJSEL2).AND.
     &        (NL3FTP(JTYPE).EQ.IJSEL3) ) IJTYPE = JTYPE
   50     CONTINUE
          IJTST(IJTP,ISTTP) = IJTYPE
  100   CONTINUE
  200  CONTINUE
  300  CONTINUE
c
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' SX-type * ST-type => ST-type '
        WRITE(LUPRI,*) ' ============================ '
        CALL IWRTMA(IJTST,9,NTYPE,9,NTYPE)
      END IF
c
      RETURN
      END
C  /* Deck nsexci */
      FUNCTION NSEXCI(NTYPE,NSTFTP,NL1FTP,NL2FTP,NL3FTP,
     &                  NORB1,NORB2,NORB3,IJTST,NTEST)
c
c Obtain total number of allowed single excitations
c
c.Input
#include <implicit.h>      
#include <priunit.h>
      DIMENSION NL1FTP(NTYPE),NL2FTP(NTYPE),NL3FTP(NTYPE)
      DIMENSION NSTFTP(NTYPE),IJTST(9,NTYPE)

c
      NELEC = NL1FTP(1) + NL2FTP(1) + NL3FTP(1)
      NN = 0
      DO 500 ISTYP=1,NTYPE
c.Off diagonal excitations
        DO 300 ITP = 1, 3
          IF(ITP.EQ.1) THEN
           IOCC = NL1FTP(ISTYP)
          ELSE IF (ITP.EQ.2) THEN
           IOCC = NL2FTP(ISTYP)
          ELSE IF (ITP.EQ.3) THEN
           IOCC = NL3FTP(ISTYP)
          END IF
          DO 200 JTP = 1, 3
            IF(JTP.EQ.1) THEN
              IUNOCC = NORB1 - NL1FTP(ISTYP)
            ELSE IF (JTP.EQ.2) THEN
              IUNOCC = NORB2 - NL2FTP(ISTYP)
            ELSE IF (JTP.EQ.3) THEN
              IUNOCC = NORB3 - NL3FTP(ISTYP)
            END IF
            JITP = (JTP-1)*3 + ITP
            IF(IJTST(JITP,ISTYP).NE.0)
     &      NN = NN + NSTFTP(ISTYP)*IOCC*IUNOCC
  200  CONTINUE
  300  CONTINUE
c. Diagonal excitations
       NN = NN + NELEC*NSTFTP(ISTYP)
  500 CONTINUE

c
      NSEXCI = NN
c
      IF(NTEST .GE. 5 ) THEN
        WRITE(LUPRI,*) ' Total number of excitations from strings ',
     &        NSEXCI
        WRITE(LUPRI,*)
     &       ' ================================================='
      END IF
c
      RETURN
      END
C  /* Deck oinstd */
      SUBROUTINE OINSTD(KORB,ISM,NFTP,NSTRIN,ISTRIN,IAB,XNDXCI)
c
c Obtain all strings of symmetry ISM with orbital KORB occupied.
c
c =====
c Input
c =====
c      KORB : orbital that should be occupied
c      ISM : Symmetry of strings
c      IAB  : 1 => alpha strings 2 => beta strings
c      XNDXCI : String info array
c
c =======
c Output
c =======
c
c      NFTP(ITYPE) : Number of strings of type ITYPE containing KORB
c      NSTRIN      : Total number of strings of sym ISM containing IORB
c      ISTRIN      : Strings containing KORB
c
#include <implicit.h>
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
C
      DIMENSION XNDXCI(*)
      DIMENSION ISTRIN(*),NFTP(*)
c
      IF(IAB.EQ.1) THEN
        KNSSO  = KNSSOA
        KISSO  = KISSOA
        NOCTP  = NOCTPA
        KSTRIN = KIASTR
        NTYPE  = NOCTPA
        NEL    = NAEL
      ELSE
        KNSSO  = KNSSOB
        KISSO  = KISSOB
        NOCTP  = NOCTPB
        KSTRIN = KIBSTR
        NTYPE  = NOCTPB
        NEL    = NBEL
      END IF
c
      CALL OINST2(XNDXCI(KSTRIN),XNDXCI(KNSSO),XNDXCI(KISSO),
     &            NTYPE,NEL,KORB,NSTRIN,ISTRIN,NFTP,ISM)

c
      RETURN
      END
C  /* Deck oinst2 */
      SUBROUTINE OINST2(STRING,NSSO,ISSO,NTYPE,NEL,IORB,LSTRIN,ISCAT,
     &                  NSCPTP,ISYM)
c
c Obtain all strings containing orbital IORB of sym ISYM
c
c =======
c  Input
c =======
c
c STRING : Array containing string occupations
c NSSO   : Number of strings per sym and type
c ISSO   : first string of given type and sym
c NEL    : Number of elecetrons in string
c IORB   : Orbital required to be occupied
c NTYPE  : Number of RAS types
c ISYM   : required symmetry of string
c
c ========
c. Output
c ========
c
c LSTRIN : NUMBER OF STRINGS WITH IORB OCCUPIED
c ISCAT(I) : EXPANDED NUMBER OF STRING I WITH IORB OCCUPIED
c NSCPTP   : Number of strings on ISCAT per type
c
c JEPPE OLSEN , SLC, DEC 1989, Revised July 1990
c
#include <implicit.h>
#include <priunit.h>
c.Input
      INTEGER STRING(NEL,*),NSSO(NTYPE,*),ISSO(NTYPE,*)
c.Output
      DIMENSION ISCAT(*),NSCPTP(*)
c
      LSTRIN = 0
      DO 300 ITYPE = 1, NTYPE
        IF(ITYPE.EQ.1) THEN
          IBASE = ISSO(1,ISYM)
        ELSE
          IBASE = IBASE + NSSO(ITYPE-1,ISYM)
        END IF
        NPTP = 0
        DO 200 ISTRIN = IBASE,IBASE + NSSO(ITYPE,ISYM)-1
         DO 100 IEL = 1, NEL
           IF(STRING(IEL,ISTRIN).EQ.IORB) THEN
             LSTRIN = LSTRIN + 1
             ISCAT(LSTRIN) = ISTRIN
             NPTP = NPTP + 1
             GOTO 200
           END IF
  100    CONTINUE
  200   CONTINUE
        NSCPTP(ITYPE) = NPTP
  300 CONTINUE
c
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(LUPRI,*) ' STRINGS CONTAINING ORBITAL ' ,IORB
        CALL IWRTMA(ISCAT,1,LSTRIN,1,LSTRIN)
        WRITE(LUPRI,*) ' Number of strings per type : '
        CALL IWRTMA(NSCPTP,1,NTYPE,1,NTYPE)

      END IF
      RETURN
      END
C  /* Deck rascg6 */
      SUBROUTINE RASCG6(C,CG,NTYPR,NFRFTP,NGRFTP,NTYPC,NFCFTP,
     &                  IBOUT,IG,SIGN,IOCOC,MXGCTP,ITPREO,IGTP,
     &                  IBIN,MXPTP)

C
C     PURPOSE: GATHER A SET OF RAS BLOCKS FOR A GIVEN SYMMETRY
C              CG(I,J) = C(IG(I),J)*SIGN(I)
C
C     =====
C     INPUT
C     =====
C     C : FULLS SYMMETRY BLOCK OF RAS VECTOR
C     NTYPR  : NUMBER OF TYPES OF ROWS
C     NFRFTP : NUMBER OF ROWS PER TYPE FOR FULL MATRIX
C     NGRFTP : NUMBER OF ROWS PER TYPE FOR GATHERED MATRIX
C     NTYPC  : NUMBER OF TYPES OF COLUMNS
C     NFCFTP : NUMBER OF COLUMNS PER SYMETRY FOR FULL MATRIX
C     IG     : GATHERING VECTOR
C     SIGN   : SIGN SHIFT ARRAY
C     IOCOC  : ARRAY TELLING WHETHER A GIVEN BLOCK IS INCLUDED OR NOT
C     ITPREO : FLAG TO REORDER ROW TYPES
C              ITPREO=0 => NO REORDERING OF TYPES
C              ITPREO=1 => REORDER ALSO ROW TYPES
C     IGTP   : TYPE REORDERING IGTP(ITP) IS REORDERED TYPE FOR ORIGINAL
C              TYPE ITP
C     IBIN   : POINTER ARRAY TO START OF GIVEN OO BLOCK IN C
C
C     ======
C     OUTPUT
C     ======
C     CG     : GATHERED MATRIX AS ABOVE
C     IBOUT  : ADRESS OF FIRST ELEMENT IN CG BELONGING TO A
C              SUBBLOCK DEFINED BY GIVEN ROW AND COLUMN TYPE
C     MXGCTP : LAST COLUMN TYPE WITH A NONVANISHING NUMBER OF ROWS
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION C(*),CG(*),NFRFTP(*),NGRFTP(*),NFCFTP(*),
     &          IBOUT(MXPTP,*),IG(*),SIGN(*),
     &          IOCOC(NTYPR,NTYPC),IGTP(*),IBIN(MXPTP,*)
C
      MXGCTP = 0
      IBSGC = 1
      IBSFC = 1
C     LOOP OVER ROW TYPE OF SCATTERED MATRIX
      DO  200 IRTP = 1, NTYPR
C       ROW TYPE OF ORIGINAL MATRIX
        IF(ITPREO.EQ.0) THEN
          IFRTP = IRTP
          IF(IRTP.EQ.1) THEN
            IBFR = 1
          ELSE
            IBFR = IBFR + NFRFTP(IFRTP-1)
          END IF
        ELSE
          IFRTP = 0
          IBFR  = 1
          DO 10 ITYP = 1, NTYPR
            IF(ITYP.EQ.1) THEN
              IBFR = 1
            ELSE
              IBFR = IBFR + NFRFTP(ITYP-1)
            END IF
            IF(IGTP(ITYP).EQ.IRTP) THEN
              IFRTP = ITYP
              GOTO 11
            END IF
   10     CONTINUE
   11     CONTINUE
          IF( IFRTP.EQ.0) GOTO 200
        END IF
C
        NFR = NFRFTP(IFRTP)
        NGR = NGRFTP(IRTP)
        IF(IRTP .EQ. 1 ) THEN
          IBGR = 1
        ELSE
          IBGR = IBGR + NGRFTP(IRTP-1)
        END IF
       DO  150 ICTP = 1, NTYPC
         IBOUT(IRTP,ICTP) = IBSGC
         IF ( IOCOC(IRTP,ICTP) .EQ. 1 ) THEN
C
            NFC = NFCFTP(ICTP)
            IF(NFC*NGR .NE. 0 ) MXGCTP = MAX(MXGCTP,ICTP)
            IF(NGR.EQ.0)GOTO  140
            IICI = IBSFC-1-NFR
            IICO = IBSGC-IBGR-NGR
C
            DO 21 IC = 1, NFC
               IICO = IICO + NGR
               IICI = IICI + NFR
               DO 20 IR = IBGR, IBGR+NGR-1
                    CG(IICO+IR ) = C(IICI + IG(IR))*SIGN(IR)
   20          CONTINUE
   21       CONTINUE
C
  140       CONTINUE
            IBSFC = IBSFC + NFR*NFC
            IBSGC = IBSGC + NGR*NFC
          END IF
C
  150   CONTINUE
  200 CONTINUE
C
      NTEST = 0
      IF( NTEST.NE.0 ) THEN
        WRITE(LUPRI,*) 'GATHERED POINTER ARRAY '
        CALL IWRTMA(IBOUT,NTYPR,NTYPC,MXPTP,MXPTP)
      END IF
C
      RETURN
      END
C  /* Deck trrss2 */
      SUBROUTINE TRRSS2(CIN,COUT,NOCTPA,NOCTPB,IOCOC,NSSOA,NSSOB,
     &                  IASM,IBSM,ICOUTP,ICINP,IWAY,MXPTP)
C
C     PURPOSE: TRANSPOSE MATRIX WITH SELECTION, I.E., A POINTER
C              VECTOR IS GENERATED IN THE FORWARD TRANSPOSITION
C              AND USED AGAIN WHEN BACK TRANSPOSING.
C
C     NOTE :   TRRSS2 DIFFERS FROM TRRSSB
C              BY HAVING LOGICAL ORDER | OF OFFSET ARRAYS
C              AND HAVING ICINP REFERRING TO A GIVEN BLOCK
C              IWAY IS ALSO AN ADDITION
C
#include <implicit.h>
C
      DIMENSION CIN(*),COUT(*)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION NSSOA(NOCTPA,*),NSSOB(NOCTPB,*)
      DIMENSION ICOUTP(MXPTP,MXPTP),ICINP(MXPTP,MXPTP)
C
C     TRANSPOSE FORWARD
C
      IF(IWAY .EQ. 1 ) THEN
        IBASE = 1
        DO  200 IBTP = 1, NOCTPB
          DO  100 IATP = 1, NOCTPA
            ICOUTP(IATP,IBTP) = IBASE
            IF(IOCOC(IATP,IBTP).EQ. 1 ) THEN
              NRIN = NSSOA(IATP,IASM)
              NCIN = NSSOB(IBTP,IBSM)
              CALL TRPMAT(CIN(ICINP(IATP,IBTP)),
     &                    NRIN,NCIN,COUT(IBASE) )
              IBASE = IBASE + NRIN * NCIN
            END IF
  100     CONTINUE
  200   CONTINUE
      ELSE
C
C     BACK TRANSPOSE
C
        DO  400 IBTP = 1, NOCTPB
          DO  300 IATP = 1, NOCTPA
            IF(IOCOC(IATP,IBTP).EQ. 1 ) THEN
              NROUT = NSSOB(IBTP,IBSM)
              NCOUT = NSSOA(IATP,IASM)
              IOUT = ICOUTP(IATP,IBTP)
              IIN = ICINP(IATP,IBTP)
              CALL TRPMAT(COUT(IOUT),NROUT,NCOUT,CIN(IIN) )
            END IF
  300     CONTINUE
  400   CONTINUE
       END IF
C
      RETURN
      END
C  /* Deck lrfexd */
      SUBROUTINE LRFEXD(IJ,JI,IRSM,L,R,SIGN,NLFTP,NRFTP,ILTFRT,
     &                  NTERM,IBRFTP,IAB,SXSTSM,XNDXCI,NSMSX,
     &                  ISEL,NFTPK,KOC)
C
C      PURPOSE: DRIVER FOR OBTAINING STRINGS THAT ARE CONNECTED
C               WITH A GIVEN EXCITATION.
C               |L(I)> = SIGN(I) E |R(I)>, WHERE E IS E(IJ) OR E(JI)
C               R-STRINGS ARE RESTRICTED TO BE OF SYMMETRY IRSM
C
C
C              If ISEL .NE. 0 only the strings in KOC are checked
C      =====
C      INPUT
C      =====
C
C      IJ,JI  : NUMBER OF TWO EXCITATIONS (IF ONLY ONE IS NEEDED, ZERO
C               IS USED FOR THE OTHER COMPONENT )
C      IRSM   : SYMMETRY OF R-STRINGS
C      SXSTSM : SX-SYM * ST-SYM => ST-SYM
C      XNDXCI : ARRAY CONTAINING STRING INFORMATION
C
C      ======
C      OUTPUT
C      ======
C      L      : L(I) ARRAY AS ABOVE
C      R      : R(I) ARRAY AS ABOVE, R(I) IS IN ASCENDING ORDER
C               L AND R ARE ENUMERATED WITH RESPECT TO START OF SO CLASS
C      SIGN   : SIGN(I) ARRAY AS ABOVE
C      NLFTP  : NUMBER OF L-STRINGS PER TYPE
C      NRFTP  : NUMBER OF R-STRINGS PER TYPE
C      ILTFRT : TYPE OF LEFT STRING FOR GIVEN TYPE OF RIGHT STRING
C      NTERM  : NUMBER OF L,R PAIRS OBTAINED
C      IBRFTP : OFFSET FOR R STRINGS OF A GIVEN TYPE
C      IAB    : = 1 => A-STRINGS, =2 => B-STRINGS.
C
#include <implicit.h>
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
C
      DIMENSION XNDXCI(*), SIGN(*)
      INTEGER SXSTSM(*)
      INTEGER L(*),R(*),NLFTP(*),NRFTP(*),ILTFRT(*),IBRFTP(*)
      INTEGER NFTPK(*),KOC(*)
C
      IF(IAB.EQ.1 ) THEN
         NEXCI = NAEXCI
         KNSO  = KNSSOA
         KISO  = KISSOA
         NTYPE = NOCTPA
         NSTR = NASTR
         KTPFS = KTPFSA
         KTIJ = KTAIJ
         KTTO = KTATO
         KTSYM = KTASYM
      ELSE
         NEXCI = NBEXCI
         KNSO  = KNSSOB
         KISO  = KISSOB
         NTYPE = NOCTPB
         NSTR = NBSTR
         KTPFS = KTPFSB
         KTIJ = KTBIJ
         KTTO = KTBTO
         KTSYM = KTBSYM
      END IF
C
      CALL LRFEX2(IJ,JI,IRSM,XNDXCI(KISSYM),L,R,SIGN,NLFTP,NRFTP,
     &            ILTFRT,XNDXCI(KTPFS),
     &            XNDXCI(KTIJ),XNDXCI(KTTO),NEXCI,NTERM,XNDXCI(KNSO),
     &            XNDXCI(KISO),
     &            NTYPE,NSTR,XNDXCI(KTSYM),NSMSX,IBRFTP,SXSTSM,
     &            ISEL,NFTPK,KOC,IAB)
C
      RETURN
      END
C  /* Deck lrfex2 */
      SUBROUTINE LRFEX2(IJ,JI,IRSM,DISSYM,L,R,SIGN,NLFTP,NRFTP,ILTFRT,
     &           ITYPST,IIJ,ITO,NEXCI,NTERM,NSSO,ISSO,NTYPE,NSTR,IEXPSM,
     &           NSMSX,IBRFTP,SXSTSM,ISEL,NFTPK,KOC,IAB)
C
C     PURPOSE: FOR TWO SYMMETRY BLOCKS OF STRINGS SET UP
C              <L| E |R> FOR EXCITATIONS E(IJ)
C              IN THE FORM OF THREE ARRAYS L,R,ISIGN
C              |L(I)> = E |R(I)> * SIGN(I)
C              THE LIST IS ORDERED ACCORDING TO R, SO ALL R
C              STRINGS OF TYPE 1 COME FIRST, ALL OF TYPE2 SECOND, ETC.
C              THE L AND R STRINGS ARE NUMBERED FROM BASE OF GIVEN
C              SYMMETRY AND TYPE.
C
C     =====
C     INPUT
C     =====
C     IJ :    SINGLE EXCITATION
C     IRSM  : SYMMETRY OF R-STRINGS
C     DISSYM: SYMMETRY OF EXCITATIONS
C     IIJ   : SINGLE EXCITATIONS FOR EACH STRING
C     ITO   : EXCITED STRINGS FOR EACH
C     NEXCI : (MAX ) NUMBER OF EXCITATIONS PER STRING
C     NSSO  : NUMBER OF STRINGS PER TYPE AND SYMMETRY
C     ISSO  : OFFSET FOR STRINGS OF GIVEN TYPE AND SYMMETRY
C     NTYPE : NUMBER OF TYPES ( RAS TYPES )
C     NSTR  : TOTAL NUMBER OF STRINGS
C     IEXPSM: BASE FOR EXCITATIONS OF GIVEN SYMMETRY FROM GIVEN STRING
C     NSMSX : NUMBER OF SYMMETRIES OF SINGLE EXCITATIONS
C
C     ======
C     OUTPUT
C     ======
C     L       : L ARRAY ABOVE
C     R       : R ARRAY ABOVE
C     SIGN    : SIGN ARRAY ABOVE
C     NLFTP(*): NUMBER OF STRINGS IN L PER TYPE OF L
C     NRFTP(*): NUMBER OF STRINGS IN R PER TYPE OF R
C     ILTFRT  : LEFT TYPE CORRESONDING TO GIVEN R-TYPE
C     NTERM   :  NUMBER OF STRINGS IN L AND R
C     IBRFTP  : FIRST STRING IN R OF GIVEN TYPE
C
#include <implicit.h>
C
#include <priunit.h>
#include <mxpdim.h>
#include <maxash.h>
      COMMON/ADDON/NL1FTP(MXOCTP,2),NL2FTP(MXOCTP,2),NL3FTP(MXOCTP,3),
     &             IJTYP(MAXASH**2),IJTST(9,MXPTP,2)
      INTEGER   R,DISSYM(*),SXSTSM(NSMSX,*)
      DIMENSION NSSO(NTYPE,*),ISSO(NTYPE,*)
      DIMENSION ITO(*),IIJ(*),ITYPST(*),
     &          IEXPSM(NSMSX+1,*)
      DIMENSION L(*),R(*),SIGN(*),NLFTP(*),NRFTP(*),IBRFTP(*),
     &          ILTFRT(*)
      DIMENSION NFTPK(*),KOC(*)
C
      IJSM = DISSYM(IJ)
      IJTP = IJTYP(IJ)
      ILSM = SXSTSM(IJSM,IRSM)
      ILSTRT = ISSO(1,ILSM)
C
      NTERM = 0
      CALL ISETVC(ILTFRT,0,NTYPE)
      CALL ISETVC(NLFTP,0,NTYPE)
      IF(ILSM.EQ.0) GOTO 401
      DO 400 IRTP = 1, NTYPE
        IBRFTP(IRTP) = NTERM + 1
C
        IF(ISEL.EQ.0) THEN
          IF( IRTP .EQ. 1 ) THEN
            IRSTR = ISSO(1,IRSM)
          ELSE
            IRSTR = IRSTR + NSSO(IRTP-1,IRSM)
          END IF
          IRSTP = IRSTR + NSSO(IRTP,IRSM)-1
        ELSE
          IF(IRTP.EQ.1) THEN
            IRSTR = 1
          ELSE
            IRSTR = IRSTR + NFTPK(IRTP-1)
          END IF
          IRSTP = IRSTR + NFTPK(IRTP)-1
        END IF
C
        NTRFTP = 0
        NTERMP = NTERM
        IROFF = ISSO(IRTP,IRSM)
        ILTP = IJTST(IJTP,IRTP,IAB)
        IF(ILTP.EQ.0) GOTO 301
        DO 300 IIR = IRSTR,IRSTP
C         IS E(IJ)|IR> OR E(JI)|IR> NONVANISHING
          IF(ISEL.EQ.0) THEN
            IR = IIR
          ELSE
            IR = KOC(IIR)
          END IF
          DO 200 IEX = IEXPSM(IJSM,IR),IEXPSM(IJSM+1,IR)-1
            IF(IIJ(IEX).EQ.IJ) THEN
              NTRFTP = NTRFTP + 1
              NTERM  = NTERM + 1
              R(NTERM) = IR- IROFF +1
              IF(ITO(IEX) .GT. NSTR) THEN
                 L(NTERM) = ITO(IEX)-NSTR+  1
     &                    - ISSO(ITYPST(ITO(IEX)-NSTR),ILSM)
                 SIGN(NTERM) = -1.0D0
                 LTYPE = ITYPST(ITO(IEX)-NSTR)
              ELSE
                 L(NTERM) = ITO(IEX) + 1
     &                    - ISSO(ITYPST(ITO(IEX)),ILSM)
                 SIGN(NTERM) =  1.0D0
                 LSTRT = ISSO(ITYPST(ITO(IEX)),ILSM)
                 LTYPE = ITYPST(ITO(IEX))
              END IF
              GOTO 201
            END IF
 200      CONTINUE
 201      CONTINUE
 300    CONTINUE
 301    CONTINUE
C
        NRFTP(IRTP) = NTRFTP
        IF(NTRFTP.NE.0) THEN
          ILTFRT(IRTP) = LTYPE
          NLFTP(ILTFRT(IRTP)) = NTRFTP
        END IF
 400  CONTINUE
 401  CONTINUE
C
      RETURN
      END
C  /* Deck smblk */
      SUBROUTINE SMBLK(NOCTPA,NOCTPB,IOCOC,NSOA,NSOB,NDTBLK,IOFF,
     &                 NDTPBL,MXPTP,MXATP,MXBTP)
C
C     PURPOSE: NUMBER OF DETERMINANT AND OFFSET VECTOR
C              FOR GIVEN SYMMETRY BLOCK
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION IOCOC(NOCTPA,NOCTPB),NSOA(NOCTPA),NSOB(NOCTPB)
      DIMENSION IOFF(MXPTP,MXPTP),NDTPBL(MXPTP,MXPTP)
C
      CALL ISETVC(NDTPBL,0,MXPTP ** 2 )
      NDTBLK = 0
      IOFF(1,1) = 1
      MXATP = 0
      MXBTP = 0
      DO 100 IATP = 1, NOCTPA
      DO 100 IBTP = 1, NOCTPB
        IOFF(IATP,IBTP) = NDTBLK + 1
        IF(IOCOC(IATP,IBTP) .EQ. 1 )  THEN
          NDTBLK = NDTBLK + NSOA(IATP)*NSOB(IBTP)
          NDTPBL(IATP,IBTP) = NSOA(IATP)*NSOB(IBTP)
          IF(NDTPBL(IATP,IBTP).NE.0) THEN
            MXATP = MAX(MXATP,IATP)
            MXBTP = MAX(MXBTP,IBTP)
          END IF
        END IF
  100 CONTINUE
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' NUMBER OF DETERMINANTS IN BLOCK ', NDTBLK
        WRITE(LUPRI,*) ' OFFSET MATRIX FOR GIVEN BLOCK '
        CALL IWRTMA(IOFF,NOCTPA,NOCTPB,MXPTP,MXPTP)
        WRITE(LUPRI,*) ' NUMBER OF DETERMINANTS PER BLOCK '
        CALL IWRTMA(NDTPBL,NOCTPA,NOCTPB,MXPTP,MXPTP)
      END IF
C
      RETURN
      END
C  /* Deck gtsxst */
      SUBROUTINE GTSXST(ISTR,ICLS,ISYM,JSYM,JCLS,NSX,IIJSX,ITOSX,
     &                  SGSX,JCLSX,XNDXCI,IAB,NSMST,IJMIN,NSMSX,SXSTSM,
     &                  MXTP,IGTJ)
C
C     PURPOSE: OBTAIN SINGLE EXCITATIONS OUT FROM STRING ISTR,ICLS,ISYM
C              INTO STRINGS OF SYMMETRY JSYM AND CLASS JCLS.
C              IF JCLS .EQ. 0 ALL EXCITATION CLASSES ARE INTRODUCED
C              ONLY EXCITATIONS LARGER THAN IJMIN ARE INCLUDED
C              LARGEST allowed excitation class is MXTP
C              IF(IGTJ.NE.0.only excitations E(ij) with i.ge.j are
C              included
C
#include <implicit.h>
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
C
      DIMENSION XNDXCI(*)
      INTEGER SXSTSM(NSMSX,*)
      DIMENSION IIJSX(*),ITOSX(*),SGSX(*),JCLSX(*)
C
      IF(IAB.EQ.1) THEN
        KTIJ = KTAIJ
        KTTO = KTATO
        KTSYM = KTASYM
        KTPFST = KTPFSA
        NSTRIN = NASTR
        NEXCI = NAEXCI
        KIOFST = KISSOA
        NOCTP = NOCTPA
      ELSE
        KTIJ = KTBIJ
        KTTO = KTBTO
        KTSYM = KTBSYM
        KTPFST = KTPFSB
        NSTRIN = NBSTR
        NEXCI = NBEXCI
        KIOFST = KISSOB
        NOCTP = NOCTPB
      END IF
C
      CALL GTSXS2(ISTR,ICLS,ISYM,JSYM,JCLS,NSX,IIJSX,ITOSX,
     &            SGSX,JCLSX,NSMST,NSTRIN,NEXCI,NOCTP,
     &            XNDXCI(KIOFST),XNDXCI(KTSYM),XNDXCI(KTIJ),
     &            XNDXCI(KTTO),XNDXCI(KTPFST),IJMIN,NSMSX,SXSTSM,MXTP,
     &            XNDXCI(KKLTP),IGTJ)
C
      RETURN
      END
C  /* Deck gtsxs2 */
      SUBROUTINE GTSXS2(ISTR,ICLS,ISYM,JSYM,JCLS,NSX,IIJSX,ITOSX,
     &                  SGSX,JCLSX,NSMST,NSTRIN,NEXCI,NOCTP,
     &                  IOFST,IOFXST,IJ,ITO,ITPFST,IJMIN,NSMSX,SXSTSM,
     &                  MXTP,KLTRP,IGTJ)
C
C     PURPOSE: OBTAIN EXCITATIONS OUT FROM STRING ISTR,ICLS,ISYM
C              TO STRINGS OF SYMMETRY JSYM AND TYPE JCLS.
C              IF JCLS .EQ. 0 NO RESTRICTIONS ON THE THE CLASS OF
C              JSTRINGS IS IMPOSED, AND THE NUMBER OF EXCITATIONS
C              PER JCLASS IS STORED IN JCLSX.
C              THE J-STRINGS ARE ENUMERATED WITH START OF GIVEN
C              SYMORB TYPE AS BASE.
C               ONLY EXCITATIONS LARGER THAN IJMIN ARE INCLUDED
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION IOFST(NOCTP,*)
      DIMENSION IOFXST(NSMSX+1,*),IJ(*),ITO(*)
      DIMENSION ITPFST(*)
      INTEGER SXSTSM(NSMSX,NSMST)
      DIMENSION IIJSX(*),ITOSX(*),SGSX(*),JCLSX(*)
      DIMENSION KLTRP(*)
C
      ISTRAB = IOFST(ICLS,ISYM)-1+ISTR
      IF(JCLS.EQ.0) CALL ISETVC(JCLSX,0,NOCTP)
C     EXCITATION SYMMETRY SO : |JSYM> = E(IJSM) |ISYM>
      IJSM = 0
      DO 10 IIJSM = 1, NSMSX
         IF(SXSTSM(IIJSM,ISYM).EQ.JSYM) IJSM = IIJSM
   10 CONTINUE
      IF(IJSM.EQ.0) GOTO 101
C
      NSX = 0
      DO 100 IIJ = IOFXST(IJSM,ISTRAB),IOFXST(IJSM+1,ISTRAB)-1
        JSTR =  ITO(IIJ)
        IJEX =  IJ (IIJ)
        IF(IGTJ.NE.0.AND.KLTRP(IJEX).GT.IJEX) GOTO 100
        IF(IJEX.GE.IJMIN) THEN
          SIGN = 1.0D0
          IF( JSTR.GT.NSTRIN) THEN
            JSTR = JSTR - NSTRIN
            SIGN = -1.0D0
          END IF
          JJCLS = ITPFST(JSTR)
          IF(JJCLS .GT. MXTP ) GOTO 101
          IF(JCLS.EQ.0) THEN
            JSTREL = JSTR - IOFST(JJCLS,JSYM)+1
            NSX = NSX + 1
            ITOSX(NSX) = JSTREL
            SGSX(NSX) = SIGN
            IIJSX(NSX) = IJ(IIJ)
            JCLSX(JJCLS) = JCLSX(JJCLS) + 1
          ELSE IF(JJCLS.EQ.JCLS) THEN
            JSTREL = JSTR - IOFST(JJCLS,JSYM)+1
            NSX = NSX + 1
            ITOSX(NSX) = JSTREL
            SGSX(NSX) = SIGN
            IIJSX(NSX) = IJ(IIJ)
          END IF
        END IF
  100 CONTINUE
  101 CONTINUE
C
      NTEST=0
      IF( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' GTSXS2 : EXCITATIONS OUT FROM STRING(ABS) ',
     &        ISTRAB
        WRITE(LUPRI,*) 
     &       ' ================================================ '
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' REQUIRED SYMMETRY OF EXCITED STRING ', JSYM
        WRITE(LUPRI,*) ' ISTR ICLS ISYM ',ISTR,ICLS,ISYM
        IF(JCLS.NE.0) WRITE(LUPRI,*)
     &  ' REQUIRED CLASS OF EXCITED STRING ', JCLS
        WRITE(LUPRI,*) ' NUMBER OF EXCITATIONS OF OBTAINED ',NSX
        IF(JCLS.EQ.0) THEN
          WRITE(LUPRI,*) ' NUMBER OF EXCITATIONS PER J CLASS '
          CALL IWRTMA(JCLSX,1,NOCTP,1,NOCTP)
        END IF
        WRITE(LUPRI,*) ' EXCITED STRING, EXCITATION AND SIGN '
        DO 1000 ISX = 1, NSX
          WRITE(LUPRI,'(2X,I5,I5,F8.3)') ITOSX(ISX),IIJSX(ISX),SGSX(ISX)
 1000   CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck sxstm */
      SUBROUTINE SXSTSM(SXSTM,NSMSX,NSMST)
C
C     PURPOSE: OBTAIN MATRIX GIVING PRODUCT SYMMETRY OF
C              EXCITATION TIMES STRING
C
#include <implicit.h>
C
      INTEGER SXSTM(NSMSX,NSMST)
C
      DO 100 ISTSM = 1, NSMST
      DO 100 ISXSM = 1, NSMSX
        SXSTM(ISXSM,ISTSM) = 1+IEOR(ISXSM-1,ISTSM-1)
  100 CONTINUE
C
      RETURN
      END
C  /* Deck zsmost */
      SUBROUTINE ZSMOST(ITOTSM,ISMOST,NSMST)
C
C     PURPOSE: A CI EXPANSION OF SYMMETRY ITOTSM IS GIVEN
C              CONSTRUCT ISMOST SO THAT IF
C              SYMMETRY OF ONE STRING IS ISM , THE SYMMETRY OF THE
C              OTHER STRING IS ISMOST(ISM),
C              WITH A ZERO INDICATING A FORBIDDEN COMBINATION
C
#include <implicit.h>
C
      DIMENSION ISMOST(*)
C
      CALL ISETVC(ISMOST,0,NSMST)
      DO 200 ISMST = 1, NSMST
        ISMOST(ISMST) = 1+IEOR(ITOTSM-1,ISMST-1)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck densit */
      SUBROUTINE DENSIT(CL,CR,LUL,LUR,ILSMOS,IRSMOS,
     &                  ILOCOC,IROCOC,
     &                  ISSOA,NSSOA,ISSOB,NSSOB,NOCTPA,NOCTPB,
     &                  DISSYM,SXSTSM,NSMST,NSMSX,XNDXCI , NORB ,
     &                  RHO1,RHO2,C2,
     &                  ICREA,IANNI,IDMTYP,PERMSM,
     &                  ITPFOB,ISPIN1,ISPIN2,ILVERV,ILTSOB,ISTLOB,
     &                  RHOW,L,R,SIGNAR,SGNAR2,IIJSX,ITOSX,VEC1,SCR1,
     &                  NSXPTP,KOC,NFTPK,RHO1A,RHO1B,RHO2A,RHO2B)
#include <implicit.h>
      REAL * 8 INPROD
      LOGICAL PERMSM
      EXTERNAL BDJEP
#include <tstjep.h>
#include <priunit.h>
c.Input vectors
      DIMENSION CL(*),CR(*)
c.Input string information
      DIMENSION ISSOA(NOCTPA,*),NSSOA(NOCTPA,*)
      DIMENSION ISSOB(NOCTPB,*),NSSOB(NOCTPB,*)
      DIMENSION ILOCOC(NOCTPA,NOCTPB),IROCOC(NOCTPA,NOCTPB)
      INTEGER   DISSYM(*),SXSTSM(NSMSX,NSMST)
      INTEGER   ICREA(*),IANNI(*),ITPFOB(*)
      INTEGER   ILSMOS(*),IRSMOS(*)
      INTEGER   ILTSOB(*),ISTLOB(*)
c
c. Detailed string information for lower level routines
      DIMENSION XNDXCI(*)
c. Output
      DIMENSION RHO1(*),RHO2(*)
#include <mxpdim.h>
#include <maxash.h>
c. Scratch blocks
      DIMENSION ILOFF(MXPTP,MXPTP),IROFF(MXPTP,MXPTP)
      DIMENSION ILOFF2(MXPTP,MXPTP),IROFF2(MXPTP,MXPTP)
      DIMENSION NLDT(MXPTP,MXPTP),NRDT(MXPTP,MXPTP)
      DIMENSION NLFTP(MXPTP),NRFTP(MXPTP)
      DIMENSION ILTFRT(MXPTP),IBRFTP(MXPTP)
      DIMENSION RHOW(*)
c. Largest number of strings of a given symmetry
      INTEGER   L(MXPST),R(MXPST)
      DIMENSION SIGNAR(MXPST),SGNAR2(MXPST)
      INTEGER   IIJSX(MXPST),ITOSX(MXPST)
      DIMENSION VEC1(MXPST),SCR1(MXPST)
      INTEGER   NSXPTP(MXPTP)
      INTEGER   KOC(2*MXPST )
      INTEGER   NFTPK(MXPTP)
c
      DIMENSION RHO1A(MAXASH ** 2),RHO1B(MAXASH ** 2)
      DIMENSION RHO2A(MAXASH ** 2),RHO2B(MAXASH ** 2)
      DIMENSION C2(*)
c
c
c This is the long awaited turbo density routine for general RAS
c Prepared by Jeppe in the spring of 1990
c ( Density matrices for the nineties )
c ======
c Input
c ======
c
c. Information about strings
c
c ISSOA(IATP,IASM) : Offset for a-strings of type IATP and symmetry IASM
c NSSOA(IATP,IASM) : Number of  a-strings of type IATP and symmetry IASM
c ISSOB(IBTP,IBSM) : Offset for b-strings of type IBTP and symmetry IBSM
c NSSOB(IBTP,IBSM) : Number of  b-strings of type IBTP and symmetry IBSM
c NOCTPA : Number of types of a-strings
c NOCTPB : Number of types of b-strings
c DISSYM : Symmetry of single excitations
c SXSTSM : String symmetry of SX symmetry times ST symmetry
c NSMST  : Number of symmetries of strings
c NSMSX : Number of symmetries of single excitations
c
c. Information about CL and CR
c
c CL : Left side vector, or buffer for it
c CR : Right side vector, or buffer for it
c LUL : Logical unit containing CL, if LUL .gt. 0
c       CL is assumed to be stored on LUL and is read in in symmetry
c       blocks
c LUR : Logical unit containing CR, if LUR .gt. 0
c       CR is assumed to be stored on LUR and is read in in symmetry
c       blocks
c ILSMOS(ISM) : If symmetry of one string in CL is of symmetry ISM,
c               the other string has symmetry ILSMOS(ISM)
c IRSMOS(ISM) : If symmetry of one string in CR is of symmetry ISM,
c               the other string has symmetry IRSMOS(ISM)
c ILOCOC(IATP,IBTP)  =  1 : this combination is allowed in CL
c IROCOC(IATP,IBTP) .ne.1 : this combination is allowed in CR
c
c XNDXCI : Array containing detailed information about strings
c          This is all this routine is allowed to know |
c NORB   : Number of orbitals
c PERMSM : Logical  .TRUE. => use alpha-beta interchange to reduce
c                             operation count
c========
c.Output
c========
c RHO1 : One-body density matrix ordered as RHO1(I,J) = <CL|E(ij)|CR>
c RHO2  : Two-body density matrix
c
c The exact form of RHO2 is no concern of this routine,
c information to obtain all types of two-body
c density matrices are produced in this routine,
c the external routine DISPDM defines
c mode of storage as well as symmetrization.
c
c========
c Flags
c========
c
c IDMTYP : -1  only rho1, spin-density matrix
c           1  only rho1, charge-density matrix
c           2  one- and two-particle density matrices
c
c
c ILVERV = 1 => Left and Right are described by the same vector
c               Take care when transposing vectors
      SIGNLL = (-1.0D0) ** (ISPIN1+ISPIN2)
c
      NORBSQ = NORB * NORB
      NTEST = 0
      IF(NTEST.NE.0) THEN
         WRITE(LUPRI,*) ' Welcome to the world of DENSIT '
         WRITE(LUPRI,*) ' ILVERV ',ILVERV
      END IF
      IF(ITSTLP.EQ.1) THEN
         WRITE(LUPRI,*)
     &'  DENSIT : Notice only alpha-alpha+beta-beta loops'
      END IF
      IF(ITSTLP.EQ.2) THEN
         WRITE(LUPRI,*)
     &'  DENSIT : Notice only alpha-beta loops'
      END IF
c
c. Loop over symmetry blocks of CL(IASM,IBSM) and CR(JASM,JBSM)
c

      IF(LUL.GT.0) CALL REWINO (LUL)
      ICLSOF = 1
      DO 7000 IASM = 1, NSMST
        IBSM = ILSMOS(IASM)
c. first element in symmetry block
        IF(LUL.GT.0.OR.IASM.EQ.1) THEN
          ILSOFF = 1
        ELSE
          ILSOFF = ILSOFF + NLDTBL
        END IF
        IF(IBSM.EQ.0) GOTO 7000
c. Number of determinants in this block and offset for occupation-
c. occupation class
        CALL SMBLK(NOCTPA,NOCTPB,ILOCOC,NSSOA(1,IASM),
     &       NSSOB(1,IBSM),NLDTBL,ILOFF,NLDT,MXPTP,IMXATP,IMXBTP)
c
        IF(NLDTBL.EQ.0) GOTO 7000
        IF(LUR.GT.0) CALL REWINO(LUR)
        DO 6000 JASM = 1, NSMST
          JBSM = IRSMOS(JASM)
c. first element in symmetry block
          IF(LUR.GT.0.OR.JASM.EQ.1) THEN
           IRSOFF = 1
          ELSE
           IRSOFF = IRSOFF + NRDTBL
          END IF
        IF(JBSM.EQ.0) GOTO 6000
c. Number of determinants in this block and offset for occupation-
c. occupation class
          CALL SMBLK(NOCTPA,NOCTPB,IROCOC,NSSOA(1,JASM),
     &       NSSOB(1,JBSM),NRDTBL,IROFF,NRDT,MXPTP,JMXATP,JMXBTP)
        IF(NRDTBL.EQ.0) GOTO 6000
          IF (NTEST .GE. 3) THEN
            WRITE(LUPRI,*) 'DENSIT: IASM, JASM, IBSM, JBSM',
     &         IASM,JASM,IBSM,JBSM
            WRITE(LUPRI,*) 'DENSIT: NLDTBL, NRDTBL',NLDTBL,NRDTBL
          END IF
c
          IF(IBSM.EQ.JBSM.AND.ITSTLP.NE.2) THEN
c****************************************************************
c                       Part 1                                  *
c                   =============                               *
c Contributions CL(IA,IB)<IA|EA(ij)EA(kl)|JA>CR(JA,IB) to RHO 2 *
c               CL(IA,IB)<IA|EA(kl)|JA>CR(JA,IB)       to RHO 1 *
c                                                               *
c****************************************************************
c
CT             CALL QENTER('DENS1',0)
               IF(.NOT.PERMSM) THEN
               CALL SETVEC(RHO2B,0.0D0,NORBSQ)
c. Transpose symmetry blocks of CL and CR to obtain
c. improved vectorization
               CALL TRRSS2(CL(ILSOFF),C2,NOCTPA,NOCTPB,ILOCOC,NSSOA,
     &              NSSOB,IASM,IBSM,ILOFF2,ILOFF,1,MXPTP)
               CALL COPVEC(C2,CL(ILSOFF),NLDTBL)
               IF(ILVERV.EQ.1.AND.IASM.EQ.JASM) THEN
                 CALL ICOPVE(ILOFF2,IROFF2,MXPTP**2)
               ELSE
                 CALL TRRSS2(CR(IRSOFF),C2,NOCTPA,NOCTPB,IROCOC,NSSOA,
     &                 NSSOB,JASM,JBSM,IROFF2,IROFF,1,MXPTP)
                 CALL COPVEC(C2,CR(IRSOFF),NRDTBL)
               END IF
c ========================================
c. Contributions to one-body density matrix
c ========================================
               CALL SETVEC(RHO1A,0.0D0,NORBSQ)
               DO 400 JATP = 1, JMXATP
               DO 390 JJASTR = 1, NSSOA(JATP,JASM)
c.All excitations out from JJSTR
                 CALL GTSXST(JJASTR,JATP,JASM,IASM,0,NSX,IIJSX,ITOSX,
     &                SGNAR2,NSXPTP,XNDXCI,1,NSMST,0,NSMSX,SXSTSM,
     &                IMXATP,0)
                 DO 380 IATP = 1, IMXATP
                   IF(IATP.EQ.1) THEN
                     IJEXOF = 1
                   ELSE
                     IJEXOF = IJEXOF + NSXPTP(IATP-1)
                   END IF
                   DO 370 IBTP = 1, IMXBTP
                     NIB = NSSOB(IBTP,IBSM)
                     IIROFF = IRSOFF - 1 + IROFF2(JATP,IBTP)
     &                     + (JJASTR-1)*NIB
                     IF(ILOCOC(IATP,IBTP).EQ.1.AND.
     &                  IROCOC(JATP,IBTP).EQ.1) THEN
                       DO 360 IIAST = 1,NSXPTP(IATP)
                         IASTR =  ITOSX( IJEXOF-1+IIAST)
                         SIGNIJ = SGNAR2(IJEXOF-1+IIAST)
                         IJ     = IIJSX( IJEXOF-1+IIAST )
                         IILOFF = ILSOFF - 1 + ILOFF2(IATP,IBTP)
     &                          + (IASTR-1)*NIB
                         RHO1A(IJ) = RHO1A(IJ)+
     &                   SIGNIJ*INPROD(CL(IILOFF),CR(IIROFF),NIB)
  360                  CONTINUE
                     END IF
  370              CONTINUE
  380            CONTINUE
  390          CONTINUE
  400        CONTINUE
             CALL VECSUM(RHO1,RHO1,RHO1A,1.0D0,1.0D0,NORBSQ)
             IF( NTEST .GE. 3  ) THEN
               WRITE(LUPRI,*) ' RHO1 after alpha-alpha  part '
               CALL WRTMAT(RHO1,NORB,NORB,NORB,NORB,0)
             END IF
c===================================
c Two-particle density contributions
c===================================
             IRB = 0
             IF(IDMTYP.EQ.2) THEN
             CALL SETVEC(RHO2B,0.0D0,NORBSQ)
             CALL MTRSP(NORB,NORB,RHO1A,NORB,RHOW,NORB)
             CALL COPVEC(RHOW,RHO1A,NORB ** 2 )
             DO 4110 LL = 1, NORB
             IRA = 0
c.Strings is sym JASM, that has orbitals L occupied
C                   OINSTD(KORB,ISM,NFTP,NSTRIN,ISTRIN,IAB,XNDXCI)
             CALL OINSTD(LL,JASM,NFTPK,NSTK,KOC,1,XNDXCI)
             IF (NSTK .GT. 2*MXPST) THEN
                WRITE (LUPRI,*) 'DENSIT.OINSTD.1 error, NSTK .gt. '//
     &               '2*MXPST'
                WRITE (LUPRI,*) 'NSTK, MXPST:',NSTK,MXPST
                CALL QUIT('DENSIT.OINSTD.1: MXPST too small')
             END IF
             DO 4100 K = 1, NORB
               KL = (K-1)*NORB + LL
               CALL SETVEC(RHO2A,0.0D0,NORBSQ)
               KLSM = DISSYM(KL)
c. Obtain all excitations |L(I> = Sign(I)E(kl)|R(I)>
               ISEL = 1
               CALL LRFEXD(KL,0,JASM,L,R,SIGNAR,NLFTP,NRFTP,ILTFRT,
     &              NTERM,IBRFTP,1,SXSTSM,XNDXCI,NSMSX,ISEL,NFTPK,KOC)
               IF(NTERM.EQ.0) GOTO 4090
               IRA = 1
               DO 800 JATP = 1, JMXATP
               DO 700 JJASTR = 1, NRFTP(JATP)
                 JASTR = R(IBRFTP(JATP)-1+JJASTR)
                 KASTR = L(IBRFTP(JATP)-1+JJASTR)
                 SIGNKL = SIGNAR(IBRFTP(JATP)-1+JJASTR)
                 KASM = SXSTSM(KLSM,JASM)
                 KATP = ILTFRT(JATP)
c. Excitations EA(ij) from KASTR
                 CALL GTSXST(KASTR,KATP,KASM,IASM,0,NSX,IIJSX,ITOSX,
     &                SGNAR2,NSXPTP,XNDXCI,1,NSMST,KL,NSMSX,SXSTSM,
     &                NOCTPA,0)
c. Sign array for double excitation
                 IF(SIGNKL.EQ.-1.0D0) CALL SCALVE(SGNAR2,-1.0D0,NSX)
                 DO 650 IATP = 1, IMXATP
                   NIACT = NSXPTP(IATP)
                   IF(IATP.EQ.1) THEN
                     IJOFF = 1
                   ELSE
                     IJOFF = IJOFF + NSXPTP(IATP-1)
                   END IF
                   IF(NIACT.EQ.0) GOTO 650
                   DO 640 IBTP = 1, IMXBTP
                 IF(ILOCOC(IATP,IBTP).EQ.1.AND.IROCOC(JATP,IBTP).EQ.1)
     &             THEN
                     ILSSO = ILOFF2(IATP,IBTP)
                     IRSSO = IROFF2(JATP,IBTP)
                     NBSSO = NSSOB(IBTP,IBSM)
C                        MINPRD(VU,A,VI,IP,NPROD,NROW)
                     CALL MINPRD(SCR1(IJOFF),CL(ILSSO+ILSOFF-1),
     &                    CR(IRSSO+IRSOFF-1+(JASTR-1)*NBSSO),
     &                    ITOSX(IJOFF),NIACT,NBSSO)
                     DO 630 IJEFF = IJOFF,IJOFF+NIACT-1
                       RHO2A(IIJSX(IJEFF)) =
     &                 RHO2A(IIJSX(IJEFF)) + SGNAR2(IJEFF)*SCR1(IJEFF)
  630                CONTINUE
                   END IF
  640              CONTINUE
  650            CONTINUE
  700          CONTINUE
  800          CONTINUE
 4090          CONTINUE
               IF(NTEST.GE.10) THEN
               WRITE(LUPRI,*) ' RHO2A after loop 1 for ',ICREA(KL),
     &                 IANNI(KL)
               CALL WRTMAT(RHO2A,NORB,NORB,NORB,NORB,0)
               END IF

               CALL DISPDM(RHO2A,RHO2B,RHO1A,RHO2,NORB,ICREA(KL),
     &              IANNI(KL),RHOW,ILTSOB,ISTLOB,ISPIN1,ISPIN2,
     &              IRA,IRB)
 4100      CONTINUE
 4110      CONTINUE
           END IF
c. Clean up by transposing relevant blocks once again
           CALL COPVEC(CL(ILSOFF),C2,NLDTBL)
           CALL TRRSS2(CL(ILSOFF),C2,NOCTPA,NOCTPB,ILOCOC,NSSOA,
     &          NSSOB,IASM,IBSM,ILOFF2,ILOFF,2,MXPTP)
           IF(.NOT.(ILVERV.EQ.1.AND.IASM.EQ.JASM)) THEN
             CALL COPVEC(CR(IRSOFF),C2,NRDTBL)
             CALL TRRSS2(CR(IRSOFF),C2,NOCTPA,NOCTPB,IROCOC,NSSOA,
     &            NSSOB,JASM,JBSM,IROFF2,IROFF,2,MXPTP)
           END IF
       END IF
CT             CALL QEXIT('DENS1',0)
       END IF
c. End of Loop 1
c
             IF(IASM.EQ.JASM.AND.ITSTLP.NE.2) THEN
CT           CALL QENTER('DENS2',0)
c****************************************************************
c                       Part 2                                  *
c                   =============                               *
c Contributions CL(IA,IB)<IB|EB(ij)EB(kl)|JB>CR(IA,JB) to RHO 2 *
c               CL(IA,IB)<IB|EB(kl)|JB>CR(IA,JB)       to RHO 1 *
c                                                               *
c****************************************************************
c
c ========================================
c. Contributions to one-body density matrix
c ========================================
               CALL SETVEC(RHO1B,0.0D0,NORBSQ)
               DO 1400 JBTP = 1, JMXBTP
               DO 1390 JJBSTR = 1, NSSOB(JBTP,JBSM)
c.All excitations out from JJBSTR
                 CALL GTSXST(JJBSTR,JBTP,JBSM,IBSM,0,NSX,IIJSX,ITOSX,
     &                SGNAR2,NSXPTP,XNDXCI,2,NSMST,0,NSMSX,SXSTSM,
     &                NOCTPB,0)
                 DO 1380 IBTP = 1, IMXBTP
                   IF(IBTP.EQ.1) THEN
                     IJEXOF = 1
                   ELSE
                     IJEXOF = IJEXOF + NSXPTP(IBTP-1)
                   END IF
                   DO 1370 IATP = 1, IMXATP
                     NIA = NSSOA(IATP,IASM)
                     IIROFF = IRSOFF - 1 + IROFF(IATP,JBTP)
     &                     + (JJBSTR-1)*NIA
                     IF(ILOCOC(IATP,IBTP).EQ.1.AND.
     &                  IROCOC(IATP,JBTP).EQ.1) THEN
                       DO 1360 IIBST = 1,NSXPTP(IBTP)
                         IBSTR =  ITOSX( IJEXOF-1+IIBST)
                         SIGNIJ = SGNAR2(IJEXOF-1+IIBST)*SIGNLL
                         IJ     = IIJSX( IJEXOF-1+IIBST )
                         IILOFF = ILSOFF - 1 + ILOFF(IATP,IBTP)
     &                          + (IBSTR-1)*NIA
                         RHO1B(IJ) = RHO1B(IJ)+
     &                   SIGNIJ*INPROD(CL(IILOFF),CR(IIROFF),NIA)
 1360                  CONTINUE
                     END IF
 1370              CONTINUE
 1380            CONTINUE
 1390          CONTINUE
 1400        CONTINUE
             IF(PERMSM) CALL SCALVE(RHO1B,2.0D0,NORBSQ)
             IF(IDMTYP.NE.-1) THEN
               CALL VECSUM(RHO1,RHO1,RHO1B,1.0D0,1.0D0,NORBSQ)
             ELSE IF( IDMTYP.EQ.-1) THEN
               IF (.NOT.PERMSM) THEN
                 CALL VECSUM(RHO1,RHO1,RHO1B,0.5D0,-0.5D0,NORBSQ)
               END IF
             END IF
c
             IF( NTEST .GE. 3  ) THEN
               WRITE(LUPRI,*) ' RHO1 after beta-beta  part '
               CALL WRTMAT(RHO1,NORB,NORB,NORB,NORB,0)
             END IF
c ========================================
c. Contributions to two-body density matrix
c ========================================
             IRB = 0
c. Obtain all excitations |L(I> = Sign(I)E(kl)|R(I)>
             IF(IDMTYP.EQ.2) THEN
             CALL SETVEC(RHO2B,0.0D0,NORBSQ)
             CALL MTRSP(NORB,NORB,RHO1B,NORB,RHOW,NORB)
             CALL COPVEC(RHOW,RHO1B,NORB ** 2 )
             DO 4210 LL = 1, NORB
c.Strings is sym JBSM, that has orbitals L occupied
CT      CALL QENTER('DOIN',0)
             CALL OINSTD(LL,JBSM,NFTPK,NSTK,KOC,2,XNDXCI)
             IF (NSTK .GT. 2*MXPST) THEN
                WRITE (LUPRI,*) 'DENSIT.OINSTD.2 error, NSTK .gt. '//
     &               '2*MXPST'
                WRITE (LUPRI,*) 'NSTK, MXPST:',NSTK,MXPST
                CALL QUIT('DENSIT.OINSTD.2: MXPST too small')
             END IF
CT      CALL QEXIT('DOIN',0)
             DO 4200 K  = 1, NORB
               IRA = 0
               KL = (K-1)*NORB + LL
c*Start Correction Jan 31, rather early in the morning
COLD           IF(ILVERV.NE.0) THEN
COLD             LK = (LL-1)*NORB + K
COLD             MAXKL = MAX(KL,LK)
COLD           ELSE
COLD             MAXKL = KL
COLD           END IF
               MAXKL = KL
c*End Correction Jan 31, rather early in the morning
               CALL SETVEC(RHO2A,0.0D0,NORBSQ)
               KLSM = DISSYM(KL)
c. Obtain all excitations |L(I> = Sign(I)E(kl)|R(I)>
               ISEL = 1
CT       CALL QENTER('DLRFE',0)
               CALL LRFEXD(KL,0,JBSM,L,R,SIGNAR,NLFTP,NRFTP,ILTFRT,
     &              NTERM,IBRFTP,2,SXSTSM,XNDXCI,NSMSX,ISEL,NFTPK,KOC)
CT       CALL QEXIT('DLRFE',0)
               IF(NTERM.EQ.0) GOTO 4190
               IRA = 1
               DO 1800 JBTP = 1, JMXBTP
c.Number of dets in this block
               NJAJB = 0
               DO 1799 JATP = 1, JMXATP
                 NJAJB = NJAJB + NRDT(JATP,JBTP)
 1799          CONTINUE
               IF(NJAJB.EQ.0) GOTO 1800
c
               DO 1700 JJBSTR = 1, NRFTP(JBTP)
                 JBSTR = R(IBRFTP(JBTP)-1+JJBSTR)
                 KBSTR = L(IBRFTP(JBTP)-1+JJBSTR)
                 SIGNKL = SIGNAR(IBRFTP(JBTP)-1+JJBSTR)*SIGNLL
                 KBSM = SXSTSM(KLSM,JBSM)
                 KBTP = ILTFRT(JBTP)
c. Excitations EB(ij) from KBSTR
                 IGTJ = 0
                 CALL GTSXST(KBSTR,KBTP,KBSM,IBSM,0,NSX,IIJSX,ITOSX,
     &                SGNAR2,NSXPTP,XNDXCI,2,NSMST,MAXKL,NSMSX,SXSTSM,
     &                IMXBTP,IGTJ)
c. Sign array for double excitation
                 IF(SIGNKL.EQ.-1.0D0) CALL SCALVE(SGNAR2,-1.0D0,NSX)
                 DO 1650 IBTP = 1, IMXBTP
                   NIACT = NSXPTP(IBTP)
                   IF(IBTP.EQ.1) THEN
                     IJOFF = 1
                   ELSE
                     IJOFF = IJOFF + NSXPTP(IBTP-1)
                   END IF
                   IF(NIACT.EQ.0) GOTO 1650
                   DO 1640 IATP = 1, IMXATP
                 IF(ILOCOC(IATP,IBTP).EQ.1.AND.IROCOC(IATP,JBTP).EQ.1)
     &             THEN
                     NASSO = NSSOA(IATP,IASM)
                     IF(NASSO.EQ.0) GOTO 1640
                     ILSSO = ILOFF(IATP,IBTP)
                     IRSSO = IROFF(IATP,JBTP)
                     CALL MINPRD(SCR1(IJOFF),CL(ILSSO+ILSOFF-1),
     &                    CR(IRSSO+IRSOFF-1+(JBSTR-1)*NASSO),
     &                    ITOSX(IJOFF),NIACT,NASSO)
                     DO 1630 IJEFF = IJOFF,IJOFF+NIACT-1
                       RHO2A(IIJSX(IJEFF)) =
     &                 RHO2A(IIJSX(IJEFF)) + SGNAR2(IJEFF)*SCR1(IJEFF)
 1630                CONTINUE
                   END IF
 1640              CONTINUE
 1650            CONTINUE
 1700          CONTINUE
 1800          CONTINUE
               IF(PERMSM) CALL SCALVE(RHO2A,2.0D0,NORBSQ)
 4190          CONTINUE
c
               IF(NTEST.GE.10) THEN
               WRITE(LUPRI,*) ' RHO2A after loop 2 for ',ICREA(KL),
     &                 IANNI(KL)
               CALL WRTMAT(RHO2A,NORB,NORB,NORB,NORB,0)
               END IF
c
               CALL DISPDM(RHO2A,RHO2B,RHO1B,RHO2,NORB,ICREA(KL),
     &                     IANNI(KL),RHOW,ILTSOB,ISTLOB,ISPIN1,ISPIN2,
     &                     IRA,IRB)

 4200        CONTINUE
 4210        CONTINUE
             END IF
CT           CALL QEXIT('DENS2',0)
               END IF
c. End of Loop 2
C1298    continue
c
c********************************************************************
c                                Part 3                             *
c                             =============                         *
c Contributions CL(IA,IB)<IB|EB(ij)|JB>*<IA|EA(IJ)|JA>*CR(JA,JB>    *
c                                                                   *
c********************************************************************
c
c. Obtain all excitations |L(I> = Sign(I)E(kl)|R(I)>
            IF(IDMTYP.EQ.2.AND.ITSTLP.NE.1) THEN
            IRA = 0
            IRB = 1
CT          CALL QENTER('DENS3',0)
            CALL SETVEC(RHO2A,0.0D0,NORBSQ)
            CALL SETVEC(RHO1A,0.0D0,NORBSQ)
            DO 4310 LL = 1,NORB
c.All strings  of symmetry JASM that has LL occupied
             CALL OINSTD(LL,JASM,NFTPK,NSTK,KOC,1,XNDXCI)
             IF (NSTK .GT. 2*MXPST) THEN
                WRITE (LUPRI,*) 'DENSIT.OINSTD.3 error, NSTK .gt. '//
     &               '2*MXPST'
                WRITE (LUPRI,*) 'NSTK, MXPST:',NSTK,MXPST
                CALL QUIT('DENSIT.OINSTD.3: MXPST too small')
             END IF
            DO 4300 K = 1, NORB
            KL = (K-1)*NORB + LL
             KLSM = DISSYM(KL)
             IF(SXSTSM(KLSM,JASM).NE.IASM) GOTO 4300
             CALL SETVEC(RHO2B,0.0D0,NORBSQ)
             ISEL = 1
             CALL LRFEXD(KL,0,JASM,L,R,SIGNAR,NLFTP,NRFTP,ILTFRT,
     &            NTERM,IBRFTP,1,SXSTSM,XNDXCI,NSMSX,ISEL,NFTPK,KOC)
c.Gather C2(I,J) = CL(L(I),J)*SIGNAR(I)
             CALL RASCG6(CL(ILSOFF),C2,NOCTPA,NSSOA(1,IASM),NLFTP,
     &                  NOCTPB,NSSOB(1,IBSM),ILOFF2,L,SIGNAR,
     &                  ILOCOC,MXGCTP,0,ILTFRT,ILOFF,MXPTP)
c( Reordering eliminated
c. The term can now be written as
c. CL(I,IB)<IB|EB(IJ)|JB> C2((R(I),JB)
             DO 2800 JBTP = 1, JMXBTP
               IF(PERMSM) THEN
                 JATPMN =  JBTP
               ELSE
                 JATPMN = 1
               END IF
             DO 2700 JJBSTR = 1, NSSOB(JBTP,JBSM)
               JBSTR = JJBSTR + ISSOB(JBTP,JBSM) - 1
c. Gather VEC1(I) = CR(R(I),JB)
               JOFF = 1
c
               DO 2650 JATP = JATPMN, JMXATP
               IROFF2(JATP,1) = JOFF
               IF(IROCOC(JATP,JBTP).EQ.1) THEN
                 DO 2640 JA = 1, NRFTP(JATP)
                   VEC1(JOFF) =
     &             CR(IRSOFF-1+IROFF(JATP,JBTP)-1+R(IBRFTP(JATP)-1+JA)
     &                +(JJBSTR-1)*NSSOA(JATP,JASM))
                   JOFF = JOFF + 1
 2640            CONTINUE
               END IF
 2650          CONTINUE
               IF(JOFF.EQ.1) GOTO 2800
c.The term is now CL(I,IB)<IB|EB(ij)|JB>CR(I,JB)
c. Excitations EB(ij) from JBSTR
               CALL GTSXST(JJBSTR,JBTP,JBSM,IBSM,0,NSX,IIJSX,ITOSX,
     &              SIGNAR,NSXPTP,XNDXCI,2,NSMST,0,NSMSX,SXSTSM,MXGCTP,
     &              0)
               DO 2550 IBTP = 1, MXGCTP
                 NIACT = NSXPTP(IBTP)
                 IF(IBTP.EQ.1) THEN
                   IJOFF = 1
                 ELSE
                   IJOFF = IJOFF + NSXPTP(IBTP-1)
                 END IF
                 IF(NIACT.EQ.0) GOTO 2550
                 DO 2540 JATP = JATPMN, JMXATP

                 IF(NRFTP(JATP).EQ.0) GOTO 2540
                 IATP = ILTFRT(JATP)
               IF(ILOCOC(IATP,IBTP).EQ.1.AND.IROCOC(JATP,JBTP).EQ.1)
     &           THEN
                   ILSSO = ILOFF2(IATP,IBTP)
                   IRSSO = IROFF2(JATP,1)
                   NASSO = NRFTP(JATP)
                   CALL MINPRD(SCR1(IJOFF),C2(ILSSO),
     &                  VEC1(IRSSO),
     &                  ITOSX(IJOFF),NIACT,NASSO)
                   IF(PERMSM.AND.JATP.NE.JBTP)
     &             CALL SCALVE(SCR1(IJOFF),2.0D0,NIACT)

                   DO 2530 IJEFF = IJOFF,IJOFF+NIACT-1
                     RHO2B(IIJSX(IJEFF)) =
     &               RHO2B(IIJSX(IJEFF)) + SIGNAR(IJEFF)*SCR1(IJEFF)
 2530              CONTINUE
                   END IF
 2540              CONTINUE
 2550            CONTINUE
 2700          CONTINUE
 2800          CONTINUE
c
               IF(NTEST.GE.10) THEN
               WRITE(LUPRI,*) ' RHO2B after loop 3 for ',ICREA(KL),
     &                 IANNI(KL)
               CALL WRTMAT(RHO2B,NORB,NORB,NORB,NORB,0)
               END IF
c
               CALL DISPDM(RHO2A,RHO2B,RHO1A,RHO2,NORB,ICREA(KL),
     &                     IANNI(KL),RHOW,ILTSOB,ISTLOB,ISPIN1,ISPIN2,
     &                     IRA,IRB)
C              CALL DISPS3(RHO2A,RHO2B,RHO1A,RHO2,RHO2AA,IDOR2A,NORB,
C    &                  ICREA(KL),IANNI(KL))
 4300 CONTINUE
 4310 CONTINUE
c. End of Loop 3
CT    CALL QEXIT('DENS3',0)
      END IF
C5000 CONTINUE
 6000 CONTINUE
c
 7000 CONTINUE
C8000 CONTINUE
c
c. Transpose and reorder Rho1 to SIRIUS order
      CALL MTRSP(NORB,NORB,RHO1,NORB,RHOW,NORB)
      CALL COPVEC(RHOW,RHO1,NORB ** 2 )
      CALL REORMT(RHOW,RHO1,NORB,NORB,ISTLOB,ISTLOB)
c
      IF(NTEST .GE.2 ) THEN
        TRACE = 0.0D0
        DO 0803 II = 1, NORB
          TRACE = TRACE + RHO1((II-1)*NORB + II )
 0803   CONTINUE
        WRITE(LUPRI,*) ' Trace of RHO1 ', TRACE
        WRITE(LUPRI,*) ' Final one-body density matrix '
        CALL WRTMAT(RHO1,NORB,NORB,NORB,NORB,0)
      END IF
c
      RETURN
      END
C  /* Deck phpvec */
      SUBROUTINE PHPVEC(NVEC,I1VEC,CVECS,NCONF,DIAG)
C
C 30-October-1990 Hans Joergen Aa. Jensen
C
C Retrieve NVEC PHP solution vectors in CVECS, starting
C with number I1VEC.
C
#include <implicit.h>
      DIMENSION CVECS(NCONF,NVEC), DIAG(LPHP)
C
C Used from common blocks :
C   PHPINF : LPHP,KPEVEC,NPCVAR,KPNTR
C
#include <phpinf.h>
C
      IF (NVEC .GT. 0) THEN
         CALL DZERO(CVECS,NVEC*NCONF)
         DO 800 IVEC = 1,NVEC
            IEVEC = I1VEC - 1 + IVEC
            JPEVEC = KPEVEC + (IEVEC-1)*NPCVAR
            CALL SCAVEC(CVECS(1,IVEC),DIAG(JPEVEC),DIAG(KPNTR),NPCVAR)
  800    CONTINUE
      END IF
      RETURN
      END
C  /* Deck bdjep */
      BLOCK DATA BDJEP
#include <implicit.h>
#include <tstjep.h>
      DATA ITSTLP/0/
      END

