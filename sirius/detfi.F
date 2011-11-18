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
C  /* Deck strpa2 */
      SUBROUTINE STRPA2(NSASOA,NSASOB,NASTR,NBSTR,
     &                  MAXSYM,NOCTPA,NOCTPB,IOCOC,
     &                  NL1MNA,NL1MNB,NEL1MN,NEL1MX,
     &                  NEL3MN,NEL3MX,NORB1,NORB2,NORB3,
     &                  IOCPTA,IOCPTB,SYMPRO,NAEL,NBEL,NEL,
     &                  REFSPC,NDTASM,MXNDT,MXVBLK,NTEST)
C
C CALCULATE NUMBER OF DETERMINANTS CONSISTENT WITH RESTRICTIONS
C ON OCCUPATION SPACE AND SYMMETRY
C
C      Mar  7 89 hjaaj: corrected loop limit for IAEL1 and IBEL1
C                       introduced REFSPC option
C      SEP 14 88 NEL1MX AND NEL3MN ADDED
C L.R. JAN 14 88 ( NDTASM AND MXNDT ADDED )
C
#include <implicit.h>
#include <priunit.h>
      INTEGER SYMPRO(8,8)
      DIMENSION NSASOA(NOCTPA,MAXSYM),NSASOB(NOCTPB,MAXSYM)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION IOCPTA(*),IOCPTB(*)
      DIMENSION NDTASM(MAXSYM)
      LOGICAL   REFSPC
C
C*** 1 : SET UP ARRAY IOCOC THAT SPECIFIES WHICH COMBINATIONS
C        OF ALPHA AND BETA TYPES ARE ALLOWED .
C
      CALL ISETVC(IOCOC,0,NOCTPA*NOCTPB)
C
C** LOOP OVER ALLOWED ALPHA TYPES
C
CHJ-910416 corrected IAEL1 and IBEL1 upper loop limits
Cold  DO 100 IAEL1 = NL1MNA, MIN(NORB1,NEL1MX)
      DO 100 IAEL1 = NL1MNA, MIN(NORB1,NEL1MX+1)
      DO 100 IAEL3 = 0     , NEL3MX
C
        IF( (IAEL1 + IAEL3 .GT. NAEL) .OR.
     &      (NAEL - IAEL1 - IAEL3 .GT. NORB2) ) GOTO 100
        IOCTPA = (  NORB1 - IAEL1 ) * ( NEL3MX + 1 ) + IAEL3 + 1
        IOCTPA = IOCPTA(IOCTPA)
        IF (IOCTPA .EQ. 0) GO TO 100
C** ALLOWED BETA TYPES
        DO 50 IBEL1 = NL1MNB,MIN(NORB1,NEL1MX+1)
        DO 50 IBEL3 = 0     ,NEL3MX
          IF( (IBEL1 + IBEL3 .GT. NBEL) .OR.
     &        (NBEL - IBEL1 - IBEL3 .GT. NORB2) ) GOTO 50
C
          IOCTPB = (  NORB1 - IBEL1 ) * ( NEL3MX + 1 ) + IBEL3 + 1
          IOCTPB = IOCPTB(IOCTPB)
          IF( IOCTPB .EQ. 0 ) GOTO 50
C
          IEL1 = IAEL1+IBEL1
          IEL3 = IAEL3+IBEL3
          IF( IEL1 .GT. 2*NORB1 .OR. IEL3 .GT. 2*NORB3 ) GOTO 50
C
          IF( NEL1MN .LE. IEL1 .AND. IEL3 .LE. NEL3MX ) THEN
             IF( IEL1 .LE. NEL1MX .AND. NEL3MN .LE. IEL3 ) THEN
C                CONGRATULATIONS , YOU ARE ALLOWED
                 IOCOC(IOCTPA,IOCTPB) = 1
             ELSE IF( REFSPC .AND.
     &                IEL1 .EQ. 2*NORB1 .AND. IEL3 .EQ. 0) THEN
C                This is a reference determinant in, e.g., MR-CID
                 IOCOC(IOCTPA,IOCTPB) = 1
             ELSE IF( IEL1 .LE. (NEL1MX+1) .AND.
     &               (NEL3MN-1) .LE. IEL3 ) THEN
C                YOU ARE NOT ALLOWED BUT I HAVE TO ACCEPT YOU AS
C                AN INTERMEDIATE STATE
                 IOCOC(IOCTPA,IOCTPB) = - 1
             END IF
           END IF
   50   CONTINUE
  100 CONTINUE
      IF ( NTEST .GE. 4 ) THEN
        WRITE(LUPRI,'(/A/)') ' Matrix giving allowed combinations '//
     &                       'of types'
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
C
C** TOTAL NUMBER OF DETERMINANTS
C
      IF( NTEST .GE. 4 ) THEN
       WRITE(lupri,'(/A)') ' Number of strings per type '//
     &                     'and sym ( alpha)'
       CALL IWRTMA(NSASOA,NOCTPA,MAXSYM,NOCTPA,MAXSYM)
       WRITE(LUPRI,'(/A)') ' Number of strings per type '//
     &                     'and sym ( beta )'
       CALL IWRTMA(NSASOB,NOCTPB,MAXSYM,NOCTPB,MAXSYM)
      END IF
C
      MXVBLK = 0
      MXNDT  = 0
      DO 300 ITOTSM = 1, MAXSYM
      NDET = 0
      DO 200 IASYM = 1,MAXSYM
        NDETSM = 0
        DO 190 IATYP = 1, NOCTPA
          IASTRI = NSASOA(IATYP,IASYM)
          IBSYM  = SYMPRO(IASYM,ITOTSM)
          DO 150 IBTYP = 1, NOCTPB
            IF(IOCOC(IATYP,IBTYP) .EQ. 1 ) THEN
              NDETSM = NDETSM + IASTRI * NSASOB(IBTYP,IBSYM)
            END IF
  150    CONTINUE
  190  CONTINUE
       NDET   = NDET + NDETSM
       MXVBLK = MAX(MXVBLK,NDETSM)
  200 CONTINUE
      NDTASM(ITOTSM) = NDET
      MXNDT = MAX(NDET,MXNDT)
  300 CONTINUE
C
      IF (NTEST .GE. 2 ) THEN
        WRITE(LUPRI,'(/A)') ' >>> OUTPUT FROM STRPA2 (ras) <<<'
C
        WRITE(LUPRI,'(/A/8I9)') ' Number of determinants per symmetry',
     &  (NDTASM(I),I=1,MAXSYM)
C
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) 'Largest determinant ci-expansion',MXNDT
        WRITE(LUPRI,*) 'mxvblk from strpa2              ',MXVBLK
C
      END IF
C
      RETURN
      END
C  /* Deck ioctyp */
      INTEGER FUNCTION IOCTYP(STRING,IAB,IFLAG)
C
C OBTAIN OCCUPATION TYPE NUMBER FOR STRING WITH
C FORBIDDEN STRINGS RETURN WITH A ZERO TYPE
C IF IFLAG .NE. 0 REDEFINED TYPE IS OBTAINED THROUGH
C OCCUPATION POINTERS.
C IAB = 1 ( 2 ) DESIGNATES ALPHA ( BETA )   STRING
C
#include <implicit.h>
#include <priunit.h>
      INTEGER  STRING(*)
C
#include <mxpdim.h>
#include <ciinfo.h>
#include <strnum.h>
C
C OCCUPATION CLASSS OF STRING
C
          IEL1 = 0
          IEL3 = 0
          IF (  IAB .EQ. 1  ) NEL = NAEL
          IF (  IAB .EQ. 2  ) NEL = NBEL
          DO 80 IEL = 1,NEL
            IF(IORB1F .LE. STRING(IEL) .AND.
     &         STRING(IEL) .LE. IORB1L   ) IEL1 = IEL1 +1
            IF(IORB3F .LE. STRING(IEL) .AND.
     &         STRING(IEL) .LE. IORB3L   ) IEL3 = IEL3 +1
   80     CONTINUE
C TYPE
      ITYP = ( NORB1 - IEL1 ) * ( NEL3MX + 1 ) + IEL3 + 1
C ALLOWED ?
      IF ( IAB .EQ. 1 ) THEN
         IF(( IEL1 + IEL3 .GT. NAEL )  .OR.
     &      ( IEL1 .LT. NL1MNA      )  .OR.
     &      ( IEL3 .GT. NEL3MX     )) ITYP = 0
      ELSE
        IF(( IEL1 + IEL3 .GT. NBEL )  .OR.
     &      ( IEL1 .LT. NL1MNB     )  .OR.
     &      ( IEL3 .GT. NEL3MX     )) ITYP = 0
      END IF
C REDEFINE
      IF ( IFLAG .NE. 0 .AND. ITYP .NE. 0 ) THEN
         IF ( IAB .EQ. 1 ) ITYP = IOCPTA(ITYP)
         IF ( IAB .EQ. 2 ) ITYP = IOCPTB(ITYP)
      END IF
C
#if defined (VAR_TR1TST    )
      NTEST = 0
      IF ( NTEST .GE.10 ) THEN
        WRITE(LUPRI,*) ' FROM IOCTYP : IEL1 IEL3 ITYP ',IEL1,IEL3,ITYP
        WRITE(LUPRI,*) ' STRING',(STRING(IEL),IEL=1,NEL)
      END IF
#endif
C
      IOCTYP = ITYP
C
      RETURN
      END
C  /* Deck numst2 */
      SUBROUTINE NUMST2(NEL,NORB,NORB1,NEL1MN,NORB3,NEL3MX,NSTAOC,
     &                  NSTRIN,NOCCTP,IPOINT,NTEST,NEL1MX,NEL3MN,
     &                  NL1FTP,NL2FTP,NL3FTP)
C
C Last revision 890819-hjaaj: error in IEL1MX and IEL3MX corrected.
C               901010-jo   : NL1FTP,NL2FTP,NL3FTP added
C CALCULATE THE NUMBER OF WAYS NEL ELECTRONS CAN BE DISTRIBUTED IN NORB
C SO THAT THE FOLLOWING CONDITIONS ARE FULFILLED
C
C NUMBER OF ELECTRONS IN THE FIRST NORB1 ORBITALS
C IS BETWEEN NEL1MN AND NEL1MX
C NUMBER OF ELECTRONS IN THE LAST NORB3 ORBITALS
C IS BETWEEN NEL3MN AND NEL3MX
C
C THE NUMBER OF STRINGS IN EACH OCCUPATION CLASS IS ALSO
C CALCULATED. THE INITIAL NUMBER OF A CLASSS IS GIVEN BY
C ( NORB1 - IEL1) *( NEL3MX + 1 ) + IEL3 + 1
C A POINTER , IPOINT, IS SET UP TO ELIMINATE CLASSES
C WITH ZERO STRINGS
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION NSTAOC(*), IPOINT(*)
      DIMENSION NL1FTP(*),NL2FTP(*),NL3FTP(*)
C
      IF (NTEST .GT. 2) THEN
         WRITE(LUPRI,*) ' NEL1MN NEL1MX NEL3MN NEL3MX (from NUMST2)'
         WRITE(LUPRI,*)   NEL1MN,NEL1MX,NEL3MN,NEL3MX
      END IF
      NSTRIN = 0
COLD  NEL1MX = MIN ( NORB1, NEL)
COLD  IEL1MX = MIN ( NORB1, NEL1MX )
CJO OCT90
      IEL1MX = MIN(NORB1,NEL1MX+1)
      NORB2  = NORB - NORB1 - NORB3
      NOCCTP = (NORB1-NEL1MN+1) * (NEL3MX + 1 )
      CALL ISETVC(NSTAOC,0,NOCCTP)
C
      DO 100 IEL1 = NEL1MN,IEL1MX
        NSTIN1 = IBION(NORB1,IEL1)
        IEL3MN = MAX ( NEL3MN,NEL-(IEL1+NORB2) )
        IEL3MX = MIN ( NEL3MX,NEL-IEL1,NORB3 )
        DO 80 IEL3 = IEL3MN, IEL3MX
         IEL2   = NEL - IEL1 - IEL3
         NSTINT = NSTIN1*IBION(NORB2,IEL2)*IBION(NORB3,IEL3)
         IOCCTP = (NORB1-IEL1)*(NEL3MX+1) + IEL3 + 1
         NSTRIN = NSTRIN + NSTINT
         NSTAOC(IOCCTP) = NSTINT
         NL1FTP(IOCCTP) = IEL1
         NL2FTP(IOCCTP) = IEL2
         NL3FTP(IOCCTP) = IEL3
  80   CONTINUE
 100  CONTINUE
C
      IF( NTEST .GE. 1 ) THEN
        WRITE(LUPRI,'(/A)')    ' <<< output from NUMST2 >>>'
        WRITE(LUPRI,'(/A,I6)') ' Number of strings generated ...', 
     &                         NSTRIN
        WRITE(LUPRI,'(/A,I6)') ' Number of occupation classes ..',
     &                         NOCCTP
        WRITE(LUPRI,'(/A)')    ' Number of strings per occupation class'
        CALL IWRTMA(NSTAOC,1,NOCCTP,1,NOCCTP)
      END IF
C
C** REMOVE OCCUPATION TYPES WITH ZERO OCCUPATION.
C   CREATE POINTER ARRAY GIVING NEW NUMBER AD FUNCTION OF OLD NUMBER
C
      IITYPE = 0
      DO 200 ITYPE = 1,NOCCTP
        IF(NSTAOC(ITYPE).NE.0) THEN
          IITYPE = IITYPE + 1
          NSTAOC(IITYPE) = NSTAOC(ITYPE)
          IPOINT( ITYPE) = IITYPE
          NL1FTP(IITYPE) = NL1FTP(ITYPE)
          NL2FTP(IITYPE) = NL2FTP(ITYPE)
          NL3FTP(IITYPE) = NL3FTP(ITYPE)
        END IF
  200 CONTINUE
      NOCCTO = NOCCTP
      NOCCTP = IITYPE
C
C
      IF( NTEST .GE. 2 .AND. NOCCTO .NE. NOCCTP ) THEN
        WRITE(LUPRI,'(/A,I6)') ' Reduced number of types .......',
     &                           NOCCTP
        WRITE(LUPRI,'(/A)')    ' Number of strings per occupation class'
        CALL IWRTMA(NSTAOC,1,NOCCTP,1,NOCCTP)
      END IF
      RETURN
      END
C  /* Deck strdim */
      SUBROUTINE STRDIM(NORB,NEL,MS2,ITPFOB,NTEST)
C
C SOME DIMENSION OF STRINGS FOR DEFINING STATIC PARTITIONING OF
C MEMORY
C
#include <implicit.h>
#include <priunit.h>
#include <mxpdim.h>
#include <maxash.h>
#include <strnum.h>
#include <ciinfo.h>
c.Added sept 90
c NL1FTP : Number of electron in ras 1 for string of given type
c NL2FTP : Number of electron in ras 2 for string of given type
c NL3FTP : Number of electron in ras 3 for string of given type
c IJTYP  : type of a given orbital excitation
C dimension: MAXASH*MAXASH
c IJTST  : Type of string arising from excitation times string type
c          (0 indicates a forbidden excitation )
      COMMON/ADDON/NL1FTP(MXOCTP,2),NL2FTP(MXOCTP,2),NL3FTP(MXOCTP,3),
     &             IJTYP(MAXASH**2),IJTST(9,MXPTP,2)
      DIMENSION ITPFOB(NORB)
C
C
      IF( MS2 .EQ. 0 ) THEN
        EQUAL = .TRUE.
      ELSE
        EQUAL = .FALSE.
      END IF
C
COLD  MS2  = MULTS - 1
      NAEL = (MS2 + NEL ) / 2
      NBEL = (NEL - MS2 ) / 2
      IF (NTEST .GT. 1) WRITE(LUPRI,*) ' MS2 NAEL NBEL ', MS2,NAEL,NBEL
C
      IF (NAEL + NBEL .NE. NEL ) THEN
         WRITE(LUPRI,'(/A/A,2I15)')
     &'STOP in STRDIM: No. of ELECTRONS and MULTIPLICITY INCONSISTENT',
     &'MS2, NEL', MS2, NEL
         CALL QUIT('NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY')
      END IF
C
      NORB1 = IORB1L - IORB1F + 1
      NORB3 = IORB3L - IORB3F + 1
      NORB2 = NORB - NORB1 - NORB3
      NL1MNA = MAX(0,NEL1MN-MIN(NBEL,NORB1) )
      NL1MNB = MAX(0,NEL1MN-MIN(NAEL,NORB1) )
      NL3MNA = MAX(0,NEL3MN-1-MIN(NBEL,NORB3) )
      NL3MNB = MAX(0,NEL3MN-1-MIN(NAEL,NORB3) )
      CALL NUMST2(NAEL,NORB,NORB1,NL1MNA,NORB3,NEL3MX,
     &            NSTAOA,NASTR,NOCTPA,IOCPTA,NTEST,NEL1MX,NL3MNA,
     &              NL1FTP(1,1),NL2FTP(1,1),NL3FTP(1,1) )
      CALL NUMST2(NBEL,NORB,NORB1,NL1MNB,NORB3,NEL3MX,
     &            NSTAOB,NBSTR,NOCTPB,IOCPTB,NTEST,NEL1MX,NL3MNB,
     &              NL1FTP(1,2),NL2FTP(1,2),NL3FTP(1,2) )
C
      CALL ISETVC(ITPFOB,0,NORB)
      IF( NORB1 .NE. 0 ) CALL ISETVC(ITPFOB(IORB1F),1,NORB1)
      IF( NORB2 .NE. 0 ) CALL ISETVC(ITPFOB(IORB1L+1),2,NORB2)
      IF( NORB3 .NE. 0 ) CALL ISETVC(ITPFOB(IORB3F),3,NORB3)
C
C
      IF(NTEST .GE. 10 ) THEN
       WRITE(LUPRI,*)'NAEL,NBEL,NASTR,NBSTR,EQUAL,MAXSYM'
       WRITE(LUPRI,*) NAEL,NBEL,NASTR,NBSTR,EQUAL,MAXSYM
      END IF
C
      CALL ZIJTYP(IJTYP,ITPFOB,NORB)
c.Type of excitation times string
      CALL SXTST(NOCTPA,NL1FTP(1,1),NL2FTP(1,1),NL3FTP(1,1),
     &     IJTST(1,1,1) )
      CALL SXTST(NOCTPB,NL1FTP(1,2),NL2FTP(1,2),NL3FTP(1,2),
     &     IJTST(1,1,2) )
c. Total number of excitations from strings
C     FUNCTION NSEXCI(NTYPE,NSTFTP,NL1FTP,NL2FTP,NL3FTP,
C    &                  NORB1,NORB2,NORB3,IJTST,NTEST)
      NAEXCI = NSEXCI(NOCTPA,NSTAOA,NL1FTP(1,1),NL2FTP(1,1),
     &                NL3FTP(1,1),NORB1,NORB2,NORB3,IJTST(1,1,1),NTEST)
      NBEXCI = NSEXCI(NOCTPB,NSTAOB,NL1FTP(1,2),NL2FTP(1,2),
     &                NL3FTP(1,2),NORB1,NORB2,NORB3,IJTST(1,1,2),NTEST)
c
      RETURN
      END
C  /* Deck numstr */
      SUBROUTINE NUMSTR(NSTRIN,NEL,NORB1,NORB2,NORB3,
     &                  NELMN1,NELMX1,NELMN3,NELMX3,
     &                  NSSO,NSTASM,MXST,
     &                  NOCTYP,MAXSYM,ORBSYM,IOC,IAB,NTEST)
c
c 15-Apr-1991 hjaaj: based on NSTRSO from Jeppe.
c Turbo routine for generating strings that fulfills RAS
c constraints
c
c Number of strings per type and symmetry is given in NSSO
c
c Jeppe Olsen May 1989
c Last revision November 1989
c
c NEL is number of active electrons
c
#include <implicit.h>
#include <priunit.h>
      DIMENSION IOC(NEL),NSSO(NOCTYP,MAXSYM), NSTASM(MAXSYM)
      INTEGER   ORBSYM(*)
c
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
c Loop over possible partitionings between RAS1,RAS2,RAS3
      CALL ISETVC(NSSO,0,MAXSYM*NOCTYP)
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
c Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD(IOC,IEL1,IORB1F,IORB1L,NONEW1,NTEST)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE. 25) THEN
           WRITE(LUPRI,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
c Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2,NTEST)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE. 25) THEN
             WRITE(LUPRI,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
c Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3,NTEST)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
c Next string has been constructed , Enlist it |.
             NSTRIN = NSTRIN + 1
c. Symmetry
C                   MSYMPR(STRING,NEL,ORBSYM)
             ISYM = MSYMPR(IOC,NEL,ORBSYM)
c. Occupation type
C                   IOCTYP(STRING,IAB,IFLAG)
             ITYP = IOCTYP(IOC,IAB,1)
             IF( NTEST .GE. 20 ) THEN
               WRITE(LUPRI,*) ' RAS 1 2 3 string; IEL1,IEL2,IEL3 =',
     *            IEL1,IEL2,IEL3
               WRITE(LUPRI,*) ' Symmetry and type of string ',ISYM,ITYP
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
c
             NSSO(ITYP,ISYM) = NSSO(ITYP,ISYM)+ 1
c
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
c
      MXST = 0
      DO 2020 ISYM = 1,MAXSYM
         NSTR = 0
         DO 2010 ITYP = 1,NOCTYP
            NSTR = NSTR + NSSO(ITYP,ISYM)
 2010    CONTINUE
         MXST = MAX(MXST,NSTR)
         NSTASM(ISYM) = NSTR
 2020 CONTINUE
      IF (NTEST .GE. 2) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' *** Output from NUMSTR ***   IAB =',IAB
        WRITE(LUPRI,*) ' Number of strings generated       ', NSTRIN
        WRITE(LUPRI,*) ' Max number of strings per symmetry',MXST
        WRITE(LUPRI,*)
        IF(NTEST .GT. 2 ) THEN
         WRITE(LUPRI,*) ' Number of strings per sym (COL) and '//
     &                  'type (ROW)'
         WRITE(LUPRI,*) '====================================='//
     &                  '==========='
         CALL IWRTMA(NSSO,NOCTYP,MAXSYM,NOCTYP,MAXSYM)
        END IF
      END IF
C
      RETURN
      END
C  /* Deck msympr */
      INTEGER FUNCTION MSYMPR(STRING,NEL,ORBSYM)
c
c 910415-hjaaj MSYMPR routine for SIRIUS
c
c Symmetry of string STRING
c
      INTEGER STRING(NEL), ORBSYM(*), SYMPRO(8,8)
c
      DATA SYMPRO/1,2,3,4,5,6,7,8,
     &            2,1,4,3,6,5,8,7,
     &            3,4,1,2,7,8,5,6,
     &            4,3,2,1,8,7,6,5,
     &            5,6,7,8,1,2,3,4,
     &            6,5,8,7,2,1,4,3,
     &            7,8,5,6,3,4,1,2,
     &            8,7,6,5,4,3,2,1 /
c
      ITOTSM = 1
      DO 10 IEL = 1, NEL
         IELSYM = ORBSYM( STRING(IEL) )
         ITOTSM = SYMPRO( ITOTSM,IELSYM )
   10 CONTINUE
c
      MSYMPR = ITOTSM
      RETURN
      END
C  /* Deck memdet */
      SUBROUTINE MEMDET(KFREE,NTEST)
C
C ALLOCATE MEMORY FOR DETERMINANT- PROGRAM
C
c L.R. October 1990, Dimensions for KIPNSA,KIPNSB changed
c                    RAS always assumed
C April 1991 hjaaj: put single excitation info. last
#include <implicit.h>
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
#include <ciinfo.h>
#include <priunit.h>
C
      NORB = NORB1 + NORB2 + NORB3
C
C
C**1 :  OCCUPATION OF STRINGS
C
      CALL MEMADD(KIASTR,NAEL*NASTR,KFREE,1)
      CALL MEMADD(KSTBAA,MAXSYM,KFREE,1)
      CALL MEMADD(KSTASA,MAXSYM,KFREE,1)
C
C.. Matrices for lexical ordering of strings
C
      CALL MEMADD(KZA,NORB*NAEL,KFREE,1)
      IF(EQUAL) THEN
        KZB = KZA
      ELSE
        CALL MEMADD(KZB,NORB*NBEL,KFREE,1)
      END IF
C.  Symmetry of orbitals in RAS order
      CALL MEMADD(KORBSM,NORB,KFREE,1)
C
C.. Arrays giving symmetry order from lexical ordering
C
      MAXSTR = NASTR
      CALL MEMADD(KIPNSA,MAXSTR,KFREE,1)
      IF(EQUAL) THEN
        KIPNSB = KIPNSA
      ELSE
        MAXSTR = NBSTR
        CALL MEMADD(KIPNSB,MAXSTR,KFREE,1)
      END IF
C
      CALL MEMADD(KTPFSA,NASTR,KFREE,1)
      IF( EQUAL ) THEN
        KTPFSB = KTPFSA
      ELSE
        CALL MEMADD(KTPFSB,NBSTR,KFREE,1)
      END IF
C
C*1A : NUMBERS AND POINTERS FOR STRINGS DIVIDED INTO SYMMETRY AND CLASS
C
C      NSSOA(ISYM,IOCTP) : NUMBER OF A STRINGS OF SYM ISYM AND OCC TYPE
C      ISSOA(ISYM,IOCTP) : FIRST ELEMENT OF SYM ISYM AND OCTYPE IOCTP
C** POINTERS FOR SYMMETRY AND OCCUPATION CLASS
       CALL MEMADD(KNSSOA,MAXSYM*NOCTPA,KFREE,1)
       CALL MEMADD(KISSOA,MAXSYM*NOCTPA,KFREE,1)
      IF( EQUAL) THEN
        KIBSTR = KIASTR
        KSTBAB = KSTBAA
        KSTASB = KSTASA
        KNSSOB = KNSSOA
        KISSOB = KISSOA
      ELSE
        CALL MEMADD(KIBSTR,NBEL*NBSTR,KFREE,1)
        CALL MEMADD(KSTBAB,MAXSYM,KFREE,1)
        CALL MEMADD(KSTASB,MAXSYM,KFREE,1)
        CALL MEMADD(KNSSOB,MAXSYM*NOCTPB,KFREE,1)
        CALL MEMADD(KISSOB,MAXSYM*NOCTPB,KFREE,1)
      END IF
C
C**2 :  SYMMETRY OFFSETS FOR CI VECTORS
C
      CALL MEMADD(KCOFF, MAXSYM,KFREE,1)
      CALL MEMADD(KHCOFF,MAXSYM,KFREE,1)
      CALL MEMADD(KIOCOC,NOCTPA*NOCTPB,KFREE,1)
      CALL MEMADD(KICSO, MAXSYM*NOCTPA,KFREE,1)
      CALL MEMADD(KIHCSO,MAXSYM*NOCTPA,KFREE,1)
      CALL MEMADD(KICOOS,MAXSYM*NOCTPA*NOCTPB,KFREE,1)
      CALL MEMADD(KIHOOS,MAXSYM*NOCTPA*NOCTPB,KFREE,1)
      CALL MEMADD(KCDTAS,MAXSYM,KFREE,1)
      CALL MEMADD(KHDTAS,MAXSYM,KFREE,1)
C
C**3:  SYMMETRY OF DISTRIBUTION AND OTHER INFO ON EXCITATIONS
C
      CALL MEMADD(KISSYM,NORB**2,KFREE,1)
      CALL MEMADD(KKLTP, NORB**2,KFREE,1)
      CALL MEMADD(KICREA,NORB**2,KFREE,1)
      CALL MEMADD(KIANNI,NORB**2,KFREE,1)
      CALL MEMADD(KKLCAN,NORB**2,KFREE,1)
C
C.. SPACE CLASSIFICATION OF ORBITALS
      CALL MEMADD(KTPFOB,NORB,KFREE,1)
C
C** SPACE FOR SWITCHING BETWEEN LUNAR AND
C   SIRIUS ORDER OF ORBITALS
C
      CALL MEMADD(KLTSOB,NORB,KFREE,1)
      CALL MEMADD(KSTLOB,NORB,KFREE,1)
C
C
C**4 : INFORMATION ABOUT SINGLE EXCITATIONS
C
      CALL MEMADD(KTAIJ, NAEXCI,KFREE,1)
      CALL MEMADD(KTATO, NAEXCI,KFREE,1)
      CALL MEMADD(KTASYM,(MAXSYM+1)*NASTR,KFREE,1)
C NUMBER OF EXCITATIONS PER STRING
C EXCITATIONS DIVIDED INTO STRING AND CLASS
C
      CALL MEMADD(KNXFSA,NASTR,KFREE,1)
      IF( EQUAL ) THEN
        KTBIJ  = KTAIJ
        KTBTO  = KTATO
        KTBSYM = KTASYM
        KNXFSB = KNXFSA
      ELSE
        CALL MEMADD(KTBIJ, NBEXCI,KFREE,1)
        CALL MEMADD(KTBTO, NBEXCI,KFREE,1)
        CALL MEMADD(KTBSYM,(MAXSYM+1)*NBSTR,KFREE,1)
        CALL MEMADD(KNXFSB,NBSTR,KFREE,1)
      END IF
C
      IF( NTEST .GE. 3 ) THEN
        WRITE(LUPRI,'(//A/)') ' >>> Output from MEMDET <<<'
        WRITE(LUPRI,*)   'EQUAL =',EQUAL
        WRITE(LUPRI,*)   'KIASTR,KSTBAA,KSTASA and KIBSTR,KSTBAB,KSTASB'
        WRITE(LUPRI,*)    KIASTR,KSTBAA,KSTASA
        WRITE(LUPRI,*)    KIBSTR,KSTBAB,KSTASB
        WRITE(LUPRI,*)   'KZA,   KZB,   KORBSM'
        WRITE(LUPRI,*)    KZA,   KZB,   KORBSM
        WRITE(LUPRI,*)   'KIPNSA,KTPFSA,KIPNSB,KTPFSB'
        WRITE(LUPRI,*)    KIPNSA,KTPFSA,KIPNSB,KTPFSB
        WRITE(LUPRI,'(A,4I8)') ' KNSSOA  KNSSOB  KISSOA  KISSOB',
     &                       KNSSOA,KNSSOB,KISSOA,KISSOB
        WRITE(LUPRI,'(A,2I8)') ' KCOFF   KHCOFF',KCOFF,KHCOFF
        WRITE(LUPRI,'(A,4I8)') ' KICSO   KIHCSO  KICOOS  KIHOOS',
     &                       KICSO,KIHCSO,KICOOS,KIHOOS
        WRITE(LUPRI,*)   'KISSYM,KTPFOB ',KISSYM,KTPFOB
        WRITE(LUPRI,*)   'KLTSOB,KSTLOB ',KLTSOB,KSTLOB
        WRITE(LUPRI,*)
     &   'KTAIJ,KTATO,KTASYM,KNXFSA and KTBIJ,KTBTO,KTBSYM,KNXFSB'
        WRITE(LUPRI,*)    KTAIJ,KTATO,KTASYM,KNXFSA
        WRITE(LUPRI,*)    KTBIJ,KTBTO,KTBSYM,KNXFSB
      END IF
C
      IF( NTEST .GE. 2 ) THEN
        WRITE(LUPRI,'(/A,I12)')
     &    ' Space needed for string information ( r*8 words ) ',KFREE-1
      END IF
      RETURN
      END
C  /* Deck detfo */
      SUBROUTINE DETFO(NEL,NORB,MULTSX,ICSYM,NTEST,
     &                 MELMN1,MELMX3,NRAS1,NRAS2,NRAS3,
     &                 NCDET,NCCSF,LCINDX,LLOCA,LLOCB,WORK,LFREE,
     &                 MELMX1,MELMN3,NOSYM,ICSF,NCNSM,REFSPC)
C
C       number of strings consistent with imposed constraints
C
C       required amount of storage for string information
C
C** L.R. JAN 12 1988 / J.O.
C
C NOSYM added Nov '88
#include <implicit.h>
      INTEGER   CLCTYP
      DIMENSION WORK(*)
      DIMENSION NRAS1(8),NRAS2(8),NRAS3(8)
      LOGICAL   REFSPC
C
#include <mxsmob.h>
#include <mxpdim.h>
#include <maxash.h>
C
#include <priunit.h>
#include <spinfo.h>
#include <ciinfo.h>
#include <detbas.h>
#include <strnum.h>
#include <mxblk.h>
C
      INTEGER SYMPRO(8,8)
      DATA SYMPRO/1,2,3,4,5,6,7,8,
     &            2,1,4,3,6,5,8,7,
     &            3,4,1,2,7,8,5,6,
     &            4,3,2,1,8,7,6,5,
     &            5,6,7,8,1,2,3,4,
     &            6,5,8,7,2,1,4,3,
     &            7,8,5,6,3,4,1,2,
     &            8,7,6,5,4,3,2,1 /
C
C     Check max number of active orbitals
C
      IF (NORB .GT. MAXASH) THEN
         WRITE (LUPRI,'(/A,2(/A,I5))')
     &      ' DETFO ERROR, this version of CI module cannot'//
     &      ' handle so many orbitals.',
     &      ' The total number of active orbitals:   ',NORB,
     &      ' is greater than fixed parameter MAXASH=',MAXASH
         CALL QUIT('DETFO: Too many active orbitals for this version')
      END IF
C
C     Transfer input to common /CIINFO/
C
      NEL1MN = MELMN1
      NEL1MX = MELMX1
      NEL3MN = MELMN3
      NEL3MX = MELMX3
C
C     Transfer input to common /SPINFO/
C
      MULTS  = MULTSX
C
C *** 1.05 : NUMBER OF ORBITALS IN CLASSES ETC.............
C
C LARGEST SYMMETRY NUMBER OCCURING
      MAXSYM = 0
      DO 20 ISYM = 1, 8
        IF((NRAS1(ISYM).NE.0) .OR.
     &     (NRAS2(ISYM).NE.0) .OR. (NRAS3(ISYM) .NE. 0 ) )
     &  MAXSYM = ISYM
   20 CONTINUE
C
      IF (MAXSYM .GT. 4 ) THEN
         MAXSYM = 8
      ELSE IF(MAXSYM .GT. 2 ) THEN
         MAXSYM = 4
      END IF
C
      NORB1 = 0
      NORB2 = 0
      NORB3 = 0
      DO 30 ISYM = 1, MAXSYM
        NORB1 = NORB1 +NRAS1(ISYM)
        NORB2 = NORB2 +NRAS2(ISYM)
        NORB3 = NORB3 +NRAS3(ISYM)
   30 CONTINUE
      IORB1F = 1
      IORB1L = NORB1
      IORB3F = NORB1+NORB2 +1
      IORB3L = NORB
C
      IF (NTEST .GE. 5) THEN
         WRITE(LUPRI,'(//A//A,I3/)') ' >>> TEST OUTPUT FROM DETFO <<<',
     &      ' MAXSYM  ',MAXSYM
         WRITE(LUPRI,*) ' NORB1 NORB2 NORB3 ',NORB1,NORB2,NORB3
         WRITE(LUPRI,*) ' IORB1F IORB1L IORB3F IORB3L '
         WRITE(LUPRI,*)   IORB1F, IORB1L, IORB3F, IORB3L
      END IF
C
C** 1.1 : NUMBER OF STRINGS , SINGLE EXCITATIONS ...
C
      CALL STRDIM(NORB,NEL,MS2,WORK(1),NTEST )
C
C** 1.2 : DYNAMIC MEMORY ALLOCATION
C
      KFREE = 1
      CALL MEMDET(KFREE,NTEST)
      LCINDX = KFREE - 1
      KFREE  = KTAIJ
C910408-hjaaj Information about single excitations not used in DETFO
      LNEED  = KFREE + 4*NORB
C     ... 4*NORB needed in CSFDIM
      IF (LNEED .GT. LFREE) THEN
         WRITE (LUPRI,'(//A,2(/A,I10))')
     *      ' *** FATAL ERROR in DETFO,'//
     *         ' insufficient memory for string information',
     *      '     Need more than',LNEED,
     *      '     Available now ',LFREE
         CALL QUIT('DETFO: insufficient memory for string information.')
      END IF
C
C** 1.3 : SYMMETRY OF ORBITALS
C         ORDER ORBITALS AFTER OCCUPATION FIRST, SYMMETRY SECOND
      CALL ZORBSM(NRAS1,NRAS2,NRAS3,MAXSYM,WORK(KORBSM),NTEST,
     &            WORK(KLTSOB),WORK(KSTLOB),NORB,NOSYM)
C
C*** 2 : STRINGS, SINGLE EXCITATIONS UND SO WEITER FOR ALPHA ELECTRONS
C
      NLMXAB = MAX(NAEL,NBEL)
      CALL MEMADD(KIOC,NLMXAB,KFREE,1)
      IEL1MX = MIN(NORB1,NEL1MX+1)
      NASTR1 = NASTR
      CALL NUMSTR(NASTR,NAEL,NORB1,NORB2,NORB3,
     &            NL1MNA,IEL1MX,NL3MNA,NEL3MX,
     &            WORK(KNSSOA),WORK(KSTASA),MXSTA,
     &            NOCTPA,MAXSYM,WORK(KORBSM),WORK(KIOC),1,NTEST)
      IF (NASTR1 .NE. NASTR) THEN
         WRITE (LUPRI,*) 'Fatal error, NASTR from NUMSTR:',NASTR
         WRITE (LUPRI,*) '             NASTR from NUMST2:',NASTR1
         CALL QUIT('ERROR: different NASTR from NUMSTR and NUMST2')
      END IF
C     CALL NUMSTR(NSTRIN,NEL,NORB1,NORB2,NORB3,
C    &            NELMN1,NELMX1,NELMN3,NELMX3,
C    &            NSSO,NSTASM,MXST,
C    &            NOCTYP,MAXSYM,ORBSYM,IOC,IAB,NTEST)
Cold  CALL NUMSTR(NASTR,NAEL,NORB,NOCTPA,NL1MNA,WORK(KORBSM),
Cold &            WORK(KNSSOA),WORK(KSTASA),MAXSYM,1,NTEST,
Cold &            MXSTA)
Cold  CALL NUMSTR(NSTRIN,NEL,NORB,NOCTP,NLMNST,ORBSYM,
Cold &            NSSO,NSTASM,MAXSYM,IAB,NTEST,MXSASM)
      IF(.NOT.EQUAL) THEN
Cold    CALL NUMSTR(NBSTR,NBEL,NORB,NOCTPB,
Cold &              NL1MNB,WORK(KORBSM),WORK(KNSSOB),
Cold &              WORK(KSTASB),MAXSYM,2,NTEST,MXSTB)
        NBSTR1 = NBSTR
        CALL NUMSTR(NBSTR,NBEL,NORB1,NORB2,NORB3,
     &              NL1MNB,IEL1MX,NL3MNB,NEL3MX,
     &              WORK(KNSSOB),WORK(KSTASB),MXSTB,
     &              NOCTPB,MAXSYM,WORK(KORBSM),WORK(KIOC),2,NTEST)
        IF (NBSTR1 .NE. NBSTR) THEN
         WRITE (LUPRI,*) 'Fatal error, NBSTR from NUMSTR:',NBSTR
         WRITE (LUPRI,*) '             NBSTR from NUMST2:',NBSTR1
         CALL QUIT('ERROR: different NBSTR from NUMSTR and NUMST2')
        END IF
        MXSASM = MAX(MXSTA,MXSTB)
      ELSE
        MXSASM = MXSTA
      END IF
C
C*** 3 :       NUMBER OF DETERMINANTS
C
         CALL STRPA2(WORK(KNSSOA),WORK(KNSSOB),NASTR,NBSTR,
     &               MAXSYM,NOCTPA,NOCTPB,WORK(KIOCOC),
     &               NL1MNA,NL1MNB,NEL1MN,NEL1MX,
     &               NEL3MN,NEL3MX,NORB1,NORB2,NORB3,
     &               IOCPTA,IOCPTB,SYMPRO,NAEL,NBEL,NEL,
     &               REFSPC,NDTASM,MXNDT ,MXVBLK,NTEST)
      NCDET = NDTASM(ICSYM)
C
C... Memory needed for csfstuff
C
C     CSFDIM(NACTOB,NACTEL,MULTP,MAXSYM,SYMPRO,ORBSYM,IWORK,
C    &       LLCSF,LCSFFO,KFREE,NCNSM,NTEST)
      IF(ICSF .NE. 0 ) THEN
        KCSF = LCINDX+1
        CALL CSFDIM(NORB,NEL,MULTS,MAXSYM,SYMPRO,WORK(KORBSM),
     &              WORK(KFREE),LLCSF,LCSFFO,KCSF,NCNSM,NTEST)
        LCINDX = KCSF-1
        NCCSF  = NCSASM(ICSYM)
      ELSE
        LLCSF  = 0
        LCSFFO = 0
        NCCSF  = NCDET
      END IF
C
C*** 4 : NEEDED  LOCAL MEMORY
C
C..4.1 : DETINFO
      LA  =  NASTR
      LB  =  NBSTR
      LDETIN = MAX(LA,LB)
C..4.2 : CI-DIAGONAL
      LCIDIA = 5*NORB + 2*NORB**2
C..4.3 : GENERATION OF SIGMA VECTOR
Cjo   LSIGMA = 3*NORB**2 + 7*MXSASM + MXVBLK
C     Revised to correspond to correct cisigd/cisgd2 allocation,
C     including csf usage/890824-hjaaj; revised for cisig9/910408-hjaaj
Cold  LSIGMA = NORB**2 + 7*MXSASM + MXVBLK + MAX(2*NBSTR,NORB**2)
      LSIGMA = NORB**2*(2+MXSMOB) + MXSASM*(1+2*MXSMOB) + MXVBLK
C     add 4 integer items to LSIGMA
      CALL MEMADD(KDUM,MXSASM*(2+2*MXSMOB)+4,LSIGMA,1)
      IF (ICSF .NE. 0 ) LSIGMA = LSIGMA + 2*MXNDT
C..4.4 : GENERATION OF DENSITY MATRICES
Cjo   LDENSI = 3*NORB**2 + 2*NBSTR
C     Revised to correspond to correct MAKTDM allocation,
C     including csf usage/890824-hjaaj; revised for densit/910408-hjaaj
Cold  LDENSI = 3*NORB**2 + MAX(NORB**2,2*NBSTR)
      LDENSI = NORB**2 + MXVBLK
C     add 3 integer items to LDENSI
      CALL MEMADD(KDUM,MAXSYM*(2+MAXSYM)+3,LDENSI,1)
      IF (ICSF .NE. 0) LDENSI = MAX(LDENSI,MXNDT) + 2*MXNDT
C..4.5 : construction of CSF stuff
C     LLCSF = LLCSF
C.. LARGEST BLOCK OF LOCAL MEMORY NEEDED ( EXCEPT IN DETINF )
      LLOCA = MAX(LDENSI,LSIGMA)
      LLOCB = 0
C
      IF (NTEST .GE. 5) WRITE(LUPRI,'(/A/6I12)')
     &' Local memory : LDETIN  LCIDIA  LSIGMA  LDENSI   LLOCA   LLCSF',
     &                 LDETIN,LCIDIA,LSIGMA,LDENSI, LLOCA, LLCSF
C
      RETURN
      END
C  /* Deck zorbsm */
      SUBROUTINE ZORBSM(NRAS1,NRAS2,NRAS3,MAXSYM,ORBSYM,NTEST,
     &                  ILTSOB,ISTLOB,NORB,NOSYM)
C
C SET ORBSYM ( SYMMETRY OF ORBITALS ARRANGED AFTER FIRST OCCUPATION TYPE,
C              SECOND AFTER SYMMETRY )
C SET ARRAYS FOR POINTING BETWEEN LUNAR AND SIRIUS ORDER OF ORBITALS
C
C     ILTSOB(K) GIVES SIRIUS NUMBER OF ORBITAL WITH LUNAR NUMBER K
C     ISTLOB(K) GIVES LUNAR  NUMBER OF ORBITAL WITH SIRIUS NUMBER K
C
#include <implicit.h>
#include <priunit.h>
      INTEGER ORBSYM(1)
      DIMENSION NRAS1(MAXSYM),NRAS2(MAXSYM),NRAS3(MAXSYM)
      DIMENSION ILTSOB(*),ISTLOB(*)
#include <cbreor.h>
C
        MINTP = 1
        MAXTP = 3
C
      IORB = 0
      DO 40 IOCCTP = MINTP,MAXTP
        ISBAS = 1
        ISOBAS = ISBAS
        DO 30 ISYM = 1, MAXSYM
          IF(ISYM.NE.1) THEN
            ISBAS = ISBAS + NRAS1(ISYM-1)+NRAS2(ISYM-1)+NRAS3(ISYM-1)
            ISOBAS = ISBAS
          END IF
          IF(IOCCTP .GE. 2 ) ISOBAS = ISOBAS + NRAS1(ISYM)
          IF(IOCCTP .EQ. 3 ) ISOBAS = ISOBAS + NRAS2(ISYM)
          IF(IOCCTP .EQ. 1 ) KMXORB = NRAS1(ISYM)
          IF(IOCCTP .EQ. 2 ) KMXORB = NRAS2(ISYM)
          IF(IOCCTP .EQ. 3 ) KMXORB = NRAS3(ISYM)
          DO 25 KORB = 1,  KMXORB
            IORB = IORB + 1
            IF( NOSYM .NE. 0 ) THEN
              ORBSYM(IORB) = 1
            ELSE
              ORBSYM(IORB) = ISYM
            END IF
            ILTSOB(IORB) = ISOBAS-1+KORB
   25     CONTINUE
   30   CONTINUE
   40 CONTINUE
c
      SLREOR = .FALSE.
      DO 60 IORB = 1,NORB
        ISTLOB(ILTSOB(IORB)) = IORB
        IF (ILTSOB(IORB) .NE. IORB) SLREOR = .TRUE.
  60  CONTINUE
C
      IF (NTEST .GE. 3) THEN
         WRITE(LUPRI,'(/A)') ' ZORBSM: ORBSYM ARRAY '
         CALL IWRTMA(ORBSYM,1,NORB,1,NORB)
         IF (SLREOR) THEN
            WRITE(LUPRI,'(/A)') ' LUNAR TO SIRIUS ORDERING '
            CALL IWRTMA(ILTSOB,1,NORB,1,NORB)
            WRITE(LUPRI,'(/A)') ' SIRIUS TO LUNAR ORDERING '
            CALL IWRTMA(ISTLOB,1,NORB,1,NORB)
         ELSE
            WRITE(LUPRI,'(/A)') ' LUNAR AND SIRIUS ORDERING ARE '//
     &                           'IDENTICAL.'
         END IF
      END IF
C
      RETURN
      END
