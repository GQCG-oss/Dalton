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
C  /* Deck setdis */
      SUBROUTINE SETDIS(SYMPRO,DISSYM,ORBSYM,NORB,MAXSYM,
     &                  KLTP,ICREA,IANNI,KLCAN)
C
C INFORMATION ABOUT DISTRIBUTIONS
C
C DISSYM(I,J) = SYMPRO( ORBSYM(I),ORBSYM(J) )
      INTEGER SYMPRO(8,8),ORBSYM(NORB),DISSYM(NORB,NORB)
      INTEGER ICREA(*),IANNI(*),KLTP(*),KLCAN(*)
C
      DO 100 I = 1, NORB
        DO 90 J = 1, NORB
          DISSYM(I,J) = SYMPRO(ORBSYM(I),ORBSYM(J) )
          IJ = (I-1)*NORB + J
          JI = (J-1)*NORB + I
          KLTP(IJ) = JI
          ICREA(IJ) = I
          IANNI(IJ) = J
          IF(I .GE. J ) THEN
           KLCAN(IJ) = I*(I-1)/2+J
          ELSE
           KLCAN(IJ) = J*(J-1)/2 + I
          END IF
   90   CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck occexc */
      SUBROUTINE OCCEXC(NEXCI,NSTRIN,TTO,TIJ,TSYM,MAXSYM,DISSYM,SYMPRO,
     &                  NEXFST,IAB,NOCTP,STRING,NEL,NTEST)
c
c Reorder single excitations according to symmetry
c
c  October 1990 version, Jeppe Olsen
c
      INTEGER   TIJ,TTO,TSYM,DISSYM,SYMPRO
      INTEGER   STRING(NEL,NSTRIN)
      DIMENSION TTO(NEXCI),TIJ(NEXCI)
      DIMENSION TSYM(MAXSYM+1,NSTRIN)
      DIMENSION SYMPRO(8,8),DISSYM(*)
      DIMENSION NEXFST(NSTRIN)
      PARAMETER (MXEXCI = 5000)
      DIMENSION I1(MXEXCI),I2(MXEXCI)
#include <priunit.h>
c
      DO 1000 ISTRIN = 1, NSTRIN
        IBASE = 1
        NNEXCI = NEXFST(ISTRIN)
        IF (NNEXCI .GT. MXEXCI) THEN
           WRITE (LUPRI,*) ' Local arrays I1 and I2 too small in OCCEXC'
           WRITE (LUPRI,*) ' Need at least ',NNEXCI,'; have',MXEXCI
           CALL QUIT('ERROR: local dimension too small in OCCEXC')
        END IF
        DO 900 ISYM = 1, MAXSYM
        DO 890 IOCTP = 1,NOCTP
          IEFF = 0
c
          DO 800 I = TSYM(ISYM,ISTRIN),TSYM(ISYM+1,ISTRIN)-1
            IF(TTO(I) .GT. NSTRIN ) THEN
             JSTRIN = TTO(I) - NSTRIN
            ELSE
             JSTRIN = TTO(I)
            END IF
c. Occupation type of string
            JOCTP = IOCTYP(STRING(1,JSTRIN),IAB,1)
            IF(JOCTP .EQ. IOCTP ) THEN
              I1(IBASE + IEFF ) = TIJ(I)
              I2(IBASE + IEFF ) = TTO(I)
              IEFF = IEFF + 1
            END IF
  800     CONTINUE
              IBASE = IBASE + IEFF
  890   CONTINUE
  900   CONTINUE
        IOFF = TSYM(1,ISTRIN)
        CALL ICOPVE(I1,TIJ(IOFF),NNEXCI)
        CALL ICOPVE(I2,TTO(IOFF),NNEXCI)
 1000 CONTINUE
c
      IF ( NTEST .GE. 10 ) THEN
        MAXPR = 60
        NPR = MIN(NSTRIN,MAXPR)
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' OUTPUT FROM OCCEXC '
        WRITE(LUPRI,*) '===================='
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' TIJ ARRAY ( ALLOWED E(IJ) ) '
        WRITE(LUPRI,*)
        IOFF = 1
        DO 1234 ISTRIN = 1,NPR
           LEXCI = NEXFST(ISTRIN)
           WRITE(LUPRI,'(A,I4,A,(T29,10I5))')
     &       ' STRING..',ISTRIN,' EXCITATIONS : ',
     &       ( TIJ(IOFF-1+IEX),IEX = 1, LEXCI)
           IOFF = IOFF + LEXCI
 1234    CONTINUE
         WRITE(LUPRI,*) ' TTO ARRAY ( EXCITED STRINGS ) '
         IOFF = 1
         DO 1235 ISTRIN = 1, NPR
           LEXCI = NEXFST(ISTRIN)
           WRITE(LUPRI,'(A,I4,A,(T29,10I5))')
     &       ' STRING..',ISTRIN,' new string..',
     &       ( TTO(IOFF-1+IEX),IEX = 1, LEXCI)
           IOFF = IOFF + LEXCI
 1235    CONTINUE
C
      END IF
C
      RETURN
      END
C  /* Deck imxsym */
      INTEGER FUNCTION IMXSYM(ORBSYM,NORB)
C
C FIND LARGEST SYMMETRY OCCURING IN ORBSYM OR PRODUCTS OF ORBSYM
C
      INTEGER ORBSYM(NORB)
#include <priunit.h>
C
      MAXSYM = 1
      DO 100 IORB = 1,NORB
        IF( ORBSYM( IORB ) .GT. MAXSYM ) MAXSYM = ORBSYM ( IORB )
  100 CONTINUE
      IF( MAXSYM .GT. 1 ) MAXSYM = MAX(MAXSYM,2)
      IF( MAXSYM .GT. 2 ) MAXSYM = MAX(MAXSYM,4)
      IF( MAXSYM .GT. 4 ) MAXSYM = MAX(MAXSYM,8)
C
C
      IMXSYM = MAXSYM
C
      RETURN
      END
C  /* Deck occstr */
      SUBROUTINE OCCSTR(STRING,STRIN2,NORB,NEL,NEWORD,NSTRIN,SYMPRO,
     &                  ORBSYM,NSTASM,ISTBAS,NEWORH,MAXSTR,NOCCTP,
     &                  NSTASO,ISTASO,IAB,MAXSYM,NTEST)
C
C   REORDER SYMMETRY ORDERED STRINGS TO
C   STRINGS THAT ARE SYMMETRY AND OCCUPATION ORDERED
C ( FIRST SYMMETRY THEN OCCUPATION CLASS )
C
C TWO POINTERS ARE INVOLVED
C    NEWORH : FROM HANDY NUMBER TO SYMMETRY ORDERED ON INPUT
C             FROM HANDY NUMBER TO SYMMETRY AND OCCUPATION ORDERED ON OUTPUT
C    NEWORS : LOCAL POINTER FROM SYMMETRY NUMBER TO OCCUPATION ORDERED NUMBER
C
      INTEGER   STRING(NEL,NSTRIN),STRIN2(NEL,NSTRIN)
      INTEGER   SYMPRO(8,8),ORBSYM(NORB)
      DIMENSION NEWORD(NSTRIN),NSTASM(MAXSYM),ISTBAS(MAXSYM)
      DIMENSION NEWORH(MAXSTR)
      DIMENSION NSTASO(NOCCTP,MAXSYM),ISTASO(NOCCTP,MAXSYM)
#include <priunit.h>
C
#include <ciinfo.h>
C
C OUTPUT : NSTASO(IOCC,ISYM) : NUMBER OF STRINGS OF SYMMETRY ISYM
C                              AND OCCUPATION IOCC
C          ISTASO(IOCC,ISYM) : OFFSET FOR STRING OF SYMMETRY ISYM
C                              AND OCCUPATION IOCC
C          NEWORH(J) : NEW NUMBER FOR STRING WITH HANDY NUMBER J
C          STRING    : STRINGS SYMMETRY AND OCCUPATION ORDERED
C
      IFLAG = 1
C
C** 1 : ARRAYS DESCRIBING NEW ORDER
C
      CALL ISETVC(NSTASO,0,MAXSYM*NOCCTP)
      CALL ISETVC(ISTASO,0,MAXSYM*NOCCTP)
      IBASE = 1
C?    WRITE(LUPRI,*) ' IORB1F,IORB1L,IORB3F,IORB3L'
C?    WRITE(LUPRI,*)   IORB1F,IORB1L,IORB3F,IORB3L
      DO 200 ISYM = 1, MAXSYM
        IF( ISYM .NE. 1 ) IBASE = IBASE + NSTASM(ISYM-1)
        DO 100 ISTRIN = IBASE, IBASE + NSTASM(ISYM) - 1
          IOCCTP =  IOCTYP(STRING(1,ISTRIN),IAB,IFLAG)
          NSTASO(IOCCTP,ISYM) = NSTASO(IOCCTP,ISYM) + 1
          NEWORD(ISTRIN) = IOCCTP
  100   CONTINUE
C
        IF ( NTEST .GE. 10 ) THEN
           WRITE(LUPRI,*) ' NUMBER OF STRINGS OF EACH OCCUPATION'
           WRITE(LUPRI,*) '======================================'
           CALL IWRTMA(NSTASO(1,ISYM),1,NOCCTP,1,NOCCTP)
        END IF
C
C** 2: OFFSET FOR EACH SYMMETRY AND OCCUPATION CLASS
C
        IF ( ISYM . EQ.  1 ) THEN
           ISTASO(1,ISYM) = 1
        ELSE
           ISTASO(1,ISYM) = ISTASO(1,ISYM-1) + NSTASM(ISYM-1)
        END IF
C
        DO 110 IOCCTP = 2,NOCCTP
          ISTASO(IOCCTP,ISYM) = ISTASO(IOCCTP-1,ISYM)
     &                        + NSTASO(IOCCTP-1,ISYM)
  110  CONTINUE
C
        IF ( NTEST .GE. 10 ) THEN
           WRITE(LUPRI,*) ' OFFSETS FOR EACH OCCUPATION  '
           CALL IWRTMA(ISTASO(1,ISYM),1,NOCCTP,1,NOCCTP)
        END IF
C** 3 : REORDERING THE STRINGS
        DO 450 IOCCTP = 1,NOCCTP
          IEFF = 0
          DO 400 ISTRIN = IBASE,IBASE + NSTASM(ISYM) - 1
          IF(NEWORD(ISTRIN) .EQ. IOCCTP ) THEN
            IEFF = IEFF + 1
            NEWORD(ISTRIN) = -(ISTASO(IOCCTP,ISYM)-1 + IEFF )
            CALL ICOPVE(STRING(1,ISTRIN),STRIN2(1,-NEWORD(ISTRIN)),NEL)
          END IF
  400   CONTINUE
  450   CONTINUE
        IF ( NSTASM( ISYM ) .NE. 0 ) THEN
          CALL ICOPVE(STRIN2(1,IBASE),STRING(1,IBASE),NSTASM(ISYM)*NEL)
          CALL ISCAVE(NEWORD(IBASE) , -1 , NSTASM( ISYM ) )
        END IF
  200 CONTINUE
C?      WRITE(LUPRI,*) ' AFTER LOOP 200 '
C
C** TWO POINTERS ARE AVAILABLE NOW :
C     NEWORH : FROM HANDY NUMBER TO SYMMETRY NUMBER
C     NEWORD : FROM REDUCED NUMBER TO OCCUPATION NUMBER
C     MERGE THESE SO
C     NEWORH : FROM HANDY NUMBER TO OCCUPATION NUMBER
C
      DO 500 ISTIN = 1,MAXSTR
        IF( NEWORH(ISTIN) .NE. 0 ) NEWORH(ISTIN) = NEWORD(NEWORH(ISTIN))
  500 CONTINUE
C?      WRITE(LUPRI,*) ' AFTER 500 '
C
      IF ( NTEST .GE. 10 ) THEN
        WRITE(LUPRI,*) 'OLD NUMBER TO NEW NUMBER '
        WRITE(LUPRI,*) '======================== '
        MAXPR = 60
        NPR = MIN(NSTRIN,MAXPR)
        CALL IWRTMA(NEWORH,1,NPR,1,MAXSTR)
        WRITE(LUPRI,*) ' REORDERED STRINGS '
        DO 1234 ISTRIN = 1,NPR
          WRITE(LUPRI,'(1X,I4,A,(T11,15I4))')
     &    ISTRIN,' : ',(STRING(IEL,ISTRIN),IEL = 1,NEL)
1234    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck iznum */
      INTEGER FUNCTION IZNUM(IOCC,NEL,Z,MAXEL,MXCORB,NEWORD,NTEST)
C
C ADRESS OF STRING IOCC
C
      INTEGER   IOCC(NEL),NEWORD(*)
      INTEGER   Z(MXCORB,MAXEL)
#include <priunit.h>
C
      IZ = 1
      DO 100 I = 1,NEL
        IZ = IZ + Z(IOCC(I),I)
  100 CONTINUE
C
      IZNUM = NEWORD(IZ)
C
      IF ( NTEST .GE. 20 ) THEN
        WRITE(LUPRI,*)
     &     ' >>> IZNUM <<< ADRESS OF STRING (Old,new) '
     &     ,IZ,IZNUM,' STRING:'
        CALL IWRTMA(IOCC,1,NEL,1,NEL)
      END IF
C
      RETURN
      END
C  /* Deck sexci */
      SUBROUTINE SEXCI(STRING,NSTRIN,NEL,NORB,TIJ,TTO,
     &                 MXCORB,MAXEL,Z,MAXTL,NEWORD,
     &                 NEXFST,IAB,NTEST)
C
C CALCULATE AND STORE INFORMATION ABOUT SINGLE EXCITATIONS OBTAINED FROM
C STRINGS
C
c October 90 : Excitations packed in a onedimensional array
#include <implicit.h>
c. Input
      INTEGER   STRING(NEL,NSTRIN),Z(MXCORB,MAXEL)
c. Output
      INTEGER   TIJ(MAXTL),TTO(MAXTL),NEWORD(NSTRIN),NEXFST(*)
C
#include <priunit.h>
C
#if defined (VAR_STACKMEM)
      INTEGER   OCC(NORB),UNOCC(NORB),ISIGN(NORB)
      INTEGER   IEXP(NORB),STRIN2(NORB)
#else
#include <maxash.h>
      INTEGER   OCC(MAXASH),UNOCC(MAXASH),ISIGN(MAXASH)
      INTEGER   IEXP(MAXASH),STRIN2(MAXASH)
#endif
C
C T*(NUM,STRING GIVES THE FOLLOWING INFORMATION :
C TIJ(NUM,STRING) : IJVALUE FOR SINGLE EXCITATION
C TTO(NUM,STRING) : NUMBER OF 'EXCITED STRING '
C
C NEXFST ( I ) : NUMBER OF ALLOWED SINGLE EXCITATIONS FOR STRING I
C
       CALL ISETVC(NEXFST,NEL,NSTRIN)
#if !defined (VAR_STACKMEM)
      IF (NORB .GT. MAXASH) THEN
         WRITE (LUPRI,'(/A,I5,A,I5,A)')
     *   ' ERROR in SEXCI, NORB =',NORB,' is greater than dimension',
     *   MAXASH,' of local variables OCC,UNOCC,IEXP,ISIGN,STRIN2.'
         CALL QUIT('SEXCI error: increase MAXASH')
      END IF
#endif
      IF ( NTEST .GE. 15 ) THEN
         WRITE(LUPRI,*) ' NEWORD ARRAY IN SEXCI '
         CALL IWRTMA(NEWORD,1,NSTRIN,1,NSTRIN)
      END IF
c
      NSXTOT = 0
      DO 1000 ISTRIN = 1,NSTRIN
c. Expand string
        CALL ISETVC(IEXP(1),0,NORB)
        DO 10 I = 1,NEL
           IEXP(STRING(I,ISTRIN) ) = 1
   10   CONTINUE
c. Occupied and unoccupied orbitals in string
        IOEFF = 0
        IUEFF = 0
        DO 20 I = 1,NORB
          IF( IEXP(I) .EQ. 1 ) THEN
            IOEFF = IOEFF + 1
            OCC(IOEFF) = I
          ELSE
            IUEFF = IUEFF + 1
            UNOCC(IUEFF) = I
          END IF
   20   CONTINUE
c
        IF ( NTEST .GE. 14 ) THEN
          WRITE(LUPRI,*)
     *      ' OCCUPIED AND UNOCCUPIED ORBITALS FOR STRING',ISTRIN
          CALL IWRTMA(OCC,1,NEL,1,NEL)
          CALL IWRTMA(UNOCC,1,NORB-NEL,1,NORB-NEL)
        END IF
c. Sign array
        CALL ISETVC(ISIGN(1),0,NORB)
        KSIGN = (-1)**NEL
        DO 30 I = 1,NEL
         ISIGN(OCC(I) ) = (-1)**(NEL-I)
   30   CONTINUE
        DO 40 I = 1,NORB
          IF(ISIGN(I) .NE. 0 ) THEN
            KSIGN = - KSIGN
          ELSE
            ISIGN(I) = KSIGN
          END IF
   40   CONTINUE
C
        IF ( NTEST .GE. 14 ) THEN
           WRITE(LUPRI,*) 'THE ISIGN ARRAY '
           CALL IWRTMA(ISIGN,1,NORB,1,NORB)
        END IF
c. Getting to business, the single excitations
c
c. Diagonal excitations
        DO 250 IOCC = 1,NEL
          NSXTOT = NSXTOT + 1
          TIJ(NSXTOT) = (OCC(IOCC)-1)*NORB + OCC(IOCC)
          TTO(NSXTOT) = ISTRIN
  250   CONTINUE
c. Offdiagonal excitations
        IJ = NEL
        DO 300 IOCC = 1,NEL
          KOCC = OCC(IOCC)
          KOSIGN = ISIGN(KOCC)
          IEXP(KOCC) = 0
          DO 200 IUNOCC = 1, (NORB-NEL)
            IJ = IJ + 1
            KUNOCC = UNOCC(IUNOCC)
            KUSIGN = ISIGN(KUNOCC)
            KL = (KUNOCC-1)*NORB + KOCC
            IF(KUNOCC.GT.KOCC) THEN
              IIISGN           =   KUSIGN * KOSIGN
            ELSE
              IIISGN           = - KUSIGN * KOSIGN
            END IF
            IEXP(KUNOCC) = 1
            IEFF = 0
            DO 190 I = 1,NORB
              IF(IEXP(I) .NE. 0 ) THEN
                IEFF = IEFF + 1
                STRIN2(IEFF) = I
              END IF
  190       CONTINUE
c. Is new string allowed  ?
            ITYPE = IOCTYP(STRIN2,IAB,1)
            IF ( ITYPE .NE. 0 ) THEN
              NSXTOT = NSXTOT + 1
              TIJ(NSXTOT) = KL
              TTO(NSXTOT) =
     *           IZNUM(STRIN2,NEL,Z,MAXEL,MXCORB,NEWORD,NTEST)
              IF (IIISGN .EQ. -1 )
     &           TTO(NSXTOT) = TTO(NSXTOT) + NSTRIN
               NEXFST(ISTRIN) = NEXFST(ISTRIN ) + 1
            ELSE
              IJ = IJ - 1
            END IF
            IEXP(KUNOCC) = 0
  200     CONTINUE
          IEXP(KOCC) = 1
  300   CONTINUE
C
 1000 CONTINUE
c
      IF(NSXTOT.GT.MAXTL) THEN
        WRITE(LUPRI,*) ' SUBROUTINE SEXCI is in error'
        WRITE(LUPRI,*) ' Number of obtained single excitations '//
     &                 'too large'
        WRITE(LUPRI,*) ' Actual and max number ',NSXTOT,MAXTL
        CALL QUIT('SEXCI ERROR: # single excit. .gt. MAXTL')
      END IF
c
      IF ( NTEST .GE. 5 ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' OUTPUT FROM SEXCI : '
        WRITE(LUPRI,*) '==================='
        WRITE(LUPRI,*) ' Number of single excitations obtained',NSXTOT
      END IF
      IF ( NTEST .GE. 10 ) THEN
        MAXPR = 60
        NPR = MIN(NSTRIN,MAXPR)
        WRITE(LUPRI,*) ' TIJ ARRAY ( ALLOWED E(IJ) )'
        IOFF = 1
        DO 1234 ISTRIN = 1,NPR
           LEXCI = NEXFST(ISTRIN)
           WRITE(LUPRI,'(A,I4,A,(T29,10I5))')
     &       ' STRING..',ISTRIN,' EXCITATIONS : ',
     &       ( TIJ(IOFF-1+IEX),IEX = 1, LEXCI)
           IOFF = IOFF + LEXCI
 1234    CONTINUE
         WRITE(LUPRI,*) ' TTO ARRAY ( EXCITED STRINGS ) '
         IOFF = 1
         DO 1235 ISTRIN = 1, NPR
           LEXCI = NEXFST(ISTRIN)
           WRITE(LUPRI,'(A,I4,A,(T29,10I5))')
     &       ' STRING..',ISTRIN,' new string..',
     &       ( TTO(IOFF-1+IEX),IEX = 1, LEXCI)
           IOFF = IOFF + LEXCI
 1235    CONTINUE
         WRITE(LUPRI,'(A)')
     &  ' NUMBER OF INTERNAL EXCITATIONS FOR EACH STRING '
         WRITE(LUPRI,'(A)')
     &  '================================================'
         CALL IWRTMA(NEXFST,1,NSTRIN,1,NSTRIN)
      END IF
C
      RETURN
      END
C  /* Deck symstr */
      SUBROUTINE SYMSTR(STRING,STRIN2,NORB,NEL,NEWORD,NSTRIN,
     &                  SYMPRO,ORBSYM,NSTASM,ISTBAS,NEWORH,
     *                  MAXSTR,MAXSYM,NTEST,MXSASM)
C
C SYMMETRY ORDER  STRINGS
C
C TWO POINTERS ARE INVOLVED
C    NEWORH : FROM HANDY NUMBER TO REDUCED NUMBER ON INPUT,
C             FROM HANDY NUMBER TO SYMMETRY ORDERED REDEUCED NUMBER
C             ON OUTPUT
C    NEWORS : LOCAL POINTER FROM REDUCED NUMBER TO SYMMETRY ORDERED
C             NUMBER
C
      INTEGER   STRING(NEL,NSTRIN),STRIN2(NEL,NSTRIN)
      INTEGER   SYMPRO(8,8),ORBSYM(NORB)
      INTEGER   NEWORD(NSTRIN),NSTASM(MAXSYM),ISTBAS(MAXSYM)
      INTEGER   NEWORH(MAXSTR)
#include <priunit.h>
C
C OUTPUT : NSTASM(J) : NUMBER OF STRINGS OF SYMMETRY J
C          ISTBAS(J) : OFFSET FOR STRING OF SYMMETRY J
C          NEWORH(J) : NEW NUMBER FOR STRING WITH HANDY NUMBER J
C          STRING    :  STRINGS SYMMETRY ORDERED
C
C** 1 : ARRAYS DESCRIBING NEW ORDER
C
C
      CALL ISETVC(NSTASM,0,MAXSYM)
      DO 100 ISTRIN = 1, NSTRIN
c. IBM compiler does not vectorize this loop correct,kill vectorization
c       IF(ISTRIN.EQ.-1) GOTO 100
C SYMMETRY OF ISTRIN
        ISYM = 1
        DO 10 IEL = 1,NEL
           ISYM = SYMPRO(ISYM,ORBSYM(STRING(IEL,ISTRIN) ) )
   10   CONTINUE
        NSTASM(ISYM) = NSTASM(ISYM) + 1
        NEWORD(ISTRIN) = ISYM
  100 CONTINUE
C
      DO 110 ISYM = 1, MAXSYM
         MXSASM = MAX(MXSASM,NSTASM(ISYM))
  110 CONTINUE
C
C
      IF (NTEST .GE. 2)
     *   WRITE(LUPRI,'(/A,I5)') ' MXSASM in SYMSTR',MXSASM
      IF ( NTEST .GE. 5 ) THEN
         WRITE(LUPRI,'(/A/A)')
     *      '  NUMBER OF STRINGS OF EACH SYMMETRY',
     *      ' ===================================='
         CALL IWRTMA(NSTASM,1,MAXSYM,1,MAXSYM)
      IF ( NTEST .GE. 15 ) THEN
         write(lupri,'(/A,I10)') 'Symmetry of each string of',NSTRIN
         CALL IWRTMA(neword,1,nstrin,1,nstrin)
      END IF
      END IF
C** 2: OFFSET FOR EACH SYMMETRY
      ISTBAS(1) = 1
      DO 200 ISYM = 2,MAXSYM
         ISTBAS(ISYM) = ISTBAS(ISYM-1) + NSTASM(ISYM-1)
  200 CONTINUE
      IF ( NTEST .GE. 10 ) THEN
         WRITE(LUPRI,'(/A)') ' OFFSETS FOR EACH SYMMETRY'
         CALL IWRTMA(ISTBAS,1,MAXSYM,1,MAXSYM)
      END IF
C** 3 : REORDERING THE STRINGS
      DO 450 ISYM= 1,MAXSYM
        IEFF = 0
        DO 400 ISTRIN = 1,NSTRIN
        IF(NEWORD(ISTRIN) .EQ. ISYM ) THEN
          IEFF = IEFF + 1
          NEWORD(ISTRIN) = -(ISTBAS(ISYM)-1 + IEFF )
          CALL ICOPVE(STRING(1,ISTRIN),STRIN2(1,-NEWORD(ISTRIN)),NEL)
        END IF
  400   CONTINUE
  450 CONTINUE
      CALL ICOPVE(STRIN2(1,1),STRING(1,1),NSTRIN*NEL)
      CALL ISCAVE(NEWORD,-1,NSTRIN)
C** TWO POINTERS ARE AVAILABLE NOW :
C     NEWORH : FROM HANDY NUMBER TO REDUCED NUMBER
C     NEWORD : FROM REDUCED NUMBER TO SYMMETRY NUMBER
C   MERGE THESE SO
C     NEWORH : FROM HANDY NUMBER TO SYMMETRY NUMBER
C
      DO 500 ISTIN = 1,MAXSTR
        IF( NEWORH(ISTIN) .NE. 0 ) NEWORH(ISTIN) = NEWORD(NEWORH(ISTIN))
  500 CONTINUE
C
      IF ( NTEST .GE. 10 ) THEN
        WRITE(LUPRI,'(/A/A)') '  OLD NUMBER TO NEW NUMBER',
     &                        ' =========================='
        MAXPR = 60
        NPR = MIN(NSTRIN,MAXPR)
        CALL IWRTMA(NEWORH,1,NPR,1,MAXSTR)
        WRITE(LUPRI,'(/A/A)') '  REORDERED STRINGS ',
     &                        '  STRING          OCCUPATION '
        DO 1234 ISTRIN = 1,NPR
          WRITE(LUPRI,'(I6,8X,(20I3))')
     &    ISTRIN,(STRING(IEL,ISTRIN),IEL = 1,NEL)
 1234   CONTINUE
        IF (NSTRIN .GT. NPR) WRITE (LUPRI,*)
     &    (NSTRIN-NPR),' STRINGS NOT PRINTED.'
      END IF
C
      RETURN
      END
C  /* Deck cioff2 */
      SUBROUTINE CIOFF2(IREFSM,MAXSYM,SYMPRO,IOCOC,NSASOA,NSASOB,
     &                  NOCTPA,NOCTPB,IOFF,IOOS,NDTAS,NTEST)
C
C GENERATE OFFSETS FOR CI VECTOR OF RAS TYPE
C
C IN OUTPUT IOFF(IOCTYP,ISYM) INDICATES WHERE ALPHA STRINGS WITH
C TYPE IOCTYP AND SYMMETRY ISYM STARTS.
C
C THE CI VECTOR IS ASSUMED ORDERED AS FOLLOWS
C         LOOP OVER SYMMETRY OF ALPHA STRING
C         ( DEFINES SYMMETRY OF BETA STRING )
C             LOOP OVER TYPES OF ALPHA STRING
C               LOOP OVER ALLOWED TYPES OF BETASTRINGS
C
C                  LOOP OVER BETA STRINGS OF GIVEN TYPA AND SYM
C                    LOOP OVER ALPHA STRINGS OF GIVEN TYPE AND SYM
C
C                    END OF LOOP OVER ALPHA STRINGS
C                  END OF LOOP OVER BETA STRINGS
C               END OF LOOP OVER TYPES OF ALPHA STRINGS
C             END OF LOOP OVER TYPES OF ALPHA STRINGS
C         END OF LOOP OVER SYMMETRY OF ALPHA STRINGS
C
C
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION NSASOA(NOCTPA,MAXSYM),NSASOB( NOCTPB,MAXSYM)
      DIMENSION IOFF(NOCTPA,MAXSYM)
      DIMENSION IOOS(NOCTPB,NOCTPA,MAXSYM), NDTAS(MAXSYM)
      INTEGER   SYMPRO(8,8)
#include <priunit.h>
C
      IOFF(1,1) = 1
      DO 100 IASYM = 1,MAXSYM
        NDETSM = 0
        DO  95 IAOCC = 1, NOCTPA
          IASTRI = NSASOA(IAOCC,IASYM)
          IBSYM = SYMPRO(IASYM,IREFSM)
          IBSTR = 0
          IOOS(1,IAOCC,IASYM) = IOFF(IAOCC,IASYM)
          DO 90 IBOCC = 1, NOCTPB
            IF( IOCOC(IAOCC,IBOCC) .EQ. 1 ) THEN
              IBSTR = IBSTR + NSASOB(IBOCC,IBSYM)
              NDETSM = NDETSM + IASTRI*NSASOB(IBOCC,IBSYM)
            END IF
            IF( IBOCC .NE. NOCTPB ) IOOS(IBOCC+1,IAOCC,IASYM) =
     &      IOOS(1,IAOCC,IASYM) + IASTRI*IBSTR
   90     CONTINUE
C
          IF ( IAOCC .NE. NOCTPA) THEN
             IOFF(IAOCC+1,IASYM) = IOFF(IAOCC,IASYM) + IASTRI*IBSTR
          ELSE
             IF( IASYM .NE. MAXSYM )
     &       IOFF(1,IASYM+1) = IOFF(NOCTPA,IASYM) + IASTRI*IBSTR
          END IF
   95   CONTINUE
        NDTAS(IASYM) = NDETSM
  100 CONTINUE
C
      IF ( NTEST .GE. 10 ) THEN
         WRITE(LUPRI,'(/A/)') ' >>> OUTPUT FROM CIOFF2 (ras) <<<'
         WRITE(LUPRI,*) ' Determinant CI vector of symmetry :',IREFSM
         WRITE(LUPRI,'(/A)') ' Offset array for type and symmetry'
         CALL IWRTMA(IOFF,NOCTPA,MAXSYM,NOCTPA,MAXSYM)
         WRITE(LUPRI,'(/A)')' Size of each symmetry block in CI vector:'
         CALL IWRTMA(NDTAS,1,MAXSYM,1,MAXSYM)
      END IF
C
      RETURN
      END
C  /* Deck symexc */
      SUBROUTINE SYMEXC(NXFST,NSTRIN,TTO,TIJ,TSYM,
     *                  MAXSYM,DISSYM,SYMPRO,NTEST,NEXCI)
C
C REORDER SINGLE EXCITATIONS IN TTO AND TIJ ACCORDING TO SYMMETRY
C
      INTEGER TIJ,TTO,TSYM,DISSYM,SYMPRO
      DIMENSION TTO(NEXCI),TIJ(NEXCI)
      DIMENSION TSYM(MAXSYM+1,NSTRIN)
      DIMENSION SYMPRO(8,8),DISSYM(*)
      PARAMETER (MXEXCI = 5000)
      DIMENSION I1(MXEXCI),I2(MXEXCI)
      DIMENSION NXFST(NSTRIN)
#include <priunit.h>
C
      CALL ISETVC(TSYM,0,(NSTRIN*(MAXSYM+1)))
      IOFF = 1
      DO 1000 ISTRIN =1 ,NSTRIN
        TSYM(1,ISTRIN) = IOFF
        LEXCI = NXFST(ISTRIN)
        IF (LEXCI .GT. MXEXCI) THEN
           WRITE (LUPRI,*) ' Local arrays I1 and I2 too small in SYMEXC'
           WRITE (LUPRI,*) ' Need',LEXCI,' have',MXEXCI
           CALL QUIT('SYMEXC : insufficient dimension of local arrays')
        END IF
        IEFF = 0
        DO  900 ISYM   = 1,MAXSYM
          DO 800 I = 1,LEXCI
            IF(DISSYM(TIJ(IOFF-1+I)) .EQ. ISYM ) THEN
             IEFF = IEFF + 1
             I1(IEFF) = TIJ(IOFF-1+I)
             I2(IEFF) = TTO(IOFF-1+I)
            END IF
  800     CONTINUE
          TSYM(ISYM+1,ISTRIN) = IEFF +TSYM(1,ISTRIN)
  900   CONTINUE
        CALL ICOPVE(I1,TIJ(IOFF),LEXCI)
        CALL ICOPVE(I2,TTO(IOFF),LEXCI)
        IOFF = IOFF + LEXCI
 1000 CONTINUE
C
      IF ( NTEST .GE. 10 ) THEN
        MAXPR = 60
        NPR = MIN(MAXPR,NSTRIN)
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' SYMMETRY-ORDERED EXCITATIONS '
        WRITE(LUPRI,*) '============================='
        NPR = MIN(NSTRIN,MAXPR)
        WRITE(LUPRI,*) ' TIJ ARRAY ( ALLOWED E(IJ) )'
        IOFF = 1
        DO 1234 ISTRIN = 1,NPR
           LEXCI = NXFST(ISTRIN)
           WRITE(LUPRI,'(A,I4,A,(T29,10I5))')
     &       ' STRING..',ISTRIN,' EXCITATIONS : ',
     &       ( TIJ(IOFF-1+IEX),IEX = 1, LEXCI)
           IOFF = IOFF + LEXCI
 1234    CONTINUE
         WRITE(LUPRI,*) ' TTO ARRAY ( EXCITED STRINGS ) '
         IOFF = 1
         DO 1235 ISTRIN = 1, NPR
           LEXCI = NXFST(ISTRIN)
           WRITE(LUPRI,'(A,I4,A,(T29,10I5))')
     &       ' STRING..',ISTRIN,' new string..',
     &       ( TTO(IOFF-1+IEX),IEX = 1, LEXCI)
           IOFF = IOFF + LEXCI
 1235    CONTINUE
      WRITE(LUPRI,*) ' TSYM array '
      WRITE(LUPRI,*) ' ========== '
      CALL IWRTMA(TSYM,MAXSYM+1,NPR,MAXSYM+1,NPR)
      END IF
C
      RETURN
      END
C  /* Deck detinf */
      SUBROUTINE DETINF(NEL,NORB,SYMPRO,MULTSX,
     *                  NCDET,NHCDET,ICSYM,IHCSYM,
     *                  NRAS1,NRAS2,NRAS3,NTEST,XNDXCI,WORK,LFREE,LUSTR,
     *                  NOSYM,ICSF,REFSPC)
C
C       GENERATE INFORMATION ABOUT STRINGS  AND DETERMINANTS
C** . APRIL 27 1987 / J.O.
C     JUNE 23 .87 : KIHOOS AND KICOOS ADDED TO CALL OF CIOFF2
C     JAN  20 .88 : ALIGNED WITH SIRIUS PACKAGE
C     Feb  89     : some changes allowing csfstuff
C     Oct  90     : New weight , CSF, ...  routines
C                   RAS always assumed
C
C MOTECC-90: This module, DETINF, and the algorithms used are
C            described in Chapter 8 Setions D.1 and D.2 of MOTECC-90
C            "RAS-CI Expansions in a CSF basis" and "Slater
C            Determinants and Strings"
C
#include <implicit.h>
      INTEGER   SYMPRO
      DIMENSION SYMPRO(8,8), NRAS1(8), NRAS2(8), NRAS3(8)
      DIMENSION XNDXCI(*), WORK(LFREE)
      LOGICAL   REFSPC
#include <priunit.h>
#include <iratdef.h>
C
#include <mxpdim.h>
#include <maxash.h>
#include <detbas.h>
#include <ciinfo.h>
#include <strnum.h>
#include <mxblk.h>
#include <spinfo.h>
#include <csfbas.h>
C
C
      IF ( NTEST .GT. 0 ) THEN
         WRITE (LUPRI,'(//A,I4,A/A/)')
     *   ' >>> Output from DETINF module at level',NTEST,' <<<',
     *   ' =============================================='
         CALL GETTIM(TTSTRT,WTSTRT)
      END IF
C
C     Check max number of active orbitals
C
      IF (NORB .GT. MAXASH) THEN
         WRITE (LUPRI,'(/A,2(/A,I5))')
     &      ' DETINF ERROR, this version of CI module cannot'//
     &      ' handle so many orbitals.',
     &      ' The total number of active orbitals:   ',NORB,
     &      ' is greater than fixed parameter MAXASH=',MAXASH
         CALL QUIT('DETINF: Too many active orbitals for this version')
      END IF
C
C     Consistency check
C
      IF (MULTS .NE. MULTSX) THEN
         WRITE (LUPRI,*) ' ERROR in DETINF'
         WRITE (LUPRI,*) '    MULTS from parameter list ',MULTSX
         WRITE (LUPRI,*) '    MULTS from common /SPINFO/',MULTS
         CALL QUIT('ERROR DETINF: MULTS in /SPINFO/ .ne. MULTS param.')
      END IF
C
C** 1.02 : LOCAL MEMORY
C
      KFREE = 1
      K0 = (NORB+IRAT-1)/IRAT
      KA = (NASTR+IRAT-1)/IRAT
      KB = (NBSTR+IRAT-1)/IRAT
      KNEED = MAX(K0,KA,KB)
C     ...890819-hjaaj: NOTE use of work in CNFORD and CSDTMT
C        is not included yet.
      IF(KNEED .GT. LFREE ) THEN
        WRITE(LUPRI,*) ' NOT ENOUGH SPACE IN DETINF '
        WRITE(LUPRI,*) ' NEEDED AND AVAILABLE ',KNEED,LFREE
        CALL QUIT(' Insufficient memory in DETINF ')
      END IF
      KFREE2 = KFREE
C
C     950623-hjaaj (bug fix for NOSYM ne 0):
C     ======================================
C     Allocate temp. space in WORK for 4 arrays which
C     after MEMDET will be copied to XNDXCI.
C     (If NOSYM .ne. 0 then MAXSYM will be 1 (.lt. NSYM) and this
C      means that allocation of KORBSM etc. from MEMDET in DETFO will
C      be moved to new place in MEMDET call below)
C
      CALL MEMADD(JORBSM,NORB,KFREE,1)
      CALL MEMADD(JLTSOB,NORB,KFREE,1)
      CALL MEMADD(JSTLOB,NORB,KFREE,1)
      CALL MEMADD(JTPFOB,NORB,KFREE,1)
C
C** 1.05 : NUMBER OF ORBITALS , SYMMETRY ETC ..
C
      CALL ZORBSM(NRAS1,NRAS2,NRAS3,8,WORK(JORBSM),NTEST,
     &            WORK(JLTSOB),WORK(JSTLOB),NORB,NOSYM)
C** 1.1 : NUMBER OF STRINGS , SINGLE EXCITATIONS ...
C
      MAXSYM = IMXSYM(WORK(JORBSM),NORB)
C     ... LARGEST SYMMETRY NUMBER OCCURING
      IF ( NTEST .GE. 5 ) WRITE(LUPRI,'(/A,I3)')
     &   ' LARGEST SYMMETRY OCCURING AMONG ACTIVE ORBITALS :',MAXSYM
      IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
      CALL STRDIM(NORB,NEL,MS2,WORK(JTPFOB),NTEST )
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in STRDIM',TEND-TSTRT
      END IF
C
C** 1.2 : DYNAMIC MEMORY ALLOCATION
C
      KFREXN = 1
      CALL MEMDET(KFREXN,NTEST)
      CALL ICOPVE(WORK(JORBSM),XNDXCI(KORBSM),NORB)
      CALL ICOPVE(WORK(JLTSOB),XNDXCI(KLTSOB),NORB)
      CALL ICOPVE(WORK(JSTLOB),XNDXCI(KSTLOB),NORB)
      CALL ICOPVE(WORK(JTPFOB),XNDXCI(KTPFOB),NORB)
C
C*** 2  : STRINGS, SINGLE EXCITATIONS UND SO WEITER FOR ALPHA ELECTRONS
C
      IF (NTEST .GE. 5) THEN
         WRITE (LUPRI,'(/A)') ' >>> ALPHA  STRINGS'
         CALL GETTIM(TSTRT,WSTRT)
         TASTRT = TSTRT
      END IF
C
C** 2.1 : SET UP ZA MATRIX ACCORDING TO HANDY'S ALGO
C
      CALL WEIGHT(XNDXCI(KZA),1,WORK(KFREE),NTEST)
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in WEIGHT',TEND-TSTRT
      END IF
C
C** 2.2 :    STRINGS : ASSOCIATE NUMBER AND OCCUPATION
C
      MAXSTR = NASTR
      CALL ISTVC2(XNDXCI(KIPNSA),0,1,MAXSTR)
      IF (NTEST .GE. 11) WRITE(LUPRI,*) ' KIASTR KSTASA ',KIASTR,KSTASA
      IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
      NASTR1 = NASTR
C          DETGN4(NEL,NORB1,NORB2,NORB3,
C    &            NELMN1,NELMX1,NELMN3,NELMX3,Z,IREORD,
C    &            MAXSTR,ISTR,IOC,NORB,IPRINT)
      CALL DETGN4(NAEL,NORB1,NORB2,NORB3,
     &     NL1MNA,NEL1MX,NL3MNA,NEL3MX,XNDXCI(KZA),XNDXCI(KIPNSA),
     &     MAXSTR,XNDXCI(KIASTR),WORK(KFREE),NORB,NTEST)


      IF (NASTR .NE. NASTR1) THEN
         WRITE (LUPRI,*) 'FATAL ERROR in DETINF:'
         WRITE (LUPRI,*) 'NASTR from NUMST2 =',NASTR1
         WRITE (LUPRI,*) 'NASTR from DETGN4 =',NASTR
         CALL QUIT('ERROR DETINF: NUMST2:NASTR .ne. DETGN4:NASTR')
      END IF
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in DETGN4',TEND-TSTRT
      END IF
C
C** 2.3 :     REORDER STRINGS ACCORDING TO SYMMETRY
C
C* A POINTER POINTING FROM REDUCED NUMBER TO SYMMETRY NUMBER
      CALL MEMADD(KEWORS,NASTR,KFREE,1)
      IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
      MXSASM = 0
      CALL SYMSTR(XNDXCI(KIASTR),XNDXCI(KTATO),NORB,NAEL,
     &            WORK(KEWORS),NASTR,SYMPRO,XNDXCI(KORBSM),
     &            XNDXCI(KSTASA),XNDXCI(KSTBAA),XNDXCI(KIPNSA),
     &            MAXSTR,MAXSYM,NTEST,MXSASM)
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in SYMSTR',TEND-TSTRT
      END IF
c
c. 2.35 : ORDER IN GIVEN SYMMETRY AFTER OCCUPATION CLASS
      IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
      CALL OCCSTR(XNDXCI(KIASTR) ,XNDXCI(KTATO),NORB,NAEL,WORK(KEWORS),
     &            NASTR,SYMPRO,XNDXCI(KORBSM),XNDXCI(KSTASA),
     &            XNDXCI(KSTBAA),XNDXCI(KIPNSA),MAXSTR,
     &            NOCTPA,XNDXCI(KNSSOA),XNDXCI(KISSOA),1,MAXSYM,NTEST)
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in OCCSTR',TEND-TSTRT
      END IF
c. OCCUPATION TYPE OF STRINGS AFTER CALL TO OCCSTR
      CALL TYPSTR(XNDXCI(KIASTR),NAEL,NASTR,XNDXCI(KTPFSA),1,NTEST)
C
C** 2.4 : TABLES OF SINGLE EXCITATIONS FROM ALPHA AND BETA STRINGS
C
C SYMMETRY OF DISTRIBUTIONS
C
      IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
      CALL SETDIS(SYMPRO,XNDXCI(KISSYM),XNDXCI(KORBSM),NORB,MAXSYM,
     &            XNDXCI(KKLTP),XNDXCI(KICREA),
     &            XNDXCI(KIANNI),XNDXCI(KKLCAN) )
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in SETDIS',TEND-TSTRT
         TSTRT = TEND
      END IF
c
      CALL SEXCI(XNDXCI(KIASTR),NASTR,NAEL,NORB,XNDXCI(KTAIJ),
     &           XNDXCI(KTATO),NORB,NAEL,XNDXCI(KZA),NAEXCI,
     &           XNDXCI(KIPNSA),XNDXCI(KNXFSA) ,1,NTEST)
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in SEXCI ',TEND-TSTRT
         TSTRT = TEND
      END IF
C
C** 2.5 :  REORDER SINGLE EXCITATIONS ACCORDING TO SYMMETRY
C
C     CALL SYMEXC(NAEXCI,NASTR,XNDXCI(KTATO),XNDXCI(KTAIJ),
C    &            XNDXCI(KTASYM),MAXSYM,XNDXCI(KISSYM),SYMPRO,NTEST)
      CALL SYMEXC(XNDXCI(KNXFSA),NASTR,XNDXCI(KTATO),XNDXCI(KTAIJ),
     &            XNDXCI(KTASYM),MAXSYM,XNDXCI(KISSYM),SYMPRO,NTEST,
     &            NAEXCI)
      IF (NTEST .GE. 6) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)') ' Time used in SYMEXC',TEND-TSTRT
         TSTRT = TEND
      END IF
C
C**2.6 : REORDER SINGLE EXCITATIONS ACCORDING TO TYPE OF EXCITED STRING
C
       CALL OCCEXC(NAEXCI,NASTR,XNDXCI(KTATO),XNDXCI(KTAIJ),
     &             XNDXCI(KTASYM),MAXSYM,XNDXCI(KISSYM),SYMPRO,
     &             XNDXCI(KNXFSA),1,NOCTPA,
     &             XNDXCI(KIASTR), NAEL, NTEST)
       IF (NTEST .GE. 6) THEN
          CALL GETTIM(TEND,WEND)
          WRITE (LUPRI,'(/A,F10.3)') ' Time used in OCCEXC',TEND-TSTRT
       END IF
C
      IF (NTEST .GE. 5) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(/A,F10.3)')
     *      ' Total CPU time used for alpha strings',TEND-TASTRT
      END IF
C
C
C
C
      IF (.NOT.EQUAL) THEN
        KFREE = KFREE2
C
C*** 3 : STRINGS, SINGLE EXCITATIONS FOR BETA ELECTRONS
C
        IF (NTEST .GE. 5) THEN
           WRITE (LUPRI,'(/A)') ' >>> BETA STRINGS'
           CALL GETTIM(TSTRT,WSTRT)
           TBSTRT = TSTRT
        END IF
C
C** 3.1 : SET UP ZB MATRIX ACCORDING TO HANDY'S ALGORITHM
C
        IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
        CALL WEIGHT(XNDXCI(KZB),2,WORK(KFREE),NTEST)
        IF (NTEST .GE. 6) THEN
           CALL GETTIM(TEND,WEND)
           WRITE (LUPRI,'(/A,F10.3)') ' Time used in WEIGHT',TEND-TSTRT
        END IF
C
C** 3.2 :    STRINGS : ASSOCIATE NUMBER AND OCCUPATION
C
        MAXSTR = NBSTR
        CALL ISTVC2(XNDXCI(KIPNSB),0,1,MAXSTR)
        IF (NTEST .GE. 6) CALL GETTIM(TSTRT,WSTRT)
        NBSTR1 = NBSTR
      CALL DETGN4(NBEL,NORB1,NORB2,NORB3,
     &     NL1MNB,NEL1MX,NL3MNB,NEL3MX,XNDXCI(KZB),XNDXCI(KIPNSB),
     &     MAXSTR,XNDXCI(KIBSTR),WORK(KFREE),NORB,NTEST)
        IF (NBSTR .NE. NBSTR1) THEN
           WRITE (LUPRI,*) 'FATAL ERROR in DETINF:'
           WRITE (LUPRI,*) 'NBSTR from NUMST2 =',NBSTR1
           WRITE (LUPRI,*) 'NBSTR from DETGN4 =',NBSTR
           CALL QUIT('ERROR DETINF, NUMST2:NBSTR .ne. DETGN4:NBSTR')
        END IF
        IF (NTEST .GE. 6) THEN
           CALL GETTIM(TEND,WEND)
           WRITE (LUPRI,'(/A,F10.3)') ' Time used in DETGN4',TEND-TSTRT
           TSTRT = TEND
        END IF
C
C** 3.3 :     REORDER STRINGS ACCORDING TO SYMMETRY
C
C* A POINTER POINTING FROM REDUCED NUMBER TO SYMMETRY NUMBER
        CALL MEMADD(KEWORS,NBSTR,KFREE,1)
        CALL SYMSTR(XNDXCI(KIBSTR),XNDXCI(KTBTO),NORB,NBEL,
     &            WORK(KEWORS),NBSTR,SYMPRO,XNDXCI(KORBSM),
     &            XNDXCI(KSTASB),XNDXCI(KSTBAB),XNDXCI(KIPNSB),
     &            MAXSTR,MAXSYM,NTEST,MXSASM)
        IF (NTEST .GE. 6) THEN
           CALL GETTIM(TEND,WEND)
           WRITE (LUPRI,'(/A,F10.3)') ' Time used in SYMSTR',TEND-TSTRT
           TSTRT = TEND
        END IF
C
C** 3.35 : ORDER IN GIVEN SYMMETRY AFTER OCCUPATION CLASS
        CALL OCCSTR(XNDXCI(KIBSTR),XNDXCI(KTBTO),NORB,NBEL,
     &              WORK(KEWORS),
     &              NBSTR,SYMPRO,XNDXCI(KORBSM),XNDXCI(KSTASB),
     &              XNDXCI(KSTBAB),XNDXCI(KIPNSB),MAXSTR,
     &              NOCTPB,XNDXCI(KNSSOB),XNDXCI(KISSOB),2,MAXSYM,
     &              NTEST)
        IF (NTEST .GE. 6) THEN
          CALL GETTIM(TEND,WEND)
          WRITE (LUPRI,'(/A,F10.3)') ' Time used in OCCSTR',TEND-TSTRT
          TSTRT = TEND
        END IF
C.. OCCUPATION TYPE OF STRINGS AFTER CALL TO OCCSTR
           CALL TYPSTR(XNDXCI(KIBSTR),NBEL,NBSTR,XNDXCI(KTPFSB),2,NTEST)
C
        CALL SEXCI(XNDXCI(KIBSTR),NBSTR,NBEL,NORB,XNDXCI(KTBIJ),
     &           XNDXCI(KTBTO),NORB,NBEL,XNDXCI(KZB),NBEXCI,
     &           XNDXCI(KIPNSB),XNDXCI(KNXFSB) ,2,NTEST)
        IF (NTEST .GE. 6) THEN
           CALL GETTIM(TEND,WEND)
           WRITE (LUPRI,'(/A,F10.3)') ' Time used in SEXCI ',TEND-TSTRT
           TSTRT = TEND
        END IF
C
C** 3.5 :  REORDER SINGLE EXCITATIONS ACCORDING TO SYMMETRY
C
C         CALL SYMEXC(NBEXCI,NBSTR,XNDXCI(KTBTO),XNDXCI(KTBIJ),
C    &               XNDXCI(KTBSYM),MAXSYM,XNDXCI(KISSYM),SYMPRO,NTEST)
          CALL SYMEXC(XNDXCI(KNXFSB),NBSTR,XNDXCI(KTBTO),XNDXCI(KTBIJ),
     &               XNDXCI(KTBSYM),MAXSYM,XNDXCI(KISSYM),SYMPRO,NTEST,
     &               NBEXCI)
         IF (NTEST .GE. 6) THEN
            CALL GETTIM(TEND,WEND)
            WRITE (LUPRI,'(/A,F10.3)') ' Time used in SYMEXC',TEND-TSTRT
         END IF
C
C** 3.6 : REORDER SINGLE EXCITATIONS ACCORDING TO TYPE OF EXCITED STRING
C
         CALL OCCEXC(NBEXCI,NBSTR,XNDXCI(KTBTO),XNDXCI(KTBIJ),
     &               XNDXCI(KTBSYM),MAXSYM,XNDXCI(KISSYM),SYMPRO,
     &               XNDXCI(KNXFSB),2,NOCTPB,
     &               XNDXCI(KIBSTR), NBEL, NTEST)
         IF (NTEST .GE. 6) THEN
           CALL GETTIM(TEND,WEND)
           WRITE (LUPRI,'(/A,F10.3)') ' Time used in OCCEXC',TEND-TSTRT
         END IF
c
         IF (NTEST .GE. 5) THEN
            CALL GETTIM(TEND,WEND)
            WRITE (LUPRI,'(/A,F10.3)')
     *      ' Total CPU time used for beta strings',TEND-TBSTRT
         END IF
c
      END IF
C
C*** 3 :       OFFSETS FOR CI AND SIGMA VECTORS
C
      IF ( NTEST .GE. 2 ) CALL GETTIM(T3STRT,W3STRT)
      CALL STRPA2(XNDXCI(KNSSOA),XNDXCI(KNSSOB),NASTR,NBSTR,
     &            MAXSYM,NOCTPA,NOCTPB,XNDXCI(KIOCOC),
     &            NL1MNA,NL1MNB,NEL1MN,NEL1MX,
     &            NEL3MN,NEL3MX,NORB1,NORB2,NORB3,
     &            IOCPTA,IOCPTB,SYMPRO,NAEL,NBEL,NEL,
     &            REFSPC,NDTASM,MXNDT ,MXVBLK,NTEST)
c
      NCDET  = NDTASM(ICSYM)
      NHCDET = NDTASM(IHCSYM)
      IF (NTEST .GT. 0) THEN
         WRITE(LUPRI,'(/A,I3/A,I10)')
     *      ' REFERENCE SYMMETRY      ',ICSYM,
     *      ' NUMBER OF DETERMINANTS GENERATED',NCDET
         WRITE(LUPRI,'(/A,I3/A,I10)')
     *      ' SIGMA VECTOR OF SYMMETRY',IHCSYM,
     *      ' NUMBER OF DETERMINANTS GENERATED',NHCDET
      END IF
C
      IF ( NTEST .GE. 2 ) CALL GETTIM(T4STRT,W4STRT)
      IF ( ICSF .NE. 0 ) THEN
         IF (ICSF .GT. 0 .AND. ICSYM .NE. IHCSYM) THEN
            WRITE (LUPRI,*)
     *      'DETINF: Sorry, CSF only implemented for ICSYM = IHCSYM'
            CALL QUIT('DETINF: CSF only implemented for ICSYM = IHCSYM')
         END IF
         IF (ICSF .EQ. -2) THEN
            JCSYM = IHCSYM
         ELSE
            JCSYM = ICSYM
         END IF
C
C.. Generate prototype determinants,upper determinants,and CSF-DET
C   matrix
C
        CALL CSDTMT(XNDXCI(KDFTP),XNDXCI(KCFTP),XNDXCI(KDTOC),
     &              WORK,NTEST)
C configurations and reordering arrays
        CALL CNFORD(XNDXCI(KICTS(1)),
     &              XNDXCI(KICONF(1)),WORK,XNDXCI(KORBSM),
     &              SYMPRO,JCSYM,NORB,XNDXCI(KDFTP),
     &              XNDXCI,NCNFTP,NEL,NTEST)
C        CNFORD(ICTSDT,ICONF,IWORK,ORBSYM,
C               SYMPRO,IREFSM,NORB,IDFTP,XNDXCI,NCNFTP,NEL,NTEST)
      END IF
      IF(LUSTR .GT. 0 ) THEN
        WRITE(LUPRI,'(/A)')
     &     'String information written on SIRIUS.STRINGINFO'//
     &     ' for use by WESTA program'
C
C WRITE STRING INFORMATION ON FILE LUSTR
C FOR USE IN CONNECTION WITH WESTA
C
        REWIND(LUSTR)
        READ(LUSTR)
        WRITE(LUSTR) NCDET,NHCDET
        WRITE(LUSTR) MAXSYM,2
        WRITE(LUSTR) NAEL,NBEL
        WRITE(LUSTR) NOCTPA,NOCTPB
        WRITE(LUSTR) NASTR,NBSTR
        WRITE(LUSTR) (XNDXCI(KSTASA-1+I),I=1,MAXSYM)
        IF(.NOT.EQUAL)
     &  WRITE(LUSTR) (XNDXCI(KSTASB-1+I),I=1,MAXSYM)
        WRITE(LUSTR) (XNDXCI(KSTBAA-1+I),I=1,MAXSYM)
        IF ( .NOT. EQUAL )
     &  WRITE(LUSTR) (XNDXCI(KSTBAB-1+I),I=1,MAXSYM)
        WRITE(LUSTR) (XNDXCI(KIASTR-1+I),I=1,NAEL*NASTR)
        IF ( .NOT. EQUAL )
     &  WRITE(LUSTR) (XNDXCI(KIBSTR-1+I),I=1,NBEL*NBSTR)
        WRITE(LUSTR) (XNDXCI(KNSSOA-1+I),I=1,MAXSYM*NOCTPA)
        WRITE(LUSTR) (XNDXCI(KISSOA-1+I),I=1,MAXSYM*NOCTPA)
        IF( .NOT.EQUAL ) THEN
          WRITE(LUSTR) (XNDXCI(KNSSOB-1+I),I=1,MAXSYM*NOCTPB)
          WRITE(LUSTR) (XNDXCI(KISSOB-1+I),I=1,MAXSYM*NOCTPB)
        END IF
C
        WRITE(LUSTR) ( XNDXCI(KCOFF-1+I),I=1,MAXSYM)
        WRITE(LUSTR) (XNDXCI(KICOOS-1+I),I=1,NOCTPB*NOCTPA*MAXSYM)
        WRITE(LUSTR) (XNDXCI(KIOCOC-1+I),I=1,NOCTPA*NOCTPB)
      END IF
C
      IF (NTEST .GT. 0) THEN
         CALL GETTIM(TEND,WEND)
         IF (NTEST .GE. 2) THEN
            WRITE (LUPRI,'(2(/A,F10.3))')
     *      ' Total CPU  time used in STRPA*',T4STRT-T3STRT,
     *      ' Total CPU  time used in CSF*  ',TEND  -T4STRT
         END IF
         WRITE (LUPRI,'(2(/A,F10.3))')
     *      ' Total CPU  time used in DETINF',TEND-TTSTRT,
     *      ' Total wall time used in DETINF',WEND-WTSTRT
      END IF
C
      RETURN
      END
