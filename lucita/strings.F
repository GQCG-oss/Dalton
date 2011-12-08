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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GENSTR_GAS(NEL,NELMN1,NELMX1,NELMN3,NELMX3,
     &                  ISTASO,IGRP,NOCTYP,NSMST,Z,LSTASO,
     &                  IREORD,STRING,IOC,IOTYP,IPRNT)
*
* Generate strings consisting of  NEL electrons fullfilling
*   1 : Between NELMN1 AND NELMX1 electrons in the first NORB1 orbitals
*   2 : Between NELMN3 AND NELMX3 electrons in the last  NORB3 orbitals
*
* In the present version the strings are directly ordered into
* symmetry and occupation type .
*
* Jeppe Olsen Winter of 1990
*
* Special GAS version, Winter of 94 All strings of group IGRP
*
* ========
* Output :
* ========
* STRING(IEL,ISTRIN) : Occupation of strings.
* IREORD             : Reordering array going from lexical
*                      order to symmetry and occupation type order.
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Input
      DIMENSION ISTASO(NSMST,*)
      INTEGER Z(NACOB,NEL)
*.Orbinp
#include "mxpdim.inc"
#include "orbinp.inc"
*
*.Output
      INTEGER STRING(NEL,*),IREORD(*)
*.Scratch arrays
      DIMENSION IOC(*),LSTASO(NOCTYP,NSMST)
*
      CALL ISETVC(LSTASO,0,NOCTYP*NSMST)
      NTESTL = 00
      NTEST = MAX(NTESTL,IPRNT)
      IF( NTEST .GE. 10) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' GENSTR speaking '
        WRITE(6,*)  ' =============== '
      END IF
*
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD_LUCI(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE. 500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD_LUCI(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE. 500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD_LUCI(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
             IF( NTEST .GE. 500 ) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry
*                   ISYMST(STRING,NEL)
             ISYM = ISYMST(IOC,NEL)
*. Occupation type
C            ITYP = IOCTP2(IOC,NEL,IOTYP)
             ITYP = 1
*
             IF(ITYP.NE.0) THEN
               LSTASO(ITYP,ISYM) = LSTASO(ITYP,ISYM)+ 1
C                      ISTRNM(IOCC,NACTOB,NEL,Z,NEWORD,IREORD)
               LEXCI = ISTRNM(IOC,NACOB,NEL,Z,IREORD,0)
               LACTU = ISTASO(ISYM,IGRP)-1+LSTASO(ITYP,ISYM)
               IREORD(LEXCI) = LACTU
               IF(NTEST.GT.10) WRITE(6,*) ' LEXCI,LACTU',
     &         LEXCI,LACTU
               CALL ICOPVE(IOC,STRING(1,LACTU),NEL)
             END IF
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
      IF(NTEST.GE.1 ) THEN
        WRITE(6,*) ' Number of strings generated   ', NSTRIN
      END IF
      IF(NTEST.GE.10)  THEN
        IF(NTEST.GE.100) THEN
          NPR = NSTRIN
        ELSE
          NPR = MIN(NSTRIN,50)
        END IF
        WRITE(6,*) ' Strings generated '
        WRITE(6,*) ' =================='
        ISTRIN = 0
        DO 100 ISYM = 1, NSMST
        DO 100 ITYP = 1,NOCTYP
          LSTRIN = MIN(LSTASO(ITYP,ISYM),NPR-ISTRIN)
          IF(LSTRIN.GT.0) THEN
            WRITE(6,*) ' Strings of type and symmetry ',ITYP,ISYM
            DO 90 KSTRIN = 1,LSTRIN
              ISTRIN = ISTRIN + 1
              WRITE(6,'(2X,I4,8X,(10I5))')
     &        ISTRIN,(STRING(IEL,ISTRIN),IEL = 1,NEL)
   90       CONTINUE
          END IF
  100   CONTINUE
*
        WRITE(6,*) ' Array giving actual place from lexical place'
        WRITE(6,*) ' ============================================'
        CALL IWRTMA(IREORD,1,NPR,1,NPR)
      END IF
 
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADAADAS1_GAS(NK,I1,XI1S,LI1,
     &                IORB,NIORB,IAC,JORB,NJORB,JAC,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,SCLFAC,NSTRI)
*
* Obtain I1(KSTR) = +/- a+/a  IORB a+/a JORB !KSTR>
* Only orbital pairs IOB .gt. JOB are included
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First I orbital to be added
* NIORB : Number of orbitals to be added : IORB to IORB-1+NIORB
*        are used. They must all be in the same TS group
* JORB : First J orbital to be added
* LORB : Number of orbitals to be added : JORB to JORB-1+NJORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
* L.R. Jan 20, 1998
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*. Local scratch, atmost 1000 orbitals in a given TS block)
      DIMENSION ISCR(1000)
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' ADADS1_GAS speaking '
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' IORB,NIORB,IAC ', IORB,NIORB,IAC
       WRITE(6,*) ' JORB,NJORB,JAC ', JORB,NJORB,JAC
*
C      WRITE(6,*) ' Kstrings in action (el,string) '
C      WRITE(6,*) ' ==============================='
C      CALL IWRTMA(KSTR,NKEL,NKSTR,NKEL,NKSTR)
       WRITE(6,*) ' 24 elements of reorder array'
       CALL IWRTMA(IREO,1,24,1,24)
*
      END IF
*
      IORBMIN = IORB
      IORBMAX = IORB + NIORB - 1
*
      JORBMIN = JORB
      JORBMAX = JORB + NJORB - 1
*
      NIJ = NIORB*NJORB
*
      KEND = MIN(NKSTR,KMAX)
      IF(KEND.LT.NKSTR) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
*
      IF(IAC.EQ.2.AND.JAC.EQ.2) THEN
*
* ==========================
* Creation- creation mapping
* ==========================
*
        DO KKSTR = KMIN,KEND
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Occupation of string ', KKSTR
           CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
         END IF
*. Loop over electrons after which JORB can be added
         DO JEL = 0, NKEL
*
           IF(JEL.EQ.0 ) THEN
             JORB1 = JORBMIN - 1
           ELSE
             JORB1 = MAX(JORBMIN-1,KSTR(JEL,KKSTR))
           END IF
           IF(JEL.EQ.NKEL) THEN
             JORB2 = JORBMAX + 1
           ELSE
             JORB2 = MIN(JORBMAX+1,KSTR(JEL+1,KKSTR))
           END IF
           IF(NTEST.GE.1000)
     &      WRITE(6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
*
           IF(JEL.GT.0.AND.JORB1.GE.JORBMIN.AND.
     &                     JORB1.LE.JORBMAX) THEN
*. vanishing for any IORB
             IJOFF = (JORB1-JORBMIN)*NIORB
             DO IIORB = 1, NIORB
               IJ = IJOFF + IIORB
               I1(KKSTR-KMIN+1,IJ) = 0
               XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
             END DO
           END IF
*
           IF(JORB1.LT.JORBMAX.AND.JORB2.GT.JORBMIN) THEN
*. Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
             SIGNJ = (-1) ** JEL * SCLFAC
*. reverse lexical number of the first JEL ELECTRONS
             ILEX0 = 1
             DO JJEL = 1, JEL
               ILEX0 = ILEX0 + IZ(KSTR(JJEL,KKSTR),JJEL)
             END DO
             DO JJORB = JORB1+1, JORB2-1
* And electron JEL + 1
               ILEX1 = ILEX0 + IZ(JJORB,JEL+1)
*. Add electron IORB
               DO IEL = JEL, NKEL
                 IF(IEL.EQ.JEL) THEN
                   IORB1 = MAX(JJORB,IORBMIN-1)
                 ELSE
                   IORB1 = MAX(IORBMIN-1,KSTR(IEL,KKSTR))
                 END IF
                   IF(IEL.EQ.NKEL) THEN
                   IORB2 = IORBMAX+1
                 ELSE
                   IORB2 = MIN(IORBMAX+1,KSTR(IEL+1,KKSTR))
                 END IF
                 IF(NTEST.GE.5000)
     &            WRITE(6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
                 IF(IEL.GT.JEL.AND.IORB1.GE.IORBMIN.AND.
     &                             IORB1.LE.IORBMAX) THEN
                   IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
                   I1(KKSTR-KMIN+1,IJ) = 0
                   XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
                 END IF
                 IF(IORB1.LT.IORBMAX.AND.IORB2.GT.IORBMIN) THEN
*. Orbitals IORB1+1 - IORB2 -1 can be added after ELECTRON IEL in KSTR
*. Reverse lexical number of the first IEL+1 electrons
                   ILEX2 = ILEX1
                   DO IIEL = JEL+1,IEL
                     ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+1)
                   END DO
*. add terms for the last electrons
                   DO IIEL = IEL+1,NKEL
                     ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+2)
                   END DO
                   IJOFF = (JJORB-JORBMIN)*NIORB
                   SIGNIJ =  SIGNJ*(-1.0D0) ** (IEL+1)
                   DO IIORB = IORB1+1, IORB2-1
                     IJ = IJOFF + IIORB - IORBMIN + 1
                     ILEX = ILEX2 + IZ(IIORB,IEL+2)
                     IACT = IREO(ILEX)
                     IF(NTEST.GE.1000) THEN
                       WRITE(6,*) 'IIORB JJORB', IIORB,JJORB
                       WRITE(6,*) ' ILEX IACT ', ILEX,IACT
                     END IF
                     I1(KKSTR-KMIN+1,IJ) = IACT
                     XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
                   END DO
                 END IF
               END DO
             END DO
           END IF
         END DO
        END DO
      ELSE IF(IAC.EQ.1.AND.JAC.EQ.1) THEN
*
* ===========================================
* annihilation - annihilation mapping (a i a j)
* ===========================================
*
        DO KKSTR = KMIN,KEND
*. Active range for electrons
         IIELMIN = 0
         IIELMAX = 0
         JJELMIN = 0
         JJELMAX = 0
*
         DO KEL = 1, NKEL
          KKORB = KSTR(KEL,KKSTR)
          IF(IIELMIN.EQ.0.AND.KKORB.GE.IORBMIN)IIELMIN = KEL
          IF(JJELMIN.EQ.0.AND.KKORB.GE.JORBMIN)JJELMIN = KEL
          IF(KKORB.LE.IORBMAX) IIELMAX = KEL
          IF(KKORB.LE.JORBMAX) JJELMAX = KEL
         END DO
         IF(IIELMIN.EQ.0) IIELMIN = NKEL  + 1
         IF(JJELMIN.EQ.0) JJELMIN = NKEL  + 1

         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Occupation of string ', KKSTR
           CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
         END IF
*. Loop over first electron to be removed
C        DO JEL = 1, NKEL
         DO JEL = JJELMIN,JJELMAX
           JJORB = KSTR(JEL,KKSTR)
*. Loop over second electron to be removed
C          DO IEL = JEL+1, NKEL
           DO IEL = MAX(JEL+1,IIELMIN),IIELMAX
             IIORB = KSTR(IEL,KKSTR)
             IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IEL JEL IORB JORB ',
     &        IEL,JEL,IORB,JORB
             END IF
             IF(IIORB.GE.IORBMIN.AND.IIORB.LE.IORBMAX.AND.
     &          JJORB.GE.JORBMIN.AND.JJORB.LE.JORBMAX     )THEN
               SIGN = (-1) ** (IEL+JEL-1) * SCLFAC
*. reverse lexical number of the double annihilated string
               ILEX = 1
               DO JJEL = 1, JEL-1
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL)
               END DO
               DO JJEL = JEL+1,IEL-1
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL-1)
               END DO
               DO JJEL = IEL+1, NKEL
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL-2)
               END DO
               IACT = IREO(ILEX)
               IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                 WRITE(6,*) ' IACT out of bounds, IACT =  ', IACT
                 Call Abend2(       ' IACT out of bounds ' )
               END IF
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' ILEX and IACT ', ILEX, IACT
               END IF
*
               IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
               I1(KKSTR-KMIN+1,IJ) = IACT
               XI1S(KKSTR-KMIN+1,IJ) = SIGN
             END IF
*            ^ End if orbitals are in correct range
           END DO
*          ^ End of loop over IEL
         END DO
*        ^ End of loop over JEL
        END DO
*       ^ End of loop over Kstrings
*
      ELSE IF(IAC.EQ.2.AND.JAC.EQ.1) THEN
*
* ===================================
* Creation-annihilation map a+ i a j
* ===================================
*
C       DO KKSTR = 1, NKSTR
        DO KKSTR = KMIN,KEND
*. Indicate where a given orbital i should be added in KKSTR
         IZERO = 0
         CALL ISETVC(ISCR(IORBMIN),IZERO,NIORB)
         IIEL = 1
         DO IIORB = IORBMIN,IORBMAX
 2810      CONTINUE
           IF(IIEL.LE.NKEL) THEN
             IF(IIORB.LT.KSTR(IIEL,KKSTR)) THEN
               ISCR(IIORB)=IIEL
             ELSE IF (IIORB.EQ.KSTR(IIEL,KKSTR)) THEN
               ISCR(IIORB) = - IIEL
               IIEL = IIEL+1
             ELSE IF (IIORB.GT.KSTR(IIEL,KKSTR)) THEN
               IIEL = IIEL + 1
               GOTO 2810
             END IF
           ELSE IF(IIEL.EQ.NKEL+1) THEN
              ISCR(IIORB) = IIEL
           END IF
         END DO
         IF(NTEST.GE.10000) THEN
           WRITE(6,*) ' ISCR from IORBMIN array for KKSTR = ', KKSTR
           WRITE(6,*) ' IORBMIN, NIORB', IORBMIN,NIORB
           CALL IWRTMA(ISCR(IORBMIN),1,NIORB,1,NIORB)
         END IF
         DO JEL = 1, NKEL
           JJORB = KSTR(JEL,KKSTR)
           IF(JJORB.GE.JORBMIN.AND.JJORB.LE.JORBMAX)THEN
             DO IIORB = IORBMIN,IORBMAX
               IEL = ISCR(IIORB)
C?             write(6,*) ' JEL IEL JJORB IIORB',JEL,IEL,JJORB,IIORB
               IACT = 0
               IF(IEL.GT.0.AND.IIORB.GT.JJORB) THEN
*. New string is  a+1 ... a+ jel-1 a+jel+1 ..a+iel-1 a+iiorb a+iel+1 ...
*. Lexical number of new string
                 ILEX = 1
                 DO KEL = 1, JEL-1
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 DO KEL = JEL+1, IEL-1
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL-1)
                 END DO
                 ILEX = ILEX + IZ(IIORB,IEL-1)
                 DO KEL = IEL, NKEL
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 IACT = IREO(ILEX)
                 IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                   WRITE(6,*) ' 1: IACT out of bounds, IACT =  ', IACT
                   WRITE(6,*) ' NSTRI = ', NSTRI
                   WRITE(6,*) 'IIORB,JJORB ',IIORB,JJORB
                   WRITE(6,*) ' Kstring : '
                   CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                   WRITE(6,*) ' ILEX = ', ILEX
                   WRITE(6,*) 'IZ matrix'
                   CALL IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
                   Call Abend2( ' IACT out of bounds' )
                 END IF
                 SIGN = (-1) ** (IEL+JEL-1) * SCLFAC
               ELSE IF(IEL.GT.0 .AND. IIORB.LT.JJORB) THEN
*. New string is  a+1 ... a+ iel-1 a+ iiorb a+iel+1 ..a+jel-1 a+jel+1 ...
                 ILEX = 1
                 DO KEL = 1, IEL-1
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 ILEX = ILEX + IZ(IIORB,IEL)
                 DO KEL = IEL,JEL-1
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL+1)
                 END DO
                 DO KEL = JEL + 1, NKEL
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
C?               write(6,*) ' ILEX =' , ILEX
                 IACT = IREO(ILEX)
                 IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                   WRITE(6,*) '2 IACT out of bounds, IACT =  ', IACT
                   Call Abend2(       ' IACT out of bounds ' )
                 END IF
C?               write(6,*) ' IACT = ', IACT
                 SIGN = (-1) ** (IEL+JEL  ) * SCLFAC
               ELSE IF (IEL.LT.0. AND. IIORB .EQ. JJORB) THEN
*. Diagonal excitation
                 SIGN = SCLFAC
                 IACT = KKSTR
               END IF
               IF(IACT.NE.0) THEN
                 IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
                 I1(KKSTR-KMIN+1,IJ) = IACT
                 XI1S(KKSTR-KMIN+1,IJ) = SIGN
               END IF
             END DO
*            ^ End of loop over IIORB
           END IF
*          ^ End of  active cases
         END DO
*        ^ End of loop over electrons to be annihilated
        END DO
*       ^ End of loop over Kstrings
      ELSE IF(IAC.EQ.1.AND.JAC.EQ.2) THEN
*
* ======================================
* Annihilation-creation  map a i a+ j
* ======================================
*
*. Diagonal excitations ?
        IF(IORBMIN.EQ.JORBMIN) THEN
         IDIAG = 1
        ELSE
         IDIAG = 0
        END IF
C       DO KKSTR = 1, NKSTR
        DO KKSTR = KMIN,KEND
*. Indicate where a given orbital j should be added in KKSTR
         IZERO = 0
         CALL ISETVC(ISCR(JORBMIN),IZERO,NJORB)
         JJEL = 1
         DO JJORB = JORBMIN,JORBMAX
 0803      CONTINUE
           IF(JJEL.LE.NKEL) THEN
             IF(JJORB.LT.KSTR(JJEL,KKSTR)) THEN
               ISCR(JJORB)=JJEL
             ELSE IF (JJORB.EQ.KSTR(JJEL,KKSTR)) THEN
               ISCR(JJORB) = - JJEL
               JJEL = JJEL+1
             ELSE IF (JJORB.GT.KSTR(JJEL,KKSTR)) THEN
               JJEL = JJEL + 1
               GOTO 0803
             END IF
           ELSE IF(JJEL.EQ.NKEL+1)THEN
              ISCR(JJORB) = JJEL
           END IF
         END DO
         IF(NTEST.GE.10000) THEN
           WRITE(6,*) ' ISCR from JORBMIN array for KKSTR = ', KKSTR
           CALL IWRTMA(ISCR(JORBMIN),1,NJORB,1,NJORB)
         END IF
         DO IEL = 1, NKEL+IDIAG
*. IEL = NKEL + 1 will be used for excitations a+j aj
           IF(IEL.LE.NKEL) IIORB = KSTR(IEL,KKSTR)
           IF((IIORB.GE.IORBMIN.AND.IIORB.LE.IORBMAX).OR.
     &         IEL.EQ.NKEL+1 )THEN
             DO JJORB = JORBMIN,JORBMAX
               JEL = ISCR(JJORB)
C?             WRITE(6,*) ' IEL IIORB JEL JJORB ',
C?   &                      IEL,IIORB,JEL,JJORB
               IACT = 0
               IF(IEL.LE.NKEL) THEN
                 IF(JEL.GT.0.AND.JJORB.GT.IIORB) THEN
*. New string is  a+1 ... a+ iel-1 a+iel+1 ..a+jel-1 a+jjorb a+jel+1 ...
*. Lexical number of new string
                   ILEX = 1
                   DO KEL = 1, IEL-1
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   DO KEL = IEL+1, JEL-1
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL-1)
                   END DO
                   ILEX = ILEX + IZ(JJORB,JEL-1)
                   DO KEL = JEL, NKEL
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   IACT = IREO(ILEX)
                   IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                     WRITE(6,*) '3 IACT out of bounds, IACT =  ', IACT
                     WRITE(6,*) ' ILEX,KKSTR ', ILEX, KKSTR
                     WRITE(6,*) ' occupation of KSTR '
                     CALL IWRTMA(KSTR,1,NKEL,1,NKEL)
                     WRITE(6,*) ' IEL JEL ', IEL,JEL
                     WRITE(6,*) ' IIORB,JJORB',IIORB,JJORB
                     Call Abend2(       ' IACT out of bounds ' )
                   END IF
                   SIGN = (-1) ** (IEL+JEL) * SCLFAC
                 ELSE IF(JEL.GT.0 .AND. JJORB.LT.IIORB) THEN
*. New string is  a+1 ... a+ jel-1 a+ jjorb a+jel+1 ..a+iel-1 a+iel+1 ...
                   ILEX = 1
                   DO KEL = 1, JEL-1
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   ILEX = ILEX + IZ(JJORB,JEL)
                   DO KEL = JEL,IEL-1
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL+1)
                   END DO
                   DO KEL = IEL + 1, NKEL
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   IACT = IREO(ILEX)
                   SIGN = (-1) ** (IEL+JEL-1) * SCLFAC
                   IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                     WRITE(6,*) '4 IACT out of bounds, IACT =  ', IACT
                     WRITE(6,*) ' NSTRI = ', NSTRI
                     WRITE(6,*) 'IIORB,JJORB ',IIORB,JJORB
                     WRITE(6,*) ' Kstring : '
                     CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                     WRITE(6,*) ' ILEX = ', ILEX
                     WRITE(6,*) 'IZ matrix'
                     CALL IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
                     Call Abend2(       ' IACT out of bounds ' )
                   END IF
                 END IF
                 IF(IACT.NE.0) THEN
                   IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
                   I1(KKSTR-KMIN+1,IJ) = IACT
                   XI1S(KKSTR-KMIN+1,IJ) = SIGN
                 END IF
               ELSE IF(IEL.EQ.NKEL+1.AND.JEL.GT.0) THEN
*. Diagonal excitations aja+j
                 JJ = (JJORB-JORB)*NJORB + JJORB-JORB+1
                 I1(KKSTR-KMIN+1,JJ) = KKSTR
                 XI1S(KKSTR-KMIN+1,JJ) = SCLFAC
               END IF
             END DO
*            ^ End of loop over JJORB
           END IF
*          ^ End of  active cases
         END DO
*        ^ End of loop over electrons to be annihilated
        END DO
*       ^ End of loop over Kstrings
      END IF
*.    ^ End of types of creation mappings
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADST1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          IJ = 0
          DO  JJORB = JORB,JORB+NJORB-1
            JJORBR = JJORB-JORB+1
            DO  IIORB = IORB, IORB + NIORB - 1
              IJ = IJ + 1
C?            WRITE(6,*) ' IJ = ', IJ
C?            IF(IIORB.GT.JJORB) THEN
                IIORBR = IIORB - IORB + 1
                WRITE(6,*)
     &          ' Info for orbitals (iorb,jorb) ', IIORB,JJORB
                WRITE(6,*) ' Excited strings and sign '
                CALL IWRTMA(I1(1,IJ),1,NK,1,NK)
                CALL WRTMT_LU(XI1S(1,IJ),1,NK,1,NK)
C?            END IF
            END DO
          END DO
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADAADAST_GAS(IOB,IOBSM,IOBTP,NIOB,IAC,
     &                       JOB,JOBSM,JOBTP,NJOB,JAC,
     &                       ISPGP,ISM,ITP,KMIN,KMAX,
     &                       I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                       SCLFAC)
*
* Obtain two-operator mappings
* a+/a IORB a+/a JORB !KSTR> = +/-!ISTR>
*
* Whether creation- or annihilation operators are in use depends
* upon IAC, JAC : 1=> Annihilation,
*                 2=> Creation
*
* In the form
* I1(KSTR) =  ISTR if a+/a IORB a+/a JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included (if the types are identical)
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95   ( adadst)
*               November 1997 : annihilation added
*
*
* ======
*. Input
* ======
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLOCSTR, KLREO, KLZ, KLZSCR
      INTEGER*8 KLLOC, KLLZ, KLLREO
!     for addressing of WORK
#include "mxpdim.inc"
*./BIGGY
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*./ORBINP/
#include "orbinp.inc"
#include "strinp.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/SSAVE_LUCI/NELIS(4), NSTRKS(4)
      COMMON/UMMAGUMMA/NSTRIA(4)
      INTEGER KELFGRP(MXPNGAS),KGRP(MXPNGAS)
      COMMON/COMJEP/MXACJ,MXACIJ,MXAADST
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADAADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP,IAC ', IOB,IOBSM,IOBTP,IAC
        WRITE(6,*) ' JOB,JOBSM,JOBTP,JAC ', JOB,JOBSM,JOBTP,JAC
        WRITE(6,*) ' I12, K12 ', I12, K12
        WRITE(6,*) ' IFRST,KFRST', IFRST,KFRST
      END IF
*
*
*. Internal affairs
*
      IF(I12.LE.4.AND.K12.LE.2) THEN
        KLLOC  = KLOCSTR(K12)
        KLLZ   = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*) ' ADST_GAS : Illegal value of I12 = ', I12
        Call Abend2(' ADST_GAS : Illegal value of I12  ')
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' KLLOC KLLREO',KLLOC,KLLREO
      END IF

*
*. Supergroup and symmetry of K strings
*
      CALL SYMCOM(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*) ' K1SM,KSM : ',  K1SM,KSM
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      IF(IAC.EQ.1) THEN
        IACADJ = 2
        IDELTA =-1
      ELSE IF(IAC.EQ.2) THEN
        IACADJ = 1
        IDELTA = 1
      END IF
      IF(JAC.EQ.1) THEN
        JACADJ = 2
        JDELTA =-1
      ELSE IF(JAC.EQ.2) THEN
        JACADJ = 1
        JDELTA = 1
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' IACADJ, JACADJ', IACADJ,JACADJ
       WRITE(6,*) ' IDELTA, JDELTA', IDELTA, JDELTA
      END IF
*. Occupation of K-strings
      IF(IOBTP.EQ.JOBTP) THEN
        IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA-JDELTA
        JEL = IEL
      ELSE
        IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA
        JEL = NELFSPGP(JOBTP,ISPGPABS)-JDELTA
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' IEL, JEL', IEL,JEL
*. Trivial zero ? (Nice, then mission is complete )
      ITRIVIAL = 0
      IF(IEL.LT.0.OR.JEL.LT.0.OR.
     &   IEL.GT.NOBPT(IOBTP).OR.JEL.GT.NOBPT(JOBTP)) THEN
*. No strings with this number of elecs - be happy : No work
        NK = 0
        KACT = 0
        KACGRP = 0
        IF(NTEST.GE.100) WRITE(6,*) ' Trivial zero excitations'
        ITRIVIAL = 1
C       GOTO 9999
      ELSE
*. Find group with IEL electrons in IOBTP, JEL in JOBTP
        IIGRP = 0
        DO IGRP = IBGPSTR(IOBTP),IBGPSTR(IOBTP)+NGPSTR(IOBTP)-1
          IF(NELFGP(IGRP).EQ.IEL) IIGRP = IGRP
        END DO
        JJGRP = 0
        DO JGRP = IBGPSTR(JOBTP),IBGPSTR(JOBTP)+NGPSTR(JOBTP)-1
          IF(NELFGP(JGRP).EQ.JEL) JJGRP = JGRP
        END DO
C?      WRITE(6,*) ' ADAADA : IIGRP, JJGRP', IIGRP,JJGRP
*
        IF(IIGRP.EQ.0.OR.JJGRP.EQ.0) THEN
          WRITE(6,*)' ADAADAST : cul de sac, active K groups not found'
          WRITE(6,*)' Active GAS spaces  ' ,IOBTP, JOBTP
          WRITE(6,*)' Number of electrons', IEL, JEL
          Call
     &       Abend2(' ADAADAST : cul de sac, active K groups not found')
        END IF
*
      END IF
*. Groups defining Kstrings
      IF(ITRIVIAL.NE.1) THEN
      CALL ICOPVE(ISPGPFTP(1,ISPGPABS),KGRP,NGAS)
      KGRP(IOBTP) = IIGRP
      KGRP(JOBTP) = JJGRP
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Groups in KGRP '
        CALL IWRTMA(KGRP,1,NGAS,1,NGAS)
      END IF
      END IF
*
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(WORK(KLLZ),NGAS,
     &                  NELFSPGP(1,ISPGPABS),
     &                  NOBPT,WORK(KLZSCR),NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,
     &                         WORK(KLLOC),NOCOB,
     &                         1,WORK(KLLZ),WORK(KLLREO))
        IF(NTEST.GE.1000) THEN
         write(6,*) ' Info on I strings generated '
         write(6,*) ' NSTRI = ', NSTRI
         WRITE(6,*) ' REORDER array '
         CALL IWRTMA(WORK(KLLREO),1,NSTRI,1,NSTRI)
       END IF
       NSTRIA(I12) = NSTRI
*
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' REORDER array for I STRINGS'
       CALL IWRTMA(WORK(KLLREO),1,NSTRI,1,NSTRI)
      END IF
*
      IF(ITRIVIAL.EQ.1) GOTO 9999
      NELK = NELIS(I12)
      IF(IAC.EQ.1) THEN
        NELK = NELK + 1
      ELSE
        NELK = NELK - 1
      END IF
      IF(JAC.EQ.1) THEN
        NELK = NELK + 1
      ELSE
        NELK = NELK - 1
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' NELK = ' , NELK
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS

       CALL GETSTR2_TOTSM_SPGP(KGRP,NGAS,KSM,NELK,NSTRK,
     &                        WORK(KLLOC),NOCOB,
     &                        0,ISUM,IDUM)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
       NSTRKS(K12) = NSTRK
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' K strings generated '
         WRITE(6,*) ' Reorder array after generation of K strings'
         CALL IWRTMA(WORK(KLLREO),1,NSTRI,1,NSTRI)
       END IF
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
      JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
*
      IZERO = 0
      ZERO = 0.0D0
      CALL ISETVC(I1  ,IZERO,LI1*NIOB*NJOB)
COLD  CALL SETVEC(XI1S,ZERO ,LI1*NIOB*NJOB)
*
      CALL ADAADAS1_GAS(NK,I1,XI1S,LI1,
     &          IIOB,NIOB,IAC,JJOB,NJOB,JAC,
     &          WORK(KLLOC),NELK,NSTRK,WORK(KLLREO),WORK(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,SCLFAC,NSTRIA(I12))
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Reorder array after ADAADAS1'
         CALL IWRTMA(WORK(KLLREO),1,NSTRI,1,NSTRI)
       END IF
 9999 CONTINUE
*
C     WRITE(6,*) ' Memcheck at end of ADAADAS1 '
C     CALL MEMCHK
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADADST_GAS(IOB,IOBSM,IOBTP,NIOB,
     &                      JOB,JOBSM,JOBTP,NJOB,
     &                      ISPGP,ISM,ITP,KMIN,KMAX,
     &                      I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                      SCLFAC)
*
*
* Obtain mappings
* a+IORB a+ JORB !KSTR> = +/-!ISTR>
* In the form
* I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95
*
* ======
*. Input
* ======
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLOCSTR, KLREO, KLZ, KLZSCR
      INTEGER*8 KLLOC, KLLZ, KLLREO

!     for addressing of WORK
#include "mxpdim.inc"
*./BIGGY
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*./ORBINP/
#include "orbinp.inc"
#include "strinp.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/SSAVE_LUCI/NELIS(4), NSTRKS(4)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP ', IOB,IOBSM,IOBTP
        WRITE(6,*) ' JOB,JOBSM,JOBTP ', JOB,JOBSM,JOBTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) 'Problemo, ADADST '
C?      WRITE(6,*) ' SCLFAC = ',SCLFAC
C?    END IF

*
*. Internal affairs
*
      IF(I12.LE.4.AND.K12.LE.2) THEN
        KLLOC  = KLOCSTR(K12)
        KLLZ   = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*) ' ADST_GAS : Illegal value of I12 = ', I12
        Call Abend2( ' ADST_GAS : Illegal value of I12  ' )
      END IF

*
*. Supergroup and symmetry of K strings
*
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP(ISPGPABS,1,IOBTP,1,K1SPGPABS)
      CALL NEWTYP(K1SPGPABS,1,JOBTP,1,KSPGPABS)
      CALL SYMCOM(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*)
     & ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',
     &   K1SM,K1SPGPABS,KSM,KSPGPABS
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(WORK(KLLZ),NGAS,
     &                  NELFSPGP(1,ISPGPABS),
     &                  NOBPT,WORK(KLZSCR),NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,
     &                         WORK(KLLOC),NOCOB,
     &                         1,WORK(KLLZ),WORK(KLLREO))
      END IF
      NELK = NELIS(I12) - 2
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS
       CALL GETSTR_TOTSM_SPGP(1,KSPGPABS,KSM,NELK,NSTRK,
     &                        WORK(KLLOC),NOCOB,
     &                        0,ISUM,IDUM)
       NSTRKS(K12) = NSTRK
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
      JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
      CALL ADADS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,JJOB,NJOB,
     &          WORK(KLLOC),NELK,NSTRK,WORK(KLLREO),WORK(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,SCLFAC)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADAST_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
     &                    IAC)
*
*
* Obtain creation or annihilation mapping
*
* IAC = 2 : Creation map
* a+IORB !KSTR> = +/-!ISTR>
*
* IAC = 1 : Annihilation map
* a IORB !KSTR> = +/-!ISTR>
*
* for orbitals of symmetry IOBSM and type IOBTP
* and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
*
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR>
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version
*               October 96   : Improved version
*               September 97 : annihilation mappings added
*                              I groups defined by IGRP
*
*
* ======
*. Input
* ======
*
*./BIGGY
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
#include "lucinp.inc"
*. Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS)
      INTEGER MXVLI(MXPNGAS),MNVLI(MXPNGAS)
      INTEGER MXVLK(MXPNGAS),MNVLK(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8)
      DIMENSION IOFFI(LOFFI)
*
#include "comjep.inc"
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension
* (NKSTR,*), Where NKSTR is the number of K-strings of
*  correct symmetry . Nk is provided by this routine.
*
      CALL QENTER('ADAST ')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADAST_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        WRITE(6,*) '  Symmetry of supergroup : ', ISPGPSM
        WRITE(6,*) ' SCLFAC = ', SCLFAC
*
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' Annihilation mapping '
        ELSE IF(IAC.EQ.2) THEN
          WRITE(6,*) ' Creation mapping '
        ELSE
          WRITE(6,*) ' Unknown IAC parameter in ADAST ',IAC
          Call Abend2(' Unknown IAC parameter in ADAST ')
        END IF
*
      END IF
*. A few preparations
      NORBTS= NOBPTS(IOBTP,IOBSM)
      NORBT= NOBPT(IOBTP)
      IACGAS = IOBTP
*. First orbital of given GASpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
*. First orbital of given GASpace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM)
*
*====================================================
*. K strings : Supergroup, symmetry and distributions
*====================================================
*
      IF(IAC.EQ.1) THEN
       IDELTA = +1
      ELSE
       IDELTA = -1
      END IF
*. Is required mapping contained within current set of maps?
*. a:) Is active GASpace included in IGRP - must be
      IACGRP = 0
      DO JGRP = 1, NIGRP
       IF(IGSFGP(IGRP(JGRP)).EQ. IACGAS) IACGRP = JGRP
      END DO
*. Note : IACGRP is not the actual active group, it is the address of the
*         active group in IGRP
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ADAST in problems '
        WRITE(6,*) ' Active GASpace not included in IGRP '
        WRITE(6,*) ' Active GASpace : ', IACGAS
        Call Abend2(' ADAST : Active GASpace not included in IGRP ')
      END IF
*. b:) active group in K strings
      NIEL = NELFGP(IGRP(IACGRP))
      NKEL = NIEL + IDELTA
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
      IF(NKEL.EQ.-1.OR.NKEL.EQ.NOBPT(IACGAS)+1) THEN
*. No strings with this number of elecs - be happy : No work
        NKSTR = 0
        KACT = 0
        KACGRP = 0
        GOTO 9999
      ELSE
*. Find group with NKEL electrons in IACGAS
        KACGRP = 0
        DO JGRP = IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
          IF(NELFGP(JGRP).EQ.NKEL) KACGRP = JGRP
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
*. KACGRP is the Active group itself
        IF(KACGRP.EQ.0) THEN
          WRITE(6,*)' ADAST : cul de sac, active K group not found'
          WRITE(6,*)' GAS space and number of electrons ',
     &               IACGAS,NKEL
          Call Abend2(' ADAST : cul de sac, active K group not found')
        END IF
      END IF
*. Okay active K group was found and is nontrivial
      CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
*. The K supergroup
      CALL ICOPVE(IGRP,KGRP,NIGRP)
      KGRP(IACGRP) = KACGRP
*. Number of strings and symmetry distributions of K strings
      CALL NST_SPGRP(NIGRP,KGRP,KSM,WORK(KNSTSGP(1)),
     &               NSMST,NKSTR,NKDIST)
C          NST_SPGRP(NGRP,IGRP,ISM_TOT,NSTSGP,NSMST,NSTRIN,NDIST)
      IF(NTEST.GE.1000) WRITE(6,*)
     & ' KSM, NKSTR : ', KSM, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*. Last active space in K strings and number of strings per group and sym
      NGASL = 1
      DO JGRP = 1, NIGRP
       IF(NELFGP(KGRP(JGRP)).GT.0) NGASL = JGRP
       CALL ICOPVE2(WORK(KNSTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &              NNSTSGP(1,JGRP))
       CALL ICOPVE2(WORK(KISTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &              IISTSGP(1,JGRP))
      END DO
C     NGASL = NIGRP
*. MIN/MAX for Kstrings
      CALL MINMAX_FOR_SYM_DIST(NIGRP,KGRP,MNVLK,MXVLK,NKDIST_TOT)
      IF(NTEST.GE.1000) THEN
        write(6,*) 'MNVLK and MXVLK '
        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
      END IF
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings
* ==============
*. Generate symmetry distributions of I strings with given symmetry
      CALL TS_SYM_PNT2(IGRP,NIGRP,MXVLI,MNVLI,ISPGPSM,
     &                 IOFFI,LOFFI)
*. Offset and dimension for active group in I strings
      CALL ICOPVE2(WORK(KISTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               IACIST)
C?    WRITE(6,*) ' IACIST for IACGRP,IGRP = ', IACGRP,IGRP(IACGRP)
C?    CALL IWRTMA(IACIST,1,NSMST,1,NSMST)
*
      CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               NACIST)
*. Number of I strings per group and sym
COLD  DO IGAS = 1, NIGRP
COLD   CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(IGAS)-1)*NSMST+1,NSMST,
COLD &              IISTSGP)
COLD  END DO
*. Last entry in IGRP with a nonvanisking number of strings
      NIGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(IGRP(JGRP)).GT.0) NIGASL = JGRP
      END DO
C?    WRITE(6,*) ' NIGASL = ', NIGASL
C     NIGASL = NIGRP

*. Number of electrons before active space
      NELB = 0
      DO JGRP = 1, IACGRP-1
        NELB = NELB + NELFGP(IGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
*
      ZERO =0.0D0
      IZERO = 0
COLD  CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
 1000 CONTINUE
*. Next distribution
        CALL NEXT_SYM_DISTR(NIGRP,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,NONEW)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution '
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
        IF(NONEW.EQ.1) GOTO 9999
        KFIRST = 0
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IACGRP)
C       CALL  SYMCOM(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        CALL  SYMCOM(3,1,IOBSM,ISMFGS(IACGRP),IACSM)
        ISMFGS(IACGRP) = IACSM
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,IOFFI,MXVLI,MNVLI)
        ISMFGS(IACGRP) = ISAVE
C?      WRITE(6,*) ' IBSTRINI ', IBSTRINI
*. Number of strings before active GAS space
        NSTB = 1
C       DO IGAS = 1, IOBTP-1
        DO IGAS = 1, IACGRP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings After active GAS space
        NSTA = 1
C       DO IGAS =  IOBTP +1, NIGRP
        DO IGAS =  IACGRP+1, NIGRP
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
C?      WRITE(6,*) ' IIAC, IACSM = ',IIAC, IACSM
*
        NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
        IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
        IF(NTEST.GE.1000) THEN
        write(6,*) ' nstb nsta niac nkac ',
     &               nstb,nsta,niac,nkac
        END IF
*. Obtain annihilation/creation mapping for all strings of this type
*. Are group mappings in expanded or compact form
        IF(IAC.EQ.1.AND.ISTAC(KACGRP,2).EQ.0) THEN
          IEC = 2
          LROW_IN = NKEL
        ELSE
          IEC = 1
C         LROW_IN = NORBTS
          LROW_IN = NORBT
        END IF
        NKACT = NSTFGP(KACGRP)
*
        MXAADST = MAX(MXAADST,NKSTR*NORBTS)
C     COMMON/COMJEP/MXACJ,MXACIJ,MXAADST
        IF(NSTA*NSTB*NIAC*NKAC.NE.0)
     &  CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                 WORK(KSTSTM(KACGRP,1)),WORK(KSTSTM(KACGRP,2)),
     &                 IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                 NKSTR,KBSTRIN,NELB,I1,XI1S,SCLFAC,IAC,
     &                 LROW_IN,IEC)
        KSTRBS = KSTRBS + NKSD
C       IF(NGASL-1.GT.0) GOTO 1000
        GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADAST_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMT_LU(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
      CALL QEXIT('ADAST ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADD_STR_GROUP(NSTADD,IOFADD,ISTADD,NSTB,NSTA,
     &                         ISTRING,IELOF,NELADD,NELTOT)
*
* Part of assembling strings in individual types to
* super group of strings
*
*. Copying strings belonging to a given type to supergroup of strings
*
* Jeppe Olsen, for once improving performance of LUCIA
*
*.Input
* =====
* NSTADD : Number of strings to be added
* IOFADD : First string to be added
* ISTADD : Strings to be added
* NSTB   : Number of strings belonging to lower gasspaces
* NSTA   : Number of strings belonging to higher gasspaces
* ISTRING: Supergroup of strings under construction
* IELOF  : Place of first electron to be added
* NELADD : Number of electrons to be added
* NELTOT : Total number of electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
C     DIMENSION ISTADD(NELADD,*)
      DIMENSION ISTADD(*)
*. Input and output
C     DIMENSION ISTRING(NELTOT,*)
      DIMENSION ISTRING(*)
*
      IF(NSTA.GT.1) THEN
        DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
          IOFFY = (IOFADD-2+IISTR)*NELADD
          IOFF1 = (IISTR-1)*NSTA + 1
          IADD2 = NSTADD*NSTA
          IOFF2 = IOFF1 - IADD2
          DO ISTB = 1, NSTB
*. Address of A(1,IISTR,ISTB)
C           IOFF2 = IOFF1 + (ISTB-1)*NSTADD*NSTA
            IOFF2 = IOFF2 + IADD2
            IOFFX = IELOF-1+(IOFF2-2)*NELTOT
            DO ISTA = 1, NSTA
              IOFFX = IOFFX + NELTOT
              DO JEL = 1, NELADD
                ISTRING(JEL+IOFFX)
     &        = ISTADD(JEL+IOFFY)
C               ISTRING(IELOF-1+JEL,IOFF2-1+ISTA)
C    &        = ISTADD(JEL,IOFADD-1+IISTR)
              END DO
            END DO
          END DO
        END DO
      ELSE IF (NSTA .EQ. 1 ) THEN
*. Address of A(1,IISTR,1)
*. A(I(after),Igas,I(before))
        DO ISTB = 1, NSTB
          IOFF0 = (ISTB-1)*NSTADD
          IOFFY  = (IOFADD-2)*NELADD
          IOFFX = IELOF-1+(IOFF0-1)*NELTOT
          DO IISTR = 1,NSTADD
*. Address of A(1,IISTR,ISTB)
C           IOFF2 = IISTR  + IOFF0
C           IOFFX  = IELOF-1+(IOFF2-1)*NELTOT
            IOFFX  = IOFFX + NELTOT
            IOFFY  = IOFFY + NELADD
            DO JEL = 1, NELADD
              ISTRING(JEL+IOFFX)
     &      = ISTADD( JEL+IOFFY)
C             ISTRING(IELOF-1+JEL+(IOFF2-1)*NELTOT)
C    &      = ISTADD(JEL+(IOFADD-1+IISTR-1)*NELADD)
C             ISTRING(IELOF-1+JEL,IOFF2)
C    &      = ISTADD(JEL,IOFADD-1+IISTR)
            END DO
          END DO
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADSTN_GAS(IOBSM,IOBTP,ISPGP,ISPGPSM,ISPGPTP,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC)
*
*
* Obtain mappings
* a+IORB !KSTR> = +/-!ISTR> for orbitals of symmetry IOBSM and type IOBTP
* and I strings belonging to supergroup ISPGP wih symmetry ISPGPSM
* and type ISPGPTP(=1=>alpha,=2=>beta)
*
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR>
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version
*               October 96   : Improved version
*
* ======
*. Input
* ======
*
*./BIGGY
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLOFFI
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
#include "lucinp.inc"
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8*8*8)
      PARAMETER(MXLNGAS=7)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension
* (NKSTR,*), Where NKSTR is the number of K-strings of
*  correct symmetry . Nk is provided by this routine.
*
CTF   DIMENSION IOFFI(LOFFI)
*     set mark for local memory
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'ADSTN ')
      call memman(KLOFFI,LOFFI,'ADDL',1,'LLOFFI')
*
      CALL QENTER('ADSTN ')
*
      IF(NGAS.GT.MXLNGAS) THEN
        WRITE(6,*) ' Ad hoc programming in ADSTN (IOFFI)'
        WRITE(6,*) ' Must be changed - or redimensioned '
        Call Abend2( 'ADST : IOFFI problem ' )
      END IF
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADSTN_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) '  ISPGP ISPGPSM ISPGPTP :  ',
     &                ISPGP,ISPGPSM,ISPGPTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Problemo : ADSTN_GAS'
C?      WRITE(6,*) ' SCLFAC .ne. 1 '
C?    END IF
*
*. Supergroup and symmetry of K strings
*
      ISPGRPABS = IBSPGPFTP(ISPGPTP)-1+ISPGP
      CALL NEWTYP(ISPGRPABS,1,IOBTP,1,KSPGRPABS)
      CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
      NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
      IF(NTEST.GE.200) WRITE(6,*)
     & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*
      NORBTS= NOBPTS(IOBTP,IOBSM)
      ZERO =0.0D0
      CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      IZERO = 0
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
*. First orbital of given GASSpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
*. First orbital of fiven GASSPace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM)
*
*. Information about I strings
* =============================
*
*. structure of group of strings defining I strings
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Number of electrons before active type
      NELB = 0
      DO IGAS = 1, IOBTP -1
        NELB = NELB + NELFGS(IGAS)
      END DO
*. Number of electrons in active space
      NACGSOB = NOBPT(IOBTP)
*. Number of strings per symmetry for each symmetry
      DO IGAS = 1, NGAS
        CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
      END DO
*. Offset and dimension for active group in I strings
      CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               IACIST)
      CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               NACIST)
C?     WRITE(6,*) ' IACIST and NACIST arrays '
C?     CALL IWRTMA(IACIST,1,NSMST,1,NSMST)
C?     CALL IWRTMA(NACIST,1,NSMST,1,NSMST)
*
*. Generate offsets for I strings with given symmetry in
*  each space
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGPSM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
*
        call i_save_r(WORK(KLOFFI),NSTRINT+1,IOFF)
CTF     IOFFI(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.200) THEN
          WRITE(6,'(A22,I8,I8,I8)')
     &    'IOFF,NSTRINT+1 NSTRII ',
     &     IOFF,IFRMR(WORK(KLOFFI),1,IOFF),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
*. Supergroup and symmetry of K strings
*
CM    CALL NEWTYP(ISPGRPABS,1,IOBTP,1,KSPGRPABS)
CM    CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
CM    NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
CM    IF(NTEST.GE.200) WRITE(6,*)
CM   & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
*
*. Gas structure of K strings
*
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,KSPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Active group of K-strings
      KACGRP = ITPFGS(IOBTP)
*. Number of strings per symmetry distribution
      DO IGAS = 1, NGAS
        CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
 1000 CONTINUE
        IF(KFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        KFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,KSM)
C?      write(6,*) ' after  SYMCOM '
C?      write(6,*) ' ngasl istsmm1 ksm',ngasl,istsmm1,ksm
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IOBTP)
        CALL  SYMCOM(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        ISMFGS(IOBTP) = IACSM
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGAS
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
        ISMFGS(IOBTP) = ISAVE
CTF     IBSTRINI = IOFFI(IOFF)
        IBSTRINI = IFRMR(WORK(KLOFFI),1,IOFF)
C?      WRITE(6,*) ' IOFF IBSTRINI ', IOFF,IBSTRINI
*. Number of strings before active GAS space
        NSTB = 1
        DO IGAS = 1, IOBTP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings before active GAS space
        NSTA = 1
        DO IGAS =  IOBTP+1, NGAS
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group
C?      write(6,*) ' IACSM = ', IACSM
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
*
        NKAC = NNSTSGP(ISMFGS(IOBTP),IOBTP)
        IKAC = IISTSGP(ISMFGS(IOBTP),IOBTP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
C?      write(6,*) ' nstb nsta niac nkac ',
C?   &               nstb,nsta,niac,nkac
*. Obtain annihilation n mapping for all strings of this type
*
        NORBTS= NOBPTS(IOBTP,IOBSM)
*
        NKACT = NSTFGP(KACGRP)
C?      write(6,*) ' KACGRP ', KACGRP
        CALL ADSTN_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                 WORK(KSTSTM(KACGRP,1)),WORK(KSTSTM(KACGRP,2)),
     &                 IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                 NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC)
        KSTRBS = KSTRBS + NKSD
        IF(NGASL-1.GT.0) GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADSTN_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMT_LU(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
*     eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'ADSTN ')

      CALL QEXIT('ADSTN ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ANNSTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NACOB,
     &                  IEC,LDIM,IPRNT)
*
* A group of strings containing NEL electrons is given
* Set up all possible ways of removing an electron
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL - 1 electrons
* NEWORD : Reordering array for N-1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL-1 strings
* IEC    : = 1 Extended map, dimension equals number of orbs
* IEC    : = 2 Compact  map, dimension equals number of elecs
* LDIM   : Row dimension ( see IEC)
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(LDIM,NSTINI),TTO(LDIM,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTESTL =  1
      NTEST = MAX(IPRNT,NTESTL)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' ANNSTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
      LUOUT = 6
*
      DO 1000 ISTRIN = 1,NSTINI
        DO 100 IEL = 1, NEL
*. String with electron removed
          DO JEL = 1, IEL-1
           STRIN2(JEL) = STRING(JEL,ISTRIN)
          END DO
          DO JEL = IEL+1, NEL
           STRIN2(JEL-1) = STRING(JEL,ISTRIN)
          END DO
          JSTRIN = ISTRNM(STRIN2,NACOB,NEL-1,Z,NEWORD,1)
C?        WRITE(6,*) ' anni-string and number '
C?        CALL IWRTMA(STRIN2,1,NEL-1,1,NEL-1)
C?        WRITE(6,*) ' JSTRIN = ', JSTRIN
*
          IORBABS = STRING(IEL,ISTRIN)
          IORB = STRING(IEL,ISTRIN)-IORBOF+1
          IF(IEC.EQ.1) THEN
            IROW = IORB
          ELSE
            IROW = IEL
          END IF
*
          TI(IROW,ISTRIN ) = -IORBABS
C         TI(IROW,ISTRIN ) = -IORB
          TTO(IROW,ISTRIN) = JSTRIN
          IIISGN = (-1)**(IEL-1)
          IF(LSGSTR.NE.0)
     &    IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
          IF(IIISGN .EQ. -1 )
     &    TTO(IROW,ISTRIN) = - TTO(IROW,ISTRIN)
  100   CONTINUE
*
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(LUOUT,*) ' Output from ANNSTR : '
        WRITE(LUOUT,*) '==================='
*
        WRITE(6,*)
        WRITE(LUOUT,*) ' Strings with an electron added or removed'
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,LDIM)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,LDIM)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE CRESTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NACOB,
     &                  IPRNT)
*
* A type of strings containing NEL electrons are given
* set up all possible ways of adding an electron to this type of strings
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL + 1 electrons
* NEWORD : Reordering array for N+1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL+1 strings
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(NORB,NSTINI),TTO(NORB,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTESTL =  1
      NTEST = MAX(IPRNT,NTESTL)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' CRESTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
C?    WRITE(6,*) ' Reorder array NEWORD '
C?    CALL IWRTMA(NEWORD,1,NSTINO,1,NSTINO)
      LUOUT = 6
*
      DO 1000 ISTRIN = 1,NSTINI
C?    write(6,*) ' Input string ',istrin,(string(i,istrin),i=1,nel)
        DO 100 IORB = IORBOF, IORBOF-1+NORB
C?      write(6,*) ' orbital ',iorb
           IPLACE = 0
           IF(NEL.EQ.0) THEN
             IPLACE = 1
             GOTO 11
           ELSE IF ( NEL .NE. 0 ) THEN
            DO 10 IEL = 1, NEL
              IF(IEL.EQ.1.AND.STRING(1,ISTRIN).GT.IORB) THEN
                IPLACE = 1
                GOTO 11
              ELSE IF( (IEL.EQ.NEL.AND.IORB.GT.STRING(IEL,ISTRIN)) .OR.
     &                 (IEL.LT.NEL.AND.IORB.GT.STRING(IEL,ISTRIN).AND.
     &                  IORB.LT.STRING(IEL+1,ISTRIN)) ) THEN
                IPLACE = IEL+1
                GOTO 11
              ELSE IF(STRING(IEL,ISTRIN).EQ.IORB) THEN
                IPLACE = 0
                GOTO 11
              END IF
   10       CONTINUE
           END IF
   11     CONTINUE
*
C?        write(6,*) ' iplace = ', iplace
          IF(IPLACE.NE.0) THEN
*. Generate next string
            DO 30 I = 1, IPLACE-1
   30       STRIN2(I) = STRING(I,ISTRIN)
            STRIN2(IPLACE) = IORB
            DO 40 I = IPLACE,NEL
   40       STRIN2(I+1) = STRING(I,ISTRIN)
C?          write(6,*) ' updated string (STRIN2) '
C?          call iwrtma(STRIN2,1,NEL+1,1,NEL+1)
            JSTRIN = ISTRNM(STRIN2,NACOB,NEL+1,Z,NEWORD,1)
C?          write(6,*) ' corresponding number ', JSTRIN
*
            TTO(IORB-IORBOF+1,ISTRIN) = JSTRIN
            IIISGN = (-1)**(IPLACE-1)
            IF(LSGSTR.NE.0)
     &      IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
            IF(IIISGN .EQ. -1 )
     &      TTO(IORB-IORBOF+1,ISTRIN) = - TTO(IORB-IORBOF+1,ISTRIN)
            TI(IORB-IORBOF+1,ISTRIN ) = IORB
          END IF
  100   CONTINUE
*
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(LUOUT,*) ' Output from CRESTR : '
        WRITE(LUOUT,*) '==================='
*
        WRITE(6,*)
        WRITE(LUOUT,*) ' Strings with an electron added  '
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,NORB)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,NORB)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
     &                              NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete
*                           orbital space, product of strings in
*                           each GAS space
*
* Compared to GETSTR2_TOTSM_SPGP : Based upon IGRP(NIGRP)
*                                  (Just a few changes in the beginning)
*
* =====
* Input
* =====
*
* IGRP :  supergroup, here as an array of GAS space
* NIGRP : Number of active groups
* ISPGRPSM : Total symmetry of superstrings
* NEL : Number of electrons
* IZ  : Reverse lexical ordering array for this supergroup (IF IDOREO.NE.0)
*
*
* ======
* Output
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0)
*
*
* Jeppe Olsen, Written  July 1995
*              Version of Dec 1997
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
      INTEGER IZ(NORBT,NEL)
      INTEGER IGRP(NIGRP)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      NTEST = 0000
      CALL QENTER('GETST')
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A)')  ' Strings to be obtained : '
        WRITE(6,'(A)') ' **************************'
        WRITE(6,'(A)')
        WRITE(6,'(A,I2)') '   Symmetry : ', ISPGRPSM
        WRITE(6,'(A,16I3)') ' Groups : ', (IGRP(I),I=1,NIGRP)
        WRITE(6,*)
      END IF
*. Absolut number of this supergroup
*. Occupation per gasspace
*. Largest occupied space
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = IGRP(IGAS)
        NELFGS(IGAS) = NELFGP(IGRP(IGAS))
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
* Largest and lowest active symmetries for each GAS space
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Obtain all strings of this symmetry
CT      CALL QENTER('GASSM')
        CALL GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1,ISTRBS),NSTR,NEL,
     &                          NNSTSGP,IISTSGP)
CT      CALL QEXIT('GASSM')
*. Reorder Info : Lexical => actual number
        IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL
              LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
            END DO
C           WRITE(6,*) ' string '
C           CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
C           WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
        END IF
*
        ISTRBS = ISTRBS + NSTR
*. ready for next symmetry distribution
        IF(NGAS-1.NE.0) GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
*
      CALL QEXIT('GETST')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETSTR_TOTSM_SPGP(ISTRTP,ISPGRP,ISPGRPSM,NEL,NSTR,ISTR,
     &                             NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete
*                           orbital space, product of strings in
*                           each GAS space
* =====
* Input
* =====
*
* ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
* ISPGRP :  supergroup number, (relative to start of this type )
* ISPGRPSM : Total symmetry of superstrings
* NEL : Number of electrons
* IZ  : Reverse lexical ordering array for this supergroup
*
*
* ======
* Output
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0)
*
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
      INTEGER IZ(NORBT,NEL)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      NTEST = 0000
      CALL QENTER('GETST')
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A,3I3)')
     & ' Strings to be obtained : Type, supergroup, symmetry ',
     &   ISTRTP,ISPGRP,ISPGRPSM
        WRITE(6,*)
      END IF
*. Absolut number of this supergroup
      ISPGRPA = IBSPGPFTP(ISTRTP) - 1 + ISPGRP
*. Occupation per gasspace
*. Largest occupied space
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPA)
        NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(WORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
* Largest and lowest active symmetries for each GAS space
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
C        NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
C        CALL NXTNUM2(ISMFGS,NGASL-1,1,MAXVAL,NONEW)
C        CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
C         SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
C?        write(6,*) ' ISTSMM1 : ', ISTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Obtain all strings of this symmetry
CT      CALL QENTER('GASSM')
        CALL GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1,ISTRBS),NSTR,NEL,
     &                          NNSTSGP,IISTSGP)
CT      CALL QEXIT('GASSM')
*. Reorder Info : Lexical => actual number
        IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL
              LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
            END DO
C           WRITE(6,*) ' string '
C           CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
C           WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
        END IF
*
        ISTRBS = ISTRBS + NSTR
*. ready for next symmetry distribution
        IF(NGAS-1.NE.0) GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
*
      CALL QEXIT('GETST')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTROC,NSTR,NEL,
     &                          NNSTSGP,IISTSGP)
*
* Obtain all superstrings containing  strings of given sym and type
*
* ( Superstring :contains electrons belonging to all gasspaces
*        string :contains electrons belonging to a given GAS space
* A super string is thus a product of NGAS strings )
*
* Jeppe Olsen, Summer of 95
*              Optimized version, october 1995
*
*. In this subroutine the ordering of strings belonging to a given type
*  is defined !!
* Currently we are using the order
* Loop over GAS 1 strings
*  Loop over GAS 2 strings
*   Loop over GAS 3 strings --
*
*     Loop over gas N strings
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
*. Specific input
      INTEGER ITPFGS(*), ISMFGS(*)
      INTEGER NNSTSGP(MXPNSMST,*), IISTSGP(MXPNSMST,*)
*. Local scratch
C     INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS), ISTRNM(MXPNGAS)
      INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS)
*. Output
      INTEGER ISTROC(NEL,*)
*. Number of strings per GAS space
C?    write(6,*) ' entering problem child '
      DO IGAS = 1, NGAS
        NSTFGS(IGAS)  = NNSTSGP(ISMFGS(IGAS),IGAS)
        IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GETSTR_GASSM_SPGP speaking '
        WRITE(6,*) '  =========================== '
        WRITE(6,*) ' ISMFGS,ITPFGS (input) '
        CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' NSTFGS, IBSTFGS ( intermediate results ) '
        CALL IWRTMA(NSTFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(IBSTFGS,1,NGAS,1,NGAS)
      END IF
*. Last gasspace with a nonvanishing number of electrons
      IGASL = 0
      DO IGAS = 1, NGAS
        IF( NELFGP(ITPFGS(IGAS)) .NE. 0 ) IGASL = IGAS
      END DO
C     WRITE(6,*) ' IGASL = ', IGASL
*
      NSTRTOT = 1
      DO IGAS = 1, NGAS
        NSTRTOT = NSTRTOT*NSTFGS(IGAS)
      END DO
C     WRITE(6,*) ' NSTRTOT = ', NSTRTOT
      IF(IGASL.EQ.0) GOTO 2810
*
      NELL = NELFGP(ITPFGS(IGASL))
      NELML = NEL - NELL
      NSTRGASL = NSTFGS(IGASL)
      IBGASL = IBSTFGS(IGASL)
*
      IF(NSTRTOT.EQ.0) GOTO 1001
*. Loop over GAS spaces
      DO IGAS = 1, IGASL
*. Number of electrons in GAS = 1, IGAS - 1
        IF(IGAS.EQ.1) THEN
          NELB = 0
        ELSE
          NELB = NELB +  NELFGP(ITPFGS(IGAS-1))
        END IF
*. Number of electron in IGAS
        NELI = NELFGP(ITPFGS(IGAS))
C?      WRITE(6,*) ' NELI and NELB ', NELI,NELB
        IF(NELI.GT.0) THEN

*. The order of strings corresponds to a matrix A(I(after),Igas,I(before))
*. where I(after) loops over strings in IGAS+1 - IGASL and
*  I(before) loop over strings in 1 - IGAS -1
          NSTA = 1
          DO JGAS = IGAS+1, IGASL
            NSTA = NSTA * NSTFGS(JGAS)
          END DO
*
          NSTB =  1
          DO JGAS = 1, IGAS-1
            NSTB = NSTB * NSTFGS(JGAS)
          END DO
*
          NSTI = NSTFGS(IGAS)

C?        write(6,*) ' before call to add_str_group '
          CALL ADD_STR_GROUP(NSTI,
     &          IBSTFGS(IGAS),
     &          WORK(KOCSTR(ITPFGS(IGAS))),
     &          NSTB,NSTA,ISTROC,NELB+1,NELI,NEL)
*. Loop over strings in IGAS
        END IF
      END DO
 1001 CONTINUE
 2810 CONTINUE
      NSTR = NSTRTOT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from  GETSTR_GASSM_SPGP '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and type strings : '
        WRITE(6,*)
        WRITE(6,*) '   AS    Sym  Type '
        WRITE(6,*) ' =================='
        DO IGAS = 1, NGAS
          WRITE(6,'(3I6)') IGAS,ISMFGS(IGAS),ITPFGS(IGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Number of strings generated : ', NSTR
        WRITE(6,*) ' Strings generated '
        CALL PRTSTR(ISTROC,NEL,NSTR)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION ISYMS1(STRING,NEL)
*
* Symmmetry of string, D2H version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "orbinp.inc"
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
*. Specific input
      INTEGER STRING(*)
*
      ISYM = 1
      DO 100 IEL = 1, NEL
        ISYM = SYMPRO(ISYM,ISMFTO(STRING(IEL)))
  100 CONTINUE
*
      ISYMS1 = ISYM
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' ISYMS1, String and symmetry '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,*) ISYM
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION ISYMS2(STRING,NEL)
*
* Symmetry of string STRING, D inf h, C inf v, O3 version
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONAB/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*
      INTEGER   STRING(NEL)
*. ML and parity of string
      MLSTR = 0
      IPARI = 1
      DO 10 IEL = 1, NEL
        IF(ISMFTO(STRING(IEL)).LE.NMLOB) THEN
          MLSTR = MLSTR + ISMFTO(STRING(IEL))-1+MNMLOB
        ELSE
          MLSTR = MLSTR + ISMFTO(STRING(IEL))-1+MNMLOB-NMLOB
          IPARI = - IPARI
        END IF
   10 CONTINUE
*
      IF(IPARI.EQ.-1) IPARI = 2
      ISYM  = (IPARI-1) * NMLST+ MLSTR - MNMLST + 1
      ISYMS2 = ISYM
*
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' STRING '
        CALL IWRTMA(STRING,1,NEL,1,NEL)
        WRITE(6,'(A,3I3)') ' MLSTR, IPARI ISYMS2 ', MLSTR,IPARI,ISYMS2
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION ISYMST(STRING,NEL)
*
* Master routine for symmetry of string
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input ( PNTGRP is used )
#include "mxpdim.inc"
#include "lucinp.inc"
*. Specific input
      INTEGER STRING(*)
      IF(PNTGRP.EQ.1) THEN
*.D2h
        ISYMST = ISYMS1(STRING,NEL)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
*.Cinfv Dinfh O3
        ISYMST = ISYMS2(STRING,NEL)
      ELSE
        WRITE(6,*) ' Sorry PNTGRP option not programmed ', PNTGRP
        WRITE(6,*) ' Enforced stop in ISYMST '
        Call Abend1( 5 )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION ISYMSTR(ISYM,NSTR)
*
* Symmetry of product of NSTR string symmetries
*
* works currently only for D2H and subgroups
*
* Jeppe Olsen, 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "multd2h.inc"
*. Input
      INTEGER ISYM(*)
*
      IF(NSTR.EQ.0) THEN
        IISYM = 1
      ELSE
        IISYM = ISYM(1)
        DO JSTR = 2, NSTR
           IISYM = MULTD2H(IISYM,ISYM(JSTR))
        END DO
      END IF
*
      ISYMSTR = IISYM
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MEMSTR_GAS(NTEST)
*
*
* Construct pointers for saving information about strings and
* their mappings
*
* GAS version
*
*========
* Input :
*========
* Number and groups of strings defined by /GASSTR/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen , Winter of 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "csm.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "stinf.inc"

*. Start of string information
      CALL MEMMAN(KSTINF,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      IF(iprint.GE.10)
     &WRITE(lupri,*) ' First word with string information',KSTINF
*
*.  Offsets for occupation and reorder array of strings
*
      DO IGRP = 1, NGRP
        NSTRIN = NSTFGP(IGRP)
        LSTRIN = NSTRIN*NELFGP(IGRP)
        CALL MEMMAN(KOCSTR(IGRP),LSTRIN,'ADDS  ',1,'OCSTR ')
        CALL MEMMAN(KSTREO(IGRP),NSTRIN,'ADDS  ',1,'STREO ')
      END DO
*
*. Number of strings per symmetry and offset for strings of given sym
*. for groups
*
      CALL MEMMAN(KNSTSGP(1),NSMST*NGRP,'ADDS  ',1,'NSTSGP')
      CALL MEMMAN(KISTSGP(1),NSMST*NGRP,'ADDS  ',1,'ISTSGP')
*
*. Number of strings per symmetry and offset for strings of given sym
*. for types
*
      DO  ITP  = 1, NSTTP
        CALL MEMMAN(KNSTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'NSTSO ')
        CALL MEMMAN(KISTSO(ITP),NSPGPFTP(ITP)*NSMST,'ADDS  ',1,'ISTSO ')
      END DO
*
**. Lexical adressing of arrays : use array indeces for complete active space
*
*. Not in use so
      DO  IGRP = 1, NGRP
        CALL MEMMAN(KZ(IGRP),NACOB*NELFGP(IGRP),'ADDS  ',1,'Zmat  ')
      END DO
*
*. Mappings between different groups
*
      DO  IGRP = 1, NGRP
        IEL = NELFGP(IGRP)
        IGAS = IGSFGP(IGRP)
        IORB = NOBPT(IGAS)
        ISTRIN = NSTFGP(IGRP)
*. IF creation is involve : Use full orbital notation
*  If only annihilation is involved, compact form will be used
        IF(ISTAC(IGRP,2).NE.0) THEN
          LENGTH = IORB*ISTRIN
          CALL MEMMAN(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          CALL MEMMAN(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
        ELSE IF(ISTAC(IGRP,1).NE.0) THEN
*. Only annihilation map so
          LENGTH = IEL*ISTRIN
          CALL MEMMAN(KSTSTM(IGRP,1),LENGTH,'ADDS  ',1,'ORBMAP')
          CALL MEMMAN(KSTSTM(IGRP,2),LENGTH,'ADDS  ',1,'STRMAP')
        ELSE
*. Neither annihilation nor creation (?!)
          KSTSTM(IGRP,1) = -1
          KSTSTM(IGRP,2) = -1
        END IF
      END DO
*
*. Symmetry of conjugated orbitals and orbital excitations
*
*     KCOBSM,KNIFSJ,KIFSJ,KIFSJO
      CALL MEMMAN(KCOBSM,NACOB,'ADDS  ',1,'Cobsm ')
      CALL MEMMAN(KNIFSJ,NACOB*NSMSX,'ADDS  ',1,'Nifsj ')
      CALL MEMMAN(KIFSJ,NACOB**2 ,'ADDS  ',1,'Ifsj  ')
      CALL MEMMAN(KIFSJO,NACOB*NSMSX,'ADDS  ',1,'Ifsjo ')
*
*. Symmetry of excitation connecting  strings of given symmetry
*
      CALL MEMMAN(KSTSTX,NSMST*NSMST,'ADDS  ',1,'Ststx ')
*
*. Occupation classes
*
      CALL MEMMAN(KIOCLS,NMXOCCLS*NGAS,'ADDS  ',1,'IOCLS ')
*. Annihilation/Creation map of supergroup types
      CALL MEMMAN(KSPGPAN,NTSPGP*NGAS,'ADDS  ',1,'SPGPAN')
      CALL MEMMAN(KSPGPCR,NTSPGP*NGAS,'ADDS  ',1,'SPGPCR')
*
*. Last word of string information
*
      CALL MEMMAN(KSTINE,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      IF(iprint.GE.10)
     &  WRITE(lupri,*) ' Last word with string information',KSTINE-1
      IF(iprint.NE.0) THEN
      WRITE(lupri,*) ' R*8 words used for storing string information ',
     &                 KSTINE - KSTINF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE NST_SPGRP(NGRP,IGRP,ISM_TOT,NSTSGP,NSMST,
     &                     NSTRIN,NDIST)
*
* Number of strings for given combination of groups and
* symmetry.
*
*. Input
*
*
*   NGRP : Number of active groups
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSTSGP  : Number of strings per symmetry and supergroup
*   NSMST   : Number of string symmetries
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific Input
      DIMENSION IGRP(NGRP)
*. General input
      DIMENSION NSTSGP(NSMST,*)
*. Scratch
#include "mxpdim.inc"
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' NST_SPGP is speaking '
        WRITE(6,*) ' ====================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NGRP)
      END IF
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST(NGRP,IGRP,MNSM,MXSM,NDISTX)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0
      NDIST = 0
 1000 CONTINUE
*. Next symmetry distribution
        CALL NEXT_SYM_DISTR(NGRP,MNSM,MXSM,ISM,ISM_TOT,IFIRST,NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = 1
          DO JGRP = 1, NGRP
            LDIST = LDIST*NSTSGP(ISM(JGRP),IGRP(JGRP))
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
      NSTRIN = LENGTH
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', LENGTH
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE NSTPTP_GAS(NGAS,ISPGRP,NSTSGP,NSMST,
     &                      NSTSSPGP,IGRP,MXNSTR,
     &                      NSMCLS,NSMCLSE,NSMCLSE1)
*
* From number of strings per group and sym to number of strings
* per supergroup and sym for given super group.
*
* Jeppe Olsen , Fall of 94
*
* Last Revision : Aug 97, NSMCLS,NSMCLSE added
*
* NSMCLS : MAX Number of symmetry classes for given supergroup,
*          i.e. number of combinations of symmetries of groups
*          containing strings
* NSMCLSE : Number of symmetry classes for given supergroup
*          obtained by restricting allowed symmetries in
*          a given group by a max and min.
* NSMCLSE1 : As NSMCLSE, but the symmetry of the last active
*            orbital space where there is more than one symmetry
*            is left out
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ISPGRP(NGAS),NSTSGP(NSMST,*)
*. Output
      DIMENSION NSTSSPGP(NSMST,IGRP)
*. Scratch
#include "mxpdim.inc"
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER MSMCLS(MXPNSMST)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NSTPTP_GAS is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Supergroup in action '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
      END IF
*
      IZERO = 0
      CALL ISETVC(NSTSSPGP(1,IGRP),IZERO,NSMST)
      CALL ISETVC(MSMCLS,IZERO,NSMST)

*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTSGP(ISYM,ISPGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. NSMCLSE
      NSMCLSE = 1
      DO IGAS = 1, NGAS
        NSMCLSE = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE
      END DO
*. NSMCLSE1
      NSMCLSE1 = 1
      DO IGAS = 1, NGASL-1
        NSMCLSE1 = (MXSM(IGAS)-MNSM(IGAS)+1)*NSMCLSE1
      END DO
*. First symmetry combination
      DO IGAS = 1, NGAS
         ISM(IGAS) = MNSM(IGAS)
      END DO
*. Loop over symmetries in each gas space
      IFIRST = 1
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
      ELSE
C       NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        ISMSPGP = ISM(1)
        NST = NSTSGP(ISM(1),ISPGRP(1))
        DO JGRP = 2, NGAS
          CALL SYMCOM(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTSGP(ISM(JGRP),ISPGRP(JGRP))
        END DO
C?      WRITE(6,*) ' Symmetry of groups '
C?      CALL IWRTMA(ISM,1,NGAS,1,NGAS)
C?      WRITE(6,*) ' Symmetry of super group and number of strings '
C?      WRITE(6,*) ISMSPGP , NST
        NSTSSPGP(ISMSPGP,IGRP) =   NSTSSPGP(ISMSPGP,IGRP) + NST
        IF(NST.NE.0) MSMCLS(ISMSPGP) = MSMCLS(ISMSPGP) + 1
        GOTO 1000
      END IF
*
      MXNSTR = -1
      NSMCLS = 0
      DO ISTRSM = 1, NSMST
        MXNSTR = MAX(MXNSTR,NSTSSPGP(ISTRSM,IGRP))
        NSMCLS = MAX(NSMCLS,MSMCLS(ISTRSM))
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     &  ' Number of strings per symmetry for supergroup',IGRP
        CALL IWRTMA10(NSTSSPGP(1,IGRP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Largest number of strings of given sym ',MXNSTR
*
        WRITE(6,'(A,3I7)') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                       NSMCLS,NSMCLSE,NSMCLSE1
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE NSTRSO_GAS(NEL,NORB1,NORB2,NORB3,
     &                  NELMN1,NELMX1,NELMN3,NELMX3,
     &                  IOC,NORB,NSTASO,ISTASO,
     &                  NOCTYP,NSMST,IOTYP,IPRNT)
*
* Number of strings per symmetry for group IOTYP
*
* Gas version, no check of type : set to 1
*
* Jeppe Olsen Winter of 1994
*
      IMPLICIT REAL*8           ( A-H,O-Z)
C     DIMENSION IOC(*),NSTASO(NOCTYP,NSMST)
      DIMENSION IOC(*),NSTASO(NSMST,*),ISTASO(NSMST,*)
*
      CALL ISETVC(NSTASO(1,IOTYP),0,NSMST)
      NTESTL = 0
      NTEST = MAX(IPRNT,NTESTL)
      NSTRIN = 0
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
      IORB2F = IORB1L + 1
      IORB2L = IORB2F+NORB2-1
      IORB3F = IORB2L + 1
      IORB3L = IORB3F+NORB3-1
* Loop over possible partitionings between RAS1,RAS2,RAS3
      DO 1001 IEL1 = NELMX1,NELMN1,-1
      DO 1003 IEL3 = NELMN3,NELMX3, 1
       IF(IEL1.GT. NORB1 ) GOTO 1001
       IF(IEL3.GT. NORB3 ) GOTO 1003
       IEL2 = NEL - IEL1-IEL3
       IF(IEL2 .LT. 0 .OR. IEL2 .GT. NORB2 ) GOTO 1003
       IFRST1 = 1
* Loop over RAS 1 occupancies
  901  CONTINUE
         IF( IEL1 .NE. 0 ) THEN
           IF(IFRST1.EQ.1) THEN
            CALL ISTVC2(IOC(1),0,1,IEL1)
            IFRST1 = 0
           ELSE
             CALL NXTORD_LUCI(IOC,IEL1,IORB1F,IORB1L,NONEW1)
             IF(NONEW1 .EQ. 1 ) GOTO 1003
           END IF
         END IF
         IF( NTEST .GE.500) THEN
           WRITE(6,*) ' RAS 1 string '
           CALL IWRTMA(IOC,1,IEL1,1,IEL1)
         END IF
         IFRST2 = 1
         IFRST3 = 1
* Loop over RAS 2 occupancies
  902    CONTINUE
           IF( IEL2 .NE. 0 ) THEN
             IF(IFRST2.EQ.1) THEN
              CALL ISTVC2(IOC(IEL1+1),IORB2F-1,1,IEL2)
              IFRST2 = 0
             ELSE
               CALL NXTORD_LUCI(IOC(IEL1+1),IEL2,IORB2F,IORB2L,NONEW2)
               IF(NONEW2 .EQ. 1 ) THEN
                 IF(IEL1 .NE. 0 ) GOTO 901
                 IF(IEL1 .EQ. 0 ) GOTO 1003
               END IF
             END IF
           END IF
           IF( NTEST .GE.500) THEN
             WRITE(6,*) ' RAS 1 2 string '
             CALL IWRTMA(IOC,1,IEL1+IEL2,1,IEL1+IEL2)
           END IF
           IFRST3 = 1
* Loop over RAS 3 occupancies
  903      CONTINUE
             IF( IEL3 .NE. 0 ) THEN
               IF(IFRST3.EQ.1) THEN
                CALL ISTVC2(IOC(IEL1+IEL2+1),IORB3F-1,1,IEL3)
                IFRST3 = 0
               ELSE
                 CALL NXTORD_LUCI(IOC(IEL1+IEL2+1),
     &           IEL3,IORB3F,IORB3L,NONEW3)
                 IF(NONEW3 .EQ. 1 ) THEN
                   IF(IEL2 .NE. 0 ) GOTO 902
                   IF(IEL1 .NE. 0 ) GOTO 901
                   GOTO 1003
                 END IF
               END IF
             END IF
             IF( NTEST .GE. 500) THEN
               WRITE(6,*) ' RAS 1 2 3 string '
               CALL IWRTMA(IOC,1,NEL,1,NEL)
             END IF
* Next string has been constructed , Enlist it !.
             NSTRIN = NSTRIN + 1
*. Symmetry of string
             ISYM = ISYMST(IOC,NEL)
C                   ISYMST(STRING,NEL)
*. occupation type of string
COLD         ITYP = IOCTP2(IOC,NEL,IOTYP)
C                   IOCTP2(STRING,NEL)
*
             NSTASO(ISYM,IOTYP) = NSTASO(ISYM,IOTYP)+ 1
*
           IF( IEL3 .NE. 0 ) GOTO 903
           IF( IEL3 .EQ. 0 .AND. IEL2 .NE. 0 ) GOTO 902
           IF( IEL3 .EQ. 0 .AND. IEL2 .EQ. 0 .AND. IEL1 .NE. 0)
     &     GOTO 901
 1003 CONTINUE
 1001 CONTINUE
*
*. The corresponding offset
*
      DO ISM = 1, NSMST
        IF(ISM.EQ.1) THEN
          ISTASO(ISM,IOTYP) = 1
        ELSE
          ISTASO(ISM,IOTYP) = ISTASO(ISM-1,IOTYP)+NSTASO(ISM-1,IOTYP)
        END IF
      END DO

      IF(NTEST.GE.5)
     &WRITE(6,*) ' Number of strings generated   ', NSTRIN
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Number of strings per sym for group = ', IOTYP
        WRITE(6,*) '================================================'
        CALL IWRTMA(NSTASO(1,IOTYP),1,NSMST,1,NSMST)
        WRITE(6,*) ' Offset for given symmetry for group = ', IOTYP
        WRITE(6,*) '================================================'
        CALL IWRTMA(ISTASO(1,IOTYP),1,NSMST,1,NSMST)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
* Testing new routines with active/passive divisions
*                          Jeppe, sept. 20 97
*
* Splitting the strings into active and passive parts allow us to
* write the a b loop  as
*
* Sigma(I_pa_a,I_pa_b,I_ac_a,I_ac_b) =
* sum(ijkl) <I_ac_a!Ea(ij)!J_ac_b><I_ac_b!Eb(kl)!I_ac_b>(ij!kl)
*           C(J_pa_a,J_pa_b,I_ac_a,I_pa_b)
*. One possibility is to use the above directly, and vectorize
*  over the passive strings.
*
* Another possibility is to use N-1 resolution in the alpha space
* as
*   1) Construct C(Ka_a,I_pa_a,I_pa_b,j,J_ac_b)
*   2) Obtain    S(Ka_a,I_pa_a,I_pa_b,I,I_ac_b) =
*      Sum(kl)   <I_ac_b!Eb(kl)!J_ac_b> Sum(ij)(ij!kl)
*      C(Ka_a,I_Pa_a,I_Pa_b,j,J_ac_b)
*
* The latter is the straightforward extension of LUCILLE's normal route.
*
* About reorganization of strings
*
* I will construct an array IREO(IA,IP)=I, that for given active
* and passive parts gives original number.
* Invoking symmetry :
* We are reordering strings of given supergroup and symmetry.
* The readressing goes thus as
* Loop over symmetry of active part => symmetry of passive part
*   Obtain offset for reordered strings
*          Number of active and passive strings of appropriate sym
*          Reorder array.
*
      SUBROUTINE REO_STR_SPGRP3(ISPGP,NIGRP,ISPGPSM,NSTR,NACTE,IACTE,
     &                          NACT,IACTGP,NAST_S,NPST_S,IBREO,IREO,
     &                          SIGNPA)
*
* Reorder strings of supergroup with given symmetry
* so passive groups are placed before active groups
*
* The division of into active and passive parts is based upon
* IACTE that give gasspaces of an operator.
*
* Output order
* Loop over Symmetry of active groups => Symmetry of passive groups
*   Loop over symmetry distributions of active groups
*   (generated by NEXTNUM2)
*     Loop over active strings of this symmetry distribution
*      Loop over symmetry distributions of passive groups
*      (generated by NEXTNUM2)
*        Loop over passive strings os this distribution
*        End of loop over passive strings
*      End of loop over symmetry distributions of passive groups
*     End of loop over active strings
*   End of loop over symmetrydistributions of active groups
* End of loop over symmetry of active group
*
* The reordered supergroup is thus organized as  a sequence of
* matrices C(p,a):
*
* Jeppe Olsen, September 1997
*
* ======
*. Input
* ======
*   ISPGP   : Groups defining supergroup
*   NISPGP  : Number of groups in supergroup
*   ISPGPSM : Symmetry of supergroup
*   NSTR    : Number of strings of this group
*   NACTE    : Number of active orbital spaces(needs not all be distinct)
*   IACTE    : The active gas orbital spaces (needs not all be distinct)
*
* ======
* Output
* ======
*
*    NAST_S : Number of active strings per symmetry
*    NPST_S : Number of passive strings per symmetry
*    IBREO  : Start of reorder loop for active strings of given sym
*    IREO   : The reorder array
*    NACT   : Number of active groups
*    IACT   : The active groups
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific Input
      DIMENSION IACTE(NACTE)
      DIMENSION ISPGP(NIGRP)
*. General Input
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
*. Output
*. Number of active strings per symmetry
      DIMENSION NAST_S(NSMST)
*. Number of passive strings per symmetry
      DIMENSION NPST_S(NSMST)
*. The active groups
      DIMENSION IACTGP(*)
*. Offsets for reordered strings with active strings of given sym
      DIMENSION IBREO(NSMST)
*. Reorder array : Note : New order to old order !
      DIMENSION IREO(NSTR)
*
*  ==============
*. Local scratch
*  ==============
*
      INTEGER ITPFGS(MXPNGAS)
*
      INTEGER IPAS(MXPNGAS)
      INTEGER IACT(MXPNGAS),IPASGP(MXPNGAS)

*
      INTEGER MXVL(MXPNGAS),MNVL(MXPNGAS)
*
      INTEGER NST_ACT(MXPNGAS),NST_PAS(MXPNGAS),NST_TOT(MXPNGAS)
*
      INTEGER N2STSMGP(MXPOBS,MXPNGAS)
C     INTEGER NSTGP(MXPNGAS),IREOGS(MXPNGAS)
      INTEGER IREOGS(MXPNGAS)
*
C     INTEGER ISMAC(MXPNGAS),ISMPA(MXPNGAS),ISMTO(MXPNGAS)
      INTEGER ISMTO(MXPNGAS)
      INTEGER ISTRAC(MXPNGAS),ISTPA(MXPNGAS),ISTTO(MXPNGAS)
*
*. Temporary solution ( for once )
      PARAMETER(LOFF=8*8*8*8)
      DIMENSION IOFF(LOFF)
      DIMENSION IDIST_ACT(MXPNGAS*LOFF), LDIST_ACT(LOFF)
      DIMENSION IDIST_PAS(MXPNGAS*LOFF), LDIST_PAS(LOFF)
*
CT    CALL QENTER('REO_PA')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================== '
        WRITE(6,*) ' REO_STR_SPGRP3 in service '
        WRITE(6,*) ' ========================== '
        WRITE(6,*)
        WRITE(6,*) '   Symmetry    Groups '
        WRITE(6,*) ' =============================================='
        WRITE(6,*)
C       WRITE(6,*)
        WRITE(6,'(A,3X,I4,3X,16(1X,I2))')
     & '     ',ISPGPSM,(ISPGP(II),II=1,NIGRP)
        WRITE(6,*) ' Number of operators ', NACTE
        WRITE(6,*) ' Gaspaces of operators ', (IACTE(I),I=1,NACTE)
      END IF
*. Info on Active operators : There can be several operators in
*  a active space. Obtain distinct active spaces
      IZERO = 0
      CALL ISETVC(IACTGP,IZERO,MXPNGAS)
      DO JACT = 1, NACTE
        IACTGP(IACTE(JACT)) = 1
      END DO
      NACT = 0
*. indeces of active groups
      DO JGRP = 1, NIGRP
        IF(IACTGP(IGSFGP(ISPGP(JGRP))).EQ.1) THEN
          NACT = NACT + 1
          IACT(NACT) = JGRP
        END IF
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of active groups ', NACT
        WRITE(6,*)
     &  ' Active gasspaces(index in IGRP) : ', (IACT(JACT),JACT=1,NACT)
      END IF
*
*. Info about supergroups
*
*. Number of strings per sym and type
      NGASL = 1
      DO IIGRP = 1, NIGRP
       IGAS = IGSFGP(ISPGP(IIGRP))
       IF( NELFGP(ISPGP(IIGRP)).GT.0) NGASL = IGAS
       CALL ICOPVE(NSTFSMGP(1,ISPGP(IIGRP)),N2STSMGP(1,IGAS),NSMST)
      END DO
*. For the moment
C     NGASL = NIGRP
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' NGASL',NGASL
       WRITE(6,*) ' N2STSMGP '
       CALL IWRTMA(N2STSMGP,NSMST,NGAS,MXPOBS,NGAS)
      END IF
*
* A supergroup contains a number of
* symmetry blocks, with each block defined by  given symmetries
* in the GAS SPACE. Obtain offsets to given symmetry combinations
*
      CALL TS_SYM_PNT2(ISPGP,NIGRP,MXVL,MNVL,ISPGPSM,IOFF,LOFF)
*. passive groups (we know these types)
      NPAS  = NIGRP - NACT
      JACT = 1
      JPAS = 0
      DO JGRP = 1, NIGRP
        IF(JGRP.EQ.IACT(JACT)) THEN
          IF (JACT.LT.NACT) JACT = JACT + 1
        ELSE
          JPAS = JPAS + 1
          IPAS(JPAS) = JGRP
        END IF
      END DO
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Number of passive types ',NPAS
        WRITE(6,*) ' Passive indeces in IGRP'
        CALL IWRTMA(IPAS,1,NPAS,1,NPAS)
      END IF
CCC*. Reorder array, old order of spaces to new order of
*. Reorder array, new order of indeces to old order
      JPAS = 1
      JACT = 1
      DO INDEX = 1, NIGRP
          IF(IACT(JACT).EQ.INDEX) THEN
            IREOGS(JACT+NPAS)  = INDEX
            IF(JACT.LT.NACT) JACT = JACT + 1
          ELSE IF(IPAS(JPAS).EQ.INDEX) THEN
            IREOGS(JPAS) = INDEX
            IF(JPAS.LT.NPAS) JPAS = JPAS + 1
          END IF
      END DO
*. Actual groups constituting active and passive parts
      JACT = 1
      JPAS = 1
      DO IIGRP = 1, NIGRP
        IF(IGSFGP(ISPGP(IIGRP)).EQ.IACT(JACT)) THEN
          IACTGP(JACT) = ISPGP(IIGRP)
          IF(JACT.LT.NACT) JACT = JACT + 1
        ELSE IF ( IGSFGP(ISPGP(IIGRP)).EQ.IPAS(JPAS) ) THEN
          IPASGP(JPAS) = ISPGP(IIGRP)
          IF(JPAS.LT.NPAS) JPAS = JPAS + 1
        END IF
      END DO
*. Sign for separating Active and passive groups, taken
*  as sign to bring passive strings in front.
      SIGNPA = 1.0D0
      DO JACTGP = 1, NACT
        JJACTGP = IACTGP(JACTGP)
        JJACTEL = NELFGP(JJACTGP)
        DO JPASGP = 1, NPAS
*
          JJPASGP = IPASGP(JPASGP)
          JJPASEL = NELFGP(JJPASGP)
          IF(JJPASGP.GT.JJACTGP) THEN
            SIGNPA = SIGNPA*(-1) **(JJPASEL*JJACTEL)
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Reordering of gasspaces, new  => old '
        CALL IWRTMA(IREOGS,1,NIGRP,1,NIGRP)
*
        WRITE(6,*) ' Groups constituting active parts '
        CALL IWRTMA(IACTGP,1,NACT,1,NACT)
*
        WRITE(6,*) ' Groups constituting passive parts '
        CALL IWRTMA(IPASGP,1,NPAS,1,NPAS)
      END IF
*. Last active and passive group with nonvanishing number of elecs
      NPASL = 1
      DO JPAS = 1, NPAS
        IF(NELFGP(IPASGP(JPAS)).NE.0) NPASL = JPAS
      END DO
*
*
      I_REO = 1
*. Loop over symmetries of active strings
      DO ISMAC_T = 1, NSMST
*. Symmetry of passive strings
        CALL SYMCOM(2,0,ISMAC_T,ISMPA_T,ISPGPSM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*)
     &    ' ISPGPSM ISMAC_T ISMPA_T',ISPGPSM,ISMAC_T,ISMPA_T
        END IF
        IBREO(ISMAC_T) = I_REO
*. Obtain symmetry distributions of active strings
        MXDIST = LOFF
        CALL SYM_DIST_FOR_SPGRP(IACTGP,NACT,ISMAC_T,
     &       NDIST_ACT,IDIST_ACT,LDIST_ACT,LEN_ACT,MXDIST)
        NAST_S(ISMAC_T) = LEN_ACT
*. And of passive strings
        IF( NAST_S(ISMAC_T) .GT. 0 ) THEN
          CALL SYM_DIST_FOR_SPGRP(IPASGP,NPASL,ISMPA_T,
     &         NDIST_PAS,IDIST_PAS,LDIST_PAS,LEN_PAS,MXDIST)
          NPST_S(ISMAC_T) = LEN_PAS
        ELSE
          NPST_S(ISMAC_T) = 0
        END IF
*. Note : NPST_S is ordered after symmetry of alpha strings
*. Loop over active strings : two loops 1)distributions 2)strings of this dist
        DO JDIST_ACT = 1, NDIST_ACT
C?        WRITE(6,*) ' Active distribution ', JDIST_ACT
*. Number of strings in each active GAS space
         DO JACT = 1, NACT
           NST_ACT(JACT) =
     &     N2STSMGP(IDIST_ACT(JACT+(JDIST_ACT-1)*NACT),IACT(JACT))
         END DO
*. Loop over active strings of given distribution
         DO ISTR_ACT = 1, LDIST_ACT(JDIST_ACT)
*. Obtain next active string
           IF(ISTR_ACT.EQ.1) THEN
            DO JACT = 1, NACT
             ISTRAC(JACT) = 1
            END DO
           ELSE
            CALL NXTNUM2_REV(ISTRAC,NACT,1,NST_ACT,NONEW_ACT)
            IF(NONEW_ACT.NE.0) THEN
              WRITE(6,*) ' Ran out of active strings !!'
            END IF
           END IF
*. Loop over distributions of passive strings
          DO JDIST_PAS = 1, NDIST_PAS
*. Number of strings in each passive GAS space
           DO JPAS = 1, NPASL
             NST_PAS(JPAS) =
     &       N2STSMGP(IDIST_PAS(JPAS+(JDIST_PAS-1)*NPASL),IPAS(JPAS))
           END DO
           DO JPAS = NPASL+1,NPAS
             NST_PAS(JPAS) = 1
           END DO
*. Offset to this distribution in original order
           DO JACT = 1, NACT
             ISMTO(IREOGS(NPAS+JACT))
     &     = IDIST_ACT(JACT+(JDIST_ACT-1)*NACT)
           END DO
           DO JPAS = 1, NPASL
             ISMTO(IREOGS(JPAS))
     &       = IDIST_PAS(JPAS+(JDIST_PAS-1)*NPASL)
           END DO
           DO JPAS = NPASL+1,NPAS
             ISMTO(IREOGS(JPAS)) = 1
           END DO
           I_ORIG_BASE = IOFF_SYM_DIST(ISMTO,NGASL,IOFF,MXVL,MNVL)
*. Number of strings, original order
           DO JGAS = 1, NGAS
            NST_TOT(JGAS) = N2STSMGP(ISMTO(JGAS),JGAS)
           END DO
C?         WRITE(6,*) ' NST_TOT  array '
C?         CALL IWRTMA(NST_TOT,1,NGAS,1,NGAS)
*. Loop over passive strings belonging to passive dist
           DO ISTR_PAS = 1, LDIST_PAS(JDIST_PAS)
C?         write(6,*) ' Passive string ', ISTR_PAS
            IF(ISTR_PAS.EQ.1) THEN
             DO JPAS = 1, NPASL
              ISTPA(JPAS) = 1
             END DO
            ELSE
             CALL NXTNUM2_REV(ISTPA,NPASL,1,NST_PAS,NONEW_PAS)
            END IF
*. Obtain original address of this string
*.. a : Original order of indeces
            DO JACT = 1, NACT
             ISTTO(IREOGS(JACT+NPAS)) = ISTRAC(JACT)
            END DO
            DO JPAS = 1, NPASL
             ISTTO(IREOGS(JPAS)) = ISTPA(JPAS)
            END DO
            DO JPAS = NPASL+1,NPAS
             ISTTO(IREOGS(JPAS)) = 1
            END DO
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Active part of string '
              CALL IWRTMA(ISTRAC,1,NACT,1,NACT)
              WRITE(6,*) ' passive part of string '
              CALL IWRTMA(ISTPA,1,NPASL,1,NPASL)
              WRITE(6,*) ' Complete string original order '
              CALL IWRTMA(ISTTO,1,NIGRP,1,NIGRP)
            END IF
*. and offset in original order
* The strings are generated as
* loop over GAS1
*   Loop over GAS2
*     ...
*       Loop over GASN
*       End of Loop over GASN
*     ...
*   End of Loop over GAS2
* End of loop over GAS1
* So address of string I1,I2, ..., In is
* (I1-1)*NI2*NI3*...NIN +
* (I2-1)*    NI3*...NIN + ...
*  NIN
            I_ORIG_REL = 0
            DO IGAS = 1, NIGRP-1
              I_ORIG_REL = (I_ORIG_REL + ISTTO(IGAS)-1)*NST_TOT(IGAS+1)
            END DO
            I_ORIG_REL = I_ORIG_REL + ISTTO(NIGRP)
            I_ORIG = I_ORIG_BASE+I_ORIG_REL - 1
*. And mapping
            IREO(I_REO) = I_ORIG
C?          WRITE(6,*) ' I_REO and I_ORIG' ,I_REO,I_ORIG
*. And update pointer, new order
            I_REO = I_REO + 1
           END DO
*          ^ End of loop over passive strings of given dist
          END DO
*         ^ End of loop over distributions of passive strings
        END DO
*       ^ End of loop over active strings of given dist
       END DO
*      ^ ENd of loop over distributions of active strings
      END DO
*     ^ End of loop over symmetries of active strings
      N_REO = I_REO-1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of elements in reordered array', N_REO
        WRITE(6,*) ' Reordering array : '
        CALL IWRTMA(IREO,1,N_REO,1,N_REO)
*
        WRITE(6,*) ' Offsets to blocks with given sym of active '
        CALL IWRTMA(IBREO,1,NSMST,1,NSMST)
        WRITE(6,*) ' Number of active per symmetry '
        CALL IWRTMA(NAST_S,1,NSMST,1,NSMST)
        WRITE(6,*) ' Number of passive per symmetry '
        CALL IWRTMA(NPST_S,1,NSMST,1,NSMST)
        WRITE(6,*) ' Sign change ', SIGNPA
*
      END IF
*
CT    CALL QEXIT('REO_PA')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SPGSPG_TO_CLASS(JOCTPA,JOCTPB,JOCCLS,NOCCLS,IOCCLS)
*
* Two supergroups of strings are given.
* Find the corresponding occupation class
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Absolute numbers
      IATP = 1
      IBTP = 2
      JJOCTPA = IBSPGPFTP(IATP) - 1 + JOCTPA
      JJOCTPB = IBSPGPFTP(IBTP) - 1 + JOCTPB
*
      JOCCLS = 0
      DO KOCCLS= 1, NOCCLS
        I_AM_OKAY = 1
        DO IGAS = 1, NGAS
          IF(NELFSPGP(IGAS,JJOCTPA)+NELFSPGP(IGAS,JJOCTPB)
     &       .NE. IOCCLS(IGAS,KOCCLS)                     )I_AM_OKAY=0
        END DO
        IF(I_AM_OKAY.EQ.1) JOCCLS = KOCCLS
      END DO
*
      IF(JOCCLS.EQ.0) THEN
        WRITE(6,*) ' Problem in SPGSPG_TO_CLASS : No match found ! '
      END IF
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info form SPGSPG_TO_CLASS'
        WRITE(6,*) ' input supergroups ',JJOCTPA,JJOCTPB
        WRITE(6,*) ' Output occupation class ', JOCCLS
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SPSPCLS(ISPSPCLS,ICLS,NCLS)
*
* Obtain mapping a-supergroup X b-supergroup => class
*
* Classes are specified by ICLS
*
* Jeppe Olsen, Jan 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "stinf.inc"
#include "strinp.inc"
*. Specific input
      INTEGER ICLS(*)
*. OUtput
      INTEGER ISPSPCLS(*)
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)

      CALL SPSPCLS_GAS(NOCTPA,NOCTPB,
     &            ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,ISPSPCLS,ICLS,NCLS,IPRDIA)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE SPSPCLS_GAS(NOCTPA,NOCTPB,IOCA,IOCB,NELFGP,
     &           MXPNGAS,NGAS,ISPSPCLS,ICLS,NCLS,IPRNT)
*
* Obtain mapping a-supergroup X b-supergroup => class
*
*
* =====
*.Input
* =====
*
* NOCTPA : Number of alpha types
* NOCTPB : Number of beta types
*
* IOCA(IGAS,ISTR) occupation of AS IGAS for alpha string type ISTR
* IOCB(IGAS,ISTR) occupation of AS IGAS for beta  string type ISTR
*
* MXPNGAS : Largest allowed number of gas spaces
* NGAS    : Actual number of gas spaces
*
*
* ======
*.Output
* ======
*
* ISPSPCLS(IATP,IBTP) => Class of this block of determinants
*                        =0 indicates unallowed(class less) combination
*
*.Input
      INTEGER IOCA(MXPNGAS,NOCTPA),IOCB(MXPNGAS,NOCTPB)
      INTEGER NELFGP(*)
      INTEGER ICLS(NGAS,NCLS)
*.Output
      INTEGER ISPSPCLS(NOCTPA,NOCTPB)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ISPSPCLS_GAS entered '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' IOCA and IOCB '
        CALL IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NGAS)
        CALL IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' ICLS '
        CALL IWRTMA(ICLS,NGAS,NCLS,NGAS,NCLS)
      END IF
*
      DO 100 IATP = 1, NOCTPA
        DO 90 IBTP = 1, NOCTPB
          IICLS = 0
          DO KCLS = 1, NCLS
            IAMOKAY = 1
            DO IGAS = 1, NGAS
              IEL = NELFGP(IOCA(IGAS,IATP))+NELFGP(IOCB(IGAS,IBTP))
              IF(IEL.NE.ICLS(IGAS,KCLS)) IAMOKAY = 0
            END DO
            IF(IAMOKAY.EQ.1) IICLS=KCLS
          END DO
          ISPSPCLS(IATP,IBTP) = IICLS
   90   CONTINUE
  100 CONTINUE
*
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Matrix giving classes for alpha-beta supergroups'
        WRITE(6,*)
        CALL IWRTMA(ISPSPCLS,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE STFDT3(NDET,IDET,IASTR2,IBSTR2,NSMST,NOCTPA,NOCTPB,
     &                  NSSOA,NSSOB,ISSOA,ISSOB,IOCOC,
     &                  ISMOST,ICOMBI)
*
* Obtain actual numbers for alpha- and beta strings corresponding to
* Given determinant numbers
*
* Jeppe Olsen, January 1989
*
*. Modified March 1993 so only string numbers are returned
*
* Changed into LUCIA form , June 1993
* Combinations enabled, september 1993
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Specific Input
      DIMENSION IDET(NDET)
*. General input
      DIMENSION NSSOA(NOCTPA,NSMST),NSSOB(NOCTPB,NSMST)
      DIMENSION ISSOA(NOCTPA,NSMST),ISSOB(NOCTPB,NSMST)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*)
*. Output
      DIMENSION IASTR2(NDET),IBSTR2(NDET)
*
COLD  IF( ICOMBI.NE.0) THEN
COLD    WRITE(6,*) ' STFDT3 does not work for COMBINATIONS'
COLD    WRITE(6,*) ' Enforced STOP'
COLD    Call Abend2('STFDT3')
COLD  END IF
*
      IDETF = 1
      DO 1000 IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBSM.LE.0) GOTO 1000
        IF(ICOMBI.NE.0. AND. IASM.LT.IBSM) GOTO 1000
        DO 999 IATP = 1,NOCTPA
          NIA = NSSOA(IATP,IASM)
          IAOFF = ISSOA(IATP,IASM)
          IF(ICOMBI.NE.0.AND.IASM.EQ.IBSM) THEN
            ISYM1 = 1
            MXBTP = IATP
          ELSE
            ISYM1 = 0
            MXBTP = NOCTPB
          END IF
          DO 901 IBTP = 1, MXBTP
            IF(IATP.EQ.IBTP.AND.ISYM1.EQ.1) THEN
              ISYM = 1
            ELSE
              ISYM = 0
            END IF
*
            IF(IOCOC(IATP,IBTP).LE.0) GOTO 901
            NIB = NSSOB(IBTP,IBSM)
            IF(NIA*NIB.EQ.0) GOTO 901
            IBOFF = ISSOB(IBTP,IBSM)
            IF(ISYM.EQ.0) THEN
              IDETL = IDETF + NIA*NIB-1
            ELSE
              IDETL = IDETF + NIA*(NIA+1)/2 - 1
            END IF
            DO 500 I = 1, NDET
              IF(IDET(I).GE.IDETF.AND.IDET(I).LE.IDETL) THEN
                IREL = IDET(I)-IDETF+1
                IF(ISYM.EQ.0) THEN
*. Find IAREL, IBREL so (IBREL-1)*NIA + IAREL = IREL
                  IBREL = (IREL-1)/NIA+1
                  IAREL = IREL - NIA*(IBREL-1)
                ELSE IF (ISYM.EQ. 1 ) THEN
*. Find IAREL, IBREL so IREL =
*                   (IBREL-1)*NIA + IAREL  - IBREL*(IBREL-1)/2
                   CALL UNPCPC(IREL,IAREL,NIA,IBREL)
                END IF
                IASTR2(I) = IAREL + IAOFF -1
                IBSTR2(I) = IBREL + IBOFF -1
              END IF
  500       CONTINUE
            IDETF = IDETL+1
  901     CONTINUE
  999   CONTINUE
 1000 CONTINUE
 1001 CONTINUE
C
      NTEST = 00
      IF( NTEST .GE. 2 ) THEN
         IF( ICOMBI .NE. 0 ) THEN
           WRITE(6,*) ' Combinations selected '
         ELSE
           WRITE(6,*) ' Determinants selected '
         END IF
C
         CALL IWRTMA(IDET,1,NDET,1,NDET)
         WRITE(6,*) ' Selected alpha and beta strings from STFDT3 '
         CALL IWRTMA(IASTR2,1,NDET,1,NDET)
         WRITE(6,*)
         CALL IWRTMA(IBSTR2,1,NDET,1,NDET)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE STRFDT_LUCI(INTSPC,JSMOST,IOCOC,NSBDET,IDET,IASTR2,
     &                       IBSTR2,ICOMBI)
*
* A CI expansion from internal subspace INTSPC and symmetry defined
* by JSMOST is considered. In this CI expansion, a set of NSBDET
*  combinations
* are given in terms of combination numbers IDET. Obtain the number of
* the corresponding alpha and betastrings
*
*
*. Revised March 22 1993 so only stringnumbers are stored
*  September 1993 : combinations enabled
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Specific input
      DIMENSION IDET(*)
*. General Input
      DIMENSION JSMOST(*),IOCOC(*)
*. Common blocks
#include "mxpdim.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "csm.inc"
*.
#include "wrkspc.inc"
*. Output
      DIMENSION IASTR2(*),IBSTR2(*)
*
      IATP = IASTFI(INTSPC)
      IBTP = IBSTFI(INTSPC)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
       CALL STFDT3(NSBDET,IDET,IASTR2,IBSTR2,NSMST,NOCTPA,NOCTPB,
     &             WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &             WORK(KISTSO(IATP)),WORK(KISTSO(IBTP)),
     &             IOCOC,JSMOST,ICOMBI)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE STRINF_GAS(IPRNT)
*
* Obtain string information for GAS expansion
*
* =====
*.Input
* =====
*
* /LUCINP/,/ORBINP/,/CSM/, /CGAS/, /GASSTR/
*
* =====
*.Output
* =====
*
* /STRINP/,/STINF/,/STRBAS/ and string information in STIN
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KFREEL
!               for addressing of WORK
#include "priunit.h"
#include "parluci.h"
*. Input
*     (and /LUCINP/ not occuring here )
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
#include "csm.inc"
#include "cstate.inc"
#include "lucinp.inc"
#include "stinf.inc"
#include "strinp.inc"
*
#include "wrkspc.inc"
*
*. A bit of scratch
C     DIMENSION IOCTYP(MXPNGAS)
*
      CALL QENTER('STRIN')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      CALL ZSTINF_GAS(NTEST)
*
**.3 : Static memory for string information
*
      CALL MEMSTR_GAS(NTEST)
*
** 4 : Info about group of strings
*
*.First free address
      CALL MEMMAN(KFREEL,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')

      DO IGRP = 1, NGRP
*. A gas group can be considered as a RAS group with 0 electrons in
*  RAS1, RAS3 !
        IGAS = IGSFGP(IGRP)
        IF(IGAS.EQ.1) THEN
          NORB1 = 0
        ELSE
          NORB1 = IELSUM(NOBPT,IGAS-1)
        END IF
        NORB2 = NOBPT(IGAS)
        NORB3 = NACOB-NORB1-NORB2
        MNRS1X = 0
        MXRS1X = 0
        MNRS3X = 0
        MXRS3X = 0
        IEL = NELFGP(IGRP)
        IOCTYPX = 1
*. Reverse lexical adresing schemes for each group of string
        CALL WEIGHT_LUCI(WORK(KZ(IGRP)),IEL,NORB1,NORB2,NORB3,
     &                   MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &                   WORK(KFREEL),NTEST)
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS(IEL,NORB1,NORB2,NORB3,
     &              MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              WORK(KFREEL),NACOB,
     &              WORK(KNSTSGP(1)),
     &              WORK(KISTSGP(1)),
     &              IOCTYPX,NSMST,IGRP,NTEST)
*. Construct the strings ordered by symmetry
        CALL GENSTR_GAS(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,
     &              WORK(KISTSGP(1)),IGRP,
     &              IOCTYPX,NSMST,WORK(KZ(IGRP)),WORK(KFREEL),
     &              WORK(KSTREO(IGRP)),WORK(KOCSTR(IGRP)),
     &              WORK(KFREEL+IOCTYPX*NSMST),IGRP,NTEST)
*
       CALL ICOPVE2(WORK(KNSTSGP(1)),1+(IGRP-1)*NSMST,NSMST,
     &              NSTFSMGP(1,IGRP))
      END DO
*
      IF(NTEST.GE.10) THEN
      WRITE(lupri,*) ' Number of strings per group and symmetry '
      CALL IWRTMA10(WORK(KNSTSGP(1)),NSMST,NGRP,NSMST,NGRP)
      WRITE(lupri,*) ' Number of strings per group and symmetry(2) '
      CALL IWRTMA10(NSTFSMGP,NSMST,NGRP,MXPNSMST,NGRP)
      END IF
*
*. Min and max of sym for each group
*
      DO IGP = 1, NGRP
       MX = 1
       DO ISM = 1, NSMST
         IF(NSTFSMGP(ISM,IGP).GT.0) MX = ISM
       END DO
*
       MN = NSMST
       DO ISM = NSMST,1,-1
         IF(NSTFSMGP(ISM,IGP).GT.0) MN = ISM
       END DO
*
       MINMAX_SM_GP(1,IGP) = MN
       MINMAX_SM_GP(2,IGP) = MX
*
      END DO
      IF(NTEST.GT.5) THEN
        WRITE(lupri,*) ' MINMAX array for sym of groups '
        WRITE(lupri,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP,1,NSMST,1,NSMST)
      END IF
*
*
* 4.5 : Creation/Annihilation mappings between different
*       types of strings
*
      DO IGRP = 1, NGRP
*
        IGAS = IGSFGP(IGRP)
        NGSOBP = NOBPT(IGAS)
*. First orbital in GAS spacce
        IGSOB = IELSUM(NOBPT,IGAS-1)+1
        IEL = NELFGP(IGRP)
        NSTINI = NSTFGP(IGRP)
*
*. Type of mapping : Only creation                  (LAC = 1)
*                    Only annihilation              (LAC = 2)
*                    Both annihilation and creation (LAC = 3)
* If only annihilation is present the string mapping arrays
* will only be over electronns
        IF(     ISTAC(IGRP,1).NE.0.AND.ISTAC(IGRP,2).NE.0) THEN
          LAC = 3
          IEC = 1
          LROW = NGSOBP
        ELSE IF(ISTAC(IGRP,1).NE.0.AND.ISTAC(IGRP,2).EQ.0) THEN
          LAC = 1
          IEC = 2
          LROW = IEL
        ELSE IF(ISTAC(IGRP,1).EQ.0.AND.ISTAC(IGRP,2).NE.0) THEN
          LAC = 2
          IEC = 0
          LROW = NGSOBP
        ELSE IF(ISTAC(IGRP,1).EQ.0.AND.ISTAC(IGRP,2).EQ.0) THEN
          LAC = 0
          IEC = 0
          LROW = 0
        END IF
*. Zero
        IF(LAC.NE.0) THEN
          CALL IZERO(WORK(KSTSTM(IGRP,1)),LROW*NSTINI)
          CALL IZERO(WORK(KSTSTM(IGRP,2)),LROW*NSTINI)
        END IF
*
        IF(ISTAC(IGRP,2).NE.0) THEN
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS(WORK(KOCSTR(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,WORK(KZ(JGRP)),
     &         WORK(KSTREO(JGRP)),0,IDUM,IDUM,
     &         WORK(KSTSTM(IGRP,1)),WORK(KSTSTM(IGRP,2)),NACOB,NTEST)
        END IF
        IF(ISTAC(IGRP,1).NE.0) THEN
          JGRP = ISTAC(IGRP,1)
          CALL ANNSTR_GAS(WORK(KOCSTR(IGRP)),NSTFGP(IGRP),
     &         NSTFGP(JGRP),IEL,NGSOBP,IGSOB,WORK(KZ(JGRP)),
     &         WORK(KSTREO(JGRP)),0,IDUM,IDUM,
     &         WORK(KSTSTM(IGRP,1)),WORK(KSTSTM(IGRP,2)),NACOB,IEC,
     /         LROW,NTEST)
        END IF
      END DO
*
*
*. Now to supergroups , i.e. strings of with given number of elecs in
*  each GAspace
*
      CALL ISETVC(NSTFSMSPGP,0,MXPNSMST*NTSPGP)
      MXNSTR = -1
      DO ITP = 1, NSTTYP
*. Loop over supergroups of given type . i.e. strings
*  with given occupation in each GAS space
        DO IGRP = 1, NSPGPFTP(ITP)
          IGRPABS = IGRP-1 + IBSPGPFTP(ITP)
          CALL NSTPTP_GAS(NGAS,ISPGPFTP(1,IGRPABS),
     &                    WORK(KNSTSGP(1)),NSMST,
     &                    WORK(KNSTSO(ITP)),IGRP,MXNSTRFSG,
     &                    NSMCLS,NSMCLSE,NSMCLSE1)
*
          MXSMCLS   = MAX(MXSMCLS,NSMCLS)
          MXSMCLSE  = MAX(MXSMCLSE,NSMCLSE)
          MXSMCLSE1 = MAX(MXSMCLSE1,NSMCLSE1)
*
          MXNSTR = MAX(MXNSTR,MXNSTRFSG)
        END DO
*
        CALL ICOPMT(WORK(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &              NSTFSMSPGP(1,IBSPGPFTP(ITP)),MXPNSMST,NSPGPFTP(ITP))
*. Corresponding offset array : Each supergroup is generated individually
*. so each supergroup starts with offset 1 !
        CALL ZSPGPIB(WORK(KNSTSO(ITP)),WORK(KISTSO(ITP)),
     &                NSPGPFTP(ITP),NSMST)
*
       IF(NTEST.GE.5) THEN
        IF(LUCI_MYPROC.EQ.LUCI_MASTER)THEN
          WRITE(lupri,*)
     &    ' Number of strings per sym (row) and supergroup(column)',
     &    ' for type = ', ITP
          CALL IWRTMA(WORK(KNSTSO(ITP)),NSMST,NSPGPFTP(ITP),
     &                NSMST,NSPGPFTP(ITP))
          WRITE(lupri,'(A,3I6)') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                         NSMCLS,NSMCLSE,NSMCLSE1
        ENDIF
       END IF
*
      END DO
*. Number of electron in each AS for each supergroup
      CALL ZNELFSPGP(NTEST)
*
* Number of holes per supergroup
      DO IISPGP = 1, NTSPGP
        NHOLE = 0
        DO IGAS = 1, NGAS
          IF(IPHGAS(IGAS).EQ.2) NHOLE = NHOLE + NELFSPGP(IGAS,IISPGP)
        END DO
        NHLFSPGP(IISPGP) = NHOLE
      END DO
      IF(NTEST.GE.10) THEN
        WRITE(lupri,*) 
     &  ' Number of electrons in hole spaces per supergroup '
        CALL IWRTMA(NHLFSPGP,1,NTSPGP,1,NTSPGP)
      END IF
*. Largest number of strings belonging to given supergroup
*. Largest Occupation block for given supergroup and sym
      MAX_STR_OC_BLK = -1
      MAX_STR_SPGP = 0
      DO ISPGP = 1, NTSPGP
        NSTR         = IELSUM(NSTFSMSPGP(1,ISPGP),NSMST)
        NEL          = IELSUM(NELFSPGP(1,ISPGP),NGAS)
        MAX_STR_SPGP = MAX(MAX_STR_SPGP,NSTR)
        DO ISTSM = 1, NSMST
          MAX_STR_OC_BLK
     &  = MAX(MAX_STR_OC_BLK,(NEL+4)*NSTFSMSPGP(ISTSM,ISPGP))
CMOD &  = MAX(MAX_STR_OC_BLK,NEL*NSTFSMSPGP(ISTSM,ISPGP))
        END DO
      END DO
*
      IF(NTEST.GE.2) THEN
      WRITE(lupri,*)
     & ' Largest number of strings of given supergroup        ',
     & MAX_STR_SPGP
      WRITE(lupri,*) ' Largest block of string occupations ',
     &              MAX_STR_OC_BLK
*
      WRITE(lupri,*)
     & ' Largest number of strings of given supergroup and sym', MXNSTR
      END IF
*
*
* Possible occupation classes
*
      CALL OCCLS(2,NMXOCCLS,WORK(KIOCLS),NACTEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),0,0)
*
* Maps creation/annihilation of given gas orb from given supergroup
* gives new supergroup.
*
      CALL IZERO(WORK(KSPGPCR),NGAS*NTSPGP)
      CALL IZERO(WORK(KSPGPAN),NGAS*NTSPGP)
*
      DO ISTTYP = 1,NSTTYP
*. Creation map from this type
        IIEL = NELFTP(ISTTYP)
*. Type of string with one elec more
        ISTTYPC = 0
        DO JSTTYP = 1, NSTTYP
          IF(MOD(ISTTYP,2).EQ.MOD(JSTTYP,2).AND.
     &       NELFTP(JSTTYP).EQ.IIEL+1           ) ISTTYPC = JSTTYP
        END DO
C?      WRITE(6,*) ' ISTTYP and ISTTYPC ',ISTTYP,ISTTYPC
        IF(NSPGPFTP(ISTTYP).GT.0) THEN
*
        IF(ISTTYPC.GE.1.AND.NSPGPFTP(ISTTYPC).GT.0) THEN
           CALL SPGP_AC(NELFSPGP(1,1),NSPGPFTP(ISTTYP),
     &                  NELFSPGP(1,1),NSPGPFTP(ISTTYPC),
     &                  NGAS,MXPNGAS,2,WORK(KSPGPCR),
     &                  IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPC))
        ELSE
        END IF
*. Annihilation maps
        ISTTYPA = 0
        DO JSTTYP = 1, NSTTYP
          IF(MOD(ISTTYP,2).EQ.MOD(JSTTYP,2).AND.
     &       NELFTP(JSTTYP).EQ.IIEL-1           ) ISTTYPA = JSTTYP
        END DO
C?      WRITE(6,*) 'ISTTYP, ISTTYPA', ISTTYP,ISTTYPA
        IF(ISTTYPA.GE.1 .AND.NSPGPFTP(ISTTYPA).GT.0) THEN
           CALL SPGP_AC(NELFSPGP(1,1),NSPGPFTP(ISTTYP),
     &                  NELFSPGP(1,1),NSPGPFTP(ISTTYPA),
     &                  NGAS,MXPNGAS,1,WORK(KSPGPAN),
     &                  IBSPGPFTP(ISTTYP),IBSPGPFTP(ISTTYPA))
        END IF
        END IF
      END DO
*
      CALL QEXIT('STRIN')
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE STRTYP_GAS(IPRNT)
*
* Find groups of strings in each GA space
*
* Output : /GASSTR/
*
* Jeppe Olsen, Oct 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "crun.inc"
*. Local scratch
      DIMENSION IOCTYP(MXPSTT),IREOSPGP(MXPSTT),ISCR(MXPSTT)
*
      CALL QENTER('STRTY')
      NTESTL = 00
      NTEST = MAX(IPRNT,NTESTL)
*. As input NCISPC GAS spaces IGSOCCX are given.
* Obtain space that cantains all these as special cases
*
C?    WRITE(6,*) ' NCISPC ', NCISPC
      DO IGAS = 1, NGAS
       MINI = IGSOCCX(IGAS,1,1)
       MAXI = IGSOCCX(IGAS,2,1)
C?     WRITE(6,*) ' MINI and MAXI for ISPC = 1 ',MINI,MAXI
       DO ICISPC = 2, NCISPC
        MINI = MIN(MINI,IGSOCCX(IGAS,1,ICISPC))
        MAXI = MAX(MAXI,IGSOCCX(IGAS,2,ICISPC))
C?     WRITE(6,*) ' MINI and MAXI for ISPC =  ',ICISPC,MINI,MAXI
       END DO
       IGSOCC(IGAS,1) = MINI
       IGSOCC(IGAS,2) = MAXI
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Compound GAS space : '
        WRITE(6,*) ' ====================='
        WRITE(6,'(A)')
        WRITE(6,'(A)') '         Min. occ    Max. occ '
        WRITE(6,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(6,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCC(IGAS,1),IGSOCC(IGAS,2)
        END DO
      END IF

*
*. Find min and max number of elecs in each subspace
*
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          MNGSOC(IGAS) = IGSOCC(IGAS,1)
          MXGSOC(IGAS) = IGSOCC(IGAS,2)
        ELSE
          MXGSOC(IGAS) = IGSOCC(IGAS,2)-IGSOCC(IGAS-1,1)
          MNGSOC(IGAS) = MAX(0,IGSOCC(IGAS,1)-IGSOCC(IGAS-1,2))
        END IF
      END DO
*
*. Particle and hole spaces  :
*  Hole spaces are always more than half occupied
*
      DO IGAS = 1, NGAS
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
C       write(6,*) ' Halffilled spaces  defined hole spaces'
        IF(IUSE_PH.EQ.1) THEN
C         IF(2*MNGSOC(IGAS).GE.NOBPT(IGAS)) THEN
          IF(2*MNGSOC(IGAS).GT.NOBPT(IGAS)) THEN
            IPHGAS(IGAS) = 2
          ELSE
             IPHGAS(IGAS) = 1
          END IF
        ELSE IF(IUSE_PH.EQ.0) THEN
          IPHGAS(IGAS) = 1
        END IF
*
      END DO
*. Min number of electrons in hole spaces
      MNHL  = 0
      DO IGAS = 1, NGAS
        IF(IPHGAS(IGAS).EQ.2) THEN
          MNHL = MNHL + MNGSOC(IGAS)
        END IF
      END DO
*
C     WRITE(6,*) ' MNHL' , MNHL
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Min and Max occupation in each GAS space: '
        WRITE(6,'(A)') ' ========================================= '
        WRITE(6,*)
        DO IGAS = 1,  NGAS
          WRITE(6,'(A,I2,4X,2I3)')
     &    '  GAS',IGAS,MNGSOC(IGAS),MXGSOC(IGAS)
        END DO
*
        WRITE(6,*) ' Particle (=1) or hole (=2) spaces '
        WRITE(6,*) ' =================================='
        CALL IWRTMA(IPHGAS,1,NGAS,1,NGAS)
       END IF
*.
*. Occupation classes corresponding to largest CI space
*
C     CALL QENTER('OCCLS')
      CALL OCCLS(1,NOCCLS,IOCCLS,NACTEL,NGAS,IGSOCC(1,1),IGSOCC(1,2),
     &           0,0)
      NMXOCCLS = NOCCLS
*
* Split into alpha and beta parts
*
*. Number of alpha. and beta electrons
*
      NAEL = (MS2 + NACTEL ) / 2
      NBEL = (NACTEL - MS2 ) / 2
*
      IF(NAEL + NBEL .NE. NACTEL ) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
        WRITE(6,*)
     &  ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
        Call
     &    Abend2(' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY ')
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I6)')   MS2,NACTEL,NAEL,NBEL
      END IF
*
      IF(NAEL + NBEL .NE. NACTEL ) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
        WRITE(6,*)
     &  ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
        Call
     &    Abend2(' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY ')
      END IF

*. Number of electrons to be subtracted or added
*
      MAXSUB = 2
      MAXADD = 2
*. electrons are only added for systems that atleast have halffilled
*. shells
      IGRP = 0
      MXAL = NAEL
      MNAL = NAEL
      MXBL = NBEL
      MNBL = NBEL
      NORBL = NTOOB
      DO IGAS = 1, NGAS
*. occupation constraints 1
       MXA1 = MIN(MXGSOC(IGAS),NOBPT(IGAS),MXAL)
       MXB1 = MIN(MXGSOC(IGAS),NOBPT(IGAS),MXBL)
       MNA1 = MAX(0,MNGSOC(IGAS)-MXA1)
       MNB1 = MAX(0,MNGSOC(IGAS)-MXB1)
*. Additional checks can be made here
       MXA = MXA1
       MXB = MXB1
       MNA = MNA1
       MNB = MNB1
*
       MXAL = MXAL - MNA
       MNAL = MAX(0,MNAL-MXA)
       MXBL = MXBL - MNB
       MNBL = MAX(0,MNBL-MXB)
*
       IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation numbers for IGAS = ', IGAS
        WRITE(6,*) ' MXAL MNAL MXBL MNBL ',MXAL,MNAL,MXBL,MNBL
        WRITE(6,*) ' MXA MNA MXB MNB ',MXA,MNA,MXB,MNB
       END IF
*
       MNAB = MIN(MNA,MNB)
       MXAB = MAX(MXA,MXB)
*. Additional holes only allowed in particle spaces
       IF(IPHGAS(IGAS).EQ.1) THEN
         MNAB = MAX(0,MNAB-MAXSUB)
       ELSE IF(IPHGAS(IGAS).EQ.2) THEN
         MNAB = MNAB
       END IF
*. Additional electrons allowed in hole spaces
       IF(IPHGAS(IGAS).EQ.2)  MXAB = MIN(MXAB + 2,NOBPT(IGAS))
*
       IF(NTEST.GE.100) WRITE(6,*) ' MNAB,MXAB',MNAB,MXAB
       NGPSTR(IGAS) = MXAB-MNAB+1
       IBGPSTR(IGAS) = IGRP + 1
       MNELFGP(IGAS) = MNAB
       MXELFGP(IGAS) = MXAB
*
       IADD = 0
       DO JGRP = IGRP+1,IGRP+NGPSTR(IGAS)
         IF(JGRP.GT.MXPSTT) THEN
           WRITE(6,*) ' Too many string groups '
           WRITE(6,*) ' Current limit ', MXPSTT
           WRITE(6,*) ' STOP : GASSTR, Too many string groups'
           Call Abend2(' GASSTR, Too many string groups')
          END IF
*
         IADD = IADD + 1
         IEL = MNAB-1+IADD
         NELFGP(JGRP) = IEL
         IGSFGP(JGRP) = IGAS
         NSTFGP(JGRP) = IBION(NOBPT(IGAS),IEL)
       END DO
       IGRP = IGRP + NGPSTR(IGAS)
      END DO
      NGRP = IGRP
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Information about Groups of strings '
        WRITE(6,'(A)') ' =================================== '
        WRITE(6,*)
        WRITE(6,*) '     GAS  MNEL  MXEL IBGRP  NGRP'
        WRITE(6,*) '    ============================'
        DO IGAS = 1, NGAS
          WRITE(6,'(5(2X,I4))') IGAS,MNELFGP(IGAS),
     &          MXELFGP(IGAS),IBGPSTR(IGAS),NGPSTR(IGAS)
        END DO
        WRITE(6,'(A,I3)')
     &  ' Total number of groups generated ', NGRP
*
        WRITE(6,'(A)') ' Information about each string group '
        WRITE(6,'(A)') ' ===================================='
        WRITE(6,*)
        IITYPE = 0
        WRITE(6,'(A)') ' GROUP  GAS   NEL      NSTR '
        WRITE(6,'(A)') ' ==========================='
        DO IGRP = 1, NGRP
          IITYPE = IITYPE + 1
          WRITE(6,'(3(2X,I4),2X,I8)')
     &    IITYPE,IGSFGP(IGRP),NELFGP(IGRP),NSTFGP(IGRP)
        END DO
      END IF
*
*. Creation-annihilation connections between groups
*
      DO IGRP = 1, NGRP
        ISTAC(IGRP,1) = 0
        ISTAC(IGRP,2) = 0
        DO JGRP = 1, NGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)-1) ISTAC(IGRP,2) = JGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)+1) ISTAC(IGRP,1) = JGRP
        END DO
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================================'
        WRITE(6,*) ' Annihilation / Creation connections'
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        CALL IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
      END IF
*.
*
* Construct number of type ( combinations of groups ) with nael and nbel strings
*
*
* Type 1 : NAEL electrons
*      2 : NBEL ELECTRONS
*      3 : NAEL -1 ELECTRONS
*      4 : NBEL -1 ELECTRONS
*      5 : NAEL -2 ELECTRONS
*      6 : NBEL -2 ELECTRONS
*
      NSTTYP = 6
      NSTTP = 6
*
      IF(IUSE_PH.EQ.1) THEN
*. allow N+1,N+2 resolution string
        NSTTYP = 10
        NSTTP  = 10
      END IF
*. alpha
      NELEC(1) = NAEL
      NELFTP(1) = NAEL
*. beta
      NELEC(2) = NBEL
      NELFTP(2) = NBEL
*. alpha -1
      NELEC(3) = NAEL-1
      NELFTP(3) = NAEL-1
*. beta  -1
      NELEC(4) = NBEL-1
      NELFTP(4) = NBEL-1
*. alpha -2
      NELEC(5) = NAEL-2
      NELFTP(5) = NAEL-2
*. beta  -2
      NELEC(6) = NBEL-2
      NELFTP(6) = NBEL-2
*
      IF(IUSE_PH.EQ.1) THEN
*. Alpha + 1
        NELEC(7) = NAEL+1
        NELFTP(7) = NAEL+1
*. beta  + 1
        NELEC(8) = NBEL+1
        NELFTP(8) = NBEL+1
*. Alpha + 2
        NELEC(9) = NAEL+2
        NELFTP(9) = NAEL+2
*. beta  + 2
        NELEC(10) = NBEL+2
        NELFTP(10) = NBEL+2
      END IF
*. Can easily be extended to relativistic case !!
      DO ITP = 1, NSTTYP
        NOCTYP(ITP) = 0
        NSPGPFTP(ITP) = 0
      END DO
*
* Loop over types, i.e.  given number of electrons
*
      IOFF = 1
      NABEL = NAEL + NBEL
      NSPGP_TOT = 0
      DO 2000 ITYP = 1, NSTTYP
*. Number of electrons in reference space ( alpha or beta )
        IF(MOD(ITYP,2).EQ.1) THEN
*. alpha type
          NELEC_REF = NELEC(1)
        ELSE
*. Beta type
          NELEC_REF = NELEC(2)
        END IF
*. If we are studying beta type, and number of alpha and beta
* electrons are identical, just refer to alpha
        IF(NAEL.EQ.NBEL.AND.MOD(ITYP,2).EQ.0) THEN
          IBSPGPFTP(ITYP) =  IBSPGPFTP(ITYP-1)
          NOCTYP(ITYP) =   NOCTYP(ITYP-1)
          NSPGPFTP(ITYP) =  NSPGPFTP(ITYP-1)
        ELSE
*
*. Number of electrons removed compared to reference
        IDEL = NELEC(ITYP) - NELEC_REF
C?      WRITE(6,*) '  GASSPC : ITYP IDEL ', ITYP,IDEL
*. Initial type of strings, relative to offset for given group
        DO IGAS = 1, NGAS
          IOCTYP(IGAS) = 1
        END DO
        NSPGP = 0
        IBSPGPFTP(ITYP) = IOFF
        IF(NELEC(ITYP).LT.0) THEN
          NOCTYP(ITYP) = 0
          NSPGPFTP(ITYP) =  0
          GOTO 2000
        END IF
*. Number of electrons in present type
*. Loop over  SUPER GROUPS with current nomenclature!
*. Temp max for loop
        MXLOOP = 10000
        IONE = 1
        NLOOP = 0
 1000   CONTINUE
*. Number of electrons in present supergroup
          NEL = 0
          DO IGAS = 1, NGAS
            NEL = NEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
          END DO
*
          IF(NEL.GT.NELEC(ITYP)) THEN
*. If the number of electrons is to large find next number that
* can be correct.
* The following uses that within a given GAS space
* the number of elecs increases as the type number increases
*
*. First integer  that can be reduced
            IRED = 0
            DO IGAS = 1, NGAS
              IF(IOCTYP(IGAS).NE.1) THEN
                IRED = IGAS
                GOTO 888
              END IF
            END DO
  888       CONTINUE
            IF(IRED.EQ.NGAS) THEN
              NONEW = 1
            ELSE IF(IRED.LT.NGAS) THEN
              IOCTYP(IRED) = 1
*. Increase remanining part
              CALL NXTNUM2(
     &        IOCTYP(IRED+1),NGAS-IRED,IONE,NGPSTR(IRED+1),NONEW)
            END IF
            GOTO 2803
          END IF

          IF(NEL.EQ.NELEC(ITYP)) THEN
*. test 1 has been passed, check additional occupation constraints
*
           I_AM_OKAY = 1
*. Number of extra holes in hole spaces
CE         IF(IUSE_PH.EQ.1) THEN
CE           IDELP = 0
CE           IDELM = 0
CE           DO IGAS = 1, NGAS
CE             IF(IPHGAS(IGAS).EQ.2) THEN
CE               NELH =  NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
CE               IF(NELH.LT.MNGSOC(IGAS))
CE    &          IDELM = IDELM +MNGSOC(IGAS)-NELH
CE               IF(NELH.GT.MXGSOC(IGAS))
CE    &          IDELP = IDELP + NELH-MXGSOC(IGAS)
CE             END IF
CE           END DO
CE           IF(IDELM.GT.0.OR.IDELP.GT.MAX(0,IDEL)) THEN
CE             I_AM_OKAY = 0
CE             WRITE(6,*) ' P/H rejected supergroup '
CE             CALL IWRTMA(IOCTYP,1,NGAS,1,NGAS)
CE             WRITE(6,*) ' IDELM, IDELP ', IDELM, IDELP
CE           END IF
CE         END IF
*
*. Check from above
*
           DO IGAS = NGAS, 1, -1
*. Number of electrons when all electrons of AS IGAS have been added
             IF(IGAS.EQ.NGAS ) THEN
               IEL = MAX(NABEL,NABEL+IDEL)
             ELSE
               IEL = IEL-NELFGP(IOCTYP(IGAS+1)+IBGPSTR(IGAS+1)-1)
               IF(IEL.LT.MAX(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &         I_AM_OKAY = 0
             END IF
           END DO
*
* Check from below
*
           IEL = 0
           IOELMX = 0
           DO IGAS = 1, NGAS
             IEL = IEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
             IOELMX = IOELMX+NOBPT(IGAS)
             IF(IEL+IOELMX.LT.MIN(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &       I_AM_OKAY = 0
           END DO
*
           IF(I_AM_OKAY.EQ.1) THEN
*. passed !!!
             NSPGP = NSPGP + 1
*. Copy supergroup to ISPGPFTP with absolute group numbers
             DO IGAS = 1, NGAS
               ISPGPFTP(IGAS,IOFF-1+NSPGP)
     &       = IOCTYP(IGAS)+IBGPSTR(IGAS)-1
             END DO
           END IF
*
          END IF
*. Next type of strings
          IONE = 1
          CALL NXTNUM2(IOCTYP,NGAS,IONE,NGPSTR,NONEW)
 2803   CONTINUE
        IF(NONEW.EQ.0) GOTO 1000
*. End of loop over possible supergroups, save information about current type
        IOFF = IOFF + NSPGP
        NOCTYP(ITYP) = NSPGP
        NSPGPFTP(ITYP) =  NSPGP
        NSPGP_TOT =  NSPGP_TOT +  NSPGP
      END IF
 2000 CONTINUE
      NTSPGP = NSPGP_TOT
*
      IF(NSPGP_TOT .GT. MXPSTT ) THEN
        WRITE(6,*) ' Too many super groups = ', NSPGP_TOT
        WRITE(6,*) ' Increase MXPSTT to this value'
        WRITE(6,*) ' See you later '
        WRITE(6,*)
        WRITE(6,*) ' STOP Increase MXPSTT '
        Call Abend2(' Increase MXPSTT')
      END IF
*
*. Reorder supergroups according to dimension
*
      DO ITYP= 1, NSTTP
        IBTYP = IBSPGPFTP(ITYP)
        NSPGP = NSPGPFTP(ITYP)
*.Dimension of supergroups
        DO ISPGP = 1, NSPGP
          IDIM = 1
          DO JGAS = 1, NGAS
            IDIM = IDIM * NSTFGP(ISPGPFTP(JGAS,ISPGP+IBTYP-1))
          END DO
          IOCTYP(ISPGP) = IDIM
        END DO
*. Reorder
C            ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
        CALL ORDINT(IOCTYP,ISCR,NSPGP,IREOSPGP,NTEST)
C?      WRITE(6,*) ' IREO array '
C?      CALL IWRTMA(IREOSPGP,1,NSPGP,1,NSPGP)
*. And reorder the definition of supergroups
        DO ISPGP = 1, NSPGP
          CALL ICOPVE(ISPGPFTP(1,ISPGP+IBTYP-1),NELFSPGP(1,ISPGP),NGAS)
        END DO
        DO ISPGP_N = 1, NSPGP
          ISPGP_O = IREOSPGP(ISPGP_N)
          CALL ICOPVE(NELFSPGP(1,ISPGP_O),ISPGPFTP(1,ISPGP_N+IBTYP-1),
     &                NGAS)
        END DO
*
      END DO
*     ^ End of loop over types
*
      IF(NTEST.GE.2) THEN
       WRITE(6,*) ' Total number of super groups ', NTSPGP
       WRITE(6,*) ' Number of alpha supergroups  ', NSPGPFTP(1)
       WRITE(6,*) ' Number of beta  supergroups  ', NSPGPFTP(2)
      END IF

       WRITE(6,*)
*
      IF (NTEST.GE.5) THEN
        WRITE(6,*) ' Information about types of strings'
        WRITE(6,*) ' =================================='
        WRITE(6,*)
        DO ITYP = 1, NSTTYP
          WRITE(6,*)
          WRITE(6,*) '      Type : ', ITYP
          WRITE(6,*) '      ==============='
          WRITE(6,*) '      Number of electrons  ',NELFTP(ITYP)
          WRITE(6,*) '      Number of super groups ', NSPGPFTP(ITYP)
          WRITE(6,*) '      Supergroups '
          DO ISPGP = 1, NSPGPFTP(ITYP)
            IOFF = IBSPGPFTP(ITYP)
            CALL IWRTMA(ISPGPFTP(1,IOFF-1+ISPGP),1,NGAS,1,NGAS)
          END DO
        END DO
      END IF
*
      CALL QEXIT('STRTY')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TS_SYM_PNT2(IGRP,NIGRP,MAXVAL,MINVAL,ISYM,
     &           IPNT,LPNT)
*
* Construct pointers to start of symmetry distributions
* for supergroup of strings with given symmetry
*
* The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
* is given as
*     1
*     +  (ISM1-MINVAL(1))
*     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
* Where L is the last group of strings with nonvanishing occupation
*
* Jeppe Olsen, September 1997
*
* Version 2 : Uses IGRP and NIGRP to define supergroup
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
#include "orbinp.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "csm.inc"
*. Specific Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
*. Output
      INTEGER MINVAL(*),MAXVAL(*),IPNT(*)
*
      NTEST = 00
*. Info on groups of strings in supergroup
      NGASL = 1
      DO IGAS = 1, NIGRP
       ITPFGS(IGAS) = IGRP(IGAS)
       IF(NELFGP(IGRP(IGAS)).GT.0) NGASL = IGAS
*. Number of strings per symmetry in each gasspace
C       CALL ICOPVE2(WORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
C    &               NNSTSGP(1,IGAS))
        CALL ICOPVE(NSTFSMGP(1,IGRP(IGAS)),NNSTSGP(1,IGAS),NSMST)
      END DO
*
C     NGASL = NIGRP
*
C     DO IGAS = 1, NIGRP
C       DO ISMST = 1, NSMST
C         IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
C       END DO
C       DO ISMST = NSMST,1,-1
C         IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
C       END DO
C     END DO
      DO IGAS = 1, NIGRP
        MINVAL(IGAS) = MINMAX_SM_GP(1,IGRP(IGAS))
        MAXVAL(IGAS) = MINMAX_SM_GP(2,IGRP(IGAS))
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' MINVAL and MAXVAL '
        CALL IWRTMA(MINVAL,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MAXVAL,1,NIGRP,1,NIGRP)
        WRITE(6,*) ' NIGRP = ', NIGRP
      END IF

*. Total number of symmetry blocks that will be generated
      NBLKS = 1
      DO IGAS = 1, NGASL-1
       NBLKS = NBLKS*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
      END DO
      IF(NBLKS.GT.LPNT) THEN
        WRITE(6,*) ' Problem in TS_SYM_PNT'
        WRITE(6,*) ' Dimension of IPNT too small'
        WRITE(6,*) ' Actual and required length',NBLKS,LPNT
        WRITE(6,*)
        WRITE(6,*) ' I will Stop and wait for instructions'
        Call Abend2( ' TS_SYM_PNT too small' )
      END IF
*. Loop over symmetry blocks in standard order
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
C       ISTSMM1 = 1
C       DO IGAS = 1, NGASL -1
C         CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
C         ISTSMM1 = JSTSMM1
C       END DO
        ISTSMM1 = ISYMSTR(ISMFGS,NGASL-1)
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISYM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL-1
          IOFF = IOFF + (ISMFGS(IGAS)-MINVAL(IGAS))*MULT
          MULT = MULT * (MAXVAL(IGAS)-MINVAL(IGAS)+1)
        END DO
*
        IPNT(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOFF, IPNT(IOFF) NSTRII ',
     &                 IOFF, IPNT(IOFF),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT'
        WRITE(6,*) ' Required total symmetry',ISYM
        WRITE(6,*) ' Number of symmetry blocks ', NBLKS
        WRITE(6,*)
        WRITE(6,*) ' Offset array  for symmetry blocks'
        CALL IWRTMA(IPNT,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE WEIGHT_LUCI(Z,NEL,NORB1,NORB2,NORB3,
     &                       MNRS1,MXRS1,MNRS3,MXRS3,ISCR,NTEST)
*
* construct vertex weights
*
* Reverse lexical ordering is used for restricted space
*
      IMPLICIT REAL*8           ( A-H,O-Z)
      INTEGER Z(*), ISCR(*)
*
      NORB = NORB1 + NORB2 + NORB3
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' >>>> WEIGHT <<<<< '
        WRITE(6,*) ' NORB1 NORB2 NORB3 ',NORB1,NORB2,NORB3
        WRITE(6,*) ' NEL MNRS1 MXRS1 MNRS3 MXRS3 '
        WRITE(6,*)   NEL,MNRS1,MXRS1,MNRS3,MXRS3
      END IF
*
      KLFREE = 1
      KLMAX  = KLFREE
      KLFREE = KLFREE + NORB
*
      KLMIN  = KLFREE
      KLFREE = KLFREE + NORB
*
      KW     = KLFREE
      KLFREE = KW + (NEL+1)*(NORB+1)
*.Max and min arrays for strings
      CALL RSMXMN_LUCI(ISCR(KLMAX),ISCR(KLMIN),NORB1,NORB2,NORB3,
     &                 NEL,MNRS1,MXRS1,MNRS3,MXRS3,NTEST)
*. Arc weights
      CALL GRAPW_LUCI(ISCR(KW),Z(1),ISCR(KLMIN),ISCR(KLMAX),NORB,
     &                NEL,NTEST)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ZSTINF_GAS(IPRNT)
*
* Set up common block /STINF/ from information in /STINP/
*
*=========
* Input
*=========
* Information in /CGAS/ and /GASSTR/
*
*======================
* Output ( in /STINF/ )
*======================
* ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
*                    or annihilating (ISTAC(ITYP,1)) an electron
*                    from a string of type  ITYP . A zero indicates
*                    that this mapping is not included
*                    Only strings belonging to the same
*                    Orbital group are mapped
*                    mapped
*. Input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(2,MXPNGAS),IGSINA,IGSDEL
C     COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),
C    &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),
C    &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS)
*. Output
#include "stinf.inc"
*./STINF/
C     COMMON/STINF/ISTAC(MXPSTT,2),NOCTYP(MXPSTT),NSTFTP(MXPSTT),
C    &             INUMAP(MXPSTT),INDMAP(MXPSTT)
*. Only the first element, i.e. ISTAC  is defined

*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
* ******************************************************************
* Mappings between strings with the same type ISTTP index , +/- 1 el
* ******************************************************************
      CALL ISETVC(ISTAC,0,2*MXPSTT)
      DO  IGAS = 1, NGAS
*. groups for a given gas spaces goes with increasing number of orbitals,
*  so the first space does not have any creation mapping
*  and the last space does not have any annihilation mapping
*
        MGRP = NGPSTR(IGAS)
        DO IGRP = 1, MGRP
          IIGRP = IGRP + IBGPSTR(IGAS) -1
          IF(IGRP.NE.1) THEN
*. Annihilation map is present : IIGRP => IIGRP - 1
            ISTAC(IIGRP,1) = IIGRP -1
          END IF
          IF(IGRP.NE.MGRP) THEN
*. Creation map is present : IIGRP => IIGRP + 1
             ISTAC(IIGRP,2) = IIGRP + 1
          END IF
        END DO
      END DO
*
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Type - type mapping array ISTAC '
        WRITE(6,*) ' =============================== '
        CALL IWRTMA(ISTAC,NGRP  ,2,MXPSTT,2)
      END IF
*
      RETURN
      END
