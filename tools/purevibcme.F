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
      PROGRAM PUREVIB
C
C     This program calculates pure vibrational contribution to nonlinear 
C     magneto-optical processes, if requested consistently to second order 
C     in perturbation theory, but a couple of lower truncation levels are
C     possible. It is though worth noticing that as soon as one goes beyond 
C     the double-harmonic approximation, doing the full second-order 
C     perturbation corrections does not cost more than the lower second-order
C     truncation level.
C
C     The program is currently restricted to pure vibrational contributions
C     to the hypermagnetizability.
C
C     Cubic and quartic force fields (in normal modes) are read from a 
C     Dalton calculation of the first- and second numerical derivatives of
C     the molecular Hessian along normal coordinates in the vibrational 
C     averaging module.
C
C     The program expects a set of files to be listed corresponding to the 
C     order of the pure vibrational contributions calculated and the number of
C     normal modes in the molecule.
C
C     Future improvements: Include directly in Dalton for automatic
C     calculation of pure vibrational contributions using the new, general
C     numerical derivative routines of Torgeir.
C
C     Should also be extended to hyperpolarizabilities at some point.
C
C     Written in San Diego, November 1999 by Kenneth Ruud.
C
C     Parameter description
C
      PARAMETER (NFREQ = 3)
      CHARACTER*90 FILEROOT
      LOGICAL PURE, ZPVA, ALFA, BETA, GAMMA, ETA
      DIMENSION Q(3,3,-NFREQ:NFREQ,-NFREQ:NFREQ), 
     &          XMU(3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          ALPHA(3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          XI(3,3,-NFREQ:NFREQ,-NFREQ:NFREQ), 
     &          G(3,3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          F(3,3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          BETAP(3,3,3,-NFREQ:NFREQ,-NFREQ:NFREQ)
      DIMENSION OMEGA(NFREQ), CUBIC(NFREQ,NFREQ,NFREQ),
     &          QUARTIC(NFREQ,NFREQ,NFREQ)
      DIMENSION ETAD(3,3,3,3,-NFREQ:NFREQ), ETAP(3,3,3,3,-NFREQ:NFREQ),
     &          ETAD2(3,3,3,3,NFREQ), ETAP2(3,3,3,3,NFREQ),
     &          ETADZ(3,3,3,3), ETAPZ(3,3,3,3)
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ), XMU3(3,NFREQ,NFREQ),
     &          F1(3,3,3,NFREQ), F2(3,3,3,NFREQ,NFREQ), 
     &          F3(3,3,3,NFREQ,NFREQ), G1(3,3,3,NFREQ), 
     &          G2(3,3,3,NFREQ,NFREQ), G3(3,3,3,NFREQ,NFREQ),
     &          ALPHA1(3,3,NFREQ), ALPHA2(3,3,NFREQ,NFREQ),
     &          ALPHA3(3,3,NFREQ,NFREQ), XIP1(3,3,NFREQ),
     &          XIP2(3,3,NFREQ,NFREQ), XIP3(3,3,NFREQ,NFREQ),
     &          Q1(3,3,NFREQ), Q2(3,3,NFREQ,NFREQ), Q3(3,3,NFREQ,NFREQ)
      DIMENSION XMUG(3,3,3,3,5), XMUQ(3,3,3,3,2)
      DIMENSION XMUF(3,3,3,3,5), AXI(3,3,3,3,5), XMUXI(3,3,3,3,2)
      DIMENSION XMU2V(3,3,5)
      DIMENSION XMUA(3,3,3,5), XMU3V(3,3,3,2)
C
      CALL DZERO(XMUG,405)
      CALL DZERO(XMUQ,405)
      CALL DZERO(XMUF,405)
      CALL DZERO(AXI,405)
      CALL DZERO(XMUXI,405)
      CALL DZERO(XMU2V,45)
      CALL DZERO(XMUA,135)
      CALL DZERO(XMU3V,54)
      CALL DZERO(Q,9*(2*NFREQ + 1)**2)
      CALL DZERO(XMU,3*(2*NFREQ + 1)**2)
      CALL DZERO(ALPHA,9*(2*NFREQ + 1)**2)
      CALL DZERO(XI,9*(2*NFREQ + 1)**2)
      CALL DZERO(G,27*(2*NFREQ + 1)**2)
      CALL DZERO(F,27*(2*NFREQ + 1)**2)
      CALL DZERO(ETAD,81*(2*NFREQ + 1))
      CALL DZERO(ETAP,81*(2*NFREQ + 1))
C
C     Read in some specific information for this calculation
C
ckr      CALL READMOD(IORDER,STEP,FILEROOT,OPTFRQ)
ckr
C      For simplicity we now just define this in the following section
C
      STEP = 0.05D0
      PURE = .TRUE.
      ZPVA = .FALSE.
      ETA  = .TRUE.
      ALFA = .FALSE.
      BETA = .FALSE.
      GAMMA = .FALSE.
C
      IORDER = 2
ckr      FILEROOT(1:9) = 'quad_hftz'
      FILEROOT(1:10) = 'quad_castz'
ckr      FILEROOT(1:7) = 'hftzeff'
ckr      FILEROOT(1:8) = 'castzeff'
ckr      FILEROOT(1:8) = 'Zquad_hf'
      LENROOT = 10
      OPTFRQ = 0.072D0
      OFRQ1 = 0.0D0
      OFRQ2 = 0.0D0
      OFRQ3 = 0.0D0
      OFRQS = -(OFRQ1 + OFRQ2 + OFRQ3)
      INFFRQ = .FALSE.
C
C     Read in properties and calculate derivatives
C
      CALL READPP(XMU,Q,ALPHA,XI,G,F,BETAP,ETAD,ETAP,FILEROOT,LENROOT,
     &            IORDER,PURE,ZPVA,OPTFRQ,NFREQ)
      IF (PURE) THEN
         CALL DERIV(IORDER,XMU,XMU1,XMU2,XMU3,NFREQ,3,STEP)
         IF (ETA .OR. BETA .OR. GAMMA) THEN
            CALL DERIV(IORDER,ALPHA,ALPHA1,ALPHA2,ALPHA3,NFREQ,9,STEP)
         END IF
         IF (ETA) THEN
            CALL DERIV(IORDER,Q,Q1,Q2,Q3,NFREQ,9,STEP)
            CALL DERIV(IORDER,XI,XIP1,XIP2,XIP3,NFREQ,9,STEP)
            CALL DERIV(IORDER,G,G1,G2,G3,NFREQ,27,STEP)
            CALL DERIV(IORDER,F,F1,F2,F3,NFREQ,27,STEP)
         END IF
         IF (GAMMA) THEN
            CALL DERIV(IORDER,BETAP,BETA1,BETA2,BETA3,NFREQ,27,STEP)
         END IF
      ELSE IF (ZPVA) THEN
         CALL DERIV2(ETAD,ETAD2,NFREQ,81,STEP)
         CALL DERIV2(ETAP,ETAP2,NFREQ,81,STEP)
      END IF
C
C     Read in vibrational frequencies
C
      CALL READFREQ(NFREQ,OMEGA,FILEROOT,LENROOT)
C
C     Read in cubic and quartic force fields if needed
C
      IF (IORDER .GE. 1) CALL FORCEF(CUBIC,QUARTIC,FILEROOT,LENROOT,
     &                               NFREQ,IORDER)
C
C     We have done all the preprocessing we need, calculate the properties
C
      IF (PURE .AND. ETA) THEN
         CALL DRIVETA(XMU1,XMU2,XMU3,F1,F2,F3,G1,G2,G3,ALPHA1,ALPHA2,
     &                ALPHA3,XIP1,XIP2,XIP3,Q1,Q2,Q3,OMEGA,CUBIC,
     &                QUARTIC,XMUG,XMUQ,XMUF,AXI,XMUXI,NFREQ,IORDER,
     &                OPTFRQ,INFFRQ)
      ELSE IF (PURE .AND. (ALFA .OR. BETA .OR. GAMMA)) THEN
         CALL DRIVPOL(XMU1,XMU2,XMU3,ALPHA1,ALPHA2,ALPHA3,OMEGA,CUBIC,
     &                QUARTIC,XMU2V,XMUA,XMU3V,NFREQ,IORDER,OFRQS,OFRQ1,
     &                OFRQ2,OFRQ3,ALFA,BETA,GAMMA)
      ELSE IF (ZPVA) THEN
         CALL ZPVADRIV(ETAD2,ETADZ,81,OMEGA,NFREQ)
         CALL ZPVADRIV(ETAP2,ETAPZ,81,OMEGA,NFREQ)
      END IF
C
C     And the data are not of much use unless we print them
C     
      IF (PURE .AND. ETA) THEN
         CALL ETAPRI(XMUG,XMUQ,XMUF,AXI,XMUXI,IORDER,OPTFRQ)
      ELSE IF (PURE .AND. (ALFA .OR. BETA .OR. GAMMA)) THEN
         CALL POLPRI(XMU2V,XMUA,XMU3V,IORDER,OFRQS,OFRQ1,OFRQ2,OFRQ3,
     &               ALFA,BETA,GAMMA)
      ELSE IF (ZPVA) THEN
         CALL ETAZPRI(ETADZ,ETAPZ,OPTFRQ)
      END IF
C
C     Now the program can go to rest a bit after having completed a nice
C     piece of work (if the user has done his part of the bargain nicely)
C
      STOP
      END
      SUBROUTINE READFREQ(NFREQ,OMEGA,FILE,LENROOT)
C
C     Reads in the vibrational frequencies from file. The name of the
C     file does not conform with usual Dalton shell-script naming.
C     
      PARAMETER (XTKAYS = 219 474.63067 D0)
      DIMENSION OMEGA(NFREQ)
      CHARACTER*90 FILE, TMPLINE
C
      FILE(LENROOT+1:LENROOT+12) = '00.freq     '
      OPEN(UNIT=8,FILE=FILE(:LENROOT+7),FORM='FORMATTED')
      REWIND (8)
 10   CONTINUE
      READ (8,'(A90)',END=900) TMPLINE
      IF (TMPLINE(:37).EQ.'   mode   irrep     cm-1     hartrees') THEN
         READ (8,'(A)') TMPLINE
         READ (8,'(A)') TMPLINE
         IMOD = 1
 20      CONTINUE
         READ (8,'(3X,I2,10X,F11.2)') NMOD,FREQ
         IF (NMOD .EQ. 0) GOTO 20
         OMEGA(NMOD) = FREQ/XTKAYS
         IMOD = IMOD + 1
         IF (IMOD .LE. NFREQ) GO TO 20
      ELSE
         GOTO 10
      END IF
      CLOSE (8,STATUS='KEEP')
      RETURN
 900  WRITE (6,'(A)') 'No frequencies found in file. Aborting'
      STOP
      END
C      
      SUBROUTINE READPP(XMU,Q,ALPHA,XI,G,F,BETAP,ETAD,ETAP,FILE,LENROOT,
     &                  IORDER,PURE,ZPVA,OPTFRQ,NFREQ)
C
C     Read in the various quantities that are needed to generate the
C     geometrical derivatives of the properties.
C
      LOGICAL PURE, ZPVA
      DIMENSION XMU(3,-NFREQ:NFREQ,-NFREQ:NFREQ), 
     &          Q(3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          ALPHA(3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          XI(3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          G(3,3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          F(3,3,3,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          ETAD(3,3,3,3,-NFREQ:NFREQ),
     &          ETAP(3,3,3,3,-NFREQ:NFREQ),
     &          BETAP(3,3,3,-NFREQ:NFREQ,-NFREQ:NFREQ)
      CHARACTER FILE*90, CHMOD*1, CHMOD2*1
C
      FILE(LENROOT+1:LENROOT+6) = '00.out'
      OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+6))
      REWIND (8)
      IF (PURE) THEN
         CALL READNUM(XMU(1,0,0),Q(1,1,0,0),ALPHA(1,1,0,0),XI(1,1,0,0),
     &                G(1,1,1,0,0),F(1,1,1,0,0))
      ELSE IF (ZPVA) THEN
         CALL READETA(ETAD(1,1,1,1,0),ETAP(1,1,1,1,0),OPTFRQ)
      END IF
      CLOSE (UNIT=8,STATUS='KEEP')
C
C     Needs to be modified when we go beyond 9 vibrational modes
C
      DO IMOD = 1, NFREQ, 1
         CHMOD = CHAR(48 + IMOD)
         FILE(LENROOT+1:LENROOT+8) = 'm'//CHMOD//'p1.out'
         OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+8))
         REWIND (8)
         IF (PURE) THEN
            CALL READNUM(XMU(1,IMOD,0),Q(1,1,IMOD,0),ALPHA(1,1,IMOD,0),
     &                   XI(1,1,IMOD,0),G(1,1,1,IMOD,0),F(1,1,1,IMOD,0))
         ELSE IF (ZPVA) THEN
            CALL READETA(ETAD(1,1,1,1,IMOD),ETAP(1,1,1,1,IMOD),
     &                   OPTFRQ)
         END IF
         CLOSE (UNIT=8,STATUS='KEEP')
         FILE(LENROOT+1:LENROOT+8)= 'm'//CHMOD//'m1.out'
         OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+8))
         REWIND (8)
         IF (PURE) THEN
            CALL READNUM(XMU(1,-IMOD,0),Q(1,1,-IMOD,0),
     &                   ALPHA(1,1,-IMOD,0),XI(1,1,-IMOD,0),
     &                   G(1,1,1,-IMOD,0),F(1,1,1,-IMOD,0))
         ELSE IF (ZPVA) THEN
            CALL READETA(ETAD(1,1,1,1,-IMOD),ETAP(1,1,1,1,-IMOD),
     &                   OPTFRQ)
         END IF
         CLOSE (UNIT=8,STATUS='KEEP')
         IF (PURE .AND. IORDER .GE. 1) THEN
            FILE(LENROOT+1:LENROOT+8) = 'm'//CHMOD//'p2.out'
            OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+8))
            REWIND (8)
            CALL READNUM(XMU(1,IMOD,IMOD),Q(1,1,IMOD,IMOD),
     &                   ALPHA(1,1,IMOD,IMOD),XI(1,1,IMOD,IMOD),
     &                   G(1,1,1,IMOD,IMOD),F(1,1,1,IMOD,IMOD))
            CLOSE (UNIT=8,STATUS='KEEP')
            FILE(LENROOT+1:LENROOT+8) = 'm'//CHMOD//'m2.out'
            OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+8))
            REWIND (8)
            CALL READNUM(XMU(1,-IMOD,-IMOD),Q(1,1,-IMOD,-IMOD),
     &                   ALPHA(1,1,-IMOD,-IMOD),XI(1,1,-IMOD,-IMOD),
     &                   G(1,1,1,-IMOD,-IMOD),F(1,1,1,-IMOD,-IMOD))
            CLOSE (UNIT=8,STATUS='KEEP')
         DO JMOD = IMOD + 1, NFREQ, 1
            CHMOD2 = CHAR(48 + JMOD)
            FILE(LENROOT+1:LENROOT+12) = 
     &           'm'//CHMOD//'m'//CHMOD2//'p1p1.out'
            OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+12))
            REWIND (8)
            CALL READNUM(XMU(1,IMOD,JMOD),Q(1,1,IMOD,JMOD),
     &                   ALPHA(1,1,IMOD,JMOD),XI(1,1,IMOD,JMOD),
     &                   G(1,1,1,IMOD,JMOD),F(1,1,1,IMOD,JMOD))
            REWIND (8)
            CALL READNUM(XMU(1,JMOD,IMOD),Q(1,1,JMOD,IMOD),
     &                   ALPHA(1,1,JMOD,IMOD),XI(1,1,JMOD,IMOD),
     &                   G(1,1,1,JMOD,IMOD),F(1,1,1,JMOD,IMOD))
            CLOSE (UNIT=8,STATUS='KEEP')
            FILE(LENROOT+1:LENROOT+12) = 
     &           'm'//CHMOD//'m'//CHMOD2//'m1p1.out'
            OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+12))
            REWIND (8)
            CALL READNUM(XMU(1,-IMOD,JMOD),Q(1,1,-IMOD,JMOD),
     &                   ALPHA(1,1,-IMOD,JMOD),XI(1,1,-IMOD,JMOD),
     &                   G(1,1,1,-IMOD,JMOD),F(1,1,1,-IMOD,JMOD))
            REWIND (8)
            CALL READNUM(XMU(1,JMOD,-IMOD),Q(1,1,JMOD,-IMOD),
     &                   ALPHA(1,1,JMOD,-IMOD),XI(1,1,JMOD,-IMOD),
     &                   G(1,1,1,JMOD,-IMOD),F(1,1,1,JMOD,-IMOD))
            CLOSE (UNIT=8,STATUS='KEEP')
            FILE(LENROOT+1:LENROOT+12) = 
     &           'm'//CHMOD//'m'//CHMOD2//'p1m1.out'
            OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+12))
            REWIND (8)
            CALL READNUM(XMU(1,IMOD,-JMOD),Q(1,1,IMOD,-JMOD),
     &                   ALPHA(1,1,IMOD,-JMOD),XI(1,1,IMOD,-JMOD),
     &                   G(1,1,1,IMOD,-JMOD),F(1,1,1,IMOD,-JMOD))
            REWIND (8)
            CALL READNUM(XMU(1,-JMOD,IMOD),Q(1,1,-JMOD,IMOD),
     &                   ALPHA(1,1,-JMOD,IMOD),XI(1,1,-JMOD,IMOD),
     &                   G(1,1,1,-JMOD,IMOD),F(1,1,1,-JMOD,IMOD))
            CLOSE (UNIT=8,STATUS='KEEP')
            FILE(LENROOT+1:LENROOT+12) = 
     &           'm'//CHMOD//'m'//CHMOD2//'m1m1.out'
            OPEN (UNIT=8,FORM='FORMATTED',NAME=FILE(:LENROOT+12))
            REWIND (8)
            CALL READNUM(XMU(1,-IMOD,-JMOD),Q(1,1,-IMOD,-JMOD),
     &                   ALPHA(1,1,-IMOD,-JMOD),XI(1,1,-IMOD,-JMOD),
     &                   G(1,1,1,-IMOD,-JMOD),F(1,1,1,-IMOD,-JMOD))
            REWIND (8)
            CALL READNUM(XMU(1,-JMOD,-IMOD),Q(1,1,-JMOD,-IMOD),
     &                   ALPHA(1,1,-JMOD,-IMOD),XI(1,1,-JMOD,-IMOD),
     &                   G(1,1,1,-JMOD,-IMOD),F(1,1,1,-JMOD,-IMOD))
            CLOSE (UNIT=8,STATUS='KEEP')
         END DO
         END IF
      END DO
      RETURN
      END
C
      SUBROUTINE READNUM(XMU,Q,ALPHA,XI,G,F)
C
C     Search in the file for the numbers we need
C
      CHARACTER TMPLINE*90, TXTOP*8
      DIMENSION XMU(3), Q(3,3), ALPHA(3,3), XI(3,3), G(3,3,3),
     &          F(3,3,3)
C
 10   CONTINUE
      READ (8,'(A90)') TMPLINE
C
C     Expectation values?
C     
      IF (TMPLINE(6:19) .EQ. 'XDIPLEN  TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') XMU(1)
      ELSE IF (TMPLINE(6:19) .EQ. 'YDIPLEN  TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') XMU(2)
      ELSE IF (TMPLINE(6:19) .EQ. 'ZDIPLEN  TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') XMU(3)
      ELSE IF (TMPLINE(6:19) .EQ. 'XXQUADRU TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') Q(1,1)
      ELSE IF (TMPLINE(6:19) .EQ. 'XYQUADRU TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') Q(1,2)
         Q(2,1) = Q(1,2)
      ELSE IF (TMPLINE(6:19) .EQ. 'XZQUADRU TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') Q(1,3)
         Q(3,1) = Q(1,3)
      ELSE IF (TMPLINE(6:19) .EQ. 'YYQUADRU TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') Q(2,2)
      ELSE IF (TMPLINE(6:19) .EQ. 'YZQUADRU TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') Q(2,3)
         Q(3,2) = Q(2,3)
      ELSE IF (TMPLINE(6:19) .EQ. 'ZZQUADRU TOTAL') THEN
         READ (TMPLINE,'(28X,F15.8)') Q(3,3)
C
C     Any linear response function values?
C
      ELSE IF (TMPLINE(35:65) .EQ. '<<XXQUADRU;XDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<XDIPLEN ;XXQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(1,1,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XXQUADRU;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;XXQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(2,1,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XXQUADRU;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;XXQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(3,1,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XYQUADRU;XDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<XDIPLEN ;XYQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(1,2,1)
         G(1,1,2) = G(1,2,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XYQUADRU;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;XYQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(2,2,1)
         G(2,1,2) = G(2,2,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XYQUADRU;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;XYQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(3,2,1)
         G(3,1,2) = G(3,2,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XZQUADRU;XDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<XDIPLEN ;XZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(1,3,1)
         G(1,1,3) = G(1,3,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XZQUADRU;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;XZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(2,3,1)
         G(2,1,3) = G(2,3,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XZQUADRU;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;XZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(3,3,1)
         G(3,1,3) = G(3,3,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YYQUADRU;XDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<XDIPLEN ;YYQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(1,2,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YYQUADRU;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;YYQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(2,2,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YYQUADRU;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;YYQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(3,2,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YZQUADRU;XDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<XDIPLEN ;YZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(1,2,3)
         G(1,3,2) = G(1,2,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YZQUADRU;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;YZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(2,2,3)
         G(2,3,2) = G(2,2,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YZQUADRU;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;YZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(3,2,3)
         G(3,3,2) = G(3,2,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<ZZQUADRU;XDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<XDIPLEN ;ZZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(1,3,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<ZZQUADRU;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;ZZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(2,3,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<ZZQUADRU;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;ZZQUADRU>> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') G(3,3,3)
C
      ELSE IF (TMPLINE(35:65) .EQ. '<<XDIPLEN ;XDIPLEN >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') ALPHA(1,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XDIPLEN ;YDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YDIPLEN ;XDIPLEN >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') ALPHA(1,2)
         ALPHA(2,1) = ALPHA(1,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XDIPLEN ;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;XDIPLEN >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') ALPHA(1,3)
         ALPHA(3,1) = ALPHA(1,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YDIPLEN ;YDIPLEN >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') ALPHA(2,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YDIPLEN ;ZDIPLEN >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZDIPLEN ;YDIPLEN >> ( 0.0000)')
     &        THEN 
         READ (TMPLINE,'(68X,F20.12)') ALPHA(3,2)
         ALPHA(2,3) = ALPHA(3,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<ZDIPLEN ;ZDIPLEN >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') ALPHA(3,3)
C
      ELSE IF (TMPLINE(35:65) .EQ. '<<XANGMOM ;XANGMOM >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') XI(1,1)
         XI(1,1) = XI(1,1)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XANGMOM ;YANGMOM >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<YANGMOM ;XANGMOM >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') XI(1,2)
         XI(1,2) = XI(1,2)
         XI(2,1) = XI(1,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<XANGMOM ;ZANGMOM >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZANGMOM ;XANGMOM >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') XI(1,3)
         XI(3,1) = XI(1,3)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YANGMOM ;YANGMOM >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') XI(2,2)
         XI(2,2) = XI(2,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<YANGMOM ;ZANGMOM >> ( 0.0000)'
     &    .OR. TMPLINE(35:65) .EQ. '<<ZANGMOM ;YANGMOM >> ( 0.0000)')
     &        THEN 
         READ (TMPLINE,'(68X,F20.12)') XI(3,2)
         XI(3,2) = XI(3,2)
         XI(2,3) = XI(3,2)
      ELSE IF (TMPLINE(35:65) .EQ. '<<ZANGMOM ;ZANGMOM >> ( 0.0000)')
     &        THEN
         READ (TMPLINE,'(68X,F20.12)') XI(3,3)
         XI(3,3) = XI(3,3)
      ELSE IF (TMPLINE(1:47) .EQ. 
     &        '@ Quadratic response function value in a.u. for') THEN
         GOTO 90
      END IF
      GOTO 10
C
C     We have done all expectation values and linear response functions.
C     Now we turn to the quadratic response functions.
C     
 90   CONTINUE
      READ (8,'(30X,A8)') TXTOP
      IF (TXTOP .EQ. 'XDIPLEN ') THEN
         IDIP = 1
      ELSE IF (TXTOP .EQ. 'YDIPLEN ') THEN
         IDIP = 2
      ELSE IF (TXTOP .EQ. 'ZDIPLEN ') THEN
         IDIP = 3
      ELSE
         GOTO 150
      END IF
      READ (8,'(30X,A8)') TXTOP
      IF (TXTOP .EQ. 'XANGMOM ') THEN
         IANG1 = 1
      ELSE IF (TXTOP .EQ. 'YANGMOM ') THEN
         IANG1 = 2
      ELSE IF (TXTOP .EQ. 'ZANGMOM ') THEN
         IANG1 = 3
      ELSE
         GOTO 150
      END IF
      READ (8,'(30X,A8)') TXTOP
      IF (TXTOP .EQ. 'XANGMOM ') THEN
         IANG2 = 1
      ELSE IF (TXTOP .EQ. 'YANGMOM ') THEN
         IANG2 = 2
      ELSE IF (TXTOP .EQ. 'ZANGMOM ') THEN
         IANG2 = 3
      ELSE
         GOTO 150
      END IF
      READ (8,'(A)') TMPLINE
      READ (8,'(29X,3F15.8)') FREQA, FREQB, VALUE
      IF (FREQA .EQ. 0.0D0 .AND. FREQB .EQ. 0.0D0) THEN
         F(IDIP,IANG1,IANG2) = VALUE
         F(IDIP,IANG2,IANG1) = VALUE
      END IF
 150  CONTINUE
      READ (8,'(A)',END=900) TMPLINE
      IF (TMPLINE(1:47) .EQ. 
     &        '@ Quadratic response function value in a.u. for') THEN
         GOTO 90
      ELSE
         GOTO 150
      END IF
 900  CONTINUE
      DO I = 1, 9
         XI(I,1) = 0.25D0*XI(I,1)
      END DO
      DO I = 1, 27
         G(I,1,1) = G(I,1,1)
         F(I,1,1) = 0.25D0*F(I,1,1)
      END DO
C
C     All input read from this file
C
      RETURN
      END
C
      SUBROUTINE READETA(ETAD,ETAP,OPTFRQ)
C
C     Search in the file for the numbers we need
C
      CHARACTER TMPLINE*90, TXTOP*8
      DIMENSION ETAD(3,3,3,3), ETAP(3,3,3,3)
C
 10   CONTINUE
      READ (8,'(A90)') TMPLINE
      IF (TMPLINE(1:47) .NE. 
     &        '@ Quadratic response function value in a.u. for') GOTO 10
C
C     Diamagnetic contribution
C     
 90   CONTINUE
      READ (8,'(30X,A8)') TXTOP
      IF (TXTOP .EQ. 'XDIPLEN ') THEN
         IDIP1 = 1
      ELSE IF (TXTOP .EQ. 'YDIPLEN ') THEN
         IDIP1 = 2
      ELSE IF (TXTOP .EQ. 'ZDIPLEN ') THEN
         IDIP1 = 3
      ELSE
         GOTO 150
      END IF
      READ (8,'(30X,A8)') TXTOP
      IF (TXTOP .EQ. 'XDIPLEN ') THEN
         IDIP2 = 1
      ELSE IF (TXTOP .EQ. 'YDIPLEN ') THEN
         IDIP2 = 2
      ELSE IF (TXTOP .EQ. 'ZDIPLEN ') THEN
         IDIP2 = 3
      ELSE
         GOTO 150
      END IF
      READ (8,'(30X,A8)') TXTOP
      IF (TXTOP .EQ. 'XXQUADRU') THEN
         IANG1 = 1
         IANG2 = 1
      ELSE IF (TXTOP .EQ. 'XYQUADRU') THEN
         IANG1 = 1
         IANG2 = 2
      ELSE IF (TXTOP .EQ. 'XZQUADRU') THEN
         IANG1 = 1
         IANG2 = 3
      ELSE IF (TXTOP .EQ. 'YYQUADRU') THEN
         IANG1 = 2
         IANG2 = 2
      ELSE IF (TXTOP .EQ. 'YZQUADRU') THEN
         IANG1 = 2
         IANG2 = 3
      ELSE IF (TXTOP .EQ. 'ZZQUADRU') THEN
         IANG1 = 3
         IANG2 = 3
      ELSE
         GOTO 150
      END IF
      READ (8,'(A)') TMPLINE
      READ (8,'(A)') TMPLINE
      IF (TMPLINE(1:9) .EQ. ' is equal') THEN
         READ (8,'(A)') TMPLINE
         READ (8,'(A)') TMPLINE
         GOTO 150
      ELSE
         READ (TMPLINE,'(29X,3F15.8)') FREQA, FREQB, VALUE
         IF (FREQA .EQ. OPTFRQ .AND. FREQB .EQ. 0.0D0) THEN
            ETAD(IDIP1,IDIP2,IANG1,IANG2) = VALUE
            ETAD(IDIP1,IDIP2,IANG2,IANG1) = VALUE
            IF (OPTFRQ .EQ. 0.0D0) THEN
               ETAD(IDIP2,IDIP1,IANG1,IANG2) = VALUE
               ETAD(IDIP2,IDIP1,IANG2,IANG1) = VALUE
            END IF
         END IF
      END IF
 150  CONTINUE
      READ (8,'(A)',END=900) TMPLINE
      IF (TMPLINE(1:47) .EQ. 
     &        '@ Quadratic response function value in a.u. for') THEN
         GOTO 90
      ELSE IF (TMPLINE(1:30) .EQ. 
     &        ' Cubic response function value') THEN
         GOTO 200
      ELSE
         GOTO 150
      END IF
C
C     Paramagnetic contribution
C
 200  CONTINUE
      READ (8,'(36X,A8,6X,F10.6)') TXTOP,FREQA
      IF (TXTOP .EQ. 'XDIPLEN ') THEN
         IDIP1 = 1
      ELSE IF (TXTOP .EQ. 'YDIPLEN ') THEN
         IDIP1 = 2
      ELSE IF (TXTOP .EQ. 'ZDIPLEN ') THEN
         IDIP1 = 3
      ELSE
         GOTO 250
      END IF
      READ (8,'(36X,A8,6X,F10.6)') TXTOP,FREQB
      IF (TXTOP .EQ. 'XDIPLEN ') THEN
         IDIP2 = 1
      ELSE IF (TXTOP .EQ. 'YDIPLEN ') THEN
         IDIP2 = 2
      ELSE IF (TXTOP .EQ. 'ZDIPLEN ') THEN
         IDIP2 = 3
      ELSE
         GOTO 250
      END IF
      READ (8,'(36X,A8,6X,F10.6)') TXTOP,FREQC
      IF (TXTOP .EQ. 'XANGMOM ') THEN
         IANG1 = 1
      ELSE IF (TXTOP .EQ. 'YANGMOM ') THEN
         IANG1 = 2
      ELSE IF (TXTOP .EQ. 'ZANGMOM ') THEN
         IANG1 = 3
      ELSE
         GOTO 250
      END IF
      READ (8,'(36X,A8,6X,F10.6)') TXTOP,FREQD
      IF (TXTOP .EQ. 'XANGMOM ') THEN
         IANG2 = 1
      ELSE IF (TXTOP .EQ. 'YANGMOM ') THEN
         IANG2 = 2
      ELSE IF (TXTOP .EQ. 'ZANGMOM ') THEN
         IANG2 = 3
      ELSE
         GOTO 250
      END IF
      READ (8,'(A)') TMPLINE
      READ (8,'(A)') TMPLINE
      IF (TMPLINE(1:9) .EQ. ' is equal') THEN
         READ (8,'(A)') TMPLINE
         READ (8,'(A)') TMPLINE
         GOTO 250
      ELSE
         READ (TMPLINE,'(21X,F20.8)') VALUE
         IF (FREQB .EQ. OPTFRQ .AND. FREQC .EQ. 0.0D0 .AND.
     &       FREQD .EQ. 0.0D0) THEN
            ETAP(IDIP1,IDIP2,IANG1,IANG2) = VALUE
            ETAP(IDIP1,IDIP2,IANG2,IANG1) = VALUE
            IF (OPTFRQ .EQ. 0.0D0) THEN
               ETAP(IDIP2,IDIP1,IANG1,IANG2) = VALUE
               ETAP(IDIP2,IDIP1,IANG2,IANG1) = VALUE
            END IF
         END IF
      END IF
 250  CONTINUE
      READ (8,'(A)',END=900) TMPLINE
      IF (TMPLINE(1:30) .EQ. 
     &        ' Cubic response function value') THEN
         GOTO 200
      ELSE
         GOTO 150
      END IF
 900  CONTINUE
C
      DO I = 1, 81
         ETAD(I,1,1,1) = - ETAD(I,1,1,1)
      END DO
      DO I = 1, 81
         ETAP(I,1,1,1) = 0.25D0*ETAP(I,1,1,1)
      END DO
C
C     All input read from this file
C
      RETURN
      END
C
      SUBROUTINE FORCEF(CUBIC,QUARTIC,FILE,LENROOT,NFREQ,IORDER)
C
C     Read in the full cubic and partial quartic force fields needed in case
C     we want to go beyond the double-harmonic approximation
C
      DIMENSION CUBIC(NFREQ,NFREQ,NFREQ), QUARTIC(NFREQ,NFREQ,NFREQ)
      CHARACTER*90 FILE, TMPLINE
C
      FILE(LENROOT+1:LENROOT+18) = '00.force          '
      OPEN (UNIT=8,FILE=FILE(:LENROOT+8),FORM='FORMATTED')
 10   CONTINUE
      READ (8,'(A90)') TMPLINE
      IF (TMPLINE(:27) .EQ. ' Anharmonic force constants') THEN
         READ (TMPLINE,'(36X,I5)') IC
 20      CONTINUE
         READ (8,'(A90)') TMPLINE
         IF (TMPLINE(16:21) .EQ. 'Column') THEN
            READ (TMPLINE,'(21X,I4)') ICOL
 30         CONTINUE
            READ (8,'(1X,I7,2X,4D15.6)') IROW,VAL1,VAL2,VAL3,VAL4
            IF (IROW .EQ. 0) GO TO 20
            CUBIC(IROW,ICOL,IC) = VAL1
            IF (NFREQ .GT. 1) THEN
               CUBIC(IROW,ICOL + 1,IC) = VAL2
               IF (NFREQ .GT. 2) THEN
                  CUBIC(IROW,ICOL + 2,IC) = VAL3
                  IF (NFREQ .GT. 3) THEN
                     CUBIC(IROW,ICOL + 3,IC) = VAL4
                  END IF
               END IF
            END IF
            GO TO 30
         ELSE IF (TMPLINE(:27) .EQ. ' Anharmonic force constants') THEN
            READ (TMPLINE,'(36X,I5)') IC
            GO TO 20
         ELSE IF (TMPLINE(2:35) .EQ. 
     &        'Numerical differentiation complete' .OR.
     &        TMPLINE(2:31) .EQ. 'Quartic force constants V_aabc') THEN
            GOTO 40
         ELSE
            GOTO 20
         END IF
      ELSE
         GOTO 10
      END IF
 40   CONTINUE
      IF (IORDER .GE. 2 .AND. TMPLINE(2:35) .EQ. 
     &    'Numerical differentiation complete') THEN
         WRITE (6,'(/2A)') 'ERROR: No quartic force field found in '//
     &        'file ',FILE(:LENROOT+8)
         STOP
      END IF
      IF (IORDER .GE. 2) THEN
         READ (TMPLINE,'(41X,I5)') IC
 50      CONTINUE
         READ (8,'(A90)') TMPLINE
         IF (TMPLINE(16:21) .EQ. 'Column') THEN
            READ (TMPLINE,'(21X,I4)') ICOL
 60         CONTINUE
            READ (8,'(1X,I7,2X,4D15.6)') IROW,VAL1,VAL2,VAL3,VAL4
            IF (IROW .EQ. 0) GO TO 50
            QUARTIC(IROW,ICOL,IC) = VAL1
            IF (NFREQ .GT. 1) THEN
               QUARTIC(IROW,ICOL + 1,IC) = VAL2
               IF (NFREQ .GT. 2) THEN
                  QUARTIC(IROW,ICOL + 2,IC) = VAL3
                  IF (NFREQ .GT. 3) THEN
                     QUARTIC(IROW,ICOL + 3,IC) = VAL4
                  END IF
               END IF
            END IF
            GO TO 60
         ELSE IF (TMPLINE(2:35) .EQ. 
     &           'Numerical differentiation complete') THEN
            GOTO 70
         ELSE IF (TMPLINE(2:31) .EQ. 'Quartic force constants V_aabc')
     &           THEN
            READ (TMPLINE,'(41X,I5)') IC
            GOTO 50
         ELSE
            GOTO 50
         END IF
      END IF
 70   CONTINUE
      DO I = 1, NFREQ*NFREQ*NFREQ
         QUARTIC(I,1,1) = QUARTIC(I,1,1)*0.50D0
      END DO
      RETURN
      END
      SUBROUTINE ZPVADRIV(PROP2,ZPVAP,NDIM,OMEGA,NFREQ)
C
C     Calculates the zero-point vibrational contributions to eta
C     according to Aastrand, Ruud and Taylor, JCP 112 (February 2000)
C
      PARAMETER (XFAMU = 1822.88852 D0)
      DIMENSION PROP2(NDIM,NFREQ), ZPVAP(NDIM), OMEGA(NFREQ)
C
      CALL DZERO(ZPVAP,NDIM)
      DO I = 1, NDIM
         DO J = 1, NFREQ
ckr           ZPVAP(I) = ZPVAP(I) + PROP2(I,J)/(4.0D0*SQRT(OMEGA(J))*XFAMU)
           ZPVAP(I) = ZPVAP(I) + PROP2(I,J)/(4.0D0*SQRT(OMEGA(J)))
        END DO
      END DO
      RETURN
      END
      SUBROUTINE DRIVPOL(XMU1,XMU2,XMU3,ALPHA1,ALPHA2,ALPHA3,OMEGA,
     &                   CUBIC,QUARTIC,XMU2V,XMUA,XMU3V,NFREQ,IORDER,
     &                   OFRQS,OFRQ1,OFRQ2,OFRQ3,ALFA,BETA,GAMMA)
C
C     Routine that drives the calculations of pure vibrational contributions
C     to (hyper)-polarizabilities
C     Order 0: Double-harmonic
C     Order 1: Second-order in perturbation theory (II')
C     Order 2: Full second-order including third derivatives as well as
C              contributions from the quartic force field
C
      LOGICAL ALFA, BETA, GAMMA
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ), XMU3(3,NFREQ,NFREQ),
     &          ALPHA1(3,3,NFREQ), ALPHA2(3,3,NFREQ,NFREQ),
     &          ALPHA3(3,3,NFREQ,NFREQ)
      DIMENSION OMEGA(NFREQ), CUBIC(NFREQ,NFREQ,NFREQ),
     &          QUARTIC(NFREQ,NFREQ,NFREQ)
      DIMENSION XMU2V(3,3,5), XMUA(3,3,3,5), XMU3V(3,3,3,2)
C
C     Polarizabilities. There are only [mu2] contributions
C
      IF (ALFA) THEN
         CALL MU200(XMU1,XMU2V(1,1,1),OMEGA,NFREQ,OFRQS)
         IF (IORDER .GE. 1) THEN
            CALL MU220(XMU1,XMU2,XMU3,XMU2V(1,1,2),
     &                 XMU2V(1,1,3),OMEGA,NFREQ,OFRQS)
            CALL MU211(XMU1,XMU2,CUBIC,XMU2V(1,1,4),OMEGA,NFREQ,OFRQS)
            IF (IORDER .EQ. 2) THEN
               CALL MU202(XMU1,CUBIC,QUARTIC,XMU2V(1,1,5),OMEGA,
     &                    NFREQ,OFRQS)
            END IF
         END IF
      END IF
C
C     First hyperpolarizability. Two contributions: [mu/alpha] and [mu3]
C
      IF (BETA) THEN
         CALL MUA00(XMU1,ALPHA1,XMUA(1,1,1,1),OMEGA,NFREQ,OFRQS,
     &              OFRQ1,OFRQ2)
         IF (IORDER .GE. 1) THEN
            CALL MUA20(XMU1,XMU2,XMU3,ALPHA1,ALPHA2,ALPHA3,
     &                 XMUA(1,1,1,2),XMUA(1,1,1,3),OMEGA,NFREQ,OFRQS,
     &                 OFRQ1,OFRQ2)
            CALL MUA11(XMU1,XMU2,ALPHA1,ALPHA2,CUBIC,XMUA(1,1,1,4),
     &                 OMEGA,NFREQ,OFRQS,OFRQ1,OFRQ2)
            CALL MU310(XMU1,XMU2,XMU3V(1,1,1,1),OMEGA,NFREQ,OFRQS,OFRQ1,
     &                 OFRQ2)
            CALL MU301(XMU1,CUBIC,XMU3V(1,1,1,2),OMEGA,NFREQ,OFRQS,
     &                 OFRQ1,OFRQ2)
            IF (IORDER .EQ. 2) THEN
               CALL MUA02(XMU1,ALPHA1,CUBIC,QUARTIC,XMUA(1,1,1,5),OMEGA,
     &                    NFREQ,OFRQS,OFRQ1,OFRQ2)
            END IF
         END IF
      END IF
C
C     Second hyperpolarizability. A lot of contributions
C
ckr      IF (GAMMA) THEN
ckr      END IF
C
      RETURN
      END
C
      SUBROUTINE DRIVETA(XMU1,XMU2,XMU3,F1,F2,F3,G1,G2,G3,ALPHA1,ALPHA2,
     &                   ALPHA3,XIP1,XIP2,XIP3,Q1,Q2,Q3,OMEGA,CUBIC,
     &                   QUARTIC,XMUG,XMUQ,XMUF,AXI,XMUXI,NFREQ,IORDER,
     &                   OPTFRQ,INFFRQ)
C
C     Routine that drives the calculations of pure vibrational contributions
C     to hypermagnetizabilities
C     Order 0: Double-harmonic
C     Order 1: Second-order in perturbation theory (II')
C     Order 2: Full second-order including third derivatives as well as
C              contributions from the quartic force field
C
      LOGICAL INFFRQ
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ), XMU3(3,NFREQ,NFREQ),
     &          F1(3,3,3,NFREQ), F2(3,3,3,NFREQ,NFREQ), 
     &          F3(3,3,3,NFREQ,NFREQ), G1(3,3,3,NFREQ), 
     &          G2(3,3,3,NFREQ,NFREQ), G3(3,3,3,NFREQ,NFREQ),
     &          ALPHA1(3,3,NFREQ), ALPHA2(3,3,NFREQ,NFREQ),
     &          ALPHA3(3,3,NFREQ,NFREQ), XIP1(3,3,NFREQ),
     &          XIP2(3,3,NFREQ,NFREQ), XIP3(3,3,NFREQ,NFREQ),
     &          Q1(3,3,NFREQ), Q2(3,3,NFREQ,NFREQ), Q3(3,3,NFREQ,NFREQ)
      DIMENSION OMEGA(NFREQ), CUBIC(NFREQ,NFREQ,NFREQ),
     &          QUARTIC(NFREQ,NFREQ,NFREQ)
      DIMENSION XMUG(3,3,3,3,5), XMUQ(3,3,3,3,2) 
      DIMENSION XMUF(3,3,3,3,5), AXI(3,3,3,3,5), XMUXI(3,3,3,3,2)
C
C     For the diamagnetic part we start with XMUG aka alphaQ
C     Calculate it as the sum of XMUF and AXI
C
      IF (.NOT. INFFRQ) THEN
         CALL MUF00(XMU1,G1,XMUG(1,1,1,1,1),OMEGA,NFREQ,OPTFRQ)
         IF (IORDER .GE. 1) THEN
            CALL MUF20(XMU1,G1,XMU2,G2,XMU3,G3,XMUG(1,1,1,1,2),
     &                 XMUG(1,1,1,1,3),OMEGA,NFREQ,OPTFRQ)
            CALL MUF11(XMU1,G1,XMU2,G2,CUBIC,XMUG(1,1,1,1,4),OMEGA,
     &                 NFREQ,OPTFRQ)
            IF (IORDER .EQ. 2) THEN
               CALL MUF02(XMU1,G1,CUBIC,QUARTIC,XMUG(1,1,1,1,5),OMEGA,
     &                    NFREQ,OPTFRQ)
            END IF
         END IF
      END IF
      CALL XAXI00(ALPHA1,Q1,XMUG(1,1,1,1,1),OMEGA,NFREQ)
      IF (IORDER .GE. 1) THEN
         CALL XAXI20(ALPHA1,Q1,ALPHA2,Q2,ALPHA3,Q3,XMUG(1,1,1,1,2),
     &              XMUG(1,1,1,1,3),OMEGA,NFREQ)
         CALL XAXI11(ALPHA1,Q1,ALPHA2,Q2,CUBIC,OMEGA,NFREQ,
     &              XMUG(1,1,1,1,4))
         IF (IORDER .EQ. 2) THEN
            CALL XAXI02(ALPHA1,Q1,CUBIC,QUARTIC,XMUG(1,1,1,1,5),
     &                 OMEGA,NFREQ)
         END IF
      END IF
C
C     Second contribution is XMUQ
C
      IF (IORDER .GE. 1 .AND. .NOT. INFFRQ) THEN
         CALL MUQ10(XMU1,Q1,XMU2,Q2,XMUQ(1,1,1,1,1),OMEGA,NFREQ,OPTFRQ)
         CALL MUQ01(XMU1,Q1,CUBIC,XMUQ(1,1,1,1,2),OMEGA,NFREQ,OPTFRQ)
      END IF
C
C     For the paramagnetic part we start with uF
C     
      IF (.NOT. INFFRQ) THEN
         CALL MUF00(XMU1,F1,XMUF(1,1,1,1,1),OMEGA,NFREQ,OPTFRQ)
         IF (IORDER .GE. 1) THEN
            CALL MUF20(XMU1,F1,XMU2,F2,XMU3,F3,XMUF(1,1,1,1,2),
     &                  XMUF(1,1,1,1,3),OMEGA,NFREQ,OPTFRQ)
            CALL MUF11(XMU1,F1,XMU2,F2,CUBIC,XMUF(1,1,1,1,4),OMEGA,
     &                 NFREQ,OPTFRQ)
            IF (IORDER .EQ. 2) THEN
               CALL MUF02(XMU1,F1,CUBIC,QUARTIC,XMUF(1,1,1,1,5),OMEGA,
     &                    NFREQ,OPTFRQ)
            END IF
         END IF
      END IF
C
C     AXI is the next contributions
C
      CALL XAXI00(ALPHA1,XIP1,AXI(1,1,1,1,1),OMEGA,NFREQ)
      IF (IORDER .GE. 1) THEN
         CALL XAXI20(ALPHA1,XIP1,ALPHA2,XIP2,ALPHA3,XIP3,AXI(1,1,1,1,2),
     &              AXI(1,1,1,1,3),OMEGA,NFREQ)
         CALL XAXI11(ALPHA1,XIP1,ALPHA2,XIP2,CUBIC,OMEGA,NFREQ,
     &              AXI(1,1,1,1,4))
         IF (IORDER .EQ. 2) THEN
            CALL XAXI02(ALPHA1,XIP1,CUBIC,QUARTIC,AXI(1,1,1,1,5),
     &                 OMEGA,NFREQ)
         END IF
      END IF
C
C     Finally, we have XMUxi
C
      IF (IORDER .GE. 1 .AND. .NOT. INFFRQ) THEN
         CALL MUXI10(XMU1,XIP1,XMU2,XIP2,XMUXI(1,1,1,1,1),OMEGA,NFREQ,
     &               OPTFRQ)
         CALL MUXI01(XMU1,XIP1,CUBIC,XMUXI(1,1,1,1,2),OMEGA,NFREQ,
     &                OPTFRQ)
      END IF
      RETURN
      END
C
      SUBROUTINE DERIV(IORDER,PROPPP,DERIV1,DERIV2,DERIV3,
     &                 NFREQ,NDIM,STEP)
C
C     Calculate derivatives of the molecular properties
C
      PARAMETER (XFAMU = 1822.88852 D0)
      DIMENSION PROPPP(NDIM,-NFREQ:NFREQ,-NFREQ:NFREQ),
     &          DERIV1(NDIM,NFREQ), DERIV2(NDIM,NFREQ,NFREQ),
     &          DERIV3(NDIM,NFREQ,NFREQ)
C
      IF (IORDER .GE. 0) THEN
         CALL DZERO(DERIV1,NDIM*NFREQ)
         DO IDIM = 1, NDIM
            DO IFREQ = 1, NFREQ
               DERIV1(IDIM,IFREQ) = (PROPPP(IDIM,IFREQ,0)
     &                            -  PROPPP(IDIM,-IFREQ,0))/
     &                              (2.0D0*STEP*SQRT(XFAMU))
            END DO
         END DO
      END IF
      IF (IORDER .GE. 1) THEN
         CALL DZERO(DERIV2,NDIM*NFREQ*NFREQ)
         DO IDIM = 1, NDIM
            DO IFREQA = 1, NFREQ
               DERIV2(IDIM,IFREQA,IFREQA) = (PROPPP(IDIM,IFREQA,0) +
     &              PROPPP(IDIM,-IFREQA,0) - 2.0D0*PROPPP(IDIM,0,0))/
     &              (2.0D0*(STEP**2)*XFAMU)
            DO IFREQB = 1, IFREQA - 1
               DERIV2(IDIM,IFREQA,IFREQB) = (PROPPP(IDIM,IFREQA,IFREQB)+
     &              PROPPP(IDIM,-IFREQA,-IFREQB) - 
     &              PROPPP(IDIM,-IFREQA,IFREQB) -
     &              PROPPP(IDIM,IFREQA,-IFREQB))/(8.0D0*XFAMU*STEP**2)
               DERIV2(IDIM,IFREQB,IFREQA) = DERIV2(IDIM,IFREQA,IFREQB)
            END DO
            END DO
         END DO
      END IF
      IF (IORDER .GE. 2) THEN
         CALL DZERO(DERIV3,NDIM*NFREQ*NFREQ)
         DO IDIM = 1, NDIM
            DO IFREQA = 1, NFREQ
            DO IFREQB = 1, IFREQA
               IF (IFREQA .EQ. IFREQB) THEN
                 DERIV3(IDIM,IFREQA,IFREQA)=(PROPPP(IDIM,IFREQA,IFREQA)-
     &                 PROPPP(IDIM,-IFREQA,-IFREQA) + 
     &                 2.0D0*(PROPPP(IDIM,-IFREQA,0) - 
     &                 PROPPP(IDIM,IFREQA,0)))/(12.0D0*XFAMU*SQRT(XFAMU)
     &                 *STEP**3)
               ELSE 
                 DERIV3(IDIM,IFREQA,IFREQB)=(PROPPP(IDIM,IFREQA,IFREQB)-
     &                 PROPPP(IDIM,-IFREQA,-IFREQB) +
     &                 PROPPP(IDIM,-IFREQA,IFREQB) -
     &                 PROPPP(IDIM,IFREQA,-IFREQB) + 2.0D0*(
     &                 PROPPP(IDIM,-IFREQB,0) - PROPPP(IDIM,IFREQB,0)))/
     &                 (12.0D0*XFAMU*SQRT(XFAMU)*STEP**3)
                 DERIV3(IDIM,IFREQB,IFREQA)=(PROPPP(IDIM,IFREQB,IFREQA)-
     &                 PROPPP(IDIM,-IFREQB,-IFREQA) +
     &                 PROPPP(IDIM,-IFREQB,IFREQA) -
     &                 PROPPP(IDIM,IFREQB,-IFREQA) + 2.0D0*(
     &                 PROPPP(IDIM,-IFREQA,0) - PROPPP(IDIM,IFREQA,0)))/
     &                 (12.0D0*XFAMU*SQRT(XFAMU)*STEP**3)
               END IF
            END DO
            END DO
         END DO
      END IF
      RETURN
      END
      SUBROUTINE DERIV2(PROPPP,PROP2,NFREQ,NDIM,STEP)
C
C     Calculate derivatives of the molecular properties
C
      PARAMETER (XFAMU = 1822.88852 D0)
      DIMENSION PROPPP(NDIM,-NFREQ:NFREQ), PROP2(NDIM,NFREQ)
C
      CALL DZERO(PROP2,NDIM*NFREQ)
      DO IDIM = 1, NDIM
         DO IFREQ = 1, NFREQ
            PROP2(IDIM,IFREQ) = (PROPPP(IDIM,IFREQ) +
     &              PROPPP(IDIM,-IFREQ) - 2.0D0*PROPPP(IDIM,0))/
     &              ((STEP**2)*XFAMU)
         END DO
      END DO
      RETURN
      END
      SUBROUTINE MUF00(XMU1,F1,XMUF00,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), F1(3,3,3,NFREQ), XMUF00(3,3,3,3),
     &          OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQ = 1, NFREQ
            PREFAC = 1.0D0/(OMEGA(IFREQ)**2 - OPTFRQ**2)
            XMUF00(IA,IB,IC,ID) = XMUF00(IA,IB,IC,ID) +
     &           (XMU1(IA,IFREQ)*F1(IB,IC,ID,IFREQ) +
     &            XMU1(IB,IFREQ)*F1(IA,IC,ID,IFREQ))*PREFAC
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUF20(XMU1,F1,XMU2,F2,XMU3,F3,XMUF20,XMUF20H,OMEGA,
     &                 NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), F1(3,3,3,NFREQ),
     &          XMU2(3,NFREQ,NFREQ), F2(3,3,3,NFREQ,NFREQ),
     &          XMU3(3,NFREQ,NFREQ), F3(3,3,3,NFREQ,NFREQ),
     &          XMUF20(3,3,3,3), XMUF20H(3,3,3,3), OMEGA(NFREQ)
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            XMUF20H(IA,IB,IC,ID) = XMUF20H(IA,IB,IC,ID) +
     &      0.25D0*(XMU2(IA,IFREQA,IFREQB)*F2(IB,IC,ID,IFREQA,IFREQB)
     &            + XMU2(IB,IFREQA,IFREQB)*F2(IA,IC,ID,IFREQA,IFREQB))*
     &              (1.0D0/OMEGAA + 1.0D0/OMEGAB)/
     &              ((OMEGAA + OMEGAB)**2-OPTFRQ**2)
            XMUF20(IA,IB,IC,ID) = XMUF20(IA,IB,IC,ID) + 0.25D0*
     &           (XMU3(IA,IFREQA,IFREQB)*F1(IB,IC,ID,IFREQB) +
     &            XMU3(IB,IFREQA,IFREQB)*F1(IA,IC,ID,IFREQB) +
     &            XMU1(IA,IFREQB)*F3(IB,IC,ID,IFREQA,IFREQB) +
     &            XMU1(IB,IFREQB)*F3(IA,IC,ID,IFREQA,IFREQB))/
     &            (OMEGAA*(OMEGAB**2-OPTFRQ**2))
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUF11(XMU1,F1,XMU2,F2,CUBIC,XMUF11,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), F1(3,3,3,NFREQ),
     &          XMU2(3,NFREQ,NFREQ), F2(3,3,3,NFREQ,NFREQ),
     &          CUBIC(NFREQ,NFREQ,NFREQ), XMUF11(3,3,3,3), OMEGA(NFREQ)
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMUF11(IA,IB,IC,ID) = XMUF11(IA,IB,IC,ID) - 0.25D0*
     &          ((XMU2(IA,IFREQA,IFREQB)*F1(IB,IC,ID,IFREQC) +
     &            XMU2(IB,IFREQA,IFREQB)*F1(IA,IC,ID,IFREQC) +
     &            F2(IB,IC,ID,IFREQA,IFREQB)*XMU1(IA,IFREQC) +
     &            F2(IA,IC,ID,IFREQA,IFREQB)*XMU1(IB,IFREQC))*
     &            CUBIC(IFREQA,IFREQB,IFREQC)*
     &            (1.0D0/OMEGAA + 1.0D0/OMEGAB)/
     &       (((OMEGAA + OMEGAB)**2 - OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2))
     &         + (XMU2(IA,IFREQA,IFREQB)*F1(IB,IC,ID,IFREQA) +
     &            XMU2(IB,IFREQA,IFREQB)*F1(IA,IC,ID,IFREQA) +
     &            F2(IB,IC,ID,IFREQA,IFREQB)*XMU1(IA,IFREQA) +
     &            F2(IA,IC,ID,IFREQA,IFREQB)*XMU1(IB,IFREQA))*
     &            CUBIC(IFREQB,IFREQC,IFREQC)/((OMEGAB**2)*OMEGAC*
     &           (OMEGAA**2 - OPTFRQ**2)))
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUF02(XMU1,F1,CUBIC,QUARTIC,XMUF02,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), F1(3,3,3,NFREQ), 
     &          CUBIC(NFREQ,NFREQ,NFREQ), QUARTIC(NFREQ,NFREQ,NFREQ),
     &          XMUF02(3,3,3,3), OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMUF02(IA,IB,IC,ID) = XMUF02(IA,IB,IC,ID) - (0.25D0/OMEGAA)*
     &             QUARTIC(IFREQA,IFREQB,IFREQC)*
     &            (XMU1(IA,IFREQB)*F1(IB,IC,ID,IFREQC) + 
     &             XMU1(IB,IFREQB)*F1(IA,IC,ID,IFREQC))/
     &           ((OMEGAB**2-OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2))
         DO IFREQD = 1, NFREQ
            OMEGAD = OMEGA(IFREQD)
            XMUF02(IA,IB,IC,ID) = XMUF02(IA,IB,IC,ID) + (0.25D0/OMEGAA)*
     &         CUBIC(IFREQA,IFREQA,IFREQB)*CUBIC(IFREQB,IFREQC,IFREQD)*
     &        (XMU1(IA,IFREQC)*F1(IB,IC,ID,IFREQD) +
     &         XMU1(IB,IFREQC)*F1(IA,IC,ID,IFREQD))/((OMEGAB**2)*
     &         (OMEGAC**2-OPTFRQ**2)*(OMEGAD**2-OPTFRQ**2))
            XMUF02(IA,IB,IC,ID) = XMUF02(IA,IB,IC,ID) + (0.50D0/OMEGAA)*
     &         CUBIC(IFREQA,IFREQB,IFREQC)*CUBIC(IFREQA,IFREQB,IFREQD)*
     &        (XMU1(IA,IFREQC)*F1(IB,IC,ID,IFREQD) + 
     &         XMU1(IB,IFREQC)*F1(IA,IC,ID,IFREQD))/
     &      (((OMEGAA + OMEGAB)**2-OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2)*
     &        (OMEGAD**2-OPTFRQ**2))
         END DO
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUQ10(XMU1,Q1,XMU2,Q2,XMUQ10,OMEGA,NFREQ,OPTFRQ)
      DIMENSION  XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ),
     &           Q1(3,3,NFREQ), Q2(3,3,NFREQ,NFREQ),
     &           XMUQ10(3,3,3,3), OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            XMUQ10(IA,IB,IC,ID) = XMUQ10(IA,IB,IC,ID) + 0.50D0*
     &         (XMU1(IA,IFREQA)*XMU2(IB,IFREQA,IFREQB)*Q1(IC,ID,IFREQB))
     &           /((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2))
            XMUQ10(IA,IB,IC,ID) = XMUQ10(IA,IB,IC,ID) + 0.50D0*
     &         (XMU1(IA,IFREQA)*Q2(IC,ID,IFREQA,IFREQB)*XMU1(IB,IFREQB))
     &           /((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2-OPTFRQ**2))
            XMUQ10(IA,IB,IC,ID) = XMUQ10(IA,IB,IC,ID) + 0.50D0*
     &         (XMU1(IB,IFREQA)*XMU2(IA,IFREQA,IFREQB)*Q1(IC,ID,IFREQB))
     &           /((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2))
            XMUQ10(IA,IB,IC,ID) = XMUQ10(IA,IB,IC,ID) + 0.50D0*
     &         (XMU1(IB,IFREQA)*Q2(IC,ID,IFREQA,IFREQB)*XMU1(IA,IFREQB))
     &           /((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2-OPTFRQ**2))
            XMUQ10(IA,IB,IC,ID) = XMUQ10(IA,IB,IC,ID) + 0.50D0*
     &         (Q1(IC,ID,IFREQA)*XMU2(IA,IFREQA,IFREQB)*XMU1(IB,IFREQB))
     &           /((OMEGAA**2)*(OMEGAB**2-OPTFRQ**2))
            XMUQ10(IA,IB,IC,ID) = XMUQ10(IA,IB,IC,ID) + 0.50D0*
     &         (Q1(IC,ID,IFREQA)*XMU2(IB,IFREQA,IFREQB)*XMU1(IA,IFREQB))
     &           /((OMEGAA**2)*(OMEGAB**2-OPTFRQ**2))
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUQ01(XMU1,Q1,CUBIC,XMUQ01,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), Q1(3,3,NFREQ), CUBIC(NFREQ,NFREQ,NFREQ),
     &          XMUQ01(3,3,3,3), OMEGA(NFREQ)
      PARAMETER (D16 = 1.0D0/6.0D0)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            PREFAC = D16*CUBIC(IFREQA,IFREQB,IFREQC)
            XMUQ01(IA,IB,IC,ID) = XMUQ01(IA,IB,IC,ID) - PREFAC*
     &           XMU1(IA,IFREQA)*XMU1(IB,IFREQB)*Q1(IC,ID,IFREQC)/
     &        ((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2-OPTFRQ**2)*(OMEGAC**2))
            XMUQ01(IA,IB,IC,ID) = XMUQ01(IA,IB,IC,ID) - PREFAC*
     &           XMU1(IA,IFREQA)*Q1(IC,ID,IFREQB)*XMU1(IB,IFREQC)/
     &        ((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2)*(OMEGAC**2-OPTFRQ**2))
            XMUQ01(IA,IB,IC,ID) = XMUQ01(IA,IB,IC,ID) - PREFAC*
     &           XMU1(IB,IFREQA)*XMU1(IA,IFREQB)*Q1(IC,ID,IFREQC)/
     &        ((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2-OPTFRQ**2)*(OMEGAC**2))
            XMUQ01(IA,IB,IC,ID) = XMUQ01(IA,IB,IC,ID) - PREFAC*
     &           XMU1(IB,IFREQA)*Q1(IC,ID,IFREQB)*XMU1(IA,IFREQC)/
     &        ((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2)*(OMEGAC**2-OPTFRQ**2))
            XMUQ01(IA,IB,IC,ID) = XMUQ01(IA,IB,IC,ID) - PREFAC*
     &           Q1(IC,ID,IFREQA)*XMU1(IA,IFREQB)*XMU1(IB,IFREQC)/
     &        ((OMEGAA**2)*(OMEGAB**2-OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2))
            XMUQ01(IA,IB,IC,ID) = XMUQ01(IA,IB,IC,ID) - PREFAC*
     &           Q1(IC,ID,IFREQA)*XMU1(IB,IFREQB)*XMU1(IA,IFREQC)/
     &        ((OMEGAA**2)*(OMEGAB**2-OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2))
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE XAXI00(ALPHA1,XIP1,AXI00,OMEGA,NFREQ)
      DIMENSION ALPHA1(3,3,NFREQ), XIP1(3,3,NFREQ), OMEGA(NFREQ),
     &          AXI00(3,3,3,3)
C
C
      DO IA = 1, 3
         DO IB = 1, 3
            DO IC = 1, 3
               DO ID = 1, 3
                  DO IFREQ = 1, NFREQ
                     OMEGAI = 1.0D0/OMEGA(IFREQ)
                     AXI00(IA,IB,IC,ID) = AXI00(IA,IB,IC,ID) +
     &                    ALPHA1(IA,IB,IFREQ)*XIP1(IC,ID,IFREQ)/
     &                    (OMEGA(IFREQ)**2)
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
      SUBROUTINE XAXI20(ALPHA1,XIP1,ALPHA2,XIP2,ALPHA3,XIP3,AXI20,
     &                  AXI20H,OMEGA,NFREQ)
      DIMENSION ALPHA1(3,3,NFREQ), XIP1(3,3,NFREQ), 
     &          ALPHA2(3,3,NFREQ,NFREQ), XIP2(3,3,NFREQ,NFREQ),
     &          ALPHA3(3,3,NFREQ,NFREQ), XIP3(3,3,NFREQ,NFREQ),
     &          OMEGA(NFREQ), AXI20(3,3,3,3), AXI20H(3,3,3,3)
C
      DO IA = 1, 3
         DO IB = 1, 3
            DO IC = 1, 3
               DO ID = 1,3
                  DO IFREQA = 1, NFREQ
                     DO IFREQB = 1, NFREQ
                        OMEGAA = OMEGA(IFREQA)
                        OMEGAB = OMEGA(IFREQB)
                        PREFAC = 0.25D0/(OMEGAA*OMEGAB)
                        AXI20H(IA,IB,IC,ID) = AXI20H(IA,IB,IC,ID) +
     &    PREFAC*ALPHA2(IA,IB,IFREQA,IFREQB)*XIP2(IC,ID,IFREQA,IFREQB)/
     &           (OMEGAA + OMEGAB)
                        AXI20(IA,IB,IC,ID) = AXI20(IA,IB,IC,ID) +
     &    PREFAC*(ALPHA3(IA,IB,IFREQA,IFREQB)*XIP1(IC,ID,IFREQB)/OMEGAB
     &         +  ALPHA1(IA,IB,IFREQB)*XIP3(IC,ID,IFREQA,IFREQB)/OMEGAB)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
      SUBROUTINE XAXI11(ALPHA1,XIP1,ALPHA2,XIP2,CUBIC,OMEGA,NFREQ,AXI11)
      DIMENSION ALPHA1(3,3,NFREQ), XIP1(3,3,NFREQ),
     &          ALPHA2(3,3,NFREQ,NFREQ), XIP2(3,3,NFREQ,NFREQ),
     &          CUBIC(NFREQ,NFREQ,NFREQ), OMEGA(NFREQ),
     &          AXI11(3,3,3,3)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            PREFAC = -0.25D0/(OMEGAA*OMEGAB*OMEGAC)
            AXI11(IA,IB,IC,ID) = AXI11(IA,IB,IC,ID) 
     &                         + PREFAC*CUBIC(IFREQA,IFREQB,IFREQC)*
     &           (ALPHA2(IA,IB,IFREQA,IFREQB)*XIP1(IC,ID,IFREQC) +
     &            ALPHA1(IA,IB,IFREQC)*XIP2(IC,ID,IFREQA,IFREQB))/
     &           (OMEGAC*(OMEGAA + OMEGAB))
     &                         + PREFAC*CUBIC(IFREQB,IFREQC,IFREQC)*
     &           (ALPHA2(IA,IB,IFREQA,IFREQB)*XIP1(IC,ID,IFREQA) +
     &            ALPHA1(IA,IB,IFREQA)*XIP2(IC,ID,IFREQA,IFREQB))/
     &           (OMEGAA*OMEGAB)
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE XAXI02(ALPHA1,XIP1,CUBIC,QUARTIC,AXI02,OMEGA,NFREQ)
      DIMENSION ALPHA1(3,3,NFREQ), XIP1(3,3,NFREQ), AXI02(3,3,3,3),
     &          CUBIC(NFREQ,NFREQ,NFREQ), QUARTIC(NFREQ,NFREQ,NFREQ),
     &          OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            AXI02(IA,IB,IC,ID) = AXI02(IA,IB,IC,ID) - (0.125D0/OMEGAA)*
     &           QUARTIC(IFREQA,IFREQB,IFREQC)*(ALPHA1(IA,IB,IFREQB)*
     &           XIP1(IC,ID,IFREQC) + ALPHA1(IA,IB,IFREQC)*
     &           XIP1(IC,ID,IFREQB))/((OMEGAB**2)*(OMEGAC**2))
         DO IFREQD = 1, NFREQ
            OMEGAD = OMEGA(IFREQD)
            AXI02(IA,IB,IC,ID) = AXI02(IA,IB,IC,ID) + (0.125D0/OMEGAA)*
     &         (CUBIC(IFREQA,IFREQA,IFREQB)*CUBIC(IFREQB,IFREQC,IFREQD)*
     &          (ALPHA1(IA,IB,IFREQC)*XIP1(IC,ID,IFREQD)
     &         + ALPHA1(IA,IB,IFREQD)*XIP1(IC,ID,IFREQC))/
     &        ((OMEGAB**2)*(OMEGAC**2)*(OMEGAD**2)) + 2.0D0*
     &          CUBIC(IFREQA,IFREQB,IFREQC)*CUBIC(IFREQA,IFREQB,IFREQD)*
     &          (ALPHA1(IA,IB,IFREQC)*XIP1(IC,ID,IFREQD) 
     &         + ALPHA1(IA,IB,IFREQD)*XIP1(IC,ID,IFREQC))/
     &        (((OMEGAA + OMEGAB)**2)*(OMEGAC**2)*(OMEGAD**2)))
         END DO
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUXI10(XMU1,XIP1,XMU2,XIP2,XMUXI10,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ), XIP1(3,3,NFREQ),
     &          XIP2(3,3,NFREQ,NFREQ), XMUXI10(3,3,3,3), OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            XMUXI10(IA,IB,IC,ID) = XMUXI10(IA,IB,IC,ID) + 0.50D0*
     &       ((XMU1(IA,IFREQA)*XMU1(IB,IFREQB)*XIP2(IC,ID,IFREQA,IFREQB)
     &       + XMU1(IB,IFREQA)*XMU1(IA,IFREQB)
     &           *XIP2(IC,ID,IFREQA,IFREQB))/
     &       ((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2-OPTFRQ**2)) + 2.0D0*
     &        (XMU1(IA,IFREQA)*XMU2(IB,IFREQA,IFREQB)*XIP1(IC,ID,IFREQB)
     &       + XMU1(IB,IFREQA)*XMU2(IA,IFREQA,IFREQB)
     &       *XIP1(IC,ID,IFREQB))/((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2)))
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUXI01(XMU1,XIP1,CUBIC,XMUXI01,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), XIP1(3,3,NFREQ),
     &          CUBIC(NFREQ,NFREQ,NFREQ), XMUXI01(3,3,3,3), OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
      DO ID = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMUXI01(IA,IB,IC,ID) = XMUXI01(IA,IB,IC,ID) - 0.50D0*
     &           (XMU1(IA,IFREQA)*XMU1(IB,IFREQB)*XIP1(IC,ID,IFREQC) +
     &            XMU1(IB,IFREQA)*XMU1(IA,IFREQB)*XIP1(IC,ID,IFREQC))*
     &            CUBIC(IFREQA,IFREQB,IFREQC)/
     &         ((OMEGAA**2-OPTFRQ**2)*(OMEGAB**2-OPTFRQ**2)*(OMEGAC**2))
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE ETAZPRI(ETADZ,ETAPZ,OPTFRQ)
      DIMENSION ETADZ(3,3,3,3), ETAPZ(3,3,3,3), ETA(3,3,3,3)
      CHARACTER*1 XI, XJ, XK, XL
C
      WRITE (6,'(3(/2X,A))')
     & '******************************************************',
     & '*** ZPVA contributions to the hypermagnetizability ***',
     & '******************************************************'
      WRITE (6,'(//,2X,A,F12.6)') 'Frequency of the optical '//
     &     'process : ',OPTFRQ
C
C     DeltaEta
C
      DELTAD = 0.0D0
      DELTAP = 0.0D0
      DO I = 1, 3
         DO J = 1, 3
            DELTAD = DELTAD+0.20D0*(ETADZ(I,J,I,J)-ETADZ(I,I,J,J)/3.0D0)
            DELTAP = DELTAP+0.20D0*(ETAPZ(I,J,I,J)-ETAPZ(I,I,J,J)/3.0D0)
         END DO
      END DO
      DELTAT = DELTAD + DELTAP
C     
      CALL HEADER('Zero-point vibrational contribution to ETA',
     &             -1)
      WRITE (6,'(12X,3A16)') '  Diamagnetic   ','  Paramagnetic  ',
     &     '     Total      '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,3(F12.6,4X))') XI,XJ,XK,XL,
     &                 ETADZ(I,J,K,L), ETAPZ(I,J,K,L), 
     &                 ETADZ(I,J,K,L) + ETAPZ(I,J,K,L)
               END DO
            END DO
         END DO
      END DO
      WRITE (6,'(/,A,3(F12.6,X))') ' Pure vibrational correction to'//
     &     ' DeltaEta (dia/para/total): ',DELTAD,DELTAP,DELTAT
C
      RETURN
      END
      SUBROUTINE ETAPRI(XMUG,XMUQ,XMUF,AXI,XMUXI,IORDER,OPTFRQ)
C
C     The results we have obtained is not of much use unless we
C     print them
C
      CHARACTER*1 XI, XJ, XK, XL
      DIMENSION XMUG(3,3,3,3,5), XMUQ(3,3,3,3,2)
      DIMENSION XMUF(3,3,3,3,5), AXI(3,3,3,3,5), XMUXI(3,3,3,3,2),
     &          ETA(3,3,3,3,3)
C
      OPEN(UNIT=6,NAME='purevib.out',STATUS='NEW',FORM='FORMATTED')
c      WRITE (6,'(3(/2X,A)')
c     & '************************************************************',
c     & '* Pure vibrational contributions to hypermagnetizabilities *',
c     & '************************************************************'
      WRITE (6,'(//,2X,A,F12.6)') 'Frequency of the optical '//
     &     'process : ',OPTFRQ
C
C     Double Harmonic approximation
C
      CALL TITLER('The double-harmonic approximation','*',-1)
      CALL HEADER('Diamagnetic contribution [muG(0,0)]',-1)
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,F12.6)') XI,XJ,XK,XL,
     &                   -XMUG(I,J,K,L,1)
               END DO
            END DO
         END DO
      END DO
C
      CALL HEADER('Paramagnetic contribution',-1)
      WRITE (6,'(12X,3A16)') '    muF(0,0)    ','  AlphaXi(0,0)  ',
     &     '  Paramagnetic  '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,3(F12.6,4X))') XI,XJ,XK,XL,
     &                   XMUF(I,J,K,L,1),AXI(I,J,K,L,1),
     &                   XMUF(I,J,K,L,1) + AXI(I,J,K,L,1)
               END DO
            END DO
         END DO
      END DO
C
      IF (IORDER .GE. 1) THEN
C
C     Corrections from first level of second-order perturbation corrections
C     
      CALL TITLER('The II-prime approximation','*',-1)
      CALL HEADER('Diamagnetic contribution',-1)
      WRITE (6,'(8X,5A14)') '   muG(2,0)   ', '   muG(1,1)   ',
     &     '  mu^2Q(1,0)  ','  mu^2Q(0,1)  ',' Paramagnetic '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------------------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,2X,5(F12.6,2X))') XI,XJ,XK,XL,
     &                   - XMUG(I,J,K,L,3), - XMUG(I,J,K,L,4),
     &                   - XMUQ(I,J,K,L,1), - XMUQ(I,J,K,L,2),
     &                   - XMUG(I,J,K,L,3) - XMUG(I,J,K,L,4) +
     &                   - XMUQ(I,J,K,L,1) - XMUQ(I,J,K,L,2)
               END DO
            END DO
         END DO
      END DO
C
      CALL HEADER('Paramagnetic contribution',-1)
      WRITE (6,'(8X,5A14)') '   muF(2,0)   ', '   muF(1,1)   ',
     &     '  AlphaXi(2,0) ','  AlphaXi(1,1) ', '     Sum      '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------------------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,5(F12.6,2X))') XI,XJ,XK,XL,
     &                   XMUF(I,J,K,L,3), XMUF(I,J,K,L,4),
     &                   AXI(I,J,K,L,3), AXI(I,J,K,L,4),
     &                   XMUF(I,J,K,L,3) + XMUF(I,J,K,L,4) +
     &                   AXI(I,J,K,L,3) + AXI(I,J,K,L,4)
               END DO
            END DO
         END DO
      END DO
      WRITE (6,'(8X,4A14)') '  mu^2Xi(1,0) ', '  mu^2Xi(0,1) ',
     &     '     Sum      ',' Paramagnetic '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '--------------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,5(F12.6,2X))') XI,XJ,XK,XL,
     &                   XMUXI(I,J,K,L,1), XMUXI(I,J,K,L,2),
     &                   XMUF(I,J,K,L,3) + XMUF(I,J,K,L,4) +
     &                   AXI(I,J,K,L,3) + AXI(I,J,K,L,4),
     &                   XMUF(I,J,K,L,3) + XMUF(I,J,K,L,4) +
     &                   AXI(I,J,K,L,3) + AXI(I,J,K,L,4) +
     &                   XMUXI(I,J,K,L,1) + XMUXI(I,J,K,L,2)
               END DO
            END DO
         END DO
      END DO
      IF (IORDER .GE. 2) THEN
C
C     Corrections from second level of second-order perturbation corrections
C     
      CALL TITLER('The (II-II-prime) contributions','*',-1)
      CALL HEADER('Diamagnetic contribution',-1)
      WRITE (6,'(8X,3A14)') '   muG(2,0)c  ', '   muG(0,2)   ',
     &     ' Diamagnetic  '
      WRITE (6,'(A)') '------------------------------------------'//
     &                '--------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,2X,3(F12.6,2X))') XI,XJ,XK,XL,
     &                   - XMUG(I,J,K,L,2),
     &                   - XMUG(I,J,K,L,5), - XMUG(I,J,K,L,2) -
     &                   XMUG(I,J,K,L,3) - XMUG(I,J,K,L,5)
               END DO
            END DO
         END DO
      END DO
C
      CALL HEADER('Paramagnetic contribution',-1)
      WRITE (6,'(8X,5A14)') '   muF(2,0)C  ', '   muF(0,2)   ',
     &     ' AlphaXi(2,0)c','  AlphaXi(0,2) ', ' Paramagnetic '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------------------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,5(F12.6,2X))') XI,XJ,XK,XL,
     &                   XMUF(I,J,K,L,2),
     &                   XMUF(I,J,K,L,5),
     &                   AXI(I,J,K,L,2),
     &                   AXI(I,J,K,L,5),
     &                   XMUF(I,J,K,L,3) + XMUF(I,J,K,L,2) +
     &                   XMUF(I,J,K,L,5) + AXI(I,J,K,L,5) +
     &                   AXI(I,J,K,L,2) + AXI(I,J,K,L,3)
               END DO
            END DO
         END DO
      END DO
      END IF
      END IF
C
C     We now summarize the results for Eta
C     
C     Double harmonic approximation
C
      DO I = 1, 3
         DO J = 1, 3
            DO K = 1, 3
               DO L = 1, 3
                  ETA(I,J,K,L,1) = - XMUG(I,J,K,L,1)
                  ETA(I,J,K,L,2) = XMUF(I,J,K,L,1) + AXI(I,J,K,L,1)
                  ETA(I,J,K,L,3) = ETA(I,J,K,L,1) + ETA(I,J,K,L,2)
               END DO
            END DO
         END DO
      END DO
C
C     DeltaEta
C
      DELTAD = 0.0D0
      DELTAP = 0.0D0
      DELTAT = 0.0D0
      DO I = 1, 3
         DO J = 1, 3
            DELTAD = DELTAD+0.20D0*(ETA(I,J,I,J,1)-ETA(I,I,J,J,1)/3.0D0)
            DELTAP = DELTAP+0.20D0*(ETA(I,J,I,J,2)-ETA(I,I,J,J,2)/3.0D0)
            DELTAT = DELTAT+0.20D0*(ETA(I,J,I,J,3)-ETA(I,I,J,J,3)/3.0D0)
         END DO
      END DO
C     
      CALL HEADER('Vibrational contribution to ETA (0,0) approximation',
     &             -1)
      WRITE (6,'(12X,3A16)') '  Diamagnetic   ','  Paramagnetic  ',
     &     '     Total      '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,3(F12.6,4X))') XI,XJ,XK,XL,
     &                 ETA(I,J,K,L,1), ETA(I,J,K,L,2), ETA(I,J,K,L,3)
               END DO
            END DO
         END DO
      END DO
      WRITE (6,'(/,A,3(F12.6,X))') ' Pure vibrational correction to'//
     &     ' DeltaEta (dia/para/total): ',DELTAD,DELTAP,DELTAT
C
C     First level of second-order corrections
C
      IF (IORDER .GE. 1) THEN
      DO I = 1, 3
         DO J = 1, 3
            DO K = 1, 3
               DO L = 1, 3
                  ETA(I,J,K,L,1) = ETA(I,J,K,L,1) - XMUG(I,J,K,L,3)
     &                           - XMUG(I,J,K,L,4) - XMUQ(I,J,K,L,1)
     &                           - XMUQ(I,J,K,L,2)
                  ETA(I,J,K,L,2) = ETA(I,J,K,L,2) + XMUF(I,J,K,L,3)
     &                           + XMUF(I,J,K,L,4) + XMUXI(I,J,K,L,1)
     &                           + XMUXI(I,J,K,L,2) + AXI(I,J,K,L,3)
     &                           + AXI(I,J,K,L,4)
                  ETA(I,J,K,L,3) = ETA(I,J,K,L,1) + ETA(I,J,K,L,2)
               END DO
            END DO
         END DO
      END DO
C
C     DeltaEta
C
      DELTAD = 0.0D0
      DELTAP = 0.0D0
      DELTAT = 0.0D0
      DO I = 1, 3
         DO J = 1, 3
            DELTAD = DELTAD+0.20D0*(ETA(I,J,I,J,1)-ETA(I,I,J,J,1)/3.0D0)
            DELTAP = DELTAP+0.20D0*(ETA(I,J,I,J,2)-ETA(I,I,J,J,2)/3.0D0)
            DELTAT = DELTAT+0.20D0*(ETA(I,J,I,J,3)-ETA(I,I,J,J,3)/3.0D0)
         END DO
      END DO
      CALL HEADER('Vibrational contribution to ETA (IIp) approximation',
     &            -1)
      WRITE (6,'(12X,3A16)') '  Diamagnetic   ','  Paramagnetic  ',
     &     '     Total      '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,3(F12.6,4X))') XI,XJ,XK,XL,
     &                 ETA(I,J,K,L,1), ETA(I,J,K,L,2), ETA(I,J,K,L,3)
               END DO
            END DO
         END DO
      END DO
      WRITE (6,'(/,A,3(F12.6,X))') ' Pure vibrational correction to'//
     &     ' DeltaEta (dia/para/total): ',DELTAD,DELTAP,DELTAT
C
C     Second level of second-order corrections
C
      IF (IORDER .GE. 2) THEN
      DO I = 1, 3
         DO J = 1, 3
            DO K = 1, 3
               DO L = 1, 3
                  ETA(I,J,K,L,1) = ETA(I,J,K,L,1) - XMUG(I,J,K,L,2)
     &                           - XMUG(I,J,K,L,5) + XMUG(I,J,K,L,3)
                  ETA(I,J,K,L,2) = ETA(I,J,K,L,2) + XMUF(I,J,K,L,2)
     &                           + XMUF(I,J,K,L,5) + AXI(I,J,K,L,2)
     &                           + AXI(I,J,K,L,5) - XMUF(I,J,K,L,3)
     &                           - AXI(I,J,K,L,3)
                  ETA(I,J,K,L,3) = ETA(I,J,K,L,1) + ETA(I,J,K,L,2)
               END DO
            END DO
         END DO
      END DO
C
C     DeltaEta
C
      DELTAD = 0.0D0
      DELTAP = 0.0D0
      DELTAT = 0.0D0
      DO I = 1, 3
         DO J = 1, 3
            DELTAD = DELTAD+0.20D0*(ETA(I,J,I,J,1)-ETA(I,I,J,J,1)/3.0D0)
            DELTAP = DELTAP+0.20D0*(ETA(I,J,I,J,2)-ETA(I,I,J,J,2)/3.0D0)
            DELTAT = DELTAT+0.20D0*(ETA(I,J,I,J,3)-ETA(I,I,J,J,3)/3.0D0)
         END DO
      END DO
      CALL HEADER('Vibrational contribution to ETA (II) approximation',
     &            -1)
      WRITE (6,'(12X,3A16)') '  Diamagnetic   ','  Paramagnetic  ',
     &     '     Total      '
      WRITE (6,'(A)') '--------------------------------------------'//
     &     '----------------'
      DO I = 1, 3
         XI = CHAR(48+I)
         DO J = 1, I
            XJ = CHAR(48+J)
            DO K = 1, 3
               XK = CHAR(48+K)
               DO L = 1, K
                  XL = CHAR(48+L)
                  WRITE (6,'(2X,4A1,6X,3(F12.6,4X))') XI,XJ,XK,XL,
     &                 ETA(I,J,K,L,1), ETA(I,J,K,L,2), ETA(I,J,K,L,3)
               END DO
            END DO
         END DO
      END DO
      WRITE (6,'(/,A,3(F12.6,X))') ' Pure vibrational correction to'//
     &     ' DeltaEta (dia/para/total): ',DELTAD,DELTAP,DELTAT
      END IF
      END IF
C
      CLOSE(6)
      RETURN
      END
C
      SUBROUTINE POLPRI(XMU2,XMUA,XMU3,IORDER,OFRQS,OFRQ1,OFRQ2,OFRQ3,
     &                  ALFA,BETA,GAMMA)
C
C     The results we have obtained is not of much use unless we
C     print them
C
      LOGICAL ALFA, BETA, GAMMA
      CHARACTER*1 XI, XJ, XK, XL
      DIMENSION XMU2(3,3,5), XMUA(3,3,3,5), XMU3(3,3,3,2)
C
      OPEN(UNIT=12,NAME='purevib.out',STATUS='NEW',FORM='FORMATTED')
      WRITE (12,'(3(/2X,A))')
     & '*************************************************************',
     & '* Pure vibrational contributions to (hyper)polarizabilities *',
     & '*************************************************************'
      WRITE (12,'(//,2X,A,4(F8.4,1X))') 'Frequencies of the optical '//
     &     'process : ',OFRQS, OFRQ1, OFRQ2, OFRQ3
      call flush(12)
C
C     Double Harmonic approximation
C
      IF (ALFA) THEN
         IF (IORDER .EQ. 0) THEN
            CALL TITLER('The double-harmonic approximation','*',-1)
            DO I = 1, 3
               XI = CHAR(48+I)
               DO J = 1, I
                  XJ = CHAR(48+J)
                  WRITE (12,'(2X,2A1,6X,F12.6)') XI,XJ,XMU2(I,J,1)
               END DO
            END DO
         ELSE IF (IORDER .EQ. 1) THEN
            CALL TITLER("The II' approximation",'*',-1)
            DO I = 1, 3
               XI = CHAR(48+I)
               DO J = 1, I
                  XJ = CHAR(48+J)
                  WRITE (12,'(2X,2A1,2X,4(F12.6,2X))')XI,XJ,XMU2(I,J,1),
     &                 XMU2(I,J,3), XMU2(I,J,4), XMU2(I,J,1) + 
     &                 XMU2(I,J,3) + XMU2(I,J,4)
               END DO
            END DO
         ELSE
            CALL TITLER('The full second-order approximation','*',-1)
            DO I = 1, 3
               XI = CHAR(48+I)
               DO J = 1, I
                  XJ = CHAR(48+J)
                  WRITE (12,'(1X,2A1,1X,6(F11.6,1X))')XI,XJ,XMU2(I,J,1),
     &                 XMU2(I,J,3), XMU2(I,J,2) - XMU2(I,J,3),
     &                 XMU2(I,J,4), XMU2(I,J,5), XMU2(I,J,1) + 
     &                 XMU2(I,J,2) + XMU2(I,J,4) + XMU2(I,J,5)
               END DO
            END DO
         END IF
      END IF
C
      IF (BETA) THEN
         IF (IORDER .GE. 0) THEN
            CALL TITLER('The double-harmonic approximation','*',-1)
      call flush(12)
            DO I = 1, 3
               XI = CHAR(48+I)
               DO J = 1, 3
                  XJ = CHAR(48+J)
               DO K = 1, 3
                  XK = CHAR(48+K)
                  WRITE (12,'(2X,3A1,6X,F12.6)') XI,XJ,XK,XMUA(I,J,K,1)
      call flush(12)
               END DO
            END DO
            END DO
         END IF
      call flush(12)
         IF (IORDER .GE. 1) THEN
            CALL TITLER("The II' approximation",'*',-1)
            WRITE (12,'(7X,5A14)') '   muA(2,0)   ', '   muA(1,1)   ',
     &           '  mu^3(1,0)   ','  mu^3(0,1)   ',"     II'      "
       WRITE (12,'(A)') '--------------------------------------------'//
     &           '----------------------------------'
            DO I = 1, 3
               XI = CHAR(48+I)
               DO J = 1, 3
                  XJ = CHAR(48+J)
               DO K = 1, 3
                  XK = CHAR(48+K)
                  WRITE (12,'(2X,3A1,2X,5(F12.6,2X))') XI,XJ,XK,
     &                 XMUA(I,J,K,3), XMUA(I,J,K,4), XMU3(I,J,K,1),
     &                 XMU3(I,J,K,2), XMUA(I,J,K,3) + XMUA(I,J,K,4) + 
     &                 XMU3(I,J,K,1) + XMU3(I,J,K,2)
               END DO
            END DO
            END DO
         END IF
      call flush(12)
         IF (IORDER .GE. 2) THEN
            CALL TITLER('The full second-order approximation','*',-1)
            WRITE (12,'(7X,4A14)') '   muA(2,0)C  ', '   muA(0,2)   ',
     &           '     II      ','  TOTAL   '
       WRITE (12,'(A)') '--------------------------------------------'//
     &           '--------------------'
            DO I = 1, 3
               XI = CHAR(48+I)
               DO J = 1, 3
                  XJ = CHAR(48+J)
               DO K = 1, 3
                  XK = CHAR(48+K)
                  WRITE (12,'(1X,3A1,1X,4(F12.6,2X))') XI,XJ,XK,
     &                 XMUA(I,J,K,2)-XMUA(I,J,K,3), XMUA(I,J,K,5),
     &                 XMUA(I,J,K,2)-XMUA(I,J,K,3) + XMUA(I,J,K,5),
     &                 XMUA(I,J,K,1)+XMUA(I,J,K,2)+XMUA(I,J,K,4) +
     &                 XMUA(I,J,K,5)+XMU3(I,J,K,1)+XMU3(I,J,K,2)
               END DO
            END DO
         END DO
         END IF
      END IF
      call flush(12)
C
      CLOSE(12)
      RETURN
      END
      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
      PARAMETER (LUPRI = 12)
C
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
      RETURN
      END
      SUBROUTINE TITLER(HEAD,A,IN)
      CHARACTER HEAD*(*), A
      PARAMETER (LUPRI = 12)
C
      LENGTH = LEN(HEAD)
      IF (IN .EQ. 200) THEN
         LENGTH = LENGTH + 2
      ELSE IF (IN .GE. 100) THEN
         MARG = IN - 100
         IF (MARG .GT. 0) MARG = MARG + 1
         LENGTH = LENGTH + 2*MARG
      END IF
      IF (IN .GE. 0 .AND. IN .LT. 100) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(72 - LENGTH)/2 + 1)
      END IF
      IF (IN .EQ. 200) THEN
         WRITE (LUPRI, '(//,80A)')
     *      (' ',I=1,INDENT), '.', ('-',I=1,LENGTH), '.'
         WRITE (LUPRI,'(80A)') (' ',I=1,INDENT),'| ', HEAD, ' |'
         WRITE (LUPRI, '(80A)')
     *      (' ',I=1,INDENT), '`', ('-',I=1,LENGTH), ''''
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD, ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
      RETURN
      END
      SUBROUTINE DZERO(DX,LENGTH)
C
C Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
C
C   Subroutine DZERO sets a real array of length *LENGTH*
C   to zero.
C...................................................................
      DIMENSION DX(*)
C
      IF (LENGTH.LE.0) RETURN
C
      DO 100 I = 1,LENGTH
  100    DX(I) = 0.0D00
C
      RETURN
      END
      SUBROUTINE MU200(XMU1,XMU200,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), XMU200(3,3), OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
         DO IFREQ = 1, NFREQ
            PREFAC = 0.50D0/(OMEGA(IFREQ)**2 - OPTFRQ**2)
            XMU200(IA,IB) = XMU200(IA,IB) +
     &           (XMU1(IA,IFREQ)*XMU1(IB,IFREQ) +
     &            XMU1(IB,IFREQ)*XMU1(IA,IFREQ))*PREFAC
         END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MU220(XMU1,XMU2,XMU3,XMU220,XMU220H,OMEGA,
     &                 NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ),
     &          XMU3(3,NFREQ,NFREQ), XMU220(3,3), 
     &          XMU220H(3,3), OMEGA(NFREQ)
      DO IA = 1, 3
      DO IB = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            XMU220H(IA,IB) = XMU220H(IA,IB) +
     &      0.125D0*(XMU2(IA,IFREQA,IFREQB)*XMU2(IB,IFREQA,IFREQB)
     &             + XMU2(IB,IFREQA,IFREQB)*XMU2(IA,IFREQA,IFREQB))*
     &              (1.0D0/OMEGAA + 1.0D0/OMEGAB)/
     &              ((OMEGAA + OMEGAB)**2-OPTFRQ**2)
            XMU220(IA,IB) = XMU220H(IA,IB) + 0.25D0*
     &           (XMU3(IA,IFREQA,IFREQB)*XMU1(IB,IFREQB) +
     &            XMU3(IB,IFREQA,IFREQB)*XMU1(IA,IFREQB))/
     &            (OMEGAA*(OMEGAB**2-OPTFRQ**2))
         END DO
         END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MU211(XMU1,XMU2,CUBIC,XMU211,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ),
     &          CUBIC(NFREQ,NFREQ,NFREQ), XMU211(3,3), OMEGA(NFREQ)
      DO IA = 1, 3
      DO IB = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMU211(IA,IB) = XMU211(IA,IB) - 0.25D0*
     &          ((XMU2(IA,IFREQA,IFREQB)*XMU1(IB,IFREQC) +
     &            XMU2(IB,IFREQA,IFREQB)*XMU1(IA,IFREQC))*
     &            CUBIC(IFREQA,IFREQB,IFREQC)*
     &            (1.0D0/OMEGAA + 1.0D0/OMEGAB)/
     &       (((OMEGAA + OMEGAB)**2 - OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2))
     &         + (XMU2(IA,IFREQA,IFREQB)*XMU1(IB,IFREQA) +
     &            XMU2(IB,IFREQA,IFREQB)*XMU1(IA,IFREQA))*
     &            CUBIC(IFREQB,IFREQC,IFREQC)/((OMEGAB**2)*OMEGAC*
     &           (OMEGAA**2 - OPTFRQ**2)))
         END DO
         END DO
         END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MU202(XMU1,CUBIC,QUARTIC,XMU202,OMEGA,NFREQ,OPTFRQ)
      DIMENSION XMU1(3,NFREQ), XMU202(3,3), OMEGA(NFREQ),
     &          CUBIC(NFREQ,NFREQ,NFREQ), QUARTIC(NFREQ,NFREQ,NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMU202(IA,IB) = XMU202(IA,IB) - (0.125D0/OMEGAA)*
     &             QUARTIC(IFREQA,IFREQB,IFREQC)*
     &            (XMU1(IA,IFREQB)*XMU1(IB,IFREQC) + 
     &             XMU1(IB,IFREQB)*XMU1(IA,IFREQC))/
     &           ((OMEGAB**2-OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2))
         DO IFREQD = 1, NFREQ
            OMEGAD = OMEGA(IFREQD)
            XMU202(IA,IB) = XMU202(IA,IB) + (0.125D0/OMEGAA)*
     &         CUBIC(IFREQA,IFREQA,IFREQB)*CUBIC(IFREQB,IFREQC,IFREQD)*
     &        (XMU1(IA,IFREQC)*XMU1(IB,IFREQD) +
     &         XMU1(IB,IFREQC)*XMU1(IA,IFREQD))/((OMEGAB**2)*
     &         (OMEGAC**2-OPTFRQ**2)*(OMEGAD**2-OPTFRQ**2))
            XMU202(IA,IB) = XMU202(IA,IB) + (0.250D0/OMEGAA)*
     &         CUBIC(IFREQA,IFREQB,IFREQC)*CUBIC(IFREQA,IFREQB,IFREQD)*
     &        (XMU1(IA,IFREQC)*XMU1(IB,IFREQD) + 
     &         XMU1(IB,IFREQC)*XMU1(IA,IFREQD))/
     &      (((OMEGAA + OMEGAB)**2-OPTFRQ**2)*(OMEGAC**2-OPTFRQ**2)*
     &        (OMEGAD**2-OPTFRQ**2))
         END DO
         END DO
         END DO
         END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUA00(XMU1,ALPHA1,XMUA00,OMEGA,NFREQ,OFRQS,OFRQ1,
     &                   OFRQ2)
      DIMENSION XMU1(3,NFREQ), ALPHA1(3,3,NFREQ), XMUA00(3,3,3), 
     &          OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
         DO IFREQ = 1, NFREQ
            OMEG = OMEGA(IFREQ)
            XMUA00(IA,IB,IC) = XMUA00(IA,IB,IC) + 0.50D0*(
     &    XMU1(IA,IFREQ)*ALPHA1(IB,IC,IFREQ)/(OMEG**2-OFRQS**2)
     &  + XMU1(IA,IFREQ)*ALPHA1(IC,IB,IFREQ)/(OMEG**2-OFRQS**2)
     &  + XMU1(IB,IFREQ)*ALPHA1(IA,IC,IFREQ)/(OMEG**2-OFRQ1**2)
     &  + XMU1(IB,IFREQ)*ALPHA1(IC,IA,IFREQ)/(OMEG**2-OFRQ1**2)
     &  + XMU1(IC,IFREQ)*ALPHA1(IA,IB,IFREQ)/(OMEG**2-OFRQ2**2)
     &  + XMU1(IC,IFREQ)*ALPHA1(IB,IA,IFREQ)/(OMEG**2-OFRQ2**2))
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUA20(XMU1,XMU2,XMU3,ALPHA1,ALPHA2,ALPHA3,XMUA20,
     &                 XMUA20H,OMEGA,NFREQ,OFRQS,OFRQ1,OFRQ2)
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ),
     &          XMU3(3,NFREQ,NFREQ), ALPHA1(3,3,NFREQ),
     &          ALPHA2(3,3,NFREQ,NFREQ), ALPHA3(3,3,NFREQ,NFREQ), 
     &          XMUA20(3,3,3), XMUA20H(3,3,3), OMEGA(NFREQ)
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            XMUA20H(IA,IB,IC) = XMUA20H(IA,IB,IC) +
     &            0.125D0*(1.0D0/OMEGAA + 1.0D0/OMEGAB)*
     &           (XMU2(IA,IFREQA,IFREQB)*ALPHA2(IB,IC,IFREQA,IFREQB)/
     &            ((OMEGAA + OMEGAB)**2-OFRQS**2) + 
     &            XMU2(IA,IFREQA,IFREQB)*ALPHA2(IC,IB,IFREQA,IFREQB)/
     &            ((OMEGAA + OMEGAB)**2-OFRQS**2) + 
     &            XMU2(IB,IFREQA,IFREQB)*ALPHA2(IA,IC,IFREQA,IFREQB)/
     &            ((OMEGAA + OMEGAB)**2-OFRQ1**2) + 
     &            XMU2(IB,IFREQA,IFREQB)*ALPHA2(IC,IA,IFREQA,IFREQB)/
     &            ((OMEGAA + OMEGAB)**2-OFRQ1**2) + 
     &            XMU2(IC,IFREQA,IFREQB)*ALPHA2(IA,IB,IFREQA,IFREQB)/
     &            ((OMEGAA + OMEGAB)**2-OFRQ2**2) + 
     &            XMU2(IC,IFREQA,IFREQB)*ALPHA2(IB,IA,IFREQA,IFREQB)/
     &            ((OMEGAA + OMEGAB)**2-OFRQ2**2))
            XMUA20(IA,IB,IC) =XMUA20H(IA,IB,IC)+(0.125D0/OMEGAA)*
     &          ((XMU3(IA,IFREQA,IFREQB)*ALPHA1(IB,IC,IFREQB) +
     &            ALPHA3(IB,IC,IFREQA,IFREQB)*XMU1(IA,IFREQB))/
     &            (OMEGAB**2-OFRQS**2) +
     &           (XMU3(IA,IFREQA,IFREQB)*ALPHA1(IC,IB,IFREQB) +
     &            ALPHA3(IC,IB,IFREQA,IFREQB)*XMU1(IA,IFREQB))/
     &            (OMEGAB**2-OFRQS**2) +
     &           (XMU3(IB,IFREQA,IFREQB)*ALPHA1(IA,IC,IFREQB) +
     &            ALPHA3(IA,IC,IFREQA,IFREQB)*XMU1(IB,IFREQB))/
     &            (OMEGAB**2-OFRQ1**2) +
     &           (XMU3(IB,IFREQA,IFREQB)*ALPHA1(IC,IA,IFREQB) +
     &            ALPHA3(IC,IA,IFREQA,IFREQB)*XMU1(IB,IFREQB))/
     &            (OMEGAB**2-OFRQ1**2) +
     &           (XMU3(IC,IFREQA,IFREQB)*ALPHA1(IA,IB,IFREQB) +
     &            ALPHA3(IA,IB,IFREQA,IFREQB)*XMU1(IC,IFREQB))/
     &            (OMEGAB**2-OFRQ2**2) +
     &           (XMU3(IC,IFREQA,IFREQB)*ALPHA1(IB,IA,IFREQB) +
     &            ALPHA3(IB,IA,IFREQA,IFREQB)*XMU1(IC,IFREQB))/
     &            (OMEGAB**2-OFRQ2**2))
         END DO
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUA11(XMU1,XMU2,ALPHA1,ALPHA2,CUBIC,XMUA11,OMEGA,NFREQ,
     &                 OFRQS,OFRQ1,OFRQ2)
      DIMENSION XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ),
     &          ALPHA1(3,3,NFREQ), ALPHA2(3,3,NFREQ,NFREQ),
     &          CUBIC(NFREQ,NFREQ,NFREQ), XMUA11(3,3,3), OMEGA(NFREQ)
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMUA11(IA,IB,IC) = XMUA11(IA,IB,IC) - 0.125D0*
     &         CUBIC(IFREQA,IFREQB,IFREQC)*(1.0D0/OMEGAA + 1.0D0/OMEGAB)
     &         *((XMU2(IA,IFREQA,IFREQB)*ALPHA1(IB,IC,IFREQC) +
     &            ALPHA2(IB,IC,IFREQA,IFREQB)*XMU1(IA,IFREQC))/
     &          (((OMEGAA + OMEGAB)**2 - OFRQS**2)*(OMEGAC**2-OFRQS**2))
     &         + (XMU2(IA,IFREQA,IFREQB)*ALPHA1(IC,IB,IFREQC) +
     &            ALPHA2(IC,IB,IFREQA,IFREQB)*XMU1(IA,IFREQC))/
     &          (((OMEGAA + OMEGAB)**2 - OFRQS**2)*(OMEGAC**2-OFRQS**2))
     &         + (XMU2(IB,IFREQA,IFREQB)*ALPHA1(IA,IC,IFREQC) +
     &            ALPHA2(IA,IC,IFREQA,IFREQB)*XMU1(IB,IFREQC))/
     &          (((OMEGAA + OMEGAB)**2 - OFRQ1**2)*(OMEGAC**2-OFRQ1**2))
     &         + (XMU2(IB,IFREQA,IFREQB)*ALPHA1(IC,IA,IFREQC) +
     &            ALPHA2(IC,IA,IFREQA,IFREQB)*XMU1(IB,IFREQC))/
     &          (((OMEGAA + OMEGAB)**2 - OFRQ1**2)*(OMEGAC**2-OFRQ1**2))
     &         + (XMU2(IC,IFREQA,IFREQB)*ALPHA1(IA,IB,IFREQC) +
     &            ALPHA2(IA,IB,IFREQA,IFREQB)*XMU1(IC,IFREQC))/
     &          (((OMEGAA + OMEGAB)**2 - OFRQ2**2)*(OMEGAC**2-OFRQ2**2))
     &         + (XMU2(IC,IFREQA,IFREQB)*ALPHA1(IB,IA,IFREQC) +
     &            ALPHA2(IB,IA,IFREQA,IFREQB)*XMU1(IC,IFREQC))/
     &         (((OMEGAA + OMEGAB)**2 - OFRQ2**2)*(OMEGAC**2-OFRQ2**2)))
            XMUA11(IA,IB,IC) = XMUA11(IA,IB,IC) - 0.125D0/
     &            ((OMEGAB**2)*OMEGAC)*CUBIC(IFREQB,IFREQC,IFREQC)*
     &          ((XMU2(IA,IFREQA,IFREQB)*ALPHA1(IB,IC,IFREQA) +
     &            ALPHA2(IB,IC,IFREQA,IFREQB)*XMU1(IA,IFREQA))/
     &           (OMEGAA**2 - OFRQS**2) +
     &           (XMU2(IA,IFREQA,IFREQB)*ALPHA1(IC,IB,IFREQA) +
     &            ALPHA2(IC,IB,IFREQA,IFREQB)*XMU1(IA,IFREQA))/
     &           (OMEGAA**2 - OFRQS**2) +
     &           (XMU2(IB,IFREQA,IFREQB)*ALPHA1(IA,IC,IFREQA) +
     &            ALPHA2(IA,IC,IFREQA,IFREQB)*XMU1(IB,IFREQA))/
     &           (OMEGAA**2 - OFRQ1**2) +
     &           (XMU2(IB,IFREQA,IFREQB)*ALPHA1(IC,IA,IFREQA) +
     &            ALPHA2(IC,IA,IFREQA,IFREQB)*XMU1(IB,IFREQA))/
     &           (OMEGAA**2 - OFRQ1**2) +
     &           (XMU2(IC,IFREQA,IFREQB)*ALPHA1(IA,IB,IFREQA) +
     &            ALPHA2(IA,IB,IFREQA,IFREQB)*XMU1(IC,IFREQA))/
     &           (OMEGAA**2 - OFRQ2**2) +
     &           (XMU2(IC,IFREQA,IFREQB)*ALPHA1(IB,IA,IFREQA) +
     &            ALPHA2(IB,IA,IFREQA,IFREQB)*XMU1(IC,IFREQA))/
     &           (OMEGAA**2 - OFRQ2**2))
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MUA02(XMU1,ALPHA1,CUBIC,QUARTIC,XMUA02,OMEGA,NFREQ,
     8                 OFRQS,OFRQ1,OFRQ2)
      DIMENSION XMU1(3,NFREQ), ALPHA1(3,3,NFREQ), XMUA02(3,3,3), 
     &          OMEGA(NFREQ), CUBIC(NFREQ,NFREQ,NFREQ), 
     &          QUARTIC(NFREQ,NFREQ,NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            XMUA02(IA,IB,IC) = XMUA02(IA,IB,IC) - (0.125D0/OMEGAA)*
     &             QUARTIC(IFREQA,IFREQB,IFREQC)*
     &             (XMU1(IA,IFREQB)*ALPHA1(IB,IC,IFREQC)/      
     &            ((OMEGAB**2-OFRQS**2)*(OMEGAC**2-OFRQS**2))
     &            + XMU1(IA,IFREQB)*ALPHA1(IC,IB,IFREQC)/      
     &            ((OMEGAB**2-OFRQS**2)*(OMEGAC**2-OFRQS**2))
     &            + XMU1(IB,IFREQB)*ALPHA1(IA,IC,IFREQC)/      
     &            ((OMEGAB**2-OFRQ1**2)*(OMEGAC**2-OFRQ1**2))
     &            + XMU1(IB,IFREQB)*ALPHA1(IC,IA,IFREQC)/      
     &            ((OMEGAB**2-OFRQ1**2)*(OMEGAC**2-OFRQ1**2))
     &            + XMU1(IC,IFREQB)*ALPHA1(IA,IB,IFREQC)/      
     &            ((OMEGAB**2-OFRQ2**2)*(OMEGAC**2-OFRQ2**2))
     &            + XMU1(IC,IFREQB)*ALPHA1(IB,IA,IFREQC)/      
     &            ((OMEGAB**2-OFRQ2**2)*(OMEGAC**2-OFRQ2**2)))
         DO IFREQD = 1, NFREQ
            OMEGAD = OMEGA(IFREQD)
            XMUA02(IA,IB,IC) = XMUA02(IA,IB,IC) +
     &           (0.125D0/(OMEGAA*(OMEGAB**2)))*
     &         CUBIC(IFREQA,IFREQA,IFREQB)*CUBIC(IFREQB,IFREQC,IFREQD)*
     &         (XMU1(IA,IFREQC)*ALPHA1(IB,IC,IFREQD)/
     &         ((OMEGAC**2-OFRQS**2)*(OMEGAD**2-OFRQS**2)) +
     &          XMU1(IA,IFREQC)*ALPHA1(IC,IB,IFREQD)/
     &         ((OMEGAC**2-OFRQS**2)*(OMEGAD**2-OFRQS**2)) +
     &          XMU1(IB,IFREQC)*ALPHA1(IA,IC,IFREQD)/
     &         ((OMEGAC**2-OFRQ1**2)*(OMEGAD**2-OFRQ1**2)) +
     &          XMU1(IB,IFREQC)*ALPHA1(IC,IA,IFREQD)/
     &         ((OMEGAC**2-OFRQ1**2)*(OMEGAD**2-OFRQ1**2)) +
     &          XMU1(IC,IFREQC)*ALPHA1(IA,IB,IFREQD)/
     &         ((OMEGAC**2-OFRQ2**2)*(OMEGAD**2-OFRQ2**2)) +
     &          XMU1(IC,IFREQC)*ALPHA1(IB,IA,IFREQD)/
     &         ((OMEGAC**2-OFRQ2**2)*(OMEGAD**2-OFRQ2**2)))
            XMUA02(IA,IB,IC) = XMUA02(IA,IB,IC) + (0.250D0/OMEGAA)*
     &         CUBIC(IFREQA,IFREQB,IFREQC)*CUBIC(IFREQA,IFREQB,IFREQD)*
     &        (XMU1(IA,IFREQC)*ALPHA1(IB,IC,IFREQD)/
     &      (((OMEGAA + OMEGAB)**2-OFRQS**2)*(OMEGAC**2-OFRQS**2)*
     &        (OMEGAD**2-OFRQS**2)) +
     &         XMU1(IA,IFREQC)*ALPHA1(IC,IB,IFREQD)/
     &      (((OMEGAA + OMEGAB)**2-OFRQS**2)*(OMEGAC**2-OFRQS**2)*
     &        (OMEGAD**2-OFRQS**2)) +
     &         XMU1(IB,IFREQC)*ALPHA1(IA,IC,IFREQD)/
     &      (((OMEGAA + OMEGAB)**2-OFRQ1**2)*(OMEGAC**2-OFRQ1**2)*
     &        (OMEGAD**2-OFRQ1**2)) +
     &         XMU1(IB,IFREQC)*ALPHA1(IC,IA,IFREQD)/
     &      (((OMEGAA + OMEGAB)**2-OFRQ1**2)*(OMEGAC**2-OFRQ1**2)*
     &        (OMEGAD**2-OFRQ1**2)) +
     &         XMU1(IC,IFREQC)*ALPHA1(IA,IB,IFREQD)/
     &      (((OMEGAA + OMEGAB)**2-OFRQ2**2)*(OMEGAC**2-OFRQ2**2)*
     &        (OMEGAD**2-OFRQ2**2)) +
     &         XMU1(IC,IFREQC)*ALPHA1(IB,IA,IFREQD)/
     &      (((OMEGAA + OMEGAB)**2-OFRQ2**2)*(OMEGAC**2-OFRQ2**2)*
     &        (OMEGAD**2-OFRQ2**2)))
         END DO
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MU310(XMU1,XMU2,XMU310,OMEGA,NFREQ,OFRQS,OFRQ1,OFRQ2)
      DIMENSION  XMU1(3,NFREQ), XMU2(3,NFREQ,NFREQ),
     &           XMU310(3,3,3), OMEGA(NFREQ)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            XMU310(IA,IB,IC) = XMU310(IA,IB,IC) + 0.50D0*
     &         (XMU1(IA,IFREQA)*XMU2(IB,IFREQA,IFREQB)*XMU1(IC,IFREQB))
     &           /((OMEGAA**2-OFRQS**2)*(OMEGAB**2-OFRQ2**2))
            XMU310(IA,IB,IC) = XMU310(IA,IB,IC) + 0.50D0*
     &         (XMU1(IA,IFREQA)*XMU2(IC,IFREQA,IFREQB)*XMU1(IB,IFREQB))
     &           /((OMEGAA**2-OFRQS**2)*(OMEGAB**2-OFRQ1**2))
            XMU310(IA,IB,IC) = XMU310(IA,IB,IC) + 0.50D0*
     &         (XMU1(IB,IFREQA)*XMU2(IA,IFREQA,IFREQB)*XMU1(IC,IFREQB))
     &           /((OMEGAA**2-OFRQ1**2)*(OMEGAB**2-OFRQ2**2))
            XMU310(IA,IB,IC) = XMU310(IA,IB,IC) + 0.50D0*
     &         (XMU1(IB,IFREQA)*XMU2(IC,IFREQA,IFREQB)*XMU1(IA,IFREQB))
     &           /((OMEGAA**2-OFRQ1**2)*(OMEGAB**2-OFRQS**2))
            XMU310(IA,IB,IC) = XMU310(IA,IB,IC) + 0.50D0*
     &         (XMU1(IC,IFREQA)*XMU2(IA,IFREQA,IFREQB)*XMU1(IB,IFREQB))
     &           /((OMEGAA**2-OFRQ2**2)*(OMEGAB**2-OFRQ1**2))
            XMU310(IA,IB,IC) = XMU310(IA,IB,IC) + 0.50D0*
     &         (XMU1(IC,IFREQA)*XMU2(IB,IFREQA,IFREQB)*XMU1(IA,IFREQB))
     &           /((OMEGAA**2-OFRQ2**2)*(OMEGAB**2-OFRQS**2))
         END DO
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE MU301(XMU1,CUBIC,XMU301,OMEGA,NFREQ,OFRQS,OFRQ1,OFRQ2)
      DIMENSION XMU1(3,NFREQ), CUBIC(NFREQ,NFREQ,NFREQ),
     &          XMU301(3,3,3), OMEGA(NFREQ)
      PARAMETER (D16 = 1.0D0/6.0D0)
C
      DO IA = 1, 3
      DO IB = 1, 3
      DO IC = 1, 3
         DO IFREQA = 1, NFREQ
         DO IFREQB = 1, NFREQ
         DO IFREQC = 1, NFREQ
            OMEGAA = OMEGA(IFREQA)
            OMEGAB = OMEGA(IFREQB)
            OMEGAC = OMEGA(IFREQC)
            PREFAC = D16*CUBIC(IFREQA,IFREQB,IFREQC)
            XMU301(IA,IB,IC) = XMU301(IA,IB,IC) - PREFAC*
     &           XMU1(IA,IFREQA)*XMU1(IB,IFREQB)*XMU1(IC,IFREQC)/
     &  ((OMEGAA**2-OFRQS**2)*(OMEGAB**2-OFRQ1**2)*(OMEGAC**2-OFRQ2**2))
            XMU301(IA,IB,IC) = XMU301(IA,IB,IC) - PREFAC*
     &           XMU1(IA,IFREQA)*XMU1(IC,IFREQB)*XMU1(IB,IFREQC)/
     &  ((OMEGAA**2-OFRQS**2)*(OMEGAB**2-OFRQ2**2)*(OMEGAC**2-OFRQ1**2))
            XMU301(IA,IB,IC) = XMU301(IA,IB,IC) - PREFAC*
     &           XMU1(IB,IFREQA)*XMU1(IA,IFREQB)*XMU1(IC,IFREQC)/
     &  ((OMEGAA**2-OFRQ1**2)*(OMEGAB**2-OFRQS**2)*(OMEGAC**2-OFRQ2**2))
            XMU301(IA,IB,IC) = XMU301(IA,IB,IC) - PREFAC*
     &           XMU1(IB,IFREQA)*XMU1(IC,IFREQB)*XMU1(IA,IFREQC)/
     &  ((OMEGAA**2-OFRQ1**2)*(OMEGAB**2-OFRQ2**2)*(OMEGAC**2-OFRQS**2))
            XMU301(IA,IB,IC) = XMU301(IA,IB,IC) - PREFAC*
     &           XMU1(IC,IFREQA)*XMU1(IA,IFREQB)*XMU1(IB,IFREQC)/
     &  ((OMEGAA**2-OFRQ2**2)*(OMEGAB**2-OFRQS**2)*(OMEGAC**2-OFRQ1**2))
            XMU301(IA,IB,IC) = XMU301(IA,IB,IC) - PREFAC*
     &           XMU1(IC,IFREQA)*XMU1(IB,IFREQB)*XMU1(IA,IFREQC)/
     &  ((OMEGAA**2-OFRQ2**2)*(OMEGAB**2-OFRQ1**2)*(OMEGAC**2-OFRQS**2))
         END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
