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
C  /* Deck nucinp */
      SUBROUTINE NUCINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 3)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <cbinuc.h>
#include <abainf.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.STOP  '/
C
      NEWDEF = (WORD .EQ. '*NUCREP')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in NUCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in NUCINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  CUT  = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in NUCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in NUCINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for NUCREP:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' NUCREP skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in NUCREP:',IPRINT
            END IF
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after NUCREP.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck nucini */
      SUBROUTINE NUCINI
C
C     Initialize /CBINUC/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbinuc.h>
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      IF (MOLHES) THEN
         MAXDIF = 2
      ELSE IF (MOLGRD) THEN
         MAXDIF = 1
      ELSE
         SKIP = .TRUE.
      END IF
      RETURN
      END
C  /* Deck nucrep */
      SUBROUTINE NUCREP(HESSNN,CSTRA,SCTRA)
C
C     Jan 1985 tuh
C     Revised Jan 9 1987 tuh - zero charges introduced
C     Revised Jun 25 1988 tuh - symmetry introduced
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0,
     *           THREE = 3.0D0)
      LOGICAL SECOND
      DIMENSION HESSNN(MXCOOR,MXCOOR), CSTRA(*), SCTRA(*)
#include <cbinuc.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
#include <energy.h>
#include <past.h>
#include <ibtfun.h>
C
C     This subroutine calculates the nuclear repulsion contributions
C     to the molecular gradient and Hessian.
C
      IF (SKIP) RETURN
      CALL QENTER('NUCREP')
      IF (PASONE .AND. IPRINT .LT. 5) IPRINT = 0
      IF (IPRINT .GT. 1) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         WRITE (LUPRI,'(A)')
     *      '  <<<<<<<<<< Output from NUCREP >>>>>>>>>> '
      END IF
C
C     **********************
C     ***** INITIALIZE *****
C     **********************
C
      SECOND = MAXDIF .EQ. 2
      NCDEP3 = 3*NUCDEP
      CALL DZERO(GRADNN(1),NCDEP3)
      IF (SECOND) CALL HESZER(HESSNN)
C
C     ******************************************
C     ***** CALCULATE GRADIENT AND HESSIAN *****
C     ******************************************
C
C     Run over symmetry-independent nuclei A
C
      DO 200 NCENTA = 1, NUCIND
         CHARGA = CHARGE(NCENTA)
         IF (ABS(CHARGA) .GT. ZERO) THEN
            NAX    = 3*NCENTA - 2
            NAY    = 3*NCENTA - 1
            NAZ    = 3*NCENTA
            CORDXA = CORD(1,NCENTA)
            CORDYA = CORD(2,NCENTA)
            CORDZA = CORD(3,NCENTA)
            MULA   = ISTBNU(NCENTA)
C
C           Run over symmetry-independent nuclei B
C
            DO 300 NCENTB =  1, NCENTA
               CHARGB = CHARGE(NCENTB)
               IF (ABS(CHARGB) .GT. ZERO) THEN
                  NBX    = 3*NCENTB - 2
                  NBY    = 3*NCENTB - 1
                  NBZ    = 3*NCENTB
                  CORBX0 = CORD(1,NCENTB)
                  CORBY0 = CORD(2,NCENTB)
                  CORBZ0 = CORD(3,NCENTB)
                  MULB   = ISTBNU(NCENTB)
C
                  MAB    = IBTOR (MULA,MULB)
                  KAB    = IBTAND(MULA,MULB)
                  HKAB   = FMULT(KAB)
                  CROSS = ONE
                  IF (NCENTA .EQ. NCENTB) THEN
                     HKAB = HALF*HKAB
                     CROSS = TWO
                  END IF
C
C                 Run over symmetry independent charge distributions
C
                  DO 400 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MAB) .EQ. 0) THEN
                     ICENTA = NUCNUM(NCENTA,1)
                     ICENTB = NUCNUM(NCENTB,ISYMOP+1)
                  IF (ICENTA .EQ. ICENTB) GO TO 400
C
                     SBX = PT(IBTAND(ISYMAX(1,1),ISYMOP))
                     SBY = PT(IBTAND(ISYMAX(2,1),ISYMOP))
                     SBZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))
                     XAB = CORDXA - SBX*CORBX0
                     YAB = CORDYA - SBY*CORBY0
                     ZAB = CORDZA - SBZ*CORBZ0
C
                     XAB2   = XAB*XAB
                     YAB2   = YAB*YAB
                     ZAB2   = ZAB*ZAB
                     RAB2   = XAB2 + YAB2 + ZAB2
                     RAB1   = SQRT(RAB2)
                     ZZR3IN = - HKAB*CHARGA*CHARGB/(RAB1*RAB2)
                     VNUCX  = XAB*ZZR3IN
                     VNUCY  = YAB*ZZR3IN
                     VNUCZ  = ZAB*ZZR3IN
C
C                    ********************
C                    ***** Gradient *****
C                    ********************
C
C                    Totally symmetric contributions only
C
                     IF (DOREPS(0)) THEN
                        IAX  = IPTCNT(NAX,0,1)
                        IAY  = IPTCNT(NAY,0,1)
                        IAZ  = IPTCNT(NAZ,0,1)
                        IBX  = IPTCNT(NBX,0,1)
                        IBY  = IPTCNT(NBY,0,1)
                        IBZ  = IPTCNT(NBZ,0,1)
                        IF (IAX.NE.0) GRADNN(IAX) =GRADNN(IAX)+VNUCX
                        IF (IAY.NE.0) GRADNN(IAY) =GRADNN(IAY)+VNUCY
                        IF (IAZ.NE.0) GRADNN(IAZ) =GRADNN(IAZ)+VNUCZ
                        IF (IBX.NE.0) GRADNN(IBX) =GRADNN(IBX)-SBX*VNUCX
                        IF (IBY.NE.0) GRADNN(IBY) =GRADNN(IBY)-SBY*VNUCY
                        IF (IBZ.NE.0) GRADNN(IBZ) =GRADNN(IBZ)-SBZ*VNUCZ
                     END IF
C
C                    *******************
C                    ***** Hessian *****
C                    *******************
C
                     IF (SECOND) THEN
                        ZZR5IN = ZZR3IN/RAB2
                        VNUCXX = (XAB2 + XAB2 - YAB2 - ZAB2)*ZZR5IN
                        VNUCXY = THREE*XAB*YAB*ZZR5IN
                        VNUCXZ = THREE*XAB*ZAB*ZZR5IN
                        VNUCYY = (YAB2 + YAB2 - ZAB2 - XAB2)*ZZR5IN
                        VNUCYZ = THREE*YAB*ZAB*ZZR5IN
                        VNUCZZ = (ZAB2 + ZAB2 - XAB2 - YAB2)*ZZR5IN
C
C                       Loop over irreps for differentiation operator
C
                        DO 500 IREPD = 0, MAXREP
                        IF (DOREPS(IREPD)) THEN
                           CRSCHI = CROSS*PT(IBTAND(ISYMOP,IREPD))
                           IAX  = IPTCNT(NAX,IREPD,1)
                           IAY  = IPTCNT(NAY,IREPD,1)
                           IAZ  = IPTCNT(NAZ,IREPD,1)
                           IBX  = IPTCNT(NBX,IREPD,1)
                           IBY  = IPTCNT(NBY,IREPD,1)
                           IBZ  = IPTCNT(NBZ,IREPD,1)
                           IF (IAX*IBX.NE.0) THEN
                              HESSNN(IAX,IBX) = HESSNN(IAX,IBX)
     *                                        + CRSCHI*SBX*VNUCXX
                           END IF
                           IF (IAX*IBY.NE.0) THEN
                              HESSNN(IAX,IBY) = HESSNN(IAX,IBY)
     *                                        + CRSCHI*SBY*VNUCXY
                           END IF
                           IF (IAX*IBZ.NE.0) THEN
                              HESSNN(IAX,IBZ) = HESSNN(IAX,IBZ)
     *                                        + CRSCHI*SBZ*VNUCXZ
                           END IF
                           IF (IAY*IBX.NE.0) THEN
                              HESSNN(IAY,IBX) = HESSNN(IAY,IBX)
     *                                        + CRSCHI*SBX*VNUCXY
                           END IF
                           IF (IAY*IBY.NE.0) THEN
                              HESSNN(IAY,IBY) = HESSNN(IAY,IBY)
     *                                        + CRSCHI*SBY*VNUCYY
                           END IF
                           IF (IAY*IBZ.NE.0) THEN
                              HESSNN(IAY,IBZ) = HESSNN(IAY,IBZ)
     *                                        + CRSCHI*SBZ*VNUCYZ
                           END IF
                           IF (IAZ*IBX.NE.0) THEN
                              HESSNN(IAZ,IBX) = HESSNN(IAZ,IBX)
     *                                        + CRSCHI*SBX*VNUCXZ
                           END IF
                           IF (IAZ*IBY.NE.0) THEN
                              HESSNN(IAZ,IBY) = HESSNN(IAZ,IBY)
     *                                        + CRSCHI*SBY*VNUCYZ
                           END IF
                           IF (IAZ*IBZ.NE.0) THEN
                              HESSNN(IAZ,IBZ) = HESSNN(IAZ,IBZ)
     *                                        + CRSCHI*SBZ*VNUCZZ
                           END IF
                           IF (IAX*IAX.NE.0) THEN
                              HESSNN(IAX,IAX) = HESSNN(IAX,IAX) - VNUCXX
                           END IF
                           IF (IAY*IAX.NE.0) THEN
                              HESSNN(IAY,IAX) = HESSNN(IAY,IAX) - VNUCXY
                           END IF
                           IF (IAY*IAY.NE.0) THEN
                              HESSNN(IAY,IAY) = HESSNN(IAY,IAY) - VNUCYY
                           END IF
                           IF (IAZ*IAX.NE.0) THEN
                              HESSNN(IAZ,IAX) = HESSNN(IAZ,IAX) - VNUCXZ
                           END IF
                           IF (IAZ*IAY.NE.0) THEN
                              HESSNN(IAZ,IAY) = HESSNN(IAZ,IAY) - VNUCYZ
                           END IF
                           IF (IAZ*IAZ.NE.0) THEN
                              HESSNN(IAZ,IAZ) = HESSNN(IAZ,IAZ) - VNUCZZ
                           END IF
                           IF (IBX*IBX.NE.0) THEN
                              HESSNN(IBX,IBX) = HESSNN(IBX,IBX) - VNUCXX
                           END IF
                           IF (IBY*IBX.NE.0) THEN
                              HESSNN(IBY,IBX) = HESSNN(IBY,IBX)
     *                                        - SBX*SBY*VNUCXY
                           END IF
                           IF (IBY*IBY.NE.0) THEN
                              HESSNN(IBY,IBY) = HESSNN(IBY,IBY) - VNUCYY
                           END IF
                           IF (IBZ*IBX.NE.0) THEN
                              HESSNN(IBZ,IBX) = HESSNN(IBZ,IBX)
     *                                        - SBX*SBZ*VNUCXZ
                           END IF
                           IF (IBZ*IBY.NE.0) THEN
                              HESSNN(IBZ,IBY) = HESSNN(IBZ,IBY)
     *                                        - SBY*SBZ*VNUCYZ
                           END IF
                           IF (IBZ*IBZ.NE.0) THEN
                              HESSNN(IBZ,IBZ) = HESSNN(IBZ,IBZ) - VNUCZZ
                           END IF
                        END IF
  500                   CONTINUE
                     END IF
                  END IF
  400             CONTINUE
               END IF
  300       CONTINUE
         END IF
  200 CONTINUE
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('Nuclear repulsion gradient',-1)
         CALL PRIGRD(GRADNN,CSTRA,SCTRA)
         IF (SECOND) THEN
            CALL HEADER('Nuclear repulsion Hessian',-1)
            CALL PRIHES(HESSNN,'CENTERS',CSTRA,SCTRA)
         END IF
      END IF
      IF (SECOND) CALL ADDHES(HESSNN)
      IF (IPRINT .GT. 1) CALL TIMER('NUCREP',TIMSTR,TIMEND)
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Program stopped in ABACUS after NUCREP as requested.'
         CALL QUIT(' ***** End of ABACUS (in NUCREP) *****')
      END IF
      CALL QEXIT('NUCREP')
      RETURN
      END
