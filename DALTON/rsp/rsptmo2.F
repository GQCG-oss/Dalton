!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
      SUBROUTINE TM1IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT ONE-INDEX
C LINEAR RESPONSE EQUATIONS THAT ARE USED IN A
C A CALCULATION OF THIRD ORDER TRANSITION MOMENTS
C
#include "implicit.h"
C
#include "priunit.h"
#include "rspprp.h"
#include "infsmo.h"
#include "indcr.h"
#include "inforb.h"
#include "infrsp.h"
#include "infspi.h"
#include "inftmo.h"
#include "infcr.h"
C
C     Keep the accumulated number of ROOTS for three-photon calcs
C
      IDOFF = 0
C
      DO 200 IDSYM = 1,NSYM
      DO 300 ICSYM = 1,NSYM
      DO 400 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NTMCNV(IDSYM).GT.0) .AND. (NCTMOP(ICSYM).GT.0) .AND.
     &       (NBTMOP(IBSYM).GT.0) .AND. (NATMOP(IASYM).GT.0) ) THEN
C
C     If three-photon absorption is specified we only compute
C     linear response functions and excitation vectors for certain 
C     frequencies
C 
            IF (THREEPHOTON) THEN
               DO IDFR = 1,NTMCNV(IDSYM)
                  BTMFR(IDFR+IDOFF) = EXCIT2(IDSYM,IDFR)/3
                  CTMFR(IDFR+IDOFF) = EXCIT2(IDSYM,IDFR)/3
                  ATMFR = EXCIT2(IDSYM,IDFR) - 
     &                 CTMFR(IDFR+IDOFF) - BTMFR(IDFR+IDOFF)
                  INUM = INCRLR('EXCITLAB',EXCIT2(IDSYM,IDFR),IDSYM)
                  DO ICOP = 1,NCTMOP(ICSYM)
                     INUM = INCRLR(CTMLB(ICSYM,ICOP),-CTMFR(IDFR+IDOFF),
     &                    ICSYM)
                  END DO
                  DO IBOP = 1,NBTMOP(IBSYM)
                     INUM = INCRLR(BTMLB(IBSYM,IBOP),-BTMFR(IDFR+IDOFF),
     &                    IBSYM)
                  END DO
                  DO IAOP = 1,NATMOP(IASYM)
                     INUM = INCRLR(ATMLB(IASYM,IAOP),ATMFR,IASYM)
                  END DO
               END DO
            ELSE
C
C     The general case
C
               DO 450 IDFR = 1,NTMCNV(IDSYM)
                  INUM = INCRLR('EXCITLAB',EXCIT2(IDSYM,IDFR),IDSYM)
 450           CONTINUE
               DO 500 ICOP = 1,NCTMOP(ICSYM)
               DO 600 ICFR = 1,NCTMFR
                  INUM = INCRLR(CTMLB(ICSYM,ICOP),-CTMFR(ICFR),ICSYM)
 600           CONTINUE
 500           CONTINUE
               DO 700 IBOP = 1,NBTMOP(IBSYM)
               DO 800 IBFR = 1,NBTMFR
                  INUM = INCRLR(BTMLB(IBSYM,IBOP),-BTMFR(IBFR),IBSYM)
 800           CONTINUE
 700           CONTINUE
               DO 900  IDFR = 1,NTMCNV(IDSYM)
               DO 1000 ICFR = 1,NCTMFR
               DO 1100 IBFR = 1,NBTMFR
C
C     If harmonic generation is specified we only accept an equal
C     number of b,c frequencies
C
                  IF (CTMOHG) THEN
                     IF (NBTMFR.NE.NCTMFR) THEN
                        WRITE (LUPRI,'(2(/A))')
     * '*** ERROR -- If harmonic generation is specified we only',
     * '             accept an equal number of b,c frequencies.'
                      CALL QUIT('wrong frequencies')
                     END IF
                     IF (IBFR.NE.ICFR) GO TO 1100
                  END IF
                  ATMFR = EXCIT2(IDSYM,IDFR) - CTMFR(ICFR) - BTMFR(IBFR)
                  DO 1200 IAOP = 1,NATMOP(IASYM)
                     INUM = INCRLR(ATMLB(IASYM,IAOP),ATMFR,IASYM)
 1200             CONTINUE
 1100          CONTINUE
 1000          CONTINUE
 900           CONTINUE
            END IF
         END IF
 400  CONTINUE
 300  CONTINUE
      IDOFF = IDOFF + NTMCNV(IDSYM)
 200  CONTINUE
      RETURN
      END
      SUBROUTINE TM2IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C TWO-INDEX LINEAR RESPONSE EQUATIONS THAT NEED TO BE SOLVED
C IN A SECOND HYPERPOLARIZABILITY CALCULATION
C
#include "implicit.h"
C
      PARAMETER ( ZEROTHR = 1.0D-10 )
C
#include "rspprp.h"
#include "infcr.h"
#include "inforb.h"
#include "infrsp.h"
#include "infpri.h"
#include "infspi.h"
#include "inftmo.h"
#include "indcr.h"
C
      CHARACTER*8 DTMLB
C
C Put label EXCITLABin list for two-index vectors
C for vectors of the type N^BX and N^CX
C
      DATA DTMLB/'EXCITLAB'/
C
      DO 300 IDSYM = 1,NSYM
      DO 200 ICSYM = 1,NSYM
      DO 100 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NTMCNV(IDSYM).GT.0) .AND. (NCTMOP(ICSYM).GT.0)
     *     .AND. (NBTMOP(IBSYM).GT.0) .AND. (NATMOP(IASYM).GT.0)) THEN
            DO 210 ICOP = 1,NCTMOP(ICSYM)
            DO 110 IBOP = 1,NBTMOP(IBSYM)
               DO 320 IDFR = 1,NTMCNV(IDSYM)
               DO 220 ICFR = 1,NCTMFR
               DO 120 IBFR = 1,NBTMFR
                  DTMFR = EXCIT2(IDSYM,IDFR)
C
C     If three photon absorption calculation only certain frequencies
C     are considered.
C
                  IF (THREEPHOTON .AND. 
     &                ( (ABS(DTMFR-3*BTMFR(IBFR)).GT.ZEROTHR) 
     &                  .OR.
     &                  (ABS(DTMFR-3*CTMFR(ICFR)).GT.ZEROTHR) ) 
     &               ) GOTO 120
C
C  If harmonic generation is specified we only need two index
C  response vectors with equal B and C frequencies.
C
                  IF (CTMOHG) THEN
                     IF (IBFR.NE.ICFR) GO TO 120
                  END IF
C
                  INUM = INCR2(BTMLB(IBSYM,IBOP),CTMLB(ICSYM,ICOP),
     *                          -BTMFR(IBFR),-CTMFR(ICFR),IBSYM,ICSYM)
                  INUM = INCR2(BTMLB(IBSYM,IBOP),DTMLB,
     *                          -BTMFR(IBFR),DTMFR,IBSYM,IDSYM)
                  INUM = INCR2(CTMLB(ICSYM,ICOP),DTMLB,
     *                          -CTMFR(ICFR),DTMFR,ICSYM,IDSYM)
120            CONTINUE
220            CONTINUE
320            CONTINUE
110         CONTINUE
210         CONTINUE
         END IF
100   CONTINUE
200   CONTINUE
300   CONTINUE
      RETURN
      END

