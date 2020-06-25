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
C  /* Deck cino */
#if defined (VAR_CINO)
CMAERKE- 900221-hjaaj: no calls of CINO yet
      SUBROUTINE CINO(INDXCI,WRK,LFRSAV)
C
C Nov 89 - PJ
C Aug 90 - HJAaJ
C PURPOSE: DETERMINE A SET OF NATURAL ORBITALS FROM A CI CALCULATION.
C          THE OCCUPATION NUMBERS OF THE NATURAL ORBITALS
C          ARE USED TO DETERMINE THE ACTIVE ORBITAL
C          SPACE OF A MCSCF AND THE NATURAL ORBITALS ARE USED
C          AS STARTING ORBITALS OF THE MCSCF CALCULATION.
C
#include "implicit.h"
C
      DIMENSION INDXCI(*),WRK(*)
C
#include "dummy.h"
      PARAMETER ( D2 = 2.0D0 )
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "infvar.h"
#include "infpri.h"
C
      CHARACTER*8 TABLE(2)
C
      DATA TABLE/'********', 'OLDORB  '/
C
      CALL QENTER('CINO  ')
      TICINO = SECOND()
      IPRCI  = IPRSIR
      LFREE  = LFRSAV
C
      IF (IPRCI .GE. 4) WRITE (LUW4,'(//A)')
     *   ' ----- OUTPUT FROM SIRIUS.CINO MODULE -----'
C
      LRCINO = LROOTS
      IF (ICHECK .NE. 0) THEN
         LRCINO = MAX(LRCINO,NROOTS+2)
         LRCINO = MIN(NCONF,LRCINO)
      END IF
C
      KCMO   = 1
      KW21   = KCMO  + NCMOT
      LW21   = LFREE - KW21
      KECI   = KW21
      KICROO = KECI    + LRCINO + 2
      KW22   = KICROO  + LRCINO + 2
      LW22   = LFREE   - KW22
      IF (LW22 .LT. 0) CALL ERRWRK(' CINO',-KW22,LFREE)
c     MAERKE 900801-hjaaj: call setsir and getcix ???
      CALL READMO(WRK(KCMO),9)
      CALL NEWIT1
      WRITE (LUIT1) TABLE(1),TABLE(1),TABLE(1),TABLE(2)
      NCMOT4 = MAX(4,NCMOT)
      CALL WRITT(LUIT1,NCMOT4,WRK(KCMO))
c     MAERKE 900801-hjaaj: Integral transformation sometimes needed !!
      IF (MAXCIT .LE. 0) THEN
         MAXCITX = 12
      ELSE
         MAXCITX = MAXCIT
      END IF
      CALL CICTL(1,LRCINO,MAXCITX,THRGRD,WRK(KCMO),INDXCI,
     *           WRK(KECI),ICONV,WRK(KICROO),WRK(KW22),LW22)
C     CALL CICTL (ICICTL,NCROOT,MAXITC,THRCIX,CMO,INDXCI,
C    *            ECI,ICONV,ICROOT,WRK,LWRK)
C
C CALL CICHCK TO REMOVE CI VECTORS THAT DO NOT HAVE THE SAME
C SYMMETRY AS ISTACI (ICHECK = 1) OR AS THE STATE OF THE LOWEST
C SYMMETRY (ICHECK = 2). (can only be done after several ci iterations
C have been performed.)  In CINO this is only of interest if this is
C an excited state calculation (ISTACI .gt. 1).
C
      IF (ICHECK .NE. 0 .AND. ISTACI .GT. 1) THEN
         IF (MAXITCX .GT. 2) THEN
            CALL CICHCK(WRK(KW21),LW21,ICHECK)
         ELSE
            WRITE(LUW4,'(//2A/A/A,I5)')
     *      ' --- WARNING CINO ---  REQUESTED SYMMETRY CHECK OF',
     *      ' CI VECTOR NOT PERFORMED',
     *      ' CURRENT ALGORITHM ONLY WORKS AFTER SEVERAL CI ITERATIONS',
     *      ' MAXIMUM CI ITERATIONS   :',MAXCITX
         END IF
      END IF
C
C READ IN CI VECTOR
C
      KCREF  = KW21
      KDV    = KCREF + NCONF
      KUDV   = KDV   + NNASHX
      KRHO1  = KUDV  + N2ASHX
      KW3    = KRHO1 + N2ORBT
      LW3    = LFREE - KW3
      REWIND LUIT1
      CALL MOLLAB('STARTVEC',LUIT1,lupri)
      DO 920 I = 1,ISTACI-1
         READ (LUIT1)
 920  CONTINUE
      CALL READT(LUIT1,NCONF,WRK(KCREF))
      CALL MAKDV(WRK(KCREF),WRK(KDV),INDXCI,WRK(KW3),LW3)
      CALL DSPTSI(NASHT,WRK(KDV),WRK(KUDV))
      CALL DZERO( WRK(KRHO1),N2ORBX)
      DO 300 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NISHI = NISH(ISYM)
         DO 400 I = IORBI+1,IORBI+NISHI
            WRK(KRHO1-1+(I-1)*NORBT+I) = D2
 400     CONTINUE
         IASHI = IASH(ISYM)
         NASHI = NASH(ISYM)
         CALL MCOPY(NASHI,NASHI,WRK(KUDV+IASHI*(NASHT+1)),NASHT,
     *              WRK(KRHO1+(IORBI+NISHI)*(NORBT+1)),NORBT)
C        CALL MCOPY(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
 300  CONTINUE
C
C DETERMINE NATURAL ORBITALS AND WRITE ORBITALS ON LUIT1
C (reuse MP2_NATORB from MP2 module here)
C
      CALL MP2_NATORB(IPRCI,WRK(KRHO1),WRK(KCMO),WRK(KW3),LW3)
      CALL NEWORB('CINOSAVE',WRK(KCMO),.TRUE.)
C
      TICINO = SECOND() - TICINO
      IF (IPRSTAT .GT. 0) WRITE (LUSTAT,1860) TICINO
      IF (IPRCI   .GE. 4) WRITE (LUW4,1860) TICINO
 1860 FORMAT (/' Time used for CI natural orbitals :',F12.3,
     *        ' CPU seconds.')
      CALL FLSHFO(LUSTAT)
      CALL FLSHFO(LUW4)
      CALL QEXIT('CINO  ')
      RETURN
      END
#endif
