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
C  /* Deck herdi1 */
      SUBROUTINE HERDI1(WORK,LWORK,IPRFCK)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
      LOGICAL   ABA, OLDDX
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include <ccinftap.h>
#include <inforb.h>
#include <cbitwo.h>
#include <nuclei.h>
#include <dorps.h>
#include <inftap.h>
#include <symmet.h>
#include <blocks.h>
#include <eribuf.h>
C
C
C     Inquire if this is an ABACUS run or not:
C
      CALL ABARUN(ABA)
C
C----------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C----------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
C
         IPRALL =  0
C
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
C
C
         IPRINT = IPRFCK
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
C
      IF (LUINTR .LE. 0) CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUINTR
C
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KJORBS = KIORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KKORBS = KJORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KLAST  = KKORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('HERDI1','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
C
      RETURN
      END
C  /* Deck herdi2 */
      SUBROUTINE HERDI2(WORK,LWORK,INDEXA,ISHELA,NAINTS,IPRFCK)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>

CCOLD LOGICAL DOINDX, RELCAL, RETUR, TKTIME, NOCONT, NOPV, NODV
      LOGICAL DOINDX, RELCAL
      INTEGER NUMDIS
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
#include <blocks.h>
#include <eribuf.h>
#include <r12int.h>
C
      IF (LUAORC(0) .LE. 0) THEN
         LBFINP = LBUF
C
#if defined (VAR_INT64)
         LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
         LRECL = LBFINP + LBFINP/IRAT + 1
#else
         IF (NBASIS .LE. 255) THEN
            LRECL = 2*LBFINP + LBFINP + 1
         ELSE
            LRECL = 2*LBFINP + 2*LBFINP + 1
         END IF
#endif
         CALL GPOPEN(LUAORC(0),'AO2DIS00','UNKNOWN','DIRECT',
     &        'UNFORMATTED',LRECL,OLDDX)
         IF (U21INT) THEN
            LU21INT = -1
            CALL GPOPEN(LU21INT,'AOTDIS00','UNKNOWN','DIRECT',
     &                  'UNFORMATTED',LRECL,OLDDX)
         END IF
      END IF
C
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KJORBS = KIORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KKORBS = KJORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KLAST  = KKORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('HERDI2','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      ITYPE  = 7
      MAXDIF = 0
      JPRINT = 0
      NBUFX(0) = 0
      IPRINT = 0
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      NDMAT  = 0
      IREPDM = 0
      IFCTYP = 0
      INDXAB = 0
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      I2TYP  = 0
      ICEDIF = 0
      IFTHRS = 0
      DOINDX = .TRUE.
      NEWDIS = .TRUE.
      RETUR  = .FALSE.
      TKTIME = .FALSE.
      NOCONT = .FALSE.
      NOPV   = .FALSE.
      NODV   = .FALSE.
      RELCAL = .FALSE.
C
      REWIND(LUINTR)
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
C
      CALL AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,WORK(KIORBS),JPRINT)
C
      NUMDIS = NAINTS
C
      IF (IPRFCK.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (LUINTA .LE. 0) CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      CALL TWOINT(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,IREPDM,
     &            IFCTYP,GMAT,INDXAB,NUMDIS,MAXDIS,
     &            ITYPE,MAXDIF,JATOM,NODV,
     &            NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &            IPRNTD,RETUR,ISHELA,I2TYP,
     &            WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL)
      IF (IPRFCK.GT.0) CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
      IF (U21INT) CALL GPCLOSE(LU21INT,'KEEP')
C
      RETURN
      END
C
C  /* Deck nerdi2 */
      SUBROUTINE NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHELA,ISHELB,NAINTS,
     &                  NBINTS,IPRFCK,SETUP,IOPT)
C
C     Option input :   IOPT = 1  :  Calculate ordenary integral distributions 
C                                   with two fixed indices.
C                           = 2  :  Calculate exchange type integrals entering
C                                   the diagonal elements.
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
#include <choles.h>

      LOGICAL DOINDX, RELCAL,SETUP
      INTEGER NUMDIS
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
#include <blocks.h>
#include <eribuf.h>
C
      IF (LUAORC(0) .LE. 0) THEN
         LBFINP = LBUF
C
#if defined (VAR_INT64)
         LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
         LRECL = LBFINP + LBFINP/IRAT + 1
#else
         IF (NBASIS .LE. 255) THEN
            LRECL = 2*LBFINP + LBFINP + 1
         ELSE
            LRECL = 2*LBFINP + 2*LBFINP + 1
         END IF
#endif
         CALL GPOPEN(LUAORC(0),'AO2DIS00','UNKNOWN','DIRECT',
     &        'UNFORMATTED',LRECL,OLDDX)
      END IF
C
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KJORBS = KIORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KKORBS = KJORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KLAST  = KKORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('NERDI2','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      ITYPE  = 7
      MAXDIF = 0
      JPRINT = 0
      NBUFX(0) = 0
      IPRINT = 0
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      NDMAT  = 0
      IREPDM = 0
      IFCTYP = 0
      INDXAB = 0
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      I2TYP  = 0
      ICEDIF = 0
      IFTHRS = 0
      DOINDX = .TRUE.
      NEWDIS = .TRUE.
      RETUR  = .FALSE.
      TKTIME = .FALSE.
      NOCONT = .FALSE.
      NOPV   = .FALSE.
      NODV   = .FALSE.
      RELCAL = .FALSE.
C
      ISHLB = ISHELB
C
      IF     (IOPT .EQ. 1) THEN
         I2TYP = 4
      ELSEIF (IOPT .EQ. 2) THEN
         I2TYP = 5
c        LSAVE = LBUF
c        LBUF  = 50 000
      ELSE
          STOP 'Wrong I2TYP in NERDI2'
      ENDIF
C
      REWIND(LUINTR)
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
C
      CALL AINDEX(ISHELB,NBINTS,INDEXB,DOINDX,WORK(KIORBS),JPRINT)
      CALL AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,WORK(KIORBS),JPRINT)
C
      IF (SETUP) RETURN
C
Casm  NUMDIS = NAINTS
C
      IF (IPRFCK.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (LUINTA .LE. 0) CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      CALL TWOINT(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,IREPDM,
     &            IFCTYP,GMAT,INDXAB,NUMDIS,MAXDIS,
     &            ITYPE,MAXDIF,JATOM,NODV,
     &            NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &            IPRNTD,RETUR,ISHELA,I2TYP,
     &            WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL)
      IF (IPRFCK.GT.0) CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
C
      RETURN
      END

