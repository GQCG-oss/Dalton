C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
#if defined (SYS_CRAY) || defined (SYS_T3D)
      NIBUF = 1
      IF (NBASIS .LE. 255) THEN
         NBITS = 8
      ELSE
         NBITS = 16
      END IF
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
#if defined (SYS_CRAY) || defined (SYS_T90) || (defined (SYS_AIX) && defined (VAR_INT64))
         LRECL = LBFINP + LBFINP + 1
#endif
#if (defined (SYS_AIX) && !defined (VAR_INT64)) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_LINUX)
         IF (NBASIS .LE. 255) THEN
            LRECL = 2*LBFINP + LBFINP + 1
         ELSE
            LRECL = 2*LBFINP + 2*LBFINP + 1
         END IF
#endif
#if defined (SYS_IRIX) || defined (SYS_DEC)
         IF (NBASIS .LE. 255) THEN
            LRECL = 2*LBFINP + LBFINP + 1
         ELSE
            LRECL = 2*LBFINP + 2*LBFINP + 1
         END IF
#endif
#if defined (SYS_NEC)
         LRECL = LBFINP + LBFINP/IRAT + 1
#endif
         CALL GPOPEN(LUAORC(0),'AO2DIS00','UNKNOWN','DIRECT',
     &        'UNFORMATTED',LRECL,OLDDX)
         IF (U21INT) THEN
            LU21INT = 0
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
