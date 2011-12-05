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
C  /* Deck sym2 */
      SUBROUTINE SYM2(SO,AO,IPNTAO,IPNTOP,ISYMR,ISYMT,ISYMTS,FACTOR,
     &                SQ12EL,SOP000,NINTX,IPRINT)
C
C     Take block of distinct AO two-electron integral derivatives and
C     generate symmetrized contributions to SO integral derivatives
C
C                                          880429   PRT & TUH
C
C     Modified for triangular looping in components tuh
C     Loop over operator irreps outside loop over components 880928 tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      LOGICAL SQ12EL, SOP000
      DIMENSION AO(NOABCD,*), SO(NOABCD,*), IPNTAO(*),
     &          IPNTOP(3,*)
#include <twocom.h>
#include <symmet.h>
#include <ibtfun.h>
      IF (IPRINT .GT. 5) THEN
        CALL TITLER('Output from SYM2','*',103)
        IF (IPRINT .GT. 10) THEN
          WRITE (LUPRI,'(2X,A, I5)')'NINTX     ',NINTX
          WRITE (LUPRI,'(2X,A,3I5)')'ISYMR/T/TS',ISYMR,ISYMT,ISYMTS
          WRITE (LUPRI,'(2X,A, L5)')'SOP000    ',SOP000
          WRITE (LUPRI,'(2X,A, I5)')'NOABCD    ',NOABCD
          WRITE (LUPRI,'(2X,A,F12.6)')'FACTOR    ',FACTOR
        END IF
      END IF
      IF (SOP000) THEN
         DO 100 I = 1, NINTX
            IAO = IPNTAO(I)
#if defined (SYS_AIX)
            ITMPR  = IBTAND(ISYMR, IPNTOP(1,I))
            ITMPT  = IBTAND(ISYMT, IPNTOP(2,I))
            ITMPTS = IBTAND(ISYMTS,IPNTOP(3,I))
            ITMP1  = IBTXOR(ITMPT,ITMPTS)
            FAC    = FACTOR*PT(IBTXOR(ITMPR,ITMP1))
#else
            FAC = FACTOR*PT(IBTXOR((IBTAND(ISYMR, IPNTOP(1,I))),
     &                     (IBTXOR (IBTAND(ISYMT, IPNTOP(2,I)),
     &                              IBTAND(ISYMTS,IPNTOP(3,I))))))
#endif
            DO 110 J = 1,NOABCD
               SO(J,I) = FAC*AO(J,IAO)
  110       CONTINUE
  100    CONTINUE
      ELSE
         DO 200 I = 1, NINTX
            IAO = IPNTAO(I)
#if defined (SYS_AIX)
            ITMPR  = IBTAND(ISYMR, IPNTOP(1,I))
            ITMPT  = IBTAND(ISYMT, IPNTOP(2,I))
            ITMPTS = IBTAND(ISYMTS,IPNTOP(3,I))
            ITMP1  = IBTXOR(ITMPT,ITMPTS)
            FAC    = FACTOR*PT(IBTXOR(ITMPR,ITMP1))
#else
            FAC = FACTOR*PT(IBTXOR((IBTAND(ISYMR, IPNTOP(1,I))),
     &                     (IBTXOR (IBTAND(ISYMT, IPNTOP(2,I)),
     &                              IBTAND(ISYMTS,IPNTOP(3,I))))))
#endif
            DO 210 J = 1,NOABCD
               SO(J,I) = SO(J,I) + FAC*AO(J,IAO)
  210       CONTINUE
  200    CONTINUE
      END IF
      IF (IPRINT .GT. 20) THEN
         DO 300 I = 1, NINTX
            WRITE (LUPRI,'(/A,I3,A,5F12.6,/,(12X,F12.6,/))')
     &         ' Unf. SO',I,':',(SO(J,I),J=1,NOABCD)
  300    CONTINUE
      END IF
      RETURN
      END
C  /* Deck drsym2 */
      SUBROUTINE DRSYM2(SOINT,AOINT,WORK,NCCINT,LWORK,HKABCD,ISYMR,
     &                  ISYMT,ISYMTS,IPNTAO,IPNTOP,MULE,SQ12EL,SOP000,
     &                  IPRINT)
C
C     This superb piece of software is a joint effort by
C     TUH and PRT, 880429
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      LOGICAL SQ12EL, SOP000
      DIMENSION AOINT(NCCINT,*), SOINT(NOABCD,*), WORK(LWORK),
     &          IPNTAO(NINTMX,*), IPNTOP(3,NINTMX,*)
#include <twocom.h>
#include <twosta.h>
#include <dirprt.h>
#include <expcom.h>
#include <symmet.h>
#include <dorps.h>
#include <doxyz.h>
#include <ibtfun.h>
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DRSYM2',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
      IF (NOATMS .GT. 1) THEN
         KLAST = NCCINT
         IF (KLAST .GT. LWORK) CALL STOPIT('DRSYM2',' ',KLAST,LWORK)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
         LWTOT  = LWTOT - KLAST
         MWDRSY = MAX(MWDRSY,KLAST)
      END IF
      ISTR   = 1
      DO 100 ICOOR = 1, 3
      IF (DOXYZ(ICOOR)) THEN
         ISYTYE = ISYMAX(ICOOR,1)
         DO 200 IREPE = 0, MAXREP
            IF (DOREPS(IREPE) .AND.
     &          (IBTAND(MULE,IBTXOR(IREPE,ISYTYE)) .EQ. 0)) THEN
               IREPX = IPTREP(IREPE,2)
               NINTX = NINTSR(IREPX)
               IF (NINTX .GT. 0) THEN
                  IF (NOATMS .EQ. 1) THEN
                     FACSYM = PT(IBTAND(ISOPDR(1),IREPE))
     &                        *SIGNDR(1,ICOOR)*HKABCD
                     CALL SYM2(SOINT(1,ISTR),AOINT(1,ICOOR),
     &                         IPNTAO(1,IREPX),IPNTOP(1,1,IREPX),
     &                         ISYMR,ISYMT,ISYMTS,FACSYM,SQ12EL,
     &                         SOP000,NINTX,IPRINT)
                  ELSE
                     FACSYM = HKABCD
                     DO 300 ICENT = 1, NOATMS
                        FACTOR = PT(IBTAND(ISOPDR(ICENT),IREPE))
     &                           *SIGNDR(ICENT,ICOOR)
                        ITYPE   = 3*(ICENT - 1) + ICOOR
                        IF (ICENT .EQ. 1) THEN
                           DO 400 INT = 1, NCCINT
                              WORK(INT) = FACTOR*AOINT(INT,ITYPE)
  400                      CONTINUE
                        ELSE
                           DO 450 INT = 1, NCCINT
                              WORK(INT) = WORK(INT)
     &                                  + FACTOR*AOINT(INT,ITYPE)
  450                      CONTINUE
                        END IF
  300                CONTINUE
                     CALL SYM2(SOINT(1,ISTR),WORK,IPNTAO(1,IREPX),
     &                         IPNTOP(1,1,IREPX),ISYMR,ISYMT,ISYMTS,
     &                         FACSYM,SQ12EL,SOP000,NINTX,IPRINT)
                  END IF
               END IF
               ISTR = ISTR + NINTX
               IF (IPRINT .GT. 20) THEN
                   WRITE (LUPRI,'(A,3I5)') ' ICOOR, IREPE, ISYTYE ',
     &                                       ICOOR, IREPE, ISYTYE
                   WRITE (LUPRI,'(A,F12.6)') ' FACSYM ', FACSYM
               END IF
            END IF
  200    CONTINUE
      END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck sposym */
      SUBROUTINE SPOSYM(SOINT,AOINT,NCCINT,HKABCD,ISYMR,ISYMT,ISYMTS,
     &                  IPNTAO,IPNTOP,SQ12EL,SOP000,IPRINT)
C
C     Driver routine for symmetrization of spin-orbit integrals
C     TUH and OV, 891111
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      LOGICAL SQ12EL, SOP000
      INTEGER X
      DIMENSION AOINT(NCCINT,*), SOINT(NOABCD,*),
     &          IPNTAO(NINTMX,*), IPNTOP(3,NINTMX,*)
#include <twocom.h>
#include <dirprt.h>
#include <expcom.h>
#include <symmet.h>
#include <ibtfun.h>
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SPOSYM',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
      ISOSTR = 1
      IRPXYZ = IBTXOR(ISYMAX(1,1),IBTXOR(ISYMAX(2,1),ISYMAX(3,1)))
      DO 100 X = 1, 3
         IREPE = IBTXOR(ISYMAX(X,1),IRPXYZ)
         IREPX = IPTREP(IREPE,2)
         NINTX = NINTSR(IREPX)
         IF (IPRINT .GT. 20) THEN
            WRITE (LUPRI,'(A,3I5)') ' X, IREPE ', X, IREPE
         END IF
         CALL SYM2(SOINT(1,ISOSTR),AOINT(1,X),IPNTAO(1,IREPX),
     &             IPNTOP(1,1,IREPX),ISYMR,ISYMT,ISYMTS,
     &             HKABCD,SQ12EL,SOP000,NINTX,IPRINT)
         ISOSTR = ISOSTR + NINTX
  100 CONTINUE
      RETURN
      END
C  /* Deck mgsym2 */
      SUBROUTINE MGSYM2(SOINT,AOINT,WORK,NCCINT,LWORK,HKABCD,ISYMR,
     &                  ISYMT,ISYMTS,IPNTAO,IPNTOP,SOP000,IPRINT)
C
C     tuh Sep 7 92
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (DP5 = 0.5D0)
      INTEGER X, Y, Z
      LOGICAL SOP000
      DIMENSION AOINT(NCCINT,3,*), SOINT(*), WORK(NCCINT,*),
     &          IPNTAO(NINTMX,*), IPNTOP(3,NINTMX,*)
      DIMENSION DIFAB(3), DIFCD(3)
#include <twocom.h>
#include <twosta.h>
#include <symmet.h>
#include <ibtfun.h>
      XAND(I) = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I) = PT(IBTAND(ISYMAX(3,1),I))
      NEXT(I) = MOD(I,3) + 1
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine MGSYM2',-1)
         WRITE (LUPRI,'(A,3I5)') ' ISYMR, ISYMT, ISYMTS ',
     &                             ISYMR, ISYMT, ISYMTS
      END IF
C
      KLAST = 2*NCCINT
      IF (KLAST .GT. LWORK) CALL STOPIT('MGSYM2',' ',KLAST,LWORK)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
      LWTOT  = LWTOT - KLAST
      MWDRSY = MAX(MWDRSY,KLAST)
C
      DIFAB(1) =             CORAX0 - XAND(ISYMR )*CORBX0
      DIFAB(2) =             CORAY0 - YAND(ISYMR )*CORBY0
      DIFAB(3) =             CORAZ0 - ZAND(ISYMR )*CORBZ0
      DIFCD(1) = XAND(ISYMT)*CORCX0 - XAND(ISYMTS)*CORDX0
      DIFCD(2) = YAND(ISYMT)*CORCY0 - YAND(ISYMTS)*CORDY0
      DIFCD(3) = ZAND(ISYMT)*CORCZ0 - ZAND(ISYMTS)*CORDZ0
      FACSYM   = DP5*HKABCD
C
      IADR = 1
      DO 100 X = 1, 3
         IREPE = ISYMAX(X,2)
         IREPX = IPTREP(IREPE,2)
         NINTX = NINTSR(IREPX)
         IF (NINTX .GT. 0) THEN
            Y   = NEXT(X)
            Z   = NEXT(Y)
            ABY = DIFAB(Y)
            ABZ = DIFAB(Z)
            CDY = DIFCD(Y)
            CDZ = DIFCD(Z)
            DO 200 I = 1, NCCINT
               AOAB = ABY*AOINT(I,Z,1) - ABZ*AOINT(I,Y,1)
               AOCD = CDY*AOINT(I,Z,2) - CDZ*AOINT(I,Y,2)
               WORK(I,1) = AOAB + AOCD
               WORK(I,2) = AOAB - AOCD
  200       CONTINUE
            CALL SYM2(SOINT(IADR),WORK(1,1),IPNTAO(1,IREPX),
     &                IPNTOP(1,1,IREPX),ISYMR,ISYMT,ISYMTS,FACSYM,
     &                .FALSE.,SOP000,NINTX,IPRINT)
            CALL SYM2(SOINT(IADR+NOABCD*NINTX),WORK(1,2),
     &                IPNTAO(1,IREPX),IPNTOP(1,1,IREPX),ISYMR,ISYMT,
     &                ISYMTS,FACSYM,.FALSE.,SOP000,NINTX,IPRINT)
         END IF
         IADR = IADR + 2*NOABCD*NINTX
  100 CONTINUE
      RETURN
      END
