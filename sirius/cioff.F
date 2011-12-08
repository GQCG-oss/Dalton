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
C  /* Deck cioff */
      SUBROUTINE CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
c Generate offset arrays for CI vector with symmetry IREFSM
c ICORHC = 1 => store in C arrays
c ICORHC = 2 => store in HC arrays
c
#include <implicit.h>
#include <detbas.h>
      DIMENSION XNDXCI(*)
c
      IF(ICORHC.EQ.1) THEN
        CALL CIOFFS(IREFSM,XNDXCI(KCOFF),XNDXCI(KICSO),XNDXCI(KICOOS),
     &       XNDXCI(KCDTAS),XNDXCI,NTEST)
      ELSE
        CALL CIOFFS(IREFSM,XNDXCI(KHCOFF),XNDXCI(KIHCSO),XNDXCI(KIHOOS),
     &       XNDXCI(KHDTAS),XNDXCI,NTEST)
      END IF
c
      RETURN
      END
C  /* Deck cioffs */
      SUBROUTINE CIOFFS(IREFSM,IOFF,ISO,IOOS,NDTAS,XNDXCI,NTEST)
c
c Obtain symmetry blocking of CI vector of symmetry IREFSM
c
c Specific input
c IREFSM : Symmetry of CI vector
c.General input
c XNDXCI : CI information vector
c.Output
c For CAS
c  IOFF : Address of first element with a-string of given symmetry
c For RAS
c  ISO  : Address of first element with a-string of given symmetry
c  IOOS : Address of first element of given SOO type
c  NDTAS: Number of dets with a-strings of given symmetry
c
#include <implicit.h>
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
#include <inforb.h>
#include <ciinfo.h>
c
      DIMENSION XNDXCI(*)
      DIMENSION IOFF(*),ISO(*),IOOS(NOCTPB,NOCTPA,MAXSYM),NDTAS(*)
c
      NTEST2 = NTEST
C            CIOFF2(IREFSM,MAXSYM,SYMPRO,IOCOC,NSASOA,NSASOB,
C    &                  NOCTPA,NOCTPB,IOFF,IOOS,NDTAS,NTEST)
        CALL CIOFF2(IREFSM,MAXSYM,MULD2H,XNDXCI(KIOCOC),XNDXCI(KNSSOA),
     &              XNDXCI(KNSSOB),NOCTPA,NOCTPB,ISO,IOOS,NDTAS,NTEST2)
c
      RETURN
      END
