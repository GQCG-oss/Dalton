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
C  /* Deck cioff */
      SUBROUTINE CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
c Generate offset arrays for CI vector with symmetry IREFSM
c ICORHC = 1 => store in C arrays
c ICORHC = 2 => store in HC arrays
c
#include "implicit.h"
#include "detbas.h"
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
#include "implicit.h"
#include "mxpdim.h"
#include "detbas.h"
#include "strnum.h"
#include "inforb.h"
#include "ciinfo.h"
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
