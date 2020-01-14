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
      REAL*8 FUNCTION GETERF(x)
#include "implicit.h"
#include "pi.h"
      DIMENSION rj(1,0:0), wvalu(1), fjws(1,0:0), indads(1,3)
      DIMENSION rexpw(1), wvals(1,3), indadr(1)
      NODS = 1
      indadr(1) = 1
      wvalu(1)   = x*x 
      rj(1,0) = 0
      JMAX = 0
      NUABCD = 1
      CALL GETGAM(NODS,INDADR,wvalu,rj,JMAX,NUABCD,fjws,indads,wvals,
     &            rexpw,0)
      GETERF=fjws(1,0)*2.d0*x/SQRTPI
      RETURN
      END
