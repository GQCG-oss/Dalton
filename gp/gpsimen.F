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
#if defined (SYS_CRAY)
      REAL
#else
      DOUBLE PRECISION
#endif
     &FUNCTION GETERF(x)
#include <implicit.h>
#include <pi.h>
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
