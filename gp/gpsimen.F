C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
C
C FUNCTION GETERF created 12.01.2005 by Simen S. Reine
C Uses the GETGAM subroutine (Boys function F_n(x) ) to calculate 
C the error function erf(x) = 2*x*F_0(x*x)/sqrt(pi)
C NOTE: Tabulated values (see GETGAM) must be set up for GETGAM to 
C return right values
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
      TABERF=fjws(1,0)*2.d0*x/SQRTPI
      RETURN
      END
