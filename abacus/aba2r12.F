C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
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
#include <single.h>
C
C  /* Deck eriher */
      SUBROUTINE ABEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,
     &                  EXP12,EXP34,NUC12,NUC34,NUABCD,
     &                  NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
C
C     Written by Wim Klopper (University of Karlsruhe, 14 November 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DM2 = -2.0D0)
      DIMENSION Q(*), R(*), W(*), A(*), EXP12(NUC12), EXP34(NUC34),
     &          IODDHR(*), INDHER(*), PQX(*), PQY(*), PQZ(*)
      IODS = 0
      DO  IOD12 = 1, NUC12
         EXPP = EXP12(IOD12)
         DO  IOD34 = 1, NUC34
            IODS = IODS + 1
            EXPQ = EXP34(IOD34)
            W(IODS) = DM2*EXPP*EXPQ/(EXPP + EXPQ)
         END DO
      END DO
      CALL WKEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,NUABCD,
     &            NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
      RETURN
      END
