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
C  /* Deck sbdorb */
      BLOCK DATA SBDORB
C     Initialize MULD2H in common block /INFORB/
#include <inforb.h>
C
C     MULTIPLICATION TABLE FOR SYMMETRIES (MULD2H is in /INFORB/)
C
      DATA MULD2H/1,2,3,4,5,6,7,8,
     *            2,1,4,3,6,5,8,7,
     *            3,4,1,2,7,8,5,6,
     *            4,3,2,1,8,7,6,5,
     *            5,6,7,8,1,2,3,4,
     *            6,5,8,7,2,1,4,3,
     *            7,8,5,6,3,4,1,2,
     *            8,7,6,5,4,3,2,1/
      END
C  /* Deck sbdtap */
      BLOCK DATA SBDTAP
C 900226-hjaaj
C Define file names for INFTAP
C Define LUSIFC label, LBSIFC, in INTAP
C To be safe, include EXTERNAL SBDTAP in main program.
#include <implicit.h>
#include <inftap.h>
      DATA LBSIFC /'SIR IPH '/
      END
C  /* Deck sbdgetd */
      BLOCK DATA SBDGETD
C 900302-hjaaj
C Define default values for CBGETD
C To be safe, include EXTERNAL SBDGETD in main program.
#include <implicit.h>
#include <cbgetdis.h>
      DATA DISTYP,IADINT,IADH2,IADH2X,IADH2D /1,4*-1/
      END
C  /* Deck sbdtra */
      BLOCK DATA SBDTRA
C 951130-hjaaj
C Initialize variables in /INFTRA/
C Initialization moved to here because NEWTRA may be set
C both in dalton and sirius input modules.
C To be safe, include EXTERNAL SBDTRA in sirctl
C
#include <implicit.h>
#include <inftra.h>
      DATA NEWTRA /.FALSE./
      END
C  /* Deck sbddim */
      BLOCK DATA SBDDIM
C 010321 hjaaj
C Initialize variables in /INFDIM/
C
#include <implicit.h>
#include <infdim.h>
      DATA MWORK /8 000 000/
C
C     Re MWORK:
C     default 8 mill. words for number of simultaneous AO distributions
C     in "old" integral transformation - defined for SORTA. The size
C     should not be bigger than resident memory, to avoid paging.
C     If MWORK too big, then SORTA may define more simult. dist. than
C     can be treated at later calls to TRACTL !!! /Mar 2001 hjaaj
C
      END
