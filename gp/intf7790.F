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
#ifdef UNDEF
/* Comdeck log */
011002-vebjornb: Interfaces between Fortran 90 code and the
                 old common blocks.
#endif
C  /* Deck gtunit */
      SUBROUTINE GTUNIT(ILPRI, ILCMD, ILERR, IUSER)
C
C     Subroutine to read unit numbers and the general print level
C
#include <implicit.h>
#include <priunit.h>
#include <gnrinf.h>
      ILPRI = LUPRI
      ILCMD = LUCMD
      ILERR = LUERR
      IUSER = IPRUSR
      RETURN
      END

C  /* Deck stunit */
      SUBROUTINE STUNIT(ILPRI, ILCMD, ILERR)
C
C     Subroutine to set unit numbers
C
#include <implicit.h>
#include <priunit.h>
      LUPRI = ILPRI
      LUCMD = ILCMD
      LUERR = ILERR
      RETURN
      END

C  /* Deck gterg */
      SUBROUTINE GTERG(ENERGY)
C
C     Returns the energy
C
#include <implicit.h>
#include <maxorb.h>
#include <gnrinf.h>
#include <infinp.h>
#include <infopt.h>
      IF ((DOMP2 .AND. .NOT. DOMC) .OR. DOCCSD .OR. ((DOCI .OR. DOCINO)
     &     .AND. .NOT. DOMC)) THEN
         ENERGY = ECORR
      ELSE
         ENERGY = EMCSCF
      END IF
      RETURN
      END

C  /* Deck stctrl */
      SUBROUTINE STCTRL(LWRI,LGRAD,LHESS,LHLFEY)
C
C     (Re)sets a number of ABACUS variables, intended to be used between
C     iterations (optimization or dynamics)
C
#include <implicit.h>
#include <mxcent.h>
#include <gnrinf.h>
#include <abainf.h>
      WRINDX = LWRI
      MOLGRD = LGRAD
      MOLHES = LHESS
      HELFEY = LHLFEY
      RETURN
      END

C  /* Deck streai */
      SUBROUTINE STREAI(RDIN,RMOL,HRIN)
C
C     (Re)sets a number of READIN variables, intended to be used between
C     iterations (optimization or dynamics)
C
#include <implicit.h>
#include <gnrinf.h>
      LOGICAL RDIN, RMOL, HRIN
      RDINPC = RDIN
      RDMLIN = RMOL
      HRINPC = HRIN
      RETURN
      END

C  /* Deck gtconv */
      SUBROUTINE GTCONV(TOJOUL,TOKAYS,DIRCON,FRMAMU,TFS2AU)
C
C     Returns conversion factors from the common block codata.h
C
#include <implicit.h>
#include <codata.h>
      TOJOUL = XTJ
      TOKAYS = XTKAYS
      DIRCON = HBAR
      FRMAMU = XFAMU
      TFS2AU = 1.D-15*TOJOUL/DIRCON
      RETURN
      END

C  /* Deck upmoli */
      SUBROUTINE UPMOLI(IATOMS,NCOORD,CRDGEO)
C
C     Update the molecular geometry (molecule input), intended to be
C     used between iterations (optimization or dynamics)
C
#include <implicit.h>
#include <mxcent.h>
#include <cbirea.h>
#include <molinp.h>
#include <nuclei.h>
      DIMENSION CRDGEO(NCOORD)
      CHARACTER*80 TMPLIN
      I = 5
      IF (BASIS) I = 6
      TMPLIN = MLINE(I)
      WRITE(TMPLIN(10:20), '(A11)') '0          '
      MLINE(I) = TMPLIN
      CALL WLKMOL(CRDGEO)
      CALL PNCMOL
      DO 100 I = 1, IATOMS
         DO 110 J = 1, 3
            CORD(J,I) = CRDGEO((I-1)*3+J)
 110     CONTINUE
 100  CONTINUE
      RETURN
      END

C  /* Deck rmvsym */
      SUBROUTINE RMVSYM
C
C     Removes (resets) the symmetry of the molecule, intended to be
C     used between iterations (optimization or dynamics)
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <gnrinf.h>
#include <huckel.h>
#include <nuclei.h>
#include <symmet.h>
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      NEWSYM = .TRUE.
      ADDSTO = .TRUE.
      CALL IZERO(NUCNUM,8*MXCENT)
      CALL IZERO(NCRREP,16)
      CALL IZERO(IPTCNT,48*MXCENT)
      CALL IZERO(NAXREP,16)
      RETURN
      END
