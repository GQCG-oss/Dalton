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
*=====================================================================*
C  /* Deck cc_bfbsort1 */
      SUBROUTINE CC_BFBSORT1(DSRHF,BSRHF,ISYRHF,LSQRAB)
*---------------------------------------------------------------------*
*
*     Purpose: presort DSRHF integral array for the BF intermediate
*              calculation in the B matrix transformation
*
*     DSRHF  : (alp bet|k delta) integrals for a fixed delta
*     BSRHF  : integrals sorted as I(alp k;bet)^del
*     ISYRHF : symmetry of the integral arrays DSRHF,BSRHF
*
*     Written by Christof Haettig July/October 1998
*     Updated by Sonia Coriani November 1999 to handle full (a b|
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <ccorb_dec.h>
#  include <maxorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <symsq_dec.h>
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <symsq.h>

      LOGICAL LSQRAB
      INTEGER ISYRHF, ISYM, ISYMAK, ISYBET, ISYMK, ISYMAB, ISYALP
      INTEGER ICOUNT, NBSRHF(8), IBSRHF(8,8)
      INTEGER NABK, NAKB, NAK, KOFF1, IJSQ

#if defined (SYS_CRAY)
      REAL DSRHF(*), BSRHF(*)
#else
      DOUBLE PRECISION DSRHF(*), BSRHF(*)
#endif
C
C     --------------------------------------
C     precalculate symmetry array for BSRHF:
C     --------------------------------------
C
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
C
C     -------------------
C     sort the integrals:
C     -------------------
C
      DO ISYMAK = 1, NSYM
      DO ISYMK  = 1, NSYM
C 
         ISYBET = MULD2H(ISYMAK,ISYRHF)
         ISYALP = MULD2H(ISYMK,ISYMAK)
         ISYMAB = MULD2H(ISYALP,ISYBET)
C
C        --------------------------------------------------------
C        get (alp k;bet) blocks out of (alp bet|k del) integrals:
C        --------------------------------------------------------
C
         DO K = 1, NRHF(ISYMK)
C
            IF (LSQRAB) THEN
              KOFF1  = IDSRHFSQ(ISYMAB,ISYMK) + N2BST(ISYMAB)*(K-1)
            ELSE
              KOFF1  = IDSRHF(ISYMAB,ISYMK) + NNBST(ISYMAB)*(K-1) 
            END IF
C
            DO A = 1, NBAS(ISYALP)
            DO B = 1, NBAS(ISYBET)
C
               IJSQ = IAODIS(ISYALP,ISYBET) + NBAS(ISYALP)*(B-1) + A
               IF (LSQRAB) THEN
                 NABK = KOFF1  + IJSQ                        !not quite sure
               ELSE
                 NABK = KOFF1  + IADRPK( I2BST(ISYMAB) + IJSQ )
               END IF
               NAK  = IT1AO(ISYALP,ISYMK)   + NBAS(ISYALP)*(K-1) + A
               NAKB = IBSRHF(ISYMAK,ISYBET) +NT1AO(ISYMAK)*(B-1) + NAK
C
               BSRHF(NAKB) = DSRHF(NABK)
C
            END DO
            END DO
C
         END DO
C
      END DO
      END DO
C
      RETURN
      END
*=====================================================================*
