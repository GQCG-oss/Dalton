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
C  /* Deck ccs_dhfao */
      SUBROUTINE CC_DHFAO(AODEN,ISYDAO,CMOP,ISYMP,CMOH,ISYMH,WORK,LWORK)
C
C     Purpose: To set up HF one electron AO-density matrix
C              allow for two different CMO vectors to handle
C              different density matrices needed for derivatives
C      
C                D_alp,bet = \sum_i  CMOP_alp,i CMOH_bet,i
C
C     Christof Haettig, spring 99, based on Asgers CCS_D1AO
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK), CMOP(*), CMOH(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KONEAI = 1
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KEND1  = KONEIA + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for work allocation '//
     &        'in CCS_D1AO')
      ENDIF
C
C--------------------------------------------------------------
C     Initialize arrays 
C--------------------------------------------------------------
C
      CALL DZERO(WORK(KONEAI),NT1AMX)
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEIA),NT1AMX)
C
C-----------------------
C     Set up MO-density.
C-----------------------
C
      DO 100 ISYM = 1,NSYM
         DO 110 I = 1,NRHF(ISYM)
C
            NII = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1) + I
C
            WORK(KONEIJ + NII - 1) = TWO
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------------
C     Transform density to AO basis.
C-----------------------------------
C
      ISYDEN = MULD2H(ISYMP,ISYMH)
C
      CALL DZERO(AODEN,N2BST(ISYDEN))
C
      IF (ISYMH.NE.1 .OR. ISYMP.NE.1) THEN
         WRITE (LUPRI,*) 'CC_DHFAO only implemented for '//
     &        'total symmetric CMO.'
         WRITE (LUPRI,*) 'ISYMH, ISYMP:',ISYMH,ISYMP
         CALL QUIT('CC_DHFAO only implemented for total symmetric CMO.')
      END IF
C
      ISYDMO = 1
      CALL CC_DENAO(AODEN,ISYDAO,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISYDMO,CMOP,ISYMP,
     *              CMOH,ISYMH,WORK(KEND1),LWRK1)
C
      RETURN
      END
