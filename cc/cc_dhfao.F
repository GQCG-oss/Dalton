C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
C     IF (ISYMH.NE.1 .OR. ISYMP.NE.1) THEN
C        WRITE (LUPRI,*) 'CC_DHFAO only implemented for '//
C    &        'total symmetric CMO.'
C        WRITE (LUPRI,*) 'ISYMH, ISYMP:',ISYMH,ISYMP
C        CALL QUIT('CC_DHFAO only implemented for total symmetric CMO.')
C     END IF
C
      ISYDMO = 1
      CALL CC_DENAO(AODEN,ISYDAO,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISYDMO,CMOP,ISYMP,
     *              CMOH,ISYMH,WORK(KEND1),LWRK1)
C
      RETURN
      END
