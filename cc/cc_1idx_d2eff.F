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
*=====================================================================*
      SUBROUTINE CC_1IDX_D2EFF(BAODEN,ICON,G,ISYMG,D,ISYMD,ISYMQ,
     &                         DHFAO,ZKABAO,
     &                         B1DHFAO,B1KABAO,B2DHFAO,B2KABAO)
*---------------------------------------------------------------------*
*
*     Purpose: Add the extra terms to the "one-index" transformed
*              2-electron density matrix which originate from the 
*              orbital relaxation
*              
*           ICON    --  analogous to ICON in CC_D2EFF
*           DHFAO   --  usual Hatree-Fock density matrix
*           ZKABAO  --  relaxation contribution to 1-electron density
*           B1DHFAO --  DHFAO  with leading index transformed
*           B1KABAO --  ZKABAO with leading index transformed
*           B2DHFAO --  DHFAO  with second  index transformed
*           B2KABAO --  ZKABAO with second  index transformed
*           ISYMQ   --  symmetry of B1DHFAO, B1KABAO, B2DHFAO, B2KABAO
*
*     Christof Haettig, March 1999, based on Asgers CC_D2EFF routine
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER ISYMG, ISYMD, ISYMQ, ICON

#if defined (SYS_CRAY)
      REAL BAODEN(*)
      REAL DHFAO(*),  B1DHFAO(*), B2DHFAO(*)
      REAL ZKABAO(*), B1KABAO(*), B2KABAO(*)
      REAL ONE, HALF, TWO, ZERO, FACI, FAC1, FAC2
#else
      DOUBLE PRECISION BAODEN(*)
      DOUBLE PRECISION DHFAO(*),  B1DHFAO(*), B2DHFAO(*)
      DOUBLE PRECISION ZKABAO(*), B1KABAO(*), B2KABAO(*)
      DOUBLE PRECISION ONE, HALF, TWO, ZERO, FACI, FAC1, FAC2
#endif
      PARAMETER(HALF=0.5D0, ONE=1.0D0, ZERO=0.0D0, TWO=2.0D0)

      INTEGER KOFFGD, KOFFAB, KOFFAD, KOFFGB, ISYMA, ISYMB

*---------------------------------------------------------------------*
*     set FACI : if ICON = 2 multiply all contributions by 0.5
*---------------------------------------------------------------------*
      FACI = ONE
      IF (ICON .EQ. 2) FACI = HALF


*---------------------------------------------------------------------*
*     Add coulomb terms:  
*---------------------------------------------------------------------*

C     ------------------------------------------------------
C     2 D^HF_alp,bet (D^zeta_gambar,del + D^Zeta_gam,delbar)
C     ------------------------------------------------------
      IF (MULD2H(ISYMG,ISYMD) .EQ. ISYMQ) THEN
         KOFFGD = IAODIS(ISYMG,ISYMD) + NBAS(ISYMG)*(D - 1) + G
         FAC1   = TWO * ( B1KABAO(KOFFGD) + B2KABAO(KOFFGD) ) * FACI
         CALL DAXPY(N2BST(ISYM0),FAC1,DHFAO,1,BAODEN,1)
      END IF

C     ------------------------------------------------------
C     2 D^zeta_alp,bet (D^HF_gambar,del + D^HF_gam,delbar)
C     ------------------------------------------------------
      IF (MULD2H(ISYMG,ISYMD) .EQ. ISYMQ) THEN
         KOFFGD = IAODIS(ISYMG,ISYMD) + NBAS(ISYMG)*(D - 1) + G
         FAC1   = TWO * ( B1DHFAO(KOFFGD) + B2DHFAO(KOFFGD) ) * FACI
         CALL DAXPY(N2BST(ISYM0),FAC1,ZKABAO,1,BAODEN,1)
      END IF

*---------------------------------------------------------------------*
*     Add exchange terms:  
*---------------------------------------------------------------------*

      ISYMB = MULD2H(ISYMG,ISYMQ)
      ISYMA = MULD2H(ISYMD,ISYMQ)

      DO B = 1, NBAS(ISYMB)

         KOFFGB = IAODIS(ISYMG,ISYMB) + NBAS(ISYMG)*(B-1) + G
         KOFFAD = IAODIS(ISYMA,ISYMD) + NBAS(ISYMA)*(D-1) + 1
         KOFFAB = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(B-1) + 1

C        --------------------------------
C        - D^HF_alp,del D^zeta_gambar,bet 
C        --------------------------------
         FAC2 = -B1KABAO(KOFFGB) * FACI
         CALL DAXPY(NBAS(ISYMA),FAC2,DHFAO(KOFFAD),1,
     &                               BAODEN(KOFFAB),1)

C        --------------------------------
C        - D^HF_alp,delbar D^zeta_gam,bet 
C        --------------------------------
         FAC2 = -ZKABAO(KOFFGB) * FACI
         CALL DAXPY(NBAS(ISYMA),FAC2,B2DHFAO(KOFFAD),1,
     &                               BAODEN(KOFFAB),1)

C        --------------------------------
C        - D^zeta_alp,del D^HF_gambar,bet 
C        --------------------------------
         FAC2 = -B1DHFAO(KOFFGB) * FACI
         CALL DAXPY(NBAS(ISYMA),FAC2,ZKABAO(KOFFAD),1,
     &                               BAODEN(KOFFAB),1)

C        --------------------------------
C        - D^zeta_alp,delbar D^HF_gam,bet 
C        --------------------------------
         FAC2 = -DHFAO(KOFFGB) * FACI
         CALL DAXPY(NBAS(ISYMA),FAC2,B2KABAO(KOFFAD),1,
     &                               BAODEN(KOFFAB),1)

      END DO

      RETURN
      END
*=====================================================================*
