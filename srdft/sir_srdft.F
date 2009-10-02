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
C File: sir_srdft.F
C
C Special rouitnes for srDFT in sirius/
C Responsible author: Hans Joergen Aa. Jensen.
C
C  /* Deck srdft_dens */
      SUBROUTINE SRDFT_DENS(DAO,UDV,PV,CMO,WRK,KFREE,LFREE)

C      Hans Joergen Aa. Jensen, 2-Jun-2009
C
C      Purpose: calculate DCAO and DVAO matrics for srDFT
C               (done in a new routine for future flexibility
C                for e.g. also  on-top pair density)
C
#include <implicit.h>
C
#include <inforb.h>
      DIMENSION DAO(N2BASX,*), UDV(NASHT,NASHT), PV(*), CMO(*), WRK(*)
C
      CALL QENTER('SRDFT_DENS')
      IF (LFREE .LT. N2BASX) CALL QUIT('Insufficient WRK dim')
      CALL GTDMSO(UDV,CMO,DAO(1,1),DAO(1,2),WRK(KFREE))
      IF (NASHT .GT. 0) THEN
          CALL DAXPY(N2BASX,1.0D0,DAO(1,2),1,DAO(1,1),1)
      END IF
      CALL QEXIT ('SRDFT_DENS')
      RETURN
      END
