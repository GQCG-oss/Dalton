C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C  /* Deck drive_fit */
      SUBROUTINE DRIVE_FIT(DENS,ADENS,LDM,INDAO,INDAOS,INDNAO,IPRINT,
     &                     CFIT)
C*****************************************************************************
C    
C     DRIVE_FIT: Driver for density fitting in Dirac.
C        
C     Written by Luuk Visscher, november 2003
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      DIMENSION DENS(LDM,*),ADENS(LDM,*),CFIT(*)
      DIMENSION INDAO(*),INDAOS(*)
      DIMENSION INDNAO(2,*)
#include "nuclei.h"
#include "symmet.h"
C
#include "ibtfun.h"
C
C     ***************************************
C     ***** Sort the Density matrix      ****
C     ***************************************
C
C     The Density matrix is already in AO-basis but we need it
C     with all functions on one center consecutively.
C
      CALL REORMT (DENS,ADENS,LDM,LDM,INDAO,INDAO)
C
C     **************************************************
C     ***** Loop over symmetry unique center pairs  ****
C     **************************************************
C
      print*,"in driver of density fitting "
C     Loop over symmetry unique centers
      DO IB = 1, NUCIND
         MULB = ISTBNU(IB)
C        Loop over symmetry unique centers
         DO IA = 1, NUCIND
            MULA = ISTBNU(IA)
            MAB  = IBTOR (MULA,MULB)
C           Run over symmetry independent charge distributions
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MAB) .NE. 0) CYCLE
               ICENTA = NUCNUM(IA,1)
               ICENTB = NUCNUM(IB,ISYMOP+1)
C              The LL part of the density matrix
               INDA = INDNAO(1,ICENTA)+1
               INDB = INDNAO(1,ICENTB)+1
C              We need a pointer array that gives the start of the fit functions for a center
               IFTA = (ICENTA-1)*2 + 1
               IFTB = (ICENTB-1)*2 + 1
               print*,"fitting atom pair :",icenta,icentb
               CALL FitPairDensity (ICENTA,ICENTB,ADENS(INDA,INDB),LDM,
     &                              CFIT(IFTA), CFIT(IFTB), IERR)
C              The SS part of the density matrix
               INDA = INDNAO(2,ICENTA)+1
               INDB = INDNAO(2,ICENTB)+1
               CALL FitPairDensity (ICENTA,ICENTB,ADENS(INDA,INDB),LDM,
     &                              CFIT(IFTA), CFIT(IFTB), IERR)
            END DO  
         END DO 
      END DO 
C
      RETURN
      END
