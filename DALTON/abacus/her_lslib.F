!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
C FILE: her_lslib.F

      SUBROUTINE DALTON_LSLIB_FCK(FMAT,DMAT,NDMAT,
     &                  ISYMDM,IFCTYP,IPRFCK,ICEDIF,IFTHRS,WORK,LWORK)
C
C     PURPOSE : Driver routine for direct calculation of the
C               two-electron part of the Fock matrices using LSLIB routines.
C               We assume the densities and fock matrices are full
C               squares and without symmetry reduction.
C
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT),
     &          ISYMDM(NDMAT), IFCTYP(NDMAT)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  DFTCOM: HFXFAC
C  INFORB: N2BASX,NBAST,?
C  GNRINF: PARCAL
C
#include "gnrinf.h"
#include "dftcom.h"
#include "inforb.h"
C defined parallel calculation types  
#include "iprtyp.h"
      IF (NDMAT .LE. 0) RETURN  
      CALL QENTER('DALTON_LSLIB_FCK')
#ifdef BUILD_LSLIB
      write(lupri,'(//A)') 'DALTON_LSLIB_FCK called :-)'
      call quit('sorry, LSLIB calls not programmed yet ...')
#else
      call quit('sorry, LSLIB not available in this Dalton build')
#endif
      CALL QEXIT('DALTON_LSLIB_FCK')
      RETURN
      END
! -- end of her_lslib.F --
