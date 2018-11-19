!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
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
