!
!...   Copyright (c) 2013 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2013 (2013), see http://daltonprogram.org"
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
#include "single.h"
C
      SUBROUTINE QFITINP(WORD)
#if defined (BUILD_QFITLIB)
      use qfit_input_readers, only: dalton_input
#endif
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      CHARACTER WORD*7
#include "abainf.h"
#if defined (BUILD_QFITLIB)
      call dalton_input(word, lucmd, lupri)
#endif
      RETURN
      END
C
      SUBROUTINE QFITINFO
#if defined (BUILD_QFITLIB)
      use qfit, only : qfit_print_info
      CALL HEADER('Charge fitting (QFITLIB) settings',-1) 
      call qfit_print_info
#endif
      RETURN
      END
C
      SUBROUTINE QFITCALC(WORK,LWORK)
#if defined (BUILD_QFITLIB)
      use qfit
#endif
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "moldip.h"
#include "magone.h"
#include "orgcom.h"
      real*8, dimension(lwork) :: work
C
      real*8, dimension(:), allocatable :: dens, dso, fso, dummat
      real*8, dimension(:), allocatable :: charges
      integer :: m, icharg
C
      allocate( dso( nnbasx ) )
      allocate( fso( nnbasx ) )
      allocate( dens( nnbasx ) )
      allocate( dummat( nnbasx ) )
      dso = 0.0d0
      fso = 0.0d0
      dens = 0.0d0
      dummat = 0.0d0
#if defined (BUILD_QFITLIB)
C 
C     OBTAIN THE AO DENSITY MATRIX
      CALL DSOFSO(dso,fso,WORK(1),LWORK,MPRINT,
     &            MNODC,MNODV)
      CALL DSYM1(dens,dummat,dso,fso,
     &           .FALSE.,NBAST,MPRINT)
C
C     MOLECULAR CHARGE
      call molchr(icharg)
C
C     INVOKE QFIT AND CALCULATE THE FITTED CHARGES
      call qfit_initialize(cord(1:3,1:natoms), charge(1:natoms),
     *                     icharg, dip0, diporg)
      call qfit_fit(dens)
#endif
C
      deallocate( dso )
      deallocate( fso )
      deallocate( dens )
      deallocate( dummat )
C
      END
C
C
      SUBROUTINE QFITRES
#if defined (BUILD_QFITLIB)
      use qfit
#endif
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "moldip.h"
#include "orgcom.h"
      real*8, dimension(:), allocatable :: charges

#if defined (BUILD_QFITLIB)
      CALL HEADER('Potential fitted charges (QFITLIB)',-1)

      allocate( charges( natoms ) )
      call qfit_get_results( charges )

      write(lupri,9000)
      write(lupri,*)

      do m = 1, natoms
          write(lupri,'(a,a4,f9.4)') '@ ',namdep(m), charges(m)
      enddo
      write(lupri,*)
      write(lupri,'(a,f9.4)') '  Sum = total charge:', sum(charges)

      call qfit_finalize
      deallocate( charges )
#endif
 9000 format(4x,'Charges have been fitted to the molecular ',
     *          'electrostatic potential.')
      END
