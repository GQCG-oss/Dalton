!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
      SUBROUTINE QFITINP(WORD)
#if defined (BUILD_QFITLIB)
      use qfit_input_readers, only: dalton_input
#endif
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      CHARACTER WORD*7
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
      SUBROUTINE QFITCALC(FDENS)
#if defined (BUILD_QFITLIB)
      use qfit
#endif
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxaqn.h"
#include "symmet.h"
#include "moldip.h"
#include "magone.h"
#include "orgcom.h"
C
      DIMENSION FDENS(N2BASX)
      integer :: icharg
C
#if defined (BUILD_QFITLIB)
      call molchr(icharg)
      call qfit_initialize(cord(1:3,1:nucdep), charge(1:nucdep),
     *                     icharg, dip0, cmxyz)
      call qfit_fit(FDENS)
      call qfitres
#endif
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

      allocate( charges( nucdep ) )
      call qfit_get_results( charges )

      write(lupri,9000)
      write(lupri,*)

      do m = 1, size(charges)
          write(lupri,'(a,a4,f12.6)') '@ ',namdep(m), charges(m)
      enddo
      write(lupri,*)
      write(lupri,'(a,f12.6)') '  Sum = total charge:', sum(charges)

      !call qfit_finalize
      deallocate( charges )
#endif
 9000 format(4x,'Charges have been fitted to the molecular ',
     *          'electrostatic potential.')
      END

      SUBROUTINE QFITGSCALC(CMO,WRK,KFRSAV,LFRSAV)
C
C     CALCULATE MEP CHARGE FITTING FOR GS
C
C     LOOSELY BASED ON CALC_MEP BY JACOB KONGSTED
C
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "priunit.h"
#include "dipole.h"
      PARAMETER ( D1 = 1.0D0 )
      DIMENSION WRK(LFRSAV)
      DIMENSION CMO(*)
      LOGICAL LOCDEB

      LOCDEB = .FALSE.

      CALL QFITINFO

      CALL QENTER('QFITGSCALC')

      KFREE  = KFRSAV
      LFREE  = LFRSAV

      CALL MEMGET('REAL',KDV   ,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENC ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENV ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENS ,NNBASX,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KDV), NNASHX)
      CALL DZERO(WRK(KDENC), N2BASX)
      CALL DZERO(WRK(KDENV), N2BASX)
      CALL DZERO(WRK(KDENS), NNBASX)

      CALL MEMGET('REAL',KXDIP ,NNBASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KYDIP ,NNBASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZDIP ,NNBASX,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KXDIP), NNBASX)
      CALL DZERO(WRK(KYDIP), NNBASX)
      CALL DZERO(WRK(KZDIP), NNBASX)

      IF (NASHT .EQ. 1) THEN
         WRK(KDV) = D1
      ELSE IF (HSROHF) THEN
         DO I = 1, NASHT
            II = I*(I+1)/2
            WRK(KDV+II-1) = D1
         END DO
      ENDIF

      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,WRK(KDV),WRK(KFREE),LFREE)
      CALL DAXPY(N2BASX,1.0D0,WRK(KDENV),1,WRK(KDENC),1)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DGEFSP(NBAST,WRK(KDENC),WRK(KDENV))
      CALL PKSYM1(WRK(KDENV),WRK(KDENS),NBAS,NSYM,1)

C     CALCULATE MOLECULAR DIPOLE
C     NUCLEAR CONTRIBUTION
      DUM = 0.0d0
      CALL DIPNUC(DUM,DUM,0,.FALSE.)
C
C     ELECTRONIC CONTRIBUTION
      LU = -1
      CALL GPOPEN(LU,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',0,
     &     .FALSE.)

      REWIND LU

      CALL MOLLAB('XDIPLEN ', LU, LUPRI)
      CALL READT(LU, NNBASX, WRK(KXDIP))
      CALL MOLLAB('YDIPLEN ', LU, LUPRI)
      CALL READT(LU, NNBASX, WRK(KYDIP))
      CALL MOLLAB('ZDIPLEN ', LU, LUPRI)
      CALL READT(LU, NNBASX, WRK(KZDIP))
      CALL GPCLOSE(LU, 'KEEP')

      DIPME(1) = -DDOT(NNBASX, WRK(KDENS), 1, WRK(KXDIP), 1)
      DIPME(2) = -DDOT(NNBASX, WRK(KDENS), 1, WRK(KYDIP), 1)
      DIPME(3) = -DDOT(NNBASX, WRK(KDENS), 1, WRK(KZDIP), 1)

      CALL DP0SUM

      IF(LOCDEB) THEN
          WRITE(LUPRI,*) 'DENSITY'
          CALL OUTPAK(WRK(KDENS), NBAST, 1, LUPRI)

          WRITE(LUPRI,*) 'XDIPLEN:'
          CALL OUTPAK(WRK(KXDIP), NBAST, 1, LUPRI)

          WRITE(LUPRI,*) 'YDIPLEN:'
          CALL OUTPAK(WRK(KYDIP), NBAST, 1, LUPRI)

          WRITE(LUPRI,*) 'ZDIPLEN:'
          CALL OUTPAK(WRK(KZDIP), NBAST, 1, LUPRI)

          WRITE(LUPRI,*) 'DIPOLE COMPONENTS'
          DO I=1,3
              WRITE(LUPRI,'(2F16.9)') DIPMN(I), DIPME(I)
          ENDDO
      ENDIF

C     CALCULATE FITTED CHARGES.
      CALL QFITCALC(WRK(KDENC))

      CALL MEMREL('QFITGSCALC',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('QFITGSCALC')

      END SUBROUTINE
