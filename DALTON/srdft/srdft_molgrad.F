

C*****************************************************************************
      SUBROUTINE SRDFT_MOLGRAD_DRIVER(WORK, LWORK, IPRINT)
C*****************************************************************************
C
C     Jan. 2019: driver routine for srDFT contributions to molecular gradient
C
C*****************************************************************************
      implicit none
! dftinf.h -> ipnt_st, ipnt_end
#include "dftinf.h"
      integer, intent(in) :: iprint
      integer :: LWORK, ipnt, ipnt_st, ipnt_end
      real*8 :: WORK(LWORK)

      CALL QUIT('SRDFT_MOLGRAD not implemented yet')
      
      do ipnt = ipnt_st, ipnt_end
      end do

      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRDFT_MOLGRAD(rho, gamma_xy, laplace, tau)
C*****************************************************************************
C
C     
C
C*****************************************************************************
      implicit none
! dftcom.h -> do_metagga, dogga
#include "dftcom.h"
      real*8, intent(in) :: rho(4), gamma_xy(6), laplace, tau(4)
      logical :: do_metagga, dogga
      integer :: isym
      
      if (isym .eq. 0) then
         if (do_metagga) then
         else if (dogga) then
         else ! lda
         end if
      else
      end if
      RETURN
      END
