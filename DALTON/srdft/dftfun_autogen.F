C File: srdft/dftfun_autogen.F


C****************************************************************************
      pure subroutine Convert_d1_alpha_beta_to_C_S(d1Ea, d1Eb,
     &                  functional_type, d1E) 
C****************************************************************************
C  Implemented by E. R. Kjellgren
C
C  Subroutine to convert the derivatives in terms of alpha and beta
C  density to be formulated in terms of charge and spin density
C 
C  functional_type = 0 is LDA
C  functional_type = 1 is GGA
C  functional_type = 2 is MGGA
C
C            drho_x, dgamma_xx, dtau_x, dlaplace_x
C  d1Ex(idx)    1,      2,        3,        4
C****************************************************************************
      implicit none

      integer, intent(in) :: functional_type
      real*8, intent(in) :: d1Ea(4), d1Eb(4)
      real*8, intent(out) :: d1E(9)
      d1E(:) = 0.0d0

      d1E(1) = 0.5d0 * (d1Ea(1) + d1Eb(1))
      d1E(2) = 0.5d0 * (d1Ea(1) - d1Eb(1))
      if (functional_type .ge. 1) then ! GGA
         d1E(3) = 0.25d0 * (d1Ea(2) + d1Eb(2))
         d1E(4) = d1E(3)
         d1E(5) = 0.5d0 * (d1Ea(2) - d1Eb(2))
      end if
      if (functional_type .ge. 2) then ! MGGA
         d1E(6) = 0.5d0*(d1Ea(3) + d1Eb(3))
         d1E(7) = 0.5d0*(d1Ea(3) - d1Eb(3))
         d1E(8) = 0.5d0*(d1Ea(4) + d1Eb(4))
         d1E(9) = 0.5d0*(d1Ea(4) - d1Eb(4))
      end if
      end subroutine


C****************************************************************************
      pure subroutine Convert_d2_alpha_beta_to_C_S(d2Ea, d2Eb,
     &                  functional_type, d2E) 
C****************************************************************************
C  Implemented by E. R. Kjellgren
C
C  Subroutine to convert the derivatives in terms of alpha and beta
C  density to be formulated in terms of charge and spin density
C
C  functional_type = 0 is LDA
C  functional_type = 1 is GGA
C  functional_type = 2 is MGGA
C
C  d2Ex(idx) 
C              drho_x, dgamma_xx, dtau_x, dlaplace_x
C drho_x         1
C dgamma_xx      2         3
C dtau_x         4         5        6
C dlaplace_x     7         8        9        10
C****************************************************************************
      implicit none

      integer, intent(in) :: functional_type
      real*8, intent(in) :: d2Ea(10), d2Eb(10)
      real*8, intent(out) :: d2E(45)
      d2E(:) = 0.0d0

      d2E(1) = 0.25d0*(d2Ea(1) + d2Eb(1))
      d2E(2) = 0.25d0*(d2Ea(1) - d2Eb(1))
      d2E(3) = d2E(1)
      if (functional_type .ge. 1) then ! GGA
         d2E(4) = 0.125d0*(d2Ea(2) + d2Eb(2))
         d2E(5) = 0.125d0*(d2Ea(2) - d2Eb(2))
         d2E(6) = 0.0625d0*(d2Ea(3) + d2Eb(3))
         d2E(7) = d2E(4)
         d2E(8) = d2E(5)
         d2E(9) = d2E(6)
         d2E(10) = d2E(6)
         d2E(11) = 2.0d0*d2E(5)
         d2E(12) = 2.0d0*d2E(4)
         d2E(13) = 0.125d0*(d2Ea(3) - d2Eb(3))
         d2E(14) = d2E(13)
         d2E(15) = 4.0d0*d2E(6)
      end if
      if (functional_type .ge. 2) then !MGGA
         d2E(16) = 0.25d0*(d2Ea(4) + d2Eb(4))
         d2E(17) = 0.25d0*(d2Ea(4) - d2Eb(4))
         d2E(18) = 0.125d0*(d2Ea(5) + d2Eb(5))
         d2E(19) = d2E(18)
         d2E(20) = 2.0d0*d2E(24)
         d2E(21) = 0.25d0*(d2Ea(6) + d2Eb(6))
         d2E(22) = 0.25d0*(d2Ea(4) - d2Eb(4))
         d2E(23) = d2E(16)
         d2E(24) = 0.125d0*(d2Ea(5) - d2Eb(5))
         d2E(25) = d2E(24)
         d2E(26) = 2.0d0*d2E(18)
         d2E(27) = 0.25d0*(d2Ea(6) - d2Eb(6))
         d2E(28) = d2E(21)
         d2E(29) = 0.25d0*(d2Ea(7) + d2Eb(7))
         d2E(30) = 0.25d0*(d2Ea(7) - d2Eb(7))
         d2E(31) = 0.125d0*(d2Ea(8) + d2Eb(8))
         d2E(32) = d2E(31)
         d2E(33) = 2.0d0*d2E(39)
         d2E(34) = 0.25d0*(d2Ea(9) + d2Eb(9))
         d2E(35) = 0.25d0*(d2Ea(9) - d2Eb(9))
         d2E(36) = 0.25d0*(d2Ea(10) + d2Eb(10))
         d2E(37) = 0.25d0*(d2Ea(7) - d2Eb(7))
         d2E(38) = d2E(29)
         d2E(39) = 0.125d0*(d2Ea(8) - d2Eb(8))
         d2E(40) = d2E(39)
         d2E(41) = 2.0d0*d2E(31)
         d2E(42) = 0.25d0*(d2Ea(9) - d2Eb(9))
         d2E(43) = d2E(34)
         d2E(44) = 0.25d0*(d2Ea(10) - d2Eb(10))
         d2E(45) = d2E(36)
      end if
      end subroutine


C****************************************************************************
      subroutine BOTH_LDAX_ERF(rho, diff_order, mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange LDA functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "pi.h" ! pi and pi2
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(4), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: alpha_eq_beta
      real*8 :: A, kF, Ea, Eb, d1Ea(4), d2Ea(10), d1Eb(4), d2Eb(10),
     &          rho_a, rho_b

      E = 0.0d0
      Ea = 0.0d0
      d1E(:) = 0.0d0
      d1Ea(:) = 0.0d0
      d2E(:) = 0.0d0
      d2Ea(:) = 0.0d0
      rho_a = rho(3)
      rho_b = rho(4)
      
      if ((rho_a + rho_b) .le. 1.d-10) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (rho_a .eq. rho_b) then
         alpha_eq_beta = .true.
      else
         alpha_eq_beta = .false.
         Eb = 0.0d0
         d1Eb(:) = 0.0d0
         d2Eb(:) = 0.0d0
      end if

      kF = (rho_a + rho_b)**(1.0d0/3.0d0)*((3.0d0*pi2)**(1.0d0/3.0d0))
      A = mu / (2.0d0*kF)

      if (diff_order .eq. 0) then
         write(lupri,*)'ESRX_LDA_ERF_case_x:'
         write(lupri,*)'WARNING: Code not covered in ctest'
         if (A .lt. 1.d-9) then
            call ESRX_LDA_ERF_case_1(rho_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_LDA_ERF_case_1(rho_b, mu, Eb)
            end if
         elseif (A .le. 1.d+2) then
            call ESRX_LDA_ERF_case_2(rho_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_LDA_ERF_case_2(rho_b, mu, Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call ESRX_LDA_ERF_case_3(rho_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_LDA_ERF_case_3(rho_b, mu, Eb)
            end if
         end if

      elseif (diff_order .eq. 1) then
         if (A .lt. 1.d-9) then
            write(lupri,*)'D1ESRX_LDA_ERF_case_1:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_LDA_ERF_case_1(rho_a, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_LDA_ERF_case_1(rho_b, mu, Eb, d1Eb)
            end if
         elseif (A .le. 1.d+2) then
            call D1ESRX_LDA_ERF_case_2(rho_a, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_LDA_ERF_case_2(rho_b, mu, Eb, d1Eb)
            end if
         elseif (A .lt. 1.d+9) then
            write(lupri,*)'D1ESRX_LDA_ERF_case_3:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_LDA_ERF_case_3(rho_b, mu, Eb, d1Eb)
            end if
         end if

      elseif (diff_order .eq. 2) then
         if (A .lt. 1.d-9) then
            write(lupri,*)'D2ESRX_LDA_ERF_case_1:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_LDA_ERF_case_1(rho_a, mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_LDA_ERF_case_1(rho_b, mu, Eb, d1Eb, d2Eb)
            end if
         elseif (A .le. 1.d+2) then
            call D2ESRX_LDA_ERF_case_2(rho_a, mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_LDA_ERF_case_2(rho_b, mu, Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 1.d+9) then
            write(lupri,*)'D2ESRX_LDA_ERF_case_3:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_LDA_ERF_case_3(rho_b, mu, Eb, d1Eb, d2Eb)
            end if
         end if

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if

C     Do the convertion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting realtion
      E = 0.5d0*Ea + 0.5d0*Eb
      if (diff_order .ge. 1) then
         call Convert_d1_alpha_beta_to_C_S(0.5d0*d1Ea, 0.5d0*d1Eb, 0,
     &                                     d1E)
      end if
      if (diff_order .ge. 2) then
         call Convert_d2_alpha_beta_to_C_S(0.5d0*d2Ea, 0.5d0*d2Eb, 0,
     &                                     d2E)
      end if
      end subroutine


C****************************************************************************
      subroutine BOTH_PBEX_ERF(rho, gamma_vec, diff_order, mu, E, d1E, 
     &                         d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange PBE functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           gamma_vec, is called gamma_vec, becuase "gamma" is an
C              intrinsic function and therefore a bad name.
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "pi.h" ! pi and pi2
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(4), gamma_vec(6), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: alpha_eq_beta
      real*8 :: A, kF, Ea, Eb, d1Ea(4), d2Ea(10), d1Eb(4), d2Eb(10),
     &          rho_a, rho_b, gamma_aa, gamma_bb

      E = 0.0d0
      Ea = 0.0d0
      d1E(:) = 0.0d0
      d1Ea(:) = 0.0d0
      d2E(:) = 0.0d0
      d2Ea(:) = 0.0d0
      rho_a = rho(3)
      rho_b = rho(4)
      gamma_aa = gamma_vec(4)
      gamma_bb = gamma_vec(5)
      
      if (((rho_a + rho_b) .le. 1.d-10) 
     &     .and. ((gamma_aa + gamma_bb) .lt. 1.0d-10)) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if ((rho_a .eq. rho_b) .and. (gamma_aa .eq. gamma_bb)) then
         alpha_eq_beta = .true.
      else
         alpha_eq_beta = .false.
         Eb = 0.0d0
         d1Eb(:) = 0.0d0
         d2Eb(:) = 0.0d0
      end if

      kF = (rho_a + rho_b)**(1.0d0/3.0d0)*((3.0d0*pi2)**(1.0d0/3.0d0))
      A = mu / (2.0d0*kF)

      if (diff_order .eq. 0) then
         write(lupri,*)'ESRX_PBE_ERF_case_x:'
         write(lupri,*)'WARNING: Code not covered in ctest'
         if (A .lt. 1.d-9) then
            call ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_ERF_case_1(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_ERF_case_2_1(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_ERF_case_2_3(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_ERF_case_2_2(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_ERF_case_3(rho_b, gamma_bb, mu, Eb)
            end if
         end if

      elseif (diff_order .eq. 1) then
         if (A .lt. 1.d-9) then
            write(lupri,*)'D1ESRX_PBE_ERF_case_1:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_ERF_case_1(rho_b, gamma_bb, mu, Eb, d1Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D1ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_ERF_case_2_1(rho_b, gamma_bb, mu, Eb,
     &                                      d1Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D1ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_ERF_case_2_3(rho_b, gamma_bb, mu, Eb,
     &                                      d1Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            write(lupri,*)'D1ESRX_PBE_ERF_case_2_2:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_ERF_case_2_2(rho_b, gamma_bb, mu, Eb,
     &                                      d1Eb)
            end if
         elseif (A .lt. 1.d+9) then
            write(lupri,*)'D1ESRX_PBE_ERF_case_3:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_ERF_case_3(rho_b, gamma_bb, mu, Eb, d1Eb)
            end if
         end if

      elseif (diff_order .eq. 2) then
         if (A .lt. 1.d-9) then
            write(lupri,*)'D2ESRX_PBE_ERF_case_1:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_ERF_case_1(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D2ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_ERF_case_2_1(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D2ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_ERF_case_2_3(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            write(lupri,*)'D2ESRX_PBE_ERF_case_2_2:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_ERF_case_2_2(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .lt. 1.d+9) then
            write(lupri,*)'D2ESRX_PBE_ERF_case_3:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_ERF_case_3(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         end if

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if

C     Do the convertion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting realtion
      E = 0.5d0*Ea + 0.5d0*Eb
      if (diff_order .ge. 1) then
         call Convert_d1_alpha_beta_to_C_S(0.5d0*d1Ea, 0.5d0*d1Eb, 1,
     &                                     d1E)
      end if
      if (diff_order .ge. 2) then
         call Convert_d2_alpha_beta_to_C_S(0.5d0*d2Ea, 0.5d0*d2Eb, 1,
     &                                     d2E)
      end if
      end subroutine


C****************************************************************************
      subroutine BOTH_TPSSX_ERF(rho, gamma_vec, tau, diff_order, mu, E,
     &                          d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange TPSS functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           gamma_vec, is called gamma_vec, becuase "gamma" is an
C              intrinsic function and therefore a bad name.
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "pi.h" ! pi and pi2
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(4), gamma_vec(6), tau(4), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: alpha_eq_beta
      real*8 :: A, kF, Ea, Eb, d1Ea(4), d2Ea(10), d1Eb(4), d2Eb(10),
     &          rho_a, rho_b, gamma_aa, gamma_bb, tau_a, tau_b

      E = 0.0d0
      Ea = 0.0d0
      d1E(:) = 0.0d0
      d1Ea(:) = 0.0d0
      d2E(:) = 0.0d0
      d2Ea(:) = 0.0d0
      rho_a = rho(3)
      rho_b = rho(4)
      gamma_aa = gamma_vec(4)
      gamma_bb = gamma_vec(5)
      tau_a = tau(3)
      tau_b = tau(4)
      
      if (((rho_a + rho_b) .le. 1.d-10) 
     &     .and. ((gamma_aa + gamma_bb) .lt. 1.0d-10)
     &     .and. ((tau_a + tau_b) .lt. 1.0d-10)) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if ((rho_a .eq. rho_b) .and. (gamma_aa .eq. gamma_bb)
     &    .and. (tau_a .eq. tau_b)) then
         alpha_eq_beta = .true.
      else
         alpha_eq_beta = .false.
         Eb = 0.0d0
         d1Eb(:) = 0.0d0
         d2Eb(:) = 0.0d0
      end if

      kF = (rho_a + rho_b)**(1.0d0/3.0d0)*((3.0d0*pi2)**(1.0d0/3.0d0))
      A = mu / (2.0d0*kF)

      if (diff_order .eq. 0) then
         write(lupri,*)'ESRX_TPSS_ERF_case_x:'
         write(lupri,*)'WARNING: Code not covered in ctest'
         if (A .lt. 1.d-9) then
            call ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_ERF_case_1(rho_b, gamma_bb, tau_b, mu, Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_ERF_case_2_1(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_ERF_case_2_3(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_ERF_case_2_2(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_ERF_case_3(rho_b, gamma_bb, tau_b, mu, Eb)
            end if
         end if

      elseif (diff_order .eq. 1) then
         if (A .lt. 1.d-9) then
            write(lupri,*)'D1ESRX_TPSS_ERF_case_1:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     &                                  Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_ERF_case_1(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb, d1Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D1ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_ERF_case_2_1(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb, d1Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D1ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_ERF_case_2_3(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb, d1Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            write(lupri,*)'D1ESRX_TPSS_ERF_case_2_2:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_ERF_case_2_2(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb, d1Eb)
            end if
         elseif (A .lt. 1.d+9) then
            write(lupri,*)'D1ESRX_TPSS_ERF_case_3:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D1ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     &                                  Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_ERF_case_3(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb, d1Eb)
            end if
         end if

      elseif (diff_order .eq. 2) then
         if (A .lt. 1.d-9) then
            write(lupri,*)'D2ESRX_TPSS_ERF_case_1:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     &                                  Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_ERF_case_1(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D2ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_ERF_case_2_1(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D2ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_ERF_case_2_3(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb, d1Eb, d2Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            write(lupri,*)'D2ESRX_TPSS_ERF_case_2_2:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_ERF_case_2_2(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 1.d+9) then
            write(lupri,*)'D2ESRX_TPSS_ERF_case_3:'
            write(lupri,*)'WARNING: Code not covered in ctest'
            call D2ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     &                                  Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_ERF_case_3(rho_b, gamma_bb, tau_b, mu,
     &                                     Eb, d1Eb, d2Eb)
            end if
         end if

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if

C     Do the convertion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting realtion
      E = 0.5d0*Ea + 0.5d0*Eb
      if (diff_order .ge. 1) then
         call Convert_d1_alpha_beta_to_C_S(0.5d0*d1Ea, 0.5d0*d1Eb, 2,
     &                                     d1E)
      end if
      if (diff_order .ge. 2) then
         call Convert_d2_alpha_beta_to_C_S(0.5d0*d2Ea, 0.5d0*d2Eb, 2,
     &                                     d2E)
      end if
      end subroutine


C****************************************************************************
      subroutine SPIN_PW92C_ERF(rho, diff_order, mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized PW92 functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: rho_c, rho_s
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)

      if (rho_c .le. 1.d-10) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (diff_order .eq. 0) then
         write(lupri,*)'ESRC_PW92_ERF:'
         write(lupri,*)'WARNING: Code not covered in ctest'
         call ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E)

      elseif (diff_order .eq. 1) then
         call D1ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E, d2E)

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if
      end subroutine


C****************************************************************************
      subroutine PW92C_ERF(rho, diff_order, mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized PW92 functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: rho_c
      
      write(lupri,*)'PW92C_ERF:'
      write(lupri,*)'WARNING: Code not covered in ctest'

      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)

      if (rho_c .le. 1.d-10) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (diff_order .eq. 0) then
         call ESRC_PW92_ERF(rho_c, mu, E)

      elseif (diff_order .eq. 1) then
         call D1ESRC_PW92_ERF(rho_c, mu, E, d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_PW92_ERF(rho_c, mu, E, d1E, d2E)

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if
      end subroutine


C****************************************************************************
      subroutine SPIN_PBEC_GWS_ERF(rho, gamma_xx, diff_order, mu, E,
     &                             d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized PW92 functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(2), gamma_xx(3), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: rho_c, rho_s, gamma_cc
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)
      gamma_cc = gamma_xx(1)

      if (rho_c .le. 1.d-10) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (diff_order .eq. 0) then
         write(lupri,*)'ESRC_SPIN_PBE_GWS_ERF:'
         write(lupri,*)'WARNING: Code not covered in ctest'
         call ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu, E)

      elseif (diff_order .eq. 1) then
         call D1ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu, E,
     &                                d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu, E,
     &                                d1E, d2E)

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if
      end subroutine


C****************************************************************************
      subroutine PBEC_GWS_ERF(rho, gamma_xx, diff_order, mu, E,
     &                             d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized PW92 functional
C
C   Input : rho_a, alpha-spin density
C           rho_b, beta-spin density
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h"
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(2), gamma_xx(3), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: rho_c, gamma_cc
      
      write(lupri,*)'PBEC_GWS_ERF:'
      write(lupri,*)'WARNING: Code not covered in ctest'

      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      gamma_cc = gamma_xx(1)

      if (rho_c .le. 1.d-10) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (diff_order .eq. 0) then
         call ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E)

      elseif (diff_order .eq. 1) then
         call D1ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E,d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E, d2E)

      else
         call quit('Requested to high order derivative: ',diff_order)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 09, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44
      E = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/rho_c
      x3 = x2**0.333333333333333d0
      x4 = 1d0/x1
      x5 = x2**0.166666666666667d0
      x6 = sqrt(x2)
      x7 = x2**0.666666666666667d0
      x8 = (0.132568899905202d0*x3 + 1.0d0)*log(4.934802200544679d0*x4/(
     & 2.2255694211506878d0*x3 + 5.9825504357710795d0*x5 + 0.80042863499
     &936356d0*x6 + 0.18970043257475591d0*x7) + 1.0d0)
      x9 = 0.202642367284676d0*x1*x8
      x10 = rho_s**4/rho_c**4
      x11 = rho_s*x2
      x12 = x11 + 1.0d0
      x13 = -x11 + 1.0d0
      x14 = 1.92366105093154d0*x12**1.33333333333333d0 + 1.9236610509315
     &4d0*x13** 1.33333333333333d0 - 3.84732210186307d0
      x15 = 0.0197517897025652d0*x14*(-x10 + 1.0d0)*(0.0690139921125583d
     &0*x3 + 1.0d0 )*log(1.0d0 + 29.608574643216677d0/(2.247591863577616
     &4d0*x3 + 8.1574147034876408d0*x5 + 0.43009724712766434d0*x6 + 0.19
     &115125951273382d0*x7))
      x16 = 0.0506605918211689d0*(0.127469618870009d0*x3 + 1.0d0)*log( 9
     &.869604401089358d0*x4/(3.8447462374472119d0*x3 + 11.12037486309468
     &5d0*x5 + 1.6447337755676092d0*x6 + 0.24058712912881922d0*x7) + 1.0
     &d0)
      x17 = 0.101321183642338d0*x8
      x18 = 2.0d0*x1*x10*x14
      x19 = x18*(x16 - x17)
      x20 = mu**2
      x21 = x2**1.66666666666667d0
      x22 = 1.0d0 - rho_s**2/rho_c**2
      x23 = x2**1.0d0
      x24 = x2**1.33333333333333d0
      x25 = (-0.00304892328706429d0*x23 + 0.000275310793616358d0*x24 + 0
     &.0128480790170175d0*x3 + 0.0315297404574782d0*x7 + 1.0d0)*exp( -0.
     &46675853320810856d0*x3)
      x26 = x22*x25
      x27 = x12**0.666666666666667d0 + x13**0.666666666666667d0
      x28 = x27**3
      x29 = x20*x3/x27**2
      x30 = 9.20897059521559d0*mu*x5/x27 + 1.0d0
      x31 = mu**3
      x32 = 1d0/x7
      x33 = 1d0/x12
      x34 = x33**0.666666666666667d0
      x35 = x3*x33**0.333333333333333d0
      x36 = 0.0244354823084343d0*x7
      x37 = 8.61370200056457d0*x12**2*(-0.0177069812999059d0*x35 + 1.0d0
     &)/(x34*(x34* x36 + 0.337569861992026d0*x35 + 1.0d0))
      x38 = 1d0/x13
      x39 = x38**0.666666666666667d0
      x40 = x3*x38**0.333333333333333d0
      x41 = 8.61370200056457d0*x13**2*(-0.0177069812999059d0*x40 + 1.0d0
     &)/(x39*(x36* x39 + 0.337569861992026d0*x40 + 1.0d0))
      x42 = x15 - x9
      x43 = -x19 + x42
      x44 = x22*(x25 - 1.0d0)
      E = -rho_c*(-x15 + x19 + x9 + (0.00316102962473761d0*mu**8*x2** 2.
     &66666666666667d0*x43 - mu**6*(-0.0533250677421794d0*x2**2.0d0*( x1
     &8*(-x16 + x17) + x42) + 0.00248266686171384d0*x21*x44) - 0.0026411
     &7541143784d0*mu**5*x21*x26 + mu**4*( -0.00020688890514282d0*x2**1.
     &0d0*(-8.61370200056457d0*x12** 2.66666666666667d0 - 8.613702000564
     &57d0*x13**2.66666666666667d0 + 23.3866625383245d0*x22*(-0.24069599
     &0468967d0*x3 + 0.26014827853397d0*x7)*exp(-0.33933171852197186d0*x
     &3) + x37 + x41 ) - 0.0418814102677522d0*x23*x44 + 0.33733820787238
     &7d0*x24*x43) - 0.0253302959105844d0*x1*x28*log((-2.481401963597600
     &3d0*x29*( -3.4485100000000002d0 + 1.227722850784289d0/(x0 - 1.0d0)
     &) + x30 + 15.312568193831799d0*x31*x6/x28)/(8.5571394854859619d0*x
     &29 + x30 )) - x31*(0.000132058770571892d0*x21*(37.6991118430775d0*
     &x22*( -3.07073492995203d0*x3 + 0.384834731559127d0*x7)*exp( -0.192
     &30865217881402d0*x3)/x23 + x32*x37 + x32*x41) + 0.0445553741830535
     &d0*x23*x26))/(0.237113969739022d0*x20*x7 + 1.0d0)**4)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 09, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/rho_c
      x3 = x2**0.333333333333333d0
      x4 = 1d0/x1
      x5 = x2**0.166666666666667d0
      x6 = sqrt(x2)
      x7 = x2**0.666666666666667d0
      x8 = 2.22556942115069d0*x3 + 5.98255043577108d0*x5 + 0.80042863499
     &9364d0*x6 + 0.189700432574756d0*x7
      x9 = 4.93480220054468d0*x4/x8 + 1.0d0
      x10 = log(x9)
      x11 = x10*(0.132568899905202d0*x3 + 1.0d0)
      x12 = 0.202642367284676d0*x1*x11
      x13 = rho_s*x2
      x14 = x13 + 1.0d0
      x15 = -x13 + 1.0d0
      x16 = 1.92366105093154d0*x14**1.33333333333333d0 + 1.9236610509315
     &4d0*x15** 1.33333333333333d0 - 3.84732210186307d0
      x17 = rho_s**4/rho_c**4
      x18 = -x17 + 1.0d0
      x19 = 2.24759186357762d0*x3 + 8.15741470348764d0*x5 + 0.4300972471
     &27664d0*x6 + 0.191151259512734d0*x7
      x20 = 1.0d0 + 29.6085746432167d0/x19
      x21 = log(x20)
      x22 = x18*x21
      x23 = x16*x22
      x24 = 0.0197517897025652d0*x23*(0.0690139921125583d0*x3 + 1.0d0)
      x25 = 3.84474623744721d0*x3 + 11.1203748630947d0*x5 + 1.6447337755
     &6761d0*x6 + 0.240587129128819d0*x7
      x26 = 1.0d0 + 9.86960440108936d0*x4/x25
      x27 = log(x26)
      x28 = 0.0506605918211689d0*x27*(0.127469618870009d0*x3 + 1.0d0)
      x29 = 0.101321183642338d0*x11
      x30 = x28 - x29
      x31 = x16*x17
      x32 = x1*x31
      x33 = 2.0d0*x32
      x34 = x30*x33
      x35 = mu**2
      x36 = x35*x7
      x37 = 0.237113969739022d0*x36 + 1.0d0
      x38 = x37**(-4)
      x39 = rho_s**2/rho_c**2
      x40 = -x39 + 1.0d0
      x41 = exp(-0.46675853320810856d0*x3)
      x42 = x2**1.0d0
      x43 = -0.00304892328706429d0*x42
      x44 = x2**1.33333333333333d0
      x45 = 0.0128480790170175d0*x3 + x43 + 0.000275310793616358d0*x44 +
     & 0.0315297404574782d0*x7 + 1.0d0
      x46 = x41*x45
      x47 = x40*x46
      x48 = x2**1.66666666666667d0
      x49 = mu**5
      x50 = x48*x49
      x51 = x47*x50
      x52 = x14**0.666666666666667d0 + x15**0.666666666666667d0
      x53 = x52**3
      x54 = x52**2
      x55 = 1d0/x54
      x56 = x35*x55
      x57 = x3*x56
      x58 = 1d0/x52
      x59 = mu*x58
      x60 = 9.20897059521559d0*x5*x59 + 1.0d0
      x61 = mu**3
      x62 = x6/x53
      x63 = -3.44851d0 + 1.22772285078429d0/(x0 - 1.0d0)
      x64 = -2.4814019635976d0*x57*x63 + x60 + 15.3125681938318d0*x61*x6
     &2
      x65 = x64/(8.55713948548596d0*x57 + x60)
      x66 = x1*log(x65)
      x67 = x40*x42
      x68 = 0.0445553741830535d0*x67
      x69 = x46*x68
      x70 = -3.07073492995203d0*x3 + 0.384834731559127d0*x7
      x71 = exp(-0.19230865217881402d0*x3)
      x72 = x40*x70*x71
      x73 = 1d0/x42
      x74 = 37.6991118430775d0*x73
      x75 = x72*x74
      x76 = 1d0/x7
      x77 = x14**2
      x78 = 1d0/x14
      x79 = x78**0.333333333333333d0
      x80 = x3*x79
      x81 = x78**0.666666666666667d0
      x82 = 0.0244354823084343d0*x7
      x83 = 0.337569861992026d0*x80 + x81*x82 + 1.0d0
      x84 = 1d0/x83
      x85 = x77*x84
      x86 = (-0.0177069812999059d0*x80 + 1.0d0)/x81
      x87 = 8.61370200056457d0*x86
      x88 = x85*x87
      x89 = x15**2
      x90 = 1d0/x15
      x91 = x90**0.333333333333333d0
      x92 = x3*x91
      x93 = x90**0.666666666666667d0
      x94 = x82*x93 + 0.337569861992026d0*x92 + 1.0d0
      x95 = 1d0/x94
      x96 = x89*x95
      x97 = (-0.0177069812999059d0*x92 + 1.0d0)/x93
      x98 = 8.61370200056457d0*x97
      x99 = x96*x98
      x100 = x48*(x75 + x76*x88 + x76*x99)
      x101 = -x12 + x24
      x102 = x101 - x34
      x103 = mu**8
      x104 = x2**2.66666666666667d0
      x105 = 0.00316102962473761d0*x103*x104
      x106 = mu**6
      x107 = x46 - 1.0d0
      x108 = x107*x40
      x109 = x108*x48
      x110 = x2**2.0d0
      x111 = -x28 + x29
      x112 = x101 + x111*x33
      x113 = x110*x112
      x114 = mu**4
      x115 = x108*x42
      x116 = x2**1.0d0
      x117 = x14**2.66666666666667d0
      x118 = 8.61370200056457d0*x117
      x119 = x15**2.66666666666667d0
      x120 = 8.61370200056457d0*x119
      x121 = exp(-0.33933171852197186d0*x3)
      x122 = x121*(-0.240695990468967d0*x3 + 0.26014827853397d0*x7)
      x123 = x122*x40
      x124 = 23.3866625383245d0*x123
      x125 = 0.337338207872387d0*x44
      x126 = x102*x105 - x106*(0.00248266686171384d0*x109 - 0.0533250677
     &421794d0*x113 ) + x114*(x102*x125 - 0.0418814102677522d0*x115 - 0.
     &00020688890514282d0*x116*(-x118 - x120 + x124 + x88 + x99)) - 0.00
     &264117541143784d0*x51 - 0.0253302959105844d0*x53*x66 - x61*( 0.000
     &132058770571892d0*x100 + x69)
      x127 = x126*x38
      x128 = x39*x46
      x129 = x41*(0.0042826930056725d0*x3 + x43 + 0.000367081058155144d0
     &*x44 + 0.0210198269716521d0*x7)
      x130 = 0.00528235082287567d0*x50
      x131 = x14**(-0.333333333333333d0)
      x132 = x15**(-0.333333333333333d0)
      x133 = x131 - x132
      x134 = x133*x66
      x135 = x3*x59
      x136 = x13*x133
      x137 = mu*x136*x3*x55
      x138 = 6.13931373014372d0*x5
      x139 = x136*x138*x58 - 1.53482843253593d0*x5
      x140 = x56*x6
      x141 = x135*x63
      x142 = 0.0506605918211689d0*x54
      x143 = mu*x1/x64
      x144 = 16.0d0*x32
      x145 = -x14**0.333333333333333d0 + x15**0.333333333333333d0
      x146 = x1*x145
      x147 = rho_s**5*x146/rho_c**5
      x148 = 10.2595256049682d0*x147
      x149 = x1*x3
      x150 = x10*x149
      x151 = 0.0179093838034102d0*x150
      x152 = (0.265137799810404d0*x3 + 2.0d0)*(0.741856473716896d0*x3 + 
     &0.99709173929518d0*x5 + 0.400214317499682d0*x6 + 0.126466955049837
     &d0*x7)/(x8**2*x9)
      x153 = x31*(0.00861024844156464d0*x149*x27 - x151 + x152 - (0.2549
     &39237740017d0 *x3 + 2.0d0)*(1.28158207914907d0*x3 + 1.853395810515
     &78d0*x5 + 0.822366887783805d0*x6 + 0.160391419419213d0*x7)/(x25**2
     &*x26))
      x154 = x23*x3
      x155 = 0.00136314985874174d0*x3 + 0.0197517897025652d0
      x156 = x155*x21
      x157 = x156*x31
      x158 = x145*x155*x22
      x159 = x155*x16*x18*(0.749197287859205d0*x3 + 1.35956911724794d0*x
     &5 + 0.215048623563832d0*x6 + 0.127434173008489d0*x7)/(x19**2*x20)
      x160 = 5.1297628024841d0*x13*x158 + x151 - x152 + x153 - 0.0009087
     &66572494496d0 *x154 + 8.0d0*x157 + 59.2171492864334d0*x159
      x161 = 0.0891107483661069d0*x42
      x162 = 75.398223686155d0*x70*x71*x73
      x163 = x13*x78
      x164 = 1d0/x3
      x165 = 0.0508408867489862d0*x164
      x166 = x165/x79
      x167 = x166*x85*(-x163 + 1.0d0)
      x168 = 5.74246800037638d0*x76
      x169 = x168*x85*x86
      x170 = x13*x90
      x171 = x165/x91
      x172 = x171*x96*(x170 + 1.0d0)
      x173 = x14*x84
      x174 = 22.9698720015055d0*x76
      x175 = x173*x174*x86
      x176 = x13*x175
      x177 = x168*x96*x97
      x178 = x15*x95
      x179 = x174*x178*x97
      x180 = x13*x179
      x181 = 0.112523287330675d0*x80
      x182 = 0.0162903215389562d0*x7
      x183 = x182*x81
      x184 = x181 + x183
      x185 = x76*x87/x83**2
      x186 = x185*x77*(-x163*x181 - x163*x183 + x184)
      x187 = 0.112523287330675d0*x92
      x188 = x182*x93
      x189 = x187 + x188
      x190 = x76*x98/x94**2
      x191 = x190*x89*(x170*x187 + x170*x188 + x189)
      x192 = 0.000132058770571892d0*x48
      x193 = x107*x39
      x194 = 0.00993066744685536d0*x48
      x195 = x41*(0.0777930888680181d0*x3*x45 - 0.00214134650283625d0*x3
     & + 0.00152446164353214d0*x42 - 0.000183540529077572d0*x44 - 0.0105
     &099134858261d0*x7)
      x196 = -x111*x144 + x111*x148 + x160
      x197 = 0.0533250677421794d0*x110
      x198 = x15**1.66666666666667d0
      x199 = x14**1.66666666666667d0
      x200 = 22.9698720015055d0*x76
      x201 = 46.773325076649d0*x122*x76
      x202 = 5.1297628024841d0*x30
      x203 = 2.56488140124205d0*x158
      x204 = x13*x46
      x205 = rho_s**3*x16/rho_c**3
      x206 = -8.0d0*x1*x205*x30 + x146*x17*x202 - 4.0d0*x156*x205 - x203
     &
      x207 = x107*x13
      x208 = x14*x184*x185 - x15*x189*x190 + x166*x173 - x171*x178 + x17
     &5 - x179
      E = -rho_c*(x12 + x127 - x24 + x34)
      d1E(1) = x101 - 0.632303919304059d0*x126*x36/x37**5 - x127 + x13*x
     &203 - x147*x202 + 0.0089546919017051d0*x150 - 0.5d0*x152 + 0.5d0*x
     &153 - 0.000454383286247248d0*x154 + 4.0d0*x157 + 29.6085746432167d
     &0* x159 + 6.0d0*x30*x32 - 0.5d0*x38*(-0.0168588246652672d0*x102*x1
     &03 *x104 + x105*(x144*x30 - x148*x30 + x160) - x106*( -0.008275556
     &2057128d0*x109 + 0.213300270968718d0*x113 + x193*x194 + x194*x195*
     &x40 - x196*x197) - x114*(0.899568554326366d0*x112*x44 - 0.08376282
     &05355044d0*x115 + 0.0006896296838094d0*x116*(x118 + x120 - x124 - 
     &x88 - x99) - x125*x196 + 0.167525641071009d0*x193* x42 + 0.1675256
     &41071009d0*x195*x67 - 0.00041377781028564d0*x48*( -23.386662538324
     &5d0*x121*x40*x76*(0.0802319968229891d0*x3 - 0.173432185689313d0*x7
     &) - 2.64527879654103d0*x123*x164 - 15.591108358883d0*x123*x76 + x1
     &3*x200*(x198 - x199) - x167 - x169 - x172 + x176 - x177 - x180 - x
     &186 - x191 - x201*x39 + 5.74246800037638d0*x76*(x117 + x119))) - 0
     &.0105647016457513d0* x128*x50 + x129*x130*x40 + 0.101321183642338d
     &0*x13*x134*x54 - x142*x143*(30.6251363876636d0*x136*x35*x62 - 3.30
     &853595146347d0* x137*x63 + x139 - 7.6562840969159d0*x140 + 0.82713
     &3987865867d0* x141 - x65*(-2.85237982849532d0*x135 + 11.4095193139
     &813d0*x137 + x139)) - 0.000821860773992031d0*x2**2.0d0*x47*x49 + 0
     &.00880391803812612d0*x51 + 2.0d0*x61*(0.000220097950953153d0* x100
     & - x128*x161 + x129*x68 - x192*(x162*x39 + x167 + x169 + x172 - x1
     &76 + x177 + x180 + x186 + x191 + 2.41662179562687d0*x2**( -0.66666
     &6666666667d0)*x72 + x40*x71*x74*(1.02357830998401d0*x3 - 0.2565564
     &87706084d0*x7) + x75) - 0.00693220036674015d0*x44*x47 + x69))
      d1E(2) = x206 - x38*(x105*x206 + x106*(x197*x206 + 0.0049653337234
     &2768d0*x207*x48 ) + x114*(x125*x206 + 0.0837628205355044d0*x207*x4
     &2 - 0.00020688890514282d0*x48*(-x13*x201 - x200*(-x198 + x199) + x
     &208 )) + x130*x204 - x134*x142 - 0.0253302959105844d0*x143*x52*(x1
     &33* (-x138 - 30.6251363876636d0*x140 + 3.30853595146347d0*x141) - 
     &x65 *(-x131 + x132)*(11.4095193139813d0*x135 + x138)) + x61*(x161*
     & x204 - x192*(-x13*x162 + x208)))
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E, d2E
     &)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 09, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399, x400, x401
     &, x402, x403, x404, x405, x406, x407, x408, x409, x410, x411, x412
     &, x413, x414, x415, x416, x417, x418, x419, x420, x421, x422, x423
     &, x424, x425, x426, x427, x428, x429, x430, x431, x432, x433, x434
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446, x447, x448, x449, x450, x451, x452, x453, x454, x455, x456
     &, x457, x458, x459, x460, x461, x462, x463, x464, x465, x466, x467
     &, x468, x469, x470, x471, x472, x473, x474, x475, x476, x477, x478
     &, x479, x480, x481, x482, x483, x484, x485, x486, x487, x488, x489
     &, x490, x491, x492, x493, x494, x495, x496, x497, x498, x499, x500
     &, x501, x502, x503, x504, x505, x506, x507, x508, x509, x510, x511
     &, x512, x513, x514, x515, x516, x517, x518, x519, x520, x521, x522
     &, x523, x524, x525, x526, x527, x528, x529, x530, x531, x532, x533
     &, x534, x535, x536, x537, x538, x539, x540, x541, x542, x543, x544
     &, x545, x546, x547, x548, x549, x550, x551, x552, x553, x554, x555
     &, x556, x557, x558, x559, x560, x561, x562, x563, x564, x565, x566
     &, x567, x568, x569, x570, x571, x572, x573, x574, x575, x576, x577
     &, x578, x579, x580, x581, x582, x583, x584, x585, x586, x587, x588
     &, x589, x590, x591, x592, x593, x594, x595, x596, x597, x598, x599
     &, x600, x601, x602, x603, x604, x605, x606, x607
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/rho_c
      x3 = x2**0.333333333333333d0
      x4 = 1d0/x1
      x5 = x2**0.166666666666667d0
      x6 = sqrt(x2)
      x7 = 0.800428634999364d0*x6
      x8 = x2**0.666666666666667d0
      x9 = 2.22556942115069d0*x3 + 5.98255043577108d0*x5 + x7 + 0.189700
     &432574756d0 *x8
      x10 = 4.93480220054468d0*x4/x9 + 1.0d0
      x11 = log(x10)
      x12 = x11*(0.132568899905202d0*x3 + 1.0d0)
      x13 = 0.202642367284676d0*x1*x12
      x14 = rho_s*x2
      x15 = x14 + 1.0d0
      x16 = -x14 + 1.0d0
      x17 = 1.92366105093154d0*x15**1.33333333333333d0 + 1.9236610509315
     &4d0*x16** 1.33333333333333d0 - 3.84732210186307d0
      x18 = rho_s**4/rho_c**4
      x19 = -x18 + 1.0d0
      x20 = 0.430097247127664d0*x6
      x21 = x20 + 2.24759186357762d0*x3 + 8.15741470348764d0*x5 + 0.1911
     &51259512734d0*x8
      x22 = 1.0d0 + 29.6085746432167d0/x21
      x23 = log(x22)
      x24 = x19*x23
      x25 = x17*x24
      x26 = 0.0197517897025652d0*x25*(0.0690139921125583d0*x3 + 1.0d0)
      x27 = 1.64473377556761d0*x6
      x28 = x27 + 3.84474623744721d0*x3 + 11.1203748630947d0*x5 + 0.2405
     &87129128819d0*x8
      x29 = 9.86960440108936d0*x4
      x30 = 1.0d0 + x29/x28
      x31 = log(x30)
      x32 = 0.0506605918211689d0*x31*(0.127469618870009d0*x3 + 1.0d0)
      x33 = 0.101321183642338d0*x12
      x34 = x32 - x33
      x35 = x17*x18
      x36 = x1*x35
      x37 = 2.0d0*x36
      x38 = x34*x37
      x39 = mu**2
      x40 = x39*x8
      x41 = 0.237113969739022d0*x40 + 1.0d0
      x42 = x41**(-4)
      x43 = mu**5
      x44 = rho_s**2/rho_c**2
      x45 = -x44 + 1.0d0
      x46 = x2**1.66666666666667d0
      x47 = exp(-0.46675853320810856d0*x3)
      x48 = x2**1.0d0
      x49 = -0.00304892328706429d0*x48
      x50 = x2**1.33333333333333d0
      x51 = 0.0128480790170175d0*x3 + x49 + 0.000275310793616358d0*x50 +
     & 0.0315297404574782d0*x8 + 1.0d0
      x52 = x47*x51
      x53 = x46*x52
      x54 = x45*x53
      x55 = x43*x54
      x56 = 0.00264117541143784d0*x55
      x57 = x15**0.666666666666667d0
      x58 = x16**0.666666666666667d0
      x59 = x57 + x58
      x60 = x59**3
      x61 = x59**2
      x62 = 1d0/x61
      x63 = x39*x62
      x64 = x3*x63
      x65 = 1d0/x59
      x66 = mu*x65
      x67 = 9.20897059521559d0*x5*x66 + 1.0d0
      x68 = 8.55713948548596d0*x64 + x67
      x69 = 1d0/x68
      x70 = mu**3
      x71 = 1d0/x60
      x72 = x6*x71
      x73 = 1.22772285078429d0/(x0 - 1.0d0)
      x74 = x73 - 3.44851d0
      x75 = -2.4814019635976d0*x64*x74 + x67 + 15.3125681938318d0*x70*x7
     &2
      x76 = x69*x75
      x77 = x1*log(x76)
      x78 = 0.0253302959105844d0*x60*x77
      x79 = x48*x52
      x80 = 0.0445553741830535d0*x45
      x81 = x79*x80
      x82 = -3.07073492995203d0*x3 + 0.384834731559127d0*x8
      x83 = exp(-0.19230865217881402d0*x3)
      x84 = x83/x48
      x85 = x82*x84
      x86 = 37.6991118430775d0*x45
      x87 = x85*x86
      x88 = 1d0/x8
      x89 = x15**2
      x90 = 8.61370200056457d0*x89
      x91 = 1d0/x15
      x92 = x91**0.333333333333333d0
      x93 = x3*x92
      x94 = 0.0177069812999059d0*x93
      x95 = -x94 + 1.0d0
      x96 = x91**0.666666666666667d0
      x97 = 0.0244354823084343d0*x8
      x98 = 0.337569861992026d0*x93 + x96*x97 + 1.0d0
      x99 = 1d0/x98
      x100 = 1d0/x96
      x101 = x100*x99
      x102 = x101*x95
      x103 = x102*x90
      x104 = x16**2
      x105 = 8.61370200056457d0*x104
      x106 = 1d0/x16
      x107 = x106**0.333333333333333d0
      x108 = x107*x3
      x109 = 0.0177069812999059d0*x108
      x110 = -x109 + 1.0d0
      x111 = x106**0.666666666666667d0
      x112 = 0.337569861992026d0*x108 + x111*x97 + 1.0d0
      x113 = 1d0/x112
      x114 = 1d0/x111
      x115 = x113*x114
      x116 = x110*x115
      x117 = x105*x116
      x118 = x46*(x103*x88 + x117*x88 + x87)
      x119 = x70*(0.000132058770571892d0*x118 + x81)
      x120 = -x13 + x26
      x121 = x120 - x38
      x122 = mu**8*x2**2.66666666666667d0
      x123 = 0.00316102962473761d0*x122
      x124 = x121*x123
      x125 = mu**6
      x126 = x52 - 1.0d0
      x127 = x126*x46
      x128 = x127*x45
      x129 = x2**2.0d0
      x130 = -x32 + x33
      x131 = x120 + x130*x37
      x132 = x129*x131
      x133 = x125*(0.00248266686171384d0*x128 - 0.0533250677421794d0*x13
     &2)
      x134 = mu**4
      x135 = x126*x48
      x136 = x135*x45
      x137 = x2**1.0d0
      x138 = x15**2.66666666666667d0
      x139 = 8.61370200056457d0*x138
      x140 = x16**2.66666666666667d0
      x141 = 8.61370200056457d0*x140
      x142 = 0.240695990468967d0*x3
      x143 = 0.26014827853397d0*x8
      x144 = exp(-0.33933171852197186d0*x3)
      x145 = x144*(-x142 + x143)
      x146 = x145*x45
      x147 = 23.3866625383245d0*x146
      x148 = x137*(x103 + x117 - x139 - x141 + x147)
      x149 = 0.337338207872387d0*x50
      x150 = x134*(x121*x149 - 0.0418814102677522d0*x136 - 0.00020688890
     &514282d0*x148 )
      x151 = -x119 + x124 - x133 + x150 - x56 - x78
      x152 = x151*x42
      x153 = 0.0105647016457513d0*x43
      x154 = x44*x53
      x155 = x153*x154
      x156 = 0.00528235082287567d0*x43
      x157 = 0.0042826930056725d0*x3 + x49 + 0.000367081058155144d0*x50 
     &+ 0.0210198269716521d0*x8
      x158 = x46*x47
      x159 = x158*x45
      x160 = x157*x159
      x161 = x156*x160
      x162 = 0.00880391803812612d0*x55
      x163 = x2**2.0d0
      x164 = x163*x43*x45
      x165 = x164*x52
      x166 = 0.000821860773992031d0*x165
      x167 = x15**(-0.333333333333333d0)
      x168 = x16**(-0.333333333333333d0)
      x169 = x167 - x168
      x170 = 0.101321183642338d0*x169
      x171 = x61*x77
      x172 = x14*x171
      x173 = x170*x172
      x174 = x121*x122
      x175 = 0.0168588246652672d0*x174
      x176 = x3*x66
      x177 = mu*x3
      x178 = x177*x62
      x179 = x14*x169
      x180 = x178*x179
      x181 = 6.13931373014372d0*x5
      x182 = x14*x65
      x183 = x169*x182
      x184 = x181*x183 - 1.53482843253593d0*x5
      x185 = -2.85237982849532d0*x176 + 11.4095193139813d0*x180 + x184
      x186 = x176*x74
      x187 = x178*x74
      x188 = x179*x187
      x189 = x6*x63
      x190 = x39*x72
      x191 = x179*x190
      x192 = x184 - 7.6562840969159d0*x189 + 30.6251363876636d0*x191
      x193 = 0.827133987865867d0*x186 - 3.30853595146347d0*x188 + x192
      x194 = x1*(-x185*x76 + x193)
      x195 = 1d0/x75
      x196 = mu*x195
      x197 = x196*x61
      x198 = 0.0506605918211689d0*x194*x197
      x199 = 16.0d0*x36
      x200 = rho_s**5/rho_c**5
      x201 = 10.2595256049682d0*x200
      x202 = x16**0.333333333333333d0
      x203 = x15**0.333333333333333d0
      x204 = x202 - x203
      x205 = x1*x204
      x206 = x201*x205
      x207 = x1*x3
      x208 = x11*x207
      x209 = 0.0179093838034102d0*x208
      x210 = x207*x31
      x211 = 0.00861024844156464d0*x210
      x212 = 0.741856473716896d0*x3 + 0.99709173929518d0*x5 + 0.40021431
     &7499682d0*x6 + 0.126466955049837d0*x8
      x213 = 1d0/x10
      x214 = x9**(-2)
      x215 = 0.265137799810404d0*x3 + 2.0d0
      x216 = x213*x214*x215
      x217 = x212*x216
      x218 = 1.28158207914907d0*x3 + 1.85339581051578d0*x5 + 0.822366887
     &783805d0*x6 + 0.160391419419213d0*x8
      x219 = 1d0/x30
      x220 = x28**(-2)
      x221 = 0.254939237740017d0*x3 + 2.0d0
      x222 = x219*x220*x221
      x223 = x218*x222
      x224 = -x209 + x211 + x217 - x223
      x225 = x224*x35
      x226 = x25*x3
      x227 = 0.00136314985874174d0*x3 + 0.0197517897025652d0
      x228 = x227*x23
      x229 = x228*x35
      x230 = -x217
      x231 = x209 + x230
      x232 = x227*x24
      x233 = x14*x232
      x234 = x17*x19*x227
      x235 = 1d0/x22
      x236 = x21**(-2)
      x237 = 0.749197287859205d0*x3 + 1.35956911724794d0*x5 + 0.21504862
     &3563832d0*x6 + 0.127434173008489d0*x8
      x238 = x235*x236*x237
      x239 = x234*x238
      x240 = 5.1297628024841d0*x204*x233 + 59.2171492864334d0*x239
      x241 = x225 - 0.000908766572494496d0*x226 + 8.0d0*x229 + x231 + x2
     &40
      x242 = x199*x34 - x206*x34 + x241
      x243 = x123*x242
      x244 = 0.0891107483661069d0*x79
      x245 = x47*x48
      x246 = x245*x80
      x247 = x45*x50
      x248 = x247*x52
      x249 = 1.02357830998401d0*x3 - 0.256556487706084d0*x8
      x250 = x84*x86
      x251 = x2**(-0.666666666666667d0)
      x252 = x45*x82*x83
      x253 = x251*x252
      x254 = 75.398223686155d0*x85
      x255 = x254*x44
      x256 = 5.74246800037638d0*x88
      x257 = x102*x89
      x258 = x256*x257
      x259 = x104*x116
      x260 = x256*x259
      x261 = 1d0/x3
      x262 = 0.0508408867489862d0*x261
      x263 = 1d0/x92
      x264 = x263*x99
      x265 = x14*x91
      x266 = -x265 + 1.0d0
      x267 = x266*x89
      x268 = x264*x267
      x269 = x262*x268
      x270 = 1d0/x107
      x271 = x113*x270
      x272 = x106*x14
      x273 = x272 + 1.0d0
      x274 = x104*x273
      x275 = x271*x274
      x276 = x262*x275
      x277 = 0.112523287330675d0*x93
      x278 = 0.0162903215389562d0*x8
      x279 = x278*x96
      x280 = x277 + x279
      x281 = -x265*x277 - x265*x279 + x280
      x282 = x98**(-2)
      x283 = x100*x282
      x284 = x283*x95
      x285 = x281*x284
      x286 = x285*x88
      x287 = x286*x90
      x288 = 0.112523287330675d0*x108
      x289 = x111*x278
      x290 = x288 + x289
      x291 = x272*x288 + x272*x289 + x290
      x292 = x112**(-2)
      x293 = x114*x292
      x294 = x110*x293
      x295 = x291*x294*x88
      x296 = x105*x295
      x297 = 22.9698720015055d0*x88
      x298 = x116*x16
      x299 = x297*x298
      x300 = x14*x299
      x301 = x102*x15
      x302 = x297*x301
      x303 = x14*x302
      x304 = x258 + x260 + x269 + x276 + x287 + x296 + x300 - x303
      x305 = x249*x250 + 2.41662179562687d0*x253 + x255 + x304 + x87
      x306 = 0.000132058770571892d0*x46
      x307 = 2.0d0*x70*(0.000220097950953153d0*x118 + x157*x246 - x244*x
     &44 - 0.00693220036674015d0*x248 - x305*x306 + x81)
      x308 = x127*x44
      x309 = x3*x51
      x310 = 0.00214134650283625d0*x3
      x311 = 0.0105099134858261d0*x8
      x312 = 0.00152446164353214d0*x48
      x313 = 0.000183540529077572d0*x50
      x314 = -x310 - x311 + x312 - x313
      x315 = 0.0777930888680181d0*x309 + x314
      x316 = x315*x45
      x317 = x158*x316
      x318 = -x130*x199 + x130*x206 + x241
      x319 = 0.0533250677421794d0*x129
      x320 = x125*(-0.0082755562057128d0*x128 + 0.213300270968718d0*x132
     & + 0.00993066744685536d0*x308 + 0.00993066744685536d0*x317 - x318*
     & x319)
      x321 = x135*x44
      x322 = 0.0837628205355044d0*x135
      x323 = x245*x316
      x324 = x88*(x138 + x140)
      x325 = 5.74246800037638d0*x324
      x326 = x16**1.66666666666667d0
      x327 = x15**1.66666666666667d0
      x328 = x326 - x327
      x329 = 22.9698720015055d0*x88
      x330 = x14*x328*x329
      x331 = x145*x88
      x332 = 46.773325076649d0*x331
      x333 = x332*x44
      x334 = x146*x88
      x335 = 15.591108358883d0*x334
      x336 = 2.64527879654103d0*x146*x261
      x337 = x144*x88
      x338 = 0.0802319968229891d0*x3 - 0.173432185689313d0*x8
      x339 = x338*x45
      x340 = x337*x339
      x341 = 23.3866625383245d0*x340
      x342 = x134*(0.899568554326366d0*x131*x50 + 0.0006896296838094d0*x
     &137*(-x103 - x117 + x139 + x141 - x147) - x149*x318 + 0.1675256410
     &71009d0*x321 - x322*x45 + 0.167525641071009d0*x323 - 0.00041377781
     &028564d0*x46 *(-x258 - x260 - x269 - x276 - x287 - x296 - x300 + x
     &303 + x325 + x330 - x333 - x335 - x336 - x341))
      x343 = -x155 + x161 + x162 - x166 + x173 - x175 - x198 + x243 + x3
     &07 - x320 - x342
      x344 = x40/x41**5
      x345 = x151*x344
      x346 = x34*x36
      x347 = x200*x34
      x348 = 5.1297628024841d0*x205
      x349 = 2.56488140124205d0*x204
      x350 = x232*x349
      x351 = x156*x53
      x352 = 0.0506605918211689d0*x169
      x353 = rho_s**3/rho_c**3
      x354 = x17*x353
      x355 = x228*x354
      x356 = -x350
      x357 = x18*x34
      x358 = 8.0d0*x1
      x359 = x34*x354
      x360 = x348*x357 - 4.0d0*x355 + x356 - x358*x359
      x361 = x196*x59
      x362 = 30.6251363876636d0*x189
      x363 = x169*(-x181 + 3.30853595146347d0*x186 - x362)
      x364 = -x167 + x168
      x365 = 11.4095193139813d0*x176 + x181
      x366 = x364*x365
      x367 = x1*(x363 - x366*x76)
      x368 = 0.0253302959105844d0*x367
      x369 = 0.00496533372342768d0*x127
      x370 = x14*x254
      x371 = x15*x264
      x372 = x16*x271
      x373 = 8.61370200056457d0*x88
      x374 = x280*x284
      x375 = x15*x374
      x376 = x290*x294
      x377 = x16*x376
      x378 = x262*x371 - x262*x372 - x299 + x302 + x373*x375 - x373*x377
     &
      x379 = -x370 + x378
      x380 = -x14*x332 - x329*(-x326 + x327) + x378
      x381 = 0.00020688890514282d0*x46
      x382 = x123*x360 + x125*(x14*x369 + x319*x360) + x134*(x14*x322 + 
     &x149*x360 - x380*x381) + x14*x351 - x171*x352 - x361*x368 + x70*(x
     &14*x244 - x306*x379)
      x383 = -x382*x42
      x384 = x157*x44
      x385 = x43*x52
      x386 = x163*x385
      x387 = x157*x47
      x388 = 0.00571025734089666d0*x3 - 0.00609784657412857d0*x48 + 0.00
     &0856522469028669d0*x50 + 0.0350330449527535d0*x8
      x389 = x169**2
      x390 = x389*x59*x77
      x391 = x15**(-1.33333333333333d0)
      x392 = x14*x391
      x393 = x16**(-1.33333333333333d0)
      x394 = x14*x393
      x395 = x392 + x394
      x396 = -6.0d0*x167 + 6.0d0*x168 + x395
      x397 = x194*x361
      x398 = 0.101321183642338d0*x59
      x399 = x194*x39
      x400 = x195*x69
      x401 = x399*x400
      x402 = -x73 + 3.44851d0
      x403 = x176*x402
      x404 = x75**(-2)
      x405 = x399*x404
      x406 = 3.80317310466043d0*x176
      x407 = x14*x396
      x408 = x389*x44
      x409 = x177*x408*x71
      x410 = 2.04643791004791d0*x5
      x411 = 8.18575164019163d0*x5
      x412 = x182*x396*x410 - x183*x410 + x408*x411*x62
      x413 = 1.10284531715449d0*x186
      x414 = x185*x66
      x415 = 2.0d0*x69
      x416 = 12.2786274602874d0*x5
      x417 = 22.8190386279626d0*x178
      x418 = x179*x417
      x419 = x75/x68**2
      x420 = x414*x419
      x421 = x15**(-0.666666666666667d0)
      x422 = 0.854960467080683d0*x14
      x423 = x16**(-0.666666666666667d0)
      x424 = x421*x422 + x422*x423
      x425 = -5.1297628024841d0*x202 + 5.1297628024841d0*x203 + x424
      x426 = x358*x425
      x427 = 164.152409679491d0*x205
      x428 = 0.0477583568090938d0*x208
      x429 = x212*x213
      x430 = x214*x3*x429
      x431 = 0.353517066413871d0*x430
      x432 = x218*x219
      x433 = x216*(0.989141964955861d0*x3 + 1.16327369584438d0*x5 + 0.60
     &0321476249523d0*x6 + 0.210778258416395d0*x8)
      x434 = 2.0d0*x433
      x435 = x212**2*x215/(x10**2*x9**4)
      x436 = x29*x435
      x437 = x215*x429*(1.48371294743379d0*x3 + 1.99418347859036d0*x5 + 
     &x7 + 0.252933910099675d0*x8)/x9**3
      x438 = 2.0d0*x437
      x439 = x35*(-0.022960662510839d0*x210 + 19.7392088021787d0*x218**2
     &*x221*x4/(x28 **4*x30**2) + 0.339918983653357d0*x220*x3*x432 - 2.0
     &d0*x221*x432* (x27 + 2.56316415829814d0*x3 + 3.70679162103156d0*x5
     & + 0.320782838838426d0*x8)/x28**3 + 2.0d0*x222*(1.70877610553209d0
     &* x3 + 2.16229511226841d0*x5 + 1.23355033167571d0*x6 + 0.267319032
     &365355d0*x8) + x428 - x431 - x434 - x436 + x438)
      x440 = x204*x224
      x441 = x204*x228
      x442 = x200*x441
      x443 = x23*x3
      x444 = x35*x443
      x445 = x233*x425
      x446 = x204*x24*x3
      x447 = x14*x446
      x448 = x19*x238
      x449 = x17*x3*x448
      x450 = x234*x235
      x451 = 118.434298572867d0*x450
      x452 = x236*(0.998929717145607d0*x3 + 1.5861639701226d0*x5 + 0.322
     &572935345748d0*x6 + 0.212390288347482d0*x8)
      x453 = x234*x237**2/(x21**4*x22**2)
      x454 = x237*(x20 + 1.49839457571841d0*x3 + 2.71913823449588d0*x5 +
     & 0.254868346016978d0*x8)/x21**3
      x455 = x227*x238
      x456 = x35*x455
      x457 = x204*x227*x448
      x458 = x14*x457
      x459 = x201*x440 - 16.0d0*x225 + 0.00242337752665199d0*x226 - 80.0
     &d0*x229 - x428 + x431 + x434 + x436 - x438 + x439 + 82.07620483974
     &55d0*x442 - 0.0145402651599119d0*x444 + 4.0d0*x445 - 0.00932351391
     &944647d0* x447 - 0.107629131579854d0*x449 - x451*x452 + x451*x454 
     &- 3506.67076961173d0*x453 + 947.474388582934d0*x456 + 607.53985935
     &7387d0*x458
      x460 = -160.0d0*x346 - x347*x426 + x347*x427 + x459
      x461 = 0.00316102962473761d0*x122
      x462 = x50*x52
      x463 = x245*x45
      x464 = 0.0891107483661069d0*x157
      x465 = x249*x45
      x466 = 75.398223686155d0*x84
      x467 = x251*x83
      x468 = 4.83324359125373d0*x467
      x469 = 0.0677878489986482d0*x261
      x470 = 1d0/x89
      x471 = 0.00786976946662482d0*x44
      x472 = x261*x264
      x473 = 1.91415600012546d0*x88
      x474 = 1d0/x104
      x475 = x261*x271
      x476 = 38.2831200025092d0*x88
      x477 = x102*x476
      x478 = x116*x476
      x479 = 0.0543010717965206d0*x8
      x480 = 0.150031049774234d0*x93
      x481 = x44*x470
      x482 = 0.0271505358982603d0*x8
      x483 = x482*x96
      x484 = -x480 - x483
      x485 = x88*x90
      x486 = 0.101681773497972d0*x261
      x487 = x263*x282
      x488 = x281*x487
      x489 = x270*x292
      x490 = x486*x489
      x491 = 0.150031049774234d0*x108
      x492 = x44*x474
      x493 = x111*x482
      x494 = x491 + x493
      x495 = x105*x88
      x496 = x261*x273*x372
      x497 = 0.271151395994593d0*x14
      x498 = x14*x88
      x499 = 15.3132480010037d0*x498
      x500 = 0.225046574661351d0*x93
      x501 = 0.0325806430779124d0*x8
      x502 = x501*x96
      x503 = x500 + x502
      x504 = x98**(-3)
      x505 = x100*x504*x95
      x506 = 0.225046574661351d0*x108
      x507 = x111*x501
      x508 = x506 + x507
      x509 = x112**(-3)
      x510 = x110*x114*x509
      x511 = x261*x371
      x512 = x266*x511
      x513 = x16*x291
      x514 = 45.939744003011d0*x88
      x515 = x14*x514
      x516 = x15*x285
      x517 = 11.4849360007528d0*x104*x295 - x105*x475*(0.015739538933249
     &6d0*x272 + x471*x474 + 0.00786976946662482d0) - x257*x473 - x259*x
     &473 + x267 *x486*x488 + x268*x469 + x274*x291*x490 + x275*x469 + x
     &281*x485* x505*(-x265*x500 - x265*x502 + x503) + x284*x485*(x265*x
     &479*x96 + 0.300062099548467d0*x265*x93 - x480*x481 - x481*x483 + x
     &484) + 11.4849360007528d0*x286*x89 + x291*x495*x510*(x272*x506 + x
     &272* x507 + x508) - x294*x495*(0.300062099548467d0*x108*x272 + x11
     &1* x272*x479 + x491*x492 + x492*x493 + x494) + x294*x513*x515 - x2
     &98 *x499 + x301*x499 + x44*x477 + x44*x478 + x472*x90*( 0.01573953
     &89332496d0*x265 - x470*x471 - 0.00786976946662482d0) + x496*x497 -
     & x497*x512 - x515*x516
      x518 = x315*x44
      x519 = 0.311172355472072d0*x3
      x520 = 0.00285512867044833d0*x3 - 0.103724118490691d0*x309 + x49 +
     & 0.000428261234514334d0*x50 + 0.0121035293512547d0*x51*x8 + 0.0175
     &165224763768d0*x8
      x521 = x130*x200
      x522 = 5.29055759308205d0*x144*x261
      x523 = x328*x88
      x524 = 5.0d0*x14
      x525 = x524*x57 + x524*x58
      x526 = 7.65662400050184d0*x88
      x527 = x337*x338
      x528 = 0.632303919304059d0*x344
      x529 = 2.0d0*x1
      x530 = 29.6085746432167d0*x450
      x531 = x14*x158
      x532 = 0.0675474557615585d0*x390
      x533 = 3.0d0*x167
      x534 = 3.0d0*x168
      x535 = x395 - x533 + x534
      x536 = 0.0168868639403896d0*x171
      x537 = x196*x367
      x538 = x364*(x181 + x362 + 3.30853595146347d0*x403)
      x539 = x169*x365
      x540 = x1*x359
      x541 = x205*x357
      x542 = -2.56488140124205d0*x202 + 2.56488140124205d0*x203 + x424
      x543 = 20.5190512099364d0*x441
      x544 = -1.28244070062102d0*x18*x440 - x18*x543 + 2.0d0*x224*x354 -
     & x232*x542 + 0.00181753314498899d0*x354*x443 - 118.434298572867d0*
     &x354*x455 + x357*x529*x542 + 0.00116543923993081d0*x446 - 75.94248
     &24196734d0* x457
      x545 = 16.0d0*x355 + 32.0d0*x540 - 41.0381024198728d0*x541 + x544
      x546 = 1.02321895502395d0*x5
      x547 = x411*x65
      x548 = x389*x547
      x549 = x190*x389
      x550 = 6.61707190292693d0*x187*x389
      x551 = 10.2083787958879d0*x189
      x552 = 22.8190386279626d0*x176 + x416
      x553 = -x392 - x394 + x533 - x534
      x554 = 0.0253302959105844d0*x1*x361
      x555 = x14*x245
      x556 = 0.135575697997296d0*x14
      x557 = x472*x556
      x558 = 0.0338939244993241d0*x261
      x559 = x371*x558
      x560 = x475*x556
      x561 = x372*x558
      x562 = 38.2831200025092d0*x498
      x563 = 7.65662400050184d0*x88
      x564 = 0.0677878489986482d0*x496
      x565 = x280*x487
      x566 = x15*x262
      x567 = x266*x565*x566
      x568 = x14*x297
      x569 = x262*x489
      x570 = x16*x273*x290*x569
      x571 = x488*x566
      x572 = x513*x569
      x573 = x265*x480
      x574 = x265*x483
      x575 = x284*x373
      x576 = x297*x513
      x577 = x272*x491 + x272*x493 + x494
      x578 = x294*x373
      x579 = x373*x503
      x580 = x505*x579
      x581 = x373*x508
      x582 = x510*x581
      x583 = x14*(x142 - x143)
      x584 = x94 - 1.0d0
      x585 = x15*x584
      x586 = x109 - 1.0d0
      x587 = x16*x586
      x588 = x283*x585
      x589 = x293*x586
      x590 = x293*x587
      x591 = x480 + x483
      x592 = x353*x543
      x593 = x421 + x423
      x594 = 0.854960467080683d0*x232*x593
      x595 = x17*x44
      x596 = 12.0d0*x228*x595
      x597 = 1.70992093416137d0*x1*x357*x593
      x598 = 41.0381024198728d0*x205*x34*x353
      x599 = 24.0d0*x1*x34*x595
      x600 = x391 + x393
      x601 = x368*x39*x65
      x602 = x592 + x594 - x596 - x597 + x598 - x599
      x603 = 0.203363546995945d0*x261
      x604 = x264*x603 + x271*x603 + x280*x580 + x290*x490 + x290*x582 +
     & x374*x514 + x376*x514 + x477 + x478 + x486*x565 - x494*x578 - x57
     &5*x591
      x605 = x410*x600
      x606 = mu*x62
      x607 = x364**2
      E = -rho_c*(x13 + x152 - x26 + x38)
      d1E(1) = x120 + x14*x350 - x152 + 0.0089546919017051d0*x208 - 0.5d
     &0*x217 + 0.5d0* x225 - 0.000454383286247248d0*x226 + 4.0d0*x229 + 
     &29.6085746432167d0*x239 - 0.5d0*x343*x42 - 0.632303919304059d0* x3
     &45 + 6.0d0*x346 - x347*x348
      d1E(2) = x360 + x383
      d2E(1) = x2*(-0.499760307959092d0*x134*x151*x50/x41**6 + x200*x224
     &*x349 + 30.7785768149046d0*x205*x347 + 0.00596979460113673d0*x208 
     &- 4.0d0 *x225 - 0.000302922190831499d0*x226 - 12.0d0*x229 + x230 +
     & x240 - x343*x528 + 1.26460783860812d0*x344*(x119 - x124 + x133 - 
     &x150 + x56 + x78) + 1.05383986550676d0*x345 - 24.0d0*x346 - x347*x
     &425* x529 - x35*(-x211 + x223 + x231) + 2.46740110027234d0*x4*x435
     & + x42*(x155 - x161 - x162 + x166 - x173 + x175 + x198 - x243 - x3
     &07 + x320 + x342) - 0.25d0*x42*(-0.101321183642338d0*x1*x197*( -1.
     &10284531715449d0*x187*x407 + 2.20569063430898d0*x188 + 11.48442614
     &53738d0*x189 + 10.2083787958879d0*x190*x407 - 30.6251363876636d0*x
     &191 - x193*x414*x415 + 81.6670303671029d0*x39 *x408*x6/x59**4 - 6.
     &61707190292693d0*x409*x74 + x412 - x413 + x420*(-5.70475965699064d
     &0*x176 + x183*x416 + x418 - 3.06965686507186d0*x5) + 1.79063317129
     &192d0*x5 - x76*( 3.80317310466043d0*x178*x407 - 7.60634620932085d0
     &*x180 + x406 + 22.8190386279626d0*x409 + x412 + 1.79063317129192d0
     &*x5)) - 0.0337176493305345d0*x122*x242 - x125*(0.0441362997638016d
     &0*x128 + 0.426600541937435d0*x129*x318 - 1.27980162581231d0*x132 +
     & 0.0794453395748429d0*x158*x518 + 0.0198613348937107d0*x159*(x314*
     & x519 + x520) - 0.125788454326835d0*x308 - 0.0662044496457024d0* x
     &317 - x319*(160.0d0*x130*x36 + x426*x521 - x427*x521 + x459)) + x1
     &34*(4.19798658685637d0*x121*x50 - 0.335051282142018d0*x136 - 0.003
     &6780249803168d0*x148 + x149*x460 - 1.79913710865273d0*x242* x50 - 
     &1.34020512856807d0*x245*x518 + 1.67525641071009d0*x321 + 0.6701025
     &64284035d0*x323 + 0.0027585187352376d0*x46*(x304 - x325 - x330 + x
     &333 + x335 + x336 + x341) - 0.00082755562057128d0*x46*( -30.626496
     &0020074d0*x14*x523 - x14*x526*(-6.0d0*x326 + 6.0d0*x327 + x525) + 
     &10.5811151861641d0*x145*x261*x44 + 0.299209d0*x146 + 1.91415600012
     &546d0*x324 - 77.9555417944151d0*x331*x44 - 5.19703611962767d0*x334
     & + 23.3866625383245d0*x337*x45*( -0.106975995763985d0*x3 + 0.28905
     &3642815522d0*x8) + x339*x522 + 31.182216717766d0*x340 + 93.5466501
     &532981d0*x44*x527 + x517) - 0.335051282142018d0*x463*(-x519*(x310 
     &+ x311 - x312 + x313) + x520)) - x153*x159*x388 + 0.13381955417951
     &7d0*x154*x43 + 0.0422588065830054d0*x158*x384*x43 - 0.035215672152
     &5045d0*x160* x43 + 0.00328744309596812d0*x164*x387 + 0.00767070055
     &725895d0* x165 + 0.0675474557615585d0*x172*x396 + 0.12363138087862
     &6d0*x174 + 0.405284734569351d0*x179*x397 - x185*x398*x401 - 0.0002
     &55740352913201d0*x2**2.33333333333333d0*x385*x45 - 0.0065748861919
     &3625d0*x386*x44 - 0.270189823046234d0*x390*x44 + x398*x405*(3.3085
     &3595146347d0*x180*x402 + x192 - 0.827133987865867d0*x403) + x460*x
     &461 - 0.0469542295366726d0*x55 + 4.0d0*x70*(-0.000586927869208408d
     &0*x118 - x244*x45 + 0.178221496732214d0*x245*x384 - x246*x388 + 0.
     &0138644007334803d0* x247*x387 + 0.0231073345558005d0*x248 + 0.0004
     &40195901906306d0* x305*x46 - x306*(0.154912426780983d0*x2**(-0.333
     &333333333333d0)* x252 + 150.79644737231d0*x249*x44*x84 + x250*(-1.
     &36477107997868d0 *x3 + 0.427594146176807d0*x8) + 1.61108119708458d
     &0*x253 - x255 + 9.66648718250747d0*x44*x467*x82 + x465*x466 + x465
     &*x468 + x517) - 0.0277288014669606d0*x44*x462 + 0.445553741830535d
     &0*x44*x79 - x463*x464 - 0.00107855455836145d0*x54)) + 0.0883792666
     &034679d0* x430 + 0.5d0*x433 - 0.5d0*x437 + 0.25d0*x439 + 20.519051
     &2099364d0 *x442 - 0.00363506628997799d0*x444 + x445 - 0.0023308784
     &7986162d0 *x447 - 0.0269072828949635d0*x449 - x452*x530 - 876.6676
     &92402934d0*x453 + x454*x530 + 236.868597145733d0*x456 + 151.884964
     &839347d0*x458)
      d2E(2) = x2*(12.0d0*x355 + x356 - x382*x528 + x383 - x42*(-0.00842
     &941233263362d0* x122*x360 + x125*(-0.0182062236525682d0*x127*x14 -
     & 0.106650135484359d0*x129*x360 + 0.00993066744685536d0*x315*x531 +
     & x319*x545) + x134*(-0.251288461606513d0*x135*x14 + x149*x545 + 0.
     &167525641071009d0*x315*x555 - 0.449784277163183d0*x360*x50 + 0.000
     &3448148419047d0*x380*x46 + x381*(x100*x281*x504*x579*x585 - x101*x
     &562*x584 - x101*x563*x585 - x114*x509*x513*x581*x586 - x115 *x562*
     &x586 + x115*x563*x587 + 46.773325076649d0*x14*x527 + x256* x280*x5
     &88 - x256*x290*x590 - x280*x283*x568*x584 + x281*x297*x588 - x290*
     &x568*x589 + 62.364433435532d0*x337*x583 + x373*x577*x590 - x373*x5
     &88*(-x573 - x574 + x591) - 0.0677878489986482d0*x512 - x522*x583 -
     & 15.3132480010037d0*x523 - x526*(-3.0d0*x326 + 3.0d0* x327 + x525)
     & + x557 - x559 + x560 + x561 + x564 - x567 + x570 - x571 + x572 - 
     &x576*x589)) + x14*x352*x537 + 0.000821860773992031d0*x14*x386 - 0.
     &0193686196838775d0*x14*x43* x53 + x14*x532 - x156*x157*x531 - x352
     &*x397 + 0.0253302959105844d0*x401*x539 + 0.0253302959105844d0*x405
     &*x538 + x461*x545 - x535*x536 - x554*(-x14*x548 - 81.6670303671029
     &d0*x14* x549 + x14*x550 + 15.3125681938318d0*x169*x189 - x169*x413
     & + x169 *x546 - x193*x366*x66*x69 - x363*x414*x69 + x364*x420*x552
     & - x410 *x535 + x413*x535 - x535*x551 - x76*(x179*x364*x547 - x364
     &*x406 + x364*x418 - x364*x546 + x406*x553 + x410*x553)) + x70*( 0.
     &0138644007334803d0*x14*x462 - 0.267332245098321d0*x14*x79 - x306*(
     &-x102*x562 - x116*x562 - x14*x249*x466 - x14*x468*x82 + x15 *x281*
     &x580 + x15*x575*(x484 + x573 + x574) + x16*x577*x578 + x256 *x375 
     &- x256*x377 - x294*x576 + x297*x516 + x298*x563 - x301*x563 + x370
     & - x374*x568 - x376*x568 + 0.0677878489986482d0*x511*(x265 - 1.0d0
     &) + 0.135575697997296d0*x512 - x513*x582 - x557 + x559 - x560 - x5
     &61 - x564 + x567 - x570 + x571 - x572) + 0.000220097950953153d0*x3
     &79*x46 - x464*x555)) + 24.0d0*x540 - 35.9083396173887d0*x541 + x54
     &4)
      d2E(3) = -x2*(x42*(x125*(x319*x602 + x369) + x134*(x149*x602 + x32
     &2 - x381*(-x332 - 38.2831200025092d0*x59*x88 + x604)) - x170*x537 
     &+ x351 + x400* x539*x601 + x404*x538*x601 + x461*x602 - x532 + x53
     &6*x600 - x554* (-x363*x366*x415*x606 + x365*x419*x552*x606*x607 - 
     &x413*x600 + x548 + 81.6670303671029d0*x549 - x550 + x551*x600 + x6
     &05 - x76*( x406*x600 + x417*x607 + x547*x607 + x605)) + x70*(x244 
     &- x306*( -x254 + x604))) - x592 - x594 + x596 + x597 - x598 + x599
     &)
      end subroutine
    

C*****************************************************************************
      pure subroutine ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu,
     &                                      E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 09, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75
      E = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/rho_c
      x3 = x2**0.333333333333333d0
      x4 = 1d0/x1
      x5 = x2**0.166666666666667d0
      x6 = sqrt(x2)
      x7 = x2**0.666666666666667d0
      x8 = (0.132568899905202d0*x3 + 1.0d0)*log(4.934802200544679d0*x4/(
     & 2.2255694211506878d0*x3 + 5.9825504357710795d0*x5 + 0.80042863499
     &936356d0*x6 + 0.18970043257475591d0*x7) + 1.0d0)
      x9 = 0.202642367284676d0*x1*x8
      x10 = rho_s**4/rho_c**4
      x11 = rho_s*x2
      x12 = x11 + 1.0d0
      x13 = -x11 + 1.0d0
      x14 = 1.92366105093154d0*x12**1.33333333333333d0 + 1.9236610509315
     &4d0*x13** 1.33333333333333d0 - 3.84732210186307d0
      x15 = 0.0197517897025652d0*x14*(-x10 + 1.0d0)*(0.0690139921125583d
     &0*x3 + 1.0d0 )*log(1.0d0 + 29.608574643216677d0/(2.247591863577616
     &4d0*x3 + 8.1574147034876408d0*x5 + 0.43009724712766434d0*x6 + 0.19
     &115125951273382d0*x7))
      x16 = 0.0506605918211689d0*(0.127469618870009d0*x3 + 1.0d0)*log( 9
     &.869604401089358d0*x4/(3.8447462374472119d0*x3 + 11.12037486309468
     &5d0*x5 + 1.6447337755676092d0*x6 + 0.24058712912881922d0*x7) + 1.0
     &d0)
      x17 = 0.101321183642338d0*x8
      x18 = 2.0d0*x1*x10*x14
      x19 = x18*(x16 - x17)
      x20 = mu**2
      x21 = (0.237113969739022d0*x20*x7 + 1.0d0)**(-4)
      x22 = rho_s**2/rho_c**2
      x23 = -x22 + 1.0d0
      x24 = x2**1.0d0
      x25 = x2**1.33333333333333d0
      x26 = (-0.00304892328706429d0*x24 + 0.000275310793616358d0*x25 + 0
     &.0128480790170175d0*x3 + 0.0315297404574782d0*x7 + 1.0d0)*exp( -0.
     &46675853320810856d0*x3)
      x27 = x23*x26
      x28 = x2**1.66666666666667d0
      x29 = 0.00264117541143784d0*mu**5*x28
      x30 = x12**0.666666666666667d0 + x13**0.666666666666667d0
      x31 = x20*x3/x30**2
      x32 = 9.20897059521559d0*mu*x5/x30 + 1.0d0
      x33 = 1d0/(8.55713948548596d0*x31 + x32)
      x34 = x0 - 1.0d0
      x35 = 1.22772285078429d0/x34
      x36 = 2.4814019635976d0*x31
      x37 = mu**3
      x38 = x30**3
      x39 = x32 + 15.3125681938318d0*x37*x6/x38
      x40 = 0.0253302959105844d0*x38
      x41 = 0.0445553741830535d0*x24
      x42 = 1d0/x7
      x43 = 1d0/x12
      x44 = x43**0.666666666666667d0
      x45 = x3*x43**0.333333333333333d0
      x46 = 0.0244354823084343d0*x7
      x47 = 8.61370200056457d0*x12**2*(-0.0177069812999059d0*x45 + 1.0d0
     &)/(x44*(x44* x46 + 0.337569861992026d0*x45 + 1.0d0))
      x48 = 1d0/x13
      x49 = x48**0.666666666666667d0
      x50 = x3*x48**0.333333333333333d0
      x51 = 8.61370200056457d0*x13**2*(-0.0177069812999059d0*x50 + 1.0d0
     &)/(x49*(x46* x49 + 0.337569861992026d0*x50 + 1.0d0))
      x52 = 0.000132058770571892d0*x28*(37.6991118430775d0*x23*(-3.07073
     &492995203d0* x3 + 0.384834731559127d0*x7)*exp(-0.19230865217881402
     &d0*x3)/x24 + x42*x47 + x42*x51)
      x53 = x15 - x9
      x54 = -x19 + x53
      x55 = 0.00316102962473761d0*mu**8*x2**2.66666666666667d0
      x56 = mu**6
      x57 = x23*(x26 - 1.0d0)
      x58 = 0.00248266686171384d0*x28*x57
      x59 = x18*(-x16 + x17) + x53
      x60 = 0.0533250677421794d0*x2**2.0d0*x59
      x61 = mu**4
      x62 = -0.0418814102677522d0*x24*x57
      x63 = 8.61370200056457d0*x12**2.66666666666667d0
      x64 = 8.61370200056457d0*x13**2.66666666666667d0
      x65 = 23.3866625383245d0*x23*(-0.240695990468967d0*x3 + 0.26014827
     &853397d0*x7) *exp(-0.33933171852197186d0*x3)
      x66 = 0.00020688890514282d0*x2**1.0d0
      x67 = 0.337338207872387d0*x25
      x68 = x21*(-x1*x40*log(x33*(-x36*(x35 - 3.4485100000000002d0) + x3
     &9)) - x27* x29 - x37*(x27*x41 + x52) + x54*x55 - x56*(x58 - x60) +
     & x61*(x54* x67 + x62 - x66*(x47 + x51 - x63 - x64 + x65)))
      x69 = x12**0.666666666666667d0 + x13**0.666666666666667d0
      x70 = x69**3
      x71 = x26*(x22 - 1.0d0)
      x72 = -1.0d0 + exp(-257.30918915441765d0*(-x21*(x29*x71 + x34*x40*
     &log(x33*(x36 *(-x35 + 3.4485100000000002d0) + x39)) + x37*(x41*x71
     & - x52) + x55*x59 + x56*(-x58 + x60) + x61*(x59*x67 + x62 + x66*(-
     &x47 - x51 + x63 + x64 - x65))) + x59)/x70)
      x73 = (x54 - x68)/x54
      x74 = 0.544841407721756d0*gamma_cc*rho_c**(-2.33333333333333d0)*x7
     &3**2.78d0/ x69**2
      x75 = 1.0d0 + x74/x72
      E = -rho_c*(-x15 + x19 + x68 - 0.003886375d0*x70*log(x74*x75/( 0.2
     &9685215956822486d0*gamma_cc**2*rho_c**(-4.666666666666667d0)* x73*
     &*5.5599999999999996d0/(x69**4*x72**2) + x75) + 1.0d0) + x9)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc,
     &                                        mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 09, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399, x400, x401
     &, x402, x403, x404, x405, x406, x407, x408, x409, x410, x411, x412
     &, x413, x414, x415, x416, x417, x418, x419, x420, x421, x422, x423
     &, x424, x425, x426, x427, x428, x429, x430, x431, x432, x433, x434
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.693147180559945d0
      x3 = -x2 + 1.0d0
      x4 = 1d0/x3
      x5 = x0**0.166666666666667d0
      x6 = sqrt(x0)
      x7 = x0**0.666666666666667d0
      x8 = 3.84474623744721d0*x1 + 11.1203748630947d0*x5 + 1.64473377556
     &761d0*x6 + 0.240587129128819d0*x7
      x9 = 9.86960440108936d0*x4/x8 + 1.0d0
      x10 = log(x9)
      x11 = 0.0506605918211689d0*x10*(0.127469618870009d0*x1 + 1.0d0)
      x12 = 2.22556942115069d0*x1 + 5.98255043577108d0*x5 + 0.8004286349
     &99364d0*x6 + 0.189700432574756d0*x7
      x13 = 1.0d0 + 4.93480220054468d0*x4/x12
      x14 = log(x13)
      x15 = x14*(0.132568899905202d0*x1 + 1.0d0)
      x16 = 0.101321183642338d0*x15
      x17 = x11 - x16
      x18 = rho_s*x0
      x19 = x18 + 1.0d0
      x20 = -x18 + 1.0d0
      x21 = 1.92366105093154d0*x19**1.33333333333333d0 + 1.9236610509315
     &4d0*x20** 1.33333333333333d0 - 3.84732210186307d0
      x22 = rho_c**(-4)
      x23 = rho_s**4
      x24 = x22*x23
      x25 = x21*x24
      x26 = x25*x3
      x27 = 2.0d0*x26
      x28 = x17*x27
      x29 = 0.202642367284676d0*x15*x3
      x30 = -x24 + 1.0d0
      x31 = 2.24759186357762d0*x1 + 8.15741470348764d0*x5 + 0.4300972471
     &27664d0*x6 + 0.191151259512734d0*x7
      x32 = 1.0d0 + 29.6085746432167d0/x31
      x33 = log(x32)
      x34 = x30*x33
      x35 = x21*x34
      x36 = 0.0197517897025652d0*x35*(0.0690139921125583d0*x1 + 1.0d0)
      x37 = -x29 + x36
      x38 = -x28 + x37
      x39 = mu**2
      x40 = x39*x7
      x41 = 0.237113969739022d0*x40 + 1.0d0
      x42 = x41**(-4)
      x43 = mu**5
      x44 = rho_c**(-2)
      x45 = rho_s**2*x44
      x46 = -x45 + 1.0d0
      x47 = exp(-0.46675853320810856d0*x1)
      x48 = x0**1.0d0
      x49 = 0.00304892328706429d0*x48
      x50 = -x49
      x51 = x0**1.33333333333333d0
      x52 = 0.0128480790170175d0*x1 + x50 + 0.000275310793616358d0*x51 +
     & 0.0315297404574782d0*x7 + 1.0d0
      x53 = x47*x52
      x54 = x46*x53
      x55 = x43*x54
      x56 = x0**1.66666666666667d0
      x57 = 0.00264117541143784d0*x56
      x58 = x55*x57
      x59 = x19**0.666666666666667d0 + x20**0.666666666666667d0
      x60 = x59**3
      x61 = 0.0253302959105844d0*x60
      x62 = x59**2
      x63 = 1d0/x62
      x64 = x39*x63
      x65 = x1*x64
      x66 = 1d0/x59
      x67 = mu*x66
      x68 = 9.20897059521559d0*x5*x67 + 1.0d0
      x69 = 1d0/(8.55713948548596d0*x65 + x68)
      x70 = x2 - 1.0d0
      x71 = 1.22772285078429d0/x70
      x72 = x71 - 3.44851d0
      x73 = 2.4814019635976d0*x65
      x74 = mu**3
      x75 = x6/x60
      x76 = x68 + 15.3125681938318d0*x74*x75
      x77 = -x72*x73 + x76
      x78 = x69*x77
      x79 = x3*log(x78)
      x80 = x61*x79
      x81 = x46*x48
      x82 = 0.0445553741830535d0*x53
      x83 = x81*x82
      x84 = -3.07073492995203d0*x1 + 0.384834731559127d0*x7
      x85 = exp(-0.19230865217881402d0*x1)
      x86 = x46*x84*x85
      x87 = 1d0/x48
      x88 = 37.6991118430775d0*x87
      x89 = x86*x88
      x90 = 1d0/x7
      x91 = x19**2
      x92 = 1d0/x19
      x93 = x92**0.333333333333333d0
      x94 = x1*x93
      x95 = x92**0.666666666666667d0
      x96 = 0.0244354823084343d0*x7
      x97 = 0.337569861992026d0*x94 + x95*x96 + 1.0d0
      x98 = 1d0/x97
      x99 = x91*x98
      x100 = (-0.0177069812999059d0*x94 + 1.0d0)/x95
      x101 = 8.61370200056457d0*x100
      x102 = x101*x99
      x103 = x20**2
      x104 = 1d0/x20
      x105 = x104**0.333333333333333d0
      x106 = x1*x105
      x107 = x104**0.666666666666667d0
      x108 = 0.337569861992026d0*x106 + x107*x96 + 1.0d0
      x109 = 1d0/x108
      x110 = x103*x109
      x111 = (-0.0177069812999059d0*x106 + 1.0d0)/x107
      x112 = 8.61370200056457d0*x111
      x113 = x110*x112
      x114 = x56*(x102*x90 + x113*x90 + x89)
      x115 = 0.000132058770571892d0*x114
      x116 = x74*(x115 + x83)
      x117 = mu**8*x0**2.66666666666667d0
      x118 = 0.00316102962473761d0*x117
      x119 = x118*x38
      x120 = mu**6
      x121 = x53 - 1.0d0
      x122 = x121*x46
      x123 = x122*x56
      x124 = 0.00248266686171384d0*x123
      x125 = x0**2.0d0
      x126 = -x11 + x16
      x127 = x126*x27
      x128 = x127 + x37
      x129 = x125*x128
      x130 = 0.0533250677421794d0*x129
      x131 = x120*(x124 - x130)
      x132 = mu**4
      x133 = x122*x48
      x134 = -0.0418814102677522d0*x133
      x135 = x19**2.66666666666667d0
      x136 = 8.61370200056457d0*x135
      x137 = x20**2.66666666666667d0
      x138 = 8.61370200056457d0*x137
      x139 = exp(-0.33933171852197186d0*x1)
      x140 = x139*(-0.240695990468967d0*x1 + 0.26014827853397d0*x7)
      x141 = x140*x46
      x142 = 23.3866625383245d0*x141
      x143 = x0**1.0d0
      x144 = 0.00020688890514282d0*x143
      x145 = 0.337338207872387d0*x51
      x146 = x132*(x134 - x144*(x102 + x113 - x136 - x138 + x142) + x145
     &*x38)
      x147 = -x116 + x119 - x131 + x146 - x58 - x80
      x148 = x147*x42
      x149 = -x148 + x38
      x150 = x149/x38
      x151 = x150**2.78d0
      x152 = 0.544841407721756d0*rho_c**(-2.33333333333333d0)
      x153 = x151*x152
      x154 = x19**0.666666666666667d0 + x20**0.666666666666667d0
      x155 = x154**2
      x156 = gamma_cc/x155
      x157 = x45 - 1.0d0
      x158 = x157*x53
      x159 = x158*x43
      x160 = -x71 + 3.44851d0
      x161 = x160*x73 + x76
      x162 = x161*x69
      x163 = x70*log(x162)
      x164 = x157*x48
      x165 = x164*x82
      x166 = -x102 - x113 + x136 + x138 - x142
      x167 = x118*x128 + x120*(-x124 + x130) + x132*(x128*x145 + x134 + 
     &x144*x166) + x159*x57 + x163*x61 + x74*(-x115 + x165)
      x168 = x167*x42
      x169 = x128 - x168
      x170 = x154**3
      x171 = 1d0/x170
      x172 = 257.309189154418d0*x171
      x173 = exp(-x169*x172)
      x174 = x173 - 1.0d0
      x175 = 1d0/x174
      x176 = x156*x175
      x177 = x153*x176 + 1.0d0
      x178 = rho_c**(-4.66666666666667d0)
      x179 = x150**5.56d0
      x180 = x178*x179
      x181 = 0.296852159568225d0/x174**2
      x182 = gamma_cc**2/x154**4
      x183 = 1d0/(x177 + x180*x181*x182)
      x184 = x177*x183
      x185 = x153*x156*x184
      x186 = x185 + 1.0d0
      x187 = log(x186)
      x188 = 0.003886375d0*x170
      x189 = x29 - x36
      x190 = x189 + x28
      x191 = x148 + x190
      x192 = x116 - x119 + x131 - x146 + x58 + x80
      x193 = x1*x3
      x194 = x14*x193
      x195 = 0.0179093838034102d0*x194
      x196 = rho_c**(-5)
      x197 = x196*x23
      x198 = x197*x21
      x199 = 0.00136314985874174d0*x1 + 0.0197517897025652d0
      x200 = x199*x33
      x201 = 8.0d0*x200
      x202 = x1*x35
      x203 = 0.000908766572494496d0*x202
      x204 = rho_s*x44
      x205 = x20**0.333333333333333d0
      x206 = x19**0.333333333333333d0
      x207 = x205 - x206
      x208 = x199*x34
      x209 = x207*x208
      x210 = 5.1297628024841d0*x209
      x211 = (0.265137799810404d0*x1 + 2.0d0)*(0.741856473716896d0*x1 + 
     &0.99709173929518d0*x5 + 0.400214317499682d0*x6 + 0.126466955049837
     &d0*x7)/(x12**2*x13)
      x212 = x199*x21*x30*(0.749197287859205d0*x1 + 1.35956911724794d0*x
     &5 + 0.215048623563832d0*x6 + 0.127434173008489d0*x7)/(x31**2*x32)
      x213 = 59.2171492864334d0*x212
      x214 = rho_s**5
      x215 = x207*x3
      x216 = 10.2595256049682d0*x214*x215
      x217 = x17*x3
      x218 = 0.00861024844156464d0*x10*x193
      x219 = (0.254939237740017d0*x1 + 2.0d0)*(1.28158207914907d0*x1 + 1
     &.85339581051578d0*x5 + 0.822366887783805d0*x6 + 0.160391419419213d
     &0*x7)/(x8**2*x9)
      x220 = x195 - x211
      x221 = x40/x41**5
      x222 = 1.26460783860812d0*x221
      x223 = x45*x53
      x224 = x43*x56
      x225 = 0.0105647016457513d0*x223*x224
      x226 = 0.0042826930056725d0*x1
      x227 = 0.0210198269716521d0*x7
      x228 = 0.000367081058155144d0*x51
      x229 = x226 + x227 + x228 + x50
      x230 = 0.00528235082287567d0*x224
      x231 = x230*x47
      x232 = x229*x231*x46
      x233 = 0.00880391803812612d0*x56
      x234 = x233*x55
      x235 = 0.000821860773992031d0*x0**2.0d0
      x236 = x235*x55
      x237 = x19**(-0.333333333333333d0)
      x238 = x20**(-0.333333333333333d0)
      x239 = x237 - x238
      x240 = x18*x239
      x241 = 0.101321183642338d0*x62
      x242 = x240*x241*x79
      x243 = 0.0168588246652672d0*x117
      x244 = x243*x38
      x245 = x1*x67
      x246 = 2.85237982849532d0*x245
      x247 = mu*x1*x240*x63
      x248 = 11.4095193139813d0*x247
      x249 = 1.53482843253593d0*x5
      x250 = 6.13931373014372d0*x5
      x251 = x240*x250*x66
      x252 = -x249 + x251
      x253 = 0.827133987865867d0*x245
      x254 = 3.30853595146347d0*x247
      x255 = x6*x64
      x256 = 30.6251363876636d0*x240*x39*x75 + x252 - 7.6562840969159d0*
     &x255
      x257 = 0.0506605918211689d0*x62
      x258 = mu*x257
      x259 = x3/x77
      x260 = x258*x259*(x253*x72 - x254*x72 + x256 - x78*(-x246 + x248 +
     & x252))
      x261 = 16.0d0*x26
      x262 = x196*x216
      x263 = -x195 + x211
      x264 = x25*(x218 - x219 + x263)
      x265 = x201*x25
      x266 = x18*x210
      x267 = -x203 + x213 + x220 + x264 + x265 + x266
      x268 = x17*x261 - x17*x262 + x267
      x269 = x118*x268
      x270 = 75.398223686155d0*x84*x85*x87
      x271 = x270*x45
      x272 = 2.41662179562687d0*x0**(-0.666666666666667d0)*x86
      x273 = x46*x85*x88*(1.02357830998401d0*x1 - 0.256556487706084d0*x7
     &)
      x274 = x18*x92
      x275 = 1d0/x1
      x276 = 0.0508408867489862d0*x275
      x277 = x276/x93
      x278 = x277*x99*(-x274 + 1.0d0)
      x279 = 5.74246800037638d0*x90
      x280 = x100*x279*x99
      x281 = x104*x18
      x282 = x276/x105
      x283 = x110*x282*(x281 + 1.0d0)
      x284 = x19*x98
      x285 = 22.9698720015055d0*x90
      x286 = x100*x284*x285
      x287 = x18*x286
      x288 = x110*x111*x279
      x289 = x109*x20
      x290 = x111*x285*x289
      x291 = x18*x290
      x292 = 0.112523287330675d0*x94
      x293 = 0.0162903215389562d0*x7
      x294 = x293*x95
      x295 = x292 + x294
      x296 = x101*x90/x97**2
      x297 = x296*x91*(-x274*x292 - x274*x294 + x295)
      x298 = 0.112523287330675d0*x106
      x299 = x107*x293
      x300 = x298 + x299
      x301 = x112*x90/x108**2
      x302 = x103*x301*(x281*x298 + x281*x299 + x300)
      x303 = 0.000132058770571892d0*x56
      x304 = 0.0445553741830535d0*x47
      x305 = 0.00693220036674015d0*x51
      x306 = 0.0891107483661069d0*x48
      x307 = 0.000220097950953153d0*x114 - x223*x306
      x308 = 2.0d0*x74
      x309 = x308*(x229*x304*x81 - x303*(x271 + x272 + x273 + x278 + x28
     &0 + x283 - x287 + x288 + x291 + x297 + x302 + x89) - x305*x54 + x3
     &07 + x83)
      x310 = x121*x45
      x311 = 0.00993066744685536d0*x56
      x312 = x310*x311
      x313 = 0.0082755562057128d0*x123
      x314 = x47*(0.0777930888680181d0*x1*x52 - 0.00214134650283625d0*x1
     & + 0.00152446164353214d0*x48 - 0.000183540529077572d0*x51 - 0.0105
     &099134858261d0*x7)
      x315 = x311*x314*x46
      x316 = 0.213300270968718d0*x129
      x317 = x126*x261
      x318 = x126*x262
      x319 = x267 - x317 + x318
      x320 = 0.0533250677421794d0*x125
      x321 = x319*x320
      x322 = x120*(x312 - x313 + x315 + x316 - x321)
      x323 = 0.167525641071009d0*x310*x48
      x324 = 0.0837628205355044d0*x133
      x325 = 0.167525641071009d0*x314*x81
      x326 = 0.0006896296838094d0*x143*x166
      x327 = 0.899568554326366d0*x128*x51
      x328 = x20**1.66666666666667d0
      x329 = x19**1.66666666666667d0
      x330 = 22.9698720015055d0*x90
      x331 = 46.773325076649d0*x140*x90
      x332 = -x278 - x280 - x283 + x287 - x288 - x291 - x297 - x302
      x333 = 0.00041377781028564d0*x56*(-23.3866625383245d0*x139*x46*x90
     &*( 0.0802319968229891d0*x1 - 0.173432185689313d0*x7) - 2.645278796
     &54103d0*x141*x275 - 15.591108358883d0*x141*x90 + x18* x330*(x328 -
     & x329) - x331*x45 + x332 + 5.74246800037638d0*x90*( x135 + x137))
      x334 = x145*x319
      x335 = x132*(x323 - x324 + x325 + x326 + x327 - x333 - x334)
      x336 = -x173 + 1.0d0
      x337 = x191/x190
      x338 = x178*x337**5.56d0
      x339 = x152*x337**2.78d0
      x340 = x156/x336
      x341 = -x339*x340 + 1.0d0
      x342 = 1d0/(0.296852159568225d0*x182*x338/x336**2 + x341)
      x343 = x339*x341
      x344 = x342*x343
      x345 = x156*x344 + 1.0d0
      x346 = log(x345)
      x347 = x20**(-0.333333333333333d0)
      x348 = x19**(-0.333333333333333d0)
      x349 = x347 - 1.0d0*x348
      x350 = 1d0/x345
      x351 = rho_c**(-3.33333333333333d0)
      x352 = x177*x351
      x353 = -x348
      x354 = x347 + x353
      x355 = 1d0/x154
      x356 = 0.726455210295675d0*x355
      x357 = rho_c**(-4.33333333333333d0)*rho_s*x356
      x358 = x354*x357
      x359 = 1d0/x149
      x360 = -x225
      x361 = x359*(-x147*x222 - x150*x319 + x268 - x42*(x232 + x234 - x2
     &36 + x242 - x244 - x260 + x269 + x309 - x322 - x335 + x360))
      x362 = x169/x128
      x363 = x362**2.78d0
      x364 = 1.27129661801743d0*x351
      x365 = 1d0/x169
      x366 = -x226 - x227 - x228 + x49
      x367 = -x237 + x238
      x368 = x70/x161
      x369 = x42*(x118*x319 + x120*(-x312 + x313 - x315 - x316 + x321) -
     & x128*x243 + x132*(-x323 + x324 - x325 - x326 - x327 + x333 + x334
     &) + x157* x231*x366 - x159*x233 + x159*x235 + x163*x18*x241*x367 +
     & x258* x368*(-x160*x253 + x160*x254 + x162*(x246 - x248 + x249 - x
     &251) + x256) + x308*(x158*x305 + x164*x304*x366 - x165 + x303*(-x2
     &71 - x272 - x273 + x332 - x89) + x307) + x360)
      x370 = 0.757329556733241d0*x351
      x371 = -x127 + x168 + x189
      x372 = 1.0d0*x347
      x373 = 1029.23675661767d0*x18*x355
      x374 = -2639.87021452335d0*x126*x196*x214*x215 + 4116.94702647068d
     &0*x126*x26 + 325.395217550575d0*x167*x221 - 1319.93510726168d0*x18
     &*x209 - 4.60824902471074d0*x194 - 2058.47351323534d0*x200*x25 + 0.
     &233833989899198d0*x202 + 257.309189154418d0*x211 - 15237.116666928
     &3d0*x212 - 257.309189154418d0*x264 + 257.309189154418d0*x369
      x375 = 0.272420703860878d0*x351
      x376 = exp(x172*x371)
      x377 = x376 - 1.0d0
      x378 = 1d0/x377
      x379 = x171*x378
      x380 = x376*x379
      x381 = x156*x378
      x382 = rho_c**(-5.66666666666667d0)
      x383 = x179*x382
      x384 = x176*x383
      x385 = x171*x175
      x386 = 0.7916057588486d0*gamma_cc*x354
      x387 = x353 + x372
      x388 = x169*x373*x387 + x374
      x389 = x173*x385
      x390 = gamma_cc/x154**5
      x391 = gamma_cc*x183
      x392 = 4.0d0*x200
      x393 = rho_s**3*x21
      x394 = x22*x393
      x395 = 2.56488140124205d0*x208
      x396 = x207*x395
      x397 = x17*x215
      x398 = 8.0d0*x217
      x399 = x18*x53
      x400 = x230*x399
      x401 = x239*x257
      x402 = x393/rho_c**3
      x403 = x392*x402
      x404 = -x403
      x405 = 5.1297628024841d0*x24
      x406 = -x396 + x397*x405 - x398*x402 + x404
      x407 = 30.6251363876636d0*x255
      x408 = 3.30853595146347d0*x245
      x409 = 11.4095193139813d0*x245 + x250
      x410 = 0.0253302959105844d0*mu*x59
      x411 = x121*x18
      x412 = 0.00496533372342768d0*x411*x56
      x413 = x306*x399
      x414 = x18*x270
      x415 = x277*x284
      x416 = x282*x289
      x417 = x19*x295*x296
      x418 = x20*x300*x301
      x419 = x286 - x290 + x415 - x416 + x417 - x418
      x420 = 0.0837628205355044d0*x411*x48
      x421 = x330*(-x328 + x329)
      x422 = x18*x331
      x423 = 0.00020688890514282d0*x56
      x424 = x42*(x118*x406 + x120*(x320*x406 + x412) + x132*(x145*x406 
     &+ x420 - x423 *(x419 - x421 - x422)) - x259*x410*(x239*(-x250 - x4
     &07 + x408*x72 ) - x367*x409*x78) + x400 - x401*x79 + x74*(-x303*(-
     &x414 + x419) + x413))
      x425 = x354*x356
      x426 = x151*x352
      x427 = -x205 + x206
      x428 = x395*x427
      x429 = x126*x3
      x430 = x427*x429
      x431 = x405*x430
      x432 = x402*x429
      x433 = 8.0d0*x432
      x434 = x404 + x428 + x431 + x433
      x435 = x150*x434
      x436 = x200*x402
      x437 = x24*x397
      x438 = x217*x402
      x439 = x349*x355
      x440 = x351*x363
      x441 = x362**5.56d0*x382
      x442 = -x286 + x290 - x415 + x416 - x417 + x418
      x443 = x42*(x118*x434 + x120*(x320*x434 + x412) + x132*(x145*x434 
     &+ x420 + x423 *(x421 + x422 + x442)) + x163*x401 + x367*x368*x410*
     &(x160*x408 - x162*x409 + x250 + x407) + x400 + x74*(x303*(x414 + x
     &442) + x413 ))
      x444 = x365*(x362*(x403 - x428 - x431 - x433) + x434 - x443)
      x445 = -659.967553630838d0*x208*x427 - 1319.93510726168d0*x24*x430
     & + 514.618378308835d0*x371*x439 - 2058.47351323534d0*x432 + 1029.2
     &3675661767d0*x436 + 257.309189154418d0*x443
      x446 = 0.003886375d0*x154
      E = -rho_c*(-x187*x188 + x191)
      d1E(1) = 0.5d0*rho_c*(x0*x192*x222 + x0*x195 - x0*x203 - x0*x211 +
     & x0*x213 + x0* x42*(x225 - x232 - x234 + x236 - x242 + x244 + x260
     & - x269 - x309 + x322 + x335) - 0.00777275d0*x151*x154*x350*x391*(
     &-x152*x176* x184*(rho_c**(-6.66666666666667d0)*rho_s*x179*x385*x38
     &6 + x151* x358 - x151*x361*x370 + x151*x364 + x151*x375*x388*x389 
     &+ x173* x181*x383*x388*x390 - 0.825249003599665d0*x361*x384 + 1.38
     &531007798505d0*x384) + x152*x363*x381*(-x357*(-x347 + x348) + x364
     & - x365*x370*(-x167*x222 + x319 + x362*(x203 - x213 + x263 - x264 
     &- x265 - x266 + x317 - x318) - x369) + x375*x380*(x371*x373* (x348
     & - x372) + x374)) + x177*x358 - 0.757329556733241d0*x352* x361 + 1
     &.27129661801743d0*x352) + 0.0155455d0*x155*x204*x346*x349 + x198*x
     &201 + 16.0d0*x198*x217 - x198*(-x218 + x219 + x220) + x204*x210 - 
     &x17*x216/rho_c**6) + x188*x346 + x192*x42 + x38
      d1E(2) = -rho_c*(0.00777275d0*x0*x155*x187*x387 + x0*x396 + x0*x42
     &4 - 5.1297628024841d0*x197*x397 + x392*x394 + x394*x398 - x391*x44
     &6*( -x185*x378*(-0.59370431913645d0*x376*x390*x441*x445/x377**2 + 
     &x379*x386*x441 - 0.544841407721756d0*x380*x440*x445 + 1.6504980071
     &9933d0*x381*x441*x444 + x425*x440 + 1.51465911346648d0*x440*x444) 
     &+ 1.51465911346648d0*x359*x426*( x406 - x424 - x435) + 0.544841407
     &721756d0*x384*(x355*( 0.726455210295675d0*x347 - 0.726455210295675
     &d0*x348) + x359*( -3.88492098935195d0*x209 - 1.51465911346648d0*x4
     &24 - 1.51465911346648d0*x435 - 6.05863645386593d0*x436 + 7.7698419
     &7870389d0*x437 - 12.1172729077319d0*x438) - x389*( -280.3854016772
     &73d0*x169*x439 + 359.577650970909d0*x209 + 140.192700838636d0*x424
     & + 560.770803354546d0*x436 - 719.155301941818d0*x437 + 1121.541606
     &70909d0*x438)) + x425*x426)/ x186)
      d1E(3) = rho_c*x342*x350*x446*(-x176*x344*(x153 + 0.59370431913645
     &d0*x176*x180) - 0.296852159568225d0*x338*x340 + x343)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu
     &, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 10, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399, x400, x401
     &, x402, x403, x404, x405, x406, x407, x408, x409, x410, x411, x412
     &, x413, x414, x415, x416, x417, x418, x419, x420, x421, x422, x423
     &, x424, x425, x426, x427, x428, x429, x430, x431, x432, x433, x434
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446, x447, x448, x449, x450, x451, x452, x453, x454, x455, x456
     &, x457, x458, x459, x460, x461, x462, x463, x464, x465, x466, x467
     &, x468, x469, x470, x471, x472, x473, x474, x475, x476, x477, x478
     &, x479, x480, x481, x482, x483, x484, x485, x486, x487, x488, x489
     &, x490, x491, x492, x493, x494, x495, x496, x497, x498, x499, x500
     &, x501, x502, x503, x504, x505, x506, x507, x508, x509, x510, x511
     &, x512, x513, x514, x515, x516, x517, x518, x519, x520, x521, x522
     &, x523, x524, x525, x526, x527, x528, x529, x530, x531, x532, x533
     &, x534, x535, x536, x537, x538, x539, x540, x541, x542, x543, x544
     &, x545, x546, x547, x548, x549, x550, x551, x552, x553, x554, x555
     &, x556, x557, x558, x559, x560, x561, x562, x563, x564, x565, x566
     &, x567, x568, x569, x570, x571, x572, x573, x574, x575, x576, x577
     &, x578, x579, x580, x581, x582, x583, x584, x585, x586, x587, x588
     &, x589, x590, x591, x592, x593, x594, x595, x596, x597, x598, x599
     &, x600, x601, x602, x603, x604, x605, x606, x607, x608, x609, x610
     &, x611, x612, x613, x614, x615, x616, x617, x618, x619, x620, x621
     &, x622, x623, x624, x625, x626, x627, x628, x629, x630, x631, x632
     &, x633, x634, x635, x636, x637, x638, x639, x640, x641, x642, x643
     &, x644, x645, x646, x647, x648, x649, x650, x651, x652, x653, x654
     &, x655, x656, x657, x658, x659, x660, x661, x662, x663, x664, x665
     &, x666, x667, x668, x669, x670, x671, x672, x673, x674, x675, x676
     &, x677, x678, x679, x680, x681, x682, x683, x684, x685, x686, x687
     &, x688, x689, x690, x691, x692, x693, x694, x695, x696, x697, x698
     &, x699, x700, x701, x702, x703, x704, x705, x706, x707, x708, x709
     &, x710, x711, x712, x713, x714, x715, x716, x717, x718, x719, x720
     &, x721, x722, x723, x724, x725, x726, x727, x728, x729, x730, x731
     &, x732, x733, x734, x735, x736, x737, x738, x739, x740, x741, x742
     &, x743, x744, x745, x746, x747, x748, x749, x750, x751, x752, x753
     &, x754, x755, x756, x757, x758, x759, x760, x761, x762, x763, x764
     &, x765, x766, x767, x768, x769, x770, x771, x772, x773, x774, x775
     &, x776, x777, x778, x779, x780, x781, x782, x783, x784, x785, x786
     &, x787, x788, x789, x790, x791, x792, x793, x794, x795, x796, x797
     &, x798, x799, x800, x801, x802, x803, x804, x805, x806, x807, x808
     &, x809, x810, x811, x812, x813, x814, x815, x816, x817, x818, x819
     &, x820, x821, x822, x823, x824, x825, x826, x827, x828, x829, x830
     &, x831, x832, x833, x834, x835, x836, x837, x838, x839, x840, x841
     &, x842, x843, x844, x845, x846, x847, x848, x849, x850, x851, x852
     &, x853, x854, x855, x856, x857, x858, x859, x860, x861, x862, x863
     &, x864, x865, x866, x867, x868, x869, x870, x871, x872, x873, x874
     &, x875, x876, x877, x878, x879, x880, x881, x882, x883, x884, x885
     &, x886, x887, x888, x889, x890, x891, x892, x893, x894, x895, x896
     &, x897, x898, x899, x900, x901, x902, x903, x904, x905, x906, x907
     &, x908, x909, x910, x911, x912, x913, x914, x915, x916, x917, x918
     &, x919, x920, x921, x922, x923, x924, x925, x926, x927, x928, x929
     &, x930, x931, x932, x933, x934, x935, x936, x937, x938, x939, x940
     &, x941, x942, x943, x944, x945, x946, x947, x948, x949, x950, x951
     &, x952, x953, x954, x955, x956, x957, x958, x959, x960, x961, x962
     &, x963, x964, x965, x966, x967, x968, x969, x970, x971, x972, x973
     &, x974, x975, x976, x977, x978, x979, x980, x981, x982, x983, x984
     &, x985, x986, x987, x988, x989, x990, x991, x992, x993, x994, x995
     &, x996, x997, x998, x999, x1000, x1001, x1002, x1003, x1004, x1005
     &, x1006, x1007, x1008, x1009, x1010, x1011, x1012, x1013, x1014, x
     &1015, x1016, x1017, x1018, x1019, x1020, x1021, x1022, x1023, x102
     &4, x1025, x1026, x1027, x1028, x1029, x1030, x1031, x1032, x1033, 
     &x1034, x1035, x1036, x1037, x1038, x1039, x1040, x1041, x1042, x10
     &43, x1044, x1045, x1046, x1047, x1048, x1049, x1050, x1051, x1052,
     & x1053, x1054, x1055, x1056, x1057, x1058, x1059, x1060, x1061, x1
     &062, x1063, x1064, x1065, x1066, x1067, x1068, x1069, x1070, x1071
     &, x1072, x1073, x1074, x1075, x1076, x1077, x1078, x1079, x1080, x
     &1081, x1082, x1083, x1084, x1085, x1086, x1087, x1088, x1089, x109
     &0, x1091, x1092, x1093, x1094, x1095, x1096, x1097, x1098, x1099, 
     &x1100, x1101, x1102, x1103, x1104, x1105, x1106, x1107, x1108, x11
     &09, x1110, x1111, x1112, x1113, x1114, x1115, x1116, x1117, x1118,
     & x1119, x1120, x1121, x1122, x1123, x1124, x1125, x1126, x1127, x1
     &128, x1129, x1130, x1131, x1132, x1133, x1134, x1135, x1136, x1137
     &, x1138, x1139, x1140, x1141, x1142, x1143, x1144, x1145, x1146, x
     &1147, x1148, x1149, x1150, x1151, x1152, x1153, x1154, x1155, x115
     &6, x1157, x1158, x1159, x1160, x1161, x1162, x1163, x1164, x1165, 
     &x1166, x1167, x1168, x1169, x1170, x1171, x1172, x1173, x1174, x11
     &75, x1176, x1177, x1178, x1179, x1180, x1181, x1182, x1183, x1184,
     & x1185, x1186, x1187, x1188, x1189, x1190, x1191, x1192, x1193, x1
     &194, x1195, x1196, x1197, x1198, x1199, x1200, x1201, x1202, x1203
     &, x1204, x1205, x1206, x1207, x1208, x1209, x1210, x1211, x1212, x
     &1213, x1214, x1215, x1216, x1217, x1218, x1219, x1220, x1221, x122
     &2, x1223, x1224, x1225, x1226, x1227, x1228, x1229, x1230, x1231, 
     &x1232, x1233, x1234, x1235, x1236, x1237, x1238, x1239, x1240, x12
     &41, x1242, x1243, x1244, x1245, x1246, x1247, x1248, x1249, x1250,
     & x1251, x1252, x1253, x1254, x1255, x1256, x1257, x1258, x1259, x1
     &260, x1261, x1262, x1263, x1264, x1265, x1266, x1267, x1268, x1269
     &, x1270, x1271, x1272, x1273, x1274, x1275, x1276, x1277, x1278, x
     &1279, x1280, x1281, x1282, x1283, x1284, x1285, x1286, x1287, x128
     &8, x1289, x1290, x1291, x1292, x1293, x1294, x1295, x1296, x1297, 
     &x1298, x1299, x1300, x1301, x1302, x1303, x1304, x1305, x1306, x13
     &07, x1308, x1309, x1310, x1311, x1312, x1313, x1314, x1315, x1316,
     & x1317, x1318, x1319, x1320, x1321, x1322, x1323, x1324, x1325, x1
     &326, x1327, x1328, x1329, x1330, x1331, x1332, x1333, x1334, x1335
     &, x1336, x1337, x1338, x1339, x1340, x1341, x1342, x1343, x1344, x
     &1345, x1346, x1347, x1348, x1349, x1350, x1351, x1352, x1353, x135
     &4, x1355, x1356, x1357, x1358, x1359, x1360, x1361, x1362, x1363, 
     &x1364, x1365, x1366, x1367, x1368, x1369, x1370, x1371, x1372, x13
     &73, x1374, x1375, x1376, x1377, x1378, x1379, x1380, x1381, x1382,
     & x1383, x1384, x1385, x1386, x1387, x1388, x1389, x1390, x1391, x1
     &392, x1393, x1394, x1395, x1396, x1397, x1398, x1399, x1400, x1401
     &, x1402, x1403, x1404, x1405, x1406, x1407, x1408, x1409, x1410, x
     &1411, x1412, x1413, x1414, x1415, x1416, x1417, x1418, x1419, x142
     &0, x1421, x1422, x1423, x1424, x1425, x1426, x1427, x1428, x1429, 
     &x1430, x1431, x1432, x1433, x1434, x1435, x1436, x1437, x1438, x14
     &39, x1440, x1441, x1442, x1443, x1444, x1445, x1446, x1447, x1448,
     & x1449, x1450, x1451, x1452, x1453, x1454, x1455, x1456, x1457, x1
     &458, x1459, x1460, x1461, x1462, x1463, x1464, x1465, x1466, x1467
     &, x1468, x1469, x1470, x1471, x1472, x1473, x1474, x1475, x1476, x
     &1477, x1478, x1479, x1480, x1481, x1482, x1483, x1484, x1485, x148
     &6, x1487, x1488, x1489, x1490, x1491, x1492, x1493, x1494, x1495, 
     &x1496, x1497, x1498, x1499, x1500, x1501, x1502, x1503, x1504, x15
     &05, x1506, x1507, x1508, x1509, x1510, x1511, x1512, x1513, x1514,
     & x1515, x1516, x1517, x1518, x1519, x1520, x1521, x1522, x1523, x1
     &524, x1525, x1526, x1527, x1528, x1529
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = rho_s*x0
      x2 = x1 + 1.0d0
      x3 = -x1 + 1.0d0
      x4 = x2**0.666666666666667d0 + x3**0.666666666666667d0
      x5 = x4**2
      x6 = 1d0/x5
      x7 = gamma_cc*x6
      x8 = gamma_cc**2
      x9 = x4**(-4)
      x10 = x8*x9
      x11 = 0.296852159568225d0*x10
      x12 = x0**0.666666666666667d0
      x13 = mu**2
      x14 = x12*x13
      x15 = 0.237113969739022d0*x14 + 1.0d0
      x16 = x15**(-4)
      x17 = rho_c**(-2)
      x18 = rho_s**2
      x19 = x17*x18
      x20 = x19 - 1.0d0
      x21 = x0**0.333333333333333d0
      x22 = exp(-0.46675853320810856d0*x21)
      x23 = x0**1.0d0
      x24 = 0.00304892328706429d0*x23
      x25 = -x24
      x26 = x0**1.33333333333333d0
      x27 = 0.0315297404574782d0*x12 + 0.0128480790170175d0*x21 + x25 + 
     &0.000275310793616358d0*x26 + 1.0d0
      x28 = x22*x27
      x29 = x20*x28
      x30 = x0**1.66666666666667d0
      x31 = mu**5
      x32 = x30*x31
      x33 = 0.00264117541143784d0*x32
      x34 = x2**0.666666666666667d0
      x35 = x3**0.666666666666667d0
      x36 = x34 + x35
      x37 = x36**3
      x38 = 0.0253302959105844d0*x37
      x39 = 0.693147180559945d0
      x40 = x39 - 1.0d0
      x41 = x36**2
      x42 = 1d0/x41
      x43 = x13*x42
      x44 = x21*x43
      x45 = x0**0.166666666666667d0
      x46 = 1d0/x36
      x47 = mu*x46
      x48 = 9.20897059521559d0*x45*x47 + 1.0d0
      x49 = 8.55713948548596d0*x44 + x48
      x50 = 1d0/x49
      x51 = 1.22772285078429d0/x40
      x52 = -x51 + 3.44851d0
      x53 = 2.4814019635976d0*x44
      x54 = mu**3
      x55 = sqrt(x0)
      x56 = 1d0/x37
      x57 = x55*x56
      x58 = x48 + 15.3125681938318d0*x54*x57
      x59 = x52*x53 + x58
      x60 = x50*x59
      x61 = x40*log(x60)
      x62 = x23*x28
      x63 = 0.0445553741830535d0*x62
      x64 = x20*x63
      x65 = 0.384834731559127d0*x12 - 3.07073492995203d0*x21
      x66 = exp(-0.19230865217881402d0*x21)
      x67 = x66/x23
      x68 = x65*x67
      x69 = -x19 + 1.0d0
      x70 = 37.6991118430775d0*x69
      x71 = x68*x70
      x72 = 1d0/x12
      x73 = x2**2
      x74 = 8.61370200056457d0*x73
      x75 = 1d0/x2
      x76 = x75**0.333333333333333d0
      x77 = x21*x76
      x78 = 0.0177069812999059d0*x77
      x79 = -x78 + 1.0d0
      x80 = x75**0.666666666666667d0
      x81 = 0.0244354823084343d0*x12
      x82 = 0.337569861992026d0*x77 + x80*x81 + 1.0d0
      x83 = 1d0/x82
      x84 = 1d0/x80
      x85 = x83*x84
      x86 = x79*x85
      x87 = x74*x86
      x88 = x3**2
      x89 = 8.61370200056457d0*x88
      x90 = 1d0/x3
      x91 = x90**0.333333333333333d0
      x92 = x21*x91
      x93 = 0.0177069812999059d0*x92
      x94 = -x93 + 1.0d0
      x95 = x90**0.666666666666667d0
      x96 = x81*x95 + 0.337569861992026d0*x92 + 1.0d0
      x97 = 1d0/x96
      x98 = 1d0/x95
      x99 = x97*x98
      x100 = x94*x99
      x101 = x100*x89
      x102 = x30*(x101*x72 + x71 + x72*x87)
      x103 = 0.000132058770571892d0*x102
      x104 = 1.64473377556761d0*x55
      x105 = x104 + 0.240587129128819d0*x12 + 3.84474623744721d0*x21 + 1
     &1.1203748630947d0*x45
      x106 = -x39 + 1.0d0
      x107 = 1d0/x106
      x108 = 9.86960440108936d0*x107
      x109 = 1.0d0 + x108/x105
      x110 = log(x109)
      x111 = 0.0506605918211689d0*x110*(0.127469618870009d0*x21 + 1.0d0)
     &
      x112 = 0.800428634999364d0*x55
      x113 = x112 + 0.189700432574756d0*x12 + 2.22556942115069d0*x21 + 5
     &.98255043577108d0*x45
      x114 = 4.93480220054468d0*x107/x113 + 1.0d0
      x115 = log(x114)
      x116 = x115*(0.132568899905202d0*x21 + 1.0d0)
      x117 = 0.101321183642338d0*x116
      x118 = -x111 + x117
      x119 = 1.92366105093154d0*x2**1.33333333333333d0 + 1.9236610509315
     &4d0*x3** 1.33333333333333d0 - 3.84732210186307d0
      x120 = rho_c**(-4)
      x121 = rho_s**4
      x122 = x120*x121
      x123 = x119*x122
      x124 = x106*x123
      x125 = 2.0d0*x124
      x126 = x118*x125
      x127 = 0.202642367284676d0*x106*x116
      x128 = -x122 + 1.0d0
      x129 = 0.430097247127664d0*x55
      x130 = 0.191151259512734d0*x12 + x129 + 2.24759186357762d0*x21 + 8
     &.15741470348764d0*x45
      x131 = 1.0d0 + 29.6085746432167d0/x130
      x132 = log(x131)
      x133 = x128*x132
      x134 = x119*x133
      x135 = 0.0197517897025652d0*x134*(0.0690139921125583d0*x21 + 1.0d0
     &)
      x136 = -x127 + x135
      x137 = x126 + x136
      x138 = mu**8*x0**2.66666666666667d0
      x139 = 0.00316102962473761d0*x138
      x140 = mu**6
      x141 = x28 - 1.0d0
      x142 = x141*x30
      x143 = x142*x69
      x144 = 0.00248266686171384d0*x143
      x145 = x0**2.0d0
      x146 = x137*x145
      x147 = 0.0533250677421794d0*x146
      x148 = mu**4
      x149 = x141*x23
      x150 = x149*x69
      x151 = -0.0418814102677522d0*x150
      x152 = x2**2.66666666666667d0
      x153 = 8.61370200056457d0*x152
      x154 = x3**2.66666666666667d0
      x155 = 8.61370200056457d0*x154
      x156 = 0.240695990468967d0*x21
      x157 = 0.26014827853397d0*x12
      x158 = exp(-0.33933171852197186d0*x21)
      x159 = x158*(-x156 + x157)
      x160 = x159*x69
      x161 = 23.3866625383245d0*x160
      x162 = -x101 + x153 + x155 - x161 - x87
      x163 = x0**1.0d0
      x164 = 0.00020688890514282d0*x163
      x165 = 0.337338207872387d0*x26
      x166 = x137*x139 + x140*(-x144 + x147) + x148*(x137*x165 + x151 + 
     &x162*x164) + x29*x33 + x38*x61 + x54*(-x103 + x64)
      x167 = x16*x166
      x168 = x137 - x167
      x169 = x4**3
      x170 = 1d0/x169
      x171 = 257.309189154418d0*x170
      x172 = exp(-x168*x171)
      x173 = x172 - 1.0d0
      x174 = x173**(-2)
      x175 = rho_c**(-4.66666666666667d0)
      x176 = x111 - x117
      x177 = x125*x176
      x178 = x136 - x177
      x179 = 1d0/x178
      x180 = x28*x69
      x181 = x180*x33
      x182 = x51 - 3.44851d0
      x183 = -x182*x53 + x58
      x184 = x183*x50
      x185 = x106*log(x184)
      x186 = x185*x38
      x187 = x63*x69
      x188 = x54*(x103 + x187)
      x189 = x139*x178
      x190 = x140*(x144 - x147)
      x191 = x101 - x153 - x155 + x161 + x87
      x192 = x148*(x151 - x164*x191 + x165*x178)
      x193 = -x181 - x186 - x188 + x189 - x190 + x192
      x194 = x16*x193
      x195 = x178 - x194
      x196 = x179*x195
      x197 = x196**5.56d0
      x198 = x175*x197
      x199 = x174*x198
      x200 = x196**2.78d0
      x201 = rho_c**(-2.33333333333333d0)
      x202 = 0.544841407721756d0*x201
      x203 = x200*x202
      x204 = 1d0/x173
      x205 = x204*x7
      x206 = x203*x205 + 1.0d0
      x207 = x11*x199 + x206
      x208 = 1d0/x207
      x209 = x203*x206
      x210 = x208*x209
      x211 = x210*x7 + 1.0d0
      x212 = log(x211)
      x213 = 0.003886375d0*x169
      x214 = x127 - x135
      x215 = x177 + x214
      x216 = x194 + x215
      x217 = x181 + x186 + x188 - x189 + x190 - x192
      x218 = x106*x21
      x219 = x115*x218
      x220 = 0.0179093838034102d0*x219
      x221 = rho_c**(-5)
      x222 = x119*x121
      x223 = x221*x222
      x224 = 0.00136314985874174d0*x21 + 0.0197517897025652d0
      x225 = x132*x224
      x226 = 8.0d0*x225
      x227 = x134*x21
      x228 = 0.000908766572494496d0*x227
      x229 = rho_s*x17
      x230 = x3**0.333333333333333d0
      x231 = x2**0.333333333333333d0
      x232 = x230 - x231
      x233 = x128*x225
      x234 = x232*x233
      x235 = 5.1297628024841d0*x234
      x236 = 0.126466955049837d0*x12 + 0.741856473716896d0*x21 + 0.99709
     &173929518d0* x45 + 0.400214317499682d0*x55
      x237 = 1d0/x114
      x238 = x113**(-2)
      x239 = 0.265137799810404d0*x21 + 2.0d0
      x240 = x237*x238*x239
      x241 = x236*x240
      x242 = x119*x128*x224
      x243 = 1d0/x131
      x244 = x130**(-2)
      x245 = 0.127434173008489d0*x12 + 0.749197287859205d0*x21 + 1.35956
     &911724794d0* x45 + 0.215048623563832d0*x55
      x246 = x243*x244*x245
      x247 = x242*x246
      x248 = 59.2171492864334d0*x247
      x249 = rho_c**(-6)
      x250 = rho_s**5
      x251 = x106*x232
      x252 = 10.2595256049682d0*x251
      x253 = x176*x252
      x254 = x106*x176
      x255 = x110*x218
      x256 = 0.00861024844156464d0*x255
      x257 = 0.160391419419213d0*x12 + 1.28158207914907d0*x21 + 1.853395
     &81051578d0* x45 + 0.822366887783805d0*x55
      x258 = 1d0/x109
      x259 = x105**(-2)
      x260 = 0.254939237740017d0*x21 + 2.0d0
      x261 = x258*x259*x260
      x262 = x257*x261
      x263 = x220 - x241
      x264 = x14/x15**5
      x265 = 1.26460783860812d0*x264
      x266 = x28*x32
      x267 = x19*x266
      x268 = 0.0105647016457513d0*x267
      x269 = 0.0042826930056725d0*x21
      x270 = 0.0210198269716521d0*x12
      x271 = 0.000367081058155144d0*x26
      x272 = x25 + x269 + x270 + x271
      x273 = x272*x69
      x274 = x22*x30
      x275 = x274*x31
      x276 = 0.00528235082287567d0*x275
      x277 = x273*x276
      x278 = 0.00880391803812612d0*x32
      x279 = x180*x278
      x280 = x0**2.0d0*x31
      x281 = 0.000821860773992031d0*x280
      x282 = x180*x281
      x283 = x2**(-0.333333333333333d0)
      x284 = x3**(-0.333333333333333d0)
      x285 = x283 - x284
      x286 = 0.101321183642338d0*x285
      x287 = x185*x41
      x288 = x1*x287
      x289 = x286*x288
      x290 = 0.0168588246652672d0*x138
      x291 = x178*x290
      x292 = mu*x41
      x293 = 0.0506605918211689d0*x292
      x294 = 1d0/x183
      x295 = x21*x47
      x296 = 2.85237982849532d0*x295
      x297 = mu*x21
      x298 = x297*x42
      x299 = x1*x285
      x300 = x298*x299
      x301 = 11.4095193139813d0*x300
      x302 = 1.53482843253593d0*x45
      x303 = 6.13931373014372d0*x45
      x304 = x1*x46
      x305 = x285*x304
      x306 = x303*x305
      x307 = -x302 + x306
      x308 = -x296 + x301 + x307
      x309 = 0.827133987865867d0*x295
      x310 = x182*x298
      x311 = x299*x310
      x312 = x43*x55
      x313 = 7.6562840969159d0*x312
      x314 = x13*x57
      x315 = x299*x314
      x316 = 30.6251363876636d0*x315
      x317 = x307 - x313 + x316
      x318 = x182*x309 - 3.30853595146347d0*x311 + x317
      x319 = x106*(-x184*x308 + x318)
      x320 = x294*x319
      x321 = x293*x320
      x322 = 16.0d0*x124
      x323 = x221*x250
      x324 = -x220 + x241
      x325 = x256 - x262 + x324
      x326 = x123*x325
      x327 = x123*x226
      x328 = x1*x235
      x329 = -x228 + x248 + x263 + x326 + x327 + x328
      x330 = x176*x322 - x253*x323 + x329
      x331 = x139*x330
      x332 = -0.256556487706084d0*x12 + 1.02357830998401d0*x21
      x333 = x67*x70
      x334 = x332*x333
      x335 = x0**(-0.666666666666667d0)
      x336 = x65*x66*x69
      x337 = x335*x336
      x338 = 2.41662179562687d0*x337
      x339 = 75.398223686155d0*x68
      x340 = x19*x339
      x341 = 5.74246800037638d0*x72
      x342 = x73*x86
      x343 = x341*x342
      x344 = x100*x88
      x345 = x341*x344
      x346 = 1d0/x21
      x347 = 0.0508408867489862d0*x346
      x348 = 1d0/x76
      x349 = x348*x83
      x350 = x1*x75
      x351 = -x350 + 1.0d0
      x352 = x351*x73
      x353 = x349*x352
      x354 = x347*x353
      x355 = 1d0/x91
      x356 = x355*x97
      x357 = x1*x90
      x358 = x357 + 1.0d0
      x359 = x358*x88
      x360 = x356*x359
      x361 = x347*x360
      x362 = 0.112523287330675d0*x77
      x363 = 0.0162903215389562d0*x12
      x364 = x363*x80
      x365 = x362 + x364
      x366 = -x350*x362 - x350*x364 + x365
      x367 = x82**(-2)
      x368 = x367*x84
      x369 = x368*x79
      x370 = x366*x369
      x371 = x370*x72
      x372 = x371*x74
      x373 = 0.112523287330675d0*x92
      x374 = x363*x95
      x375 = x373 + x374
      x376 = x357*x373 + x357*x374 + x375
      x377 = x96**(-2)
      x378 = x377*x98
      x379 = x378*x94
      x380 = x376*x379*x72
      x381 = x380*x89
      x382 = 22.9698720015055d0*x72
      x383 = x100*x3
      x384 = x382*x383
      x385 = x1*x384
      x386 = x2*x86
      x387 = x382*x386
      x388 = x1*x387
      x389 = x343 + x345 + x354 + x361 + x372 + x381 + x385 - x388
      x390 = x334 + x338 + x340 + x389 + x71
      x391 = 0.000132058770571892d0*x30
      x392 = x22*x23
      x393 = 0.0445553741830535d0*x392
      x394 = 0.00693220036674015d0*x26
      x395 = 0.0891107483661069d0*x62
      x396 = 0.000220097950953153d0*x102 - x19*x395
      x397 = 2.0d0*x54
      x398 = x397*(-x180*x394 + x187 + x273*x393 - x390*x391 + x396)
      x399 = x142*x19
      x400 = 0.00993066744685536d0*x399
      x401 = 0.0082755562057128d0*x143
      x402 = x21*x27
      x403 = 0.00214134650283625d0*x21
      x404 = 0.0105099134858261d0*x12
      x405 = 0.00152446164353214d0*x23
      x406 = 0.000183540529077572d0*x26
      x407 = -x403 - x404 + x405 - x406
      x408 = 0.0777930888680181d0*x402 + x407
      x409 = x408*x69
      x410 = x274*x409
      x411 = 0.00993066744685536d0*x410
      x412 = 0.213300270968718d0*x146
      x413 = x118*x322
      x414 = x118*x323
      x415 = x252*x414
      x416 = x329 - x413 + x415
      x417 = 0.0533250677421794d0*x145
      x418 = x416*x417
      x419 = x140*(x400 - x401 + x411 + x412 - x418)
      x420 = x149*x19
      x421 = 0.167525641071009d0*x420
      x422 = 0.0837628205355044d0*x149
      x423 = x422*x69
      x424 = x392*x409
      x425 = 0.167525641071009d0*x424
      x426 = x162*x163
      x427 = 0.0006896296838094d0*x426
      x428 = x137*x26
      x429 = 0.899568554326366d0*x428
      x430 = x72*(x152 + x154)
      x431 = 5.74246800037638d0*x430
      x432 = x3**1.66666666666667d0
      x433 = x2**1.66666666666667d0
      x434 = x72*(x432 - x433)
      x435 = x1*x434
      x436 = 22.9698720015055d0*x435
      x437 = 2.64527879654103d0*x160*x346
      x438 = x160*x72
      x439 = 15.591108358883d0*x438
      x440 = x158*x72
      x441 = -0.173432185689313d0*x12 + 0.0802319968229891d0*x21
      x442 = x441*x69
      x443 = x440*x442
      x444 = 23.3866625383245d0*x443
      x445 = x159*x72
      x446 = 46.773325076649d0*x445
      x447 = x19*x446
      x448 = -x343 - x345 - x354 - x361 - x372 - x381 - x385 + x388
      x449 = x30*(x431 + x436 - x437 - x439 - x444 - x447 + x448)
      x450 = 0.00041377781028564d0*x449
      x451 = x165*x416
      x452 = x148*(x421 - x423 + x425 + x427 + x429 - x450 - x451)
      x453 = -x172 + 1.0d0
      x454 = x216/x215
      x455 = x175*x454**5.56d0
      x456 = x202*x454**2.78d0
      x457 = x7/x453
      x458 = -x456*x457 + 1.0d0
      x459 = 1d0/(x11*x455/x453**2 + x458)
      x460 = x456*x458
      x461 = x459*x460
      x462 = x461*x7 + 1.0d0
      x463 = log(x462)
      x464 = x3**(-0.333333333333333d0)
      x465 = x2**(-0.333333333333333d0)
      x466 = x464 - 1.0d0*x465
      x467 = x229*x466
      x468 = gamma_cc*x4
      x469 = x200*x208
      x470 = 1d0/x462
      x471 = rho_c**(-3.33333333333333d0)
      x472 = 1.27129661801743d0*x471
      x473 = x206*x472
      x474 = rho_c**(-4.33333333333333d0)
      x475 = x206*x474
      x476 = 0.726455210295675d0*rho_s
      x477 = 1d0/x4
      x478 = -x465
      x479 = x464 + x478
      x480 = x477*x479
      x481 = x476*x480
      x482 = x475*x481
      x483 = 0.757329556733241d0*x471
      x484 = 1d0/x195
      x485 = -x268
      x486 = x277 + x279 - x282 + x289 - x291 - x321 + x331 + x398 - x41
     &9 - x452 + x485
      x487 = -x16*x486 - x193*x265 + x330
      x488 = -x196*x416 + x487
      x489 = x484*x488
      x490 = x206*x489
      x491 = x483*x490
      x492 = 1.27129661801743d0*x471
      x493 = x474*x477
      x494 = -x464 + x465
      x495 = x476*x494
      x496 = x493*x495
      x497 = 1d0/x168
      x498 = x228 - x248 + x324 - x326 - x327 - x328 + x413 - x415
      x499 = 1d0/x137
      x500 = x168*x499
      x501 = x24 - x269 - x270 - x271
      x502 = x20*x501
      x503 = -x283 + x284
      x504 = x41*x61
      x505 = x1*x504
      x506 = x302 - x306
      x507 = x296 - x301 + x506
      x508 = x309*x52
      x509 = x300*x52
      x510 = 3.30853595146347d0*x509
      x511 = x317 - x508 + x510
      x512 = x507*x60 + x511
      x513 = x40/x59
      x514 = x512*x513
      x515 = -x340
      x516 = -x334 - x338 + x448 + x515 - x71
      x517 = x20*x393
      x518 = -x137*x290 + x139*x416 + x140*(-x400 + x401 - x411 - x412 +
     & x418) + x148 *(-x421 + x423 - x425 - x427 - x429 + x450 + x451) +
     & x276*x502 - x278*x29 + x281*x29 + x293*x514 + x397*(x29*x394 + x3
     &91*x516 + x396 + x501*x517 - x64) + x485 + 0.101321183642338d0*x50
     &3*x505
      x519 = x16*x518
      x520 = x166*x265
      x521 = x416 - x519 - x520
      x522 = x498*x500 + x521
      x523 = x497*x522
      x524 = x483*x523
      x525 = -x126 + x167 + x214
      x526 = exp(x171*x525)
      x527 = 0.272420703860878d0*x526
      x528 = 1029.23675661767d0*x477
      x529 = x525*x528
      x530 = 1.0d0*x464
      x531 = x465 - x530
      x532 = x1*x531
      x533 = x166*x264
      x534 = x118*x124
      x535 = x251*x414
      x536 = x123*x225
      x537 = -1319.93510726168d0*x1*x234 - 4.60824902471074d0*x219 + 0.2
     &33833989899198d0*x227 + 257.309189154418d0*x241 - 15237.1166669283
     &d0*x247 - 257.309189154418d0*x326 + 257.309189154418d0*x519 + 325.
     &395217550575d0*x533 + 4116.94702647068d0*x534 - 2639.87021452335d0
     &*x535 - 2058.47351323534d0*x536
      x538 = x529*x532 + x537
      x539 = x526 - 1.0d0
      x540 = 1d0/x539
      x541 = x170*x540
      x542 = x538*x541
      x543 = x527*x542
      x544 = x471*x543
      x545 = x492 - x496 - x524 + x544
      x546 = x500**2.78d0
      x547 = x202*x546
      x548 = x540*x7
      x549 = x547*x548
      x550 = x545*x549
      x551 = rho_c**(-5.66666666666667d0)
      x552 = 1.38531007798505d0*x551
      x553 = x197*x205
      x554 = 0.7916057588486d0*rho_s
      x555 = rho_c**(-6.66666666666667d0)
      x556 = gamma_cc*x479
      x557 = x170*x556
      x558 = x197*x204
      x559 = x557*x558
      x560 = x555*x559
      x561 = x4**(-5)
      x562 = gamma_cc*x561
      x563 = 0.296852159568225d0*x562
      x564 = x168*x528
      x565 = x478 + x530
      x566 = x1*x565
      x567 = x537 + x564*x566
      x568 = x197*x551
      x569 = x174*x568
      x570 = x172*x567*x569
      x571 = x551*x553
      x572 = x489*x571
      x573 = x200*x471
      x574 = x489*x573
      x575 = x479*x493
      x576 = x200*x575
      x577 = x172*x204
      x578 = x170*x577
      x579 = x567*x578
      x580 = x573*x579
      x581 = x200*x492 + x476*x576 - 0.757329556733241d0*x574 + 0.272420
     &703860878d0* x580
      x582 = x552*x553 + x554*x560 + x563*x570 - 0.825249003599665d0*x57
     &2 + x581
      x583 = x208*x582
      x584 = x205*x583
      x585 = x202*x206*x584
      x586 = x470*(x473 + x482 - x491 + x550 - x585)
      x587 = x0*x16*(x268 - x277 - x279 + x282 - x289 + x291 + x321 - x3
     &31 - x398 + x419 + x452) + x0*x217*x265 + x0*x220 - x0*x228 - x0*x
     &241 + x0* x248 + x223*x226 + 16.0d0*x223*x254 - x223*(-x256 + x262
     & + x263) + x229*x235 - x249*x250*x253 + 0.0155455d0*x463*x467*x5 -
     & 0.00777275d0*x468*x469*x586
      x588 = 4.0d0*x225
      x589 = rho_s**3
      x590 = x119*x589
      x591 = x120*x590
      x592 = x588*x591
      x593 = 2.56488140124205d0*x233
      x594 = x232*x593
      x595 = x0*x594
      x596 = x176*x251
      x597 = x221*x596
      x598 = 5.1297628024841d0*x121*x597
      x599 = 8.0d0*x254
      x600 = x591*x599
      x601 = 0.00528235082287567d0*x266
      x602 = x1*x601
      x603 = 0.0506605918211689d0*x285
      x604 = rho_c**(-3)
      x605 = x590*x604
      x606 = x588*x605
      x607 = -x606
      x608 = 5.1297628024841d0*x122
      x609 = x596*x608
      x610 = x599*x605
      x611 = -x594 + x607 + x609 - x610
      x612 = mu*x36
      x613 = 0.0253302959105844d0*x612
      x614 = 30.6251363876636d0*x312
      x615 = 3.30853595146347d0*x295
      x616 = x285*(x182*x615 - x303 - x614)
      x617 = 11.4095193139813d0*x295 + x303
      x618 = x503*x617
      x619 = -x184*x618 + x616
      x620 = x106*x294
      x621 = x619*x620
      x622 = 0.00496533372342768d0*x142
      x623 = x1*x622
      x624 = x1*x395
      x625 = x1*x339
      x626 = -x625
      x627 = x2*x349
      x628 = x347*x627
      x629 = x3*x356
      x630 = x347*x629
      x631 = 8.61370200056457d0*x72
      x632 = x365*x369
      x633 = x2*x632
      x634 = x631*x633
      x635 = x375*x379
      x636 = x3*x635
      x637 = x631*x636
      x638 = -x384 + x387 + x628 - x630 + x634 - x637
      x639 = x626 + x638
      x640 = x1*x422
      x641 = x72*(-x432 + x433)
      x642 = 22.9698720015055d0*x641
      x643 = x1*x446
      x644 = x638 - x642 - x643
      x645 = 0.00020688890514282d0*x30
      x646 = x139*x611 + x140*(x417*x611 + x623) + x148*(x165*x611 + x64
     &0 - x644*x645 ) - x287*x603 + x54*(-x391*x639 + x624) + x602 - x61
     &3*x621
      x647 = x16*x646
      x648 = x0*x647
      x649 = x0*x565
      x650 = 0.00777275d0*x649
      x651 = x212*x5
      x652 = x650*x651
      x653 = x477*(0.726455210295675d0*x464 - 0.726455210295675d0*x465)
      x654 = x225*x605
      x655 = -6.05863645386593d0*x654
      x656 = x122*x596
      x657 = x254*x605
      x658 = -x230 + x231
      x659 = x593*x658
      x660 = x106*x118
      x661 = x658*x660
      x662 = x608*x661
      x663 = x605*x660
      x664 = 8.0d0*x663
      x665 = x607 + x659 + x662 + x664
      x666 = x196*x665
      x667 = x466*x477
      x668 = x168*x667
      x669 = x484*(-3.88492098935195d0*x234 - 1.51465911346648d0*x647 + 
     &x655 + 7.76984197870389d0*x656 - 12.1172729077319d0*x657 - 1.51465
     &911346648d0*x666) - x578*(359.577650970909d0*x234 + 140.1927008386
     &36d0*x647 + 560.770803354546d0*x654 - 719.155301941818d0*x656 + 11
     &21.54160670909d0*x657 - 280.385401677273d0*x668) + x653
      x670 = x571*x669
      x671 = 0.726455210295675d0*x480
      x672 = x206*x573
      x673 = x611 - x647
      x674 = -x666 + x673
      x675 = 1.51465911346648d0*x674
      x676 = x484*x675
      x677 = x672*x676
      x678 = gamma_cc*x541
      x679 = x479*x678
      x680 = x500**5.56d0
      x681 = x551*x680
      x682 = x679*x681
      x683 = 0.7916057588486d0*x682
      x684 = x606 - x659 - x662 - x664
      x685 = x500*x684
      x686 = x513*x613
      x687 = x303 + x52*x615 + x614
      x688 = -x60*x617 + x687
      x689 = x503*x688
      x690 = x384 - x387 - x628 + x630 - x634 + x637
      x691 = x625 + x690
      x692 = x642 + x643 + x690
      x693 = x139*x665 + x140*(x417*x665 + x623) + x148*(x165*x665 + x64
     &0 + x645*x692 ) + x504*x603 + x54*(x391*x691 + x624) + x602 + x686
     &*x689
      x694 = x16*x693
      x695 = x665 - x694
      x696 = x685 + x695
      x697 = x497*x696
      x698 = x548*x680
      x699 = x551*x698
      x700 = x697*x699
      x701 = 1.65049800719933d0*x700
      x702 = 0.59370431913645d0*x562
      x703 = 1029.23675661767d0*x654
      x704 = 659.967553630838d0*x658
      x705 = x122*x661
      x706 = 514.618378308835d0*x667
      x707 = -x233*x704 + x525*x706 - 2058.47351323534d0*x663 + 257.3091
     &89154418d0* x694 + x703 - 1319.93510726168d0*x705
      x708 = x526*x707
      x709 = x539**(-2)
      x710 = x681*x709
      x711 = x708*x710
      x712 = x702*x711
      x713 = x471*x546
      x714 = x671*x713
      x715 = x697*x713
      x716 = 1.51465911346648d0*x715
      x717 = x541*x708
      x718 = 0.544841407721756d0*x717
      x719 = x713*x718
      x720 = x714 + x716 - x719
      x721 = x683 + x701 - x712 + x720
      x722 = x210*x548
      x723 = x671*x672 + x677 - x721*x722
      x724 = 0.544841407721756d0*x670 + x723
      x725 = gamma_cc*x724
      x726 = 1d0/x211
      x727 = 0.003886375d0*x4
      x728 = x208*x727
      x729 = x726*x728
      x730 = x725*x729
      x731 = x198*x205
      x732 = x459*(-x205*x461*(x203 + 0.59370431913645d0*x731) - 0.29685
     &2159568225d0* x455*x457 + x460)
      x733 = x727*x732
      x734 = 0.0477583568090938d0*x219
      x735 = x222*x249
      x736 = x132*x21
      x737 = 0.0145402651599119d0*x736
      x738 = x250/rho_c**7
      x739 = x225*x232
      x740 = 82.0762048397455d0*x739
      x741 = x133*x21
      x742 = x232*x741
      x743 = rho_s*x604
      x744 = 0.00242337752665199d0*x227
      x745 = x236*x237
      x746 = x21*x238*x745
      x747 = 0.353517066413871d0*x746
      x748 = x2**(-0.666666666666667d0)
      x749 = 0.854960467080683d0*x1
      x750 = x748*x749
      x751 = x3**(-0.666666666666667d0)
      x752 = x749*x751
      x753 = x750 + x752
      x754 = -5.1297628024841d0*x230 + 5.1297628024841d0*x231 + x753
      x755 = 4.0d0*x754
      x756 = x240*(0.210778258416395d0*x12 + 0.989141964955861d0*x21 + 1
     &.16327369584438d0*x45 + 0.600321476249523d0*x55)
      x757 = 2.0d0*x756
      x758 = x236**2*x239/(x113**4*x114**2)
      x759 = x108*x758
      x760 = x224*x246
      x761 = x128*x246
      x762 = x119*x21*x761
      x763 = 0.107629131579854d0*x762
      x764 = x224*x761
      x765 = x232*x764
      x766 = x242*x243
      x767 = 118.434298572867d0*x766
      x768 = x244*(0.212390288347482d0*x12 + 0.998929717145607d0*x21 + 1
     &.5861639701226d0*x45 + 0.322572935345748d0*x55)
      x769 = x767*x768
      x770 = 164.152409679491d0*x596
      x771 = x239*x745*(x112 + 0.252933910099675d0*x12 + 1.4837129474337
     &9d0*x21 + 1.99418347859036d0*x45)/x113**3
      x772 = 2.0d0*x771
      x773 = x242*x245**2/(x130**4*x131**2)
      x774 = 3506.67076961173d0*x773
      x775 = x245*(0.254868346016978d0*x12 + x129 + 1.49839457571841d0*x
     &21 + 2.71913823449588d0*x45)/x130**3
      x776 = x767*x775
      x777 = x599*x754
      x778 = x232*x325
      x779 = 10.2595256049682d0*x778
      x780 = x257*x258
      x781 = x734 - x747 - x757 - x759 + x772
      x782 = 0.339918983653357d0*x21*x259*x780 - 0.022960662510839d0*x25
     &5 + 2.0d0* x261*(0.267319032365355d0*x12 + 1.70877610553209d0*x21 
     &+ 2.16229511226841d0*x45 + 1.23355033167571d0*x55) + x781 - 2.0d0*
     & x260*x780*(x104 + 0.320782838838426d0*x12 + 2.56316415829814d0* x
     &21 + 3.70679162103156d0*x45)/x105**3 + 19.7392088021787d0*x107* x2
     &57**2*x260/(x105**4*x109**2)
      x783 = x193*x264
      x784 = x148*x26/x15**6
      x785 = x466*x479
      x786 = x120*x18
      x787 = x212*x4
      x788 = x2**(-1.33333333333333d0)
      x789 = 0.333333333333333d0*x1
      x790 = x788*x789
      x791 = x3**(-1.33333333333333d0)
      x792 = x789*x791
      x793 = x790 + x792
      x794 = x124*x176
      x795 = x1*x233
      x796 = x1*x742
      x797 = x1*x232
      x798 = x764*x797
      x799 = x123*x782
      x800 = 16.0d0*x326
      x801 = x323*x779
      x802 = 80.0d0*x536
      x803 = x323*x740
      x804 = x123*x737
      x805 = x123*x760
      x806 = 947.474388582934d0*x805
      x807 = -x734 + x744 + x747 + x757 + x759 - x763 - x769 - x772 - x7
     &74 + x776 + x799 - x800 + x801 - x802 + x803 - x804 + x806
      x808 = x755*x795 - 0.00932351391944647d0*x796 + 607.539859357387d0
     &*x798 + x807
      x809 = x323*x770 - x323*x777 - 160.0d0*x794 + x808
      x810 = 1.79913710865273d0*x26
      x811 = 1.91415600012546d0*x430
      x812 = 0.299209d0*x160
      x813 = 5.29055759308205d0*x346
      x814 = x158*x813
      x815 = x442*x814
      x816 = 31.182216717766d0*x443
      x817 = 23.3866625383245d0*x440*x69*(0.289053642815522d0*x12 - 0.10
     &6975995763985d0*x21)
      x818 = 5.19703611962767d0*x438
      x819 = 30.6264960020074d0*x435
      x820 = 5.0d0*x1
      x821 = x34*x820
      x822 = x35*x820
      x823 = x821 + x822
      x824 = 7.65662400050184d0*x72
      x825 = x1*x824*(-6.0d0*x432 + 6.0d0*x433 + x823)
      x826 = 10.5811151861641d0*x159*x19*x346
      x827 = x440*x441
      x828 = 93.5466501532981d0*x19*x827
      x829 = 77.9555417944151d0*x19*x445
      x830 = 0.0677878489986482d0*x346
      x831 = x353*x830
      x832 = x360*x830
      x833 = 1d0/x73
      x834 = 0.00786976946662482d0*x19
      x835 = x346*x349
      x836 = x74*x835*(0.0157395389332496d0*x350 - x833*x834 - 0.0078697
     &6946662482d0)
      x837 = 1.91415600012546d0*x72
      x838 = x342*x837
      x839 = x344*x837
      x840 = 1d0/x88
      x841 = x346*x356
      x842 = x841*x89*(0.0157395389332496d0*x357 + x834*x840 + 0.0078697
     &6946662482d0)
      x843 = 38.2831200025092d0*x72
      x844 = x843*x86
      x845 = x100*x843
      x846 = 11.4849360007528d0*x371*x73
      x847 = 11.4849360007528d0*x380*x88
      x848 = 0.0543010717965206d0*x12
      x849 = 0.150031049774234d0*x77
      x850 = x19*x833
      x851 = 0.0271505358982603d0*x12
      x852 = x80*x851
      x853 = -x849 - x852
      x854 = x72*x74
      x855 = x369*x854*(0.300062099548467d0*x350*x77 + x350*x80*x848 - x
     &849*x850 - x850*x852 + x853)
      x856 = 0.101681773497972d0*x346
      x857 = x348*x367
      x858 = x366*x857
      x859 = x352*x856*x858
      x860 = x355*x377
      x861 = x856*x860
      x862 = x359*x376*x861
      x863 = 0.150031049774234d0*x92
      x864 = x19*x840
      x865 = x851*x95
      x866 = x863 + x865
      x867 = x72*x89
      x868 = x379*x867*(x357*x848*x95 + 0.300062099548467d0*x357*x92 + x
     &863*x864 + x864*x865 + x866)
      x869 = x346*x358*x629
      x870 = 0.271151395994593d0*x1
      x871 = x869*x870
      x872 = 15.3132480010037d0*x1*x72
      x873 = x386*x872
      x874 = 0.225046574661351d0*x77
      x875 = 0.0325806430779124d0*x12
      x876 = x80*x875
      x877 = x874 + x876
      x878 = x82**(-3)
      x879 = x79*x84*x878
      x880 = x366*x854*x879*(-x350*x874 - x350*x876 + x877)
      x881 = 0.225046574661351d0*x92
      x882 = x875*x95
      x883 = x881 + x882
      x884 = x96**(-3)
      x885 = x884*x94*x98
      x886 = x376*x867*x885*(x357*x881 + x357*x882 + x883)
      x887 = x346*x627
      x888 = x351*x887
      x889 = x870*x888
      x890 = x383*x872
      x891 = x3*x376
      x892 = 45.939744003011d0*x72
      x893 = x1*x892
      x894 = x379*x891*x893
      x895 = x2*x370
      x896 = x893*x895
      x897 = x19*x844 + x19*x845 + x831 + x832 + x836 - x838 - x839 - x8
     &42 + x846 + x847 + x855 + x859 + x862 - x868 + x871 + x873 + x880 
     &+ x886 - x889 - x890 + x894 - x896
      x898 = 0.00082755562057128d0*x30
      x899 = 0.311172355472072d0*x21
      x900 = 0.0121035293512547d0*x12*x27 + 0.0175165224763768d0*x12 + 0
     &.00285512867044833d0*x21 + x25 + 0.000428261234514334d0*x26 - 0.10
     &3724118490691d0*x402
      x901 = 0.335051282142018d0*x392*x69
      x902 = x19*x408
      x903 = -0.335051282142018d0*x150 - 1.34020512856807d0*x392*x902 + 
     &1.67525641071009d0*x420 + 0.670102564284035d0*x424
      x904 = 0.125788454326835d0*x399
      x905 = 0.0441362997638016d0*x143
      x906 = 0.0794453395748429d0*x274*x902
      x907 = 0.0662044496457024d0*x410
      x908 = x407*x899 + x900
      x909 = x274*x69
      x910 = 0.0198613348937107d0*x908*x909
      x911 = 1.27980162581231d0*x146
      x912 = 0.426600541937435d0*x145*x416
      x913 = x106*x414*x754
      x914 = 8.0d0*x913
      x915 = 160.0d0*x534
      x916 = 164.152409679491d0*x535
      x917 = x914 + x915 - x916
      x918 = x808 + x917
      x919 = x417*x918
      x920 = 0.000440195901906306d0*x30
      x921 = x333*(0.427594146176807d0*x12 - 1.36477107997868d0*x21)
      x922 = x332*x69
      x923 = 75.398223686155d0*x67
      x924 = x922*x923
      x925 = 0.154912426780983d0*x0**(-0.333333333333333d0)*x336
      x926 = x335*x66
      x927 = 4.83324359125373d0*x926
      x928 = x922*x927
      x929 = 1.61108119708458d0*x337
      x930 = 150.79644737231d0*x19*x332*x67
      x931 = 9.66648718250747d0*x19*x65*x926
      x932 = 0.0231073345558005d0*x26
      x933 = x22*x273
      x934 = 0.0138644007334803d0*x26
      x935 = 0.00107855455836145d0*x30
      x936 = 0.0350330449527535d0*x12 + 0.00571025734089666d0*x21 - 0.00
     &609784657412857d0*x23 + 0.000856522469028669d0*x26
      x937 = 0.0891107483661069d0*x392
      x938 = x19*x272
      x939 = 0.178221496732214d0*x392
      x940 = x19*x28
      x941 = -0.000586927869208408d0*x102 + 0.445553741830535d0*x19*x62 
     &- 0.0277288014669606d0*x26*x940
      x942 = 4.0d0*x54
      x943 = 0.123631380878626d0*x138
      x944 = 0.00316102962473761d0*x138
      x945 = 0.0337176493305345d0*x138
      x946 = 0.00328744309596812d0*x280
      x947 = 0.00767070055725895d0*x280
      x948 = 1.79063317129192d0*x45
      x949 = 3.80317310466043d0*x295
      x950 = 6.0d0*x283
      x951 = 6.0d0*x284
      x952 = x2**(-1.33333333333333d0)
      x953 = x1*x952
      x954 = x3**(-1.33333333333333d0)
      x955 = x1*x954
      x956 = x953 + x955
      x957 = -x950 + x951 + x956
      x958 = x1*x957
      x959 = x298*x958
      x960 = 3.80317310466043d0*x959
      x961 = x285**2
      x962 = x19*x961
      x963 = x297*x56
      x964 = 22.8190386279626d0*x962*x963
      x965 = 7.60634620932085d0*x300
      x966 = 2.04643791004791d0*x45
      x967 = x304*x957*x966
      x968 = 8.18575164019163d0*x45
      x969 = x961*x968
      x970 = x19*x42*x969
      x971 = x305*x966
      x972 = x967 + x970 - x971
      x973 = 1.10284531715449d0*x295
      x974 = x182*x973
      x975 = x308*x47
      x976 = 2.0d0*x50
      x977 = 3.06965686507186d0*x45
      x978 = 5.70475965699064d0*x295
      x979 = 12.2786274602874d0*x45
      x980 = x305*x979
      x981 = 22.8190386279626d0*x300
      x982 = x49**(-2)
      x983 = x183*x982
      x984 = x975*x983
      x985 = 6.61707190292693d0*x961
      x986 = x19*x963*x985
      x987 = 81.6670303671029d0*x13*x55*x962/x36**4 + 11.4844261453738d0
     &*x312 + 10.2083787958879d0*x314*x958 - 30.6251363876636d0*x315 + 1
     &.79063317129192d0*x45 + x972
      x988 = 0.101321183642338d0*x292
      x989 = 0.000255740352913201d0*x0**2.33333333333333d0*x31
      x990 = 0.0105647016457513d0*x936
      x991 = 0.0352156721525045d0*x275
      x992 = 0.0469542295366726d0*x32
      x993 = 0.101321183642338d0*x36
      x994 = x13/x183**2
      x995 = x319*x994
      x996 = 0.0422588065830054d0*x275
      x997 = 0.270189823046234d0*x36
      x998 = x320*x50
      x999 = x13*x993
      x1000 = x308*x999
      x1001 = 0.405284734569351d0*x612
      x1002 = 0.133819554179517d0*x267 - 0.00657488619193625d0*x280*x940
     &
      x1003 = x16*(-x1000*x998 + x1001*x299*x320 + x1002 - x140*(-x904 +
     & x905 + x906 - x907 + x910 - x911 + x912 - x919) + x148*(-0.003678
     &0249803168d0* x163*x191 + x165*x809 + 4.19798658685637d0*x178*x26 
     &+ 0.0027585187352376d0*x30*(x389 - x431 - x436 + x437 + x439 + x44
     &4 + x447) - x330*x810 - x898*(x811 + x812 + x815 + x816 + x817 - x
     &818 - x819 - x825 + x826 + x828 - x829 + x897) - x901*(-x899*( x40
     &3 + x404 - x405 + x406) + x900) + x903) + x178*x943 + x180* x947 -
     & x180*x989 - x180*x992 - x185*x962*x997 - x273*x991 + 0.0675474557
     &615585d0*x288*x957 - x31*x909*x990 - x330*x945 + x511 *x993*x995 -
     & x620*x988*(-x182*x986 - x184*(x948 + x949 + x960 + x964 - x965 + 
     &x972) - 1.10284531715449d0*x310*x958 + 2.20569063430898d0*x311 - x
     &318*x975*x976 - x974 + x984*(-x977 - x978 + x980 + x981) + x987) +
     & x809*x944 + x933*x946 + x938*x996 + x942*(x180*x932 - x180*x935 -
     & x273*x937 + x390*x920 - x391*(x515 + x897 + x921 + x924 + x925 + 
     &x928 + x929 + x930 + x931) - x393* x69*x936 - x395*x69 + x933*x934
     & + x938*x939 + x941))
      x1004 = -x473 - x482 + x491 - x550 + x585
      x1005 = x208*x726
      x1006 = gamma_cc*x200
      x1007 = x1004*x1005*x1006
      x1008 = x175*x680
      x1009 = x549 + 1.0d0
      x1010 = x1008*x11*x709 + x1009
      x1011 = 1d0/x1010
      x1012 = x1009*x546
      x1013 = x555*x680
      x1014 = x1013*x678
      x1015 = x526*x538
      x1016 = x1009*x547
      x1017 = x1011*x548
      x1018 = x1016*x1017
      x1019 = x8/x211**2
      x1020 = x1004*x469
      x1021 = x200*x475
      x1022 = rho_c**(-5.33333333333333d0)
      x1023 = rho_s*x1022
      x1024 = x1023*x480
      x1025 = x200*x206
      x1026 = x1024*x1025
      x1027 = x479**2
      x1028 = x1027*x6
      x1029 = 2.9058208411827d0*rho_c**(-6.33333333333333d0)*x18
      x1030 = x1028*x1029
      x1031 = 2.66666666666667d0*x464
      x1032 = 2.66666666666667d0*x465
      x1033 = 0.444444444444444d0*x1
      x1034 = x1033*x788
      x1035 = x1033*x791
      x1036 = x1034 + x1035
      x1037 = x1023*x477
      x1038 = x1037*(x1031 - x1032 + x1036)
      x1039 = 1.08968281544351d0*x1038
      x1040 = x1021*x489
      x1041 = x488**2
      x1042 = x195**(-2)
      x1043 = x1021*x1042
      x1044 = 2.10537616771841d0*x1043
      x1045 = 4.03909096924395d0*x1024
      x1046 = x179*x330
      x1047 = 0.757329556733241d0*x1040
      x1048 = 0.757329556733241d0*x488
      x1049 = x545*x548
      x1050 = x546*x573
      x1051 = x1049*x1050
      x1052 = x474*x546
      x1053 = x545*x679
      x1054 = 2.9058208411827d0*rho_s
      x1055 = x195/x178**2
      x1056 = 5.56d0*x1055
      x1057 = x122 - 1.0d0
      x1058 = x1057*x225
      x1059 = x1*x1058*x755
      x1060 = x1057*x797
      x1061 = 0.00466175695972324d0*x1060*x736
      x1062 = 0.00466175695972323d0*x796
      x1063 = 303.769929678693d0*x798
      x1064 = 303.769929678693d0*x1060*x760
      x1065 = 2.78d0*x196
      x1066 = x179*x487
      x1067 = 5.5573346245051d0*x784
      x1068 = 7.03121958266113d0*x264
      x1069 = 22.24d0*x754
      x1070 = x323*x778
      x1071 = x323*x739
      x1072 = x107*x758
      x1073 = x123*x736
      x1074 = x754*x795
      x1075 = x766*x768
      x1076 = x766*x775
      x1077 = 28.5214811818116d0*x1070 + 228.171849454493d0*x1071 + 27.4
     &375002350284d0 *x1072 - 0.0404219371445552d0*x1073 + 11.12d0*x1074
     & - 329.247350032569d0*x1075 + 329.247350032569d0*x1076 - 0.1327682
     &31929281d0*x219 + 0.00673698952409253d0*x227 - 44.48d0* x326 - 222
     &.4d0*x536 + 0.982777444630562d0*x746 + 5.56d0*x756 - 0.29920898579
     &1994d0*x762 - 5.56d0*x771 - 9748.54473952062d0*x773 - 0.0259193686
     &960612d0*x796 + 1688.96080901354d0*x798 + 2.78d0* x799 + 2633.9788
     &0026056d0*x805
      x1078 = -2.78d0*x1003 + x1056*x416**2 - x1065*(-x1059 + x1061 - x1
     &062 + x1063 - x1064 + x807 + x917) - 5.56d0*x1066*x416 - x1067*x19
     &3 - x1068* x486 - x1069*x254*x323 + x1077 + 456.343698908985d0*x32
     &3*x596 + 11.7186993044352d0*x783 - 444.8d0*x794
      x1079 = x1021*x484
      x1080 = 0.272420703860878d0*x1079
      x1081 = x204*x557*x583
      x1082 = x490*x573
      x1083 = x201*x582
      x1084 = x545*x546
      x1085 = x204*x540
      x1086 = 1.5832115176972d0*rho_s*x560 + x570*x702 + 2.7706201559701
     &d0*x571 - 1.65049800719933d0*x572
      x1087 = 1.08968281544351d0*x174
      x1088 = x207**(-2)
      x1089 = x10*x1088
      x1090 = 8.47531078678287d0*x474
      x1091 = x474*x484
      x1092 = x1091*x488
      x1093 = x1042*x474
      x1094 = 2.10537616771841d0*x1093
      x1095 = x474*x578
      x1096 = x474*x567**2
      x1097 = 514.618378308835d0*x170
      x1098 = exp(-x1097*x168)
      x1099 = x4**(-6)
      x1100 = 0.544841407721756d0*x1099
      x1101 = x1098*x1100*x174
      x1102 = 0.272420703860878d0*x1099
      x1103 = x1102*x577
      x1104 = 1.45291042059135d0*x479
      x1105 = x577*x9
      x1106 = x1104*x1105
      x1107 = 0.757329556733241d0*x1092
      x1108 = x1048*x1093
      x1109 = 2058.47351323534d0*x477
      x1110 = x1*x1109*(2.0d0*x464 - 2.0d0*x465 + x793)
      x1111 = x479*x565
      x1112 = x168*x6
      x1113 = 5489.26270196091d0*x19
      x1114 = 38.2831200025092d0*x72
      x1115 = x1114*x86
      x1116 = x100*x1114
      x1117 = -x1115*x19 - x1116*x19 - x831 - x832 - x836 + x838 + x839 
     &+ x842 - x846 - x847 - x855 - x859 - x862 + x868 - x871 - x873 - x
     &880 - x886 + x889 + x890 - x894 + x896
      x1118 = x22*x502
      x1119 = x19*x501
      x1120 = x52*x973
      x1121 = x47*x507
      x1122 = x59*x982
      x1123 = x1121*x1122
      x1124 = -x953 - x955
      x1125 = x503**2
      x1126 = x40/x59**2
      x1127 = x1126*x512
      x1128 = x16*(x1*x1001*x503*x514 + x1000*x50*x514 + x1002 + x1118*x
     &946 - x1119* x996 + x1125*x19*x61*x997 + x1127*x999*(x313 - x316 +
     & x506 + x508 - x510) + x137*x943 + x140*(x904 - x905 - x906 + x907
     & - x910 + x911 - x912 + x919) + x148*(x165*x918 - x416*x810 + 0.00
     &36780249803168d0*x426 + 4.19798658685637d0*x428 - 0.00275851873523
     &76d0*x449 + x898*(x1117 - x811 - x812 - x815 - x816 - x817 + x818 
     &+ x819 + x825 - x826 - x828 + x829) - x901* x908 + x903) + x20*x27
     &5*x990 - x29*x947 + x29*x989 + x29*x992 - x416*x945 - x502*x991 + 
     &0.0675474557615585d0*x505*(x1124 + x950 - x951) + x513*x988*(x1120
     & + x1121*x511*x976 + x1123*(x977 + x978 - x980 - x981) - 2.2056906
     &3430898d0*x509 + 1.10284531715449d0*x52* x959 + x52*x986 + x60*(-x
     &948 - x949 - x960 - x964 + x965 - x967 - x970 + x971) + x987) + x9
     &18*x944 + x942*(x1118*x934 - x1119*x939 + x20*x395 - x29*x932 + x2
     &9*x935 + x391*(x1117 + x340 - x921 - x924 - x925 - x928 - x929 - x
     &930 - x931) - x502*x937 - x516*x920 + x517*x936 + x941))
      x1129 = x264*x518
      x1130 = x166*x784
      x1131 = -2639.87021452335d0*x1070 - 21118.9617161868d0*x1071 - 253
     &9.53990571917d0*x1072 + 3.74134383838717d0*x1073 - 1029.2367566176
     &7d0*x1074 + 30474.2333338565d0*x1075 - 30474.2333338565d0*x1076 + 
     &257.309189154418d0*x1128 + 650.790435101151d0*x1129 + 514.37167845
     &0064d0*x1130 + 12.2886640658953d0*x219 - 0.623557306397862d0*x227 
     &+ 4116.94702647068d0*x326 - 1084.65072516858d0*x533 - 41169.470264
     &7068d0*x534 + 42237.9234323736d0*x535 + 20584.7351323534d0*x536 - 
     &90.9631897112017d0*x746 - 514.618378308835d0*x756 + 27.69396457620
     &63d0*x762 + 514.618378308835d0*x771 + 902298.612360293d0*x773 + 2.
     &3990258066827d0*x796 - 156325.588590238d0*x798 - 257.309189154418d
     &0*x799 - 243793.866670852d0*x805 - 2058.47351323534d0*x913
      x1132 = 0.272420703860878d0*x1095
      x1133 = 0.272420703860878d0*x1091
      x1134 = x202*x553
      x1135 = x494*x546
      x1136 = x1037*x1135
      x1137 = x494**2
      x1138 = x546*x6
      x1139 = 1.08968281544351d0*x546
      x1140 = x1013*x548
      x1141 = rho_c**(-7.66666666666667d0)*rho_s*x680
      x1142 = x1141*x678
      x1143 = x1142*x494
      x1144 = x540*x9
      x1145 = 0.888888888888889d0*x1
      x1146 = x1145*x788
      x1147 = x1145*x791
      x1148 = x1052*x523
      x1149 = x168**(-2)
      x1150 = x1149*x522**2
      x1151 = 2.10537616771841d0*x1052
      x1152 = x1052*x526
      x1153 = gamma_cc*x1013
      x1154 = x1153*x709
      x1155 = x1154*x561
      x1156 = x1015*x1155
      x1157 = x538**2
      x1158 = x1052*x1157
      x1159 = exp(x1097*x525)
      x1160 = x1159*x709
      x1161 = x1100*x1160
      x1162 = x4**(-8)
      x1163 = x1157*x1162
      x1164 = x1153*x1159/x539**3
      x1165 = 0.890556478704675d0*x1164
      x1166 = 0.296852159568225d0*x526
      x1167 = x1015*x1099
      x1168 = x1141*x494
      x1169 = gamma_cc*x1168*x709
      x1170 = 0.757329556733241d0*x1148
      x1171 = x416*x499
      x1172 = x1140*x523
      x1173 = 4.58838446001414d0*x1140
      x1174 = 0.825249003599665d0*x1172
      x1175 = x498 + x519 + x520
      x1176 = x1175*x522
      x1177 = x1052*x1149
      x1178 = 0.757329556733241d0*x1177
      x1179 = x1140*x1149
      x1180 = x526*x542
      x1181 = x525*x6
      x1182 = x1109*x1175*x532 + x1110*x525 + x1113*x1181*x494*x531 + x1
     &131
      x1183 = x1052*x527*x541
      x1184 = x1155*x1166
      x1185 = x500*(x1059 - x1061 + x1062 - x1063 + x1064 - x744 + x763 
     &+ x769 + x774 - x776 + x781 - x799 + x800 - x801 + x802 - x803 + x
     &804 - x806 - x914 - x915 + x916)
      x1186 = x499*x521
      x1187 = 5.56d0*x498
      x1188 = x168/x137**2
      x1189 = x1188*x498**2
      x1190 = x1052*x497
      x1191 = 0.272420703860878d0*x1190
      x1192 = 11.12d0*x498
      x1193 = x1140*x497
      x1194 = 0.148426079784112d0*x1193
      x1195 = 0.00777275d0*x1005
      x1196 = x221*x590
      x1197 = 0.00181753314498899d0*x736
      x1198 = x121*x249
      x1199 = 20.5190512099364d0*x739
      x1200 = 0.00116543923993081d0*x741
      x1201 = x1200*x232
      x1202 = 2.56488140124205d0*x231
      x1203 = 2.56488140124205d0*x230
      x1204 = x1202 - x1203 + x753
      x1205 = x1204*x233
      x1206 = 75.9424824196734d0*x764
      x1207 = x1206*x232
      x1208 = 41.0381024198728d0*x596
      x1209 = x1204*x254
      x1210 = 1.28244070062102d0*x778
      x1211 = 2.0d0*x325
      x1212 = x264*x646
      x1213 = 0.0103636666666667d0*x787
      x1214 = x17*x651
      x1215 = 0.000220097950953153d0*x30
      x1216 = x1*x65*x927
      x1217 = x1*x332*x923
      x1218 = 0.135575697997296d0*x1
      x1219 = x1218*x835
      x1220 = 0.0338939244993241d0*x346
      x1221 = x1220*x627
      x1222 = x1218*x841
      x1223 = x1220*x629
      x1224 = 0.0677878489986482d0*x887*(x350 - 1.0d0)
      x1225 = 0.135575697997296d0*x888
      x1226 = 7.65662400050184d0*x72
      x1227 = x1226*x386
      x1228 = 0.0677878489986482d0*x869
      x1229 = x1226*x383
      x1230 = x365*x857
      x1231 = x2*x347
      x1232 = x1230*x1231*x351
      x1233 = x1*x382
      x1234 = x1233*x632
      x1235 = x341*x633
      x1236 = x347*x860
      x1237 = x1236*x3*x358*x375
      x1238 = x1233*x635
      x1239 = x341*x636
      x1240 = x1231*x858
      x1241 = x1236*x891
      x1242 = x350*x849
      x1243 = x350*x852
      x1244 = x369*x631
      x1245 = x1244*x2*(x1242 + x1243 + x853)
      x1246 = x382*x895
      x1247 = x382*x891
      x1248 = x1247*x379
      x1249 = x357*x863 + x357*x865 + x866
      x1250 = x379*x631
      x1251 = x1249*x1250*x3
      x1252 = x631*x877
      x1253 = x1252*x879
      x1254 = x1253*x2*x366
      x1255 = x631*x883
      x1256 = x1255*x885
      x1257 = x1256*x891
      x1258 = x1*x937
      x1259 = x1*x28
      x1260 = -0.267332245098321d0*x1*x62 + x1259*x934
      x1261 = 3.0d0*x432
      x1262 = 3.0d0*x433
      x1263 = x1*(x156 - x157)
      x1264 = x93 - 1.0d0
      x1265 = x1264*x3
      x1266 = 46.773325076649d0*x1*x827
      x1267 = x78 - 1.0d0
      x1268 = x1267*x2
      x1269 = x1268*x368
      x1270 = x1265*x378
      x1271 = x1264*x378
      x1272 = x849 + x852
      x1273 = x1*x1114
      x1274 = x1219 - x1221 + x1222 + x1223 + x1228 - x1232 + x1237 - x1
     &240 + x1241
      x1275 = 2.0d0*x122
      x1276 = 20.5190512099364d0*x122
      x1277 = x1211*x605
      x1278 = 16.0d0*x654
      x1279 = x1197*x605
      x1280 = x605*x760
      x1281 = 118.434298572867d0*x1280
      x1282 = x1277 + x1278 + x1279 - x1281
      x1283 = x1201 - x1205 - x1207 + x1209*x1275 - x1210*x122 - x1276*x
     &739 + x1282 - 41.0381024198728d0*x656 + 32.0d0*x657
      x1284 = 0.0003448148419047d0*x30
      x1285 = 0.449784277163183d0*x26
      x1286 = x1*x408
      x1287 = -0.251288461606513d0*x1*x149 + 0.167525641071009d0*x1286*x
     &392
      x1288 = 0.106650135484359d0*x145
      x1289 = -0.0182062236525682d0*x1*x142 + 0.00993066744685536d0*x128
     &6*x274
      x1290 = 0.00842941233263362d0*x138
      x1291 = 3.0d0*x283
      x1292 = 3.0d0*x284
      x1293 = -x1291 + x1292 + x956
      x1294 = 0.0168868639403896d0*x1293
      x1295 = 1.02321895502395d0*x45
      x1296 = 15.3125681938318d0*x312
      x1297 = x46*x969
      x1298 = x310*x985
      x1299 = 10.2083787958879d0*x312
      x1300 = 22.8190386279626d0*x295 + x979
      x1301 = x1300*x503
      x1302 = x50*x618
      x1303 = x1302*x47
      x1304 = x1124 + x1291 - x1292
      x1305 = x1304*x949
      x1306 = x503*x949
      x1307 = x503*x981
      x1308 = x1304*x966
      x1309 = x1295*x503
      x1310 = x305*x503*x968
      x1311 = x1308 - x1309 + x1310
      x1312 = x613*x620
      x1313 = 0.0675474557615585d0*x36
      x1314 = x1313*x961
      x1315 = x1314*x185
      x1316 = x503*x687
      x1317 = 0.0253302959105844d0*x1316
      x1318 = x603*x612
      x1319 = x1*x276
      x1320 = 0.0253302959105844d0*x13
      x1321 = x1320*x617
      x1322 = x1321*x285
      x1323 = mu*x621
      x1324 = -0.0193686196838775d0*x1*x266 + x1259*x281
      x1325 = x16*(x1*x1315 + x1*x1323*x603 + x1283*x944 - x1290*x611 - 
     &x1294*x287 - x1312*(-x1*x1297 + x1*x1298 - 81.6670303671029d0*x1*x
     &314*x961 - x1293*x1299 - x1293*x966 + x1293*x974 + x1295*x285 + x1
     &296*x285 + x1301*x984 - x1303*x318 - x184*(x1305 - x1306 + x1307 +
     & x1311) - x285*x974 - x50*x616*x975) + x1317*x995 - x1318*x320 - x
     &1319*x272 + x1322*x998 + x1324 + x140*(x1283*x417 - x1288*x611 + x
     &1289) + x148*(x1283*x165 + x1284*x644 - x1285*x611 + x1287 + x645*
     &(x1226* x1265*x99 - x1226*x1268*x85 - x1233*x1267*x365*x368 - x123
     &3*x1271 *x375 - x1247*x1271 + x1249*x1270*x631 + x1252*x1268*x366*
     &x84* x878 - x1255*x1264*x884*x891*x98 + 62.364433435532d0*x1263*x4
     &40 - x1263*x814 - x1264*x1273*x99 + x1266 - x1267*x1273*x85 + x126
     &9* x341*x365 + x1269*x366*x382 - x1269*x631*(-x1242 - x1243 + x127
     &2 ) - x1270*x341*x375 + x1274 - 15.3132480010037d0*x434 - x824*( -
     &x1261 + x1262 + x823) - 0.0677878489986482d0*x888)) + x54*(x1215 *
     &x639 - x1258*x272 + x1260 - x391*(-x1*x1115 - x1*x1116 - x1216 - x
     &1217 - x1219 + x1221 - x1222 - x1223 + x1224 + x1225 - x1227 - x12
     &28 + x1229 + x1232 - x1234 + x1235 - x1237 - x1238 - x1239 + x1240
     & - x1241 + x1245 + x1246 - x1248 + x1251 + x1254 - x1257 + x625)))
     &
      x1326 = x1012*x471
      x1327 = 0.003886375d0*x477
      x1328 = x1011*x1327
      x1329 = x1019*x1328*(x1018*(-x683 - x701 + x712 - x714 - x716 + x7
     &19) + x1326* x671 + 1.51465911346648d0*x1326*x697 + 0.544841407721
     &756d0*x699*( x497*(3.88492098935195d0*x233*x658 + x655 + 12.117272
     &9077319d0* x663 + 1.51465911346648d0*x685 - 1.51465911346648d0*x69
     &4 + 7.76984197870389d0*x705) + x653 - x718))
      x1330 = 1.45291042059135d0*x1028
      x1331 = x1023*x1330
      x1332 = 1.33333333333333d0*x464
      x1333 = 1.33333333333333d0*x465
      x1334 = x1036 + x1332 - x1333
      x1335 = x484*x674
      x1336 = 2.01954548462197d0*x1335
      x1337 = x555*x669
      x1338 = x1337*x553
      x1339 = 1.00977274231099d0*x480
      x1340 = x179*x611
      x1341 = x1042*x695
      x1342 = x488*x674
      x1343 = x179*x673
      x1344 = -x1202 + x1203 - x750 - x752
      x1345 = x1058*x1344
      x1346 = x122*x325
      x1347 = x1346*x658
      x1348 = 1.28244070062102d0*x1347
      x1349 = x1200*x658
      x1350 = 32.0d0*x663
      x1351 = x1344*x660
      x1352 = x1275*x1351
      x1353 = 10.2595256049682d0*x122
      x1354 = x225*x658
      x1355 = x1353*x1354
      x1356 = x118*x251
      x1357 = x1276*x1356
      x1358 = x1353*x739
      x1359 = 20.5190512099364d0*x705
      x1360 = x1206*x658
      x1361 = x1282 + x1348 - x1349 - x1350 + x1352 + x1355 + x1357 - x1
     &358 - x1359 + x1360
      x1362 = x122*x778
      x1363 = x1209*x122
      x1364 = x122*x739
      x1365 = -114.085924727246d0*x1364
      x1366 = x325*x605
      x1367 = x605*x736
      x1368 = -658.494700065139d0*x1280 + 11.12d0*x1366 + 0.010105484286
     &1388d0*x1367 + 88.96d0*x654
      x1369 = x1056*x416*x665 - 2.78d0*x1066*x665 - 5.56d0*x1205 - 3.515
     &60979133057d0* x1212 - 5.56d0*x1325 - 2.78d0*x1343*x416 - 7.130370
     &29545289d0* x1362 + 11.12d0*x1363 + x1365 + x1368 - 5.56d0*x196*(-
     &x1345 + x1361) - 228.171849454493d0*x656 + 177.92d0*x657 + 0.00647
     &98421740153d0*x742 - 422.240202253384d0*x765
      x1370 = x208*x721
      x1371 = x200*x473
      x1372 = x1371*x548
      x1373 = x1021*x476
      x1374 = 0.757329556733241d0*x1082
      x1375 = x1374*x548
      x1376 = x10*x1370
      x1377 = x1084*x203*x709
      x1378 = 0.726455210295675d0*x672
      x1379 = x10*x583
      x1380 = 0.544841407721756d0*x669
      x1381 = -1.1874086382729d0*x562*x711 + 1.5832115176972d0*x682 + 3.
     &30099601439866d0*x700
      x1382 = x1089*x209
      x1383 = x1382*(-x1139*x471*x717 + x1381 + 1.45291042059135d0*x480*
     &x713 + 3.02931822693296d0*x715)
      x1384 = x1085*x582
      x1385 = 1.69506215735657d0*x575
      x1386 = 0.544841407721756d0*x493
      x1387 = x1091*x674
      x1388 = 659.967553630838d0*x234 + 257.309189154418d0*x647 - 1319.9
     &3510726168d0* x656 + 2058.47351323534d0*x657 - 514.618378308835d0*
     &x668 + x703
      x1389 = x1105*x479
      x1390 = x474*x567
      x1391 = x1388*x1390
      x1392 = x1388*x578
      x1393 = x594 + x606 - x609 + x610 + x647
      x1394 = 2744.63135098045d0*x479
      x1395 = 60948.4666677131d0*x1280 - 1029.23675661767d0*x1366 - 0.93
     &5335959596793d0*x1367 - 8233.89405294136d0*x654
      x1396 = x1013*x679
      x1397 = 2.63868586282867d0*x1144
      x1398 = x1052*x697
      x1399 = x1155*x708
      x1400 = x1144*x708
      x1401 = x1099*x708
      x1402 = x1140*x697
      x1403 = x546*x575
      x1404 = x1154*x479
      x1405 = x499*x665
      x1406 = x684 + x694
      x1407 = x1149*x522
      x1408 = x1407*x696
      x1409 = x1052*x538
      x1410 = x1162*x538
      x1411 = x1*x844 + x1*x845 - x1224 - x1225 + x1227 - x1229 + x1234 
     &- x1235 + x1238 + x1239 - x1245 - x1246 + x1248 - x1251 - x1254 + 
     &x1257 + x1274
      x1412 = x1344*x233
      x1413 = x1361 + x1412
      x1414 = 6.61707190292693d0*x52
      x1415 = mu*x1125
      x1416 = x513*x688
      x1417 = x16*(0.0506605918211689d0*x1*x1415*x1416 - x1127*x13*x1317
     & - x1290*x665 + x1294*x504 + x1302*x1320*x514 + x1313*x299*x503*x6
     &1 + x1318* x514 + x1319*x501 + x1324 + x140*(-x1288*x665 + x1289 +
     & x1413* x417) + x1413*x944 + x148*(-x1284*x692 - x1285*x665 + x128
     &7 + x1413*x165 + x645*(x1*x159*x813 - 62.364433435532d0*x1*x445 + 
     &x1266 + x1411 + 15.3132480010037d0*x641 + x824*(x1261 - x1262 - x8
     &21 - x822))) + x54*(-x1215*x691 + x1258*x501 + x1260 + x391*( x121
     &6 + x1217 + x1411 + x626)) + x686*(x1120*x1304 - x1120*x503 + x112
     &1*x1316*x50 - x1123*x1301 - x1296*x503 + x1299*x1304 - x1303* x511
     & + x1311 + x1414*x300*x503 + 81.6670303671029d0*x315*x503 + x60*(-
     &x1305 + x1306 - x1307 - x1308 + x1309 - x1310)))
      x1418 = x658*x741
      x1419 = x264*x693
      x1420 = x122*x1354
      x1421 = x122*x1356
      x1422 = x122*x1351
      x1423 = x658*x764
      x1424 = x1175*x706 + x1181*x1394*x532 - x1346*x704 + 5279.74042904
     &67d0*x1364 + x1395 + x1406*x528*x532 - 514.618378308835d0*x1412 + 
     &514.618378308835d0*x1417 + 0.599756451670674d0*x1418 + 325.3952175
     &50575d0*x1419 - 5279.7404290467d0*x1420 - 10559.4808580934d0*x1421
     & - 1029.23675661767d0*x1422 - 39081.3971475595d0*x1423 + x529*(x53
     &1 - x790 - x792) + 16467.7881058827d0*x663 + 10559.4808580934d0*x7
     &05
      x1425 = x1186*x684
      x1426 = x499*x695
      x1427 = x500*(-x1277 - x1278 - x1279 + x1281 + x1345 - x1348 + x13
     &49 + x1350 - x1352 - x1355 - x1357 + x1358 + x1359 - x1360)
      x1428 = x1188*x684
      x1429 = gamma_cc*x729
      x1430 = 12.0d0*x225
      x1431 = x119*x786
      x1432 = x748 + x751
      x1433 = x1432*x233
      x1434 = 0.854960467080683d0*x1433
      x1435 = 1.70992093416137d0*x254
      x1436 = 24.0d0*x254
      x1437 = x952 + x954
      x1438 = 0.0168868639403896d0*x1437
      x1439 = x119*x19
      x1440 = x1436*x1439
      x1441 = x589*x604
      x1442 = x1199*x1441
      x1443 = x1208*x1441
      x1444 = x122*x1432
      x1445 = x1435*x1444
      x1446 = x1430*x1439
      x1447 = x1434 - x1446
      x1448 = -x1440 + x1442 + x1443 - x1445 + x1447
      x1449 = 0.203363546995945d0*x346
      x1450 = x1449*x349
      x1451 = x1449*x356
      x1452 = x1230*x856
      x1453 = x375*x861
      x1454 = x632*x892
      x1455 = x635*x892
      x1456 = x1244*x1272
      x1457 = x1250*x866
      x1458 = x1253*x365
      x1459 = x1256*x375
      x1460 = x1450 + x1451 + x1452 + x1453 + x1454 + x1455 - x1456 - x1
     &457 + x1458 + x1459 + x844 + x845
      x1461 = x1437*x966
      x1462 = x1299*x1437
      x1463 = 81.6670303671029d0*x314
      x1464 = x1415*x42*x617
      x1465 = x1300*x1464
      x1466 = x1125*x298
      x1467 = x1125*x46
      x1468 = x1461 + x1467*x968
      x1469 = x1437*x949 + 22.8190386279626d0*x1466 + x1468
      x1470 = 38.2831200025092d0*x36*x72
      x1471 = x16*(x106*x1317*x46*x619*x994 - x1312*(-mu*x42*x616*x618*x
     &976 + x1297 - x1298 - x1437*x974 + x1461 + x1462 + x1463*x961 + x1
     &465*x983 - x1469*x184) - x1315 + x1322*x46*x50*x621 - x1323*x286 +
     & x140*( x1448*x417 + x622) + x1438*x287 + x1448*x944 + x148*(x1448
     &*x165 + x422 - x645*(x1460 - x1470 - x446)) + x54*(-x391*(x1460 - 
     &x339) + x395) + x601)
      x1472 = x788 + x791
      x1473 = x1472*x477
      x1474 = 4.21075233543682d0*x1042*x674**2
      x1475 = 4.03909096924395d0*x1335*x480
      x1476 = x1340*x676
      x1477 = x1341*x675
      x1478 = -x1434 + x1440 - x1442 - x1443 + x1445 + x1446
      x1479 = x665**2
      x1480 = x1343*x665
      x1481 = x1441*x739
      x1482 = x1441*x596
      x1483 = x1439*x254
      x1484 = x1444*x254
      x1485 = x1439*x225
      x1486 = 2.3767900984843d0*x1433 - 33.36d0*x1485
      x1487 = 3.02931822693296d0*x1335
      x1488 = x1370*x672
      x1489 = x204*x568
      x1490 = x1376*x540
      x1491 = x709*x721
      x1492 = x1099*x1388**2
      x1493 = x1354*x1441
      x1494 = x1444*x660
      x1495 = x1441*x661
      x1496 = x1439*x660
      x1497 = x6*x785
      x1498 = 24.0d0*x1496
      x1499 = 41.0381024198728d0*x1495
      x1500 = 1.70992093416137d0*x1494
      x1501 = 20.5190512099364d0*x1493
      x1502 = x1447 + x1498 + x1499 + x1500 - x1501
      x1503 = -x1115 - x1116 - x1450 - x1451 - x1452 - x1453 - x1454 - x
     &1455 + x1456 + x1457 - x1458 - x1459
      x1504 = x16*(mu*x286*x513*x689 - x1126*x1320*x1467*x687*x688 + x13
     &14*x61 + x1321 *x1416*x1467*x50 + x140*(x1502*x417 + x622) - x1438
     &*x504 + x148*( x1502*x165 + x422 + x645*(x1470 + x1503 + x446)) + 
     &x1502*x944 + x54*(x391*(x1503 + x339) + x395) + x601 + x686*(x1120
     &*x1437 + x1122*x1465 + x1125*x1463 + x1414*x1466 + x1462 - x1464*x
     &687*x976 + x1468 - x1469*x60))
      x1505 = x696**2
      x1506 = x707**2
      x1507 = x1162*x1506
      x1508 = 0.59370431913645d0*x526
      x1509 = x1406*x696
      x1510 = 1029.23675661767d0*x1406*x667 - 219.989184543613d0*x1433 +
     & 171.539459436278d0*x1473*x525 + 3087.71026985301d0*x1485 + 5279.7
     &404290467d0*x1493 - 439.978369087225d0*x1494 - 10559.4808580934d0*
     &x1495 - 6175.42053970602d0*x1496 + 1372.31567549023d0*x1497*x525 +
     & 257.309189154418d0*x1504
      x1511 = x500*(0.854960467080683d0*x1058*x1432 + x1446 - x1498 - x1
     &499 - x1500 + x1501)
      x1512 = x1426*x684
      x1513 = x1188*x684**2
      x1514 = x1008*x548
      x1515 = 0.59370431913645d0*x1514 + x547
      x1516 = -x1515*x722 + x209 + 0.296852159568225d0*x731
      x1517 = 0.3958028794243d0*x557
      x1518 = x1515*x208
      x1519 = x1518*x679
      x1520 = x1139*x201
      x1521 = x10*x1518
      x1522 = x205*x210
      x1523 = 1.1874086382729d0*x1514 + x1520
      x1524 = x1382*x1523
      x1525 = x1011*x1515
      x1526 = x1525*x548
      x1527 = -x1016*x1526 + x1016 + 0.296852159568225d0*x1514
      x1528 = x1328*x726
      x1529 = 0.593704319136449d0*x1008
      E = -rho_c*(-x212*x213 + x216)
      d1E(1) = 0.5d0*rho_c*x587 + x16*x217 + x178 + x213*x463
      d1E(2) = -rho_c*(x592 + x595 - x598 + x600 + x648 + x652 - x730)
      d1E(3) = rho_c*x470*x733
      d2E(1) = -0.25d0*rho_c*(x1003*x17 - 0.062182d0*x1007*x467 + 0.0155
     &455d0*x1011* x1019*x1020*x477*(-x1012*x472 + x1012*x496 + x1012*x5
     &24 + x1018*( -x1014*x494*x554 + x1015*x563*x710 + x492*x546 - x496
     &*x546 - 0.825249003599665d0*x523*x699 - x524*x546 + x543*x713 + x5
     &52*x698 ) + x202*x698*(-x492 + x496 + x524 - x544)) - x1195*x468*(
     & -2.17936563088702d0*x10*x1083*x1084*x1085*x469 - x1021*x1054* x10
     &81 + 8.47531078678287d0*x1021 + x1025*x1030 + x1025*x1039 + x1025*
     &x1083*x1087*x1089*(1.45291042059135d0*rho_s*x576 + x1086 + 2.54259
     &323603486d0*x573 - 1.51465911346648d0*x574 + 0.544841407721756d0*x
     &580) + 6.7802486294263d0*x1026 - 7.06840919617691d0*x1040 + x1041*
     &x1044 - x1043*x1048*x521 - x1045 *x200*x490 + x1046*x1047 - 3.0293
     &1822693296d0*x1049*x546*x574 + 5.08518647206972d0*x1051 + x1052*x1
     &053*x1054*x200 + x1078*x1080 + 3.02931822693296d0*x1082*x584 + x11
     &34*(x1023*x1106*x567 + 6.7802486294263d0*x1024 + x1030 + x1039 + x
     &1041*x1094 - x1045* x489 + x1046*x1107 + x1078*x1133 + x1090 - 1.5
     &1465911346648d0* x1092*x579 - 7.06840919617691d0*x1092 + 2.5425932
     &3603486d0*x1095* x567 + x1096*x1101 - x1096*x1103 - x1108*x521 - x
     &1132*(x1109*x521 *x566 - x1110*x168 - x1111*x1112*x1113 + x1131)) 
     &- 5.08518647206972d0*x584*x672 - x722*(5.27737172565733d0*gamma_cc
     &* rho_c**(-8.66666666666667d0)*x1137*x1144*x18*x680 - 1.4529104205
     &9135d0*x1015*x1023*x1135*x1144 + x1029*x1137*x1138 + x1038*x1139 +
     & x1090*x546 - x1099*x1158*x527*x540 + 4.03909096924395d0*x1136*x52
     &3 - 6.7802486294263d0*x1136 + 15.7001808838306d0*x1140 + 0.5937043
     &1913645d0*x1142*(x1146 + x1147 + 5.33333333333333d0*x464 - 5.33333
     &333333333d0*x465) + 8.80265603839643d0*x1143*x523 - 14.77664083184
     &05d0*x1143 - 1.51465911346648d0*x1148*x1180 - 7.06840919617691d0*x
     &1148 + x1150 *x1151 + x1150*x1173 + 2.54259323603486d0*x1152*x542 
     &- x1154* x1163*x1166 - 3.30099601439866d0*x1156*x523 + 5.541240311
     &9402d0* x1156 + x1158*x1161 + x1163*x1165 - 3.1664230353944d0*x116
     &7*x1169 + x1170*x1171 + x1171*x1174 - 15.4046480671938d0*x1172 + x
     &1176* x1178 + 0.825249003599665d0*x1176*x1179 - x1182*x1183 - x118
     &2* x1184 + x1191*(-x1067*x166 - x1068*x518 + x1077 - 2.78d0*x1128 
     &+ 2.78d0*x1185 + x1186*x1187 + 5.56d0*x1189 + 11.7186993044352d0* 
     &x533 + 444.8d0*x534 - 456.343698908985d0*x535 + 22.24d0*x913) + x1
     &194*(x1069*x795 + 57.0429623636231d0*x1070 + 456.343698908985d0 *x
     &1071 + 54.8750004700568d0*x1072 - 0.0808438742891104d0*x1073 - 658
     &.494700065139d0*x1075 + 658.494700065139d0*x1076 - 5.56d0* x1128 -
     & 14.0624391653223d0*x1129 - 11.1146692490102d0*x1130 + 5.56d0*x118
     &5 + x1186*x1192 + 11.12d0*x1189 - 0.265536463858562d0* x219 + 0.01
     &34739790481851d0*x227 - 88.96d0*x326 + 23.4373986088704d0*x533 + 8
     &89.6d0*x534 - 912.68739781797d0*x535 - 444.8d0*x536 + 1.9655548892
     &6112d0*x746 + 11.12d0*x756 - 0.598417971583988d0*x762 - 11.12d0*x7
     &71 - 19497.0894790412d0*x773 - 0.0518387373921224d0*x796 + 3377.92
     &161802707d0*x798 + 5.56d0* x799 + 5267.95760052111d0*x805 + 44.48d
     &0*x913))) + 1.99904123183637d0*x17*x193*x784 + 2.52921567721623d0*
     &x17*x264* x486 + x17*x734 - x17*x744 - x17*x747 - x17*x757 - x17*x
     &759 + x17 *x763 + x17*x769 + x17*x772 + x17*x774 - x17*x776 - 4.21
     &535946202706d0*x17*x783 + 80.0d0*x225*x735 - x233*x743*x755 + 160.
     &0d0*x254*x735 + 16.0d0*x325*x735 + 0.031091d0*x651*x743*( 2.0d0*x4
     &64 - 2.0d0*x465 + x793) + x735*x737 - 947.474388582934d0* x735*x76
     &0 - x735*x782 - x738*x740 - x738*x770 + x738*x777 - x738* x779 + 0
     &.00932351391944647d0*x742*x743 - 607.539859357387d0*x743* x765 - 0
     &.0414546666666667d0*x785*x786*x787) + x587
      d2E(2) = rho_c*(-x1007*x650 - x1020*x1329 + x1195*x467*x725 + x119
     &6*x1197 + x1196 *x1211 + 16.0d0*x1196*x225 + 32.0d0*x1196*x254 - 1
     &18.434298572867d0*x1196*x760 - x1198*x1199 - x1198*x1208 + 2.0d0 *
     &x1198*x1209 - x1198*x1210 + x1201*x17 - x1205*x17 - x1207*x17 - 0.
     &632303919304059d0*x1212*x17 - x1213*x743*x785 + 0.00777275d0* x121
     &4*(x466 + x793) - x1325*x17 + x1429*(-rho_c**( -7.66666666666667d0
     &)*x476*x559*x669 - x1021*x1048*x1341 - 0.544841407721756d0*x1021*x
     &1334*x477 - 3.53420459808846d0*x1021* x1335 - 1.69506215735657d0*x
     &1021*x480 - x1025*x1331 - x1026*x1336 + x1040*x1339 + x1044*x1342 
     &+ x1047*x1340 - 0.726455210295675d0* x1050*x1053 - x1051*x676 + x1
     &080*x1369 + x1081*x1378 + x1134*( x1022*x1388*x1389*x476 - x1024*x
     &1336 + x1092*x1339 + x1094*x1342 + 1.27129661801743d0*x1095*x1388 
     &+ x1101*x1391 - x1103*x1391 + x1107*x1340 - x1107*x1392 - x1108*x6
     &95 - x1132*(x1112*x1394*x566 + 514.618378308835d0*x1205 + 325.3952
     &17550575d0*x1212 + 514.618378308835d0*x1325 + 659.967553630838d0*x
     &1362 - 1029.23675661767d0*x1363 + 10559.4808580934d0*x1364 - x1393
     &*x528* x566 + x1395 - x521*x706 + x564*(x565 + x793) + 21118.96171
     &61868d0*x656 - 16467.7881058827d0*x657 - 0.599756451670674d0*x742 
     &+ 39081.3971475595d0*x765) + x1133*x1369 - x1331 - x1334*x1386 - x
     &1385 - 0.757329556733241d0*x1387*x579 - 3.53420459808846d0*x1387 -
     & 0.363227605147837d0*x1389*x1390) + 0.757329556733241d0*x1338*x489
     & - 1.27129661801743d0*x1338 + x1370 *x1372 + x1370*x1373*x679 - x1
     &370*x1375 + x1376*x1377 + x1379* x1380*x569 - x1383*x1384 + x584*x
     &677 - x722*(0.296852159568225d0* x1014*(-x1031 + x1032 - x1146 - x
     &1147) - x1022*x1400*x495*x546 + x1023*x1104*x1138*x494 + 1.2712966
     &1801743d0*x1052*x717 - x1102* x1409*x540*x708 + 2.01954548462197d0
     &*x1136*x697 + 0.825249003599665d0*x1140*x1406*x1407 + 4.4013280191
     &9821d0*x1143* x697 - 0.363227605147837d0*x1144*x1152*x479*x538 + x
     &1151*x1408 - 0.296852159568225d0*x1154*x1410*x708 - 1.650498007199
     &33d0*x1156* x697 + x1161*x1409*x707 + x1165*x1410*x707 - 0.7916057
     &588486d0* x1167*x1404 + x1168*x1397*x556 - 1.5832115176972d0*x1169
     &*x1401 + x1170*x1405 - x1170*x717 + x1173*x1408 + x1174*x1405 + x1
     &178* x1406*x522 - 0.757329556733241d0*x1180*x1398 - x1183*x1424 - 
     &x1184*x1424 + x1191*(x1187*x1428 + 7.13037029545289d0*x1347 - 57.0
     &429623636231d0*x1364 + x1368 + 5.56d0*x1412 - 5.56d0*x1417 - 0.006
     &4798421740153d0*x1418 - 3.51560979133057d0*x1419 + 57.042962363623
     &1d0*x1420 + 114.085924727246d0*x1421 + 11.12d0* x1422 + 422.240202
     &253384d0*x1423 + 2.78d0*x1425 + 2.78d0*x1426* x498 + 5.56d0*x1427 
     &- 177.92d0*x663 - 114.085924727246d0*x705) + x1194*(-x1068*x693 + 
     &x1187*x1426 + x1192*x1428 - 1316.98940013028d0*x1280 + 14.26074059
     &09058d0*x1347 + x1365 + 22.24d0*x1366 + 0.0202109685722776d0*x1367
     & + 11.12d0*x1412 - 11.12d0*x1417 - 0.0129596843480306d0*x1418 + 11
     &4.085924727246d0* x1420 + 228.171849454493d0*x1421 + 22.24d0*x1422
     & + 844.480404506768d0*x1423 + 5.56d0*x1425 + 11.12d0*x1427 + 177.9
     &2d0*x654 - 355.84d0*x663 - 228.171849454493d0*x705) - x1385* x546 
     &+ x1386*x546*(-x1034 - x1035 - x1332 + x1333) + 2.20066400959911d0
     &*x1396*x523 - 3.69416020796013d0*x1396 - 3.53420459808846d0*x1398 
     &- 1.65049800719933d0*x1399*x523 + 2.7706201559701d0*x1399 - 7.7023
     &2403359688d0*x1402 + 1.00977274231099d0*x1403*x523))) - x592 - x59
     &5 + x598 - x600 - x648 - x652 + x730
      d2E(3) = -rho_c*(0.0155455d0*x1005*x649*x725 - x1111*x1213*x17 + x
     &1198*x1432* x1435 - x1199*x221*x589 + 0.00259091666666667d0*x1214*
     &x1472 + x1329*x208*x724 - x1429*(gamma_cc*x1104*x1337*x170*x558 + 
     &0.544841407721756d0*rho_c**(-6.66666666666667d0)*x553*(x1087* x109
     &8*x1492 - x1106*x1388 + x1330 - x1392*x1487 + x1474 + x1475 + x147
     &6 - x1477 - 0.544841407721756d0*x1492*x577 + x477*( 0.242151736765
     &225d0*x788 + 0.242151736765225d0*x791) + x484*( 3.02931822693296d0
     &*x1055*x1479 + 1.29497366311732d0*x1433 - 1.51465911346648d0*x1471
     & + 1.51465911346648d0*x1478*x196 - 3.02931822693296d0*x1480 + 31.0
     &793679148156d0*x1481 + 62.1587358296312d0*x1482 - 36.3518187231956
     &d0*x1483 - 2.58994732623463d0*x1484 - 18.1759093615978d0*x1485) - 
     &x578*( 560.770803354546d0*x1393*x667 - 119.859216990303d0*x1433 - 
     &93.461800559091d0*x1473*x168 + 1682.31241006364d0*x1485 + 2876.621
     &20776727d0*x1493 - 239.718433980606d0*x1494 - 5753.24241553455d0*x
     &1495 - 3364.62482012728d0*x1496 - 747.694404472728d0*x1497*x168 + 
     &140.192700838636d0*x1504)) + x1021*x1330 + 0.242151736765225d0*x10
     &21*x1473 + x1021*x1474 + x1021*x1475 + x1021*x1476 - x1021*x1477 +
     & 0.544841407721756d0* x1079*(x1056*x1479 + x1065*x1478 - 2.78d0*x1
     &471 - 5.56d0*x1480 + 57.0429623636231d0*x1481 + 114.085924727246d0
     &*x1482 - 66.72d0* x1483 - 4.75358019696859d0*x1484 + x1486) - x110
     &4*x1488*x678 + x1338*x1487 + x1383*x1491 - x1487*x1488*x548 - 1.08
     &968281544351d0 *x1489*x1490*x669 - x722*(0.263868586282867d0*x1014
     &*x1472 + x1027 *x1153*x1397 - x1052*x1104*x1400 + x1052*x1330 + x1
     &099*x1139* x1160*x1506*x474 - x1100*x1152*x1506*x540 - 0.544841407
     &721756d0* x1152*x1510*x541 - x1154*x1507*x1508 - x1155*x1508*x1510
     & + 1.78111295740935d0*x1164*x1507 + 4.21075233543682d0*x1177*x1505
     & + 1.51465911346648d0*x1177*x1509 + 9.17676892002827d0*x1179*x1505
     & + 1.65049800719933d0*x1179*x1509 + 0.544841407721756d0*x1190*(x14
     &86 - 57.0429623636231d0*x1493 + 4.75358019696859d0*x1494 + 114.085
     &924727246d0*x1495 + 66.72d0*x1496 - 2.78d0*x1504 + 2.78d0* x1511 +
     & 5.56d0*x1512 + 5.56d0*x1513) + 0.296852159568225d0*x1193* (4.7535
     &8019696859d0*x1433 - 66.72d0*x1485 - 114.085924727246d0* x1493 + 9
     &.50716039393719d0*x1494 + 228.171849454493d0*x1495 + 133.44d0*x149
     &6 - 5.56d0*x1504 + 5.56d0*x1511 + 11.12d0*x1512 + 11.12d0*x1513) +
     & 8.80265603839643d0*x1396*x697 + 1.51465911346648d0*x1398*x1405 - 
     &3.02931822693296d0*x1398*x717 - 6.60199202879732d0*x1399*x697 - 3.
     &1664230353944d0*x1401*x1404 + 1.65049800719933d0*x1402*x1405 + 4.0
     &3909096924395d0*x1403*x697 + 0.242151736765225d0*x1472*x493*x546))
     & + x1430*x1431 + x1431*x1436 - x1434*x17 + x1471*x17 - 41.03810241
     &98728d0*x589*x597)
      d2E(4) = x470*(rho_c*(x1006*x1088*x1327*x1516*x586 + 0.00777275d0*
     &x467*x732 - x728*(rho_s*x1517*x555*x558 + x1021*x481 + x1049*x1520
     &*x200 + x1371 - x1372*x1518 - x1373*x1519 - x1374 + x1375*x1518 - 
     &x1377* x1521 - 0.296852159568225d0*x1379*x199 + x1384*x1524 - x152
     &2*x582 - x1522*(x1086 + x581) + 0.692655038992524d0*x571 - 0.41262
     &4501799832d0*x572)) + x733)
      d2E(5) = -rho_c*x1005*(x1516*x650 + x1527*x1528*x725 - x727*(0.825
     &249003599665d0* x1335*x571 - x1378*x1519 - x1380*x1489*x1521*x540 
     &+ x1489*x1517 - 0.296852159568225d0*x1490*x198*x204 + x1491*x1524 
     &- x1518*x548* x677 + 1.08968281544351d0*x670 - x722*(x1381 + x720)
     & + x723))
      d2E(6) = -rho_c*x1528*(x1005*x1516*x1527 + x540*(0.32347469700879d
     &0*rho_c**( -7.0d0)*x1009*x1017*x500**8.34d0 + x1009*x1520*x1525 + 
     &x1526* x1529 - x1529 - x1016*x1515*x1523*x548/x1010**2))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_LDA_ERF_case_1(rho_a, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea
      Ea = 0.0d0
      Ea = -1.8610514726982d0*rho_a**1.33333333333333d0
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_LDA_ERF_case_2(rho_a, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4
      Ea = 0.0d0
      x0 = rho_a**0.333333333333333d0
      x1 = 1d0/x0
      x2 = mu*x1
      x3 = 0.00844343197019482d0*1d0/rho_a
      x4 = mu**2
      Ea = 4.9628039271952d0*rho_a**1.33333333333333d0*(0.12827824385304
     &2d0*x2*(mu **3*x3 - mu*(-0.256556487706084d0*x1 + x3*x4)*exp( -15.
     &192666241151986d0*rho_a**0.66666666666666663d0/x4) - 0.38483473155
     &9127d0*x2 + 1.77245385090552d0*erf( 3.8977770897207535d0*x0/mu)) -
     & 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_LDA_ERF_case_3(rho_a, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea
      Ea = 0.0d0
      Ea = -3.14159265358979d0*rho_a**2.0d0/mu**2
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_LDA_ERF_case_1(rho_a, mu, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      Ea = -1.8610514726982d0*rho_a**1.33333333333333d0
      d1Ea(1) = -2.4814019635976d0*rho_a**0.333333333333333d0
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_LDA_ERF_case_2(rho_a, mu, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**0.333333333333333d0
      x1 = 1d0/x0
      x2 = 0.256556487706084d0*x1
      x3 = mu**2
      x4 = 0.00844343197019482d0*1d0/rho_a
      x5 = x3*x4
      x6 = -x2 + x5
      x7 = exp(-15.192666241151986d0*rho_a**0.66666666666666663d0/x3)
      x8 = mu*x7
      x9 = 1d0/mu
      x10 = mu**3
      x11 = mu*x1
      x12 = x10*x4 - 0.384834731559127d0*x11 + 1.77245385090552d0*erf( 3
     &.8977770897207535d0*x0*x9)
      x13 = x12 - x6*x8
      x14 = 0.128278243853042d0*x11
      x15 = rho_a**1.33333333333333d0
      x16 = 4.9628039271952d0*x15
      x17 = 1d0/x15
      x18 = 0.00844343197019482d0*rho_a**(-2.0d0)
      x19 = x7*x9
      Ea = x16*(x13*x14 - 0.375d0)
      d1Ea(1) = -mu*x16*(-0.128278243853042d0*x1*(0.128278243853042d0*mu
     &*x17 + 2.59851805981384d0*rho_a**(-0.666666666666667d0)*x19 + 10.1
     &28444160768d0*rho_a**(-0.333333333333333d0)*x19*x6 - x10*x18 + x8*
     &(-0.0855188292353615d0*x17 + x18*x3)) + 0.0427594146176807d0 *x13*
     &x17) + 6.61707190292693d0*rho_a**0.333333333333333d0*(x14*( x12 + 
     &x8*(x2 - x5)) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = mu**(-2)
      Ea = -3.14159265358979d0*rho_a**2.0d0*x0
      d1Ea(1) = -6.28318530717958d0*rho_a**1.0d0*x0
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_LDA_ERF_case_1(rho_a, mu, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      Ea = -1.8610514726982d0*rho_a**1.33333333333333d0
      d1Ea(1) = -2.4814019635976d0*rho_a**0.333333333333333d0
      d2Ea(1) = -0.827133987865866d0*rho_a**(-0.666666666666667d0)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_LDA_ERF_case_2(rho_a, mu, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**0.333333333333333d0
      x1 = 1d0/x0
      x2 = 0.256556487706084d0*x1
      x3 = mu**2
      x4 = 1d0/rho_a
      x5 = 0.00844343197019482d0*x4
      x6 = x3*x5
      x7 = -x2 + x6
      x8 = exp(-15.192666241151986d0*rho_a**0.66666666666666663d0/x3)
      x9 = mu*x8
      x10 = 1d0/mu
      x11 = mu**3
      x12 = mu*x1
      x13 = x11*x5 - 0.384834731559127d0*x12 + 1.77245385090552d0*erf( 3
     &.8977770897207535d0*x0*x10)
      x14 = x13 - x7*x9
      x15 = 0.128278243853042d0*x12
      x16 = rho_a**1.33333333333333d0
      x17 = 4.9628039271952d0*x16
      x18 = rho_a**0.333333333333333d0
      x19 = x2 - x6
      x20 = x13 + x19*x9
      x21 = x15*x20 - 0.375d0
      x22 = 1d0/x16
      x23 = x10*x8
      x24 = rho_a**(-0.333333333333333d0)*x23
      x25 = 10.128444160768d0*x24
      x26 = 0.00844343197019482d0*rho_a**(-2.0d0)
      x27 = 0.0855188292353615d0*x22
      x28 = x26*x3 - x27
      x29 = rho_a**(-0.666666666666667d0)
      x30 = 0.564189583547756d0
      x31 = x23*x30
      x32 = 0.128278243853042d0*mu*x22 - x11*x26 + x28*x9 + 4.6057533417
     &6456d0*x29* x31
      x33 = 0.128278243853042d0*x1
      x34 = 0.0427594146176807d0*x14*x22 - x33*(x25*x7 + x32)
      x35 = mu*x17
      x36 = rho_a**(-2.33333333333333d0)
      x37 = 0.0168868639403896d0*rho_a**(-3.0d0)
      x38 = x8/x11
      Ea = x17*(x14*x15 - 0.375d0)
      d1Ea(1) = 6.61707190292693d0*x18*x21 - x34*x35
      d2Ea(1) = -13.2341438058539d0*mu*x18*x34 + 2.20569063430898d0*x21*
     &x29 + x35*( 0.0570125528235743d0*x20*x36 - x27*(-x19*x25 + x32) + 
     &x33*( -0.171037658470723d0*mu*x36 - 3.07050222784304d0*rho_a**( -1
     &.66666666666667d0)*x31 + 3.37614805358933d0*rho_a**( -1.3333333333
     &3333d0)*x19*x23 + x11*x37 + 102.585381117795d0*x19* x29*x38 - 20.2
     &56888321536d0*x24*x28 - 46.6491155403329d0*x30*x38* x4 + x9*(-x3*x
     &37 + 0.114025105647149d0*x36)))
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = mu**(-2)
      x1 = 6.28318530717958d0*x0
      Ea = -3.14159265358979d0*rho_a**2.0d0*x0
      d1Ea(1) = -rho_a**1.0d0*x1
      d2Ea(1) = -x1
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu**2*rho_a**(-0.666666666666667d0)
      Ea = -1.8610514726982d0*x0*(1.804d0 - 0.804d0/(0.0519881082712438d
     &0*gamma_aa* rho_a**(-2.66666666666667d0)*(1.18478216491348d0*x1 - 
     &7.0d0)*exp( -0.31265084907439067d0*x1)/(0.233951463049512d0*mu**4/
     &x0 + 49.1113297378198d0*mu*rho_a**(-0.333333333333333d0) - 10.6630
     &394842213d0*x1 - 81.0d0) + 1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu**2
      x2 = rho_a**0.666666666666667d0
      x3 = x1/x2
      x4 = rho_a**0.333333333333333d0
      x5 = 1d0/x4
      x6 = mu*x5
      x7 = 0.00844343197019482d0*1d0/rho_a
      Ea = 4.9628039271952d0*x0*(1.804d0 - 0.804d0/(0.0519881082712438d0
     &*gamma_aa* rho_a**(-2.66666666666667d0)*(1.18478216491348d0*x3 - 7
     &.0d0)*exp( -0.31265084907439067d0*x3)/(0.233951463049512d0*mu**4/x
     &0 - 10.6630394842213d0*x3 + 49.1113297378198d0*x6 - 81.0d0) + 1.0d
     &0)) *(0.128278243853042d0*x6*(mu**3*x7 - mu*(x1*x7 - 0.25655648770
     &6084d0*x5)*exp(-15.192666241151986d0*x2/x1) - 0.384834731559127d0*
     &x6 + 1.77245385090552d0*erf( 3.8977770897207535d0*x4/mu)) - 0.375d
     &0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu**2
      x2 = 1d0/x1
      x3 = rho_a**0.666666666666667d0
      x4 = rho_a**0.333333333333333d0
      x5 = 1d0/x4
      x6 = mu*x5
      x7 = 0.00844343197019482d0*1d0/rho_a
      Ea = 4.9628039271952d0*x0*(1.804d0 + 0.804d0/(0.0519881082712438d0
     &*gamma_aa* rho_a**(-2.66666666666667d0)*x2*(0.213719543995407d0*x0
     &*x2 - 0.844037013397333d0*x3 + 14.4026298885989d0*rho_a**2.0d0/mu*
     &*4)* exp(-0.31265084907439067d0*x1/x3) - 1.0d0))*(0.12827824385304
     &2d0* x6*(mu**3*x7 - mu*(x1*x7 - 0.256556487706084d0*x5)*exp( -15.1
     &92666241151986d0*x2*x3) - 0.384834731559127d0*x6 + 1.7724538509055
     &2d0*erf(3.8977770897207535d0*x4/mu)) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**0.333333333333333d0
      x2 = 1d0/x1
      x3 = mu*x2
      x4 = 0.00844343197019482d0*1d0/rho_a
      x5 = erf(3.8977770897207535d0*x1/mu)
      x6 = mu**2
      x7 = 1d0/x6
      x8 = rho_a**0.666666666666667d0
      x9 = 15.192666241152d0*x7*x8
      x10 = x6/x8
      x11 = exp(x9)
      x12 = mu**4*(x11 - 1.0d0)/x0
      x13 = 3.0d0*x11
      Ea = 4.9628039271952d0*x0*(1.804d0 + 0.804d0/(3.15935190989552d0*g
     &amma_aa* rho_a**(-2.0d0)*x7*(0.00060945584614891d0*x10*(7.0d0*x11 
     &+ 11.0d0 ) - 0.000722072416819482d0*x12 + 0.0185185185185185d0)*ex
     &p( -0.31265084907439067d0*x10)/(-0.131642462768164d0*x10*(x13 - 2.
     &0d0) + 1.81893813843777d0*x11*x3*x5 + 0.00866486900183379d0*x12 - 
     &x13) - 1.0d0))*(0.128278243853042d0*x3*(mu**3*x4 - mu*( -0.2565564
     &87706084d0*x2 + x4*x6)*exp(-x9) - 0.384834731559127d0* x3 + 1.7724
     &5385090552d0*x5) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3
      Ea = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**2.0d0
      x3 = rho_a**0.666666666666667d0
      Ea = -3.14159265358979d0*x1*x2*(1.804d0 + 0.804d0/(0.0519881082712
     &438d0* gamma_aa*rho_a**(-2.66666666666667d0)*x1*(0.213719543995407
     &d0* rho_a**1.33333333333333d0*x1 - 0.844037013397333d0*x3 + 14.402
     &6298885989d0*x2/mu**4)*exp(-0.31265084907439067d0*x0/x3) - 1.0d0))
     &
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea, d1E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu**2
      x2 = rho_a**(-0.666666666666667d0)*x1
      x3 = 1.18478216491348d0*x2 - 7.0d0
      x4 = 1d0/x0
      x5 = 1d0/(0.233951463049512d0*mu**4*x4 + 49.1113297378198d0*mu*rho
     &_a**( -0.333333333333333d0) - 10.6630394842213d0*x2 - 81.0d0)
      x6 = exp(-0.31265084907439067d0*x2)
      x7 = x3*x5*x6
      x8 = 0.0519881082712438d0*rho_a**(-2.66666666666667d0)
      x9 = gamma_aa*x7*x8 + 1.0d0
      x10 = 1.804d0 - 0.804d0/x9
      x11 = x9**(-2)
      x12 = rho_a**(-4.33333333333333d0)*x1
      Ea = -1.8610514726982d0*x0*x10
      d1Ea(1) = -1.49628538404935d0*gamma_aa*x0*x11*x5*x6*(mu*x3*x5*x8*(
     & 0.311935284066016d0*mu**3*rho_a**(-2.33333333333333d0) - 7.108692
     &98948088d0*mu*rho_a**(-1.66666666666667d0) + 16.3704432459399d0*x4
     &) - 0.138634955389984d0*rho_a**( -3.66666666666667d0)*x3 + 0.01083
     &60841285172d0*x12*x3 - 0.0410630556449071d0*x12) - 2.4814019635976
     &d0*rho_a** 0.333333333333333d0*x10
      d1Ea(2) = -0.0777890465506374d0*x11*x4*x7
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea, d
     &1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**0.333333333333333d0
      x1 = 1d0/x0
      x2 = 0.256556487706084d0*x1
      x3 = mu**2
      x4 = 0.00844343197019482d0*1d0/rho_a
      x5 = x3*x4
      x6 = -x2 + x5
      x7 = rho_a**0.666666666666667d0
      x8 = exp(-15.192666241151986d0*x7/x3)
      x9 = mu*x8
      x10 = 1d0/mu
      x11 = mu**3
      x12 = mu*x1
      x13 = x11*x4 - 0.384834731559127d0*x12 + 1.77245385090552d0*erf( 3
     &.8977770897207535d0*x0*x10)
      x14 = x13 - x6*x9
      x15 = 0.128278243853042d0*x12
      x16 = x14*x15 - 0.375d0
      x17 = rho_a**1.33333333333333d0
      x18 = x3/x7
      x19 = 1.18478216491348d0*x18 - 7.0d0
      x20 = 1d0/x17
      x21 = 1d0/(0.233951463049512d0*mu**4*x20 + 49.1113297378198d0*x12 
     &- 10.6630394842213d0*x18 - 81.0d0)
      x22 = exp(-0.31265084907439067d0*x18)
      x23 = x19*x21*x22
      x24 = 0.0519881082712438d0*rho_a**(-2.66666666666667d0)
      x25 = gamma_aa*x23*x24 + 1.0d0
      x26 = 1.804d0 - 0.804d0/x25
      x27 = 4.9628039271952d0*x17*x26
      x28 = x15*(x13 + x9*(x2 - x5)) - 0.375d0
      x29 = 0.00844343197019482d0*rho_a**(-2.0d0)
      x30 = x10*x8
      x31 = x25**(-2)
      x32 = rho_a**(-4.33333333333333d0)*x3
      Ea = x16*x27
      d1Ea(1) = 3.99009435746494d0*gamma_aa*x17*x21*x22*x28*x31*(mu*x19*
     &x21*x24*( -7.10869298948088d0*mu*rho_a**(-1.66666666666667d0) + 0.
     &311935284066016d0*rho_a**(-2.33333333333333d0)*x11 + 16.3704432459
     &399d0*x20) - 0.138634955389984d0*rho_a**( -3.66666666666667d0)*x19
     & + 0.0108360841285172d0*x19*x32 - 0.0410630556449071d0*x32) - mu*x
     &27*(-0.128278243853042d0*x1*( 0.128278243853042d0*mu*x20 + 2.59851
     &805981384d0*rho_a**( -0.666666666666667d0)*x30 + 10.128444160768d0
     &*rho_a**( -0.333333333333333d0)*x30*x6 - x11*x29 + x9*( -0.0855188
     &292353615d0*x20 + x29*x3)) + 0.0427594146176807d0*x14* x20) + 6.61
     &707190292693d0*rho_a**0.333333333333333d0*x26*x28
      d1Ea(2) = 0.207437457468366d0*x16*x20*x23*x31
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea, d
     &1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**0.666666666666667d0
      x1 = 0.844037013397333d0*x0
      x2 = mu**(-4)
      x3 = rho_a**2.0d0
      x4 = 14.4026298885989d0*x2*x3
      x5 = mu**2
      x6 = 1d0/x5
      x7 = rho_a**1.33333333333333d0
      x8 = 0.213719543995407d0*x6*x7
      x9 = -x1 + x4 + x8
      x10 = 0.0519881082712438d0*rho_a**(-2.66666666666667d0)*x6
      x11 = exp(-0.31265084907439067d0*x5/x0)
      x12 = gamma_aa*x11
      x13 = x10*x12*x9
      x14 = x13 - 1.0d0
      x15 = rho_a**0.333333333333333d0
      x16 = 1d0/x15
      x17 = 0.256556487706084d0*x16
      x18 = rho_a**1.0d0
      x19 = 0.00844343197019482d0/x18
      x20 = x19*x5
      x21 = -x17 + x20
      x22 = exp(-15.192666241151986d0*x0*x6)
      x23 = mu*x22
      x24 = 1d0/mu
      x25 = mu**3
      x26 = mu*x16
      x27 = x19*x25 - 0.384834731559127d0*x26 + 1.77245385090552d0*erf( 
     &3.8977770897207535d0*x15*x24)
      x28 = -x21*x23 + x27
      x29 = 0.128278243853042d0*x26
      x30 = x28*x29 - 0.375d0
      x31 = 4.9628039271952d0*x7
      x32 = rho_a**0.333333333333333d0
      x33 = -x13 + 1.0d0
      x34 = 1.804d0 - 0.804d0/x33
      x35 = x29*(x23*(x17 - x20) + x27) - 0.375d0
      x36 = x1 - x4 - x8
      x37 = rho_a**(-0.333333333333333d0)
      x38 = 1d0/x7
      x39 = 0.00844343197019482d0/x3
      x40 = x22*x24
      Ea = x30*x31*(1.804d0 + 0.804d0/x14)
      d1Ea(1) = -mu*x31*x34*(-0.128278243853042d0*x16*(0.128278243853042
     &d0*mu*x38 + 2.59851805981384d0*rho_a**(-0.666666666666667d0)*x40 +
     & 10.128444160768d0*x21*x37*x40 + x23*(-0.0855188292353615d0*x38 + 
     &x39*x5) - x25*x39) + 0.0427594146176807d0*x28*x38) + 3.99009435746
     &494d0*x12*x35*x7*(0.0108360841285172d0*rho_a**( -4.33333333333333d
     &0)*x36 - 0.138634955389984d0*rho_a**( -3.66666666666667d0)*x36*x6 
     &+ x10*(-28.8052597771978d0*x18*x2 - 0.284959391993876d0*x32*x6 + 0
     &.562691342264888d0*x37))/x33**2 + 6.61707190292693d0*x32*x34*x35
      d1Ea(2) = -0.207437457468366d0*x11*x30*x38*x6*x9/x14**2
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea, d
     &1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 1d0/x0
      x2 = mu**2
      x3 = 1d0/x2
      x4 = rho_a**0.666666666666667d0
      x5 = 15.192666241152d0*x3*x4
      x6 = exp(x5)
      x7 = mu**4*(x6 - 1.0d0)
      x8 = x1*x7
      x9 = 0.000722072416819482d0*x8
      x10 = 7.0d0*x6 + 11.0d0
      x11 = x2/x4
      x12 = 0.00060945584614891d0*x10*x11
      x13 = x12 - x9 + 0.0185185185185185d0
      x14 = 3.0d0*x6
      x15 = x14 - 2.0d0
      x16 = rho_a**0.333333333333333d0
      x17 = 1d0/x16
      x18 = mu*x17
      x19 = 1d0/mu
      x20 = erf(3.8977770897207535d0*x16*x19)
      x21 = x20*x6
      x22 = 1d0/(-0.131642462768164d0*x11*x15 - x14 + 1.81893813843777d0
     &*x18*x21 + 0.00866486900183379d0*x8)
      x23 = rho_a**(-2.0d0)
      x24 = 3.15935190989552d0*x23*x3
      x25 = x22*x24
      x26 = exp(-0.31265084907439067d0*x11)
      x27 = gamma_aa*x26
      x28 = x13*x25*x27
      x29 = x28 - 1.0d0
      x30 = 0.256556487706084d0*x17
      x31 = 1d0/rho_a
      x32 = 0.00844343197019482d0*x31
      x33 = x2*x32
      x34 = -x30 + x33
      x35 = exp(-x5)
      x36 = mu*x35
      x37 = mu**3
      x38 = -0.384834731559127d0*x18 + 1.77245385090552d0*x20 + x32*x37
      x39 = -x34*x36 + x38
      x40 = 0.128278243853042d0*x18
      x41 = x39*x40 - 0.375d0
      x42 = 4.9628039271952d0*x0
      x43 = x40*(x36*(x30 - x33) + x38) - 0.375d0
      x44 = -x28 + 1.0d0
      x45 = 1.804d0 - 0.804d0/x44
      x46 = mu*x1
      x47 = 0.00844343197019482d0*x23
      x48 = 0.564189583547756d0
      x49 = rho_a**(-0.666666666666667d0)
      x50 = x19*x35
      x51 = rho_a**(-0.333333333333333d0)
      x52 = -x12 + x9 - 0.0185185185185185d0
      x53 = x31*x6
      x54 = rho_a**(-1.66666666666667d0)*x2
      x55 = x54*x6
      x56 = rho_a**(-2.33333333333333d0)*x7
      Ea = x41*x42*(1.804d0 + 0.804d0/x29)
      d1Ea(1) = -mu*x42*x45*(0.0427594146176807d0*x1*x39 - 0.12827824385
     &3042d0*x17*( 10.128444160768d0*x34*x50*x51 + x36*(-0.0855188292353
     &615d0*x1 + x2*x47) - x37*x47 + 0.128278243853042d0*x46 + 4.6057533
     &4176456d0* x48*x49*x50)) + 6.61707190292693d0*rho_a**0.33333333333
     &3333d0*x43 *x45 + 3.99009435746494d0*x0*x22*x27*x43*(0.65851603810
     &2422d0* rho_a**(-3.66666666666667d0)*x52 - 6.31870381979104d0*rho_
     &a**( -3.0d0)*x3*x52 + x24*(0.000406303897432606d0*x10*x54 - 0.0432
     &098765432099d0*x53 + 0.00731347015378692d0*x55 - 0.000962763222425
     &976d0*x56) + x25*x52*(-0.087761641845443d0*x15* x54 - 18.423013367
     &0582d0*x19*x21*x49 + 0.606312712812591d0*x21* x46 + 30.38533248230
     &4d0*x3*x51*x6 - 4.72654360241471d0*x31*x48 + 4.0d0*x53 - 0.0877616
     &41845443d0*x55 + 0.0115531586691117d0*x56))/ x44**2
      d1Ea(2) = -12.6061122289202d0*x13*x22*x26*x3*x41*x49/x29**2
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea, d1E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**2.0d0
      x3 = rho_a**0.666666666666667d0
      x4 = 0.844037013397333d0*x3
      x5 = mu**(-4)
      x6 = 14.4026298885989d0*x2*x5
      x7 = 0.213719543995407d0*rho_a**1.33333333333333d0*x1
      x8 = -x4 + x6 + x7
      x9 = 0.0519881082712438d0*rho_a**(-2.66666666666667d0)*x1
      x10 = exp(-0.31265084907439067d0*x0/x3)
      x11 = gamma_aa*x10
      x12 = x11*x8*x9
      x13 = x12 - 1.0d0
      x14 = rho_a**1.0d0
      x15 = -x12 + 1.0d0
      x16 = x4 - x6 - x7
      Ea = -3.14159265358979d0*x1*x2*(1.804d0 + 0.804d0/x13)
      d1Ea(1) = -x1*(2.52584049348619d0*x11*x2*(0.0108360841285172d0*rho
     &_a**( -4.33333333333333d0)*x16 - 0.138634955389984d0*rho_a**( -3.6
     &6666666666667d0)*x1*x16 + x9*(0.562691342264888d0*rho_a**( -0.3333
     &33333333333d0) - 0.284959391993876d0*rho_a** 0.333333333333333d0*x
     &1 - 28.8052597771978d0*x14*x5))/x15**2 + 6.28318530717958d0*x14*(1
     &.804d0 - 0.804d0/x15))
      d1Ea(2) = 0.131313669051252d0*rho_a**(-0.666666666666667d0)*x10*x5
     &*x8/x13**2
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea, d1E
     &a, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = mu**2
      x3 = rho_a**(-0.666666666666667d0)*x2
      x4 = 1.18478216491348d0*x3 - 7.0d0
      x5 = exp(-0.31265084907439067d0*x3)
      x6 = x4*x5
      x7 = mu**4
      x8 = 1d0/x0
      x9 = 49.1113297378198d0*mu*rho_a**(-0.333333333333333d0) - 10.6630
     &394842213d0 *x3 + 0.233951463049512d0*x7*x8 - 81.0d0
      x10 = 1d0/x9
      x11 = 0.0519881082712438d0*x10
      x12 = gamma_aa*x1*x11*x6 + 1.0d0
      x13 = 1d0/x12
      x14 = -0.804d0*x13 + 1.804d0
      x15 = rho_a**0.333333333333333d0
      x16 = rho_a**(-4.33333333333333d0)
      x17 = x16*x2
      x18 = rho_a**(-3.66666666666667d0)*x4
      x19 = x17*x4
      x20 = mu*rho_a**(-1.66666666666667d0)
      x21 = rho_a**(-2.33333333333333d0)
      x22 = mu**3
      x23 = x21*x22
      x24 = -7.10869298948088d0*x20 + 0.311935284066016d0*x23 + 16.37044
     &32459399d0* x8
      x25 = x10*x24
      x26 = mu*x1
      x27 = x26*x4
      x28 = x25*x27
      x29 = -0.0410630556449071d0*x17 - 0.138634955389984d0*x18 + 0.0108
     &360841285172d0*x19 + 0.0519881082712438d0*x28
      x30 = x29*x5
      x31 = 1.49628538404935d0*x0
      x32 = x30*x31
      x33 = x10/x12**2
      x34 = gamma_aa*x33
      x35 = 0.277269910779967d0*x18
      x36 = x9**(-2)
      x37 = exp(-0.62530169814878134d0*x3)
      x38 = x36*x37/x12**3
      x39 = rho_a**(-5.33333333333333d0)*x2
      x40 = rho_a**(-6.0d0)*x7
      x41 = x16*x22*x25
      Ea = -1.8610514726982d0*x0*x14
      d1Ea(1) = -2.4814019635976d0*x14*x15 - x32*x34
      d1Ea(2) = -0.0777890465506374d0*x33*x6*x8
      d2Ea(1) = gamma_aa**2*x29*x31*x38*(-0.0821261112898142d0*x17 + 0.0
     &216721682570343d0*x19 + 0.103976216542488d0*x28 - x35) - 0.8271339
     &87865866d0*rho_a**(-0.666666666666667d0)*x14 - 3.99009435746494d0*
     &x15*x30*x34 - x31*x34*x5*(-mu*x25*x35 + 0.508328169763273d0*rho_a*
     &*(-4.66666666666667d0)*x4 + 0.0519881082712438d0*x1*x2*x24*x36*x4*
     &(-14.2173859789618d0*x20 + 0.623870568132033d0*x23 + 32.7408864918
     &799d0*x8) + x11*x27*( -0.727848996154038d0*rho_a**(-3.333333333333
     &33d0)*x22 - 21.8272576612533d0*x21 + 11.8478216491348d0*x26) - 0.0
     &758525888996201d0*x39*x4 + 0.28744138951435d0*x39 + 0.002258607268
     &94828d0*x4*x40 + 0.0216721682570343d0*x4*x41 - 0.0171178656172922d
     &0*x40 - 0.0821261112898142d0*x41)
      d2Ea(2) = x33*(0.155578093101275d0*gamma_aa*x10*x13*x29*x37*x4*x8 
     &- 0.103718728734183d0*x21*x6 - x32)
      d2Ea(3) = 0.00808821074878273d0*rho_a**(-4.0d0)*x38*x4**2
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea, d
     &1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**(-2.66666666666667d0)
      x1 = mu**2
      x2 = rho_a**0.666666666666667d0
      x3 = x1/x2
      x4 = exp(-0.31265084907439067d0*x3)
      x5 = gamma_aa*x4
      x6 = rho_a**0.333333333333333d0
      x7 = 1d0/x6
      x8 = mu*x7
      x9 = mu**4
      x10 = rho_a**1.33333333333333d0
      x11 = 1d0/x10
      x12 = 0.233951463049512d0*x11*x9 - 10.6630394842213d0*x3 + 49.1113
     &297378198d0* x8 - 81.0d0
      x13 = 1d0/x12
      x14 = 1.18478216491348d0*x3 - 7.0d0
      x15 = 0.0519881082712438d0*x14
      x16 = x13*x15
      x17 = x0*x16*x5 + 1.0d0
      x18 = 1d0/x17
      x19 = -0.804d0*x18 + 1.804d0
      x20 = 0.256556487706084d0*x7
      x21 = 1d0/rho_a
      x22 = 0.00844343197019482d0*x21
      x23 = x1*x22
      x24 = -x20 + x23
      x25 = exp(-15.192666241151986d0*x2/x1)
      x26 = mu*x25
      x27 = 1d0/mu
      x28 = mu**3
      x29 = x22*x28 - 0.384834731559127d0*x8 + 1.77245385090552d0*erf( 3
     &.8977770897207535d0*x27*x6)
      x30 = -x24*x26 + x29
      x31 = 0.128278243853042d0*x8
      x32 = x30*x31 - 0.375d0
      x33 = x19*x32
      x34 = 4.9628039271952d0*x10
      x35 = x31*(x26*(x20 - x23) + x29) - 0.375d0
      x36 = rho_a**0.333333333333333d0
      x37 = x19*x36
      x38 = mu*x11
      x39 = 0.00844343197019482d0*rho_a**(-2.0d0)
      x40 = rho_a**(-0.666666666666667d0)
      x41 = 0.564189583547756d0
      x42 = x25*x27
      x43 = x41*x42
      x44 = 0.0855188292353615d0*x11
      x45 = x1*x39 - x44
      x46 = rho_a**(-0.333333333333333d0)*x42
      x47 = 10.128444160768d0*x24*x46 + x26*x45 - x28*x39 + 0.1282782438
     &53042d0*x38 + 4.60575334176456d0*x40*x43
      x48 = 0.128278243853042d0*x7
      x49 = 0.0427594146176807d0*x11*x30 - x47*x48
      x50 = mu*x49
      x51 = x19*x34
      x52 = x13/x17**2
      x53 = x35*x4
      x54 = rho_a**(-4.33333333333333d0)
      x55 = x1*x54
      x56 = rho_a**(-3.66666666666667d0)*x14
      x57 = x14*x55
      x58 = mu*x0
      x59 = mu*rho_a**(-1.66666666666667d0)
      x60 = rho_a**(-2.33333333333333d0)
      x61 = x28*x60
      x62 = 16.3704432459399d0*x11 - 7.10869298948088d0*x59 + 0.31193528
     &4066016d0* x61
      x63 = x13*x62
      x64 = x58*x63
      x65 = x15*x64 - 0.0410630556449071d0*x55 - 0.138634955389984d0*x56
     & + 0.0108360841285172d0*x57
      x66 = 3.99009435746494d0*x10
      x67 = x65*x66
      x68 = x53*x67
      x69 = 0.207437457468366d0*x14*x4
      x70 = x32*x52
      x71 = x5*x70
      x72 = 0.277269910779967d0*x56
      x73 = x12**(-2)
      x74 = exp(-0.62530169814878134d0*x3)
      x75 = x32*x73*x74/x17**3
      x76 = 0.0168868639403896d0*rho_a**(-3.0d0)
      x77 = x25/x28
      x78 = rho_a**(-5.33333333333333d0)*x1
      x79 = rho_a**(-6.0d0)*x9
      x80 = x28*x54*x63
      Ea = x33*x34
      d1Ea(1) = gamma_aa*x52*x68 + 6.61707190292693d0*x35*x37 - x50*x51
      d1Ea(2) = x11*x69*x70
      d2Ea(1) = -gamma_aa**2*x67*x75*(0.103976216542488d0*x14*x64 - 0.08
     &21261112898142d0 *x55 + 0.0216721682570343d0*x57 - x72) + mu*x51*(
     & 0.0570125528235743d0*x30*x60 - x44*x47 + x48*( -0.171037658470723
     &d0*mu*x60 - 3.07050222784304d0*rho_a**( -1.66666666666667d0)*x43 -
     & 3.37614805358933d0*rho_a**( -1.33333333333333d0)*x24*x42 - 46.649
     &1155403329d0*x21*x41*x77 - 102.585381117795d0*x24*x40*x77 - x26*(x
     &1*x76 - 0.114025105647149d0*x60) + x28*x76 - 20.256888321536d0*x45
     &*x46)) - 7.98018871492988d0*x10*x5*x50*x52*x65 + 2.20569063430898d
     &0*x33* x40 + 10.6402516199065d0*x36*x65*x71 - 13.2341438058539d0*x
     &37*x50 + x66*x71*(-mu*x63*x72 + 0.508328169763273d0*rho_a**( -4.66
     &666666666667d0)*x14 + x0*x1*x15*x62*x73*(32.7408864918799d0* x11 -
     & 14.2173859789618d0*x59 + 0.623870568132033d0*x61) - 0.07585258889
     &96201d0*x14*x78 + 0.00225860726894828d0*x14*x79 + 0.02167216825703
     &43d0*x14*x80 - x16*x58*(0.727848996154038d0*rho_a **(-3.3333333333
     &3333d0)*x28 - 11.8478216491348d0*x58 + 21.8272576612533d0*x60) + 0
     &.28744138951435d0*x78 - 0.0171178656172922d0*x79 - 0.0821261112898
     &142d0*x80)
      d2Ea(2) = x52*(-0.414874914936733d0*gamma_aa*x11*x13*x14*x18*x35*x
     &65*x74 + 0.276583276624489d0*x14*x53*x60 - x38*x49*x69 + x68)
      d2Ea(3) = -0.0215685619967539d0*rho_a**(-4.0d0)*x14**2*x75
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea, d
     &1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**(-2.66666666666667d0)*x1
      x3 = 0.0519881082712438d0*x2
      x4 = rho_a**0.666666666666667d0
      x5 = 0.844037013397333d0*x4
      x6 = mu**(-4)
      x7 = rho_a**2.0d0
      x8 = 14.4026298885989d0*x6*x7
      x9 = rho_a**1.33333333333333d0
      x10 = 0.213719543995407d0*x1*x9
      x11 = x10 - x5 + x8
      x12 = x0/x4
      x13 = exp(-0.31265084907439067d0*x12)
      x14 = x11*x13
      x15 = gamma_aa*x14*x3
      x16 = x15 - 1.0d0
      x17 = 1.804d0 + 0.804d0/x16
      x18 = rho_a**0.333333333333333d0
      x19 = 1d0/x18
      x20 = 0.256556487706084d0*x19
      x21 = rho_a**1.0d0
      x22 = 1d0/x21
      x23 = 0.00844343197019482d0*x22
      x24 = x0*x23
      x25 = -x20 + x24
      x26 = exp(-15.192666241151986d0*x1*x4)
      x27 = mu*x26
      x28 = 1d0/mu
      x29 = mu**3
      x30 = mu*x19
      x31 = x23*x29 - 0.384834731559127d0*x30 + 1.77245385090552d0*erf( 
     &3.8977770897207535d0*x18*x28)
      x32 = -x25*x27 + x31
      x33 = 0.128278243853042d0*x30
      x34 = x32*x33 - 0.375d0
      x35 = x17*x34
      x36 = 4.9628039271952d0*x9
      x37 = rho_a**0.333333333333333d0
      x38 = -x15 + 1.0d0
      x39 = 1d0/x38
      x40 = -0.804d0*x39 + 1.804d0
      x41 = x33*(x27*(x20 - x24) + x31) - 0.375d0
      x42 = x38**(-2)
      x43 = -x10 + x5 - x8
      x44 = rho_a**(-4.33333333333333d0)
      x45 = 0.0108360841285172d0*x44
      x46 = rho_a**(-3.66666666666667d0)*x1
      x47 = 0.138634955389984d0*x46
      x48 = rho_a**(-0.333333333333333d0)
      x49 = 0.562691342264888d0*x48
      x50 = 28.8052597771978d0*x6
      x51 = x21*x50
      x52 = 0.284959391993876d0*x1*x37
      x53 = x3*(x49 - x51 - x52) + x43*x45 - x43*x47
      x54 = x13*x41
      x55 = 3.99009435746494d0*x9
      x56 = x53*x54*x55
      x57 = 1d0/x9
      x58 = 0.00844343197019482d0/x7
      x59 = rho_a**(-0.666666666666667d0)
      x60 = 0.564189583547756d0
      x61 = x26*x28
      x62 = x60*x61
      x63 = 0.0855188292353615d0*x57
      x64 = x0*x58 - x63
      x65 = x48*x61
      x66 = 0.128278243853042d0*mu*x57 + 10.128444160768d0*x25*x65 + x27
     &*x64 - x29* x58 + 4.60575334176456d0*x59*x62
      x67 = 0.128278243853042d0*x19
      x68 = 0.0427594146176807d0*x32*x57 - x66*x67
      x69 = mu*x36
      x70 = x1*x11
      x71 = 0.207437457468366d0*x57
      x72 = x13/x16**2
      x73 = x34*x72
      x74 = -x49 + x51 + x52
      x75 = -x11*x45 + x11*x47 - x3*x74
      x76 = gamma_aa*x73
      x77 = mu*x68
      x78 = 0.0216721682570343d0*x44
      x79 = rho_a**(-1.33333333333333d0)
      x80 = 0.277269910779967d0*x46
      x81 = exp(-0.62530169814878134d0*x12)
      x82 = x34*x81/x16**3
      x83 = rho_a**(-2.33333333333333d0)
      x84 = 0.0168868639403896d0*rho_a**(-3.0d0)
      x85 = x26/x29
      Ea = x35*x36
      d1Ea(1) = gamma_aa*x42*x56 + 6.61707190292693d0*x37*x40*x41 - x40*
     &x68*x69
      d1Ea(2) = -x70*x71*x73
      d2Ea(1) = gamma_aa**2*x55*x75*x82*(-x11*x78 + x11*x80 - 0.10397621
     &6542488d0*x2*x74 ) - 7.98018871492988d0*gamma_aa*x72*x75*x77*x9 - 
     &13.2341438058539d0*x17*x37*x77 + x17*x69*(0.0570125528235743d0* x3
     &2*x83 - x63*x66 + x67*(-0.171037658470723d0*mu*x83 - 3.07050222784
     &304d0*rho_a**(-1.66666666666667d0)*x62 - 46.6491155403329d0*x22*x6
     &0*x85 - 102.585381117795d0*x25*x59*x85 - 3.37614805358933d0*x25*x6
     &1*x79 - x27*(x0*x84 - 0.114025105647149d0*x83) + x29*x84 - 20.2568
     &88321536d0*x64*x65)) + 2.20569063430898d0*x35*x59 + 10.64025161990
     &65d0*x37*x75*x76 + x55*x76*(-0.00225860726894828d0*rho_a**(-6.0d0)
     &*x0*x11 + 0.0758525888996201d0*rho_a**(-5.33333333333333d0)*x11 - 
     &0.508328169763273d0*rho_a**(-4.66666666666667d0)*x70 - x3*( 0.0949
     &864639979585d0*x1*x59 + x50 + 0.187563780754963d0*x79) - x74*x78 +
     & x74*x80)
      d2Ea(2) = x42*(0.414874914936733d0*gamma_aa*x39*x41*x53*x57*x70*x8
     &1 + x14*x28*x68* x71 - 0.276583276624489d0*x54*x70*x83 + x56)
      d2Ea(3) = 0.0215685619967539d0*rho_a**(-4.0d0)*x11**2*x6*x82
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea, d
     &1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**0.666666666666667d0
      x3 = 15.192666241152d0*x1*x2
      x4 = exp(x3)
      x5 = 3.0d0*x4
      x6 = rho_a**1.33333333333333d0
      x7 = 1d0/x6
      x8 = mu**4
      x9 = x8*(x4 - 1.0d0)
      x10 = x7*x9
      x11 = x5 - 2.0d0
      x12 = x0/x2
      x13 = rho_a**0.333333333333333d0
      x14 = 1d0/x13
      x15 = mu*x14
      x16 = 1d0/mu
      x17 = erf(3.8977770897207535d0*x13*x16)
      x18 = x17*x4
      x19 = 0.00866486900183379d0*x10 - 0.131642462768164d0*x11*x12 + 1.
     &81893813843777d0*x15*x18 - x5
      x20 = 1d0/x19
      x21 = 0.000722072416819482d0*x10
      x22 = 7.0d0*x4 + 11.0d0
      x23 = 0.00060945584614891d0*x12*x22
      x24 = -x21 + x23 + 0.0185185185185185d0
      x25 = rho_a**(-2.0d0)
      x26 = x1*x25
      x27 = 3.15935190989552d0*x26
      x28 = x24*x27
      x29 = x20*x28
      x30 = exp(-0.31265084907439067d0*x12)
      x31 = gamma_aa*x30
      x32 = x29*x31
      x33 = x32 - 1.0d0
      x34 = 1.804d0 + 0.804d0/x33
      x35 = 0.256556487706084d0*x14
      x36 = 1d0/rho_a
      x37 = 0.00844343197019482d0*x36
      x38 = x0*x37
      x39 = -x35 + x38
      x40 = exp(-x3)
      x41 = mu*x40
      x42 = mu**3
      x43 = -0.384834731559127d0*x15 + 1.77245385090552d0*x17 + x37*x42
      x44 = -x39*x41 + x43
      x45 = 0.128278243853042d0*x15
      x46 = x44*x45 - 0.375d0
      x47 = x34*x46
      x48 = 4.9628039271952d0*x6
      x49 = rho_a**0.333333333333333d0
      x50 = x45*(x41*(x35 - x38) + x43) - 0.375d0
      x51 = -x32 + 1.0d0
      x52 = 1d0/x51
      x53 = -0.804d0*x52 + 1.804d0
      x54 = mu*x7
      x55 = 0.00844343197019482d0*x25
      x56 = rho_a**(-0.666666666666667d0)
      x57 = 0.564189583547756d0
      x58 = x16*x40
      x59 = x57*x58
      x60 = 0.0855188292353615d0*x7
      x61 = x0*x55 - x60
      x62 = rho_a**(-0.333333333333333d0)
      x63 = x58*x62
      x64 = 10.128444160768d0*x39*x63 + x41*x61 - x42*x55 + 0.1282782438
     &53042d0*x54 + 4.60575334176456d0*x56*x59
      x65 = 0.128278243853042d0*x14
      x66 = 0.0427594146176807d0*x44*x7 - x64*x65
      x67 = mu*x48
      x68 = x20/x51**2
      x69 = x21 - x23 - 0.0185185185185185d0
      x70 = rho_a**(-3.66666666666667d0)
      x71 = 0.658516038102422d0*x70
      x72 = rho_a**(-3.0d0)
      x73 = x1*x72
      x74 = 6.31870381979104d0*x73
      x75 = x36*x4
      x76 = 0.0432098765432099d0*x75
      x77 = rho_a**(-1.66666666666667d0)
      x78 = x0*x77
      x79 = x4*x78
      x80 = 0.00731347015378692d0*x79
      x81 = rho_a**(-2.33333333333333d0)
      x82 = x81*x9
      x83 = 0.000962763222425976d0*x82
      x84 = 0.000406303897432606d0*x22*x78
      x85 = x36*x57
      x86 = x1*x4
      x87 = 30.385332482304d0*x86
      x88 = x11*x78
      x89 = x18*x54
      x90 = x4*x56
      x91 = x16*x17*x90
      x92 = x62*x87 + 4.0d0*x75 - 0.087761641845443d0*x79 + 0.0115531586
     &691117d0*x82 - 4.72654360241471d0*x85 - 0.087761641845443d0*x88 + 
     &0.606312712812591d0*x89 - 18.4230133670582d0*x91
      x93 = x20*x92
      x94 = x27*x93
      x95 = x27*(-x76 + x80 - x83 + x84) + x69*x71 - x69*x74 + x69*x94
      x96 = x30*x50
      x97 = 3.99009435746494d0*x6
      x98 = x95*x96*x97
      x99 = x1*x24
      x100 = 12.6061122289202d0*x30*x56
      x101 = x20/x33**2
      x102 = x101*x46
      x103 = mu*x66
      x104 = mu*x81
      x105 = 0.0168868639403896d0*x72
      x106 = 1d0/x42
      x107 = x106*x40
      x108 = rho_a**(-1.66666666666667d0)
      x109 = rho_a**(-1.33333333333333d0)
      x110 = x76 - x80 + x83 - x84
      x111 = -x110*x27 - x24*x71 + x24*x74 - x24*x94
      x112 = x102*x31
      x113 = 1.31703207620484d0*x70
      x114 = x113*x24
      x115 = 12.6374076395821d0*x73
      x116 = x115*x24
      x117 = 6.31870381979104d0*x26
      x118 = x110*x117
      x119 = x19**(-2)
      x120 = exp(-0.62530169814878134d0*x12)
      x121 = x119*x120*x46/x33**3
      x122 = x25*x4
      x123 = rho_a**(-2.66666666666667d0)*x0
      x124 = x123*x4
      x125 = rho_a**(-3.33333333333333d0)*x9
      x126 = 1d0/x8
      x127 = x16*x18
      Ea = x47*x48
      d1Ea(1) = gamma_aa*x68*x98 + 6.61707190292693d0*x49*x50*x53 - x53*
     &x66*x67
      d1Ea(2) = -x100*x102*x99
      d2Ea(1) = gamma_aa**2*x111*x121*x97*(-x114 + x116 - x117*x24*x93 -
     & x118) - 7.98018871492988d0*x101*x103*x111*x31*x6 - 13.23414380585
     &39d0* x103*x34*x49 + 10.6402516199065d0*x111*x112*x49 + x112*x97*(
     & -0.137257065627884d0*rho_a**(-5.33333333333333d0)*x0*x24 + 3.7315
     &9088258039d0*rho_a**(-4.66666666666667d0)*x24 - 18.9561114593731d0
     &*rho_a**(-4.0d0)*x99 - x110*x113 + x110*x115 - x114*x93 + x116*x93
     & - x118*x93 - x119*x28*x92*(60.7706649646079d0 *x62*x86 + 8.0d0*x7
     &5 - 0.175523283690886d0*x79 + 0.0231063173382234d0*x82 - 9.4530872
     &0482942d0*x85 - 0.175523283690886d0*x88 + 1.21262542562518d0*x89 -
     & 36.8460267341165d0*x91) - x27*(0.43764882176158d0*x109*x86 - 0.14
     &6090534979424d0*x122 + 0.000677173162387677d0*x123*x22 + 0.0219404
     &104613607d0*x124 - 0.00224644751899394d0*x125) - x29*( -47.8725329
     &504926d0*x1*x109*x57 - 0.808416950416787d0*x104*x18 - 186.59646216
     &1332d0*x106*x17*x75 + 12.2820089113722d0*x108*x127 + x109*x87 + 0.
     &146269403075738d0*x11*x123 - 7.55555555555556d0*x122 + 0.263284925
     &536329d0*x124 - 0.0269573702279273d0*x125 + 307.756143353386d0*x12
     &6*x90 + 6.14100445568608d0*x127*x77 + 6.30205813655294d0*x25*x57))
     & + x34*x67*(0.0570125528235743d0*x44* x81 - x60*x64 + x65*(-0.1710
     &37658470723d0*x104 + x105*x42 - 102.585381117795d0*x107*x39*x56 - 
     &46.6491155403329d0*x107*x85 - 3.07050222784304d0*x108*x59 - 3.3761
     &4805358933d0*x109*x39*x58 - x41*(x0*x105 - 0.114025105647149d0*x81
     &) - 20.256888321536d0*x61* x63)) + 2.20569063430898d0*x47*x56
      d2Ea(2) = x68*(25.2122244578404d0*gamma_aa*x120*x20*x50*x52*x56*x9
     &5*x99 + x100*x16 *x24*x66 - 16.8081496385603d0*x108*x96*x99 + x98)
     &
      d2Ea(3) = 79.6542894935927d0*rho_a**(-2.66666666666667d0)*x121*x12
     &6*x24**2
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea, d1E
     &a, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**2.0d0
      x3 = rho_a**(-2.66666666666667d0)
      x4 = 0.0519881082712438d0*x3
      x5 = x1*x4
      x6 = rho_a**0.666666666666667d0
      x7 = x0/x6
      x8 = exp(-0.31265084907439067d0*x7)
      x9 = 0.844037013397333d0*x6
      x10 = mu**(-4)
      x11 = 14.4026298885989d0*x10*x2
      x12 = 0.213719543995407d0*rho_a**1.33333333333333d0*x1
      x13 = x11 + x12 - x9
      x14 = x13*x8
      x15 = gamma_aa*x14*x5
      x16 = x15 - 1.0d0
      x17 = rho_a**1.0d0
      x18 = -x15 + 1.0d0
      x19 = 1d0/x18
      x20 = rho_a**(-4.33333333333333d0)
      x21 = -x11 - x12 + x9
      x22 = x20*x21
      x23 = rho_a**(-3.66666666666667d0)
      x24 = x1*x21
      x25 = x23*x24
      x26 = 28.8052597771978d0*x10
      x27 = 0.562691342264888d0*rho_a**(-0.333333333333333d0) - 0.284959
     &391993876d0* rho_a**0.333333333333333d0*x1 - x17*x26
      x28 = x1*x27
      x29 = 0.0108360841285172d0*x22 - 0.138634955389984d0*x25 + x28*x4
      x30 = x29*x8
      x31 = 2.52584049348619d0*x2
      x32 = x30*x31
      x33 = x18**(-2)
      x34 = gamma_aa*x33
      x35 = rho_a**(-0.666666666666667d0)
      x36 = exp(-0.62530169814878134d0*x7)
      x37 = x29*x36
      Ea = -3.14159265358979d0*x1*x2*(1.804d0 + 0.804d0/x16)
      d1Ea(1) = -x1*(6.28318530717958d0*x17*(-0.804d0*x19 + 1.804d0) + x
     &32*x34)
      d1Ea(2) = 0.131313669051252d0*x10*x14*x35/x16**2
      d2Ea(1) = x1*(gamma_aa**2*x31*x37*(0.0216721682570343d0*x22 - 0.27
     &7269910779967d0* x25 + 0.103976216542488d0*x28*x3)/x18**3 - 10.103
     &3619739448d0*x17 *x30*x34 + 5.05168098697239d0*x19 - x31*x34*x8*( 
     &0.00225860726894828d0*rho_a**(-6.0d0)*x0*x21 - 0.0758525888996201d
     &0*rho_a**(-5.33333333333333d0)*x21 + 0.508328169763273d0*rho_a**(-
     &4.66666666666667d0)*x24 + 0.0216721682570343d0*x20*x27 - 0.2772699
     &10779967d0*x23*x28 - x5*( 0.187563780754963d0*rho_a**(-1.333333333
     &33333d0) + 0.0949864639979585d0*rho_a**(-0.666666666666667d0)*x1 +
     & x26)) - 11.334866294152d0)
      d2Ea(2) = x1*x33*(-0.262627338102504d0*gamma_aa*x1*x13*x19*x35*x37
     & + 0.262627338102504d0*rho_a**(-1.66666666666667d0)*x1*x14 - x32)
      d2Ea(3) = -0.0136534984882616d0*rho_a**(-3.33333333333333d0)*x13**
     &2*x36/(mu**6*x16 **3)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_PW92_ERF(rho_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18
      E = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/rho_c
      x3 = x2**0.333333333333333d0
      x4 = x2**0.166666666666667d0
      x5 = sqrt(x2)
      x6 = x2**0.666666666666667d0
      x7 = x1*(0.132568899905202d0*x3 + 1.0d0)*log(1.0d0 + 4.93480220054
     &4679d0/(x1* (2.2255694211506878d0*x3 + 5.9825504357710795d0*x4 + 0
     &.80042863499936356d0*x5 + 0.18970043257475591d0*x6)))
      x8 = mu**2
      x9 = x2**1.66666666666667d0
      x10 = x2**1.0d0
      x11 = x2**1.33333333333333d0
      x12 = (-0.00304892328706429d0*x10 + 0.000275310793616358d0*x11 + 0
     &.0128480790170175d0*x3 + 0.0315297404574782d0*x6 + 1.0d0)*exp( -0.
     &46675853320810856d0*x3)
      x13 = x0 - 1.0d0
      x14 = x3*x8
      x15 = 4.60448529760779d0*mu*x4 + 1.0d0
      x16 = mu**3
      x17 = x12 - 1.0d0
      x18 = 1d0/(0.337569861992026d0*x3 + 0.0244354823084343d0*x6 + 1.0d
     &0)
      E = -rho_c*(0.202642367284676d0*x7 - (0.000640558526213818d0*mu**8
     &*x2** 2.66666666666667d0*x7 + mu**6*(0.00248266686171384d0*x17*x9 
     &+ 0.0108059179628909d0*x2**2.0d0*x7) + 0.00264117541143784d0*mu**5
     &* x12*x9 - mu**4*(-0.0418814102677522d0*x10*x17 - 0.06835901301883
     &05d0*x11*x7 + 0.00020688890514282d0*x2**1.0d0*( -x18*(-0.305045320
     &493917d0*x3 + 17.2274040011291d0) - ( -5.62907590342551d0*x3 + 6.0
     &84d0*x6)*exp(-0.33933171852197186d0* x3) + 17.2274040011291d0)) - 
     &0.202642367284676d0*x13*log(( 0.62035049089940009d0*x14*(3.4485100
     &000000002d0 - 1.227722850784289d0/x13) + x15 + 1.9140710242289749d
     &0*x16*x5)/( 2.1392848713714905d0*x14 + x15)) + x16*(0.044555374183
     &0535d0*x10* x12 + 0.000132058770571892d0*x9*(17.2274040011291d0*x1
     &8*( -0.0177069812999059d0*x3 + 1.0d0)/x6 + 37.6991118430775d0*( -3
     &.07073492995203d0*x3 + 0.384834731559127d0*x6)*exp( -0.19230865217
     &881402d0*x3)/x10)))/(0.237113969739022d0*x6*x8 + 1.0d0)**4)
      end subroutine



C*****************************************************************************
      pure subroutine D1ESRC_PW92_ERF(rho_c, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = 0.693147180559945d0
      x4 = -x3 + 1.0d0
      x5 = x0**0.166666666666667d0
      x6 = sqrt(x0)
      x7 = x0**0.666666666666667d0
      x8 = 2.22556942115069d0*x1 + 5.98255043577108d0*x5 + 0.80042863499
     &9364d0*x6 + 0.189700432574756d0*x7
      x9 = 1.0d0 + 4.93480220054468d0/(x4*x8)
      x10 = x4*log(x9)
      x11 = x10*x2
      x12 = 0.202642367284676d0*x11
      x13 = mu**2
      x14 = x13*x7
      x15 = 0.237113969739022d0*x14 + 1.0d0
      x16 = x15**(-4)
      x17 = mu**5
      x18 = x0**1.66666666666667d0
      x19 = exp(-0.46675853320810856d0*x1)
      x20 = x0**1.0d0
      x21 = 0.00304892328706429d0*x20
      x22 = x0**1.33333333333333d0
      x23 = 0.0128480790170175d0*x1 - x21 + 0.000275310793616358d0*x22 +
     & 0.0315297404574782d0*x7 + 1.0d0
      x24 = x19*x23
      x25 = x18*x24
      x26 = 0.00264117541143784d0*x17*x25
      x27 = x0**2.66666666666667d0
      x28 = x11*x27
      x29 = 0.000640558526213818d0*mu**8*x28
      x30 = x1*x13
      x31 = 4.60448529760779d0*mu*x5 + 1.0d0
      x32 = mu**3
      x33 = x3 - 1.0d0
      x34 = 3.44851d0 - 1.22772285078429d0/x33
      x35 = 0.6203504908994d0*x30*x34 + x31 + 1.91407102422897d0*x32*x6
      x36 = x35/(2.13928487137149d0*x30 + x31)
      x37 = 0.202642367284676d0*x33
      x38 = x37*log(x36)
      x39 = x24 - 1.0d0
      x40 = x18*x39
      x41 = x0**2.0d0
      x42 = x11*x41
      x43 = mu**6*(0.00248266686171384d0*x40 + 0.0108059179628909d0*x42)
     &
      x44 = 0.0445553741830535d0*x20
      x45 = x24*x44
      x46 = exp(-0.19230865217881402d0*x1)
      x47 = x46*(-3.07073492995203d0*x1 + 0.384834731559127d0*x7)
      x48 = 37.6991118430775d0/x20
      x49 = x47*x48
      x50 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x7 + 1.0d0
      x51 = 1d0/x50
      x52 = 1d0/x7
      x53 = -0.0177069812999059d0*x1 + 1.0d0
      x54 = x51*x52*x53
      x55 = x18*(x49 + 17.2274040011291d0*x54)
      x56 = x32*(x45 + 0.000132058770571892d0*x55)
      x57 = mu**4
      x58 = 0.0418814102677522d0*x20*x39
      x59 = exp(-0.33933171852197186d0*x1)
      x60 = x0**1.0d0*(-x51*(-0.305045320493917d0*x1 + 17.2274040011291d
     &0) - x59*( -5.62907590342551d0*x1 + 6.084d0*x7) + 17.2274040011291
     &d0)
      x61 = x11*x22
      x62 = x57*(-x58 + 0.00020688890514282d0*x60 - 0.0683590130188305d0
     &*x61)
      x63 = x16*(x26 + x29 - x38 + x43 + x56 - x62)
      x64 = (0.741856473716896d0*x1 + 0.99709173929518d0*x5 + 0.40021431
     &7499682d0*x6 + 0.126466955049837d0*x7)/(x8**2*x9)
      x65 = x19*(-0.0042826930056725d0*x1 + x21 - 0.000367081058155144d0
     &*x22 - 0.0210198269716521d0*x7)
      x66 = mu**7
      x67 = x2*x64
      x68 = 0.767414216267966d0*x5
      x69 = mu*x1
      x70 = x19*(0.0777930888680181d0*x1*x23 - 0.00214134650283625d0*x1 
     &+ 0.00152446164353214d0*x20 - 0.000183540529077572d0*x22 - 0.01050
     &99134858261d0*x7)
      x71 = 1d0/x1
      x72 = 0.101681773497972d0*x51*x71
      x73 = 11.4849360007528d0*x54
      x74 = 17.2274040011291d0*x52*x53*(0.112523287330675d0*x1 + 0.01629
     &03215389562d0*x7)/x50**2
      x75 = x59*(-0.240695990468967d0*x1 + 0.26014827853397d0*x7)
      E = -rho_c*(x12 - x63)
      d1E(1) = -mu*x16*(-0.000410930386996015d0*x0**2.0d0*x24*x57 + 2.83
     &060463816878d-5 *x0**3.0d0*x10*x66 + x13*(-0.000132058770571892d0*
     &x18*( 2.41662179562687d0*x0**(-0.666666666666667d0)*x47 + x46*x48*
     &( 1.02357830998401d0*x1 - 0.256556487706084d0*x7) + x49 + x72 + x7
     &3 + x74) - 0.00693220036674015d0*x22*x24 - x44*x65 + x45 + 0.00022
     &0097950953153d0*x55) + x17*(0.00047750955226877d0*x0** 2.333333333
     &33333d0*x10 - 0.00496533372342768d0*x18*x70 + 0.0041377781028564d0
     &*x40 - 0.0533250677421794d0*x41*x67 + 0.0216118359257818d0*x42) - 
     &0.00264117541143784d0*x18*x57*x65 + 0.00440195901906306d0*x25*x57 
     &- 0.00316102962473761d0*x27*x66*x67 + 0.00170815606990352d0*x28*x6
     &6 + x32*(0.00302075971817058d0*x10* x18 + 0.00020688890514282d0*x1
     &8*(-23.3866625383245d0*x52*x59*( 0.0802319968229891d0*x1 - 0.17343
     &2185689313d0*x7) - 15.591108358883d0*x52*x75 + 11.4849360007528d0*
     &x52 - 2.64527879654103d0*x71*x75 - x72 - x73 - x74) - 0.0837628205
     &355044d0*x20*x70 - 0.337338207872387d0*x22*x67 + x58 - 0.000344814
     &8419047d0*x60 + 0.091145350691774d0*x61) + x37*( -0.95703551211448
     &7d0*x13*x6 - 0.206783496966467d0*x34*x69 + x36*( x68 + 0.713094957
     &12383d0*x69) - x68)/x35) + 0.0089546919017051d0* x1*x10 - x12 - 0.
     &632303919304059d0*x14*(-x26 - x29 + x38 - x43 - x56 + x62)/x15**5 
     &+ x63 - 0.5d0*x64*(0.265137799810404d0*x1 + 2.0d0)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_PW92_ERF(rho_c, mu, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = 0.693147180559945d0
      x4 = -x3 + 1.0d0
      x5 = 1d0/x4
      x6 = x0**0.166666666666667d0
      x7 = sqrt(x0)
      x8 = 0.800428634999364d0*x7
      x9 = x0**0.666666666666667d0
      x10 = 2.22556942115069d0*x1 + 5.98255043577108d0*x6 + x8 + 0.18970
     &0432574756d0 *x9
      x11 = 1.0d0 + 4.93480220054468d0*x5/x10
      x12 = x4*log(x11)
      x13 = x12*x2
      x14 = 0.202642367284676d0*x13
      x15 = mu**2
      x16 = x15*x9
      x17 = 0.237113969739022d0*x16 + 1.0d0
      x18 = x17**(-4)
      x19 = mu**5
      x20 = x0**1.66666666666667d0
      x21 = exp(-0.46675853320810856d0*x1)
      x22 = x0**1.0d0
      x23 = 0.00304892328706429d0*x22
      x24 = -x23
      x25 = x0**1.33333333333333d0
      x26 = 0.0128480790170175d0*x1 + x24 + 0.000275310793616358d0*x25 +
     & 0.0315297404574782d0*x9 + 1.0d0
      x27 = x21*x26
      x28 = x20*x27
      x29 = 0.00264117541143784d0*x19*x28
      x30 = x0**2.66666666666667d0
      x31 = x13*x30
      x32 = 0.000640558526213818d0*mu**8*x31
      x33 = x1*x15
      x34 = 4.60448529760779d0*mu*x6 + 1.0d0
      x35 = 2.13928487137149d0*x33 + x34
      x36 = 1d0/x35
      x37 = x3 - 1.0d0
      x38 = 1.22772285078429d0/x37
      x39 = -x38 + 3.44851d0
      x40 = 0.6203504908994d0*x33
      x41 = mu**3
      x42 = x34 + 1.91407102422897d0*x41*x7
      x43 = x39*x40 + x42
      x44 = x36*x43
      x45 = 0.202642367284676d0*x37
      x46 = x45*log(x44)
      x47 = x27 - 1.0d0
      x48 = x20*x47
      x49 = x0**2.0d0
      x50 = x13*x49
      x51 = mu**6*(0.00248266686171384d0*x48 + 0.0108059179628909d0*x50)
     &
      x52 = 0.0445553741830535d0*x22
      x53 = x27*x52
      x54 = exp(-0.19230865217881402d0*x1)
      x55 = x54*(-3.07073492995203d0*x1 + 0.384834731559127d0*x9)
      x56 = 1d0/x22
      x57 = 37.6991118430775d0*x56
      x58 = x55*x57
      x59 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x9 + 1.0d0
      x60 = 1d0/x59
      x61 = 1d0/x9
      x62 = x61*(-0.0177069812999059d0*x1 + 1.0d0)
      x63 = x60*x62
      x64 = x20*(x58 + 17.2274040011291d0*x63)
      x65 = x41*(x53 + 0.000132058770571892d0*x64)
      x66 = mu**4
      x67 = x22*x47
      x68 = 0.0418814102677522d0*x67
      x69 = x0**1.0d0
      x70 = exp(-0.33933171852197186d0*x1)
      x71 = x70*(-5.62907590342551d0*x1 + 6.084d0*x9)
      x72 = x60*(-0.305045320493917d0*x1 + 17.2274040011291d0)
      x73 = x69*(-x71 - x72 + 17.2274040011291d0)
      x74 = x13*x25
      x75 = x66*(-x68 + 0.00020688890514282d0*x73 - 0.0683590130188305d0
     &*x74)
      x76 = x29 + x32 - x46 + x51 + x65 - x75
      x77 = x18*x76
      x78 = x1*x12
      x79 = 0.265137799810404d0*x1 + 2.0d0
      x80 = 0.741856473716896d0*x1 + 0.99709173929518d0*x6 + 0.400214317
     &499682d0*x7 + 0.126466955049837d0*x9
      x81 = 1d0/x11
      x82 = x81/x10**2
      x83 = x80*x82
      x84 = x79*x83
      x85 = -x29 - x32 + x46 - x51 - x65 + x75
      x86 = x17**(-5)
      x87 = x16*x86
      x88 = x85*x87
      x89 = 0.0042826930056725d0*x1
      x90 = 0.0210198269716521d0*x9
      x91 = 0.000367081058155144d0*x25
      x92 = x21*(x23 - x89 - x90 - x91)
      x93 = x52*x92
      x94 = x25*x27
      x95 = 0.00693220036674015d0*x94
      x96 = 0.000220097950953153d0*x64
      x97 = 1.02357830998401d0*x1 - 0.256556487706084d0*x9
      x98 = x54*x57
      x99 = x0**(-0.666666666666667d0)
      x100 = x55*x99
      x101 = 1d0/x1
      x102 = 0.101681773497972d0*x101*x60
      x103 = 11.4849360007528d0*x63
      x104 = x59**(-2)
      x105 = 0.112523287330675d0*x1 + 0.0162903215389562d0*x9
      x106 = x104*x105
      x107 = 17.2274040011291d0*x62
      x108 = x106*x107
      x109 = x102 + x103 + x108
      x110 = 2.41662179562687d0*x100 + x109 + x58 + x97*x98
      x111 = 0.000132058770571892d0*x20
      x112 = x110*x111
      x113 = 11.4849360007528d0*x61
      x114 = x70*(-0.240695990468967d0*x1 + 0.26014827853397d0*x9)
      x115 = x114*x61
      x116 = 15.591108358883d0*x115
      x117 = 2.64527879654103d0*x101*x114
      x118 = x70*(0.0802319968229891d0*x1 - 0.173432185689313d0*x9)
      x119 = x118*x61
      x120 = 23.3866625383245d0*x119
      x121 = -x102 - x103 - x108 + x113 - x116 - x117 - x120
      x122 = 0.00020688890514282d0*x20
      x123 = x12*x20
      x124 = x1*x26
      x125 = -0.00214134650283625d0*x1 + 0.00152446164353214d0*x22 - 0.0
     &00183540529077572d0*x25 - 0.0105099134858261d0*x9
      x126 = x21*(0.0777930888680181d0*x124 + x125)
      x127 = x126*x22
      x128 = x2*x83
      x129 = x128*x25
      x130 = 0.00302075971817058d0*x123 - 0.0837628205355044d0*x127 - 0.
     &337338207872387d0*x129 + x68 + 0.091145350691774d0*x74
      x131 = 0.767414216267966d0*x6
      x132 = mu*x1
      x133 = x131 + 0.71309495712383d0*x132
      x134 = x15*x7
      x135 = 0.957035512114487d0*x134
      x136 = x132*x39
      x137 = 0.206783496966467d0*x136
      x138 = -x131 - x135 - x137
      x139 = x133*x44 + x138
      x140 = x45/x43
      x141 = x139*x140
      x142 = x66*x92
      x143 = x142*x20
      x144 = x0**2.33333333333333d0
      x145 = x12*x144
      x146 = x126*x20
      x147 = 0.0533250677421794d0*x49
      x148 = mu**7
      x149 = x0**3.0d0*x148
      x150 = x12*x149
      x151 = x28*x66
      x152 = x0**2.0d0
      x153 = x27*x66
      x154 = x152*x153
      x155 = x148*x31
      x156 = x148*x30
      x157 = x128*x156
      x158 = 2.83060463816878d-5*x150 + 0.00440195901906306d0*x151 - 0.0
     &00410930386996015d0*x154 + 0.00170815606990352d0*x155 - 0.00316102
     &962473761d0*x157 + x19*(-x128*x147 + 0.00047750955226877d0*x145 - 
     &0.00496533372342768d0*x146 + 0.0041377781028564d0*x48 + 0.02161183
     &59257818d0*x50)
      x159 = x141 - 0.00264117541143784d0*x143 + x15*(-x112 + x53 - x93 
     &- x95 + x96) + x158 + x41*(x121*x122 + x130 - 0.0003448148419047d0
     &*x73)
      x160 = mu*x18
      x161 = x82*(0.989141964955861d0*x1 + 1.16327369584438d0*x6 + 0.600
     &321476249523d0*x7 + 0.210778258416395d0*x9)
      x162 = 0.5d0*x79
      x163 = x5*x80**2/(x10**4*x11**2)
      x164 = x80*x81*(1.48371294743379d0*x1 + 1.99418347859036d0*x6 + x8
     & + 0.252933910099675d0*x9)/x10**3
      x165 = -x40*(x38 - 3.44851d0) + x42
      x166 = x133*x36
      x167 = 0.00264117541143784d0*x20*x66
      x168 = x21*(0.00571025734089666d0*x1 - 0.00609784657412857d0*x22 +
     & 0.000856522469028669d0*x25 + 0.0350330449527535d0*x9)
      x169 = x161*x2
      x170 = 0.00316102962473761d0*x156
      x171 = x163*x2
      x172 = x164*x2
      x173 = mu*(x131 + x135 + x137)
      x174 = 0.89531658564596d0*x6
      x175 = 0.0891107483661069d0*x22
      x176 = x54*x97
      x177 = -0.203363546995945d0*x101*x106 + x104*x107*(0.1500310497742
     &34d0*x1 + 0.0271505358982603d0*x9) - x105*x107*(0.225046574661351d
     &0*x1 + 0.0325806430779124d0*x9)/x59**3 - 22.9698720015055d0*x106*x
     &62 + 3.82831200025092d0*x63
      x178 = x21*(0.311172355472072d0*x1*x125 + 0.00285512867044833d0*x1
     & - 0.103724118490691d0*x124 + x24 + 0.000428261234514334d0*x25 + 0
     &.0121035293512547d0*x26*x9 + 0.0175165224763768d0*x9)
      x179 = 0.337338207872387d0*x25
      E = -rho_c*(x14 - x77)
      d1E(1) = -x14 - x159*x160 + x77 + 0.0089546919017051d0*x78 - 0.5d0
     &*x84 - 0.632303919304059d0*x88
      d2E(1) = x0*(0.0883792666034679d0*x1*x83 - 1.26460783860812d0*x159
     &*x41*x86*x9 - 2.0d0*x160*(-x15*(x112 - x53 + x93 + x95 - x96) + x1
     &58 + x167*x21 *(x24 + x89 + x90 + x91) + x41*(-x122*(x109 - x113 +
     & x116 + x117 + x120) + x130 + 0.0003448148419047d0*x69*(x71 + x72 
     &- 17.2274040011291d0)) - 0.202642367284676d0*x4*(x138 + x165*x166)
     &/ x165) - x160*(-mu*x141*x166 - 6.39350882283002d-5*x0** 2.3333333
     &3333333d0*x153 + x139*x173*x45/x43**2 + x140*(mu*x133* x43*(1.4261
     &8991424766d0*x132 + 1.53482843253593d0*x6)/x35**2 + 1.435553268171
     &73d0*x134 + 0.275711329288622d0*x136 - 2.0d0*x166* x173 + x174 - x
     &44*(0.950793276165107d0*x132 + x174)) - 0.000821860773992031d0*x14
     &2*x152 + 0.00880391803812612d0*x143 + 0.000279369479946145d0*x149*
     &x83 + x15*(0.000440195901906306d0* x110*x20 + x111*(-0.15491242678
     &0983d0*x0**(-0.333333333333333d0)* x55 - 1.61108119708458d0*x100 -
     & 75.398223686155d0*x176*x56 - 4.83324359125373d0*x176*x99 + x177 -
     & x98*(-1.36477107997868d0*x1 + 0.427594146176807d0*x9)) - x168*x52
     & - x175*x27 + x175*x92 - 0.0138644007334803d0*x25*x92 - 0.00107855
     &455836145d0*x28 - 0.000586927869208408d0*x64 + 0.0231073345558005d
     &0*x94) - 0.000188706975877918d0*x150 - 0.0117385573841682d0*x151 +
     & 0.00191767513931474d0*x154 - 0.00626323892297956d0*x155 + 0.01559
     &90559481421d0*x156*x171 + 0.0168588246652672d0*x157 - x167 *x168 +
     & x169*x170 - x170*x172 + x19*(0.213300270968718d0*x128*x49 + 0.004
     &71283037863406d0*x144*x83 - 0.0025467176121001d0*x145 + 0.01655111
     &24114256d0*x146 + x147*x169 - x147*x172 + 0.263148661638301d0*x171
     &*x49 - 0.00496533372342768d0*x178*x20 - 0.0110340749409504d0*x48 -
     & 0.0648355077773455d0*x50) + x41*( -0.0006896296838094d0*x121*x20 
     &+ x122*(-5.29055759308205d0*x101* x118 + 5.19703611962767d0*x115 -
     & 31.182216717766d0*x119 + x177 - 23.3866625383245d0*x61*x70*(-0.10
     &6975995763985d0*x1 + 0.289053642815522d0*x9) - 3.82831200025092d0*
     &x61 - x70*( -0.0720184066122292d0*x1 + 0.0778387062718705d0*x9)) -
     & 0.0120830388726823d0*x123 + 0.167525641071009d0*x127 + 0.89956855
     &4326366d0*x129 + x169*x179 + 1.66469733053645d0*x171* x25 - x172*x
     &179 - 0.0837628205355044d0*x178*x22 + 0.0298137034090898d0*x20*x83
     & - 0.0837628205355044d0*x67 + 0.0009195062450792d0*x73 - 0.2126724
     &84947473d0*x74)) + x161*x162 - x162*x164 + 2.46740110027234d0*x163
     &*x79 + 1.26460783860812d0* x76*x87 + 0.00596979460113673d0*x78 - x
     &84 + 1.05383986550676d0* x88 - 0.499760307959092d0*x25*x66*x85/x17
     &**6)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35
      E = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/rho_c
      x3 = x2**0.333333333333333d0
      x4 = 0.132568899905202d0*x3 + 1.0d0
      x5 = x2**0.166666666666667d0
      x6 = sqrt(x2)
      x7 = x2**0.666666666666667d0
      x8 = log(1.0d0 + 4.934802200544679d0/(x1*(2.2255694211506878d0*x3 
     &+ 5.9825504357710795d0*x5 + 0.80042863499936356d0*x6 + 0.189700432
     &57475591d0*x7)))
      x9 = x1*x4*x8
      x10 = 0.202642367284676d0*x9
      x11 = mu**2
      x12 = (0.237113969739022d0*x11*x7 + 1.0d0)**(-4)
      x13 = x2**1.66666666666667d0
      x14 = x2**1.0d0
      x15 = x2**1.33333333333333d0
      x16 = (-0.00304892328706429d0*x14 + 0.000275310793616358d0*x15 + 0
     &.0128480790170175d0*x3 + 0.0315297404574782d0*x7 + 1.0d0)*exp( -0.
     &46675853320810856d0*x3)
      x17 = mu**5*x13*x16
      x18 = 0.00264117541143784d0*x17
      x19 = mu**8*x2**2.66666666666667d0*x9
      x20 = 0.000640558526213818d0*x19
      x21 = x0 - 1.0d0
      x22 = x11*x3
      x23 = 4.60448529760779d0*mu*x5 + 1.0d0
      x24 = mu**3
      x25 = x21*log((0.62035049089940009d0*x22*(3.4485100000000002d0 - 1
     &.227722850784289d0/x21) + x23 + 1.9140710242289749d0*x24*x6)/( 2.1
     &392848713714905d0*x22 + x23))
      x26 = 0.202642367284676d0*x25
      x27 = x16 - 1.0d0
      x28 = mu**6*(0.00248266686171384d0*x13*x27 + 0.0108059179628909d0*
     &x2**2.0d0*x9 )
      x29 = 1d0/(0.337569861992026d0*x3 + 0.0244354823084343d0*x7 + 1.0d
     &0)
      x30 = x24*(0.000132058770571892d0*x13*(17.2274040011291d0*x29*( -0
     &.0177069812999059d0*x3 + 1.0d0)/x7 + 37.6991118430775d0*( -3.07073
     &492995203d0*x3 + 0.384834731559127d0*x7)*exp( -0.19230865217881402
     &d0*x3)/x14) + 0.0445553741830535d0*x14*x16)
      x31 = mu**4*(-0.0418814102677522d0*x14*x27 - 0.0683590130188305d0*
     &x15*x9 + 0.00020688890514282d0*x2**1.0d0*(-x29*(-0.305045320493917
     &d0*x3 + 17.2274040011291d0) - (-5.62907590342551d0*x3 + 6.084d0*x7
     &)*exp( -0.33933171852197186d0*x3) + 17.2274040011291d0))
      x32 = exp(-x12*(0.084949837941456843d0*x17 + 0.020602699373253298d
     &0*x19 - 6.5177179017939446d0*x25 + 32.163648644302206d0*x28 + 32.1
     &63648644302206d0*x30 - 32.163648644302206d0*x31) + 6.5177179017939
     &446d0*x9) - 1.0d0
      x33 = (x10 + x12*(-x18 - x20 + x26 - x28 - x30 + x31))/(x4*x8)
      x34 = 11.5212611093509d0*gamma_cc*rho_c**(-2.33333333333333d0)*x1*
     &*(-2.78d0)* x33**2.78d0
      x35 = 1.0d0 + x34/x32
      E = -rho_c*(x10 - x12*(x18 + x20 - x26 + x28 + x30 - x31) - 0.0310
     &91d0*log( x34*x35/(132.73945754984086d0*gamma_cc**2*rho_c**( -4.66
     &6666666666667d0)*x1**(-5.5599999999999996d0)*x33** 5.5599999999999
     &996d0/x32**2 + x35) + 1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 6.51771790179394d0*x1
      x3 = 1d0/rho_c
      x4 = x3**0.333333333333333d0
      x5 = 0.132568899905202d0*x4 + 1.0d0
      x6 = 1d0/x1
      x7 = x3**0.166666666666667d0
      x8 = sqrt(x3)
      x9 = x3**0.666666666666667d0
      x10 = 2.22556942115069d0*x4 + 5.98255043577108d0*x7 + 0.8004286349
     &99364d0*x8 + 0.189700432574756d0*x9
      x11 = 1.0d0 + 4.93480220054468d0*x6/x10
      x12 = log(x11)
      x13 = x12*x5
      x14 = x13*x2
      x15 = mu**2
      x16 = x15*x9
      x17 = 0.237113969739022d0*x16 + 1.0d0
      x18 = x17**(-4)
      x19 = mu**4
      x20 = x3**1.0d0
      x21 = exp(-0.46675853320810856d0*x4)
      x22 = 0.00304892328706429d0*x20
      x23 = x3**1.33333333333333d0
      x24 = -x22 + 0.000275310793616358d0*x23 + 0.0128480790170175d0*x4 
     &+ 0.0315297404574782d0*x9 + 1.0d0
      x25 = x21*x24
      x26 = x25 - 1.0d0
      x27 = 0.0418814102677522d0*x20*x26
      x28 = exp(-0.33933171852197186d0*x4)
      x29 = 0.337569861992026d0*x4 + 0.0244354823084343d0*x9 + 1.0d0
      x30 = 1d0/x29
      x31 = x3**1.0d0*(-x28*(-5.62907590342551d0*x4 + 6.084d0*x9) - x30*
     &( -0.305045320493917d0*x4 + 17.2274040011291d0) + 17.2274040011291
     &d0)
      x32 = 0.00020688890514282d0*x31
      x33 = x1*x13
      x34 = x23*x33
      x35 = 0.0683590130188305d0*x34
      x36 = x19*(-x27 + x32 - x35)
      x37 = 32.1636486443022d0*x36
      x38 = x0 - 1.0d0
      x39 = x15*x4
      x40 = 4.60448529760779d0*mu*x7 + 1.0d0
      x41 = 1d0/(2.13928487137149d0*x39 + x40)
      x42 = 1.22772285078429d0/x38
      x43 = -x42 + 3.44851d0
      x44 = 0.6203504908994d0*x39
      x45 = mu**3
      x46 = x40 + 1.91407102422897d0*x45*x8
      x47 = x43*x44 + x46
      x48 = x41*x47
      x49 = x38*log(x48)
      x50 = 6.51771790179394d0*x49
      x51 = 0.0445553741830535d0*x20
      x52 = x25*x51
      x53 = x3**1.66666666666667d0
      x54 = exp(-0.19230865217881402d0*x4)
      x55 = x54*(-3.07073492995203d0*x4 + 0.384834731559127d0*x9)
      x56 = 37.6991118430775d0/x20
      x57 = x55*x56
      x58 = 1d0/x9
      x59 = -0.0177069812999059d0*x4 + 1.0d0
      x60 = x30*x58*x59
      x61 = x53*(x57 + 17.2274040011291d0*x60)
      x62 = x45*(x52 + 0.000132058770571892d0*x61)
      x63 = 32.1636486443022d0*x62
      x64 = x26*x53
      x65 = x3**2.0d0
      x66 = x33*x65
      x67 = mu**6*(0.00248266686171384d0*x64 + 0.0108059179628909d0*x66)
     &
      x68 = 32.1636486443022d0*x67
      x69 = mu**5
      x70 = x25*x53
      x71 = x69*x70
      x72 = 0.0849498379414568d0*x71
      x73 = x3**2.66666666666667d0
      x74 = x33*x73
      x75 = mu**8*x74
      x76 = 0.0206026993732533d0*x75
      x77 = x63 + x68 + x72 + x76
      x78 = exp(x14 - x18*(-x37 - x50 + x77)) - 1.0d0
      x79 = rho_c**(-4.66666666666667d0)
      x80 = 0.202642367284676d0*x33
      x81 = 0.00264117541143784d0*x71
      x82 = 0.000640558526213818d0*x75
      x83 = 0.202642367284676d0*x49
      x84 = x36 - x62 - x67 - x81 - x82 + x83
      x85 = x18*x84 + x80
      x86 = 1d0/x5
      x87 = 1d0/x12
      x88 = x86*x87
      x89 = x85*x88
      x90 = x1**(-5.56d0)
      x91 = x89**5.56d0*x90
      x92 = x79*x91
      x93 = 132.739457549841d0*gamma_cc**2
      x94 = x92*x93
      x95 = 1d0/x78
      x96 = x1**(-2.78d0)
      x97 = x89**2.78d0*x96
      x98 = 11.5212611093509d0*rho_c**(-2.33333333333333d0)
      x99 = x97*x98
      x100 = gamma_cc*x99
      x101 = x100*x95 + 1.0d0
      x102 = 0.031091d0*log(x100*x101/(x101 + x94/x78**2) + 1.0d0)
      x103 = x18*(-x36 + x62 + x67 + x81 + x82 - x83)
      x104 = -x103 + x80
      x105 = x1*x12
      x106 = x105*x4
      x107 = (0.741856473716896d0*x4 + 0.99709173929518d0*x7 + 0.4002143
     &17499682d0*x8 + 0.126466955049837d0*x9)/(x10**2*x11)
      x108 = x107*(0.265137799810404d0*x4 + 2.0d0)
      x109 = x16*x84/x17**5
      x110 = x21*(x22 - 0.000367081058155144d0*x23 - 0.0042826930056725d
     &0*x4 - 0.0210198269716521d0*x9)
      x111 = mu**7
      x112 = x107*x5
      x113 = 0.767414216267966d0*x7
      x114 = mu*x4
      x115 = x21*(0.00152446164353214d0*x20 - 0.000183540529077572d0*x23
     & + 0.0777930888680181d0*x24*x4 - 0.00214134650283625d0*x4 - 0.0105
     &099134858261d0*x9)
      x116 = 1d0/x4
      x117 = 0.101681773497972d0*x116*x30
      x118 = 11.4849360007528d0*x60
      x119 = 17.2274040011291d0*x58*x59*(0.112523287330675d0*x4 + 0.0162
     &903215389562d0*x9)/x29**2
      x120 = x28*(-0.240695990468967d0*x4 + 0.26014827853397d0*x9)
      x121 = mu*x18*(2.83060463816878d-5*x105*x111*x3**3.0d0 - 0.0026411
     &7541143784d0* x110*x19*x53 - 0.00316102962473761d0*x111*x112*x73 +
     & 0.00170815606990352d0*x111*x74 + x15*(-x110*x51 - 0.0069322003667
     &4015d0*x23*x25 + x52 - 0.000132058770571892d0*x53* (x117 + x118 + 
     &x119 + 2.41662179562687d0*x3**( -0.666666666666667d0)*x55 + x54*x5
     &6*(1.02357830998401d0*x4 - 0.256556487706084d0*x9) + x57) + 0.0002
     &20097950953153d0*x61) - 0.000410930386996015d0*x19*x25*x3**2.0d0 +
     & 0.00440195901906306d0* x19*x70 + 0.202642367284676d0*x38*(-x113 -
     & 0.206783496966467d0* x114*x43 - 0.957035512114487d0*x15*x8 + x48*
     &(x113 + 0.71309495712383d0*x114))/x47 + x45*(0.00302075971817058d0
     &*x105* x53 - 0.337338207872387d0*x112*x23 - 0.0837628205355044d0*x
     &115* x20 + x27 - 0.0003448148419047d0*x31 + 0.091145350691774d0*x3
     &4 + 0.00020688890514282d0*x53*(-2.64527879654103d0*x116*x120 - x11
     &7 - x118 - x119 - 15.591108358883d0*x120*x58 - 23.3866625383245d0*
     &x28 *x58*(0.0802319968229891d0*x4 - 0.173432185689313d0*x9) + 11.4
     &849360007528d0*x58)) + x69*(0.00047750955226877d0*x105*x3** 2.3333
     &3333333333d0 - 0.0533250677421794d0*x112*x65 - 0.00496533372342768
     &d0*x115*x53 + 0.0041377781028564d0*x64 + 0.0216118359257818d0*x66)
     &)
      x122 = exp(x14 + x18*(x37 + x50 - x63 - x68 - x72 - x76))
      x123 = x122 - 1.0d0
      x124 = 1d0/x123
      x125 = x100*x124 + 1.0d0
      x126 = x123**(-2)
      x127 = gamma_cc/(x125 + x126*x94)
      x128 = x125*x127*x99
      x129 = rho_c**(-3.33333333333333d0)*x97
      x130 = 26.8829425884854d0*x129
      x131 = 1d0/x85
      x132 = x4*x85*x86
      x133 = x107*x6*x85*x87
      x134 = 5.76063055467543d0*x129*x131*(-0.0497880869734803d0*x106 + 
     &2.78d0*x108 + 3.51560979133057d0*x109 + 5.56d0*x121 + 0.2456943611
     &57641d0*x132 - 27.4375002350284d0*x133)
      x135 = x122*(-0.576031128088842d0*x106 + x107*(8.52779903542516d0*
     &x4 + 64.3272972886044d0) + 40.6744021938219d0*x109 + 64.3272972886
     &044d0*x121)
      x136 = 5.76063055467543d0*x124*x135
      x137 = gamma_cc*rho_c**(-5.66666666666667d0)*x91
      x138 = x124*x137
      x139 = exp(x14 - x18*(32.163648644302206d0*x19*(x27 - x32 + x35) +
     & x2*log(x41*( -x44*(x42 - 3.4485100000000002d0) + x46)) + x77)) - 
     &1.0d0
      x140 = x104*x88
      x141 = x140**5.56d0*x79*x90
      x142 = x140**2.78d0*x96*x98
      x143 = gamma_cc/x139
      x144 = x142*x143 + 1.0d0
      x145 = 1d0/(x144 + x141*x93/x139**2)
      x146 = x142*x144
      x147 = x145*x146
      x148 = gamma_cc*x95
      E = -rho_c*(-x102 + x104)
      d1E(1) = -0.5d0*rho_c*(-0.0179093838034102d0*x106*x3 + x108*x3 + 1
     &.26460783860812d0*x109*x3 + 2.0d0*x121*x3 - 0.062182d0*x127*( x124
     &*x128*(132.739457549841d0*x126*x135*x137 + x129*x136 + x130 - 66.3
     &697287749204d0*x131*x138*(-0.0995761739469607d0*x106 + 5.56d0 *x10
     &8 + 7.03121958266113d0*x109 + 11.12d0*x121 + 0.491388722315281d0*x
     &132 - 54.8750004700568d0*x133) - x134 + 619.450801899257d0*x138) -
     & x125*x130 + x125*x134 + 11.5212611093509d0*x138*(x131*(-0.2868107
     &75078269d0*x106 + 16.0145529419977d0*x108 + 20.252129182255d0*x109
     & + 32.0291058839954d0*x121 + 1.41535444399617d0*x132 - 158.0573021
     &97819d0*x133) - x136 - 26.8829425884854d0))/(x128 + 1.0d0)) + x102
     & + x103 - x80
      d1E(3) = 0.031091d0*rho_c*x145*(132.739457549841d0*x141*x143 + x14
     &6 - x147*x148*( 265.478915099682d0*x148*x92 + x99))/(gamma_cc*x147
     & + 1.0d0)
      end subroutine



C*****************************************************************************
      pure subroutine D2ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E, d2
     &E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 12, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399, x400, x401
     &, x402, x403, x404, x405, x406, x407, x408, x409, x410, x411, x412
     &, x413, x414, x415, x416, x417, x418, x419, x420, x421, x422, x423
     &, x424, x425, x426, x427, x428, x429, x430, x431, x432, x433, x434
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446, x447, x448, x449, x450, x451, x452, x453, x454, x455, x456
     &, x457, x458, x459, x460, x461, x462, x463, x464, x465, x466, x467
     &, x468, x469, x470, x471, x472, x473, x474, x475, x476, x477, x478
     &, x479, x480, x481, x482, x483, x484, x485, x486, x487, x488, x489
     &, x490, x491, x492, x493, x494, x495, x496, x497, x498, x499, x500
     &, x501, x502, x503
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 6.51771790179394d0*x1
      x3 = 1d0/rho_c
      x4 = x3**0.333333333333333d0
      x5 = 0.132568899905202d0*x4 + 1.0d0
      x6 = x3**0.166666666666667d0
      x7 = sqrt(x3)
      x8 = 0.800428634999364d0*x7
      x9 = x3**0.666666666666667d0
      x10 = 2.22556942115069d0*x4 + 5.98255043577108d0*x6 + x8 + 0.18970
     &0432574756d0 *x9
      x11 = 1d0/x1
      x12 = 4.93480220054468d0*x11
      x13 = 1.0d0 + x12/x10
      x14 = log(x13)
      x15 = x14*x5
      x16 = x15*x2
      x17 = mu**2
      x18 = x17*x9
      x19 = 0.237113969739022d0*x18 + 1.0d0
      x20 = x19**(-4)
      x21 = mu**4
      x22 = x3**1.0d0
      x23 = exp(-0.46675853320810856d0*x4)
      x24 = 0.00304892328706429d0*x22
      x25 = -x24
      x26 = x3**1.33333333333333d0
      x27 = x25 + 0.000275310793616358d0*x26 + 0.0128480790170175d0*x4 +
     & 0.0315297404574782d0*x9 + 1.0d0
      x28 = x23*x27
      x29 = x28 - 1.0d0
      x30 = x22*x29
      x31 = 0.0418814102677522d0*x30
      x32 = x3**1.0d0
      x33 = exp(-0.33933171852197186d0*x4)
      x34 = x33*(-5.62907590342551d0*x4 + 6.084d0*x9)
      x35 = 0.337569861992026d0*x4 + 0.0244354823084343d0*x9 + 1.0d0
      x36 = 1d0/x35
      x37 = x36*(-0.305045320493917d0*x4 + 17.2274040011291d0)
      x38 = x32*(-x34 - x37 + 17.2274040011291d0)
      x39 = 0.00020688890514282d0*x38
      x40 = x1*x15
      x41 = x26*x40
      x42 = 0.0683590130188305d0*x41
      x43 = x21*(-x31 + x39 - x42)
      x44 = 32.1636486443022d0*x43
      x45 = x0 - 1.0d0
      x46 = x17*x4
      x47 = 4.60448529760779d0*mu*x6 + 1.0d0
      x48 = 2.13928487137149d0*x46 + x47
      x49 = 1d0/x48
      x50 = 1.22772285078429d0/x45
      x51 = -x50 + 3.44851d0
      x52 = 0.6203504908994d0*x46
      x53 = mu**3
      x54 = x47 + 1.91407102422897d0*x53*x7
      x55 = x51*x52 + x54
      x56 = x49*x55
      x57 = x45*log(x56)
      x58 = 6.51771790179394d0*x57
      x59 = 0.0445553741830535d0*x22
      x60 = x28*x59
      x61 = x3**1.66666666666667d0
      x62 = exp(-0.19230865217881402d0*x4)
      x63 = x62*(-3.07073492995203d0*x4 + 0.384834731559127d0*x9)
      x64 = 1d0/x22
      x65 = 37.6991118430775d0*x64
      x66 = x63*x65
      x67 = 1d0/x9
      x68 = x67*(-0.0177069812999059d0*x4 + 1.0d0)
      x69 = x36*x68
      x70 = x61*(x66 + 17.2274040011291d0*x69)
      x71 = x53*(x60 + 0.000132058770571892d0*x70)
      x72 = 32.1636486443022d0*x71
      x73 = x29*x61
      x74 = x3**2.0d0
      x75 = x40*x74
      x76 = mu**6*(0.00248266686171384d0*x73 + 0.0108059179628909d0*x75)
     &
      x77 = 32.1636486443022d0*x76
      x78 = mu**5
      x79 = x28*x61
      x80 = x78*x79
      x81 = 0.0849498379414568d0*x80
      x82 = x3**2.66666666666667d0
      x83 = x40*x82
      x84 = mu**8*x83
      x85 = 0.0206026993732533d0*x84
      x86 = x72 + x77 + x81 + x85
      x87 = exp(x16 - x20*(-x44 - x58 + x86))
      x88 = x87 - 1.0d0
      x89 = 1d0/x88
      x90 = rho_c**(-2.33333333333333d0)
      x91 = x1**(-2.78d0)
      x92 = 0.202642367284676d0*x1
      x93 = x15*x92
      x94 = 0.00264117541143784d0*x80
      x95 = 0.000640558526213818d0*x84
      x96 = 0.202642367284676d0*x57
      x97 = x43 - x71 - x76 - x94 - x95 + x96
      x98 = x20*x97 + x93
      x99 = 1d0/x5
      x100 = 1d0/x14
      x101 = x100*x99
      x102 = x101*x98
      x103 = x102**2.78d0
      x104 = x103*x91
      x105 = x104*x90
      x106 = 11.5212611093509d0*x105
      x107 = gamma_cc*x106
      x108 = x107*x89 + 1.0d0
      x109 = x106*x108
      x110 = x88**(-2)
      x111 = rho_c**(-4.66666666666667d0)
      x112 = x102**5.56d0
      x113 = x1**(-5.56d0)
      x114 = x112*x113
      x115 = x111*x114
      x116 = gamma_cc**2
      x117 = 132.739457549841d0*x116
      x118 = x115*x117
      x119 = x108 + x110*x118
      x120 = 1d0/x119
      x121 = gamma_cc*x120
      x122 = x109*x121
      x123 = x122 + 1.0d0
      x124 = 0.031091d0*log(x123)
      x125 = -x43 + x71 + x76 + x94 + x95 - x96
      x126 = x125*x20
      x127 = -x126 + x93
      x128 = x1*x14
      x129 = x128*x4
      x130 = 0.0179093838034102d0*x129
      x131 = x130*x3
      x132 = 0.265137799810404d0*x4 + 2.0d0
      x133 = 0.741856473716896d0*x4 + 0.99709173929518d0*x6 + 0.40021431
     &7499682d0*x7 + 0.126466955049837d0*x9
      x134 = 1d0/x13
      x135 = x134/x10**2
      x136 = x133*x135
      x137 = x132*x136
      x138 = x137*x3
      x139 = x19**(-5)
      x140 = x139*x18
      x141 = x140*x97
      x142 = 1.26460783860812d0*x141
      x143 = 0.0042826930056725d0*x4
      x144 = 0.0210198269716521d0*x9
      x145 = 0.000367081058155144d0*x26
      x146 = -x143 - x144 - x145 + x24
      x147 = x23*x59
      x148 = x146*x147
      x149 = x26*x28
      x150 = 0.00693220036674015d0*x149
      x151 = 0.000220097950953153d0*x70
      x152 = 1.02357830998401d0*x4 - 0.256556487706084d0*x9
      x153 = x62*x65
      x154 = x3**(-0.666666666666667d0)
      x155 = x154*x63
      x156 = 1d0/x4
      x157 = 0.101681773497972d0*x156*x36
      x158 = 11.4849360007528d0*x69
      x159 = x35**(-2)
      x160 = 0.112523287330675d0*x4 + 0.0162903215389562d0*x9
      x161 = x159*x160
      x162 = x161*x68
      x163 = 17.2274040011291d0*x162
      x164 = x157 + x158 + x163
      x165 = x152*x153 + 2.41662179562687d0*x155 + x164 + x66
      x166 = 0.000132058770571892d0*x61
      x167 = x165*x166
      x168 = 11.4849360007528d0*x67
      x169 = x33*(-0.240695990468967d0*x4 + 0.26014827853397d0*x9)
      x170 = x169*x67
      x171 = 15.591108358883d0*x170
      x172 = 2.64527879654103d0*x156*x169
      x173 = x33*(0.0802319968229891d0*x4 - 0.173432185689313d0*x9)
      x174 = x173*x67
      x175 = 23.3866625383245d0*x174
      x176 = -x157 - x158 - x163 + x168 - x171 - x172 - x175
      x177 = 0.00020688890514282d0*x61
      x178 = x128*x61
      x179 = x27*x4
      x180 = 0.00214134650283625d0*x4
      x181 = 0.0105099134858261d0*x9
      x182 = 0.00152446164353214d0*x22
      x183 = 0.000183540529077572d0*x26
      x184 = -x180 - x181 + x182 - x183
      x185 = x23*(0.0777930888680181d0*x179 + x184)
      x186 = x185*x22
      x187 = x136*x5
      x188 = x187*x26
      x189 = 0.00302075971817058d0*x178 - 0.0837628205355044d0*x186 - 0.
     &337338207872387d0*x188 + x31 + 0.091145350691774d0*x41
      x190 = 0.767414216267966d0*x6
      x191 = mu*x4
      x192 = x190 + 0.71309495712383d0*x191
      x193 = x17*x7
      x194 = 0.957035512114487d0*x193
      x195 = 0.206783496966467d0*x191
      x196 = x195*x51
      x197 = -x190 - x194 - x196
      x198 = x192*x56 + x197
      x199 = 0.202642367284676d0*x45/x55
      x200 = x198*x199
      x201 = x146*x23
      x202 = x21*x61
      x203 = 0.00264117541143784d0*x202
      x204 = x3**2.33333333333333d0
      x205 = x128*x204
      x206 = 0.00496533372342768d0*x61
      x207 = 0.0533250677421794d0*x74
      x208 = mu**7
      x209 = x208*x3**3.0d0
      x210 = x128*x209
      x211 = x21*x79
      x212 = x21*x3**2.0d0
      x213 = x212*x28
      x214 = x208*x83
      x215 = x208*x82
      x216 = x187*x215
      x217 = 2.83060463816878d-5*x210 + 0.00440195901906306d0*x211 - 0.0
     &00410930386996015d0*x213 + 0.00170815606990352d0*x214 - 0.00316102
     &962473761d0*x216 + x78*(-x185*x206 - x187*x207 + 0.000477509552268
     &77d0*x205 + 0.0041377781028564d0*x73 + 0.0216118359257818d0*x75)
      x218 = x17*(-x148 - x150 + x151 - x167 + x60) + x200 - x201*x203 +
     & x217 + x53*( x176*x177 + x189 - 0.0003448148419047d0*x38)
      x219 = mu*x20
      x220 = x218*x219
      x221 = 2.0d0*x220
      x222 = exp(x16 + x20*(x44 + x58 - x72 - x77 - x81 - x85))
      x223 = x222 - 1.0d0
      x224 = x223**(-2)
      x225 = 1d0/x223
      x226 = x107*x225 + 1.0d0
      x227 = x118*x224 + x226
      x228 = 1d0/x227
      x229 = gamma_cc*x228
      x230 = x106*x226
      x231 = x229*x230
      x232 = rho_c**(-3.33333333333333d0)*x104
      x233 = 26.8829425884854d0*x232
      x234 = x226*x233
      x235 = 1d0/x98
      x236 = 0.0497880869734803d0*x129
      x237 = 2.78d0*x137
      x238 = 3.51560979133057d0*x141
      x239 = x4*x99
      x240 = x239*x98
      x241 = 0.245694361157641d0*x240
      x242 = x100*x98
      x243 = x11*x242
      x244 = x136*x243
      x245 = 27.4375002350284d0*x244
      x246 = 5.56d0*x220
      x247 = -x236 + x237 + x238 + x241 - x245 + x246
      x248 = x235*x247
      x249 = 5.76063055467543d0*x232
      x250 = x248*x249
      x251 = x226*x250
      x252 = 8.52779903542516d0*x4 + 64.3272972886044d0
      x253 = -0.576031128088842d0*x129 + x136*x252
      x254 = 40.6744021938219d0*x141 + 64.3272972886044d0*x220 + x253
      x255 = x222*x254
      x256 = 5.76063055467543d0*x225
      x257 = x255*x256
      x258 = x235*(-0.286810775078269d0*x129 + 16.0145529419977d0*x137 +
     & 20.252129182255d0*x141 + 32.0291058839954d0*x220 + 1.415354443996
     &17d0*x240 - 158.057302197819d0*x244)
      x259 = gamma_cc*x225
      x260 = rho_c**(-5.66666666666667d0)*x114
      x261 = x259*x260
      x262 = 11.5212611093509d0*x261
      x263 = x262*(-x257 + x258 - 26.8829425884854d0)
      x264 = 619.450801899257d0*x260
      x265 = gamma_cc*x260
      x266 = 132.739457549841d0*x265
      x267 = x224*x255
      x268 = 0.0995761739469607d0*x129
      x269 = 5.56d0*x137
      x270 = 7.03121958266113d0*x141
      x271 = 0.491388722315281d0*x240
      x272 = 54.8750004700568d0*x244
      x273 = 11.12d0*x220
      x274 = -x268 + x269 + x270 + x271 - x272 + x273
      x275 = 66.3697287749204d0*x235
      x276 = x225*x231
      x277 = x276*(x232*x257 + x233 - x250 + x259*x264 - x261*x274*x275 
     &+ x266*x267)
      x278 = 0.062182d0*x229
      x279 = 0.5d0*rho_c
      x280 = x50 - 3.44851d0
      x281 = -x280*x52 + x54
      x282 = x281*x49
      x283 = exp(x16 - x20*(x2*log(x282) + 32.163648644302206d0*x21*(x31
     & - x39 + x42 ) + x86)) - 1.0d0
      x284 = x101*x127
      x285 = x111*x113*x284**5.56d0
      x286 = 11.5212611093509d0*x91
      x287 = x284**2.78d0*x286*x90
      x288 = gamma_cc/x283
      x289 = x287*x288 + 1.0d0
      x290 = 1d0/(x117*x285/x283**2 + x289)
      x291 = x287*x289
      x292 = gamma_cc*x89
      x293 = 265.478915099682d0*x115
      x294 = x290*x291
      x295 = x290*(132.739457549841d0*x285*x288 + x291 - x292*x294*(x106
     & + x292*x293 ))
      x296 = 1d0/(gamma_cc*x294 + 1.0d0)
      x297 = 0.031091d0*x296
      x298 = x125*x140
      x299 = 1.26460783860812d0*x298
      x300 = x32*(x34 + x37 - 17.2274040011291d0)
      x301 = x164 - x168 + x171 + x172 + x175
      x302 = x192*x282 + x197
      x303 = x92/x281
      x304 = x302*x303
      x305 = x23*(x143 + x144 + x145 + x25)
      x306 = -x17*(x148 + x150 - x151 + x167 - x60) + x203*x305 + x217 -
     & x304 + x53*( -x177*x301 + x189 + 0.0003448148419047d0*x300)
      x307 = x219*x306
      x308 = 2.0d0*x307
      x309 = 1d0/x123
      x310 = x234 - x251 - x263 - x277
      x311 = rho_c**(-2)
      x312 = x136*x4
      x313 = x135*(0.989141964955861d0*x4 + 1.16327369584438d0*x6 + 0.60
     &0321476249523d0*x7 + 0.210778258416395d0*x9)
      x314 = x132*x313
      x315 = x133**2/(x10**4*x13**2)
      x316 = x133*x134*(1.48371294743379d0*x4 + 1.99418347859036d0*x6 + 
     &x8 + 0.252933910099675d0*x9)/x10**3
      x317 = x132*x316
      x318 = x21*x26/x19**6
      x319 = x125*x318
      x320 = x139*x53*x9
      x321 = x306*x320
      x322 = 0.0006896296838094d0*x61
      x323 = 3.82831200025092d0*x67
      x324 = x33*(-0.0720184066122292d0*x4 + 0.0778387062718705d0*x9)
      x325 = 5.19703611962767d0*x170
      x326 = 23.3866625383245d0*x33*x67*(-0.106975995763985d0*x4 + 0.289
     &053642815522d0*x9)
      x327 = 31.182216717766d0*x174
      x328 = 5.29055759308205d0*x156*x173
      x329 = 3.82831200025092d0*x69
      x330 = 0.203363546995945d0*x156*x161
      x331 = 22.9698720015055d0*x162
      x332 = 17.2274040011291d0*x68
      x333 = x159*x332*(0.150031049774234d0*x4 + 0.0271505358982603d0*x9
     &)
      x334 = x160*x332*(0.225046574661351d0*x4 + 0.0325806430779124d0*x9
     &)/x35**3
      x335 = 0.311172355472072d0*x4
      x336 = -0.103724118490691d0*x179 + x25 + 0.000428261234514334d0*x2
     &6 + 0.0121035293512547d0*x27*x9 + 0.00285512867044833d0*x4 + 0.017
     &5165224763768d0*x9
      x337 = x22*x23
      x338 = 0.0837628205355044d0*x337
      x339 = x313*x5
      x340 = 0.337338207872387d0*x26
      x341 = x11*x315
      x342 = x341*x5
      x343 = x316*x5
      x344 = 0.0298137034090898d0*x136*x61 - 0.0120830388726823d0*x178 +
     & 0.167525641071009d0*x186 + 0.899568554326366d0*x188 + 1.664697330
     &53645d0*x26*x342 - 0.0837628205355044d0*x30 + x339* x340 - x340*x3
     &43 - 0.212672484947473d0*x41
      x345 = 0.0891107483661069d0*x146*x337
      x346 = 0.0138644007334803d0*x201*x26
      x347 = -0.00609784657412857d0*x22 + 0.000856522469028669d0*x26 + 0
     &.00571025734089666d0*x4 + 0.0350330449527535d0*x9
      x348 = x147*x347
      x349 = 0.0891107483661069d0*x22*x28
      x350 = 0.0231073345558005d0*x149
      x351 = 0.00107855455836145d0*x79
      x352 = 0.000586927869208408d0*x70
      x353 = 0.000440195901906306d0*x165*x61
      x354 = x152*x62
      x355 = x329 - x330 - x331 + x333 - x334
      x356 = x166*(-x153*(-1.36477107997868d0*x4 + 0.427594146176807d0*x
     &9) - 4.83324359125373d0*x154*x354 - 1.61108119708458d0*x155 - 0.15
     &4912426780983d0*x3**(-0.333333333333333d0)*x63 - 75.398223686155d0
     &*x354*x64 + x355)
      x357 = 0.0110340749409504d0*x73
      x358 = 0.0025467176121001d0*x205
      x359 = 0.0648355077773455d0*x75
      x360 = 0.0165511124114256d0*x185*x61
      x361 = 0.00471283037863406d0*x136*x204
      x362 = 0.213300270968718d0*x187*x74
      x363 = x207*x339
      x364 = 0.263148661638301d0*x342*x74
      x365 = x207*x343
      x366 = x184*x335 + x336
      x367 = x206*x23*x366
      x368 = 0.89531658564596d0*x6
      x369 = 0.950793276165107d0*x191 + x368
      x370 = 0.275711329288622d0*x191
      x371 = mu*x192
      x372 = x371*(1.42618991424766d0*x191 + 1.53482843253593d0*x6)/x48*
     &*2
      x373 = x190 + x194
      x374 = x196 + x373
      x375 = x371*x49
      x376 = 1.43555326817173d0*x193 + x368 - 2.0d0*x374*x375
      x377 = 0.000821860773992031d0*x212
      x378 = 0.00880391803812612d0*x202
      x379 = 0.00316102962473761d0*x215
      x380 = 0.000279369479946145d0*x136*x209 - x203*x23*x347 - 6.393508
     &82283002d-5* x21*x28*x3**2.33333333333333d0 - 0.000188706975877918
     &d0*x210 - 0.0117385573841682d0*x211 + 0.00191767513931474d0*x213 -
     & 0.00626323892297956d0*x214 + 0.0155990559481421d0*x215*x342 + 0.0
     &168588246652672d0*x216 + x339*x379 - x343*x379
      x381 = x219*(-mu*x302*x92*(-x195*x280 + x373)/x281**2 - x17*(-x345
     & + x346 + x348 + x349 - x350 + x351 + x352 - x353 - x356) - x303*(
     &-x280* x370 + x281*x372 - x282*x369 + x376) + x304*x375 + x305*x37
     &7 - x305*x378 + x380 + x53*(-x177*(x323 + x324 - x325 + x326 + x32
     &7 + x328 - x329 + x330 + x331 - x333 + x334) - 0.0009195062450792d
     &0* x300 + x301*x322 - x338*(-x335*(x180 + x181 - x182 + x183) + x3
     &36 ) + x344) - x78*(x357 + x358 + x359 - x360 - x361 - x362 - x363
     & - x364 + x365 + x367))
      x382 = x108*x233
      x383 = x236 - x237 - x238 - x241 + x245 - x246
      x384 = x235*x383
      x385 = x249*x384
      x386 = x108*x385
      x387 = x257 - x258 + 26.8829425884854d0
      x388 = x260*x292
      x389 = x388*(x268 - x269 - x270 - x271 + x272 - x273)
      x390 = x253 - 40.6744021938219d0*x298 + 64.3272972886044d0*x307
      x391 = x110*x390*x87
      x392 = x87*x89
      x393 = x390*x392
      x394 = x233 + x249*x393 + x385
      x395 = x264*x292 + x266*x391 + x275*x389 + x394
      x396 = x395*x89
      x397 = x122*x396
      x398 = x262*x387 + x382 + x386 - x397
      x399 = x116*x120
      x400 = rho_c**(-4.33333333333333d0)
      x401 = x104*x400
      x402 = x108*x401
      x403 = x384*x402
      x404 = x383**2
      x405 = x98**(-2)
      x406 = 5.76063055467543d0*x405
      x407 = x402*x406
      x408 = 0.509120303595743d0*x239
      x409 = x1**(-3.78d0)
      x410 = x103*x400
      x411 = x100*x136
      x412 = -x130 + x137
      x413 = x383*(-x299 + x308 + x412)
      x414 = rho_c**(-6.66666666666667d0)
      x415 = x114*x414
      x416 = x259*x415
      x417 = x387*x416
      x418 = 23.0425222187017d0*x384
      x419 = 0.0663841159646404d0*x129
      x420 = 0.491388722315281d0*x312
      x421 = 2.78d0*x132
      x422 = x313*x421
      x423 = x132*x341
      x424 = 13.7187501175142d0*x423
      x425 = x316*x421
      x426 = 5.85934965221761d0*x141
      x427 = x318*x97
      x428 = 2.77866731225255d0*x427
      x429 = x9*x98/x5**2
      x430 = 0.0217142874477198d0*x429
      x431 = 0.327592481543521d0*x240
      x432 = x243*x313
      x433 = 27.4375002350284d0*x432
      x434 = x315/x1**2
      x435 = x434*x98/x14**2
      x436 = 270.797273074527d0*x435
      x437 = x242*x434
      x438 = 135.398636537263d0*x437
      x439 = x102*x11*x312
      x440 = 2.42490614820429d0*x439
      x441 = x243*x316
      x442 = 27.4375002350284d0*x441
      x443 = x218*x320
      x444 = 7.03121958266113d0*x443
      x445 = x130 - x137 - x142 - x221
      x446 = x239*x445
      x447 = 0.245694361157641d0*x446
      x448 = x11*x411*x445
      x449 = 27.4375002350284d0*x448
      x450 = x219*(0.202642367284676d0*mu*x198*x374*x45/x55**2 + x17*(x3
     &45 - x346 - x348 - x349 + x350 - x351 - x352 + x353 + x356) + x199
     &*(-x369*x56 + x370*x51 + x372*x55 + x376) - x200*x375 - x201*x377 
     &+ x201*x378 + x380 + x53*(-x176*x322 + x177*(-x323 - x324 + x325 -
     & x326 - x327 - x328 + x355) - x338*x366 + x344 + 0.000919506245079
     &2d0*x38 ) + x78*(-x357 - x358 - x359 + x360 + x361 + x362 + x363 +
     & x364 - x365 - x367))
      x451 = 5.56d0*x450
      x452 = -11.5212611093509d0*x419 + 11.5212611093509d0*x420 + 11.521
     &2611093509d0* x422 + 11.5212611093509d0*x424 - 11.5212611093509d0*
     &x425 + 11.5212611093509d0*x426 - 11.5212611093509d0*x428 - 11.5212
     &611093509d0*x430 + 11.5212611093509d0*x431 - 11.5212611093509d0*x4
     &33 - 11.5212611093509d0*x436 - 11.5212611093509d0*x438 + 11.521261
     &1093509d0*x440 + 11.5212611093509d0*x442 - 11.5212611093509d0*x444
     & + 11.5212611093509d0*x447 - 11.5212611093509d0*x449 - 11.52126110
     &93509d0*x451
      x453 = x108*x121*x232
      x454 = x396*x453
      x455 = x260*x387*x399
      x456 = x225*x396
      x457 = 11.5212611093509d0*x232
      x458 = 132.739457549841d0*x235*x389 + 265.478915099682d0*x265*x391
     & + 1238.90160379851d0*x388
      x459 = 23.0425222187017d0*x105
      x460 = x110*x395
      x461 = x116/x119**2
      x462 = 179.219617256569d0*x91
      x463 = x390**2
      x464 = 13.0354358035879d0*x40
      x465 = 0.169899675882914d0*x80
      x466 = 0.0412053987465066d0*x84
      x467 = 13.0354358035879d0*x57
      x468 = 64.3272972886044d0*x76
      x469 = 64.3272972886044d0*x71
      x470 = 64.3272972886044d0*x43
      x471 = 5.76063055467543d0*x91
      x472 = x392*x471
      x473 = x235*x91
      x474 = x383*x473
      x475 = x405*x471
      x476 = 56.8551446754745d0*x409*x411
      x477 = 2.0d0*x252
      x478 = 1.53608300823691d0*x129 - 9.86960440108936d0*x252*x341 - 11
     &.3703987139002d0*x312 - x313*x477 + x316*x477
      x479 = 53.7658851769707d0*x401
      x480 = x225*x255
      x481 = gamma_cc*x415
      x482 = x267*x481
      x483 = 11.5212611093509d0*x401
      x484 = x254**2
      x485 = x484*exp(x20*(-x465 - x466 + x467 - x468 - x469 + x470) + x
     &464)
      x486 = x222*x484
      x487 = x256*x401
      x488 = 132.739457549841d0*x224*x481
      x489 = x401*x406
      x490 = x235*x416
      x491 = x274*x490
      x492 = 66.3697287749204d0*x405*x416
      x493 = x235*x274
      x494 = x142 + x221 + x412
      x495 = x222*(-135.581340646073d0*x141 + 64.296459806258d0*x427 + 1
     &62.697608775288d0*x443 + 128.654594577209d0*x450 + x478)
      x496 = x122*x225
      x497 = x120*x309
      x498 = x106 + x259*x293
      x499 = 132.739457549841d0*x115
      x500 = 530.957830199363d0*x115*x259 + x459
      x501 = x109 + x292*x499 - x496*x498
      x502 = 265.478915099682d0*x115
      x503 = x225*x229
      E = -rho_c*(-x124 + x127)
      d1E(1) = x124 + x126 - x279*(-x131 + x138 + x142*x3 + x221*x3 - x2
     &78*(-x234 + x251 + x263 + x277)/(x231 + 1.0d0)) - x93
      d1E(3) = rho_c*x295*x297
      d2E(1) = x131 - x138 - x278*x309*x310 + x279*(0.031091d0*gamma_cc*
     &x497*(rho_c**( -6.66666666666667d0)*x112*x286*x292*(x110*x286*x463
     &*exp(-x20*( x465 + x466 - x467 + x468 + x469 - x470) + x464) - x38
     &4*x476 + 11.5212611093509d0*x393*x474 + 53.7658851769707d0*x393*x9
     &1 + x404 *x475 + x408*x474 + x413*x475 - x452*x473 + x462 - x463*x
     &472 - x472*(135.581340646073d0*x298 - 64.296459806258d0*x319 + 162
     &.697608775288d0*x321 + 128.654594577209d0*x381 + x478) + 53.765885
     &1769707d0*x474) - 56.8551446754745d0*x108*x384*x409*x410 *x411 + x
     &108*x459*x460*x461*(53.7658851769707d0*x232 + x384*x457 + x393*x45
     &7 + x458) - x235*x402*x452 + 179.219617256569d0*x402 + x403*x408 +
     & 53.7658851769707d0*x403 + x404*x407 + x407*x413 + x417*x418 + 107
     &.531770353941d0*x417 - x418*x454 - 107.531770353941d0*x454 - 46.08
     &50444374035d0*x455*x456 - x496*( 655.042967216062d0*x1**(-6.56d0)*
     &x112*x259*x411*x414*x493 + x224* x483*x485 + x235*x483*(x419 - x42
     &0 - x422 - x424 + x425 - x426 + x428 + x430 - x431 + x433 + x436 +
     & x438 - x440 - x442 + x444 - x447 + x449 + x451) - 5.8657079537985
     &5d0*x239*x491 + x247**2*x489 - x247*x489*x494 - x248*x401*x408 + x
     &248*x410*x476 - x248*x479 - x248*x480*x483 + x274**2*x492 - x274*x
     &492*x494 + x410*x462 + 7020.44242152492d0*x416 + x479*x480 - 265.4
     &78915099682d0*x482* x493 + 2477.80320759703d0*x482 - x486*x487 - x
     &486*x488 - x487* x495 - x488*x495 + 132.739457549841d0*x490*(0.132
     &768231929281d0* x129 - 11.7186993044352d0*x141 - 0.655184963087042
     &d0*x240 - 0.982777444630562d0*x312 - 5.56d0*x314 + 5.56d0*x317 - 2
     &7.4375002350284d0*x423 + 5.5573346245051d0*x427 + 0.04342857489543
     &97d0*x429 + 54.8750004700568d0*x432 + 541.594546149053d0*x435 + 27
     &0.797273074527d0*x437 - 4.84981229640858d0*x439 - 54.8750004700568
     &d0*x441 + 14.0624391653223d0*x443 - 0.491388722315281d0*x446 + 54.
     &8750004700568d0*x448 + 11.12d0*x450) - 1238.90160379851d0*x491 + 3
     &98.218372649523d0*x481*x485/x223**3)) + x12*x132*x311*x315 - 0.023
     &8791784045469d0*x129*x311 - 2.10767973101353d0*x298*x311 + 0.17675
     &8533206936d0*x311*x312 + x311*x314 - x311*x317 + 0.999520615918184
     &d0*x311*x319 - 2.52921567721623d0*x311*x321 - 2.0d0*x311*x381 - 0.
     &062182d0*x228*x310*x398*x399/x123**2) + x299* x3 - x3*x308
      d2E(4) = x297*(rho_c*x120*(-x109*x456*x461*x500 + x121*x225*x382*x
     &498 + x121*x296 *x398*x501 + x122*x89*(x394 + x458) + 11.521261109
     &3509d0*x224* x455*x498 + x256*x384*x453*x498 - 23.0425222187017d0*
     &x261*x387 - x382 - 66.3697287749204d0*x384*x388 - x386 - 309.72540
     &0949629d0* x388 + x397 + x399*x460*x499) + x295)
      d2E(6) = -0.031091d0*rho_c*x228*x309*(x225*(3058.65189989063d0*rho
     &_c**(-7.0d0)*x1 **(-8.34d0)*x102**8.34d0*x226*x503 + x226*x228*x45
     &9*x498 + x498* x502*x503 - x502 - x230*x259*x498*x500/x227**2) + x
     &497*x501*(x230 + x259*x499 - x276*x498))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu, E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu*rho_a**(-0.333333333333333d0)
      x2 = mu**2*rho_a**(-0.666666666666667d0)
      x3 = gamma_aa*rho_a**(-2.66666666666667d0)
      x4 = gamma_aa**2
      x5 = rho_a**(-5.33333333333333d0)
      x6 = 1/(rho_a**2*tau_a**2)
      x7 = x4*x6
      x8 = 0.0109702052306804d0*x3
      x9 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x10 = x3*x9
      x11 = (x8*x9*(0.0274255130767009d0*x10 - 1.0d0) + 1.0d0)**(-0.5d0)
     &
      Ea = -1.8610514726982d0*x0*((0.804d0/(0.00449266744437808d0*x3 + 1
     &.0d0) - 0.804d0/(1.0d0 + (1.86977850568312d-6*gamma_aa**3*rho_a**(
     &-8.0d0 ) + 0.0204668008035082d0*x3*(0.02485875d0*x7/(0.015625d0*x7
     & + 1.0d0)**2 + 0.123456790123457d0) + 6.38454405015215d-6*x4*x5 + 
     &0.0021416494625995d0*x7 - 1.24378109452736d0*sqrt( 0.0001353885781
     &53653d0*x5 + 0.0028125d0*x6)*(x11*( 0.00222451383844352d0*x10 - 0.
     &0811111111111111d0) + 0.00197734563417202d0*x3)*abs(gamma_aa) + 0.
     &0896750813832074d0*( x11*(0.0123414808845154d0*x10 - 0.45d0) + x8)
     &**2)/( 0.0204006002177468d0*x3 + 1.0d0)**2))*exp(-1.92417365779563
     &32d0* x1) + 1.804d0 - 0.804d0/(0.0519881082712438d0*x3*( 1.1847821
     &6491348d0*x2 - 7.0d0)*exp(-0.31265084907439067d0*x2)/( 0.233951463
     &049512d0*mu**4/x0 + 49.1113297378198d0*x1 - 10.6630394842213d0*x2 
     &- 81.0d0) + 1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, mu,
     & Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**0.333333333333333d0
      x2 = 1d0/x1
      x3 = mu*x2
      x4 = mu**2
      x5 = rho_a**0.666666666666667d0
      x6 = x4/x5
      x7 = gamma_aa*rho_a**(-2.66666666666667d0)
      x8 = 0.00844343197019482d0*1d0/rho_a
      x9 = gamma_aa**2
      x10 = rho_a**(-5.33333333333333d0)
      x11 = 1/(rho_a**2*tau_a**2)
      x12 = x11*x9
      x13 = 0.0109702052306804d0*x7
      x14 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x15 = x14*x7
      x16 = (x13*x14*(0.0274255130767009d0*x15 - 1.0d0) + 1.0d0)**(-0.5d
     &0)
      Ea = -x0*(-4.9628039271952d0*(1.804d0 - 0.804d0/(0.051988108271243
     &8d0*x7*( 1.18478216491348d0*x6 - 7.0d0)*exp(-0.31265084907439067d0
     &*x6)/( 0.233951463049512d0*mu**4/x0 + 49.1113297378198d0*x3 - 10.6
     &630394842213d0*x6 - 81.0d0) + 1.0d0))*(0.128278243853042d0*x3 *(mu
     &**3*x8 + mu*(0.256556487706084d0*x2 - x4*x8)*exp( -15.192666241151
     &986d0*x5/x4) - 0.384834731559127d0*x3 + 1.77245385090552d0*erf(3.8
     &977770897207535d0*x1/mu)) - 0.375d0) + (1.49628538404935d0/(0.0044
     &9266744437808d0*x7 + 1.0d0) - 1.49628538404935d0/(1.0d0 + (1.86977
     &850568312d-6*gamma_aa**3* rho_a**(-8.0d0) + 6.38454405015215d-6*x1
     &0*x9 + 0.0021416494625995d0*x12 + 0.0204668008035082d0*x7*(0.02485
     &875d0* x12/(0.015625d0*x12 + 1.0d0)**2 + 0.123456790123457d0) - 1.
     &24378109452736d0*sqrt(0.000135388578153653d0*x10 + 0.0028125d0* x1
     &1)*(x16*(0.00222451383844352d0*x15 - 0.0811111111111111d0) + 0.001
     &97734563417202d0*x7)*abs(gamma_aa) + 0.0896750813832074d0*( x13 + 
     &x16*(0.0123414808845154d0*x15 - 0.45d0))**2)/( 0.0204006002177468d
     &0*x7 + 1.0d0)**2))*exp(-1.9241736577956332d0* x3))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, mu,
     & Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu**2
      x2 = 1d0/x1
      x3 = rho_a**0.666666666666667d0
      x4 = gamma_aa*rho_a**(-2.66666666666667d0)
      x5 = rho_a**0.333333333333333d0
      x6 = 1d0/x5
      x7 = mu*x6
      x8 = 0.00844343197019482d0*1d0/rho_a
      x9 = gamma_aa**2
      x10 = rho_a**(-5.33333333333333d0)
      x11 = 1/(rho_a**2*tau_a**2)
      x12 = x11*x9
      x13 = 0.0109702052306804d0*x4
      x14 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x15 = x14*x4
      x16 = (x13*x14*(0.0274255130767009d0*x15 - 1.0d0) + 1.0d0)**(-0.5d
     &0)
      Ea = -x0*(-4.9628039271952d0*(1.804d0 - 0.804d0/(-0.05198810827124
     &38d0*x2*x4* (0.213719543995407d0*x0*x2 - 0.844037013397333d0*x3 + 
     &14.4026298885989d0*rho_a**2.0d0/mu**4)*exp(-0.31265084907439067d0 
     &*x1/x3) + 1.0d0))*(0.128278243853042d0*x7*(mu**3*x8 + mu*(-x1*x8 +
     & 0.256556487706084d0*x6)*exp(-15.192666241151986d0*x2*x3) - 0.3848
     &34731559127d0*x7 + 1.77245385090552d0*erf( 3.8977770897207535d0*x5
     &/mu)) - 0.375d0) + (1.49628538404935d0/( 0.00449266744437808d0*x4 
     &+ 1.0d0) - 1.49628538404935d0/(1.0d0 + ( 1.86977850568312d-6*gamma
     &_aa**3*rho_a**(-8.0d0) + 6.38454405015215d-6*x10*x9 + 0.0021416494
     &625995d0*x12 + 0.0204668008035082d0*x4*(0.02485875d0*x12/(0.015625
     &d0*x12 + 1.0d0 )**2 + 0.123456790123457d0) - 1.24378109452736d0*sq
     &rt( 0.000135388578153653d0*x10 + 0.0028125d0*x11)*(x16*( 0.0022245
     &1383844352d0*x15 - 0.0811111111111111d0) + 0.00197734563417202d0*x
     &4)*abs(gamma_aa) + 0.0896750813832074d0*( x13 + x16*(0.01234148088
     &45154d0*x15 - 0.45d0))**2)/( 0.0204006002177468d0*x4 + 1.0d0)**2))
     &*exp(-1.9241736577956332d0* x7))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, mu,
     & Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**0.333333333333333d0
      x2 = 1d0/x1
      x3 = mu*x2
      x4 = 0.00844343197019482d0*1d0/rho_a
      x5 = erf(3.8977770897207535d0*x1/mu)
      x6 = mu**2
      x7 = 1d0/x6
      x8 = rho_a**0.666666666666667d0
      x9 = 15.192666241152d0*x7*x8
      x10 = x6/x8
      x11 = exp(x9)
      x12 = mu**4*(x11 - 1.0d0)/x0
      x13 = 3.0d0*x11
      x14 = gamma_aa*rho_a**(-2.66666666666667d0)
      x15 = gamma_aa**2
      x16 = rho_a**(-5.33333333333333d0)
      x17 = 1/(rho_a**2*tau_a**2)
      x18 = x15*x17
      x19 = 0.0109702052306804d0*x14
      x20 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x21 = x14*x20
      x22 = (x19*x20*(0.0274255130767009d0*x21 - 1.0d0) + 1.0d0)**(-0.5d
     &0)
      Ea = -x0*(-4.9628039271952d0*(1.804d0 - 0.804d0/(-3.15935190989552
     &d0*gamma_aa *rho_a**(-2.0d0)*x7*(0.00060945584614891d0*x10*(7.0d0*
     &x11 + 11.0d0) - 0.000722072416819482d0*x12 + 0.0185185185185185d0)
     &*exp( -0.31265084907439067d0*x10)/(-0.131642462768164d0*x10*(x13 -
     & 2.0d0) + 1.81893813843777d0*x11*x3*x5 + 0.00866486900183379d0*x12
     & - x13) + 1.0d0))*(0.128278243853042d0*x3*(mu**3*x4 + mu*( 0.25655
     &6487706084d0*x2 - x4*x6)*exp(-x9) - 0.384834731559127d0*x3 + 1.772
     &45385090552d0*x5) - 0.375d0) + (1.49628538404935d0/( 0.00449266744
     &437808d0*x14 + 1.0d0) - 1.49628538404935d0/(1.0d0 + (1.86977850568
     &312d-6*gamma_aa**3*rho_a**(-8.0d0) + 0.0204668008035082d0*x14*(0.0
     &2485875d0*x18/(0.015625d0*x18 + 1.0d0)**2 + 0.123456790123457d0) +
     & 6.38454405015215d-6*x15*x16 + 0.0021416494625995d0*x18 - 1.243781
     &09452736d0*( 0.00197734563417202d0*x14 + x22*(0.00222451383844352d
     &0*x21 - 0.0811111111111111d0))*sqrt(0.000135388578153653d0*x16 + 0
     &.0028125d0*x17)*abs(gamma_aa) + 0.0896750813832074d0*(x19 + x22* (
     &0.0123414808845154d0*x21 - 0.45d0))**2)/(0.0204006002177468d0* x14
     & + 1.0d0)**2))*exp(-1.9241736577956332d0*x3))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu, E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12
      Ea = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**0.666666666666667d0
      x3 = gamma_aa*rho_a**(-2.66666666666667d0)
      x4 = rho_a**0.333333333333333d0
      x5 = gamma_aa**2
      x6 = rho_a**(-5.33333333333333d0)
      x7 = 1/(rho_a**2*tau_a**2)
      x8 = x5*x7
      x9 = 0.0109702052306804d0*x3
      x10 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x11 = x10*x3
      x12 = (x10*x9*(0.0274255130767009d0*x11 - 1.0d0) + 1.0d0)**(-0.5d0
     &)
      Ea = -rho_a*(3.14159265358979d0*rho_a**1.0d0*x1*(1.804d0 - 0.804d0
     &/( -0.0519881082712438d0*x1*x3*(0.213719543995407d0*rho_a** 1.3333
     &3333333333d0*x1 - 0.844037013397333d0*x2 + 14.4026298885989d0*rho_
     &a**2.0d0/mu**4)*exp(-0.31265084907439067d0 *x0/x2) + 1.0d0)) + 1.8
     &610514726982d0*x4*(0.804d0/( 0.00449266744437808d0*x3 + 1.0d0) - 0
     &.804d0/(1.0d0 + ( 1.86977850568312d-6*gamma_aa**3*rho_a**(-8.0d0) 
     &+ 0.0204668008035082d0*x3*(0.02485875d0*x8/(0.015625d0*x8 + 1.0d0)
     & **2 + 0.123456790123457d0) + 6.38454405015215d-6*x5*x6 + 0.002141
     &6494625995d0*x8 - 1.24378109452736d0*sqrt( 0.000135388578153653d0*
     &x6 + 0.0028125d0*x7)*(x12*( 0.00222451383844352d0*x11 - 0.08111111
     &11111111d0) + 0.00197734563417202d0*x3)*abs(gamma_aa) + 0.08967508
     &13832074d0*( x12*(0.0123414808845154d0*x11 - 0.45d0) + x9)**2)/( 0
     &.0204006002177468d0*x3 + 1.0d0)**2))*exp(-1.9241736577956332d0* mu
     &/x4))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**0.333333333333333d0
      x1 = mu/x0
      x2 = exp(-1.9241736577956332d0*x1)
      x3 = rho_a**(-2.66666666666667d0)
      x4 = gamma_aa*x3
      x5 = 0.00449266744437808d0*x4 + 1.0d0
      x6 = 0.804d0/x5
      x7 = 0.0204006002177468d0*x4 + 1.0d0
      x8 = x7**(-2)
      x9 = 0.0109702052306804d0*x4
      x10 = tau_a/gamma_aa
      x11 = 8.0d0*rho_a*x10
      x12 = x11 - 1.0d0
      x13 = x12*x4
      x14 = 0.0123414808845154d0*x13 - 0.45d0
      x15 = 0.0274255130767009d0*x12
      x16 = x15*x4 - 1.0d0
      x17 = x12*x16
      x18 = x17*x9 + 1.0d0
      x19 = x18**(-0.5d0)
      x20 = x14*x19 + x9
      x21 = x20**2
      x22 = 0.00197734563417202d0*x3
      x23 = gamma_aa*x22
      x24 = x19*(0.00222451383844352d0*x13 - 0.0811111111111111d0) + x23
     &
      x25 = rho_a**(-5.33333333333333d0)
      x26 = rho_a**(-2)
      x27 = tau_a**(-2)
      x28 = 0.0028125d0*x27
      x29 = sqrt(0.000135388578153653d0*x25 + x26*x28)
      x30 = abs(gamma_aa)
      x31 = x29*x30
      x32 = x24*x31
      x33 = gamma_aa**3
      x34 = rho_a**(-8.0d0)
      x35 = x33*x34
      x36 = gamma_aa**2
      x37 = x25*x36
      x38 = x27*x36
      x39 = x26*x38
      x40 = 0.015625d0*x39 + 1.0d0
      x41 = x40**(-2)
      x42 = 0.02485875d0*x39*x41 + 0.123456790123457d0
      x43 = 0.0204668008035082d0*x42
      x44 = 1.86977850568312d-6*x35 + 6.38454405015215d-6*x37 + 0.002141
     &6494625995d0 *x39 + x4*x43
      x45 = x8*(0.0896750813832074d0*x21 - 1.24378109452736d0*x32 + x44)
     & + 1.0d0
      x46 = -0.804d0/x45
      x47 = x2*(x46 + x6)
      x48 = mu**2
      x49 = rho_a**(-0.666666666666667d0)*x48
      x50 = rho_a**1.33333333333333d0
      x51 = 1d0/x50
      x52 = 1d0/(0.233951463049512d0*mu**4*x51 + 49.1113297378198d0*x1 -
     & 10.6630394842213d0*x49 - 81.0d0)
      x53 = x52*exp(-0.31265084907439067d0*x49)
      x54 = 1.18478216491348d0*x49 - 7.0d0
      x55 = 0.0519881082712438d0*x54
      x56 = x4*x53*x55 + 1.0d0
      x57 = 1.804d0 - 0.804d0/x56
      x58 = 1.8610514726982d0*x50
      x59 = 1.8610514726982d0*x0
      x60 = 0.6203504908994d0*rho_a**(-0.666666666666667d0)
      x61 = rho_a**(-4.33333333333333d0)
      x62 = x48*x61
      x63 = rho_a**(-3.66666666666667d0)
      x64 = rho_a**(-1.66666666666667d0)
      x65 = rho_a**(-2.33333333333333d0)
      x66 = x53/x56**2
      x67 = x5**(-2)
      x68 = rho_a**(-6.33333333333333d0)
      x69 = rho_a**(-3)
      x70 = gamma_aa*x63
      x71 = 0.001553671875d0*x39/x40 - 0.0497175d0
      x72 = 0.0204668008035082d0*x41*x71
      x73 = 1d0/x29
      x74 = 1d0/x30
      x75 = 1.24378109452736d0*x74
      x76 = 1d0/x7
      x77 = 1.50330191856923d-6*x35 + 5.13317341632233d-6*x37 + 0.001721
     &88616792999d0*x39 + 0.0164553078460206d0*x4*x42
      x78 = tau_a*x3
      x79 = x12*x70
      x80 = 0.00548510261534019d0*x13
      x81 = x14*x18**(-1.5d0)*(-0.0438808209227215d0*x16*x78 + 0.0146269
     &403075738d0* x16*x79 - x80*(0.219404104613608d0*x78 - 0.0731347015
     &378692d0*x79 ))
      x82 = 1.24378109452736d0*x31
      x83 = x0*x8
      x84 = -x11 + 1.0d0
      x85 = x4*x84
      x86 = 0.0123414808845154d0*x85 + 0.45d0
      x87 = 0.0274255130767009d0*x85 + 1.0d0
      x88 = x84*x87*x9 + 1.0d0
      x89 = x88**(-0.5d0)
      x90 = -x86*x89 + x9
      x91 = x90**2
      x92 = x23 - x89*(0.00222451383844352d0*x85 + 0.0811111111111111d0)
     &
      x93 = x31*x92
      x94 = (x8*(x44 + 0.0896750813832074d0*x91 - 1.24378109452736d0*x93
     &) + 1.0d0)** (-2)
      x95 = rho_a**(-4.66666666666667d0)
      x96 = 0.0177961107075482d0*x64
      x97 = x3*x84
      x98 = 0.0438808209227215d0*x64
      x99 = x10*x64
      x100 = x88**(-1.5d0)*(-x10*x16*x98 + 0.00548510261534019d0*x17*x3 
     &+ x80*(x15*x3 - 0.219404104613608d0*x99))
      x101 = 0.180246913580247d0*x86
      x102 = 2.0d0*x86
      x103 = tau_a**(-3)
      x104 = x103*x26*x36
      x105 = (0.00120345402803247d0*gamma_aa*x61*x84 + x87*x98)/x88
      Ea = -x58*(x47 + x57)
      d1Ea(1) = -rho_a*(1.49628538404935d0*gamma_aa*x0*x66*(mu*x3*x52*x5
     &5*( 0.311935284066016d0*mu**3*x65 - 7.10869298948088d0*mu*x64 + 16
     &.3704432459399d0*x51) + 0.0108360841285172d0*x54*x62 - 0.138634955
     &389984d0*x54*x63 - 0.0410630556449071d0*x62) + 1.19366207318922d0*
     &mu*1d0/rho_a*x47 + x2*(0.0179261670197794d0* gamma_aa*rho_a**(-3.3
     &3333333333333d0)*x67 + x60*(x46 + 1.804d0) - x60*(-x6 + 1.804d0) +
     & 1.49628538404935d0*x83*( -1.49582280454649d-5*rho_a**(-9.0d0)*x33
     & + rho_a**( -5.66666666666667d0)*x27*x33*x72 + 0.0896750813832074d
     &0*x20*( 2.0d0*x19*(0.0987318470761234d0*x78 - 0.0329106156920411d0
     &*x79) - 0.0585077612302953d0*x70 + 2.0d0*x81) + x24*x36*x73*x75*(x
     &28*x69 + 0.000361036208409741d0*x68) - 3.40509016008115d-5*x36*x68
     & - 0.00428329892519899d0*x38*x69 - 0.0545781354760218d0*x42*x70 + 
     &0.135327364628503d0*x70*x76*(0.0720987654320988d0*x21 - x32 + x77 
     &) - x82*(x19*(0.0177961107075482d0*x78 - 0.00593203690251605d0* x7
     &9) - 0.00527292169112538d0*x70 + 0.180246913580247d0*x81))/x45 **2
     &) + x57*x60) - x47*x59 - x57*x59
      d1Ea(2) = -rho_a*(-2.0d0*x2*(0.00336115631620864d0*x65*x67 + 0.748
     &142692024676d0* x83*x94*(-1.27690881003043d-5*gamma_aa*x25 - 0.004
     &28329892519899d0*gamma_aa*x26*x27 + gamma_aa*x29*x75*x92 - x3*x43 
     &+ 0.0507477617356885d0*x3*x76*(x77 + 0.0720987654320988d0* x91 - x
     &93) - 5.60933551704935d-6*x34*x36 + x38*x72*x95 + x82*( x100*x101 
     &+ x22 - x89*(x10*x96 + 0.00222451383844352d0*x97)) - 0.08967508138
     &32074d0*x90*(x100*x102 + 0.0219404104613608d0*x3 - 2.0d0*x89*(0.01
     &23414808845154d0*x97 + 0.0987318470761234d0*x99 )))) + 0.077789046
     &5506374d0*x54*x65*x66)
      d1Ea(3) = -x2*x58*x8*x94*(0.0164553078460206d0*x103*x33*x41*x71*x9
     &5 + 0.0028125d0* x104*x73*x74*x92 - 0.00344377233585999d0*x104 - x
     &31*x89*(-x101* x105 + x96) + 0.0720987654320988d0*x89*x90*(-x102*x
     &105 + 0.197463694152247d0*x64))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = gamma_aa*x1
      x3 = rho_a**0.333333333333333d0
      x4 = 1d0/x3
      x5 = mu*x4
      x6 = mu**2
      x7 = rho_a**0.666666666666667d0
      x8 = x6/x7
      x9 = 1d0/x0
      x10 = 1d0/(0.233951463049512d0*mu**4*x9 + 49.1113297378198d0*x5 - 
     &10.6630394842213d0*x8 - 81.0d0)
      x11 = x10*exp(-0.31265084907439067d0*x8)
      x12 = 1.18478216491348d0*x8 - 7.0d0
      x13 = 0.0519881082712438d0*x12
      x14 = x11*x13*x2 + 1.0d0
      x15 = 1.804d0 - 0.804d0/x14
      x16 = 0.256556487706084d0*x4
      x17 = 1d0/rho_a
      x18 = 0.00844343197019482d0*x17
      x19 = x18*x6
      x20 = exp(-15.192666241151986d0*x7/x6)
      x21 = mu*x20
      x22 = 1d0/mu
      x23 = mu**3
      x24 = x18*x23 - 0.384834731559127d0*x5 + 1.77245385090552d0*erf( 3
     &.8977770897207535d0*x22*x3)
      x25 = 0.128278243853042d0*x5
      x26 = x25*(x21*(x16 - x19) + x24) - 0.375d0
      x27 = x15*x26
      x28 = exp(-1.9241736577956332d0*x5)
      x29 = 0.00449266744437808d0*x2 + 1.0d0
      x30 = 1d0/x29
      x31 = 0.0204006002177468d0*x2 + 1.0d0
      x32 = x31**(-2)
      x33 = 0.0109702052306804d0*x2
      x34 = tau_a/gamma_aa
      x35 = 8.0d0*rho_a*x34
      x36 = x35 - 1.0d0
      x37 = x2*x36
      x38 = 0.0123414808845154d0*x37 - 0.45d0
      x39 = 0.0274255130767009d0*x36
      x40 = x2*x39 - 1.0d0
      x41 = x36*x40
      x42 = x33*x41 + 1.0d0
      x43 = x42**(-0.5d0)
      x44 = x33 + x38*x43
      x45 = x44**2
      x46 = 0.00197734563417202d0*x1
      x47 = gamma_aa*x46
      x48 = x43*(0.00222451383844352d0*x37 - 0.0811111111111111d0) + x47
     &
      x49 = rho_a**(-5.33333333333333d0)
      x50 = rho_a**(-2)
      x51 = tau_a**(-2)
      x52 = 0.0028125d0*x51
      x53 = sqrt(0.000135388578153653d0*x49 + x50*x52)
      x54 = abs(gamma_aa)
      x55 = x53*x54
      x56 = x48*x55
      x57 = gamma_aa**3
      x58 = rho_a**(-8.0d0)
      x59 = x57*x58
      x60 = gamma_aa**2
      x61 = x49*x60
      x62 = x51*x60
      x63 = x50*x62
      x64 = 0.015625d0*x63 + 1.0d0
      x65 = x64**(-2)
      x66 = 0.02485875d0*x63*x65 + 0.123456790123457d0
      x67 = 0.0204668008035082d0*x66
      x68 = x2*x67 + 1.86977850568312d-6*x59 + 6.38454405015215d-6*x61 +
     & 0.0021416494625995d0*x63
      x69 = x32*(0.0896750813832074d0*x45 - 1.24378109452736d0*x56 + x68
     &) + 1.0d0
      x70 = 1d0/x69
      x71 = -x16 + x19
      x72 = -x21*x71 + x24
      x73 = 4.9628039271952d0*x15*x3
      x74 = 0.804d0*x30
      x75 = -x35 + 1.0d0
      x76 = x2*x75
      x77 = 0.0123414808845154d0*x76 + 0.45d0
      x78 = 0.0274255130767009d0*x76 + 1.0d0
      x79 = x33*x75*x78 + 1.0d0
      x80 = x79**(-0.5d0)
      x81 = x33 - x77*x80
      x82 = x81**2
      x83 = x47 - x80*(0.00222451383844352d0*x76 + 0.0811111111111111d0)
     &
      x84 = x55*x83
      x85 = x32*(x68 + 0.0896750813832074d0*x82 - 1.24378109452736d0*x84
     &) + 1.0d0
      x86 = 1.8610514726982d0*x28
      x87 = rho_a**(-0.666666666666667d0)
      x88 = 0.00844343197019482d0*rho_a**(-2.0d0)
      x89 = x20*x22
      x90 = -0.804d0*x70
      x91 = rho_a**(-4.33333333333333d0)
      x92 = x6*x91
      x93 = rho_a**(-3.66666666666667d0)
      x94 = rho_a**(-1.66666666666667d0)
      x95 = rho_a**(-2.33333333333333d0)
      x96 = x11*x26/x14**2
      x97 = x29**(-2)
      x98 = 0.6203504908994d0*x87
      x99 = rho_a**(-6.33333333333333d0)
      x100 = rho_a**(-3)
      x101 = gamma_aa*x93
      x102 = 0.001553671875d0*x63/x64 - 0.0497175d0
      x103 = 0.0204668008035082d0*x102*x65
      x104 = 1d0/x53
      x105 = 1d0/x54
      x106 = 1.24378109452736d0*x105
      x107 = 1d0/x31
      x108 = 0.0164553078460206d0*x2*x66 + 1.50330191856923d-6*x59 + 5.1
     &3317341632233d-6*x61 + 0.00172188616792999d0*x63
      x109 = tau_a*x1
      x110 = x101*x36
      x111 = 0.00548510261534019d0*x37
      x112 = x38*x42**(-1.5d0)*(-0.0438808209227215d0*x109*x40 + 0.01462
     &69403075738d0 *x110*x40 - x111*(0.219404104613608d0*x109 - 0.07313
     &47015378692d0 *x110))
      x113 = 1.24378109452736d0*x55
      x114 = x3*x32
      x115 = x85**(-2)
      x116 = rho_a**(-4.66666666666667d0)
      x117 = 0.0177961107075482d0*x94
      x118 = x1*x75
      x119 = 0.0438808209227215d0*x94
      x120 = x34*x94
      x121 = x79**(-1.5d0)*(0.00548510261534019d0*x1*x41 + x111*(x1*x39 
     &- 0.219404104613608d0*x120) - x119*x34*x40)
      x122 = 0.180246913580247d0*x77
      x123 = 2.0d0*x77
      x124 = tau_a**(-3)
      x125 = x124*x50*x60
      x126 = (0.00120345402803247d0*gamma_aa*x75*x91 + x119*x78)/x79
      Ea = -x0*(-4.9628039271952d0*x27 + x28*(1.49628538404935d0*x30 - 1
     &.49628538404935d0*x70))
      d1Ea(1) = -rho_a*(-3.99009435746494d0*gamma_aa*x3*x96*(mu*x1*x10*x
     &13*( -7.10869298948088d0*mu*x94 + 0.311935284066016d0*x23*x95 + 16
     &.3704432459399d0*x9) + 0.0108360841285172d0*x12*x92 - 0.1386349553
     &89984d0*x12*x93 - 0.0410630556449071d0*x92) + 1.19366207318922d0*m
     &u*x17*x28*(x74 + x90) + mu*x73*( -0.128278243853042d0*x4*(0.128278
     &243853042d0*mu*x9 + 10.128444160768d0*rho_a**(-0.333333333333333d0
     &)*x71*x89 + x21*(x6 *x88 - 0.0855188292353615d0*x9) - x23*x88 + 2.
     &59851805981384d0* x87*x89) + 0.0427594146176807d0*x72*x9) - 1.6542
     &6797573173d0*x27* x87 + x28*(0.0179261670197794d0*gamma_aa*rho_a**
     &( -3.33333333333333d0)*x97 + 1.49628538404935d0*x114*( -1.49582280
     &454649d-5*rho_a**(-9.0d0)*x57 + rho_a**( -5.66666666666667d0)*x103
     &*x51*x57 - 0.00428329892519899d0*x100* x62 + 0.135327364628503d0*x
     &101*x107*(x108 + 0.0720987654320988d0* x45 - x56) - 0.054578135476
     &0218d0*x101*x66 + x104*x106*x48*x60*( x100*x52 + 0.000361036208409
     &741d0*x99) - x113*( -0.00527292169112538d0*x101 + 0.18024691358024
     &7d0*x112 + x43*( 0.0177961107075482d0*x109 - 0.00593203690251605d0
     &*x110)) + 0.0896750813832074d0*x44*(-0.0585077612302953d0*x101 + 2
     &.0d0*x112 + 2.0d0*x43*(0.0987318470761234d0*x109 - 0.0329106156920
     &411d0* x110)) - 3.40509016008115d-5*x60*x99)/x69**2 - x98*(-x74 + 
     &1.804d0) + x98*(x90 + 1.804d0))) - x3*x86*(x74 - 0.804d0/x85) + x7
     &3*(x25*x72 - 0.375d0)
      d1Ea(2) = rho_a*(0.207437457468366d0*x12*x95*x96 + 2.0d0*x28*(0.74
     &8142692024676d0* x114*x115*(gamma_aa*x106*x53*x83 - 1.276908810030
     &43d-5*gamma_aa* x49 - 0.00428329892519899d0*gamma_aa*x50*x51 + 0.0
     &507477617356885d0*x1*x107*(x108 + 0.0720987654320988d0*x82 - x84) 
     &- x1*x67 + x103*x116*x62 + x113*(x121*x122 + x46 - x80*(x117 *x34 
     &+ 0.00222451383844352d0*x118)) - 5.60933551704935d-6*x58*x60 - 0.0
     &896750813832074d0*x81*(0.0219404104613608d0*x1 + x121*x123 - 2.0d0
     &*x80*(0.0123414808845154d0*x118 + 0.0987318470761234d0*x120 ))) + 
     &0.00336115631620864d0*x95*x97))
      d1Ea(3) = -x0*x115*x32*x86*(0.0164553078460206d0*x102*x116*x124*x5
     &7*x65 + 0.0028125d0*x104*x105*x125*x83 - 0.00344377233585999d0*x12
     &5 - x55 *x80*(x117 - x122*x126) + 0.0720987654320988d0*x80*x81*(-x
     &123* x126 + 0.197463694152247d0*x94))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = gamma_aa*x1
      x3 = mu**2
      x4 = 1d0/x3
      x5 = 0.0519881082712438d0*x4
      x6 = rho_a**0.666666666666667d0
      x7 = 0.844037013397333d0*x6
      x8 = mu**(-4)
      x9 = rho_a**2.0d0
      x10 = 14.4026298885989d0*x8*x9
      x11 = 0.213719543995407d0*x0*x4
      x12 = exp(-0.31265084907439067d0*x3/x6)
      x13 = x12*(x10 + x11 - x7)
      x14 = x13*x2*x5
      x15 = -x14 + 1.0d0
      x16 = 1.804d0 - 0.804d0/x15
      x17 = rho_a**0.333333333333333d0
      x18 = 1d0/x17
      x19 = 0.256556487706084d0*x18
      x20 = rho_a**1.0d0
      x21 = 1d0/x20
      x22 = 0.00844343197019482d0*x21
      x23 = x22*x3
      x24 = exp(-15.192666241151986d0*x4*x6)
      x25 = mu*x24
      x26 = 1d0/mu
      x27 = mu**3
      x28 = mu*x18
      x29 = x22*x27 - 0.384834731559127d0*x28 + 1.77245385090552d0*erf( 
     &3.8977770897207535d0*x17*x26)
      x30 = 0.128278243853042d0*x28
      x31 = x30*(x25*(x19 - x23) + x29) - 0.375d0
      x32 = x16*x31
      x33 = exp(-1.9241736577956332d0*x28)
      x34 = 0.00449266744437808d0*x2 + 1.0d0
      x35 = 1d0/x34
      x36 = 0.0204006002177468d0*x2 + 1.0d0
      x37 = x36**(-2)
      x38 = 0.0109702052306804d0*x2
      x39 = tau_a/gamma_aa
      x40 = 8.0d0*rho_a*x39
      x41 = x40 - 1.0d0
      x42 = x2*x41
      x43 = 0.0123414808845154d0*x42 - 0.45d0
      x44 = 0.0274255130767009d0*x41
      x45 = x2*x44 - 1.0d0
      x46 = x41*x45
      x47 = x38*x46 + 1.0d0
      x48 = x47**(-0.5d0)
      x49 = x38 + x43*x48
      x50 = x49**2
      x51 = 0.00197734563417202d0*x1
      x52 = gamma_aa*x51
      x53 = x48*(0.00222451383844352d0*x42 - 0.0811111111111111d0) + x52
     &
      x54 = rho_a**(-5.33333333333333d0)
      x55 = rho_a**(-2)
      x56 = tau_a**(-2)
      x57 = 0.0028125d0*x56
      x58 = sqrt(0.000135388578153653d0*x54 + x55*x57)
      x59 = abs(gamma_aa)
      x60 = x58*x59
      x61 = x53*x60
      x62 = gamma_aa**3
      x63 = rho_a**(-8.0d0)
      x64 = x62*x63
      x65 = gamma_aa**2
      x66 = x54*x65
      x67 = x56*x65
      x68 = x55*x67
      x69 = 0.015625d0*x68 + 1.0d0
      x70 = x69**(-2)
      x71 = 0.02485875d0*x68*x70 + 0.123456790123457d0
      x72 = 0.0204668008035082d0*x71
      x73 = x2*x72 + 1.86977850568312d-6*x64 + 6.38454405015215d-6*x66 +
     & 0.0021416494625995d0*x68
      x74 = x37*(0.0896750813832074d0*x50 - 1.24378109452736d0*x61 + x73
     &) + 1.0d0
      x75 = 1d0/x74
      x76 = -x19 + x23
      x77 = -x25*x76 + x29
      x78 = 4.9628039271952d0*x17
      x79 = 0.804d0*x35
      x80 = -x40 + 1.0d0
      x81 = x2*x80
      x82 = 0.0123414808845154d0*x81 + 0.45d0
      x83 = 0.0274255130767009d0*x81 + 1.0d0
      x84 = x38*x80*x83 + 1.0d0
      x85 = x84**(-0.5d0)
      x86 = x38 - x82*x85
      x87 = x86**2
      x88 = x52 - x85*(0.00222451383844352d0*x81 + 0.0811111111111111d0)
     &
      x89 = x60*x88
      x90 = x37*(x73 + 0.0896750813832074d0*x87 - 1.24378109452736d0*x89
     &) + 1.0d0
      x91 = 1.8610514726982d0*x33
      x92 = rho_a**(-0.666666666666667d0)
      x93 = rho_a**(-4.33333333333333d0)
      x94 = -x10 - x11 + x7
      x95 = rho_a**(-3.66666666666667d0)
      x96 = rho_a**(-0.333333333333333d0)
      x97 = x31/x15**2
      x98 = 1d0/x0
      x99 = 0.00844343197019482d0/x9
      x100 = x24*x26
      x101 = -0.804d0*x75
      x102 = x34**(-2)
      x103 = 0.6203504908994d0*x92
      x104 = rho_a**(-6.33333333333333d0)
      x105 = rho_a**(-3)
      x106 = gamma_aa*x95
      x107 = 0.001553671875d0*x68/x69 - 0.0497175d0
      x108 = 0.0204668008035082d0*x107*x70
      x109 = 1d0/x58
      x110 = 1d0/x59
      x111 = 1.24378109452736d0*x110
      x112 = 1d0/x36
      x113 = 0.0164553078460206d0*x2*x71 + 1.50330191856923d-6*x64 + 5.1
     &3317341632233d-6*x66 + 0.00172188616792999d0*x68
      x114 = tau_a*x1
      x115 = x106*x41
      x116 = 0.00548510261534019d0*x42
      x117 = x43*x47**(-1.5d0)*(-0.0438808209227215d0*x114*x45 + 0.01462
     &69403075738d0 *x115*x45 - x116*(0.219404104613608d0*x114 - 0.07313
     &47015378692d0 *x115))
      x118 = 1.24378109452736d0*x60
      x119 = x17*x37
      x120 = rho_a**(-2.33333333333333d0)
      x121 = x90**(-2)
      x122 = rho_a**(-4.66666666666667d0)
      x123 = rho_a**(-1.66666666666667d0)
      x124 = 0.0177961107075482d0*x123
      x125 = x1*x80
      x126 = 0.0438808209227215d0*x123
      x127 = x123*x39
      x128 = x84**(-1.5d0)*(0.00548510261534019d0*x1*x46 + x116*(x1*x44 
     &- 0.219404104613608d0*x127) - x126*x39*x45)
      x129 = 0.180246913580247d0*x82
      x130 = 2.0d0*x82
      x131 = tau_a**(-3)
      x132 = x131*x55*x65
      x133 = (0.00120345402803247d0*gamma_aa*x80*x93 + x126*x83)/x84
      Ea = -x0*(-4.9628039271952d0*x32 + x33*(1.49628538404935d0*x35 - 1
     &.49628538404935d0*x75))
      d1Ea(1) = -rho_a*(-3.99009435746494d0*gamma_aa*x12*x17*x97*(x1*x5*
     &( -0.284959391993876d0*rho_a**0.333333333333333d0*x4 - 28.80525977
     &71978d0*x20*x8 + 0.562691342264888d0*x96) - 0.138634955389984d0*x4
     &*x94*x95 + 0.0108360841285172d0*x93*x94) + mu*x16*x78*(-0.12827824
     &3853042d0*x18*(0.128278243853042d0*mu*x98 + 10.128444160768d0*x100
     &*x76*x96 + 2.59851805981384d0*x100*x92 + x25*(x3*x99 - 0.085518829
     &2353615d0*x98) - x27*x99) + 0.0427594146176807d0*x77*x98) + 1.1936
     &6207318922d0*mu*x21*x33*( x101 + x79) - 1.65426797573173d0*x32*x92
     & + x33*( 0.0179261670197794d0*gamma_aa*rho_a**(-3.33333333333333d0
     &)*x102 + x103*(x101 + 1.804d0) - x103*(-x79 + 1.804d0) + 1.4962853
     &8404935d0*x119*(-1.49582280454649d-5*rho_a**(-9.0d0)*x62 + rho_a**
     &(-5.66666666666667d0)*x108*x56*x62 - 3.40509016008115d-5 *x104*x65
     & - 0.00428329892519899d0*x105*x67 + 0.135327364628503d0* x106*x112
     &*(x113 + 0.0720987654320988d0*x50 - x61) - 0.0545781354760218d0*x1
     &06*x71 + x109*x111*x53*x65*( 0.000361036208409741d0*x104 + x105*x5
     &7) - x118*( -0.00527292169112538d0*x106 + 0.180246913580247d0*x117
     & + x48*( 0.0177961107075482d0*x114 - 0.00593203690251605d0*x115)) 
     &+ 0.0896750813832074d0*x49*(-0.0585077612302953d0*x106 + 2.0d0*x11
     &7 + 2.0d0*x48*(0.0987318470761234d0*x114 - 0.0329106156920411d0* x
     &115)))/x74**2)) - x17*x91*(x79 - 0.804d0/x90) + x78*(1.804d0 + 0.8
     &04d0/(x14 - 1.0d0))*(x30*x77 - 0.375d0)
      d1Ea(2) = -rho_a*(0.207437457468366d0*x120*x13*x4*x97 - 2.0d0*x33*
     &( 0.00336115631620864d0*x102*x120 + 0.748142692024676d0*x119*x121*
     &( gamma_aa*x111*x58*x88 - 1.27690881003043d-5*gamma_aa*x54 - 0.004
     &28329892519899d0*gamma_aa*x55*x56 + 0.0507477617356885d0*x1* x112*
     &(x113 + 0.0720987654320988d0*x87 - x89) - x1*x72 + x108*x122 *x67 
     &+ x118*(x128*x129 + x51 - x85*(x124*x39 + 0.00222451383844352d0*x1
     &25)) - 5.60933551704935d-6*x63*x65 - 0.0896750813832074d0*x86*(0.0
     &219404104613608d0*x1 + x128*x130 - 2.0d0*x85*(0.0123414808845154d0
     &*x125 + 0.0987318470761234d0*x127 )))))
      d1Ea(3) = -x0*x121*x37*x91*(0.0164553078460206d0*x107*x122*x131*x6
     &2*x70 + 0.0028125d0*x109*x110*x132*x88 - 0.00344377233585999d0*x13
     &2 - x60 *x85*(x124 - x129*x133) + 0.0720987654320988d0*x85*x86*( 0
     &.197463694152247d0*x123 - x130*x133))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 1d0/x0
      x2 = mu**2
      x3 = 1d0/x2
      x4 = rho_a**0.666666666666667d0
      x5 = 15.192666241152d0*x3*x4
      x6 = exp(x5)
      x7 = mu**4*(x6 - 1.0d0)
      x8 = x1*x7
      x9 = 0.000722072416819482d0*x8
      x10 = 7.0d0*x6 + 11.0d0
      x11 = x2/x4
      x12 = 0.00060945584614891d0*x10*x11
      x13 = x12 - x9 + 0.0185185185185185d0
      x14 = 3.0d0*x6
      x15 = x14 - 2.0d0
      x16 = rho_a**0.333333333333333d0
      x17 = 1d0/x16
      x18 = mu*x17
      x19 = 1d0/mu
      x20 = erf(3.8977770897207535d0*x16*x19)
      x21 = x20*x6
      x22 = 1d0/(-0.131642462768164d0*x11*x15 - x14 + 1.81893813843777d0
     &*x18*x21 + 0.00866486900183379d0*x8)
      x23 = rho_a**(-2.0d0)
      x24 = 3.15935190989552d0*x23*x3
      x25 = x22*x24
      x26 = exp(-0.31265084907439067d0*x11)
      x27 = gamma_aa*x26
      x28 = x13*x25*x27
      x29 = -x28 + 1.0d0
      x30 = 1.804d0 - 0.804d0/x29
      x31 = 0.256556487706084d0*x17
      x32 = 1d0/rho_a
      x33 = 0.00844343197019482d0*x32
      x34 = x2*x33
      x35 = exp(-x5)
      x36 = mu*x35
      x37 = mu**3
      x38 = -0.384834731559127d0*x18 + 1.77245385090552d0*x20 + x33*x37
      x39 = 0.128278243853042d0*x18
      x40 = x39*(x36*(x31 - x34) + x38) - 0.375d0
      x41 = x30*x40
      x42 = exp(-1.9241736577956332d0*x18)
      x43 = rho_a**(-2.66666666666667d0)
      x44 = gamma_aa*x43
      x45 = 0.00449266744437808d0*x44 + 1.0d0
      x46 = 1d0/x45
      x47 = 0.0204006002177468d0*x44 + 1.0d0
      x48 = x47**(-2)
      x49 = 0.0109702052306804d0*x44
      x50 = tau_a/gamma_aa
      x51 = 8.0d0*rho_a*x50
      x52 = x51 - 1.0d0
      x53 = x44*x52
      x54 = 0.0123414808845154d0*x53 - 0.45d0
      x55 = 0.0274255130767009d0*x53 - 1.0d0
      x56 = x49*x52*x55 + 1.0d0
      x57 = x56**(-0.5d0)
      x58 = x49 + x54*x57
      x59 = x58**2
      x60 = 0.00197734563417202d0*x43
      x61 = gamma_aa*x60
      x62 = x57*(0.00222451383844352d0*x53 - 0.0811111111111111d0) + x61
     &
      x63 = rho_a**(-5.33333333333333d0)
      x64 = rho_a**(-2)
      x65 = tau_a**(-2)
      x66 = 0.0028125d0*x65
      x67 = sqrt(0.000135388578153653d0*x63 + x64*x66)
      x68 = abs(gamma_aa)
      x69 = x67*x68
      x70 = x62*x69
      x71 = gamma_aa**3
      x72 = rho_a**(-8.0d0)
      x73 = x71*x72
      x74 = gamma_aa**2
      x75 = x63*x74
      x76 = x65*x74
      x77 = x64*x76
      x78 = 0.015625d0*x77 + 1.0d0
      x79 = x78**(-2)
      x80 = 0.02485875d0*x77*x79 + 0.123456790123457d0
      x81 = 0.0204668008035082d0*x80
      x82 = x44*x81 + 1.86977850568312d-6*x73 + 6.38454405015215d-6*x75 
     &+ 0.0021416494625995d0*x77
      x83 = x48*(0.0896750813832074d0*x59 - 1.24378109452736d0*x70 + x82
     &) + 1.0d0
      x84 = 1d0/x83
      x85 = x28 - 1.0d0
      x86 = -x31 + x34
      x87 = -x36*x86 + x38
      x88 = x39*x87 - 0.375d0
      x89 = 4.9628039271952d0*x16
      x90 = 0.804d0*x46
      x91 = -x51 + 1.0d0
      x92 = x44*x91
      x93 = 0.0123414808845154d0*x92 + 0.45d0
      x94 = 0.0274255130767009d0*x92 + 1.0d0
      x95 = x49*x91*x94 + 1.0d0
      x96 = x95**(-0.5d0)
      x97 = x49 - x93*x96
      x98 = x97**2
      x99 = x61 - x96*(0.00222451383844352d0*x92 + 0.0811111111111111d0)
     &
      x100 = x69*x99
      x101 = x48*(-1.24378109452736d0*x100 + x82 + 0.0896750813832074d0*
     &x98) + 1.0d0
      x102 = 1.8610514726982d0*x42
      x103 = rho_a**(-0.666666666666667d0)
      x104 = -0.804d0*x84
      x105 = mu*x1
      x106 = 0.00844343197019482d0*x23
      x107 = 0.564189583547756d0
      x108 = x19*x35
      x109 = rho_a**(-0.333333333333333d0)
      x110 = rho_a**(-3.66666666666667d0)
      x111 = -x12 + x9 - 0.0185185185185185d0
      x112 = x32*x6
      x113 = rho_a**(-1.66666666666667d0)
      x114 = x113*x2
      x115 = x114*x6
      x116 = rho_a**(-2.33333333333333d0)
      x117 = x116*x7
      x118 = x45**(-2)
      x119 = 0.6203504908994d0*x103
      x120 = rho_a**(-6.33333333333333d0)
      x121 = rho_a**(-3)
      x122 = gamma_aa*x110
      x123 = 0.001553671875d0*x77/x78 - 0.0497175d0
      x124 = 0.0204668008035082d0*x123*x79
      x125 = 1d0/x67
      x126 = 1d0/x68
      x127 = 1.24378109452736d0*x126
      x128 = 1d0/x47
      x129 = 0.0164553078460206d0*x44*x80 + 1.50330191856923d-6*x73 + 5.
     &13317341632233d-6*x75 + 0.00172188616792999d0*x77
      x130 = tau_a*x43
      x131 = x122*x52
      x132 = x54*x56**(-1.5d0)*(-0.0438808209227215d0*x130*x55 + 0.01462
     &69403075738d0 *x131*x55 - 0.00548510261534019d0*x53*(0.21940410461
     &3608d0*x130 - 0.0731347015378692d0*x131))
      x133 = 1.24378109452736d0*x69
      x134 = x16*x48
      x135 = x101**(-2)
      x136 = rho_a**(-4.66666666666667d0)
      x137 = 0.0177961107075482d0*x113
      x138 = x43*x91
      x139 = 0.0438808209227215d0*x113*x94
      x140 = x113*x50
      x141 = x95**(-1.5d0)*(0.00548510261534019d0*x138*x94 + x139*x50 + 
     &0.00548510261534019d0*x92*(0.0274255130767009d0*x138 + 0.219404104
     &613608d0*x140))
      x142 = 0.180246913580247d0*x93
      x143 = 2.0d0*x93
      x144 = tau_a**(-3)
      x145 = x144*x64*x74
      x146 = (0.00120345402803247d0*gamma_aa*rho_a**(-4.33333333333333d0
     &)*x91 + x139) /x95
      Ea = -x0*(-4.9628039271952d0*x41 + x42*(1.49628538404935d0*x46 - 1
     &.49628538404935d0*x84))
      d1Ea(1) = -rho_a*(mu*x30*x89*(0.0427594146176807d0*x1*x87 - 0.1282
     &78243853042d0* x17*(4.60575334176456d0*x103*x107*x108 + 0.12827824
     &3853042d0*x105 - x106*x37 + 10.128444160768d0*x108*x109*x86 + x36*
     &( -0.0855188292353615d0*x1 + x106*x2))) + 1.19366207318922d0*mu*x3
     &2 *x42*(x104 + x90) - 1.65426797573173d0*x103*x41 - 3.990094357464
     &94d0*x16*x22*x27*x40*(-6.31870381979104d0*rho_a**( -3.0d0)*x111*x3
     & + 0.658516038102422d0*x110*x111 + x111*x25*( -18.4230133670582d0*
     &x103*x19*x21 + 0.606312712812591d0*x105*x21 - 4.72654360241471d0*x
     &107*x32 + 30.385332482304d0*x109*x3*x6 + 4.0d0*x112 - 0.0877616418
     &45443d0*x114*x15 - 0.087761641845443d0* x115 + 0.0115531586691117d
     &0*x117) + x24*(0.000406303897432606d0* x10*x114 - 0.04320987654320
     &99d0*x112 + 0.00731347015378692d0*x115 - 0.000962763222425976d0*x1
     &17))/x29**2 + x42*( 0.0179261670197794d0*gamma_aa*rho_a**(-3.33333
     &333333333d0)*x118 + x119*(x104 + 1.804d0) - x119*(-x90 + 1.804d0) 
     &+ 1.49628538404935d0*x134*(-1.49582280454649d-5*rho_a**(-9.0d0)*x7
     &1 + rho_a**(-5.66666666666667d0)*x124*x65*x71 - 3.40509016008115d-
     &5 *x120*x74 - 0.00428329892519899d0*x121*x76 + 0.135327364628503d0
     &* x122*x128*(x129 + 0.0720987654320988d0*x59 - x70) - 0.0545781354
     &760218d0*x122*x80 + x125*x127*x62*x74*( 0.000361036208409741d0*x12
     &0 + x121*x66) - x133*( -0.00527292169112538d0*x122 + 0.18024691358
     &0247d0*x132 + x57*( 0.0177961107075482d0*x130 - 0.0059320369025160
     &5d0*x131)) + 0.0896750813832074d0*x58*(-0.0585077612302953d0*x122 
     &+ 2.0d0*x132 + 2.0d0*x57*(0.0987318470761234d0*x130 - 0.0329106156
     &920411d0* x131)))/x83**2)) - x102*x16*(x90 - 0.804d0/x101) + x88*x
     &89*( 1.804d0 + 0.804d0/x85)
      d1Ea(2) = rho_a*(-12.6061122289202d0*rho_a**(-1.66666666666667d0)*
     &x13*x22*x26*x3* x88/x85**2 + 2.0d0*x42*(0.00336115631620864d0*x116
     &*x118 + 0.748142692024676d0*x134*x135*(gamma_aa*x127*x67*x99 - 1.2
     &7690881003043d-5*gamma_aa*x63 - 0.00428329892519899d0*gamma_aa *x6
     &4*x65 + x124*x136*x76 + 0.0507477617356885d0*x128*x43*(-x100 + x12
     &9 + 0.0720987654320988d0*x98) + x133*(x141*x142 + x60 - x96*( x137
     &*x50 + 0.00222451383844352d0*x138)) - x43*x81 - 5.60933551704935d-
     &6*x72*x74 - 0.0896750813832074d0*x97*(x141*x143 + 0.02194041046136
     &08d0*x43 - 2.0d0*x96*(0.0123414808845154d0*x138 + 0.09873184707612
     &34d0*x140)))))
      d1Ea(3) = -x0*x102*x135*x48*(0.0164553078460206d0*x123*x136*x144*x
     &71*x79 + 0.0028125d0*x125*x126*x145*x99 - 0.00344377233585999d0*x1
     &45 - x69 *x96*(x137 - x142*x146) + 0.0720987654320988d0*x96*x97*( 
     &0.197463694152247d0*x113 - x143*x146))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**(-2.66666666666667d0)
      x1 = gamma_aa*x0
      x2 = mu**2
      x3 = 1d0/x2
      x4 = 0.0519881082712438d0*x3
      x5 = rho_a**0.666666666666667d0
      x6 = 0.844037013397333d0*x5
      x7 = mu**(-4)
      x8 = 14.4026298885989d0*rho_a**2.0d0*x7
      x9 = rho_a**1.33333333333333d0
      x10 = 0.213719543995407d0*x3*x9
      x11 = exp(-0.31265084907439067d0*x2/x5)
      x12 = x11*(x10 - x6 + x8)
      x13 = x1*x12*x4
      x14 = -x13 + 1.0d0
      x15 = 1d0/x14
      x16 = rho_a**1.0d0
      x17 = x16*x3
      x18 = rho_a**0.333333333333333d0
      x19 = exp(-1.9241736577956332d0*mu/x18)
      x20 = 0.00449266744437808d0*x1 + 1.0d0
      x21 = 0.804d0/x20
      x22 = 0.0204006002177468d0*x1 + 1.0d0
      x23 = x22**(-2)
      x24 = 0.0109702052306804d0*x1
      x25 = tau_a/gamma_aa
      x26 = 8.0d0*rho_a*x25
      x27 = x26 - 1.0d0
      x28 = x1*x27
      x29 = 0.0123414808845154d0*x28 - 0.45d0
      x30 = 0.0274255130767009d0*x28 - 1.0d0
      x31 = x24*x27*x30 + 1.0d0
      x32 = x31**(-0.5d0)
      x33 = x24 + x29*x32
      x34 = x33**2
      x35 = 0.00197734563417202d0*x0
      x36 = gamma_aa*x35
      x37 = x32*(0.00222451383844352d0*x28 - 0.0811111111111111d0) + x36
     &
      x38 = rho_a**(-5.33333333333333d0)
      x39 = rho_a**(-2)
      x40 = tau_a**(-2)
      x41 = 0.0028125d0*x40
      x42 = sqrt(0.000135388578153653d0*x38 + x39*x41)
      x43 = abs(gamma_aa)
      x44 = x42*x43
      x45 = x37*x44
      x46 = gamma_aa**3
      x47 = rho_a**(-8.0d0)
      x48 = x46*x47
      x49 = gamma_aa**2
      x50 = x38*x49
      x51 = x40*x49
      x52 = x39*x51
      x53 = 0.015625d0*x52 + 1.0d0
      x54 = x53**(-2)
      x55 = 0.02485875d0*x52*x54 + 0.123456790123457d0
      x56 = 0.0204668008035082d0*x55
      x57 = x1*x56 + 1.86977850568312d-6*x48 + 6.38454405015215d-6*x50 +
     & 0.0021416494625995d0*x52
      x58 = x23*(0.0896750813832074d0*x34 - 1.24378109452736d0*x45 + x57
     &) + 1.0d0
      x59 = -0.804d0/x58
      x60 = x19*(x21 + x59)
      x61 = 3.14159265358979d0*x17*(-0.804d0*x15 + 1.804d0) + 1.86105147
     &26982d0*x18* x60
      x62 = rho_a**(-4.33333333333333d0)
      x63 = -x10 + x6 - x8
      x64 = rho_a**(-3.66666666666667d0)
      x65 = x20**(-2)
      x66 = 0.6203504908994d0*rho_a**(-0.666666666666667d0)
      x67 = rho_a**(-6.33333333333333d0)
      x68 = rho_a**(-3)
      x69 = gamma_aa*x64
      x70 = 0.001553671875d0*x52/x53 - 0.0497175d0
      x71 = 0.0204668008035082d0*x54*x70
      x72 = 1d0/x42
      x73 = 1d0/x43
      x74 = 1.24378109452736d0*x73
      x75 = 1d0/x22
      x76 = 0.0164553078460206d0*x1*x55 + 1.50330191856923d-6*x48 + 5.13
     &317341632233d-6*x50 + 0.00172188616792999d0*x52
      x77 = tau_a*x0
      x78 = x27*x69
      x79 = x29*x31**(-1.5d0)*(-0.00548510261534019d0*x28*(0.21940410461
     &3608d0*x77 - 0.0731347015378692d0*x78) - 0.0438808209227215d0*x30*
     &x77 + 0.0146269403075738d0*x30*x78)
      x80 = 1.24378109452736d0*x44
      x81 = x18*x23
      x82 = rho_a**(-1.66666666666667d0)
      x83 = -x26 + 1.0d0
      x84 = x1*x83
      x85 = 0.0123414808845154d0*x84 + 0.45d0
      x86 = 0.0274255130767009d0*x84 + 1.0d0
      x87 = x24*x83*x86 + 1.0d0
      x88 = x87**(-0.5d0)
      x89 = x24 - x85*x88
      x90 = x89**2
      x91 = x36 - x88*(0.00222451383844352d0*x84 + 0.0811111111111111d0)
     &
      x92 = x44*x91
      x93 = (x23*(x57 + 0.0896750813832074d0*x90 - 1.24378109452736d0*x9
     &2) + 1.0d0) **(-2)
      x94 = rho_a**(-4.66666666666667d0)
      x95 = 0.0177961107075482d0*x82
      x96 = x0*x83
      x97 = 0.0438808209227215d0*x82*x86
      x98 = x25*x82
      x99 = x87**(-1.5d0)*(x25*x97 + 0.00548510261534019d0*x84*(0.027425
     &5130767009d0 *x96 + 0.219404104613608d0*x98) + 0.00548510261534019
     &d0*x86*x96)
      x100 = 0.180246913580247d0*x85
      x101 = 2.0d0*x85
      x102 = tau_a**(-3)
      x103 = x102*x39*x49
      x104 = (0.00120345402803247d0*gamma_aa*x62*x83 + x97)/x87
      Ea = -rho_a*x61
      d1Ea(1) = -rho_a*(2.52584049348619d0*gamma_aa*x11*x17*(x0*x4*(0.56
     &2691342264888d0* rho_a**(-0.333333333333333d0) - 0.284959391993876
     &d0*rho_a** 0.333333333333333d0*x3 - 28.8052597771978d0*x16*x7) - 0
     &.138634955389984d0*x3*x63*x64 + 0.0108360841285172d0*x62*x63)/ x14
     &**2 + 1.19366207318922d0*mu*x60/x16 + x19*( 0.0179261670197794d0*g
     &amma_aa*rho_a**(-3.33333333333333d0)*x65 - x66*(-x21 + 1.804d0) + 
     &x66*(x59 + 1.804d0) + 1.49628538404935d0* x81*(-1.49582280454649d-
     &5*rho_a**(-9.0d0)*x46 + rho_a**( -5.66666666666667d0)*x40*x46*x71 
     &+ 0.0896750813832074d0*x33*( 2.0d0*x32*(0.0987318470761234d0*x77 -
     & 0.0329106156920411d0*x78) - 0.0585077612302953d0*x69 + 2.0d0*x79)
     & + x37*x49*x72*x74*(x41*x68 + 0.000361036208409741d0*x67) - 3.4050
     &9016008115d-5*x49*x67 - 0.00428329892519899d0*x51*x68 - 0.05457813
     &54760218d0*x55*x69 + 0.135327364628503d0*x69*x75*(0.07209876543209
     &88d0*x34 - x45 + x76 ) - x80*(x32*(0.0177961107075482d0*x77 - 0.00
     &593203690251605d0* x78) - 0.00527292169112538d0*x69 + 0.1802469135
     &80247d0*x79))/x58 **2) + x3*(-2.52584049348619d0*x15 + 5.667433147
     &07599d0)) - x61
      d1Ea(2) = rho_a*(0.131313669051252d0*x12*x7*x82/(x13 - 1.0d0)**2 +
     & 2.0d0*x19*( 0.00336115631620864d0*rho_a**(-2.33333333333333d0)*x6
     &5 + 0.748142692024676d0*x81*x93*(-1.27690881003043d-5*gamma_aa*x38
     & - 0.00428329892519899d0*gamma_aa*x39*x40 + gamma_aa*x42*x74*x91 -
     & x0*x56 + 0.0507477617356885d0*x0*x75*(x76 + 0.0720987654320988d0*
     & x90 - x92) - 5.60933551704935d-6*x47*x49 + x51*x71*x94 + x80*( x1
     &00*x99 + x35 - x88*(x25*x95 + 0.00222451383844352d0*x96)) - 0.0896
     &750813832074d0*x89*(0.0219404104613608d0*x0 + x101*x99 - 2.0d0*x88
     &*(0.0123414808845154d0*x96 + 0.0987318470761234d0*x98 )))))
      d1Ea(3) = -1.8610514726982d0*x19*x23*x9*x93*(0.0164553078460206d0*
     &x102*x46*x54*x70 *x94 + 0.0028125d0*x103*x72*x73*x91 - 0.003443772
     &33585999d0*x103 - x44*x88*(-x100*x104 + x95) + 0.0720987654320988d
     &0*x88*x89*( -x101*x104 + 0.197463694152247d0*x82))
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**0.333333333333333d0
      x1 = mu/x0
      x2 = exp(-1.9241736577956332d0*x1)
      x3 = rho_a**(-2.66666666666667d0)
      x4 = gamma_aa*x3
      x5 = 0.00449266744437808d0*x4 + 1.0d0
      x6 = 0.804d0/x5
      x7 = 0.0204006002177468d0*x4 + 1.0d0
      x8 = x7**(-2)
      x9 = 0.0109702052306804d0*x4
      x10 = 1d0/gamma_aa
      x11 = tau_a*x10
      x12 = 8.0d0*rho_a*x11
      x13 = x12 - 1.0d0
      x14 = x13*x4
      x15 = 0.0123414808845154d0*x14 - 0.45d0
      x16 = 0.0274255130767009d0*x14 - 1.0d0
      x17 = x13*x16*x9 + 1.0d0
      x18 = x17**(-0.5d0)
      x19 = x15*x18 + x9
      x20 = x19**2
      x21 = 0.00197734563417202d0*x3
      x22 = gamma_aa*x21
      x23 = x18*(0.00222451383844352d0*x14 - 0.0811111111111111d0) + x22
     &
      x24 = rho_a**(-5.33333333333333d0)
      x25 = rho_a**(-2)
      x26 = tau_a**(-2)
      x27 = 0.0028125d0*x26
      x28 = 0.000135388578153653d0*x24 + x25*x27
      x29 = sqrt(x28)
      x30 = abs(gamma_aa)
      x31 = x29*x30
      x32 = x23*x31
      x33 = 1.24378109452736d0*x32
      x34 = gamma_aa**3
      x35 = rho_a**(-8.0d0)
      x36 = x34*x35
      x37 = gamma_aa**2
      x38 = x24*x37
      x39 = x25*x26
      x40 = x37*x39
      x41 = 0.015625d0*x40 + 1.0d0
      x42 = x41**(-2)
      x43 = 0.02485875d0*x40*x42 + 0.123456790123457d0
      x44 = x4*x43
      x45 = 1.86977850568312d-6*x36 + 6.38454405015215d-6*x38 + 0.002141
     &6494625995d0 *x40 + 0.0204668008035082d0*x44
      x46 = x8*(0.0896750813832074d0*x20 - x33 + x45) + 1.0d0
      x47 = 1d0/x46
      x48 = -0.804d0*x47
      x49 = x2*(x48 + x6)
      x50 = mu**2
      x51 = rho_a**(-0.666666666666667d0)*x50
      x52 = 1.18478216491348d0*x51 - 7.0d0
      x53 = 0.0519881082712438d0*x52
      x54 = mu**4
      x55 = rho_a**1.33333333333333d0
      x56 = 1d0/x55
      x57 = 49.1113297378198d0*x1 - 10.6630394842213d0*x51 + 0.233951463
     &049512d0*x54 *x56 - 81.0d0
      x58 = 1d0/x57
      x59 = x58*exp(-0.31265084907439067d0*x51)
      x60 = x4*x53*x59 + 1.0d0
      x61 = 1.804d0 - 0.804d0/x60
      x62 = 1.8610514726982d0*x55
      x63 = 1.8610514726982d0*x0
      x64 = rho_a**(-0.666666666666667d0)
      x65 = 0.6203504908994d0*x64
      x66 = rho_a**(-4.33333333333333d0)
      x67 = x50*x66
      x68 = rho_a**(-3.66666666666667d0)
      x69 = x52*x68
      x70 = x52*x67
      x71 = rho_a**(-1.66666666666667d0)
      x72 = mu*x71
      x73 = rho_a**(-2.33333333333333d0)
      x74 = mu**3
      x75 = x73*x74
      x76 = 16.3704432459399d0*x56 - 7.10869298948088d0*x72 + 0.31193528
     &4066016d0* x75
      x77 = x58*x76
      x78 = mu*x3
      x79 = x53*x78
      x80 = -0.0410630556449071d0*x67 - 0.138634955389984d0*x69 + 0.0108
     &360841285172d0*x70 + x77*x79
      x81 = 1.49628538404935d0*x0
      x82 = x59/x60**2
      x83 = x81*x82
      x84 = x80*x83
      x85 = mu*1d0/rho_a
      x86 = 1.19366207318922d0*x85
      x87 = x5**(-2)
      x88 = rho_a**(-3.33333333333333d0)
      x89 = x87*x88
      x90 = -x6 + 1.804d0
      x91 = x46**(-2)
      x92 = rho_a**(-6.33333333333333d0)
      x93 = x37*x92
      x94 = rho_a**(-9.0d0)
      x95 = x34*x94
      x96 = rho_a**(-3)
      x97 = x26*x96
      x98 = x37*x97
      x99 = gamma_aa*x68
      x100 = 0.0545781354760218d0*x43
      x101 = x40/x41
      x102 = 0.001553671875d0*x101
      x103 = x102 - 0.0497175d0
      x104 = x103*x34
      x105 = rho_a**(-5.66666666666667d0)
      x106 = x26*x42
      x107 = x105*x106
      x108 = 0.0204668008035082d0*x107
      x109 = 1d0/x29
      x110 = 1d0/x30
      x111 = x110*x23
      x112 = x109*x111
      x113 = x27*x96 + 0.000361036208409741d0*x92
      x114 = x113*x37
      x115 = x112*x114
      x116 = 1d0/x7
      x117 = 1.50330191856923d-6*x36 + 5.13317341632233d-6*x38 + 0.00172
     &188616792999d0*x40 + 0.0164553078460206d0*x44
      x118 = x117 + 0.0720987654320988d0*x20 - x32
      x119 = x116*x118
      x120 = 0.135327364628503d0*x119
      x121 = 0.00527292169112538d0*x68
      x122 = tau_a*x3
      x123 = x13*x99
      x124 = x17**(-1.5d0)
      x125 = 0.0438808209227215d0*x122
      x126 = 0.0146269403075738d0*x16
      x127 = 0.219404104613608d0*x122 - 0.0731347015378692d0*x123
      x128 = 0.00548510261534019d0*x14
      x129 = x123*x126 - x125*x16 - x127*x128
      x130 = x124*x129
      x131 = x130*x15
      x132 = -gamma_aa*x121 + 0.180246913580247d0*x131 + x18*(0.01779611
     &07075482d0* x122 - 0.00593203690251605d0*x123)
      x133 = x132*x31
      x134 = 0.0585077612302953d0*x68
      x135 = 0.0987318470761234d0*x122 - 0.0329106156920411d0*x123
      x136 = 2.0d0*x18
      x137 = 2.0d0*x15
      x138 = -gamma_aa*x134 + x130*x137 + x135*x136
      x139 = x138*x19
      x140 = -x100*x99 + x104*x108 + 1.24378109452736d0*x115 + x120*x99 
     &- 1.24378109452736d0*x133 + 0.0896750813832074d0*x139 - 3.40509016
     &008115d-5*x93 - 1.49582280454649d-5*x95 - 0.00428329892519899d0*x9
     &8
      x141 = x8*x81
      x142 = x2*(0.0179261670197794d0*gamma_aa*x89 + x140*x141*x91 - x65
     &*x90 + x65*( x48 + 1.804d0))
      x143 = x52*x73
      x144 = 0.0777890465506374d0*x143*x82
      x145 = x73*x87
      x146 = 0.00336115631620864d0*x145
      x147 = -x12 + 1.0d0
      x148 = x147*x4
      x149 = 0.0123414808845154d0*x148 + 0.45d0
      x150 = 0.0274255130767009d0*x148 + 1.0d0
      x151 = x147*x150*x9 + 1.0d0
      x152 = x151**(-0.5d0)
      x153 = -x149*x152 + x9
      x154 = 2.0d0*x149
      x155 = 0.0438808209227215d0*x71
      x156 = x155*x16
      x157 = x11*x156
      x158 = x13*x3
      x159 = 0.00548510261534019d0*x158
      x160 = x159*x16
      x161 = x11*x71
      x162 = 0.219404104613608d0*x161
      x163 = 0.0274255130767009d0*x158 - x162
      x164 = x128*x163
      x165 = -x157 + x160 + x164
      x166 = x151**(-1.5d0)
      x167 = x165*x166
      x168 = 0.0219404104613608d0*x3
      x169 = 0.0987318470761234d0*x161
      x170 = x147*x3
      x171 = x169 + 0.0123414808845154d0*x170
      x172 = 2.0d0*x171
      x173 = -x152*x172 + x168
      x174 = x154*x167 + x173
      x175 = x153*x174
      x176 = 0.180246913580247d0*x149
      x177 = 0.0177961107075482d0*x71
      x178 = x11*x177
      x179 = -x152*(0.00222451383844352d0*x170 + x178) + x21
      x180 = x167*x176 + x179
      x181 = 1.24378109452736d0*x31
      x182 = x35*x37
      x183 = x3*x43
      x184 = 1.27690881003043d-5*x24
      x185 = x153**2
      x186 = -x152*(0.00222451383844352d0*x148 + 0.0811111111111111d0) +
     & x22
      x187 = x186*x31
      x188 = x117 + 0.0720987654320988d0*x185 - x187
      x189 = x116*x3
      x190 = 0.0507477617356885d0*x189
      x191 = 0.00428329892519899d0*x39
      x192 = x110*x186
      x193 = 1.24378109452736d0*gamma_aa
      x194 = x193*x29
      x195 = 0.0204668008035082d0*x106
      x196 = rho_a**(-4.66666666666667d0)
      x197 = x103*x37
      x198 = x196*x197
      x199 = -gamma_aa*x184 - gamma_aa*x191 - 5.60933551704935d-6*x182 -
     & 0.0204668008035082d0*x183 + x188*x190 + x192*x194 + x195*x198
      x200 = -0.0896750813832074d0*x175 + x180*x181 + x199
      x201 = x8*(0.0896750813832074d0*x185 - 1.24378109452736d0*x187 + x
     &45) + 1.0d0
      x202 = x201**(-2)
      x203 = x202*x8
      x204 = x200*x203
      x205 = 0.748142692024676d0*x0
      x206 = 2.0d0*x2
      x207 = 0.197463694152247d0*x71
      x208 = 1d0/x151
      x209 = 0.00120345402803247d0*x66
      x210 = x150*x155
      x211 = x208*(gamma_aa*x147*x209 + x210)
      x212 = 0.0720987654320988d0*x152*x153
      x213 = x152*x31
      x214 = tau_a**(-3)
      x215 = x214*x37
      x216 = x215*x25
      x217 = 0.0028125d0*x216
      x218 = x109*x192
      x219 = 0.0164553078460206d0*x104
      x220 = x214*x42
      x221 = x196*x219*x220 - 0.00344377233585999d0*x216
      x222 = x217*x218 + x221
      x223 = x212*(-x154*x211 + x207) - x213*(-x176*x211 + x177) + x222
      x224 = x2*x8
      x225 = x202*x224
      x226 = gamma_aa*x80
      x227 = x226*x82
      x228 = 0.413566993932933d0*rho_a**(-1.66666666666667d0)
      x229 = 0.997523589366235d0*x64
      x230 = 0.277269910779967d0*x69
      x231 = x57**(-2)
      x232 = x231*exp(-0.62530169814878134d0*x51)/x60**3
      x233 = x24*x50
      x234 = rho_a**(-6.0d0)
      x235 = x234*x54
      x236 = x66*x74*x77
      x237 = mu*x56
      x238 = gamma_aa*x66
      x239 = x5**(-3)
      x240 = 1d0/x201
      x241 = x43*x99
      x242 = 0.270654729257005d0*x99
      x243 = 1/(x201**3*x7**4)
      x244 = x140*x243*x81
      x245 = rho_a**(-7.33333333333333d0)
      x246 = x245*x37
      x247 = rho_a**(-4)
      x248 = x247*x26
      x249 = gamma_aa*x196
      x250 = rho_a**(-6.66666666666667d0)
      x251 = tau_a**(-4)
      x252 = x251*x42
      x253 = 0.00014565673828125d0*gamma_aa**4*x247*x252
      x254 = x34*(-0.010875703125d0*x101 + x253 + 0.1491525d0)
      x255 = x28**(-1.5d0)*x30
      x256 = 1.24378109452736d0*x23
      x257 = x118*x8
      x258 = x110*x132
      x259 = x109*x114
      x260 = 1.24378109452736d0*x138
      x261 = tau_a*x68
      x262 = x13*x249
      x263 = x130*x135
      x264 = x124*(-0.087761641845443d0*x122*x127 + 0.0292538806151477d0
     &*x123*x127 - x128*(-1.17015522460591d0*x261 + 0.26816057230552d0*x
     &262) + 0.234031044921181d0*x16*x261 - 0.0536321144611041d0*x16*x26
     &2)
      x265 = 0.180246913580247d0*x15
      x266 = 0.131642462768165d0*x16
      x267 = 0.0164553078460206d0*x14
      x268 = x129*x17**(-2.5d0)
      x269 = x268*(-x122*x266 + 0.0438808209227215d0*x123*x16 - x127*x26
     &7)
      x270 = 0.000790938253668807d0*x3
      x271 = gamma_aa*x270
      x272 = x18*(0.000889805535377408d0*x14 - 0.0324444444444444d0) + x
     &271
      x273 = 1.24378109452736d0*x272
      x274 = x107*x219 + x115 - x133 + 0.0720987654320988d0*x139 - 0.043
     &8808209227215d0*x241 - 2.73769248870524d-5*x93 - 1.20264153485538d
     &-5*x95 - 0.00344377233585999d0*x98
      x275 = x141*x202
      x276 = x166*(x11*x210 + 0.00548510261534019d0*x148*(x162 + 0.02742
     &55130767009d0 *x170) + 0.00548510261534019d0*x150*x170)
      x277 = -0.0896750813832074d0*x153*(x154*x276 + x173) + x181*(x176*
     &x276 + x179) + x199
      x278 = x203*x205
      x279 = gamma_aa*x24
      x280 = 0.00856659785039798d0*gamma_aa
      x281 = 0.0409336016070163d0*x106
      x282 = gamma_aa*x29
      x283 = 2.48756218905473d0*x282
      x284 = 0.101495523471377d0*x189
      x285 = -0.179350162766415d0*x175 + 2.48756218905473d0*x180*x31 - 1
     &.12186710340987d-5*x182 - 0.0409336016070163d0*x183 + x188*x284 + 
     &x192*x283 + x198*x281 - 2.55381762006086d-5*x279 - x280*x39
      x286 = gamma_aa*x92
      x287 = -x102 + 0.0497175d0
      x288 = x287*x37
      x289 = 0.00932203125d0*x101
      x290 = x10*x109*x30
      x291 = x157 - x160 - x164
      x292 = x124*x291
      x293 = x18*(0.00222451383844352d0*x158 - x178) + x21 + x265*x292
      x294 = x110*x293
      x295 = 1.24378109452736d0*x259
      x296 = 0.00711844428301927d0*x71
      x297 = 0.0720987654320988d0*x15
      x298 = x18*(-x11*x296 + 0.000889805535377408d0*x158) + x270 + x292
     &*x297
      x299 = x111*x29
      x300 = 0.0123414808845154d0*x158 - x169
      x301 = x136*x300 + x137*x292 + x168
      x302 = 0.0720987654320988d0*x19
      x303 = -gamma_aa*x299 + 0.00344377233585999d0*gamma_aa*x39 + 0.016
     &4553078460206d0*x106*x196*x288 + 4.50990575570768d-6*x182 + 0.0164
     &553078460206d0*x183 + 1.02663468326447d-5*x279 - x293*x31 + x301*x
     &302
      x304 = 0.135327364628503d0*x116*x99
      x305 = x10*x122
      x306 = x13*x68
      x307 = x130*x300
      x308 = x266*x71
      x309 = x11*x308 - 0.0164553078460206d0*x158*x16 - x163*x267
      x310 = x268*x309
      x311 = x127*x155
      x312 = x124*(x11*x311 + 0.0146269403075738d0*x123*x163 - x125*x163
     & + x126*x306 - x127*x159 - x128*(0.585077612302953d0*x305 - 0.0731
     &347015378692d0*x306) - 0.117015522460591d0*x16*x305)
      x313 = 2.0d0*x135
      x314 = rho_a**(-5.0d0)
      x315 = x15*x163*(-0.0109702052306804d0*x158 + 0.087761641845443d0*
     &x161)
      x316 = x291*x300
      x317 = 1d0/x17
      x318 = x291*x309
      x319 = x15*x317
      x320 = x13*x209
      x321 = gamma_aa*x320 + x156
      x322 = x208*x321
      x323 = x154*x322 + x207
      x324 = 0.180246913580247d0*x321
      x325 = x149*x208*x324 + x177
      x326 = x212*x323 - x213*x325 + x222
      x327 = x317*x321
      x328 = x137*x327
      x329 = x319*x324
      x330 = x18*(x177 - x329)
      x331 = x112*x217 + x18*x302*(x207 - x328) + x221 - x31*x330
      x332 = x215*x96
      x333 = x105*x220
      x334 = x253 - x289 + 0.099435d0
      x335 = x214*x25
      x336 = x109*x216
      x337 = 0.00962763222425977d0*tau_a
      x338 = x317*(0.00521496745480737d0*x13*x279 + 0.0731347015378692d0
     &*x16*x3 - x311 - x337*x66)
      x339 = x129*x317
      x340 = 0.00361036208409741d0*x13*x238 + x308
      x341 = x129*x340/x17**2
      x342 = gamma_aa*x335
      x343 = 0.0153120874502694d0*x196*x215*x42
      x344 = 0.9305257363491d0*x152
      x345 = x110*x325
      x346 = 0.0720987654320988d0*x149*x322 + x296
      x347 = x149*(-x10*x337*x88 + x155*x163 + x320)
      x348 = x165*x208*x340
      x349 = 0.9305257363491d0*x166
      x350 = -x152*(0.000889805535377408d0*x148 + 0.0324444444444444d0) 
     &+ x271
      x351 = x240*x8
      x352 = x25*x251*x37
      x353 = x15*x88
      x354 = x321*x71
      x355 = 1.8610514726982d0*x124
      Ea = -x62*(x49 + x61)
      d1Ea(1) = -rho_a*(gamma_aa*x84 + x142 + x49*x86 + x61*x65) - x49*x
     &63 - x61*x63
      d1Ea(2) = -rho_a*(x144 - x206*(x146 + x204*x205))
      d1Ea(3) = -x223*x225*x62
      d2Ea(1) = rho_a*(-gamma_aa*x83*(-mu*x230*x77 + 0.508328169763273d0
     &*x196*x52 + x231 *x3*x50*x53*x76*(32.7408864918799d0*x56 - 14.2173
     &859789618d0*x72 + 0.623870568132033d0*x75) - 0.0758525888996201d0*
     &x233*x52 + 0.28744138951435d0*x233 + 0.00225860726894828d0*x235*x5
     &2 - 0.0171178656172922d0*x235 + 0.0216721682570343d0*x236*x52 - 0.
     &0821261112898142d0*x236 + x58*x79*(-21.8272576612533d0*x73 - 0.727
     &848996154038d0*x74*x88 + 11.8478216491348d0*x78)) + 1.591549430918
     &95d0*mu*rho_a**(-2.0d0)*x49 - 1.28278243853042d0* x142*x237 + x2*(
     &-0.00042952697051865d0*rho_a**(-7.0d0)*x239*x37 - x140*x203*x229 -
     & x228*x90 + x228*(-0.804d0*x240 + 1.804d0) + 0.0537785010593382d0*
     &x238*x87 + x244*(0.0409336016070163d0*x104* x107 + 2.4875621890547
     &3d0*x115 + x119*x242 - 2.48756218905473d0* x133 + 0.17935016276641
     &5d0*x139 - 0.109156270952044d0*x241 - 6.81018032016229d-5*x93 - 2.
     &99164560909299d-5*x95 - 0.00856659785039798d0*x98) - x275*(0.00013
     &4624052409184d0*rho_a** (-10.0d0)*x34 - 0.109156270952044d0*x104*x
     &106*x250 - 1.24378109452736d0*x112*x37*(0.00228656265326169d0*x245
     & + 0.0084375d0*x248) + x113**2*x255*x256 + x116*x242*x274 - 0.4962
     &00336971176d0*x119*x249 - x181*(x18*(-0.0949125904402569d0* x261 +
     & 0.0217508019758922d0*x262) + 0.0193340462007931d0*x249 + 0.360493
     &827160494d0*x263 + x264*x265 + x265*x269) + x195*x250* x254 + 0.02
     &20860757144586d0*x246*x257 + 0.000215655710138473d0* x246 + 0.0128
     &49896775597d0*x248*x37 + 0.200119830078746d0*x249* x43 + 2.4875621
     &8905473d0*x258*x259 + x260*(0.0720987654320988d0* x131 + x18*(0.00
     &711844428301927d0*x122 - 0.00237281476100642d0* x123) - 0.00210916
     &867645015d0*x99) + x273*(x136*( -0.526569851072658d0*x261 + 0.1206
     &72257537484d0*x262) + x137*x264 + x137*x269 + 0.214528457844416d0*
     &x249 + 4.0d0*x263))) - x227* x229 + x228*x61 + x232*x37*x80*x81*(-
     &x230 + 0.103976216542488d0* x52*x77*x78 - 0.0821261112898142d0*x67
     & + 0.0216721682570343d0*x70 ) - 0.76560437251347d0*x49*x50*x73) - 
     &2.99257076809871d0*x0*x227 - 2.0d0*x142 - 2.38732414637843d0*x49*x
     &85 - 1.2407009817988d0*x61 *x64
      d2Ea(2) = rho_a*(0.155578093101275d0*x143*x226*x232 + 0.6413912192
     &65211d0*x2*x237* (0.00672231263241728d0*x145 + x275*x277) - x2*( -
     &0.000161072613944494d0*gamma_aa*x234*x239 - 0.498761794683118d0* x
     &204*x64 + x244*x285 + x275*(-x100*x68 - 0.0545781354760218d0* x107
     &*x288 + x108*x197 + x108*x37*(-x253 + x289 - 0.099435d0) + x112*x1
     &93*(0.000722072416819483d0*x92 + 0.005625d0*x97) - x113* x256*x290
     & + x120*x68 - x181*(-x121 + 0.180246913580247d0*x135* x292 + x18*(
     &0.0474562952201284d0*x305 - 0.00593203690251605d0* x306) + x265*x3
     &10 + x265*x312 + 0.180246913580247d0*x307) - x190* x274 - x194*x25
     &8 - 0.00828227839292197d0*x257*x286 + x260*x298 + x273*(-x134 + x1
     &36*(0.263284925536329d0*x305 - 0.0329106156920411d0*x306) + x137*x
     &310 + x137*x312 + x292*x313 + 2.0d0*x307) - x280*x97 - 6.810180320
     &16229d-5*x286 + x294*x295 + x303*x304 - 4.48746841363948d-5*x37*x9
     &4) + 0.015685396142307d0* x89) - 0.0259296821835458d0*x52*x82*x88 
     &- x84) - x144 + x206*( x146 + x277*x278)
      d2Ea(3) = rho_a*(-x206*(-x200*x205*x243*x285 + 3.02011151145926d-5
     &*x239*x314 + x278*(1.12186710340987d-5*gamma_aa*x35 - x124*x181*(x
     &265*x317* x318 + 0.180246913580247d0*x315 + 0.360493827160494d0*x3
     &16) + 2.48756218905473d0*x124*x272*(x315 + 2.0d0*x316 + x318*x319)
     & + x184 + x191 + x195*x249*(-0.007768359375d0*x101 + x253 + 0.0497
     &175d0) + 0.00310585439734574d0*x24*x257 + x249*x281*x287 - x283*x2
     &94 - x284*x303 + 1.24378109452736d0*x298*x301 - 1.24378109452736d0
     &*x299 + x33/x37)) + 0.00808821074878273d0*x232* x314*x52**2)
      d2Ea(4) = -x224*(rho_a*x91*(-3.7221029453964d0*x0*x140*x331*x47*x8
     & + x331*x65 + x331*x86 + x81*(-0.0545781354760218d0*x104*x333 + x1
     &10*x295*x330 - 0.00699626865671642d0*x112*x332 + 0.003498134328358
     &21d0*x113* x23*x255*x335 - x18*x181*(-x135*x317*x324 + x177*x339 +
     & x265*x338 - x265*x341 - 0.0296601845125803d0*x3) + x18*x260*(x296
     & - x297* x327) + x18*x273*(x137*x338 - x137*x341 + x207*x339 - 0.3
     &29106156920411d0*x3 - x313*x327) + 0.00349813432835821d0*x258* x33
     &6 + x304*x331 + 0.00856659785039798d0*x332 + 0.0204668008035082d0*
     &x333*x334*x34)) + x202*x326*x63)
      d2Ea(5) = x203*x206*x55*(-x103*x343 - 0.00261710363348184d0*x110*x
     &180*x336 - x174* x344*x346 + 0.00261710363348184d0*x186*x290*x335 
     &+ 0.0379665670791648d0*x189*x223 - 0.00523420726696369d0*x218*x342
     & - 1.8610514726982d0*x223*x277*x351 + x282*x344*x345 + x31*x349*( 
     &-x165*x177 + x171*x324 - x176*x348 + 0.180246913580247d0*x347) + x
     &334*x343 + 0.00640903757728955d0*x342 - x349*x350*(-x154*x348 - x1
     &65*x207 + x172*x321 + 2.0d0*x347))
      d2Ea(6) = x225*x55*(-0.0104684145339274d0*x152*x336*x345 - 0.03062
     &41749005388d0* x196*x252*x254 - 1.8610514726982d0*x208*x323*x346 +
     & 0.0157026218008911d0*x218*x352 - x31*x355*(-x329*x340 + 0.0034707
     &019870171d0*x353 + 0.0355922214150963d0*x354) + 4.62948127536866d0
     &*x326**2*x351 + x350*x355*(-x328*x340 + 0.0385105288970391d0*x353 
     &+ 0.394927388304494d0*x354) - 0.0192271127318687d0*x352 - 1.472120
     &79383354d-5*x186*x247*x255/ tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = gamma_aa*x1
      x3 = mu**2
      x4 = rho_a**0.666666666666667d0
      x5 = x3/x4
      x6 = 1.18478216491348d0*x5 - 7.0d0
      x7 = 0.0519881082712438d0*x6
      x8 = rho_a**0.333333333333333d0
      x9 = 1d0/x8
      x10 = mu*x9
      x11 = mu**4
      x12 = 1d0/x0
      x13 = 49.1113297378198d0*x10 + 0.233951463049512d0*x11*x12 - 10.66
     &30394842213d0*x5 - 81.0d0
      x14 = 1d0/x13
      x15 = x14*exp(-0.31265084907439067d0*x5)
      x16 = x15*x2*x7 + 1.0d0
      x17 = 1.804d0 - 0.804d0/x16
      x18 = 0.256556487706084d0*x9
      x19 = 1d0/rho_a
      x20 = 0.00844343197019482d0*x19
      x21 = x20*x3
      x22 = x18 - x21
      x23 = exp(-15.192666241151986d0*x4/x3)
      x24 = mu*x23
      x25 = 1d0/mu
      x26 = mu**3
      x27 = -0.384834731559127d0*x10 + x20*x26 + 1.77245385090552d0*erf(
     & 3.8977770897207535d0*x25*x8)
      x28 = x22*x24 + x27
      x29 = 0.128278243853042d0*x10
      x30 = x28*x29 - 0.375d0
      x31 = x17*x30
      x32 = exp(-1.9241736577956332d0*x10)
      x33 = 0.00449266744437808d0*x2 + 1.0d0
      x34 = 1d0/x33
      x35 = 0.0204006002177468d0*x2 + 1.0d0
      x36 = x35**(-2)
      x37 = 0.0109702052306804d0*x2
      x38 = 1d0/gamma_aa
      x39 = tau_a*x38
      x40 = 8.0d0*rho_a*x39
      x41 = x40 - 1.0d0
      x42 = x2*x41
      x43 = 0.0123414808845154d0*x42 - 0.45d0
      x44 = 0.0274255130767009d0*x42 - 1.0d0
      x45 = x37*x41*x44 + 1.0d0
      x46 = x45**(-0.5d0)
      x47 = x37 + x43*x46
      x48 = x47**2
      x49 = 0.00197734563417202d0*x1
      x50 = gamma_aa*x49
      x51 = x46*(0.00222451383844352d0*x42 - 0.0811111111111111d0) + x50
     &
      x52 = rho_a**(-5.33333333333333d0)
      x53 = rho_a**(-2)
      x54 = tau_a**(-2)
      x55 = 0.0028125d0*x54
      x56 = 0.000135388578153653d0*x52 + x53*x55
      x57 = sqrt(x56)
      x58 = abs(gamma_aa)
      x59 = x57*x58
      x60 = x51*x59
      x61 = 1.24378109452736d0*x60
      x62 = gamma_aa**3
      x63 = rho_a**(-8.0d0)
      x64 = x62*x63
      x65 = gamma_aa**2
      x66 = x52*x65
      x67 = x53*x54
      x68 = x65*x67
      x69 = 0.015625d0*x68 + 1.0d0
      x70 = x69**(-2)
      x71 = 0.02485875d0*x68*x70 + 0.123456790123457d0
      x72 = x2*x71
      x73 = 1.86977850568312d-6*x64 + 6.38454405015215d-6*x66 + 0.002141
     &6494625995d0 *x68 + 0.0204668008035082d0*x72
      x74 = x36*(0.0896750813832074d0*x48 - x61 + x73) + 1.0d0
      x75 = 1d0/x74
      x76 = -x18 + x21
      x77 = -x24*x76 + x27
      x78 = x29*x77 - 0.375d0
      x79 = x17*x8
      x80 = 4.9628039271952d0*x79
      x81 = 0.804d0*x34
      x82 = -x40 + 1.0d0
      x83 = x2*x82
      x84 = 0.0123414808845154d0*x83 + 0.45d0
      x85 = 0.0274255130767009d0*x83 + 1.0d0
      x86 = x37*x82*x85 + 1.0d0
      x87 = x86**(-0.5d0)
      x88 = x37 - x84*x87
      x89 = x88**2
      x90 = x50 - x87*(0.00222451383844352d0*x83 + 0.0811111111111111d0)
     &
      x91 = x59*x90
      x92 = x36*(x73 + 0.0896750813832074d0*x89 - 1.24378109452736d0*x91
     &) + 1.0d0
      x93 = 1d0/x92
      x94 = -0.804d0*x93
      x95 = 1.8610514726982d0*x8
      x96 = rho_a**(-0.666666666666667d0)
      x97 = x31*x96
      x98 = x23*x25
      x99 = rho_a**(-0.333333333333333d0)*x98
      x100 = 10.128444160768d0*x99
      x101 = mu*x12
      x102 = 0.00844343197019482d0*rho_a**(-2.0d0)
      x103 = 0.0855188292353615d0*x12
      x104 = x102*x3 - x103
      x105 = 0.564189583547756d0
      x106 = x105*x98
      x107 = 0.128278243853042d0*x101 - x102*x26 + x104*x24 + 4.60575334
     &176456d0*x106 *x96
      x108 = 0.128278243853042d0*x9
      x109 = -x108*(x100*x76 + x107) + 0.0427594146176807d0*x12*x77
      x110 = mu*x80
      x111 = mu*x19
      x112 = 1.19366207318922d0*x111
      x113 = -0.804d0*x75
      x114 = x32*(x113 + x81)
      x115 = rho_a**(-4.33333333333333d0)
      x116 = x115*x3
      x117 = rho_a**(-3.66666666666667d0)
      x118 = x117*x6
      x119 = x116*x6
      x120 = rho_a**(-1.66666666666667d0)
      x121 = mu*x120
      x122 = rho_a**(-2.33333333333333d0)
      x123 = x122*x26
      x124 = 16.3704432459399d0*x12 - 7.10869298948088d0*x121 + 0.311935
     &284066016d0* x123
      x125 = x124*x14
      x126 = mu*x1
      x127 = x126*x7
      x128 = -0.0410630556449071d0*x116 - 0.138634955389984d0*x118 + 0.0
     &108360841285172d0*x119 + x125*x127
      x129 = 3.99009435746494d0*x8
      x130 = x128*x129
      x131 = x15/x16**2
      x132 = x131*x30
      x133 = gamma_aa*x132
      x134 = x33**(-2)
      x135 = rho_a**(-3.33333333333333d0)
      x136 = x134*x135
      x137 = -x81 + 1.804d0
      x138 = 0.6203504908994d0*x96
      x139 = x74**(-2)
      x140 = rho_a**(-6.33333333333333d0)
      x141 = x140*x65
      x142 = rho_a**(-9.0d0)
      x143 = x142*x62
      x144 = rho_a**(-3)
      x145 = x144*x54
      x146 = x145*x65
      x147 = gamma_aa*x117
      x148 = 0.0545781354760218d0*x71
      x149 = x68/x69
      x150 = 0.001553671875d0*x149
      x151 = x150 - 0.0497175d0
      x152 = x151*x62
      x153 = rho_a**(-5.66666666666667d0)
      x154 = x54*x70
      x155 = x153*x154
      x156 = 0.0204668008035082d0*x155
      x157 = 1d0/x57
      x158 = 1d0/x58
      x159 = x158*x51
      x160 = x157*x159
      x161 = 0.000361036208409741d0*x140 + x144*x55
      x162 = x161*x65
      x163 = x160*x162
      x164 = 1d0/x35
      x165 = 1.50330191856923d-6*x64 + 5.13317341632233d-6*x66 + 0.00172
     &188616792999d0*x68 + 0.0164553078460206d0*x72
      x166 = x165 + 0.0720987654320988d0*x48 - x60
      x167 = x164*x166
      x168 = 0.135327364628503d0*x167
      x169 = 0.00527292169112538d0*x117
      x170 = tau_a*x1
      x171 = x147*x41
      x172 = x45**(-1.5d0)
      x173 = 0.0438808209227215d0*x170
      x174 = 0.0146269403075738d0*x44
      x175 = 0.219404104613608d0*x170 - 0.0731347015378692d0*x171
      x176 = 0.00548510261534019d0*x42
      x177 = x171*x174 - x173*x44 - x175*x176
      x178 = x172*x177
      x179 = x178*x43
      x180 = -gamma_aa*x169 + 0.180246913580247d0*x179 + x46*(0.01779611
     &07075482d0* x170 - 0.00593203690251605d0*x171)
      x181 = x180*x59
      x182 = 0.0585077612302953d0*x117
      x183 = 0.0987318470761234d0*x170 - 0.0329106156920411d0*x171
      x184 = 2.0d0*x46
      x185 = 2.0d0*x43
      x186 = -gamma_aa*x182 + x178*x185 + x183*x184
      x187 = x186*x47
      x188 = -3.40509016008115d-5*x141 - 1.49582280454649d-5*x143 - 0.00
     &428329892519899d0*x146 - x147*x148 + x147*x168 + x152*x156 + 1.243
     &78109452736d0*x163 - 1.24378109452736d0*x181 + 0.0896750813832074d
     &0*x187
      x189 = 1.49628538404935d0*x8
      x190 = x189*x36
      x191 = x32*(0.0179261670197794d0*gamma_aa*x136 - x137*x138 + x138*
     &(x113 + 1.804d0) + x139*x188*x190)
      x192 = x122*x6
      x193 = 0.207437457468366d0*x192
      x194 = x122*x134
      x195 = 0.00336115631620864d0*x194
      x196 = 2.0d0*x84
      x197 = 0.0438808209227215d0*x120
      x198 = x197*x44
      x199 = x198*x39
      x200 = x1*x41
      x201 = 0.00548510261534019d0*x200
      x202 = x201*x44
      x203 = x120*x39
      x204 = 0.219404104613608d0*x203
      x205 = 0.0274255130767009d0*x200 - x204
      x206 = x176*x205
      x207 = -x199 + x202 + x206
      x208 = x86**(-1.5d0)
      x209 = x207*x208
      x210 = 0.0219404104613608d0*x1
      x211 = 0.0987318470761234d0*x203
      x212 = x1*x82
      x213 = x211 + 0.0123414808845154d0*x212
      x214 = 2.0d0*x213
      x215 = x210 - x214*x87
      x216 = x196*x209 + x215
      x217 = x216*x88
      x218 = 0.180246913580247d0*x84
      x219 = 0.0177961107075482d0*x120
      x220 = x219*x39
      x221 = x49 - x87*(0.00222451383844352d0*x212 + x220)
      x222 = x209*x218 + x221
      x223 = 1.24378109452736d0*x59
      x224 = x63*x65
      x225 = x1*x71
      x226 = 1.27690881003043d-5*x52
      x227 = x165 + 0.0720987654320988d0*x89 - x91
      x228 = x1*x164
      x229 = 0.0507477617356885d0*x228
      x230 = 0.00428329892519899d0*x67
      x231 = x158*x90
      x232 = 1.24378109452736d0*gamma_aa
      x233 = x232*x57
      x234 = 0.0204668008035082d0*x154
      x235 = rho_a**(-4.66666666666667d0)
      x236 = x151*x65
      x237 = x235*x236
      x238 = -gamma_aa*x226 - gamma_aa*x230 - 5.60933551704935d-6*x224 -
     & 0.0204668008035082d0*x225 + x227*x229 + x231*x233 + x234*x237
      x239 = -0.0896750813832074d0*x217 + x222*x223 + x238
      x240 = x92**(-2)
      x241 = x240*x36
      x242 = x239*x241
      x243 = 0.748142692024676d0*x8
      x244 = 2.0d0*x32
      x245 = 0.197463694152247d0*x120
      x246 = 1d0/x86
      x247 = 0.00120345402803247d0*x115
      x248 = x197*x85
      x249 = x246*(gamma_aa*x247*x82 + x248)
      x250 = 0.0720987654320988d0*x87*x88
      x251 = x59*x87
      x252 = tau_a**(-3)
      x253 = x252*x65
      x254 = x253*x53
      x255 = 0.0028125d0*x254
      x256 = x157*x231
      x257 = 0.0164553078460206d0*x152
      x258 = x252*x70
      x259 = x235*x257*x258 - 0.00344377233585999d0*x254
      x260 = x255*x256 + x259
      x261 = x250*(-x196*x249 + x245) - x251*(-x218*x249 + x219) + x260
      x262 = 1.8610514726982d0*x261
      x263 = x32*x36
      x264 = x0*x240*x263
      x265 = mu*x109
      x266 = x128*x133
      x267 = 7.98018871492988d0*x8
      x268 = rho_a**(-1.66666666666667d0)
      x269 = gamma_aa*x128
      x270 = 0.277269910779967d0*x118
      x271 = x13**(-2)
      x272 = x271*exp(-0.62530169814878134d0*x5)/x16**3
      x273 = x272*x30
      x274 = mu*x122
      x275 = 0.0168868639403896d0*rho_a**(-3.0d0)
      x276 = x23/x26
      x277 = x3*x52
      x278 = rho_a**(-6.0d0)
      x279 = x11*x278
      x280 = x115*x125*x26
      x281 = gamma_aa*x115
      x282 = x33**(-3)
      x283 = 0.413566993932933d0*x268
      x284 = x147*x71
      x285 = 0.270654729257005d0*x147
      x286 = 1/(x35**4*x92**3)
      x287 = x188*x189*x286
      x288 = rho_a**(-7.33333333333333d0)
      x289 = x288*x65
      x290 = rho_a**(-4)
      x291 = x290*x54
      x292 = gamma_aa*x235
      x293 = rho_a**(-6.66666666666667d0)
      x294 = tau_a**(-4)
      x295 = x294*x70
      x296 = 0.00014565673828125d0*gamma_aa**4*x290*x295
      x297 = x62*(-0.010875703125d0*x149 + x296 + 0.1491525d0)
      x298 = x56**(-1.5d0)*x58
      x299 = 1.24378109452736d0*x51
      x300 = x166*x36
      x301 = x158*x180
      x302 = x157*x162
      x303 = 1.24378109452736d0*x186
      x304 = tau_a*x117
      x305 = x292*x41
      x306 = x178*x183
      x307 = x172*(-0.087761641845443d0*x170*x175 + 0.0292538806151477d0
     &*x171*x175 - x176*(-1.17015522460591d0*x304 + 0.26816057230552d0*x
     &305) + 0.234031044921181d0*x304*x44 - 0.0536321144611041d0*x305*x4
     &4)
      x308 = 0.180246913580247d0*x43
      x309 = 0.131642462768165d0*x44
      x310 = 0.0164553078460206d0*x42
      x311 = x177*x45**(-2.5d0)
      x312 = x311*(-x170*x309 + 0.0438808209227215d0*x171*x44 - x175*x31
     &0)
      x313 = 0.000790938253668807d0*x1
      x314 = gamma_aa*x313
      x315 = x314 + x46*(0.000889805535377408d0*x42 - 0.0324444444444444
     &d0)
      x316 = 1.24378109452736d0*x315
      x317 = -2.73769248870524d-5*x141 - 1.20264153485538d-5*x143 - 0.00
     &344377233585999d0*x146 + x155*x257 + x163 - x181 + 0.0720987654320
     &988d0*x187 - 0.0438808209227215d0*x284
      x318 = x190*x240
      x319 = x131*x78
      x320 = x208*(0.00548510261534019d0*x212*x85 + x248*x39 + 0.0054851
     &0261534019d0* x83*(x204 + 0.0274255130767009d0*x212))
      x321 = x223*(x218*x320 + x221) + x238 - 0.0896750813832074d0*x88*(
     &x196*x320 + x215)
      x322 = x241*x243
      x323 = gamma_aa*x52
      x324 = 0.00856659785039798d0*gamma_aa
      x325 = 0.0409336016070163d0*x154
      x326 = gamma_aa*x57
      x327 = 2.48756218905473d0*x326
      x328 = 0.101495523471377d0*x228
      x329 = -0.179350162766415d0*x217 + 2.48756218905473d0*x222*x59 - 1
     &.12186710340987d-5*x224 - 0.0409336016070163d0*x225 + x227*x328 + 
     &x231*x327 + x237*x325 - 2.55381762006086d-5*x323 - x324*x67
      x330 = gamma_aa*x140
      x331 = -x150 + 0.0497175d0
      x332 = x331*x65
      x333 = 0.00932203125d0*x149
      x334 = x157*x38*x58
      x335 = x199 - x202 - x206
      x336 = x172*x335
      x337 = x308*x336 + x46*(0.00222451383844352d0*x200 - x220) + x49
      x338 = x158*x337
      x339 = 1.24378109452736d0*x302
      x340 = 0.00711844428301927d0*x120
      x341 = 0.0720987654320988d0*x43
      x342 = x313 + x336*x341 + x46*(0.000889805535377408d0*x200 - x340*
     &x39)
      x343 = x159*x57
      x344 = 0.0123414808845154d0*x200 - x211
      x345 = x184*x344 + x185*x336 + x210
      x346 = 0.0720987654320988d0*x47
      x347 = -gamma_aa*x343 + 0.00344377233585999d0*gamma_aa*x67 + 0.016
     &4553078460206d0*x154*x235*x332 + 4.50990575570768d-6*x224 + 0.0164
     &553078460206d0*x225 + 1.02663468326447d-5*x323 - x337*x59 + x345*x
     &346
      x348 = 0.135327364628503d0*x147*x164
      x349 = x170*x38
      x350 = x117*x41
      x351 = x178*x344
      x352 = x120*x309
      x353 = -0.0164553078460206d0*x200*x44 - x205*x310 + x352*x39
      x354 = x311*x353
      x355 = x175*x197
      x356 = x172*(0.0146269403075738d0*x171*x205 - x173*x205 + x174*x35
     &0 - x175*x201 - x176*(0.585077612302953d0*x349 - 0.073134701537869
     &2d0*x350) - 0.117015522460591d0*x349*x44 + x355*x39)
      x357 = 2.0d0*x183
      x358 = rho_a**(-5.0d0)
      x359 = x205*x43*(-0.0109702052306804d0*x200 + 0.087761641845443d0*
     &x203)
      x360 = x335*x344
      x361 = 1d0/x45
      x362 = x335*x353
      x363 = x361*x43
      x364 = x247*x41
      x365 = gamma_aa*x364 + x198
      x366 = x246*x365
      x367 = x196*x366 + x245
      x368 = 0.180246913580247d0*x365
      x369 = x219 + x246*x368*x84
      x370 = x250*x367 - x251*x369 + x260
      x371 = x361*x365
      x372 = x185*x371
      x373 = x363*x368
      x374 = x46*(x219 - x373)
      x375 = x160*x255 + x259 + x346*x46*(x245 - x372) - x374*x59
      x376 = x144*x253
      x377 = x153*x258
      x378 = x296 - x333 + 0.099435d0
      x379 = x252*x53
      x380 = x157*x254
      x381 = 0.00962763222425977d0*tau_a
      x382 = x361*(0.0731347015378692d0*x1*x44 - x115*x381 + 0.005214967
     &45480737d0* x323*x41 - x355)
      x383 = x177*x361
      x384 = 0.00361036208409741d0*x281*x41 + x352
      x385 = x177*x384/x45**2
      x386 = gamma_aa*x379
      x387 = 0.0153120874502694d0*x235*x253*x70
      x388 = 0.9305257363491d0*x87
      x389 = x158*x369
      x390 = x340 + 0.0720987654320988d0*x366*x84
      x391 = x84*(-x135*x38*x381 + x197*x205 + x364)
      x392 = x207*x246*x384
      x393 = 0.9305257363491d0*x208
      x394 = x314 - x87*(0.000889805535377408d0*x83 + 0.0324444444444444
     &d0)
      x395 = x36*x93
      x396 = x294*x53*x65
      x397 = x135*x43
      x398 = x120*x365
      x399 = 1.8610514726982d0*x172
      Ea = -x0*(-4.9628039271952d0*x31 + x32*(1.49628538404935d0*x34 - 1
     &.49628538404935d0*x75))
      d1Ea(1) = -rho_a*(x109*x110 + x112*x114 - x130*x133 + x191 - 1.654
     &26797573173d0* x97) - x32*x95*(x81 + x94) + x78*x80
      d1Ea(2) = rho_a*(x132*x193 + x244*(x195 + x242*x243))
      d1Ea(3) = -x262*x264
      d2Ea(1) = rho_a*(1.59154943091895d0*mu*rho_a**(-2.0d0)*x114 - 1.28
     &278243853042d0* x101*x191 + x110*(-x103*(-x100*x22 + x107) + x108*
     &( 3.37614805358933d0*rho_a**(-1.33333333333333d0)*x22*x98 - 20.256
     &888321536d0*x104*x99 - 46.6491155403329d0*x105*x19*x276 - 3.070502
     &22784304d0*x106*x268 + 102.585381117795d0*x22*x276*x96 + x24*(0.11
     &4025105647149d0*x122 - x275*x3) + x26*x275 - 0.171037658470723d0*x
     &274) + 0.0570125528235743d0*x122*x28) - 0.76560437251347d0*x114*x1
     &22*x3 + x129*x133*(-mu*x125*x270 + x1* x124*x271*x3*x7*(32.7408864
     &918799d0*x12 - 14.2173859789618d0*x121 + 0.623870568132033d0*x123)
     & + x127*x14*(-21.8272576612533d0*x122 + 11.8478216491348d0*x126 - 
     &0.727848996154038d0*x135*x26) + 0.508328169763273d0*x235*x6 - 0.07
     &58525888996201d0*x277*x6 + 0.28744138951435d0*x277 + 0.00225860726
     &894828d0*x279*x6 - 0.0171178656172922d0*x279 + 0.0216721682570343d
     &0*x280*x6 - 0.0821261112898142d0*x280) - x130*x273*x65*(-0.0821261
     &112898142d0 *x116 + 0.0216721682570343d0*x119 + 0.103976216542488d
     &0*x125*x126 *x6 - x270) - x131*x265*x267*x269 - 3.30853595146347d0
     &*x17*x265* x96 + 2.66006290497663d0*x266*x96 - 1.10284531715449d0*
     &x268*x31 + x32*(-0.00042952697051865d0*rho_a**(-7.0d0)*x282*x65 + 
     &0.0537785010593382d0*x134*x281 - x137*x283 - 0.997523589366235d0* 
     &x188*x241*x96 + x283*(x94 + 1.804d0) + x287*(-6.81018032016229d-5 
     &*x141 - 2.99164560909299d-5*x143 - 0.00856659785039798d0*x146 + 0.
     &0409336016070163d0*x152*x155 + 2.48756218905473d0*x163 + x167* x28
     &5 - 2.48756218905473d0*x181 + 0.179350162766415d0*x187 - 0.1091562
     &70952044d0*x284) - x318*(0.000134624052409184d0*rho_a**( -10.0d0)*
     &x62 - 0.109156270952044d0*x152*x154*x293 - 1.24378109452736d0*x160
     &*x65*(0.00228656265326169d0*x288 + 0.0084375d0*x291) + x161**2*x29
     &8*x299 + x164*x285*x317 - 0.496200336971176d0*x167*x292 - x223*(0.
     &0193340462007931d0*x292 + 0.360493827160494d0*x306 + x307*x308 + x
     &308*x312 + x46*( -0.0949125904402569d0*x304 + 0.0217508019758922d0
     &*x305)) + x234* x293*x297 + 0.0220860757144586d0*x289*x300 + 0.000
     &215655710138473d0*x289 + 0.012849896775597d0*x291*x65 + 0.20011983
     &0078746d0*x292*x71 + 2.48756218905473d0*x301*x302 + x303*(-0.00210
     &916867645015d0*x147 + 0.0720987654320988d0*x179 + x46*(0.007118444
     &28301927d0*x170 - 0.00237281476100642d0*x171)) + x316*(x184*(-0.52
     &6569851072658d0*x304 + 0.120672257537484d0*x305 ) + x185*x307 + x1
     &85*x312 + 0.214528457844416d0*x292 + 4.0d0*x306 )))) - 2.387324146
     &37843d0*x111*x114 - 2.0d0*x191 - 9.9256078543904d0*x265*x79 + x266
     &*x267 + 3.30853595146347d0*x97
      d2Ea(2) = rho_a*(0.641391219265211d0*x101*x32*(0.00672231263241728
     &d0*x194 + x318* x321) - 0.207437457468366d0*x109*x131*x274*x6 + x1
     &30*x319 + 0.0691458191561221d0*x135*x319*x6 - 0.414874914936733d0*
     &x192*x269 *x272*x78 - x32*(-0.000161072613944494d0*gamma_aa*x278*x
     &282 + 0.015685396142307d0*x136 - 0.498761794683118d0*x242*x96 + x2
     &87* x329 + x318*(-x117*x148 + x117*x168 - 4.48746841363948d-5*x142
     &* x65 - x145*x324 - 0.0545781354760218d0*x155*x332 + x156*x236 + x
     &156*x65*(-x296 + x333 - 0.099435d0) + x160*x232*( 0.00072207241681
     &9483d0*x140 + 0.005625d0*x145) - x161*x299*x334 - x223*(-x169 + 0.
     &180246913580247d0*x183*x336 + x308*x354 + x308* x356 + 0.180246913
     &580247d0*x351 + x46*(0.0474562952201284d0*x349 - 0.005932036902516
     &05d0*x350)) - x229*x317 - x233*x301 - 0.00828227839292197d0*x300*x
     &330 + x303*x342 + x316*(-x182 + x184* (0.263284925536329d0*x349 - 
     &0.0329106156920411d0*x350) + x185* x354 + x185*x356 + x336*x357 + 
     &2.0d0*x351) - 6.81018032016229d-5* x330 + x338*x339 + x347*x348)))
     & + x193*x319 + x244*(x195 + x321* x322)
      d2Ea(3) = -rho_a*(-x244*(x239*x243*x286*x329 - 3.02011151145926d-5
     &*x282*x358 - x322*(1.12186710340987d-5*gamma_aa*x63 - x172*x223*(x
     &308*x361* x362 + 0.180246913580247d0*x359 + 0.360493827160494d0*x3
     &60) + 2.48756218905473d0*x172*x315*(x359 + 2.0d0*x360 + x362*x363)
     & + x226 + x230 + x234*x292*(-0.007768359375d0*x149 + x296 + 0.0497
     &175d0) + x292*x325*x331 + 0.00310585439734574d0*x300*x52 - x327*x3
     &38 - x328*x347 + 1.24378109452736d0*x342*x345 - 1.24378109452736d0
     &*x343 + x61/x65)) + 0.0215685619967539d0*x273* x358*x6**2)
      d2Ea(4) = -x263*(rho_a*x139*(x112*x375 + x138*x375 - 3.72210294539
     &64d0*x188*x36* x375*x75*x8 + x189*(-0.0545781354760218d0*x152*x377
     & + x158*x339* x374 - 0.00699626865671642d0*x160*x376 + 0.003498134
     &32835821d0* x161*x298*x379*x51 - x223*x46*(-0.0296601845125803d0*x
     &1 - x183* x361*x368 + x219*x383 + x308*x382 - x308*x385) + 0.00349
     &813432835821d0*x301*x380 + x303*x46*(x340 - x341*x371) + x316*x46*
     &(-0.329106156920411d0*x1 + x185*x382 - x185*x385 + x245* x383 - x3
     &57*x371) + x348*x375 + 0.00856659785039798d0*x376 + 0.020466800803
     &5082d0*x377*x378*x62)) + x240*x370*x95)
      d2Ea(5) = x0*x241*x244*(-x151*x387 - 0.00261710363348184d0*x158*x2
     &22*x380 - x216* x388*x390 + 0.0379665670791648d0*x228*x261 - 0.005
     &23420726696369d0*x256*x386 - x262*x321*x395 + x326*x388*x389 + 0.0
     &0261710363348184d0*x334*x379*x90 + x378*x387 + 0.00640903757728955
     &d0*x386 - x393*x394*(-x196*x392 - x207*x245 + x214*x365 + 2.0d0*x3
     &91) + x393*x59*(-x207*x219 + x213*x368 - x218 *x392 + 0.1802469135
     &80247d0*x391))
      d2Ea(6) = x264*(-0.0306241749005388d0*x235*x295*x297 - 1.861051472
     &6982d0*x246*x367 *x390 + 0.0157026218008911d0*x256*x396 + 4.629481
     &27536866d0*x370 **2*x395 - 0.0104684145339274d0*x380*x389*x87 + x3
     &94*x399*(-x372* x384 + 0.0385105288970391d0*x397 + 0.3949273883044
     &94d0*x398) - 0.0192271127318687d0*x396 - x399*x59*(-x373*x384 + 0.
     &0034707019870171d0*x397 + 0.0355922214150963d0*x398) - 1.472120793
     &83354d-5*x290*x298*x90/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399, x400, x401
     &, x402, x403, x404, x405, x406
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = gamma_aa*x1
      x3 = mu**2
      x4 = rho_a**0.666666666666667d0
      x5 = x3/x4
      x6 = exp(-0.31265084907439067d0*x5)
      x7 = 1d0/x3
      x8 = 0.844037013397333d0*x4
      x9 = mu**(-4)
      x10 = rho_a**2.0d0
      x11 = 14.4026298885989d0*x10*x9
      x12 = 0.213719543995407d0*x0*x7
      x13 = x11 + x12 - x8
      x14 = x13*x7
      x15 = x14*x6
      x16 = 0.0519881082712438d0*x15*x2
      x17 = -x16 + 1.0d0
      x18 = 1.804d0 - 0.804d0/x17
      x19 = rho_a**0.333333333333333d0
      x20 = 1d0/x19
      x21 = 0.256556487706084d0*x20
      x22 = rho_a**1.0d0
      x23 = 1d0/x22
      x24 = 0.00844343197019482d0*x23
      x25 = x24*x3
      x26 = x21 - x25
      x27 = exp(-15.192666241151986d0*x4*x7)
      x28 = mu*x27
      x29 = 1d0/mu
      x30 = mu**3
      x31 = mu*x20
      x32 = x24*x30 - 0.384834731559127d0*x31 + 1.77245385090552d0*erf( 
     &3.8977770897207535d0*x19*x29)
      x33 = x26*x28 + x32
      x34 = 0.128278243853042d0*x31
      x35 = x33*x34 - 0.375d0
      x36 = x18*x35
      x37 = exp(-1.9241736577956332d0*x31)
      x38 = 0.00449266744437808d0*x2 + 1.0d0
      x39 = 1d0/x38
      x40 = 0.0204006002177468d0*x2 + 1.0d0
      x41 = x40**(-2)
      x42 = 0.0109702052306804d0*x2
      x43 = 1d0/gamma_aa
      x44 = tau_a*x43
      x45 = 8.0d0*rho_a*x44
      x46 = x45 - 1.0d0
      x47 = x2*x46
      x48 = 0.0123414808845154d0*x47 - 0.45d0
      x49 = 0.0274255130767009d0*x47 - 1.0d0
      x50 = x42*x46*x49 + 1.0d0
      x51 = x50**(-0.5d0)
      x52 = x42 + x48*x51
      x53 = x52**2
      x54 = 0.00197734563417202d0*x1
      x55 = gamma_aa*x54
      x56 = x51*(0.00222451383844352d0*x47 - 0.0811111111111111d0) + x55
     &
      x57 = rho_a**(-5.33333333333333d0)
      x58 = rho_a**(-2)
      x59 = tau_a**(-2)
      x60 = 0.0028125d0*x59
      x61 = 0.000135388578153653d0*x57 + x58*x60
      x62 = sqrt(x61)
      x63 = abs(gamma_aa)
      x64 = x62*x63
      x65 = x56*x64
      x66 = 1.24378109452736d0*x65
      x67 = gamma_aa**3
      x68 = rho_a**(-8.0d0)
      x69 = x67*x68
      x70 = gamma_aa**2
      x71 = x57*x70
      x72 = x58*x59
      x73 = x70*x72
      x74 = 0.015625d0*x73 + 1.0d0
      x75 = x74**(-2)
      x76 = 0.02485875d0*x73*x75 + 0.123456790123457d0
      x77 = x2*x76
      x78 = 1.86977850568312d-6*x69 + 6.38454405015215d-6*x71 + 0.002141
     &6494625995d0 *x73 + 0.0204668008035082d0*x77
      x79 = x41*(0.0896750813832074d0*x53 - x66 + x78) + 1.0d0
      x80 = 1d0/x79
      x81 = x16 - 1.0d0
      x82 = -x21 + x25
      x83 = -x28*x82 + x32
      x84 = x34*x83 - 0.375d0
      x85 = 4.9628039271952d0*x19
      x86 = 0.804d0*x39
      x87 = -x45 + 1.0d0
      x88 = x2*x87
      x89 = 0.0123414808845154d0*x88 + 0.45d0
      x90 = 0.0274255130767009d0*x88 + 1.0d0
      x91 = x42*x87*x90 + 1.0d0
      x92 = x91**(-0.5d0)
      x93 = x42 - x89*x92
      x94 = x93**2
      x95 = x55 - x92*(0.00222451383844352d0*x88 + 0.0811111111111111d0)
     &
      x96 = x64*x95
      x97 = x41*(x78 + 0.0896750813832074d0*x94 - 1.24378109452736d0*x96
     &) + 1.0d0
      x98 = 1d0/x97
      x99 = -0.804d0*x98
      x100 = 1.8610514726982d0*x19
      x101 = rho_a**(-0.666666666666667d0)
      x102 = x101*x36
      x103 = x19*x6
      x104 = 3.99009435746494d0*x103
      x105 = x17**(-2)
      x106 = x105*x35
      x107 = rho_a**(-4.33333333333333d0)
      x108 = -x11 - x12 + x8
      x109 = x107*x108
      x110 = rho_a**(-3.66666666666667d0)
      x111 = x108*x7
      x112 = x110*x111
      x113 = rho_a**(-0.333333333333333d0)
      x114 = 0.562691342264888d0*x113
      x115 = 28.8052597771978d0*x9
      x116 = x115*x22
      x117 = 0.284959391993876d0*rho_a**0.333333333333333d0*x7
      x118 = x114 - x116 - x117
      x119 = x118*x7
      x120 = 0.0519881082712438d0*x1
      x121 = 0.0108360841285172d0*x109 - 0.138634955389984d0*x112 + x119
     &*x120
      x122 = gamma_aa*x121
      x123 = x106*x122
      x124 = 1d0/x0
      x125 = x27*x29
      x126 = x113*x125
      x127 = 10.128444160768d0*x126
      x128 = mu*x124
      x129 = 0.00844343197019482d0/x10
      x130 = 0.0855188292353615d0*x124
      x131 = x129*x3 - x130
      x132 = 0.564189583547756d0
      x133 = x125*x132
      x134 = 4.60575334176456d0*x101*x133 + 0.128278243853042d0*x128 - x
     &129*x30 + x131*x28
      x135 = 0.128278243853042d0*x20
      x136 = 0.0427594146176807d0*x124*x83 - x135*(x127*x82 + x134)
      x137 = mu*x136
      x138 = x18*x85
      x139 = mu*x23
      x140 = 1.19366207318922d0*x139
      x141 = -0.804d0*x80
      x142 = x37*(x141 + x86)
      x143 = x38**(-2)
      x144 = rho_a**(-3.33333333333333d0)
      x145 = x143*x144
      x146 = -x86 + 1.804d0
      x147 = 0.6203504908994d0*x101
      x148 = x79**(-2)
      x149 = rho_a**(-6.33333333333333d0)
      x150 = x149*x70
      x151 = rho_a**(-9.0d0)
      x152 = x151*x67
      x153 = rho_a**(-3)
      x154 = x153*x59
      x155 = x154*x70
      x156 = gamma_aa*x110
      x157 = 0.0545781354760218d0*x76
      x158 = x73/x74
      x159 = 0.001553671875d0*x158
      x160 = x159 - 0.0497175d0
      x161 = x160*x67
      x162 = rho_a**(-5.66666666666667d0)
      x163 = x59*x75
      x164 = x162*x163
      x165 = 0.0204668008035082d0*x164
      x166 = 1d0/x62
      x167 = 1d0/x63
      x168 = x167*x56
      x169 = x166*x168
      x170 = 0.000361036208409741d0*x149 + x153*x60
      x171 = x170*x70
      x172 = x169*x171
      x173 = 1d0/x40
      x174 = 1.50330191856923d-6*x69 + 5.13317341632233d-6*x71 + 0.00172
     &188616792999d0*x73 + 0.0164553078460206d0*x77
      x175 = x174 + 0.0720987654320988d0*x53 - x65
      x176 = x173*x175
      x177 = 0.135327364628503d0*x176
      x178 = 0.00527292169112538d0*x110
      x179 = tau_a*x1
      x180 = x156*x46
      x181 = x50**(-1.5d0)
      x182 = 0.0438808209227215d0*x179
      x183 = 0.0146269403075738d0*x49
      x184 = 0.219404104613608d0*x179 - 0.0731347015378692d0*x180
      x185 = 0.00548510261534019d0*x47
      x186 = x180*x183 - x182*x49 - x184*x185
      x187 = x181*x186
      x188 = x187*x48
      x189 = -gamma_aa*x178 + 0.180246913580247d0*x188 + x51*(0.01779611
     &07075482d0* x179 - 0.00593203690251605d0*x180)
      x190 = x189*x64
      x191 = 0.0585077612302953d0*x110
      x192 = 0.0987318470761234d0*x179 - 0.0329106156920411d0*x180
      x193 = 2.0d0*x51
      x194 = 2.0d0*x48
      x195 = -gamma_aa*x191 + x187*x194 + x192*x193
      x196 = x195*x52
      x197 = -3.40509016008115d-5*x150 - 1.49582280454649d-5*x152 - 0.00
     &428329892519899d0*x155 - x156*x157 + x156*x177 + x161*x165 + 1.243
     &78109452736d0*x172 - 1.24378109452736d0*x190 + 0.0896750813832074d
     &0*x196
      x198 = 1.49628538404935d0*x19
      x199 = x198*x41
      x200 = x37*(0.0179261670197794d0*gamma_aa*x145 - x146*x147 + x147*
     &(x141 + 1.804d0) + x148*x197*x199)
      x201 = rho_a**(-2.33333333333333d0)
      x202 = 0.207437457468366d0*x201
      x203 = x15*x202
      x204 = x143*x201
      x205 = 0.00336115631620864d0*x204
      x206 = 2.0d0*x89
      x207 = rho_a**(-1.66666666666667d0)
      x208 = 0.0438808209227215d0*x207
      x209 = x208*x49
      x210 = x209*x44
      x211 = x1*x46
      x212 = 0.00548510261534019d0*x211
      x213 = x212*x49
      x214 = x207*x44
      x215 = 0.219404104613608d0*x214
      x216 = 0.0274255130767009d0*x211 - x215
      x217 = x185*x216
      x218 = -x210 + x213 + x217
      x219 = x91**(-1.5d0)
      x220 = x218*x219
      x221 = 0.0219404104613608d0*x1
      x222 = 0.0987318470761234d0*x214
      x223 = x1*x87
      x224 = x222 + 0.0123414808845154d0*x223
      x225 = 2.0d0*x224
      x226 = x221 - x225*x92
      x227 = x206*x220 + x226
      x228 = x227*x93
      x229 = 0.180246913580247d0*x89
      x230 = 0.0177961107075482d0*x207
      x231 = x230*x44
      x232 = x54 - x92*(0.00222451383844352d0*x223 + x231)
      x233 = x220*x229 + x232
      x234 = 1.24378109452736d0*x64
      x235 = x68*x70
      x236 = x1*x76
      x237 = 1.27690881003043d-5*x57
      x238 = x174 + 0.0720987654320988d0*x94 - x96
      x239 = x1*x173
      x240 = 0.0507477617356885d0*x239
      x241 = 0.00428329892519899d0*x72
      x242 = x167*x95
      x243 = 1.24378109452736d0*gamma_aa
      x244 = x243*x62
      x245 = 0.0204668008035082d0*x163
      x246 = rho_a**(-4.66666666666667d0)
      x247 = x160*x70
      x248 = x246*x247
      x249 = -gamma_aa*x237 - gamma_aa*x241 - 5.60933551704935d-6*x235 -
     & 0.0204668008035082d0*x236 + x238*x240 + x242*x244 + x245*x248
      x250 = -0.0896750813832074d0*x228 + x233*x234 + x249
      x251 = x97**(-2)
      x252 = x251*x41
      x253 = x250*x252
      x254 = 0.748142692024676d0*x19
      x255 = 2.0d0*x37
      x256 = 0.197463694152247d0*x207
      x257 = 1d0/x91
      x258 = 0.00120345402803247d0*x107
      x259 = x208*x90
      x260 = x257*(gamma_aa*x258*x87 + x259)
      x261 = 0.0720987654320988d0*x92*x93
      x262 = x64*x92
      x263 = tau_a**(-3)
      x264 = x263*x70
      x265 = x264*x58
      x266 = 0.0028125d0*x265
      x267 = x166*x242
      x268 = 0.0164553078460206d0*x161
      x269 = x263*x75
      x270 = x246*x268*x269 - 0.00344377233585999d0*x265
      x271 = x266*x267 + x270
      x272 = x261*(-x206*x260 + x256) - x262*(-x229*x260 + x230) + x271
      x273 = 1.8610514726982d0*x272
      x274 = x37*x41
      x275 = x0*x251*x274
      x276 = 7.98018871492988d0*x103
      x277 = x137*x18
      x278 = rho_a**(-1.66666666666667d0)
      x279 = rho_a**(-1.33333333333333d0)
      x280 = x120*x7
      x281 = rho_a**(-6.0d0)
      x282 = exp(-0.62530169814878134d0*x5)
      x283 = x282*x35/x17**3
      x284 = 0.0168868639403896d0*rho_a**(-3.0d0)
      x285 = x27/x30
      x286 = gamma_aa*x107
      x287 = x38**(-3)
      x288 = 0.413566993932933d0*x278
      x289 = x156*x76
      x290 = 0.270654729257005d0*x156
      x291 = 1/(x40**4*x97**3)
      x292 = x197*x198*x291
      x293 = rho_a**(-7.33333333333333d0)
      x294 = x293*x70
      x295 = rho_a**(-4)
      x296 = x295*x59
      x297 = gamma_aa*x246
      x298 = rho_a**(-6.66666666666667d0)
      x299 = tau_a**(-4)
      x300 = x299*x75
      x301 = 0.00014565673828125d0*gamma_aa**4*x295*x300
      x302 = x67*(-0.010875703125d0*x158 + x301 + 0.1491525d0)
      x303 = x61**(-1.5d0)*x63
      x304 = 1.24378109452736d0*x56
      x305 = x175*x41
      x306 = x167*x189
      x307 = x166*x171
      x308 = 1.24378109452736d0*x195
      x309 = tau_a*x110
      x310 = x297*x46
      x311 = x187*x192
      x312 = x181*(-0.087761641845443d0*x179*x184 + 0.0292538806151477d0
     &*x180*x184 - x185*(-1.17015522460591d0*x309 + 0.26816057230552d0*x
     &310) + 0.234031044921181d0*x309*x49 - 0.0536321144611041d0*x310*x4
     &9)
      x313 = 0.180246913580247d0*x48
      x314 = 0.131642462768165d0*x49
      x315 = 0.0164553078460206d0*x47
      x316 = x186*x50**(-2.5d0)
      x317 = x316*(-x179*x314 + 0.0438808209227215d0*x180*x49 - x184*x31
     &5)
      x318 = 0.000790938253668807d0*x1
      x319 = gamma_aa*x318
      x320 = x319 + x51*(0.000889805535377408d0*x47 - 0.0324444444444444
     &d0)
      x321 = 1.24378109452736d0*x320
      x322 = -2.73769248870524d-5*x150 - 1.20264153485538d-5*x152 - 0.00
     &344377233585999d0*x155 + x164*x268 + x172 - x190 + 0.0720987654320
     &988d0*x196 - 0.0438808209227215d0*x289
      x323 = x199*x251
      x324 = x81**(-2)
      x325 = x324*x84
      x326 = x219*(0.00548510261534019d0*x223*x90 + x259*x44 + 0.0054851
     &0261534019d0* x88*(x215 + 0.0274255130767009d0*x223))
      x327 = x234*(x229*x326 + x232) + x249 - 0.0896750813832074d0*x93*(
     &x206*x326 + x226)
      x328 = x252*x254
      x329 = -0.0108360841285172d0*x107*x13 + 0.138634955389984d0*x110*x
     &14 - x280*( -x114 + x116 + x117)
      x330 = gamma_aa*x57
      x331 = 0.00856659785039798d0*gamma_aa
      x332 = 0.0409336016070163d0*x163
      x333 = gamma_aa*x62
      x334 = 2.48756218905473d0*x333
      x335 = 0.101495523471377d0*x239
      x336 = -0.179350162766415d0*x228 + 2.48756218905473d0*x233*x64 - 1
     &.12186710340987d-5*x235 - 0.0409336016070163d0*x236 + x238*x335 + 
     &x242*x334 + x248*x332 - 2.55381762006086d-5*x330 - x331*x72
      x337 = gamma_aa*x149
      x338 = -x159 + 0.0497175d0
      x339 = x338*x70
      x340 = 0.00932203125d0*x158
      x341 = x166*x43*x63
      x342 = x210 - x213 - x217
      x343 = x181*x342
      x344 = x313*x343 + x51*(0.00222451383844352d0*x211 - x231) + x54
      x345 = x167*x344
      x346 = 1.24378109452736d0*x307
      x347 = 0.00711844428301927d0*x207
      x348 = 0.0720987654320988d0*x48
      x349 = x318 + x343*x348 + x51*(0.000889805535377408d0*x211 - x347*
     &x44)
      x350 = x168*x62
      x351 = 0.0123414808845154d0*x211 - x222
      x352 = x193*x351 + x194*x343 + x221
      x353 = 0.0720987654320988d0*x52
      x354 = -gamma_aa*x350 + 0.00344377233585999d0*gamma_aa*x72 + 0.016
     &4553078460206d0*x163*x246*x339 + 4.50990575570768d-6*x235 + 0.0164
     &553078460206d0*x236 + 1.02663468326447d-5*x330 - x344*x64 + x352*x
     &353
      x355 = 0.135327364628503d0*x156*x173
      x356 = x179*x43
      x357 = x110*x46
      x358 = x187*x351
      x359 = x207*x314
      x360 = -0.0164553078460206d0*x211*x49 - x216*x315 + x359*x44
      x361 = x316*x360
      x362 = x184*x208
      x363 = x181*(0.0146269403075738d0*x180*x216 - x182*x216 + x183*x35
     &7 - x184*x212 - x185*(0.585077612302953d0*x356 - 0.073134701537869
     &2d0*x357) - 0.117015522460591d0*x356*x49 + x362*x44)
      x364 = 2.0d0*x192
      x365 = rho_a**(-5.0d0)
      x366 = x216*x48*(-0.0109702052306804d0*x211 + 0.087761641845443d0*
     &x214)
      x367 = x342*x351
      x368 = 1d0/x50
      x369 = x342*x360
      x370 = x368*x48
      x371 = x258*x46
      x372 = gamma_aa*x371 + x209
      x373 = x257*x372
      x374 = x206*x373 + x256
      x375 = 0.180246913580247d0*x372
      x376 = x230 + x257*x375*x89
      x377 = x261*x374 - x262*x376 + x271
      x378 = x368*x372
      x379 = x194*x378
      x380 = x370*x375
      x381 = x51*(x230 - x380)
      x382 = x169*x266 + x270 + x353*x51*(x256 - x379) - x381*x64
      x383 = x153*x264
      x384 = x162*x269
      x385 = x301 - x340 + 0.099435d0
      x386 = x263*x58
      x387 = x166*x265
      x388 = 0.00962763222425977d0*tau_a
      x389 = x368*(0.0731347015378692d0*x1*x49 - x107*x388 + 0.005214967
     &45480737d0* x330*x46 - x362)
      x390 = x186*x368
      x391 = 0.00361036208409741d0*x286*x46 + x359
      x392 = x186*x391/x50**2
      x393 = gamma_aa*x386
      x394 = 0.0153120874502694d0*x246*x264*x75
      x395 = 0.9305257363491d0*x92
      x396 = x167*x376
      x397 = x347 + 0.0720987654320988d0*x373*x89
      x398 = x89*(-x144*x388*x43 + x208*x216 + x371)
      x399 = x218*x257*x391
      x400 = 0.9305257363491d0*x219
      x401 = x319 - x92*(0.000889805535377408d0*x88 + 0.0324444444444444
     &d0)
      x402 = x41*x98
      x403 = x299*x58*x70
      x404 = x144*x48
      x405 = x207*x372
      x406 = 1.8610514726982d0*x181
      Ea = -x0*(-4.9628039271952d0*x36 + x37*(1.49628538404935d0*x39 - 1
     &.49628538404935d0*x80))
      d1Ea(1) = -rho_a*(-1.65426797573173d0*x102 - x104*x123 + x137*x138
     & + x140*x142 + x200) - x100*x37*(x86 + x99) + x84*x85*(1.804d0 + 0
     &.804d0/x81)
      d1Ea(2) = -rho_a*(x106*x203 - x255*(x205 + x253*x254))
      d1Ea(3) = -x273*x275
      d2Ea(1) = rho_a*(gamma_aa*x104*x106*(0.0216721682570343d0*x107*x11
     &8 + 0.00225860726894828d0*x108*x281*x3 - 0.0758525888996201d0*x108
     &* x57 - 0.277269910779967d0*x110*x119 + 0.508328169763273d0*x111* 
     &x246 - x280*(0.0949864639979585d0*x101*x7 + x115 + 0.1875637807549
     &63d0*x279)) + 1.59154943091895d0*mu*rho_a**(-2.0d0 )*x142 + mu*x13
     &8*(-x130*(-x127*x26 + x134) + x135*( -0.171037658470723d0*mu*x201 
     &+ 102.585381117795d0*x101*x26*x285 + 3.37614805358933d0*x125*x26*x
     &279 - 20.256888321536d0*x126*x131 - 46.6491155403329d0*x132*x23*x2
     &85 - 3.07050222784304d0*x133*x278 + x28*(0.114025105647149d0*x201 
     &- x284*x3) + x284*x30) + 0.0570125528235743d0*x201*x33) + 2.660062
     &90497663d0*x101*x123*x6 - 3.30853595146347d0*x101*x277 - x105*x122
     &*x137*x276 - 3.99009435746494d0*x121*x19*x283*x70*(0.1039762165424
     &88d0*x1*x119 + 0.0216721682570343d0*x109 - 0.277269910779967d0*x11
     &2) - 1.28278243853042d0*x128*x200 - 0.76560437251347d0*x142*x201*x
     &3 - 1.10284531715449d0*x278*x36 + x37*(-0.00042952697051865d0*rho_
     &a** (-7.0d0)*x287*x70 - 0.997523589366235d0*x101*x197*x252 + 0.053
     &7785010593382d0*x143*x286 - x146*x288 + x288*(x99 + 1.804d0 ) + x2
     &92*(-6.81018032016229d-5*x150 - 2.99164560909299d-5*x152 - 0.00856
     &659785039798d0*x155 + 0.0409336016070163d0*x161*x164 + 2.487562189
     &05473d0*x172 + x176*x290 - 2.48756218905473d0*x190 + 0.17935016276
     &6415d0*x196 - 0.109156270952044d0*x289) - x323*( 0.000134624052409
     &184d0*rho_a**(-10.0d0)*x67 - 0.109156270952044d0 *x161*x163*x298 -
     & 1.24378109452736d0*x169*x70*( 0.00228656265326169d0*x293 + 0.0084
     &375d0*x296) + x170**2*x303* x304 + x173*x290*x322 - 0.496200336971
     &176d0*x176*x297 - x234*( 0.0193340462007931d0*x297 + 0.36049382716
     &0494d0*x311 + x312*x313 + x313*x317 + x51*(-0.0949125904402569d0*x
     &309 + 0.0217508019758922d0*x310)) + x245*x298*x302 + 0.02208607571
     &44586d0*x294*x305 + 0.000215655710138473d0*x294 + 0.01284989677559
     &7d0*x296*x70 + 0.200119830078746d0*x297*x76 + 2.48756218905473d0*x
     &306*x307 + x308*(-0.00210916867645015d0*x156 + 0.0720987654320988d
     &0*x188 + x51*(0.00711844428301927d0*x179 - 0.00237281476100642d0*x
     &180)) + x321*(x193*(-0.526569851072658d0* x309 + 0.120672257537484
     &d0*x310) + x194*x312 + x194*x317 + 0.214528457844416d0*x297 + 4.0d
     &0*x311)))) + 3.30853595146347d0* x102 + x123*x276 - 2.387324146378
     &43d0*x139*x142 - 9.9256078543904d0*x19*x277 - 2.0d0*x200
      d2Ea(2) = rho_a*(-0.414874914936733d0*gamma_aa*x14*x201*x282*x329*
     &x84/x81**3 + x104*x325*x329 + 0.641391219265211d0*x128*x37*( 0.006
     &72231263241728d0*x204 + x323*x327) + x13*x136*x202*x29*x324* x6 - 
     &0.0691458191561221d0*x144*x15*x325 - x37*( -0.000161072613944494d0
     &*gamma_aa*x281*x287 - 0.498761794683118d0* x101*x253 + 0.015685396
     &142307d0*x145 + x292*x336 + x323*(-x110* x157 + x110*x177 - 4.4874
     &6841363948d-5*x151*x70 - x154*x331 - 0.0545781354760218d0*x164*x33
     &9 + x165*x247 + x165*x70*(-x301 + x340 - 0.099435d0) + x169*x243*(
     &0.000722072416819483d0*x149 + 0.005625d0*x154) - x170*x304*x341 - 
     &x234*(-x178 + 0.180246913580247d0*x192*x343 + x313*x361 + x313*x36
     &3 + 0.180246913580247d0*x358 + x51*(0.0474562952201284d0*x356 - 0.
     &00593203690251605d0*x357)) - x240*x322 - x244*x306 - 0.00828227839
     &292197d0*x305*x337 + x308*x349 + x321*(-x191 + x193* (0.2632849255
     &36329d0*x356 - 0.0329106156920411d0*x357) + x194* x361 + x194*x363
     & + x343*x364 + 2.0d0*x358) - 6.81018032016229d-5* x337 + x345*x346
     & + x354*x355))) - x203*x325 + x255*(x205 + x327* x328)
      d2Ea(3) = -rho_a*(0.0215685619967539d0*x13**2*x283*x365*x9 - x255*
     &(x250*x254*x291* x336 - 3.02011151145926d-5*x287*x365 - x328*(1.12
     &186710340987d-5* gamma_aa*x68 - x181*x234*(x313*x368*x369 + 0.1802
     &46913580247d0* x366 + 0.360493827160494d0*x367) + 2.48756218905473
     &d0*x181*x320*( x366 + 2.0d0*x367 + x369*x370) + x237 + x241 + x245
     &*x297*( -0.007768359375d0*x158 + x301 + 0.0497175d0) + x297*x332*x
     &338 + 0.00310585439734574d0*x305*x57 - x334*x345 - x335*x354 + 1.2
     &4378109452736d0*x349*x352 - 1.24378109452736d0*x350 + x66/x70 )))
      d2Ea(4) = -x274*(rho_a*x148*(x140*x382 + x147*x382 - 3.72210294539
     &64d0*x19*x197* x382*x41*x80 + x198*(-0.0545781354760218d0*x161*x38
     &4 + x167*x346* x381 - 0.00699626865671642d0*x169*x383 + 0.00349813
     &432835821d0* x170*x303*x386*x56 - x234*x51*(-0.0296601845125803d0*
     &x1 - x192* x368*x375 + x230*x390 + x313*x389 - x313*x392) + 0.0034
     &9813432835821d0*x306*x387 + x308*x51*(x347 - x348*x378) + x321*x51
     &*(-0.329106156920411d0*x1 + x194*x389 - x194*x392 + x256* x390 - x
     &364*x378) + x355*x382 + 0.00856659785039798d0*x383 + 0.02046680080
     &35082d0*x384*x385*x67)) + x100*x251*x377)
      d2Ea(5) = x0*x252*x255*(-x160*x394 - 0.00261710363348184d0*x167*x2
     &33*x387 - x227* x395*x397 + 0.0379665670791648d0*x239*x272 - 0.005
     &23420726696369d0*x267*x393 - x273*x327*x402 + x333*x395*x396 + 0.0
     &0261710363348184d0*x341*x386*x95 + x385*x394 + 0.00640903757728955
     &d0*x393 - x400*x401*(-x206*x399 - x218*x256 + x225*x372 + 2.0d0*x3
     &98) + x400*x64*(-x218*x230 + x224*x375 - x229 *x399 + 0.1802469135
     &80247d0*x398))
      d2Ea(6) = x275*(-0.0306241749005388d0*x246*x300*x302 - 1.861051472
     &6982d0*x257*x374 *x397 + 0.0157026218008911d0*x267*x403 + 4.629481
     &27536866d0*x377 **2*x402 - 0.0104684145339274d0*x387*x396*x92 + x4
     &01*x406*(-x379* x391 + 0.0385105288970391d0*x404 + 0.3949273883044
     &94d0*x405) - 0.0192271127318687d0*x403 - x406*x64*(-x380*x391 + 0.
     &0034707019870171d0*x404 + 0.0355922214150963d0*x405) - 1.472120793
     &83354d-5*x295*x303*x95/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358, x359, x360, x361, x362, x363, x364, x365, x366, x367, x368
     &, x369, x370, x371, x372, x373, x374, x375, x376, x377, x378, x379
     &, x380, x381, x382, x383, x384, x385, x386, x387, x388, x389, x390
     &, x391, x392, x393, x394, x395, x396, x397, x398, x399, x400, x401
     &, x402, x403, x404, x405, x406, x407, x408, x409, x410, x411, x412
     &, x413, x414, x415, x416, x417, x418, x419, x420, x421, x422, x423
     &, x424, x425, x426, x427, x428, x429, x430, x431, x432, x433, x434
     &, x435, x436, x437, x438, x439, x440
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 1d0/x0
      x2 = mu**4
      x3 = mu**2
      x4 = 1d0/x3
      x5 = rho_a**0.666666666666667d0
      x6 = 15.192666241152d0*x4*x5
      x7 = exp(x6)
      x8 = x2*(x7 - 1.0d0)
      x9 = x1*x8
      x10 = 0.000722072416819482d0*x9
      x11 = 7.0d0*x7 + 11.0d0
      x12 = x3/x5
      x13 = 0.00060945584614891d0*x11*x12
      x14 = -x10 + x13 + 0.0185185185185185d0
      x15 = rho_a**(-2.0d0)
      x16 = 3.15935190989552d0*x15
      x17 = x16*x4
      x18 = x14*x17
      x19 = 3.0d0*x7
      x20 = x19 - 2.0d0
      x21 = rho_a**0.333333333333333d0
      x22 = 1d0/x21
      x23 = mu*x22
      x24 = 1d0/mu
      x25 = erf(3.8977770897207535d0*x21*x24)
      x26 = x25*x7
      x27 = -0.131642462768164d0*x12*x20 - x19 + 1.81893813843777d0*x23*
     &x26 + 0.00866486900183379d0*x9
      x28 = 1d0/x27
      x29 = x28*exp(-0.31265084907439067d0*x12)
      x30 = gamma_aa*x29
      x31 = x18*x30
      x32 = -x31 + 1.0d0
      x33 = 1.804d0 - 0.804d0/x32
      x34 = 0.256556487706084d0*x22
      x35 = 1d0/rho_a
      x36 = 0.00844343197019482d0*x35
      x37 = x3*x36
      x38 = x34 - x37
      x39 = exp(-x6)
      x40 = mu*x39
      x41 = mu**3
      x42 = -0.384834731559127d0*x23 + 1.77245385090552d0*x25 + x36*x41
      x43 = x38*x40 + x42
      x44 = 0.128278243853042d0*x23
      x45 = x43*x44 - 0.375d0
      x46 = x33*x45
      x47 = exp(-1.9241736577956332d0*x23)
      x48 = rho_a**(-2.66666666666667d0)
      x49 = gamma_aa*x48
      x50 = 0.00449266744437808d0*x49 + 1.0d0
      x51 = 1d0/x50
      x52 = 0.0204006002177468d0*x49 + 1.0d0
      x53 = x52**(-2)
      x54 = 0.0109702052306804d0*x49
      x55 = 1d0/gamma_aa
      x56 = tau_a*x55
      x57 = 8.0d0*rho_a*x56
      x58 = x57 - 1.0d0
      x59 = x49*x58
      x60 = 0.0123414808845154d0*x59 - 0.45d0
      x61 = 0.0274255130767009d0*x59 - 1.0d0
      x62 = x54*x58*x61 + 1.0d0
      x63 = x62**(-0.5d0)
      x64 = x54 + x60*x63
      x65 = x64**2
      x66 = 0.00197734563417202d0*x48
      x67 = gamma_aa*x66
      x68 = x63*(0.00222451383844352d0*x59 - 0.0811111111111111d0) + x67
     &
      x69 = rho_a**(-5.33333333333333d0)
      x70 = rho_a**(-2)
      x71 = tau_a**(-2)
      x72 = 0.0028125d0*x71
      x73 = 0.000135388578153653d0*x69 + x70*x72
      x74 = sqrt(x73)
      x75 = abs(gamma_aa)
      x76 = x74*x75
      x77 = x68*x76
      x78 = 1.24378109452736d0*x77
      x79 = gamma_aa**3
      x80 = rho_a**(-8.0d0)
      x81 = x79*x80
      x82 = gamma_aa**2
      x83 = x69*x82
      x84 = x70*x71
      x85 = x82*x84
      x86 = 0.015625d0*x85 + 1.0d0
      x87 = x86**(-2)
      x88 = 0.02485875d0*x85*x87 + 0.123456790123457d0
      x89 = x49*x88
      x90 = 1.86977850568312d-6*x81 + 6.38454405015215d-6*x83 + 0.002141
     &6494625995d0 *x85 + 0.0204668008035082d0*x89
      x91 = x53*(0.0896750813832074d0*x65 - x78 + x90) + 1.0d0
      x92 = 1d0/x91
      x93 = x31 - 1.0d0
      x94 = -x34 + x37
      x95 = -x40*x94 + x42
      x96 = x44*x95 - 0.375d0
      x97 = 4.9628039271952d0*x21
      x98 = 0.804d0*x51
      x99 = -x57 + 1.0d0
      x100 = x49*x99
      x101 = 0.0123414808845154d0*x100 + 0.45d0
      x102 = 0.0274255130767009d0*x100 + 1.0d0
      x103 = x102*x54*x99 + 1.0d0
      x104 = x103**(-0.5d0)
      x105 = -x101*x104 + x54
      x106 = x105**2
      x107 = -x104*(0.00222451383844352d0*x100 + 0.0811111111111111d0) +
     & x67
      x108 = x107*x76
      x109 = x53*(0.0896750813832074d0*x106 - 1.24378109452736d0*x108 + 
     &x90) + 1.0d0
      x110 = 1d0/x109
      x111 = -0.804d0*x110
      x112 = 1.8610514726982d0*x21
      x113 = rho_a**(-0.666666666666667d0)
      x114 = x113*x46
      x115 = mu*x35
      x116 = 1.19366207318922d0*x115
      x117 = -0.804d0*x92
      x118 = x47*(x117 + x98)
      x119 = rho_a**(-0.333333333333333d0)
      x120 = x24*x39
      x121 = x119*x120
      x122 = 10.128444160768d0*x121
      x123 = mu*x1
      x124 = 0.00844343197019482d0*x15
      x125 = 0.0855188292353615d0*x1
      x126 = x124*x3 - x125
      x127 = 0.564189583547756d0
      x128 = x120*x127
      x129 = 4.60575334176456d0*x113*x128 + 0.128278243853042d0*x123 - x
     &124*x41 + x126*x40
      x130 = 0.128278243853042d0*x22
      x131 = 0.0427594146176807d0*x1*x95 - x130*(x122*x94 + x129)
      x132 = mu*x131
      x133 = x33*x97
      x134 = rho_a**(-3.66666666666667d0)
      x135 = x10 - x13 - 0.0185185185185185d0
      x136 = x134*x135
      x137 = rho_a**(-3.0d0)
      x138 = x135*x4
      x139 = x137*x138
      x140 = x35*x7
      x141 = 0.0432098765432099d0*x140
      x142 = rho_a**(-1.66666666666667d0)
      x143 = x142*x3
      x144 = x143*x7
      x145 = 0.00731347015378692d0*x144
      x146 = rho_a**(-2.33333333333333d0)
      x147 = x146*x8
      x148 = 0.000962763222425976d0*x147
      x149 = 0.000406303897432606d0*x11*x143
      x150 = -x141 + x145 - x148 + x149
      x151 = x150*x4
      x152 = x127*x35
      x153 = x4*x7
      x154 = 30.385332482304d0*x153
      x155 = x143*x20
      x156 = x123*x26
      x157 = x113*x7
      x158 = x157*x24*x25
      x159 = x119*x154 + 4.0d0*x140 - 0.087761641845443d0*x144 + 0.01155
     &31586691117d0 *x147 - 4.72654360241471d0*x152 - 0.087761641845443d
     &0*x155 + 0.606312712812591d0*x156 - 18.4230133670582d0*x158
      x160 = x159*x28
      x161 = x138*x16
      x162 = 0.658516038102422d0*x136 - 6.31870381979104d0*x139 + x151*x
     &16 + x160* x161
      x163 = 3.99009435746494d0*x21
      x164 = x162*x163*x45
      x165 = x30/x32**2
      x166 = x50**(-2)
      x167 = rho_a**(-3.33333333333333d0)
      x168 = x166*x167
      x169 = -x98 + 1.804d0
      x170 = 0.6203504908994d0*x113
      x171 = x91**(-2)
      x172 = rho_a**(-6.33333333333333d0)
      x173 = x172*x82
      x174 = rho_a**(-9.0d0)
      x175 = x174*x79
      x176 = rho_a**(-3)
      x177 = x176*x71
      x178 = x177*x82
      x179 = gamma_aa*x134
      x180 = 0.0545781354760218d0*x88
      x181 = x85/x86
      x182 = 0.001553671875d0*x181
      x183 = x182 - 0.0497175d0
      x184 = x183*x79
      x185 = rho_a**(-5.66666666666667d0)
      x186 = x71*x87
      x187 = x185*x186
      x188 = 0.0204668008035082d0*x187
      x189 = 1d0/x74
      x190 = 1d0/x75
      x191 = x190*x68
      x192 = x189*x191
      x193 = 0.000361036208409741d0*x172 + x176*x72
      x194 = x193*x82
      x195 = x192*x194
      x196 = 1d0/x52
      x197 = 1.50330191856923d-6*x81 + 5.13317341632233d-6*x83 + 0.00172
     &188616792999d0*x85 + 0.0164553078460206d0*x89
      x198 = x197 + 0.0720987654320988d0*x65 - x77
      x199 = x196*x198
      x200 = 0.135327364628503d0*x199
      x201 = 0.00527292169112538d0*x134
      x202 = tau_a*x48
      x203 = x179*x58
      x204 = x62**(-1.5d0)
      x205 = 0.0438808209227215d0*x202
      x206 = 0.0146269403075738d0*x61
      x207 = 0.219404104613608d0*x202 - 0.0731347015378692d0*x203
      x208 = 0.00548510261534019d0*x59
      x209 = x203*x206 - x205*x61 - x207*x208
      x210 = x204*x209
      x211 = x210*x60
      x212 = -gamma_aa*x201 + 0.180246913580247d0*x211 + x63*(0.01779611
     &07075482d0* x202 - 0.00593203690251605d0*x203)
      x213 = x212*x76
      x214 = 0.0585077612302953d0*x134
      x215 = 0.0987318470761234d0*x202 - 0.0329106156920411d0*x203
      x216 = 2.0d0*x63
      x217 = 2.0d0*x60
      x218 = -gamma_aa*x214 + x210*x217 + x215*x216
      x219 = x218*x64
      x220 = -3.40509016008115d-5*x173 - 1.49582280454649d-5*x175 - 0.00
     &428329892519899d0*x178 - x179*x180 + x179*x200 + x184*x188 + 1.243
     &78109452736d0*x195 - 1.24378109452736d0*x213 + 0.0896750813832074d
     &0*x219
      x221 = 1.49628538404935d0*x21
      x222 = x221*x53
      x223 = x47*(0.0179261670197794d0*gamma_aa*x168 - x169*x170 + x170*
     &(x117 + 1.804d0) + x171*x220*x222)
      x224 = rho_a**(-1.66666666666667d0)
      x225 = x14*x4
      x226 = x224*x225
      x227 = x29/x93**2
      x228 = 12.6061122289202d0*x227
      x229 = x146*x166
      x230 = 2.0d0*x101
      x231 = x103**(-1.5d0)
      x232 = 0.0438808209227215d0*x142
      x233 = x102*x232
      x234 = x48*x99
      x235 = x142*x56
      x236 = 0.219404104613608d0*x235
      x237 = x231*(0.00548510261534019d0*x100*(0.0274255130767009d0*x234
     & + x236) + 0.00548510261534019d0*x102*x234 + x233*x56)
      x238 = 0.0219404104613608d0*x48
      x239 = 0.0987318470761234d0*x235
      x240 = 0.0123414808845154d0*x234 + x239
      x241 = 2.0d0*x240
      x242 = -x104*x241 + x238
      x243 = 0.0896750813832074d0*x105
      x244 = 0.180246913580247d0*x101
      x245 = 0.0177961107075482d0*x142
      x246 = x245*x56
      x247 = -x104*(0.00222451383844352d0*x234 + x246) + x66
      x248 = 1.24378109452736d0*x76
      x249 = x80*x82
      x250 = x48*x88
      x251 = 1.27690881003043d-5*x69
      x252 = 0.0720987654320988d0*x106 - x108 + x197
      x253 = x196*x48
      x254 = 0.0507477617356885d0*x253
      x255 = 0.00428329892519899d0*x84
      x256 = x107*x190
      x257 = 1.24378109452736d0*gamma_aa
      x258 = x257*x74
      x259 = 0.0204668008035082d0*x186
      x260 = rho_a**(-4.66666666666667d0)
      x261 = x183*x82
      x262 = x260*x261
      x263 = -gamma_aa*x251 - gamma_aa*x255 - 5.60933551704935d-6*x249 -
     & 0.0204668008035082d0*x250 + x252*x254 + x256*x258 + x259*x262
      x264 = -x243*(x230*x237 + x242) + x248*(x237*x244 + x247) + x263
      x265 = x109**(-2)
      x266 = x265*x53
      x267 = 0.748142692024676d0*x21
      x268 = x266*x267
      x269 = 2.0d0*x47
      x270 = -x226*x228*x96 + x269*(0.00336115631620864d0*x229 + x264*x2
     &68)
      x271 = 0.197463694152247d0*x142
      x272 = 1d0/x103
      x273 = rho_a**(-4.33333333333333d0)
      x274 = 0.00120345402803247d0*x273
      x275 = x272*(gamma_aa*x274*x99 + x233)
      x276 = 0.0720987654320988d0*x104*x105
      x277 = x104*x76
      x278 = tau_a**(-3)
      x279 = x278*x82
      x280 = x279*x70
      x281 = 0.0028125d0*x280
      x282 = x189*x256
      x283 = 0.0164553078460206d0*x184
      x284 = x278*x87
      x285 = x260*x283*x284 - 0.00344377233585999d0*x280
      x286 = x281*x282 + x285
      x287 = x276*(-x230*x275 + x271) - x277*(-x244*x275 + x245) + x286
      x288 = 1.8610514726982d0*x287
      x289 = x47*x53
      x290 = x0*x265*x289
      x291 = x132*x33
      x292 = x165*x45
      x293 = x162*x292
      x294 = 7.98018871492988d0*x21
      x295 = mu*x146
      x296 = 0.0168868639403896d0*x137
      x297 = 1d0/x41
      x298 = x297*x39
      x299 = rho_a**(-1.33333333333333d0)
      x300 = 1.31703207620484d0*x136
      x301 = 12.6374076395821d0*x139
      x302 = 6.31870381979104d0*x15
      x303 = x151*x302
      x304 = x27**(-2)
      x305 = x304*exp(-0.62530169814878134d0*x12)
      x306 = x15*x7
      x307 = x3*x48
      x308 = x307*x7
      x309 = x167*x8
      x310 = 1d0/x2
      x311 = x24*x26
      x312 = gamma_aa*x273
      x313 = x50**(-3)
      x314 = 0.413566993932933d0*x224
      x315 = x113*x266
      x316 = x179*x88
      x317 = 0.270654729257005d0*x179
      x318 = 1/(x109**3*x52**4)
      x319 = x220*x221*x318
      x320 = rho_a**(-7.33333333333333d0)
      x321 = x320*x82
      x322 = rho_a**(-4)
      x323 = x322*x71
      x324 = gamma_aa*x260
      x325 = rho_a**(-6.66666666666667d0)
      x326 = tau_a**(-4)
      x327 = x326*x87
      x328 = 0.00014565673828125d0*gamma_aa**4*x322*x327
      x329 = x79*(-0.010875703125d0*x181 + x328 + 0.1491525d0)
      x330 = x73**(-1.5d0)*x75
      x331 = 1.24378109452736d0*x68
      x332 = x198*x53
      x333 = x190*x212
      x334 = x189*x194
      x335 = 1.24378109452736d0*x218
      x336 = tau_a*x134
      x337 = x324*x58
      x338 = x210*x215
      x339 = x204*(-0.087761641845443d0*x202*x207 + 0.0292538806151477d0
     &*x203*x207 - x208*(-1.17015522460591d0*x336 + 0.26816057230552d0*x
     &337) + 0.234031044921181d0*x336*x61 - 0.0536321144611041d0*x337*x6
     &1)
      x340 = 0.180246913580247d0*x60
      x341 = 0.131642462768165d0*x61
      x342 = 0.0164553078460206d0*x59
      x343 = x209*x62**(-2.5d0)
      x344 = x343*(-x202*x341 + 0.0438808209227215d0*x203*x61 - x207*x34
     &2)
      x345 = 0.000790938253668807d0*x48
      x346 = gamma_aa*x345
      x347 = x346 + x63*(0.000889805535377408d0*x59 - 0.0324444444444444
     &d0)
      x348 = 1.24378109452736d0*x347
      x349 = -2.73769248870524d-5*x173 - 1.20264153485538d-5*x175 - 0.00
     &344377233585999d0*x178 + x187*x283 + x195 - x213 + 0.0720987654320
     &988d0*x219 - 0.0438808209227215d0*x316
      x350 = x222*x265
      x351 = x96*(-0.658516038102422d0*x134*x14 + 6.31870381979104d0*x13
     &7*x225 - x160 *x18 - x17*(x141 - x145 + x148 - x149))
      x352 = x305/x93**3
      x353 = x232*x61
      x354 = x353*x56
      x355 = x48*x58
      x356 = 0.00548510261534019d0*x355
      x357 = x356*x61
      x358 = -x236 + 0.0274255130767009d0*x355
      x359 = x208*x358
      x360 = -x354 + x357 + x359
      x361 = x231*x360
      x362 = x230*x361 + x242
      x363 = x244*x361 + x247
      x364 = -x243*x362 + x248*x363 + x263
      x365 = gamma_aa*x69
      x366 = 0.00856659785039798d0*gamma_aa
      x367 = 0.0409336016070163d0*x186
      x368 = gamma_aa*x74
      x369 = 2.48756218905473d0*x368
      x370 = 0.101495523471377d0*x253
      x371 = -0.179350162766415d0*x105*x362 - 1.12186710340987d-5*x249 -
     & 0.0409336016070163d0*x250 + x252*x370 + x256*x369 + x262*x367 + 2
     &.48756218905473d0*x363*x76 - 2.55381762006086d-5*x365 - x366*x84
      x372 = gamma_aa*x172
      x373 = -x182 + 0.0497175d0
      x374 = x373*x82
      x375 = 0.00932203125d0*x181
      x376 = x189*x55*x75
      x377 = x354 - x357 - x359
      x378 = x204*x377
      x379 = x340*x378 + x63*(-x246 + 0.00222451383844352d0*x355) + x66
      x380 = x190*x379
      x381 = 1.24378109452736d0*x334
      x382 = 0.00711844428301927d0*x142
      x383 = 0.0720987654320988d0*x60
      x384 = x345 + x378*x383 + x63*(0.000889805535377408d0*x355 - x382*
     &x56)
      x385 = x191*x74
      x386 = -x239 + 0.0123414808845154d0*x355
      x387 = x216*x386 + x217*x378 + x238
      x388 = 0.0720987654320988d0*x64
      x389 = -gamma_aa*x385 + 0.00344377233585999d0*gamma_aa*x84 + 0.016
     &4553078460206d0*x186*x260*x374 + 4.50990575570768d-6*x249 + 0.0164
     &553078460206d0*x250 + 1.02663468326447d-5*x365 - x379*x76 + x387*x
     &388
      x390 = 0.135327364628503d0*x179*x196
      x391 = x202*x55
      x392 = x134*x58
      x393 = x210*x386
      x394 = x142*x341
      x395 = -x342*x358 - 0.0164553078460206d0*x355*x61 + x394*x56
      x396 = x343*x395
      x397 = x207*x232
      x398 = x204*(0.0146269403075738d0*x203*x358 - x205*x358 + x206*x39
     &2 - x207*x356 - x208*(0.585077612302953d0*x391 - 0.073134701537869
     &2d0*x392) - 0.117015522460591d0*x391*x61 + x397*x56)
      x399 = 2.0d0*x215
      x400 = x358*x60*(0.087761641845443d0*x235 - 0.0109702052306804d0*x
     &355)
      x401 = x377*x386
      x402 = 1d0/x62
      x403 = x377*x395
      x404 = x402*x60
      x405 = x274*x58
      x406 = gamma_aa*x405 + x353
      x407 = x272*x406
      x408 = x230*x407 + x271
      x409 = 0.180246913580247d0*x406
      x410 = x101*x272*x409 + x245
      x411 = x276*x408 - x277*x410 + x286
      x412 = x402*x406
      x413 = x217*x412
      x414 = x404*x409
      x415 = x63*(x245 - x414)
      x416 = x192*x281 + x285 + x388*x63*(x271 - x413) - x415*x76
      x417 = x176*x279
      x418 = x185*x284
      x419 = x328 - x375 + 0.099435d0
      x420 = x278*x70
      x421 = x189*x280
      x422 = 0.00962763222425977d0*tau_a
      x423 = x402*(-x273*x422 + 0.00521496745480737d0*x365*x58 - x397 + 
     &0.0731347015378692d0*x48*x61)
      x424 = x209*x402
      x425 = 0.00361036208409741d0*x312*x58 + x394
      x426 = x209*x425/x62**2
      x427 = gamma_aa*x420
      x428 = 0.0153120874502694d0*x260*x279*x87
      x429 = 0.9305257363491d0*x104
      x430 = x190*x410
      x431 = 0.0720987654320988d0*x101*x407 + x382
      x432 = x101*(-x167*x422*x55 + x232*x358 + x405)
      x433 = x272*x360*x425
      x434 = 0.9305257363491d0*x231
      x435 = -x104*(0.000889805535377408d0*x100 + 0.0324444444444444d0) 
     &+ x346
      x436 = x110*x53
      x437 = x326*x70*x82
      x438 = x167*x60
      x439 = x142*x406
      x440 = 1.8610514726982d0*x204
      Ea = -x0*(-4.9628039271952d0*x46 + x47*(1.49628538404935d0*x51 - 1
     &.49628538404935d0*x92))
      d1Ea(1) = -rho_a*(-1.65426797573173d0*x114 + x116*x118 + x132*x133
     & - x164*x165 + x223) - x112*x47*(x111 + x98) + x96*x97*(1.804d0 + 
     &0.804d0/x93)
      d1Ea(2) = rho_a*x270
      d1Ea(3) = -x288*x290
      d2Ea(1) = rho_a*(1.59154943091895d0*mu*rho_a**(-2.0d0)*x118 + mu*x
     &133*(-x125*( -x122*x38 + x129) + x130*(102.585381117795d0*x113*x29
     &8*x38 + 3.37614805358933d0*x120*x299*x38 - 20.256888321536d0*x121*
     &x126 - 3.07050222784304d0*x128*x224 - 46.6491155403329d0*x152*x298
     & - 0.171037658470723d0*x295 + x296*x41 + x40*(0.114025105647149d0*
     & x146 - x296*x3)) + 0.0570125528235743d0*x146*x43) - 3.30853595146
     &347d0*x113*x291 + 2.66006290497663d0*x113*x293 - 0.76560437251347d
     &0*x118*x146*x3 - 1.28278243853042d0*x123*x223 - x132*x162*x165*x29
     &4 + x163*x292*(18.9561114593731d0*rho_a**( -4.0d0)*x138 + 1.317032
     &07620484d0*x134*x150 - 3.73159088258039d0* x135*x260 + 0.137257065
     &627884d0*x135*x3*x69 - 12.6374076395821d0* x137*x151 + x159*x161*x
     &304*(60.7706649646079d0*x119*x153 + 8.0d0* x140 - 0.17552328369088
     &6d0*x144 + 0.0231063173382234d0*x147 - 9.45308720482942d0*x152 - 0
     &.175523283690886d0*x155 + 1.21262542562518d0*x156 - 36.84602673411
     &65d0*x158) + x160*x300 - x160*x301 + x160*x303 + x161*x28*(6.30205
     &813655294d0*x127*x15 - 47.8725329504926d0*x127*x299*x4 - 186.59646
     &2161332d0*x140*x25* x297 + 6.14100445568608d0*x142*x311 + x154*x29
     &9 + 307.756143353386d0*x157*x310 + 0.146269403075738d0*x20*x307 + 
     &12.2820089113722d0*x224*x311 - 0.808416950416787d0*x26*x295 - 7.55
     &555555555556d0*x306 + 0.263284925536329d0*x308 - 0.026957370227927
     &3d0*x309) + x17*(-0.000677173162387677d0*x11* x307 - 0.43764882176
     &158d0*x153*x299 + 0.146090534979424d0*x306 - 0.0219404104613607d0*
     &x308 + 0.00224644751899394d0*x309)) - x164* x305*x82*(x138*x160*x3
     &02 + x300 - x301 + x303)/x32**3 - 1.10284531715449d0*x224*x46 + x4
     &7*(-0.00042952697051865d0*rho_a** (-7.0d0)*x313*x82 + 0.0537785010
     &593382d0*x166*x312 - x169*x314 - 0.997523589366235d0*x220*x315 + x
     &314*(x111 + 1.804d0) + x319*( -6.81018032016229d-5*x173 - 2.991645
     &60909299d-5*x175 - 0.00856659785039798d0*x178 + 0.0409336016070163
     &d0*x184*x187 + 2.48756218905473d0*x195 + x199*x317 - 2.48756218905
     &473d0*x213 + 0.179350162766415d0*x219 - 0.109156270952044d0*x316) 
     &- x350*( 0.000134624052409184d0*rho_a**(-10.0d0)*x79 - 0.109156270
     &952044d0 *x184*x186*x325 - 1.24378109452736d0*x192*x82*( 0.0022865
     &6265326169d0*x320 + 0.0084375d0*x323) + x193**2*x330* x331 + x196*
     &x317*x349 - 0.496200336971176d0*x199*x324 - x248*( 0.0193340462007
     &931d0*x324 + 0.360493827160494d0*x338 + x339*x340 + x340*x344 + x6
     &3*(-0.0949125904402569d0*x336 + 0.0217508019758922d0*x337)) + x259
     &*x325*x329 + 0.0220860757144586d0*x321*x332 + 0.000215655710138473
     &d0*x321 + 0.012849896775597d0*x323*x82 + 0.200119830078746d0*x324*
     &x88 + 2.48756218905473d0*x333*x334 + x335*(-0.00210916867645015d0*
     &x179 + 0.0720987654320988d0*x211 + x63*(0.00711844428301927d0*x202
     & - 0.00237281476100642d0*x203)) + x348*(x216*(-0.526569851072658d0
     &* x336 + 0.120672257537484d0*x337) + x217*x339 + x217*x344 + 0.214
     &528457844416d0*x324 + 4.0d0*x338)))) + 3.30853595146347d0* x114 - 
     &2.38732414637843d0*x115*x118 - 9.9256078543904d0*x21*x291 - 2.0d0*
     &x223 + x293*x294
      d2Ea(2) = rho_a*(-25.2122244578404d0*gamma_aa*x226*x351*x352 - 4.2
     &0203740964007d0* rho_a**(-2.66666666666667d0)*x225*x227*x96 + 0.64
     &1391219265211d0* x123*x47*(0.00672231263241728d0*x229 + x264*x350)
     & + x131*x14*x224 *x228*x24 + x163*x227*x351 - x47*(-0.000161072613
     &944494d0* gamma_aa*rho_a**(-6.0d0)*x313 + 0.015685396142307d0*x168
     & - 0.498761794683118d0*x315*x364 + x319*x371 + x350*(-x134*x180 + 
     &x134*x200 - 4.48746841363948d-5*x174*x82 - x177*x366 - 0.054578135
     &4760218d0*x187*x374 + x188*x261 + x188*x82*(-x328 + x375 - 0.09943
     &5d0) + x192*x257*(0.000722072416819483d0*x172 + 0.005625d0*x177) -
     & x193*x331*x376 - x248*(-x201 + 0.180246913580247d0*x215*x378 + x3
     &40*x396 + x340*x398 + 0.180246913580247d0*x393 + x63*(0.0474562952
     &201284d0*x391 - 0.00593203690251605d0*x392)) - x254*x349 - x258*x3
     &33 - 0.00828227839292197d0*x332*x372 + x335*x384 + x348*(-x214 + x
     &216* (0.263284925536329d0*x391 - 0.0329106156920411d0*x392) + x217
     &* x396 + x217*x398 + x378*x399 + 2.0d0*x393) - 6.81018032016229d-5
     &* x372 + x380*x381 + x389*x390))) + x270
      d2Ea(3) = rho_a*(79.6542894935927d0*rho_a**(-3.66666666666667d0)*x
     &14**2*x310*x352* x96 - x269*(3.02011151145926d-5*rho_a**(-5.0d0)*x
     &313 - x267*x318* x364*x371 + x268*(1.12186710340987d-5*gamma_aa*x8
     &0 - x204*x248*( x340*x402*x403 + 0.180246913580247d0*x400 + 0.3604
     &93827160494d0* x401) + 2.48756218905473d0*x204*x347*(x400 + 2.0d0*
     &x401 + x403* x404) + x251 + x255 + x259*x324*(-0.007768359375d0*x1
     &81 + x328 + 0.0497175d0) + x324*x367*x373 + 0.00310585439734574d0*
     &x332*x69 - x369*x380 - x370*x389 + 1.24378109452736d0*x384*x387 - 
     &1.24378109452736d0*x385 + x78/x82)))
      d2Ea(4) = -x289*(rho_a*x171*(x116*x416 + x170*x416 - 3.72210294539
     &64d0*x21*x220* x416*x53*x92 + x221*(-0.0545781354760218d0*x184*x41
     &8 + x190*x381* x415 - 0.00699626865671642d0*x192*x417 + 0.00349813
     &432835821d0* x193*x330*x420*x68 - x248*x63*(-x215*x402*x409 + x245
     &*x424 + x340 *x423 - x340*x426 - 0.0296601845125803d0*x48) + 0.003
     &49813432835821d0*x333*x421 + x335*x63*(x382 - x383*x412) + x348*x6
     &3*(x217*x423 - x217*x426 + x271*x424 - x399*x412 - 0.3291061569204
     &11d0*x48) + x390*x416 + 0.00856659785039798d0*x417 + 0.02046680080
     &35082d0*x418*x419*x79)) + x112*x265*x411)
      d2Ea(5) = x0*x266*x269*(0.00261710363348184d0*x107*x376*x420 - x18
     &3*x428 - 0.00261710363348184d0*x190*x363*x421 + 0.0379665670791648
     &d0*x253* x287 - x264*x288*x436 - 0.00523420726696369d0*x282*x427 -
     & x362* x429*x431 + x368*x429*x430 + x419*x428 + 0.0064090375772895
     &5d0* x427 - x434*x435*(-x230*x433 + x241*x406 - x271*x360 + 2.0d0*
     &x432 ) + x434*x76*(x240*x409 - x244*x433 - x245*x360 + 0.180246913
     &580247d0*x432))
      d2Ea(6) = x290*(-0.0104684145339274d0*x104*x421*x430 - 0.030624174
     &9005388d0*x260* x327*x329 - 1.8610514726982d0*x272*x408*x431 + 0.0
     &157026218008911d0*x282*x437 + 4.62948127536866d0*x411**2*x436 + x4
     &35*x440*(-x413*x425 + 0.0385105288970391d0*x438 + 0.39492738830449
     &4d0*x439) - 0.0192271127318687d0*x437 - x440*x76* (-x414*x425 + 0.
     &0034707019870171d0*x438 + 0.0355922214150963d0* x439) - 1.47212079
     &383354d-5*x107*x322*x330/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.2
C   Generated: September 13, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102, x103, x104
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113, x114, x115
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125, x126
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227, x228, x229, x230, x231, x232, x233, x234, x235, x236
     &, x237, x238, x239, x240, x241, x242, x243, x244, x245, x246, x247
     &, x248, x249, x250, x251, x252, x253, x254, x255, x256, x257, x258
     &, x259, x260, x261, x262, x263, x264, x265, x266, x267, x268, x269
     &, x270, x271, x272, x273, x274, x275, x276, x277, x278, x279, x280
     &, x281, x282, x283, x284, x285, x286, x287, x288, x289, x290, x291
     &, x292, x293, x294, x295, x296, x297, x298, x299, x300, x301, x302
     &, x303, x304, x305, x306, x307, x308, x309, x310, x311, x312, x313
     &, x314, x315, x316, x317, x318, x319, x320, x321, x322, x323, x324
     &, x325, x326, x327, x328, x329, x330, x331, x332, x333, x334, x335
     &, x336, x337, x338, x339, x340, x341, x342, x343, x344, x345, x346
     &, x347, x348, x349, x350, x351, x352, x353, x354, x355, x356, x357
     &, x358
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = mu**2
      x1 = rho_a**0.666666666666667d0
      x2 = x0/x1
      x3 = exp(-0.31265084907439067d0*x2)
      x4 = rho_a**(-2.66666666666667d0)
      x5 = gamma_aa*x4
      x6 = 1d0/x0
      x7 = 0.844037013397333d0*x1
      x8 = mu**(-4)
      x9 = 14.4026298885989d0*rho_a**2.0d0*x8
      x10 = rho_a**1.33333333333333d0
      x11 = 0.213719543995407d0*x10*x6
      x12 = x11 - x7 + x9
      x13 = x12*x6
      x14 = 0.0519881082712438d0*x13*x3*x5
      x15 = -x14 + 1.0d0
      x16 = 1d0/x15
      x17 = rho_a**1.0d0
      x18 = x17*x6
      x19 = rho_a**0.333333333333333d0
      x20 = 1.8610514726982d0*x19
      x21 = exp(-1.9241736577956332d0*mu/x19)
      x22 = 0.00449266744437808d0*x5 + 1.0d0
      x23 = 0.804d0/x22
      x24 = 0.0204006002177468d0*x5 + 1.0d0
      x25 = x24**(-2)
      x26 = 0.0109702052306804d0*x5
      x27 = 1d0/gamma_aa
      x28 = tau_a*x27
      x29 = 8.0d0*rho_a*x28
      x30 = x29 - 1.0d0
      x31 = x30*x5
      x32 = 0.0123414808845154d0*x31 - 0.45d0
      x33 = 0.0274255130767009d0*x31 - 1.0d0
      x34 = x26*x30*x33 + 1.0d0
      x35 = x34**(-0.5d0)
      x36 = x26 + x32*x35
      x37 = x36**2
      x38 = 0.00197734563417202d0*x4
      x39 = gamma_aa*x38
      x40 = x35*(0.00222451383844352d0*x31 - 0.0811111111111111d0) + x39
     &
      x41 = rho_a**(-5.33333333333333d0)
      x42 = rho_a**(-2)
      x43 = tau_a**(-2)
      x44 = 0.0028125d0*x43
      x45 = 0.000135388578153653d0*x41 + x42*x44
      x46 = sqrt(x45)
      x47 = abs(gamma_aa)
      x48 = x46*x47
      x49 = x40*x48
      x50 = 1.24378109452736d0*x49
      x51 = gamma_aa**3
      x52 = rho_a**(-8.0d0)
      x53 = x51*x52
      x54 = gamma_aa**2
      x55 = x41*x54
      x56 = x42*x43
      x57 = x54*x56
      x58 = 0.015625d0*x57 + 1.0d0
      x59 = x58**(-2)
      x60 = 0.02485875d0*x57*x59 + 0.123456790123457d0
      x61 = x5*x60
      x62 = 1.86977850568312d-6*x53 + 6.38454405015215d-6*x55 + 0.002141
     &6494625995d0 *x57 + 0.0204668008035082d0*x61
      x63 = x25*(0.0896750813832074d0*x37 - x50 + x62) + 1.0d0
      x64 = 1d0/x63
      x65 = -0.804d0*x64
      x66 = x21*(x23 + x65)
      x67 = 3.14159265358979d0*x18*(-0.804d0*x16 + 1.804d0) + x20*x66
      x68 = rho_a**(-4.33333333333333d0)
      x69 = -x11 + x7 - x9
      x70 = x68*x69
      x71 = rho_a**(-3.66666666666667d0)
      x72 = x6*x69
      x73 = x71*x72
      x74 = 0.562691342264888d0*rho_a**(-0.333333333333333d0)
      x75 = 28.8052597771978d0*x8
      x76 = x17*x75
      x77 = 0.284959391993876d0*rho_a**0.333333333333333d0*x6
      x78 = x74 - x76 - x77
      x79 = x6*x78
      x80 = 0.0519881082712438d0*x4
      x81 = 0.0108360841285172d0*x70 - 0.138634955389984d0*x73 + x79*x80
     &
      x82 = gamma_aa*x3/x15**2
      x83 = x81*x82
      x84 = 2.52584049348619d0*x18
      x85 = mu/x17
      x86 = 1.19366207318922d0*x85
      x87 = x22**(-2)
      x88 = rho_a**(-3.33333333333333d0)
      x89 = x87*x88
      x90 = -x23 + 1.804d0
      x91 = rho_a**(-0.666666666666667d0)
      x92 = 0.6203504908994d0*x91
      x93 = x63**(-2)
      x94 = rho_a**(-6.33333333333333d0)
      x95 = x54*x94
      x96 = rho_a**(-9.0d0)
      x97 = x51*x96
      x98 = rho_a**(-3)
      x99 = x43*x98
      x100 = x54*x99
      x101 = gamma_aa*x71
      x102 = 0.0545781354760218d0*x60
      x103 = x57/x58
      x104 = 0.001553671875d0*x103
      x105 = x104 - 0.0497175d0
      x106 = x105*x51
      x107 = rho_a**(-5.66666666666667d0)
      x108 = x43*x59
      x109 = x107*x108
      x110 = 0.0204668008035082d0*x109
      x111 = 1d0/x46
      x112 = 1d0/x47
      x113 = x112*x40
      x114 = x111*x113
      x115 = x44*x98 + 0.000361036208409741d0*x94
      x116 = x115*x54
      x117 = x114*x116
      x118 = 1d0/x24
      x119 = 1.50330191856923d-6*x53 + 5.13317341632233d-6*x55 + 0.00172
     &188616792999d0*x57 + 0.0164553078460206d0*x61
      x120 = x119 + 0.0720987654320988d0*x37 - x49
      x121 = x118*x120
      x122 = 0.135327364628503d0*x121
      x123 = 0.00527292169112538d0*x71
      x124 = tau_a*x4
      x125 = x101*x30
      x126 = x34**(-1.5d0)
      x127 = 0.0438808209227215d0*x124
      x128 = 0.0146269403075738d0*x33
      x129 = 0.219404104613608d0*x124 - 0.0731347015378692d0*x125
      x130 = 0.00548510261534019d0*x31
      x131 = x125*x128 - x127*x33 - x129*x130
      x132 = x126*x131
      x133 = x132*x32
      x134 = -gamma_aa*x123 + 0.180246913580247d0*x133 + x35*(0.01779611
     &07075482d0* x124 - 0.00593203690251605d0*x125)
      x135 = x134*x48
      x136 = 0.0585077612302953d0*x71
      x137 = 0.0987318470761234d0*x124 - 0.0329106156920411d0*x125
      x138 = 2.0d0*x35
      x139 = 2.0d0*x32
      x140 = -gamma_aa*x136 + x132*x139 + x137*x138
      x141 = x140*x36
      x142 = -0.00428329892519899d0*x100 - x101*x102 + x101*x122 + x106*
     &x110 + 1.24378109452736d0*x117 - 1.24378109452736d0*x135 + 0.08967
     &50813832074d0*x141 - 3.40509016008115d-5*x95 - 1.49582280454649d-5
     &*x97
      x143 = 1.49628538404935d0*x19
      x144 = x143*x25
      x145 = x21*(0.0179261670197794d0*gamma_aa*x89 + x142*x144*x93 - x9
     &0*x92 + x92*( x65 + 1.804d0))
      x146 = rho_a**(-1.66666666666667d0)
      x147 = x12*x146*x8
      x148 = x14 - 1.0d0
      x149 = x3/x148**2
      x150 = 0.131313669051252d0*x149
      x151 = rho_a**(-2.33333333333333d0)
      x152 = x151*x87
      x153 = -x29 + 1.0d0
      x154 = x153*x5
      x155 = 0.0123414808845154d0*x154 + 0.45d0
      x156 = 2.0d0*x155
      x157 = 0.0274255130767009d0*x154 + 1.0d0
      x158 = x153*x157*x26 + 1.0d0
      x159 = x158**(-1.5d0)
      x160 = 0.0438808209227215d0*x146
      x161 = x157*x160
      x162 = x153*x4
      x163 = x146*x28
      x164 = 0.219404104613608d0*x163
      x165 = x159*(0.00548510261534019d0*x154*(0.0274255130767009d0*x162
     & + x164) + 0.00548510261534019d0*x157*x162 + x161*x28)
      x166 = 0.0219404104613608d0*x4
      x167 = x158**(-0.5d0)
      x168 = 0.0987318470761234d0*x163
      x169 = 0.0123414808845154d0*x162 + x168
      x170 = 2.0d0*x169
      x171 = x166 - x167*x170
      x172 = -x155*x167 + x26
      x173 = 0.0896750813832074d0*x172
      x174 = 0.180246913580247d0*x155
      x175 = 0.0177961107075482d0*x146
      x176 = x175*x28
      x177 = -x167*(0.00222451383844352d0*x162 + x176) + x38
      x178 = 1.24378109452736d0*x48
      x179 = x52*x54
      x180 = x4*x60
      x181 = 1.27690881003043d-5*x41
      x182 = x172**2
      x183 = -x167*(0.00222451383844352d0*x154 + 0.0811111111111111d0) +
     & x39
      x184 = x183*x48
      x185 = x119 + 0.0720987654320988d0*x182 - x184
      x186 = x118*x4
      x187 = 0.0507477617356885d0*x186
      x188 = 0.00428329892519899d0*x56
      x189 = x112*x183
      x190 = 1.24378109452736d0*gamma_aa
      x191 = x190*x46
      x192 = 0.0204668008035082d0*x108
      x193 = rho_a**(-4.66666666666667d0)
      x194 = x105*x54
      x195 = x193*x194
      x196 = -gamma_aa*x181 - gamma_aa*x188 - 5.60933551704935d-6*x179 -
     & 0.0204668008035082d0*x180 + x185*x187 + x189*x191 + x192*x195
      x197 = -x173*(x156*x165 + x171) + x178*(x165*x174 + x177) + x196
      x198 = x25*(0.0896750813832074d0*x182 - 1.24378109452736d0*x184 + 
     &x62) + 1.0d0
      x199 = x198**(-2)
      x200 = x199*x25
      x201 = 0.748142692024676d0*x19
      x202 = x200*x201
      x203 = 2.0d0*x21
      x204 = x147*x150 + x203*(0.00336115631620864d0*x152 + x197*x202)
      x205 = 0.197463694152247d0*x146
      x206 = 1d0/x158
      x207 = 0.00120345402803247d0*x68
      x208 = x206*(gamma_aa*x153*x207 + x161)
      x209 = 0.0720987654320988d0*x167*x172
      x210 = x167*x48
      x211 = tau_a**(-3)
      x212 = x211*x54
      x213 = x212*x42
      x214 = 0.0028125d0*x213
      x215 = x111*x189
      x216 = 0.0164553078460206d0*x106
      x217 = x211*x59
      x218 = x193*x216*x217 - 0.00344377233585999d0*x213
      x219 = x214*x215 + x218
      x220 = x209*(-x156*x208 + x205) - x210*(-x174*x208 + x175) + x219
      x221 = 1.8610514726982d0*x220
      x222 = x21*x25
      x223 = x10*x199*x222
      x224 = 5.05168098697239d0*x83
      x225 = x6*x80
      x226 = rho_a**(-6.0d0)
      x227 = exp(-0.62530169814878134d0*x2)
      x228 = mu/x10
      x229 = gamma_aa*x68
      x230 = x22**(-3)
      x231 = 0.413566993932933d0*rho_a**(-1.66666666666667d0)
      x232 = 1d0/x198
      x233 = x200*x91
      x234 = x101*x60
      x235 = 0.270654729257005d0*x101
      x236 = 1/(x198**3*x24**4)
      x237 = x142*x143*x236
      x238 = rho_a**(-7.33333333333333d0)
      x239 = x238*x54
      x240 = rho_a**(-4)
      x241 = x240*x43
      x242 = gamma_aa*x193
      x243 = rho_a**(-6.66666666666667d0)
      x244 = tau_a**(-4)
      x245 = x244*x59
      x246 = 0.00014565673828125d0*gamma_aa**4*x240*x245
      x247 = x51*(-0.010875703125d0*x103 + x246 + 0.1491525d0)
      x248 = x45**(-1.5d0)*x47
      x249 = 1.24378109452736d0*x40
      x250 = x120*x25
      x251 = x112*x134
      x252 = x111*x116
      x253 = 1.24378109452736d0*x140
      x254 = tau_a*x71
      x255 = x242*x30
      x256 = x132*x137
      x257 = x126*(-0.087761641845443d0*x124*x129 + 0.0292538806151477d0
     &*x125*x129 - x130*(-1.17015522460591d0*x254 + 0.26816057230552d0*x
     &255) + 0.234031044921181d0*x254*x33 - 0.0536321144611041d0*x255*x3
     &3)
      x258 = 0.180246913580247d0*x32
      x259 = 0.131642462768165d0*x33
      x260 = 0.0164553078460206d0*x31
      x261 = x131*x34**(-2.5d0)
      x262 = x261*(-x124*x259 + 0.0438808209227215d0*x125*x33 - x129*x26
     &0)
      x263 = 0.000790938253668807d0*x4
      x264 = gamma_aa*x263
      x265 = x264 + x35*(0.000889805535377408d0*x31 - 0.0324444444444444
     &d0)
      x266 = 1.24378109452736d0*x265
      x267 = -0.00344377233585999d0*x100 + x109*x216 + x117 - x135 + 0.0
     &720987654320988d0*x141 - 0.0438808209227215d0*x234 - 2.73769248870
     &524d-5*x95 - 1.20264153485538d-5*x97
      x268 = x144*x199
      x269 = -0.0108360841285172d0*x12*x68 + 0.138634955389984d0*x13*x71
     & - x225*(-x74 + x76 + x77)
      x270 = x227/x148**3
      x271 = x160*x33
      x272 = x271*x28
      x273 = x30*x4
      x274 = 0.00548510261534019d0*x273
      x275 = x274*x33
      x276 = -x164 + 0.0274255130767009d0*x273
      x277 = x130*x276
      x278 = -x272 + x275 + x277
      x279 = x159*x278
      x280 = x156*x279 + x171
      x281 = x174*x279 + x177
      x282 = -x173*x280 + x178*x281 + x196
      x283 = gamma_aa*x41
      x284 = 0.00856659785039798d0*gamma_aa
      x285 = 0.0409336016070163d0*x108
      x286 = gamma_aa*x46
      x287 = 2.48756218905473d0*x286
      x288 = 0.101495523471377d0*x186
      x289 = -0.179350162766415d0*x172*x280 - 1.12186710340987d-5*x179 -
     & 0.0409336016070163d0*x180 + x185*x288 + x189*x287 + x195*x285 + 2
     &.48756218905473d0*x281*x48 - 2.55381762006086d-5*x283 - x284*x56
      x290 = gamma_aa*x94
      x291 = -x104 + 0.0497175d0
      x292 = x291*x54
      x293 = 0.00932203125d0*x103
      x294 = x111*x27*x47
      x295 = x272 - x275 - x277
      x296 = x126*x295
      x297 = x258*x296 + x35*(-x176 + 0.00222451383844352d0*x273) + x38
      x298 = x112*x297
      x299 = 1.24378109452736d0*x252
      x300 = 0.00711844428301927d0*x146
      x301 = 0.0720987654320988d0*x32
      x302 = x263 + x296*x301 + x35*(0.000889805535377408d0*x273 - x28*x
     &300)
      x303 = x113*x46
      x304 = -x168 + 0.0123414808845154d0*x273
      x305 = x138*x304 + x139*x296 + x166
      x306 = 0.0720987654320988d0*x36
      x307 = -gamma_aa*x303 + 0.00344377233585999d0*gamma_aa*x56 + 0.016
     &4553078460206d0*x108*x193*x292 + 4.50990575570768d-6*x179 + 0.0164
     &553078460206d0*x180 + 1.02663468326447d-5*x283 - x297*x48 + x305*x
     &306
      x308 = 0.135327364628503d0*x101*x118
      x309 = x124*x27
      x310 = x30*x71
      x311 = x132*x304
      x312 = x146*x259
      x313 = -x260*x276 - 0.0164553078460206d0*x273*x33 + x28*x312
      x314 = x261*x313
      x315 = x129*x160
      x316 = x126*(0.0146269403075738d0*x125*x276 - x127*x276 + x128*x31
     &0 - x129*x274 - x130*(0.585077612302953d0*x309 - 0.073134701537869
     &2d0*x310) + x28*x315 - 0.117015522460591d0*x309*x33)
      x317 = 2.0d0*x137
      x318 = x276*x32*(0.087761641845443d0*x163 - 0.0109702052306804d0*x
     &273)
      x319 = x295*x304
      x320 = 1d0/x34
      x321 = x295*x313
      x322 = x32*x320
      x323 = x207*x30
      x324 = gamma_aa*x323 + x271
      x325 = x206*x324
      x326 = x156*x325 + x205
      x327 = 0.180246913580247d0*x324
      x328 = x155*x206*x327 + x175
      x329 = x209*x326 - x210*x328 + x219
      x330 = x320*x324
      x331 = x139*x330
      x332 = x322*x327
      x333 = x35*(x175 - x332)
      x334 = x114*x214 + x218 + x306*x35*(x205 - x331) - x333*x48
      x335 = x212*x98
      x336 = x107*x217
      x337 = x246 - x293 + 0.099435d0
      x338 = x211*x42
      x339 = x111*x213
      x340 = 0.00962763222425977d0*tau_a
      x341 = x320*(0.00521496745480737d0*x283*x30 - x315 + 0.07313470153
     &78692d0*x33* x4 - x340*x68)
      x342 = x131*x320
      x343 = 0.00361036208409741d0*x229*x30 + x312
      x344 = x131*x343/x34**2
      x345 = gamma_aa*x338
      x346 = 0.0153120874502694d0*x193*x212*x59
      x347 = 0.9305257363491d0*x167
      x348 = x112*x328
      x349 = 0.0720987654320988d0*x155*x325 + x300
      x350 = x155*(x160*x276 - x27*x340*x88 + x323)
      x351 = x206*x278*x343
      x352 = 0.9305257363491d0*x159
      x353 = -x167*(0.000889805535377408d0*x154 + 0.0324444444444444d0) 
     &+ x264
      x354 = x232*x25
      x355 = x244*x42*x54
      x356 = x32*x88
      x357 = x146*x324
      x358 = 1.8610514726982d0*x126
      Ea = -rho_a*x67
      d1Ea(1) = -rho_a*(x145 + x6*(-2.52584049348619d0*x16 + 5.667433147
     &07599d0) + x66* x86 + x83*x84) - x67
      d1Ea(2) = rho_a*x204
      d1Ea(3) = -x221*x223
      d2Ea(1) = rho_a*(1.59154943091895d0*mu*rho_a**(-2.0d0)*x66 - 0.765
     &60437251347d0*x0 *x151*x66 - 1.28278243853042d0*x145*x228 + x21*( 
     &-0.00042952697051865d0*rho_a**(-7.0d0)*x230*x54 - 0.99752358936623
     &5d0*x142*x233 + 0.0537785010593382d0*x229*x87 - x231*x90 + x231*(-
     &0.804d0*x232 + 1.804d0) + x237*( -0.00856659785039798d0*x100 + 0.0
     &409336016070163d0*x106*x109 + 2.48756218905473d0*x117 + x121*x235 
     &- 2.48756218905473d0*x135 + 0.179350162766415d0*x141 - 0.109156270
     &952044d0*x234 - 6.81018032016229d-5*x95 - 2.99164560909299d-5*x97)
     & - x268*( 0.000134624052409184d0*rho_a**(-10.0d0)*x51 - 0.10915627
     &0952044d0 *x106*x108*x243 - 1.24378109452736d0*x114*x54*( 0.002286
     &56265326169d0*x238 + 0.0084375d0*x241) + x115**2*x248* x249 + x118
     &*x235*x267 - 0.496200336971176d0*x121*x242 - x178*( 0.019334046200
     &7931d0*x242 + 0.360493827160494d0*x256 + x257*x258 + x258*x262 + x
     &35*(-0.0949125904402569d0*x254 + 0.0217508019758922d0*x255)) + x19
     &2*x243*x247 + 0.0220860757144586d0*x239*x250 + 0.00021565571013847
     &3d0*x239 + 0.012849896775597d0*x241*x54 + 0.200119830078746d0*x242
     &*x60 + 2.48756218905473d0*x251*x252 + x253*(-0.00210916867645015d0
     &*x101 + 0.0720987654320988d0*x133 + x35*(0.00711844428301927d0*x12
     &4 - 0.00237281476100642d0*x125)) + x266*(x138*(-0.526569851072658d
     &0* x254 + 0.120672257537484d0*x255) + x139*x257 + x139*x262 + 0.21
     &4528457844416d0*x242 + 4.0d0*x256))) - x224*x6 - x82*x84*( 0.00225
     &860726894828d0*x0*x226*x69 + 0.508328169763273d0*x193*x72 - x225*(
     &0.187563780754963d0*rho_a**(-1.33333333333333d0) + 0.0949864639979
     &585d0*x6*x91 + x75) - 0.0758525888996201d0*x41*x69 + 0.02167216825
     &70343d0*x68*x78 - 0.277269910779967d0*x71*x79) + x227*x54*x81*x84*
     &(0.103976216542488d0*x4*x79 + 0.0216721682570343d0*x70 - 0.2772699
     &10779967d0*x73)/x15**3) - 2.0d0*x145 - x18*x224 - x6*(-5.051680986
     &97239d0*x16 + 11.334866294152d0) - 2.38732414637843d0*x66*x85
      d2Ea(2) = rho_a*(0.262627338102504d0*gamma_aa*x147*x269*x270 + x12
     &*x150*x4*x8 - x149*x269*x84 + 0.641391219265211d0*x21*x228*( 0.006
     &72231263241728d0*x152 + x197*x268) - x21*( -0.000161072613944494d0
     &*gamma_aa*x226*x230 - 0.498761794683118d0* x233*x282 + x237*x289 +
     & x268*(-x102*x71 - 0.0545781354760218d0* x109*x292 + x110*x194 + x
     &110*x54*(-x246 + x293 - 0.099435d0) + x114*x190*(0.000722072416819
     &483d0*x94 + 0.005625d0*x99) - x115* x249*x294 + x122*x71 - x178*(-
     &x123 + 0.180246913580247d0*x137* x296 + x258*x314 + x258*x316 + 0.
     &180246913580247d0*x311 + x35*( 0.0474562952201284d0*x309 - 0.00593
     &203690251605d0*x310)) - x187* x267 - x191*x251 - 0.008282278392921
     &97d0*x250*x290 + x253*x302 + x266*(-x136 + x138*(0.263284925536329
     &d0*x309 - 0.0329106156920411d0*x310) + x139*x314 + x139*x316 + x29
     &6*x317 + 2.0d0*x311) - x284*x99 - 6.81018032016229d-5*x290 + x298*
     &x299 + x307*x308 - 4.48746841363948d-5*x54*x96) + 0.01568539614230
     &7d0* x89)) + x204
      d2Ea(3) = -rho_a*(x203*(3.02011151145926d-5*rho_a**(-5.0d0)*x230 -
     & x201*x236*x282* x289 + x202*(1.12186710340987d-5*gamma_aa*x52 - x
     &126*x178*(x258* x320*x321 + 0.180246913580247d0*x318 + 0.360493827
     &160494d0*x319) + 2.48756218905473d0*x126*x265*(x318 + 2.0d0*x319 +
     & x321*x322) + x181 + x188 + x192*x242*(-0.007768359375d0*x103 + x2
     &46 + 0.0497175d0) + x242*x285*x291 + 0.00310585439734574d0*x250*x4
     &1 - x287*x298 - x288*x307 + 1.24378109452736d0*x302*x305 - 1.24378
     &109452736d0*x303 + x50/x54)) + 0.0136534984882616d0*x12**2 *x270*x
     &68/mu**6)
      d2Ea(4) = -x222*(rho_a*x93*(-3.7221029453964d0*x142*x19*x25*x334*x
     &64 + x143*( -0.0545781354760218d0*x106*x336 + x112*x299*x333 - 0.0
     &0699626865671642d0*x114*x335 + 0.00349813432835821d0*x115*x248 *x3
     &38*x40 - x178*x35*(-x137*x320*x327 + x175*x342 + x258*x341 - x258*
     &x344 - 0.0296601845125803d0*x4) + 0.00349813432835821d0*x251 *x339
     & + x253*x35*(x300 - x301*x330) + x266*x35*(x139*x341 - x139* x344 
     &+ x205*x342 - x317*x330 - 0.329106156920411d0*x4) + x308* x334 + 0
     &.00856659785039798d0*x335 + 0.0204668008035082d0*x336* x337*x51) +
     & x334*x86 + x334*x92) + x199*x20*x329)
      d2Ea(5) = x10*x200*x203*(-x105*x346 - 0.00261710363348184d0*x112*x
     &281*x339 + 0.00261710363348184d0*x183*x294*x338 + 0.03796656707916
     &48d0*x186* x220 - x197*x221*x354 - 0.00523420726696369d0*x215*x345
     & - x280* x347*x349 + x286*x347*x348 + x337*x346 + 0.00640903757728
     &955d0* x345 - x352*x353*(-x156*x351 + x170*x324 - x205*x278 + 2.0d
     &0*x350 ) + x352*x48*(x169*x327 - x174*x351 - x175*x278 + 0.1802469
     &13580247d0*x350))
      d2Ea(6) = x223*(-0.0104684145339274d0*x167*x339*x348 - 0.030624174
     &9005388d0*x193* x245*x247 - 1.8610514726982d0*x206*x326*x349 + 0.0
     &157026218008911d0*x215*x355 + 4.62948127536866d0*x329**2*x354 + x3
     &53*x358*(-x331*x343 + 0.0385105288970391d0*x356 + 0.39492738830449
     &4d0*x357) - 0.0192271127318687d0*x355 - x358*x48* (-x332*x343 + 0.
     &0034707019870171d0*x356 + 0.0355922214150963d0* x357) - 1.47212079
     &383354d-5*x183*x240*x248/tau_a**6)
      end subroutine



