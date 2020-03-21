C File: srdft/srdftfun.F


C*****************************************************************************
      pure subroutine Convert_d1_alpha_beta_to_C_S(d1Ea, d1Eb,
     &                  functional_type, d1E) 
C*****************************************************************************
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
C*****************************************************************************
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


C*****************************************************************************
      pure subroutine Convert_d2_alpha_beta_to_C_S(d2Ea, d2Eb,
     &                  functional_type, d2E) 
C*****************************************************************************
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
C*****************************************************************************
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


C*****************************************************************************
      subroutine LDAX_ERF(rho, diff_order, mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange LDA functional
C
C   Input : rho, vector of densities
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
!#include "pi.h" --" pi, pi2
#include "pi.h" 
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
            call D1ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_LDA_ERF_case_3(rho_b, mu, Eb, d1Eb)
            end if
         end if

      elseif (diff_order .eq. 2) then
         if (A .lt. 1.d-9) then
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
            call D2ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_LDA_ERF_case_3(rho_b, mu, Eb, d1Eb, d2Eb)
            end if
         end if

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if

C     Do the conversion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting relation
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


C*****************************************************************************
      subroutine PBEX_GWS_ERF(rho, gamma_vec, diff_order, mu, E, d1E, 
     &                         d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange PBE functional
C
C   Input : rho, vector of densities
C           gamma_vec, is called gamma_vec, becuase "gamma" is an
C              intrinsic function and therefore a bad name.
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
!#include "pi.h" --" pi, pi2
#include "pi.h" 
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
         if (A .lt. 1.d-9) then
            call ESRX_PBE_GWS_ERF_case_1(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_GWS_ERF_case_1(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call ESRX_PBE_GWS_ERF_case_2_1(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_GWS_ERF_case_2_1(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call ESRX_PBE_GWS_ERF_case_2_3(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_GWS_ERF_case_2_3(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call ESRX_PBE_GWS_ERF_case_2_2(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_GWS_ERF_case_2_2(rho_b, gamma_bb, mu, Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call ESRX_PBE_GWS_ERF_case_3(rho_a, gamma_aa, mu, Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_PBE_GWS_ERF_case_3(rho_b, gamma_bb, mu, Eb)
            end if
         end if

      elseif (diff_order .eq. 1) then
         if (A .lt. 1.d-9) then
            call D1ESRX_PBE_GWS_ERF_case_1(rho_a, gamma_aa, mu, Ea, 
     &                                     d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_GWS_ERF_case_1(rho_b, gamma_bb, mu, Eb, 
     &                                        d1Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D1ESRX_PBE_GWS_ERF_case_2_1(rho_a, gamma_aa, mu, Ea, 
     &                                       d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_GWS_ERF_case_2_1(rho_b, gamma_bb, mu, Eb,
     &                                      d1Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D1ESRX_PBE_GWS_ERF_case_2_3(rho_a, gamma_aa, mu, Ea,
     &                                       d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_GWS_ERF_case_2_3(rho_b, gamma_bb, mu, Eb,
     &                                      d1Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call D1ESRX_PBE_GWS_ERF_case_2_2(rho_a, gamma_aa, mu, Ea, 
     &                                       d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_GWS_ERF_case_2_2(rho_b, gamma_bb, mu, Eb,
     &                                      d1Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call D1ESRX_PBE_GWS_ERF_case_3(rho_a, gamma_aa, mu, Ea,
     &                                     d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_PBE_GWS_ERF_case_3(rho_b, gamma_bb, mu, Eb, 
     &                                        d1Eb)
            end if
         end if

      elseif (diff_order .eq. 2) then
         if (A .lt. 1.d-9) then
            call D2ESRX_PBE_GWS_ERF_case_1(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_GWS_ERF_case_1(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D2ESRX_PBE_GWS_ERF_case_2_1(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_GWS_ERF_case_2_1(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D2ESRX_PBE_GWS_ERF_case_2_3(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_GWS_ERF_case_2_3(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call D2ESRX_PBE_GWS_ERF_case_2_2(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_GWS_ERF_case_2_2(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call D2ESRX_PBE_GWS_ERF_case_3(rho_a, gamma_aa, mu, Ea,
     &                                   d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_PBE_GWS_ERF_case_3(rho_b, gamma_bb, mu, Eb,
     &                                    d1Eb, d2Eb)
            end if
         end if

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if

C     Do the conversion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting relation
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


C*****************************************************************************
      subroutine TPSSX_GWS_ERF(rho, gamma_vec, tau, diff_order, mu, E,
     &                          d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange TPSS functional
C
C   Input : rho, vector of densities
C           gamma_vec, is called gamma_vec, becuase "gamma" is an
C              intrinsic function and therefore a bad name.
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
!#include "pi.h" --" pi, pi2
#include "pi.h" 
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
         if (A .lt. 1.d-9) then
            call ESRX_TPSS_GWS_ERF_case_1(rho_a, gamma_aa, tau_a, mu, 
     &                                    Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_GWS_ERF_case_1(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call ESRX_TPSS_GWS_ERF_case_2_1(rho_a, gamma_aa, tau_a, mu,
     &                                      Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_GWS_ERF_case_2_1(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call ESRX_TPSS_GWS_ERF_case_2_3(rho_a, gamma_aa, tau_a, mu,
     &                                      Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_GWS_ERF_case_2_3(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call ESRX_TPSS_GWS_ERF_case_2_2(rho_a, gamma_aa, tau_a, mu,
     &                                      Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_GWS_ERF_case_2_2(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call ESRX_TPSS_GWS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     &                                    Ea)
            if (.not.alpha_eq_beta) then
               call ESRX_TPSS_GWS_ERF_case_3(rho_b, gamma_bb, tau_b, mu,
     &                                       Eb)
            end if
         end if

      elseif (diff_order .eq. 1) then
         if (A .lt. 1.d-9) then
            call D1ESRX_TPSS_GWS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     &                                  Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_GWS_ERF_case_1(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb, d1Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D1ESRX_TPSS_GWS_ERF_case_2_1(rho_a, gamma_aa, tau_a,
     &                                        mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_GWS_ERF_case_2_1(rho_b, gamma_bb, tau_b,
     &                                           mu, Eb, d1Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D1ESRX_TPSS_GWS_ERF_case_2_3(rho_a, gamma_aa, tau_a,
     &                                        mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_GWS_ERF_case_2_3(rho_b, gamma_bb, tau_b,
     &                                           mu, Eb, d1Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call D1ESRX_TPSS_GWS_ERF_case_2_2(rho_a, gamma_aa, tau_a,
     &                                        mu, Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_GWS_ERF_case_2_2(rho_b, gamma_bb, tau_b,
     &                                           mu, Eb, d1Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call D1ESRX_TPSS_GWS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     &                                      Ea, d1Ea)
            if (.not.alpha_eq_beta) then
               call D1ESRX_TPSS_GWS_ERF_case_3(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb, d1Eb)
            end if
         end if

      elseif (diff_order .eq. 2) then
         if (A .lt. 1.d-9) then
            call D2ESRX_TPSS_GWS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     &                                      Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_GWS_ERF_case_1(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 0.075d0) then ! case_2_1
            call D2ESRX_TPSS_GWS_ERF_case_2_1(rho_a, gamma_aa, tau_a,
     &                                        mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_GWS_ERF_case_2_1(rho_b, gamma_bb, tau_b,
     &                                           mu, Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 50d0) then ! case_2_3
            call D2ESRX_TPSS_GWS_ERF_case_2_3(rho_a, gamma_aa, tau_a,
     &                                        mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_GWS_ERF_case_2_3(rho_b, gamma_bb, tau_b,
     &                                           mu, Eb, d1Eb, d2Eb)
            end if
         elseif (A .le. 1.d+2) then ! case_2_2
            call D2ESRX_TPSS_GWS_ERF_case_2_2(rho_a, gamma_aa, tau_a,
     &                                        mu, Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_GWS_ERF_case_2_2(rho_b, gamma_bb, tau_b,
     &                                           Eb, d1Eb, d2Eb)
            end if
         elseif (A .lt. 1.d+9) then
            call D2ESRX_TPSS_GWS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     &                                      Ea, d1Ea, d2Ea)
            if (.not.alpha_eq_beta) then
               call D2ESRX_TPSS_GWS_ERF_case_3(rho_b, gamma_bb, tau_b,
     &                                         mu, Eb, d1Eb, d2Eb)
            end if
         end if

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if

C     Do the conversion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting relation
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


C*****************************************************************************
      subroutine PW92C_ERF(rho, diff_order, special_case, mu, 
     &                     E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized PW92 functional
C
C   Input : rho, vector of densities
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h" 
      integer, intent(in) :: diff_order, special_case
      real*8, intent(in) :: rho(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: spin_polarized
      real*8 :: rho_c, rho_s
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)
      
      if (rho_c .le. 1.d-10) then
         return
      end if

      if (abs(rho_s) .le. 1.d-20) then
         spin_polarized = .false.
      else
         spin_polarized = .true.
      end if
      
      if (spin_polarized) then
         if (diff_order .eq. 0) then
            call ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            call D2ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E, d2E)

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      else
         if (diff_order .eq. 0) then
            call ESRC_PW92_ERF(rho_c, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_PW92_ERF(rho_c, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            if (special_case .eq. 1) then
               call D2ESRC_PW92_ERF_singletref_triplet(rho_c, mu, 
     &                                     E, d1E, d2E)
            else
               call D2ESRC_PW92_ERF(rho_c, mu, E, d1E, d2E)
            end if

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      end if
      end subroutine


C*****************************************************************************
      subroutine VWN5C_ERF(rho, diff_order, special_case, mu, 
     &                     E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized VWN5 functional
C
C   Input : rho, vector of densities
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h" 
      integer, intent(in) :: diff_order, special_case
      real*8, intent(in) :: rho(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: spin_polarized
      real*8 :: rho_c, rho_s
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)

      if (rho_c .le. 1.d-10) then
         return
      end if

      if (abs(rho_s) .le. 1.d-20) then
         spin_polarized = .false.
      else
         spin_polarized = .true.
      end if

      if (spin_polarized) then
         if (diff_order .eq. 0) then
            call ESRC_SPIN_VWN5_ERF(rho_c, rho_s, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_SPIN_VWN5_ERF(rho_c, rho_s, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            call D2ESRC_SPIN_VWN5_ERF(rho_c, rho_s, mu, E, d1E, d2E)

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      else
         if (diff_order .eq. 0) then
            call ESRC_VWN5_ERF(rho_c, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_VWN5_ERF(rho_c, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            if (special_case .eq. 1) then
               call D2ESRC_VWN5_ERF_singletref_triplet(rho_c, mu, 
     &                                     E, d1E, d2E)
            else
               call D2ESRC_VWN5_ERF(rho_c, mu, E, d1E, d2E)
            end if

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      end if
      end subroutine


C*****************************************************************************
      subroutine PBEC_GWS_ERF(rho, gamma_xx, diff_order, special_case,
     &                        mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to control calling of correlation functional PBEGWS
C
C   Input : rho, vector of densities
C           gamma_xx, vector of density gradient norms
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C           special_case, keyword to invoke special call.
C                         special_case = 1, singlet-reference triplet
C                                           hessian.
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h" 
      integer, intent(in) :: diff_order, special_case
      real*8, intent(in) :: rho(2), gamma_xx(3), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: spin_polarized
      real*8 :: rho_c, rho_s, gamma_cc
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)
      gamma_cc = gamma_xx(1)

      if (rho_c .le. 1.d-10) then
         return
      end if

      if (abs(rho_s) .le. 1.d-20) then
         spin_polarized = .false.
      else
         spin_polarized = .true.
      end if

      if (spin_polarized) then
         if (diff_order .eq. 0) then
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
      else
         if (diff_order .eq. 0) then
            call ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            if (special_case .eq. 1) then
               call D2ESRC_PBE_GWS_ERF_singletref_triplet(rho_c,
     &                    gamma_cc, mu, E, d1E, d2E)
            else
               call D2ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E, d2E)
            end if
         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      end if
      end subroutine


C*****************************************************************************
      subroutine TPSSC_GWS_ERF(rho, gamma_xx, tau, diff_order,
     &                         special_case ,mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation TPSS functional
C
C   Input : rho, vector of densities
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h" 
      integer, intent(in) :: diff_order, special_case
      real*8, intent(in) :: rho(2), gamma_xx(3), tau(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: spin_polarized
      real*8 :: rho_c, rho_s, gamma_cc, gamma_ss, gamma_cs, tau_c
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)
      gamma_cc = gamma_xx(1)
      gamma_ss = gamma_xx(2)
      gamma_cs = gamma_xx(3)
      tau_c = tau(1)

      if (rho_c .le. 1.d-10) then
         return
      end if

      if ((abs(rho_s) .le. 1.d-20) .and. 
     &    (abs(gamma_ss) .le. 1.d-20) .and. 
     &    (abs(gamma_cs) .le. 1.d-20)) then
         spin_polarized = .false.
      else
         spin_polarized = .true.
      end if

      if (spin_polarized) then
         if (diff_order .eq. 0) then
            call ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc,
     &                                 gamma_ss, gamma_cs, tau_c, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc,
     &                            gamma_ss, gamma_cs, tau_c, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            call D2ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc,
     &                       gamma_ss, gamma_cs, tau_c, mu, E, d1E, d2E)

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      else
         if (diff_order .eq. 0) then
            call ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E, d1E)

         elseif (diff_order .eq. 2) then
            if (special_case .eq. 1) then
               call D2ESRC_TPSS_GWS_ERF_singletref_triplet(rho_c,
     &                    gamma_cc, tau_c, mu, E, d1E, d2E)
            else
               call D2ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu,
     &                                  E, d1E, d2E)
            end if

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      end if
      end subroutine


C*****************************************************************************
      subroutine PBEC_nomu(rho, gamma_xx, diff_order, special_case,
     &                     E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized PBE functional
C This functional is just the normal PBE, no mu dependency.
C
C   Input : rho, vector of densities
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h" 
      integer, intent(in) :: diff_order, special_case
      real*8, intent(in) :: rho(2), gamma_xx(3)
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: spin_polarized
      real*8 :: rho_c, rho_s, gamma_cc
      
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      rho_s = rho(2)
      gamma_cc = gamma_xx(1)

      if (rho_c .le. 1.d-10) then
         return
      end if

      if (abs(rho_s) .le. 1.d-20) then
         spin_polarized = .false.
      else
         spin_polarized = .true.
      end if

      if (spin_polarized) then
         if (diff_order .eq. 0) then
            call ESRC_SPIN_PBE(rho_c, rho_s, gamma_cc, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_SPIN_PBE(rho_c, rho_s, gamma_cc, E, d1E)

         elseif (diff_order .eq. 2) then
            call D2ESRC_SPIN_PBE(rho_c, rho_s, gamma_cc, E, d1E, d2E)

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      else
         if (diff_order .eq. 0) then
            call ESRC_PBE(rho_c, gamma_cc, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_PBE(rho_c, gamma_cc, E, d1E)

         elseif (diff_order .eq. 2) then
            if (special_case .eq. 1) then
               call D2ESRC_PBE_singletref_triplet(rho_c,
     &                    gamma_cc, E, d1E, d2E)
            else
               call D2ESRC_PBE(rho_c, gamma_cc, E, d1E, d2E)
            end if

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      end if
      end subroutine


C*****************************************************************************
      subroutine wPBEX(rho, gamma_vec, diff_order, mu, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of Exchange wPBE functional
C
C   Input : rho, vector of densities
C           gamma_vec, is called gamma_vec, becuase "gamma" is an
C              intrinsic function and therefore a bad name.
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
!#include "pi.h" --" pi, pi2
#include "pi.h" 
#include "priunit.h" 
      integer, intent(in) :: diff_order
      real*8, intent(in) :: rho(4), gamma_vec(6), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: alpha_eq_beta
      real*8 :: Ea, Eb, d1Ea(4), d2Ea(10), d1Eb(4), d2Eb(10),
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


      if (diff_order .eq. 0) then
         call ESRX_wPBE(rho_a, gamma_aa, mu, Ea)
         if (.not.alpha_eq_beta) then
            call ESRX_wPBE(rho_b, gamma_bb, mu, Eb)
         end if

      elseif (diff_order .eq. 1) then
         call D1ESRX_wPBE(rho_a, gamma_aa, mu, Ea, d1Ea)
         if (.not.alpha_eq_beta) then
            call D1ESRX_wPBE(rho_b, gamma_bb, mu, Eb, d1Eb)
         end if

      elseif (diff_order .eq. 2) then
         call D2ESRX_wPBE(rho_a, gamma_aa, mu, Ea, d1Ea, d2Ea)
         if (.not.alpha_eq_beta) then
            call D2ESRX_wPBE(rho_b, gamma_bb, mu, Eb, d1Eb, d2Eb)
         end if

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if

C     Do the conversion from alpha, beta to charge, spin
      if (alpha_eq_beta) then
         Eb = Ea
         if (diff_order .ge. 1) then
            d1Eb(:) = d1Ea(:)
         end if
         if (diff_order .ge. 2) then
            d2Eb(:) = d2Ea(:)
         end if
      end if
      ! Spin splitting relation
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


C*****************************************************************************
      subroutine VWN5C_nomu(rho, diff_order, special_case, E, d1E, d2E)
C*****************************************************************************
C Implemented by E. R. Kjellgren
C
C Subroutine to to control calling of correlation spin-polarized VWN5 functional
C This functional is just the normal VWN5, no mu dependency.
C
C   Input : rho, vector of densities
C           diff_order, max order of differentation
C           mu, range-seperation parameter
C   Output : E, energy
C            d1E, vector of first derivatives of energy
C            d2E, vector of second derivatives of energy
C*****************************************************************************
      implicit none
#include "priunit.h" 
      integer, intent(in) :: diff_order, special_case
      real*8, intent(in) :: rho(2)
      real*8, intent(out) :: E, d1E(9), d2E(45)
      logical :: spin_polarized
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

      if (abs(rho_s) .le. 1.d-20) then
         spin_polarized = .false.
      else
         spin_polarized = .true.
      end if

      if (spin_polarized) then
         if (diff_order .eq. 0) then
            call ESRC_SPIN_VWN5(rho_c, rho_s, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_SPIN_VWN5(rho_c, rho_s, E, d1E)

         elseif (diff_order .eq. 2) then
            call D2ESRC_SPIN_VWN5(rho_c, rho_s, E, d1E, d2E)

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      else
         if (diff_order .eq. 0) then
            call ESRC_VWN5(rho_c, E)

         elseif (diff_order .eq. 1) then
            call D1ESRC_VWN5(rho_c, E, d1E)

         elseif (diff_order .eq. 2) then
            if (special_case .eq. 1) then
               call D2ESRC_VWN5_singletref_triplet(rho_c,
     &                           E, d1E, d2E)
            else
               call D2ESRC_VWN5(rho_c, E, d1E, d2E)
            end if

         else 
            call quit('Requested to high order derivative: ',diff_order)
         end if
      end if
      end subroutine
