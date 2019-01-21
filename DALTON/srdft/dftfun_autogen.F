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
!#include "pi.h" --> pi, pi2
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
!#include "pi.h" --> pi, pi2
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
!#include "pi.h" --> pi, pi2
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
         call D1ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E, d2E)

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if
      end subroutine


C****************************************************************************
      subroutine SPIN_TPSSC_GWS_ERF(rho, gamma_xx, tau, diff_order, mu,
     &                              E, d1E, d2E)
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
      real*8, intent(in) :: rho(2), gamma_xx(3), tau(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
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
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (diff_order .eq. 0) then
         write(lupri,*)'ESRC_SPIN_TPSS_GWS_ERF:' 
         write(lupri,*)'WARNING: Code not covered in ctest' 
         call ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, gamma_ss,
     &                               gamma_cs, tau_c, mu, E)

      elseif (diff_order .eq. 1) then
         call D1ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, gamma_ss,
     &                                 gamma_cs, tau_c, mu, E, d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, gamma_ss,
     &                                 gamma_cs, tau_c, mu, E, d1E, d2E)

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if
      end subroutine


C****************************************************************************
      subroutine TPSSC_GWS_ERF(rho, gamma_xx, tau, diff_order, mu, E,
     &                         d1E, d2E)
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
      real*8, intent(in) :: rho(2), gamma_xx(3), tau(2), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: rho_c, gamma_cc, tau_c
      
      write(lupri,*)'PBEC_GWS_ERF:' 
      write(lupri,*)'WARNING: Code not covered in ctest' 

      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      rho_c = rho(1)
      gamma_cc = gamma_xx(1)
      tau_c = tau(1)

      if (rho_c .le. 1.d-10) then
         ! Are we completely sure this check is valid
         !  for derivatives? - E. R. K
         return
      end if

      if (diff_order .eq. 0) then
         call ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E)

      elseif (diff_order .eq. 1) then
         call D1ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E, d1E)

      elseif (diff_order .eq. 2) then
         call D2ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E, d1E,
     &                           d2E)

      else 
         call quit('Requested to high order derivative: ',diff_order) 
      end if
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_LDA_ERF_case_1(rho_a, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5
      Ea = 0.0d0
      x0 = 2.14502939711103d0
      x1 = rho_a**0.333333333333333d0
      x2 = 1/(x0*x1)
      x3 = mu*x2
      x4 = 0.00844343197019482d0*1d0/rho_a
      x5 = mu**2
      Ea = 4.9628039271952d0*rho_a**1.33333333333333d0*(0.27516060407455
     &2d0*x3*(mu **3*x4 - mu*(-0.550321208149104d0*x2 + x4*x5)*exp( -15.
     &192666241151989d0*rho_a**0.66666666666666663d0/x5) - 0.82548181222
     &3657d0*x3 + 1.77245385090552d0*erf( 1.8171205928321397d0*x0*x1/mu)
     &) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_LDA_ERF_case_3(rho_a, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = 0.682784063255296d0
      Ea = -2.72568088924821d0*rho_a**1.33333333333333d0*x0
      d1Ea(1) = -3.63424118566428d0*rho_a**0.333333333333333d0*x0
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_LDA_ERF_case_2(rho_a, mu, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = 0.682784063255296d0
      x1 = 2.14502939711103d0
      x2 = 1d0/x1
      x3 = rho_a**0.333333333333333d0
      x4 = 1d0/x3
      x5 = x2*x4
      x6 = 0.550321208149104d0*x5
      x7 = mu**2
      x8 = 0.00844343197019481d0
      x9 = 1d0/rho_a*x8
      x10 = x7*x9
      x11 = x10 - x6
      x12 = 4.60115111447049d0
      x13 = exp(-3.3019272488946267d0*rho_a**0.66666666666666663d0*x12/x
     &7)
      x14 = mu*x13
      x15 = 1d0/mu
      x16 = x1*x15
      x17 = mu**3
      x18 = mu*x5
      x19 = x17*x9 - 0.825481812223657d0*x18 + 1.77245385090552d0*erf( 1
     &.8171205928321397d0*x16*x3)
      x20 = -x11*x14 + x19
      x21 = 0.275160604074552d0*x18
      x22 = rho_a**1.33333333333333d0
      x23 = 7.26848237132856d0*x22
      x24 = 1d0/x22
      x25 = x2*x24
      x26 = rho_a**(-2.0d0)*x8
      Ea = x0*x23*(x20*x21 - 0.375d0)
      d1Ea(1) = -0.318309886183791d0*mu*x23*(0.0917202013581841d0*x20*x2
     &4 - 0.275160604074552d0*x4*(0.275160604074552d0*mu*x25 + 1.2114137
     &2855476d0*rho_a**(-0.666666666666667d0)*x13*x16 + 2.20128483259642
     &d0*rho_a**(-0.333333333333333d0)*x11*x12*x13*x15 + x14*(-0.1834404
     &02716368d0*x25 + x26*x7) - x17*x26)) + 9.69130982843808d0*rho_a**0
     &.333333333333333d0*x0*(x21*(x14*(-x10 + x6) + x19) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = 3.14159265358979d0/mu**2
      Ea = -1.0d0*rho_a**2.0d0*x0
      d1Ea(1) = -2.0d0*rho_a**1.0d0*x0
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_LDA_ERF_case_1(rho_a, mu, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = 0.682784063255296d0
      Ea = -2.72568088924821d0*rho_a**1.33333333333333d0*x0
      d1Ea(1) = -3.63424118566428d0*rho_a**0.333333333333333d0*x0
      d2Ea(1) = -1.21141372855476d0*rho_a**(-0.666666666666667d0)*x0
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_LDA_ERF_case_2(rho_a, mu, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = 0.682784063255296d0
      x1 = 2.14502939711103d0
      x2 = 1d0/x1
      x3 = rho_a**0.333333333333333d0
      x4 = 1d0/x3
      x5 = x2*x4
      x6 = 0.550321208149104d0*x5
      x7 = mu**2
      x8 = 1d0/rho_a
      x9 = 9.86960440108936d0
      x10 = 1d0/x9
      x11 = 0.0833333333333333d0*x10
      x12 = x11*x8
      x13 = x12*x7
      x14 = x13 - x6
      x15 = 4.60115111447049d0
      x16 = exp(-3.3019272488946267d0*rho_a**0.66666666666666663d0*x15/x
     &7)
      x17 = mu*x16
      x18 = 1d0/mu
      x19 = x1*x18
      x20 = mu**3
      x21 = mu*x5
      x22 = x12*x20 - 0.825481812223657d0*x21 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x19*x3)
      x23 = -x14*x17 + x22
      x24 = 0.275160604074552d0*x21
      x25 = rho_a**1.33333333333333d0
      x26 = 7.26848237132856d0*x25
      x27 = rho_a**0.333333333333333d0
      x28 = -x13 + x6
      x29 = x17*x28 + x22
      x30 = x0*(x24*x29 - 0.375d0)
      x31 = 1d0/x25
      x32 = x15*x16*x18
      x33 = rho_a**(-0.333333333333333d0)*x32
      x34 = 2.20128483259642d0*x33
      x35 = 0.183440402716368d0*x31
      x36 = rho_a**(-2.0d0)*x11
      x37 = -x2*x35 + x36*x7
      x38 = rho_a**(-0.666666666666667d0)
      x39 = x16*x19
      x40 = 0.275160604074552d0*mu*x2*x31 + x17*x37 - x20*x36 + 1.211413
     &72855476d0* x38*x39
      x41 = 0.275160604074552d0*x4
      x42 = 0.0917202013581841d0*x23*x31 - x41*(x14*x34 + x40)
      x43 = 0.318309886183791d0
      x44 = mu*x26*x43
      x45 = rho_a**(-2.33333333333333d0)
      x46 = x2*x45
      x47 = 0.166666666666667d0*rho_a**(-3.0d0)*x10
      x48 = x16/x20
      Ea = x0*x26*(x23*x24 - 0.375d0)
      d1Ea(1) = 9.69130982843808d0*x27*x30 - x42*x44
      d2Ea(1) = -19.3826196568762d0*mu*x27*x42*x43 + 3.23043660947936d0*
     &x30*x38 + x44*( 0.122293601810912d0*x29*x45 - x35*(-x28*x34 + x40)
     & + x41*( -0.366880805432736d0*mu*x46 - 0.80760915236984d0*rho_a**(
     & -1.66666666666667d0)*x39 + 0.733761610865473d0*rho_a**( -1.333333
     &33333333d0)*x28*x32 + x17*(0.244587203621824d0*x46 - x47 *x7) + x2
     &0*x47 + 102.585381117795d0*x28*x38*x48 - 4.40256966519284d0*x33*x3
     &7 - 2.66666666666667d0*x48*x8*x9))
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_LDA_ERF_case_3(rho_a, mu, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = 3.14159265358979d0/mu**2
      x1 = 2.0d0*x0
      Ea = -1.0d0*rho_a**2.0d0*x0
      d1Ea(1) = -rho_a**1.0d0*x1
      d2Ea(1) = -x1
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.217336917462899d0
      x2 = mu**2*rho_a**(-0.666666666666667d0)*x1
      Ea = -1.8610514726982d0*x0*(1.804d0 - 0.804d0/(0.239205142311446d0
     &*gamma_aa* rho_a**(-2.66666666666667d0)*x1*(5.45136177849642d0*x2 
     &- 7.0d0)* exp(-1.4385538026587772d0*x2)/(0.233951463049512d0*mu**4
     &/x0 + 49.1113297378199d0*mu*rho_a**(-0.333333333333333d0) - 49.062
     &2560064678d0*x2 - 81.0d0) + 1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 4.60115111447049d0
      x2 = 1d0/x1
      x3 = mu**2
      x4 = rho_a**0.666666666666667d0
      x5 = x2*x3/x4
      x6 = rho_a**0.333333333333333d0
      x7 = 1d0/x6
      x8 = 2.14502939711103d0
      x9 = x7/x8
      x10 = mu*x9
      x11 = 0.00844343197019482d0*1d0/rho_a
      Ea = 4.9628039271952d0*x0*(1.804d0 - 0.804d0/(0.239205142311446d0*
     &gamma_aa* rho_a**(-2.66666666666667d0)*x2*(5.45136177849642d0*x5 -
     & 7.0d0)* exp(-1.4385538026587772d0*x5)/(0.233951463049512d0*mu**4/
     &x0 + 49.1113297378199d0*mu*x7 - 49.0622560064678d0*x5 - 81.0d0) + 
     &1.0d0))*(0.275160604074552d0*x10*(mu**3*x11 - mu*(x11*x3 - 0.55032
     &1208149104d0*x9)*exp(-3.3019272488946267d0*x1*x4/x3) - 0.825481812
     &223657d0*x10 + 1.77245385090552d0*erf( 1.8171205928321397d0*x6*x8/
     &mu)) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 4.60115111447049d0
      x2 = 1d0/x1
      x3 = mu**2
      x4 = 1d0/x3
      x5 = rho_a**0.666666666666667d0
      x6 = x1*x5
      x7 = 2.14502939711103d0
      x8 = rho_a**0.333333333333333d0
      x9 = 1/(x7*x8)
      x10 = mu*x9
      x11 = 0.00844343197019482d0*1d0/rho_a
      Ea = 4.9628039271952d0*x0*(1.804d0 + 0.804d0/(0.239205142311446d0*
     &gamma_aa* rho_a**(-2.66666666666667d0)*x2*x4*(0.213719543995407d0*
     &x0*x4 - 0.183440402716368d0*x6 + 14.4026298885989d0*rho_a**2.0d0/m
     &u**4)* exp(-1.4385538026587772d0*x2*x3/x5) - 1.0d0))*( 0.275160604
     &074552d0*x10*(mu**3*x11 - mu*(x11*x3 - 0.550321208149104d0*x9)*exp
     &(-3.3019272488946267d0*x4*x6) - 0.825481812223657d0*x10 + 1.772453
     &85090552d0*erf( 1.8171205928321397d0*x7*x8/mu)) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 2.14502939711103d0
      x2 = rho_a**0.333333333333333d0
      x3 = 1/(x1*x2)
      x4 = mu*x3
      x5 = 0.00844343197019482d0*1d0/rho_a
      x6 = erf(1.8171205928321397d0*x1*x2/mu)
      x7 = 4.60115111447049d0
      x8 = mu**2
      x9 = 1d0/x8
      x10 = rho_a**0.666666666666667d0
      x11 = 3.30192724889463d0*x10*x7*x9
      x12 = x8/(x10*x7)
      x13 = exp(x11)
      x14 = 0.0472353356922751d0*mu**4*(x13 - 1.0d0)/x0
      x15 = 3.0d0*x13
      Ea = 4.9628039271952d0*x0*(1.804d0 + 0.804d0/(3.15935190989552d0*g
     &amma_aa* rho_a**(-2.0d0)*x9*(0.00280419844572861d0*x12*(7.0d0*x13 
     &+ 11.0d0 ) - 0.015286700226364d0*x14 + 0.0185185185185185d0)*exp( 
     &-1.4385538026587772d0*x12)/(-0.60570686427738d0*x12*(x15 - 2.0d0 )
     & + 3.90167577847542d0*x13*x4*x6 + 0.183440402716368d0*x14 - x15 ) 
     &- 1.0d0))*(0.275160604074552d0*x4*(mu**3*x5 - mu*( -0.550321208149
     &104d0*x3 + x5*x8)*exp(-x11) - 0.825481812223657d0* x4 + 1.77245385
     &090552d0*x6) - 0.375d0)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5
      Ea = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = rho_a**2.0d0
      x3 = 4.60115111447049d0
      x4 = 1d0/x3
      x5 = rho_a**0.666666666666667d0
      Ea = -3.14159265358979d0*x1*x2*(1.804d0 + 0.804d0/(0.2392051423114
     &46d0* gamma_aa*rho_a**(-2.66666666666667d0)*x1*x4*(0.2137195439954
     &07d0* rho_a**1.33333333333333d0*x1 - 0.183440402716368d0*x3*x5 + 1
     &4.4026298885989d0*x2/mu**4)*exp(-1.4385538026587772d0*x0*x4/x5) - 
     &1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea, d1E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = 0.682784063255296d0
      x1 = rho_a**1.33333333333333d0
      x2 = 0.826307487110758d0
      x3 = 0.0472353356922751d0
      x4 = 1d0/x1
      x5 = 0.217336917462899d0
      x6 = mu**2
      x7 = rho_a**(-0.666666666666667d0)*x5*x6
      x8 = 1d0/(4.95289087334194d0*mu**4*x3*x4 + 59.4346904801033d0*mu*r
     &ho_a**( -0.333333333333333d0)*x2 - 49.0622560064678d0*x7 - 81.0d0)
     &
      x9 = x8*exp(-1.4385538026587772d0*x7)
      x10 = gamma_aa*x9
      x11 = 5.45136177849642d0*x7 - 7.0d0
      x12 = x11*x5
      x13 = 0.239205142311446d0*rho_a**(-2.66666666666667d0)*x12
      x14 = x10*x13 + 1.0d0
      x15 = 1.804d0 - 0.804d0/x14
      x16 = x14**(-2)
      x17 = rho_a**(-4.33333333333333d0)*x3*x6
      Ea = -2.72568088924821d0*x0*x1*x15
      d1Ea(1) = -x0*(3.63424118566428d0*rho_a**0.333333333333333d0*x15 +
     & 2.19144743495556d0*x1*x10*x16*(mu*x13*x8*(6.60385449778925d0*mu**
     & 3*rho_a**(-2.33333333333333d0)*x3 - 32.7081706709785d0*mu*rho_a**
     & (-1.66666666666667d0)*x5 + 19.8115634933678d0*x2*x4) - 0.63788037
     &9497189d0*rho_a**(-3.66666666666667d0)*x12 + 0.229406311391776d0*x
     &11*x17 - 0.869329180010942d0*x17))
      d1Ea(2) = -0.0777890465506374d0*x11*x16*x4*x9
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea, d
     &1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = 2.14502939711103d0
      x1 = 1d0/x0
      x2 = rho_a**0.333333333333333d0
      x3 = 1d0/x2
      x4 = x1*x3
      x5 = 0.550321208149104d0*x4
      x6 = mu**2
      x7 = 0.00844343197019481d0
      x8 = 1d0/rho_a*x7
      x9 = x6*x8
      x10 = -x5 + x9
      x11 = 4.60115111447049d0
      x12 = rho_a**0.666666666666667d0
      x13 = exp(-3.3019272488946267d0*x11*x12/x6)
      x14 = mu*x13
      x15 = 1d0/mu
      x16 = x0*x15
      x17 = mu**3
      x18 = mu*x4
      x19 = x17*x8 - 0.825481812223657d0*x18 + 1.77245385090552d0*erf( 1
     &.8171205928321397d0*x16*x2)
      x20 = -x10*x14 + x19
      x21 = 0.275160604074552d0*x18
      x22 = x20*x21 - 0.375d0
      x23 = 0.682784063255296d0
      x24 = 0.826307487110758d0
      x25 = 0.0472353356922751d0
      x26 = rho_a**1.33333333333333d0
      x27 = 1d0/x26
      x28 = 1d0/x11
      x29 = x28*x6/x12
      x30 = 1d0/(4.95289087334194d0*mu**4*x25*x27 + 59.4346904801033d0*m
     &u*x24*x3 - 49.0622560064678d0*x29 - 81.0d0)
      x31 = x30*exp(-1.4385538026587772d0*x29)
      x32 = gamma_aa*x31
      x33 = 5.45136177849642d0*x29 - 7.0d0
      x34 = x28*x33
      x35 = 0.239205142311446d0*rho_a**(-2.66666666666667d0)*x34
      x36 = x32*x35 + 1.0d0
      x37 = 1.804d0 - 0.804d0/x36
      x38 = x23*x37
      x39 = 7.26848237132856d0*x26
      x40 = x21*(x14*(x5 - x9) + x19) - 0.375d0
      x41 = x1*x27
      x42 = rho_a**(-2.0d0)*x7
      x43 = x36**(-2)
      x44 = rho_a**(-4.33333333333333d0)*x25*x6
      Ea = x22*x38*x39
      d1Ea(1) = -0.318309886183791d0*mu*x37*x39*(0.0917202013581841d0*x2
     &0*x27 - 0.275160604074552d0*x3*(0.275160604074552d0*mu*x41 + 1.211
     &41372855476d0*rho_a**(-0.666666666666667d0)*x13*x16 + 2.2012848325
     &9642d0*rho_a**(-0.333333333333333d0)*x10*x11*x13*x15 + x14*(-0.183
     &440402716368d0*x41 + x42*x6) - x17*x42)) + 9.69130982843808d0*rho_
     &a**0.333333333333333d0*x38*x40 + 5.84385982654816d0*x23*x26*x32*x4
     &0*x43*(mu*x30*x35*( -32.7081706709785d0*mu*rho_a**(-1.666666666666
     &67d0)*x28 + 6.60385449778925d0*rho_a**(-2.33333333333333d0)*x17*x2
     &5 + 19.8115634933678d0*x24*x27) - 0.637880379497189d0*rho_a**( -3.
     &66666666666667d0)*x34 + 0.229406311391776d0*x33*x44 - 0.8693291800
     &10942d0*x44)
      d1Ea(2) = 0.207437457468366d0*x22*x27*x31*x33*x43
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea, d
     &1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = 0.682784063255296d0
      x1 = 4.60115111447049d0
      x2 = rho_a**0.666666666666667d0
      x3 = x1*x2
      x4 = 0.183440402716368d0*x3
      x5 = rho_a**1.33333333333333d0
      x6 = 21.1705915781930d0
      x7 = mu**2
      x8 = 1d0/x7
      x9 = x6*x8
      x10 = 0.010095114404623d0*x5*x9
      x11 = rho_a**2.0d0
      x12 = 97.4090910340024d0/mu**4
      x13 = 0.147857142857143d0*x11*x12
      x14 = x10 + x13 - x4
      x15 = 1d0/x1
      x16 = x15*x8
      x17 = 0.239205142311446d0*rho_a**(-2.66666666666667d0)*x16
      x18 = exp(-1.4385538026587772d0*x15*x7/x2)
      x19 = gamma_aa*x18
      x20 = x14*x17*x19
      x21 = x20 - 1.0d0
      x22 = 2.14502939711103d0
      x23 = 1d0/x22
      x24 = rho_a**0.333333333333333d0
      x25 = 1d0/x24
      x26 = x23*x25
      x27 = 0.550321208149104d0*x26
      x28 = rho_a**1.0d0
      x29 = 0.00844343197019481d0
      x30 = x29/x28
      x31 = x30*x7
      x32 = -x27 + x31
      x33 = exp(-3.3019272488946267d0*x3*x8)
      x34 = mu*x33
      x35 = 1d0/mu
      x36 = x22*x35
      x37 = mu**3
      x38 = mu*x26
      x39 = x30*x37 - 0.825481812223657d0*x38 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x24*x36)
      x40 = -x32*x34 + x39
      x41 = 0.275160604074552d0*x38
      x42 = x40*x41 - 0.375d0
      x43 = 7.26848237132856d0*x5
      x44 = rho_a**0.333333333333333d0
      x45 = -x20 + 1.0d0
      x46 = 1.804d0 - 0.804d0/x45
      x47 = x0*(x41*(x34*(x27 - x31) + x39) - 0.375d0)
      x48 = 1d0/x5
      x49 = x23*x48
      x50 = x29/x11
      x51 = rho_a**(-0.333333333333333d0)*x1
      x52 = -x10 - x13 + x4
      Ea = x0*x42*x43*(1.804d0 + 0.804d0/x21)
      d1Ea(1) = -0.318309886183791d0*mu*x43*x46*(-0.275160604074552d0*x2
     &5*( 0.275160604074552d0*mu*x49 + 1.21141372855476d0*rho_a**( -0.66
     &6666666666667d0)*x33*x36 + 2.20128483259642d0*x32*x33*x35* x51 + x
     &34*(-0.183440402716368d0*x49 + x50*x7) - x37*x50) + 0.091720201358
     &1841d0*x40*x48) + 5.84385982654816d0*x19*x47*x5*( 0.22940631139177
     &6d0*rho_a**(-4.33333333333333d0)*x52/x6 - 0.637880379497189d0*rho_
     &a**(-3.66666666666667d0)*x16*x52 + x17*( -0.295714285714286d0*x12*
     &x28 - 0.0134601525394973d0*x44*x9 + 0.122293601810912d0*x51))/x45*
     &*2 + 9.69130982843808d0*x44*x46*x47
      d1Ea(2) = -0.207437457468366d0*x14*x18*x42*x48*x8/x21**2
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea, d
     &1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 1d0/x0
      x2 = rho_a**0.666666666666667d0
      x3 = 4.60115111447049d0
      x4 = mu**2
      x5 = 1d0/x4
      x6 = x3*x5
      x7 = 3.30192724889463d0*x2*x6
      x8 = exp(x7)
      x9 = 0.0472353356922751d0*mu**4*(x8 - 1.0d0)
      x10 = x1*x9
      x11 = 0.015286700226364d0*x10
      x12 = 7.0d0*x8 + 11.0d0
      x13 = 1d0/x3
      x14 = x13*x4
      x15 = x14/x2
      x16 = 0.00280419844572861d0*x12*x15
      x17 = -x11 + x16 + 0.0185185185185185d0
      x18 = 3.0d0*x8
      x19 = x18 - 2.0d0
      x20 = 2.14502939711103d0
      x21 = 1d0/x20
      x22 = rho_a**0.333333333333333d0
      x23 = 1d0/x22
      x24 = x21*x23
      x25 = mu*x24
      x26 = 1d0/mu
      x27 = x20*x26
      x28 = erf(1.8171205928321397d0*x22*x27)
      x29 = x28*x8
      x30 = 1d0/(0.183440402716368d0*x10 - 0.60570686427738d0*x15*x19 - 
     &x18 + 3.90167577847542d0*x25*x29)
      x31 = rho_a**(-2.0d0)
      x32 = 3.15935190989552d0*x31*x5
      x33 = x30*x32
      x34 = exp(-1.4385538026587772d0*x15)
      x35 = gamma_aa*x34
      x36 = x17*x33*x35
      x37 = x36 - 1.0d0
      x38 = 7.26848237132856d0*x0
      x39 = 0.682784063255296d0
      x40 = 0.550321208149104d0*x24
      x41 = 1d0/rho_a
      x42 = 0.00844343197019481d0
      x43 = x41*x42
      x44 = x4*x43
      x45 = -x40 + x44
      x46 = exp(-x7)
      x47 = mu*x46
      x48 = mu**3
      x49 = -0.825481812223657d0*x25 + 1.77245385090552d0*x28 + x43*x48
      x50 = -x45*x47 + x49
      x51 = 0.275160604074552d0*x25
      x52 = x39*(x50*x51 - 0.375d0)
      x53 = -x36 + 1.0d0
      x54 = 1.804d0 - 0.804d0/x53
      x55 = x39*(x51*(x47*(x40 - x44) + x49) - 0.375d0)
      x56 = x1*x21
      x57 = mu*x56
      x58 = x31*x42
      x59 = rho_a**(-0.666666666666667d0)
      x60 = x27*x59
      x61 = rho_a**(-0.333333333333333d0)
      x62 = x11 - x16 - 0.0185185185185185d0
      x63 = x41*x8
      x64 = rho_a**(-1.66666666666667d0)*x14
      x65 = x64*x8
      x66 = rho_a**(-2.33333333333333d0)*x9
      Ea = x38*x52*(1.804d0 + 0.804d0/x37)
      d1Ea(1) = -0.318309886183791d0*mu*x38*x54*(0.0917202013581841d0*x1
     &*x50 - 0.275160604074552d0*x23*(2.20128483259642d0*x26*x3*x45*x46*
     &x61 + 1.21141372855476d0*x46*x60 + x47*(x4*x58 - 0.183440402716368
     &d0* x56) - x48*x58 + 0.275160604074552d0*x57)) + 9.69130982843808d
     &0* rho_a**0.333333333333333d0*x54*x55 + 5.84385982654816d0*x0*x30*
     & x35*x55*(3.02993180261165d0*rho_a**(-3.66666666666667d0)*x13*x62 
     &- 6.31870381979104d0*rho_a**(-3.0d0)*x5*x62 + x32*( 0.001869465630
     &48574d0*x12*x64 - 0.0432098765432099d0*x63 + 0.0336503813487433d0*
     &x65 - 0.0203822669684853d0*x66) + x33*x62*( -0.40380457618492d0*x1
     &9*x64 + 1.30055859282514d0*x29*x57 - 8.58869971286677d0*x29*x60 - 
     &2.66666666666667d0*x41 + 6.60385449778925d0*x6*x61*x8 + 4.0d0*x63 
     &- 0.40380457618492d0*x65 + 0.244587203621824d0*x66))/x53**2
      d1Ea(2) = -18.4628097041666d0*x17*x30*x34*x5*x52*x59/x37**2
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea, d1E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**2.0d0
      x1 = 4.60115111447049d0
      x2 = rho_a**0.666666666666667d0
      x3 = 0.183440402716368d0*x1*x2
      x4 = mu**2
      x5 = 1d0/x4
      x6 = 21.1705915781930d0
      x7 = x5*x6
      x8 = 0.010095114404623d0*rho_a**1.33333333333333d0*x7
      x9 = mu**(-4)
      x10 = 97.4090910340024d0*x9
      x11 = 0.147857142857143d0*x0*x10
      x12 = x11 - x3 + x8
      x13 = 1d0/x1
      x14 = x13*x5
      x15 = 0.239205142311446d0*rho_a**(-2.66666666666667d0)*x14
      x16 = exp(-1.4385538026587772d0*x13*x4/x2)
      x17 = gamma_aa*x16
      x18 = x12*x15*x17
      x19 = x18 - 1.0d0
      x20 = 3.14159265358979d0*x5
      x21 = rho_a**1.0d0
      x22 = -x18 + 1.0d0
      x23 = -x11 + x3 - x8
      Ea = -1.0d0*x0*x20*(1.804d0 + 0.804d0/x19)
      d1Ea(1) = -x20*(0.804d0*x0*x17*(0.229406311391776d0*rho_a**(-4.333
     &33333333333d0)* x23/x6 - 0.637880379497189d0*rho_a**(-3.6666666666
     &6667d0)*x14*x23 + x15*(0.122293601810912d0*rho_a**(-0.333333333333
     &333d0)*x1 - 0.0134601525394973d0*rho_a**0.333333333333333d0*x7 - 0
     &.295714285714286d0*x10*x21))/x22**2 + 2.0d0*x21*(1.804d0 - 0.804d0
     &/x22))
      d1Ea(2) = 0.131313669051252d0*rho_a**(-0.666666666666667d0)*x12*x1
     &6*x9/x19**2
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_1(rho_a, gamma_aa, mu, Ea, d1E
     &a, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = 0.217336917462899d0
      x1 = mu**2
      x2 = x0*x1
      x3 = rho_a**(-0.666666666666667d0)*x2
      x4 = exp(-1.4385538026587772d0*x3)
      x5 = gamma_aa*x4
      x6 = 0.826307487110758d0
      x7 = 0.0472353356922751d0
      x8 = mu**4
      x9 = rho_a**1.33333333333333d0
      x10 = 1d0/x9
      x11 = 59.4346904801033d0*mu*rho_a**(-0.333333333333333d0)*x6 + 4.9
     &5289087334194d0*x10*x7*x8 - 49.0622560064678d0*x3 - 81.0d0
      x12 = 1d0/x11
      x13 = rho_a**(-2.66666666666667d0)
      x14 = 5.45136177849642d0*x3 - 7.0d0
      x15 = x0*x14
      x16 = x13*x15
      x17 = 0.239205142311446d0*x16
      x18 = x12*x17
      x19 = x18*x5 + 1.0d0
      x20 = 1d0/x19
      x21 = -0.804d0*x20 + 1.804d0
      x22 = 0.682784063255296d0
      x23 = x22*x9
      x24 = rho_a**0.333333333333333d0
      x25 = rho_a**(-4.33333333333333d0)
      x26 = x1*x7
      x27 = x25*x26
      x28 = rho_a**(-3.66666666666667d0)*x15
      x29 = x14*x27
      x30 = x10*x6
      x31 = mu*x0
      x32 = rho_a**(-1.66666666666667d0)*x31
      x33 = rho_a**(-2.33333333333333d0)
      x34 = mu**3*x7
      x35 = x33*x34
      x36 = 19.8115634933678d0*x30 - 32.7081706709785d0*x32 + 6.60385449
     &778925d0*x35
      x37 = x12*x36
      x38 = mu*x37
      x39 = x17*x38 - 0.869329180010942d0*x27 - 0.637880379497189d0*x28 
     &+ 0.229406311391776d0*x29
      x40 = x39*x4
      x41 = 2.19144743495556d0*x40
      x42 = x12/x19**2
      x43 = gamma_aa*x42
      x44 = 0.148394183600699d0*x14
      x45 = x4*x44
      x46 = 1.27576075899438d0*x28
      x47 = 2.19144743495556d0*x9
      x48 = x11**(-2)
      x49 = exp(-2.8771076053175544d0*x3)
      x50 = x48*x49/x19**3
      x51 = 0.0102659822546843d0*rho_a**(-6.0d0)*x8
      x52 = rho_a**(-5.33333333333333d0)*x26
      x53 = x25*x34*x37
      Ea = -2.72568088924821d0*x21*x23
      d1Ea(1) = -x22*(3.63424118566428d0*x21*x24 + x41*x43*x9)
      d1Ea(2) = -0.524205495546598d0*x10*x42*x45
      d2Ea(1) = x22*(gamma_aa**2*x39*x47*x50*(0.478410284622892d0*x16*x3
     &8 - 1.73865836002188d0*x27 + 0.458812622783553d0*x29 - x46) - 1.21
     &141372855476d0*rho_a**(-0.666666666666667d0)*x21 - 5.8438598265481
     &6d0*x24*x40*x43 - x42*x47*x5*(mu*x18*( -15.4089938281749d0*rho_a**
     &(-3.33333333333333d0)*x34 + 54.5136177849642d0*x13*x31 - 26.415417
     &991157d0*x33*x6) + 2.33889472482303d0*rho_a**(-4.66666666666667d0)
     &*x15 + 0.239205142311446d0*x13*x14*x2*x36*x48*(39.6231269867355d0*
     &x30 - 65.416341341957d0*x32 + 13.2077089955785d0*x35) + 0.22000888
     &1071042d0*x14*x51 - 1.60584417974243d0*x14*x52 + 0.458812622783553
     &d0*x14*x53 - x38*x46 - 1.66743573022264d0*x51 + 6.08530426007659d0
     &*x52 - 1.73865836002188d0*x53))
      d2Ea(2) = x42*(1.0484109910932d0*gamma_aa*x10*x12*x20*x39*x44*x49 
     &- x23*x41 - 0.698940660728797d0*x33*x45)
      d2Ea(3) = 0.00808821074878272d0*rho_a**(-4.0d0)*x14**2*x50
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_2_1(rho_a, gamma_aa, mu, Ea, d
     &1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &2, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = 2.14502939711103d0
      x1 = 1d0/x0
      x2 = rho_a**0.333333333333333d0
      x3 = 1d0/x2
      x4 = x1*x3
      x5 = 0.550321208149104d0*x4
      x6 = mu**2
      x7 = 1d0/rho_a
      x8 = 9.86960440108936d0
      x9 = 1d0/x8
      x10 = 0.0833333333333333d0*x9
      x11 = x10*x7
      x12 = x11*x6
      x13 = x12 - x5
      x14 = 4.60115111447049d0
      x15 = rho_a**0.666666666666667d0
      x16 = exp(-3.3019272488946267d0*x14*x15/x6)
      x17 = mu*x16
      x18 = 1d0/mu
      x19 = x0*x18
      x20 = mu**3
      x21 = mu*x4
      x22 = x11*x20 - 0.825481812223657d0*x21 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x19*x2)
      x23 = -x13*x17 + x22
      x24 = 0.275160604074552d0*x21
      x25 = x23*x24 - 0.375d0
      x26 = 0.682784063255296d0
      x27 = 1d0/x14
      x28 = x27*x6
      x29 = x28/x15
      x30 = exp(-1.4385538026587772d0*x29)
      x31 = gamma_aa*x30
      x32 = 0.826307487110758d0
      x33 = 21.1705915781930d0
      x34 = 1d0/x33
      x35 = mu**4
      x36 = rho_a**1.33333333333333d0
      x37 = 1d0/x36
      x38 = 59.4346904801033d0*mu*x3*x32 - 49.0622560064678d0*x29 + 4.95
     &289087334194d0*x34*x35*x37 - 81.0d0
      x39 = 1d0/x38
      x40 = rho_a**(-2.66666666666667d0)
      x41 = 5.45136177849642d0*x29 - 7.0d0
      x42 = x27*x41
      x43 = x40*x42
      x44 = 0.239205142311446d0*x43
      x45 = x39*x44
      x46 = x31*x45 + 1.0d0
      x47 = 1d0/x46
      x48 = -0.804d0*x47 + 1.804d0
      x49 = x26*x48
      x50 = x25*x49
      x51 = 7.26848237132856d0*x36
      x52 = rho_a**0.333333333333333d0
      x53 = x24*(x17*(-x12 + x5) + x22) - 0.375d0
      x54 = 0.318309886183791d0
      x55 = mu*x37
      x56 = rho_a**(-2.0d0)*x10
      x57 = rho_a**(-0.666666666666667d0)
      x58 = x16*x19
      x59 = 0.183440402716368d0*x37
      x60 = -x1*x59 + x56*x6
      x61 = rho_a**(-0.333333333333333d0)*x14*x16*x18
      x62 = 0.275160604074552d0*x1*x55 + 2.20128483259642d0*x13*x61 + x1
     &7*x60 - x20* x56 + 1.21141372855476d0*x57*x58
      x63 = 0.275160604074552d0*x3
      x64 = 0.0917202013581841d0*x23*x37 - x62*x63
      x65 = mu*x54*x64
      x66 = x48*x51
      x67 = x39/x46**2
      x68 = x30*x53
      x69 = rho_a**(-4.33333333333333d0)
      x70 = x34*x6
      x71 = x69*x70
      x72 = rho_a**(-3.66666666666667d0)*x42
      x73 = x41*x71
      x74 = x32*x37
      x75 = mu*x27
      x76 = rho_a**(-1.66666666666667d0)*x75
      x77 = rho_a**(-2.33333333333333d0)
      x78 = x20*x34
      x79 = x77*x78
      x80 = 19.8115634933678d0*x74 - 32.7081706709785d0*x76 + 6.60385449
     &778925d0*x79
      x81 = x39*x80
      x82 = mu*x81
      x83 = x44*x82 - 0.869329180010942d0*x71 - 0.637880379497189d0*x72 
     &+ 0.229406311391776d0*x73
      x84 = x26*x83
      x85 = 5.84385982654816d0*x36
      x86 = x84*x85
      x87 = x68*x86
      x88 = 1.39788132145759d0*x30
      x89 = x25*x67
      x90 = 0.148394183600699d0*x41
      x91 = x37*x90
      x92 = x31*x89
      x93 = 1.27576075899438d0*x72
      x94 = x38**(-2)
      x95 = exp(-2.8771076053175544d0*x29)
      x96 = x25*x94*x95/x46**3
      x97 = x1*x77
      x98 = 0.166666666666667d0*rho_a**(-3.0d0)*x9
      x99 = x16/x20
      x100 = 0.0102659822546843d0*rho_a**(-6.0d0)*x35
      x101 = rho_a**(-5.33333333333333d0)*x70
      x102 = x69*x78*x81
      Ea = x50*x51
      d1Ea(1) = gamma_aa*x67*x87 + 9.69130982843808d0*x49*x52*x53 - x65*
     &x66
      d1Ea(2) = x88*x89*x91
      d2Ea(1) = -gamma_aa**2*x86*x96*(0.478410284622892d0*x43*x82 - 1.73
     &865836002188d0* x71 + 0.458812622783553d0*x73 - x93) + mu*x54*x66*
     &( 0.122293601810912d0*x23*x77 - x59*x62 + x63*(-0.366880805432736d
     &0 *mu*x97 - 0.80760915236984d0*rho_a**(-1.66666666666667d0)*x58 - 
     &0.733761610865473d0*rho_a**(-1.33333333333333d0)*x13*x14*x16*x18 -
     & 4.84565491421904d0*x13*x33*x57*x99 - x17*(x6*x98 - 0.244587203621
     &824d0*x97) + x20*x98 - 4.40256966519284d0*x60*x61 - 2.666666666666
     &67d0*x7*x8*x99)) + x26*x85*x92*(-mu*x45*( 15.4089938281749d0*rho_a
     &**(-3.33333333333333d0)*x78 + 26.415417991157d0*x32*x77 - 54.51361
     &77849642d0*x40*x75) + 2.33889472482303d0*rho_a**(-4.66666666666667
     &d0)*x42 + 0.220008881071042d0*x100*x41 - 1.66743573022264d0*x100 -
     & 1.60584417974243d0*x101*x41 + 6.08530426007659d0*x101 + 0.4588126
     &22783553d0*x102*x41 - 1.73865836002188d0*x102 + 0.239205142311446d
     &0*x28*x40*x41*x80*x94*(39.6231269867355d0*x74 - 65.416341341957d0*
     &x76 + 13.2077089955785d0*x79) - x82*x93) - 11.6877196530963d0*x31*
     &x36*x65*x67*x83 - 19.3826196568762d0*x48* x52*x65 + 3.230436609479
     &36d0*x50*x57 + 15.5836262041284d0*x52*x84 *x92
      d2Ea(2) = x67*(-2.79576264291519d0*gamma_aa*x39*x47*x53*x83*x91*x9
     &5 - 0.0691804894611514d0*x41*x55*x64*x88 + 1.86384176194346d0*x68*
     &x77 *x90 + x87)
      d2Ea(3) = -0.0215685619967539d0*rho_a**(-4.0d0)*x41**2*x96
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_2_2(rho_a, gamma_aa, mu, Ea, d
     &1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x105, x106
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 7.26848237132856d0*x0
      x2 = mu**2
      x3 = 1d0/x2
      x4 = 4.60115111447049d0
      x5 = 1d0/x4
      x6 = x3*x5
      x7 = rho_a**(-2.66666666666667d0)*x6
      x8 = 0.239205142311446d0*x7
      x9 = rho_a**0.666666666666667d0
      x10 = x4*x9
      x11 = 0.183440402716368d0*x10
      x12 = 21.1705915781930d0
      x13 = x12*x3
      x14 = 0.010095114404623d0*x0*x13
      x15 = rho_a**2.0d0
      x16 = 97.4090910340024d0
      x17 = mu**(-4)
      x18 = x16*x17
      x19 = 0.147857142857143d0*x15*x18
      x20 = -x11 + x14 + x19
      x21 = x2*x5/x9
      x22 = exp(-1.4385538026587772d0*x21)
      x23 = x20*x22
      x24 = gamma_aa*x23*x8
      x25 = x24 - 1.0d0
      x26 = 1.804d0 + 0.804d0/x25
      x27 = 0.682784063255296d0
      x28 = 2.14502939711103d0
      x29 = 1d0/x28
      x30 = rho_a**0.333333333333333d0
      x31 = 1d0/x30
      x32 = x29*x31
      x33 = 0.550321208149104d0*x32
      x34 = rho_a**1.0d0
      x35 = 1d0/x34
      x36 = 9.86960440108936d0
      x37 = 1d0/x36
      x38 = 0.0833333333333333d0*x37
      x39 = x35*x38
      x40 = x2*x39
      x41 = -x33 + x40
      x42 = exp(-3.3019272488946267d0*x10*x3)
      x43 = mu*x42
      x44 = 1d0/mu
      x45 = x28*x44
      x46 = mu**3
      x47 = mu*x32
      x48 = x39*x46 - 0.825481812223657d0*x47 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x30*x45)
      x49 = -x41*x43 + x48
      x50 = 0.275160604074552d0*x47
      x51 = x49*x50 - 0.375d0
      x52 = x27*x51
      x53 = x26*x52
      x54 = rho_a**0.333333333333333d0
      x55 = -x24 + 1.0d0
      x56 = 1d0/x55
      x57 = -0.804d0*x56 + 1.804d0
      x58 = x50*(x43*(x33 - x40) + x48) - 0.375d0
      x59 = x27*x58
      x60 = 1d0/x0
      x61 = x38/x15
      x62 = rho_a**(-0.666666666666667d0)
      x63 = x42*x45
      x64 = 0.183440402716368d0*x60
      x65 = x2*x61 - x29*x64
      x66 = rho_a**(-0.333333333333333d0)*x4
      x67 = x42*x44*x66
      x68 = 0.275160604074552d0*mu*x29*x60 + 2.20128483259642d0*x41*x67 
     &+ x43*x65 - x46*x61 + 1.21141372855476d0*x62*x63
      x69 = 0.275160604074552d0*x31
      x70 = 0.0917202013581841d0*x49*x60 - x68*x69
      x71 = 0.318309886183791d0
      x72 = mu*x1*x71
      x73 = x55**(-2)
      x74 = x11 - x14 - x19
      x75 = 1d0/x12
      x76 = rho_a**(-4.33333333333333d0)*x75
      x77 = 0.229406311391776d0*x76
      x78 = rho_a**(-3.66666666666667d0)*x6
      x79 = 0.637880379497189d0*x78
      x80 = 0.122293601810912d0*x66
      x81 = 0.0134601525394973d0*x13*x54
      x82 = 0.295714285714286d0*x18
      x83 = x34*x82
      x84 = x74*x77 - x74*x79 + x8*(x80 - x81 - x83)
      x85 = 5.84385982654816d0*x0
      x86 = x22*x85
      x87 = x59*x84*x86
      x88 = x25**(-2)
      x89 = 1.39788132145759d0*x23*x60
      x90 = 0.148394183600699d0*x3
      x91 = mu*x70*x71
      x92 = -x80 + x81 + x83
      x93 = -x20*x77 + x20*x79 - x8*x92
      x94 = x52*x93
      x95 = gamma_aa*x88
      x96 = x22*x95
      x97 = 0.458812622783553d0*x76
      x98 = rho_a**(-1.33333333333333d0)*x4
      x99 = 1.27576075899438d0*x78
      x100 = exp(-2.8771076053175544d0*x21)
      x101 = x100/x25**3
      x102 = rho_a**(-2.33333333333333d0)
      x103 = x102*x29
      x104 = 0.166666666666667d0*rho_a**(-3.0d0)*x37
      x105 = x42/x46
      x106 = x58*x90
      Ea = x1*x53
      d1Ea(1) = gamma_aa*x73*x87 + 9.69130982843808d0*x54*x57*x59 - x57*
     &x70*x72
      d1Ea(2) = -x51*x88*x89*x90
      d2Ea(1) = gamma_aa**2*x101*x85*x94*(-x20*x97 + x20*x99 - 0.4784102
     &84622892d0*x7* x92) - 11.6877196530963d0*x0*x91*x93*x96 - 19.38261
     &96568762d0*x26 *x54*x91 + x26*x72*(0.122293601810912d0*x102*x49 - 
     &x64*x68 + x69* (-0.366880805432736d0*mu*x103 - 0.80760915236984d0*
     &rho_a**( -1.66666666666667d0)*x63 + x104*x46 - 4.84565491421904d0*
     &x105*x12 *x41*x62 - 2.66666666666667d0*x105*x35*x36 - 0.7337616108
     &65473d0* x41*x42*x44*x98 - x43*(-0.244587203621824d0*x103 + x104*x
     &2) - 4.40256966519284d0*x65*x67)) + x52*x86*x95*(-0.22000888107104
     &2d0* rho_a**(-6.0d0)*x2*x20/x16 + 1.60584417974243d0*rho_a**( -5.3
     &3333333333333d0)*x20*x75 - 2.33889472482303d0*rho_a**( -4.66666666
     &666667d0)*x20*x6 - x8*(0.00448671751316577d0*x13*x62 + x82 + 0.040
     &7645339369707d0*x98) - x92*x97 + x92*x99) + 3.23043660947936d0*x53
     &*x62 + 15.5836262041284d0*x54*x94*x96
      d2Ea(2) = x73*(2.79576264291519d0*gamma_aa*x100*x106*x20*x56*x60*x
     &84 - 1.86384176194346d0*x102*x106*x23 + 0.0691804894611514d0*x44*x
     &70* x89 + x87)
      d2Ea(3) = 0.0215685619967539d0*rho_a**(-4.0d0)*x101*x17*x20**2*x51
     &
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_2_3(rho_a, gamma_aa, mu, Ea, d
     &1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x127, x128, x129, x130, x131, x132, x133, x134, x135, x136, x137
     &, x138, x139, x140, x141, x142
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 7.26848237132856d0*x0
      x2 = rho_a**0.666666666666667d0
      x3 = 4.60115111447049d0
      x4 = mu**2
      x5 = 1d0/x4
      x6 = x3*x5
      x7 = 3.30192724889463d0*x2*x6
      x8 = exp(x7)
      x9 = 3.0d0*x8
      x10 = 1d0/x0
      x11 = 21.1705915781930d0
      x12 = 1d0/x11
      x13 = mu**4
      x14 = x12*x13*(x8 - 1.0d0)
      x15 = x10*x14
      x16 = x9 - 2.0d0
      x17 = 1d0/x3
      x18 = x17*x4
      x19 = x18/x2
      x20 = 2.14502939711103d0
      x21 = 1d0/x20
      x22 = rho_a**0.333333333333333d0
      x23 = 1d0/x22
      x24 = x21*x23
      x25 = mu*x24
      x26 = 1d0/mu
      x27 = x20*x26
      x28 = erf(1.8171205928321397d0*x22*x27)
      x29 = x28*x8
      x30 = 0.183440402716368d0*x15 - 0.60570686427738d0*x16*x19 + 3.901
     &67577847542d0*x25*x29 - x9
      x31 = 1d0/x30
      x32 = 0.015286700226364d0*x15
      x33 = 7.0d0*x8 + 11.0d0
      x34 = 0.00280419844572861d0*x19*x33
      x35 = -x32 + x34 + 0.0185185185185185d0
      x36 = rho_a**(-2.0d0)
      x37 = x36*x5
      x38 = 3.15935190989552d0*x37
      x39 = x35*x38
      x40 = x31*x39
      x41 = exp(-1.4385538026587772d0*x19)
      x42 = gamma_aa*x41
      x43 = x40*x42
      x44 = x43 - 1.0d0
      x45 = 1.804d0 + 0.804d0/x44
      x46 = 0.682784063255296d0
      x47 = 0.550321208149104d0*x24
      x48 = 1d0/rho_a
      x49 = 9.86960440108936d0
      x50 = 1d0/x49
      x51 = 0.0833333333333333d0*x50
      x52 = x48*x51
      x53 = x4*x52
      x54 = -x47 + x53
      x55 = exp(-x7)
      x56 = mu*x55
      x57 = mu**3
      x58 = -0.825481812223657d0*x25 + 1.77245385090552d0*x28 + x52*x57
      x59 = -x54*x56 + x58
      x60 = 0.275160604074552d0*x25
      x61 = x46*(x59*x60 - 0.375d0)
      x62 = x45*x61
      x63 = rho_a**0.333333333333333d0
      x64 = -x43 + 1.0d0
      x65 = 1d0/x64
      x66 = -0.804d0*x65 + 1.804d0
      x67 = x46*(x60*(x56*(x47 - x53) + x58) - 0.375d0)
      x68 = 0.318309886183791d0
      x69 = mu*x10*x21
      x70 = x36*x51
      x71 = rho_a**(-0.666666666666667d0)
      x72 = x27*x55
      x73 = 0.183440402716368d0*x10
      x74 = -x21*x73 + x4*x70
      x75 = rho_a**(-0.333333333333333d0)
      x76 = x26*x3*x55*x75
      x77 = 2.20128483259642d0*x54*x76 + x56*x74 - x57*x70 + 0.275160604
     &074552d0*x69 + 1.21141372855476d0*x71*x72
      x78 = 0.275160604074552d0*x23
      x79 = x68*(0.0917202013581841d0*x10*x59 - x77*x78)
      x80 = mu*x1
      x81 = x31/x64**2
      x82 = x32 - x34 - 0.0185185185185185d0
      x83 = rho_a**(-3.66666666666667d0)*x17
      x84 = 3.02993180261165d0*x83
      x85 = rho_a**(-3.0d0)
      x86 = x5*x85
      x87 = 6.31870381979104d0*x86
      x88 = x48*x8
      x89 = 0.0432098765432099d0*x88
      x90 = rho_a**(-1.66666666666667d0)
      x91 = x18*x90
      x92 = x8*x91
      x93 = 0.0336503813487433d0*x92
      x94 = rho_a**(-2.33333333333333d0)
      x95 = x14*x94
      x96 = 0.0203822669684853d0*x95
      x97 = 0.00186946563048574d0*x33*x91
      x98 = 0.564189583547756d0
      x99 = x48*x98
      x100 = x6*x8
      x101 = 6.60385449778925d0*x100
      x102 = x16*x91
      x103 = x29*x69
      x104 = x71*x8
      x105 = x104*x27*x28
      x106 = x101*x75 - 0.40380457618492d0*x102 + 1.30055859282514d0*x10
     &3 - 8.58869971286677d0*x105 + 4.0d0*x88 - 0.40380457618492d0*x92 +
     & 0.244587203621824d0*x95 - 4.72654360241471d0*x99
      x107 = x106*x31
      x108 = x107*x38
      x109 = x108*x82 + x38*(-x89 + x93 - x96 + x97) + x82*x84 - x82*x87
     &
      x110 = x41*x67
      x111 = 5.84385982654816d0*x0
      x112 = x109*x110*x111
      x113 = x35*x5
      x114 = 18.4628097041666d0*x41*x71
      x115 = x31/x44**2
      x116 = x115*x61
      x117 = mu*x79
      x118 = x21*x94
      x119 = mu*x118
      x120 = 0.166666666666667d0*x50*x85
      x121 = rho_a**(-1.66666666666667d0)
      x122 = 1d0/x57
      x123 = x122*x55
      x124 = rho_a**(-1.33333333333333d0)
      x125 = x89 - x93 + x96 - x97
      x126 = -x108*x35 - x125*x38 - x35*x84 + x35*x87
      x127 = x116*x42
      x128 = 6.0598636052233d0*x83
      x129 = x128*x35
      x130 = 12.6374076395821d0*x86
      x131 = x130*x35
      x132 = 6.31870381979104d0*x37
      x133 = x125*x132
      x134 = x30**(-2)
      x135 = exp(-2.8771076053175544d0*x19)
      x136 = x134*x135*x61/x44**3
      x137 = x36*x8
      x138 = rho_a**(-2.66666666666667d0)*x18
      x139 = x138*x8
      x140 = rho_a**(-3.33333333333333d0)*x14
      x141 = 1d0/x13
      x142 = x27*x29
      Ea = x1*x62
      d1Ea(1) = gamma_aa*x112*x81 + 9.69130982843808d0*x63*x66*x67 - x66
     &*x79*x80
      d1Ea(2) = -x113*x114*x116
      d2Ea(1) = gamma_aa**2*x111*x126*x136*(-x107*x132*x35 - x129 + x131
     & - x133) - 11.6877196530963d0*x0*x115*x117*x126*x42 + x111*x127*( 
     &-2.90581327762917d0*rho_a**(-5.33333333333333d0)*x12*x35*x4 + 17.1
     &696135481327d0*rho_a**(-4.66666666666667d0)*x17*x35 - 18.956111459
     &3731d0*rho_a**(-4.0d0)*x113 - x106*x134*x39*( 13.2077089955785d0*x
     &100*x75 - 0.80760915236984d0*x102 + 2.60111718565028d0*x103 - 17.1
     &773994257335d0*x105 + 8.0d0*x88 - 0.80760915236984d0*x92 + 0.48917
     &4407243648d0*x95 - 9.45308720482942d0*x99) - x107*x129 + x107*x131
     & - x107*x133 - x125*x128 + x125*x130 - x38*(0.0951172458529316d0*x
     &100*x124 - 0.146090534979424d0*x137 + 0.00311577605080957d0*x138*x
     &33 + 0.10095114404623d0*x139 - 0.0475586229264658d0*x140) - x40*(x
     &101* x124 + 14.5369647426571d0*x104*x11*x141 - 1.73407812376686d0*
     &x119 *x29 + 5.72579980857785d0*x121*x142 - 18.9061744096588d0*x122
     &*x28 *x49*x88 - 27.0091844287146d0*x124*x5 - 7.55555555555556d0*x1
     &37 + 0.673007626974866d0*x138*x16 + 1.21141372855476d0*x139 - 0.57
     &070347511759d0*x140 + 2.86289990428892d0*x142*x90 + 6.302058136552
     &94d0*x36*x98)) - 19.3826196568762d0*x117*x45*x63 + 15.583626204128
     &4d0*x126*x127*x63 + x45*x68*x80*( 0.122293601810912d0*x59*x94 - x7
     &3*x77 + x78*(-4.84565491421904d0* x11*x123*x54*x71 - 0.36688080543
     &2736d0*x119 + x120*x57 - 0.80760915236984d0*x121*x72 - 2.666666666
     &66667d0*x123*x48*x49 - 0.733761610865473d0*x124*x26*x3*x54*x55 - x
     &56*( -0.244587203621824d0*x118 + x120*x4) - 4.40256966519284d0*x74
     &*x76 )) + 3.23043660947936d0*x62*x71
      d2Ea(2) = x81*(36.9256194083333d0*gamma_aa*x109*x113*x135*x31*x65*
     &x67*x71 - 24.6170796055555d0*x110*x113*x121 + x112 + x114*x26*x35*
     &x79)
      d2Ea(3) = 116.661026201793d0*rho_a**(-2.66666666666667d0)*x136*x14
     &1*x35**2
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_PBE_ERF_case_3(rho_a, gamma_aa, mu, Ea, d1E
     &a, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, mu
      real*8, intent(out) :: Ea, d1Ea(4), d2Ea(10)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**2.0d0
      x1 = mu**2
      x2 = 1d0/x1
      x3 = 4.60115111447049d0
      x4 = 1d0/x3
      x5 = x2*x4
      x6 = rho_a**(-2.66666666666667d0)
      x7 = 0.239205142311446d0*x6
      x8 = x5*x7
      x9 = rho_a**0.666666666666667d0
      x10 = x1*x4/x9
      x11 = exp(-1.4385538026587772d0*x10)
      x12 = 0.183440402716368d0*x3*x9
      x13 = 21.1705915781930d0
      x14 = x13*x2
      x15 = 0.010095114404623d0*rho_a**1.33333333333333d0*x14
      x16 = 97.4090910340024d0
      x17 = mu**(-4)
      x18 = x16*x17
      x19 = 0.147857142857143d0*x0*x18
      x20 = -x12 + x15 + x19
      x21 = x11*x20
      x22 = gamma_aa*x21*x8
      x23 = x22 - 1.0d0
      x24 = 3.14159265358979d0
      x25 = x2*x24
      x26 = rho_a**1.0d0
      x27 = -x22 + 1.0d0
      x28 = 1d0/x27
      x29 = x12 - x15 - x19
      x30 = 1d0/x13
      x31 = rho_a**(-4.33333333333333d0)*x30
      x32 = x29*x31
      x33 = rho_a**(-3.66666666666667d0)
      x34 = x29*x5
      x35 = x33*x34
      x36 = 0.295714285714286d0*x18
      x37 = 0.122293601810912d0*rho_a**(-0.333333333333333d0)*x3 - 0.013
     &4601525394973d0*rho_a**0.333333333333333d0*x14 - x26*x36
      x38 = x37*x5
      x39 = 0.229406311391776d0*x32 - 0.637880379497189d0*x35 + x38*x7
      x40 = x11*x39
      x41 = 0.804d0*x0
      x42 = x40*x41
      x43 = x27**(-2)
      x44 = gamma_aa*x43
      x45 = rho_a**(-0.666666666666667d0)
      x46 = 0.682784063255296d0
      x47 = x21*x46
      x48 = exp(-2.8771076053175544d0*x10)
      x49 = x39*x48
      x50 = 0.384641868836805d0*x2
      Ea = -1.0d0*x0*x25*(1.804d0 + 0.804d0/x23)
      d1Ea(1) = -x25*(2.0d0*x26*(-0.804d0*x28 + 1.804d0) + x42*x44)
      d1Ea(2) = 0.192320934418403d0*x17*x45*x47/x23**2
      d2Ea(1) = x25*(gamma_aa**2*x41*x49*(0.458812622783553d0*x32 - 1.27
     &576075899438d0* x35 + 0.478410284622892d0*x38*x6)/x27**3 - x11*x41
     &*x44*( 0.220008881071042d0*rho_a**(-6.0d0)*x1*x29/x16 - 1.60584417
     &974243d0*rho_a**(-5.33333333333333d0)*x29*x30 + 2.33889472482303d0
     &*rho_a**(-4.66666666666667d0)*x34 + 0.458812622783553d0*x31*x37 - 
     &1.27576075899438d0*x33*x38 - x8*( 0.0407645339369707d0*rho_a**(-1.
     &33333333333333d0)*x3 + 0.00448671751316577d0*rho_a**(-0.6666666666
     &66667d0)*x14 + x36)) - 3.216d0*x26*x40*x44 + 1.608d0*x28 - 3.608d0
     &)
      d2Ea(2) = x2*x43*(-gamma_aa*x20*x28*x45*x46*x49*x50 + rho_a**(-1.6
     &6666666666667d0) *x47*x50 - x24*x42)
      d2Ea(3) = -0.0136534984882616d0*rho_a**(-3.33333333333333d0)*x20**
     &2*x48/(mu**6*x23 **3)
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu, E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu*rho_a**(-0.333333333333333d0)
      x2 = 0.0472353356922751d0
      x3 = 0.217336917462899d0
      x4 = mu**2*rho_a**(-0.666666666666667d0)*x3
      x5 = gamma_aa*rho_a**(-2.66666666666667d0)*x3
      x6 = gamma_aa**2
      x7 = rho_a**(-5.33333333333333d0)*x2
      x8 = 1/(rho_a**2*tau_a**2)
      x9 = x6*x8
      x10 = 0.050475572023115d0*x5
      x11 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x12 = x11*x5
      x13 = (x10*x11*(0.126188930057787d0*x12 - 1.0d0) + 1.0d0)**(-0.5d0
     &)
      Ea = -x0*(1.8610514726982d0*(0.804d0/(0.0206714418186455d0*x5 + 1.
     &0d0) - 0.804d0/(1.0d0 + (1.86977850568311d-6*gamma_aa**3*rho_a**(-
     &8.0d0 ) + 0.0941708433267071d0*x5*(0.02485875d0*x9/(0.015625d0*x9 
     &+ 1.0d0)**2 + 0.123456790123457d0) + 0.000135164574498753d0*x6*x7 
     &+ 0.0021416494625995d0*x9 + 0.0896750813832074d0*(x10 + x13*( 0.05
     &67850185260044d0*x12 - 0.45d0))**2 - 1.24378109452736d0*sqrt( 0.00
     &286625629244325d0*x7 + 0.0028125d0*x8)*(x13*( 0.0102353243269094d0
     &*x12 - 0.0811111111111111d0) + 0.00909806606836394d0*x5)*abs(gamma
     &_aa))/(0.0938662444277524d0*x5 + 1.0d0)**2))*exp(-1.92417365779563
     &32d0*x1) + 3.35733685674755d0 - 1.49628538404935d0/(0.239205142311
     &446d0*x5*(5.45136177849642d0* x4 - 7.0d0)*exp(-1.4385538026587772d
     &0*x4)/(4.95289087334194d0*mu **4*x2/x0 + 49.1113297378199d0*x1 - 4
     &9.0622560064678d0*x4 - 81.0d0) + 1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, mu,
     & Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**0.333333333333333d0
      x2 = 1d0/x1
      x3 = 0.0472353356922751d0
      x4 = 4.60115111447049d0
      x5 = 1d0/x4
      x6 = mu**2
      x7 = rho_a**0.666666666666667d0
      x8 = x5*x6/x7
      x9 = gamma_aa*rho_a**(-2.66666666666667d0)*x5
      x10 = 2.14502939711103d0
      x11 = x2/x10
      x12 = mu*x11
      x13 = 0.00844343197019482d0*1d0/rho_a
      x14 = gamma_aa**2
      x15 = rho_a**(-5.33333333333333d0)*x3
      x16 = 1/(rho_a**2*tau_a**2)
      x17 = x14*x16
      x18 = 0.050475572023115d0*x9
      x19 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x20 = x19*x9
      x21 = (x18*x19*(0.126188930057787d0*x20 - 1.0d0) + 1.0d0)**(-0.5d0
     &)
      Ea = -x0*(-4.9628039271952d0*(1.804d0 - 0.804d0/(0.239205142311446
     &d0*x9*( 5.45136177849642d0*x8 - 7.0d0)*exp(-1.4385538026587772d0*x
     &8)/( 4.95289087334194d0*mu**4*x3/x0 + 49.1113297378199d0*mu*x2 - 4
     &9.0622560064678d0*x8 - 81.0d0) + 1.0d0))*(0.275160604074552d0* x12
     &*(mu**3*x13 + mu*(0.550321208149104d0*x11 - x13*x6)*exp( -3.301927
     &2488946267d0*x4*x7/x6) - 0.825481812223657d0*x12 + 1.7724538509055
     &2d0*erf(1.8171205928321397d0*x1*x10/mu)) - 0.375d0 ) + 1.861051472
     &6982d0*(0.804d0/(0.0206714418186455d0*x9 + 1.0d0) - 0.804d0/(1.0d0
     & + (1.86977850568311d-6*gamma_aa**3*rho_a**( -8.0d0) + 0.000135164
     &574498753d0*x14*x15 + 0.0021416494625995d0* x17 + 0.09417084332670
     &71d0*x9*(0.02485875d0*x17/(0.015625d0*x17 + 1.0d0)**2 + 0.12345679
     &0123457d0) - 1.24378109452736d0*sqrt( 0.00286625629244325d0*x15 + 
     &0.0028125d0*x16)*(x21*( 0.0102353243269094d0*x20 - 0.0811111111111
     &111d0) + 0.00909806606836394d0*x9)*abs(gamma_aa) + 0.0896750813832
     &074d0*( x18 + x21*(0.0567850185260044d0*x20 - 0.45d0))**2)/( 0.093
     &8662444277524d0*x9 + 1.0d0)**2))*exp(-4.1274090611182839d0* x12))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, mu,
     & Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = mu**2
      x2 = 1d0/x1
      x3 = 4.60115111447049d0
      x4 = rho_a**0.666666666666667d0
      x5 = x3*x4
      x6 = 21.1705915781930d0
      x7 = 97.4090910340024d0
      x8 = 1d0/x3
      x9 = gamma_aa*rho_a**(-2.66666666666667d0)*x8
      x10 = 2.14502939711103d0
      x11 = rho_a**0.333333333333333d0
      x12 = 1/(x10*x11)
      x13 = mu*x12
      x14 = 0.00844343197019482d0*1d0/rho_a
      x15 = gamma_aa**2
      x16 = rho_a**(-5.33333333333333d0)/x6
      x17 = 1/(rho_a**2*tau_a**2)
      x18 = x15*x17
      x19 = 0.050475572023115d0*x9
      x20 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x21 = x20*x9
      x22 = (x19*x20*(0.126188930057787d0*x21 - 1.0d0) + 1.0d0)**(-0.5d0
     &)
      Ea = -x0*(-4.9628039271952d0*(1.804d0 - 0.804d0/(-0.23920514231144
     &6d0*x2*x9*( 0.010095114404623d0*x0*x2*x6 - 0.183440402716368d0*x5 
     &+ 0.147857142857143d0*rho_a**2.0d0*x7/mu**4)*exp( -1.4385538026587
     &772d0*x1*x8/x4) + 1.0d0))*(0.275160604074552d0* x13*(mu**3*x14 + m
     &u*(-x1*x14 + 0.550321208149104d0*x12)*exp( -3.3019272488946267d0*x
     &2*x5) - 0.825481812223657d0*x13 + 1.77245385090552d0*erf(1.8171205
     &928321397d0*x10*x11/mu)) - 0.375d0) + 1.8610514726982d0*(0.804d0/(
     &0.0206714418186455d0*x9 + 1.0d0) - 0.804d0/(1.0d0 + (0.00018213342
     &4673508d0*gamma_aa**3* rho_a**(-8.0d0)/x7 + 0.000135164574498753d0
     &*x15*x16 + 0.0021416494625995d0*x18 + 0.0941708433267071d0*x9*(0.0
     &2485875d0* x18/(0.015625d0*x18 + 1.0d0)**2 + 0.123456790123457d0) 
     &- 1.24378109452736d0*sqrt(0.00286625629244325d0*x16 + 0.0028125d0*
     & x17)*(x22*(0.0102353243269094d0*x21 - 0.0811111111111111d0) + 0.0
     &0909806606836394d0*x9)*abs(gamma_aa) + 0.0896750813832074d0*( x19 
     &+ x22*(0.0567850185260044d0*x21 - 0.45d0))**2)/( 0.093866244427752
     &4d0*x9 + 1.0d0)**2))*exp(-4.1274090611182839d0* x13))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, mu,
     & Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6
      Ea = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 2.14502939711103d0
      x2 = rho_a**0.333333333333333d0
      x3 = 1/(x1*x2)
      x4 = mu*x3
      x5 = 0.00844343197019482d0*1d0/rho_a
      x6 = erf(1.8171205928321397d0*x1*x2/mu)
      x7 = 4.60115111447049d0
      x8 = mu**2
      x9 = 1d0/x8
      x10 = rho_a**0.666666666666667d0
      x11 = 3.30192724889463d0*x10*x7*x9
      x12 = 1d0/x7
      x13 = x12*x8/x10
      x14 = 0.0472353356922751d0
      x15 = exp(x11)
      x16 = mu**4*x14*(x15 - 1.0d0)/x0
      x17 = 3.0d0*x15
      x18 = gamma_aa*rho_a**(-2.66666666666667d0)*x12
      x19 = gamma_aa**2
      x20 = rho_a**(-5.33333333333333d0)*x14
      x21 = 1/(rho_a**2*tau_a**2)
      x22 = x19*x21
      x23 = 0.050475572023115d0*x18
      x24 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x25 = x18*x24
      x26 = (x23*x24*(0.126188930057787d0*x25 - 1.0d0) + 1.0d0)**(-0.5d0
     &)
      Ea = -x0*(-4.9628039271952d0*(1.804d0 - 0.804d0/(-3.15935190989552
     &d0*gamma_aa *rho_a**(-2.0d0)*x9*(0.00280419844572861d0*x13*(7.0d0*
     &x15 + 11.0d0) - 0.015286700226364d0*x16 + 0.0185185185185185d0)*ex
     &p( -1.4385538026587772d0*x13)/(-0.60570686427738d0*x13*(x17 - 2.0d
     &0 ) + 3.90167577847542d0*x15*x4*x6 + 0.183440402716368d0*x16 - x17
     & ) + 1.0d0))*(0.275160604074552d0*x4*(mu**3*x5 + mu*( 0.5503212081
     &49104d0*x3 - x5*x8)*exp(-x11) - 0.825481812223657d0* x4 + 1.772453
     &85090552d0*x6) - 0.375d0) + 1.8610514726982d0*( 0.804d0/(0.0206714
     &418186455d0*x18 + 1.0d0) - 0.804d0/(1.0d0 + ( 1.86977850568311d-6*
     &gamma_aa**3*rho_a**(-8.0d0) + 0.0941708433267071d0*x18*(0.02485875
     &d0*x22/(0.015625d0*x22 + 1.0d0)**2 + 0.123456790123457d0) + 0.0001
     &35164574498753d0*x19*x20 + 0.0021416494625995d0*x22 - 1.2437810945
     &2736d0*( 0.00909806606836394d0*x18 + x26*(0.0102353243269094d0*x25
     & - 0.0811111111111111d0))*sqrt(0.00286625629244325d0*x20 + 0.00281
     &25d0*x21)*abs(gamma_aa) + 0.0896750813832074d0*(x23 + x26* (0.0567
     &850185260044d0*x25 - 0.45d0))**2)/(0.0938662444277524d0* x18 + 1.0
     &d0)**2))*exp(-4.1274090611182839d0*x4))
      end subroutine


C*****************************************************************************
      pure subroutine ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu, E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_a, gamma_aa, tau_a, mu
      real*8, intent(out) :: Ea
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16
      Ea = 0.0d0
      x0 = mu**2
      x1 = 1d0/x0
      x2 = 4.60115111447049d0
      x3 = rho_a**0.666666666666667d0
      x4 = 21.1705915781930d0
      x5 = 97.4090910340024d0
      x6 = 1d0/x2
      x7 = gamma_aa*rho_a**(-2.66666666666667d0)*x6
      x8 = rho_a**0.333333333333333d0
      x9 = gamma_aa**2
      x10 = rho_a**(-5.33333333333333d0)/x4
      x11 = 1/(rho_a**2*tau_a**2)
      x12 = x11*x9
      x13 = 0.050475572023115d0*x7
      x14 = -1.0d0 + 8.0d0*rho_a*tau_a/gamma_aa
      x15 = x14*x7
      x16 = (x13*x14*(0.126188930057787d0*x15 - 1.0d0) + 1.0d0)**(-0.5d0
     &)
      Ea = -rho_a*(3.14159265358979d0*rho_a**1.0d0*x1*(1.804d0 - 0.804d0
     &/( -0.239205142311446d0*x1*x7*(0.010095114404623d0*rho_a** 1.33333
     &333333333d0*x1*x4 - 0.183440402716368d0*x2*x3 + 0.147857142857143d
     &0*rho_a**2.0d0*x5/mu**4)*exp( -1.4385538026587772d0*x0*x6/x3) + 1.
     &0d0)) + 1.8610514726982d0*x8* (0.804d0/(0.0206714418186455d0*x7 + 
     &1.0d0) - 0.804d0/(1.0d0 + ( 0.000182133424673508d0*gamma_aa**3*rho
     &_a**(-8.0d0)/x5 + 0.000135164574498753d0*x10*x9 + 0.00214164946259
     &95d0*x12 + 0.0941708433267071d0*x7*(0.02485875d0*x12/(0.015625d0*x
     &12 + 1.0d0 )**2 + 0.123456790123457d0) - 1.24378109452736d0*sqrt( 
     &0.00286625629244325d0*x10 + 0.0028125d0*x11)*(x16*( 0.010235324326
     &9094d0*x15 - 0.0811111111111111d0) + 0.00909806606836394d0*x7)*abs
     &(gamma_aa) + 0.0896750813832074d0*( x13 + x16*(0.0567850185260044d
     &0*x15 - 0.45d0))**2)/( 0.0938662444277524d0*x7 + 1.0d0)**2))*exp(-
     &1.9241736577956332d0* mu/x8))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x116, x117, x118, x119, x120, x121, x122
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.682784063255296d0
      x2 = rho_a**(-2.66666666666667d0)
      x3 = 0.217336917462899d0
      x4 = x2*x3
      x5 = gamma_aa*x4
      x6 = 0.826307487110758d0
      x7 = rho_a**0.333333333333333d0
      x8 = mu/x7
      x9 = 0.0472353356922751d0
      x10 = 1d0/x0
      x11 = mu**2
      x12 = rho_a**(-0.666666666666667d0)*x11*x3
      x13 = 1d0/(4.95289087334194d0*mu**4*x10*x9 - 49.0622560064678d0*x1
     &2 + 59.4346904801033d0*x6*x8 - 81.0d0)
      x14 = x13*exp(-1.4385538026587772d0*x12)
      x15 = 5.45136177849642d0*x12 - 7.0d0
      x16 = 0.239205142311446d0*x15
      x17 = x14*x16*x5 + 1.0d0
      x18 = x1*(1.804d0 - 0.804d0/x17)
      x19 = 2.72568088924821d0*x18
      x20 = exp(-1.9241736577956332d0*x8)
      x21 = 0.0206714418186455d0*x5 + 1.0d0
      x22 = 0.804d0/x21
      x23 = 0.0938662444277524d0*x5 + 1.0d0
      x24 = x23**(-2)
      x25 = 0.050475572023115d0*x5
      x26 = tau_a/gamma_aa
      x27 = 8.0d0*rho_a*x26
      x28 = x27 - 1.0d0
      x29 = x28*x5
      x30 = 0.0567850185260044d0*x29 - 0.45d0
      x31 = 0.126188930057787d0*x28
      x32 = x31*x5 - 1.0d0
      x33 = x28*x32
      x34 = x25*x33 + 1.0d0
      x35 = x34**(-0.5d0)
      x36 = x25 + x30*x35
      x37 = x36**2
      x38 = 0.00909806606836394d0*x4
      x39 = gamma_aa*x38
      x40 = x35*(0.0102353243269094d0*x29 - 0.0811111111111111d0) + x39
      x41 = rho_a**(-5.33333333333333d0)*x9
      x42 = rho_a**(-2)
      x43 = tau_a**(-2)
      x44 = 0.0028125d0*x43
      x45 = sqrt(0.00286625629244325d0*x41 + x42*x44)
      x46 = abs(gamma_aa)
      x47 = x45*x46
      x48 = x40*x47
      x49 = gamma_aa**2
      x50 = x41*x49
      x51 = rho_a**(-8.0d0)
      x52 = 0.0102659822546843d0
      x53 = gamma_aa**3
      x54 = x52*x53
      x55 = x51*x54
      x56 = x43*x49
      x57 = x42*x56
      x58 = 0.015625d0*x57 + 1.0d0
      x59 = x58**(-2)
      x60 = 0.02485875d0*x57*x59 + 0.123456790123457d0
      x61 = 0.0941708433267071d0*x60
      x62 = x5*x61 + 0.000135164574498753d0*x50 + 0.000182133424673508d0
     &*x55 + 0.0021416494625995d0*x57
      x63 = x24*(0.0896750813832074d0*x37 - 1.24378109452736d0*x48 + x62
     &) + 1.0d0
      x64 = -0.804d0/x63
      x65 = x20*(x22 + x64)
      x66 = 0.682784063255296d0
      x67 = 2.72568088924821d0*x66
      x68 = x65*x67
      x69 = rho_a**(-0.666666666666667d0)
      x70 = rho_a**(-4.33333333333333d0)*x9
      x71 = x11*x70
      x72 = rho_a**(-3.66666666666667d0)
      x73 = rho_a**(-1.66666666666667d0)
      x74 = x3*x73
      x75 = rho_a**(-2.33333333333333d0)
      x76 = x14/x17**2
      x77 = 0.148394183600699d0
      x78 = x77/x21**2
      x79 = 0.90856029641607d0*x66*x69
      x80 = rho_a**(-6.33333333333333d0)*x9
      x81 = rho_a**(-3)
      x82 = gamma_aa*x72
      x83 = x3*x82
      x84 = 0.001553671875d0*x57/x58 - 0.0497175d0
      x85 = 0.0941708433267071d0*x3*x59*x84
      x86 = 1d0/x45
      x87 = 1d0/x46
      x88 = 1.24378109452736d0*x87
      x89 = 1d0/x23
      x90 = 0.0757133580346725d0*x5*x60 + 0.000108672317896998d0*x50 + 0
     &.0001464352734375d0*x55 + 0.00172188616792999d0*x57
      x91 = tau_a*x2
      x92 = x28*x82
      x93 = x35*(0.454280148208035d0*x91 - 0.151426716069345d0*x92)
      x94 = 0.180246913580247d0*x3
      x95 = x3*x32
      x96 = 0.0252377860115575d0*gamma_aa*x2*x28*x9
      x97 = x30*x34**(-1.5d0)*(-0.20190228809246d0*x91*x95 + 0.067300762
     &6974867d0* x92*x95 - x96*(1.0095114404623d0*x91 - 0.33650381348743
     &3d0*x92))
      x98 = 1.24378109452736d0*x47
      x99 = 2.0d0*x3
      x100 = x24*x66*x7
      x101 = -x27 + 1.0d0
      x102 = 0.0567850185260044d0*x101
      x103 = x102*x5 + 0.45d0
      x104 = x101*x5
      x105 = 0.126188930057787d0*x104 + 1.0d0
      x106 = x101*x105*x25 + 1.0d0
      x107 = x106**(-0.5d0)
      x108 = -x103*x107 + x25
      x109 = x108**2
      x110 = -x107*(0.0102353243269094d0*x104 + 0.0811111111111111d0) + 
     &x39
      x111 = x110*x47
      x112 = (x24*(0.0896750813832074d0*x109 - 1.24378109452736d0*x111 +
     & x62) + 1.0d0 )**(-2)
      x113 = rho_a**(-4.66666666666667d0)
      x114 = x26*x73
      x115 = x107*(x102*x2 + 0.454280148208035d0*x114)
      x116 = 0.20190228809246d0*x74
      x117 = x106**(-1.5d0)*(-x116*x26*x32 + 0.0252377860115575d0*x33*x4
     & + x96*( -1.0095114404623d0*x114 + x2*x31))
      x118 = 0.180246913580247d0*x103
      x119 = 2.0d0*x103
      x120 = tau_a**(-3)
      x121 = x120*x42*x49
      x122 = (0.0254778337106067d0*gamma_aa*x101*x70 + x105*x116)/x106
      Ea = -x0*(x19 + x68)
      d1Ea(1) = -rho_a*(2.19144743495556d0*gamma_aa*x1*x7*x76*(mu*x13*x1
     &6*x4*( 6.60385449778925d0*mu**3*x75*x9 - 32.7081706709785d0*mu*x74
     & + 19.8115634933678d0*x10*x6) - 0.637880379497189d0*x15*x3*x72 + 0
     &.229406311391776d0*x15*x71 - 0.869329180010942d0*x71) + 1.19366207
     &318922d0*mu*1d0/rho_a*x65 + 0.90856029641607d0*x18*x69 + x20*(0.12
     &080100840081d0*gamma_aa*rho_a**(-3.33333333333333d0)* x78 + 2.1914
     &4743495556d0*x100*(-0.00145706739738806d0*rho_a**( -9.0d0)*x54 + r
     &ho_a**(-5.66666666666667d0)*x43*x53*x85 + 0.0896750813832074d0*x36
     &*(-0.269203050789947d0*x83 + x93*x99 + 2.0d0*x97) + x40*x49*x86*x8
     &8*(x44*x81 + 0.00764335011318201d0*x80 ) - 0.000720877730660018d0*
     &x49*x80 - 0.00428329892519899d0*x56* x81 - 0.251122248871219d0*x60
     &*x83 + 0.622661654578789d0*x83*x89*( 0.0720987654320988d0*x37 - x4
     &8 + x90) - x98*( -0.0242615095156372d0*x83 + x93*x94 + 0.180246913
     &580247d0*x97))/ x63**2 - x79*(-x22 + 1.804d0) + x79*(x64 + 1.804d0
     &))) - x19*x7 - x68*x7
      d1Ea(2) = -rho_a*(0.524205495546598d0*x15*x75*x76*x77 - 2.0d0*x20*
     &( 1.09572371747778d0*x100*x112*(gamma_aa*x110*x45*x88 - 0.00027032
     &9148997507d0*gamma_aa*x41 - 0.00428329892519899d0* gamma_aa*x42*x4
     &3 - 0.0896750813832074d0*x108*(-x115*x99 + x117* x119 + 0.10095114
     &404623d0*x4) + x113*x56*x85 - x4*x61 + 0.233498120467046d0*x4*x89*
     &(0.0720987654320988d0*x109 - x111 + x90) - 0.000546400274020523d0*
     &x49*x51*x52 + x98*(-x115*x94 + x117 *x118 + x38)) + 0.022650189075
     &1519d0*x75*x78))
      d1Ea(3) = -x0*x112*x20*x24*x67*(0.0720987654320988d0*x107*x108*(-x
     &119*x122 + 0.90856029641607d0*x74) - x107*x47*(-x118*x122 + 0.0818
     &825946152754d0*x74) + 0.0028125d0*x110*x121*x86*x87 + 0.0757133580
     &346725d0*x113*x120*x3*x53*x59*x84 - 0.00344377233585999d0*x121)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = 4.60115111447049d0
      x3 = 1d0/x2
      x4 = x1*x3
      x5 = gamma_aa*x4
      x6 = 0.826307487110758d0
      x7 = rho_a**0.333333333333333d0
      x8 = 1d0/x7
      x9 = 0.0472353356922751d0
      x10 = 1d0/x0
      x11 = mu**2
      x12 = rho_a**0.666666666666667d0
      x13 = x11*x3/x12
      x14 = 1d0/(4.95289087334194d0*mu**4*x10*x9 + 59.4346904801033d0*mu
     &*x6*x8 - 49.0622560064678d0*x13 - 81.0d0)
      x15 = x14*exp(-1.4385538026587772d0*x13)
      x16 = 5.45136177849642d0*x13 - 7.0d0
      x17 = 0.239205142311446d0*x16
      x18 = x15*x17*x5 + 1.0d0
      x19 = 1.804d0 - 0.804d0/x18
      x20 = 7.26848237132856d0*x19
      x21 = 0.682784063255296d0
      x22 = 2.14502939711103d0
      x23 = 1d0/x22
      x24 = x23*x8
      x25 = 0.550321208149104d0*x24
      x26 = 1d0/rho_a
      x27 = 0.00844343197019481d0
      x28 = x26*x27
      x29 = x11*x28
      x30 = exp(-3.3019272488946267d0*x12*x2/x11)
      x31 = mu*x30
      x32 = 1d0/mu
      x33 = x22*x32
      x34 = mu**3
      x35 = mu*x24
      x36 = x28*x34 - 0.825481812223657d0*x35 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x33*x7)
      x37 = 0.275160604074552d0*x35
      x38 = x37*(x31*(x25 - x29) + x36) - 0.375d0
      x39 = x21*x38
      x40 = exp(-4.1274090611182839d0*x35)
      x41 = 0.0206714418186455d0*x5 + 1.0d0
      x42 = 0.804d0/x41
      x43 = 0.0938662444277524d0*x5 + 1.0d0
      x44 = x43**(-2)
      x45 = 0.050475572023115d0*x5
      x46 = tau_a/gamma_aa
      x47 = 8.0d0*rho_a*x46
      x48 = x47 - 1.0d0
      x49 = x48*x5
      x50 = 0.0567850185260044d0*x49 - 0.45d0
      x51 = 0.126188930057787d0*x48
      x52 = x5*x51 - 1.0d0
      x53 = x48*x52
      x54 = x45*x53 + 1.0d0
      x55 = x54**(-0.5d0)
      x56 = x45 + x50*x55
      x57 = x56**2
      x58 = 0.00909806606836394d0*x4
      x59 = gamma_aa*x58
      x60 = x55*(0.0102353243269094d0*x49 - 0.0811111111111111d0) + x59
      x61 = rho_a**(-5.33333333333333d0)*x9
      x62 = rho_a**(-2)
      x63 = tau_a**(-2)
      x64 = 0.0028125d0*x63
      x65 = sqrt(0.00286625629244325d0*x61 + x62*x64)
      x66 = abs(gamma_aa)
      x67 = x65*x66
      x68 = x60*x67
      x69 = gamma_aa**2
      x70 = x61*x69
      x71 = rho_a**(-8.0d0)
      x72 = 0.0102659822546843d0
      x73 = gamma_aa**3
      x74 = x72*x73
      x75 = x71*x74
      x76 = x63*x69
      x77 = x62*x76
      x78 = 0.015625d0*x77 + 1.0d0
      x79 = x78**(-2)
      x80 = 0.02485875d0*x77*x79 + 0.123456790123457d0
      x81 = 0.0941708433267071d0*x80
      x82 = x5*x81 + 0.000135164574498753d0*x70 + 0.000182133424673508d0
     &*x75 + 0.0021416494625995d0*x77
      x83 = x44*(0.0896750813832074d0*x57 - 1.24378109452736d0*x68 + x82
     &) + 1.0d0
      x84 = -0.804d0/x83
      x85 = x40*(x42 + x84)
      x86 = 0.682784063255296d0
      x87 = 2.72568088924821d0*x86
      x88 = -x25 + x29
      x89 = -x31*x88 + x36
      x90 = x20*x7
      x91 = -x47 + 1.0d0
      x92 = 0.0567850185260044d0*x91
      x93 = x5*x92 + 0.45d0
      x94 = x5*x91
      x95 = 0.126188930057787d0*x94 + 1.0d0
      x96 = x45*x91*x95 + 1.0d0
      x97 = x96**(-0.5d0)
      x98 = x45 - x93*x97
      x99 = x98**2
      x100 = x59 - x97*(0.0102353243269094d0*x94 + 0.0811111111111111d0)
     &
      x101 = x100*x67
      x102 = x44*(-1.24378109452736d0*x101 + x82 + 0.0896750813832074d0*
     &x99) + 1.0d0
      x103 = x40*x87
      x104 = rho_a**(-0.666666666666667d0)
      x105 = x10*x23
      x106 = rho_a**(-2.0d0)*x27
      x107 = 0.318309886183791d0*mu
      x108 = rho_a**(-4.33333333333333d0)*x9
      x109 = x108*x11
      x110 = rho_a**(-3.66666666666667d0)
      x111 = rho_a**(-1.66666666666667d0)
      x112 = x111*x3
      x113 = rho_a**(-2.33333333333333d0)
      x114 = x15/x18**2
      x115 = 0.148394183600699d0
      x116 = x115/x41**2
      x117 = 0.90856029641607d0*x104*x86
      x118 = rho_a**(-6.33333333333333d0)*x9
      x119 = rho_a**(-3)
      x120 = gamma_aa*x110
      x121 = x120*x3
      x122 = 0.001553671875d0*x77/x78 - 0.0497175d0
      x123 = 0.0941708433267071d0*x122*x3*x79
      x124 = 1d0/x65
      x125 = 1d0/x66
      x126 = 1.24378109452736d0*x125
      x127 = 1d0/x43
      x128 = 0.0757133580346725d0*x5*x80 + 0.000108672317896998d0*x70 + 
     &0.0001464352734375d0*x75 + 0.00172188616792999d0*x77
      x129 = tau_a*x1
      x130 = x120*x48
      x131 = x55*(0.454280148208035d0*x129 - 0.151426716069345d0*x130)
      x132 = 0.180246913580247d0*x3
      x133 = x3*x52
      x134 = 0.0252377860115575d0*gamma_aa*x1*x48*x9
      x135 = x50*x54**(-1.5d0)*(-0.20190228809246d0*x129*x133 + 0.067300
     &7626974867d0* x130*x133 - x134*(1.0095114404623d0*x129 - 0.3365038
     &13487433d0* x130))
      x136 = 1.24378109452736d0*x67
      x137 = 2.0d0*x3
      x138 = x44*x7*x86
      x139 = x102**(-2)
      x140 = rho_a**(-4.66666666666667d0)
      x141 = x111*x46
      x142 = x97*(x1*x92 + 0.454280148208035d0*x141)
      x143 = 0.20190228809246d0*x112
      x144 = x96**(-1.5d0)*(x134*(x1*x51 - 1.0095114404623d0*x141) - x14
     &3*x46*x52 + 0.0252377860115575d0*x4*x53)
      x145 = 0.180246913580247d0*x93
      x146 = 2.0d0*x93
      x147 = tau_a**(-3)
      x148 = x147*x62*x69
      x149 = (0.0254778337106067d0*gamma_aa*x108*x91 + x143*x95)/x96
      Ea = -x0*(-x20*x39 + x85*x87)
      d1Ea(1) = -rho_a*(-5.84385982654816d0*gamma_aa*x114*x39*x7*(mu*x14
     &*x17*x4*( -32.7081706709785d0*mu*x112 + 19.8115634933678d0*x10*x6 
     &+ 6.60385449778925d0*x113*x34*x9) + 0.229406311391776d0*x109*x16 -
     & 0.869329180010942d0*x109 - 0.637880379497189d0*x110*x16*x3) - 2.4
     &2282745710952d0*x104*x19*x39 + 3.75d0*x107*x26*x85 + x107*x90* (0.
     &0917202013581841d0*x10*x89 - 0.275160604074552d0*x8*( 0.2751606040
     &74552d0*mu*x105 + 2.20128483259642d0*rho_a**( -0.333333333333333d0
     &)*x2*x30*x32*x88 + 1.21141372855476d0*x104* x30*x33 - x106*x34 + x
     &31*(-0.183440402716368d0*x105 + x106*x11 ))) + x40*(0.120801008400
     &81d0*gamma_aa*rho_a**( -3.33333333333333d0)*x116 - x117*(-x42 + 1.
     &804d0) + x117*(x84 + 1.804d0) + 2.19144743495556d0*x138*(-0.001457
     &06739738806d0*rho_a **(-9.0d0)*x74 + rho_a**(-5.66666666666667d0)*
     &x123*x63*x73 - 0.000720877730660018d0*x118*x69 - 0.004283298925198
     &99d0*x119*x76 + 0.622661654578789d0*x121*x127*(x128 + 0.0720987654
     &320988d0*x57 - x68) - 0.251122248871219d0*x121*x80 + x124*x126*x60
     &*x69*( 0.00764335011318201d0*x118 + x119*x64) - x136*( -0.02426150
     &95156372d0*x121 + x131*x132 + 0.180246913580247d0*x135 ) + 0.08967
     &50813832074d0*x56*(-0.269203050789947d0*x121 + x131* x137 + 2.0d0*
     &x135))/x83**2)) - x103*x7*(x42 - 0.804d0/x102) + x21 *x90*(x37*x89
     & - 0.375d0)
      d1Ea(2) = rho_a*(1.39788132145759d0*x113*x114*x115*x16*x38 + 2.0d0
     &*x40*( 0.0226501890751519d0*x113*x116 + 1.09572371747778d0*x138*x1
     &39*( gamma_aa*x100*x126*x65 - 0.000270329148997507d0*gamma_aa*x61 
     &- 0.00428329892519899d0*gamma_aa*x62*x63 + x123*x140*x76 + 0.23349
     &8120467046d0*x127*x4*(-x101 + x128 + 0.0720987654320988d0* x99) + 
     &x136*(-x132*x142 + x144*x145 + x58) - x4*x81 - 0.00054640027402052
     &3d0*x69*x71*x72 - 0.0896750813832074d0*x98*( -x137*x142 + x144*x14
     &6 + 0.10095114404623d0*x4))))
      d1Ea(3) = -x0*x103*x139*x44*(0.0028125d0*x100*x124*x125*x148 + 0.0
     &757133580346725d0*x122*x140*x147*x3*x73*x79 - 0.00344377233585999d
     &0*x148 - x67*x97*(0.0818825946152754d0*x112 - x145*x149) + 0.07209
     &87654320988d0*x97*x98*(0.90856029641607d0* x112 - x146*x149))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.682784063255296d0
      x2 = 7.26848237132856d0*x1
      x3 = rho_a**(-2.66666666666667d0)
      x4 = 4.60115111447049d0
      x5 = 1d0/x4
      x6 = x3*x5
      x7 = gamma_aa*x6
      x8 = mu**2
      x9 = 1d0/x8
      x10 = 0.239205142311446d0*x9
      x11 = rho_a**0.666666666666667d0
      x12 = x11*x4
      x13 = 0.183440402716368d0*x12
      x14 = 21.1705915781930d0
      x15 = x14*x9
      x16 = 0.010095114404623d0*x0*x15
      x17 = rho_a**2.0d0
      x18 = 97.4090910340024d0
      x19 = x18/mu**4
      x20 = 0.147857142857143d0*x17*x19
      x21 = exp(-1.4385538026587772d0*x5*x8/x11)
      x22 = x21*(-x13 + x16 + x20)
      x23 = x10*x22*x7
      x24 = -x23 + 1.0d0
      x25 = 1.804d0 - 0.804d0/x24
      x26 = 2.14502939711103d0
      x27 = 1d0/x26
      x28 = rho_a**0.333333333333333d0
      x29 = 1d0/x28
      x30 = x27*x29
      x31 = 0.550321208149104d0*x30
      x32 = rho_a**1.0d0
      x33 = 1d0/x32
      x34 = 0.00844343197019481d0
      x35 = x33*x34
      x36 = x35*x8
      x37 = exp(-3.3019272488946267d0*x12*x9)
      x38 = mu*x37
      x39 = 1d0/mu
      x40 = x26*x39
      x41 = mu**3
      x42 = mu*x30
      x43 = x35*x41 - 0.825481812223657d0*x42 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x28*x40)
      x44 = 0.275160604074552d0*x42
      x45 = x44*(x38*(x31 - x36) + x43) - 0.375d0
      x46 = x25*x45
      x47 = exp(-4.1274090611182839d0*x42)
      x48 = 0.0206714418186455d0*x7 + 1.0d0
      x49 = 0.804d0/x48
      x50 = 0.0938662444277524d0*x7 + 1.0d0
      x51 = x50**(-2)
      x52 = 0.050475572023115d0*x7
      x53 = tau_a/gamma_aa
      x54 = 8.0d0*rho_a*x53
      x55 = x54 - 1.0d0
      x56 = x55*x7
      x57 = 0.0567850185260044d0*x56 - 0.45d0
      x58 = 0.126188930057787d0*x55
      x59 = x58*x7 - 1.0d0
      x60 = x55*x59
      x61 = x52*x60 + 1.0d0
      x62 = x61**(-0.5d0)
      x63 = x52 + x57*x62
      x64 = x63**2
      x65 = 0.00909806606836394d0*x6
      x66 = gamma_aa*x65
      x67 = x62*(0.0102353243269094d0*x56 - 0.0811111111111111d0) + x66
      x68 = 1d0/x14
      x69 = rho_a**(-5.33333333333333d0)*x68
      x70 = rho_a**(-2)
      x71 = tau_a**(-2)
      x72 = 0.0028125d0*x71
      x73 = sqrt(0.00286625629244325d0*x69 + x70*x72)
      x74 = abs(gamma_aa)
      x75 = x73*x74
      x76 = x67*x75
      x77 = gamma_aa**2
      x78 = x69*x77
      x79 = rho_a**(-8.0d0)
      x80 = 1d0/x18
      x81 = gamma_aa**3
      x82 = x80*x81
      x83 = x79*x82
      x84 = x71*x77
      x85 = x70*x84
      x86 = 0.015625d0*x85 + 1.0d0
      x87 = x86**(-2)
      x88 = 0.02485875d0*x85*x87 + 0.123456790123457d0
      x89 = 0.0941708433267071d0*x88
      x90 = x7*x89 + 0.000135164574498753d0*x78 + 0.000182133424673508d0
     &*x83 + 0.0021416494625995d0*x85
      x91 = x51*(0.0896750813832074d0*x64 - 1.24378109452736d0*x76 + x90
     &) + 1.0d0
      x92 = -0.804d0/x91
      x93 = x47*(x49 + x92)
      x94 = 0.682784063255296d0
      x95 = 2.72568088924821d0*x94
      x96 = -x31 + x36
      x97 = -x38*x96 + x43
      x98 = -x54 + 1.0d0
      x99 = 0.0567850185260044d0*x98
      x100 = x7*x99 + 0.45d0
      x101 = x7*x98
      x102 = 0.126188930057787d0*x101 + 1.0d0
      x103 = x102*x52*x98 + 1.0d0
      x104 = x103**(-0.5d0)
      x105 = -x100*x104 + x52
      x106 = x105**2
      x107 = -x104*(0.0102353243269094d0*x101 + 0.0811111111111111d0) + 
     &x66
      x108 = x107*x75
      x109 = x51*(0.0896750813832074d0*x106 - 1.24378109452736d0*x108 + 
     &x90) + 1.0d0
      x110 = x47*x95
      x111 = rho_a**(-0.666666666666667d0)
      x112 = 1d0/x0
      x113 = x112*x27
      x114 = x34/x17
      x115 = rho_a**(-0.333333333333333d0)*x4
      x116 = 0.318309886183791d0*mu
      x117 = x13 - x16 - x20
      x118 = rho_a**(-4.33333333333333d0)*x68
      x119 = rho_a**(-3.66666666666667d0)
      x120 = x45/x24**2
      x121 = 0.148394183600699d0
      x122 = x121/x48**2
      x123 = 0.90856029641607d0*x111*x94
      x124 = rho_a**(-6.33333333333333d0)*x68
      x125 = rho_a**(-3)
      x126 = gamma_aa*x119
      x127 = x126*x5
      x128 = 0.001553671875d0*x85/x86 - 0.0497175d0
      x129 = 0.0941708433267071d0*x128*x5*x87
      x130 = 1d0/x73
      x131 = 1d0/x74
      x132 = 1.24378109452736d0*x131
      x133 = 1d0/x50
      x134 = 0.0757133580346725d0*x7*x88 + 0.000108672317896998d0*x78 + 
     &0.0001464352734375d0*x83 + 0.00172188616792999d0*x85
      x135 = tau_a*x3
      x136 = x126*x55
      x137 = x62*(0.454280148208035d0*x135 - 0.151426716069345d0*x136)
      x138 = 0.180246913580247d0*x5
      x139 = x5*x59
      x140 = 0.0252377860115575d0*gamma_aa*x3*x55*x68
      x141 = x57*x61**(-1.5d0)*(-0.20190228809246d0*x135*x139 + 0.067300
     &7626974867d0* x136*x139 - x140*(1.0095114404623d0*x135 - 0.3365038
     &13487433d0* x136))
      x142 = 1.24378109452736d0*x75
      x143 = 2.0d0*x5
      x144 = x28*x51*x94
      x145 = rho_a**(-2.33333333333333d0)
      x146 = x109**(-2)
      x147 = rho_a**(-4.66666666666667d0)
      x148 = rho_a**(-1.66666666666667d0)
      x149 = x148*x53
      x150 = x104*(0.454280148208035d0*x149 + x3*x99)
      x151 = x148*x5
      x152 = 0.20190228809246d0*x151
      x153 = x103**(-1.5d0)*(x140*(-1.0095114404623d0*x149 + x3*x58) - x
     &152*x53*x59 + 0.0252377860115575d0*x6*x60)
      x154 = 0.180246913580247d0*x100
      x155 = 2.0d0*x100
      x156 = tau_a**(-3)
      x157 = x156*x70*x77
      x158 = (0.0254778337106067d0*gamma_aa*x118*x98 + x102*x152)/x103
      Ea = -x0*(-x2*x46 + x93*x95)
      d1Ea(1) = -rho_a*(-5.84385982654816d0*gamma_aa*x1*x120*x21*x28*(x1
     &0*x6*( -0.0134601525394973d0*rho_a**0.333333333333333d0*x15 + 0.12
     &2293601810912d0*x115 - 0.295714285714286d0*x19*x32) + 0.2294063113
     &91776d0*x117*x118 - 0.637880379497189d0*x117*x119*x5* x9) - 2.4228
     &2745710952d0*x1*x111*x46 + 7.26848237132856d0*x116* x25*x28*(0.091
     &7202013581841d0*x112*x97 - 0.275160604074552d0*x29* (0.27516060407
     &4552d0*mu*x113 + 1.21141372855476d0*x111*x37*x40 - x114*x41 + 2.20
     &128483259642d0*x115*x37*x39*x96 + x38*( -0.183440402716368d0*x113 
     &+ x114*x8))) + 3.75d0*x116*x33*x93 + x47*(0.12080100840081d0*gamma
     &_aa*rho_a**(-3.33333333333333d0)* x122 - x123*(-x49 + 1.804d0) + x
     &123*(x92 + 1.804d0) + 2.19144743495556d0*x144*(-0.0014570673973880
     &6d0*rho_a**(-9.0d0)* x82 + rho_a**(-5.66666666666667d0)*x129*x71*x
     &81 - 0.000720877730660018d0*x124*x77 - 0.00428329892519899d0*x125*
     &x84 + 0.622661654578789d0*x127*x133*(x134 + 0.0720987654320988d0*x
     &64 - x76) - 0.251122248871219d0*x127*x88 + x130*x132*x67*x77*( 0.0
     &0764335011318201d0*x124 + x125*x72) - x142*( -0.0242615095156372d0
     &*x127 + x137*x138 + 0.180246913580247d0*x141 ) + 0.089675081383207
     &4d0*x63*(-0.269203050789947d0*x127 + x137* x143 + 2.0d0*x141))/x91
     &**2)) - x110*x28*(x49 - 0.804d0/x109) + x2 *x28*(1.804d0 + 0.804d0
     &/(x23 - 1.0d0))*(x44*x97 - 0.375d0)
      d1Ea(2) = -rho_a*(1.39788132145759d0*x120*x121*x145*x22*x9 - 2.0d0
     &*x47*( 0.0226501890751519d0*x122*x145 + 1.09572371747778d0*x144*x1
     &46*( gamma_aa*x107*x132*x73 - 0.000270329148997507d0*gamma_aa*x69 
     &- 0.00428329892519899d0*gamma_aa*x70*x71 - 0.0896750813832074d0* x
     &105*(-x143*x150 + x153*x155 + 0.10095114404623d0*x6) + x129*x147 *
     &x84 + 0.233498120467046d0*x133*x6*(0.0720987654320988d0*x106 - x10
     &8 + x134) + x142*(-x138*x150 + x153*x154 + x65) - x6*x89 - 0.00054
     &6400274020523d0*x77*x79*x80)))
      d1Ea(3) = -x0*x110*x146*x51*(0.0720987654320988d0*x104*x105*(0.908
     &56029641607d0* x151 - x155*x158) - x104*x75*(0.0818825946152754d0*
     &x151 - x154* x158) + 0.0028125d0*x107*x130*x131*x157 + 0.075713358
     &0346725d0* x128*x147*x156*x5*x81*x87 - 0.00344377233585999d0*x157)
     &
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x138, x139, x140, x141, x142, x143, x144, x145, x146, x147, x148
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164, x165, x166, x167, x168, x169, x170
     &, x171
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.682784063255296d0
      x2 = 7.26848237132856d0*x1
      x3 = 1d0/x0
      x4 = 0.0472353356922751d0
      x5 = rho_a**0.666666666666667d0
      x6 = 4.60115111447049d0
      x7 = mu**2
      x8 = 1d0/x7
      x9 = x6*x8
      x10 = 3.30192724889463d0*x5*x9
      x11 = exp(x10)
      x12 = mu**4*x4*(x11 - 1.0d0)
      x13 = x12*x3
      x14 = 0.015286700226364d0*x13
      x15 = 7.0d0*x11 + 11.0d0
      x16 = 1d0/x6
      x17 = x16*x7/x5
      x18 = 0.00280419844572861d0*x15*x17
      x19 = -x14 + x18 + 0.0185185185185185d0
      x20 = 3.0d0*x11
      x21 = x20 - 2.0d0
      x22 = 2.14502939711103d0
      x23 = 1d0/x22
      x24 = rho_a**0.333333333333333d0
      x25 = 1d0/x24
      x26 = x23*x25
      x27 = mu*x26
      x28 = 1d0/mu
      x29 = x22*x28
      x30 = erf(1.8171205928321397d0*x24*x29)
      x31 = x11*x30
      x32 = 1d0/(0.183440402716368d0*x13 - 0.60570686427738d0*x17*x21 - 
     &x20 + 3.90167577847542d0*x27*x31)
      x33 = rho_a**(-2.0d0)
      x34 = 3.15935190989552d0*x33*x8
      x35 = x32*x34
      x36 = exp(-1.4385538026587772d0*x17)
      x37 = gamma_aa*x36
      x38 = x19*x35*x37
      x39 = -x38 + 1.0d0
      x40 = 1.804d0 - 0.804d0/x39
      x41 = 0.550321208149104d0*x26
      x42 = 1d0/rho_a
      x43 = 0.00844343197019481d0
      x44 = x42*x43
      x45 = x44*x7
      x46 = exp(-x10)
      x47 = mu*x46
      x48 = mu**3
      x49 = -0.825481812223657d0*x27 + 1.77245385090552d0*x30 + x44*x48
      x50 = 0.275160604074552d0*x27
      x51 = x50*(x47*(x41 - x45) + x49) - 0.375d0
      x52 = x40*x51
      x53 = exp(-4.1274090611182839d0*x27)
      x54 = rho_a**(-2.66666666666667d0)
      x55 = x16*x54
      x56 = gamma_aa*x55
      x57 = 0.0206714418186455d0*x56 + 1.0d0
      x58 = 0.804d0/x57
      x59 = 0.0938662444277524d0*x56 + 1.0d0
      x60 = x59**(-2)
      x61 = 0.050475572023115d0*x56
      x62 = tau_a/gamma_aa
      x63 = 8.0d0*rho_a*x62
      x64 = x63 - 1.0d0
      x65 = x56*x64
      x66 = 0.0567850185260044d0*x65 - 0.45d0
      x67 = 0.126188930057787d0*x65 - 1.0d0
      x68 = x61*x64*x67 + 1.0d0
      x69 = x68**(-0.5d0)
      x70 = x61 + x66*x69
      x71 = x70**2
      x72 = 0.00909806606836394d0*x55
      x73 = gamma_aa*x72
      x74 = x69*(0.0102353243269094d0*x65 - 0.0811111111111111d0) + x73
      x75 = rho_a**(-5.33333333333333d0)*x4
      x76 = rho_a**(-2)
      x77 = tau_a**(-2)
      x78 = 0.0028125d0*x77
      x79 = sqrt(0.00286625629244325d0*x75 + x76*x78)
      x80 = abs(gamma_aa)
      x81 = x79*x80
      x82 = x74*x81
      x83 = gamma_aa**2
      x84 = x75*x83
      x85 = rho_a**(-8.0d0)
      x86 = 0.0102659822546843d0
      x87 = gamma_aa**3
      x88 = x86*x87
      x89 = x85*x88
      x90 = x77*x83
      x91 = x76*x90
      x92 = 0.015625d0*x91 + 1.0d0
      x93 = x92**(-2)
      x94 = 0.02485875d0*x91*x93 + 0.123456790123457d0
      x95 = 0.0941708433267071d0*x94
      x96 = x56*x95 + 0.000135164574498753d0*x84 + 0.000182133424673508d
     &0*x89 + 0.0021416494625995d0*x91
      x97 = x60*(0.0896750813832074d0*x71 - 1.24378109452736d0*x82 + x96
     &) + 1.0d0
      x98 = -0.804d0/x97
      x99 = x53*(x58 + x98)
      x100 = 0.682784063255296d0
      x101 = 2.72568088924821d0*x100
      x102 = x38 - 1.0d0
      x103 = -x41 + x45
      x104 = -x103*x47 + x49
      x105 = x104*x50 - 0.375d0
      x106 = -x63 + 1.0d0
      x107 = x106*x56
      x108 = 0.0567850185260044d0*x107 + 0.45d0
      x109 = 0.126188930057787d0*x107 + 1.0d0
      x110 = x106*x109
      x111 = x110*x61 + 1.0d0
      x112 = x111**(-0.5d0)
      x113 = -x108*x112 + x61
      x114 = x113**2
      x115 = -x112*(0.0102353243269094d0*x107 + 0.0811111111111111d0) + 
     &x73
      x116 = x115*x81
      x117 = x60*(0.0896750813832074d0*x114 - 1.24378109452736d0*x116 + 
     &x96) + 1.0d0
      x118 = x101*x53
      x119 = rho_a**(-0.666666666666667d0)
      x120 = 0.318309886183791d0*mu
      x121 = x23*x3
      x122 = mu*x121
      x123 = x33*x43
      x124 = x119*x29
      x125 = rho_a**(-0.333333333333333d0)
      x126 = rho_a**(-3.66666666666667d0)
      x127 = x14 - x18 - 0.0185185185185185d0
      x128 = x11*x42
      x129 = rho_a**(-1.66666666666667d0)
      x130 = x129*x16
      x131 = x130*x7
      x132 = x11*x131
      x133 = rho_a**(-2.33333333333333d0)
      x134 = x12*x133
      x135 = x1*x32
      x136 = 0.148394183600699d0/x57**2
      x137 = 0.90856029641607d0*x100*x119
      x138 = rho_a**(-6.33333333333333d0)*x4
      x139 = rho_a**(-3)
      x140 = gamma_aa*x126
      x141 = x140*x16
      x142 = 0.001553671875d0*x91/x92 - 0.0497175d0
      x143 = 0.0941708433267071d0*x142*x16*x93
      x144 = 1d0/x79
      x145 = 1d0/x80
      x146 = 1.24378109452736d0*x145
      x147 = 1d0/x59
      x148 = 0.0757133580346725d0*x56*x94 + 0.000108672317896998d0*x84 +
     & 0.0001464352734375d0*x89 + 0.00172188616792999d0*x91
      x149 = tau_a*x54
      x150 = x140*x64
      x151 = x69*(0.454280148208035d0*x149 - 0.151426716069345d0*x150)
      x152 = 0.180246913580247d0*x16
      x153 = x16*x67
      x154 = gamma_aa*x4
      x155 = 0.0252377860115575d0*x154
      x156 = x66*x68**(-1.5d0)*(-0.20190228809246d0*x149*x153 + 0.067300
     &7626974867d0* x150*x153 - x155*x54*x64*(1.0095114404623d0*x149 - 0
     &.336503813487433d0*x150))
      x157 = 1.24378109452736d0*x81
      x158 = 2.0d0*x16
      x159 = x100*x24*x60
      x160 = x117**(-2)
      x161 = rho_a**(-4.66666666666667d0)
      x162 = x129*x62
      x163 = x106*x54
      x164 = x112*(0.454280148208035d0*x162 + 0.0567850185260044d0*x163)
     &
      x165 = 0.20190228809246d0*x109*x130
      x166 = x111**(-1.5d0)*(0.0252377860115575d0*x110*x55 + x155*x163*(
     & 1.0095114404623d0*x162 + 0.126188930057787d0*x163) + x165*x62)
      x167 = 0.180246913580247d0*x108
      x168 = 2.0d0*x108
      x169 = tau_a**(-3)
      x170 = x169*x76*x83
      x171 = (0.0254778337106067d0*rho_a**(-4.33333333333333d0)*x106*x15
     &4 + x165)/ x111
      Ea = -x0*(x101*x99 - x2*x52)
      d1Ea(1) = -rho_a*(-2.42282745710952d0*x1*x119*x52 + 7.268482371328
     &56d0*x120*x24* x40*(0.0917202013581841d0*x104*x3 - 0.2751606040745
     &52d0*x25*( 2.20128483259642d0*x103*x125*x28*x46*x6 + 0.27516060407
     &4552d0* x122 - x123*x48 + 1.21141372855476d0*x124*x46 + x47*( -0.1
     &83440402716368d0*x121 + x123*x7))) + 3.75d0*x120*x42*x99 - 5.84385
     &982654816d0*x135*x24*x37*x51*(-6.31870381979104d0*rho_a**( -3.0d0)
     &*x127*x8 + 3.02993180261165d0*x126*x127*x16 + x127*x35*( 6.6038544
     &9778925d0*x11*x125*x9 + 1.30055859282514d0*x122*x31 - 8.5886997128
     &6677d0*x124*x31 + 4.0d0*x128 - 0.40380457618492d0* x131*x21 - 0.40
     &380457618492d0*x132 + 0.244587203621824d0*x134 - 2.66666666666667d
     &0*x42) + x34*(-0.0432098765432099d0*x128 + 0.00186946563048574d0*x
     &131*x15 + 0.0336503813487433d0*x132 - 0.0203822669684853d0*x134))/
     &x39**2 + x53*(0.12080100840081d0* gamma_aa*rho_a**(-3.333333333333
     &33d0)*x136 - x137*(-x58 + 1.804d0 ) + x137*(x98 + 1.804d0) + 2.191
     &44743495556d0*x159*( -0.00145706739738806d0*rho_a**(-9.0d0)*x88 + 
     &rho_a**( -5.66666666666667d0)*x143*x77*x87 - 0.000720877730660018d
     &0*x138* x83 - 0.00428329892519899d0*x139*x90 + 0.622661654578789d0
     &*x141* x147*(x148 + 0.0720987654320988d0*x71 - x82) - 0.2511222488
     &71219d0*x141*x94 + x144*x146*x74*x83*( 0.00764335011318201d0*x138 
     &+ x139*x78) - x157*( -0.0242615095156372d0*x141 + x151*x152 + 0.18
     &0246913580247d0*x156 ) + 0.0896750813832074d0*x70*(-0.269203050789
     &947d0*x141 + x151* x158 + 2.0d0*x156))/x97**2)) + x105*x2*x24*(1.8
     &04d0 + 0.804d0/ x102) - x118*x24*(x58 - 0.804d0/x117)
      d1Ea(2) = rho_a*(-18.4628097041666d0*rho_a**(-1.66666666666667d0)*
     &x105*x135*x19* x36*x8/x102**2 + 2.0d0*x53*(0.0226501890751519d0*x1
     &33*x136 + 1.09572371747778d0*x159*x160*(gamma_aa*x115*x146*x79 - 0
     &.000270329148997507d0*gamma_aa*x75 - 0.00428329892519899d0* gamma_
     &aa*x76*x77 - 0.0896750813832074d0*x113*(-x158*x164 + x166* x168 + 
     &0.10095114404623d0*x55) + x143*x161*x90 + 0.233498120467046d0*x147
     &*x55*(0.0720987654320988d0*x114 - x116 + x148) + x157*(-x152*x164 
     &+ x166*x167 + x72) - x55*x95 - 0.000546400274020523d0*x83*x85*x86)
     &))
      d1Ea(3) = -x0*x118*x160*x60*(0.0720987654320988d0*x112*x113*(0.908
     &56029641607d0* x130 - x168*x171) - x112*x81*(0.0818825946152754d0*
     &x130 - x167* x171) + 0.0028125d0*x115*x144*x145*x170 + 0.075713358
     &0346725d0* x142*x16*x161*x169*x87*x93 - 0.00344377233585999d0*x170
     &)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x116, x117, x118, x119, x120, x121, x122, x123, x124, x125
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      x0 = rho_a**1.0d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = 4.60115111447049d0
      x3 = 1d0/x2
      x4 = x1*x3
      x5 = gamma_aa*x4
      x6 = mu**2
      x7 = 1d0/x6
      x8 = 0.239205142311446d0*x7
      x9 = rho_a**0.666666666666667d0
      x10 = 0.183440402716368d0*x2*x9
      x11 = rho_a**1.33333333333333d0
      x12 = 21.1705915781930d0
      x13 = x12*x7
      x14 = 0.010095114404623d0*x11*x13
      x15 = 97.4090910340024d0
      x16 = mu**(-4)
      x17 = x15*x16
      x18 = 0.147857142857143d0*rho_a**2.0d0*x17
      x19 = exp(-1.4385538026587772d0*x3*x6/x9)
      x20 = x19*(-x10 + x14 + x18)
      x21 = x20*x5*x8
      x22 = -x21 + 1.0d0
      x23 = 3.14159265358979d0
      x24 = x23*x7
      x25 = 1.0d0*x24*(1.804d0 - 0.804d0/x22)
      x26 = 0.682784063255296d0
      x27 = rho_a**0.333333333333333d0
      x28 = x26*x27
      x29 = exp(-1.9241736577956332d0*mu/x27)
      x30 = 0.0206714418186455d0*x5 + 1.0d0
      x31 = 0.804d0/x30
      x32 = 0.0938662444277524d0*x5 + 1.0d0
      x33 = x32**(-2)
      x34 = 0.050475572023115d0*x5
      x35 = tau_a/gamma_aa
      x36 = 8.0d0*rho_a*x35
      x37 = x36 - 1.0d0
      x38 = x37*x5
      x39 = 0.0567850185260044d0*x38 - 0.45d0
      x40 = 0.126188930057787d0*x38 - 1.0d0
      x41 = x34*x37*x40 + 1.0d0
      x42 = x41**(-0.5d0)
      x43 = x34 + x39*x42
      x44 = x43**2
      x45 = 0.00909806606836394d0*x4
      x46 = gamma_aa*x45
      x47 = x42*(0.0102353243269094d0*x38 - 0.0811111111111111d0) + x46
      x48 = 1d0/x12
      x49 = rho_a**(-5.33333333333333d0)*x48
      x50 = rho_a**(-2)
      x51 = tau_a**(-2)
      x52 = 0.0028125d0*x51
      x53 = sqrt(0.00286625629244325d0*x49 + x50*x52)
      x54 = abs(gamma_aa)
      x55 = x53*x54
      x56 = x47*x55
      x57 = gamma_aa**2
      x58 = x49*x57
      x59 = rho_a**(-8.0d0)
      x60 = 1d0/x15
      x61 = gamma_aa**3
      x62 = x60*x61
      x63 = x59*x62
      x64 = x51*x57
      x65 = x50*x64
      x66 = 0.015625d0*x65 + 1.0d0
      x67 = x66**(-2)
      x68 = 0.02485875d0*x65*x67 + 0.123456790123457d0
      x69 = 0.0941708433267071d0*x68
      x70 = x5*x69 + 0.000135164574498753d0*x58 + 0.000182133424673508d0
     &*x63 + 0.0021416494625995d0*x65
      x71 = x33*(0.0896750813832074d0*x44 - 1.24378109452736d0*x56 + x70
     &) + 1.0d0
      x72 = -0.804d0/x71
      x73 = x29*(x31 + x72)
      x74 = x0*x25 + 2.72568088924821d0*x28*x73
      x75 = x10 - x14 - x18
      x76 = rho_a**(-4.33333333333333d0)*x48
      x77 = rho_a**(-3.66666666666667d0)
      x78 = 0.148394183600699d0/x30**2
      x79 = 0.90856029641607d0*rho_a**(-0.666666666666667d0)*x26
      x80 = rho_a**(-6.33333333333333d0)*x48
      x81 = rho_a**(-3)
      x82 = gamma_aa*x77
      x83 = x3*x82
      x84 = 0.001553671875d0*x65/x66 - 0.0497175d0
      x85 = 0.0941708433267071d0*x3*x67*x84
      x86 = 1d0/x53
      x87 = 1d0/x54
      x88 = 1.24378109452736d0*x87
      x89 = 1d0/x32
      x90 = 0.0757133580346725d0*x5*x68 + 0.000108672317896998d0*x58 + 0
     &.0001464352734375d0*x63 + 0.00172188616792999d0*x65
      x91 = tau_a*x1
      x92 = x37*x82
      x93 = x42*(0.454280148208035d0*x91 - 0.151426716069345d0*x92)
      x94 = 0.180246913580247d0*x3
      x95 = x3*x40
      x96 = 0.0252377860115575d0*gamma_aa*x48
      x97 = x39*x41**(-1.5d0)*(-x1*x37*x96*(1.0095114404623d0*x91 - 0.33
     &6503813487433d0*x92) - 0.20190228809246d0*x91*x95 + 0.067300762697
     &4867d0*x92*x95)
      x98 = 1.24378109452736d0*x55
      x99 = 2.0d0*x3
      x100 = x28*x33
      x101 = rho_a**(-1.66666666666667d0)
      x102 = -x36 + 1.0d0
      x103 = x102*x5
      x104 = 0.0567850185260044d0*x103 + 0.45d0
      x105 = 0.126188930057787d0*x103 + 1.0d0
      x106 = x102*x105
      x107 = x106*x34 + 1.0d0
      x108 = x107**(-0.5d0)
      x109 = -x104*x108 + x34
      x110 = x109**2
      x111 = -x108*(0.0102353243269094d0*x103 + 0.0811111111111111d0) + 
     &x46
      x112 = x111*x55
      x113 = (x33*(0.0896750813832074d0*x110 - 1.24378109452736d0*x112 +
     & x70) + 1.0d0 )**(-2)
      x114 = rho_a**(-4.66666666666667d0)
      x115 = x101*x35
      x116 = x1*x102
      x117 = x108*(0.454280148208035d0*x115 + 0.0567850185260044d0*x116)
     &
      x118 = x101*x3
      x119 = 0.20190228809246d0*x105*x118
      x120 = x107**(-1.5d0)*(0.0252377860115575d0*x106*x4 + x116*x96*( 1
     &.0095114404623d0*x115 + 0.126188930057787d0*x116) + x119*x35)
      x121 = 0.180246913580247d0*x104
      x122 = 2.0d0*x104
      x123 = tau_a**(-3)
      x124 = x123*x50*x57
      x125 = (0.0254778337106067d0*gamma_aa*x102*x76 + x119)/x107
      Ea = -rho_a*x74
      d1Ea(1) = -rho_a*(0.804d0*gamma_aa*x0*x19*x24*(-0.637880379497189d
     &0*x3*x7*x75*x77 + x4*x8*(0.122293601810912d0*rho_a**(-0.3333333333
     &33333d0)*x2 - 0.0134601525394973d0*rho_a**0.333333333333333d0*x13 
     &- 0.295714285714286d0*x0*x17) + 0.229406311391776d0*x75*x76)/x22**
     &2 + 3.75d0*mu*x73/(x0*x23) + x25 + x29*(0.12080100840081d0*gamma_a
     &a *rho_a**(-3.33333333333333d0)*x78 + 2.19144743495556d0*x100*( -0
     &.00145706739738806d0*rho_a**(-9.0d0)*x62 + rho_a**( -5.66666666666
     &667d0)*x51*x61*x85 + 0.0896750813832074d0*x43*( -0.269203050789947
     &d0*x83 + x93*x99 + 2.0d0*x97) + x47*x57*x86*x88 *(x52*x81 + 0.0076
     &4335011318201d0*x80) - 0.000720877730660018d0* x57*x80 - 0.0042832
     &9892519899d0*x64*x81 - 0.251122248871219d0*x68 *x83 + 0.6226616545
     &78789d0*x83*x89*(0.0720987654320988d0*x44 - x56 + x90) - x98*(-0.0
     &242615095156372d0*x83 + x93*x94 + 0.180246913580247d0*x97))/x71**2
     & - x79*(-x31 + 1.804d0) + x79*( x72 + 1.804d0))) - x74
      d1Ea(2) = rho_a*(0.131313669051252d0*x101*x16*x20/(x21 - 1.0d0)**2
     & + 2.0d0*x29*( 0.0226501890751519d0*rho_a**(-2.33333333333333d0)*x
     &78 + 1.09572371747778d0*x100*x113*(gamma_aa*x111*x53*x88 - 0.00027
     &0329148997507d0*gamma_aa*x49 - 0.00428329892519899d0* gamma_aa*x50
     &*x51 - 0.0896750813832074d0*x109*(-x117*x99 + x120* x122 + 0.10095
     &114404623d0*x4) + x114*x64*x85 - x4*x69 + 0.233498120467046d0*x4*x
     &89*(0.0720987654320988d0*x110 - x112 + x90) - 0.000546400274020523
     &d0*x57*x59*x60 + x98*(-x117*x94 + x120 *x121 + x45))))
      d1Ea(3) = -2.72568088924821d0*x11*x113*x26*x29*x33*(0.072098765432
     &0988d0*x108*x109 *(0.90856029641607d0*x118 - x122*x125) - x108*x55
     &*( 0.0818825946152754d0*x118 - x121*x125) + 0.0028125d0*x111*x124*
     & x86*x87 + 0.0757133580346725d0*x114*x123*x3*x61*x67*x84 - 0.00344
     &377233585999d0*x124)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x402, x403, x404
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.682784063255296d0
      x2 = rho_a**(-2.66666666666667d0)
      x3 = 0.217336917462899d0
      x4 = x2*x3
      x5 = gamma_aa*x4
      x6 = mu**2
      x7 = rho_a**(-0.666666666666667d0)*x3*x6
      x8 = 5.45136177849642d0*x7 - 7.0d0
      x9 = 0.239205142311446d0*x8
      x10 = 0.826307487110758d0
      x11 = rho_a**0.333333333333333d0
      x12 = mu/x11
      x13 = 0.0472353356922751d0
      x14 = mu**4
      x15 = 1d0/x0
      x16 = 59.4346904801033d0*x10*x12 + 4.95289087334194d0*x13*x14*x15 
     &- 49.0622560064678d0*x7 - 81.0d0
      x17 = 1d0/x16
      x18 = x17*exp(-1.4385538026587772d0*x7)
      x19 = x18*x5*x9 + 1.0d0
      x20 = x1*(1.804d0 - 0.804d0/x19)
      x21 = 2.72568088924821d0*x20
      x22 = 0.682784063255296d0
      x23 = 2.72568088924821d0*x22
      x24 = 0.466194077035412d0
      x25 = exp(-4.1274090611182839d0*x12*x24)
      x26 = 0.0206714418186455d0*x5 + 1.0d0
      x27 = 0.804d0/x26
      x28 = 0.0938662444277524d0*x5 + 1.0d0
      x29 = x28**(-2)
      x30 = 0.050475572023115d0*x5
      x31 = 1d0/gamma_aa
      x32 = tau_a*x31
      x33 = 8.0d0*rho_a*x32
      x34 = x33 - 1.0d0
      x35 = x34*x5
      x36 = 0.0567850185260044d0*x35 - 0.45d0
      x37 = 0.126188930057787d0*x35 - 1.0d0
      x38 = x34*x37
      x39 = x30*x38 + 1.0d0
      x40 = x39**(-0.5d0)
      x41 = x30 + x36*x40
      x42 = x41**2
      x43 = 0.00909806606836394d0*x4
      x44 = gamma_aa*x43
      x45 = x40*(0.0102353243269094d0*x35 - 0.0811111111111111d0) + x44
      x46 = rho_a**(-5.33333333333333d0)*x13
      x47 = rho_a**(-2)
      x48 = tau_a**(-2)
      x49 = 0.0028125d0*x48
      x50 = 0.00286625629244325d0*x46 + x47*x49
      x51 = sqrt(x50)
      x52 = abs(gamma_aa)
      x53 = x51*x52
      x54 = x45*x53
      x55 = 1.24378109452736d0*x54
      x56 = gamma_aa**2
      x57 = x46*x56
      x58 = rho_a**(-8.0d0)
      x59 = 0.0102659822546843d0
      x60 = gamma_aa**3
      x61 = x59*x60
      x62 = x58*x61
      x63 = x47*x48
      x64 = x56*x63
      x65 = 0.015625d0*x64 + 1.0d0
      x66 = x65**(-2)
      x67 = 0.02485875d0*x64*x66 + 0.123456790123457d0
      x68 = x5*x67
      x69 = 0.000135164574498753d0*x57 + 0.000182133424673508d0*x62 + 0.
     &0021416494625995d0*x64 + 0.0941708433267071d0*x68
      x70 = x29*(0.0896750813832074d0*x42 - x55 + x69) + 1.0d0
      x71 = 1d0/x70
      x72 = -0.804d0*x71
      x73 = x25*(x27 + x72)
      x74 = x23*x73
      x75 = rho_a**(-0.666666666666667d0)
      x76 = x20*x75
      x77 = rho_a**(-4.33333333333333d0)
      x78 = x13*x77
      x79 = x6*x78
      x80 = rho_a**(-3.66666666666667d0)
      x81 = x3*x80
      x82 = x8*x81
      x83 = x79*x8
      x84 = x10*x15
      x85 = rho_a**(-1.66666666666667d0)
      x86 = x3*x85
      x87 = mu*x86
      x88 = rho_a**(-2.33333333333333d0)
      x89 = mu**3
      x90 = x13*x89
      x91 = x88*x90
      x92 = 19.8115634933678d0*x84 - 32.7081706709785d0*x87 + 6.60385449
     &778925d0*x91
      x93 = x17*x92
      x94 = mu*x4
      x95 = x9*x94
      x96 = -0.869329180010942d0*x79 - 0.637880379497189d0*x82 + 0.22940
     &6311391776d0 *x83 + x93*x95
      x97 = x18/x19**2
      x98 = x1*x96*x97
      x99 = 2.19144743495556d0*x11
      x100 = x98*x99
      x101 = 0.318309886183791d0*mu
      x102 = 1d0/rho_a*x101
      x103 = 3.75d0*x102
      x104 = rho_a**(-3.33333333333333d0)
      x105 = 0.148394183600699d0
      x106 = x105/x26**2
      x107 = x104*x106
      x108 = -x27 + 1.804d0
      x109 = 0.90856029641607d0*x75
      x110 = x109*x22
      x111 = x70**(-2)
      x112 = 2.19144743495556d0*x11
      x113 = rho_a**(-9.0d0)
      x114 = x113*x61
      x115 = rho_a**(-6.33333333333333d0)*x13
      x116 = x115*x56
      x117 = rho_a**(-3)
      x118 = x117*x48
      x119 = x118*x56
      x120 = gamma_aa*x81
      x121 = 0.251122248871219d0*x67
      x122 = x64/x65
      x123 = 0.001553671875d0*x122
      x124 = x123 - 0.0497175d0
      x125 = x124*x60
      x126 = x48*x66
      x127 = 0.0941708433267071d0*x126
      x128 = rho_a**(-5.66666666666667d0)*x3
      x129 = x127*x128
      x130 = 1d0/x51
      x131 = 1d0/x52
      x132 = x131*x45
      x133 = x130*x132
      x134 = 0.00764335011318201d0*x115 + x117*x49
      x135 = x134*x56
      x136 = x133*x135
      x137 = 1d0/x28
      x138 = 0.000108672317896998d0*x57 + 0.0001464352734375d0*x62 + 0.0
     &0172188616792999d0*x64 + 0.0757133580346725d0*x68
      x139 = x138 + 0.0720987654320988d0*x42 - x54
      x140 = x137*x139
      x141 = 0.622661654578789d0*x140
      x142 = 0.0242615095156372d0*x81
      x143 = tau_a*x2
      x144 = x34*x80
      x145 = 0.151426716069345d0*x144
      x146 = -gamma_aa*x145 + 0.454280148208035d0*x143
      x147 = x146*x40
      x148 = 0.180246913580247d0*x3
      x149 = x3*x37
      x150 = 0.20190228809246d0*x143
      x151 = 0.0673007626974867d0*x38
      x152 = 0.336503813487433d0*x144
      x153 = x13*(-gamma_aa*x152 + 1.0095114404623d0*x143)
      x154 = x2*x34
      x155 = 0.0252377860115575d0*x154
      x156 = x153*x155
      x157 = -gamma_aa*x156 + x120*x151 - x149*x150
      x158 = x39**(-1.5d0)
      x159 = x158*x36
      x160 = x157*x159
      x161 = -gamma_aa*x142 + x147*x148 + 0.180246913580247d0*x160
      x162 = x161*x53
      x163 = 0.269203050789947d0*x81
      x164 = 2.0d0*x3
      x165 = -gamma_aa*x163 + x147*x164 + 2.0d0*x160
      x166 = x165*x41
      x167 = -0.00145706739738806d0*x114 - 0.000720877730660018d0*x116 -
     & 0.00428329892519899d0*x119 - x120*x121 + x120*x141 + x125*x129 + 
     &1.24378109452736d0*x136 - 1.24378109452736d0*x162 + 0.089675081383
     &2074d0*x166
      x168 = x22*x29
      x169 = x167*x168
      x170 = x25*(0.12080100840081d0*gamma_aa*x107 - x108*x110 + x110*(x
     &72 + 1.804d0 ) + x111*x112*x169)
      x171 = x105*x88
      x172 = x171*x8
      x173 = 0.524205495546598d0*x172*x97
      x174 = x106*x88
      x175 = 0.0226501890751519d0*x174
      x176 = -x33 + 1.0d0
      x177 = x176*x5
      x178 = 0.0567850185260044d0*x177 + 0.45d0
      x179 = 0.126188930057787d0*x177 + 1.0d0
      x180 = x176*x179
      x181 = x180*x30 + 1.0d0
      x182 = x181**(-0.5d0)
      x183 = -x178*x182 + x30
      x184 = 2.0d0*x178
      x185 = x181**(-1.5d0)
      x186 = 0.20190228809246d0*x86
      x187 = x186*x37
      x188 = x187*x32
      x189 = x37*x4
      x190 = x189*x34
      x191 = 0.0252377860115575d0*x190
      x192 = x32*x85
      x193 = 1.0095114404623d0*x192
      x194 = 0.126188930057787d0*x2
      x195 = -x193 + x194*x34
      x196 = gamma_aa*x13
      x197 = x155*x196
      x198 = x195*x197
      x199 = -x188 + x191 + x198
      x200 = x185*x199
      x201 = 0.10095114404623d0*x4
      x202 = 0.454280148208035d0*x192
      x203 = 0.0567850185260044d0*x2
      x204 = x176*x203 + x202
      x205 = x182*x204
      x206 = -x164*x205 + x201
      x207 = x184*x200 + x206
      x208 = x183*x207
      x209 = 0.180246913580247d0*x178
      x210 = -x148*x205 + x43
      x211 = x200*x209 + x210
      x212 = 1.24378109452736d0*x53
      x213 = 0.000270329148997507d0*x46
      x214 = x4*x67
      x215 = 0.00428329892519899d0*x63
      x216 = x56*x59
      x217 = x216*x58
      x218 = -x182*(0.0102353243269094d0*x177 + 0.0811111111111111d0) + 
     &x44
      x219 = x131*x218
      x220 = 1.24378109452736d0*gamma_aa
      x221 = x220*x51
      x222 = x183**2
      x223 = x218*x53
      x224 = x138 + 0.0720987654320988d0*x222 - x223
      x225 = x137*x4
      x226 = 0.233498120467046d0*x225
      x227 = rho_a**(-4.66666666666667d0)
      x228 = x227*x3
      x229 = x228*x56
      x230 = x124*x229
      x231 = -gamma_aa*x213 - gamma_aa*x215 + x127*x230 - 0.094170843326
     &7071d0*x214 - 0.000546400274020523d0*x217 + x219*x221 + x224*x226
      x232 = -0.0896750813832074d0*x208 + x211*x212 + x231
      x233 = x168*x232
      x234 = x29*(0.0896750813832074d0*x222 - 1.24378109452736d0*x223 + 
     &x69) + 1.0d0
      x235 = x234**(-2)
      x236 = x11*x235
      x237 = 1.09572371747778d0*x236
      x238 = 2.0d0*x25
      x239 = 0.90856029641607d0*x86
      x240 = 1d0/x181
      x241 = 0.0254778337106067d0*x77
      x242 = x176*x196
      x243 = x179*x186
      x244 = x240*(x241*x242 + x243)
      x245 = 0.0720987654320988d0*x182*x183
      x246 = 0.0818825946152754d0*x86
      x247 = x182*x53
      x248 = tau_a**(-3)
      x249 = x248*x56
      x250 = x249*x47
      x251 = 0.0028125d0*x250
      x252 = x130*x219
      x253 = 0.0757133580346725d0*x125
      x254 = x248*x66
      x255 = x228*x253*x254 - 0.00344377233585999d0*x250
      x256 = x251*x252 + x255
      x257 = x245*(-x184*x244 + x239) - x247*(-x209*x244 + x246) + x256
      x258 = x25*x29
      x259 = x0*x235*x258
      x260 = gamma_aa*x98
      x261 = 0.60570686427738d0*rho_a**(-1.66666666666667d0)
      x262 = 1.27576075899438d0*x82
      x263 = x1*x99
      x264 = x16**(-2)
      x265 = x264*exp(-2.8771076053175544d0*x7)/x19**3
      x266 = x265*x96
      x267 = rho_a**(-6.0d0)
      x268 = x14*x267*x59
      x269 = x46*x6
      x270 = x78*x89*x93
      x271 = mu*x15*x24
      x272 = 0.0322515344331995d0
      x273 = x272/x26**3
      x274 = x22*x261
      x275 = 1d0/x234
      x276 = x235*x75
      x277 = x120*x67
      x278 = x126*x128
      x279 = 1.24532330915758d0*x120
      x280 = x22/(x234**3*x28**4)
      x281 = x112*x167*x280
      x282 = rho_a**(-7.33333333333333d0)*x13
      x283 = x282*x56
      x284 = rho_a**(-4)
      x285 = x284*x48
      x286 = gamma_aa*x227
      x287 = x286*x3
      x288 = rho_a**(-6.66666666666667d0)*x3
      x289 = tau_a**(-4)
      x290 = x289*x66
      x291 = 0.00014565673828125d0*gamma_aa**4*x284*x290
      x292 = x60*(-0.010875703125d0*x122 + x291 + 0.1491525d0)
      x293 = x50**(-1.5d0)*x52
      x294 = 1.24378109452736d0*x45
      x295 = x139*x29
      x296 = x131*x161
      x297 = x130*x135
      x298 = 0.0720987654320988d0*x3
      x299 = 1.24378109452736d0*x165
      x300 = tau_a*x80
      x301 = x286*x34
      x302 = x40*(-2.42282745710952d0*x300 + 0.555231292254265d0*x301)
      x303 = 0.360493827160494d0*x3
      x304 = x157*x158
      x305 = x146*x304
      x306 = gamma_aa*x153
      x307 = x159*(-0.40380457618492d0*x143*x153 + 0.134601525394973d0*x
     &144*x306 + 1.07681220315979d0*x149*x300 - 0.246769463224118d0*x149
     &*x301 - x197*(-5.38406101579893d0*x300 + 1.23384731612059d0*x301))
     &
      x308 = 0.60570686427738d0*x37
      x309 = 0.0757133580346725d0*x154
      x310 = x157*x36
      x311 = x310*x39**(-2.5d0)
      x312 = x311*(0.20190228809246d0*x120*x38 - x143*x3*x308 - x306*x30
     &9)
      x313 = 0.00363922642734557d0*x4
      x314 = gamma_aa*x313
      x315 = x314 + x40*(0.00409412973076377d0*x35 - 0.0324444444444444d
     &0)
      x316 = 1.24378109452736d0*x315
      x317 = -0.0011714821875d0*x114 - 0.000579585695450655d0*x116 - 0.0
     &0344377233585999d0*x119 + x136 - x162 + 0.0720987654320988d0* x166
     & + x253*x278 - 0.20190228809246d0*x277
      x318 = x168*x235
      x319 = x112*x318
      x320 = x185*(0.0252377860115575d0*x180*x4 + 0.0252377860115575d0*x
     &2*x242*(x176* x194 + x193) + x243*x32)
      x321 = -0.0896750813832074d0*x183*(x184*x320 + x206) + x212*(x209*
     &x320 + x210) + x231
      x322 = x168*x237
      x323 = gamma_aa*x46
      x324 = 0.00856659785039798d0*gamma_aa
      x325 = gamma_aa*x51
      x326 = 2.48756218905473d0*x325
      x327 = 0.466996240934092d0*x225
      x328 = 0.188341686653414d0*x126*x230 - 0.179350162766415d0*x208 + 
     &2.48756218905473d0*x211*x53 - 0.188341686653414d0*x214 - 0.0010928
     &0054804105d0*x217 + x219*x326 + x224*x327 - 0.000540658297995013d0
     &*x323 - x324*x63
      x329 = gamma_aa*x115
      x330 = -x123 + 0.0497175d0
      x331 = x129*x56
      x332 = 0.00932203125d0*x122
      x333 = x130*x31*x52
      x334 = -x202 + x203*x34
      x335 = x334*x40
      x336 = x188 - x191 - x198
      x337 = x159*x336
      x338 = x148*x335 + 0.180246913580247d0*x337 + x43
      x339 = x131*x338
      x340 = 1.24378109452736d0*x297
      x341 = x298*x335 + x313 + 0.0720987654320988d0*x337
      x342 = x126*x330
      x343 = x132*x51
      x344 = x164*x335 + x201 + 2.0d0*x337
      x345 = 0.0720987654320988d0*x41
      x346 = -gamma_aa*x343 + 0.00344377233585999d0*gamma_aa*x63 + 0.075
     &7133580346725d0*x214 + 0.0004393058203125d0*x217 + 0.0757133580346
     &725d0*x229*x342 + 0.000217344635793995d0*x323 - x338*x53 + x344*x3
     &45
      x347 = 0.622661654578789d0*x120*x137
      x348 = x143*x31
      x349 = x40*(-x145 + 1.21141372855476d0*x348)
      x350 = x304*x334
      x351 = x158*x336
      x352 = x146*x148
      x353 = x308*x86
      x354 = x195*x196
      x355 = -0.0757133580346725d0*x190 - x309*x354 + x32*x353
      x356 = x311*x355
      x357 = 0.20190228809246d0*x85
      x358 = x153*x357
      x359 = x159*(-x13*x150*x195 + 0.0673007626974867d0*x144*x354 - 0.5
     &38406101579893d0*x149*x348 + x151*x81 - x156 - x197*(-x152 + 2.692
     &03050789947d0*x348) + x32*x358)
      x360 = x146*x164
      x361 = rho_a**(-5.0d0)
      x362 = x13*x36
      x363 = x195*x362*(-0.050475572023115d0*x154 + 0.40380457618492d0*x
     &192)
      x364 = x334*x336
      x365 = 1d0/x39
      x366 = x36*x365
      x367 = x336*x355*x366
      x368 = x241*x34
      x369 = x187 + x196*x368
      x370 = x240*x369
      x371 = x184*x370 + x239
      x372 = x209*x370 + x246
      x373 = x245*x371 - x247*x372 + x256
      x374 = x366*x369
      x375 = 2.0d0*x374
      x376 = 0.180246913580247d0*x374
      x377 = x40*(x246 - x376)
      x378 = x133*x251 + x255 + x345*x40*(x239 - x375) - x377*x53
      x379 = x117*x249
      x380 = x128*x254
      x381 = x291 - x332 + 0.099435d0
      x382 = x248*x47
      x383 = x130*x250
      x384 = 0.0327530378461102d0*x86
      x385 = x365*x369
      x386 = x366*(-0.203822669684853d0*tau_a*x78 + 0.336503813487433d0*
     &x189 + 0.110403946079296d0*x323*x34 - x358)
      x387 = x157*x365
      x388 = 0.0764335011318201d0*gamma_aa*x34*x78 + x353
      x389 = x310*x388/x39**2
      x390 = gamma_aa*x382
      x391 = 0.103185226527957d0*x228*x249*x66
      x392 = 1.3628404446241d0*x182
      x393 = x131*x372
      x394 = 0.0720987654320988d0*x178*x370 + x384
      x395 = x204*x369
      x396 = x13*(-0.203822669684853d0*x104*x32 + x195*x357 + x368)
      x397 = x199*x240*x388
      x398 = 1.3628404446241d0*x185
      x399 = -x182*(0.00409412973076377d0*x177 + 0.0324444444444444d0) +
     & x314
      x400 = x275*x29
      x401 = x289*x47*x56
      x402 = x104*x362
      x403 = x369*x86
      x404 = 2.72568088924821d0*x158
      Ea = -x0*(x21 + x74)
      d1Ea(1) = -rho_a*(gamma_aa*x100 + x103*x73 + x170 + 0.908560296416
     &07d0*x76) - x11* x21 - x11*x74
      d1Ea(2) = -rho_a*(x173 - x238*(x175 + x233*x237))
      d1Ea(3) = -x23*x257*x259
      d2Ea(1) = rho_a*(-gamma_aa*x263*x97*(-mu*x262*x93 + x17*x95*(-26.4
     &15417991157d0* x10*x88 - 15.4089938281749d0*x104*x90 + 54.51361778
     &49642d0*x94) + 2.33889472482303d0*x228*x8 + x264*x4*x6*x9*x92*( 39
     &.6231269867355d0*x84 - 65.416341341957d0*x87 + 13.2077089955785d0*
     &x91) + 0.220008881071042d0*x268*x8 - 1.66743573022264d0*x268 - 1.6
     &0584417974243d0*x269*x8 + 6.08530426007659d0*x269 + 0.458812622783
     &553d0*x270*x8 - 1.73865836002188d0*x270) + 5.0d0*rho_a**(-2.0d0)*x
     &101*x73 - 2.75160604074552d0*x170*x271 - 5.15926132639785d0*x171*x
     &6*x73 + x20*x261 + x25*(0.36240302520243d0*gamma_aa*x106*x77 - 0.0
     &133180320895522d0*rho_a**(-7.0d0)*x273*x56 - x108*x274 - 1.4609649
     &5663704d0*x169*x276 + x274*(-0.804d0*x275 + 1.804d0) + x281*(-0.00
     &291413479477612d0*x114 - 0.00144175546132004d0*x116 - 0.0085665978
     &5039798d0*x119 + 0.188341686653414d0*x125*x278 + 2.48756218905473d
     &0*x136 + x140*x279 - 2.48756218905473d0*x162 + 0.179350162766415d0
     &*x166 - 0.502244497742438d0*x277) - x319*( 0.0131136065764925d0*rh
     &o_a**(-10.0d0)*x61 - 0.502244497742438d0* x125*x126*x288 + x127*x2
     &88*x292 - 1.24378109452736d0*x133*x56*( 0.0484078840501527d0*x282 
     &+ 0.0084375d0*x285) + x134**2*x293*x294 + x137*x279*x317 - 2.28309
     &273345556d0*x140*x287 - x212*(x148*x302 + 0.0889588682240029d0*x28
     &7 + x303*x305 + 0.180246913580247d0* x307 + 0.180246913580247d0*x3
     &12) + 0.46757528851585d0*x283*x295 + 0.00456555896084678d0*x283 + 
     &0.012849896775597d0*x285*x56 + 0.920781579194469d0*x287*x67 + 2.48
     &756218905473d0*x296*x297 + x299*(-0.00970460380625487d0*x120 + x14
     &7*x298 + 0.0720987654320988d0*x160) + x316*(x164*x302 + 0.98707785
     &2896471d0*x287 + 4.0d0*x3*x305 + 2.0d0*x307 + 2.0d0* x312))) - 1.4
     &6096495663704d0*x260*x75 + x263*x266*x56*(-x262 - 1.73865836002188
     &d0*x79 + 0.478410284622892d0*x8*x93*x94 + 0.458812622783553d0*x83)
     &) - 7.5d0*x102*x73 - 4.38289486991112d0* x11*x260 - 2.0d0*x170 - 1
     &.81712059283214d0*x76
      d2Ea(2) = rho_a*(1.0484109910932d0*gamma_aa*x172*x266 - x100 - 0.1
     &74735165182199d0 *x104*x105*x8*x97 + 1.37580302037276d0*x25*x271*(
     & 0.0453003781503037d0*x174 + x319*x321) - x25*( -0.004994262033582
     &09d0*gamma_aa*x267*x273 + 0.105700882350709d0* x107 - 0.7304824783
     &1852d0*x233*x276 + x281*x328 + x319*( -0.00437120219216418d0*x113*
     &x216 - x118*x324 - x121*x81 + x124* x331 + x133*x220*(0.0152867002
     &26364d0*x115 + 0.005625d0*x118) - x134*x294*x333 + x141*x81 - x212
     &*(-x142 + x148*x349 + x148*x350 + x351*x352 + 0.180246913580247d0*
     &x356 + 0.180246913580247d0*x359) - x221*x296 - x226*x317 - 0.25112
     &2248871219d0*x278*x330*x56 - 0.175340733193444d0*x295*x329 + x299*
     &x341 + x316*(-x163 + x164* x349 + x164*x350 + x351*x360 + 2.0d0*x3
     &56 + 2.0d0*x359) - 0.00144175546132004d0*x329 + x331*(-x291 + x332
     & - 0.099435d0) + x339*x340 + x346*x347))) - x173 + x238*(x175 + x3
     &21*x322)
      d2Ea(3) = rho_a*(-x238*(-1.09572371747778d0*x11*x232*x280*x328 + 0
     &.000936424131296642d0*x273*x361 + x322*(0.00109280054804105d0* gam
     &ma_aa*x58*x59 + x127*x287*(-0.007768359375d0*x122 + x291 + 0.04971
     &75d0) - x158*x212*(x303*x364 + 0.180246913580247d0*x363 + 0.180246
     &913580247d0*x367) + 2.48756218905473d0*x158*x315*(x164* x364 + x36
     &3 + x367) + x213 + x215 + 0.188341686653414d0*x287*x342 + 0.065752
     &7749475415d0*x295*x46 - x326*x339 - x327*x346 + 1.24378109452736d0
     &*x341*x344 - 1.24378109452736d0*x343 + x55/x56 )) + 0.250785300325
     &332d0*x265*x272*x361*x8**2)
      d2Ea(4) = -x258*(rho_a*x111*(x103*x378 - x22*(-x109*x378 + 5.45136
     &177849642d0*x11* x167*x29*x378*x71 - x112*(-0.251122248871219d0*x1
     &25*x380 + x131* x340*x377 - 0.00699626865671642d0*x133*x379 + 0.00
     &349813432835821d0*x134*x293*x382*x45 - x212*x40*(x246*x387 - x352*
     &x385 + 0.180246913580247d0*x386 - 0.180246913580247d0*x389 - 0.136
     &470991025459d0*x4) + 0.00349813432835821d0*x296*x383 + x299* x40*(
     &-0.0720987654320988d0*x374 + x384) + x316*x40*(x239*x387 - x360*x3
     &85 + 2.0d0*x386 - 2.0d0*x389 - 1.51426716069345d0*x4) + x347*x378 
     &+ 0.00856659785039798d0*x379 + 0.0941708433267071d0* x380*x381*x60
     &))) + x23*x236*x373)
      d2Ea(5) = x0*x238*x318*(-x124*x391 - 0.00383298875050529d0*x131*x2
     &11*x383 - x207* x392*x394 + 0.00383298875050529d0*x218*x333*x382 +
     & 0.255849428582226d0*x225*x257 - 0.00766597750101059d0*x252*x390 -
     & 2.72568088924821d0*x257*x321*x400 + x325*x392*x393 + x381*x391 + 
     &0.00938662444277524d0*x390 - x398*x399*(x164*x395 + x184*x396 - x1
     &84*x397 - x199*x239) + x398*x53*(x148*x395 - x199*x246 + x209* x39
     &6 - x209*x397))
      d2Ea(6) = x22*x259*(-0.0153319550020212d0*x182*x383*x393 - 0.20637
     &0453055914d0* x228*x290*x292 - 2.72568088924821d0*x240*x371*x394 +
     & 0.0229979325030318d0*x252*x401 + 6.78030071952291d0*x373**2*x400 
     &+ x399*x404*(-x375*x388 + 0.815290678739414d0*x402 + 1.81712059283
     &214d0*x403) - 0.0281598733283257d0*x401 - x404*x53*( -x376*x388 + 
     &0.073476814256762d0*x402 + 0.163765189230551d0*x403 ) - 2.15605617
     &215923d-5*x218*x284*x293/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_2_1(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446, x447, x448, x449, x450
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = rho_a**(-2.66666666666667d0)
      x2 = 4.60115111447049d0
      x3 = 1d0/x2
      x4 = x1*x3
      x5 = gamma_aa*x4
      x6 = mu**2
      x7 = rho_a**0.666666666666667d0
      x8 = x3*x6/x7
      x9 = 5.45136177849642d0*x8 - 7.0d0
      x10 = 0.239205142311446d0*x9
      x11 = 0.826307487110758d0
      x12 = rho_a**0.333333333333333d0
      x13 = 1d0/x12
      x14 = 21.1705915781930d0
      x15 = 1d0/x14
      x16 = mu**4
      x17 = 1d0/x0
      x18 = 59.4346904801033d0*mu*x11*x13 + 4.95289087334194d0*x15*x16*x
     &17 - 49.0622560064678d0*x8 - 81.0d0
      x19 = 1d0/x18
      x20 = x19*exp(-1.4385538026587772d0*x8)
      x21 = x10*x20*x5 + 1.0d0
      x22 = 1.804d0 - 0.804d0/x21
      x23 = 7.26848237132856d0*x22
      x24 = 0.682784063255296d0
      x25 = 2.14502939711103d0
      x26 = 1d0/x25
      x27 = x13*x26
      x28 = 0.550321208149104d0*x27
      x29 = 1d0/rho_a
      x30 = 9.86960440108936d0
      x31 = 1d0/x30
      x32 = 0.0833333333333333d0*x31
      x33 = x29*x32
      x34 = x33*x6
      x35 = x28 - x34
      x36 = exp(-3.3019272488946267d0*x2*x7/x6)
      x37 = mu*x36
      x38 = 1d0/mu
      x39 = x25*x38
      x40 = mu**3
      x41 = mu*x27
      x42 = x33*x40 - 0.825481812223657d0*x41 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x12*x39)
      x43 = x35*x37 + x42
      x44 = 0.275160604074552d0*x41
      x45 = x43*x44 - 0.375d0
      x46 = x24*x45
      x47 = 0.682784063255296d0
      x48 = 2.72568088924821d0*x47
      x49 = exp(-4.1274090611182839d0*x41)
      x50 = 0.0206714418186455d0*x5 + 1.0d0
      x51 = 0.804d0/x50
      x52 = 0.0938662444277524d0*x5 + 1.0d0
      x53 = x52**(-2)
      x54 = 0.050475572023115d0*x5
      x55 = 1d0/gamma_aa
      x56 = tau_a*x55
      x57 = 8.0d0*rho_a*x56
      x58 = x57 - 1.0d0
      x59 = x5*x58
      x60 = 0.0567850185260044d0*x59 - 0.45d0
      x61 = 0.126188930057787d0*x59 - 1.0d0
      x62 = x58*x61
      x63 = x54*x62 + 1.0d0
      x64 = x63**(-0.5d0)
      x65 = x54 + x60*x64
      x66 = x65**2
      x67 = 0.00909806606836394d0*x4
      x68 = gamma_aa*x67
      x69 = x64*(0.0102353243269094d0*x59 - 0.0811111111111111d0) + x68
      x70 = rho_a**(-5.33333333333333d0)*x15
      x71 = rho_a**(-2)
      x72 = tau_a**(-2)
      x73 = 0.0028125d0*x72
      x74 = 0.00286625629244325d0*x70 + x71*x73
      x75 = sqrt(x74)
      x76 = abs(gamma_aa)
      x77 = x75*x76
      x78 = x69*x77
      x79 = 1.24378109452736d0*x78
      x80 = gamma_aa**2
      x81 = x70*x80
      x82 = rho_a**(-8.0d0)
      x83 = 0.0102659822546843d0
      x84 = gamma_aa**3
      x85 = x83*x84
      x86 = x82*x85
      x87 = x71*x72
      x88 = x80*x87
      x89 = 0.015625d0*x88 + 1.0d0
      x90 = x89**(-2)
      x91 = 0.02485875d0*x88*x90 + 0.123456790123457d0
      x92 = x5*x91
      x93 = 0.000135164574498753d0*x81 + 0.000182133424673508d0*x86 + 0.
     &0021416494625995d0*x88 + 0.0941708433267071d0*x92
      x94 = x53*(0.0896750813832074d0*x66 - x79 + x93) + 1.0d0
      x95 = 1d0/x94
      x96 = -0.804d0*x95
      x97 = x49*(x51 + x96)
      x98 = x12*x23
      x99 = -x28 + x34
      x100 = -x37*x99 + x42
      x101 = x100*x44 - 0.375d0
      x102 = x101*x24
      x103 = -x57 + 1.0d0
      x104 = x103*x5
      x105 = 0.0567850185260044d0*x104 + 0.45d0
      x106 = 0.126188930057787d0*x104 + 1.0d0
      x107 = x103*x106
      x108 = x107*x54 + 1.0d0
      x109 = x108**(-0.5d0)
      x110 = -x105*x109 + x54
      x111 = x110**2
      x112 = -x109*(0.0102353243269094d0*x104 + 0.0811111111111111d0) + 
     &x68
      x113 = x112*x77
      x114 = x53*(0.0896750813832074d0*x111 - 1.24378109452736d0*x113 + 
     &x93) + 1.0d0
      x115 = 1d0/x114
      x116 = -0.804d0*x115
      x117 = x12*x48
      x118 = rho_a**(-0.666666666666667d0)
      x119 = x22*x46
      x120 = x118*x119
      x121 = x2*x36*x38
      x122 = rho_a**(-0.333333333333333d0)*x121
      x123 = 2.20128483259642d0*x122
      x124 = 0.183440402716368d0*x17
      x125 = rho_a**(-2.0d0)*x32
      x126 = -x124*x26 + x125*x6
      x127 = mu*x17*x26
      x128 = x36*x39
      x129 = 1.21141372855476d0*x118*x128 - x125*x40 + x126*x37 + 0.2751
     &60604074552d0 *x127
      x130 = 0.275160604074552d0*x13
      x131 = 0.0917202013581841d0*x100*x17 - x130*(x123*x99 + x129)
      x132 = 0.318309886183791d0*mu
      x133 = x132*x98
      x134 = x132*x29
      x135 = 3.75d0*x134
      x136 = rho_a**(-4.33333333333333d0)
      x137 = x136*x15
      x138 = x137*x6
      x139 = rho_a**(-3.66666666666667d0)
      x140 = x139*x3
      x141 = x140*x9
      x142 = x138*x9
      x143 = x11*x17
      x144 = rho_a**(-1.66666666666667d0)
      x145 = x144*x3
      x146 = mu*x145
      x147 = rho_a**(-2.33333333333333d0)
      x148 = x15*x40
      x149 = x147*x148
      x150 = 19.8115634933678d0*x143 - 32.7081706709785d0*x146 + 6.60385
     &449778925d0* x149
      x151 = x150*x19
      x152 = mu*x4
      x153 = x10*x152
      x154 = -0.869329180010942d0*x138 - 0.637880379497189d0*x141 + 0.22
     &9406311391776d0*x142 + x151*x153
      x155 = x20/x21**2
      x156 = x154*x155
      x157 = 5.84385982654816d0*x12
      x158 = x156*x157
      x159 = gamma_aa*x46
      x160 = rho_a**(-3.33333333333333d0)
      x161 = 0.148394183600699d0
      x162 = x161/x50**2
      x163 = x160*x162
      x164 = -x51 + 1.804d0
      x165 = 0.90856029641607d0*x118
      x166 = x165*x47
      x167 = x94**(-2)
      x168 = 2.19144743495556d0*x12
      x169 = rho_a**(-9.0d0)
      x170 = x169*x85
      x171 = rho_a**(-6.33333333333333d0)*x15
      x172 = x171*x80
      x173 = rho_a**(-3)
      x174 = x173*x72
      x175 = x174*x80
      x176 = gamma_aa*x140
      x177 = 0.251122248871219d0*x91
      x178 = x88/x89
      x179 = 0.001553671875d0*x178
      x180 = x179 - 0.0497175d0
      x181 = x180*x84
      x182 = x72*x90
      x183 = 0.0941708433267071d0*x182
      x184 = rho_a**(-5.66666666666667d0)*x3
      x185 = x183*x184
      x186 = 1d0/x75
      x187 = 1d0/x76
      x188 = x187*x69
      x189 = x186*x188
      x190 = 0.00764335011318201d0*x171 + x173*x73
      x191 = x190*x80
      x192 = x189*x191
      x193 = 1d0/x52
      x194 = 0.000108672317896998d0*x81 + 0.0001464352734375d0*x86 + 0.0
     &0172188616792999d0*x88 + 0.0757133580346725d0*x92
      x195 = x194 + 0.0720987654320988d0*x66 - x78
      x196 = x193*x195
      x197 = 0.622661654578789d0*x196
      x198 = 0.0242615095156372d0*x140
      x199 = tau_a*x1
      x200 = x139*x58
      x201 = 0.151426716069345d0*x200
      x202 = -gamma_aa*x201 + 0.454280148208035d0*x199
      x203 = x202*x64
      x204 = 0.180246913580247d0*x3
      x205 = x3*x61
      x206 = 0.20190228809246d0*x199
      x207 = 0.0673007626974867d0*x62
      x208 = 0.336503813487433d0*x200
      x209 = x15*(-gamma_aa*x208 + 1.0095114404623d0*x199)
      x210 = x1*x58
      x211 = 0.0252377860115575d0*x210
      x212 = x209*x211
      x213 = -gamma_aa*x212 + x176*x207 - x205*x206
      x214 = x63**(-1.5d0)
      x215 = x214*x60
      x216 = x213*x215
      x217 = -gamma_aa*x198 + x203*x204 + 0.180246913580247d0*x216
      x218 = x217*x77
      x219 = 0.269203050789947d0*x140
      x220 = 2.0d0*x3
      x221 = -gamma_aa*x219 + x203*x220 + 2.0d0*x216
      x222 = x221*x65
      x223 = -0.00145706739738806d0*x170 - 0.000720877730660018d0*x172 -
     & 0.00428329892519899d0*x175 - x176*x177 + x176*x197 + x181*x185 + 
     &1.24378109452736d0*x192 - 1.24378109452736d0*x218 + 0.089675081383
     &2074d0*x222
      x224 = x47*x53
      x225 = x223*x224
      x226 = x49*(0.12080100840081d0*gamma_aa*x163 - x164*x166 + x166*(x
     &96 + 1.804d0 ) + x167*x168*x225)
      x227 = x147*x161
      x228 = x227*x9
      x229 = 1.39788132145759d0*x155
      x230 = x228*x229
      x231 = x147*x162
      x232 = 0.0226501890751519d0*x231
      x233 = 2.0d0*x105
      x234 = x108**(-1.5d0)
      x235 = 0.20190228809246d0*x145
      x236 = x235*x61
      x237 = x236*x56
      x238 = x4*x61
      x239 = x238*x58
      x240 = 0.0252377860115575d0*x239
      x241 = x144*x56
      x242 = 1.0095114404623d0*x241
      x243 = 0.126188930057787d0*x1
      x244 = -x242 + x243*x58
      x245 = gamma_aa*x15
      x246 = x211*x245
      x247 = x244*x246
      x248 = -x237 + x240 + x247
      x249 = x234*x248
      x250 = 0.10095114404623d0*x4
      x251 = 0.454280148208035d0*x241
      x252 = 0.0567850185260044d0*x1
      x253 = x103*x252 + x251
      x254 = x109*x253
      x255 = -x220*x254 + x250
      x256 = x233*x249 + x255
      x257 = x110*x256
      x258 = 0.180246913580247d0*x105
      x259 = -x204*x254 + x67
      x260 = x249*x258 + x259
      x261 = 1.24378109452736d0*x77
      x262 = 0.000270329148997507d0*x70
      x263 = x4*x91
      x264 = 0.00428329892519899d0*x87
      x265 = x80*x83
      x266 = x265*x82
      x267 = x112*x187
      x268 = 1.24378109452736d0*gamma_aa
      x269 = x268*x75
      x270 = 0.0720987654320988d0*x111 - x113 + x194
      x271 = x193*x4
      x272 = 0.233498120467046d0*x271
      x273 = rho_a**(-4.66666666666667d0)
      x274 = x273*x3
      x275 = x274*x80
      x276 = x180*x275
      x277 = -gamma_aa*x262 - gamma_aa*x264 + x183*x276 - 0.094170843326
     &7071d0*x263 - 0.000546400274020523d0*x266 + x267*x269 + x270*x272
      x278 = -0.0896750813832074d0*x257 + x260*x261 + x277
      x279 = x114**(-2)
      x280 = x224*x279
      x281 = x278*x280
      x282 = 1.09572371747778d0*x12
      x283 = 2.0d0*x49
      x284 = 0.90856029641607d0*x145
      x285 = 1d0/x108
      x286 = 0.0254778337106067d0*x136
      x287 = x103*x245
      x288 = x106*x235
      x289 = x285*(x286*x287 + x288)
      x290 = 0.0720987654320988d0*x109*x110
      x291 = 0.0818825946152754d0*x145
      x292 = x109*x77
      x293 = tau_a**(-3)
      x294 = x293*x80
      x295 = x294*x71
      x296 = 0.0028125d0*x295
      x297 = x186*x267
      x298 = 0.0757133580346725d0*x181
      x299 = x293*x90
      x300 = x274*x298*x299 - 0.00344377233585999d0*x295
      x301 = x296*x297 + x300
      x302 = x290*(-x233*x289 + x284) - x292*(-x258*x289 + x291) + x301
      x303 = x49*x53
      x304 = x0*x279*x303
      x305 = x131*x132
      x306 = x22*x305
      x307 = x156*x159
      x308 = 11.6877196530963d0*x12
      x309 = rho_a**(-1.66666666666667d0)
      x310 = 1.27576075899438d0*x141
      x311 = x18**(-2)
      x312 = x311*exp(-2.8771076053175544d0*x8)/x21**3
      x313 = x154*x312
      x314 = x147*x26
      x315 = 0.166666666666667d0*rho_a**(-3.0d0)*x31
      x316 = x36/x40
      x317 = rho_a**(-6.0d0)
      x318 = x16*x317*x83
      x319 = x6*x70
      x320 = x137*x151*x40
      x321 = 0.0322515344331995d0
      x322 = x321/x50**3
      x323 = 0.60570686427738d0*x309*x47
      x324 = x176*x91
      x325 = x182*x184
      x326 = 1.24532330915758d0*x176
      x327 = x47/(x114**3*x52**4)
      x328 = x168*x223*x327
      x329 = rho_a**(-7.33333333333333d0)*x15
      x330 = x329*x80
      x331 = rho_a**(-4)
      x332 = x331*x72
      x333 = gamma_aa*x273
      x334 = x3*x333
      x335 = rho_a**(-6.66666666666667d0)*x3
      x336 = tau_a**(-4)
      x337 = x336*x90
      x338 = 0.00014565673828125d0*gamma_aa**4*x331*x337
      x339 = x84*(-0.010875703125d0*x178 + x338 + 0.1491525d0)
      x340 = x74**(-1.5d0)*x76
      x341 = 1.24378109452736d0*x69
      x342 = x195*x53
      x343 = x187*x217
      x344 = x186*x191
      x345 = 0.0720987654320988d0*x3
      x346 = 1.24378109452736d0*x221
      x347 = tau_a*x139
      x348 = x333*x58
      x349 = x64*(-2.42282745710952d0*x347 + 0.555231292254265d0*x348)
      x350 = 0.360493827160494d0*x3
      x351 = x213*x214
      x352 = x202*x351
      x353 = gamma_aa*x209
      x354 = x215*(-0.40380457618492d0*x199*x209 + 0.134601525394973d0*x
     &200*x353 + 1.07681220315979d0*x205*x347 - 0.246769463224118d0*x205
     &*x348 - x246*(-5.38406101579893d0*x347 + 1.23384731612059d0*x348))
     &
      x355 = 0.60570686427738d0*x61
      x356 = 0.0757133580346725d0*x210
      x357 = x213*x60
      x358 = x357*x63**(-2.5d0)
      x359 = x358*(0.20190228809246d0*x176*x62 - x199*x3*x355 - x353*x35
     &6)
      x360 = 0.00363922642734557d0*x4
      x361 = gamma_aa*x360
      x362 = x361 + x64*(0.00409412973076377d0*x59 - 0.0324444444444444d
     &0)
      x363 = 1.24378109452736d0*x362
      x364 = -0.0011714821875d0*x170 - 0.000579585695450655d0*x172 - 0.0
     &0344377233585999d0*x175 + x192 - x218 + 0.0720987654320988d0* x222
     & + x298*x325 - 0.20190228809246d0*x324
      x365 = x168*x280
      x366 = x234*(0.0252377860115575d0*x1*x287*(x103*x243 + x242) + 0.0
     &252377860115575d0*x107*x4 + x288*x56)
      x367 = -0.0896750813832074d0*x110*(x233*x366 + x255) + x261*(x258*
     &x366 + x259) + x277
      x368 = x280*x282
      x369 = gamma_aa*x70
      x370 = 0.00856659785039798d0*gamma_aa
      x371 = gamma_aa*x75
      x372 = 2.48756218905473d0*x371
      x373 = 0.466996240934092d0*x271
      x374 = 0.188341686653414d0*x182*x276 - 0.179350162766415d0*x257 + 
     &2.48756218905473d0*x260*x77 - 0.188341686653414d0*x263 - 0.0010928
     &0054804105d0*x266 + x267*x372 + x270*x373 - 0.000540658297995013d0
     &*x369 - x370*x87
      x375 = gamma_aa*x171
      x376 = -x179 + 0.0497175d0
      x377 = x185*x80
      x378 = 0.00932203125d0*x178
      x379 = x186*x55*x76
      x380 = -x251 + x252*x58
      x381 = x380*x64
      x382 = x237 - x240 - x247
      x383 = x215*x382
      x384 = x204*x381 + 0.180246913580247d0*x383 + x67
      x385 = x187*x384
      x386 = 1.24378109452736d0*x344
      x387 = x345*x381 + x360 + 0.0720987654320988d0*x383
      x388 = x182*x376
      x389 = x188*x75
      x390 = x220*x381 + x250 + 2.0d0*x383
      x391 = 0.0720987654320988d0*x65
      x392 = -gamma_aa*x389 + 0.00344377233585999d0*gamma_aa*x87 + 0.075
     &7133580346725d0*x263 + 0.0004393058203125d0*x266 + 0.0757133580346
     &725d0*x275*x388 + 0.000217344635793995d0*x369 - x384*x77 + x390*x3
     &91
      x393 = 0.622661654578789d0*x176*x193
      x394 = x199*x55
      x395 = x64*(-x201 + 1.21141372855476d0*x394)
      x396 = x351*x380
      x397 = x214*x382
      x398 = x202*x204
      x399 = x145*x355
      x400 = x244*x245
      x401 = -0.0757133580346725d0*x239 - x356*x400 + x399*x56
      x402 = x358*x401
      x403 = 0.20190228809246d0*x144
      x404 = x209*x403
      x405 = x215*(x140*x207 - x15*x206*x244 + 0.0673007626974867d0*x200
     &*x400 - 0.538406101579893d0*x205*x394 - x212 - x246*(-x208 + 2.692
     &03050789947d0*x394) + x404*x56)
      x406 = x202*x220
      x407 = rho_a**(-5.0d0)
      x408 = x15*x60
      x409 = x244*x408*(-0.050475572023115d0*x210 + 0.40380457618492d0*x
     &241)
      x410 = x380*x382
      x411 = 1d0/x63
      x412 = x411*x60
      x413 = x382*x401*x412
      x414 = x286*x58
      x415 = x236 + x245*x414
      x416 = x285*x415
      x417 = x233*x416 + x284
      x418 = x258*x416 + x291
      x419 = x290*x417 - x292*x418 + x301
      x420 = x412*x415
      x421 = 2.0d0*x420
      x422 = 0.180246913580247d0*x420
      x423 = x64*(x291 - x422)
      x424 = x189*x296 + x300 + x391*x64*(x284 - x421) - x423*x77
      x425 = x173*x294
      x426 = x184*x299
      x427 = x338 - x378 + 0.099435d0
      x428 = x293*x71
      x429 = x186*x295
      x430 = 0.0327530378461102d0*x145
      x431 = x411*x415
      x432 = x412*(-0.203822669684853d0*tau_a*x137 + 0.336503813487433d0
     &*x238 + 0.110403946079296d0*x369*x58 - x404)
      x433 = x213*x411
      x434 = 0.0764335011318201d0*gamma_aa*x137*x58 + x399
      x435 = x357*x434/x63**2
      x436 = gamma_aa*x428
      x437 = 0.103185226527957d0*x274*x294*x90
      x438 = 1.3628404446241d0*x109
      x439 = x187*x418
      x440 = 0.0720987654320988d0*x105*x416 + x430
      x441 = x253*x415
      x442 = x15*(-0.203822669684853d0*x160*x56 + x244*x403 + x414)
      x443 = x248*x285*x434
      x444 = 1.3628404446241d0*x234
      x445 = -x109*(0.00409412973076377d0*x104 + 0.0324444444444444d0) +
     & x361
      x446 = x115*x53
      x447 = x336*x71*x80
      x448 = x160*x408
      x449 = x145*x415
      x450 = 2.72568088924821d0*x214
      Ea = -x0*(-x23*x46 + x48*x97)
      d1Ea(1) = -rho_a*(-2.42282745710952d0*x120 + x131*x133 + x135*x97 
     &- x158*x159 + x226) + x102*x98 - x117*x49*(x116 + x51)
      d1Ea(2) = rho_a*(x230*x45 + x283*(x232 + x281*x282))
      d1Ea(3) = -x302*x304*x48
      d2Ea(1) = rho_a*(-gamma_aa*x156*x305*x308 + 5.0d0*rho_a**(-2.0d0)*
     &x132*x97 - 4.84565491421904d0*x118*x306 + 3.89590655103211d0*x118*
     &x307 - 1.61521830473968d0*x119*x309 - 2.75160604074552d0*x127*x226
     & + x133*(-x124*(-x123*x35 + x129) + x130*(-0.366880805432736d0*mu*
     & x314 + 0.733761610865473d0*rho_a**(-1.33333333333333d0)*x121*x35 
     &+ 4.84565491421904d0*x118*x14*x316*x35 - 4.40256966519284d0*x122* 
     &x126 - 0.80760915236984d0*x128*x309 - 2.66666666666667d0*x29*x30* 
     &x316 + x315*x40 + x37*(0.244587203621824d0*x314 - x315*x6)) + 0.12
     &2293601810912d0*x147*x43) + x155*x157*x159*(-mu*x151*x310 + x10*x1
     &50*x311*x4*x6*(39.6231269867355d0*x143 - 65.416341341957d0* x146 +
     & 13.2077089955785d0*x149) + x153*x19*(-26.415417991157d0* x11*x147
     & - 15.4089938281749d0*x148*x160 + 54.5136177849642d0*x152 ) + 2.33
     &889472482303d0*x274*x9 + 0.220008881071042d0*x318*x9 - 1.667435730
     &22264d0*x318 - 1.60584417974243d0*x319*x9 + 6.08530426007659d0*x31
     &9 + 0.458812622783553d0*x320*x9 - 1.73865836002188d0*x320) - x157*
     &x313*x46*x80*(-1.73865836002188d0 *x138 + 0.458812622783553d0*x142
     & + 0.478410284622892d0*x151*x152* x9 - x310) - 5.15926132639785d0*
     &x227*x6*x97 + x49*( 0.36240302520243d0*gamma_aa*x136*x162 - 0.0133
     &180320895522d0* rho_a**(-7.0d0)*x322*x80 - 1.46096495663704d0*x118
     &*x225*x279 - x164*x323 + x323*(x116 + 1.804d0) + x328*(-0.00291413
     &479477612d0* x170 - 0.00144175546132004d0*x172 - 0.008566597850397
     &98d0*x175 + 0.188341686653414d0*x181*x325 + 2.48756218905473d0*x19
     &2 + x196* x326 - 2.48756218905473d0*x218 + 0.179350162766415d0*x22
     &2 - 0.502244497742438d0*x324) - x365*(0.0131136065764925d0*rho_a**
     &( -10.0d0)*x85 - 0.502244497742438d0*x181*x182*x335 + x183*x335* x
     &339 - 1.24378109452736d0*x189*x80*(0.0484078840501527d0*x329 + 0.0
     &084375d0*x332) + x190**2*x340*x341 + x193*x326*x364 - 2.2830927334
     &5556d0*x196*x334 - x261*(x204*x349 + 0.0889588682240029d0*x334 + x
     &350*x352 + 0.180246913580247d0*x354 + 0.180246913580247d0*x359) + 
     &0.46757528851585d0*x330*x342 + 0.00456555896084678d0*x330 + 0.0128
     &49896775597d0*x332*x80 + 0.920781579194469d0*x334*x91 + 2.48756218
     &905473d0*x343*x344 + x346*(-0.00970460380625487d0*x176 + x203*x345
     & + 0.0720987654320988d0*x216) + x363*(x220*x349 + 4.0d0*x3*x352 + 
     &0.987077852896471d0*x334 + 2.0d0*x354 + 2.0d0*x359)))) - 14.536964
     &7426571d0*x12*x306 + 4.84565491421904d0*x120 - 7.5d0* x134*x97 - 2
     &.0d0*x226 + x307*x308
      d2Ea(2) = rho_a*(-2.79576264291519d0*gamma_aa*x101*x228*x313 - 0.0
     &691804894611514d0*mu*x131*x147*x229*x9 + 0.465960440485865d0* x101
     &*x155*x160*x161*x9 + x102*x158 + 1.37580302037276d0*x127*x49* (0.0
     &453003781503037d0*x231 + x365*x367) - x49*( -0.00499426203358209d0
     &*gamma_aa*x317*x322 - 0.73048247831852d0* x118*x281 + 0.1057008823
     &50709d0*x163 + x328*x374 + x365*(-x140* x177 + x140*x197 - 0.00437
     &120219216418d0*x169*x265 - x174*x370 + x180*x377 + x189*x268*(0.01
     &5286700226364d0*x171 + 0.005625d0*x174 ) - x190*x341*x379 - x261*(
     &-x198 + x204*x395 + x204*x396 + x397* x398 + 0.180246913580247d0*x
     &402 + 0.180246913580247d0*x405) - x269*x343 - x272*x364 - 0.251122
     &248871219d0*x325*x376*x80 - 0.175340733193444d0*x342*x375 + x346*x
     &387 + x363*(-x219 + x220* x395 + x220*x396 + x397*x406 + 2.0d0*x40
     &2 + 2.0d0*x405) - 0.00144175546132004d0*x375 + x377*(-x338 + x378 
     &- 0.099435d0) + x385*x386 + x392*x393))) + x101*x230 + x283*(x232 
     &+ x367*x368)
      d2Ea(3) = -rho_a*(-x283*(x278*x282*x327*x374 - 0.00093642413129664
     &2d0*x322*x407 - x368*(0.00109280054804105d0*gamma_aa*x82*x83 + x18
     &3*x334*( -0.007768359375d0*x178 + x338 + 0.0497175d0) - x214*x261*
     &(x350* x410 + 0.180246913580247d0*x409 + 0.180246913580247d0*x413)
     & + 2.48756218905473d0*x214*x362*(x220*x410 + x409 + x413) + x262 +
     & x264 + 0.188341686653414d0*x334*x388 + 0.0657527749475415d0*x342*
     & x70 - x372*x385 - x373*x392 + 1.24378109452736d0*x387*x390 - 1.24
     &378109452736d0*x389 + x79/x80)) + 0.668760800867552d0*x312* x321*x
     &407*x45*x9**2)
      d2Ea(4) = -x303*(rho_a*x167*(x135*x424 - x47*(5.45136177849642d0*x
     &12*x223*x424*x53 *x95 - x165*x424 - x168*(-0.251122248871219d0*x18
     &1*x426 + x187* x386*x423 - 0.00699626865671642d0*x189*x425 + 0.003
     &49813432835821d0*x190*x340*x428*x69 - x261*x64*(x291*x433 - x398*x
     &431 - 0.136470991025459d0*x4 + 0.180246913580247d0*x432 - 0.180246
     &913580247d0*x435) + 0.00349813432835821d0*x343*x429 + x346*x64*(-0
     &.0720987654320988d0*x420 + x430) + x363*x64*(x284* x433 - 1.514267
     &16069345d0*x4 - x406*x431 + 2.0d0*x432 - 2.0d0* x435) + x393*x424 
     &+ 0.00856659785039798d0*x425 + 0.0941708433267071d0*x426*x427*x84)
     &)) + x117*x279*x419)
      d2Ea(5) = x0*x280*x283*(0.00383298875050529d0*x112*x379*x428 - x18
     &0*x437 - 0.00383298875050529d0*x187*x260*x429 - x256*x438*x440 + 0
     &.255849428582226d0*x271*x302 - 0.00766597750101059d0*x297*x436 - 2
     &.72568088924821d0*x302*x367*x446 + x371*x438*x439 + x427*x437 + 0.
     &00938662444277524d0*x436 - x444*x445*(x220*x441 + x233*x442 - x233
     &*x443 - x248*x284) + x444*x77*(x204*x441 - x248*x291 + x258* x442 
     &- x258*x443))
      d2Ea(6) = x304*x47*(-0.0153319550020212d0*x109*x429*x439 - 0.20637
     &0453055914d0* x274*x337*x339 - 2.72568088924821d0*x285*x417*x440 +
     & 0.0229979325030318d0*x297*x447 + 6.78030071952291d0*x419**2*x446 
     &+ x445*x450*(-x421*x434 + 0.815290678739414d0*x448 + 1.81712059283
     &214d0*x449) - 0.0281598733283257d0*x447 - x450*x77*( -x422*x434 + 
     &0.073476814256762d0*x448 + 0.163765189230551d0*x449 ) - 2.15605617
     &215923d-5*x112*x331*x340/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_2_2(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446, x447, x448, x449, x450, x451, x452, x453, x454, x455, x456
     &, x457, x458, x459, x460
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.682784063255296d0
      x2 = 7.26848237132856d0*x1
      x3 = rho_a**(-2.66666666666667d0)
      x4 = 4.60115111447049d0
      x5 = 1d0/x4
      x6 = x3*x5
      x7 = gamma_aa*x6
      x8 = mu**2
      x9 = rho_a**0.666666666666667d0
      x10 = x5*x8/x9
      x11 = exp(-1.4385538026587772d0*x10)
      x12 = 1d0/x8
      x13 = x4*x9
      x14 = 0.183440402716368d0*x13
      x15 = 21.1705915781930d0
      x16 = x12*x15
      x17 = 0.010095114404623d0*x0*x16
      x18 = rho_a**2.0d0
      x19 = 97.4090910340024d0
      x20 = mu**(-4)
      x21 = x19*x20
      x22 = 0.147857142857143d0*x18*x21
      x23 = -x14 + x17 + x22
      x24 = x12*x23
      x25 = x11*x24
      x26 = 0.239205142311446d0*x25*x7
      x27 = -x26 + 1.0d0
      x28 = 1.804d0 - 0.804d0/x27
      x29 = 2.14502939711103d0
      x30 = 1d0/x29
      x31 = rho_a**0.333333333333333d0
      x32 = 1d0/x31
      x33 = x30*x32
      x34 = 0.550321208149104d0*x33
      x35 = rho_a**1.0d0
      x36 = 1d0/x35
      x37 = 9.86960440108936d0
      x38 = 1d0/x37
      x39 = 0.0833333333333333d0*x38
      x40 = x36*x39
      x41 = x40*x8
      x42 = x34 - x41
      x43 = exp(-3.3019272488946267d0*x12*x13)
      x44 = mu*x43
      x45 = 1d0/mu
      x46 = x29*x45
      x47 = mu**3
      x48 = mu*x33
      x49 = x40*x47 - 0.825481812223657d0*x48 + 1.77245385090552d0*erf( 
     &1.8171205928321397d0*x31*x46)
      x50 = x42*x44 + x49
      x51 = 0.275160604074552d0*x48
      x52 = x50*x51 - 0.375d0
      x53 = x28*x52
      x54 = 0.682784063255296d0
      x55 = 2.72568088924821d0*x54
      x56 = exp(-4.1274090611182839d0*x48)
      x57 = 0.0206714418186455d0*x7 + 1.0d0
      x58 = 0.804d0/x57
      x59 = 0.0938662444277524d0*x7 + 1.0d0
      x60 = x59**(-2)
      x61 = 0.050475572023115d0*x7
      x62 = 1d0/gamma_aa
      x63 = tau_a*x62
      x64 = 8.0d0*rho_a*x63
      x65 = x64 - 1.0d0
      x66 = x65*x7
      x67 = 0.0567850185260044d0*x66 - 0.45d0
      x68 = 0.126188930057787d0*x66 - 1.0d0
      x69 = x65*x68
      x70 = x61*x69 + 1.0d0
      x71 = x70**(-0.5d0)
      x72 = x61 + x67*x71
      x73 = x72**2
      x74 = 0.00909806606836394d0*x6
      x75 = gamma_aa*x74
      x76 = x71*(0.0102353243269094d0*x66 - 0.0811111111111111d0) + x75
      x77 = 1d0/x15
      x78 = rho_a**(-5.33333333333333d0)*x77
      x79 = rho_a**(-2)
      x80 = tau_a**(-2)
      x81 = 0.0028125d0*x80
      x82 = 0.00286625629244325d0*x78 + x79*x81
      x83 = sqrt(x82)
      x84 = abs(gamma_aa)
      x85 = x83*x84
      x86 = x76*x85
      x87 = 1.24378109452736d0*x86
      x88 = gamma_aa**2
      x89 = x78*x88
      x90 = rho_a**(-8.0d0)
      x91 = 1d0/x19
      x92 = gamma_aa**3
      x93 = x91*x92
      x94 = x90*x93
      x95 = x79*x80
      x96 = x88*x95
      x97 = 0.015625d0*x96 + 1.0d0
      x98 = x97**(-2)
      x99 = 0.02485875d0*x96*x98 + 0.123456790123457d0
      x100 = x7*x99
      x101 = 0.0941708433267071d0*x100 + 0.000135164574498753d0*x89 + 0.
     &000182133424673508d0*x94 + 0.0021416494625995d0*x96
      x102 = x60*(x101 + 0.0896750813832074d0*x73 - x87) + 1.0d0
      x103 = 1d0/x102
      x104 = -0.804d0*x103
      x105 = x56*(x104 + x58)
      x106 = x26 - 1.0d0
      x107 = -x34 + x41
      x108 = -x107*x44 + x49
      x109 = x108*x51 - 0.375d0
      x110 = x109*x31
      x111 = -x64 + 1.0d0
      x112 = x111*x7
      x113 = 0.0567850185260044d0*x112 + 0.45d0
      x114 = 0.126188930057787d0*x112 + 1.0d0
      x115 = x111*x114
      x116 = x115*x61 + 1.0d0
      x117 = x116**(-0.5d0)
      x118 = -x113*x117 + x61
      x119 = x118**2
      x120 = -x117*(0.0102353243269094d0*x112 + 0.0811111111111111d0) + 
     &x75
      x121 = x120*x85
      x122 = x60*(x101 + 0.0896750813832074d0*x119 - 1.24378109452736d0*
     &x121) + 1.0d0
      x123 = 1d0/x122
      x124 = -0.804d0*x123
      x125 = x31*x55
      x126 = rho_a**(-0.666666666666667d0)
      x127 = x1*x53
      x128 = x126*x127
      x129 = 1d0/x0
      x130 = rho_a**(-0.333333333333333d0)*x4
      x131 = x43*x45
      x132 = x130*x131
      x133 = 2.20128483259642d0*x132
      x134 = 0.183440402716368d0*x129
      x135 = x39/x18
      x136 = -x134*x30 + x135*x8
      x137 = mu*x129*x30
      x138 = x43*x46
      x139 = 1.21141372855476d0*x126*x138 - x135*x47 + x136*x44 + 0.2751
     &60604074552d0 *x137
      x140 = 0.275160604074552d0*x32
      x141 = 0.0917202013581841d0*x108*x129 - x140*(x107*x133 + x139)
      x142 = 0.318309886183791d0*mu
      x143 = x141*x142*x28
      x144 = 7.26848237132856d0*x31
      x145 = x14 - x17 - x22
      x146 = rho_a**(-4.33333333333333d0)
      x147 = x146*x77
      x148 = x145*x147
      x149 = rho_a**(-3.66666666666667d0)
      x150 = x149*x5
      x151 = x12*x145
      x152 = x150*x151
      x153 = 0.122293601810912d0*x130
      x154 = 0.0134601525394973d0*rho_a**0.333333333333333d0*x16
      x155 = 0.295714285714286d0*x21
      x156 = x155*x35
      x157 = x153 - x154 - x156
      x158 = x12*x157
      x159 = 0.239205142311446d0*x6
      x160 = 0.229406311391776d0*x148 - 0.637880379497189d0*x152 + x158*
     &x159
      x161 = x160*x31
      x162 = x1*x11
      x163 = 5.84385982654816d0*x162
      x164 = x27**(-2)
      x165 = x164*x52
      x166 = gamma_aa*x165
      x167 = x163*x166
      x168 = x142*x36
      x169 = 3.75d0*x168
      x170 = rho_a**(-3.33333333333333d0)
      x171 = 0.148394183600699d0
      x172 = x171/x57**2
      x173 = x170*x172
      x174 = -x58 + 1.804d0
      x175 = 0.90856029641607d0*x126
      x176 = x175*x54
      x177 = x102**(-2)
      x178 = 2.19144743495556d0*x31
      x179 = rho_a**(-9.0d0)
      x180 = x179*x93
      x181 = rho_a**(-6.33333333333333d0)*x77
      x182 = x181*x88
      x183 = rho_a**(-3)
      x184 = x183*x80
      x185 = x184*x88
      x186 = gamma_aa*x150
      x187 = 0.251122248871219d0*x99
      x188 = x96/x97
      x189 = 0.001553671875d0*x188
      x190 = x189 - 0.0497175d0
      x191 = x190*x92
      x192 = x80*x98
      x193 = 0.0941708433267071d0*x192
      x194 = rho_a**(-5.66666666666667d0)*x5
      x195 = x193*x194
      x196 = 1d0/x83
      x197 = 1d0/x84
      x198 = x197*x76
      x199 = x196*x198
      x200 = 0.00764335011318201d0*x181 + x183*x81
      x201 = x200*x88
      x202 = x199*x201
      x203 = 1d0/x59
      x204 = 0.0757133580346725d0*x100 + 0.000108672317896998d0*x89 + 0.
     &0001464352734375d0*x94 + 0.00172188616792999d0*x96
      x205 = x204 + 0.0720987654320988d0*x73 - x86
      x206 = x203*x205
      x207 = 0.622661654578789d0*x206
      x208 = 0.0242615095156372d0*x150
      x209 = tau_a*x3
      x210 = x149*x65
      x211 = 0.151426716069345d0*x210
      x212 = -gamma_aa*x211 + 0.454280148208035d0*x209
      x213 = x212*x71
      x214 = 0.180246913580247d0*x5
      x215 = x5*x68
      x216 = 0.20190228809246d0*x209
      x217 = 0.0673007626974867d0*x69
      x218 = 0.336503813487433d0*x210
      x219 = x77*(-gamma_aa*x218 + 1.0095114404623d0*x209)
      x220 = x3*x65
      x221 = 0.0252377860115575d0*x220
      x222 = x219*x221
      x223 = -gamma_aa*x222 + x186*x217 - x215*x216
      x224 = x70**(-1.5d0)
      x225 = x224*x67
      x226 = x223*x225
      x227 = -gamma_aa*x208 + x213*x214 + 0.180246913580247d0*x226
      x228 = x227*x85
      x229 = 0.269203050789947d0*x150
      x230 = 2.0d0*x5
      x231 = -gamma_aa*x229 + x213*x230 + 2.0d0*x226
      x232 = x231*x72
      x233 = -0.00145706739738806d0*x180 - 0.000720877730660018d0*x182 -
     & 0.00428329892519899d0*x185 - x186*x187 + x186*x207 + x191*x195 + 
     &1.24378109452736d0*x202 - 1.24378109452736d0*x228 + 0.089675081383
     &2074d0*x232
      x234 = x54*x60
      x235 = x233*x234
      x236 = x56*(0.12080100840081d0*gamma_aa*x173 - x174*x176 + x176*(x
     &104 + 1.804d0 ) + x177*x178*x235)
      x237 = rho_a**(-2.33333333333333d0)
      x238 = x171*x237
      x239 = 1.39788132145759d0*x238*x25
      x240 = x172*x237
      x241 = 0.0226501890751519d0*x240
      x242 = 2.0d0*x113
      x243 = x116**(-1.5d0)
      x244 = rho_a**(-1.66666666666667d0)
      x245 = x244*x5
      x246 = 0.20190228809246d0*x245
      x247 = x246*x68
      x248 = x247*x63
      x249 = x6*x68
      x250 = x249*x65
      x251 = 0.0252377860115575d0*x250
      x252 = x244*x63
      x253 = 1.0095114404623d0*x252
      x254 = 0.126188930057787d0*x3
      x255 = -x253 + x254*x65
      x256 = gamma_aa*x77
      x257 = x221*x256
      x258 = x255*x257
      x259 = -x248 + x251 + x258
      x260 = x243*x259
      x261 = 0.10095114404623d0*x6
      x262 = 0.454280148208035d0*x252
      x263 = 0.0567850185260044d0*x3
      x264 = x111*x263 + x262
      x265 = x117*x264
      x266 = -x230*x265 + x261
      x267 = x242*x260 + x266
      x268 = x118*x267
      x269 = 0.180246913580247d0*x113
      x270 = -x214*x265 + x74
      x271 = x260*x269 + x270
      x272 = 1.24378109452736d0*x85
      x273 = 0.000270329148997507d0*x78
      x274 = x6*x99
      x275 = 0.00428329892519899d0*x95
      x276 = x88*x91
      x277 = x276*x90
      x278 = x120*x197
      x279 = 1.24378109452736d0*gamma_aa
      x280 = x279*x83
      x281 = 0.0720987654320988d0*x119 - x121 + x204
      x282 = x203*x6
      x283 = 0.233498120467046d0*x282
      x284 = rho_a**(-4.66666666666667d0)
      x285 = x284*x5
      x286 = x285*x88
      x287 = x190*x286
      x288 = -gamma_aa*x273 - gamma_aa*x275 + x193*x287 - 0.094170843326
     &7071d0*x274 - 0.000546400274020523d0*x277 + x278*x280 + x281*x283
      x289 = -0.0896750813832074d0*x268 + x271*x272 + x288
      x290 = x122**(-2)
      x291 = x234*x290
      x292 = x289*x291
      x293 = 1.09572371747778d0*x31
      x294 = 2.0d0*x56
      x295 = 0.90856029641607d0*x245
      x296 = 1d0/x116
      x297 = 0.0254778337106067d0*x146
      x298 = x111*x256
      x299 = x114*x246
      x300 = x296*(x297*x298 + x299)
      x301 = 0.0720987654320988d0*x117*x118
      x302 = 0.0818825946152754d0*x245
      x303 = x117*x85
      x304 = tau_a**(-3)
      x305 = x304*x88
      x306 = x305*x79
      x307 = 0.0028125d0*x306
      x308 = x196*x278
      x309 = 0.0757133580346725d0*x191
      x310 = x304*x98
      x311 = x285*x309*x310 - 0.00344377233585999d0*x306
      x312 = x307*x308 + x311
      x313 = x301*(-x242*x300 + x295) - x303*(-x269*x300 + x302) + x312
      x314 = x56*x60
      x315 = x0*x290*x314
      x316 = x162*x166
      x317 = 11.6877196530963d0*x161
      x318 = rho_a**(-1.66666666666667d0)
      x319 = rho_a**(-1.33333333333333d0)*x4
      x320 = x12*x159
      x321 = rho_a**(-6.0d0)
      x322 = x11*x141
      x323 = exp(-2.8771076053175544d0*x10)
      x324 = x323*x52/x27**3
      x325 = x237*x30
      x326 = 0.166666666666667d0*rho_a**(-3.0d0)*x38
      x327 = x43/x47
      x328 = 0.0322515344331995d0
      x329 = x328/x57**3
      x330 = 0.60570686427738d0*x318*x54
      x331 = x186*x99
      x332 = x192*x194
      x333 = 1.24532330915758d0*x186
      x334 = x54/(x122**3*x59**4)
      x335 = x178*x233*x334
      x336 = rho_a**(-7.33333333333333d0)*x77
      x337 = x336*x88
      x338 = rho_a**(-4)
      x339 = x338*x80
      x340 = gamma_aa*x284
      x341 = x340*x5
      x342 = rho_a**(-6.66666666666667d0)*x5
      x343 = tau_a**(-4)
      x344 = x343*x98
      x345 = 0.00014565673828125d0*gamma_aa**4*x338*x344
      x346 = x92*(-0.010875703125d0*x188 + x345 + 0.1491525d0)
      x347 = x82**(-1.5d0)*x84
      x348 = 1.24378109452736d0*x76
      x349 = x205*x60
      x350 = x197*x227
      x351 = x196*x201
      x352 = 0.0720987654320988d0*x5
      x353 = 1.24378109452736d0*x231
      x354 = tau_a*x149
      x355 = x340*x65
      x356 = x71*(-2.42282745710952d0*x354 + 0.555231292254265d0*x355)
      x357 = 0.360493827160494d0*x5
      x358 = x223*x224
      x359 = x212*x358
      x360 = gamma_aa*x219
      x361 = x225*(-0.40380457618492d0*x209*x219 + 0.134601525394973d0*x
     &210*x360 + 1.07681220315979d0*x215*x354 - 0.246769463224118d0*x215
     &*x355 - x257*(-5.38406101579893d0*x354 + 1.23384731612059d0*x355))
     &
      x362 = 0.60570686427738d0*x68
      x363 = 0.0757133580346725d0*x220
      x364 = x223*x67
      x365 = x364*x70**(-2.5d0)
      x366 = x365*(0.20190228809246d0*x186*x69 - x209*x362*x5 - x360*x36
     &3)
      x367 = 0.00363922642734557d0*x6
      x368 = gamma_aa*x367
      x369 = x368 + x71*(0.00409412973076377d0*x66 - 0.0324444444444444d
     &0)
      x370 = 1.24378109452736d0*x369
      x371 = -0.0011714821875d0*x180 - 0.000579585695450655d0*x182 - 0.0
     &0344377233585999d0*x185 + x202 - x228 + 0.0720987654320988d0* x232
     & + x309*x332 - 0.20190228809246d0*x331
      x372 = x178*x291
      x373 = x106**(-2)
      x374 = x109*x373
      x375 = x243*(0.0252377860115575d0*x115*x6 + 0.0252377860115575d0*x
     &298*x3*(x111* x254 + x253) + x299*x63)
      x376 = -0.0896750813832074d0*x118*(x242*x375 + x266) + x272*(x269*
     &x375 + x270) + x288
      x377 = x291*x293
      x378 = -0.229406311391776d0*x147*x23 + 0.637880379497189d0*x150*x2
     &4 - x320*( -x153 + x154 + x156)
      x379 = gamma_aa*x78
      x380 = 0.00856659785039798d0*gamma_aa
      x381 = gamma_aa*x83
      x382 = 2.48756218905473d0*x381
      x383 = 0.466996240934092d0*x282
      x384 = 0.188341686653414d0*x192*x287 - 0.179350162766415d0*x268 + 
     &2.48756218905473d0*x271*x85 - 0.188341686653414d0*x274 - 0.0010928
     &0054804105d0*x277 + x278*x382 + x281*x383 - 0.000540658297995013d0
     &*x379 - x380*x95
      x385 = gamma_aa*x181
      x386 = -x189 + 0.0497175d0
      x387 = x195*x88
      x388 = 0.00932203125d0*x188
      x389 = x196*x62*x84
      x390 = -x262 + x263*x65
      x391 = x390*x71
      x392 = x248 - x251 - x258
      x393 = x225*x392
      x394 = x214*x391 + 0.180246913580247d0*x393 + x74
      x395 = x197*x394
      x396 = 1.24378109452736d0*x351
      x397 = x352*x391 + x367 + 0.0720987654320988d0*x393
      x398 = x192*x386
      x399 = x198*x83
      x400 = x230*x391 + x261 + 2.0d0*x393
      x401 = 0.0720987654320988d0*x72
      x402 = -gamma_aa*x399 + 0.00344377233585999d0*gamma_aa*x95 + 0.075
     &7133580346725d0*x274 + 0.0004393058203125d0*x277 + 0.0757133580346
     &725d0*x286*x398 + 0.000217344635793995d0*x379 - x394*x85 + x400*x4
     &01
      x403 = 0.622661654578789d0*x186*x203
      x404 = x209*x62
      x405 = x71*(-x211 + 1.21141372855476d0*x404)
      x406 = x358*x390
      x407 = x224*x392
      x408 = x212*x214
      x409 = x245*x362
      x410 = x255*x256
      x411 = -0.0757133580346725d0*x250 - x363*x410 + x409*x63
      x412 = x365*x411
      x413 = 0.20190228809246d0*x244
      x414 = x219*x413
      x415 = x225*(x150*x217 + 0.0673007626974867d0*x210*x410 - 0.538406
     &101579893d0* x215*x404 - x216*x255*x77 - x222 - x257*(-x218 + 2.69
     &203050789947d0*x404) + x414*x63)
      x416 = x212*x230
      x417 = rho_a**(-5.0d0)
      x418 = x67*x77
      x419 = x255*x418*(-0.050475572023115d0*x220 + 0.40380457618492d0*x
     &252)
      x420 = x390*x392
      x421 = 1d0/x70
      x422 = x421*x67
      x423 = x392*x411*x422
      x424 = x297*x65
      x425 = x247 + x256*x424
      x426 = x296*x425
      x427 = x242*x426 + x295
      x428 = x269*x426 + x302
      x429 = x301*x427 - x303*x428 + x312
      x430 = x422*x425
      x431 = 2.0d0*x430
      x432 = 0.180246913580247d0*x430
      x433 = x71*(x302 - x432)
      x434 = x199*x307 + x311 + x401*x71*(x295 - x431) - x433*x85
      x435 = x183*x305
      x436 = x194*x310
      x437 = x345 - x388 + 0.099435d0
      x438 = x304*x79
      x439 = x196*x306
      x440 = 0.0327530378461102d0*x245
      x441 = x421*x425
      x442 = x422*(-0.203822669684853d0*tau_a*x147 + 0.336503813487433d0
     &*x249 + 0.110403946079296d0*x379*x65 - x414)
      x443 = x223*x421
      x444 = 0.0764335011318201d0*gamma_aa*x147*x65 + x409
      x445 = x364*x444/x70**2
      x446 = gamma_aa*x438
      x447 = 0.103185226527957d0*x285*x305*x98
      x448 = 1.3628404446241d0*x117
      x449 = x197*x428
      x450 = 0.0720987654320988d0*x113*x426 + x440
      x451 = x264*x425
      x452 = x77*(-0.203822669684853d0*x170*x63 + x255*x413 + x424)
      x453 = x259*x296*x444
      x454 = 1.3628404446241d0*x243
      x455 = -x117*(0.00409412973076377d0*x112 + 0.0324444444444444d0) +
     & x368
      x456 = x123*x60
      x457 = x343*x79*x88
      x458 = x170*x418
      x459 = x245*x425
      x460 = 2.72568088924821d0*x224
      Ea = -x0*(x105*x55 - x2*x53)
      d1Ea(1) = -rho_a*(x105*x169 - 2.42282745710952d0*x128 + x143*x144 
     &- x161*x167 + x236) + x110*x2*(1.804d0 + 0.804d0/x106) - x125*x56*
     &(x124 + x58)
      d1Ea(2) = -rho_a*(x165*x239 - x294*(x241 + x292*x293))
      d1Ea(3) = -x313*x315*x55
      d2Ea(1) = rho_a*(-gamma_aa*x142*x164*x317*x322 + 5.0d0*rho_a**(-2.
     &0d0)*x105*x142 - 5.84385982654816d0*x1*x161*x324*x88*(0.4588126227
     &83553d0*x148 - 1.27576075899438d0*x152 + 0.478410284622892d0*x158*
     &x6) - 5.15926132639785d0*x105*x238*x8 - 4.84565491421904d0*x126*x1
     &43 + 3.89590655103211d0*x126*x160*x316 - 1.61521830473968d0*x127*x
     &318 - 2.75160604074552d0*x137*x236 + x142*x144*x28*(-x134*(-x133*x
     &42 + x139) + x140*(-0.366880805432736d0*mu*x325 + 4.84565491421904
     &d0 *x126*x15*x327*x42 + 0.733761610865473d0*x131*x319*x42 - 4.4025
     &6966519284d0*x132*x136 - 0.80760915236984d0*x138*x318 + x326*x47 -
     & 2.66666666666667d0*x327*x36*x37 + x44*( 0.244587203621824d0*x325 
     &- x326*x8)) + 0.122293601810912d0*x237* x50) + x167*x31*(0.2200088
     &81071042d0*x145*x321*x8*x91 - 1.60584417974243d0*x145*x78 + 0.4588
     &12622783553d0*x147*x157 - 1.27576075899438d0*x150*x158 + 2.3388947
     &2482303d0*x151*x285 - x320*(0.00448671751316577d0*x126*x16 + x155 
     &+ 0.0407645339369707d0*x319)) + x56*(0.36240302520243d0*gamma_aa* 
     &x146*x172 - 0.0133180320895522d0*rho_a**(-7.0d0)*x329*x88 - 1.4609
     &6495663704d0*x126*x235*x290 - x174*x330 + x330*(x124 + 1.804d0) + 
     &x335*(-0.00291413479477612d0*x180 - 0.00144175546132004d0*x182 - 0
     &.00856659785039798d0*x185 + 0.188341686653414d0*x191*x332 + 2.4875
     &6218905473d0*x202 + x206* x333 - 2.48756218905473d0*x228 + 0.17935
     &0162766415d0*x232 - 0.502244497742438d0*x331) - x372*(0.0131136065
     &764925d0*rho_a**( -10.0d0)*x93 - 0.502244497742438d0*x191*x192*x34
     &2 + x193*x342* x346 - 1.24378109452736d0*x199*x88*(0.0484078840501
     &527d0*x336 + 0.0084375d0*x339) + x200**2*x347*x348 + x203*x333*x37
     &1 - 2.28309273345556d0*x206*x341 - x272*(x214*x356 + 0.08895886822
     &40029d0*x341 + x357*x359 + 0.180246913580247d0*x361 + 0.1802469135
     &80247d0*x366) + 0.46757528851585d0*x337*x349 + 0.00456555896084678
     &d0*x337 + 0.012849896775597d0*x339*x88 + 0.920781579194469d0*x341*
     &x99 + 2.48756218905473d0*x350*x351 + x353*(-0.00970460380625487d0*
     &x186 + x213*x352 + 0.0720987654320988d0*x226) + x370*(x230*x356 + 
     &0.987077852896471d0*x341 + 4.0d0*x359*x5 + 2.0d0*x361 + 2.0d0* x36
     &6)))) - 7.5d0*x105*x168 + 4.84565491421904d0*x128 - 14.53696474265
     &71d0*x143*x31 - 2.0d0*x236 + x316*x317
      d2Ea(2) = rho_a*(-2.79576264291519d0*gamma_aa*x109*x238*x24*x323*x
     &378/x106**3 + x110*x163*x373*x378 + 1.37580302037276d0*x137*x56*( 
     &0.0453003781503037d0*x240 + x372*x376) - 0.465960440485865d0*x170 
     &*x171*x25*x374 + 0.0967061140270375d0*x23*x237*x322*x373*x45 - x56
     &*(-0.00499426203358209d0*gamma_aa*x321*x329 - 0.73048247831852d0*x
     &126*x292 + 0.105700882350709d0*x173 + x335* x384 + x372*(-x150*x18
     &7 + x150*x207 - 0.00437120219216418d0*x179* x276 - x184*x380 + x19
     &0*x387 + x199*x279*(0.015286700226364d0* x181 + 0.005625d0*x184) -
     & x200*x348*x389 - x272*(-x208 + x214* x405 + x214*x406 + x407*x408
     & + 0.180246913580247d0*x412 + 0.180246913580247d0*x415) - x280*x35
     &0 - x283*x371 - 0.251122248871219d0*x332*x386*x88 - 0.175340733193
     &444d0*x349*x385 + x353*x397 + x370*(-x229 + x230*x405 + x230*x406 
     &+ x407*x416 + 2.0d0*x412 + 2.0d0*x415) - 0.00144175546132004d0*x38
     &5 + x387*( -x345 + x388 - 0.099435d0) + x395*x396 + x402*x403))) -
     & x239*x374 + x294*(x241 + x376*x377)
      d2Ea(3) = -rho_a*(0.668760800867552d0*x20*x23**2*x324*x328*x417 - 
     &x294*(x289*x293* x334*x384 - 0.000936424131296642d0*x329*x417 - x3
     &77*( 0.00109280054804105d0*gamma_aa*x90*x91 + x193*x341*( -0.00776
     &8359375d0*x188 + x345 + 0.0497175d0) - x224*x272*(x357* x420 + 0.1
     &80246913580247d0*x419 + 0.180246913580247d0*x423) + 2.487562189054
     &73d0*x224*x369*(x230*x420 + x419 + x423) + x273 + x275 + 0.1883416
     &86653414d0*x341*x398 + 0.0657527749475415d0*x349* x78 - x382*x395 
     &- x383*x402 + 1.24378109452736d0*x397*x400 - 1.24378109452736d0*x3
     &99 + x87/x88)))
      d2Ea(4) = -x314*(rho_a*x177*(x169*x434 - x54*(5.45136177849642d0*x
     &103*x233*x31* x434*x60 - x175*x434 - x178*(-0.251122248871219d0*x1
     &91*x436 + x197*x396*x433 - 0.00699626865671642d0*x199*x435 + 0.003
     &49813432835821d0*x200*x347*x438*x76 - x272*x71*(x302*x443 - x408*x
     &441 + 0.180246913580247d0*x442 - 0.180246913580247d0*x445 - 0.1364
     &70991025459d0*x6) + 0.00349813432835821d0*x350*x439 + x353* x71*(-
     &0.0720987654320988d0*x430 + x440) + x370*x71*(x295*x443 - x416*x44
     &1 + 2.0d0*x442 - 2.0d0*x445 - 1.51426716069345d0*x6) + x403*x434 +
     & 0.00856659785039798d0*x435 + 0.0941708433267071d0* x436*x437*x92)
     &)) + x125*x290*x429)
      d2Ea(5) = x0*x291*x294*(0.00383298875050529d0*x120*x389*x438 - x19
     &0*x447 - 0.00383298875050529d0*x197*x271*x439 - x267*x448*x450 + 0
     &.255849428582226d0*x282*x313 - 0.00766597750101059d0*x308*x446 - 2
     &.72568088924821d0*x313*x376*x456 + x381*x448*x449 + x437*x447 + 0.
     &00938662444277524d0*x446 - x454*x455*(x230*x451 + x242*x452 - x242
     &*x453 - x259*x295) + x454*x85*(x214*x451 - x259*x302 + x269* x452 
     &- x269*x453))
      d2Ea(6) = x315*x54*(-0.0153319550020212d0*x117*x439*x449 - 0.20637
     &0453055914d0* x285*x344*x346 - 2.72568088924821d0*x296*x427*x450 +
     & 0.0229979325030318d0*x308*x457 + 6.78030071952291d0*x429**2*x456 
     &+ x455*x460*(-x431*x444 + 0.815290678739414d0*x458 + 1.81712059283
     &214d0*x459) - 0.0281598733283257d0*x457 - x460*x85*( -x432*x444 + 
     &0.073476814256762d0*x458 + 0.163765189230551d0*x459 ) - 2.15605617
     &215923d-5*x120*x338*x347/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_2_3(rho_a, gamma_aa, tau_a, m
     &u, Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x435, x436, x437, x438, x439, x440, x441, x442, x443, x444, x445
     &, x446, x447, x448, x449, x450, x451, x452, x453, x454, x455, x456
     &, x457, x458, x459, x460, x461, x462, x463, x464, x465, x466, x467
     &, x468, x469, x470, x471, x472, x473, x474, x475, x476, x477, x478
     &, x479, x480, x481, x482, x483, x484, x485, x486, x487, x488, x489
     &
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.33333333333333d0
      x1 = 0.682784063255296d0
      x2 = 7.26848237132856d0*x1
      x3 = 1d0/x0
      x4 = 21.1705915781930d0
      x5 = 1d0/x4
      x6 = mu**4
      x7 = rho_a**0.666666666666667d0
      x8 = 4.60115111447049d0
      x9 = mu**2
      x10 = 1d0/x9
      x11 = x10*x8
      x12 = 3.30192724889463d0*x11*x7
      x13 = exp(x12)
      x14 = x5*x6*(x13 - 1.0d0)
      x15 = x14*x3
      x16 = 0.015286700226364d0*x15
      x17 = 7.0d0*x13 + 11.0d0
      x18 = 1d0/x8
      x19 = x18*x9/x7
      x20 = 0.00280419844572861d0*x17*x19
      x21 = -x16 + x20 + 0.0185185185185185d0
      x22 = rho_a**(-2.0d0)
      x23 = 3.15935190989552d0*x22
      x24 = x10*x23
      x25 = x21*x24
      x26 = 3.0d0*x13
      x27 = x26 - 2.0d0
      x28 = 2.14502939711103d0
      x29 = 1d0/x28
      x30 = rho_a**0.333333333333333d0
      x31 = 1d0/x30
      x32 = x29*x31
      x33 = mu*x32
      x34 = 1d0/mu
      x35 = x28*x34
      x36 = erf(1.8171205928321397d0*x30*x35)
      x37 = x13*x36
      x38 = 0.183440402716368d0*x15 - 0.60570686427738d0*x19*x27 - x26 +
     & 3.90167577847542d0*x33*x37
      x39 = 1d0/x38
      x40 = x39*exp(-1.4385538026587772d0*x19)
      x41 = gamma_aa*x40
      x42 = x25*x41
      x43 = -x42 + 1.0d0
      x44 = 1.804d0 - 0.804d0/x43
      x45 = 0.550321208149104d0*x32
      x46 = 1d0/rho_a
      x47 = 9.86960440108936d0
      x48 = 1d0/x47
      x49 = 0.0833333333333333d0*x48
      x50 = x46*x49
      x51 = x50*x9
      x52 = x45 - x51
      x53 = exp(-x12)
      x54 = mu*x53
      x55 = mu**3
      x56 = -0.825481812223657d0*x33 + 1.77245385090552d0*x36 + x50*x55
      x57 = x52*x54 + x56
      x58 = 0.275160604074552d0*x33
      x59 = x57*x58 - 0.375d0
      x60 = x44*x59
      x61 = 0.682784063255296d0
      x62 = 2.72568088924821d0*x61
      x63 = exp(-4.1274090611182839d0*x33)
      x64 = rho_a**(-2.66666666666667d0)
      x65 = x18*x64
      x66 = gamma_aa*x65
      x67 = 0.0206714418186455d0*x66 + 1.0d0
      x68 = 0.804d0/x67
      x69 = 0.0938662444277524d0*x66 + 1.0d0
      x70 = x69**(-2)
      x71 = 0.050475572023115d0*x66
      x72 = 1d0/gamma_aa
      x73 = tau_a*x72
      x74 = 8.0d0*rho_a*x73
      x75 = x74 - 1.0d0
      x76 = x66*x75
      x77 = 0.0567850185260044d0*x76 - 0.45d0
      x78 = 0.126188930057787d0*x76 - 1.0d0
      x79 = x75*x78
      x80 = x71*x79 + 1.0d0
      x81 = x80**(-0.5d0)
      x82 = x71 + x77*x81
      x83 = x82**2
      x84 = 0.00909806606836394d0*x65
      x85 = gamma_aa*x84
      x86 = x81*(0.0102353243269094d0*x76 - 0.0811111111111111d0) + x85
      x87 = rho_a**(-5.33333333333333d0)*x5
      x88 = rho_a**(-2)
      x89 = tau_a**(-2)
      x90 = 0.0028125d0*x89
      x91 = 0.00286625629244325d0*x87 + x88*x90
      x92 = sqrt(x91)
      x93 = abs(gamma_aa)
      x94 = x92*x93
      x95 = x86*x94
      x96 = 1.24378109452736d0*x95
      x97 = gamma_aa**2
      x98 = x87*x97
      x99 = rho_a**(-8.0d0)
      x100 = 0.0102659822546843d0
      x101 = gamma_aa**3
      x102 = x100*x101
      x103 = x102*x99
      x104 = x88*x89
      x105 = x104*x97
      x106 = 0.015625d0*x105 + 1.0d0
      x107 = x106**(-2)
      x108 = 0.02485875d0*x105*x107 + 0.123456790123457d0
      x109 = x108*x66
      x110 = 0.000182133424673508d0*x103 + 0.0021416494625995d0*x105 + 0
     &.0941708433267071d0*x109 + 0.000135164574498753d0*x98
      x111 = x70*(x110 + 0.0896750813832074d0*x83 - x96) + 1.0d0
      x112 = 1d0/x111
      x113 = -0.804d0*x112
      x114 = x63*(x113 + x68)
      x115 = x42 - 1.0d0
      x116 = -x45 + x51
      x117 = -x116*x54 + x56
      x118 = x117*x58 - 0.375d0
      x119 = x118*x30
      x120 = -x74 + 1.0d0
      x121 = x120*x66
      x122 = 0.0567850185260044d0*x121 + 0.45d0
      x123 = 0.126188930057787d0*x121 + 1.0d0
      x124 = x120*x123
      x125 = x124*x71 + 1.0d0
      x126 = x125**(-0.5d0)
      x127 = -x122*x126 + x71
      x128 = x127**2
      x129 = -x126*(0.0102353243269094d0*x121 + 0.0811111111111111d0) + 
     &x85
      x130 = x129*x94
      x131 = x70*(x110 + 0.0896750813832074d0*x128 - 1.24378109452736d0*
     &x130) + 1.0d0
      x132 = 1d0/x131
      x133 = -0.804d0*x132
      x134 = x30*x62
      x135 = rho_a**(-0.666666666666667d0)
      x136 = x1*x60
      x137 = x135*x136
      x138 = 0.318309886183791d0
      x139 = mu*x138
      x140 = x139*x46
      x141 = 3.75d0*x140
      x142 = rho_a**(-0.333333333333333d0)
      x143 = x34*x53*x8
      x144 = x142*x143
      x145 = 2.20128483259642d0*x144
      x146 = 0.183440402716368d0*x3
      x147 = x22*x49
      x148 = -x146*x29 + x147*x9
      x149 = mu*x29*x3
      x150 = x35*x53
      x151 = 1.21141372855476d0*x135*x150 - x147*x55 + x148*x54 + 0.2751
     &60604074552d0 *x149
      x152 = 0.275160604074552d0*x31
      x153 = 0.0917202013581841d0*x117*x3 - x152*(x116*x145 + x151)
      x154 = x139*x153*x44
      x155 = 7.26848237132856d0*x30
      x156 = x43**(-2)
      x157 = x16 - x20 - 0.0185185185185185d0
      x158 = rho_a**(-3.66666666666667d0)
      x159 = x158*x18
      x160 = x157*x159
      x161 = rho_a**(-3.0d0)
      x162 = x10*x157
      x163 = x161*x162
      x164 = x13*x46
      x165 = 0.0432098765432099d0*x164
      x166 = rho_a**(-1.66666666666667d0)
      x167 = x166*x18
      x168 = x167*x9
      x169 = x13*x168
      x170 = 0.0336503813487433d0*x169
      x171 = rho_a**(-2.33333333333333d0)
      x172 = x14*x171
      x173 = 0.0203822669684853d0*x172
      x174 = 0.00186946563048574d0*x168*x17
      x175 = -x165 + x170 - x173 + x174
      x176 = x10*x175
      x177 = 0.564189583547756d0
      x178 = x177*x46
      x179 = x11*x13
      x180 = 6.60385449778925d0*x179
      x181 = x168*x27
      x182 = x149*x37
      x183 = x13*x135
      x184 = x183*x35*x36
      x185 = x142*x180 + 4.0d0*x164 - 0.40380457618492d0*x169 + 0.244587
     &203621824d0* x172 - 4.72654360241471d0*x178 - 0.40380457618492d0*x
     &181 + 1.30055859282514d0*x182 - 8.58869971286677d0*x184
      x186 = x185*x39
      x187 = x162*x23
      x188 = 3.02993180261165d0*x160 - 6.31870381979104d0*x163 + x176*x2
     &3 + x186*x187
      x189 = x156*x188*x41*x59
      x190 = 5.84385982654816d0*x1
      x191 = x190*x30
      x192 = rho_a**(-3.33333333333333d0)
      x193 = 0.148394183600699d0
      x194 = x193/x67**2
      x195 = x192*x194
      x196 = -x68 + 1.804d0
      x197 = 0.90856029641607d0*x135
      x198 = x197*x61
      x199 = x111**(-2)
      x200 = 2.19144743495556d0*x30
      x201 = rho_a**(-9.0d0)
      x202 = x102*x201
      x203 = rho_a**(-6.33333333333333d0)*x5
      x204 = x203*x97
      x205 = rho_a**(-3)
      x206 = x205*x89
      x207 = x206*x97
      x208 = gamma_aa*x159
      x209 = 0.251122248871219d0*x108
      x210 = x105/x106
      x211 = 0.001553671875d0*x210
      x212 = x211 - 0.0497175d0
      x213 = x101*x212
      x214 = x107*x89
      x215 = 0.0941708433267071d0*x214
      x216 = rho_a**(-5.66666666666667d0)*x18
      x217 = x215*x216
      x218 = 1d0/x92
      x219 = 1d0/x93
      x220 = x219*x86
      x221 = x218*x220
      x222 = 0.00764335011318201d0*x203 + x205*x90
      x223 = x222*x97
      x224 = x221*x223
      x225 = 1d0/x69
      x226 = 0.0001464352734375d0*x103 + 0.00172188616792999d0*x105 + 0.
     &0757133580346725d0*x109 + 0.000108672317896998d0*x98
      x227 = x226 + 0.0720987654320988d0*x83 - x95
      x228 = x225*x227
      x229 = 0.622661654578789d0*x228
      x230 = 0.0242615095156372d0*x159
      x231 = tau_a*x64
      x232 = x158*x75
      x233 = 0.151426716069345d0*x232
      x234 = -gamma_aa*x233 + 0.454280148208035d0*x231
      x235 = x234*x81
      x236 = 0.180246913580247d0*x18
      x237 = x18*x78
      x238 = 0.20190228809246d0*x231
      x239 = 0.0673007626974867d0*x79
      x240 = 0.336503813487433d0*x232
      x241 = x5*(-gamma_aa*x240 + 1.0095114404623d0*x231)
      x242 = x64*x75
      x243 = 0.0252377860115575d0*x242
      x244 = x241*x243
      x245 = -gamma_aa*x244 + x208*x239 - x237*x238
      x246 = x80**(-1.5d0)
      x247 = x246*x77
      x248 = x245*x247
      x249 = -gamma_aa*x230 + x235*x236 + 0.180246913580247d0*x248
      x250 = x249*x94
      x251 = 0.269203050789947d0*x159
      x252 = 2.0d0*x18
      x253 = -gamma_aa*x251 + x235*x252 + 2.0d0*x248
      x254 = x253*x82
      x255 = -0.00145706739738806d0*x202 - 0.000720877730660018d0*x204 -
     & 0.00428329892519899d0*x207 - x208*x209 + x208*x229 + x213*x217 + 
     &1.24378109452736d0*x224 - 1.24378109452736d0*x250 + 0.089675081383
     &2074d0*x254
      x256 = x61*x70
      x257 = x255*x256
      x258 = x63*(0.12080100840081d0*gamma_aa*x195 - x196*x198 + x198*(x
     &113 + 1.804d0 ) + x199*x200*x257)
      x259 = rho_a**(-1.66666666666667d0)
      x260 = x40/x115**2
      x261 = 18.4628097041666d0*x259*x260
      x262 = x10*x21
      x263 = x1*x118
      x264 = x262*x263
      x265 = x171*x194
      x266 = 1.09572371747778d0*x30
      x267 = 2.0d0*x122
      x268 = x125**(-1.5d0)
      x269 = 0.20190228809246d0*x167
      x270 = x123*x269
      x271 = x166*x73
      x272 = 1.0095114404623d0*x271
      x273 = x120*x64
      x274 = gamma_aa*x5
      x275 = x268*(0.0252377860115575d0*x124*x65 + x270*x73 + 0.02523778
     &60115575d0* x273*x274*(x272 + 0.126188930057787d0*x273))
      x276 = 0.10095114404623d0*x65
      x277 = 0.454280148208035d0*x271
      x278 = 0.0567850185260044d0*x273 + x277
      x279 = x126*x278
      x280 = -x252*x279 + x276
      x281 = 0.0896750813832074d0*x127
      x282 = 0.180246913580247d0*x122
      x283 = -x236*x279 + x84
      x284 = 1.24378109452736d0*x94
      x285 = 0.000270329148997507d0*x87
      x286 = x108*x65
      x287 = 0.00428329892519899d0*x104
      x288 = x100*x97
      x289 = x288*x99
      x290 = x129*x219
      x291 = 1.24378109452736d0*gamma_aa
      x292 = x291*x92
      x293 = 0.0720987654320988d0*x128 - x130 + x226
      x294 = x225*x65
      x295 = 0.233498120467046d0*x294
      x296 = rho_a**(-4.66666666666667d0)
      x297 = x18*x296
      x298 = x297*x97
      x299 = x212*x298
      x300 = -gamma_aa*x285 - gamma_aa*x287 + x215*x299 - 0.094170843326
     &7071d0*x286 - 0.000546400274020523d0*x289 + x290*x292 + x293*x295
      x301 = -x281*(x267*x275 + x280) + x284*(x275*x282 + x283) + x300
      x302 = x131**(-2)
      x303 = x256*x302
      x304 = x301*x303
      x305 = 2.0d0*x63
      x306 = -x261*x264 + x305*(0.0226501890751519d0*x265 + x266*x304)
      x307 = 0.90856029641607d0*x167
      x308 = 1d0/x125
      x309 = rho_a**(-4.33333333333333d0)
      x310 = 0.0254778337106067d0*x309
      x311 = x308*(x120*x274*x310 + x270)
      x312 = 0.0720987654320988d0*x126*x127
      x313 = 0.0818825946152754d0*x167
      x314 = x126*x94
      x315 = tau_a**(-3)
      x316 = x315*x97
      x317 = x316*x88
      x318 = 0.0028125d0*x317
      x319 = x218*x290
      x320 = 0.0757133580346725d0*x213
      x321 = x107*x315
      x322 = x297*x320*x321 - 0.00344377233585999d0*x317
      x323 = x318*x319 + x322
      x324 = x312*(-x267*x311 + x307) - x314*(-x282*x311 + x313) + x323
      x325 = x63*x70
      x326 = x0*x302*x325
      x327 = x1*x189
      x328 = 11.6877196530963d0*x30
      x329 = x171*x29
      x330 = mu*x329
      x331 = 0.166666666666667d0*x161*x48
      x332 = 1d0/x55
      x333 = x332*x53
      x334 = rho_a**(-1.33333333333333d0)
      x335 = x156*x41
      x336 = 6.0598636052233d0*x160
      x337 = 12.6374076395821d0*x163
      x338 = 6.31870381979104d0*x22
      x339 = x176*x338
      x340 = x191*x59
      x341 = x38**(-2)
      x342 = x341*exp(-2.8771076053175544d0*x19)
      x343 = x13*x22
      x344 = x65*x9
      x345 = x13*x344
      x346 = x14*x192
      x347 = 1d0/x6
      x348 = x35*x37
      x349 = 0.0322515344331995d0/x67**3
      x350 = 0.60570686427738d0*x259*x61
      x351 = x108*x208
      x352 = x214*x216
      x353 = 1.24532330915758d0*x208
      x354 = x61/(x131**3*x69**4)
      x355 = x200*x255*x354
      x356 = rho_a**(-7.33333333333333d0)*x5
      x357 = x356*x97
      x358 = rho_a**(-4)
      x359 = x358*x89
      x360 = gamma_aa*x296
      x361 = x18*x360
      x362 = rho_a**(-6.66666666666667d0)*x18
      x363 = tau_a**(-4)
      x364 = x107*x363
      x365 = 0.00014565673828125d0*gamma_aa**4*x358*x364
      x366 = x101*(-0.010875703125d0*x210 + x365 + 0.1491525d0)
      x367 = x91**(-1.5d0)*x93
      x368 = 1.24378109452736d0*x86
      x369 = x227*x70
      x370 = x219*x249
      x371 = x218*x223
      x372 = 0.0720987654320988d0*x18
      x373 = 1.24378109452736d0*x253
      x374 = tau_a*x158
      x375 = x360*x75
      x376 = x81*(-2.42282745710952d0*x374 + 0.555231292254265d0*x375)
      x377 = 0.360493827160494d0*x18
      x378 = x245*x246
      x379 = x234*x378
      x380 = gamma_aa*x241
      x381 = x243*x274
      x382 = x247*(-0.40380457618492d0*x231*x241 + 0.134601525394973d0*x
     &232*x380 + 1.07681220315979d0*x237*x374 - 0.246769463224118d0*x237
     &*x375 - x381*(-5.38406101579893d0*x374 + 1.23384731612059d0*x375))
     &
      x383 = 0.60570686427738d0*x78
      x384 = 0.0757133580346725d0*x242
      x385 = x245*x77
      x386 = x385*x80**(-2.5d0)
      x387 = x386*(-x18*x231*x383 + 0.20190228809246d0*x208*x79 - x380*x
     &384)
      x388 = 0.00363922642734557d0*x65
      x389 = gamma_aa*x388
      x390 = x389 + x81*(0.00409412973076377d0*x76 - 0.0324444444444444d
     &0)
      x391 = 1.24378109452736d0*x390
      x392 = -0.0011714821875d0*x202 - 0.000579585695450655d0*x204 - 0.0
     &0344377233585999d0*x207 + x224 - x250 + 0.0720987654320988d0* x254
     & + x320*x352 - 0.20190228809246d0*x351
      x393 = x200*x303
      x394 = -3.02993180261165d0*x159*x21 + 6.31870381979104d0*x161*x262
     & - x186*x25 - x24*(x165 - x170 + x173 - x174)
      x395 = x342/x115**3
      x396 = x269*x78
      x397 = x396*x73
      x398 = x65*x78
      x399 = x398*x75
      x400 = 0.0252377860115575d0*x399
      x401 = 0.126188930057787d0*x242 - x272
      x402 = x381*x401
      x403 = -x397 + x400 + x402
      x404 = x268*x403
      x405 = x267*x404 + x280
      x406 = x282*x404 + x283
      x407 = -x281*x405 + x284*x406 + x300
      x408 = gamma_aa*x87
      x409 = 0.00856659785039798d0*gamma_aa
      x410 = gamma_aa*x92
      x411 = 2.48756218905473d0*x410
      x412 = 0.466996240934092d0*x294
      x413 = -x104*x409 - 0.179350162766415d0*x127*x405 + 0.188341686653
     &414d0*x214* x299 - 0.188341686653414d0*x286 - 0.00109280054804105d
     &0*x289 + x290*x411 + x293*x412 + 2.48756218905473d0*x406*x94 - 0.0
     &00540658297995013d0*x408
      x414 = gamma_aa*x203
      x415 = -x211 + 0.0497175d0
      x416 = x217*x97
      x417 = 0.00932203125d0*x210
      x418 = x218*x72*x93
      x419 = 0.0567850185260044d0*x242 - x277
      x420 = x419*x81
      x421 = x397 - x400 - x402
      x422 = x247*x421
      x423 = x236*x420 + 0.180246913580247d0*x422 + x84
      x424 = x219*x423
      x425 = 1.24378109452736d0*x371
      x426 = x372*x420 + x388 + 0.0720987654320988d0*x422
      x427 = x214*x415
      x428 = x220*x92
      x429 = x252*x420 + x276 + 2.0d0*x422
      x430 = 0.0720987654320988d0*x82
      x431 = 0.00344377233585999d0*gamma_aa*x104 - gamma_aa*x428 + 0.075
     &7133580346725d0*x286 + 0.0004393058203125d0*x289 + 0.0757133580346
     &725d0*x298*x427 + 0.000217344635793995d0*x408 - x423*x94 + x429*x4
     &30
      x432 = 0.622661654578789d0*x208*x225
      x433 = x231*x72
      x434 = x81*(-x233 + 1.21141372855476d0*x433)
      x435 = x378*x419
      x436 = x246*x421
      x437 = x234*x236
      x438 = x167*x383
      x439 = x274*x401
      x440 = -x384*x439 - 0.0757133580346725d0*x399 + x438*x73
      x441 = x386*x440
      x442 = 0.20190228809246d0*x166
      x443 = x241*x442
      x444 = x247*(x159*x239 + 0.0673007626974867d0*x232*x439 - 0.538406
     &101579893d0* x237*x433 - x238*x401*x5 - x244 - x381*(-x240 + 2.692
     &03050789947d0*x433) + x443*x73)
      x445 = x234*x252
      x446 = x5*x77
      x447 = x401*x446*(-0.050475572023115d0*x242 + 0.40380457618492d0*x
     &271)
      x448 = x419*x421
      x449 = 1d0/x80
      x450 = x449*x77
      x451 = x421*x440*x450
      x452 = x310*x75
      x453 = x274*x452 + x396
      x454 = x308*x453
      x455 = x267*x454 + x307
      x456 = x282*x454 + x313
      x457 = x312*x455 - x314*x456 + x323
      x458 = x450*x453
      x459 = 2.0d0*x458
      x460 = 0.180246913580247d0*x458
      x461 = x81*(x313 - x460)
      x462 = x221*x318 + x322 + x430*x81*(x307 - x459) - x461*x94
      x463 = x205*x316
      x464 = x216*x321
      x465 = x365 - x417 + 0.099435d0
      x466 = x315*x88
      x467 = x218*x317
      x468 = 0.0327530378461102d0*x167
      x469 = x449*x453
      x470 = x309*x5
      x471 = x450*(-0.203822669684853d0*tau_a*x470 + 0.336503813487433d0
     &*x398 + 0.110403946079296d0*x408*x75 - x443)
      x472 = x245*x449
      x473 = 0.0764335011318201d0*gamma_aa*x470*x75 + x438
      x474 = x385*x473/x80**2
      x475 = gamma_aa*x466
      x476 = 0.103185226527957d0*x107*x297*x316
      x477 = 1.3628404446241d0*x126
      x478 = x219*x456
      x479 = 0.0720987654320988d0*x122*x454 + x468
      x480 = x278*x453
      x481 = x5*(-0.203822669684853d0*x192*x73 + x401*x442 + x452)
      x482 = x308*x403*x473
      x483 = 1.3628404446241d0*x268
      x484 = -x126*(0.00409412973076377d0*x121 + 0.0324444444444444d0) +
     & x389
      x485 = x132*x70
      x486 = x363*x88*x97
      x487 = x192*x446
      x488 = x167*x453
      x489 = 2.72568088924821d0*x246
      Ea = -x0*(x114*x62 - x2*x60)
      d1Ea(1) = -rho_a*(x114*x141 - 2.42282745710952d0*x137 + x154*x155 
     &- x189*x191 + x258) + x119*x2*(1.804d0 + 0.804d0/x115) - x134*x63*
     &(x133 + x68)
      d1Ea(2) = rho_a*x306
      d1Ea(3) = -x324*x326*x62
      d2Ea(1) = rho_a*(5.0d0*rho_a**(-2.0d0)*x114*x139 - 5.1592613263978
     &5d0*x114*x171* x193*x9 - 4.84565491421904d0*x135*x154 + 3.89590655
     &103211d0*x135* x327 - 1.61521830473968d0*x136*x259 - x139*x153*x18
     &8*x328*x335 + x139*x155*x44*(-x146*(-x145*x52 + x151) + x152*( 4.8
     &4565491421904d0*x135*x333*x4*x52 + 0.733761610865473d0*x143* x334*
     &x52 - 4.40256966519284d0*x144*x148 - 0.80760915236984d0*x150 *x259
     & - 0.366880805432736d0*x330 + x331*x55 - 2.66666666666667d0* x333*
     &x46*x47 + x54*(0.244587203621824d0*x329 - x331*x9)) + 0.1222936018
     &10912d0*x171*x57) - 2.75160604074552d0*x149*x258 - x188*x340*x342*
     &x97*(x162*x186*x338 + x336 - x337 + x339)/x43**3 + x335*x340*(18.9
     &561114593731d0*rho_a**(-4.0d0)*x162 - 17.1696135481327d0*x157*x297
     & + 2.90581327762917d0*x157*x87*x9 + 6.0598636052233d0*x159*x175 - 
     &12.6374076395821d0*x161*x176 + x185 *x187*x341*(13.2077089955785d0
     &*x142*x179 + 8.0d0*x164 - 0.80760915236984d0*x169 + 0.489174407243
     &648d0*x172 - 9.45308720482942d0*x178 - 0.80760915236984d0*x181 + 2
     &.60111718565028d0*x182 - 17.1773994257335d0*x184) + x186*x336 - x1
     &86*x337 + x186*x339 + x187*x39*(-27.0091844287146d0*x10*x334 - 18.
     &9061744096588d0*x164*x332*x36*x47 + 2.86289990428892d0*x166* x348 
     &+ 6.30205813655294d0*x177*x22 + x180*x334 + 14.5369647426571d0*x18
     &3*x347*x4 + 5.72579980857785d0*x259*x348 + 0.673007626974866d0*x27
     &*x344 - 1.73407812376686d0*x330*x37 - 7.55555555555556d0*x343 + 1.
     &21141372855476d0*x345 - 0.57070347511759d0*x346) + x24*(-0.0031157
     &7605080957d0*x17*x344 - 0.0951172458529316d0*x179*x334 + 0.1460905
     &34979424d0*x343 - 0.10095114404623d0*x345 + 0.0475586229264658d0*x
     &346)) + x63*( 0.36240302520243d0*gamma_aa*x194*x309 - 0.0133180320
     &895522d0* rho_a**(-7.0d0)*x349*x97 - 1.46096495663704d0*x135*x257*
     &x302 - x196*x350 + x350*(x133 + 1.804d0) + x355*(-0.00291413479477
     &612d0* x202 - 0.00144175546132004d0*x204 - 0.00856659785039798d0*x
     &207 + 0.188341686653414d0*x213*x352 + 2.48756218905473d0*x224 + x2
     &28* x353 - 2.48756218905473d0*x250 + 0.179350162766415d0*x254 - 0.
     &502244497742438d0*x351) - x393*(0.0131136065764925d0*rho_a**( -10.
     &0d0)*x102 + 0.920781579194469d0*x108*x361 - 0.502244497742438d0*x2
     &13*x214*x362 + x215*x362*x366 - 1.24378109452736d0*x221*x97*(0.048
     &4078840501527d0*x356 + 0.0084375d0*x359) + x222**2*x367*x368 + x22
     &5*x353*x392 - 2.28309273345556d0*x228*x361 - x284*(x236*x376 + 0.0
     &889588682240029d0*x361 + x377*x379 + 0.180246913580247d0*x382 + 0.
     &180246913580247d0*x387) + 0.46757528851585d0*x357*x369 + 0.0045655
     &5896084678d0*x357 + 0.012849896775597d0*x359*x97 + 2.4875621890547
     &3d0*x370*x371 + x373*(-0.00970460380625487d0*x208 + x235*x372 + 0.
     &0720987654320988d0*x248) + x391*(4.0d0*x18*x379 + x252*x376 + 0.98
     &7077852896471d0*x361 + 2.0d0*x382 + 2.0d0*x387 )))) - 7.5d0*x114*x
     &140 + 4.84565491421904d0*x137 - 14.5369647426571d0*x154*x30 - 2.0d
     &0*x258 + x327*x328
      d2Ea(2) = rho_a*(-36.9256194083333d0*gamma_aa*x259*x264*x394*x395 
     &- 6.15426990138888d0*rho_a**(-2.66666666666667d0)*x260*x264 + x119
     &* x190*x260*x394 + x138*x153*x21*x261*x34 + 1.37580302037276d0*x14
     &9 *x63*(x200*x304 + 0.0453003781503037d0*x265) - x63*( -0.00499426
     &203358209d0*gamma_aa*rho_a**(-6.0d0)*x349 - 0.73048247831852d0*x13
     &5*x303*x407 + 0.105700882350709d0*x195 + x355*x413 + x393*(-x159*x
     &209 + x159*x229 - 0.00437120219216418d0* x201*x288 - x206*x409 + x
     &212*x416 + x221*x291*( 0.015286700226364d0*x203 + 0.005625d0*x206)
     & - x222*x368*x418 - x284*(-x230 + x236*x434 + x236*x435 + x436*x43
     &7 + 0.180246913580247d0*x441 + 0.180246913580247d0*x444) - x292*x3
     &70 - x295*x392 - 0.251122248871219d0*x352*x415*x97 - 0.17534073319
     &3444d0*x369*x414 + x373*x426 + x391*(-x251 + x252* x434 + x252*x43
     &5 + x436*x445 + 2.0d0*x441 + 2.0d0*x444) - 0.00144175546132004d0*x
     &414 + x416*(-x365 + x417 - 0.099435d0) + x424*x425 + x431*x432))) 
     &+ x306
      d2Ea(3) = rho_a*(116.661026201793d0*rho_a**(-3.66666666666667d0)*x
     &21**2*x263*x347* x395 - x305*(0.000936424131296642d0*rho_a**(-5.0d
     &0)*x349 + x266* x303*(0.00109280054804105d0*gamma_aa*x100*x99 + x2
     &15*x361*( -0.007768359375d0*x210 + x365 + 0.0497175d0) - x246*x284
     &*(x377* x448 + 0.180246913580247d0*x447 + 0.180246913580247d0*x451
     &) + 2.48756218905473d0*x246*x390*(x252*x448 + x447 + x451) + x285 
     &+ x287 + 0.188341686653414d0*x361*x427 + 0.0657527749475415d0*x369
     &* x87 - x411*x424 - x412*x431 + 1.24378109452736d0*x426*x429 - 1.2
     &4378109452736d0*x428 + x96/x97) - x266*x354*x407*x413))
      d2Ea(4) = -x325*(rho_a*x199*(x141*x462 - x61*(5.45136177849642d0*x
     &112*x255*x30* x462*x70 - x197*x462 - x200*(0.0941708433267071d0*x1
     &01*x464*x465 - 0.251122248871219d0*x213*x464 + x219*x425*x461 - 0.
     &00699626865671642d0*x221*x463 + 0.00349813432835821d0*x222*x367 *x
     &466*x86 - x284*x81*(x313*x472 - x437*x469 + 0.180246913580247d0 *x
     &471 - 0.180246913580247d0*x474 - 0.136470991025459d0*x65) + 0.0034
     &9813432835821d0*x370*x467 + x373*x81*(-0.0720987654320988d0 *x458 
     &+ x468) + x391*x81*(x307*x472 - x445*x469 + 2.0d0*x471 - 2.0d0*x47
     &4 - 1.51426716069345d0*x65) + x432*x462 + 0.00856659785039798d0*x4
     &63))) + x134*x302*x457)
      d2Ea(5) = x0*x303*x305*(0.00383298875050529d0*x129*x418*x466 - x21
     &2*x476 - 0.00383298875050529d0*x219*x406*x467 + 0.255849428582226d
     &0*x294* x324 - 2.72568088924821d0*x301*x324*x485 - 0.0076659775010
     &1059d0* x319*x475 - x405*x477*x479 + x410*x477*x478 + x465*x476 + 
     &0.00938662444277524d0*x475 - x483*x484*(x252*x480 + x267*x481 - x2
     &67*x482 - x307*x403) + x483*x94*(x236*x480 + x282*x481 - x282* x48
     &2 - x313*x403))
      d2Ea(6) = x326*x61*(-0.0153319550020212d0*x126*x467*x478 - 0.20637
     &0453055914d0* x297*x364*x366 - 2.72568088924821d0*x308*x455*x479 +
     & 0.0229979325030318d0*x319*x486 + 6.78030071952291d0*x457**2*x485 
     &+ x484*x489*(-x459*x473 + 0.815290678739414d0*x487 + 1.81712059283
     &214d0*x488) - 0.0281598733283257d0*x486 - x489*x94*( -x460*x473 + 
     &0.073476814256762d0*x487 + 0.163765189230551d0*x488 ) - 2.15605617
     &215923d-5*x129*x358*x367/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRX_TPSS_ERF_case_3(rho_a, gamma_aa, tau_a, mu,
     & Ea, d1Ea, d2Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x402, x403, x404, x405
      Ea = 0.0d0
      d1Ea(:) = 0.0d0
      d2Ea(:) = 0.0d0
      x0 = rho_a**1.0d0
      x1 = 4.60115111447049d0
      x2 = 1d0/x1
      x3 = mu**2
      x4 = rho_a**0.666666666666667d0
      x5 = x2*x3/x4
      x6 = exp(-1.4385538026587772d0*x5)
      x7 = rho_a**(-2.66666666666667d0)
      x8 = x2*x7
      x9 = gamma_aa*x8
      x10 = 1d0/x3
      x11 = 0.183440402716368d0*x1*x4
      x12 = rho_a**1.33333333333333d0
      x13 = 21.1705915781930d0
      x14 = x10*x13
      x15 = 0.010095114404623d0*x12*x14
      x16 = 97.4090910340024d0
      x17 = mu**(-4)
      x18 = x16*x17
      x19 = 0.147857142857143d0*rho_a**2.0d0*x18
      x20 = -x11 + x15 + x19
      x21 = x10*x20
      x22 = 0.239205142311446d0*x21*x6*x9
      x23 = -x22 + 1.0d0
      x24 = 3.14159265358979d0
      x25 = x10*x24
      x26 = x25*(1.804d0 - 0.804d0/x23)
      x27 = 1.0d0*x26
      x28 = rho_a**0.333333333333333d0
      x29 = 0.466194077035412d0*mu
      x30 = exp(-4.1274090611182839d0*x29/x28)
      x31 = 0.0206714418186455d0*x9 + 1.0d0
      x32 = 0.804d0/x31
      x33 = 0.0938662444277524d0*x9 + 1.0d0
      x34 = x33**(-2)
      x35 = 0.050475572023115d0*x9
      x36 = 1d0/gamma_aa
      x37 = tau_a*x36
      x38 = 8.0d0*rho_a*x37
      x39 = x38 - 1.0d0
      x40 = x39*x9
      x41 = 0.0567850185260044d0*x40 - 0.45d0
      x42 = 0.126188930057787d0*x40 - 1.0d0
      x43 = x39*x42
      x44 = x35*x43 + 1.0d0
      x45 = x44**(-0.5d0)
      x46 = x35 + x41*x45
      x47 = x46**2
      x48 = 0.00909806606836394d0*x8
      x49 = gamma_aa*x48
      x50 = x45*(0.0102353243269094d0*x40 - 0.0811111111111111d0) + x49
      x51 = 1d0/x13
      x52 = rho_a**(-5.33333333333333d0)*x51
      x53 = rho_a**(-2)
      x54 = tau_a**(-2)
      x55 = 0.0028125d0*x54
      x56 = 0.00286625629244325d0*x52 + x53*x55
      x57 = sqrt(x56)
      x58 = abs(gamma_aa)
      x59 = x57*x58
      x60 = x50*x59
      x61 = 1.24378109452736d0*x60
      x62 = gamma_aa**2
      x63 = x52*x62
      x64 = rho_a**(-8.0d0)
      x65 = 1d0/x16
      x66 = gamma_aa**3
      x67 = x65*x66
      x68 = x64*x67
      x69 = x53*x54
      x70 = x62*x69
      x71 = 0.015625d0*x70 + 1.0d0
      x72 = x71**(-2)
      x73 = 0.02485875d0*x70*x72 + 0.123456790123457d0
      x74 = x73*x9
      x75 = 0.000135164574498753d0*x63 + 0.000182133424673508d0*x68 + 0.
     &0021416494625995d0*x70 + 0.0941708433267071d0*x74
      x76 = x34*(0.0896750813832074d0*x47 - x61 + x75) + 1.0d0
      x77 = 1d0/x76
      x78 = -0.804d0*x77
      x79 = x30*(x32 + x78)
      x80 = 0.682784063255296d0
      x81 = x28*x80
      x82 = 2.72568088924821d0*x81
      x83 = x0*x27 + x79*x82
      x84 = x11 - x15 - x19
      x85 = rho_a**(-4.33333333333333d0)
      x86 = x51*x85
      x87 = x84*x86
      x88 = rho_a**(-3.66666666666667d0)
      x89 = x2*x88
      x90 = x10*x84
      x91 = x89*x90
      x92 = 0.122293601810912d0*rho_a**(-0.333333333333333d0)*x1
      x93 = 0.0134601525394973d0*rho_a**0.333333333333333d0*x14
      x94 = 0.295714285714286d0*x18
      x95 = x0*x94
      x96 = x92 - x93 - x95
      x97 = x10*x96
      x98 = 0.239205142311446d0*x8
      x99 = x25*(0.229406311391776d0*x87 - 0.637880379497189d0*x91 + x97
     &*x98)
      x100 = gamma_aa*x6/x23**2
      x101 = x100*x99
      x102 = 0.804d0*x0
      x103 = mu/x24
      x104 = x103/x0
      x105 = 3.75d0*x104
      x106 = rho_a**(-3.33333333333333d0)
      x107 = 0.148394183600699d0
      x108 = x107/x31**2
      x109 = x106*x108
      x110 = -x32 + 1.804d0
      x111 = rho_a**(-0.666666666666667d0)
      x112 = 0.90856029641607d0*x111
      x113 = x112*x80
      x114 = x76**(-2)
      x115 = 2.19144743495556d0*x28
      x116 = rho_a**(-9.0d0)
      x117 = x116*x67
      x118 = rho_a**(-6.33333333333333d0)*x51
      x119 = x118*x62
      x120 = rho_a**(-3)
      x121 = x120*x54
      x122 = x121*x62
      x123 = gamma_aa*x89
      x124 = 0.251122248871219d0*x73
      x125 = x70/x71
      x126 = 0.001553671875d0*x125
      x127 = x126 - 0.0497175d0
      x128 = x127*x66
      x129 = x54*x72
      x130 = 0.0941708433267071d0*x129
      x131 = rho_a**(-5.66666666666667d0)*x2
      x132 = x130*x131
      x133 = 1d0/x57
      x134 = 1d0/x58
      x135 = x134*x50
      x136 = x133*x135
      x137 = 0.00764335011318201d0*x118 + x120*x55
      x138 = x137*x62
      x139 = x136*x138
      x140 = 1d0/x33
      x141 = 0.000108672317896998d0*x63 + 0.0001464352734375d0*x68 + 0.0
     &0172188616792999d0*x70 + 0.0757133580346725d0*x74
      x142 = x141 + 0.0720987654320988d0*x47 - x60
      x143 = x140*x142
      x144 = 0.622661654578789d0*x143
      x145 = 0.0242615095156372d0*x89
      x146 = tau_a*x7
      x147 = x39*x88
      x148 = 0.151426716069345d0*x147
      x149 = -gamma_aa*x148 + 0.454280148208035d0*x146
      x150 = x149*x45
      x151 = 0.180246913580247d0*x2
      x152 = x2*x42
      x153 = 0.20190228809246d0*x146
      x154 = 0.0673007626974867d0*x43
      x155 = 0.336503813487433d0*x147
      x156 = x51*(-gamma_aa*x155 + 1.0095114404623d0*x146)
      x157 = x39*x7
      x158 = 0.0252377860115575d0*x157
      x159 = x156*x158
      x160 = -gamma_aa*x159 + x123*x154 - x152*x153
      x161 = x44**(-1.5d0)
      x162 = x161*x41
      x163 = x160*x162
      x164 = -gamma_aa*x145 + x150*x151 + 0.180246913580247d0*x163
      x165 = x164*x59
      x166 = 0.269203050789947d0*x89
      x167 = 2.0d0*x2
      x168 = -gamma_aa*x166 + x150*x167 + 2.0d0*x163
      x169 = x168*x46
      x170 = -0.00145706739738806d0*x117 - 0.000720877730660018d0*x119 -
     & 0.00428329892519899d0*x122 - x123*x124 + x123*x144 + x128*x132 + 
     &1.24378109452736d0*x139 - 1.24378109452736d0*x165 + 0.089675081383
     &2074d0*x169
      x171 = x34*x80
      x172 = x170*x171
      x173 = x30*(0.12080100840081d0*gamma_aa*x109 - x110*x113 + x113*(x
     &78 + 1.804d0 ) + x114*x115*x172)
      x174 = rho_a**(-1.66666666666667d0)
      x175 = 0.682784063255296d0
      x176 = x17*x174*x175*x20
      x177 = x22 - 1.0d0
      x178 = x6/x177**2
      x179 = 0.192320934418403d0*x178
      x180 = rho_a**(-2.33333333333333d0)
      x181 = x108*x180
      x182 = -x38 + 1.0d0
      x183 = x182*x9
      x184 = 0.0567850185260044d0*x183 + 0.45d0
      x185 = 0.126188930057787d0*x183 + 1.0d0
      x186 = x182*x185
      x187 = x186*x35 + 1.0d0
      x188 = x187**(-0.5d0)
      x189 = -x184*x188 + x35
      x190 = x189**2
      x191 = -x188*(0.0102353243269094d0*x183 + 0.0811111111111111d0) + 
     &x49
      x192 = x191*x59
      x193 = x34*(0.0896750813832074d0*x190 - 1.24378109452736d0*x192 + 
     &x75) + 1.0d0
      x194 = x193**(-2)
      x195 = 2.0d0*x184
      x196 = x187**(-1.5d0)
      x197 = x174*x2
      x198 = 0.20190228809246d0*x197
      x199 = x185*x198
      x200 = x174*x37
      x201 = 1.0095114404623d0*x200
      x202 = x182*x7
      x203 = gamma_aa*x51
      x204 = x196*(0.0252377860115575d0*x186*x8 + x199*x37 + 0.025237786
     &0115575d0* x202*x203*(x201 + 0.126188930057787d0*x202))
      x205 = 0.10095114404623d0*x8
      x206 = 0.454280148208035d0*x200
      x207 = 0.0567850185260044d0*x202 + x206
      x208 = x188*x207
      x209 = -x167*x208 + x205
      x210 = 0.0896750813832074d0*x189
      x211 = 0.180246913580247d0*x184
      x212 = -x151*x208 + x48
      x213 = 1.24378109452736d0*x59
      x214 = 0.000270329148997507d0*x52
      x215 = x73*x8
      x216 = 0.00428329892519899d0*x69
      x217 = x62*x65
      x218 = x217*x64
      x219 = x134*x191
      x220 = 1.24378109452736d0*gamma_aa
      x221 = x220*x57
      x222 = x141 + 0.0720987654320988d0*x190 - x192
      x223 = x140*x8
      x224 = 0.233498120467046d0*x223
      x225 = rho_a**(-4.66666666666667d0)
      x226 = x2*x225
      x227 = x226*x62
      x228 = x127*x227
      x229 = -gamma_aa*x214 - gamma_aa*x216 + x130*x228 - 0.094170843326
     &7071d0*x215 - 0.000546400274020523d0*x218 + x219*x221 + x222*x224
      x230 = -x210*(x195*x204 + x209) + x213*(x204*x211 + x212) + x229
      x231 = x194*x230
      x232 = 1.09572371747778d0*x81
      x233 = x232*x34
      x234 = 2.0d0*x30
      x235 = x176*x179 + x234*(0.0226501890751519d0*x181 + x231*x233)
      x236 = 0.90856029641607d0*x197
      x237 = 1d0/x187
      x238 = 0.0254778337106067d0*x85
      x239 = x237*(x182*x203*x238 + x199)
      x240 = 0.0720987654320988d0*x188*x189
      x241 = 0.0818825946152754d0*x197
      x242 = x188*x59
      x243 = tau_a**(-3)
      x244 = x243*x62
      x245 = x244*x53
      x246 = 0.0028125d0*x245
      x247 = x133*x219
      x248 = 0.0757133580346725d0*x128
      x249 = x243*x72
      x250 = x226*x248*x249 - 0.00344377233585999d0*x245
      x251 = x246*x247 + x250
      x252 = x240*(-x195*x239 + x236) - x242*(-x211*x239 + x241) + x251
      x253 = 2.72568088924821d0*x252
      x254 = x30*x34
      x255 = x12*x194
      x256 = x254*x255*x80
      x257 = 1.608d0*x101
      x258 = x10*x98
      x259 = rho_a**(-6.0d0)
      x260 = x102*x25
      x261 = exp(-2.8771076053175544d0*x5)
      x262 = x29/x12
      x263 = 0.0322515344331995d0/x31**3
      x264 = 0.60570686427738d0*rho_a**(-1.66666666666667d0)*x80
      x265 = 1d0/x193
      x266 = x111*x194
      x267 = x123*x73
      x268 = x129*x131
      x269 = 1.24532330915758d0*x123
      x270 = 1/(x193**3*x33**4)
      x271 = x115*x170*x270*x80
      x272 = rho_a**(-7.33333333333333d0)*x51
      x273 = x272*x62
      x274 = rho_a**(-4)
      x275 = x274*x54
      x276 = gamma_aa*x225
      x277 = x2*x276
      x278 = rho_a**(-6.66666666666667d0)*x2
      x279 = tau_a**(-4)
      x280 = x279*x72
      x281 = 0.00014565673828125d0*gamma_aa**4*x274*x280
      x282 = x66*(-0.010875703125d0*x125 + x281 + 0.1491525d0)
      x283 = x56**(-1.5d0)*x58
      x284 = 1.24378109452736d0*x50
      x285 = x142*x34
      x286 = x134*x164
      x287 = x133*x138
      x288 = 0.0720987654320988d0*x2
      x289 = 1.24378109452736d0*x168
      x290 = tau_a*x88
      x291 = x276*x39
      x292 = x45*(-2.42282745710952d0*x290 + 0.555231292254265d0*x291)
      x293 = 0.360493827160494d0*x2
      x294 = x160*x161
      x295 = x149*x294
      x296 = gamma_aa*x156
      x297 = x158*x203
      x298 = x162*(-0.40380457618492d0*x146*x156 + 0.134601525394973d0*x
     &147*x296 + 1.07681220315979d0*x152*x290 - 0.246769463224118d0*x152
     &*x291 - x297*(-5.38406101579893d0*x290 + 1.23384731612059d0*x291))
     &
      x299 = 0.60570686427738d0*x42
      x300 = 0.0757133580346725d0*x157
      x301 = x160*x41
      x302 = x301*x44**(-2.5d0)
      x303 = x302*(0.20190228809246d0*x123*x43 - x146*x2*x299 - x296*x30
     &0)
      x304 = 0.00363922642734557d0*x8
      x305 = gamma_aa*x304
      x306 = x305 + x45*(0.00409412973076377d0*x40 - 0.0324444444444444d
     &0)
      x307 = 1.24378109452736d0*x306
      x308 = -0.0011714821875d0*x117 - 0.000579585695450655d0*x119 - 0.0
     &0344377233585999d0*x122 + x139 - x165 + 0.0720987654320988d0* x169
     & + x248*x268 - 0.20190228809246d0*x267
      x309 = x115*x171
      x310 = x194*x309
      x311 = -0.229406311391776d0*x20*x86 + 0.637880379497189d0*x21*x89 
     &- x258*(-x92 + x93 + x95)
      x312 = x261/x177**3
      x313 = x198*x42
      x314 = x313*x37
      x315 = x42*x8
      x316 = x315*x39
      x317 = 0.0252377860115575d0*x316
      x318 = 0.126188930057787d0*x157 - x201
      x319 = x297*x318
      x320 = -x314 + x317 + x319
      x321 = x196*x320
      x322 = x195*x321 + x209
      x323 = x211*x321 + x212
      x324 = -x210*x322 + x213*x323 + x229
      x325 = gamma_aa*x52
      x326 = 0.00856659785039798d0*gamma_aa
      x327 = gamma_aa*x57
      x328 = 2.48756218905473d0*x327
      x329 = 0.466996240934092d0*x223
      x330 = 0.188341686653414d0*x129*x228 - 0.179350162766415d0*x189*x3
     &22 - 0.188341686653414d0*x215 - 0.00109280054804105d0*x218 + x219*
     &x328 + x222*x329 + 2.48756218905473d0*x323*x59 - 0.000540658297995
     &013d0*x325 - x326*x69
      x331 = gamma_aa*x118
      x332 = -x126 + 0.0497175d0
      x333 = x132*x62
      x334 = 0.00932203125d0*x125
      x335 = x133*x36*x58
      x336 = 0.0567850185260044d0*x157 - x206
      x337 = x336*x45
      x338 = x314 - x317 - x319
      x339 = x162*x338
      x340 = x151*x337 + 0.180246913580247d0*x339 + x48
      x341 = x134*x340
      x342 = 1.24378109452736d0*x287
      x343 = x288*x337 + x304 + 0.0720987654320988d0*x339
      x344 = x129*x332
      x345 = x135*x57
      x346 = x167*x337 + x205 + 2.0d0*x339
      x347 = 0.0720987654320988d0*x46
      x348 = -gamma_aa*x345 + 0.00344377233585999d0*gamma_aa*x69 + 0.075
     &7133580346725d0*x215 + 0.0004393058203125d0*x218 + 0.0757133580346
     &725d0*x227*x344 + 0.000217344635793995d0*x325 - x340*x59 + x346*x3
     &47
      x349 = 0.622661654578789d0*x123*x140
      x350 = x146*x36
      x351 = x45*(-x148 + 1.21141372855476d0*x350)
      x352 = x294*x336
      x353 = x161*x338
      x354 = x149*x151
      x355 = x197*x299
      x356 = x203*x318
      x357 = -x300*x356 - 0.0757133580346725d0*x316 + x355*x37
      x358 = x302*x357
      x359 = 0.20190228809246d0*x174
      x360 = x156*x359
      x361 = x162*(0.0673007626974867d0*x147*x356 - 0.538406101579893d0*
     &x152*x350 - x153*x318*x51 + x154*x89 - x159 - x297*(-x155 + 2.6920
     &3050789947d0*x350) + x360*x37)
      x362 = x149*x167
      x363 = x41*x51
      x364 = x318*x363*(-0.050475572023115d0*x157 + 0.40380457618492d0*x
     &200)
      x365 = x336*x338
      x366 = 1d0/x44
      x367 = x366*x41
      x368 = x338*x357*x367
      x369 = x238*x39
      x370 = x203*x369 + x313
      x371 = x237*x370
      x372 = x195*x371 + x236
      x373 = x211*x371 + x241
      x374 = x240*x372 - x242*x373 + x251
      x375 = x367*x370
      x376 = 2.0d0*x375
      x377 = 0.180246913580247d0*x375
      x378 = x45*(x241 - x377)
      x379 = x136*x246 + x250 + x347*x45*(x236 - x376) - x378*x59
      x380 = x120*x244
      x381 = x131*x249
      x382 = x281 - x334 + 0.099435d0
      x383 = x243*x53
      x384 = x133*x245
      x385 = 0.0327530378461102d0*x197
      x386 = x366*x370
      x387 = x367*(-0.203822669684853d0*tau_a*x86 + 0.336503813487433d0*
     &x315 + 0.110403946079296d0*x325*x39 - x360)
      x388 = x160*x366
      x389 = 0.0764335011318201d0*gamma_aa*x39*x86 + x355
      x390 = x301*x389/x44**2
      x391 = gamma_aa*x383
      x392 = 0.103185226527957d0*x226*x244*x72
      x393 = 1.3628404446241d0*x188
      x394 = x134*x373
      x395 = 0.0720987654320988d0*x184*x371 + x385
      x396 = x207*x370
      x397 = x51*(-0.203822669684853d0*x106*x37 + x318*x359 + x369)
      x398 = x237*x320*x389
      x399 = 1.3628404446241d0*x196
      x400 = -x188*(0.00409412973076377d0*x183 + 0.0324444444444444d0) +
     & x305
      x401 = x265*x34
      x402 = x279*x53*x62
      x403 = x106*x363
      x404 = x197*x370
      x405 = 2.72568088924821d0*x161
      Ea = -rho_a*x83
      d1Ea(1) = -rho_a*(x101*x102 + x105*x79 + x173 + x27) - x83
      d1Ea(2) = rho_a*x235
      d1Ea(3) = -x253*x256
      d2Ea(1) = rho_a*(5.0d0*rho_a**(-2.0d0)*x103*x79 - x100*x260*(2.338
     &89472482303d0* x226*x90 - x258*(0.0407645339369707d0*rho_a**(-1.33
     &333333333333d0 )*x1 + 0.00448671751316577d0*x111*x14 + x94) + 0.22
     &0008881071042d0*x259*x3*x65*x84 - 1.60584417974243d0*x52*x84 + 0.4
     &58812622783553d0*x86*x96 - 1.27576075899438d0*x89*x97) + x102*x261
     &*x62*x99*(0.478410284622892d0*x8*x97 + 0.458812622783553d0*x87 - 1
     &.27576075899438d0*x91)/x23**3 - 5.15926132639785d0*x107*x180*x3*x7
     &9 - 2.75160604074552d0*x173* x262 - x257 + x30*(0.36240302520243d0
     &*gamma_aa*x108*x85 - 0.0133180320895522d0*rho_a**(-7.0d0)*x263*x62
     & - x110*x264 - 1.46096495663704d0*x172*x266 + x264*(-0.804d0*x265 
     &+ 1.804d0) + x271*(-0.00291413479477612d0*x117 - 0.001441755461320
     &04d0*x119 - 0.00856659785039798d0*x122 + 0.188341686653414d0*x128*
     &x268 + 2.48756218905473d0*x139 + x143*x269 - 2.48756218905473d0*x1
     &65 + 0.179350162766415d0*x169 - 0.502244497742438d0*x267) - x310*(
     & 0.0131136065764925d0*rho_a**(-10.0d0)*x67 - 0.502244497742438d0* 
     &x128*x129*x278 + x130*x278*x282 - 1.24378109452736d0*x136*x62*( 0.
     &0484078840501527d0*x272 + 0.0084375d0*x275) + x137**2*x283*x284 + 
     &x140*x269*x308 - 2.28309273345556d0*x143*x277 - x213*(x151*x292 + 
     &0.0889588682240029d0*x277 + x293*x295 + 0.180246913580247d0* x298 
     &+ 0.180246913580247d0*x303) + 0.46757528851585d0*x273*x285 + 0.004
     &56555896084678d0*x273 + 0.012849896775597d0*x275*x62 + 0.920781579
     &194469d0*x277*x73 + 2.48756218905473d0*x286*x287 + x289*(-0.009704
     &60380625487d0*x123 + x150*x288 + 0.0720987654320988d0*x163) + x307
     &*(x167*x292 + 4.0d0*x2*x295 + 0.987077852896471d0*x277 + 2.0d0*x29
     &8 + 2.0d0*x303)))) - x0*x257 - 7.5d0*x104*x79 - 2.0d0*x173 - 2.0d0
     &*x26
      d2Ea(2) = rho_a*(0.384641868836805d0*gamma_aa*x176*x311*x312 + x17
     &*x175*x179*x20* x7 - x178*x260*x311 + 1.37580302037276d0*x262*x30*
     &( 0.0453003781503037d0*x181 + x231*x309) - x30*( -0.00499426203358
     &209d0*gamma_aa*x259*x263 + 0.105700882350709d0* x109 - 0.730482478
     &31852d0*x171*x266*x324 + x271*x330 + x310*( -0.00437120219216418d0
     &*x116*x217 - x121*x326 - x124*x89 + x127* x333 + x136*x220*(0.0152
     &86700226364d0*x118 + 0.005625d0*x121) - x137*x284*x335 + x144*x89 
     &- x213*(-x145 + x151*x351 + x151*x352 + x353*x354 + 0.180246913580
     &247d0*x358 + 0.180246913580247d0*x361) - x221*x286 - x224*x308 - 0
     &.251122248871219d0*x268*x332*x62 - 0.175340733193444d0*x285*x331 +
     & x289*x343 + x307*(-x166 + x167* x351 + x167*x352 + x353*x362 + 2.
     &0d0*x358 + 2.0d0*x361) - 0.00144175546132004d0*x331 + x333*(-x281 
     &+ x334 - 0.099435d0) + x341*x342 + x348*x349))) + x235
      d2Ea(3) = -rho_a*(x234*(0.000936424131296642d0*rho_a**(-5.0d0)*x26
     &3 + x194*x233*( 0.00109280054804105d0*gamma_aa*x64*x65 + x130*x277
     &*( -0.007768359375d0*x125 + x281 + 0.0497175d0) - x161*x213*(x293*
     & x365 + 0.180246913580247d0*x364 + 0.180246913580247d0*x368) + 2.4
     &8756218905473d0*x161*x306*(x167*x365 + x364 + x368) + x214 + x216 
     &+ 0.188341686653414d0*x277*x344 + 0.0657527749475415d0*x285* x52 -
     & x328*x341 - x329*x348 + 1.24378109452736d0*x343*x346 - 1.24378109
     &452736d0*x345 + x61/x62) - x232*x270*x324*x330) + 0.09200831297404
     &85d0*x107*x20**2*x312*x85/mu**6)
      d2Ea(4) = -x254*(rho_a*x114*(x105*x379 - x80*(-x112*x379 - x115*( 
     &-0.251122248871219d0*x128*x381 + x134*x342*x378 - 0.00699626865671
     &642d0*x136*x380 + 0.00349813432835821d0*x137*x283 *x383*x50 - x213
     &*x45*(x241*x388 - x354*x386 + 0.180246913580247d0 *x387 - 0.180246
     &913580247d0*x390 - 0.136470991025459d0*x8) + 0.00349813432835821d0
     &*x286*x384 + x289*x45*(-0.0720987654320988d0 *x375 + x385) + x307*
     &x45*(x236*x388 - x362*x386 + 2.0d0*x387 - 2.0d0*x390 - 1.514267160
     &69345d0*x8) + x349*x379 + 0.00856659785039798d0*x380 + 0.094170843
     &3267071d0*x381*x382*x66) + 5.45136177849642d0*x170*x28*x34*x379*x7
     &7)) + x194*x374*x82)
      d2Ea(5) = x171*x234*x255*(-x127*x392 - 0.00383298875050529d0*x134*
     &x323*x384 + 0.00383298875050529d0*x191*x335*x383 + 0.2558494285822
     &26d0*x223* x252 - x230*x253*x401 - 0.00766597750101059d0*x247*x391
     & - x322* x393*x395 + x327*x393*x394 + x382*x392 + 0.00938662444277
     &524d0* x391 - x399*x400*(x167*x396 + x195*x397 - x195*x398 - x236*
     &x320) + x399*x59*(x151*x396 + x211*x397 - x211*x398 - x241*x320))
      d2Ea(6) = x256*(-0.0153319550020212d0*x188*x384*x394 - 0.206370453
     &055914d0*x226* x280*x282 - 2.72568088924821d0*x237*x372*x395 + 0.0
     &229979325030318d0*x247*x402 + 6.78030071952291d0*x374**2*x401 + x4
     &00*x405*(-x376*x389 + 0.815290678739414d0*x403 + 1.81712059283214d
     &0*x404) - 0.0281598733283257d0*x402 - x405*x59*( -x377*x389 + 0.07
     &3476814256762d0*x403 + 0.163765189230551d0*x404 ) - 2.156056172159
     &23d-5*x191*x274*x283/tau_a**6)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_PW92_ERF(rho_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28
      E = 0.0d0
      x0 = 9.86960440108936d0
      x1 = 2.0d0/x0
      x2 = 0.693147180559945d0
      x3 = -x2 + 1.0d0
      x4 = 1d0/rho_c
      x5 = 0.682784063255296d0*x4**0.333333333333333d0
      x6 = 2.14502939711103d0
      x7 = x4**0.666666666666667d0
      x8 = x7/x6
      x9 = 0.564189583547756d0*sqrt(x4)
      x10 = 0.826307487110758d0*x4**0.166666666666667d0
      x11 = x3*(0.194159335344114d0*x5 + 1.0d0)*log(0.5d0*x0/(x3*( 7.240
     &1019343168302d0*x10 + 3.2595509194222925d0*x5 + 0.4069130045175293
     &2d0*x8 + 1.4187228164796675d0*x9)) + 1.0d0)
      x12 = mu**2
      x13 = 1.41421356237310d0
      x14 = 4.60115111447049d0
      x15 = x4**1.33333333333333d0
      x16 = 3.14159265358979d0
      x17 = x4**1.0d0
      x18 = x17/x16
      x19 = (-0.009578475d0*x18 + 0.0188171922990732d0*x5 + 0.0676322201
     &645715d0*x8 + 1.0d0 + 0.00126674656487366d0*x15/x14)*exp( -0.68361
     &07611867116d0*x5)
      x20 = x4**1.66666666666667d0
      x21 = 0.0837224526465878d0*x20
      x22 = x12*x5
      x23 = 5.57236303607474d0*mu*x10 + 1.0d0
      x24 = mu**3
      x25 = 0.148394183600699d0
      x26 = x19 - 1.0d0
      x27 = 2.28942848510666d0*x14
      x28 = 1.63540853354893d0*x27*(-0.0259335011650457d0*x5 + 1.0d0)/( 
     &0.494402081358784d0*x5 + 0.0524148278841779d0*x8 + 1.0d0)
      E = -rho_c*(x1*x11 - (0.000640558526213818d0*mu**8*x11*x4** 2.6666
     &6666666667d0 + mu**6*(0.0108059179628909d0*x11*x4**2.0d0 + 0.01673
     &02167879722d0*x20*x25*x26) + 0.0223069557172962d0*mu**5* x13*x19*x
     &21 - mu**4*(-0.0683590130188305d0*x11*x15 - 0.13157433081915d0*x18
     &*x26 + 0.00139418473233101d0*x25*x4**1.0d0* (1.63540853354893d0*x2
     &7 - x28 - 10.9027235569928d0*x6*( -0.352521395009435d0*x5 + 0.5580
     &25705063192d0*x8)*exp( -0.49698248213959023d0*x5))) - x1*(x2 - 1.0
     &d0)*log(( 6.7683429898050367d0*x22 + x23 + 3.3926025578013106d0*x2
     &4*x9)/( 3.1331792677937806d0*x22 + x23)) + x13*x24*(0.031505407223
     &1411d0* x17*x19 + 0.00111534778586481d0*x21*(12.0d0*x16*( -4.49737
     &346725955d0*x5 + 0.825481812223657d0*x8)*exp( -0.28165369188898165
     &d0*x5)/x17 + x28/x7)))/(0.508616435555896d0* x12*x8 + 1.0d0)**4)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_PW92_ERF(rho_c, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x9
     &2, x93, x94, x95
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 9.86960440108936d0
      x1 = 2.0d0/x0
      x2 = 1.46459188756152d0
      x3 = 1d0/rho_c
      x4 = x3**0.333333333333333d0
      x5 = x4/x2
      x6 = 0.194159335344114d0*x5 + 1.0d0
      x7 = 0.693147180559945d0
      x8 = -x7 + 1.0d0
      x9 = 2.14502939711103d0
      x10 = x3**0.666666666666667d0
      x11 = x10/x9
      x12 = 0.564189583547756d0*sqrt(x3)
      x13 = 0.826307487110758d0*x3**0.166666666666667d0
      x14 = 0.406913004517529d0*x11 + 1.41872281647967d0*x12 + 7.2401019
     &3431683d0* x13 + 3.25955091942229d0*x5
      x15 = x0/(x14*x8)
      x16 = x8*log(0.5d0*x15 + 1.0d0)
      x17 = x16*x6
      x18 = x1*x17
      x19 = mu**2
      x20 = x11*x19
      x21 = 0.508616435555896d0*x20 + 1.0d0
      x22 = x21**(-4)
      x23 = mu**5
      x24 = x3**1.66666666666667d0
      x25 = 0.0837224526465878d0*x24
      x26 = 1.41421356237310d0
      x27 = exp(-0.6836107611867116d0*x5)
      x28 = 4.60115111447049d0
      x29 = x3**1.33333333333333d0
      x30 = x29/x28
      x31 = 3.14159265358979d0
      x32 = x3**1.0d0
      x33 = x32/x31
      x34 = 0.009578475d0*x33
      x35 = 0.0676322201645715d0*x11 + 0.00126674656487366d0*x30 - x34 +
     & 0.0188171922990732d0*x5 + 1.0d0
      x36 = x27*x35
      x37 = x26*x36
      x38 = x25*x37
      x39 = 0.0223069557172962d0*x23*x38
      x40 = x19*x5
      x41 = 5.57236303607474d0*mu*x13 + 1.0d0
      x42 = mu**3
      x43 = 3.39260255780131d0*x12*x42 + 6.76834298980504d0*x40 + x41
      x44 = x43/(3.13317926779378d0*x40 + x41)
      x45 = x1*(x7 - 1.0d0)
      x46 = x45*log(x44)
      x47 = x3**2.66666666666667d0
      x48 = 0.133841734303777d0*x47
      x49 = 0.00478594012208448d0*x17
      x50 = mu**8*x48*x49
      x51 = x36 - 1.0d0
      x52 = 0.148394183600699d0
      x53 = x24*x52
      x54 = x51*x53
      x55 = x3**2.0d0
      x56 = 1.0525946465532d0*x55
      x57 = 0.0102659822546843d0*x17
      x58 = mu**6*(0.0167302167879722d0*x54 + x56*x57)
      x59 = 0.0315054072231411d0*x32
      x60 = x36*x59
      x61 = exp(-0.28165369188898165d0*x5)
      x62 = x61*(0.825481812223657d0*x11 - 4.49737346725955d0*x5)
      x63 = 12.0d0*x31/x32
      x64 = x62*x63
      x65 = 1d0/x10
      x66 = 2.28942848510666d0
      x67 = x28*x66
      x68 = x65*x67
      x69 = 0.0524148278841779d0*x11 + 0.494402081358784d0*x5 + 1.0d0
      x70 = 1d0/x69
      x71 = -0.0259335011650457d0*x5 + 1.0d0
      x72 = x70*x71
      x73 = 1.63540853354893d0*x72
      x74 = x25*(x64 + x68*x73)
      x75 = x26*x42*(x60 + 0.00111534778586481d0*x74)
      x76 = mu**4
      x77 = 0.13157433081915d0*x33*x51
      x78 = 3.10428814221102d0*x29
      x79 = 0.022020833726518d0*x17
      x80 = exp(-0.49698248213959023d0*x5)
      x81 = x80*(0.558025705063192d0*x11 - 0.352521395009435d0*x5)
      x82 = x81*x9
      x83 = x3**1.0d0*x52*(-x67*x73 + 1.63540853354893d0*x67 - 10.902723
     &5569928d0* x82)
      x84 = x76*(-x77 - x78*x79 + 0.00139418473233101d0*x83)
      x85 = x22*(x39 - x46 + x50 + x58 + x75 - x84)
      x86 = x6*(0.27127533634502d0*x11 + 0.709361408239834d0*x12 + 1.206
     &68365571947d0*x13 + 1.08651697314076d0*x5)/(x14**2*(x15 + 2.0d0))
      x87 = x27*(-0.045088146776381d0*x11 - 0.00168899541983155d0*x30 + 
     &x34 - 0.00627239743302441d0*x5)
      x88 = mu**7
      x89 = 0.928727172679123d0*x13
      x90 = mu*x5
      x91 = x27*(-0.0225440733881905d0*x11 - 0.000844497709915773d0*x30 
     &+ 0.0047892375d0*x33 + 0.113935126864452d0*x35*x5 - 0.003136198716
     &5122d0*x5)
      x92 = 1d0/x4
      x93 = 0.0141372897033723d0*x31*x66*x70*x92
      x94 = 1.09027235569928d0*x68*x72
      x95 = 1.63540853354893d0*x68*x71*(0.034943218589452d0*x11 + 0.1648
     &00693786261d0*x5)/x69**2
      E = -rho_c*(x18 - x85)
      d1E(1) = -mu*x22*(2.83060463816878d-5*x16*x3**3.0d0*x88 + x19*x26*
     &( -0.00111534778586481d0*x25*(2.41662179562687d0*x3**( -0.66666666
     &6666667d0)*x62 + x61*x63*(-0.550321208149104d0*x11 + 1.49912448908
     &652d0*x5) + x64 + x93 + x94 + x95) - 0.00490180588786583d0*x29*x36
     & - x59*x87 + x60 + 0.00185891297644135d0*x74) + x23*(0.00047750955
     &226877d0*x16*x3** 2.33333333333333d0 - 0.0334604335759443d0*x53*x9
     &1 + 0.0278836946466203d0*x54 + 2.10518929310641d0*x55*x57 - 0.1013
     &21183642338d0*x56*x86) - 0.0223069557172962d0*x25*x26*x76* x87 - 0
     &.000290571663240495d0*x3**2.0d0*x37*x76 + 0.037178259528827d0*x38*
     &x76 + x42*(0.00302075971817057d0*x16*x24 + 4.13905085628136d0*x29*
     &x79 - 0.263148661638301d0*x33*x91 + 0.00139418473233101d0*x53*(-1.
     &80615420514536d0*x2*x81*x92 - 10.9027235569928d0*x65*x80*x9*(-0.37
     &2017136708795d0*x11 + 0.117507131669812d0*x5) - 7.26848237132856d0
     &*x65*x82 + 1.09027235569928d0*x68 - x93 - x94 - x95) + x77 - 0.217
     &336917462899d0*x78*x86 - 0.00232364122055169d0*x83) + 0.3569112914
     &76739d0*x47*x49*x88 - 0.0472353356922751d0*x48*x86* x88 + x45*(-1.
     &69630127890066d0*x12*x19 + x44*(x89 + 1.04439308926459d0*x90) - x8
     &9 - 2.25611432993501d0*x90)/x43) + 0.00895469190170509d0*x16*x4 - 
     &x18 - 1.35631049481572d0*x20*(-x39 + x46 - x50 - x58 - x75 + x84)/
     &x21**5 + x85 - 2.0d0*x86
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_PW92_ERF(rho_c, mu, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x171, x172, x173, x174, x175, x176, x177, x178, x179, x180, x181
     &, x182, x183, x184, x185, x186, x187, x188, x189, x190, x191, x192
     &, x193, x194, x195, x196, x197, x198, x199, x200, x201, x202, x203
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 9.86960440108936d0
      x1 = 2.0d0/x0
      x2 = 1.46459188756152d0
      x3 = 1d0/rho_c
      x4 = x3**0.333333333333333d0
      x5 = x4/x2
      x6 = 0.194159335344114d0*x5 + 1.0d0
      x7 = 0.693147180559945d0
      x8 = -x7 + 1.0d0
      x9 = 2.14502939711103d0
      x10 = x3**0.666666666666667d0
      x11 = x10/x9
      x12 = 0.564189583547756d0*sqrt(x3)
      x13 = 1.41872281647967d0*x12
      x14 = 0.826307487110758d0*x3**0.166666666666667d0
      x15 = 0.406913004517529d0*x11 + x13 + 7.24010193431683d0*x14 + 3.2
     &5955091942229d0*x5
      x16 = 1d0/x8
      x17 = x0*x16
      x18 = x17/x15
      x19 = x8*log(0.5d0*x18 + 1.0d0)
      x20 = x19*x6
      x21 = x1*x20
      x22 = mu**2
      x23 = x11*x22
      x24 = 0.508616435555896d0*x23 + 1.0d0
      x25 = x24**(-4)
      x26 = mu**5
      x27 = 1.41421356237310d0
      x28 = exp(-0.6836107611867116d0*x5)
      x29 = 4.60115111447049d0
      x30 = x3**1.33333333333333d0
      x31 = x30/x29
      x32 = 3.14159265358979d0
      x33 = x3**1.0d0
      x34 = x33/x32
      x35 = 0.009578475d0*x34
      x36 = -x35
      x37 = 0.0676322201645715d0*x11 + 0.00126674656487366d0*x31 + x36 +
     & 0.0188171922990732d0*x5 + 1.0d0
      x38 = x28*x37
      x39 = x3**1.66666666666667d0
      x40 = 0.0837224526465878d0*x39
      x41 = x38*x40
      x42 = x27*x41
      x43 = 0.0223069557172962d0*x26*x42
      x44 = x22*x5
      x45 = 5.57236303607474d0*mu*x14 + 1.0d0
      x46 = 3.13317926779378d0*x44 + x45
      x47 = 1d0/x46
      x48 = mu**3
      x49 = 3.39260255780131d0*x12*x48 + 6.76834298980504d0*x44 + x45
      x50 = x47*x49
      x51 = x1*(x7 - 1.0d0)
      x52 = x51*log(x50)
      x53 = x3**2.66666666666667d0
      x54 = 0.133841734303777d0*x53
      x55 = 0.00478594012208448d0*x20
      x56 = mu**8*x54*x55
      x57 = x38 - 1.0d0
      x58 = 0.148394183600699d0
      x59 = x39*x58
      x60 = x57*x59
      x61 = x3**2.0d0
      x62 = 1.0525946465532d0*x61
      x63 = 0.0102659822546843d0*x20
      x64 = mu**6*(0.0167302167879722d0*x60 + x62*x63)
      x65 = 0.179587122125167d0*x33
      x66 = 0.175432441092201d0*x65
      x67 = x38*x66
      x68 = exp(-0.28165369188898165d0*x5)
      x69 = x68*(0.825481812223657d0*x11 - 4.49737346725955d0*x5)
      x70 = x32/x33
      x71 = 12.0d0*x70
      x72 = x69*x71
      x73 = 1d0/x10
      x74 = 2.28942848510666d0
      x75 = x29*x74
      x76 = x73*x75
      x77 = 0.0524148278841779d0*x11 + 0.494402081358784d0*x5 + 1.0d0
      x78 = 1d0/x77
      x79 = -0.0259335011650457d0*x5 + 1.0d0
      x80 = x78*x79
      x81 = 1.63540853354893d0*x80
      x82 = x40*(x72 + x76*x81)
      x83 = x27*x48*(x67 + 0.00111534778586481d0*x82)
      x84 = mu**4
      x85 = x34*x57
      x86 = 0.13157433081915d0*x85
      x87 = 3.10428814221102d0*x30
      x88 = 0.022020833726518d0*x20
      x89 = 1.63540853354893d0*x75
      x90 = exp(-0.49698248213959023d0*x5)
      x91 = x90*(0.558025705063192d0*x11 - 0.352521395009435d0*x5)
      x92 = x9*x91
      x93 = 10.9027235569928d0*x92
      x94 = x75*x81
      x95 = x3**1.0d0*x58
      x96 = x95*(x89 - x93 - x94)
      x97 = x84*(-x86 - x87*x88 + 0.00139418473233101d0*x96)
      x98 = x43 - x52 + x56 + x64 + x83 - x97
      x99 = x25*x98
      x100 = 0.129439556896076d0*x4
      x101 = 0.0691804894611514d0
      x102 = x101*x19
      x103 = 2.0d0*x6
      x104 = x18 + 2.0d0
      x105 = 1d0/x104
      x106 = x15**(-2)
      x107 = 0.27127533634502d0*x11 + 0.709361408239834d0*x12 + 1.206683
     &65571947d0* x14 + 1.08651697314076d0*x5
      x108 = x105*x106*x107
      x109 = -x43 + x52 - x56 - x64 - x83 + x97
      x110 = x24**(-5)
      x111 = x110*x23
      x112 = x109*x111
      x113 = 1.09027235569928d0*x76
      x114 = 1d0/x4
      x115 = x114*x32*x74
      x116 = 0.0141372897033723d0*x115*x78
      x117 = x114*x2
      x118 = x117*x91
      x119 = 1.80615420514536d0*x118
      x120 = x73*x92
      x121 = 7.26848237132856d0*x120
      x122 = x90*(-0.372017136708795d0*x11 + 0.117507131669812d0*x5)
      x123 = x122*x73*x9
      x124 = 10.9027235569928d0*x123
      x125 = x76*x80
      x126 = 1.09027235569928d0*x125
      x127 = x77**(-2)
      x128 = 0.034943218589452d0*x11 + 0.164800693786261d0*x5
      x129 = x127*x128
      x130 = 1.63540853354893d0*x76*x79
      x131 = x129*x130
      x132 = x113 - x116 - x119 - x121 - x124 - x126 - x131
      x133 = 0.00139418473233101d0*x59
      x134 = 0.200908840802769d0*x39
      x135 = 0.0150354743280612d0*x19
      x136 = x37*x5
      x137 = -0.0225440733881905d0*x11 - 0.000844497709915773d0*x31 + 0.
     &0047892375d0* x34 - 0.0031361987165122d0*x5
      x138 = x28*(0.113935126864452d0*x136 + x137)
      x139 = x138*x34
      x140 = x30*x88
      x141 = 0.217336917462899d0
      x142 = x108*x6
      x143 = x141*x142
      x144 = x134*x135 - 0.263148661638301d0*x139 + 4.13905085628136d0*x
     &140 - x143* x87 + x86
      x145 = 0.00168899541983155d0*x31
      x146 = 0.045088146776381d0*x11
      x147 = 0.00627239743302441d0*x5
      x148 = x28*(-x145 - x146 - x147 + x35)
      x149 = x148*x66
      x150 = 0.122619224952946d0*x30
      x151 = x150*x38
      x152 = 0.0399758348639607d0*x151
      x153 = 0.00185891297644135d0*x82
      x154 = -0.550321208149104d0*x11 + 1.49912448908652d0*x5
      x155 = x68*x71
      x156 = 2.14502939711103d0*x3**(-0.666666666666667d0)
      x157 = x156*x69
      x158 = x116 + x126 + x131
      x159 = x40*(x154*x155 + 1.12661476755593d0*x157 + x158 + x72)
      x160 = 0.00111534778586481d0*x159
      x161 = x22*x27
      x162 = 0.928727172679123d0*x14
      x163 = mu*x5
      x164 = 2.25611432993501d0*x163
      x165 = x12*x22
      x166 = 1.69630127890066d0*x165
      x167 = x162 + 1.04439308926459d0*x163
      x168 = -x162 - x164 - x166 + x167*x50
      x169 = 1d0/x49
      x170 = x169*x51
      x171 = x168*x170
      x172 = x27*x84
      x173 = x148*x172
      x174 = x173*x40
      x175 = x3**2.33333333333333d0
      x176 = 0.0681236923205142d0*x175
      x177 = 0.00700944907716014d0*x19
      x178 = x61*x63
      x179 = x138*x59
      x180 = 0.101321183642338d0
      x181 = x180*x62
      x182 = mu**7
      x183 = x3**3.0d0
      x184 = 0.00866220739124163d0*x182*x183
      x185 = 0.00326776364305339d0*x19
      x186 = x182*x53
      x187 = x186*x55
      x188 = x42*x84
      x189 = x172*x38
      x190 = 0.0571643564037363d0*x3**2.0d0
      x191 = x189*x190
      x192 = 0.0472353356922751d0
      x193 = x142*x192
      x194 = x182*x54
      x195 = x184*x185 + 0.356911291476739d0*x187 + 0.037178259528827d0*
     &x188 - 0.00508309165921971d0*x191 - x193*x194 + x26*(-x142*x181 + 
     &x176* x177 + 2.10518929310641d0*x178 - 0.0334604335759443d0*x179 +
     & 0.0278836946466203d0*x60)
      x196 = x161*(-x149 - x152 + x153 - x160 + x67) + x171 - 0.02230695
     &57172962d0* x174 + x195 + x48*(x132*x133 + x144 - 0.00232364122055
     &169d0*x96)
      x197 = mu*x25
      x198 = x103*x105
      x199 = x106*(0.452125560575033d0*x11 + 1.06404211235975d0*x12 + 1.
     &40779759833938d0*x14 + 1.44868929752102d0*x5)
      x200 = x107**2/(x104**2*x15**4)
      x201 = x107*(0.542550672690039d0*x11 + x13 + 2.41336731143894d0*x1
     &4 + 2.17303394628153d0*x5)/x15**3
      x202 = 0.0223069557172962d0*x172*x40
      x203 = x28*(0.0751469112939684d0*x11 + 0.00394098931294027d0*x31 -
     & 0.01915695d0 *x34 + 0.00836319657736588d0*x5)
      x204 = x108*x184
      x205 = x105*x199
      x206 = x192*x194*x6
      x207 = x16*x200
      x208 = x207*x6
      x209 = x105*x201
      x210 = x162 + x164 + x166
      x211 = mu*x167
      x212 = x211*x47
      x213 = 1.08351503479231d0*x14
      x214 = 0.350864882184401d0*x65
      x215 = x154*x68
      x216 = -0.0282745794067445d0*x115*x129 + 0.363424118566428d0*x125 
     &+ x127*x130*( 0.0582386976490866d0*x11 + 0.219734258381682d0*x5) -
     & x128*x130*( 0.0698864371789039d0*x11 + 0.329601387572523d0*x5)/x7
     &7**3 - 2.18054471139857d0*x129*x76*x79
      x217 = x108*x176
      x218 = x181*x6
      x219 = x28*(0.0259624262672375d0*x11*x37 + 0.0375734556469842d0*x1
     &1 - 0.151913502485936d0*x136 + 0.455740507457808d0*x137*x5 + 0.001
     &97049465647014d0*x31 + x36 + 0.00418159828868294d0*x5)
      x220 = x6*x87
      x221 = x141*x220
      E = -rho_c*(x21 - x99)
      d1E(1) = x100*x102 - x103*x108 - 1.35631049481572d0*x112 - x196*x1
     &97 - x21 + x99
      d2E(1) = x3*(0.682784063255296d0*x100*x108 + 0.0862930379307174d0*
     &x102*x4 + x103* x17*x200 + 0.129439556896076d0*x108*x5 - 2.2994726
     &9793409d0*x109* x31*x84/x24**6 - 2.71262098963145d0*x11*x110*x196*
     &x48 + 2.71262098963145d0*x111*x98 + 2.26051749135954d0*x112 - 4.0d
     &0* x142 - 2.0d0*x197*(-x1*x168*x169*x8 - x161*(x149 + x152 - x153 
     &+ x160 - x67) + x195 + x202*x28*(x145 + x146 + x147 + x36) + x48*(
     & -x133*(-x113 + x119 + x121 + x124 + x158) + x144 + 0.002323641220
     &55169d0*x95*(-x89 + x93 + x94))) - x197*(mu*x168* x210*x51/x49**2 
     &+ x161*(-0.0799516697279215d0*x148*x150 + x148* x214 + 0.133252782
     &879869d0*x151 + 0.0037178259528827d0*x159 - x203*x66 - x214*x38 + 
     &0.00111534778586481d0*x40*(-x155*( 0.917202013581841d0*x11 - 1.998
     &83265211535d0*x5) - 2.25322953511185d0*x156*x215 - 0.7510765117039
     &51d0*x157 - 24.0d0* x215*x70 + x216 - 0.154912426780983d0*x3**(-0.
     &333333333333333d0)* x69) - 0.00910930363347549d0*x41 - 0.004957101
     &27051027d0*x82) + x170*(3.00815243991335d0*x163 + 2.54445191835098
     &d0*x165 - 2.0d0* x210*x212 + x211*x49*(1.85745434535825d0*x14 + 2.
     &08878617852919d0 *x163)/x46**2 + x213 - x50*(1.39252411901946d0*x1
     &63 + x213)) - x171*x212 - 0.0101661833184394d0*x173*x190 + 0.07435
     &6519057654d0* x174 - 0.0577480492749442d0*x182*x183*x185 + 0.71382
     &2582953479d0* x186*x193 - 1.30867473541471d0*x187 - 0.099142025410
     &2054d0*x188 - 4.52089344419912d-5*x189*x3**2.33333333333333d0 + 0.
     &023721094409692d0*x191 + 0.466194077035412d0*x194*x208 - x202* x20
     &3 + 0.064503068866399d0*x204 + x205*x206 - x206*x209 + x26*( x101*
     &x217 + 4.21037858621281d0*x142*x180*x61 - 0.363326359042743d0*x175
     &*x177 - 6.31556787931922d0*x178 + 0.111534778586481d0*x179 + x205*
     &x218 + x208*x62 - x209*x218 + 0.0691804894611514d0*x217 - 0.033460
     &4335759443d0*x219*x59 - 0.074356519057654d0*x60) + x48*(0.14839418
     &3600699d0*x108*x134 + 0.200908840802769d0*x108*x59 - 0.00464728244
     &110338d0*x132*x59 + x133*(-3.61230841029072d0*x117*x122 + 2.220446
     &04925031d-16*x118 + 2.42282745710952d0*x120 - 14.5369647426571d0*x
     &123 + x216 - 10.9027235569928d0*x73*x9*x90*(0.620028561181324d0*x1
     &1 - 0.156676175559749d0*x5) - 0.363424118566428d0*x76 - x90*( 0.16
     &6966313186253d0*x11 - 0.105477574079378d0*x5)) - 0.803635363211076
     &d0*x135*x39 + 0.526297323276602d0*x139 - 9.65778533132318d0*x140 +
     & 8.27810171256273d0*x143*x30 + x205*x221 + 2.14502939711103d0*x207
     &*x220 - x209*x221 - 0.263148661638301d0* x219*x34 - 0.263148661638
     &301d0*x85 + 0.00619637658813783d0*x96)) + x198*x199 - x198*x201)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = 0.682784063255296d0*x0**0.333333333333333d0
      x2 = 2.14502939711103d0
      x3 = x0**0.666666666666667d0
      x4 = x3/x2
      x5 = 0.564189583547756d0*sqrt(x0)
      x6 = 0.826307487110758d0*x0**0.166666666666667d0
      x7 = 9.86960440108936d0
      x8 = 0.306852819440055d0
      x9 = x7/x8
      x10 = 2.0d0*(0.194159335344114d0*x1 + 1.0d0)*log(0.5d0*x9/( 3.2595
     &509194222925d0*x1 + 0.40691300451752932d0*x4 + 1.4187228164796675d
     &0*x5 + 7.2401019343168302d0*x6) + 1.0d0)
      x11 = x8/x7
      x12 = x10*x11
      x13 = rho_s**4/rho_c**4
      x14 = rho_s*x0
      x15 = x14 + 1.0d0
      x16 = -x14 + 1.0d0
      x17 = 1.92366105093154d0*x15**1.33333333333333d0 + 1.9236610509315
     &4d0*x16** 1.33333333333333d0 - 3.84732210186307d0
      x18 = 0.0197517897025652d0*x17*(0.101077332976288d0*x1 + 1.0d0)*(-
     &x13 + 1.0d0) *log(1.0d0 + 29.608574643216677d0/(3.2918048099450625
     &d0*x1 + 0.41002507094961249d0*x4 + 0.76232752193529008d0*x5 + 9.87
     &21297225692712d0*x6))
      x19 = (0.186690969707574d0*x1 + 1.0d0)*log(x9/(5.630984149097876d0
     &*x1 + 0.51606646454786342d0*x4 + 2.9152147142191778d0*x5 + 13.4579
     &13714394447d0*x6) + 1.0d0)
      x20 = x11*x13*x17
      x21 = x20*(-x10 + x19)
      x22 = mu**2
      x23 = 1.41421356237310d0
      x24 = x0**1.66666666666667d0
      x25 = 0.0837224526465878d0*x24
      x26 = 1.0d0 - rho_s**2/rho_c**2
      x27 = 4.60115111447049d0
      x28 = x0**1.33333333333333d0/x27
      x29 = 3.14159265358979d0
      x30 = x0**1.0d0
      x31 = x30/x29
      x32 = (0.0188171922990732d0*x1 + 0.00126674656487366d0*x28 - 0.009
     &578475d0*x31 + 0.0676322201645715d0*x4 + 1.0d0)*exp(-0.68361076118
     &67116d0*x1)
      x33 = x26*x32
      x34 = x15**0.666666666666667d0 + x16**0.666666666666667d0
      x35 = x34**3
      x36 = x1*x22/x34**2
      x37 = 11.1447260721495d0*mu*x6/x34 + 1.0d0
      x38 = mu**3
      x39 = 1d0/x3
      x40 = 1d0/x15
      x41 = x40**0.666666666666667d0
      x42 = 0.0524148278841779d0*x4
      x43 = x1*x40**0.333333333333333d0
      x44 = 2.28942848510666d0*x27
      x45 = 0.817704266774463d0*x44
      x46 = x15**2*x45*(-0.0259335011650457d0*x43 + 1.0d0)/(x41*(x41*x42
     & + 0.494402081358784d0*x43 + 1.0d0))
      x47 = 1d0/x16
      x48 = x47**0.666666666666667d0
      x49 = x1*x47**0.333333333333333d0
      x50 = x16**2*x45*(-0.0259335011650457d0*x49 + 1.0d0)/(x48*(x42*x48
     & + 0.494402081358784d0*x49 + 1.0d0))
      x51 = -x12 + x18
      x52 = -x21 + x51
      x53 = 0.148394183600699d0
      x54 = x26*(x32 - 1.0d0)
      E = -rho_c*(x12 - x18 + x21 + 0.25d0*(0.0126441184989504d0*mu**8*x
     &0** 2.66666666666667d0*x52 - 4.0d0*mu**6*(-0.0533250677421794d0*x0
     &** 2.0d0*(x20*(x10 - x19) + x51) + 0.0167302167879722d0*x24*x53*x5
     &4 ) - 0.0892278228691848d0*mu**5*x23*x25*x33 + 4.0d0*mu**4*( -0.00
     &139418473233101d0*x0**1.0d0*x53*(10.9027235569928d0*x2*x26*( -0.35
     &2521395009435d0*x1 + 0.558025705063192d0*x4)*exp( -0.4969824821395
     &9023d0*x1) - 0.817704266774463d0*x44*(x15** 2.66666666666667d0 + x
     &16**2.66666666666667d0) + x46 + x50) + 1.55214407110551d0*x28*x52 
     &- 0.13157433081915d0*x31*x54) - x11* x35*log((27.073371959220147d0
     &*x36 + x37 + 27.140820462410485d0* x38*x5/x35)/(12.532717071175123
     &d0*x36 + x37)) - 4.0d0*x23*x38*( 0.00111534778586481d0*x25*(12.0d0
     &*x26*x29*(-4.49737346725955d0*x1 + 0.825481812223657d0*x4)*exp(-0.
     &28165369188898165d0*x1)/x30 + x39*x46 + x39*x50) + 0.0315054072231
     &411d0*x30*x33))/( 0.508616435555896d0*x22*x4 + 1.0d0)**4)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x204, x205, x206, x207, x208, x209, x210, x211, x212, x213, x214
     &, x215, x216, x217, x218, x219, x220, x221, x222, x223, x224, x225
     &, x226, x227
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 2.14502939711103d0
      x1 = 1d0/rho_c
      x2 = x1**0.666666666666667d0
      x3 = x2/x0
      x4 = 0.564189583547756d0*sqrt(x1)
      x5 = 1.46459188756152d0
      x6 = x1**0.333333333333333d0
      x7 = x6/x5
      x8 = 0.826307487110758d0*x1**0.166666666666667d0
      x9 = 0.406913004517529d0*x3 + 1.41872281647967d0*x4 + 3.2595509194
     &2229d0*x7 + 7.24010193431683d0*x8
      x10 = 9.86960440108936d0
      x11 = 0.306852819440055d0
      x12 = x10/x11
      x13 = x12/x9
      x14 = log(0.5d0*x13 + 1.0d0)
      x15 = 0.388318670688228d0*x7 + 2.0d0
      x16 = x14*x15
      x17 = x11/x10
      x18 = x16*x17
      x19 = rho_s*x1
      x20 = x19 + 1.0d0
      x21 = -x19 + 1.0d0
      x22 = 1.92366105093154d0*x20**1.33333333333333d0 + 1.9236610509315
     &4d0*x21** 1.33333333333333d0 - 3.84732210186307d0
      x23 = rho_s**4/rho_c**4
      x24 = -x23 + 1.0d0
      x25 = 0.101077332976288d0*x7 + 1.0d0
      x26 = 0.410025070949612d0*x3 + 0.76232752193529d0*x4 + 3.291804809
     &94506d0*x7 + 9.87212972256927d0*x8
      x27 = 1.0d0 + 29.6085746432167d0/x26
      x28 = log(x27)
      x29 = x25*x28
      x30 = x24*x29
      x31 = 0.0197517897025652d0*x22*x30
      x32 = 0.186690969707574d0*x7 + 1.0d0
      x33 = 0.516066464547863d0*x3 + 2.91521471421918d0*x4 + 5.630984149
     &09788d0*x7 + 13.4579137143944d0*x8
      x34 = x12/x33 + 1.0d0
      x35 = log(x34)
      x36 = x32*x35
      x37 = -x16 + x36
      x38 = x22*x23
      x39 = x17*x38
      x40 = x37*x39
      x41 = 1.41421356237310d0
      x42 = mu**5
      x43 = rho_s**2/rho_c**2
      x44 = -x43 + 1.0d0
      x45 = exp(-0.6836107611867116d0*x7)
      x46 = 4.60115111447049d0
      x47 = x1**1.33333333333333d0
      x48 = x47/x46
      x49 = 3.14159265358979d0
      x50 = x1**1.0d0
      x51 = x50/x49
      x52 = -0.009578475d0*x51
      x53 = 0.0676322201645715d0*x3 + 0.00126674656487366d0*x48 + x52 + 
     &0.0188171922990732d0*x7 + 1.0d0
      x54 = x45*x53
      x55 = x44*x54
      x56 = x41*x42*x55
      x57 = x1**1.66666666666667d0
      x58 = 0.0837224526465878d0*x57
      x59 = 0.0892278228691848d0*x58
      x60 = x20**0.666666666666667d0 + x21**0.666666666666667d0
      x61 = x60**3
      x62 = mu**2
      x63 = x60**2
      x64 = 1d0/x63
      x65 = x62*x64
      x66 = x65*x7
      x67 = 1d0/x60
      x68 = mu*x67
      x69 = 11.1447260721495d0*x68*x8 + 1.0d0
      x70 = mu**3
      x71 = x4/x61
      x72 = 27.0733719592201d0*x66 + x69 + 27.1408204624105d0*x70*x71
      x73 = x72/(12.5327170711751d0*x66 + x69)
      x74 = x17*log(x73)
      x75 = 0.179587122125167d0*x50
      x76 = 0.175432441092201d0*x75
      x77 = x55*x76
      x78 = 0.825481812223657d0*x3 - 4.49737346725955d0*x7
      x79 = exp(-0.28165369188898165d0*x7)
      x80 = x44*x78*x79
      x81 = x49/x50
      x82 = 12.0d0*x81
      x83 = x80*x82
      x84 = x20**2
      x85 = 1d0/x20
      x86 = x85**0.666666666666667d0
      x87 = 0.0524148278841779d0*x3
      x88 = x85**0.333333333333333d0
      x89 = x7*x88
      x90 = x86*x87 + 0.494402081358784d0*x89 + 1.0d0
      x91 = 1d0/x90
      x92 = x84*x91
      x93 = (-0.0259335011650457d0*x89 + 1.0d0)/x86
      x94 = x92*x93
      x95 = 1d0/x2
      x96 = 2.28942848510666d0
      x97 = x46*x96
      x98 = x95*x97
      x99 = 0.817704266774463d0*x98
      x100 = x21**2
      x101 = 1d0/x21
      x102 = x101**0.666666666666667d0
      x103 = x101**0.333333333333333d0
      x104 = x103*x7
      x105 = x102*x87 + 0.494402081358784d0*x104 + 1.0d0
      x106 = 1d0/x105
      x107 = x100*x106
      x108 = (-0.0259335011650457d0*x104 + 1.0d0)/x102
      x109 = x107*x108
      x110 = x58*(x109*x99 + x83 + x94*x99)
      x111 = 4.0d0*x41*x70
      x112 = -x18 + x31
      x113 = x112 - x40
      x114 = 0.0472353356922751d0
      x115 = mu**8
      x116 = x1**2.66666666666667d0
      x117 = 0.267683468607554d0*x114*x115*x116
      x118 = 0.148394183600699d0
      x119 = x118*x57
      x120 = x54 - 1.0d0
      x121 = x120*x44
      x122 = x119*x121
      x123 = 0.101321183642338d0
      x124 = x1**2.0d0
      x125 = x16 - x36
      x126 = x125*x39
      x127 = x112 + x126
      x128 = x123*x124*x127
      x129 = 4.0d0*mu**6
      x130 = -0.13157433081915d0*x121*x51
      x131 = x97*(x20**2.66666666666667d0 + x21**2.66666666666667d0)
      x132 = 0.817704266774463d0*x131
      x133 = 0.558025705063192d0*x3 - 0.352521395009435d0*x7
      x134 = exp(-0.49698248213959023d0*x7)
      x135 = x133*x134*x44
      x136 = x0*x135
      x137 = 10.9027235569928d0*x136
      x138 = 0.817704266774463d0*x97
      x139 = x138*x94
      x140 = x109*x138
      x141 = x1**1.0d0*x118
      x142 = 1.55214407110551d0*x48
      x143 = 4.0d0*mu**4
      x144 = -x111*(0.00111534778586481d0*x110 + x77) + x113*x117 - x129
     &*( 0.0167302167879722d0*x122 - 0.526297323276602d0*x128) + x143*( 
     &x113*x142 + x130 - 0.00139418473233101d0*x141*(-x132 + x137 + x139
     & + x140)) - x56*x59 - x61*x74
      x145 = x3*x62
      x146 = 0.508616435555896d0*x145 + 1.0d0
      x147 = 0.25d0/x146**4
      x148 = x144*x147
      x149 = x43*x54
      x150 = x41*x42
      x151 = 0.17845564573837d0*x150*x58
      x152 = x44*x45
      x153 = x152*(0.045088146776381d0*x3 + 0.00168899541983155d0*x48 + 
     &x52 + 0.00627239743302441d0*x7)
      x154 = x20**(-0.333333333333333d0)
      x155 = x21**(-0.333333333333333d0)
      x156 = x154 - x155
      x157 = 2.0d0*x156*x63*x74
      x158 = x68*x7
      x159 = x156*x19
      x160 = mu*x159*x64*x7
      x161 = 7.42981738143299d0*x8
      x162 = x159*x161*x67 - 1.85745434535825d0*x8
      x163 = x4*x65
      x164 = mu*x17/x72
      x165 = x32*(0.344044309698576d0*x3 + 1.45760735710959d0*x4 + 1.876
     &99471636596d0 *x7 + 2.24298561906574d0*x8)/(x33**2*x34)
      x166 = 0.0691804894611514d0*x11*x6
      x167 = 0.062230323235858d0*x166*x35
      x168 = x15*(0.27127533634502d0*x3 + 0.709361408239834d0*x4 + 1.086
     &51697314076d0 *x7 + 1.20668365571947d0*x8)/(x9**2*(x13 + 2.0d0))
      x169 = 0.129439556896076d0*x14*x166
      x170 = -x168 + x169
      x171 = -x20**0.333333333333333d0 + x21**0.333333333333333d0
      x172 = 2.56488140124205d0*x17*x171
      x173 = rho_s**5*x172/rho_c**5
      x174 = 0.0790071588102608d0*x29
      x175 = x22*x24
      x176 = 0.0506609980493537d0*x171*x30
      x177 = x170 + x174*x38 + 0.58482233974552d0*x175*x25*(0.2733500472
     &99742d0*x3 + 0.381163760967645d0*x4 + 1.09726826998169d0*x7 + 1.64
     &535495376155d0*x8)/(x26**2*x27) - 0.000665486074881265d0*x175 *x28
     &*x7 + x176*x19
      x178 = -x173*x37 + x177 - x38*(x165 - x167 + x170)
      x179 = 0.350864882184401d0*x75
      x180 = 24.0d0*x78*x79*x81
      x181 = x19*x85
      x182 = 1d0/x6
      x183 = 0.00706864485168613d0*x182*x49*x96
      x184 = x183/x88
      x185 = x184*x92*(-x181 + 1.0d0)
      x186 = x101*x19
      x187 = x183/x103
      x188 = x107*x187*(x186 + 1.0d0)
      x189 = 0.545136177849642d0*x98
      x190 = x189*x94
      x191 = x20*x91
      x192 = 2.18054471139857d0*x98
      x193 = x191*x192*x93
      x194 = x19*x193
      x195 = x109*x189
      x196 = x106*x21
      x197 = x108*x192*x196
      x198 = x19*x197
      x199 = 0.034943218589452d0*x3
      x200 = x199*x86
      x201 = 0.164800693786261d0*x89
      x202 = x200 + x201
      x203 = x93*x99/x90**2
      x204 = x203*x84*(-x181*x200 - x181*x201 + x202)
      x205 = x102*x199
      x206 = 0.164800693786261d0*x104
      x207 = x205 + x206
      x208 = x108*x99/x105**2
      x209 = x100*x208*(x186*x205 + x186*x206 + x207)
      x210 = 0.00111534778586481d0*x58
      x211 = x120*x43
      x212 = 0.0334604335759443d0*x119
      x213 = x152*(-0.0225440733881905d0*x3 - 0.000844497709915773d0*x48
     & + 0.0047892375d0*x51 + 0.113935126864452d0*x53*x7 - 0.00313619871
     &65122d0*x7)
      x214 = x125*x173 - 4.0d0*x126 + x177 + x38*(-x165 + x167 + x168 - 
     &x169)
      x215 = 0.526297323276602d0*x123*x124
      x216 = 0.263148661638301d0*x51
      x217 = x21**1.66666666666667d0
      x218 = x20**1.66666666666667d0
      x219 = 2.18054471139857d0*x98
      x220 = x0*x134*x95
      x221 = 21.8054471139857d0*x133*x220
      x222 = 0.00139418473233101d0*x119
      x223 = x19*x54
      x224 = rho_s**3*x22/rho_c**3
      x225 = -4.0d0*x17*x224*x37 + x172*x23*x37 - x174*x224 - x176
      x226 = x120*x19
      x227 = x184*x191 - x187*x196 + x193 - x197 + x20*x202*x203 - x207*
     &x208*x21
      E = -rho_c*(x148 + x18 - x31 + x40)
      d1E(1) = x112 - 0.339077623703931d0*x144*x145/x146**5 - x147*( -0.
     &00116228665296198d0*x1**2.0d0*x56 + x111*( 0.00185891297644135d0*x
     &110 - x149*x179 + x153*x76 - x210*( 2.41662179562687d0*x1**(-0.666
     &666666666667d0)*x80 + x180*x43 + x185 + x188 + x190 - x194 + x195 
     &+ x198 + x204 + x209 + x44*x79* x82*(-0.550321208149104d0*x3 + 1.4
     &9912448908652d0*x7) + x83) - 0.00490180588786583d0*x47*x55 + x77) 
     &- 0.713822582953479d0*x113* x114*x115*x116 + x117*(x178 + 4.0d0*x4
     &0) - x129*( -0.0278836946466203d0*x122 + 1.0525946465532d0*x128 + 
     &x211*x212 + x212*x213 - x214*x215) - x143*(2.06952542814068d0*x127
     &*x48 + x130 + 0.00232364122055169d0*x141*(x132 - x137 - x139 - x14
     &0) - x142* x214 + x211*x216 + x213*x216 - x222*(0.545136177849642d
     &0*x131*x95 - 1.80615420514536d0*x135*x182*x5 - 7.26848237132856d0*
     &x136*x95 - x185 - x188 + x19*x219*(x217 - x218) - x190 + x194 - x1
     &95 - x198 - x204 - x209 - 10.9027235569928d0*x220*x44*(-0.37201713
     &6708795d0 *x3 + 0.117507131669812d0*x7) - x221*x43)) - x149*x151 +
     & x150* x153*x59 + x157*x19 - x164*x63*(-9.02445731974005d0*x158 + 
     &54.281640924821d0*x159*x62*x71 + 36.0978292789602d0*x160 + x162 - 
     &13.5704102312052d0*x163 - x73*(-4.17757235705837d0*x158 + 16.71028
     &94282335d0*x160 + x162)) + 0.148713038115308d0*x56*x58) - x148 + x
     &178 + 3.0d0*x40
      d1E(2) = -x147*(x111*(x179*x223 - x210*(-x180*x19 + x227)) + x117*
     &x225 + x129*( x212*x226 + x215*x225) + x143*(x142*x225 + x216*x226
     & - x222*(-x19 *x221 - x219*(-x217 + x218) + x227)) + x151*x223 - x
     &157 + x164* x60*(x156*(36.0978292789602d0*x158 + x161 + 54.2816409
     &24821d0* x163) + x73*(-x154 + x155)*(16.7102894282335d0*x158 + x16
     &1))) + x225
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E, d2E
     &)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
      x0 = 2.14502939711103d0
      x1 = 1d0/rho_c
      x2 = x1**0.666666666666667d0
      x3 = x2/x0
      x4 = 0.564189583547756d0*sqrt(x1)
      x5 = 1.41872281647967d0*x4
      x6 = 1.46459188756152d0
      x7 = x1**0.333333333333333d0
      x8 = x7/x6
      x9 = 0.826307487110758d0*x1**0.166666666666667d0
      x10 = 0.406913004517529d0*x3 + x5 + 3.25955091942229d0*x8 + 7.2401
     &0193431683d0 *x9
      x11 = 9.86960440108936d0
      x12 = 0.306852819440055d0
      x13 = x11/x12
      x14 = x13/x10
      x15 = log(0.5d0*x14 + 1.0d0)
      x16 = 0.194159335344114d0*x8 + 1.0d0
      x17 = 2.0d0*x16
      x18 = x15*x17
      x19 = x12/x11
      x20 = x18*x19
      x21 = rho_s*x1
      x22 = x21 + 1.0d0
      x23 = -x21 + 1.0d0
      x24 = 1.92366105093154d0*x22**1.33333333333333d0 + 1.9236610509315
     &4d0*x23** 1.33333333333333d0 - 3.84732210186307d0
      x25 = rho_s**4/rho_c**4
      x26 = -x25 + 1.0d0
      x27 = x24*x26
      x28 = 0.101077332976288d0*x8 + 1.0d0
      x29 = 0.76232752193529d0*x4
      x30 = x29 + 0.410025070949612d0*x3 + 3.29180480994506d0*x8 + 9.872
     &12972256927d0*x9
      x31 = 1.0d0 + 29.6085746432167d0/x30
      x32 = log(x31)
      x33 = x28*x32
      x34 = 0.0197517897025652d0*x33
      x35 = x27*x34
      x36 = 0.186690969707574d0*x8 + 1.0d0
      x37 = 2.91521471421918d0*x4
      x38 = 0.516066464547863d0*x3 + x37 + 5.63098414909788d0*x8 + 13.45
     &79137143944d0*x9
      x39 = x13/x38 + 1.0d0
      x40 = log(x39)
      x41 = x36*x40
      x42 = -x18 + x41
      x43 = x24*x25
      x44 = x19*x43
      x45 = x42*x44
      x46 = x1**1.66666666666667d0
      x47 = 0.0837224526465878d0*x46
      x48 = rho_s**2/rho_c**2
      x49 = -x48 + 1.0d0
      x50 = exp(-0.6836107611867116d0*x8)
      x51 = 4.60115111447049d0
      x52 = x1**1.33333333333333d0
      x53 = x52/x51
      x54 = 3.14159265358979d0
      x55 = x1**1.0d0
      x56 = x55/x54
      x57 = -0.009578475d0*x56
      x58 = 0.0676322201645715d0*x3 + 0.00126674656487366d0*x53 + x57 + 
     &0.0188171922990732d0*x8 + 1.0d0
      x59 = x50*x58
      x60 = x49*x59
      x61 = x47*x60
      x62 = 1.41421356237310d0
      x63 = mu**5*x62
      x64 = x61*x63
      x65 = x22**0.666666666666667d0
      x66 = x23**0.666666666666667d0
      x67 = x65 + x66
      x68 = x67**3
      x69 = mu**2
      x70 = x67**2
      x71 = 1d0/x70
      x72 = x69*x71
      x73 = x72*x8
      x74 = 1d0/x67
      x75 = mu*x74
      x76 = 11.1447260721495d0*x75*x9 + 1.0d0
      x77 = 12.5327170711751d0*x73 + x76
      x78 = 1d0/x77
      x79 = mu**3
      x80 = 1d0/x68
      x81 = x4*x80
      x82 = 27.0733719592201d0*x73 + x76 + 27.1408204624105d0*x79*x81
      x83 = x78*x82
      x84 = x19*log(x83)
      x85 = 0.179587122125167d0*x55
      x86 = x59*x85
      x87 = 0.175432441092201d0*x49
      x88 = x86*x87
      x89 = exp(-0.28165369188898165d0*x8)
      x90 = x89*(0.825481812223657d0*x3 - 4.49737346725955d0*x8)
      x91 = x54/x55
      x92 = x90*x91
      x93 = 12.0d0*x49
      x94 = x92*x93
      x95 = x22**2
      x96 = 1d0/x22
      x97 = x96**0.333333333333333d0
      x98 = x8*x97
      x99 = 0.0259335011650457d0*x98
      x100 = -x99 + 1.0d0
      x101 = x96**0.666666666666667d0
      x102 = 0.0524148278841779d0*x3
      x103 = x101*x102 + 0.494402081358784d0*x98 + 1.0d0
      x104 = 1d0/x103
      x105 = 1d0/x101
      x106 = x104*x105
      x107 = x100*x106
      x108 = x107*x95
      x109 = 1d0/x2
      x110 = 2.28942848510666d0
      x111 = x110*x51
      x112 = x109*x111
      x113 = 0.817704266774463d0*x112
      x114 = x23**2
      x115 = 1d0/x23
      x116 = x115**0.333333333333333d0
      x117 = x116*x8
      x118 = 0.0259335011650457d0*x117
      x119 = -x118 + 1.0d0
      x120 = x115**0.666666666666667d0
      x121 = x102*x120 + 0.494402081358784d0*x117 + 1.0d0
      x122 = 1d0/x121
      x123 = 1d0/x120
      x124 = x122*x123
      x125 = x119*x124
      x126 = x114*x125
      x127 = x47*(x108*x113 + x113*x126 + x94)
      x128 = x62*x79
      x129 = 4.0d0*x128
      x130 = -x20 + x35
      x131 = x130 - x45
      x132 = 0.0472353356922751d0*mu**8*x1**2.66666666666667d0
      x133 = 0.267683468607554d0*x132
      x134 = x59 - 1.0d0
      x135 = 0.148394183600699d0
      x136 = x135*x46
      x137 = x134*x136
      x138 = x137*x49
      x139 = x18 - x41
      x140 = x139*x44
      x141 = x130 + x140
      x142 = 0.101321183642338d0*x1**2.0d0
      x143 = x141*x142
      x144 = mu**6
      x145 = 4.0d0*x144
      x146 = x134*x56
      x147 = -0.13157433081915d0*x146*x49
      x148 = x111*(x22**2.66666666666667d0 + x23**2.66666666666667d0)
      x149 = 0.817704266774463d0*x148
      x150 = 0.558025705063192d0*x3
      x151 = 0.352521395009435d0*x8
      x152 = exp(-0.49698248213959023d0*x8)
      x153 = x152*(x150 - x151)
      x154 = x153*x49
      x155 = x0*x154
      x156 = 10.9027235569928d0*x155
      x157 = 0.817704266774463d0*x111
      x158 = x108*x157
      x159 = x126*x157
      x160 = x1**1.0d0*x135
      x161 = x160*(-x149 + x156 + x158 + x159)
      x162 = 1.55214407110551d0*x53
      x163 = mu**4
      x164 = 4.0d0*x163
      x165 = -x129*(0.00111534778586481d0*x127 + x88) + x131*x133 - x145
     &*( 0.0167302167879722d0*x138 - 0.526297323276602d0*x143) + x164*( 
     &x131*x162 + x147 - 0.00139418473233101d0*x161) - 0.089227822869184
     &8d0*x64 - x68*x84
      x166 = x3*x69
      x167 = 0.508616435555896d0*x166 + 1.0d0
      x168 = x167**(-4)
      x169 = 0.25d0*x168
      x170 = x165*x169
      x171 = x48*x59
      x172 = x47*x63
      x173 = 0.17845564573837d0*x172
      x174 = 0.045088146776381d0*x3 + 0.00168899541983155d0*x53 + x57 + 
     &0.00627239743302441d0*x8
      x175 = x49*x50
      x176 = x174*x175
      x177 = x172*x176
      x178 = x60*x63
      x179 = 0.0571643564037363d0*x1**2.0d0
      x180 = x178*x179
      x181 = x22**(-0.333333333333333d0)
      x182 = x23**(-0.333333333333333d0)
      x183 = x181 - x182
      x184 = 2.0d0*x70*x84
      x185 = x183*x184
      x186 = x131*x132
      x187 = x75*x8
      x188 = mu*x8
      x189 = x188*x71
      x190 = x183*x21
      x191 = x189*x190
      x192 = 7.42981738143299d0*x9
      x193 = x21*x74
      x194 = x183*x193
      x195 = x192*x194 - 1.85745434535825d0*x9
      x196 = -4.17757235705837d0*x187 + 16.7102894282335d0*x191 + x195
      x197 = x4*x72
      x198 = x69*x81
      x199 = x190*x198
      x200 = -9.02445731974005d0*x187 + 36.0978292789602d0*x191 + x195 -
     & 13.5704102312052d0*x197 + 54.281640924821d0*x199
      x201 = x19*(-x196*x83 + x200)
      x202 = 1d0/x82
      x203 = mu*x202
      x204 = x203*x70
      x205 = 0.344044309698576d0*x3 + 1.45760735710959d0*x4 + 1.87699471
     &636596d0*x8 + 2.24298561906574d0*x9
      x206 = 1d0/x39
      x207 = x38**(-2)
      x208 = x206*x207*x36
      x209 = x205*x208
      x210 = 0.062230323235858d0*x7
      x211 = 0.0691804894611514d0*x12
      x212 = x211*x40
      x213 = x210*x212
      x214 = 0.27127533634502d0*x3 + 0.709361408239834d0*x4 + 1.08651697
     &314076d0*x8 + 1.20668365571947d0*x9
      x215 = x14 + 2.0d0
      x216 = 1d0/x215
      x217 = x10**(-2)
      x218 = x17*x216*x217
      x219 = x214*x218
      x220 = 0.129439556896076d0*x7
      x221 = x15*x211
      x222 = x220*x221
      x223 = -x219 + x222
      x224 = x209 - x213 + x223
      x225 = x224*x43
      x226 = x19*x42
      x227 = x23**0.333333333333333d0
      x228 = x22**0.333333333333333d0
      x229 = x227 - x228
      x230 = rho_s**5/rho_c**5
      x231 = x229*x230
      x232 = 2.56488140124205d0*x231
      x233 = 0.0790071588102608d0*x33
      x234 = x32*x8
      x235 = x234*x27
      x236 = 0.273350047299742d0*x3 + 0.381163760967645d0*x4 + 1.0972682
     &6998169d0*x8 + 1.64535495376155d0*x9
      x237 = 1d0/x31
      x238 = x30**(-2)
      x239 = x27*x28
      x240 = x237*x238*x239
      x241 = 0.58482233974552d0*x240
      x242 = x26*x33
      x243 = x229*x242
      x244 = 0.0506609980493537d0*x243
      x245 = x21*x244 + x223 + x233*x43 - 0.000665486074881265d0*x235 + 
     &x236*x241
      x246 = -x225 - x226*x232 + x245
      x247 = x246 + 4.0d0*x45
      x248 = 0.350864882184401d0*x86
      x249 = x174*x85
      x250 = x50*x87
      x251 = 0.122619224952946d0*x52
      x252 = x251*x60
      x253 = x89*(-0.550321208149104d0*x3 + 1.49912448908652d0*x8)
      x254 = x91*x93
      x255 = x49*x90
      x256 = 2.14502939711103d0*x1**(-0.666666666666667d0)
      x257 = x255*x256
      x258 = 24.0d0*x92
      x259 = x258*x48
      x260 = 1d0/x97
      x261 = x104*x260
      x262 = x261*x95
      x263 = x21*x96
      x264 = -x263 + 1.0d0
      x265 = 1d0/x7
      x266 = x110*x265*x54
      x267 = 0.00706864485168613d0*x266
      x268 = x264*x267
      x269 = x262*x268
      x270 = 1d0/x116
      x271 = x122*x270
      x272 = x114*x271
      x273 = x115*x21
      x274 = x273 + 1.0d0
      x275 = x267*x274
      x276 = x272*x275
      x277 = 0.545136177849642d0*x112
      x278 = x108*x277
      x279 = x126*x277
      x280 = x103**(-2)
      x281 = x105*x280
      x282 = x100*x281
      x283 = x113*x282
      x284 = 0.034943218589452d0*x3
      x285 = x101*x284
      x286 = 0.164800693786261d0*x98
      x287 = x285 + x286
      x288 = -x263*x285 - x263*x286 + x287
      x289 = x288*x95
      x290 = x283*x289
      x291 = x121**(-2)
      x292 = x123*x291
      x293 = x119*x292
      x294 = x113*x293
      x295 = x120*x284
      x296 = 0.164800693786261d0*x117
      x297 = x295 + x296
      x298 = x273*x295 + x273*x296 + x297
      x299 = x114*x298
      x300 = x294*x299
      x301 = 2.18054471139857d0*x112
      x302 = x125*x23
      x303 = x301*x302
      x304 = x21*x303
      x305 = x107*x22
      x306 = x301*x305
      x307 = x21*x306
      x308 = x269 + x276 + x278 + x279 + x290 + x300 + x304 - x307
      x309 = x253*x254 + 1.12661476755593d0*x257 + x259 + x308 + x94
      x310 = 0.00111534778586481d0*x47
      x311 = 0.0334604335759443d0*x137
      x312 = x58*x8
      x313 = 0.000844497709915773d0*x53
      x314 = 0.0047892375d0*x56
      x315 = 0.0225440733881905d0*x3
      x316 = 0.0031361987165122d0*x8
      x317 = -x313 + x314 - x315 - x316
      x318 = 0.113935126864452d0*x312 + x317
      x319 = x175*x318
      x320 = 0.0334604335759443d0*x136
      x321 = -x209 + x213 + x219 - x222
      x322 = x321*x43
      x323 = x139*x19
      x324 = -4.0d0*x140 + x232*x323 + x245 + x322
      x325 = 0.526297323276602d0*x142
      x326 = 0.263148661638301d0*x146
      x327 = x319*x56
      x328 = 2.06952542814068d0*x53
      x329 = x109*x148
      x330 = 0.545136177849642d0*x329
      x331 = x23**1.66666666666667d0
      x332 = x22**1.66666666666667d0
      x333 = x331 - x332
      x334 = 2.18054471139857d0*x112
      x335 = x21*x333*x334
      x336 = x0*x109
      x337 = x153*x336
      x338 = 21.8054471139857d0*x337
      x339 = x338*x48
      x340 = x265*x6
      x341 = x154*x340
      x342 = 1.80615420514536d0*x341
      x343 = x109*x155
      x344 = 7.26848237132856d0*x343
      x345 = x152*x336
      x346 = -0.372017136708795d0*x3 + 0.117507131669812d0*x8
      x347 = x346*x49
      x348 = x345*x347
      x349 = 10.9027235569928d0*x348
      x350 = 0.00139418473233101d0*x136
      x351 = x129*(0.00185891297644135d0*x127 - x248*x48 + x249*x250 - 0
     &.0399758348639607d0*x252 - x309*x310 + x88) + x133*x247 - x145*( -
     &0.0278836946466203d0*x138 + 1.0525946465532d0*x143 + x311*x48 + x3
     &19*x320 - x324*x325) - x164*(x141*x328 + x147 + 0.0023236412205516
     &9d0*x160*(x149 - x156 - x158 - x159) - x162* x324 + x326*x48 + 0.2
     &63148661638301d0*x327 - x350*(-x269 - x276 - x278 - x279 - x290 - 
     &x300 - x304 + x307 + x330 + x335 - x339 - x342 - x344 - x349)) - x
     &171*x173 + 0.0892278228691848d0*x177 - 0.0203323666368788d0*x180 +
     & x185*x21 - 0.713822582953479d0*x186 - x201*x204 + 0.1487130381153
     &08d0*x64
      x352 = x166/x167**5
      x353 = x165*x352
      x354 = x21*x59
      x355 = x203*x67
      x356 = 36.0978292789602d0*x187 + x192 + 54.281640924821d0*x197
      x357 = x183*x356
      x358 = -x181 + x182
      x359 = 16.7102894282335d0*x187 + x192
      x360 = x358*x359
      x361 = x19*(x357 + x360*x83)
      x362 = rho_s**3/rho_c**3
      x363 = x24*x362
      x364 = -x244
      x365 = x229*x25
      x366 = 2.56488140124205d0*x365
      x367 = 4.0d0*x363
      x368 = x226*x366 - x226*x367 - x233*x363 + x364
      x369 = x21*x258
      x370 = x22*x261
      x371 = x23*x271
      x372 = x282*x287
      x373 = x22*x372
      x374 = x293*x297
      x375 = x23*x374
      x376 = x113*x373 - x113*x375 + x267*x370 - x267*x371 - x303 + x306
     &
      x377 = -x369 + x376
      x378 = -x21*x338 - x334*(-x331 + x332) + x376
      x379 = x129*(x21*x248 - x310*x377) + x133*x368 + x145*(x21*x311 + 
     &x325*x368) + x164*(x162*x368 + x21*x326 - x350*x378) + x173*x354 -
     & x185 + x355 *x361
      x380 = -x169*x379
      x381 = x48*x50
      x382 = x179*x63
      x383 = 0.121994199821273d0*x382
      x384 = 0.0751469112939684d0*x3 + 0.00394098931294027d0*x53 - 0.019
     &15695d0*x56 + 0.00836319657736588d0*x8
      x385 = x183**2
      x386 = 8.0d0*x385*x67*x84
      x387 = x22**(-1.33333333333333d0)
      x388 = x21*x387
      x389 = x23**(-1.33333333333333d0)
      x390 = x21*x389
      x391 = x388 + x390
      x392 = -6.0d0*x181 + 6.0d0*x182 + x391
      x393 = x21*x392
      x394 = 12.0d0*x183
      x395 = x201*x355
      x396 = 3.0d0*x69
      x397 = x201*x396
      x398 = x202*x78
      x399 = x397*x398
      x400 = x82**(-2)
      x401 = x397*x400
      x402 = 5.57009647607783d0*x187
      x403 = x385*x48
      x404 = x188*x80
      x405 = x189*x393
      x406 = 2.476605793811d0*x9
      x407 = 9.90642317524398d0*x9
      x408 = x193*x392*x406 - x194*x406 + x403*x407*x71 + 2.167030069584
     &62d0*x9
      x409 = 12.0326097596534d0*x187
      x410 = x196*x75
      x411 = 2.0d0*x78
      x412 = 14.859634762866d0*x9
      x413 = 33.420578856467d0*x189
      x414 = x190*x413
      x415 = x82/x77**2
      x416 = x410*x415
      x417 = 144.751042466189d0*x385
      x418 = 72.1956585579204d0*x385
      x419 = 3.0d0*x19
      x420 = 20.5190512099364d0*x231
      x421 = x33*x43
      x422 = 0.455583412166236d0*x3 + 0.571745641451468d0*x4 + 1.4630243
     &5997558d0*x8 + 1.91958077938847d0*x9
      x423 = x236*x237
      x424 = x239*x423*(x29 + 0.546700094599483d0*x3 + 2.19453653996337d
     &0*x8 + 3.29070990752309d0*x9)/x30**3
      x425 = x238*x423
      x426 = x28*x425
      x427 = x426*x43
      x428 = x229*x26
      x429 = x234*x428
      x430 = x21*x429
      x431 = x426*x428
      x432 = x21*x431
      x433 = x221*x7
      x434 = 0.172586075861435d0*x433
      x435 = x218*(0.452125560575033d0*x3 + 1.06404211235975d0*x4 + 1.44
     &868929752102d0*x8 + 1.40779759833938d0*x9)
      x436 = x214*x216
      x437 = x17*x436*(0.542550672690039d0*x3 + x5 + 2.17303394628153d0*
     &x8 + 2.41336731143894d0*x9)/x10**3
      x438 = 0.682784063255296d0
      x439 = x217*x436
      x440 = x220*x438*x439
      x441 = 0.129439556896076d0*x439*x8
      x442 = x13*x17*x214**2/(x10**4*x215**2)
      x443 = x435 - x437 + x440 + x441 + x442
      x444 = -x434 + x443
      x445 = 0.40528798439483d0*x33
      x446 = x22**(-0.666666666666667d0)
      x447 = 0.854960467080683d0*x21
      x448 = x23**(-0.666666666666667d0)
      x449 = x446*x447 + x447*x448
      x450 = -5.1297628024841d0*x227 + 5.1297628024841d0*x228 + x449
      x451 = x26*x34
      x452 = x21*x450*x451 + x231*x445 - 0.00532388859905012d0*x234*x43 
     &- 17.3157558993759d0*x236**2*x239/(x30**4*x31**2) - 0.039408188244
     &2864d0*x27*x425*x8
      x453 = 0.000887314766508353d0*x235 - x241*x422 - 0.395035794051304
     &d0*x421 + 0.58482233974552d0*x424 + 4.67857871796416d0*x427 - 0.00
     &341378571249706d0*x430 + 2.99999988448829d0*x432 + x444 + x452
      x454 = x208*(0.573407182830959d0*x3 + 2.18641103566438d0*x4 + 2.50
     &265962182128d0*x8 + 2.6168165555767d0*x9)
      x455 = 0.0829737643144773d0*x212*x7
      x456 = x205*x206
      x457 = x36*x456*(0.688088619397151d0*x3 + x37 + 3.75398943273192d0
     &*x8 + 4.48597123813148d0*x9)/x38**3
      x458 = x207*x456
      x459 = x210*x438*x458
      x460 = 0.062230323235858d0*x458*x8
      x461 = x13*x205**2*x36/(x38**4*x39**2)
      x462 = 5.1297628024841d0*x231
      x463 = x230*x450
      x464 = -x224*x462 - x226*x463 - x43*(x444 - x454 + x455 + x457 - x
     &459 - x460 - x461)
      x465 = 8.0d0*x225 + x226*x420 - 20.0d0*x45 + x453 + x464
      x466 = 0.803050405822664d0*x132
      x467 = 0.0799516697279215d0*x251
      x468 = 0.350864882184401d0*x249
      x469 = x253*x49
      x470 = 24.0d0*x91
      x471 = 2.25322953511185d0*x256
      x472 = 0.00942485980224818d0*x266
      x473 = x264*x472
      x474 = x274*x472
      x475 = 1d0/x95
      x476 = 0.0115260005177981d0*x48
      x477 = 0.817704266774463d0*x266
      x478 = 0.181712059283214d0*x112
      x479 = 1d0/x114
      x480 = 0.0141372897033723d0*x266
      x481 = x260*x280
      x482 = x288*x481
      x483 = x270*x291
      x484 = x480*x483
      x485 = 3.63424118566428d0*x112
      x486 = x107*x485
      x487 = x125*x485
      x488 = 1.09027235569928d0*x112
      x489 = 0.116477395298173d0*x3
      x490 = 0.0582386976490866d0*x3
      x491 = x101*x490
      x492 = x475*x48
      x493 = 0.219734258381682d0*x98
      x494 = -x491 - x493
      x495 = x120*x490
      x496 = x479*x48
      x497 = 0.219734258381682d0*x117
      x498 = x495 + x497
      x499 = 0.0376994392089927d0*x21*x266
      x500 = 1.45369647426571d0*x112*x21
      x501 = 0.0698864371789039d0*x3
      x502 = x101*x501
      x503 = 0.329601387572523d0*x98
      x504 = x502 + x503
      x505 = x103**(-3)
      x506 = x100*x105*x113*x505
      x507 = x120*x501
      x508 = 0.329601387572523d0*x117
      x509 = x507 + x508
      x510 = x121**(-3)
      x511 = x113*x119*x123*x510
      x512 = x264*x370
      x513 = x23*x293
      x514 = 4.36108942279713d0*x112
      x515 = x21*x514
      x516 = x22*x282
      x517 = -x108*x478 - x114*x294*(0.439468516763364d0*x117*x273 + x12
     &0*x273*x489 + x495*x496 + x496*x497 + x498) - x126*x478 + x262*x47
     &3 + x262*x477 *(0.0230520010355962d0*x263 - x475*x476 - 0.01152600
     &05177981d0) + x264*x480*x482*x95 + x272*x474 - x272*x477*(0.023052
     &0010355962d0* x273 + x476*x479 + 0.0115260005177981d0) + x274*x299
     &*x484 + x274* x371*x499 + x282*x289*x488 + x283*x95*(x101*x263*x48
     &9 + 0.439468516763364d0*x263*x98 - x491*x492 - x492*x493 + x494) -
     & x288*x515*x516 + x289*x506*(-x263*x502 - x263*x503 + x504) + x293
     & *x299*x488 + x298*x513*x515 + x299*x511*(x273*x507 + x273*x508 + 
     &x509) - x302*x500 + x305*x500 + x48*x486 + x48*x487 - x499*x512
      x518 = 12.0d0*x128
      x519 = x318*x381
      x520 = 0.455740507457808d0*x8
      x521 = 0.0259624262672375d0*x3*x58 + 0.0375734556469842d0*x3 - 0.1
     &51913502485936d0*x312 + 0.00197049465647014d0*x53 + x57 + 0.004181
     &59828868294d0*x8
      x522 = 12.0d0*x144
      x523 = 3.61230841029072d0*x152*x340
      x524 = x345*x346
      x525 = 5.0d0*x21
      x526 = x525*x65 + x525*x66
      x527 = 0.726848237132856d0*x112
      x528 = x112*x333
      x529 = 12.0d0*x163
      x530 = 0.0833333333333333d0*x168
      x531 = x21*x50
      x532 = 0.535366937215109d0*x172
      x533 = 3.0d0*x181
      x534 = 3.0d0*x182
      x535 = x391 - x533 + x534
      x536 = x203*x361
      x537 = x356*x358
      x538 = x183*x359
      x539 = x33*x363
      x540 = x226*x363
      x541 = x226*x365
      x542 = x363*x426
      x543 = -2.56488140124205d0*x227 + 2.56488140124205d0*x228 + x449
      x544 = x226*x25
      x545 = x224*x366 - x224*x367 + 0.00266194429952506d0*x234*x363 - x
     &365*x445 - x451*x543 + x543*x544
      x546 = 0.00170689285624853d0*x429 - 1.49999994224414d0*x431 + 0.31
     &6028635241043d0*x539 + 16.0d0*x540 - 20.5190512099364d0*x541 - 2.3
     &3928935898208d0*x542 + x545
      x547 = 1.2383028969055d0*x9
      x548 = x407*x74
      x549 = x385*x548
      x550 = x189*x418
      x551 = x198*x417
      x552 = 12.0326097596534d0*x187
      x553 = 18.0938803082737d0*x197
      x554 = 33.420578856467d0*x187 + x412
      x555 = -x388 - x390 + x533 - x534
      x556 = 5.57009647607783d0*x187
      x557 = x355*x419
      x558 = x318*x531
      x559 = 0.0188497196044964d0*x266
      x560 = x21*x559
      x561 = x261*x560
      x562 = 0.00471242990112409d0*x266
      x563 = x370*x562
      x564 = x271*x560
      x565 = x371*x562
      x566 = x371*x474
      x567 = 0.726848237132856d0*x112
      x568 = x287*x481
      x569 = x22*x268*x568
      x570 = x21*x301
      x571 = x23*x483
      x572 = x275*x297*x571
      x573 = x22*x267*x482
      x574 = x267*x298*x571
      x575 = x263*x491
      x576 = x263*x493
      x577 = x288*x301
      x578 = x298*x301
      x579 = x273*x495 + x273*x497 + x498
      x580 = x504*x506
      x581 = x509*x511
      x582 = x21*(-x150 + x151)
      x583 = x99 - 1.0d0
      x584 = x21*x485
      x585 = x22*x583
      x586 = x118 - 1.0d0
      x587 = x23*x586
      x588 = x281*x585
      x589 = x292*x587
      x590 = x491 + x493
      x591 = x229*x362
      x592 = x445*x591
      x593 = x446 + x448
      x594 = 0.0168869993497846d0*x242*x593
      x595 = x24*x48
      x596 = x33*x595
      x597 = x544*x593
      x598 = 20.5190512099364d0*x226*x591
      x599 = 12.0d0*x226*x595
      x600 = x387 + x389
      x601 = x361*x396*x74
      x602 = x592 + x594 - 0.237021476430782d0*x596 - 0.854960467080683d
     &0*x597 + x598 - x599
      x603 = x406*x600
      x604 = mu*x71
      x605 = x358**2
      x606 = 0.0282745794067445d0*x266
      x607 = x261*x606 + x271*x606 - x283*x590 + x287*x580 - x294*x498 +
     & x297*x484 + x297*x581 + x372*x514 + x374*x514 + x480*x568 + x486 
     &+ x487
      E = -rho_c*(x170 + x20 - x35 + x45)
      d1E(1) = x130 - x169*x351 - x170 + x246 - 0.339077623703931d0*x353
     & + 3.0d0*x45
      d1E(2) = x368 + x380
      d2E(1) = x1*(-4.0d0*x16*x439 - 0.574868174483523d0*x163*x165*x53/x
     &167**6 - 0.5d0* x168*x351 + 0.101321996098707d0*x21*x243 + 6.0d0*x
     &225 + 15.3892884074523d0*x226*x231 - 0.000443657383254177d0*x235 +
     & 1.16964467949104d0*x236*x240 - 0.58482233974552d0*x240*x422 - 0.6
     &78155247407862d0*x351*x352 - 0.113025874567977d0*x353 - 0.23702147
     &6430782d0*x421 + 0.58482233974552d0*x424 + 4.67857871796416d0*x427
     & - 0.00341378571249706d0*x430 + 2.99999988448829d0*x432 + 0.086293
     &0379307174d0*x433 + x443 - 12.0d0*x45 + x452 + x464 - x530*(-0.000
     &542507213303895d0*x1** 2.33333333333333d0*x178 - 4.28293549772087d
     &0*x132*x247 + 3.39065726902902d0*x171*x172 - 0.243988399642546d0*x
     &171*x382 + 1.07073387443022d0*x172*x174*x381 - 0.267683468607554d0
     &*x172*x175 *x384 + x176*x383 - 0.892278228691848d0*x177 + 0.284653
     &132916304d0*x180 + x184*x393 + 7.85204841248826d0*x186 - x196*x399
     &*x67 + x200*x401*x67 - x204*x419*(-24.0652195193068d0* x191 + 20.3
     &556153468079d0*x197 + 18.0938803082737d0*x198*x393 - 54.2816409248
     &21d0*x199 - x200*x410*x411 + x4*x417*x48*x69/x67**4 + x404*x418*x4
     &8 + 12.0326097596534d0*x405 + x408 + x409 + x416*( -8.355144714116
     &75d0*x187 + x194*x412 + x414 - 3.71490869071649d0* x9) - x83*(-11.
     &1401929521557d0*x191 + x402 + 33.420578856467d0* x403*x404 + 5.570
     &09647607783d0*x405 + x408)) + x21*x394*x395 - x386*x48 + x465*x466
     & + x518*(-0.00495710127051027d0*x127 - 0.159903339455843d0*x171*x2
     &51 - x175*x468 + x176*x467 - x248*x49 + 0.701729764368802d0*x249*x
     &381 - x250*x384*x85 + 0.133252782879869d0*x252 + 0.003717825952882
     &7d0*x309*x47 - x310*( 0.154912426780983d0*x1**(-0.333333333333333d
     &0)*x255 + 48.0d0*x253 *x48*x91 + x254*x89*(0.917202013581841d0*x3 
     &- 1.99883265211535d0* x8) + 4.50645907022371d0*x256*x48*x90 + 0.75
     &1076511703951d0*x257 - x259 + x469*x470 + x469*x471 + x517) + 1.75
     &432441092201d0*x48* x86 - 0.00910930363347549d0*x61) - x522*(-0.11
     &1534778586481d0* x136*x319 + 0.133841734303777d0*x136*x519 - 0.211
     &916079314314d0* x137*x48 + 0.074356519057654d0*x138 + 2.1051892931
     &0641d0*x142* x324 - 3.15778393965961d0*x143 + x175*x320*(x317*x520
     & + x521) - x325*(20.0d0*x140 + x321*x462 - 8.0d0*x322 - x323*x420 
     &+ x323* x463 + x43*(x434 - x435 + x437 - x440 - x441 - x442 + x454
     & - x455 - x457 + x459 + x460 + x461) + x453)) + x529*(4.8288926656
     &6159d0* x131*x53 + 0.00464728244110338d0*x136*(x308 - x330 - x335 
     &+ x339 + x342 + x344 + x349) + 1.3157433081915d0*x146*x48 - 0.0061
     &9637658813783d0*x161 + x162*x465 - 0.263148661638301d0*x175 *x56*(
     &-x520*(x313 - x314 + x315 + x316) + x521) - 4.13905085628136d0*x24
     &7*x53 - x326*x49 + 0.526297323276602d0*x327 - x350*(7.224616820581
     &44d0*x153*x340*x48 + 0.299209d0*x154 - x21* x527*(-6.0d0*x331 + 6.
     &0d0*x332 + x526) - 2.90739294853142d0*x21* x528 + 0.18171205928321
     &4d0*x329 - 36.3424118566428d0*x337*x48 - 2.22044604925031d-16*x341
     & - 2.42282745710952d0*x343 + 10.9027235569928d0*x345*x49*(0.620028
     &561181324d0*x3 - 0.156676175559749d0*x8) + x347*x523 + 14.53696474
     &26571d0*x348 + 43.6108942279714d0*x48*x524 + x517) - 1.05259464655
     &32d0*x519*x56 ) - 1.18970430492246d0*x64))
      d2E(2) = x1*(-0.339077623703931d0*x352*x379 + x364 + x380 + 0.0017
     &0689285624853d0 *x429 - 1.49999994224414d0*x431 - x530*(-2.1414677
     &4886044d0*x132* x368 - 1.96301210312207d0*x172*x354 - x174*x531*x5
     &32 - 6.0d0*x183 *x395 - x184*x535 - 6.0d0*x190*x536 + x21*x386 + x
     &354*x383 + x399 *x538 + x401*x537 + x466*x546 + x518*(-1.052594646
     &5532d0*x21*x86 - x310*(-x21*x253*x470 - x21*x471*x90 - x21*x486 - 
     &x21*x487 + x22 *x283*(x494 + x575 + x576) + x22*x288*x580 + x23*x2
     &94*x579 - x23* x298*x581 + x277*x373 - x277*x375 + x302*x567 - x30
     &5*x567 + x369 + x370*x472*(x263 - 1.0d0) - x372*x570 - x374*x570 +
     & x512*x559 - x513*x578 + x516*x577 - x561 + x563 - x564 - x565 - x
     &566 + x569 - x572 + x573 - x574) + x354*x467 + 0.00185891297644135
     &d0*x377*x47 - x468*x531) + x522*(0.0669208671518886d0*x136*x558 - 
     &0.122688256445129d0*x137*x21 - 1.0525946465532d0*x142*x368 + x325 
     &*x546) + x529*(0.00232364122055169d0*x136*x378 - 0.789445984914903
     &d0*x146*x21 + x162*x546 - x328*x368 + x350*(x105 *x113*x288*x504*x
     &505*x585 - x106*x567*x585 - x106*x583*x584 - x113*x123*x298*x509*x
     &510*x587 + x113*x579*x589 - x113*x588*(-x575 - x576 + x590) + x124
     &*x567*x587 - x124*x584*x586 + 21.8054471139857d0*x21*x524 + x277*x
     &287*x588 - x277*x297*x589 - x281*x287*x570*x583 - x292*x297*x570*x
     &586 + 29.0739294853142d0* x345*x582 - x370*x473 - x523*x582 - x527
     &*(-3.0d0*x331 + 3.0d0* x332 + x526) - 1.45369647426571d0*x528 + x5
     &61 - x563 + x564 + x565 + x566 - x569 + x572 - x573 + x574 + x577*
     &x588 - x578*x589) + 0.526297323276602d0*x558*x56) - x557*(27.14082
     &04624105d0*x183* x197 + x183*x409 + x183*x547 - x200*x360*x75*x78 
     &- x21*x549 - x21 *x550 - x21*x551 + x357*x410*x78 + x358*x416*x554
     & - x406*x535 - x535*x552 - x535*x553 - x83*(x190*x358*x548 - x358*
     &x402 + x358* x414 - x358*x547 + x406*x555 + x555*x556))) + 0.23702
     &1476430782d0 *x539 + 12.0d0*x540 - 17.9541698086943d0*x541 - 2.339
     &28935898208d0*x542 + x545)
      d2E(3) = -x1*(x530*(x184*x600 - x386 + x394*x536 - x398*x538*x601 
     &- x400*x537* x601 + x466*x602 + x518*(x248 - x310*(-x258 + x607)) 
     &+ x522*(x311 + x325*x602) + x529*(x162*x602 + x326 - x350*(-x338 -
     & x485*x67 + x607)) + x532*x59 - x557*(x357*x360*x411*x604 + x359*x
     &415*x554* x604*x605 + x549 + x550 + x551 + x552*x600 + x553*x600 +
     & x603 - x83*(x413*x605 + x548*x605 + x556*x600 + x603))) - x592 - 
     &x594 + 0.237021476430782d0*x596 + 0.854960467080682d0*x597 - x598 
     &+ x599 )
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46
      E = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 9.86960440108936d0
      x3 = 1d0/x2
      x4 = x1*x3
      x5 = 1d0/rho_c
      x6 = 0.682784063255296d0*x5**0.333333333333333d0
      x7 = 0.194159335344114d0*x6 + 1.0d0
      x8 = 2.14502939711103d0
      x9 = x5**0.666666666666667d0
      x10 = x9/x8
      x11 = 0.564189583547756d0*sqrt(x5)
      x12 = 0.826307487110758d0*x5**0.166666666666667d0
      x13 = x2/x1
      x14 = log(0.5d0*x13/(0.40691300451752932d0*x10 + 1.418722816479667
     &5d0*x11 + 7.2401019343168302d0*x12 + 3.2595509194222925d0*x6) + 1.
     &0d0)
      x15 = x14*x7
      x16 = 2.0d0*x15*x4
      x17 = mu**2
      x18 = (0.508616435555896d0*x10*x17 + 1.0d0)**(-4)
      x19 = 1.41421356237310d0
      x20 = 4.60115111447049d0
      x21 = x5**1.33333333333333d0
      x22 = 3.14159265358979d0
      x23 = x5**1.0d0
      x24 = x23/x22
      x25 = (0.0676322201645715d0*x10 - 0.009578475d0*x24 + 0.0188171922
     &990732d0*x6 + 1.0d0 + 0.00126674656487366d0*x21/x20)*exp( -0.68361
     &07611867116d0*x6)
      x26 = x5**1.66666666666667d0
      x27 = 0.0837224526465878d0*x26
      x28 = 0.0223069557172962d0*mu**5*x19*x25*x27
      x29 = x17*x6
      x30 = 5.57236303607474d0*mu*x12 + 1.0d0
      x31 = mu**3
      x32 = 2.0d0*x3*(x0 - 1.0d0)*log((3.3926025578013106d0*x11*x31 + 6.
     &7683429898050367d0*x29 + x30)/(3.1331792677937806d0*x29 + x30))
      x33 = x1*x15
      x34 = 0.000640558526213818d0*mu**8*x33*x5**2.66666666666667d0
      x35 = 0.148394183600699d0
      x36 = x25 - 1.0d0
      x37 = mu**6*(0.0167302167879722d0*x26*x35*x36 + 0.0108059179628909
     &d0*x33*x5** 2.0d0)
      x38 = 2.28942848510666d0*x20
      x39 = 1.63540853354893d0*x38*(-0.0259335011650457d0*x6 + 1.0d0)/( 
     &0.0524148278841779d0*x10 + 0.494402081358784d0*x6 + 1.0d0)
      x40 = x19*x31*(0.0315054072231411d0*x23*x25 + 0.00111534778586481d
     &0*x27*( 12.0d0*x22*(0.825481812223657d0*x10 - 4.49737346725955d0*x
     &6)*exp( -0.28165369188898165d0*x6)/x23 + x39/x9))
      x41 = mu**4*(-0.0683590130188305d0*x21*x33 - 0.13157433081915d0*x2
     &4*x36 + 0.00139418473233101d0*x35*x5**1.0d0*(1.63540853354893d0*x3
     &8 - x39 - 10.9027235569928d0*x8*(0.558025705063192d0*x10 - 0.35252
     &1395009435d0*x6)*exp(-0.49698248213959023d0*x6)))
      x42 = x16 - x18*(x28 - x32 + x34 + x37 + x40 - x41)
      x43 = exp(x13*x42) - 1.0d0
      x44 = (x16 + x18*(-x28 + x32 - x34 - x37 - x40 + x41))/(x14*x7)
      x45 = 3.53536660621036d0*gamma_cc*rho_c**(-2.33333333333333d0)*x1*
     &*(-3.78d0)* x44**2.78d0
      x46 = 1.0d0 + x45/x43
      E = -rho_c*(-x4*log(x45*x46/(12.498817040307328d0*gamma_cc**2*rho_
     &c**( -4.666666666666667d0)*x1**(-7.5599999999999996d0)*x44** 5.559
     &9999999999996d0/x43**2 + x46) + 1.0d0) + x42)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x149, x150, x151, x152, x153, x154, x155, x156, x157, x158, x159
     &, x160, x161, x162, x163, x164
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 2.14502939711103d0
      x1 = 1d0/rho_c
      x2 = x1**0.666666666666667d0
      x3 = x2/x0
      x4 = 0.564189583547756d0*sqrt(x1)
      x5 = 1.46459188756152d0
      x6 = x1**0.333333333333333d0
      x7 = x6/x5
      x8 = 0.826307487110758d0*x1**0.166666666666667d0
      x9 = 0.406913004517529d0*x3 + 1.41872281647967d0*x4 + 3.2595509194
     &2229d0*x7 + 7.24010193431683d0*x8
      x10 = 9.86960440108936d0
      x11 = 0.693147180559945d0
      x12 = -x11 + 1.0d0
      x13 = x10/x12
      x14 = x13/x9
      x15 = log(0.5d0*x14 + 1.0d0)
      x16 = 1d0/x10
      x17 = x12*x16
      x18 = 0.194159335344114d0*x7 + 1.0d0
      x19 = 2.0d0*x18
      x20 = x15*x17*x19
      x21 = mu**2
      x22 = x21*x3
      x23 = 0.508616435555896d0*x22 + 1.0d0
      x24 = x23**(-4)
      x25 = mu**4
      x26 = exp(-0.6836107611867116d0*x7)
      x27 = 4.60115111447049d0
      x28 = x1**1.33333333333333d0
      x29 = x28/x27
      x30 = 3.14159265358979d0
      x31 = x1**1.0d0
      x32 = x31/x30
      x33 = 0.009578475d0*x32
      x34 = 0.00126674656487366d0*x29 + 0.0676322201645715d0*x3 - x33 + 
     &0.0188171922990732d0*x7 + 1.0d0
      x35 = x26*x34
      x36 = x35 - 1.0d0
      x37 = 0.13157433081915d0*x32*x36
      x38 = 3.10428814221102d0*x28
      x39 = x12*x15
      x40 = x18*x39
      x41 = 0.022020833726518d0*x40
      x42 = x38*x41
      x43 = 0.148394183600699d0
      x44 = 2.28942848510666d0
      x45 = x27*x44
      x46 = exp(-0.49698248213959023d0*x7)
      x47 = x46*(0.558025705063192d0*x3 - 0.352521395009435d0*x7)
      x48 = x0*x47
      x49 = 0.0524148278841779d0*x3 + 0.494402081358784d0*x7 + 1.0d0
      x50 = 1d0/x49
      x51 = -0.0259335011650457d0*x7 + 1.0d0
      x52 = x50*x51
      x53 = 1.63540853354893d0*x52
      x54 = x1**1.0d0*x43*(-x45*x53 + 1.63540853354893d0*x45 - 10.902723
     &5569928d0* x48)
      x55 = 0.00139418473233101d0*x54
      x56 = x25*(-x37 - x42 + x55)
      x57 = x21*x7
      x58 = 5.57236303607474d0*mu*x8 + 1.0d0
      x59 = mu**3
      x60 = 3.39260255780131d0*x4*x59 + 6.76834298980504d0*x57 + x58
      x61 = x60/(3.13317926779378d0*x57 + x58)
      x62 = 2.0d0*log(x61)
      x63 = x16*(x11 - 1.0d0)
      x64 = x62*x63
      x65 = x1**1.66666666666667d0
      x66 = x43*x65
      x67 = x36*x66
      x68 = x1**2.0d0
      x69 = 1.0525946465532d0*x68
      x70 = 0.0102659822546843d0*x40
      x71 = mu**6*(0.0167302167879722d0*x67 + x69*x70)
      x72 = 1.41421356237310d0
      x73 = 0.0315054072231411d0*x31
      x74 = x35*x73
      x75 = exp(-0.28165369188898165d0*x7)
      x76 = x75*(0.825481812223657d0*x3 - 4.49737346725955d0*x7)
      x77 = 12.0d0*x30/x31
      x78 = x76*x77
      x79 = 1d0/x2
      x80 = x45*x79
      x81 = 0.0837224526465878d0*x65
      x82 = x81*(x53*x80 + x78)
      x83 = x59*x72*(x74 + 0.00111534778586481d0*x82)
      x84 = mu**5
      x85 = x35*x72
      x86 = x81*x85
      x87 = 0.0223069557172962d0*x84*x86
      x88 = x1**2.66666666666667d0
      x89 = 0.133841734303777d0*x88
      x90 = 0.00478594012208448d0*x40
      x91 = mu**8*x89*x90
      x92 = x71 + x83 + x87 + x91
      x93 = x24*(-x56 - x64 + x92)
      x94 = x20 - x93
      x95 = exp(x13*x94) - 1.0d0
      x96 = x56 + x64 - x71 - x83 - x87 - x91
      x97 = x20 + x24*x96
      x98 = 1d0/x18
      x99 = 1d0/x15
      x100 = x98*x99
      x101 = x100*x97
      x102 = x101**5.56d0
      x103 = 70524387.5934282d0*x12**(-7.56d0)
      x104 = rho_c**(-4.66666666666667d0)*x103
      x105 = 1.77226878060434d-7*gamma_cc**2*x104
      x106 = x102*x105
      x107 = 1d0/x95
      x108 = x101**2.78d0
      x109 = rho_c**(-2.33333333333333d0)
      x110 = 8397.87994635719d0*x12**(-3.78d0)
      x111 = 0.000420983227766184d0*x110
      x112 = x109*x111
      x113 = x108*x112
      x114 = gamma_cc*x113
      x115 = x107*x114 + 1.0d0
      x116 = x17*log(x114*x115/(x106/x95**2 + x115) + 1.0d0)
      x117 = 0.0691804894611514d0*x39*x6
      x118 = 0.129439556896076d0*x117
      x119 = (0.27127533634502d0*x3 + 0.709361408239834d0*x4 + 1.0865169
     &7314076d0*x7 + 1.20668365571947d0*x8)/(x9**2*(x14 + 2.0d0))
      x120 = x119*x19
      x121 = x22*x96/x23**5
      x122 = 1.35631049481572d0*x121
      x123 = x26*(-0.00168899541983155d0*x29 - 0.045088146776381d0*x3 + 
     &x33 - 0.00627239743302441d0*x7)
      x124 = mu**7
      x125 = x119*x18
      x126 = 0.928727172679123d0*x8
      x127 = mu*x7
      x128 = x26*(-0.000844497709915773d0*x29 - 0.0225440733881905d0*x3 
     &+ 0.0047892375d0*x32 + 0.113935126864452d0*x34*x7 - 0.003136198716
     &5122d0*x7)
      x129 = 1d0/x6
      x130 = 0.0141372897033723d0*x129*x30*x44*x50
      x131 = 1.09027235569928d0*x52*x80
      x132 = 1.63540853354893d0*x51*x80*(0.034943218589452d0*x3 + 0.1648
     &00693786261d0 *x7)/x49**2
      x133 = mu*x24*(-0.000290571663240495d0*x1**2.0d0*x25*x85 + 2.83060
     &463816878d-5* x1**3.0d0*x124*x39 - 0.0223069557172962d0*x123*x25*x
     &72*x81 - 0.0472353356922751d0*x124*x125*x89 + 0.356911291476739d0*
     &x124*x88 *x90 + x21*x72*(-x123*x73 - 0.00490180588786583d0*x28*x35
     & + x74 - 0.00111534778586481d0*x81*(2.41662179562687d0*x1**( -0.66
     &6666666666667d0)*x76 + x130 + x131 + x132 + x75*x77*( -0.550321208
     &149104d0*x3 + 1.49912448908652d0*x7) + x78) + 0.00185891297644135d
     &0*x82) + 0.037178259528827d0*x25*x86 + x59*( -0.217336917462899d0*
     &x125*x38 - 0.263148661638301d0*x128*x32 + 4.13905085628136d0*x28*x
     &41 + x37 + 0.00302075971817057d0*x39*x65 - 0.00232364122055169d0*x
     &54 + 0.00139418473233101d0*x66*( -10.9027235569928d0*x0*x46*x79*(-
     &0.372017136708795d0*x3 + 0.117507131669812d0*x7) - 1.8061542051453
     &6d0*x129*x47*x5 - x130 - x131 - x132 - 7.26848237132856d0*x48*x79 
     &+ 1.09027235569928d0*x80 )) + x84*(0.00047750955226877d0*x1**2.333
     &33333333333d0*x39 - 0.101321183642338d0*x125*x69 - 0.0334604335759
     &443d0*x128*x66 + 0.0278836946466203d0*x67 + 2.10518929310641d0*x68
     &*x70) + 2.0d0* x63*(-x126 - 2.25611432993501d0*x127 - 1.6963012789
     &0066d0*x21*x4 + x61*(x126 + 1.04439308926459d0*x127))/x60)
      x134 = x13*x97
      x135 = exp(x134)
      x136 = x135 - 1.0d0
      x137 = 1d0/x136
      x138 = x114*x137 + 1.0d0
      x139 = x136**(-2)
      x140 = gamma_cc/(x106*x139 + x138)
      x141 = x138*x140
      x142 = rho_c**(-3.33333333333333d0)*x108
      x143 = x138*x142
      x144 = 1d0/x97
      x145 = x7*x97*x98
      x146 = x119*x134*x99
      x147 = x144*(-0.359841968171091d0*x117 + 3.77054317558771d0*x121 +
     & 5.56d0*x125 + 2.78d0*x133 + 0.179920984085546d0*x145 - 2.78d0*x14
     &6)
      x148 = 0.000982294198121095d0*x110
      x149 = 0.000420983227766184d0*x137
      x150 = x135*(-x118 + x120 + x122 + x133)
      x151 = 82883.752878387d0*x12**(-4.78d0)*x149*x150
      x152 = x111*x147
      x153 = gamma_cc*x102
      x154 = rho_c**(-5.66666666666667d0)*x153
      x155 = x103*x137*x154
      x156 = x100*x94
      x157 = x156**5.56d0
      x158 = exp(x13*(x20 - x24*(x17*x62 + x25*(x37 + x42 - x55) + x92))
     &) - 1.0d0
      x159 = x112*x156**2.78d0
      x160 = gamma_cc/x158
      x161 = x159*x160 + 1.0d0
      x162 = 1d0/(x105*x157/x158**2 + x161)
      x163 = x159*x161
      x164 = gamma_cc*x162*x163
      E = -rho_c*(-x116 + x94)
      d1E(1) = -rho_c*(-x1*x118 + x1*x120 + x1*x122 + x1*x133 - 850.8831
     &36251163d0*x12 **(-2.78d0)*x140*(x108*x109*x141*x149*(246.71675933
     &8856d0*x12**( -8.56d0)*x139*x150*x154 + x142*x148 + x142*x151 - x1
     &42*x152 - 1.77226878060434d-7*x144*x155*(-0.719683936342183d0*x117
     & + 7.54108635117542d0*x121 + 11.12d0*x125 + 5.56d0*x133 + 0.359841
     &968171091d0*x145 - 5.56d0*x146) + 8.27058764282027d-7* x155) + 0.0
     &00420983227766184d0*x143*x147 - 0.000982294198121095d0 *x143 + x14
     &9*x154*(-x148 - x151 + x152))/(x113*x141 + 1.0d0)) + x116 - x20 + 
     &x93
      d1E(3) = rho_c*x162*x17*(1.77226878060434d-7*x104*x157*x160 - x107
     &*x164*( 3.54453756120869d-7*x104*x107*x153 + x113) + x163)/(x164 +
     & 1.0d0)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E, d2
     &E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x501, x502, x503, x504, x505, x506, x507, x508, x509, x510, x511
     &, x512, x513, x514, x515, x516, x517, x518, x519, x520
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 2.14502939711103d0
      x1 = 1d0/rho_c
      x2 = x1**0.666666666666667d0
      x3 = x2/x0
      x4 = 0.564189583547756d0*sqrt(x1)
      x5 = 1.41872281647967d0*x4
      x6 = 1.46459188756152d0
      x7 = x1**0.333333333333333d0
      x8 = x7/x6
      x9 = 0.826307487110758d0*x1**0.166666666666667d0
      x10 = 0.406913004517529d0*x3 + x5 + 3.25955091942229d0*x8 + 7.2401
     &0193431683d0 *x9
      x11 = 9.86960440108936d0
      x12 = 0.693147180559945d0
      x13 = -x12 + 1.0d0
      x14 = 1d0/x13
      x15 = x11*x14
      x16 = x15/x10
      x17 = log(0.5d0*x16 + 1.0d0)
      x18 = 1d0/x11
      x19 = x13*x18
      x20 = 0.194159335344114d0*x8 + 1.0d0
      x21 = 2.0d0*x20
      x22 = x17*x19*x21
      x23 = mu**2
      x24 = x23*x3
      x25 = 0.508616435555896d0*x24 + 1.0d0
      x26 = x25**(-4)
      x27 = mu**4
      x28 = exp(-0.6836107611867116d0*x8)
      x29 = 4.60115111447049d0
      x30 = x1**1.33333333333333d0
      x31 = x30/x29
      x32 = 3.14159265358979d0
      x33 = x1**1.0d0
      x34 = x33/x32
      x35 = 0.009578475d0*x34
      x36 = -x35
      x37 = 0.0676322201645715d0*x3 + 0.00126674656487366d0*x31 + x36 + 
     &0.0188171922990732d0*x8 + 1.0d0
      x38 = x28*x37
      x39 = x38 - 1.0d0
      x40 = x34*x39
      x41 = 0.13157433081915d0*x40
      x42 = 3.10428814221102d0*x30
      x43 = x13*x17
      x44 = x20*x43
      x45 = 0.022020833726518d0*x44
      x46 = x42*x45
      x47 = 2.28942848510666d0
      x48 = x29*x47
      x49 = 1.63540853354893d0*x48
      x50 = exp(-0.49698248213959023d0*x8)
      x51 = x50*(0.558025705063192d0*x3 - 0.352521395009435d0*x8)
      x52 = x0*x51
      x53 = 10.9027235569928d0*x52
      x54 = 0.0524148278841779d0*x3 + 0.494402081358784d0*x8 + 1.0d0
      x55 = 1d0/x54
      x56 = -0.0259335011650457d0*x8 + 1.0d0
      x57 = x55*x56
      x58 = 1.63540853354893d0*x57
      x59 = x48*x58
      x60 = 0.148394183600699d0
      x61 = x1**1.0d0*x60
      x62 = x61*(x49 - x53 - x59)
      x63 = 0.00139418473233101d0*x62
      x64 = x27*(-x41 - x46 + x63)
      x65 = x23*x8
      x66 = 5.57236303607474d0*mu*x9 + 1.0d0
      x67 = 3.13317926779378d0*x65 + x66
      x68 = 1d0/x67
      x69 = mu**3
      x70 = 3.39260255780131d0*x4*x69 + 6.76834298980504d0*x65 + x66
      x71 = x68*x70
      x72 = 2.0d0*log(x71)
      x73 = x18*(x12 - 1.0d0)
      x74 = x72*x73
      x75 = x1**1.66666666666667d0
      x76 = x60*x75
      x77 = x39*x76
      x78 = x1**2.0d0
      x79 = 1.0525946465532d0*x78
      x80 = 0.0102659822546843d0*x44
      x81 = mu**6*(0.0167302167879722d0*x77 + x79*x80)
      x82 = 1.41421356237310d0
      x83 = 0.179587122125167d0*x33
      x84 = 0.175432441092201d0*x83
      x85 = x38*x84
      x86 = exp(-0.28165369188898165d0*x8)
      x87 = x86*(0.825481812223657d0*x3 - 4.49737346725955d0*x8)
      x88 = x32/x33
      x89 = 12.0d0*x88
      x90 = x87*x89
      x91 = 1d0/x2
      x92 = x48*x91
      x93 = 0.0837224526465878d0*x75
      x94 = x93*(x58*x92 + x90)
      x95 = x69*x82*(x85 + 0.00111534778586481d0*x94)
      x96 = mu**5
      x97 = x38*x93
      x98 = x82*x97
      x99 = 0.0223069557172962d0*x96*x98
      x100 = x1**2.66666666666667d0
      x101 = 0.133841734303777d0*x100
      x102 = 0.00478594012208448d0*x44
      x103 = mu**8*x101*x102
      x104 = x103 + x81 + x95 + x99
      x105 = x104 - x64 - x74
      x106 = x105*x26
      x107 = -x106 + x22
      x108 = x107*x15
      x109 = exp(x108)
      x110 = x109 - 1.0d0
      x111 = 1d0/x110
      x112 = x13**(-3.78d0)
      x113 = 8397.87994635719d0*x112
      x114 = 0.000420983227766184d0*x113
      x115 = -x103 + x64 + x74 - x81 - x95 - x99
      x116 = x115*x26 + x22
      x117 = 1d0/x20
      x118 = 1d0/x17
      x119 = x117*x118
      x120 = x116*x119
      x121 = x120**2.78d0
      x122 = rho_c**(-2.33333333333333d0)
      x123 = x121*x122
      x124 = x114*x123
      x125 = gamma_cc*x124
      x126 = x111*x125 + 1.0d0
      x127 = x124*x126
      x128 = x110**(-2)
      x129 = rho_c**(-4.66666666666667d0)
      x130 = 70524387.5934283d0
      x131 = x120**5.56d0
      x132 = x13**(-7.56d0)
      x133 = x131*x132
      x134 = x130*x133
      x135 = x129*x134
      x136 = gamma_cc**2
      x137 = 1.77226878060434d-7*x136
      x138 = x135*x137
      x139 = x126 + x128*x138
      x140 = 1d0/x139
      x141 = gamma_cc*x140
      x142 = x127*x141
      x143 = x142 + 1.0d0
      x144 = x19*log(x143)
      x145 = 0.0691804894611514d0
      x146 = x145*x43
      x147 = x146*x7
      x148 = 0.129439556896076d0*x147
      x149 = 0.27127533634502d0*x3 + 0.709361408239834d0*x4 + 1.08651697
     &314076d0*x8 + 1.20668365571947d0*x9
      x150 = x16 + 2.0d0
      x151 = 1d0/x150
      x152 = x151/x10**2
      x153 = x149*x152
      x154 = x153*x21
      x155 = x25**(-5)
      x156 = x155*x24
      x157 = x115*x156
      x158 = 1.35631049481572d0*x157
      x159 = 1.09027235569928d0*x92
      x160 = 1d0/x7
      x161 = x160*x32*x47
      x162 = 0.0141372897033723d0*x161*x55
      x163 = x160*x6
      x164 = x163*x51
      x165 = 1.80615420514536d0*x164
      x166 = x52*x91
      x167 = 7.26848237132856d0*x166
      x168 = x50*(-0.372017136708795d0*x3 + 0.117507131669812d0*x8)
      x169 = x0*x168*x91
      x170 = 10.9027235569928d0*x169
      x171 = x57*x92
      x172 = 1.09027235569928d0*x171
      x173 = x54**(-2)
      x174 = 0.034943218589452d0*x3 + 0.164800693786261d0*x8
      x175 = x173*x174
      x176 = x175*x56*x92
      x177 = 1.63540853354893d0*x176
      x178 = x159 - x162 - x165 - x167 - x170 - x172 - x177
      x179 = 0.00139418473233101d0*x76
      x180 = 0.200908840802769d0*x75
      x181 = 0.0150354743280612d0*x43
      x182 = x37*x8
      x183 = 0.000844497709915773d0*x31
      x184 = 0.0047892375d0*x34
      x185 = 0.0225440733881905d0*x3
      x186 = 0.0031361987165122d0*x8
      x187 = -x183 + x184 - x185 - x186
      x188 = x28*(0.113935126864452d0*x182 + x187)
      x189 = x188*x34
      x190 = x30*x45
      x191 = 0.217336917462899d0
      x192 = x153*x20
      x193 = x191*x192
      x194 = x180*x181 - 0.263148661638301d0*x189 + 4.13905085628136d0*x
     &190 - x193* x42 + x41
      x195 = 0.00168899541983155d0*x31
      x196 = 0.045088146776381d0*x3
      x197 = 0.00627239743302441d0*x8
      x198 = x28*(-x195 - x196 - x197 + x35)
      x199 = x198*x84
      x200 = 0.122619224952946d0*x30
      x201 = x200*x38
      x202 = 0.0399758348639607d0*x201
      x203 = 0.00185891297644135d0*x94
      x204 = -0.550321208149104d0*x3 + 1.49912448908652d0*x8
      x205 = x86*x89
      x206 = 2.14502939711103d0*x1**(-0.666666666666667d0)
      x207 = x206*x87
      x208 = x162 + x172 + x177
      x209 = x93*(x204*x205 + 1.12661476755593d0*x207 + x208 + x90)
      x210 = 0.00111534778586481d0*x209
      x211 = x23*x82
      x212 = 0.928727172679123d0*x9
      x213 = mu*x8
      x214 = 2.25611432993501d0*x213
      x215 = x23*x4
      x216 = 1.69630127890066d0*x215
      x217 = x212 + 1.04439308926459d0*x213
      x218 = -x212 - x214 - x216 + x217*x71
      x219 = 2.0d0/x70
      x220 = x219*x73
      x221 = x218*x220
      x222 = x27*x82
      x223 = x222*x93
      x224 = 0.0223069557172962d0*x223
      x225 = x1**2.33333333333333d0
      x226 = 0.0681236923205142d0*x225
      x227 = 0.00700944907716014d0*x43
      x228 = x78*x80
      x229 = x188*x76
      x230 = 0.101321183642338d0
      x231 = x230*x79
      x232 = mu**7
      x233 = x1**3.0d0
      x234 = 0.00866220739124163d0*x232*x233
      x235 = 0.00326776364305339d0*x43
      x236 = x100*x232
      x237 = x102*x236
      x238 = x27*x98
      x239 = x222*x38
      x240 = 0.0571643564037363d0*x1**2.0d0
      x241 = x239*x240
      x242 = 0.0472353356922751d0
      x243 = x192*x242
      x244 = x101*x232
      x245 = x234*x235 + 0.356911291476739d0*x237 + 0.037178259528827d0*
     &x238 - 0.00508309165921971d0*x241 - x243*x244 + x96*(-x192*x231 + 
     &x226* x227 + 2.10518929310641d0*x228 - 0.0334604335759443d0*x229 +
     & 0.0278836946466203d0*x77)
      x246 = -x198*x224 + x211*(-x199 - x202 + x203 - x210 + x85) + x221
     & + x245 + x69 *(x178*x179 + x194 - 0.00232364122055169d0*x62)
      x247 = mu*x26
      x248 = x246*x247
      x249 = x116*x15
      x250 = exp(x249)
      x251 = x250 - 1.0d0
      x252 = x251**(-2)
      x253 = 1d0/x251
      x254 = x125*x253 + 1.0d0
      x255 = x138*x252 + x254
      x256 = 1d0/x255
      x257 = gamma_cc*x256
      x258 = x124*x254
      x259 = x257*x258
      x260 = rho_c**(-3.33333333333333d0)*x121
      x261 = 0.000982294198121095d0*x260
      x262 = x254*x261
      x263 = 1d0/x116
      x264 = 0.359841968171091d0*x7
      x265 = x146*x264
      x266 = 5.56d0*x20
      x267 = x153*x266
      x268 = 3.77054317558771d0*x157
      x269 = x117*x8
      x270 = x116*x269
      x271 = 0.179920984085546d0*x270
      x272 = 2.78d0*x118
      x273 = x153*x249
      x274 = x272*x273
      x275 = 2.78d0*x248
      x276 = -x265 + x267 + x268 + x271 - x274 + x275
      x277 = x263*x276
      x278 = 0.000420983227766184d0*x260
      x279 = x254*x277*x278
      x280 = 0.000982294198121095d0*x113
      x281 = 0.000420983227766184d0*x253
      x282 = 82883.752878387d0*x13**(-4.78d0)
      x283 = -x148 + x154
      x284 = x158 + x248 + x283
      x285 = x250*x284
      x286 = x282*x285
      x287 = x281*x286
      x288 = x114*x277
      x289 = -x280 - x287 + x288
      x290 = rho_c**(-5.66666666666667d0)
      x291 = gamma_cc*x290
      x292 = x131*x291
      x293 = x281*x292
      x294 = x289*x293
      x295 = x134*x253
      x296 = 8.27058764282027d-7*x291
      x297 = 3.54453756120869d-7*x252
      x298 = 696047806.176231d0*x13**(-8.56d0)
      x299 = x292*x298
      x300 = x285*x299
      x301 = 1.77226878060434d-7*x263
      x302 = 0.719683936342183d0*x147
      x303 = 11.12d0*x20
      x304 = x153*x303
      x305 = 7.54108635117542d0*x157
      x306 = 0.359841968171091d0*x270
      x307 = 5.56d0*x118
      x308 = x273*x307
      x309 = 5.56d0*x248
      x310 = -x302 + x304 + x305 + x306 - x308 + x309
      x311 = x291*x295
      x312 = x310*x311
      x313 = x260*x280
      x314 = x260*x288
      x315 = x260*x287 + x313 - x314
      x316 = x295*x296 + x297*x300 - x301*x312 + x315
      x317 = x254*x257
      x318 = x123*x281*x316*x317
      x319 = 850.883136251163d0*x13**(-2.78d0)
      x320 = exp(x15*(x22 - x26*(x104 + x19*x72 + x27*(x41 + x46 - x63))
     &)) - 1.0d0
      x321 = x107*x119
      x322 = x129*x130*x132*x321**5.56d0
      x323 = x114*x122*x321**2.78d0
      x324 = gamma_cc/x320
      x325 = x323*x324 + 1.0d0
      x326 = 1d0/(x137*x322/x320**2 + x325)
      x327 = x323*x325
      x328 = gamma_cc*x326*x327
      x329 = gamma_cc*x135
      x330 = 3.54453756120869d-7*x111
      x331 = x19*x326*(-x111*x328*(x124 + x329*x330) + 1.77226878060434d
     &-7*x322*x324 + x327)/(x328 + 1.0d0)
      x332 = x105*x156
      x333 = x61*(-x49 + x53 + x59)
      x334 = -x159 + x165 + x167 + x170 + x208
      x335 = x19*x219
      x336 = x218*x335
      x337 = x28*(x195 + x196 + x197 + x36)
      x338 = -x211*(x199 + x202 - x203 + x210 - x85) + x224*x337 + x245 
     &- x336 + x69* (-x179*x334 + x194 + 0.00232364122055169d0*x333)
      x339 = x247*x338
      x340 = x262 - x279 - x294 - x318
      x341 = 1d0/x143
      x342 = x256*x341
      x343 = gamma_cc*x319*x340*x342
      x344 = rho_c**(-2)
      x345 = 0.172586075861435d0*x147
      x346 = 0.682784063255296d0*x153
      x347 = x346*x7
      x348 = 0.129439556896076d0*x347
      x349 = x153*x8
      x350 = 0.129439556896076d0*x349
      x351 = x152*(0.452125560575033d0*x3 + 1.06404211235975d0*x4 + 1.44
     &868929752102d0*x8 + 1.40779759833938d0*x9)
      x352 = x21*x351
      x353 = x149**2/(x10**4*x150**2)
      x354 = x15*x353
      x355 = x21*x354
      x356 = x149*x151*(0.542550672690039d0*x3 + x5 + 2.17303394628153d0
     &*x8 + 2.41336731143894d0*x9)/x10**3
      x357 = x21*x356
      x358 = x27*x31/x25**6
      x359 = 2.29947269793409d0*x105*x358
      x360 = 2.26051749135954d0*x332
      x361 = x155*x3*x69
      x362 = 2.71262098963145d0*x338*x361
      x363 = 0.00464728244110338d0*x76
      x364 = 0.363424118566428d0*x92
      x365 = x50*(0.166966313186253d0*x3 - 0.105477574079378d0*x8)
      x366 = 2.22044604925031d-16*x164
      x367 = 3.61230841029072d0*x163*x168
      x368 = 2.42282745710952d0*x166
      x369 = 10.9027235569928d0*x0*x50*x91*(0.620028561181324d0*x3 - 0.1
     &56676175559749d0*x8)
      x370 = 14.5369647426571d0*x169
      x371 = 0.363424118566428d0*x171
      x372 = 0.0282745794067445d0*x161*x175
      x373 = 2.18054471139857d0*x176
      x374 = 1.63540853354893d0*x56*x92
      x375 = x173*x374*(0.0582386976490866d0*x3 + 0.219734258381682d0*x8
     &)
      x376 = x174*x374*(0.0698864371789039d0*x3 + 0.329601387572523d0*x8
     &)/x54**3
      x377 = 0.455740507457808d0*x8
      x378 = -0.151913502485936d0*x182 + 0.0259624262672375d0*x3*x37 + 0
     &.0375734556469842d0*x3 + 0.00197049465647014d0*x31 + x36 + 0.00418
     &159828868294d0*x8
      x379 = x28*(-x377*(x183 - x184 + x185 + x186) + x378)
      x380 = 0.263148661638301d0*x34
      x381 = x20*x351
      x382 = x191*x42
      x383 = x14*x20*x353
      x384 = x20*x356
      x385 = 0.148394183600699d0*x153*x180 + 0.200908840802769d0*x153*x7
     &6 - 0.803635363211076d0*x181*x75 + 0.526297323276602d0*x189 - 9.65
     &778533132318d0*x190 + 8.27810171256273d0*x193*x30 + x381*x382 - x3
     &82*x384 + 2.14502939711103d0*x383*x42 - 0.263148661638301d0* x40
      x386 = 0.0334604335759443d0*x76
      x387 = x153*x226
      x388 = x145*x387 + 4.21037858621281d0*x192*x230*x78 - 0.3633263590
     &42743d0*x225* x227 - 6.31556787931922d0*x228 + 0.111534778586481d0
     &*x229 + x231* x381 - x231*x384 + x383*x79 + 0.0691804894611514d0*x
     &387 - 0.074356519057654d0*x77
      x389 = 0.0799516697279215d0*x198*x200
      x390 = 0.350864882184401d0*x83
      x391 = x198*x390
      x392 = x28*(0.0751469112939684d0*x3 + 0.00394098931294027d0*x31 - 
     &0.01915695d0* x34 + 0.00836319657736588d0*x8)
      x393 = x392*x84
      x394 = 0.00910930363347549d0*x97
      x395 = 0.133252782879869d0*x201
      x396 = x38*x390
      x397 = 0.00495710127051027d0*x94
      x398 = 0.0037178259528827d0*x209
      x399 = x204*x86
      x400 = x371 - x372 - x373 + x375 - x376
      x401 = 0.00111534778586481d0*x93*(-0.154912426780983d0*x1**( -0.33
     &3333333333333d0)*x87 - x205*(0.917202013581841d0*x3 - 1.9988326521
     &1535d0*x8) - 2.25322953511185d0*x206*x399 - 0.751076511703951d0*x2
     &07 - 24.0d0*x399*x88 + x400)
      x402 = 1.08351503479231d0*x9
      x403 = 2.0d0*x212 + 2.0d0*x214 + 2.0d0*x216
      x404 = mu*x217
      x405 = x404*x68
      x406 = 3.00815243991335d0*x213 + 2.54445191835098d0*x215 + x402 - 
     &x403*x405 + x404*x70*(2.08878617852919d0*x213 + 1.85745434535825d0
     &*x9)/x67**2 - x71*(1.39252411901946d0*x213 + x402)
      x407 = mu*x218*x403/x70**2
      x408 = 0.0101661833184394d0*x222*x240
      x409 = 0.074356519057654d0*x223
      x410 = x153*x234
      x411 = x242*x244
      x412 = -4.52089344419912d-5*x1**2.33333333333333d0*x239 - x224*x39
     &2 - 0.0577480492749442d0*x232*x233*x235 + 0.713822582953479d0*x236
     &* x243 - 1.30867473541471d0*x237 - 0.0991420254102054d0*x238 + 0.0
     &23721094409692d0*x241 + 0.466194077035412d0*x244*x383 + x381* x411
     & - x384*x411 + 0.064503068866399d0*x410
      x413 = x247*(-x19*x407 - x211*(x389 - x391 + x393 + x394 - x395 + 
     &x396 + x397 - x398 - x401) - x335*x406 + x336*x405 + x337*x408 - x
     &337*x409 + x412 + x69*(-x179*(x364 + x365 - x366 + x367 - x368 + x
     &369 + x370 - x371 + x372 + x373 - x375 + x376) - 0.006196376588137
     &83d0*x333 + x334*x363 - x379*x380 + x385) + x96*(-x379*x386 + x388
     &))
      x414 = x265 - x267 - x268 - x271 + x274 - x275
      x415 = x263*x414
      x416 = x126*x415
      x417 = x280 + x287 - x288
      x418 = x111*x134
      x419 = x283 - 1.35631049481572d0*x332 + x339
      x420 = x109*x111
      x421 = x282*x420
      x422 = x419*x421
      x423 = x109*x128*x299*x419
      x424 = x260*x415
      x425 = x291*(x302 - x304 - x305 - x306 + x308 - x309)
      x426 = x114*x424 + x278*x422 + x296*x418 + x301*x418*x425 + x313 +
     & 3.54453756120869d-7*x423
      x427 = x111*x426
      x428 = x141*x427
      x429 = x126*x428
      x430 = x136*x256
      x431 = 0.00327431399373698d0*x113
      x432 = rho_c**(-4.33333333333333d0)*x121
      x433 = x431*x432
      x434 = 0.00196458839624219d0*x113
      x435 = x415*x434
      x436 = x126*x432
      x437 = x116**(-2)
      x438 = x114*x437
      x439 = x414**2*x438
      x440 = x117*x7
      x441 = 0.156226553915631d0*x112*x440
      x442 = x282*x432
      x443 = x118*x153
      x444 = 0.000420983227766184d0*x443
      x445 = x414*x419*x438
      x446 = gamma_cc*rho_c**(-6.66666666666667d0)
      x447 = x131*x446
      x448 = x253*x417
      x449 = x447*x448
      x450 = 0.000841966455532367d0*x113
      x451 = 0.479789290894789d0*x147
      x452 = x264*x346
      x453 = 0.359841968171091d0*x349
      x454 = x266*x351
      x455 = x266*x354
      x456 = x266*x356
      x457 = x115*x358
      x458 = 6.39253410025677d0*x457
      x459 = 6.28423862597952d0*x157
      x460 = x116*x3/x20**2
      x461 = 0.023288892456339d0*x460
      x462 = 0.239894645447394d0*x270
      x463 = x249*x351
      x464 = x272*x463
      x465 = 97.4090910340024d0*x116*x353/x13**2
      x466 = x465/x17**2
      x467 = 5.56d0*x466
      x468 = x272*x465
      x469 = 6.73880859569814d0*x120*x14*x153
      x470 = x264*x469
      x471 = x249*x356
      x472 = x272*x471
      x473 = x246*x361
      x474 = 7.54108635117542d0*x473
      x475 = x148 - x154 - x158 - x248
      x476 = x269*x475
      x477 = 0.359841968171091d0*x476
      x478 = x15*x153*x475
      x479 = x307*x478
      x480 = x28*(x187*x377 + x378)
      x481 = x247*(-x198*x408 + x198*x409 + x211*(-x389 + x391 - x393 - 
     &x394 + x395 - x396 - x397 + x398 + x401) + x220*x406 - x221*x405 +
     & x407*x73 + x412 + x69*(-x178*x363 + x179*(-x364 - x365 + x366 - x
     &367 + x368 - x369 - x370 + x400) - x380*x480 + x385 + 0.0061963765
     &8813783d0* x62) + x96*(-x386*x480 + x388))
      x482 = 2.78d0*x481
      x483 = x114*x263
      x484 = x483*(-x451 + x452 + x453 + x454 + x455 - x456 - x458 + x45
     &9 - x461 + x462 - x464 - x467 - x468 + x470 + x472 - x474 + x477 -
     & x479 - x482)
      x485 = x260*x434
      x486 = 1.65411752856405d-6*x291
      x487 = 0.000841966455532367d0*x422
      x488 = 818029.852187331d0*x13**(-5.78d0)
      x489 = x419**2*x488
      x490 = x282*x444
      x491 = x345 - x348 - x350 - x352 - x355 + x357
      x492 = x114*x131
      x493 = x295*x446
      x494 = x253*x286
      x495 = x298*x447
      x496 = x252*x285*x495
      x497 = x284**2
      x498 = x432*x488
      x499 = x250*x281
      x500 = x497*exp(2.0d0*x249)
      x501 = 6869716491.20552d0*x13**(-9.56d0)*x447
      x502 = x250*x297
      x503 = x277*x432
      x504 = x432*x438
      x505 = x263*x310
      x506 = 1.77226878060434d-7*x437*x493
      x507 = -2.26051749135954d0*x157 + 2.29947269793409d0*x457 + 2.7126
     &2098963145d0* x473 + x481 + x491
      x508 = x140*x341
      x509 = x253*x329
      x510 = x124 + 3.54453756120869d-7*x509
      x511 = x253*x510
      x512 = x258 - x259*x511 + 1.77226878060434d-7*x509
      x513 = x254*x313
      x514 = x257*x511
      x515 = x252*x430
      x516 = x253*x259
      x517 = x123*x450
      x518 = x258*(7.08907512241738d-7*x509 + x517)/x255**2
      x519 = rho_c*x342
      x520 = 3.54453756120869d-7*x135
      E = -rho_c*(x107 - x144)
      d1E(1) = -rho_c*(-x1*x148 + x1*x154 + x1*x158 + x1*x248 - x257*x31
     &9*(-x262 + x279 + x294 + x318)/(x259 + 1.0d0)) + x106 + x144 - x22
     &
      d1E(3) = rho_c*x331
      d2E(1) = rho_c*(gamma_cc*x19*x508*(gamma_cc*rho_c**(-6.66666666666
     &667d0)*x111* x492*(0.000841966455532367d0*x128*x489*exp(2.0d0*x108
     &) + x415* x441 + x415*x487 - x415*x490 - 0.000420983227766184d0*x4
     &20*x489 - 0.000420983227766184d0*x421*(-x359 + x360 + x362 + x413 
     &+ x491) + 0.00196458839624219d0*x422 + x431 + x435 + x439 + x445 -
     & x484) + x126*x433 + x127*x128*x136*x426*(x134*x263*x330*x425 + x2
     &60*x487 + x418*x486 + 7.08907512241738d-7*x423 + x424*x450 + x485)
     &/x139** 2 - x131*x136*x140*x290*x427*x448*x450 - x142*x253*( -0.55
     &2318141716643d0*x133*x253*x440*x446*x505 + 0.000841966455532367d0*
     &x252*x498*x500 + x253*x301*x310*x443*x495 + x276**2*x504 - x276*x2
     &84*x504 - x284*x310*x506 + x301*x493*( 11.12d0*x118*x478 + 0.95957
     &8581789577d0*x147 - 12.568477251959d0* x157 - 0.479789290894789d0*
     &x270 - x303*x351 - x303*x354 + x303* x356 + x307*x463 + x307*x465 
     &- x307*x471 - 0.719683936342183d0* x347 - 0.719683936342183d0*x349
     & + 12.7850682005135d0*x457 + 0.046577784912678d0*x460 + 11.12d0*x4
     &66 - 0.719683936342183d0* x469*x7 + 15.0821727023508d0*x473 - 0.71
     &9683936342183d0*x476 + 5.56d0*x481) + x310**2*x506 + x432*x483*(x4
     &51 - x452 - x453 - x454 - x455 + x456 + x458 - x459 + x461 - x462 
     &+ x464 + x467 + x468 - x470 - x472 + x474 - x477 + x479 + x482) + 
     &0.00196458839624219d0*x432*x494 + x433 - x434*x503 - x441*x503 - x
     &442*x499*x507 + x490*x503 - 1.65411752856405d-6*x493*x505 + 4.6866
     &6633093149d-6*x493 - 0.000841966455532367d0*x494*x503 - x495*x502*
     &x507 - 7.08907512241738d-7*x496*x505 + 3.30823505712811d-6*x496 - 
     &x497*x498*x499 - x497*x501*x502 + 1.06336126836261d-6*x500*x501/x2
     &51**3) - x260*x416*x428*x450 + x415*x449*x450 + x416*x432*x441 - x
     &416*x442*x444 - x429*x485 + x434*x449 + x435*x436 + x436*x439 + x4
     &36*x445 - x436*x484) - 7145614.42661715d0*x13**(-6.56d0)*x140*x340
     &*x430*( -0.000420983227766184d0*x123*x429 + x126*x261 + x278*x416 
     &+ x293* x417)/x143**2 - x344*x345 + x344*x348 + x344*x350 + x344*x
     &352 + x344*x355 - x344*x357 + x344*x359 - x344*x360 - x344*x362 - 
     &x344* x413) + 0.258879113792152d0*x1*x147 - 4.0d0*x1*x192 + 2.7126
     &2098963145d0*x1*x332 - 2.0d0*x1*x339 - 2.0d0*x343
      d2E(4) = x331 - x519*(x19*(-1.77226878060434d-7*x135*x316*x515 + x
     &136*x252*x316* x518 - x253*x289*x292*x450 - x254*x314 - 1.77226878
     &060434d-7*x277 *x311 + x289*x290*x492*x510*x515 + 4.13529382141014
     &d-7*x311 + x314*x317*x511 - x316*x516 - x513*x514 + x513 - x516*( 
     &7.08907512241738d-7*x252*x300 - 3.54453756120869d-7*x263*x312 + x2
     &95*x486 + x315)) - x343*x512)
      d2E(6) = -x19*x519*(x253*(-gamma_cc*x511*x518 + 88.3758007628711d0
     &*rho_c**(-7.0d0 )*x120**8.34d0*x13**(-11.34d0)*x253*x317 + x254*x2
     &56*x510*x517 + x514*x520 - x520) + x508*x512*(1.77226878060434d-7*
     &x111*x329 + x127 - x142*x511))
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu, 
     &E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x
     &79, x80, x81, x82, x83
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = 0.682784063255296d0*x0**0.333333333333333d0
      x2 = 2.14502939711103d0
      x3 = x0**0.666666666666667d0
      x4 = x3/x2
      x5 = 0.564189583547756d0*sqrt(x0)
      x6 = 0.826307487110758d0*x0**0.166666666666667d0
      x7 = 9.86960440108936d0
      x8 = 0.693147180559945d0
      x9 = -x8 + 1.0d0
      x10 = 1d0/x9
      x11 = x10*x7
      x12 = 2.0d0*(0.194159335344114d0*x1 + 1.0d0)*log(0.5d0*x11/( 3.259
     &5509194222925d0*x1 + 0.40691300451752932d0*x4 + 1.4187228164796675
     &d0*x5 + 7.2401019343168302d0*x6) + 1.0d0)
      x13 = 1d0/x7
      x14 = x13*x9
      x15 = x12*x14
      x16 = rho_s**4/rho_c**4
      x17 = rho_s*x0
      x18 = x17 + 1.0d0
      x19 = -x17 + 1.0d0
      x20 = 1.92366105093154d0*x18**1.33333333333333d0 + 1.9236610509315
     &4d0*x19** 1.33333333333333d0 - 3.84732210186307d0
      x21 = 0.0197517897025652d0*x20*(0.101077332976288d0*x1 + 1.0d0)*(-
     &x16 + 1.0d0) *log(1.0d0 + 29.608574643216677d0/(3.2918048099450625
     &d0*x1 + 0.41002507094961249d0*x4 + 0.76232752193529008d0*x5 + 9.87
     &21297225692712d0*x6))
      x22 = (0.186690969707574d0*x1 + 1.0d0)*log(x11/(5.630984149097876d
     &0*x1 + 0.51606646454786342d0*x4 + 2.9152147142191778d0*x5 + 13.457
     &913714394447d0*x6) + 1.0d0)
      x23 = x14*x16*x20
      x24 = x23*(-x12 + x22)
      x25 = rho_s**2/rho_c**2
      x26 = -x25 + 1.0d0
      x27 = 4.60115111447049d0
      x28 = x0**1.33333333333333d0/x27
      x29 = 3.14159265358979d0
      x30 = x0**1.0d0
      x31 = x30/x29
      x32 = (0.0188171922990732d0*x1 + 0.00126674656487366d0*x28 - 0.009
     &578475d0*x31 + 0.0676322201645715d0*x4 + 1.0d0)*exp(-0.68361076118
     &67116d0*x1)
      x33 = x26*x32
      x34 = 1.41421356237310d0
      x35 = x0**1.66666666666667d0
      x36 = 0.0837224526465878d0*x35
      x37 = 0.0892278228691848d0*mu**5*x34*x36
      x38 = x18**0.666666666666667d0 + x19**0.666666666666667d0
      x39 = x38**3
      x40 = mu**2
      x41 = x1*x40/x38**2
      x42 = 11.1447260721495d0*mu*x6/x38 + 1.0d0
      x43 = mu**3
      x44 = x39*log((27.073371959220147d0*x41 + x42 + 27.140820462410485
     &d0*x43*x5/ x39)/(12.532717071175123d0*x41 + x42))
      x45 = 0.0315054072231411d0*x30
      x46 = 1d0/x3
      x47 = 1d0/x18
      x48 = x47**0.666666666666667d0
      x49 = 0.0524148278841779d0*x4
      x50 = x1*x47**0.333333333333333d0
      x51 = 2.28942848510666d0*x27
      x52 = 0.817704266774463d0*x51
      x53 = x18**2*x52*(-0.0259335011650457d0*x50 + 1.0d0)/(x48*(x48*x49
     & + 0.494402081358784d0*x50 + 1.0d0))
      x54 = 1d0/x19
      x55 = x54**0.666666666666667d0
      x56 = x1*x54**0.333333333333333d0
      x57 = x19**2*x52*(-0.0259335011650457d0*x56 + 1.0d0)/(x55*(x49*x55
     & + 0.494402081358784d0*x56 + 1.0d0))
      x58 = 0.00111534778586481d0*x36*(12.0d0*x26*x29*(-4.49737346725955
     &d0*x1 + 0.825481812223657d0*x4)*exp(-0.28165369188898165d0*x1)/x30
     & + x46* x53 + x46*x57)
      x59 = 4.0d0*x34*x43
      x60 = -x15 + x21
      x61 = -x24 + x60
      x62 = 0.0126441184989504d0*mu**8*x0**2.66666666666667d0
      x63 = 0.148394183600699d0
      x64 = x26*(x32 - 1.0d0)
      x65 = 0.0167302167879722d0*x35*x63*x64
      x66 = x23*(x12 - x22) + x60
      x67 = 0.0533250677421794d0*x0**2.0d0*x66
      x68 = 4.0d0*mu**6
      x69 = -0.13157433081915d0*x31*x64
      x70 = 0.817704266774463d0*x51*(x18**2.66666666666667d0 + x19** 2.6
     &6666666666667d0)
      x71 = 10.9027235569928d0*x2*x26*(-0.352521395009435d0*x1 + 0.55802
     &5705063192d0 *x4)*exp(-0.49698248213959023d0*x1)
      x72 = 0.00139418473233101d0*x0**1.0d0*x63
      x73 = 1.55214407110551d0*x28
      x74 = 4.0d0*mu**4
      x75 = 0.25d0/(0.508616435555896d0*x4*x40 + 1.0d0)**4
      x76 = x75*(-x14*x44 - x33*x37 - x59*(x33*x45 + x58) + x61*x62 - x6
     &8*(x65 - x67 ) + x74*(x61*x73 + x69 - x72*(x53 + x57 - x70 + x71))
     &)
      x77 = x18**0.666666666666667d0 + x19**0.666666666666667d0
      x78 = x77**3
      x79 = x32*(x25 - 1.0d0)
      x80 = -1.0d0 + exp(-8.0d0*x11*(x66 - x75*(x13*x44*(x8 - 1.0d0) + x
     &37*x79 + x59 *(x45*x79 - x58) + x62*x66 + x68*(-x65 + x67) + x74*(
     &x66*x73 + x69 + x72*(-x53 - x57 + x70 - x71))))/x78)
      x81 = (x61 - x76)/x61
      x82 = 0.167187784415092d0*gamma_cc*rho_c**(-2.33333333333333d0)*x1
     &0*x81** 2.78d0/x77**2
      x83 = 1.0d0 + x82/x80
      E = -rho_c*(-0.125d0*x14*x78*log(x82*x83/(0.027951755257627309d0*g
     &amma_cc**2 *rho_c**(-4.666666666666667d0)*x81**5.5599999999999996d
     &0/(x77**4* x80**2*x9**2) + x83) + 1.0d0) + x15 - x21 + x24 + x76)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu
     &, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &, x446, x447, x448, x449, x450, x451, x452, x453, x454, x455, x456
     &, x457, x458, x459, x460, x461, x462, x463, x464, x465, x466, x467
     &, x468, x469, x470, x471, x472, x473, x474, x475, x476, x477, x478
     &, x479, x480, x481, x482, x483, x484, x485, x486, x487, x488, x489
     &, x490, x491, x492, x493, x494, x495, x496, x497, x498, x499, x500
     &, x501, x502, x503, x504, x505, x506, x507, x508, x509
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/x1
      x3 = 1.46459188756152d0
      x4 = 1d0/rho_c
      x5 = x4**0.333333333333333d0
      x6 = x5/x3
      x7 = 0.186690969707574d0*x6 + 1.0d0
      x8 = 2.14502939711103d0
      x9 = x4**0.666666666666667d0
      x10 = x9/x8
      x11 = 0.564189583547756d0*sqrt(x4)
      x12 = 0.826307487110758d0*x4**0.166666666666667d0
      x13 = 0.516066464547863d0*x10 + 2.91521471421918d0*x11 + 13.457913
     &7143944d0* x12 + 5.63098414909788d0*x6
      x14 = 9.86960440108936d0
      x15 = x14*x2
      x16 = 1.0d0 + x15/x13
      x17 = log(x16)
      x18 = x17*x7
      x19 = 0.406913004517529d0*x10 + 1.41872281647967d0*x11 + 7.2401019
     &3431683d0* x12 + 3.25955091942229d0*x6
      x20 = x15/x19
      x21 = log(0.5d0*x20 + 1.0d0)
      x22 = 0.194159335344114d0*x6 + 1.0d0
      x23 = 2.0d0*x22
      x24 = x21*x23
      x25 = x18 - x24
      x26 = 1d0/x14
      x27 = x1*x26
      x28 = rho_s*x4
      x29 = x28 + 1.0d0
      x30 = -x28 + 1.0d0
      x31 = 1.92366105093154d0*x29**1.33333333333333d0 + 1.9236610509315
     &4d0*x30** 1.33333333333333d0 - 3.84732210186307d0
      x32 = rho_c**(-4)
      x33 = rho_s**4
      x34 = x32*x33
      x35 = x31*x34
      x36 = x27*x35
      x37 = x25*x36
      x38 = x24*x27
      x39 = -x34 + 1.0d0
      x40 = 0.101077332976288d0*x6 + 1.0d0
      x41 = 0.410025070949612d0*x10 + 0.76232752193529d0*x11 + 9.8721297
     &2256927d0* x12 + 3.29180480994506d0*x6
      x42 = 1.0d0 + 29.6085746432167d0/x41
      x43 = log(x42)
      x44 = x40*x43
      x45 = x39*x44
      x46 = x31*x45
      x47 = 0.0197517897025652d0*x46
      x48 = -x38 + x47
      x49 = -x37 + x48
      x50 = 1d0/x49
      x51 = mu**2
      x52 = x10*x51
      x53 = 0.508616435555896d0*x52 + 1.0d0
      x54 = x53**(-4)
      x55 = rho_c**(-2)
      x56 = rho_s**2*x55
      x57 = -x56 + 1.0d0
      x58 = exp(-0.6836107611867116d0*x6)
      x59 = 4.60115111447049d0
      x60 = x4**1.33333333333333d0
      x61 = x60/x59
      x62 = 3.14159265358979d0
      x63 = x4**1.0d0
      x64 = x63/x62
      x65 = 0.009578475d0*x64
      x66 = -x65
      x67 = 0.0676322201645715d0*x10 + 0.0188171922990732d0*x6 + 0.00126
     &674656487366d0*x61 + x66 + 1.0d0
      x68 = x58*x67
      x69 = x57*x68
      x70 = 1.41421356237310d0
      x71 = mu**5*x70
      x72 = x69*x71
      x73 = x4**1.66666666666667d0
      x74 = 0.0837224526465878d0*x73
      x75 = 0.0892278228691848d0*x74
      x76 = x72*x75
      x77 = x29**0.666666666666667d0 + x30**0.666666666666667d0
      x78 = x77**3
      x79 = x77**2
      x80 = 1d0/x79
      x81 = x51*x80
      x82 = x6*x81
      x83 = 1d0/x77
      x84 = mu*x83
      x85 = 11.1447260721495d0*x12*x84 + 1.0d0
      x86 = mu**3
      x87 = x11/x78
      x88 = 27.0733719592201d0*x82 + x85 + 27.1408204624105d0*x86*x87
      x89 = x88/(12.5327170711751d0*x82 + x85)
      x90 = log(x89)
      x91 = x78*x90
      x92 = x27*x91
      x93 = 0.179587122125167d0
      x94 = 0.175432441092201d0*x63*x93
      x95 = x69*x94
      x96 = 0.825481812223657d0*x10 - 4.49737346725955d0*x6
      x97 = exp(-0.28165369188898165d0*x6)
      x98 = x57*x96*x97
      x99 = x62/x63
      x100 = 12.0d0*x99
      x101 = x100*x98
      x102 = x29**2
      x103 = 1d0/x29
      x104 = x103**0.666666666666667d0
      x105 = 0.0524148278841779d0*x10
      x106 = x103**0.333333333333333d0
      x107 = x106*x6
      x108 = x104*x105 + 0.494402081358784d0*x107 + 1.0d0
      x109 = 1d0/x108
      x110 = x102*x109
      x111 = (-0.0259335011650457d0*x107 + 1.0d0)/x104
      x112 = x110*x111
      x113 = 1d0/x9
      x114 = 2.28942848510666d0
      x115 = x114*x59
      x116 = x113*x115
      x117 = 0.817704266774463d0*x116
      x118 = x30**2
      x119 = 1d0/x30
      x120 = x119**0.666666666666667d0
      x121 = x119**0.333333333333333d0
      x122 = x121*x6
      x123 = x105*x120 + 0.494402081358784d0*x122 + 1.0d0
      x124 = 1d0/x123
      x125 = x118*x124
      x126 = (-0.0259335011650457d0*x122 + 1.0d0)/x120
      x127 = x125*x126
      x128 = x74*(x101 + x112*x117 + x117*x127)
      x129 = 0.00111534778586481d0*x128
      x130 = 4.0d0*x70*x86
      x131 = x130*(x129 + x95)
      x132 = 0.0472353356922751d0*mu**8*x4**2.66666666666667d0
      x133 = 0.267683468607554d0*x132
      x134 = x133*x49
      x135 = 0.148394183600699d0
      x136 = x135*x73
      x137 = x68 - 1.0d0
      x138 = x137*x57
      x139 = x136*x138
      x140 = 0.0167302167879722d0*x139
      x141 = 0.101321183642338d0
      x142 = x4**2.0d0
      x143 = -x18 + x24
      x144 = x143*x36
      x145 = x144 + x48
      x146 = x141*x142*x145
      x147 = 0.526297323276602d0*x146
      x148 = 4.0d0*mu**6
      x149 = x148*(x140 - x147)
      x150 = 0.13157433081915d0*x138*x64
      x151 = -x150
      x152 = x115*(x29**2.66666666666667d0 + x30**2.66666666666667d0)
      x153 = 0.817704266774463d0*x152
      x154 = 0.558025705063192d0*x10 - 0.352521395009435d0*x6
      x155 = exp(-0.49698248213959023d0*x6)
      x156 = x154*x155*x57
      x157 = x156*x8
      x158 = 10.9027235569928d0*x157
      x159 = 0.817704266774463d0*x115
      x160 = x112*x159
      x161 = x127*x159
      x162 = x135*x4**1.0d0
      x163 = 0.00139418473233101d0*x162
      x164 = 1.55214407110551d0*x61
      x165 = 4.0d0*mu**4
      x166 = x165*(x151 - x163*(-x153 + x158 + x160 + x161) + x164*x49)
      x167 = -x131 + x134 - x149 + x166 - x76 - x92
      x168 = x167*x54
      x169 = 0.25d0*x168
      x170 = x50*(-x169 + x49)
      x171 = x170**2.78d0
      x172 = rho_c**(-2.33333333333333d0)
      x173 = x171*x172
      x174 = 14.4549425392770d0
      x175 = 0.0115661327577615d0*x174
      x176 = x173*x175
      x177 = x56 - 1.0d0
      x178 = x177*x68
      x179 = x71*x75
      x180 = x26*(x0 - 1.0d0)
      x181 = x178*x94
      x182 = x153 - x158 - x160 - x161
      x183 = x130*(-x129 + x181) + x133*x145 + x148*(-x140 + x147) + x16
     &5*(x145*x164 + x151 + x163*x182) + x178*x179 + x180*x91
      x184 = x183*x54
      x185 = 0.25d0*x184
      x186 = x145 - x185
      x187 = x29**0.666666666666667d0 + x30**0.666666666666667d0
      x188 = x187**3
      x189 = 1d0/x188
      x190 = 8.0d0*x15*x189
      x191 = exp(-x186*x190)
      x192 = x191 - 1.0d0
      x193 = x2/x192
      x194 = x187**2
      x195 = 1d0/x194
      x196 = gamma_cc*x195
      x197 = x193*x196
      x198 = x176*x197 + 1.0d0
      x199 = 208.945363813799d0
      x200 = x170**5.56d0
      x201 = x199*x200
      x202 = x1**(-2)
      x203 = x202/x192**2
      x204 = rho_c**(-4.66666666666667d0)
      x205 = 0.000133775426970164d0*x204
      x206 = gamma_cc**2*x205/x187**4
      x207 = gamma_cc/(x198 + x201*x203*x206)
      x208 = x195*x198*x207
      x209 = x176*x2*x208 + 1.0d0
      x210 = x27*log(x209)
      x211 = 0.125d0*x188
      x212 = x38 - x47
      x213 = x212 + x37
      x214 = x169 + x213
      x215 = x131 - x134 + x149 - x166 + x76 + x92
      x216 = x1*x5/x174
      x217 = x21*x216
      x218 = rho_c**(-5)
      x219 = x218*x33
      x220 = x219*x31
      x221 = x31*x39
      x222 = x43*x6
      x223 = x221*x222
      x224 = x30**0.333333333333333d0
      x225 = x29**0.333333333333333d0
      x226 = x224 - x225
      x227 = rho_s*x55
      x228 = (0.27127533634502d0*x10 + 0.709361408239834d0*x11 + 1.20668
     &365571947d0* x12 + 1.08651697314076d0*x6)/(x19**2*(x20 + 2.0d0))
      x229 = rho_s**5*x226
      x230 = x25*x27
      x231 = x40*(0.273350047299742d0*x10 + 0.381163760967645d0*x11 + 1.
     &64535495376155d0*x12 + 1.09726826998169d0*x6)/(x41**2*x42)
      x232 = x221*x231
      x233 = x7*(0.344044309698576d0*x10 + 1.45760735710959d0*x11 + 2.24
     &298561906574d0*x12 + 1.87699471636596d0*x6)/(x13**2*x16)
      x234 = 0.062230323235858d0*x17*x216
      x235 = x228*x23
      x236 = 0.129439556896076d0*x217
      x237 = -x235 + x236
      x238 = x233 - x234 + x237
      x239 = x52/x53**5
      x240 = x56*x68
      x241 = 0.17845564573837d0*x71*x74
      x242 = x240*x241
      x243 = 0.00168899541983155d0*x61
      x244 = 0.045088146776381d0*x10
      x245 = 0.00627239743302441d0*x6
      x246 = x57*x58
      x247 = x246*(x243 + x244 + x245 + x66)
      x248 = x179*x247
      x249 = 0.00116228665296198d0*x4**2.0d0
      x250 = x249*x72
      x251 = 0.148713038115308d0*x74
      x252 = x251*x72
      x253 = x27*x79
      x254 = x29**(-0.333333333333333d0)
      x255 = x30**(-0.333333333333333d0)
      x256 = x254 - x255
      x257 = 2.0d0*x90
      x258 = x256*x257
      x259 = x253*x258
      x260 = x259*x28
      x261 = 0.713822582953479d0*x132
      x262 = x261*x49
      x263 = x6*x84
      x264 = 4.17757235705837d0*x263
      x265 = x256*x28
      x266 = mu*x265*x6*x80
      x267 = 16.7102894282335d0*x266
      x268 = 1.85745434535825d0*x12
      x269 = 7.42981738143299d0*x12
      x270 = x265*x269*x83
      x271 = -x268 + x270
      x272 = x11*x81
      x273 = -9.02445731974005d0*x263 + 54.281640924821d0*x265*x51*x87 +
     & 36.0978292789602d0*x266 + x271 - 13.5704102312052d0*x272
      x274 = mu/x88
      x275 = x253*x274*(x273 - x89*(-x264 + x267 + x271))
      x276 = x238*x35
      x277 = 4.0d0*x37
      x278 = x218*x229
      x279 = x230*x278
      x280 = x35*x44
      x281 = 0.0790071588102608d0*x280
      x282 = 0.000665486074881265d0*x223
      x283 = 0.58482233974552d0*x232
      x284 = 0.0506609980493537d0*x45
      x285 = x226*x284
      x286 = x28*x285
      x287 = x237 + x281 - x282 + x283 + x286
      x288 = x133*(-x276 + x277 - 2.56488140124205d0*x279 + x287)
      x289 = 24.0d0*x96*x97*x99
      x290 = x289*x56
      x291 = 2.41662179562687d0*x4**(-0.666666666666667d0)*x98
      x292 = x100*x57*x97*(-0.550321208149104d0*x10 + 1.49912448908652d0
     &*x6)
      x293 = x103*x28
      x294 = 1d0/x5
      x295 = 0.00706864485168613d0*x114*x294*x62
      x296 = x295/x106
      x297 = x110*x296*(-x293 + 1.0d0)
      x298 = x119*x28
      x299 = x295/x121
      x300 = x125*x299*(x298 + 1.0d0)
      x301 = 0.545136177849642d0*x116
      x302 = x112*x301
      x303 = x109*x29
      x304 = 2.18054471139857d0*x116
      x305 = x111*x303*x304
      x306 = x28*x305
      x307 = x127*x301
      x308 = x124*x30
      x309 = x126*x304*x308
      x310 = x28*x309
      x311 = 0.034943218589452d0*x10
      x312 = x104*x311
      x313 = 0.164800693786261d0*x107
      x314 = x312 + x313
      x315 = x111*x117/x108**2
      x316 = x102*x315*(-x293*x312 - x293*x313 + x314)
      x317 = x120*x311
      x318 = 0.164800693786261d0*x122
      x319 = x317 + x318
      x320 = x117*x126/x123**2
      x321 = x118*x320*(x298*x317 + x298*x318 + x319)
      x322 = 0.00111534778586481d0*x74
      x323 = 0.00490180588786583d0*x60
      x324 = 0.350864882184401d0*x63*x93
      x325 = 0.00185891297644135d0*x128 - x240*x324
      x326 = x130*(x247*x94 - x322*(x101 + x290 + x291 + x292 + x297 + x
     &300 + x302 - x306 + x307 + x310 + x316 + x321) - x323*x69 + x325 +
     & x95)
      x327 = x137*x56
      x328 = 0.0334604335759443d0*x136
      x329 = x327*x328
      x330 = 0.0278836946466203d0*x139
      x331 = x246*(-0.0225440733881905d0*x10 + 0.113935126864452d0*x6*x6
     &7 - 0.0031361987165122d0*x6 - 0.000844497709915773d0*x61 + 0.00478
     &92375d0*x64)
      x332 = x328*x331
      x333 = 1.0525946465532d0*x146
      x334 = 4.0d0*x144
      x335 = -x334
      x336 = x235 - x236
      x337 = x35*(-x233 + x234 + x336)
      x338 = x143*x27
      x339 = 2.56488140124205d0*x338
      x340 = x278*x339
      x341 = x287 + x335 + x337 + x340
      x342 = 0.526297323276602d0*x141*x142
      x343 = x341*x342
      x344 = x148*(x329 - x330 + x332 + x333 - x343)
      x345 = 0.263148661638301d0*x64
      x346 = x327*x345
      x347 = x331*x345
      x348 = 0.00232364122055169d0*x162*x182
      x349 = 2.06952542814068d0*x145*x61
      x350 = x113*x155*x8
      x351 = 21.8054471139857d0*x154*x350
      x352 = x30**1.66666666666667d0
      x353 = x29**1.66666666666667d0
      x354 = 2.18054471139857d0*x116
      x355 = -x297 - x300 - x302 + x306 - x307 - x310 - x316 - x321
      x356 = 0.00139418473233101d0*x136
      x357 = x356*(0.545136177849642d0*x113*x152 - 7.26848237132856d0*x1
     &13*x157 - 1.80615420514536d0*x156*x294*x3 + x28*x354*(x352 - x353)
     & - 10.9027235569928d0*x350*x57*(-0.372017136708795d0*x10 + 0.11750
     &7131669812d0*x6) - x351*x56 + x355)
      x358 = x164*x341
      x359 = x165*(x151 + x346 + x347 + x348 + x349 - x357 - x358)
      x360 = x30**(-0.333333333333333d0)
      x361 = x29**(-0.333333333333333d0)
      x362 = x360 - 1.0d0*x361
      x363 = x196*x2
      x364 = -x191 + 1.0d0
      x365 = x214/x213
      x366 = x199*x365**5.56d0
      x367 = x172*x175*x365**2.78d0
      x368 = x363/x364
      x369 = -x367*x368 + 1.0d0
      x370 = 1d0/(x202*x206*x366/x364**2 + x369)
      x371 = x367*x369
      x372 = x370*x371
      x373 = x363*x372 + 1.0d0
      x374 = x27*log(x373)
      x375 = rho_c**(-3.33333333333333d0)
      x376 = x198*x375
      x377 = rho_c**(-4.33333333333333d0)*rho_s
      x378 = 1d0/x187
      x379 = -x361
      x380 = x360 + x379
      x381 = x378*x380
      x382 = 0.015421510343682d0*x381
      x383 = 8.0d0*x22
      x384 = x21*x27*x383
      x385 = 0.0790071588102608d0*x46
      x386 = -x384 + x385
      x387 = -x168 - x277 + x386
      x388 = 1d0/x387
      x389 = -x242
      x390 = x54*(x248 - x250 + x252 + x260 - x262 - x275 + x288 + x326 
     &- x344 - x359 + x389)
      x391 = x387*x50
      x392 = 4.0d0*x276
      x393 = 1.35631049481572d0*x239
      x394 = x167*x393
      x395 = 16.0d0*x37
      x396 = 10.2595256049682d0*x279
      x397 = x228*x383
      x398 = 0.517758227584304d0*x217
      x399 = 0.316028635241043d0*x280
      x400 = 0.00266194429952506d0*x223
      x401 = 2.33928935898208d0*x232
      x402 = 0.202643992197415d0*x45
      x403 = x226*x402
      x404 = x28*x403
      x405 = -x397 + x398 + x399 - x400 + x401 + x404
      x406 = -x341*x391 - x390 - x392 - x394 + x395 - x396 + x405
      x407 = x388*x406
      x408 = 0.032153849066577d0*x376
      x409 = 1d0/x145
      x410 = x186*x409
      x411 = x410**2.78d0
      x412 = x174*x375
      x413 = 0.107950572405774d0*x412
      x414 = 0.0616860413747281d0*x174
      x415 = x177*x58*(-x243 - x244 - x245 + x65)
      x416 = x178*x71
      x417 = -x254 + x255
      x418 = x180*x79
      x419 = x54*(x130*(x178*x323 - x181 + x322*(-x101 - x290 - x291 - x
     &292 + x355) + x325 + x415*x94) + x133*x341 - x145*x261 + x148*(-x3
     &29 + x330 - x332 - x333 + x343) + x165*(x150 - x346 - x347 - x348 
     &- x349 + x357 + x358) + x179*x415 + x249*x416 - x251*x416 + x257*x
     &28*x417* x418 + x274*x418*(x273 + x89*(x264 - x267 + x268 - x270))
     & + x389)
      x420 = -x184 + x334 + x386
      x421 = x409*x420
      x422 = 4.0d0*x337
      x423 = x183*x393
      x424 = 16.0d0*x144
      x425 = 10.2595256049682d0*x338
      x426 = x278*x425
      x427 = 1d0/x420
      x428 = 0.128615396266308d0*x412*x427
      x429 = x184 + x335 + x384 - x385
      x430 = 1.0d0*x360
      x431 = 2.0d0*x378
      x432 = x28*x431
      x433 = x34 - 1.0d0
      x434 = x31*x433
      x435 = 0.202643992197415d0*x433*x44
      x436 = x397 - x398 - x399
      x437 = exp(x190*(-x144 + x185 + x212))
      x438 = 13.2006183429726d0
      x439 = x437 - 1.0d0
      x440 = 1d0/x439
      x441 = x2*x440
      x442 = x189*x441
      x443 = x375*x437*x438*x442
      x444 = 0.00289153318944038d0*x172
      x445 = x363*x440
      x446 = x381*x414
      x447 = rho_c**(-5.66666666666667d0)
      x448 = x197*x447
      x449 = x201*x448
      x450 = x189*x193
      x451 = 0.00142693788768175d0*x380
      x452 = 0.128615396266308d0*x388
      x453 = x379 + x430
      x454 = x390 + x392 + x394 - x395 + x396 + x400 - x401 - x404 + x42
     &0*x432*x453 + x436
      x455 = x191*x438*x450
      x456 = gamma_cc*x447
      x457 = 4.41396426734162d0*x456/x187**5
      x458 = x187/x373
      x459 = 1.46459188756152d0*x207
      x460 = 0.125d0*rho_c
      x461 = 0.0790071588102608d0*x44
      x462 = rho_s**3*x31
      x463 = x32*x462
      x464 = 4.0d0*x230
      x465 = x28*x68
      x466 = x241*x465
      x467 = 36.0978292789602d0*x263 + x269 + 54.281640924821d0*x272
      x468 = x89*(16.7102894282335d0*x263 + x269)
      x469 = x274*x77
      x470 = x462/rho_c**3
      x471 = x461*x470
      x472 = -x471
      x473 = x226*x230*x34
      x474 = -x285 - x464*x470 + x472 + 2.56488140124205d0*x473
      x475 = x137*x28
      x476 = x328*x475
      x477 = x324*x465
      x478 = x28*x289
      x479 = x296*x303
      x480 = x299*x308
      x481 = x29*x314*x315
      x482 = x30*x319*x320
      x483 = x305 - x309 + x479 - x480 + x481 - x482
      x484 = x345*x475
      x485 = x354*(-x352 + x353)
      x486 = x28*x351
      x487 = x54*(x130*(-x322*(-x478 + x483) + x477) + x133*x474 + x148*
     &(x342*x474 + x476) + x165*(x164*x474 - x356*(x483 - x485 - x486) +
     & x484) - x259 + x27*x469*(x256*x467 + x417*x468) + x466)
      x488 = 0.25d0*x4
      x489 = 0.316028635241043d0*x44*x470
      x490 = -x489
      x491 = 10.2595256049682d0*x473
      x492 = 16.0d0*x230*x470
      x493 = -x224 + x225
      x494 = x284*x493
      x495 = x34*x493
      x496 = x339*x495
      x497 = x338*x470
      x498 = 4.0d0*x497
      x499 = x472 + x494 + x496 + x498
      x500 = -x391*x499 - x403 - x487 + x490 + x491 - x492
      x501 = x362*x431
      x502 = x410**5.56d0
      x503 = x199*x502
      x504 = x425*x495
      x505 = 16.0d0*x497
      x506 = -x305 + x309 - x479 + x480 - x481 + x482
      x507 = x54*(x130*(x322*(x478 + x506) + x477) + x133*x499 + x148*(x
     &342*x499 + x476) + x165*(x164*x499 + x356*(x485 + x486 + x506) + x
     &484) + x180*x417*x469*(x467 - x468) + x258*x418 + x466)
      x508 = x402*x493 + x421*(x471 - x494 - x496 - x498) + x490 + x504 
     &+ x505 - x507
      x509 = x429*x501 + x435*x493 + x489 - x504 - x505 + x507
      E = -rho_c*(-x210*x211 + x214)
      d1E(1) = x211*x374 + 0.25d0*x215*x54 + x460*(-x171*x458*x459*(-x19
     &3*x208*x444*( gamma_cc*rho_c**(-6.66666666666667d0)*rho_s*x201*x45
     &0*x451 + x171 *x375*x454*x455 + x171*x377*x446 - x171*x406*x412*x4
     &52 + x171* x413 + x191*x200*x203*x454*x457 - 0.00297516549581645d0
     &*x407*x449 + 0.00249714130344307d0*x449) + x198*x377*x382 + 0.0269
     &876431014436d0*x376 - x407*x408 + x411*x444*x445*(-x377* x378*x414
     &*(-x360 + x361) + x413 - x428*(x405 - x419 + x421*(-x281 + x282 - 
     &x283 - x286 + x334 + x336 - x337 - x340) + x422 - x423 - x424 + x4
     &26) + x443*(-0.00266194429952506d0*x222*x434 + x226*x28* x435 + 2.
     &33928935898208d0*x231*x434 + x419 - x422 + x423 + x424 - x426 + x4
     &29*x432*(x361 - x430) + x436))) + 2.0d0*x194*x227*x362* x374 + 2.7
     &1262098963145d0*x215*x239*x4 + 1.03551645516861d0*x217* x4 - 16.0d
     &0*x22*x228*x4 + 32.0d0*x220*x230 - 8.0d0*x220*x238 + 0.63205727048
     &2087d0*x220*x44 - 0.00532388859905012d0*x223*x4 + 0.40528798439483
     &d0*x226*x227*x45 + 4.67857871796416d0*x232*x4 + 2.0d0*x4*x54*(x242
     & - x248 + x250 - x252 - x260 + x262 + x275 - x288 - x326 + x344 + 
     &x359) - 20.5190512099364d0*x229*x230/rho_c** 6) + x49
      d1E(2) = -rho_c*(-0.125d0*x187*x459*(x171*x376*x382 + x171*x388*x4
     &08*x500 - 0.00289153318944038d0*x173*x208*x441*(-x202*x437*x457*x5
     &02*x509/ x439**2 + 0.0616860413747281d0*x381*x411*x412 + x411*x428
     &*x508 - x411*x443*x509 + 0.00297516549581645d0*x427*x445*x447*x503
     &*x508 + x442*x451*x456*x503) + 0.00289153318944038d0*x200*x448*(x1
     &74*x452 *x500 + x446 - x455*(x403 - x420*x501 + x487 + x489 - x491
     & + x492 )))/x209 + x194*x210*x453*x488 - 2.56488140124205d0*x219*x
     &226* x230 + x285*x4 + x461*x463 + x463*x464 + x487*x488)
      d1E(3) = x26*x370*x458*x460*(-x197*x372*(x176 + 0.0002675508539403
     &28d0*x197*x201* x204) - x205*x366*x368 + x371)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu
     &, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
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
     &524, x1525, x1526, x1527, x1528, x1529, x1530, x1531, x1532, x1533
     &, x1534, x1535, x1536, x1537, x1538, x1539, x1540, x1541, x1542, x
     &1543, x1544, x1545, x1546, x1547, x1548, x1549, x1550, x1551, x155
     &2, x1553, x1554, x1555, x1556, x1557, x1558, x1559, x1560, x1561, 
     &x1562, x1563, x1564, x1565, x1566, x1567, x1568, x1569, x1570, x15
     &71, x1572, x1573, x1574, x1575, x1576, x1577, x1578, x1579, x1580,
     & x1581, x1582, x1583, x1584, x1585, x1586, x1587, x1588, x1589, x1
     &590, x1591, x1592, x1593, x1594, x1595, x1596, x1597, x1598, x1599
     &, x1600, x1601, x1602, x1603, x1604, x1605, x1606, x1607, x1608, x
     &1609, x1610, x1611, x1612, x1613, x1614, x1615, x1616, x1617, x161
     &8, x1619, x1620, x1621, x1622, x1623, x1624, x1625, x1626, x1627, 
     &x1628, x1629, x1630, x1631, x1632, x1633, x1634, x1635, x1636, x16
     &37, x1638, x1639, x1640, x1641, x1642, x1643, x1644, x1645, x1646,
     & x1647, x1648, x1649, x1650, x1651, x1652, x1653, x1654, x1655, x1
     &656, x1657, x1658, x1659, x1660, x1661, x1662, x1663, x1664, x1665
     &, x1666, x1667, x1668, x1669, x1670, x1671, x1672, x1673, x1674
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
      x7 = 0.693147180559945d0
      x8 = -x7 + 1.0d0
      x9 = 1d0/x8
      x10 = gamma_cc*x9
      x11 = x10*x6
      x12 = 1.46459188756152d0
      x13 = x0**0.333333333333333d0
      x14 = x13/x12
      x15 = 0.186690969707574d0*x14 + 1.0d0
      x16 = 2.14502939711103d0
      x17 = x0**0.666666666666667d0
      x18 = x17/x16
      x19 = 0.564189583547756d0*sqrt(x0)
      x20 = 2.91521471421918d0*x19
      x21 = 0.826307487110758d0*x0**0.166666666666667d0
      x22 = 5.63098414909788d0*x14 + 0.516066464547863d0*x18 + x20 + 13.
     &4579137143944d0*x21
      x23 = 9.86960440108936d0
      x24 = x23*x9
      x25 = 1.0d0 + x24/x22
      x26 = log(x25)
      x27 = x15*x26
      x28 = 1.41872281647967d0*x19
      x29 = 3.25955091942229d0*x14 + 0.406913004517529d0*x18 + 7.2401019
     &3431683d0* x21 + x28
      x30 = x24/x29
      x31 = log(0.5d0*x30 + 1.0d0)
      x32 = 0.194159335344114d0*x14 + 1.0d0
      x33 = 2.0d0*x32
      x34 = x31*x33
      x35 = x27 - x34
      x36 = 1d0/x23
      x37 = x36*x8
      x38 = 1.92366105093154d0*x2**1.33333333333333d0 + 1.92366105093154
     &d0*x3** 1.33333333333333d0 - 3.84732210186307d0
      x39 = rho_c**(-4)
      x40 = rho_s**4
      x41 = x39*x40
      x42 = x38*x41
      x43 = x37*x42
      x44 = x35*x43
      x45 = x34*x37
      x46 = 0.101077332976288d0*x14 + 1.0d0
      x47 = 0.76232752193529d0*x19
      x48 = 3.29180480994506d0*x14 + 0.410025070949612d0*x18 + 9.8721297
     &2256927d0* x21 + x47
      x49 = 1.0d0 + 29.6085746432167d0/x48
      x50 = log(x49)
      x51 = x46*x50
      x52 = -x41 + 1.0d0
      x53 = x38*x52
      x54 = x51*x53
      x55 = 0.0197517897025652d0*x54
      x56 = -x45 + x55
      x57 = -x44 + x56
      x58 = 1d0/x57
      x59 = mu**2
      x60 = x18*x59
      x61 = 0.508616435555896d0*x60 + 1.0d0
      x62 = x61**(-4)
      x63 = rho_c**(-2)
      x64 = rho_s**2
      x65 = x63*x64
      x66 = -x65 + 1.0d0
      x67 = exp(-0.6836107611867116d0*x14)
      x68 = 4.60115111447049d0
      x69 = x0**1.33333333333333d0
      x70 = x69/x68
      x71 = 3.14159265358979d0
      x72 = x0**1.0d0
      x73 = x72/x71
      x74 = 0.009578475d0*x73
      x75 = -x74
      x76 = 0.0188171922990732d0*x14 + 0.0676322201645715d0*x18 + 0.0012
     &6674656487366d0*x70 + x75 + 1.0d0
      x77 = x67*x76
      x78 = x66*x77
      x79 = x0**1.66666666666667d0
      x80 = 0.0837224526465878d0*x79
      x81 = 1.41421356237310d0
      x82 = mu**5*x81
      x83 = x80*x82
      x84 = 0.0892278228691848d0*x83
      x85 = x78*x84
      x86 = x2**0.666666666666667d0
      x87 = x3**0.666666666666667d0
      x88 = x86 + x87
      x89 = x88**3
      x90 = x88**2
      x91 = 1d0/x90
      x92 = x59*x91
      x93 = x14*x92
      x94 = 1d0/x88
      x95 = mu*x94
      x96 = 11.1447260721495d0*x21*x95 + 1.0d0
      x97 = 12.5327170711751d0*x93 + x96
      x98 = 1d0/x97
      x99 = mu**3
      x100 = 1d0/x89
      x101 = x100*x19
      x102 = 27.1408204624105d0*x101*x99 + 27.0733719592201d0*x93 + x96
      x103 = x102*x98
      x104 = log(x103)
      x105 = x104*x89
      x106 = x105*x37
      x107 = 0.179587122125167d0*x72
      x108 = x107*x77
      x109 = 0.175432441092201d0*x108
      x110 = x109*x66
      x111 = exp(-0.28165369188898165d0*x14)
      x112 = x111*(-4.49737346725955d0*x14 + 0.825481812223657d0*x18)
      x113 = x71/x72
      x114 = x112*x113
      x115 = 12.0d0*x66
      x116 = x114*x115
      x117 = x2**2
      x118 = 1d0/x2
      x119 = x118**0.333333333333333d0
      x120 = x119*x14
      x121 = 0.0259335011650457d0*x120
      x122 = -x121 + 1.0d0
      x123 = x118**0.666666666666667d0
      x124 = 0.0524148278841779d0*x18
      x125 = 0.494402081358784d0*x120 + x123*x124 + 1.0d0
      x126 = 1d0/x125
      x127 = 1d0/x123
      x128 = x126*x127
      x129 = x122*x128
      x130 = x117*x129
      x131 = 1d0/x17
      x132 = 2.28942848510666d0
      x133 = x132*x68
      x134 = x131*x133
      x135 = 0.817704266774463d0*x134
      x136 = x3**2
      x137 = 1d0/x3
      x138 = x137**0.333333333333333d0
      x139 = x138*x14
      x140 = 0.0259335011650457d0*x139
      x141 = -x140 + 1.0d0
      x142 = x137**0.666666666666667d0
      x143 = x124*x142 + 0.494402081358784d0*x139 + 1.0d0
      x144 = 1d0/x143
      x145 = 1d0/x142
      x146 = x144*x145
      x147 = x141*x146
      x148 = x136*x147
      x149 = x80*(x116 + x130*x135 + x135*x148)
      x150 = 0.00111534778586481d0*x149
      x151 = x81*x99
      x152 = 4.0d0*x151
      x153 = x152*(x110 + x150)
      x154 = 0.0472353356922751d0*mu**8*x0**2.66666666666667d0
      x155 = 0.267683468607554d0*x154
      x156 = x155*x57
      x157 = x77 - 1.0d0
      x158 = 0.148394183600699d0
      x159 = x158*x79
      x160 = x157*x159
      x161 = x160*x66
      x162 = 0.0167302167879722d0*x161
      x163 = -x27 + x34
      x164 = x163*x43
      x165 = x164 + x56
      x166 = 0.101321183642338d0*x0**2.0d0
      x167 = x165*x166
      x168 = 0.526297323276602d0*x167
      x169 = mu**6
      x170 = 4.0d0*x169
      x171 = x170*(x162 - x168)
      x172 = x157*x73
      x173 = 0.13157433081915d0*x172*x66
      x174 = -x173
      x175 = x133*(x2**2.66666666666667d0 + x3**2.66666666666667d0)
      x176 = 0.817704266774463d0*x175
      x177 = 0.558025705063192d0*x18
      x178 = 0.352521395009435d0*x14
      x179 = exp(-0.49698248213959023d0*x14)
      x180 = x179*(x177 - x178)
      x181 = x180*x66
      x182 = x16*x181
      x183 = 10.9027235569928d0*x182
      x184 = 0.817704266774463d0*x133
      x185 = x130*x184
      x186 = x148*x184
      x187 = -x176 + x183 + x185 + x186
      x188 = x0**1.0d0*x158
      x189 = 0.00139418473233101d0*x188
      x190 = 1.55214407110551d0*x70
      x191 = mu**4
      x192 = 4.0d0*x191
      x193 = x192*(x174 - x187*x189 + x190*x57)
      x194 = -x106 - x153 + x156 - x171 + x193 - x85
      x195 = x194*x62
      x196 = 0.25d0*x195
      x197 = x58*(-x196 + x57)
      x198 = x197**5.56d0
      x199 = x8**(-2)
      x200 = x65 - 1.0d0
      x201 = x200*x77
      x202 = x36*(x7 - 1.0d0)
      x203 = x109*x200
      x204 = x176 - x183 - x185 - x186
      x205 = x105*x202 + x152*(-x150 + x203) + x155*x165 + x170*(-x162 +
     & x168) + x192 *(x165*x190 + x174 + x189*x204) + x201*x84
      x206 = x205*x62
      x207 = 0.25d0*x206
      x208 = x165 - x207
      x209 = x4**3
      x210 = 1d0/x209
      x211 = x210*x24
      x212 = 8.0d0*x211
      x213 = exp(-x208*x212)
      x214 = x213 - 1.0d0
      x215 = x199/x214**2
      x216 = x198*x215
      x217 = 208.945363813799d0
      x218 = rho_c**(-4.66666666666667d0)*x217
      x219 = 0.000133775426970164d0*x218
      x220 = gamma_cc**2
      x221 = x4**(-4)
      x222 = x220*x221
      x223 = x219*x222
      x224 = 14.4549425392770d0
      x225 = rho_c**(-2.33333333333333d0)
      x226 = x224*x225
      x227 = x197**2.78d0
      x228 = 0.0115661327577615d0*x227
      x229 = x226*x228
      x230 = 1d0/x214
      x231 = x11*x230
      x232 = x229*x231 + 1.0d0
      x233 = x216*x223 + x232
      x234 = 1d0/x233
      x235 = x229*x232
      x236 = x234*x235
      x237 = x11*x236 + 1.0d0
      x238 = log(x237)
      x239 = x238*x37
      x240 = 0.125d0*x209
      x241 = x45 - x55
      x242 = x241 + x44
      x243 = x196 + x242
      x244 = x106 + x153 - x156 + x171 - x193 + x85
      x245 = x8/x224
      x246 = x245*x31
      x247 = x13*x246
      x248 = x38*x40
      x249 = rho_c**(-5)
      x250 = x249*x51
      x251 = x248*x250
      x252 = x14*x50
      x253 = x252*x53
      x254 = x3**0.333333333333333d0
      x255 = x2**0.333333333333333d0
      x256 = x254 - x255
      x257 = 0.40528798439483d0*x256
      x258 = x51*x52
      x259 = x258*x63
      x260 = 1.08651697314076d0*x14 + 0.27127533634502d0*x18 + 0.7093614
     &08239834d0* x19 + 1.20668365571947d0*x21
      x261 = x30 + 2.0d0
      x262 = 1d0/x261
      x263 = x262/x29**2
      x264 = x260*x263
      x265 = x264*x32
      x266 = rho_s**5
      x267 = x256*x266
      x268 = 20.5190512099364d0*x267
      x269 = rho_c**(-6)
      x270 = x35*x37
      x271 = x269*x270
      x272 = x46*x53
      x273 = 1.09726826998169d0*x14 + 0.273350047299742d0*x18 + 0.381163
     &760967645d0* x19 + 1.64535495376155d0*x21
      x274 = 1d0/x49
      x275 = x274/x48**2
      x276 = x273*x275
      x277 = x272*x276
      x278 = x0*x277
      x279 = x249*x270
      x280 = x248*x279
      x281 = 1.87699471636596d0*x14 + 0.344044309698576d0*x18 + 1.457607
     &35710959d0* x19 + 2.24298561906574d0*x21
      x282 = 1d0/x25
      x283 = x22**(-2)
      x284 = x15*x282*x283
      x285 = x281*x284
      x286 = 0.062230323235858d0*x13
      x287 = x245*x26
      x288 = x286*x287
      x289 = x264*x33
      x290 = 0.129439556896076d0*x13
      x291 = x246*x290
      x292 = -x289 + x291
      x293 = x285 - x288 + x292
      x294 = x249*x293
      x295 = x248*x294
      x296 = x60/x61**5
      x297 = x0*x244*x296
      x298 = 2.0d0*x0
      x299 = x65*x77
      x300 = 0.17845564573837d0*x83
      x301 = x299*x300
      x302 = 0.00168899541983155d0*x70
      x303 = 0.045088146776381d0*x18
      x304 = 0.00627239743302441d0*x14
      x305 = x302 + x303 + x304 + x75
      x306 = x66*x67
      x307 = x305*x306
      x308 = x307*x84
      x309 = x78*x82
      x310 = 0.0571643564037363d0*x0**2.0d0
      x311 = 0.0203323666368788d0*x310
      x312 = x309*x311
      x313 = 0.148713038115308d0*x83
      x314 = x313*x78
      x315 = x2**(-0.333333333333333d0)
      x316 = x3**(-0.333333333333333d0)
      x317 = x315 - x316
      x318 = x37*x90
      x319 = 2.0d0*x1
      x320 = x104*x319
      x321 = x318*x320
      x322 = x317*x321
      x323 = 0.713822582953479d0*x154
      x324 = x323*x57
      x325 = x14*x95
      x326 = 4.17757235705837d0*x325
      x327 = mu*x14
      x328 = x327*x91
      x329 = x1*x317
      x330 = x328*x329
      x331 = 16.7102894282335d0*x330
      x332 = 1.85745434535825d0*x21
      x333 = 7.42981738143299d0*x21
      x334 = x1*x94
      x335 = x317*x334
      x336 = x333*x335
      x337 = -x332 + x336
      x338 = -x326 + x331 + x337
      x339 = 9.02445731974005d0*x325
      x340 = x19*x92
      x341 = 13.5704102312052d0*x340
      x342 = 36.0978292789602d0*x330
      x343 = x101*x59
      x344 = x329*x343
      x345 = 54.281640924821d0*x344
      x346 = x337 - x339 - x341 + x342 + x345
      x347 = -x103*x338 + x346
      x348 = 1d0/x102
      x349 = mu*x348
      x350 = x318*x349
      x351 = x347*x350
      x352 = x293*x42
      x353 = 4.0d0*x44
      x354 = 2.56488140124205d0*x267
      x355 = x42*x51
      x356 = 0.0790071588102608d0*x355
      x357 = 0.000665486074881265d0*x253
      x358 = 0.58482233974552d0*x277
      x359 = 0.0506609980493537d0*x258
      x360 = x256*x359
      x361 = x1*x360
      x362 = x292 + x356 - x357 + x358 + x361
      x363 = -x279*x354 - x352 + x353 + x362
      x364 = x155*x363
      x365 = x111*(1.49912448908652d0*x14 - 0.550321208149104d0*x18)
      x366 = x113*x115
      x367 = x365*x366
      x368 = x112*x66
      x369 = 2.14502939711103d0*x0**(-0.666666666666667d0)
      x370 = x368*x369
      x371 = 1.12661476755593d0*x370
      x372 = 24.0d0*x114
      x373 = x372*x65
      x374 = 1d0/x119
      x375 = x126*x374
      x376 = x117*x375
      x377 = x1*x118
      x378 = -x377 + 1.0d0
      x379 = 1d0/x13
      x380 = x132*x379*x71
      x381 = 0.00706864485168613d0*x380
      x382 = x378*x381
      x383 = x376*x382
      x384 = 1d0/x138
      x385 = x144*x384
      x386 = x136*x385
      x387 = x1*x137
      x388 = x387 + 1.0d0
      x389 = x381*x388
      x390 = x386*x389
      x391 = 0.545136177849642d0*x134
      x392 = x130*x391
      x393 = x148*x391
      x394 = x125**(-2)
      x395 = x127*x394
      x396 = x122*x395
      x397 = x135*x396
      x398 = 0.034943218589452d0*x18
      x399 = x123*x398
      x400 = 0.164800693786261d0*x120
      x401 = x399 + x400
      x402 = -x377*x399 - x377*x400 + x401
      x403 = x117*x402
      x404 = x397*x403
      x405 = x143**(-2)
      x406 = x145*x405
      x407 = x141*x406
      x408 = x135*x407
      x409 = x142*x398
      x410 = 0.164800693786261d0*x139
      x411 = x409 + x410
      x412 = x387*x409 + x387*x410 + x411
      x413 = x136*x412
      x414 = x408*x413
      x415 = 2.18054471139857d0*x134
      x416 = x147*x3
      x417 = x415*x416
      x418 = x1*x417
      x419 = x129*x2
      x420 = x415*x419
      x421 = x1*x420
      x422 = x383 + x390 + x392 + x393 + x404 + x414 + x418 - x421
      x423 = x116 + x367 + x371 + x373 + x422
      x424 = 0.00111534778586481d0*x80
      x425 = 0.175432441092201d0*x107
      x426 = 0.122619224952946d0*x69
      x427 = 0.0399758348639607d0*x426
      x428 = 0.350864882184401d0*x108
      x429 = 0.00185891297644135d0*x149 - x428*x65
      x430 = x152*(x110 + x307*x425 - x423*x424 - x427*x78 + x429)
      x431 = 0.0334604335759443d0*x160
      x432 = x431*x65
      x433 = 0.0278836946466203d0*x161
      x434 = x14*x76
      x435 = 0.000844497709915773d0*x70
      x436 = 0.0047892375d0*x73
      x437 = 0.0225440733881905d0*x18
      x438 = 0.0031361987165122d0*x14
      x439 = -x435 + x436 - x437 - x438
      x440 = 0.113935126864452d0*x434 + x439
      x441 = x306*x440
      x442 = 0.0334604335759443d0*x159
      x443 = x441*x442
      x444 = 1.0525946465532d0*x167
      x445 = 4.0d0*x164
      x446 = -x445
      x447 = x289 - x291
      x448 = -x285 + x288 + x447
      x449 = x42*x448
      x450 = x163*x37
      x451 = x249*x450
      x452 = x354*x451
      x453 = x362 + x446 + x449 + x452
      x454 = 0.526297323276602d0*x166
      x455 = x453*x454
      x456 = x170*(x432 - x433 + x443 + x444 - x455)
      x457 = 0.263148661638301d0*x172
      x458 = x457*x65
      x459 = x441*x73
      x460 = 0.263148661638301d0*x459
      x461 = x188*x204
      x462 = 0.00232364122055169d0*x461
      x463 = x165*x70
      x464 = 2.06952542814068d0*x463
      x465 = x131*x175
      x466 = 0.545136177849642d0*x465
      x467 = x131*x182
      x468 = 7.26848237132856d0*x467
      x469 = x131*x16
      x470 = x179*x469
      x471 = 0.117507131669812d0*x14 - 0.372017136708795d0*x18
      x472 = x471*x66
      x473 = x470*x472
      x474 = 10.9027235569928d0*x473
      x475 = x12*x379
      x476 = x181*x475
      x477 = 1.80615420514536d0*x476
      x478 = x180*x469
      x479 = 21.8054471139857d0*x478
      x480 = x479*x65
      x481 = x3**1.66666666666667d0
      x482 = x2**1.66666666666667d0
      x483 = x134*(x481 - x482)
      x484 = x1*x483
      x485 = 2.18054471139857d0*x484
      x486 = -x383 - x390 - x392 - x393 - x404 - x414 - x418 + x421
      x487 = x466 - x468 - x474 - x477 - x480 + x485 + x486
      x488 = 0.00139418473233101d0*x159
      x489 = x487*x488
      x490 = x190*x453
      x491 = x192*(x174 + x458 + x460 + x462 + x464 - x489 - x490)
      x492 = x62*(x301 - x308 + x312 - x314 - x322 + x324 + x351 - x364 
     &- x430 + x456 + x491)
      x493 = 2.0d0*x63
      x494 = x3**(-0.333333333333333d0)
      x495 = x2**(-0.333333333333333d0)
      x496 = x494 - 1.0d0*x495
      x497 = rho_s*x496
      x498 = x493*x497
      x499 = x243/x242
      x500 = x499**5.56d0
      x501 = -x213 + 1.0d0
      x502 = x199*x223
      x503 = 0.0115661327577615d0*x226
      x504 = x499**2.78d0*x503
      x505 = x11/x501
      x506 = -x504*x505 + 1.0d0
      x507 = 1d0/(x500*x502/x501**2 + x506)
      x508 = x504*x506
      x509 = x507*x508
      x510 = x11*x509 + 1.0d0
      x511 = x37*log(x510)
      x512 = x5*x511
      x513 = x227*x234
      x514 = 1.46459188756152d0
      x515 = gamma_cc*x514
      x516 = x4*x515
      x517 = 1d0/x510
      x518 = rho_c**(-3.33333333333333d0)
      x519 = 0.0269876431014436d0*x518
      x520 = x232*x519
      x521 = 1d0/x4
      x522 = rho_c**(-4.33333333333333d0)
      x523 = x232*x522
      x524 = x521*x523
      x525 = -x495
      x526 = x494 + x525
      x527 = rho_s*x526
      x528 = 0.015421510343682d0*x527
      x529 = x524*x528
      x530 = 8.0d0*x32
      x531 = x31*x37*x530
      x532 = 0.0790071588102608d0*x54
      x533 = -x531 + x532
      x534 = -x195 - x353 + x533
      x535 = 1d0/x534
      x536 = x534*x58
      x537 = -x301
      x538 = x308 - x312 + x314 + x322 - x324 - x351 + x364 + x430 - x45
     &6 - x491 + x537
      x539 = x538*x62
      x540 = 4.0d0*x352
      x541 = 1.35631049481572d0*x296
      x542 = x194*x541
      x543 = 16.0d0*x44
      x544 = 10.2595256049682d0*x267
      x545 = x279*x544
      x546 = x264*x530
      x547 = 0.517758227584304d0*x247
      x548 = 0.316028635241043d0*x355
      x549 = 0.00266194429952506d0*x253
      x550 = 2.33928935898208d0*x277
      x551 = 0.202643992197415d0*x258
      x552 = x256*x551
      x553 = x1*x552
      x554 = -x546 + x547 + x548 - x549 + x550 + x553
      x555 = -x539 - x540 - x542 + x543 - x545 + x554
      x556 = -x453*x536 + x555
      x557 = x535*x556
      x558 = x232*x518
      x559 = 0.032153849066577d0*x558
      x560 = x557*x559
      x561 = -x164 + x207 + x241
      x562 = exp(x212*x561)
      x563 = x562 - 1.0d0
      x564 = 1d0/x563
      x565 = 0.00289153318944038d0*x564
      x566 = x224*x518
      x567 = 0.107950572405774d0*x566
      x568 = -x494 + x495
      x569 = x224*x522
      x570 = x521*x569
      x571 = rho_s*x570
      x572 = 0.0616860413747281d0*x568*x571
      x573 = 0.128615396266308d0*x566
      x574 = -x206 + x445 + x533
      x575 = 1d0/x574
      x576 = -x356 + x357 - x358 - x361 + x445 + x447 - x449 - x452
      x577 = 1d0/x165
      x578 = x574*x577
      x579 = x67*(-x302 - x303 - x304 + x74)
      x580 = x200*x579
      x581 = x201*x82
      x582 = -x315 + x316
      x583 = x202*x90
      x584 = x320*x583
      x585 = x332 - x336
      x586 = x326 - x331 + x585
      x587 = x103*x586 + x346
      x588 = x349*x583
      x589 = -x373
      x590 = -x116 - x367 - x371 + x486 + x589
      x591 = x152*(x201*x427 - x203 + x424*x590 + x425*x580 + x429) + x1
     &55*x453 - x165*x323 + x170*(-x432 + x433 - x443 - x444 + x455) + x
     &192*(x173 - x458 - x460 - x462 - x464 + x489 + x490) - x201*x313 +
     & x311* x581 + x537 + x580*x84 + x582*x584 + x587*x588
      x592 = x591*x62
      x593 = 4.0d0*x449
      x594 = x205*x541
      x595 = 16.0d0*x164
      x596 = x451*x544
      x597 = x554 - x592 + x593 - x594 - x595 + x596
      x598 = x576*x578 + x597
      x599 = x575*x598
      x600 = x573*x599
      x601 = 1.0d0*x494
      x602 = x495 - x601
      x603 = x206 + x446 + x531 - x532
      x604 = x521*x603
      x605 = x41 - 1.0d0
      x606 = x252*x605
      x607 = x38*x606
      x608 = x38*x605
      x609 = x276*x46
      x610 = 2.33928935898208d0*x609
      x611 = x51*x605
      x612 = x1*x256
      x613 = x546 - x547 - x548
      x614 = x592 - x593 + x594 + x595 - x596
      x615 = -0.00266194429952506d0*x607 + x608*x610 + 0.202643992197415
     &d0*x611*x612 + x613 + x614
      x616 = x319*x602*x604 + x615
      x617 = 142.664564503142d0
      x618 = x210*x9
      x619 = x617*x618
      x620 = 0.0925290620620922d0*x619
      x621 = x562*x564
      x622 = x620*x621
      x623 = x518*x616*x622
      x624 = x567 - x572 - x600 + x623
      x625 = x225*x624
      x626 = x208*x577
      x627 = x626**2.78d0
      x628 = x6*x627
      x629 = x10*x628
      x630 = x565*x625*x629
      x631 = x225*x232
      x632 = x198*x231
      x633 = rho_c**(-5.66666666666667d0)
      x634 = x217*x633
      x635 = 0.00249714130344307d0*x634
      x636 = 2062.20808228389d0*gamma_cc
      x637 = x636/x4**5
      x638 = 0.00214040683152263d0*x637
      x639 = x525 + x601
      x640 = x521*x574
      x641 = x549 - x550 - x553 + x613
      x642 = x539 + x540 + x542 - x543 + x545 + x641
      x643 = x319*x639*x640 + x642
      x644 = x216*x633
      x645 = x213*x643*x644
      x646 = x198*x230
      x647 = rho_c**(-6.66666666666667d0)
      x648 = x217*x647
      x649 = gamma_cc*x618
      x650 = x648*x649
      x651 = x527*x650
      x652 = x646*x651
      x653 = x632*x634
      x654 = x557*x653
      x655 = x227*x566
      x656 = x557*x655
      x657 = 0.0616860413747281d0*x526
      x658 = x227*x518
      x659 = x213*x230
      x660 = x643*x659
      x661 = x227*x567 + x227*x571*x657 + x620*x658*x660 - 0.12861539626
     &6308d0*x656
      x662 = x632*x635 + x638*x645 + 0.00142693788768175d0*x652 - 0.0029
     &7516549581645d0*x654 + x661
      x663 = x231*x662
      x664 = x234*x663
      x665 = 0.00289153318944038d0*x631*x664
      x666 = x517*(x520 + x529 - x560 + x630 - x665)
      x667 = x513*x516*x666
      x668 = 0.125d0*rho_c
      x669 = 0.0790071588102608d0*x51
      x670 = rho_s**3
      x671 = x38*x670
      x672 = x39*x671
      x673 = x669*x672
      x674 = x0*x360
      x675 = x256*x279
      x676 = 2.56488140124205d0*x40*x675
      x677 = 4.0d0*x270*x672
      x678 = x1*x77
      x679 = x300*x678
      x680 = 2.0d0*x104
      x681 = x317*x680
      x682 = 36.0978292789602d0*x325 + x333 + 54.281640924821d0*x340
      x683 = x317*x682
      x684 = 16.7102894282335d0*x325 + x333
      x685 = x103*x684
      x686 = x582*x685 + x683
      x687 = x37*x88
      x688 = x349*x687
      x689 = rho_c**(-3)
      x690 = x671*x689
      x691 = x669*x690
      x692 = -x691
      x693 = 2.56488140124205d0*x41
      x694 = x256*x270
      x695 = 4.0d0*x690
      x696 = -x270*x695 - x360 + x692 + x693*x694
      x697 = x1*x431
      x698 = x1*x428
      x699 = x1*x372
      x700 = -x699
      x701 = x2*x375
      x702 = x381*x701
      x703 = x3*x385
      x704 = x381*x703
      x705 = x396*x401
      x706 = x2*x705
      x707 = x135*x706
      x708 = x407*x411
      x709 = x3*x708
      x710 = x135*x709
      x711 = -x417 + x420 + x702 - x704 + x707 - x710
      x712 = x700 + x711
      x713 = x1*x457
      x714 = x134*(-x481 + x482)
      x715 = 2.18054471139857d0*x714
      x716 = x1*x479
      x717 = x711 - x715 - x716
      x718 = x152*(-x424*x712 + x698) + x155*x696 + x170*(x454*x696 + x6
     &97) + x192*( x190*x696 - x488*x717 + x713) - x318*x681 + x679 + x6
     &86*x688
      x719 = x62*x718
      x720 = 0.25d0*x0
      x721 = x719*x720
      x722 = x239*x5
      x723 = x639*x720*x722
      x724 = 1d0/x237
      x725 = x234*x724
      x726 = x516*x725
      x727 = x521*x526
      x728 = 0.015421510343682d0*x227
      x729 = -x254 + x255
      x730 = x359*x729
      x731 = x450*x729
      x732 = x693*x731
      x733 = x450*x695
      x734 = x692 + x730 + x732 + x733
      x735 = x51*x690
      x736 = 0.316028635241043d0*x735
      x737 = -x736
      x738 = x270*x41
      x739 = x256*x738
      x740 = 10.2595256049682d0*x739
      x741 = 16.0d0*x690
      x742 = x270*x741
      x743 = -x552 - x719 + x737 + x740 - x742
      x744 = -x536*x734 + x743
      x745 = x535*x744
      x746 = x227*x745
      x747 = x224*x521
      x748 = x657*x747
      x749 = 0.128615396266308d0*x224
      x750 = x496*x521
      x751 = 2.0d0*x750
      x752 = x736 + x742
      x753 = x552 + x719 - x740 + x752
      x754 = -x574*x751 + x753
      x755 = x659*x754
      x756 = -x620*x755 + x745*x749 + x748
      x757 = x632*x756
      x758 = x633*x757
      x759 = x513*x631
      x760 = x526*x649
      x761 = x634*x760
      x762 = x626**5.56d0
      x763 = x564*x762
      x764 = 0.00142693788768175d0*x763
      x765 = x761*x764
      x766 = x691 - x730 - x732 - x733
      x767 = x202*x88
      x768 = x349*x767
      x769 = x682 - x685
      x770 = x582*x769
      x771 = x417 - x420 - x702 + x704 - x707 + x710
      x772 = x699 + x771
      x773 = x715 + x716 + x771
      x774 = x152*(x424*x772 + x698) + x155*x734 + x170*(x454*x734 + x69
     &7) + x192*( x190*x734 + x488*x773 + x713) + x583*x681 + x679 + x76
     &8*x770
      x775 = x62*x774
      x776 = x551*x729
      x777 = x450*x741
      x778 = x41*x450
      x779 = 10.2595256049682d0*x778
      x780 = x729*x779
      x781 = x737 + x777 + x780
      x782 = -x775 + x776 + x781
      x783 = x578*x766 + x782
      x784 = x575*x783
      x785 = x11*x564
      x786 = x762*x785
      x787 = 0.00297516549581645d0*x634*x786
      x788 = x784*x787
      x789 = 0.202643992197415d0*x729
      x790 = x736 - x777 - x780
      x791 = x775 + x790
      x792 = x611*x789 + x791
      x793 = x603*x751 + x792
      x794 = x563**(-2)
      x795 = x762*x794
      x796 = x562*x795
      x797 = x793*x796
      x798 = x199*x633
      x799 = x638*x798
      x800 = x797*x799
      x801 = x566*x627
      x802 = x521*x657*x801
      x803 = x627*x784
      x804 = x573*x803
      x805 = x622*x793
      x806 = x518*x627
      x807 = x805*x806
      x808 = x802 + x804 - x807
      x809 = x765 + x788 - x800 + x808
      x810 = x11*x565
      x811 = x809*x810
      x812 = x558*x727*x728 + x559*x746 + 0.00289153318944038d0*x758 - x
     &759*x811
      x813 = 0.125d0*x812
      x814 = x726*x813
      x815 = x36*x4
      x816 = 0.000267550853940328d0*x218
      x817 = x507*(-x219*x500*x505 - x231*x509*(x229 + x632*x816) + x508
     &)
      x818 = x815*x817
      x819 = 0.258879113792152d0*x247
      x820 = 0.00133097214976253d0*x253
      x821 = 0.101321996098707d0*x258
      x822 = x256*x63
      x823 = 4.0d0*x265
      x824 = 5.1297628024841d0*x267
      x825 = 0.5d0*x0
      x826 = x497*x63
      x827 = x248*x269
      x828 = rho_c**(-7)
      x829 = x267*x828
      x830 = 9.72691162547591d0*x51
      x831 = x252*x52
      x832 = rho_s*x689
      x833 = x256*x832
      x834 = x2**(-0.666666666666667d0)
      x835 = 0.854960467080683d0*x1
      x836 = x834*x835
      x837 = x3**(-0.666666666666667d0)
      x838 = x835*x837
      x839 = x836 + x838
      x840 = -5.1297628024841d0*x254 + 5.1297628024841d0*x255 + x839
      x841 = 1d0/x514
      x842 = x264*x841
      x843 = x13*x842
      x844 = 3.10654936550583d0*x63
      x845 = x14*x264
      x846 = x263*(1.44868929752102d0*x14 + 0.452125560575033d0*x18 + 1.
     &06404211235975d0*x19 + 1.40779759833938d0*x21)
      x847 = x32*x846
      x848 = 48.0d0*x63
      x849 = x24*x260**2/(x261**2*x29**4)
      x850 = x32*x849
      x851 = x14*x276
      x852 = x53*x851
      x853 = x52*x609
      x854 = 492.457229038473d0*x270
      x855 = x275*(1.46302435997558d0*x14 + 0.455583412166236d0*x18 + 0.
     &571745641451468d0*x19 + 1.91958077938847d0*x21)
      x856 = x272*x855
      x857 = 14.0357361538925d0*x63
      x858 = x273**2/(x48**4*x49**2)
      x859 = x272*x858
      x860 = x260*x262*(2.17303394628153d0*x14 + 0.542550672690039d0*x18
     & + 2.41336731143894d0*x21 + x28)/x29**3
      x861 = x32*x860
      x862 = x266*x840
      x863 = 24.0d0*x270
      x864 = x273*x274*(2.19453653996337d0*x14 + 0.546700094599483d0*x18
     & + 3.29070990752309d0*x21 + x47)/x48**3
      x865 = x272*x864
      x866 = x284*(2.50265962182128d0*x14 + 0.573407182830959d0*x18 + 2.
     &18641103566438d0*x19 + 2.6168165555767d0*x21)
      x867 = 0.0829737643144773d0*x13*x287
      x868 = x281*x282
      x869 = x15*x868*(3.75398943273192d0*x14 + 0.688088619397151d0*x18 
     &+ x20 + 4.48597123813148d0*x21)/x22**3
      x870 = x283*x868
      x871 = x286*x841*x870
      x872 = 0.062230323235858d0*x14*x870
      x873 = x15*x24*x281**2/(x22**4*x25**2)
      x874 = x33*x846
      x875 = 0.172586075861435d0*x247
      x876 = x33*x860
      x877 = x290*x842
      x878 = 0.129439556896076d0*x845
      x879 = x33*x849
      x880 = x874 - x875 - x876 + x877 + x878 + x879
      x881 = -x866 + x867 + x869 - x871 - x872 - x873 + x880
      x882 = x191*x70/x61**6
      x883 = x194*x882
      x884 = x194*x296
      x885 = x296*x538
      x886 = x39*x64
      x887 = x496*x526
      x888 = x238*x8*x815
      x889 = 8.0d0*x887*x888
      x890 = x2**(-1.33333333333333d0)
      x891 = 0.333333333333333d0*x1
      x892 = x890*x891
      x893 = x3**(-1.33333333333333d0)
      x894 = x891*x893
      x895 = x892 + x894
      x896 = 0.211916079314314d0*x160*x65
      x897 = 0.074356519057654d0*x161
      x898 = x65*x67
      x899 = x440*x898
      x900 = 0.133841734303777d0*x159*x899
      x901 = 0.111534778586481d0*x159*x441
      x902 = 0.455740507457808d0*x14
      x903 = 0.00418159828868294d0*x14 + 0.0259624262672375d0*x18*x76 + 
     &0.0375734556469842d0*x18 - 0.151913502485936d0*x434 + 0.0019704946
     &5647014d0*x70 + x75
      x904 = x306*(x439*x902 + x903)
      x905 = x442*x904
      x906 = 3.15778393965961d0*x167
      x907 = 2.10518929310641d0*x166*x453
      x908 = -x874 + x875 + x876 - x877 - x878 - x879
      x909 = x42*(x866 - x867 - x869 + x871 + x872 + x873 + x908)
      x910 = 8.0d0*x449
      x911 = x249*x448
      x912 = x824*x911
      x913 = x451*x862
      x914 = 20.0d0*x164
      x915 = x268*x451
      x916 = x250*x267
      x917 = 0.40528798439483d0*x916
      x918 = 0.000887314766508353d0*x253
      x919 = 0.395035794051304d0*x355
      x920 = 0.0197517897025652d0*x258
      x921 = x1*x840
      x922 = x920*x921
      x923 = x252*x42
      x924 = 0.00532388859905012d0*x923
      x925 = 17.3157558993759d0*x859
      x926 = 0.58482233974552d0*x272
      x927 = x855*x926
      x928 = x864*x926
      x929 = x42*x609
      x930 = 4.67857871796416d0*x929
      x931 = 0.0394081882442864d0*x852
      x932 = x612*x831
      x933 = 0.00341378571249706d0*x932
      x934 = x612*x853
      x935 = 2.99999988448829d0*x934
      x936 = x880 + x917 + x918 - x919 + x922 - x924 - x925 - x927 + x92
     &8 + x930 - x931 - x933 + x935
      x937 = x909 - x910 + x912 + x913 + x914 - x915 + x936
      x938 = x454*x937
      x939 = 12.0d0*x169
      x940 = 0.00464728244110338d0*x159
      x941 = x42*x881
      x942 = x279*x862
      x943 = x268*x279 - x294*x824 + 8.0d0*x352 - 20.0d0*x44 + x936 - x9
     &41 - x942
      x944 = 4.13905085628136d0*x70
      x945 = 0.299209d0*x181
      x946 = 0.181712059283214d0*x465
      x947 = 14.5369647426571d0*x473
      x948 = 10.9027235569928d0*x470*x66*(-0.156676175559749d0*x14 + 0.6
     &20028561181324d0*x18)
      x949 = 3.61230841029072d0*x475
      x950 = x179*x949
      x951 = x472*x950
      x952 = 2.22044604925031d-16*x476
      x953 = 2.42282745710952d0*x467
      x954 = x470*x471
      x955 = 43.6108942279714d0*x65*x954
      x956 = 7.22461682058144d0*x180*x475*x65
      x957 = 36.3424118566428d0*x478*x65
      x958 = 5.0d0*x1
      x959 = x86*x958
      x960 = x87*x958
      x961 = x959 + x960
      x962 = 0.726848237132856d0*x134
      x963 = x1*x962*(-6.0d0*x481 + 6.0d0*x482 + x961)
      x964 = 2.90739294853142d0*x484
      x965 = 0.00942485980224818d0*x380
      x966 = x378*x965
      x967 = x376*x966
      x968 = x388*x965
      x969 = x386*x968
      x970 = 1d0/x117
      x971 = 0.0115260005177981d0*x65
      x972 = 0.817704266774463d0*x380
      x973 = x376*x972*(0.0230520010355962d0*x377 - x970*x971 - 0.011526
     &0005177981d0)
      x974 = 0.181712059283214d0*x134
      x975 = x130*x974
      x976 = x148*x974
      x977 = 1d0/x136
      x978 = x386*x972*(0.0230520010355962d0*x387 + x971*x977 + 0.011526
     &0005177981d0)
      x979 = 0.0141372897033723d0*x380
      x980 = x374*x394
      x981 = x402*x980
      x982 = x117*x378*x979*x981
      x983 = x384*x405
      x984 = x979*x983
      x985 = x388*x413*x984
      x986 = 3.63424118566428d0*x134
      x987 = x129*x986
      x988 = x65*x987
      x989 = x147*x986
      x990 = x65*x989
      x991 = 1.09027235569928d0*x134
      x992 = x396*x403*x991
      x993 = x407*x413*x991
      x994 = 0.116477395298173d0*x18
      x995 = 0.0582386976490866d0*x18
      x996 = x123*x995
      x997 = x65*x970
      x998 = 0.219734258381682d0*x120
      x999 = -x996 - x998
      x1000 = x117*x397*(0.439468516763364d0*x120*x377 + x123*x377*x994 
     &- x996*x997 - x997*x998 + x999)
      x1001 = x142*x995
      x1002 = x65*x977
      x1003 = 0.219734258381682d0*x139
      x1004 = x1001 + x1003
      x1005 = x136*x408*(x1001*x1002 + x1002*x1003 + x1004 + 0.439468516
     &763364d0*x139* x387 + x142*x387*x994)
      x1006 = 0.0376994392089927d0*x1*x380
      x1007 = x1006*x388*x703
      x1008 = 1.45369647426571d0*x1*x134
      x1009 = x1008*x419
      x1010 = 0.0698864371789039d0*x18
      x1011 = x1010*x123
      x1012 = 0.329601387572523d0*x120
      x1013 = x1011 + x1012
      x1014 = x125**(-3)
      x1015 = x1014*x122*x127*x135
      x1016 = x1015*x403*(-x1011*x377 - x1012*x377 + x1013)
      x1017 = x1010*x142
      x1018 = 0.329601387572523d0*x139
      x1019 = x1017 + x1018
      x1020 = x143**(-3)
      x1021 = x1020*x135*x141*x145
      x1022 = x1021*x413*(x1017*x387 + x1018*x387 + x1019)
      x1023 = x378*x701
      x1024 = x1006*x1023
      x1025 = x1008*x416
      x1026 = x3*x407
      x1027 = 4.36108942279713d0*x134
      x1028 = x1*x1027
      x1029 = x1026*x1028*x412
      x1030 = x2*x396
      x1031 = x1028*x1030*x402
      x1032 = x1000 - x1005 + x1007 + x1009 + x1016 + x1022 - x1024 - x1
     &025 + x1029 - x1031 + x967 + x969 + x973 - x975 - x976 - x978 + x9
     &82 + x985 + x988 + x990 + x992 + x993
      x1033 = 0.263148661638301d0*x73
      x1034 = 1.3157433081915d0*x172*x65 - x457*x66 + 0.526297323276602d
     &0*x459 - 1.0525946465532d0*x73*x899
      x1035 = 12.0d0*x191
      x1036 = 0.0037178259528827d0*x80
      x1037 = x111*x366*(-1.99883265211535d0*x14 + 0.917202013581841d0*x
     &18)
      x1038 = x365*x66
      x1039 = 24.0d0*x113
      x1040 = x1038*x1039
      x1041 = 0.154912426780983d0*x0**(-0.333333333333333d0)*x368
      x1042 = 2.25322953511185d0*x369
      x1043 = x1038*x1042
      x1044 = 0.751076511703951d0*x370
      x1045 = 48.0d0*x113*x365*x65
      x1046 = 4.50645907022371d0*x112*x369*x65
      x1047 = 0.133252782879869d0*x426
      x1048 = 0.0799516697279215d0*x426
      x1049 = 0.00910930363347549d0*x80
      x1050 = 0.00836319657736588d0*x14 + 0.0751469112939684d0*x18 + 0.0
     &0394098931294027d0*x70 - 0.01915695d0*x73
      x1051 = x1050*x425
      x1052 = 0.350864882184401d0*x107
      x1053 = 0.701729764368802d0*x107
      x1054 = x305*x898
      x1055 = 1.75432441092201d0*x108*x65 - 0.00495710127051027d0*x149 -
     & 0.159903339455843d0*x299*x426
      x1056 = 12.0d0*x151
      x1057 = 0.803050405822664d0*x154
      x1058 = 7.85204841248826d0*x154
      x1059 = 4.28293549772087d0*x154
      x1060 = 5.57009647607783d0*x325
      x1061 = x317**2
      x1062 = x1061*x65
      x1063 = x100*x327
      x1064 = 33.420578856467d0*x1062*x1063
      x1065 = 6.0d0*x315
      x1066 = 6.0d0*x316
      x1067 = x2**(-1.33333333333333d0)
      x1068 = x1*x1067
      x1069 = x3**(-1.33333333333333d0)
      x1070 = x1*x1069
      x1071 = x1068 + x1070
      x1072 = -x1065 + x1066 + x1071
      x1073 = x1*x1072
      x1074 = x1073*x328
      x1075 = 5.57009647607783d0*x1074
      x1076 = 11.1401929521557d0*x330
      x1077 = 2.16703006958462d0*x21
      x1078 = 2.476605793811d0*x21
      x1079 = x1072*x1078*x334
      x1080 = 9.90642317524398d0*x21
      x1081 = x1061*x1080
      x1082 = x1081*x65*x91
      x1083 = x1078*x335
      x1084 = x1077 + x1079 + x1082 - x1083
      x1085 = x338*x95
      x1086 = 2.0d0*x98
      x1087 = x1086*x346
      x1088 = 3.71490869071649d0*x21
      x1089 = 8.35514471411675d0*x325
      x1090 = 14.859634762866d0*x21
      x1091 = x1090*x335
      x1092 = 33.420578856467d0*x330
      x1093 = x102/x97**2
      x1094 = x1085*x1093
      x1095 = 12.0326097596534d0*x325
      x1096 = 72.1956585579204d0*x1061
      x1097 = 144.751042466189d0*x1062*x19*x59/x88**4 + x1063*x1096*x65 
     &+ 18.0938803082737d0*x1073*x343 + 12.0326097596534d0*x1074 + x1084
     & + x1095 - 24.0652195193068d0*x330 + 20.3556153468079d0*x340 - 54.
     &281640924821d0*x344
      x1098 = 8.0d0*x104
      x1099 = x1061*x1098
      x1100 = x1099*x687
      x1101 = 3.0d0*x59
      x1102 = x1101/x102**2
      x1103 = x1102*x347
      x1104 = x310*x82
      x1105 = 0.121994199821273d0*x1104
      x1106 = 0.284653132916304d0*x310
      x1107 = 0.267683468607554d0*x1050*x83
      x1108 = 0.892278228691848d0*x83
      x1109 = 0.000542507213303895d0*x0**2.33333333333333d0
      x1110 = 1.18970430492246d0*x83
      x1111 = x1101*x348
      x1112 = x1111*x98
      x1113 = x1112*x338
      x1114 = 1.07073387443022d0*x83
      x1115 = x347*x688
      x1116 = -0.243988399642546d0*x1104*x299 + 3.39065726902902d0*x299*
     &x83
      x1117 = x62*(x1035*(-x1033*x306*(-x902*(x435 - x436 + x437 + x438)
     & + x903) + x1034 - 0.00619637658813783d0*x187*x188 + x190*x943 - x
     &363*x944 - x488*(x1032 + x945 + x946 + x947 + x948 + x951 - x952 -
     & x953 + x955 + x956 - x957 - x963 - x964) + 4.82889266566159d0*x57
     &*x70 + x940*(x422 - x466 + x468 + x474 + x477 + x480 - x485)) + x1
     &054* x1114 + x1056*(x1036*x423 + x1047*x78 + x1048*x307 - x1049*x7
     &8 - x1051*x306 - x1052*x307 + x1053*x1054 + x1055 - x424*(x1032 + 
     &x1037 + x1040 + x1041 + x1043 + x1044 + x1045 + x1046 + x589) - x4
     &28*x66) + x1057*x943 + x1058*x57 - x1059*x363 + x1072*x321 - x1100
     &*x65 + x1103*x346*x687 + x1105*x307 + x1106*x309 - x1107* x306 - x
     &1108*x307 - x1109*x309 - x1110*x78 - x1113*x347*x687 + 12.0d0*x111
     &5*x329 + x1116 - 3.0d0*x350*(-x103*(x1060 + x1064 + x1075 - x1076 
     &+ x1084) - x1085*x1087 + x1094*(-x1088 - x1089 + x1091 + x1092) + 
     &x1097) - x939*(-x896 + x897 + x900 - x901 + x905 - x906 + x907 - x
     &938))
      x1118 = -x520 - x529 + x560 - x630 + x665
      x1119 = x227*x515
      x1120 = x1118*x1119
      x1121 = x725*x826
      x1122 = x226*x627
      x1123 = 0.0115661327577615d0*x1122
      x1124 = x1123*x785 + 1.0d0
      x1125 = x1124*x627
      x1126 = rho_s*x568
      x1127 = x1125*x518
      x1128 = 0.032153849066577d0*x1127
      x1129 = x762*x810
      x1130 = x616*x796
      x1131 = x1124 + x502*x795
      x1132 = 1d0/x1131
      x1133 = x1125*x1132
      x1134 = x1133*x225*x810
      x1135 = x521*x9
      x1136 = x1132*x1135
      x1137 = 21.170591578193d0*x1136*x220/x237**2
      x1138 = 3.0d0*x1118*x513
      x1139 = x227*x523
      x1140 = x227*x232
      x1141 = rho_c**(-5.33333333333333d0)
      x1142 = x1141*x521
      x1143 = x1142*x527
      x1144 = x1140*x1143
      x1145 = rho_c**(-6.33333333333333d0)*x64
      x1146 = x526**2
      x1147 = x1146*x6
      x1148 = 0.0308430206873641d0*x1147
      x1149 = x1140*x1148
      x1150 = 2.66666666666667d0*x494
      x1151 = 2.66666666666667d0*x495
      x1152 = 0.444444444444444d0*x1
      x1153 = x1152*x890
      x1154 = x1152*x893
      x1155 = x1153 + x1154
      x1156 = rho_s*(x1150 - x1151 + x1155)
      x1157 = x1139*x557
      x1158 = x556**2
      x1159 = x534**(-2)
      x1160 = x1139*x1159
      x1161 = 0.0893877004050842d0*x1160
      x1162 = 0.032153849066577d0*x58
      x1163 = x1157*x1162
      x1164 = 0.032153849066577d0*x1160
      x1165 = x1164*x556
      x1166 = x564*x624
      x1167 = x1166*x658
      x1168 = x1167*x629
      x1169 = x527*x649
      x1170 = 0.00771075517184102d0*x1169
      x1171 = 0.0160769245332885d0*x557
      x1172 = 8.34d0*x536
      x1173 = x555*x58
      x1174 = 0.158014317620522d0*x355
      x1175 = x612*x821
      x1176 = 1.16964467949104d0*x277
      x1177 = x451*x824
      x1178 = 8.0d0*x164
      x1179 = 2.0d0*x449
      x1180 = x534/x57**2
      x1181 = 8.34d0*x1180
      x1182 = x267*x294
      x1183 = x267*x279
      x1184 = x258*x921
      x1185 = 0.658919704477575d0*x1184 - 5.75747149073746d0*x247 + 0.02
     &96008206107187d0*x253 - 13.1783940895515d0*x355 + 4.3181036180531d
     &0*x843 + 4.3181036180531d0*x845 + 66.72d0*x847 + 66.72d0*x850 - 1.
     &3146571598294d0*x852 - 19.5096732539106d0*x856 - 577.653616803179d
     &0*x859 - 66.72d0*x861 + 19.5096732539106d0*x865 + 13.5204071594115
     &d0*x916 - 0.177604923664312d0*x923 + 156.077386031284d0*x929 - 0.1
     &13883891368902d0*x932 + 100.079996146529d0*x934
      x1186 = -2.78d0*x1117 - x1172*x937 - 16.68d0*x1173*x453 + x1181*x4
     &53*(x1174 + x1175 + x1176 + x1177 - x1178 + x1179 + x819 - x820 - 
     &x823) - 171.128887090869d0*x1182 + 684.515548363478d0*x1183 + x118
     &5 + 266.88d0*x352 - 667.2d0*x44 - 19.1776023007703d0*x883 + 18.852
     &7158779386d0*x884 - 22.6232590535263d0*x885 - 33.36d0*x941 - 33.36
     &d0*x942
      x1187 = 0.00385537758592051d0*x1139
      x1188 = x1187*x535
      x1189 = x513*x558
      x1190 = x1189*x663
      x1191 = x234*x662
      x1192 = x230*x662
      x1193 = x513*x625
      x1194 = 0.00144576659472019d0*x222
      x1195 = x199*x564
      x1196 = x526*x570
      x1197 = rho_s*x1196
      x1198 = x619*x659
      x1199 = x1198*x643
      x1200 = 0.00428081366304526d0*x637
      x1201 = x1200*x645 + 0.0028538757753635d0*x652 + 0.004994282606886
     &13d0*x653 - 0.0059503309916329d0*x654
      x1202 = x233**(-2)
      x1203 = x1140*x1202
      x1204 = x1203*x662
      x1205 = 0.000722883297360095d0*x222
      x1206 = x1205*x225
      x1207 = x1204*x1206
      x1208 = 4.31802289623097d0*x569
      x1209 = x1141*x747
      x1210 = 3.45441831698478d0*x1209
      x1211 = 1.48046499299348d0*x1147
      x1212 = x1211*x224
      x1213 = 0.555174372372553d0*x1156*x1209
      x1214 = 7.20246219091326d0*x569
      x1215 = 4.29060961944404d0*x1159
      x1216 = x1215*x569
      x1217 = x1209*x527
      x1218 = 5.18162747547716d0*x522
      x1219 = 16.0d0*x211
      x1220 = 4.44139497898043d0*x215*exp(-x1219*x208)
      x1221 = x4**(-6)
      x1222 = 1408.04281369971d0*x1221
      x1223 = x1222*x522
      x1224 = x1223*x643**2
      x1225 = 2.22069748949021d0*x199
      x1226 = x1225*x659
      x1227 = x1141*x527
      x1228 = x221*x617*x9
      x1229 = x1227*x1228
      x1230 = 1.5433847551957d0*x58
      x1231 = x557*x569
      x1232 = 1.5433847551957d0*x1159
      x1233 = x556*x569
      x1234 = x1199*x522
      x1235 = x1*x639
      x1236 = 12.0d0*x521
      x1237 = 6.0d0*x521
      x1238 = x1237*x574
      x1239 = x1*(2.0d0*x494 - 2.0d0*x495 + x895)
      x1240 = 16.0d0*x6
      x1241 = x1240*x65
      x1242 = x526*x639
      x1243 = 0.237021476430782d0*x258
      x1244 = 0.0409654285499647d0*x831
      x1245 = 35.9999986138594d0*x853
      x1246 = 2.07103291033722d0*x247 + 4.74042952861565d0*x355 - 1.5532
     &7468275291d0* x843 - 1.55327468275291d0*x845 - 24.0d0*x847 - 24.0d
     &0*x850 + 24.0d0*x861 - 4.86345581273796d0*x916 + 0.063886663188601
     &5d0*x923 - 56.14294461557d0*x929
      x1247 = 0.370116248248369d0*x1198
      x1248 = x1247*x522
      x1249 = 0.185058124124184d0*x535
      x1250 = 0.000240961099120032d0*x632
      x1251 = x1250*x225
      x1252 = x568**2
      x1253 = 1.48046499299348d0*x224*x628
      x1254 = x648*x786
      x1255 = rho_c**(-7.66666666666667d0)
      x1256 = x217*x763
      x1257 = x1255*x1256
      x1258 = x1257*x649
      x1259 = x1126*x1258
      x1260 = 0.0570775155072701d0*x10*x221
      x1261 = 0.888888888888889d0*x1
      x1262 = x1261*x890
      x1263 = x1261*x893
      x1264 = x599*x627
      x1265 = x574**(-2)
      x1266 = x1265*x598**2
      x1267 = x569*x627
      x1268 = 4.29060961944404d0*x1267
      x1269 = x1126*x1209
      x1270 = x621*x627
      x1271 = x1270*x616
      x1272 = x647*x796
      x1273 = x199*x637
      x1274 = x1272*x1273
      x1275 = x1274*x616
      x1276 = x616**2
      x1277 = x1223*x1276
      x1278 = exp(x1219*x561)
      x1279 = x199*x794
      x1280 = x1279*x627
      x1281 = 4.44139497898043d0*x1278*x1280
      x1282 = x1225*x1270
      x1283 = 20353.1779648712d0*gamma_cc/(x4**8*x8**3)
      x1284 = x1276*x1283
      x1285 = 0.154109291869629d0*x1278*x647*x762/x563**3
      x1286 = 0.0513697639565431d0*x1272
      x1287 = rho_s*x1141
      x1288 = x1228*x1287*x568
      x1289 = x1221*x199*x636
      x1290 = x1126*x1255*x1289
      x1291 = 1.5433847551957d0*x1267
      x1292 = x577*x599
      x1293 = x1291*x1292
      x1294 = 0.198503041880874d0*x1254
      x1295 = 0.0357019859497974d0*x1254
      x1296 = x1292*x1295
      x1297 = x1265*x598
      x1298 = x1297*(x614 + x641)
      x1299 = x616*x621
      x1300 = x522*x619
      x1301 = 6.17353902078279d0*x1300
      x1302 = -x1000 + x1005 - x1007 - x1009 - x1016 - x1022 + x1024 + x
     &1025 - x1029 + x1031 - x967 - x969 - x973 + x975 + x976 + x978 - x
     &982 - x985 - x988 - x990 - x992 - x993
      x1303 = x200*x67
      x1304 = x579*x65
      x1305 = x586*x95
      x1306 = x1093*x1305
      x1307 = -x1068 - x1070
      x1308 = x1102*x587
      x1309 = x582**2
      x1310 = x1098*x767
      x1311 = x587*x768
      x1312 = x62*(12.0d0*x1*x1311*x582 + x1035*(-x1033*x904 + x1034 + x
     &190*x937 - x453*x944 + 0.00619637658813783d0*x461 + 4.828892665661
     &59d0*x463 - x487*x940 + x488*(x1302 - x945 - x946 - x947 - x948 - 
     &x951 + x952 + x953 - x955 - x956 + x957 + x963 + x964)) + x1056*(-
     &x1036* x590 - x1047*x201 + x1048*x580 + x1049*x201 + x1051*x1303 -
     & x1052 *x580 - x1053*x1304 + x1055 + x200*x428 + x424*(-x1037 - x1
     &040 - x1041 - x1043 - x1044 - x1045 - x1046 + x1302 + x373)) + x10
     &57* x937 + x1058*x165 - x1059*x453 + x1105*x580 - x1106*x581 + x11
     &07* x1303 - x1108*x580 + x1109*x581 + x1110*x201 + x1113*x587*x767
     & - x1114*x1304 + x1116 + x1308*x767*(x339 + x341 - x342 - x345 + x
     &585) + x1309*x1310*x65 + x584*(x1065 - x1066 + x1307) + 3.0d0* x58
     &8*(x103*(-x1060 - x1064 - x1075 + x1076 - x1077 - x1079 - x1082 + 
     &x1083) + x1087*x1305 + x1097 + x1306*(x1088 + x1089 - x1091 - x109
     &2)) + x939*(x896 - x897 - x900 + x901 - x905 + x906 - x907 + x938)
     &)
      x1313 = x267*x911
      x1314 = x1237*x603
      x1315 = x1*x602
      x1316 = x205*x882
      x1317 = x296*x591
      x1318 = x205*x296
      x1319 = 0.237021476430782d0*x611
      x1320 = x267*x451
      x1321 = x46*x608
      x1322 = 7.01786807694624d0*x1321
      x1323 = x605*x609
      x1324 = x1236*x1315*x615 + x1239*x1314 + x1241*x568*x602*x603 + x1
     &246 + x1312 - 61.5571536298092d0*x1313 + 6.89841809380227d0*x1316 
     &+ 8.13786296889434d0*x1317 - 6.78155247407862d0*x1318 + x1319*x921
     & + 246.228614519237d0*x1320 - 207.78907079251d0*x1321*x858 - x1322
     & *x855 + x1322*x864 + 35.9999986138594d0*x1323*x612 - 240.0d0*x164
     & + 96.0d0*x449 - 0.0409654285499647d0*x606*x612 + 0.01064777719810
     &02d0*x607 - 0.472898258931437d0*x608*x851 - 12.0d0*x909 - 12.0d0*x
     &913
      x1325 = 0.370116248248369d0*x1270*x1300
      x1326 = 0.00856162732609051d0*x1274
      x1327 = x578*(x908 - x909 + x910 - x912 - x913 - x914 + x915 - x91
     &7 - x918 + x919 - x922 + x924 + x925 + x927 - x928 - x930 + x931 +
     & x933 - x935)
      x1328 = x577*x597
      x1329 = 16.68d0*x576
      x1330 = x574/x165**2
      x1331 = x1330*(-x1174 - x1175 - x1176 - x1177 + x1178 - x1179 - x8
     &19 + x820 + x823)
      x1332 = 0.185058124124184d0*x1267*x575
      x1333 = 0.00214040683152263d0*x1254*x575
      x1334 = 0.000240961099120032d0*x759*x785
      x1335 = 3.0d0*x726
      x1336 = 0.0416666666666667d0*rho_c
      x1337 = x249*x671
      x1338 = x269*x40
      x1339 = x1338*x256
      x1340 = 2.56488140124205d0*x255
      x1341 = 2.56488140124205d0*x254
      x1342 = x1340 - x1341 + x839
      x1343 = x256*x293
      x1344 = x63*x722
      x1345 = 3.0d0*x481
      x1346 = 3.0d0*x482
      x1347 = 21.8054471139857d0*x1*x954
      x1348 = x1*(-x177 + x178)
      x1349 = 0.726848237132856d0*x134
      x1350 = x140 - 1.0d0
      x1351 = x1350*x3
      x1352 = x121 - 1.0d0
      x1353 = x1352*x2
      x1354 = x1353*x395
      x1355 = x402*x415
      x1356 = x1001*x387 + x1003*x387 + x1004
      x1357 = x1351*x406
      x1358 = x1*x986
      x1359 = x412*x415
      x1360 = x377*x996
      x1361 = x377*x998
      x1362 = x996 + x998
      x1363 = x1*x415
      x1364 = 0.00471242990112409d0*x380
      x1365 = x1364*x703
      x1366 = x1364*x701
      x1367 = x703*x968
      x1368 = 0.0188497196044964d0*x380
      x1369 = x1*x1368
      x1370 = x1369*x375
      x1371 = x1369*x385
      x1372 = x3*x983
      x1373 = x1372*x381*x412
      x1374 = x2*x381*x981
      x1375 = x1372*x389*x411
      x1376 = x401*x980
      x1377 = x1376*x2*x382
      x1378 = x1365 - x1366 + x1367 + x1370 + x1371 + x1373 - x1374 + x1
     &375 - x1377
      x1379 = 0.00232364122055169d0*x159
      x1380 = 0.00170689285624853d0*x831
      x1381 = x41*x51
      x1382 = x1342*x738
      x1383 = 1.49999994224414d0*x853
      x1384 = x252*x690
      x1385 = 0.00266194429952506d0*x1384
      x1386 = x610*x690
      x1387 = x1385 - x1386
      x1388 = -x1342*x920 + x1343*x693 + x1380*x256 - x1381*x257 + x1382
     & - x1383*x256 + x1387 - x293*x695 - 20.5190512099364d0*x739 + x752
     &
      x1389 = 2.06952542814068d0*x70
      x1390 = x1*x67
      x1391 = x1390*x440
      x1392 = -0.789445984914903d0*x1*x172 + 0.526297323276602d0*x1391*x
     &73
      x1393 = 1.0525946465532d0*x166
      x1394 = -0.122688256445129d0*x1*x160 + 0.0669208671518886d0*x1391*
     &x159
      x1395 = 0.00185891297644135d0*x80
      x1396 = x1*x1042*x112
      x1397 = x1*x1039*x365
      x1398 = x701*x965*(x377 - 1.0d0)
      x1399 = x1023*x1368
      x1400 = x1*x987
      x1401 = x1349*x419
      x1402 = x1*x989
      x1403 = x1349*x416
      x1404 = x1363*x705
      x1405 = x391*x706
      x1406 = x1363*x708
      x1407 = x391*x709
      x1408 = x2*x397*(x1360 + x1361 + x999)
      x1409 = x1030*x1355
      x1410 = x1026*x1359
      x1411 = x1356*x3*x408
      x1412 = x1013*x1015
      x1413 = x1412*x2*x402
      x1414 = x1019*x1021
      x1415 = x1414*x3*x412
      x1416 = x1390*x305
      x1417 = -1.0525946465532d0*x1*x108 + x1048*x678
      x1418 = 2.14146774886044d0*x154
      x1419 = 3.0d0*x315
      x1420 = 3.0d0*x316
      x1421 = x1071 - x1419 + x1420
      x1422 = x1421*x680
      x1423 = 1.2383028969055d0*x21
      x1424 = 27.1408204624105d0*x340
      x1425 = x1081*x94
      x1426 = x1096*x328
      x1427 = 144.751042466189d0*x343
      x1428 = x1061*x1427
      x1429 = 12.0326097596534d0*x325
      x1430 = 18.0938803082737d0*x340
      x1431 = x582*x98
      x1432 = x1431*x684
      x1433 = -x1432*x346*x95
      x1434 = x1090 + 33.420578856467d0*x325
      x1435 = x1434*x582
      x1436 = x1307 + x1419 - x1420
      x1437 = x1436*x325
      x1438 = 5.57009647607783d0*x1437
      x1439 = x1060*x582
      x1440 = x1092*x582
      x1441 = x1078*x1436
      x1442 = x1423*x582
      x1443 = x1080*x335*x582
      x1444 = x1441 - x1442 + x1443
      x1445 = 3.0d0*x688
      x1446 = 6.0d0*x317
      x1447 = x582*x682
      x1448 = x1447*x37
      x1449 = x317*x37
      x1450 = x1449*x686
      x1451 = 6.0d0*x1*x349
      x1452 = x1112*x684
      x1453 = 0.535366937215109d0*x83
      x1454 = x1105*x678 - 1.96301210312207d0*x678*x83
      x1455 = x62*(x1*x1100 + x1035*(x1379*x717 + x1388*x190 - x1389*x69
     &6 + x1392 + x488*(x1013*x1014*x127*x135*x1353*x402 - x1019*x1020*x
     &135*x1351* x145*x412 - x128*x1349*x1353 - x128*x1352*x1358 + x1347
     & + 29.0739294853142d0*x1348*x470 - x1348*x950 + x1349*x1351*x146 -
     & x135*x1354*(-x1360 - x1361 + x1362) + x135*x1356*x1357 - x1350* x
     &1358*x146 - x1350*x1363*x406*x411 - x1352*x1363*x395*x401 + x1354*
     &x1355 + x1354*x391*x401 - x1357*x1359 - x1357*x391*x411 + x1378 - 
     &1.45369647426571d0*x483 - x701*x966 - x962*(-x1345 + x1346 + x961)
     &)) + x1056*(-x1052*x1416 + x1395*x712 + x1417 - x424 *(-x1365 + x1
     &366 - x1367 - x1370 - x1371 - x1373 + x1374 - x1375 + x1377 - x139
     &6 - x1397 + x1398 + x1399 - x1400 - x1401 - x1402 + x1403 - x1404 
     &+ x1405 - x1406 - x1407 + x1408 + x1409 - x1410 + x1411 + x1413 - 
     &x1415 + x699)) + x1057*x1388 + x1103*x1448 - x1115*x1446 - x1416*x
     &1453 - x1418*x696 - x1422*x318 - x1445*(-x1* x1425 - x1*x1426 - x1
     &*x1428 - x103*(x1438 - x1439 + x1440 + x1444 ) - x1078*x1421 + x10
     &85*x683*x98 + x1094*x1435 + x1095*x317 - x1421*x1429 - x1421*x1430
     & + x1423*x317 + x1424*x317 + x1433) + x1449*x1452*x347 - x1450*x14
     &51 + x1454 + x939*(x1388*x454 - x1393 *x696 + x1394))
      x1456 = x515*x812
      x1457 = x639*x725
      x1458 = 0.015421510343682d0*x727
      x1459 = x1137*(x1127*x1458 + x1128*x784 + x1129*x633*(x748 + x749*
     &x784 - x805) + x1134*(-x765 - x788 + x800 - x802 - x804 + x807))
      x1460 = 1.33333333333333d0*x494
      x1461 = 1.33333333333333d0*x495
      x1462 = x1155 + x1460 - x1461
      x1463 = x1139*x745
      x1464 = x647*x757
      x1465 = x646*x756
      x1466 = x556*x744
      x1467 = 0.00803846226664426d0*x557
      x1468 = 0.00385537758592051d0*x760
      x1469 = 0.00803846226664426d0*x746
      x1470 = x448*x695
      x1471 = -x1340 + x1341 - x836 - x838
      x1472 = x1471*x920
      x1473 = x448*x729
      x1474 = x1473*x693
      x1475 = x1381*x789
      x1476 = x1381*x256
      x1477 = 0.202643992197415d0*x1476
      x1478 = x1380*x729
      x1479 = x1471*x778
      x1480 = x256*x779
      x1481 = x1383*x729
      x1482 = x1387 + x1470 + x1472 + x1474 + x1475 - x1477 - x1478 + x1
     &479 + x1480 + x1481 + x790
      x1483 = x58*x743
      x1484 = 3.0d0*x453
      x1485 = 48.0d0*x690
      x1486 = x1485*x293
      x1487 = 0.158014317620522d0*x735
      x1488 = x729*x821
      x1489 = x729*x778
      x1490 = 5.1297628024841d0*x1489
      x1491 = x450*x690
      x1492 = 8.0d0*x1491
      x1493 = -x1487 + x1488 + x1490 + x1492
      x1494 = 30.7785768149046d0*x41
      x1495 = x1343*x1494
      x1496 = x1243*x1342
      x1497 = 0.0204827142749824d0*x831
      x1498 = x1497*x256
      x1499 = 4.06893148444717d0*x296
      x1500 = x1499*x718
      x1501 = 4.86345581273796d0*x1476
      x1502 = 12.0d0*x1382
      x1503 = 192.0d0*x270*x690
      x1504 = 17.9999993069297d0*x853
      x1505 = x1504*x256
      x1506 = 246.228614519237d0*x739
      x1507 = 3.79234362289252d0*x735
      x1508 = 0.0319433315943007d0*x1384
      x1509 = 28.071472307785d0*x609*x690
      x1510 = x1507 + x1508 - x1509
      x1511 = x535*(-3.0d0*x1173*x734 + x1180*x1484*x1493 - x1455 - 3.0d
     &0*x1482*x536 - x1483*x1484 - x1486 + x1495 - x1496 + x1498 - x1500
     & - x1501 + x1502 + x1503 - x1505 - x1506 + x1510)
      x1512 = x1189*x809
      x1513 = x1512*x785
      x1514 = x234*x809
      x1515 = x1514*x564
      x1516 = x621*x793
      x1517 = -x1200*x797*x798 + 0.0059503309916329d0*x634*x784*x786 + 0
     &.0028538757753635d0*x761*x763
      x1518 = -0.185058124124184d0*x1516*x619*x806 + x1517 + 0.123372082
     &749456d0*x727* x801 + 0.257230792532616d0*x784*x801
      x1519 = 1.72720915849239d0*x1196
      x1520 = 0.555174372372553d0*x570
      x1521 = 3.60123109545663d0*x569
      x1522 = 2.59081373773858d0*x1300
      x1523 = 2.05784634026093d0*x1196
      x1524 = x1228*x526
      x1525 = 1.48046499299348d0*x1524*x522
      x1526 = x1230*x696
      x1527 = x1232*x782
      x1528 = x1223*x643
      x1529 = 3.0867695103914d0*x1300
      x1530 = 6.0d0*x750
      x1531 = x1240*x526
      x1532 = -x1507 - x1508 + x1509
      x1533 = 0.514461585065233d0*x569
      x1534 = x648*x763
      x1535 = x1534*x760
      x1536 = x1534*x649
      x1537 = x1270*x793
      x1538 = x647*x797
      x1539 = x1273*x1538
      x1540 = x1254*x784
      x1541 = x1289*x526
      x1542 = -x776 + x791
      x1543 = x1291*x1542
      x1544 = x1295*x1542
      x1545 = x1223*x793
      x1546 = x1283*x616
      x1547 = x1265*x1268
      x1548 = x1378 - x1398 - x1399 + x1400 + x1401 + x1402 - x1403 + x1
     &404 - x1405 + x1406 + x1407 - x1408 - x1409 + x1410 - x1411 - x141
     &3 + x1415
      x1549 = x1*x579
      x1550 = x329*x582
      x1551 = 3.0d0*x768
      x1552 = x202*x769
      x1553 = x62*(x1035*(-x1379*x773 - x1389*x734 + x1392 + x1482*x190 
     &+ x488*(x1* x180*x949 - 29.0739294853142d0*x1*x478 + x1347 + x1548
     & + 1.45369647426571d0*x714 + x962*(x1345 - x1346 - x959 - x960))) 
     &+ x1056*(x1052*x1549 - x1395*x772 + x1417 + x424*(x1396 + x1397 + 
     &x1548 + x700)) + x1057*x1482 + x1111*x1432*x202*x587 - x1308* x144
     &7*x202 + x1309*x1451*x1552 + x1310*x1550 + x1311*x1446 - x1418*x73
     &4 + x1422*x583 + x1453*x1549 + x1454 + x1551*(x103*( -x1438 + x143
     &9 - x1440 - x1441 + x1442 - x1443) - x1095*x582 + x1305*x1431*x682
     & - x1306*x1435 - x1424*x582 + x1427*x1550 + x1430 *x1436 + x1433 +
     & 12.0326097596534d0*x1437 + x1444 + 72.1956585579204d0*x330*x582) 
     &+ x939*(-x1393*x734 + x1394 + x1482 *x454))
      x1554 = x1485*x448
      x1555 = x1473*x1494
      x1556 = x1499*x774
      x1557 = 2.43172790636898d0*x1476
      x1558 = 2.43172790636898d0*x1381*x729
      x1559 = 12.0d0*x1479
      x1560 = 192.0d0*x1491
      x1561 = 123.114307259618d0*x1489
      x1562 = 123.114307259618d0*x256*x778
      x1563 = x1237*x1315*x792 + x1314*(x602 - x892 - x894) + x1315*x153
     &1*x603 + x1319 *x1471 + 17.9999993069297d0*x1323*x729 + x1530*x615
     & + x1532 + x1553 - x1554 - x1555 + x1556 + x1557 - x1558 - x1559 +
     & x1560 + x1561 - x1562 - 0.0204827142749824d0*x606*x729
      x1564 = x577*x782
      x1565 = 3.0d0*x576
      x1566 = x1330*(x1487 - x1488 - x1490 - x1492)
      x1567 = x575*(x1243*x1471 + 3.0d0*x1328*x766 - x1497*x729 + x1504*
     &x729 + x1510 - x1553 + x1554 + x1555 - x1556 - x1557 + x1558 + x15
     &59 - x1560 - x1561 + x1562 + x1564*x1565 + x1565*x1566 + 3.0d0*x57
     &8*(-x1385 + x1386 - x1470 - x1472 - x1474 - x1475 + x1477 + x1478 
     &- x1479 - x1480 - x1481 + x781))
      x1568 = x257*x670
      x1569 = x38*x51
      x1570 = x834 + x837
      x1571 = x1570*x258
      x1572 = 0.0168869993497846d0*x1571
      x1573 = 12.0d0*x270
      x1574 = x1453*x77
      x1575 = x1067 + x1069
      x1576 = x1575*x680
      x1577 = 12.0d0*x349
      x1578 = x51*x689
      x1579 = x38*x65
      x1580 = x670*x689
      x1581 = 20.5190512099364d0*x1580
      x1582 = 0.854960467080683d0*x1570
      x1583 = x1569*x65
      x1584 = 0.237021476430782d0*x1583
      x1585 = x1572 - x1584
      x1586 = x1568*x1578 - x1573*x1579 + x1581*x694 - x1582*x738 + x158
     &5
      x1587 = x1078*x1575
      x1588 = mu*x684*x91
      x1589 = x1309*x328
      x1590 = x1309*x94
      x1591 = x1080*x1590 + x1587
      x1592 = x1309*x1588
      x1593 = -x103*(5.57009647607783d0*x1575*x325 + 33.420578856467d0*x
     &1589 + x1591) + x1093*x1434*x1592 + x1429*x1575 + x1430*x1575
      x1594 = 0.0282745794067445d0*x380
      x1595 = x1594*x375
      x1596 = x1594*x385
      x1597 = x1376*x979
      x1598 = x411*x984
      x1599 = x1027*x705
      x1600 = x1027*x708
      x1601 = x1362*x397
      x1602 = x1004*x408
      x1603 = x1412*x401
      x1604 = x1414*x411
      x1605 = x1595 + x1596 + x1597 + x1598 + x1599 + x1600 - x1601 - x1
     &602 + x1603 + x1604 + x987 + x989
      x1606 = x88*x986
      x1607 = x62*(x1035*(x1586*x190 + x457 - x488*(x1605 - x1606 - x479
     &)) + x1056*( -x424*(x1605 - x372) + x428) + x1057*x1586 - x1100 - 
     &x1102*x1448* x686*x94 - x1445*(x1086*x1588*x582*x683 + x1425 + x14
     &26 + x1428 + x1587 + x1593) - x1450*x1452*x94 + x1450*x1577 + x157
     &4 + x1576* x318 + x939*(x1586*x454 + x431))
      x1608 = x890 + x893
      x1609 = x744**2
      x1610 = 0.00771075517184102d0*x760
      x1611 = x1578*x670
      x1612 = x1611*x729
      x1613 = 0.40528798439483d0*x1612
      x1614 = x1579*x450
      x1615 = 12.0d0*x1614
      x1616 = x1581*x731
      x1617 = x1582*x778
      x1618 = x1585 - x1613 + x1615 + x1616 + x1617
      x1619 = x1611*x256
      x1620 = x1580*x694
      x1621 = x1579*x270
      x1622 = x1570*x738
      x1623 = 0.563350298308813d0*x1571 - 7.9070364537309d0*x1583
      x1624 = -x1172*x1618 + x1181*x1493*x734 - 16.68d0*x1483*x734 - 2.7
     &8d0*x1607 + 13.5204071594115d0*x1619 + 684.515548363478d0*x1620 - 
     &400.32d0* x1621 - 28.5214811818116d0*x1622 + x1623
      x1625 = x1465*x798
      x1626 = x1203*x1279*x809
      x1627 = 0.246744165498913d0*x1608
      x1628 = x1222*x754**2
      x1629 = 2.96092998598695d0*x1524
      x1630 = 2.84425771716939d0*x1583
      x1631 = 2.0d0*x1608
      x1632 = x1240*x887
      x1633 = 12.0d0*x750
      x1634 = x783**2
      x1635 = x577*x734
      x1636 = x793**2
      x1637 = x1223*x1636
      x1638 = x1283*x1636
      x1639 = x1265*x783
      x1640 = x1570*x778
      x1641 = x1580*x731
      x1642 = x1552*x1590
      x1643 = -x1595 - x1596 - x1597 - x1598 - x1599 - x1600 + x1601 + x
     &1602 - x1603 - x1604 - x987 - x989
      x1644 = x62*(x1035*(x1618*x190 + x457 + x488*(x1606 + x1643 + x479
     &)) + x1056*( x424*(x1643 + x372) + x428) + x1057*x1618 + x1099*x76
     &7 - x1102* x1642*x682 + x1452*x1642 + x1551*(-x1086*x1592*x682 + x
     &1309*x1427 + 72.1956585579204d0*x1589 + x1591 + x1593) + x1574 - x
     &1576*x583 + x1577*x202*x317*x770 + x939*(x1618*x454 + x431))
      x1645 = 0.202643992197415d0*x1570*x611 + 4.86345581273796d0*x1612 
     &- 144.0d0* x1614 + x1630 + x1631*x604 + x1632*x603 + x1633*x792 - 
     &10.2595256049682d0*x1640 - 246.228614519237d0*x1641 + x1644
      x1646 = x578*(-x1572 + x1584 + x1613 - x1615 - x1616 - x1617)
      x1647 = x1564*x766
      x1648 = x1566*x766
      x1649 = x762*x816
      x1650 = x1649*x785
      x1651 = x1123 + x1650
      x1652 = x1651*x785
      x1653 = -x1652*x236 + x219*x632 + x235
      x1654 = 0.000178367235960219d0*x646
      x1655 = x224*x227*x520
      x1656 = 0.032153849066577d0*x232
      x1657 = x1656*x656
      x1658 = x1652*x234
      x1659 = x1651*x564
      x1660 = x1122*x624
      x1661 = 0.00289153318944038d0*x513
      x1662 = x1651*x222
      x1663 = 3.34438567425411d-5*x218*x222
      x1664 = x226*x232
      x1665 = x1661*x1664
      x1666 = 0.0231322655155231d0*x1122 + 0.000535101707880657d0*x218*x
     &786
      x1667 = x1666*x222*x226
      x1668 = x1125*x503
      x1669 = x1132*x1651
      x1670 = -x1668*x1669*x785 + x1668 + x219*x786
      x1671 = x1136*x724
      x1672 = x232*x655
      x1673 = x1656*x655*x745
      x1674 = x1664*x513
      E = -rho_c*(-x239*x240 + x243)
      d1E(1) = x240*x511 + 0.25d0*x244*x62 + x57 + x668*(rho_s*x257*x259
     & + 1.03551645516861d0*x0*x247 - 0.00532388859905012d0*x0*x253 - 16
     &.0d0*x0*x265 + 0.632057270482087d0*x251 - x268*x271 + 4.6785787179
     &6416d0*x278 + 32.0d0*x280 - 8.0d0*x295 + 2.71262098963145d0*x297 +
     & x298*x492 + x498*x512 - x667)
      d1E(2) = -rho_c*(x673 + x674 - x676 + x677 + x721 + x723 - x814)
      d1E(3) = x517*x668*x818
      d2E(1) = rho_s*x821*x822 + x0*x819 - x0*x820 - x0*x823 - x1336*(x1
     &117*x493 - 12.0d0*x1120*x1121 + x1137*x1138*(0.015421510343682d0*x
     &1125*x1126 *x521*x522 - x1125*x519 + x1128*x599 + x1129*x225*(-x56
     &7 + x572 + x600 - x623) + x1134*(-x1126*x650*x764 + x1130*x799 + x
     &567*x627 - x572*x627 - x599*x787 - x600*x627 + x623*x627 + x635*x7
     &86)) - x1335*(-x1139*x1170*x1191*x230 + 0.0899588103381452d0*x1139
     & + x1142*x1156*x228*x232 - 0.0857435975108721d0*x1144*x557 + 0.071
     &9670482705162d0*x1144 + x1145*x1149 - 0.150051295644026d0* x1157 +
     & x1158*x1161 + x1163*x363 - x1165*x597 + x1166*x1170*x227* x522*x6
     &27 - x1168*x1171 + 0.0134938215507218d0*x1168 + x1171* x1190 + x11
     &86*x1188 - 0.0134938215507218d0*x1190 - x1192*x1193* x1194*x1195*x
     &627 + x1207*x215*(0.123372082749456d0*x1197*x227 + 0.1850581241241
     &84d0*x1199*x658 + x1201 + 0.215901144811548d0*x655 - 0.25723079253
     &2616d0*x656) + x1251*(x1145*x1212 + x1158*x1216 + x1186*x1249*x569
     & + x1199*x1218 + x1208 + x1210*x527 + x1213 - x1214*x557 - 4.11569
     &268052186d0*x1217*x557 + x1220*x1224 - x1224* x1226 + 2.9609299859
     &8695d0*x1229*x660 + x1230*x1231*x363 - x1232* x1233*x597 - 6.17353
     &902078279d0*x1234*x557 - x1248*(x1117 + 61.5571536298092d0*x1182 -
     & 246.228614519237d0*x1183 - x1235*x1236 *x642 - x1238*x1239 - x124
     &1*x1242*x574 - x1243*x921 + x1244*x612 - x1245*x612 + x1246 - 0.01
     &06477771981002d0*x253 - 96.0d0*x352 + 240.0d0*x44 + 0.472898258931
     &437d0*x852 + 7.01786807694624d0*x856 + 207.78907079251d0*x859 - 7.
     &01786807694624d0*x865 + 6.89841809380227d0*x883 - 6.78155247407862
     &d0*x884 + 8.13786296889434d0*x885 + 12.0d0*x941 + 12.0d0*x942)) - 
     &x1334*( rho_c**(-8.66666666666667d0)*x1252*x1256*x1260*x64 + 0.006
     &42122049456788d0*rho_s*x1258*(x1262 + x1263 + 5.33333333333333d0*x
     &494 - 5.33333333333333d0*x495) - x1126*x1210* x627 - 0.13698603721
     &7448d0*x1130*x1290 + x1145*x1252*x1253 + x1208*x627 + x1213*x627 -
     & x1214*x1264 + x1218*x1271*x619 - 0.333218535531443d0*x1254*x599 +
     & 0.169805608634128d0*x1254 + 0.190410591732253d0*x1259*x599 - 0.15
     &9817043420356d0*x1259 + 4.11569268052186d0*x1264*x1269 - x1264*x12
     &99*x1301 + x1266*x1268 + x1266*x1294 - 2.96092998598695d0*x1271*x1
     &288 - 0.285615887598379d0*x1275*x599 + 0.239725565130534d0*x1275 +
     & x1277*x1281 - x1277*x1282 + x1284*x1285 - x1284*x1286 + x1291* x1
     &298 + x1293*x453 + x1295*x1298 + x1296*x453 - x1324*x1325 - x1324*
     &x1326 + x1332*(x1185 - 2.78d0*x1312 + 171.128887090869d0* x1313 - 
     &19.1776023007703d0*x1316 - 22.6232590535263d0*x1317 + 18.852715877
     &9386d0*x1318 - 684.515548363478d0*x1320 + 8.34d0* x1327 + x1328*x1
     &329 + 8.34d0*x1331*x576 + 667.2d0*x164 - 266.88d0 *x449 + 33.36d0*
     &x909 + 33.36d0*x913) + x1333*(1.31783940895515d0* x1184 - 5.56d0*x
     &1312 + 342.257774181739d0*x1313 - 38.3552046015406d0*x1316 - 45.24
     &65181070525d0*x1317 + 37.7054317558771d0*x1318 - 1369.03109672696d
     &0*x1320 + 16.68d0* x1327 + 33.36d0*x1328*x576 + x1329*x1331 + 1334
     &.4d0*x164 - 11.5149429814749d0*x247 + 0.0592016412214373d0*x253 - 
     &26.356788179103d0*x355 - 533.76d0*x449 + 8.63620723610619d0*x843 +
     & 8.63620723610619d0*x845 + 133.44d0*x847 + 133.44d0*x850 - 2.62931
     &431965879d0*x852 - 39.0193465078211d0*x856 - 1155.30723360636d0*x8
     &59 - 133.44d0*x861 + 39.0193465078211d0*x865 + 66.72d0*x909 + 66.7
     &2d0*x913 + 27.040814318823d0*x916 - 0.355209847328624d0*x923 + 312
     &.154772062569d0*x929 - 0.227767782737804d0*x932 + 200.159992293059
     &d0*x934))) + 4.14206582067443d0*x247*x63 + 0.127773326377203d0*x25
     &2*x827 - 0.0212955543962005d0*x253*x63 - 0.474042952861565d0*x258*
     &x832* x840 + 480.0d0*x270*x827 - 192.0d0*x293*x827 + 123.114307259
     &618d0 *x293*x829 + 9.4808590572313d0*x51*x827 - 112.28588923114d0*
     &x609* x827 + 0.945796517862875d0*x63*x852 + 415.578141585021d0*x63
     &*x859 + 13.7968361876045d0*x63*x883 - 13.5631049481572d0*x63*x884 
     &+ 16.2757259377887d0*x63*x885 + 6.0d0*x722*x832*(2.0d0*x494 - 2.0d
     &0 *x495 + x895) + 24.0d0*x827*x881 + x828*x862*x863 - x829*x830 - 
     &x829*x854 + 0.0819308570999294d0*x831*x833 - 71.9999972277189d0* x
     &833*x853 - x843*x844 - x844*x845 - x847*x848 - x848*x850 + x848* x
     &861 + x856*x857 - x857*x865 - x886*x889) + 0.158014317620522d0* x2
     &51 - x271*x824 + 1.16964467949104d0*x278 + 8.0d0*x280 - 2.0d0* x29
     &5 + 0.678155247407862d0*x297 + x492*x825 + 0.5d0*x512*x826 - 0.25d
     &0*x667
      d2E(2) = x1336*(-6.0d0*x0*x1120*x1457 + 6.0d0*x1121*x1456 - x1138*
     &x1459 + x1244* x822 - x1245*x822 + x1335*(-0.00385537758592051d0*r
     &ho_c**( -7.66666666666667d0)*x1169*x1465 + 0.010717949688859d0*x11
     &39* x1511 - 0.0359835241352581d0*x1139*x727 - 0.042871798755436d0*
     & x1143*x232*x746 - x1149*x1287 + 0.042871798755436d0*x1157*x727 + 
     &x1161*x1466 + x1163*x696 - x1165*x782 - x1166*x1469*x518*x629 - x1
     &167*x1468*x627 + x1169*x1187*x1515 + x1189*x1192*x1468 + x1191* x1
     &205*x644*x756 + x1193*x1205*x1280*x809 - x1195*x1207*x1518*x230 + 
     &x1251*(-x1212*x1287 + x1216*x1466 - 2.05784634026093d0*x1217* x745
     & + x1220*x1528*x754 - x1225*x1528*x755 + 1.48046499299348d0* x1229
     &*x755 + x1231*x1526 - x1233*x1527 - 3.0867695103914d0*x1234* x745 
     &- x1248*(-x1235*x1237*x753 + x1235*x1531*x574 + x1238*(x639 + x895
     &) + x1455 + x1486 - x1495 + x1496 - x1498 + x1500 + x1501 - x1502 
     &- x1503 + x1505 + x1506 + x1530*x642 + x1532) - x1462*x1520 + x151
     &1*x1533 - x1519 - x1521*x745 + x1522*x755 + x1523*x557 - x1525*x66
     &0 - x1529*x557*x755) - x1334*(-x1225*x1271*x1545 + x1227 *x1253*x5
     &68 + 0.0119006619832658d0*x1254*x1567 + x1257*x1260*x527 *x568 + 0
     &.0952052958661265d0*x1259*x784 - x1264*x1516*x1529 + x1264*x1523 +
     & 2.05784634026093d0*x1269*x803 - x1271*x1525 - 0.0684930186087241d
     &0*x1272*x1541*x616 - 0.14280794379919d0*x1275* x784 + x1281*x1545*
     &x616 + x1285*x1546*x793 - 1.48046499299348d0* x1288*x1537 - 0.0684
     &930186087241d0*x1290*x797 + x1293*x734 + x1294*x1297*x783 + x1296*
     &x734 + x1297*x1543 + x1297*x1544 - x1299 *x1529*x803 - x1325*x1563
     & - x1326*x1563 - x1519*x627 + x1520*x627 *(-x1153 - x1154 - x1460 
     &+ x1461) - x1521*x803 + x1522*x1537 + x1533*x1567*x627 + 0.0952052
     &958661265d0*x1535*x599 - 0.0799085217101781d0*x1535 + 0.0064212204
     &9456788d0*x1536*(-x1150 + x1151 - x1262 - x1263) - 0.0513697639565
     &431d0*x1538*x1546 - 0.14280794379919d0*x1539*x599 + 0.119862782565
     &267d0*x1539 - 0.166609267765721d0*x1540 + x1547*x598*x783) - x1462
     &*x228*x524 - 0.0750256478220131d0*x1463 + x1464*x1467 - 0.00674691
     &077536089d0* x1464 - x1467*x1513 + x1469*x558*x664 + 0.00674691077
     &536089d0* x1513) + 0.0638866631886015d0*x1337*x252 - 56.1429446155
     &7d0*x1337 *x609 + x1338*x1342*x863 + 61.5571536298092d0*x1338*x134
     &3 - x1339 *x830 - x1339*x854 - 0.474042952861565d0*x1342*x259 + 6.
     &0d0*x1344 *(x496 + x895) - x1455*x493 + 7.58468724578504d0*x250*x6
     &71 + 384.0d0*x279*x671 - 96.0d0*x294*x671 - 8.13786296889434d0*x29
     &6* x63*x718 - x832*x889) - x673 - x674 + x676 - x677 - x721 - x723
     & + x814
      d2E(3) = -rho_c*(-0.333333333333333d0*x1242*x63*x888 + 0.854960467
     &080682d0*x1338* x1570*x270 + 0.0833333333333333d0*x1344*x1608 + x1
     &456*x1457*x825 + x1459*x234*x813 - x1568*x250 + 0.237021476430782d
     &0*x1569*x886 - x1572*x63 + x1573*x38*x886 + 0.0833333333333333d0*x
     &1607*x63 - 20.5190512099364d0*x670*x675 - 0.125d0*x726*(rho_c**( -
     &6.66666666666667d0)*x1250*(x1212 + x1215*x1609*x224 + x1220* x1628
     & - x1226*x1628 - x1247*(-0.202643992197415d0*x1571 + x1607 - 4.863
     &45581273796d0*x1619 - 246.228614519237d0*x1620 + 144.0d0* x1621 + 
     &10.2595256049682d0*x1622 + x1630 - x1631*x640 - x1632* x574 + x163
     &3*x753) + x1249*x1624*x224 + x1526*x224*x745 - x1527* x224*x744 + 
     &x1627*x747 - x1629*x755 + 4.11569268052186d0*x526* x745*x747 - 6.1
     &7353902078279d0*x619*x745*x755) + x1139*x1148 + 0.0051405034478940
     &1d0*x1139*x1608*x521 + x1161*x1609 + x1162* x1463*x696 - x1164*x74
     &4*x782 + x1188*x1624 - x1194*x1515*x1625 + x1206*x1518*x1626 - x13
     &34*(x1146*x1260*x1534 + 4.11569268052186d0 *x1196*x803 + x1211*x12
     &67 + x1265*x1294*x1634 + x1281*x1637 - x1282*x1637 + x1285*x1638 -
     & x1286*x1638 - x1301*x1516*x803 - x1325*x1645 - x1326*x1645 + x133
     &2*(-13.5204071594115d0*x1612 + 400.32d0*x1614 + x1623 + 28.5214811
     &818116d0*x1640 + 684.515548363478d0*x1641 - 2.78d0*x1644 + 8.34d0*
     &x1646 + 16.68d0* x1647 + 8.34d0*x1648) + x1333*(1.12670059661763d0
     &*x1571 - 15.8140729074618d0*x1583 - 27.040814318823d0*x1612 + 800.
     &64d0* x1614 + 57.0429623636231d0*x1640 + 1369.03109672696d0*x1641 
     &- 5.56d0*x1644 + 16.68d0*x1646 + 33.36d0*x1647 + 16.68d0*x1648) + 
     &0.190410591732253d0*x1535*x784 + 0.00570775155072701d0*x1536* x160
     &8 - x1537*x1629*x522 - 0.136986037217448d0*x1538*x1541 - 0.2856158
     &87598379d0*x1539*x784 + 0.0357019859497974d0*x1540*x1635 + x1543*x
     &1639 + x1544*x1639 + x1547*x1634 + x1627*x570*x627 + 1.54338475519
     &57d0*x1635*x569*x803) + 0.0857435975108721d0*x1463* x727 + 0.01607
     &69245332885d0*x1464*x745 + x1465*x1610*x647 - x1512 *x1610*x564 - 
     &0.0160769245332885d0*x1514*x558*x746*x785))
      d2E(4) = 0.125d0*x517*(rho_c*(x1119*x1135*x1202*x1653*x666 - x234*
     &x815*(-x1191* x1663*x216 + x1197*x232*x728 + x1204*x1667*x199*x230
     &*x565 - x1279 *x1660*x1661*x1662 + x1654*x651 - x1655*x1658 + x165
     &5 + x1657* x1658 - x1657 - x1659*x232*x513*x528*x569*x649 + 0.0057
     &8306637888076d0*x1660*x227*x785 - x1665*x231*(x1201 + x661 ) - x16
     &65*x663 + 0.000312142662930383d0*x653 - 0.000371895686977057d0*x65
     &4) + x36*x498*x817) + x818)
      d2E(5) = -x668*x725*(x1456*x1670*x1671 + x1653*x298*x36*x639 - x81
     &5*(x1458*x1672 - x1515*x1663*x199*x646 - x1625*x1662*x224*x234*x56
     &5 + 0.00289153318944038d0*x1626*x1667 + x1654*x761 - x1658*x1673 -
     & 0.015421510343682d0*x1659*x1672*x234*x760 + x1673 - x1674*x810*( 
     &x1517 + x808) - x1674*x811 + 0.00578306637888076d0*x224*x758 + 0.0
     &00371895686977057d0*x653*x745))
      d2E(6) = -x1671*x36*x668*(x1653*x1670*x725 + x564*(0.0093463840640
     &7122d0*rho_c**( -7.0d0)*x1124*x1132*x626**8.34d0*x785 + 0.02313226
     &55155231d0* x1133*x1651*x226 - x1649 + x1650*x1669 - x1652*x1666*x
     &1668/x1131 **2))
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, tau_c, mu
      real*8, intent(out) :: E
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = 0.682784063255296d0*x0**0.333333333333333d0
      x2 = 0.194159335344114d0*x1 + 1.0d0
      x3 = 2.14502939711103d0
      x4 = 1d0/x3
      x5 = x0**0.666666666666667d0
      x6 = x4*x5
      x7 = 0.564189583547756d0*sqrt(x0)
      x8 = 0.826307487110758d0*x0**0.166666666666667d0
      x9 = 9.86960440108936d0
      x10 = 0.693147180559945d0
      x11 = -x10 + 1.0d0
      x12 = 1d0/x11
      x13 = x12*x9
      x14 = log(0.5d0*x13/(3.2595509194222925d0*x1 + 0.40691300451752932
     &d0*x6 + 1.4187228164796675d0*x7 + 7.2401019343168302d0*x8) + 1.0d0
     &)
      x15 = x14*x2
      x16 = 2.0d0*x15
      x17 = exp(x16) - 1.0d0
      x18 = gamma_cc**2
      x19 = rho_c**(-4.66666666666667d0)*x18
      x20 = 208.945363813799d0*x19/x11**2
      x21 = gamma_cc*rho_c**(-2.33333333333333d0)
      x22 = 14.454942539277d0*x12*x21
      x23 = 0.00289153318944038d0*x22
      x24 = 1.0d0 + x23/x17
      x25 = log(x23*x24/(x24 + 8.3609641856352646d-6*x20/x17**2) + 1.0d0
     &)
      x26 = -x16 + 1.0d0*x25
      x27 = 1d0/x9
      x28 = x11*x27
      x29 = x26*x28
      x30 = (0.235215882559859d0*x1 + 1.0d0)*log(x13/(7.0945954610727853
     &d0*x1 + 0.81920444870878784d0*x6 + 4.1227361860783676d0*x7 + 15.10
     &5997393829993d0*x8) + 1.0d0)
      x31 = exp(2.0000000000000004d0*x30) - 1.0d0
      x32 = 0.00578306637888076d0*x22
      x33 = 1.0d0 + x32/x31
      x34 = -1.0d0*x30 + 0.5d0*log(x32*x33/(3.3443856742541065d-5*x20/x3
     &1**2 + x33) + 1.0d0)
      x35 = x28*x34
      x36 = x18/(rho_c**2*tau_c**2)
      x37 = x34*x36
      x38 = x26*(0.00828125d0*x36 + 1.0d0)
      x39 = -x38
      x40 = -2.8d0*x38
      x41 = gamma_cc**3/(rho_c**3*tau_c**3)
      x42 = x28*x41
      x43 = 1024.0d0*x15 - 512.0d0*x25
      x44 = 0.001953125d0*x28*exp(-1.0053738478084204d0*mu*rho_c**( -0.3
     &3333333333333331d0)*x4)
      x45 = x16*x28
      x46 = mu**2
      x47 = (0.508616435555896d0*x46*x6 + 1.0d0)**(-4)
      x48 = 1.41421356237310d0
      x49 = 4.60115111447049d0
      x50 = x0**1.33333333333333d0
      x51 = 3.14159265358979d0
      x52 = x0**1.0d0
      x53 = x52/x51
      x54 = (0.0188171922990732d0*x1 - 0.009578475d0*x53 + 0.06763222016
     &45715d0*x6 + 1.0d0 + 0.00126674656487366d0*x50/x49)*exp(-0.6836107
     &611867116d0* x1)
      x55 = x0**1.66666666666667d0
      x56 = 0.0837224526465878d0*x55
      x57 = 0.0223069557172962d0*mu**5*x48*x54*x56
      x58 = x1*x46
      x59 = 5.57236303607474d0*mu*x8 + 1.0d0
      x60 = mu**3
      x61 = 2.0d0*x27*(x10 - 1.0d0)*log((6.7683429898050367d0*x58 + x59 
     &+ 3.3926025578013106d0*x60*x7)/(3.1331792677937806d0*x58 + x59))
      x62 = x11*x15
      x63 = 0.000640558526213818d0*mu**8*x0**2.66666666666667d0*x62
      x64 = 0.148394183600699d0
      x65 = x54 - 1.0d0
      x66 = mu**6*(0.0108059179628909d0*x0**2.0d0*x62 + 0.01673021678797
     &22d0*x55*x64 *x65)
      x67 = 2.28942848510666d0*x49
      x68 = 1.63540853354893d0*x67*(-0.0259335011650457d0*x1 + 1.0d0)/( 
     &0.494402081358784d0*x1 + 0.0524148278841779d0*x6 + 1.0d0)
      x69 = x48*x60*(0.0315054072231411d0*x52*x54 + 0.00111534778586481d
     &0*x56*( 12.0d0*x51*(-4.49737346725955d0*x1 + 0.825481812223657d0*x
     &6)*exp( -0.28165369188898165d0*x1)/x52 + x68/x5))
      x70 = mu**4*(0.00139418473233101d0*x0**1.0d0*x64*(-10.902723556992
     &8d0*x3*( -0.352521395009435d0*x1 + 0.558025705063192d0*x6)*exp( -0
     &.49698248213959023d0*x1) + 1.63540853354893d0*x67 - x68) - 0.06835
     &90130188305d0*x50*x62 - 0.13157433081915d0*x53*x65)
      x71 = (x45 + x47*(-x57 + x61 - x63 - x66 - x69 + x70))/(x14*x2)
      x72 = x45 - x47*(x57 - x61 + x63 + x66 + x69 - x70)
      x73 = exp(x13*x72) - 1.0d0
      x74 = 3.53536660621036d0*x11**(-3.78d0)*x21*x71**2.78d0
      x75 = 1.0d0 + x74/x73
      x76 = -x28*log(x74*x75/(12.498817040307328d0*x11**(-7.559999999999
     &9996d0)*x19* x71**5.5599999999999996d0/x73**2 + x75) + 1.0d0) + x7
     &2
      x77 = x36*(x26 + x34)
      x78 = -rho_c*(-x44*(x43 - (x39 + 0.011953125d0*x77)*(-x42*(x40 + 0
     &.03346875d0* x77) + 512.0d0)) + x76)
      PBE_value = x29
      PBE_beta_value = x35
      PBE_alpha_value = x35
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(-x44*(x43 - (0.02390625d0*x37 + x39)*(-x42*(0.0669375d
     &0*x37 + x40) + 512.0d0)) + x76)
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = x78
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = x78
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(-x44*(-x26*(0.015625d0*x36 - 1.0d0)*(-x29*x41*(0.04375
     &d0*x36 - 2.8d0) + 512.0d0) + x43) + x76)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E,
     & d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, tau_c, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
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
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1.46459188756152d0
      x1 = 1d0/rho_c
      x2 = x1**0.333333333333333d0
      x3 = x2/x0
      x4 = 0.194159335344114d0*x3 + 1.0d0
      x5 = 2.14502939711103d0
      x6 = 1d0/x5
      x7 = x1**0.666666666666667d0
      x8 = x6*x7
      x9 = 0.564189583547756d0*sqrt(x1)
      x10 = 0.826307487110758d0*x1**0.166666666666667d0
      x11 = 7.24010193431683d0*x10 + 3.25955091942229d0*x3 + 0.406913004
     &517529d0*x8 + 1.41872281647967d0*x9
      x12 = 9.86960440108936d0
      x13 = 0.693147180559945d0
      x14 = -x13 + 1.0d0
      x15 = 1d0/x14
      x16 = x12*x15
      x17 = x16/x11
      x18 = log(0.5d0*x17 + 1.0d0)
      x19 = x18*x4
      x20 = 2.0d0*x19
      x21 = gamma_cc*x15
      x22 = exp(x20)
      x23 = x22 - 1.0d0
      x24 = x23**(-2)
      x25 = 208.945363813799d0
      x26 = rho_c**(-4.66666666666667d0)
      x27 = x25*x26
      x28 = 8.36096418563526d-6*x27
      x29 = gamma_cc**2
      x30 = x29/x14**2
      x31 = 14.4549425392770d0
      x32 = rho_c**(-2.33333333333333d0)
      x33 = x31*x32
      x34 = 0.00289153318944038d0*x33
      x35 = 1d0/x23
      x36 = x21*x35
      x37 = x34*x36 + 1.0d0
      x38 = 1d0/(x24*x28*x30 + x37)
      x39 = x34*x37
      x40 = x38*x39
      x41 = x21*x40 + 1.0d0
      x42 = log(x41)
      x43 = 1.0d0*x42
      x44 = -x20 + x43
      x45 = 1d0/x12
      x46 = x14*x45
      x47 = x44*x46
      x48 = 0.235215882559859d0*x3 + 1.0d0
      x49 = 15.10599739383d0*x10 + 7.09459546107279d0*x3 + 0.81920444870
     &8788d0*x8 + 4.12273618607837d0*x9
      x50 = x16/x49 + 1.0d0
      x51 = log(x50)
      x52 = x48*x51
      x53 = exp(2.0000000000000004d0*x52)
      x54 = x53 - 1.0d0
      x55 = x54**(-2)
      x56 = 3.34438567425411d-5*x27
      x57 = 0.00578306637888076d0*x33
      x58 = 1d0/x54
      x59 = x21*x58
      x60 = x57*x59 + 1.0d0
      x61 = 1d0/(x30*x55*x56 + x60)
      x62 = x57*x60
      x63 = x61*x62
      x64 = x21*x63 + 1.0d0
      x65 = log(x64)
      x66 = -1.0d0*x52 + 0.5d0*x65
      x67 = x46*x66
      x68 = rho_c**(-2)
      x69 = tau_c**(-2)
      x70 = x29*x69
      x71 = x68*x70
      x72 = x66*x71
      x73 = 0.00828125d0*x71 + 1.0d0
      x74 = x44*x73
      x75 = -x74
      x76 = 0.02390625d0*x72 + x75
      x77 = -2.8d0*x74
      x78 = 0.0669375d0*x72 + x77
      x79 = x46*x78
      x80 = rho_c**(-3)
      x81 = tau_c**(-3)
      x82 = x80*x81
      x83 = gamma_cc**3*x82
      x84 = -x79*x83 + 512.0d0
      x85 = x76*x84
      x86 = 1024.0d0*x19
      x87 = -512.0d0*x42 + x86
      x88 = exp(-1.0053738478084204d0*mu*rho_c**(-0.33333333333333331d0)
     &*x6)
      x89 = 0.001953125d0*x88
      x90 = x46*x89
      x91 = x20*x46
      x92 = mu**2
      x93 = x8*x92
      x94 = 0.508616435555896d0*x93 + 1.0d0
      x95 = x94**(-4)
      x96 = mu**4
      x97 = exp(-0.6836107611867116d0*x3)
      x98 = 4.60115111447049d0
      x99 = x1**1.33333333333333d0
      x100 = x99/x98
      x101 = 3.14159265358979d0
      x102 = x1**1.0d0
      x103 = x102/x101
      x104 = 0.009578475d0*x103
      x105 = -x104
      x106 = 0.00126674656487366d0*x100 + x105 + 0.0188171922990732d0*x3
     & + 0.0676322201645715d0*x8 + 1.0d0
      x107 = x106*x97
      x108 = x107 - 1.0d0
      x109 = 0.13157433081915d0*x103*x108
      x110 = 3.10428814221102d0*x99
      x111 = x14*x19
      x112 = 0.022020833726518d0*x111
      x113 = x110*x112
      x114 = 2.28942848510666d0
      x115 = x114*x98
      x116 = 1.63540853354893d0*x115
      x117 = 0.558025705063192d0*x8
      x118 = 0.352521395009435d0*x3
      x119 = x117 - x118
      x120 = exp(-0.49698248213959023d0*x3)
      x121 = x120*x5
      x122 = 10.9027235569928d0*x121
      x123 = x119*x122
      x124 = 0.494402081358784d0*x3 + 0.0524148278841779d0*x8 + 1.0d0
      x125 = 1d0/x124
      x126 = 0.0259335011650457d0*x3
      x127 = -x126 + 1.0d0
      x128 = x125*x127
      x129 = 1.63540853354893d0*x115
      x130 = x128*x129
      x131 = 0.148394183600699d0
      x132 = x1**1.0d0*x131
      x133 = x132*(x116 - x123 - x130)
      x134 = 0.00139418473233101d0*x133
      x135 = x96*(-x109 - x113 + x134)
      x136 = x3*x92
      x137 = 5.57236303607474d0*mu*x10 + 1.0d0
      x138 = mu**3
      x139 = 6.76834298980504d0*x136 + x137 + 3.39260255780131d0*x138*x9
     &
      x140 = x139/(3.13317926779378d0*x136 + x137)
      x141 = 2.0d0*log(x140)
      x142 = x45*(x13 - 1.0d0)
      x143 = x141*x142
      x144 = x1**1.66666666666667d0
      x145 = x131*x144
      x146 = x108*x145
      x147 = x1**2.0d0
      x148 = 1.0525946465532d0*x147
      x149 = 0.0102659822546843d0*x111
      x150 = mu**6*(0.0167302167879722d0*x146 + x148*x149)
      x151 = 1.41421356237310d0
      x152 = 0.0315054072231411d0*x102
      x153 = x107*x152
      x154 = exp(-0.28165369188898165d0*x3)
      x155 = x154*(-4.49737346725955d0*x3 + 0.825481812223657d0*x8)
      x156 = 12.0d0*x101/x102
      x157 = x155*x156
      x158 = 1d0/x7
      x159 = x115*x158
      x160 = x128*x159
      x161 = x157 + 1.63540853354893d0*x160
      x162 = 0.0837224526465878d0*x144
      x163 = 0.00111534778586481d0*x162
      x164 = x138*x151*(x153 + x161*x163)
      x165 = mu**5
      x166 = x107*x151
      x167 = 0.0223069557172962d0*x162
      x168 = x165*x166*x167
      x169 = x1**2.66666666666667d0
      x170 = 0.133841734303777d0*x169
      x171 = 0.00478594012208448d0*x111
      x172 = mu**8*x170*x171
      x173 = x150 + x164 + x168 + x172
      x174 = -x135 - x143 + x173
      x175 = -x174*x95 + x91
      x176 = exp(x16*x175)
      x177 = x176 - 1.0d0
      x178 = 1d0/x177
      x179 = x135 + x143 - x150 - x164 - x168 - x172
      x180 = x179*x95 + x91
      x181 = 1d0/x4
      x182 = 1d0/x18
      x183 = x181*x182
      x184 = x180*x183
      x185 = x184**2.78d0
      x186 = 8397.87994635719d0*x14**(-3.78d0)
      x187 = 0.000420983227766184d0*x186
      x188 = x187*x32
      x189 = x185*x188
      x190 = gamma_cc*x189
      x191 = x178*x190 + 1.0d0
      x192 = x184**5.56d0
      x193 = x177**(-2)
      x194 = 70524387.5934282d0*x14**(-7.56d0)
      x195 = x194*x26
      x196 = 1.77226878060434d-7*x195*x29
      x197 = 1d0/(x191 + x192*x193*x196)
      x198 = x175 - x46*log(x190*x191*x197 + 1.0d0)
      x199 = x2/x31
      x200 = (2.517666232305d0*x10 + 2.36486515369093d0*x3 + 0.546136299
     &139192d0*x8 + 2.06136809303918d0*x9)/(x49**2*x50)
      x201 = rho_c**(-3.33333333333333d0)
      x202 = x201*x31
      x203 = 0.0134938215507218d0*x202
      x204 = 0.00578306637888076d0*x58
      x205 = x53*(x16*x200*(0.470431765119717d0*x3 + 2.0d0) - 0.15681058
     &8373239d0*x3* x51)
      x206 = x204*x205
      x207 = rho_c**(-5.66666666666667d0)
      x208 = x207*x25
      x209 = x208*x21
      x210 = x59*x63
      x211 = x61/x64
      x212 = gamma_cc*x45
      x213 = 0.0784052941866195d0*x1*x14*x199*x51 - 1.0d0*x1*x200*x48 + 
     &0.5d0*x211* x212*(-x203*x60 - x204*x209*(x206 + 0.0134938215507218
     &d0) + x210* (x202*x206 + x203 + 6.68877134850821d-5*x205*x209*x55 
     &+ 0.000156071331465192d0*x208*x59))
      x214 = x213*x71
      x215 = x70*x80
      x216 = x215*x67
      x217 = 0.00674691077536089d0*x202
      x218 = x217*x37
      x219 = (1.20668365571947d0*x10 + 1.08651697314076d0*x3 + 0.2712753
     &3634502d0*x8 + 0.709361408239834d0*x9)/(x11**2*(x17 + 2.0d0))
      x220 = x22*(2.0d0*x16*x219*(0.388318670688228d0*x3 + 2.0d0) - 0.25
     &8879113792152d0*x18*x3)
      x221 = 0.00144576659472019d0*x35
      x222 = x220*x221
      x223 = x208*x36
      x224 = 0.00289153318944038d0*x223*(x222 + 0.00674691077536089d0)
      x225 = x36*x40
      x226 = x225*(x202*x222 + 8.36096418563526d-6*x209*x220*x24 + x217 
     &+ 3.90178328662979d-5*x223)
      x227 = x38/x41
      x228 = x212*x227
      x229 = 1.0d0*x228
      x230 = x219*x4
      x231 = 2.0d0*x230
      x232 = x14*x18
      x233 = x199*x232
      x234 = 0.129439556896076d0*x233
      x235 = -x1*x231 + x1*x234
      x236 = x229*(-x218 - x224 + x226) + x235
      x237 = x215*x47
      x238 = -x236*x73 + 0.0165625d0*x237
      x239 = 3.0d0*x1
      x240 = 0.0231875d0*x71 + 2.8d0
      x241 = x236*x240 - 0.046375d0*x237
      x242 = x46*x76
      x243 = x218 + x224 - x226
      x244 = x228*x243
      x245 = 1024.0d0*x1*x230
      x246 = 66.2730531307909d0*x1*x233
      x247 = -x245 + x246
      x248 = -512.0d0*x244 + x247
      x249 = 512.0d0*x42 - x86
      x250 = x249 + x85
      x251 = 0.0158297237324883d0*mu*rho_c**(-1.33333333333333d0)*x14*x8
     &8
      x252 = 1d0/x2
      x253 = 0.0141372897033723d0*x101*x114*x125*x252
      x254 = -x117 + x118
      x255 = 7.26848237132856d0*x121*x158
      x256 = 1.80615420514536d0*x0*x120*x252
      x257 = x126 - 1.0d0
      x258 = x125*x257
      x259 = 1.63540853354893d0*x159*(0.164800693786261d0*x3 + 0.0349432
     &18589452d0*x8 )/x124**2
      x260 = x122*x158*(0.117507131669812d0*x3 - 0.372017136708795d0*x8)
     &
      x261 = 1.09027235569928d0*x159
      x262 = x260 - x261
      x263 = 0.00139418473233101d0*x145
      x264 = 0.00232364122055169d0*x132
      x265 = 0.000844497709915773d0*x100
      x266 = 0.0047892375d0*x103
      x267 = 0.0225440733881905d0*x8
      x268 = 0.0031361987165122d0*x3
      x269 = 0.113935126864452d0*x106*x3
      x270 = x97*(x265 - x266 + x267 + x268 - x269)
      x271 = 0.263148661638301d0*x103
      x272 = x109 - 0.217336917462899d0*x110*x230 + 4.13905085628136d0*x
     &112*x99 + 0.00302075971817057d0*x144*x232
      x273 = 0.0334604335759443d0*x145
      x274 = 0.00047750955226877d0*x1**2.33333333333333d0*x232 + 0.02788
     &36946466203d0 *x146 + 2.10518929310641d0*x147*x149 - 0.10132118364
     &2338d0*x148* x230
      x275 = 0.00168899541983155d0*x100
      x276 = 0.045088146776381d0*x8
      x277 = 0.00627239743302441d0*x3
      x278 = x105 + x275 + x276 + x277
      x279 = x152*x97
      x280 = 1.09027235569928d0*x160
      x281 = x127*x259
      x282 = x253 + x280 + x281
      x283 = x163*(2.41662179562687d0*x1**(-0.666666666666667d0)*x155 + 
     &x154*x156*( 1.49912448908652d0*x3 - 0.550321208149104d0*x8) + x157
     & + x282)
      x284 = 0.00185891297644135d0*x161*x162
      x285 = 0.00490180588786583d0*x107*x99
      x286 = -x153 + x283 - x284 + x285
      x287 = x151*x92
      x288 = 0.928727172679123d0*x10
      x289 = mu*x3
      x290 = 2.0d0*(x140*(x288 + 1.04439308926459d0*x289) - x288 - 2.256
     &11432993501d0 *x289 - 1.69630127890066d0*x9*x92)/x139
      x291 = x151*x167*x96*x97
      x292 = mu**7
      x293 = x166*x96
      x294 = -0.000290571663240495d0*x1**2.0d0*x293 + 2.83060463816878d-
     &5*x1**3.0d0* x232*x292 + 0.037178259528827d0*x162*x293 + 0.3569112
     &91476739d0* x169*x171*x292 - 0.0472353356922751d0*x170*x230*x292
      x295 = x278*x291 - x290*x46 + x294
      x296 = mu*x95
      x297 = x175*x183
      x298 = x297**5.56d0
      x299 = x141*x46 + x173 + x96*(x109 + x113 - x134)
      x300 = x299*x95
      x301 = exp(x16*(-x300 + x91)) - 1.0d0
      x302 = x188*x297**2.78d0
      x303 = gamma_cc/x301
      x304 = x302*x303 + 1.0d0
      x305 = 1d0/(x196*x298/x301**2 + x304)
      x306 = x302*x304
      x307 = gamma_cc*x305*x306
      x308 = x307 + 1.0d0
      x309 = 1d0/x308
      x310 = x185*x201
      x311 = x191*x310
      x312 = 1d0/x180
      x313 = 0.359841968171091d0*x233
      x314 = 5.56d0*x230
      x315 = x93/x94**5
      x316 = x179*x315
      x317 = 3.77054317558771d0*x316
      x318 = x180*x181*x3
      x319 = 0.179920984085546d0*x318
      x320 = x16*x180
      x321 = x182*x219*x320
      x322 = 2.78d0*x321
      x323 = x119*x256
      x324 = x119*x255
      x325 = x97*(-x265 + x266 - x267 - x268 + x269)
      x326 = -x271*x325 + x272
      x327 = x165*(-x273*x325 + x274)
      x328 = x104 - x275 - x276 - x277
      x329 = x279*x328
      x330 = x296*(x138*(-0.00232364122055169d0*x133 + x263*(-x253 - x26
     &0 + x261 - x280 - x281 - x323 - x324) + x326) + x142*x290 + x287*(
     &x153 - x283 + x284 - x285 - x329) - x291*x328 + x294 + x327)
      x331 = 2.78d0*x330
      x332 = x313 - x314 - x317 - x319 + x322 - x331
      x333 = 0.000982294198121095d0*x186
      x334 = x231 - x234
      x335 = exp(x320)
      x336 = 0.000420983227766184d0/(x335 - 1.0d0)
      x337 = 82883.752878387d0*x14**(-4.78d0)
      x338 = x187*x312
      x339 = gamma_cc*x192
      x340 = x207*x339
      x341 = x178*x194*x340
      x342 = 0.000420983227766184d0*x178
      x343 = x176*(-1.35631049481572d0*x174*x315 + x296*(x138*(-x263*(x2
     &62 + x282 + x323 + x324) + x264*(-x116 + x123 + x130) + x326) - x2
     &87*(x286 + x329) + x295 + x327) + x334)
      x344 = gamma_cc*x197
      x345 = -512.0d0*x1*x296*(x138*(-x263*(-1.09027235569928d0*x159*x25
     &8 + x253 - x254*x255 - x254*x256 - x257*x259 + x262) - x264*(x116 
     &+ x122* x254 + x129*x258) + x270*x271 + x272) + x165*(x270*x273 + 
     &x274) - x287*(-x278*x279 + x286) + x295) + 694.43097334565d0*x1*x2
     &99*x315 - 435652.165760595d0*x14**(-2.78d0)*x309*x344*(-x185*x191*
     &x32* x342*x344*(246.716759338856d0*x14**(-8.56d0)*x193*x340*x343 +
     & x310*x332*x338 + x310*x333 + x310*x337*x342*x343 + 1.772268780604
     &34d-7*x312*x341*(-11.12d0*x230 + 0.719683936342183d0*x233 - 7.5410
     &8635117542d0*x316 - 0.359841968171091d0*x318 + 5.56d0*x321 - 5.56d
     &0*x330) + 8.27058764282027d-7*x341) + 0.000420983227766184d0*x311*
     &x312*x332 + 0.000982294198121095d0*x311 + x336*x340*(x333 + x335*x
     &336*x337* (1.35631049481572d0*x316 + x330 + x334) - x338*(-x313 + 
     &x314 + x317 + x319 - x322 + x331))) + x247
      x346 = 0.001953125d0*rho_c
      x347 = x300 + x46*log(x308) - x91
      x348 = x27*x36
      x349 = x34 + 1.67219283712705d-5*x348
      x350 = x227*(-x225*x349 + x28*x36 + x39)
      x351 = 512.0d0*x350
      x352 = gamma_cc*x14*x68*x69
      x353 = x352*x66
      x354 = 6.68877134850821d-5*x27*x59 + x57
      x355 = x211*x71*(-x210*x354 + x56*x59 + x62)
      x356 = 1.0d0*x350
      x357 = -x356*x73
      x358 = x352*x44
      x359 = x357 - 0.0165625d0*x358
      x360 = 3.0d0*x14
      x361 = -x240*x356
      x362 = -0.046375d0*x358 + x361
      x363 = x29*x82
      x364 = x14*x305*x309*(-x178*x307*(3.54453756120869d-7*x178*x195*x3
     &39 + x189) + 1.77226878060434d-7*x195*x298*x303 + x306)
      x365 = rho_c*x45
      x366 = -8.4d0*x74
      x367 = gamma_cc*x1/tau_c
      x368 = x1*x29*x81*x90
      x369 = x44 + x66
      x370 = x369*x71
      x371 = 0.011953125d0*x370 + x75
      x372 = 0.03346875d0*x370 + x77
      x373 = x372*x46
      x374 = -x373*x83 + 512.0d0
      x375 = x371*x374
      x376 = -rho_c*(x198 - x90*(-x375 + x87))
      x377 = 0.5d0*x46
      x378 = 0.5d0*rho_c
      x379 = x215*(-1.0d0*x19*x46 + x213*x378 + x236*x378 + x377*x42 - x
     &377*x52 + 0.25d0*x46*x65)
      x380 = x215*x46
      x381 = x369*x380
      x382 = x371*x46
      x383 = x249 + x375
      x384 = x346*(-x251*x383 + x345 - x88*(x248 + x374*(x238 + 0.023906
     &25d0*x379 - 0.035859375d0*x381) + x382*x83*(x239*x373 + x241 - 0.0
     &669375d0* x379 + 0.10040625d0*x381))) + x347 - x383*x90
      x385 = 0.00144576659472019d0*x33
      x386 = x385*x60
      x387 = x71*(x211*(x28*x59 - x354*x386*x59*x61 + x386) + x227*(-x21
     &*x221*x33* x349*x37*x38 + 4.18048209281763d-6*x348 + x37*x385))
      x388 = x352*x369
      x389 = x365*(x364 - x89*(x351 - x363*x382*(gamma_cc*(x362 + 0.0669
     &375d0*x387 + 0.0669375d0*x388) + x360*x372) + x374*(x359 + 0.02390
     &625d0*x387 + 0.02390625d0*x388)))
      x390 = -x368*(x367*x382*(x366 + 0.10040625d0*x370 - x71*(0.041125d
     &0*x19 - 0.0205625d0*x42 + 0.0669375d0*x52 - 0.03346875d0*x65)) + x
     &374*( 0.0146875d0*x19 - 0.00734375d0*x42 + 0.02390625d0*x52 - 0.01
     &1953125d0*x65))
      x391 = 0.015625d0*x71
      x392 = x391 - 1.0d0
      x393 = 0.04375d0*x71
      x394 = x393 - 2.8d0
      x395 = x394*x47*x83
      x396 = -x395 + 512.0d0
      x397 = x392*x396*x44
      x398 = x20 - x43
      x399 = 0.03125d0*x398
      x400 = -x229*x243 + x235
      x401 = x400*x71
      x402 = x398*x46
      x403 = x249 + x397
      x404 = x350*x71
      x405 = x394*x44
      x406 = x392*x47
      PBE_value = x47
      PBE_beta_value = x67
      PBE_alpha_value = x67
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x198 - x90*(-x85 + x87))
      d1E(1) = -x250*x90 + x346*(-x250*x251 + x345 - x88*(x242*x83*(-0.0
     &669375d0*x214 + 0.133875d0*x216 + x239*x79 + x241) + x248 + x84*(0
     &.02390625d0* x214 - 0.0478125d0*x216 + x238))) + x347
      d1E(3) = x365*(x364 - x89*(-x242*x363*(gamma_cc*(0.133875d0*x353 +
     & 0.03346875d0* x355 + x362) + x360*x78) + x351 + x84*(0.0478125d0*
     &x353 + 0.011953125d0*x355 + x359)))
      d1E(6) = -x368*(x242*x367*(x366 - x71*(-0.09275d0*x19 + 0.046375d0
     &*x42 + 0.133875d0*x52 - 0.0669375d0*x65) + 0.2008125d0*x72) + x84*
     &( -0.033125d0*x19 + 0.0165625d0*x42 + 0.0478125d0*x52 - 0.02390625
     &d0*x65))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = x376
      d1E(1) = x384
      d1E(3) = x389
      d1E(6) = x390
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = x376
      d1E(1) = x384
      d1E(3) = x389
      d1E(6) = x390
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x198 - x90*(-x397 + x87))
      d1E(1) = x346*(-x251*x403 + x345 + x88*(512.0d0*x244 + x245 - x246
     & - x402*x83*( -x391 + 1.0d0)*(-0.0875d0*x215*x402 + x239*x402*(-x3
     &93 + 2.8d0) + x240*x400 - 0.0669375d0*x401) + (x395 - 512.0d0)*(x3
     &80*x399 - x400*x73 + 0.02390625d0*x401))) + x347 - x403*x90
      d1E(3) = x365*(x364 - x89*(x351 - x363*x406*(gamma_cc*(0.0875d0*x3
     &58 + x361 + 0.0669375d0*x404) + x360*x405) + x396*(x357 + 0.03125d
     &0*x358 + 0.02390625d0*x404)))
      d1E(6) = -x368*(x367*x406*(3.0d0*x405 - x71*(0.175d0*x19 - 0.0875d
     &0*x42)) + x396* x399)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E,
     & d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, tau_c, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
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
     &75
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1.46459188756152d0
      x1 = 1d0/rho_c
      x2 = x1**0.333333333333333d0
      x3 = x2/x0
      x4 = 0.194159335344114d0*x3 + 1.0d0
      x5 = 2.14502939711103d0
      x6 = 1d0/x5
      x7 = x1**0.666666666666667d0
      x8 = x6*x7
      x9 = 0.564189583547756d0*sqrt(x1)
      x10 = 1.41872281647967d0*x9
      x11 = 0.826307487110758d0*x1**0.166666666666667d0
      x12 = x10 + 7.24010193431683d0*x11 + 3.25955091942229d0*x3 + 0.406
     &913004517529d0*x8
      x13 = 9.86960440108936d0
      x14 = 0.693147180559945d0
      x15 = -x14 + 1.0d0
      x16 = 1d0/x15
      x17 = x13*x16
      x18 = x17/x12
      x19 = log(0.5d0*x18 + 1.0d0)
      x20 = x19*x4
      x21 = 2.0d0*x20
      x22 = gamma_cc*x16
      x23 = 208.945363813799d0
      x24 = rho_c**(-4.66666666666667d0)
      x25 = x23*x24
      x26 = 8.36096418563526d-6*x25
      x27 = exp(x21)
      x28 = x27 - 1.0d0
      x29 = x28**(-2)
      x30 = gamma_cc**2
      x31 = x15**(-2)
      x32 = x30*x31
      x33 = x29*x32
      x34 = x26*x33
      x35 = 14.4549425392770d0
      x36 = rho_c**(-2.33333333333333d0)
      x37 = x35*x36
      x38 = 0.00289153318944038d0*x37
      x39 = 1d0/x28
      x40 = x22*x39
      x41 = x38*x40 + 1.0d0
      x42 = x34 + x41
      x43 = 1d0/x42
      x44 = x38*x41
      x45 = x43*x44
      x46 = x22*x45 + 1.0d0
      x47 = log(x46)
      x48 = 1.0d0*x47
      x49 = -x21 + x48
      x50 = 1d0/x13
      x51 = x15*x50
      x52 = x49*x51
      x53 = 0.235215882559859d0*x3 + 1.0d0
      x54 = 4.12273618607837d0*x9
      x55 = 15.10599739383d0*x11 + 7.09459546107279d0*x3 + x54 + 0.81920
     &4448708788d0 *x8
      x56 = x17/x55 + 1.0d0
      x57 = log(x56)
      x58 = x53*x57
      x59 = 1.0d0*x58
      x60 = 3.34438567425411d-5*x25
      x61 = exp(2.0000000000000004d0*x58)
      x62 = x61 - 1.0d0
      x63 = x62**(-2)
      x64 = x32*x63
      x65 = x60*x64
      x66 = 0.00578306637888076d0*x37
      x67 = 1d0/x62
      x68 = x22*x67
      x69 = x66*x68 + 1.0d0
      x70 = x65 + x69
      x71 = 1d0/x70
      x72 = x66*x69
      x73 = x71*x72
      x74 = x22*x73 + 1.0d0
      x75 = log(x74)
      x76 = 0.5d0*x75
      x77 = -x59 + x76
      x78 = x51*x77
      x79 = rho_c**(-2)
      x80 = tau_c**(-2)
      x81 = x30*x80
      x82 = x79*x81
      x83 = x77*x82
      x84 = 0.00828125d0*x82 + 1.0d0
      x85 = x49*x84
      x86 = -x85
      x87 = 0.02390625d0*x83 + x86
      x88 = -2.8d0*x85
      x89 = 0.0669375d0*x83 + x88
      x90 = rho_c**(-3)
      x91 = x51*x90
      x92 = tau_c**(-3)
      x93 = gamma_cc**3*x92
      x94 = x91*x93
      x95 = -x89*x94 + 512.0d0
      x96 = x87*x95
      x97 = 1024.0d0*x20
      x98 = -512.0d0*x47 + x97
      x99 = mu*x6
      x100 = exp(-1.0053738478084204d0*rho_c**(-0.33333333333333331d0)*x
     &99)
      x101 = 0.001953125d0*x100
      x102 = x101*x51
      x103 = x21*x51
      x104 = mu**2
      x105 = x104*x8
      x106 = 0.508616435555896d0*x105 + 1.0d0
      x107 = x106**(-4)
      x108 = mu**5
      x109 = exp(-0.6836107611867116d0*x3)
      x110 = 4.60115111447049d0
      x111 = x1**1.33333333333333d0
      x112 = x111/x110
      x113 = 3.14159265358979d0
      x114 = x1**1.0d0
      x115 = x114/x113
      x116 = 0.009578475d0*x115
      x117 = -x116
      x118 = 0.00126674656487366d0*x112 + x117 + 0.0188171922990732d0*x3
     & + 0.0676322201645715d0*x8 + 1.0d0
      x119 = x109*x118
      x120 = x1**1.66666666666667d0
      x121 = 0.0837224526465878d0*x120
      x122 = x119*x121
      x123 = 1.41421356237310d0
      x124 = 0.0223069557172962d0*x123
      x125 = x108*x122*x124
      x126 = x104*x3
      x127 = 5.57236303607474d0*mu*x11 + 1.0d0
      x128 = 3.13317926779378d0*x126 + x127
      x129 = 1d0/x128
      x130 = mu**3
      x131 = 6.76834298980504d0*x126 + x127 + 3.39260255780131d0*x130*x9
     &
      x132 = x129*x131
      x133 = 2.0d0*log(x132)
      x134 = x50*(x14 - 1.0d0)
      x135 = x133*x134
      x136 = x1**2.66666666666667d0
      x137 = 0.133841734303777d0*x136
      x138 = x15*x20
      x139 = 0.00478594012208448d0*x138
      x140 = mu**8*x137*x139
      x141 = x119 - 1.0d0
      x142 = 0.148394183600699d0
      x143 = x120*x142
      x144 = x141*x143
      x145 = x1**2.0d0
      x146 = 1.0525946465532d0*x145
      x147 = 0.0102659822546843d0*x138
      x148 = mu**6*(0.0167302167879722d0*x144 + x146*x147)
      x149 = 0.179587122125167d0*x114
      x150 = 0.175432441092201d0*x149
      x151 = x119*x150
      x152 = exp(-0.28165369188898165d0*x3)
      x153 = x152*(-4.49737346725955d0*x3 + 0.825481812223657d0*x8)
      x154 = x113/x114
      x155 = 12.0d0*x154
      x156 = x153*x155
      x157 = 1d0/x7
      x158 = 2.28942848510666d0
      x159 = x110*x158
      x160 = x157*x159
      x161 = 0.494402081358784d0*x3 + 0.0524148278841779d0*x8 + 1.0d0
      x162 = 1d0/x161
      x163 = 0.0259335011650457d0*x3
      x164 = -x163 + 1.0d0
      x165 = x162*x164
      x166 = x160*x165
      x167 = x156 + 1.63540853354893d0*x166
      x168 = 0.00111534778586481d0*x121
      x169 = x123*x130*(x151 + x167*x168)
      x170 = mu**4
      x171 = x115*x141
      x172 = 0.13157433081915d0*x171
      x173 = 3.10428814221102d0*x111
      x174 = 0.0220208337265180d0
      x175 = x138*x174
      x176 = x173*x175
      x177 = 1.63540853354893d0*x159
      x178 = 0.558025705063192d0*x8
      x179 = 0.352521395009435d0*x3
      x180 = x178 - x179
      x181 = exp(-0.49698248213959023d0*x3)
      x182 = x181*x5
      x183 = 10.9027235569928d0*x182
      x184 = x180*x183
      x185 = 1.63540853354893d0*x159
      x186 = x165*x185
      x187 = x1**1.0d0*x142
      x188 = x187*(x177 - x184 - x186)
      x189 = 0.00139418473233101d0*x188
      x190 = x170*(-x172 - x176 + x189)
      x191 = -x125 + x135 - x140 - x148 - x169 + x190
      x192 = x103 + x107*x191
      x193 = 1d0/x4
      x194 = 1d0/x19
      x195 = x193*x194
      x196 = x192*x195
      x197 = x196**2.78d0
      x198 = x15**(-3.78d0)
      x199 = 8397.87994635719d0*x198
      x200 = 0.000420983227766184d0*x199
      x201 = x197*x36
      x202 = x200*x201
      x203 = x125 + x140 + x148 + x169
      x204 = -x135 - x190 + x203
      x205 = x103 - x107*x204
      x206 = x17*x205
      x207 = exp(x206)
      x208 = x207 - 1.0d0
      x209 = 1d0/x208
      x210 = gamma_cc*x209
      x211 = x202*x210 + 1.0d0
      x212 = x197*x211
      x213 = x200*x36
      x214 = x212*x213
      x215 = x208**(-2)
      x216 = x196**5.56d0
      x217 = x15**(-7.56d0)
      x218 = 70524387.5934282d0*x217
      x219 = x218*x24
      x220 = x216*x219
      x221 = 1.77226878060434d-7*x30
      x222 = x220*x221
      x223 = x211 + x215*x222
      x224 = 1d0/x223
      x225 = gamma_cc*x224
      x226 = x214*x225 + 1.0d0
      x227 = x205 - x51*log(x226)
      x228 = 1d0/x35
      x229 = x2*x228
      x230 = x15*x229*x57
      x231 = x1*x230
      x232 = 1.0d0*x53
      x233 = 2.517666232305d0*x11 + 2.36486515369093d0*x3 + 0.5461362991
     &39192d0*x8 + 2.06136809303918d0*x9
      x234 = 1d0/x56
      x235 = x234/x55**2
      x236 = x233*x235
      x237 = x1*x236
      x238 = 0.00578306637888076d0*x67
      x239 = x22*x238
      x240 = x3*x57
      x241 = 0.470431765119717d0*x3 + 2.0d0
      x242 = x17*x241
      x243 = x236*x242 - 0.156810588373239d0*x240
      x244 = x243*x61
      x245 = x238*x244
      x246 = x245 + 0.0134938215507218d0
      x247 = rho_c**(-5.66666666666667d0)
      x248 = x23*x247
      x249 = x246*x248
      x250 = x239*x249
      x251 = rho_c**(-3.33333333333333d0)
      x252 = x251*x35
      x253 = 0.0134938215507218d0*x252
      x254 = x253*x69
      x255 = x248*x68
      x256 = 6.68877134850821d-5*x63
      x257 = x22*x248
      x258 = x244*x257
      x259 = x245*x252 + x253
      x260 = 0.000156071331465192d0*x255 + x256*x258 + x259
      x261 = x68*x73
      x262 = x260*x261
      x263 = -x254 + x262
      x264 = -x250 + x263
      x265 = x71/x74
      x266 = gamma_cc*x265
      x267 = x264*x266
      x268 = 0.5d0*x50
      x269 = 0.0784052941866195d0*x231 - x232*x237 + x267*x268
      x270 = x269*x82
      x271 = x81*x90
      x272 = x271*x78
      x273 = x19*x3
      x274 = 1.20668365571947d0*x11 + 1.08651697314076d0*x3 + 0.27127533
     &634502d0*x8 + 0.709361408239834d0*x9
      x275 = x18 + 2.0d0
      x276 = 1d0/x275
      x277 = x276/x12**2
      x278 = x274*x277
      x279 = 0.388318670688228d0*x3 + 2.0d0
      x280 = x17*x279
      x281 = -0.258879113792152d0*x273 + 2.0d0*x278*x280
      x282 = x27*x39
      x283 = x281*x282
      x284 = 0.00144576659472019d0*x283
      x285 = x284 + 0.00674691077536089d0
      x286 = x248*x40
      x287 = x285*x286
      x288 = 0.00289153318944038d0*x287
      x289 = 0.00674691077536089d0*x252
      x290 = x289*x41
      x291 = x27*x29
      x292 = x281*x291
      x293 = x257*x292
      x294 = x252*x284 + x289
      x295 = 3.90178328662979d-5*x286 + 8.36096418563526d-6*x293 + x294
      x296 = x40*x45
      x297 = x295*x296
      x298 = -x290 + x297
      x299 = -x288 + x298
      x300 = x43/x46
      x301 = x300*x50
      x302 = gamma_cc*x301
      x303 = x299*x302
      x304 = 1.0d0*x303
      x305 = x278*x4
      x306 = 2.0d0*x305
      x307 = x15*x19
      x308 = x229*x307
      x309 = 0.129439556896076d0*x308
      x310 = x1*x309
      x311 = -x1*x306 + x310
      x312 = x304 + x311
      x313 = x271*x52
      x314 = -x312*x84 + 0.0165625d0*x313
      x315 = 0.02390625d0*x270 - 0.0478125d0*x272 + x314
      x316 = x1*x51
      x317 = x316*x89
      x318 = 0.0669375d0*x270
      x319 = 0.133875d0*x272
      x320 = 0.0231875d0*x82 + 2.8d0
      x321 = x312*x320
      x322 = 0.046375d0*x313
      x323 = x321 - x322
      x324 = 3.0d0*x317 - x318 + x319 + x323
      x325 = x87*x94
      x326 = x290 - x297
      x327 = x288 + x326
      x328 = x302*x327
      x329 = x1*x305
      x330 = 1024.0d0*x329
      x331 = x1*x308
      x332 = 66.2730531307909d0*x331
      x333 = -x330 + x332
      x334 = -512.0d0*x328 + x333
      x335 = x100*(x315*x95 + x324*x325 + x334)
      x336 = 0.0472353356922751d0
      x337 = 0.33512461593614d0*mu*x336
      x338 = rho_c**(-1.33333333333333d0)
      x339 = x100*x338
      x340 = x337*x339
      x341 = 512.0d0*x47 - x97
      x342 = x341 + x96
      x343 = x15*x342
      x344 = mu*x107
      x345 = 512.0d0*x344
      x346 = 1d0/x2
      x347 = x113*x158*x346
      x348 = 0.0141372897033723d0*x162*x347
      x349 = -x178 + x179
      x350 = x157*x182
      x351 = 7.26848237132856d0*x350
      x352 = x0*x181*x346
      x353 = 1.80615420514536d0*x352
      x354 = x163 - 1.0d0
      x355 = x162*x354
      x356 = x161**(-2)
      x357 = 0.164800693786261d0*x3 + 0.034943218589452d0*x8
      x358 = x356*x357
      x359 = x160*x358
      x360 = 1.63540853354893d0*x359
      x361 = 0.117507131669812d0*x3 - 0.372017136708795d0*x8
      x362 = x157*x361
      x363 = x183*x362
      x364 = 1.09027235569928d0*x160
      x365 = x363 - x364
      x366 = 0.00139418473233101d0*x143
      x367 = 0.00232364122055169d0*x187
      x368 = x118*x3
      x369 = 0.113935126864452d0*x368
      x370 = 0.000844497709915773d0*x112
      x371 = 0.0047892375d0*x115
      x372 = 0.0225440733881905d0*x8
      x373 = 0.0031361987165122d0*x3
      x374 = x370 - x371 + x372 + x373
      x375 = x109*(-x369 + x374)
      x376 = 0.263148661638301d0*x115
      x377 = 0.200908840802769d0*x120
      x378 = 0.0150354743280612d0*x307
      x379 = x111*x175
      x380 = 0.217336917462899d0
      x381 = x305*x380
      x382 = x172 - x173*x381 + x377*x378 + 4.13905085628136d0*x379
      x383 = 0.0334604335759443d0*x143
      x384 = x1**2.33333333333333d0
      x385 = 0.0681236923205142d0*x384
      x386 = 0.00700944907716014d0*x307
      x387 = x145*x147
      x388 = 0.101321183642338d0
      x389 = x146*x388
      x390 = 0.0278836946466203d0*x144 - x305*x389 + x385*x386 + 2.10518
     &929310641d0* x387
      x391 = 0.00168899541983155d0*x112
      x392 = 0.045088146776381d0*x8
      x393 = 0.00627239743302441d0*x3
      x394 = x117 + x391 + x392 + x393
      x395 = x109*x150
      x396 = 1.49912448908652d0*x3 - 0.550321208149104d0*x8
      x397 = x152*x155
      x398 = 2.14502939711103d0*x1**(-0.666666666666667d0)
      x399 = x153*x398
      x400 = 1.09027235569928d0*x166
      x401 = x164*x360
      x402 = x348 + x400 + x401
      x403 = x156 + x396*x397 + 1.12661476755593d0*x399 + x402
      x404 = x168*x403
      x405 = x121*x167
      x406 = 0.00185891297644135d0*x405
      x407 = 0.122619224952946d0*x111
      x408 = x119*x407
      x409 = 0.0399758348639607d0*x408
      x410 = -x151 + x404 - x406 + x409
      x411 = x104*x123
      x412 = 0.928727172679123d0*x11
      x413 = mu*x3
      x414 = 2.25611432993501d0*x413
      x415 = x104*x9
      x416 = 1.69630127890066d0*x415
      x417 = x412 + 1.04439308926459d0*x413
      x418 = x132*x417 - x412 - x414 - x416
      x419 = 2.0d0/x131
      x420 = x418*x419
      x421 = x420*x51
      x422 = x109*x394
      x423 = x170*x422
      x424 = x121*x124
      x425 = mu**7
      x426 = x1**3.0d0
      x427 = 0.00866220739124163d0*x425*x426
      x428 = 0.00326776364305339d0*x307
      x429 = x136*x425
      x430 = x139*x429
      x431 = x123*x170
      x432 = x122*x431
      x433 = 0.0571643564037363d0*x1**2.0d0
      x434 = x119*x431
      x435 = x433*x434
      x436 = 0.0472353356922751d0
      x437 = x305*x436
      x438 = x137*x425
      x439 = x427*x428 + 0.356911291476739d0*x430 + 0.037178259528827d0*
     &x432 - 0.00508309165921971d0*x435 - x437*x438
      x440 = -x421 + x423*x424 + x439
      x441 = x108*(x375*x383 + x390) + x130*(-x366*(-1.09027235569928d0*
     &x160*x355 + x348 - x349*x351 - x349*x353 - x354*x360 + x365) - x36
     &7*(x177 + x183*x349 + x185*x355) + x375*x376 + x382) - x411*(-x394
     &*x395 + x410) + x440
      x442 = x1*x441
      x443 = 0.000982294198121095d0*x251
      x444 = 0.000420983227766184d0*x251
      x445 = 1d0/x192
      x446 = 0.359841968171091d0*x308
      x447 = 5.56d0*x4
      x448 = x278*x447
      x449 = x106**(-5)
      x450 = x105*x449
      x451 = x191*x450
      x452 = 3.77054317558771d0*x451
      x453 = x193*x3
      x454 = x192*x453
      x455 = 0.179920984085546d0*x454
      x456 = x17*x192
      x457 = x194*x278
      x458 = x456*x457
      x459 = 2.78d0*x458
      x460 = x180*x353
      x461 = x180*x351
      x462 = -x348 - x363 + x364 - x400 - x401 - x460 - x461
      x463 = -x370 + x371 - x372 - x373
      x464 = x109*(x369 + x463)
      x465 = -x376*x464 + x382
      x466 = x108*(-x383*x464 + x390)
      x467 = x116 - x391 - x392 - x393
      x468 = x395*x467
      x469 = x134*x420
      x470 = x109*x467
      x471 = x170*x424
      x472 = x130*(-0.00232364122055169d0*x188 + x366*x462 + x465) + x41
     &1*(x151 - x404 + x406 - x409 - x468) + x439 + x466 + x469 - x470*x
     &471
      x473 = x344*x472
      x474 = 2.78d0*x473
      x475 = x446 - x448 - x452 - x455 + x459 - x474
      x476 = x445*x475
      x477 = exp(x456)
      x478 = x477 - 1.0d0
      x479 = 1d0/x478
      x480 = 0.000420983227766184d0*x479
      x481 = 0.000982294198121095d0*x199
      x482 = 1.35631049481572d0*x451
      x483 = x306 - x309
      x484 = x473 + x482 + x483
      x485 = 82883.752878387d0*x15**(-4.78d0)
      x486 = x200*x445
      x487 = x477*x480*x484*x485 + x481 - x486*(-x446 + x448 + x452 + x4
     &55 - x459 + x474)
      x488 = gamma_cc*x247
      x489 = x216*x488
      x490 = x487*x489
      x491 = 8.27058764282027d-7*x218
      x492 = x209*x489
      x493 = x218*x492
      x494 = -11.12d0*x305 + 0.719683936342183d0*x308 - 7.54108635117542
     &d0*x451 - 0.359841968171091d0*x454 + 5.56d0*x458 - 5.56d0*x473
      x495 = x445*x494
      x496 = x493*x495
      x497 = 3.54453756120869d-7*x207*x215
      x498 = -x177 + x184 + x186
      x499 = x365 + x402 + x460 + x461
      x500 = x130*(-x366*x499 + x367*x498 + x465) - x411*(x410 + x468) +
     & x440 + x466
      x501 = x344*x500
      x502 = x204*x450
      x503 = 1.35631049481572d0*x502
      x504 = x483 + x501 - x503
      x505 = 696047806.176231d0*x15**(-8.56d0)
      x506 = x504*x505
      x507 = x489*x506
      x508 = x251*x481
      x509 = x475*x486
      x510 = x251*x509
      x511 = x207*x209
      x512 = x504*x511
      x513 = x485*x512
      x514 = x197*x444*x513 + x197*x508 + x197*x510
      x515 = x491*x492 + 1.77226878060434d-7*x496 + x497*x507 + x514
      x516 = 0.000420983227766184d0*x209
      x517 = -x212*x225*x36*x515*x516 + x212*x443 + x212*x444*x476 + x48
     &0*x490
      x518 = 512.0d0*x517
      x519 = 850.883136251163d0
      x520 = x15**(-2.78d0)
      x521 = x133*x51 + x170*(x172 + x176 - x189) + x203
      x522 = x107*x521
      x523 = x17*(x103 - x522)
      x524 = exp(x523)
      x525 = x524 - 1.0d0
      x526 = x525**(-2)
      x527 = x195*x205
      x528 = x527**5.56d0
      x529 = x219*x528
      x530 = x527**2.78d0
      x531 = 1d0/x525
      x532 = gamma_cc*x531
      x533 = x213*x530*x532 + 1.0d0
      x534 = x221*x526*x529 + x533
      x535 = 1d0/x534
      x536 = x530*x533
      x537 = x213*x536
      x538 = x535*x537
      x539 = gamma_cc*x538 + 1.0d0
      x540 = 1d0/x539
      x541 = x224*x540
      x542 = gamma_cc*x519*x520*x541
      x543 = x450*x521
      x544 = x1*x543
      x545 = x333 - x345*x442 - x518*x542 + 694.43097334565d0*x544
      x546 = 0.001953125d0*rho_c
      x547 = -x103 + x51*log(x539) + x522
      x548 = x26*x40
      x549 = 1.67219283712705d-5*x25
      x550 = x40*x549
      x551 = x38 + x550
      x552 = x43*x551
      x553 = x44*x552
      x554 = x40*x553
      x555 = x44 + x548 - x554
      x556 = x300*x555
      x557 = 512.0d0*x556
      x558 = x79*x80
      x559 = gamma_cc*x558
      x560 = x15*x77
      x561 = x559*x560
      x562 = x60*x68
      x563 = 6.68877134850821d-5*x25
      x564 = x563*x68
      x565 = x564 + x66
      x566 = x261*x565
      x567 = x562 - x566 + x72
      x568 = x265*x567
      x569 = 0.011953125d0*x82
      x570 = 1.0d0*x556
      x571 = -x570*x84
      x572 = x15*x49
      x573 = x559*x572
      x574 = x571 - 0.0165625d0*x573
      x575 = 0.0478125d0*x561 + x568*x569 + x574
      x576 = x15*x89
      x577 = 0.03346875d0*x82
      x578 = -x320*x570
      x579 = -0.046375d0*x573 + x578
      x580 = gamma_cc*(0.133875d0*x561 + x568*x577 + x579)
      x581 = 3.0d0*x576 + x580
      x582 = x30*x92
      x583 = x581*x582
      x584 = x87*x91
      x585 = x557 + x575*x95 - x583*x584
      x586 = x101*x585
      x587 = 3.54453756120869d-7*x220
      x588 = x210*x538
      x589 = x535*x540
      x590 = x589*(1.77226878060434d-7*x529*x532 + x537 - x588*(x202 + x
     &210*x587))
      x591 = x15*x590
      x592 = rho_c*x50
      x593 = -0.033125d0*x20 + 0.0165625d0*x47 + 0.0478125d0*x58 - 0.023
     &90625d0*x75
      x594 = -8.4d0*x85
      x595 = -0.09275d0*x20 + 0.046375d0*x47 + 0.133875d0*x58 - 0.066937
     &5d0*x75
      x596 = x595*x82
      x597 = x594 - x596 + 0.2008125d0*x83
      x598 = 1d0/tau_c
      x599 = x316*x598
      x600 = x597*x599
      x601 = gamma_cc*x87
      x602 = x593*x95 + x600*x601
      x603 = x1*x102
      x604 = x582*x603
      x605 = 0.682784063255296d0*x2
      x606 = x605*x79
      x607 = 0.0784052941866195d0*x236
      x608 = x3*x79
      x609 = x235*(2.9372772710225d0*x11 + 3.15315353825457d0*x3 + 0.910
     &227165231986d0*x8 + 3.09205213955878d0*x9)
      x610 = x232*x79
      x611 = x233**2/(x55**4*x56**2)
      x612 = x233*x234*(5.03533246461d0*x11 + 4.72973030738186d0*x3 + x5
     &4 + 1.09227259827838d0*x8)/x55**3
      x613 = x70**(-2)
      x614 = x613/x74**2
      x615 = rho_c**(-4.33333333333333d0)
      x616 = x35*x615
      x617 = 0.0449794051690726d0*x616
      x618 = rho_c**(-6.66666666666667d0)
      x619 = x23*x618
      x620 = x619*x68
      x621 = 0.0269876431014436d0*x252
      x622 = x69*x71
      x623 = x622*x68
      x624 = x260*x71
      x625 = x249*x64
      x626 = x244*x67
      x627 = 0.000312142662930383d0*x255 + 0.000133775426970164d0*x258*x
     &63
      x628 = x613*x72
      x629 = x260*x628*x64
      x630 = 0.0269876431014436d0*x626
      x631 = exp(4.0000000000000009d0*x58)
      x632 = x243**2
      x633 = x63*x632
      x634 = 0.0115661327577615d0*x631*x633
      x635 = x238*x61
      x636 = x632*x635
      x637 = 6.73880859569814d0*x16*x2
      x638 = 97.4090910340024d0*x31
      x639 = -0.313621176746478d0*x236*x637 + 0.209080784497652d0*x240 -
     & x241*x611* x638 - x242*x609 + x242*x612
      x640 = x635*x639
      x641 = rho_c**(-6.66666666666667d0)
      x642 = x23*x641
      x643 = x22*x619
      x644 = x61*x643
      x645 = x16*x264*x268*x30*x614*(x250 + x254 - x262) + x17*x610*x611
     & - 0.104540392248826d0*x230*x79 + x266*x268*(x239*x642*(x630 + x63
     &4 - x636 - x640 + 0.0449794051690726d0) + 0.0269876431014436d0*x24
     &6 *x620 - x260*x621*x623 + x261*(-0.000624285325860767d0*x244*x63*
     & x643 + x256*x639*x644 - x616*x630 - x616*x634 + x616*x636 + x616*
     & x640 - x617 - 0.000884404211636086d0*x620 + 6.68877134850821d-5* 
     &x633*x644 - 0.000200663140455246d0*x631*x632*x643/x62**3) + x617* 
     &x69 - 0.0115661327577615d0*x624*x625 + x629*(0.0115661327577615d0 
     &*x252*x626 + x621 + x627)) + x606*x607 + x607*x608 + x609*x610 - x
     &610*x612
      x646 = x645*x82
      x647 = x269*x271
      x648 = x81/rho_c**4
      x649 = x648*x78
      x650 = 0.172586075861435d0*x308
      x651 = 0.129439556896076d0*x278
      x652 = x605*x651
      x653 = x3*x651
      x654 = x277*(1.40779759833938d0*x11 + 1.44868929752102d0*x3 + 0.45
     &2125560575033d0*x8 + 1.06404211235975d0*x9)
      x655 = x4*x654
      x656 = 2.0d0*x655
      x657 = x274**2/(x12**4*x275**2)
      x658 = x4*x657
      x659 = x17*x658
      x660 = 2.0d0*x659
      x661 = x274*x276*(x10 + 2.41336731143894d0*x11 + 2.17303394628153d
     &0*x3 + 0.542550672690039d0*x8)/x12**3
      x662 = x4*x661
      x663 = 2.0d0*x662
      x664 = x16*x327
      x665 = x42**(-2)
      x666 = x665/x46**2
      x667 = x30*x50*x666
      x668 = 0.690344303445739d0*x273
      x669 = 1.03551645516861d0*x278*x637
      x670 = 4.0d0*x280
      x671 = x654*x670
      x672 = x638*x657
      x673 = 4.0d0*x279*x672
      x674 = x661*x670
      x675 = x668 - x669 - x671 - x673 + x674
      x676 = 0.000722883297360095d0*x282
      x677 = x675*x676
      x678 = 0.00674691077536089d0*x283
      x679 = exp(4.0d0*x20)
      x680 = x281**2
      x681 = x29*x680
      x682 = 0.00144576659472019d0*x679*x681
      x683 = x676*x680
      x684 = x678 + x682 - x683 + 0.0224897025845363d0
      x685 = 0.00289153318944038d0*x40*x642
      x686 = 0.0224897025845363d0*x616
      x687 = x40*x619
      x688 = 0.000221101052909021d0*x687
      x689 = x616*x678
      x690 = x616*x682
      x691 = x616*x683
      x692 = 7.80356657325958d-5*x292*x643
      x693 = 1.25414462784529d-5*x643*x679*x680/x28**3
      x694 = 4.18048209281763d-6*x643
      x695 = x27*x681*x694
      x696 = x616*x677
      x697 = x291*x675*x694
      x698 = 0.0134938215507218d0*x252
      x699 = x295*x43
      x700 = x40*x41
      x701 = x248*x285*x33
      x702 = 7.80356657325958d-5*x286 + 1.67219283712705d-5*x293
      x703 = x44*x665
      x704 = x295*x33*x703
      x705 = 0.0134938215507218d0*x285*x687 + x41*x686 - x698*x699*x700 
     &- 0.00578306637888076d0*x699*x701 + x704*(0.00289153318944038d0* x
     &252*x283 + x698 + x702)
      x706 = x299*x664*x667 + x302*(x296*(-x686 - x688 - x689 - x690 + x
     &691 - x692 - x693 + x695 + x696 + x697) + x685*(-x677 + x684) + x7
     &05) - x650* x79 + x652*x79 + x653*x79 + x656*x79 + x660*x79 - x663
     &*x79
      x707 = -x706*x84
      x708 = x271*x312
      x709 = x52*x648
      x710 = x707 + 0.033125d0*x708 - 0.0496875d0*x709
      x711 = x315*x90
      x712 = x324*x93
      x713 = -x321 + x322
      x714 = 6.0d0*x1
      x715 = x51*x79
      x716 = 12.0d0*x715
      x717 = x320*x706
      x718 = -0.09275d0*x708 + 0.139125d0*x709 + x717
      x719 = 1024.0d0*x79
      x720 = 66.2730531307909d0*x278
      x721 = -88.3640708410546d0*x308*x79 + x606*x720 + x608*x720 + x655
     &*x719 + x659* x719 - x662*x719
      x722 = -512.0d0*x16*x327**2*x667 + 512.0d0*x302*(-x296*(x686 + x68
     &8 + x689 + x690 - x691 + x692 + x693 - x695 - x696 - x697) + x685*
     &(x676*( -x668 + x669 + x671 + x673 - x674) + x684) + x705) + x721
      x723 = 0.67024923187228d0*x338*x99
      x724 = x100*x343
      x725 = mu*x336
      x726 = 0.446832821248187d0*rho_c**(-2.33333333333333d0)*x725
      x727 = 0.112308508206345d0*rho_c**(-2.66666666666667d0)*x104*x174
      x728 = 0.111534778586481d0*x143
      x729 = 0.455740507457808d0*x3
      x730 = 0.00197049465647014d0*x112 + x117 + 0.0259624262672375d0*x1
     &18*x8 + 0.00418159828868294d0*x3 - 0.151913502485936d0*x368 + 0.03
     &75734556469842d0*x8
      x731 = x109*(-x374*x729 + x730)
      x732 = x278*x385
      x733 = x16*x658
      x734 = -0.074356519057654d0*x144 + 4.21037858621281d0*x145*x305*x3
     &88 + x146* x733 + x228*x732 - 0.363326359042743d0*x384*x386 - 6.31
     &556787931922d0*x387 + x389*x655 - x389*x662 + 0.0691804894611514d0
     &*x732
      x735 = -x383*x731 + x734
      x736 = 0.263148661638301d0*x171
      x737 = 0.803635363211076d0*x120*x378
      x738 = 9.65778533132318d0*x379
      x739 = 0.00619637658813783d0*x187*x498
      x740 = 0.526297323276602d0*x115*x464
      x741 = 0.200908840802769d0*x143*x278
      x742 = 0.148394183600699d0*x278*x377
      x743 = 8.27810171256273d0*x111*x381
      x744 = x173*x380
      x745 = x655*x744
      x746 = 2.14502939711103d0*x173*x733
      x747 = x662*x744
      x748 = x376*x731
      x749 = 0.00464728244110338d0*x143
      x750 = x499*x749
      x751 = x181*(-0.105477574079378d0*x3 + 0.166966313186253d0*x8)
      x752 = 0.363424118566428d0*x160
      x753 = 14.5369647426571d0*x182*x362
      x754 = x157*x183*(-0.156676175559749d0*x3 + 0.620028561181324d0*x8
     &)
      x755 = 3.61230841029072d0*x352*x361
      x756 = 2.22044604925031d-16*x180*x352
      x757 = 2.42282745710952d0*x180*x350
      x758 = 0.0282745794067445d0*x347*x358
      x759 = 0.363424118566428d0*x166
      x760 = 2.18054471139857d0*x164*x359
      x761 = 1.63540853354893d0*x160*x164
      x762 = x356*x761*(0.219734258381682d0*x3 + 0.0582386976490866d0*x8
     &)
      x763 = x357*x761*(0.329601387572523d0*x3 + 0.0698864371789039d0*x8
     &)/x161**3
      x764 = x758 - x759 + x760 - x762 + x763
      x765 = x366*(x751 + x752 + x753 + x754 + x755 - x756 - x757 + x764
     &)
      x766 = x397*(-1.99883265211535d0*x3 + 0.917202013581841d0*x8)
      x767 = x152*x396
      x768 = 24.0d0*x154*x767
      x769 = 0.154912426780983d0*x1**(-0.333333333333333d0)*x153
      x770 = 2.25322953511185d0*x398*x767
      x771 = 0.751076511703951d0*x399
      x772 = 0.350864882184401d0*x149
      x773 = 0.0799516697279215d0*x407
      x774 = 0.00495710127051027d0*x405
      x775 = 0.0037178259528827d0*x121*x403
      x776 = 0.00910930363347549d0*x122
      x777 = 0.00394098931294027d0*x112 - 0.01915695d0*x115 + 0.00836319
     &657736588d0* x3 + 0.0751469112939684d0*x8
      x778 = x395*x777
      x779 = x119*x772
      x780 = 0.133252782879869d0*x408
      x781 = x774 - x775 + x776 + x778 + x779 - x780
      x782 = 1.08351503479231d0*x11
      x783 = 3.00815243991335d0*x413
      x784 = 2.54445191835098d0*x415
      x785 = 2.0d0*x412 + 2.0d0*x414 + 2.0d0*x416
      x786 = mu*x417
      x787 = x129*x786
      x788 = x785*x787
      x789 = x132*(1.39252411901946d0*x413 + x782)
      x790 = x131*x786*(1.85745434535825d0*x11 + 2.08878617852919d0*x413
     &)/x128**2
      x791 = x419*x51
      x792 = mu*x418*x785/x131**2
      x793 = x123*x423
      x794 = 0.0101661833184394d0*x433
      x795 = 0.074356519057654d0*x121
      x796 = x278*x427
      x797 = x436*x438
      x798 = -4.52089344419912d-5*x1**2.33333333333333d0*x434 - x109*x47
     &1*x777 - 0.0577480492749442d0*x425*x426*x428 + 0.713822582953479d0
     &*x429* x437 - 1.30867473541471d0*x430 - 0.0991420254102054d0*x432 
     &+ 0.023721094409692d0*x435 + 0.466194077035412d0*x438*x733 + x655*
     & x797 - x662*x797 + 0.064503068866399d0*x796
      x799 = x421*x787 - x51*x792 + x793*x794 - x793*x795 + x798
      x800 = 0.00327431399373698d0*x199
      x801 = x615*x800
      x802 = 0.00196458839624219d0*x199
      x803 = 1d0/x205
      x804 = x475*x803
      x805 = x802*x804
      x806 = x536*x615
      x807 = x475**2
      x808 = x200/x205**2
      x809 = x807*x808
      x810 = x193*x2
      x811 = 0.156226553915631d0*x198*x810
      x812 = x533*x804
      x813 = x530*x615*x812
      x814 = 0.000420983227766184d0*x485
      x815 = x457*x814
      x816 = -x306 + x309
      x817 = -x501 + x503 + x816
      x818 = x475*x808*x817
      x819 = x210*x618
      x820 = x197*(x481 + x509 + x512*x814)
      x821 = x530*x820
      x822 = x819*x821
      x823 = 0.000841966455532367d0*x199
      x824 = 0.359841968171091d0*x453
      x825 = x17*x457
      x826 = 5.56d0*x825
      x827 = -x758 + x759 - x760 + x762 - x763
      x828 = x109*(x463*x729 + x730)
      x829 = -x736 - x737 - x738 + x740 + x741 + x742 + x743 + x745 + x7
     &46 - x747
      x830 = x464*x728
      x831 = x470*x773
      x832 = x470*x772
      x833 = x168*(-x766 - x768 - x769 - x770 - x771 + x827)
      x834 = x782 + x783 + x784 - x788 - x789 + x790
      x835 = x431*x470
      x836 = x344*(x108*(-x383*x828 + x734 + x830) + x130*(0.00619637658
     &813783d0*x188 + x366*(-x751 - x752 - x753 - x754 - x755 + x756 + x
     &757 + x827) - x376*x828 - x462*x749 + x829) + x134*x419*x834 + x13
     &4*x792 + x411 *(-x774 + x775 - x776 - x778 - x779 + x780 - x831 + 
     &x832 + x833) - x469*x787 - x794*x835 + x795*x835 + x798)
      x837 = x192*x8/x4**2
      x838 = 0.359841968171091d0*x278
      x839 = x112*x170/x106**6
      x840 = x191*x839
      x841 = x130*x449*x8
      x842 = x472*x841
      x843 = 2.78d0*x194
      x844 = x456*x654
      x845 = x192*x672
      x846 = x19**(-2)
      x847 = 5.56d0*x845
      x848 = x456*x661
      x849 = x196*x637
      x850 = x17*x447*x657 + x3*x838 - 0.479789290894789d0*x308 + x447*x
     &654 - x447* x661 + 6.28423862597952d0*x451 + 0.239894645447394d0*x
     &454 + x605* x838 - 2.78d0*x836 - 0.023288892456339d0*x837 + x838*x
     &849 - 6.39253410025677d0*x840 - 7.54108635117542d0*x842 - x843*x84
     &4 - x843*x845 + x843*x848 - x846*x847
      x851 = x200*x803*(-x484*x824 + x484*x826 + x850)
      x852 = x251*x530
      x853 = x802*x852
      x854 = x488*x528
      x855 = x531*x854
      x856 = x524*x531
      x857 = x814*x856
      x858 = x505*x524*x526*x817*x854
      x859 = x804*x852
      x860 = x218*x494*x803*x855
      x861 = x200*x859 + x491*x855 + x508*x530 - x817*x852*x857 - 3.5445
     &3756120869d-7 *x858 + 1.77226878060434d-7*x860
      x862 = x535*x861
      x863 = x532*x862
      x864 = x533*x863
      x865 = x209*x30
      x866 = 1.65411752856405d-6*x218
      x867 = x817*x856
      x868 = 0.000841966455532367d0*x485*x867
      x869 = 0.00196458839624219d0*x485
      x870 = 818029.852187331d0*x15**(-5.78d0)
      x871 = x817**2*x870
      x872 = 2.29947269793409d0*x204*x839
      x873 = 2.26051749135954d0*x502
      x874 = 2.71262098963145d0*x500*x841
      x875 = x344*(x108*(x735 + x830) + x130*(-x739 - x748 + x750 - x765
     & + x829) - x411*(x781 + x831 - x832 - x833) - x791*x834 + x799)
      x876 = x216*x819
      x877 = x218*x876
      x878 = x197*x615
      x879 = x207*x215
      x880 = gamma_cc*x216*x618
      x881 = x506*x879*x880
      x882 = x504**2
      x883 = x870*x878
      x884 = x882*exp(2.0d0*x206)
      x885 = 6869716491.20552d0*x15**(-9.56d0)*x880
      x886 = x476*x878
      x887 = x192**(-2)
      x888 = x200*x878*x887
      x889 = 1.77226878060434d-7*x877
      x890 = x887*x889
      x891 = x495*x876
      x892 = x650 - x652 - x653 - x656 - x660 + x663 - x872 + x873 + x87
     &4 + x875
      x893 = -x473 - x482 + x816
      x894 = 0.719683936342183d0*x278
      x895 = 5.56d0*x194
      x896 = 512.0d0*x51
      x897 = x518/x539**2
      x898 = x224*x30
      x899 = gamma_cc*x589*x896*(x200*x528*x532*x641*(0.0008419664555323
     &67d0*x526* x871*exp(2.0d0*x523) + x800 + x804*x811 - x804*x815 - x
     &804*x868 + x805 + x809 - x818 - x851 - 0.000420983227766184d0*x856
     &*x871 + x857*(-x650 + x652 + x653 + x656 + x660 - x663 + x872 - x8
     &73 - x874 - x875) - x867*x869) - x247*x531*x821*x823*x862*x865 + x
     &30* x526*x537*x861*(x823*x859 - x852*x868 + x853 + x855*x866 - 7.0
     &8907512241738d-7*x858 + 3.54453756120869d-7*x860)/x534**2 + x536*x
     &801 - x588*(x197*x801 + 0.000841966455532367d0*x215*x883* x884 + 0
     &.552318141716643d0*x217*x810*x891 - x445*x889*(-x194*x847 + x3*x89
     &4 - 0.959578581789577d0*x308 + 12.568477251959d0*x451 + 0.71968393
     &6342183d0*x453*x893 + 0.479789290894789d0*x454 + x605* x894 + 11.1
     &2d0*x655 + 11.12d0*x659 - 11.12d0*x662 - 11.12d0*x825* x893 - 5.56
     &d0*x836 - 0.046577784912678d0*x837 - 12.7850682005135d0*x840 - 15.
     &0821727023508d0*x842 - x844*x895 - 11.12d0*x845*x846 + x848*x895 +
     & x849*x894) - 1.77226878060434d-7* x457*x505*x891 + x475*x504*x888
     & - x486*x878*(x824*x893 - x826* x893 + x850) + x494**2*x890 + x494
     &*x504*x890 + 1.65411752856405d-6*x495*x877 + 7.08907512241738d-7*x
     &495*x881 - x497*x505*x880*x892 - x497*x882*x885 - x511*x814*x878*x
     &892 - 0.000420983227766184d0*x511*x882*x883 + x512*x869*x878 + 0.0
     &00841966455532367d0*x513*x886 + x802*x886 + x807*x888 + x811* x886
     & - x815*x886 + 4.68666633093149d-6*x877 + 3.30823505712811d-6 *x88
     &1 + 1.06336126836261d-6*x884*x885/x208**3) + x802*x822 + x804 *x82
     &2*x823 + x805*x806 + x806*x809 - x806*x818 - x806*x851 + x811 *x81
     &3 - x812*x823*x852*x863 - x813*x815 - x853*x864) - 7145614.4266171
     &5d0*x15**(-6.56d0)*x530*x535*x897*x898*( -0.000420983227766184d0*x
     &36*x864 + x443*x533 + x444*x812 + x488* x516*x820) - x345*x79*(x10
     &8*(-x375*x728 + x735) - x130*(x736 + x737 + x738 + x739 - x740 - x
     &741 - x742 - x743 - x745 - x746 + x747 + x748 - x750 + x765) - x41
     &1*(x168*(x764 + x766 + x768 + x769 + x770 + x771) + x422*x772 - x4
     &22*x773 + x781) + x791*(-x782 - x783 - x784 + x788 + x789 - x790) 
     &+ x799) - 1388.8619466913d0* x441*x79*x841 + 1177.33002134225d0*x5
     &21*x79*x839 - 1157.38495557608d0*x543*x79 + x721
      x900 = 0.00130908053100055d0*x339*x725
      x901 = -4.0d0*x329 + 0.258879113792152d0*x331 - 2.0d0*x344*x442 - 
     &2.0d0*x517* x542 + 2.71262098963145d0*x544
      x902 = x51*x590
      x903 = x269*x559
      x904 = gamma_cc*x90
      x905 = x80*x904
      x906 = x78*x905
      x907 = -x562 + x566 - x72
      x908 = x565*x71
      x909 = x68*x908
      x910 = 0.0115661327577615d0*x37
      x911 = 0.000133775426970164d0*x25*x68 + x910
      x912 = x265*(x16*x267*x907 - 0.0115661327577615d0*x246*x255 + x254
     &*x909 - 7.80356657325958d-5*x255 + x261*(x259 + x627) + x263 + x62
     &4*x65 + 0.00578306637888076d0*x625*x908 - x629*x911)
      x913 = x50*x912
      x914 = x271*x50
      x915 = x568*x914
      x916 = x312*x559
      x917 = -x44 - x548 + x554
      x918 = 1.9508916433149d-5*x286
      x919 = 0.00578306637888076d0*x287
      x920 = x34*x699
      x921 = x40*x552
      x922 = x290*x921
      x923 = x296*(x294 + x702)
      x924 = 0.00289153318944038d0*x552*x701
      x925 = 0.00578306637888076d0*x37
      x926 = x25*x40
      x927 = x925 + 3.34438567425411d-5*x926
      x928 = x704*x927
      x929 = x298 - x918 - x919 + x920 + x922 + x923 + x924 - x928
      x930 = x22*x299*x300*x917 + x929
      x931 = 1.0d0*x301
      x932 = x930*x931
      x933 = x52*x905
      x934 = x556*x914
      x935 = -x84*x932 - 0.0165625d0*x916 + 0.033125d0*x933 + 0.0165625d
     &0*x934
      x936 = x50*x90
      x937 = x575*x936
      x938 = -x320*x932 - 0.046375d0*x916 + 0.09275d0*x933 + 0.046375d0*
     &x934
      x939 = 3.0d0*x580
      x940 = x1*x50
      x941 = 9.0d0*x317
      x942 = 0.2008125d0*x270
      x943 = 0.401625d0*x272
      x944 = 3.0d0*x321
      x945 = 0.139125d0*x313
      x946 = x944 - x945
      x947 = 512.0d0*x16
      x948 = x666*x947
      x949 = gamma_cc*x327*x50*x555*x948 - 512.0d0*x301*(x326 + x918 + x
     &919 - x920 - x922 - x923 - x924 + x928)
      x950 = x212*x508
      x951 = x212*x510
      x952 = gamma_cc*x479
      x953 = x202 + x587*x952
      x954 = x952*x953
      x955 = x224*x954
      x956 = x478**(-2)
      x957 = 1.77226878060434d-7*x220
      x958 = x214*x224
      x959 = x210*x958
      x960 = x223**(-2)
      x961 = x201*x823
      x962 = x220*x952
      x963 = x961 + 7.08907512241738d-7*x962
      x964 = x210*x957 + x214 - x954*x958
      x965 = gamma_cc*x519*x520*x897*x960*x964 - x541*x896*(-x200*x216*x
     &247*x487*x898 *x953*x956 + x214*x479*x515*x865*x960*x963 - x215*x5
     &15*x898*x957 + 1.77226878060434d-7*x476*x493 + x479*x490*x823 + 4.
     &13529382141014d-7*x493 - x515*x959 - x950*x955 + x950 - x951* x955
     & + x951 - x959*(x492*x866 + 3.54453756120869d-7*x496 + 7.089075122
     &41738d-7*x507*x879 + x514))
      x966 = x558*x560
      x967 = x559*x568
      x968 = x16*x569
      x969 = rho_c**(-7.0d0)
      x970 = 3020.29322777679d0*x969
      x971 = x623*x970
      x972 = x565*x622
      x973 = x565*x68*x911
      x974 = x265*x67
      x975 = x974*(x563 - x564*x908 + x628*x973 - x910*x972 - 3.86816087
     &015788d-7* x971)
      x976 = x567*x614*x907
      x977 = x558*x572
      x978 = x556*x559
      x979 = 1.0d0*x16
      x980 = x84*x979
      x981 = x43*x700*x970
      x982 = 4.83520108769735d-8*x981
      x983 = x550*x552
      x984 = x41*x552*x925
      x985 = x40*x551*x927
      x986 = x703*x985
      x987 = x300*x39
      x988 = x987*(x549 - x982 - x983 - x984 + x986)
      x989 = x555*x666*x917
      x990 = -x980*x988 - x980*x989
      x991 = -0.0165625d0*x977 - 0.033125d0*x978 + x990
      x992 = x16*x577
      x993 = x320*x979
      x994 = -x988*x993 - x989*x993
      x995 = -0.046375d0*x977 - 0.09275d0*x978 + x994
      x996 = x904*x92
      x997 = x51*x996
      x998 = x555**2*x948 + x947*x987*(-x549 + x982 + x983 + x984 - x986
     &)
      x999 = 3.54453756120869d-7*x220
      x1000 = x202*x952 + 1.0d0
      x1001 = x1000 + x222*x956
      x1002 = 1d0/x1001
      x1003 = x1002*x953
      x1004 = x1003*x952
      x1005 = x1000*x202
      x1006 = x1002*x15
      x1007 = x1006*x224*x964*(-x1004*x1005 + x1005 + 1.77226878060434d-
     &7*x962)/x226** 2 + x1006*x479*(88.3758007628711d0*x1000*x1002*x15*
     &*(-11.34d0)* x196**8.34d0*x952*x969 + x1000*x1003*x961 + x1004*x99
     &9 - x999 - x1005*x954*x963/x1001**2)/x226
      x1008 = rho_c**(-3.33333333333333d0)*x15*x337
      x1009 = x1*x78
      x1010 = x237*x53
      x1011 = x267*x50
      x1012 = x1*x52
      x1013 = -0.033125d0*x1012 + 0.0165625d0*x303 - 0.033125d0*x329 + 0
     &.00214384266109126d0*x331
      x1014 = x593*x94
      x1015 = -0.09275d0*x1012 + 0.046375d0*x303 - 0.09275d0*x329 + 0.00
     &600275945105553d0*x331
      x1016 = x271*x51
      x1017 = 3.0d0*x1016
      x1018 = -x944 + x945
      x1019 = x599*x601
      x1020 = x101*x582
      x1021 = x15*(x59 - x76)
      x1022 = x21 - x48
      x1023 = x1022*x15
      x1024 = x266*x567
      x1025 = gamma_cc*x556
      x1026 = 0.0165625d0*x1025
      x1027 = 3.0d0*x15
      x1028 = 0.046375d0*x1025
      x1029 = gamma_cc*x101*x92*x940
      x1030 = 33.6d0*x85
      x1031 = x30*x603/tau_c**4
      x1032 = x49 + x77
      x1033 = x1032*x82
      x1034 = 0.011953125d0*x1033 + x86
      x1035 = 0.03346875d0*x1033 + x88
      x1036 = -x1035*x94 + 512.0d0
      x1037 = x1034*x1036
      x1038 = -rho_c*(-x102*(-x1037 + x98) + x227)
      x1039 = x20*x51
      x1040 = x47*x51
      x1041 = 0.5d0*rho_c
      x1042 = -1.0d0*x1039 + 0.5d0*x1040 + x1041*x269 + x1041*x312 - 0.5
     &d0*x51*x58 + 0.25d0*x51*x75
      x1043 = x1042*x271
      x1044 = x1016*x1032
      x1045 = 0.02390625d0*x1043 - 0.035859375d0*x1044 + x314
      x1046 = x1035*x316
      x1047 = 0.0669375d0*x1043
      x1048 = 0.10040625d0*x1044
      x1049 = 3.0d0*x1046 - x1047 + x1048 + x323
      x1050 = x1034*x94
      x1051 = x100*(x1036*x1045 + x1049*x1050 + x334)
      x1052 = x1037 + x341
      x1053 = x1052*x15
      x1054 = -x102*x1052 + x546*(-x1051 - x1053*x340 + x545) + x547
      x1055 = 0.00144576659472019d0*x37
      x1056 = x1055*x69
      x1057 = x1055*x41
      x1058 = x265*(-x1056*x909 + x1056 + x26*x68) + x300*(-x1057*x921 +
     & x1057 + 4.18048209281763d-6*x926)
      x1059 = x1058*x82
      x1060 = x1032*x15
      x1061 = x1060*x559
      x1062 = 0.02390625d0*x1059 + 0.02390625d0*x1061 + x574
      x1063 = x1035*x15
      x1064 = gamma_cc*(0.0669375d0*x1059 + 0.0669375d0*x1061 + x579)
      x1065 = 3.0d0*x1063 + x1064
      x1066 = x1065*x582
      x1067 = x1034*x91
      x1068 = x1036*x1062 - x1066*x1067 + x557
      x1069 = x101*x1068
      x1070 = x592*(-x1069 + x591)
      x1071 = 0.02390625d0*x58
      x1072 = 0.011953125d0*x75
      x1073 = x1071 - x1072 + 0.0146875d0*x20 - 0.00734375d0*x47
      x1074 = 0.0669375d0*x58
      x1075 = 0.03346875d0*x75
      x1076 = x1074 - x1075 + 0.041125d0*x20 - 0.0205625d0*x47
      x1077 = x1076*x82
      x1078 = 0.10040625d0*x1033 - x1077 + x594
      x1079 = x1078*x599
      x1080 = gamma_cc*x1034
      x1081 = x1036*x1073 + x1079*x1080
      x1082 = -x1081*x604
      x1083 = x271*(x1041*x645 + x1041*x706 + x269 + x304 + x310 - 2.0d0
     &*x329)
      x1084 = x1042*x648
      x1085 = x1032*x51
      x1086 = x1085*x648
      x1087 = x1049*x93
      x1088 = 2.0d0*x90
      x1089 = x100*x1053
      x1090 = -0.00390625d0*x1051 - x1053*x900 + x546*(-x100*(x1036*(0.0
     &2390625d0* x1083 - 0.1434375d0*x1084 + 0.1434375d0*x1086 + x710) +
     & x1045* x1087*x1088 + x1050*(-x1035*x716 - 0.0669375d0*x1083 + 0.4
     &01625d0 *x1084 - 0.401625d0*x1086 + x714*(x1047 - x1048 + x713) + 
     &x718) + x722) - x1051*x723 + x1089*x726 - x1089*x727 + x899) + x90
     &1
      x1091 = x1042*x905
      x1092 = x914*(0.25d0*rho_c*x912 + x1041*x300*x930 + x1058)
      x1093 = x1058*x914
      x1094 = x1085*x905
      x1095 = x1062*x936
      x1096 = x582*x936
      x1097 = 3.0d0*x1064
      x1098 = 9.0d0*x1046
      x1099 = 0.2008125d0*x1043
      x1100 = 0.30121875d0*x1044
      x1101 = -x1069*x50 + x546*(-x100*(x1036*(0.0478125d0*x1091 + 0.023
     &90625d0*x1092 - 0.07171875d0*x1093 - 0.07171875d0*x1094 + x935) - 
     &x1045*x1065* x1096 + x1067*x582*(-gamma_cc*(0.133875d0*x1091 + 0.0
     &669375d0* x1092 - 0.2008125d0*x1093 - 0.2008125d0*x1094 + x938) + 
     &x1097* x940 + x1098 - x1099 + x1100 + x946) + x1087*x1095 + x949) 
     &- x1068*x340 + x965) + x902
      x1102 = x1060*x558
      x1103 = gamma_cc*x1058
      x1104 = x1103*x558
      x1105 = 0.02390625d0*x82
      x1106 = x16*(x974*(x1056*x613*x973 - x38*x972 - x549*x909 + x549 -
     & 9.6704021753947d-8*x971) + 0.25d0*x976 + x987*(x1057*x665*x985 + 
     &x26 - x548*x552 - x553 - 2.41760054384867d-8*x981) + 0.5d0*x989)
      x1107 = 0.0669375d0*x82
      x1108 = -x592*(x1007 - x101*(x1034*x997*(6.0d0*x1063 + 6.0d0*x1064
     & + x30*( 0.0669375d0*x1102 + 0.26775d0*x1104 + x1106*x1107 + x995)
     &) - x1036*(0.02390625d0*x1102 + 0.095625d0*x1104 + x1105*x1106 + x
     &991 ) + 2.0d0*x1066*x1095 + x998))
      x1109 = 0.02390625d0*rho_c
      x1110 = x1032*x316
      x1111 = x1073*x94
      x1112 = 0.0669375d0*rho_c
      x1113 = x1080*x599
      x1114 = -x1020*(rho_c*(x1008*x1081 + x79*(gamma_cc*x1045*x1079 + x
     &1036*(-x1*( -0.0478125d0*x1039 + 0.02390625d0*x1040 - x1071*x51 + 
     &x1072*x51 + x1109*x269 + x1109*x312) + x1013 + 0.07171875d0*x1110)
     & + x1049* x1111 + x1113*(x1017*x1076 + x1018 - x1098 + x1099 - x11
     &00 - x82* (-x1*(-0.133875d0*x1039 + 0.0669375d0*x1040 - x1074*x51 
     &+ x1075* x51 + x1112*x269 + x1112*x312) + x1015 + 0.2008125d0*x111
     &0)))) + x1081*x715)
      x1115 = -x1029*(x1036*(x1026 - 0.0478125d0*x1060 - 0.0478125d0*x11
     &03 + 0.033125d0*x572) + x1062*x1079*x30 - x1065*x1111 + x1113*(-x1
     &027* x1077 + 9.0d0*x1063 + x1097 - x82*(x1028 - 0.133875d0*x1060 -
     & 0.133875d0*x1103 + 0.09275d0*x572)))
      x1116 = -x1031*(-x1036*(0.0440625d0*x20 - 0.02203125d0*x47 + 0.071
     &71875d0*x58 - 0.035859375d0*x75) + 2.0d0*x1078*x1111 + x1113*(x103
     &0 - 0.401625d0*x1033 + 6.0d0*x1077 + x82*(0.123375d0*x20 - 0.06168
     &75d0*x47 + 0.2008125d0*x58 - 0.10040625d0*x75)))
      x1117 = 0.015625d0*x82
      x1118 = x1117 - 1.0d0
      x1119 = 0.04375d0*x82
      x1120 = x1119 - 2.8d0
      x1121 = x1120*x52
      x1122 = x90*x93
      x1123 = x1121*x1122
      x1124 = -x1123 + 512.0d0
      x1125 = x1118*x1124*x49
      x1126 = x1016*x1022
      x1127 = 0.03125d0*x1126
      x1128 = x311 - 1.0d0*x328
      x1129 = x1128*x82
      x1130 = 0.02390625d0*x1129
      x1131 = x1128*x84
      x1132 = x1022*x316
      x1133 = x1128*x320
      x1134 = -0.0875d0*x1126 - 0.0669375d0*x1129 + x1133
      x1135 = x1022*x94
      x1136 = x100*(-x1135*(-x1117 + 1.0d0)*(3.0d0*x1132*(-x1119 + 2.8d0
     &) + x1134) + 512.0d0*x328 + x330 - x332 + (x1123 - 512.0d0)*(x1127
     & + x1130 - x1131))
      x1137 = x1125 + x341
      x1138 = x1137*x15
      x1139 = x1105*x556 + x571 + 0.03125d0*x573
      x1140 = x1120*x572
      x1141 = gamma_cc*(x1107*x556 + 0.0875d0*x573 + x578)
      x1142 = 3.0d0*x1140 + x1141
      x1143 = x1118*x52
      x1144 = x1143*x582*x90
      x1145 = x1124*x1139 - x1142*x1144 + x557
      x1146 = x101*x1145
      x1147 = x1022*x1124
      x1148 = 0.175d0*x20 - 0.0875d0*x47
      x1149 = x1148*x82
      x1150 = x1120*x49
      x1151 = -x1149 + 3.0d0*x1150
      x1152 = gamma_cc*x1012*x1118*x598
      x1153 = 0.03125d0*x1147 + x1151*x1152
      x1154 = x706*x82
      x1155 = -x1127 - x1130 + x1131
      x1156 = x1012*x1120
      x1157 = x1134 + 3.0d0*x1156
      x1158 = x1157*x93
      x1159 = x100*x1138
      x1160 = 0.0625d0*x1022
      x1161 = x51*x905
      x1162 = x1128*x559
      x1163 = x1025*x664 + x929
      x1164 = x1163*x301
      x1165 = x1163*x931
      x1166 = x1096*x1142
      x1167 = 0.2625d0*x1126
      x1168 = 9.0d0*x1156
      x1169 = 0.2008125d0*x1129
      x1170 = 3.0d0*x1133
      x1171 = 3.0d0*x1141
      x1172 = x1105*x16
      x1173 = x1107*x16
      x1174 = x1151*x599
      x1175 = 0.03125d0*x1135
      PBE_value = x52
      PBE_beta_value = x78
      PBE_alpha_value = x78
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(-x102*(-x96 + x98) + x227)
      d1E(1) = -x102*x342 + x546*(-x335 - x340*x343 + x545) + x547
      d1E(3) = x592*(-x586 + x591)
      d1E(6) = -x602*x604
      d2E(1) = -0.00390625d0*x335 - x343*x900 + x546*(-x100*(x325*(-0.06
     &69375d0*x646 + 0.26775d0*x647 - 0.401625d0*x649 + x714*(x318 - x31
     &9 + x713) - x716*x89 + x718) + 2.0d0*x711*x712 + x722 + x95*(0.023
     &90625d0* x646 - 0.095625d0*x647 + 0.1434375d0*x649 + x710)) - x335
     &*x723 + x724*x726 - x724*x727 + x899) + x901
      d2E(4) = -x50*x586 + x546*(-x100*(-x50*x583*x711 + x582*x584*(-gam
     &ma_cc*(x577* x913 + 0.133875d0*x903 - 0.26775d0*x906 - 0.0669375d0
     &*x915 + x938 ) + x939*x940 + x941 - x942 + x943 + x946) + x712*x93
     &7 + x949 + x95*(x569*x913 + 0.0478125d0*x903 - 0.095625d0*x906 - 0
     &.02390625d0*x915 + x935)) - x340*x585 + x965) + x902
      d2E(6) = -x592*(x1007 - x101*(2.0d0*x583*x937 + x87*x997*(x30*(0.1
     &33875d0*x966 + 0.133875d0*x967 + x975*x992 + x976*x992 + x995) + 6
     &.0d0*x576 + 6.0d0*x580) - x95*(0.0478125d0*x966 + 0.0478125d0*x967
     & + x968* x975 + x968*x976 + x991) + x998))
      d2E(16) = -x1020*(rho_c*(x1008*x602 + x79*(gamma_cc*x315*x600 + x1
     &014*x324 + x1019 *(x1017*x595 + x1018 - x82*(0.26775d0*x1009 + 0.1
     &33875d0*x1010 - 0.0669375d0*x1011 + x1015 - 0.0104965087592337d0*x
     &231) - x941 + x942 - x943) + x95*(0.095625d0*x1009 + 0.0478125d0*x
     &1010 - 0.02390625d0*x1011 + x1013 - 0.00374875312829775d0*x231))) 
     &+ x602 *x715)
      d2E(18) = -x1029*(-x1014*x581 + x1019*(-x1027*x596 + 9.0d0*x576 + 
     &x82*(-0.26775d0* x1021 + 0.09275d0*x1023 + 0.0669375d0*x1024 - x10
     &28) + x939) + x30*x575*x600 - x95*(-0.095625d0*x1021 + 0.033125d0*
     &x1023 + 0.02390625d0*x1024 - x1026))
      d2E(21) = -x1031*(2.0d0*x1014*x597 + x1019*(x1030 + 6.0d0*x596 + x
     &82*(-0.27825d0* x20 + 0.139125d0*x47 + 0.401625d0*x58 - 0.2008125d
     &0*x75) - 0.80325d0*x83) - x95*(-0.099375d0*x20 + 0.0496875d0*x47 +
     & 0.1434375d0*x58 - 0.07171875d0*x75))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = x1038
      d1E(1) = x1054
      d1E(3) = x1070
      d1E(6) = x1082
      d2E(1) = x1090
      d2E(4) = x1101
      d2E(6) = x1108
      d2E(16) = x1114
      d2E(18) = x1115
      d2E(21) = x1116
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = x1038
      d1E(1) = x1054
      d1E(3) = x1070
      d1E(6) = x1082
      d2E(1) = x1090
      d2E(4) = x1101
      d2E(6) = x1108
      d2E(16) = x1114
      d2E(18) = x1115
      d2E(21) = x1116
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(-x102*(-x1125 + x98) + x227)
      d1E(1) = -x102*x1137 + x546*(x1136 - x1138*x340 + x545) + x547
      d1E(3) = x592*(-x1146 + x591)
      d1E(6) = -x1153*x604
      d2E(1) = 0.00390625d0*x1136 - x1138*x900 + x546*(-x100*(-x1088*x11
     &55*x1158 - x1122*x1143*(12.0d0*x1121*x79 + x1134*x714 + 0.0669375d
     &0*x1154 - 0.175d0*x708 + 0.2625d0*x709 - x717) + x1124*(0.02390625
     &d0*x1154 + x707 - 0.0625d0*x708 + 0.09375d0*x709) + x722) + x1136*
     &x723 + x1159*x726 - x1159*x727 + x899) + x901
      d2E(4) = -x1146*x50 + x546*(-x100*(-x1124*(-x1105*x1164 - x1160*x1
     &161 - 0.03125d0 *x1162 + x1165*x84 + 0.03125d0*x934) + x1139*x1158
     &*x936 + x1144*( gamma_cc*(-0.175d0*x1022*x1161 - x1107*x1164 - 0.0
     &875d0*x1162 + x1165*x320 + 0.0875d0*x934) - x1167 + x1168 - x1169 
     &+ x1170 + x1171*x940) + x1155*x1166 + x949) - x1145*x340 + x965) +
     & x902
      d2E(6) = -x592*(x1007 - x101*(-x1124*(x1172*x988 + x1172*x989 + 0.
     &03125d0*x977 + 0.0625d0*x978 + x990) + 2.0d0*x1139*x1166 + x1143*x
     &996*(6.0d0* x1140 + 6.0d0*x1141 + x30*(x1173*x988 + x1173*x989 + 0
     &.0875d0* x977 + 0.175d0*x978 + x994)) + x998))
      d2E(16) = -x1020*(rho_c*(x1008*x1153 + x79*(-gamma_cc*x1155*x1174 
     &- x1124*( 0.0625d0*x1132 - 0.03125d0*x328 - 0.0625d0*x329 + 0.0040
     &4498615300238d0*x331) + x1152*(x1017*x1148 + x1167 - x1168 + x1169
     & - x1170 + x82*(0.175d0*x1132 - 0.0875d0*x328 - 0.175d0* x329 + 0.
     &0113259612284067d0*x331)) + x1157*x1175)) + x1153*x715)
      d2E(18) = -x1029*(-x1124*(-0.0625d0*x1023 + 0.03125d0*x1025) + x11
     &39*x1174*x30 - x1142*x1175 + x1152*(-x1027*x1149 + 9.0d0*x1140 + x
     &1171 + x82*( -0.175d0*x1023 + 0.0875d0*x1025)))
      d2E(21) = -x1031*(-0.09375d0*x1147 + x1151*x1160*x94 + x1152*(6.0d
     &0*x1149 - 12.0d0 *x1150 + x82*(0.525d0*x20 - 0.2625d0*x47)))
      end if
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, gam
     &ma_ss, gamma_cs, tau_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, gamma_ss, gamma_cs, 
     &tau_c, mu
      real*8, intent(out) :: E
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
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
     &
      E = 0.0d0
      x0 = 0.682784063255296d0
      x1 = 1d0/rho_c
      x2 = x0*x1**0.333333333333333d0
      x3 = 2.14502939711103d0
      x4 = 1d0/x3
      x5 = x1**0.666666666666667d0
      x6 = x4*x5
      x7 = 0.564189583547756d0
      x8 = sqrt(x1)*x7
      x9 = 0.826307487110758d0
      x10 = x1**0.166666666666667d0*x9
      x11 = 9.86960440108936d0
      x12 = 0.693147180559945d0
      x13 = -x12 + 1.0d0
      x14 = 1d0/x13
      x15 = x11*x14
      x16 = (0.186690969707574d0*x2 + 1.0d0)*log(x15/(13.457913714394447
     &d0*x10 + 5.630984149097876d0*x2 + 0.51606646454786342d0*x6 + 2.915
     &2147142191778d0*x8) + 1.0d0)
      x17 = (0.194159335344114d0*x2 + 1.0d0)*log(0.5d0*x15/(7.2401019343
     &168302d0*x10 + 3.2595509194222925d0*x2 + 0.40691300451752932d0*x6 
     &+ 1.4187228164796675d0*x8) + 1.0d0)
      x18 = 2.0d0*x17
      x19 = -x16 + x18
      x20 = rho_s*x1
      x21 = x20 + 1.0d0
      x22 = x21**1.33333333333333d0
      x23 = -x20 + 1.0d0
      x24 = x23**1.33333333333333d0
      x25 = 1.92366105093154d0*x22 + 1.92366105093154d0*x24 - 3.84732210
     &186307d0
      x26 = rho_s**4/rho_c**4
      x27 = 1d0/x11
      x28 = x13*x27
      x29 = x25*x26*x28
      x30 = x19*x29
      x31 = x18*x28
      x32 = x25*(0.101077332976288d0*x2 + 1.0d0)*(-x26 + 1.0d0)*log(1.0d
     &0 + 29.608574643216677d0/(9.8721297225692712d0*x10 + 3.29180480994
     &50625d0*x2 + 0.41002507094961249d0*x6 + 0.76232752193529008d0*x8))
     &
      x33 = 0.0197517897025652d0*x32
      x34 = -x31 + x33
      x35 = x30 + x34
      x36 = x21**0.666666666666667d0 + x23**0.666666666666667d0
      x37 = x36**3
      x38 = 8.0d0*x15/x37
      x39 = 1.0d0 - exp(-x35*x38)
      x40 = rho_c**(-4.66666666666667d0)
      x41 = gamma_cc**2
      x42 = 208.945363813799d0/x13**2
      x43 = 0.000133775426970164d0*x40*x41*x42/x36**4
      x44 = 14.454942539277d0*x14
      x45 = 0.0115661327577615d0*gamma_cc*rho_c**(-2.33333333333333d0)*x
     &44/x36**2
      x46 = 1.0d0 - x45/x39
      x47 = x28*x37
      x48 = 0.125d0*x47
      x49 = x29*(x16 - x18)
      x50 = x34 - x49
      x51 = rho_c + rho_s
      x52 = 1d0/x51
      x53 = x0*x52**0.333333333333333d0
      x54 = 0.819204448708788d0*x4
      x55 = 4.12273618607837d0*x7
      x56 = 15.10599739383d0*x9
      x57 = (0.235215882559859d0*x53 + 1.0d0)*log(x15/(x52**0.1666666666
     &6666666d0* x56 + sqrt(x52)*x55 + x52**0.66666666666666663d0*x54 + 
     &7.0945954610727853d0*x53) + 1.0d0)
      x58 = exp(2.0000000000000004d0*x57) - 1.0d0
      x59 = 2.0d0*gamma_cs
      x60 = gamma_cc + gamma_ss
      x61 = x59 + x60
      x62 = 3.34438567425411d-5*x42
      x63 = 0.00578306637888076d0*x44
      x64 = x51**(-2.33333333333333d0)*x61*x63
      x65 = 1.0d0 + x64/x58
      x66 = -1.0d0*x57 + 0.5d0*log(x64*x65/(x51**(-4.666666666666667d0)*
     &x61**2*x62/ x58**2 + x65) + 1.0d0)
      x67 = rho_c - rho_s
      x68 = 1d0/x67
      x69 = x0*x68**0.333333333333333d0
      x70 = (0.235215882559859d0*x69 + 1.0d0)*log(x15/(x54*x68** 0.66666
     &666666666663d0 + x55*sqrt(x68) + x56*x68** 0.16666666666666666d0 +
     & 7.0945954610727853d0*x69) + 1.0d0)
      x71 = exp(2.0000000000000004d0*x70) - 1.0d0
      x72 = -x59 + x60
      x73 = x63*x67**(-2.33333333333333d0)*x72
      x74 = 1.0d0 + x73/x71
      x75 = -1.0d0*x70 + 0.5d0*log(x73*x74/(x62*x67**(-4.666666666666667
     &d0)*x72**2/ x71**2 + x74) + 1.0d0)
      x76 = x51*x66
      x77 = x67*x75
      x78 = rho_c**(-3)
      x79 = 0.5d0*x26
      x80 = rho_c**2
      x81 = 1d0/x80
      x82 = rho_s**2
      x83 = x81*x82
      x84 = 4.60115111447049d0
      x85 = 1d0/x84
      x86 = (x79 + 0.87d0*x83 + 0.53d0 + 2.26d0*rho_s**6/rho_c**6)/( 0.0
     &60093732096142d0*x40*x85*(1d0/x24 + 1d0/x22)*(gamma_cc*x82 + gamma
     &_ss*x80 - rho_c*rho_s*x59) + 1.0d0)**4
      x87 = x86 + 1.0d0
      x88 = x41/tau_c**2
      x89 = x78*x87*x88
      x90 = x27*x89*(-0.5d0*x12 + 0.5d0)*(x76 + x77)
      x91 = x81*x88
      x92 = x31 - x33
      x93 = exp(x38*(-x30 + x92)) - 1.0d0
      x94 = x45/x93 + 1.0d0
      x95 = log(x45*x94/(x43/x93**2 + x94) + 1.0d0)
      x96 = x35 + x48*x95
      x97 = x96*(x86*x91 + 64.0d0)
      x98 = 2.8d0*x97
      x99 = gamma_cc**3*x78/tau_c**3
      x100 = x17*x28
      x101 = x47*x95
      x102 = -4194304.0d0*x100 + 262144.0d0*x101 + 2097152.0d0*x30 + 414
     &22.505278314d0*x32
      x103 = 4.76837158203125d-7*exp(-1.0053738478084204d0*mu*rho_c**( -
     &0.33333333333333331d0)*x4)
      x104 = -x83 + 1.0d0
      x105 = x1**1.33333333333333d0*x85
      x106 = 3.14159265358979d0
      x107 = x1**1.0d0
      x108 = x107/x106
      x109 = (0.00126674656487366d0*x105 - 0.009578475d0*x108 + 0.018817
     &1922990732d0* x2 + 0.0676322201645715d0*x6 + 1.0d0)*exp(-0.6836107
     &611867116d0* x2)
      x110 = x104*x109
      x111 = 1.41421356237310d0
      x112 = x1**1.66666666666667d0
      x113 = 0.0837224526465878d0*x112
      x114 = 0.0892278228691848d0*mu**5*x111*x113
      x115 = x21**0.666666666666667d0 + x23**0.666666666666667d0
      x116 = x115**3
      x117 = mu**2
      x118 = x117*x2/x115**2
      x119 = 11.1447260721495d0*mu*x10/x115 + 1.0d0
      x120 = mu**3
      x121 = x116*log((27.073371959220147d0*x118 + x119 + 27.14082046241
     &0485d0*x120* x8/x116)/(12.532717071175123d0*x118 + x119))
      x122 = 0.0315054072231411d0*x107
      x123 = 1d0/x5
      x124 = 1d0/x21
      x125 = x124**0.666666666666667d0
      x126 = 0.0524148278841779d0*x6
      x127 = x124**0.333333333333333d0*x2
      x128 = 2.28942848510666d0*x84
      x129 = 0.817704266774463d0*x128
      x130 = x129*x21**2*(-0.0259335011650457d0*x127 + 1.0d0)/(x125*(x12
     &5*x126 + 0.494402081358784d0*x127 + 1.0d0))
      x131 = 1d0/x23
      x132 = x131**0.666666666666667d0
      x133 = x131**0.333333333333333d0*x2
      x134 = x129*x23**2*(-0.0259335011650457d0*x133 + 1.0d0)/(x132*(x12
     &6*x132 + 0.494402081358784d0*x133 + 1.0d0))
      x135 = 0.00111534778586481d0*x113*(12.0d0*x104*x106*(-4.4973734672
     &5955d0*x2 + 0.825481812223657d0*x6)*exp(-0.28165369188898165d0*x2)
     &/x107 + x123*x130 + x123*x134)
      x136 = 4.0d0*x111*x120
      x137 = 0.0126441184989504d0*mu**8*x1**2.66666666666667d0
      x138 = 0.148394183600699d0
      x139 = x104*(x109 - 1.0d0)
      x140 = 0.0167302167879722d0*x112*x138*x139
      x141 = 0.0533250677421794d0*x1**2.0d0*x35
      x142 = 4.0d0*mu**6
      x143 = -0.13157433081915d0*x108*x139
      x144 = 0.817704266774463d0*x128*(x21**2.66666666666667d0 + x23** 2
     &.66666666666667d0)
      x145 = 10.9027235569928d0*x104*x3*(-0.352521395009435d0*x2 + 0.558
     &025705063192d0*x6)*exp(-0.49698248213959023d0*x2)
      x146 = 0.00139418473233101d0*x1**1.0d0*x138
      x147 = 1.55214407110551d0*x105
      x148 = 4.0d0*mu**4
      x149 = 0.25d0/(0.508616435555896d0*x117*x6 + 1.0d0)**4
      x150 = x149*(-x110*x114 - x121*x28 - x136*(x110*x122 + x135) + x13
     &7*x50 - x142* (x140 - x141) + x148*(x143 - x146*(x130 + x134 - x14
     &4 + x145) + x147*x50))
      x151 = (-x150 + x50)/x50
      x152 = x109*(x83 - 1.0d0)
      x153 = -1.0d0 + exp(-x38*(-x149*(x114*x152 + x121*x27*(x12 - 1.0d0
     &) + x136*( x122*x152 - x135) + x137*x35 + x142*(-x140 + x141) + x1
     &48*(x143 + x146*(-x130 - x134 + x144 - x145) + x147*x35)) + x35))
      x154 = x151**2.78d0*x45
      x155 = 1.0d0 + x154/x153
      x156 = x150 - 0.125d0*x47*log(x154*x155/(x151**5.5599999999999996d
     &0*x43/x153**2 + x155) + 1.0d0) + x49 + x92
      x157 = x89*(x28*x76 + x67*x96)
      x158 = x89*(x28*x77 + x51*x96)
      x159 = x87*x91*(-1.0d0*x100 + 0.0625d0*x101 + x19*x25*x28*x79 + 0.
     &0098758948512826d0*x32)
      PBE_value = x48*log(x45*x46/(x46 + x43/x39**2) + 1.0d0) + x50
      PBE_beta_value = x28*x66
      PBE_alpha_value = x28*x75
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x103*(x102 - (-x90 + x97)*(x99*(-2.8d0*x90 + x98) + 32
     &768.0d0)) + x156)
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x103*(x102 - (-0.5d0*x157 + x97)*(x99*(-1.4d0*x157 + x
     &98) + 32768.0d0)) + x156)
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x103*(x102 - (-0.5d0*x158 + x97)*(x99*(-1.4d0*x158 + x
     &98) + 32768.0d0)) + x156)
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x103*(x102 - (-2.0d0*x159 + x97)*(x99*(-5.6d0*x159 + x
     &98) + 32768.0d0)) + x156)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, g
     &amma_ss, gamma_cs, tau_c, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 18, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, gamma_ss, gamma_cs, 
     &tau_c, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
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
     &, x985, x986, x987, x988, x989, x990
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = rho_s*x0
      x2 = x1 + 1.0d0
      x3 = x2**1.33333333333333d0
      x4 = -x1 + 1.0d0
      x5 = x4**1.33333333333333d0
      x6 = 1.92366105093154d0*x3 + 1.92366105093154d0*x5 - 3.84732210186
     &307d0
      x7 = rho_c**(-4)
      x8 = rho_s**4
      x9 = x7*x8
      x10 = x6*x9
      x11 = 1.46459188756152d0
      x12 = 1d0/x11
      x13 = x0**0.333333333333333d0
      x14 = x12*x13
      x15 = 0.186690969707574d0*x14 + 1.0d0
      x16 = 2.14502939711103d0
      x17 = 1d0/x16
      x18 = x0**0.666666666666667d0
      x19 = x17*x18
      x20 = 0.564189583547756d0
      x21 = sqrt(x0)*x20
      x22 = 0.826307487110758d0
      x23 = x0**0.166666666666667d0*x22
      x24 = 5.63098414909788d0*x14 + 0.516066464547863d0*x19 + 2.9152147
     &1421918d0* x21 + 13.4579137143944d0*x23
      x25 = 9.86960440108936d0
      x26 = 0.693147180559945d0
      x27 = -x26 + 1.0d0
      x28 = 1d0/x27
      x29 = x25*x28
      x30 = 1.0d0 + x29/x24
      x31 = log(x30)
      x32 = x15*x31
      x33 = 0.194159335344114d0*x14 + 1.0d0
      x34 = 3.25955091942229d0*x14 + 0.406913004517529d0*x19 + 1.4187228
     &1647967d0* x21 + 7.24010193431683d0*x23
      x35 = x29/x34
      x36 = log(0.5d0*x35 + 1.0d0)
      x37 = x33*x36
      x38 = 2.0d0*x37
      x39 = 1d0/x25
      x40 = x27*x39
      x41 = x40*(-x32 + x38)
      x42 = x10*x41
      x43 = x38*x40
      x44 = -x9 + 1.0d0
      x45 = 0.101077332976288d0*x14 + 1.0d0
      x46 = 3.29180480994506d0*x14 + 0.410025070949612d0*x19 + 0.7623275
     &2193529d0* x21 + 9.87212972256927d0*x23
      x47 = 1.0d0 + 29.6085746432167d0/x46
      x48 = log(x47)
      x49 = x45*x48
      x50 = x44*x49
      x51 = x50*x6
      x52 = 0.0197517897025652d0*x51
      x53 = -x43 + x52
      x54 = x42 + x53
      x55 = x2**0.666666666666667d0 + x4**0.666666666666667d0
      x56 = x55**3
      x57 = 1d0/x56
      x58 = 8.0d0*x29*x57
      x59 = exp(-x54*x58)
      x60 = -x59 + 1.0d0
      x61 = gamma_cc**2
      x62 = x27**(-2)
      x63 = x62/x55**4
      x64 = rho_c**(-4.66666666666667d0)
      x65 = 208.945363813799d0
      x66 = x64*x65
      x67 = 0.000133775426970164d0*x66
      x68 = x61*x63*x67
      x69 = 14.4549425392770d0
      x70 = rho_c**(-2.33333333333333d0)
      x71 = 0.0115661327577615d0*x70
      x72 = x69*x71
      x73 = x55**2
      x74 = x28/x73
      x75 = gamma_cc*x74
      x76 = x72*x75
      x77 = 1.0d0 - x76/x60
      x78 = x40*x56
      x79 = 0.125d0*x78
      x80 = x40*(x32 - x38)
      x81 = x10*x80
      x82 = x53 - x81
      x83 = rho_c + rho_s
      x84 = 1d0/x83
      x85 = x84**0.333333333333333d0
      x86 = x12*x85
      x87 = 0.235215882559859d0*x86 + 1.0d0
      x88 = x84**0.666666666666667d0
      x89 = 0.819204448708788d0*x17
      x90 = sqrt(x84)
      x91 = 4.12273618607837d0*x20
      x92 = x84**0.166666666666667d0
      x93 = 15.10599739383d0*x22
      x94 = 7.09459546107279d0*x86 + x88*x89 + x90*x91 + x92*x93
      x95 = x29/x94 + 1.0d0
      x96 = log(x95)
      x97 = x87*x96
      x98 = 1.0d0*x97
      x99 = 2.0d0*gamma_cs
      x100 = gamma_cc + gamma_ss
      x101 = x100 + x99
      x102 = x101*x28
      x103 = exp(2.0000000000000004d0*x97)
      x104 = x103 - 1.0d0
      x105 = x104**(-2)
      x106 = x65*x83**(-4.66666666666667d0)
      x107 = 3.34438567425411d-5*x62
      x108 = x83**(-2.33333333333333d0)
      x109 = 0.00578306637888076d0*x69
      x110 = x108*x109
      x111 = 1d0/x104
      x112 = x102*x111
      x113 = x110*x112 + 1.0d0
      x114 = 1d0/(x101**2*x105*x106*x107 + x113)
      x115 = x110*x113
      x116 = x114*x115
      x117 = x102*x116 + 1.0d0
      x118 = log(x117)
      x119 = 0.5d0*x118
      x120 = x119 - x98
      x121 = rho_c - rho_s
      x122 = 1d0/x121
      x123 = x122**0.333333333333333d0
      x124 = x12*x123
      x125 = 0.235215882559859d0*x124 + 1.0d0
      x126 = x122**0.666666666666667d0
      x127 = sqrt(x122)
      x128 = x122**0.166666666666667d0
      x129 = 7.09459546107279d0*x124 + x126*x89 + x127*x91 + x128*x93
      x130 = 1.0d0 + x29/x129
      x131 = log(x130)
      x132 = x125*x131
      x133 = 1.0d0*x132
      x134 = x100 - x99
      x135 = x134*x28
      x136 = exp(2.0000000000000004d0*x132)
      x137 = x136 - 1.0d0
      x138 = x137**(-2)
      x139 = x121**(-4.66666666666667d0)*x65
      x140 = x121**(-2.33333333333333d0)
      x141 = x109*x140
      x142 = 1d0/x137
      x143 = x142*x28
      x144 = x134*x143
      x145 = x141*x144 + 1.0d0
      x146 = 1d0/(x107*x134**2*x138*x139 + x145)
      x147 = x141*x145
      x148 = x146*x147
      x149 = x135*x148 + 1.0d0
      x150 = log(x149)
      x151 = 0.5d0*x150
      x152 = -x133 + x151
      x153 = tau_c**(-2)
      x154 = x153*x61
      x155 = rho_c**(-3)
      x156 = rho_c**(-6)
      x157 = rho_s**6*x156
      x158 = 0.5d0*x9
      x159 = rho_c**2
      x160 = 1d0/x159
      x161 = rho_s**2
      x162 = x160*x161
      x163 = 2.26d0*x157 + x158 + 0.87d0*x162 + 0.53d0
      x164 = gamma_cc*x161
      x165 = gamma_ss*x159 - rho_c*rho_s*x99 + x164
      x166 = 4.60115111447049d0
      x167 = 1d0/x166
      x168 = 1d0/x5 + 1d0/x3
      x169 = x168*x64
      x170 = x167*x169
      x171 = 0.060093732096142d0*x165*x170 + 1.0d0
      x172 = x171**(-4)
      x173 = x163*x172
      x174 = x173 + 1.0d0
      x175 = x155*x174
      x176 = x154*x175
      x177 = x120*x83
      x178 = x121*x152
      x179 = -0.5d0*x26 + 0.5d0
      x180 = x179*x39
      x181 = x180*(x177 + x178)
      x182 = x176*x181
      x183 = x154*x160
      x184 = x173*x183
      x185 = x184 + 64.0d0
      x186 = x43 - x52
      x187 = x186 - x42
      x188 = exp(x187*x58)
      x189 = x188 - 1.0d0
      x190 = 1d0/x189
      x191 = x190*x76 + 1.0d0
      x192 = x189**(-2)
      x193 = 1d0/(x191 + x192*x68)
      x194 = x191*x193
      x195 = x194*x76 + 1.0d0
      x196 = log(x195)
      x197 = x196*x79
      x198 = x197 + x54
      x199 = x185*x198
      x200 = -x182 + x199
      x201 = 2.8d0*x182
      x202 = 2.8d0*x199
      x203 = -x201 + x202
      x204 = tau_c**(-3)
      x205 = x155*x204
      x206 = gamma_cc**3*x205
      x207 = x203*x206 + 32768.0d0
      x208 = x37*x40
      x209 = x196*x78
      x210 = -4194304.0d0*x208 + 262144.0d0*x209 + 2097152.0d0*x42 + 414
     &22.505278314d0*x51
      x211 = -x200*x207 + x210
      x212 = mu*x17
      x213 = exp(-1.0053738478084204d0*rho_c**(-0.33333333333333331d0)*x
     &212)
      x214 = 4.76837158203125d-7*x213
      x215 = x211*x214
      x216 = 1d0/x82
      x217 = mu**2
      x218 = x19*x217
      x219 = 0.508616435555896d0*x218 + 1.0d0
      x220 = x219**(-4)
      x221 = -x162 + 1.0d0
      x222 = exp(-0.6836107611867116d0*x14)
      x223 = x0**1.33333333333333d0
      x224 = x167*x223
      x225 = 3.14159265358979d0
      x226 = x0**1.0d0
      x227 = x226/x225
      x228 = 0.009578475d0*x227
      x229 = -x228
      x230 = 0.0188171922990732d0*x14 + 0.0676322201645715d0*x19 + 0.001
     &26674656487366d0*x224 + x229 + 1.0d0
      x231 = x222*x230
      x232 = x221*x231
      x233 = 1.41421356237310d0
      x234 = mu**5*x233
      x235 = x232*x234
      x236 = x0**1.66666666666667d0
      x237 = 0.0837224526465878d0*x236
      x238 = 0.0892278228691848d0*x237
      x239 = x2**0.666666666666667d0 + x4**0.666666666666667d0
      x240 = x239**3
      x241 = x239**2
      x242 = 1d0/x241
      x243 = x217*x242
      x244 = x14*x243
      x245 = 1d0/x239
      x246 = mu*x245
      x247 = 11.1447260721495d0*x23*x246 + 1.0d0
      x248 = mu**3
      x249 = x21/x240
      x250 = 27.0733719592201d0*x244 + x247 + 27.1408204624105d0*x248*x2
     &49
      x251 = x250/(12.5327170711751d0*x244 + x247)
      x252 = log(x251)
      x253 = x240*x252
      x254 = 0.179587122125167d0
      x255 = 0.175432441092201d0*x226*x254
      x256 = x232*x255
      x257 = -4.49737346725955d0*x14 + 0.825481812223657d0*x19
      x258 = exp(-0.28165369188898165d0*x14)
      x259 = x221*x257*x258
      x260 = x225/x226
      x261 = 12.0d0*x260
      x262 = x259*x261
      x263 = x2**2
      x264 = 1d0/x2
      x265 = x264**0.666666666666667d0
      x266 = 0.0524148278841779d0*x19
      x267 = x264**0.333333333333333d0
      x268 = x14*x267
      x269 = x265*x266 + 0.494402081358784d0*x268 + 1.0d0
      x270 = 1d0/x269
      x271 = x263*x270
      x272 = (-0.0259335011650457d0*x268 + 1.0d0)/x265
      x273 = x271*x272
      x274 = 1d0/x18
      x275 = 2.28942848510666d0
      x276 = x166*x275
      x277 = x274*x276
      x278 = 0.817704266774463d0*x277
      x279 = x4**2
      x280 = 1d0/x4
      x281 = x280**0.666666666666667d0
      x282 = x280**0.333333333333333d0
      x283 = x14*x282
      x284 = x266*x281 + 0.494402081358784d0*x283 + 1.0d0
      x285 = 1d0/x284
      x286 = x279*x285
      x287 = (-0.0259335011650457d0*x283 + 1.0d0)/x281
      x288 = x286*x287
      x289 = x237*(x262 + x273*x278 + x278*x288)
      x290 = 0.00111534778586481d0*x289
      x291 = 4.0d0*x233*x248
      x292 = 0.0472353356922751d0*mu**8*x0**2.66666666666667d0
      x293 = 0.267683468607554d0*x292
      x294 = 0.148394183600699d0
      x295 = x236*x294
      x296 = x231 - 1.0d0
      x297 = x221*x296
      x298 = x295*x297
      x299 = 0.0167302167879722d0*x298
      x300 = 0.101321183642338d0
      x301 = x0**2.0d0
      x302 = x300*x301*x54
      x303 = 0.526297323276602d0*x302
      x304 = 4.0d0*mu**6
      x305 = 0.13157433081915d0*x227*x297
      x306 = -x305
      x307 = x276*(x2**2.66666666666667d0 + x4**2.66666666666667d0)
      x308 = 0.817704266774463d0*x307
      x309 = -0.352521395009435d0*x14 + 0.558025705063192d0*x19
      x310 = exp(-0.49698248213959023d0*x14)
      x311 = x221*x309*x310
      x312 = x16*x311
      x313 = 10.9027235569928d0*x312
      x314 = 0.817704266774463d0*x276
      x315 = x273*x314
      x316 = x288*x314
      x317 = x0**1.0d0*x294
      x318 = 0.00139418473233101d0*x317
      x319 = 1.55214407110551d0*x224
      x320 = 4.0d0*mu**4
      x321 = -x235*x238 - x253*x40 - x291*(x256 + x290) + x293*x82 - x30
     &4*(x299 - x303) + x320*(x306 - x318*(-x308 + x313 + x315 + x316) +
     & x319*x82 )
      x322 = x220*x321
      x323 = 0.25d0*x322
      x324 = -x323 + x82
      x325 = x216*x324
      x326 = x325**2.78d0
      x327 = x326*x72
      x328 = x162 - 1.0d0
      x329 = x231*x328
      x330 = x234*x238
      x331 = x39*(x26 - 1.0d0)
      x332 = x255*x329
      x333 = x308 - x313 - x315 - x316
      x334 = x253*x331 + x291*(-x290 + x332) + x293*x54 + x304*(-x299 + 
     &x303) + x320* (x306 + x318*x333 + x319*x54) + x329*x330
      x335 = x220*x334
      x336 = 0.25d0*x335
      x337 = -x336 + x54
      x338 = exp(-x337*x58)
      x339 = x338 - 1.0d0
      x340 = 1d0/x339
      x341 = x340*x75
      x342 = x327*x341 + 1.0d0
      x343 = x325**5.56d0
      x344 = x343/x339**2
      x345 = 1d0/(x342 + x344*x68)
      x346 = x342*x345
      x347 = x327*x346*x75 + 1.0d0
      x348 = log(x347)
      x349 = 0.125d0*x348*x78
      x350 = x186 + x81
      x351 = x323 + x350
      x352 = -x349 + x351
      x353 = x27/x69
      x354 = 0.0784052941866196d0*x353
      x355 = 0.546136299139192d0*x17
      x356 = 2.06136809303918d0*x20
      x357 = 2.517666232305d0*x22
      x358 = (x355*x88 + x356*x90 + x357*x92 + 2.36486515369093d0*x86)/(
     &x94**2*x95)
      x359 = x83**(-3.33333333333333d0)
      x360 = 0.0115661327577615d0*x359
      x361 = x360*x69
      x362 = 0.00192768879296025d0*x83**(-3.33333333333333d0)
      x363 = x362*x69
      x364 = 0.00578306637888076d0*x111
      x365 = x103*(x29*x358*(0.470431765119717d0*x86 + 2.0d0) - 0.156810
     &588373239d0* x86*x96)
      x366 = x359*x365
      x367 = x102*x65
      x368 = x83**(-5.66666666666667d0)
      x369 = x112*x65
      x370 = x112*x116
      x371 = x114/x117
      x372 = 0.5d0*x39
      x373 = 0.5d0*x83
      x374 = 0.25d0*x40
      x375 = 0.5d0*x40
      x376 = x118*x374 + x373*(x101*x371*x372*(-x108*x364*x367*(x360 + x
     &362 + x364* x366) - x113*x361 - x113*x363 + x370*(6.68877134850821
     &d-5*x105* x365*x367*x368 + x109*x111*x366 + x361 + x363 + 0.000133
     &775426970164d0*x368*x369 + 2.22959044950274d-5*x369*x83** (-5.6666
     &6666666667d0))) + x354*x84*x85*x96 - 1.0d0*x358*x84*x87) - x375*x9
     &7
      x377 = x122*x123*x131*x354
      x378 = (2.36486515369093d0*x124 + x126*x355 + x127*x356 + x128*x35
     &7)/(x129**2* x130)
      x379 = 1.0d0*x122*x125*x378
      x380 = x121**(-3.33333333333333d0)
      x381 = 0.0115661327577615d0*x380
      x382 = x381*x69
      x383 = x145*x382
      x384 = 0.00192768879296025d0*x121**(-3.33333333333333d0)
      x385 = x384*x69
      x386 = x145*x385
      x387 = 0.156810588373239d0*x124*x131
      x388 = x29*x378*(0.470431765119717d0*x124 + 2.0d0)
      x389 = x136*(-x387 + x388)
      x390 = 0.00578306637888076d0*x142
      x391 = x380*x390
      x392 = x381 + x384
      x393 = x135*x65
      x394 = x140*x390*x393
      x395 = x121**(-5.66666666666667d0)
      x396 = x144*x65
      x397 = 0.000133775426970164d0*x395*x396
      x398 = 2.22959044950274d-5*x121**(-5.66666666666667d0)*x396
      x399 = x109*x142*x380
      x400 = 6.68877134850821d-5*x138*x393*x395
      x401 = x144*x148
      x402 = x146/x149
      x403 = x134*x372*x402
      x404 = 0.5d0*x121
      x405 = x150*x374
      x406 = x132*x375
      x407 = x404*(x377 - x379 + x403*(-x383 - x386 - x394*(x389*x391 + 
     &x392) + x401* (x382 + x385 + x389*x399 + x389*x400 + x397 + x398))
     &) + x405 - x406
      x408 = x176*(x376 + x407)
      x409 = x154*x174
      x410 = x409*x7
      x411 = x181*x410
      x412 = 2.0d0*x162 + 13.56d0*x9 + 1.74d0
      x413 = x155*x161*x412
      x414 = 0.480749856769136d0*x169
      x415 = rho_c**(-5.66666666666667d0)
      x416 = x165*x415
      x417 = x2**(-2.33333333333333d0)
      x418 = x4**(-2.33333333333333d0)
      x419 = rho_c**(-6.66666666666667d0)
      x420 = rho_s*x419
      x421 = x167/x171
      x422 = x163*x421
      x423 = x422*(-0.320499904512757d0*x165*x420*(x417 - x418) + 1.1217
     &4966579465d0* x168*x416 - x414*(-gamma_cs*rho_s + gamma_ss*rho_c))
     &
      x424 = -x413 + x423
      x425 = x154*x155
      x426 = x172*x425
      x427 = x424*x426
      x428 = x181*x427
      x429 = x13*x353
      x430 = x36*x429
      x431 = 0.129439556896076d0*x430
      x432 = 0.0790071588102608d0*x49
      x433 = rho_c**(-5)
      x434 = x433*x8
      x435 = x434*x6
      x436 = x14*x48
      x437 = 0.000665486074881265d0*x436
      x438 = x44*x6
      x439 = x437*x438
      x440 = rho_s*x160
      x441 = x4**0.333333333333333d0
      x442 = x2**0.333333333333333d0
      x443 = x441 - x442
      x444 = x443*x50
      x445 = 0.0506609980493537d0*x444
      x446 = 2.0d0*x0
      x447 = x33*(1.08651697314076d0*x14 + 0.27127533634502d0*x19 + 0.70
     &9361408239834d0*x21 + 1.20668365571947d0*x23)/(x34**2*(x35 + 2.0d0
     &))
      x448 = 2.56488140124205d0*x41
      x449 = rho_s**5*x443
      x450 = x156*x449
      x451 = x45*(1.09726826998169d0*x14 + 0.273350047299742d0*x19 + 0.3
     &81163760967645d0*x21 + 1.64535495376155d0*x23)/(x46**2*x47)
      x452 = 0.58482233974552d0*x451
      x453 = x438*x452
      x454 = x41*x435
      x455 = x15*(1.87699471636596d0*x14 + 0.344044309698576d0*x19 + 1.4
     &5760735710959d0*x21 + 2.24298561906574d0*x23)/(x24**2*x30)
      x456 = 0.062230323235858d0*x31*x429
      x457 = 2.0d0*x447
      x458 = -x431 + x457
      x459 = x6*(-x455 + x456 + x458)
      x460 = x434*x459
      x461 = 0.25d0*x196
      x462 = x40*x73
      x463 = x4**(-0.333333333333333d0)
      x464 = x2**(-0.333333333333333d0)
      x465 = x463 - 1.0d0*x464
      x466 = x440*x465
      x467 = x462*x466
      x468 = rho_c**(-3.33333333333333d0)
      x469 = 0.0269876431014436d0*x468
      x470 = 1d0/x55
      x471 = rho_c**(-4.33333333333333d0)*rho_s
      x472 = x470*x471
      x473 = -x463 + x464
      x474 = 0.015421510343682d0*x473
      x475 = x469*x69
      x476 = x470*x69
      x477 = x471*x474*x476
      x478 = 13.2006183429726d0
      x479 = x28*x57
      x480 = x478*x479
      x481 = x468*x480
      x482 = 1.0d0*x463
      x483 = x464 - x482
      x484 = 2.0d0*x470
      x485 = x1*x484
      x486 = x483*x485
      x487 = x9 - 1.0d0
      x488 = x487*x6
      x489 = x487*x49
      x490 = x443*x489
      x491 = x1*x490
      x492 = 4.0d0*x42
      x493 = x459*x9
      x494 = x10*x432
      x495 = x433*x449
      x496 = x448*x495
      x497 = x458 + x492 - x493 - x494 - x496
      x498 = x188*(x187*x486 - x437*x488 + x452*x488 + 0.050660998049353
     &7d0*x491 + x497)
      x499 = x190*x481*x498
      x500 = x190*x75
      x501 = x500*x71
      x502 = x415*x500
      x503 = x502*x65
      x504 = gamma_cc*rho_s
      x505 = x419*x479*x504*x65
      x506 = 2062.20808228389d0
      x507 = 0.00214040683152263d0*gamma_cc*x192*x415*x506
      x508 = x62/x55**5
      x509 = x194*x501
      x510 = 1.46459188756152d0*gamma_cc
      x511 = x193*x510/x195
      x512 = x511*x55*(-x191*x469 + x191*x472*x474 + x501*(-x475 + x477 
     &- x499) + x509*(-0.000356734471920438d0*x190*x473*x505 + x475 - x4
     &77 + x498 *x507*x508 + x499 + 0.000624285325860766d0*x503))
      x513 = x0*x431 - x0*x439 + x0*x453 + x432*x435 + x440*x445 - x446*
     &x447 + x448* x450 - 4.0d0*x454 + x460 + x461*x467 + 0.125d0*x512
      x514 = x172*x183
      x515 = x198*x514
      x516 = x185*x513 + x515*(-x163*x446 + x424)
      x517 = 3.0d0*x0
      x518 = 0.04375d0*x184 + 2.8d0
      x519 = 64.0d0*x513*x518 + 64.0d0*x515*(-x0*(0.19775d0*x157 + 0.076
     &125d0*x162 + 0.04375d0*x9 + 0.046375d0) - 0.04375d0*x413 + 0.04375
     &d0*x423)
      x520 = x0*x488
      x521 = 1395.62545290939d0*x436
      x522 = 5378946.16037756d0*x450
      x523 = 1226461.339442d0*x451
      x524 = 165690.021113256d0*x49
      x525 = 271454.42562372d0*x0*x430 - 4194304.0d0*x0*x447 + x435*x524
     &
      x526 = 524288.0d0*x196*x467 + x41*x522 - 106243.813381198d0*x440*x
     &490 - 8388608.0d0*x454 + 2097152.0d0*x460 + 262144.0d0*x512 + x520
     &*x521 - x520*x523 + x525
      x527 = 0.33512461593614d0*rho_c**(-1.33333333333333d0)*x212*x213
      x528 = x162*x231
      x529 = 0.17845564573837d0*x234*x237
      x530 = -x528*x529
      x531 = 0.00168899541983155d0*x224
      x532 = 0.045088146776381d0*x19
      x533 = 0.00627239743302441d0*x14
      x534 = x221*x222
      x535 = x534*(x229 + x531 + x532 + x533)
      x536 = 0.00116228665296198d0*x0**2.0d0
      x537 = 0.148713038115308d0*x237
      x538 = x241*x40
      x539 = x2**(-0.333333333333333d0)
      x540 = x4**(-0.333333333333333d0)
      x541 = x539 - x540
      x542 = 2.0d0*x252
      x543 = x541*x542
      x544 = x538*x543
      x545 = 0.713822582953479d0*x292
      x546 = x14*x246
      x547 = 4.17757235705837d0*x546
      x548 = x1*x541
      x549 = mu*x14*x242*x548
      x550 = 16.7102894282335d0*x549
      x551 = 1.85745434535825d0*x23
      x552 = 7.42981738143299d0*x23
      x553 = x245*x548*x552
      x554 = -x551 + x553
      x555 = x21*x243
      x556 = 54.281640924821d0*x217*x249*x548 - 9.02445731974005d0*x546 
     &+ 36.0978292789602d0*x549 + x554 - 13.5704102312052d0*x555
      x557 = mu/x250
      x558 = x431 - x457
      x559 = x6*(x455 - x456 + x558)
      x560 = x559*x9
      x561 = 4.0d0*x81
      x562 = x495*x80
      x563 = x1*x445
      x564 = -x439 + x453 + x494 + x558 + x563
      x565 = 24.0d0*x257*x258*x260
      x566 = x162*x565
      x567 = 2.41662179562687d0*x0**(-0.666666666666667d0)*x259
      x568 = x221*x258*x261*(1.49912448908652d0*x14 - 0.550321208149104d
     &0*x19)
      x569 = x1*x264
      x570 = 1d0/x13
      x571 = 0.00706864485168613d0*x225*x275*x570
      x572 = x571/x267
      x573 = x271*x572*(-x569 + 1.0d0)
      x574 = x1*x280
      x575 = x571/x282
      x576 = x286*x575*(x574 + 1.0d0)
      x577 = 0.545136177849642d0*x277
      x578 = x273*x577
      x579 = x2*x270
      x580 = 2.18054471139857d0*x277
      x581 = x272*x579*x580
      x582 = x1*x581
      x583 = x288*x577
      x584 = x285*x4
      x585 = x287*x580*x584
      x586 = x1*x585
      x587 = 0.034943218589452d0*x19
      x588 = x265*x587
      x589 = 0.164800693786261d0*x268
      x590 = x588 + x589
      x591 = x272*x278/x269**2
      x592 = x263*x591*(-x569*x588 - x569*x589 + x590)
      x593 = x281*x587
      x594 = 0.164800693786261d0*x283
      x595 = x593 + x594
      x596 = x278*x287/x284**2
      x597 = x279*x596*(x574*x593 + x574*x594 + x595)
      x598 = 0.00111534778586481d0*x237
      x599 = 0.00490180588786583d0*x223
      x600 = 0.350864882184401d0*x226*x254
      x601 = 0.00185891297644135d0*x289 - x528*x600
      x602 = x162*x296
      x603 = 0.0334604335759443d0*x295
      x604 = x602*x603
      x605 = 0.0278836946466203d0*x298
      x606 = x534*(0.113935126864452d0*x14*x230 - 0.0031361987165122d0*x
     &14 - 0.0225440733881905d0*x19 - 0.000844497709915773d0*x224 + 0.00
     &47892375d0*x227)
      x607 = x603*x606
      x608 = 1.0525946465532d0*x302
      x609 = -x492
      x610 = x493 + x496 + x564 + x609
      x611 = 0.526297323276602d0*x300*x301
      x612 = x610*x611
      x613 = 0.263148661638301d0*x227
      x614 = x602*x613
      x615 = x606*x613
      x616 = 0.00232364122055169d0*x317*x333
      x617 = 2.06952542814068d0*x224*x54
      x618 = x16*x274*x310
      x619 = 21.8054471139857d0*x309*x618
      x620 = x4**1.66666666666667d0
      x621 = x2**1.66666666666667d0
      x622 = 2.18054471139857d0*x277
      x623 = -x573 - x576 - x578 + x582 - x583 - x586 - x592 - x597
      x624 = 0.00139418473233101d0*x295
      x625 = x624*(x1*x622*(x620 - x621) - 1.80615420514536d0*x11*x311*x
     &570 - x162* x619 - 10.9027235569928d0*x221*x618*(0.117507131669812
     &d0*x14 - 0.372017136708795d0*x19) + 0.545136177849642d0*x274*x307 
     &- 7.26848237132856d0*x274*x312 + x623)
      x626 = x319*x610
      x627 = x220*(x1*x544 - x235*x536 + x235*x537 + x291*(-x232*x599 + 
     &x255*x535 + x256 - x598*(x262 + x566 + x567 + x568 + x573 + x576 +
     & x578 - x582 + x583 + x586 + x592 + x597) + x601) + x293*(-x560 + 
     &x561 - 2.56488140124205d0*x562 + x564) - x304*(x604 - x605 + x607 
     &+ x608 - x612) - x320*(x306 + x614 + x615 + x616 + x617 - x625 - x
     &626) + x330*x535 + x530 - x538*x557*(-x251*(-x547 + x550 + x554) +
     & x556 ) - x545*x82)
      x628 = 8388608.0d0*x80
      x629 = 106243.813381198d0*x444
      x630 = x0*x438
      x631 = x219**(-5)
      x632 = -x464
      x633 = x463 + x632
      x634 = 0.015421510343682d0*x633
      x635 = 8.0d0*x208
      x636 = 0.0790071588102608d0*x51
      x637 = -x635 + x636
      x638 = -x322 - x561 + x637
      x639 = 1d0/x638
      x640 = x216*x638
      x641 = 4.0d0*x560
      x642 = 1.35631049481572d0*x218*x631
      x643 = x321*x642
      x644 = 16.0d0*x81
      x645 = 10.2595256049682d0*x562
      x646 = 8.0d0*x447
      x647 = 0.517758227584304d0*x430
      x648 = 0.316028635241043d0*x49
      x649 = x10*x648
      x650 = 0.00266194429952506d0*x436
      x651 = x438*x650
      x652 = 2.33928935898208d0*x451
      x653 = x438*x652
      x654 = 0.202643992197415d0*x444
      x655 = x1*x654
      x656 = -x646 + x647 + x649 - x651 + x653 + x655
      x657 = -x610*x640 - x627 - x641 - x643 + x644 - x645 + x656
      x658 = x639*x657
      x659 = 0.032153849066577d0*x342*x468
      x660 = 1d0/x54
      x661 = x337*x660
      x662 = x661**2.78d0
      x663 = x468*x69
      x664 = 0.107950572405774d0*x663
      x665 = 0.0616860413747281d0*x476
      x666 = x222*x328*(x228 - x531 - x532 - x533)
      x667 = x234*x329
      x668 = -x539 + x540
      x669 = x241*x331
      x670 = x220*(x1*x542*x668*x669 + x291*(x255*x666 + x329*x599 - x33
     &2 + x598*( -x262 - x566 - x567 - x568 + x623) + x601) + x293*x610 
     &+ x304*( -x604 + x605 - x607 - x608 + x612) + x320*(x305 - x614 - 
     &x615 - x616 - x617 + x625 + x626) + x330*x666 + x530 + x536*x667 -
     & x537* x667 - x54*x545 + x557*x669*(x251*(x547 - x550 + x551 - x55
     &3) + x556))
      x671 = -x335 + x492 + x637
      x672 = x660*x671
      x673 = 4.0d0*x493
      x674 = x334*x642
      x675 = 16.0d0*x42
      x676 = 10.2595256049682d0*x41
      x677 = x495*x676
      x678 = 1d0/x671
      x679 = 0.128615396266308d0*x663*x678
      x680 = x335 + x609 + x635 - x636
      x681 = x646 - x647 - x649
      x682 = exp(x58*(x187 + x336))
      x683 = x682 - 1.0d0
      x684 = 1d0/x683
      x685 = x481*x682*x684
      x686 = 0.00289153318944038d0*x70
      x687 = x684*x75
      x688 = x686*x687
      x689 = x633*x665
      x690 = x341*x343
      x691 = x415*x65
      x692 = x690*x691
      x693 = 0.00142693788768175d0*x633
      x694 = 0.128615396266308d0*x639
      x695 = x482 + x632
      x696 = x485*x671*x695 + x627 + x641 + x643 - x644 + x645 + x651 - 
     &x653 - x655 + x681
      x697 = x338*x340*x480
      x698 = 0.00214040683152263d0*gamma_cc*x415*x506*x508
      x699 = x345*x510*x55/x347
      x700 = x348*x462
      x701 = -711097.316705946d0*x0*x218*x321*x631 - 524288.0d0*x0*x627 
     &+ 262144.0d0* x326*x699*(x341*x346*x686*(x326*x468*x696*x697 + x32
     &6*x471*x689 - x326*x657*x663*x694 + x326*x664 + x338*x344*x696*x69
     &8 + x340*x343 *x505*x693 - 0.00297516549581645d0*x658*x692 + 0.002
     &49714130344307d0*x692) - x342*x469 - x342*x472*x634 + x658* x659 -
     & x662*x688*(-x471*x473*x665 + x664 - x679*(x656 - x670 + x672*(x43
     &9 - x453 + x497 - x563) + x673 - x674 - x675 + x677) + x685*(x486*
     &x680 - x488*x650 + x488*x652 + 0.202643992197415d0* x491 + x670 - 
     &x673 + x674 + x675 - x677 + x681))) - 2097152.0d0* x434*x559 + x43
     &5*x628 + x440*x629 + 524288.0d0*x466*x700 - x521* x630 - x522*x80 
     &+ x523*x630 + x525
      x702 = 4.76837158203125d-7*rho_c
      x703 = x324 + x349
      x704 = -x202
      x705 = x206*(x201 + x704)
      x706 = -x705 + 32768.0d0
      x707 = x136*(x387 - x388)
      x708 = x404*(-x377 + x379 + x403*(x383 + x386 + x394*(-x391*x707 +
     & x392) + x401 *(-x382 - x385 - x397 - x398 + x399*x707 + x400*x707
     &))) - x405 + x406
      x709 = x176*(x376 + x708)
      x710 = x412*x440 - x422*(x414*(-gamma_cs*rho_c + x504) + 0.3204999
     &04512757d0* x416*(-x417 + x418))
      x711 = x426*x710
      x712 = x181*x711
      x713 = rho_s**3*x6
      x714 = x432*x713
      x715 = x7*x714
      x716 = -x441 + x442
      x717 = 0.0506609980493537d0*x716
      x718 = x50*x717
      x719 = x448*x716
      x720 = x7*x713
      x721 = 4.0d0*x41
      x722 = x0*x462
      x723 = x483*x722
      x724 = x468*x634
      x725 = x634*x69
      x726 = x465*x484
      x727 = x155*x714
      x728 = x155*x713
      x729 = x721*x728
      x730 = x719*x9
      x731 = x727 - x729 - x730
      x732 = x188*(x187*x726 + x489*x717 + x731)
      x733 = x478*x74
      x734 = x190*x732*x733
      x735 = x634*x663
      x736 = 0.000356734471920438d0*x503*x633
      x737 = x468*x734
      x738 = x507*x63*x732
      x739 = x511*(x191*x724 + 0.0115661327577615d0*x502*(x725 - x734) +
     & x509*(-x735 - x736 + x737 + x738))
      x740 = x0*x718 + x434*x719 + x461*x723 - x715 + x720*x721 + 0.125d
     &0*x739
      x741 = x515*x710
      x742 = -x185*x740 - x741
      x743 = -x518*x740 - 0.04375d0*x741
      x744 = -x199
      x745 = x182 + x744
      x746 = 64.0d0*x745
      x747 = x59 - 1.0d0
      x748 = 1d0/x747
      x749 = x748*x76 + 1.0d0
      x750 = x443*x80
      x751 = x750*x9
      x752 = 2.56488140124205d0*x751
      x753 = 4.0d0*x80
      x754 = x728*x753
      x755 = x748*x75
      x756 = 1d0/(x68/x747**2 + x749)
      x757 = x72*x749
      x758 = x756*x757
      x759 = x75*x758 + 1.0d0
      x760 = x756/x759
      x761 = 262144.0d0*x760
      x762 = log(x759)
      x763 = x434*x750
      x764 = x0*x629 - x510*x761*(0.0115661327577615d0*x415*x755*(-x59*x
     &733*x748*( x445 - x54*x726 + x727 - x752 + x754) + x725) - x501*x7
     &49*x756*( x735 + x736 - x737 - x738) + x724*x749) + x524*x720 + x6
     &28*x720 + 524288.0d0*x695*x722*x762 - 5378946.16037756d0*x763
      x765 = x1*x231
      x766 = x529*x765
      x767 = 36.0978292789602d0*x546 + x552 + 54.281640924821d0*x555
      x768 = x251*(16.7102894282335d0*x546 + x552)
      x769 = x239*x557
      x770 = -x727
      x771 = -x445 + x752 - x754 + x770
      x772 = x1*x296
      x773 = x603*x772
      x774 = x600*x765
      x775 = x1*x565
      x776 = x572*x579
      x777 = x575*x584
      x778 = x2*x590*x591
      x779 = x4*x595*x596
      x780 = x581 - x585 + x776 - x777 + x778 - x779
      x781 = x613*x772
      x782 = x622*(-x620 + x621)
      x783 = x1*x619
      x784 = x220*(x291*(-x598*(-x775 + x780) + x774) + x293*x771 + x304
     &*(x611*x771 + x773) + x320*(x319*x771 - x624*(x780 - x782 - x783) 
     &+ x781) + x40 *x769*(x541*x767 + x668*x768) - x544 + x766)
      x785 = 0.25d0*x0
      x786 = x326*x342
      x787 = x648*x728
      x788 = -x787
      x789 = 10.2595256049682d0*x751
      x790 = 16.0d0*x728
      x791 = x790*x80
      x792 = x718 + x729 + x730 + x770
      x793 = -x640*x792 - x654 - x784 + x788 + x789 - x791
      x794 = x661**5.56d0
      x795 = x691*x794
      x796 = 0.202643992197415d0*x716
      x797 = x676*x716*x9
      x798 = x41*x790
      x799 = -x581 + x585 - x776 + x777 - x778 + x779
      x800 = x220*(x291*(x598*(x775 + x799) + x774) + x293*x792 + x304*(
     &x611*x792 + x773) + x320*(x319*x792 + x624*(x782 + x783 + x799) + 
     &x781) + x331*x668*x769*(x767 - x768) + x543*x669 + x766)
      x801 = x50*x796 + x672*(-x718 + x731) + x788 + x797 + x798 - x800
      x802 = x489*x796 + x680*x726 + x787 - x797 - x798 + x800
      x803 = x0*x445 + x695*x700*x785 - 0.125d0*x699*(x326*x639*x659*x79
     &3 - x345*x688 *x786*(gamma_cc*x479*x684*x693*x795 + 0.061686041374
     &7281d0*x470* x633*x662*x663 + x662*x679*x801 - x662*x685*x802 + 0.
     &00297516549581645d0*x678*x687*x795*x801 - x682*x698*x794*x802/ x68
     &3**2) + 0.00289153318944038d0*x415*x690*(x689 + x69*x694*x793 - x6
     &97*(x654 - x671*x726 + x784 + x787 - x789 + x791)) + x470* x724*x7
     &86) + x715 + x720*x753 - 2.56488140124205d0*x763 + x784* x785
      x804 = 0.000267550853940328d0*x66
      x805 = -x500*x758*(x500*x804 + x72) + x67*x755 + x757
      x806 = x39*x55
      x807 = x805*x806
      x808 = x761*x807
      x809 = x106*x112
      x810 = 6.68877134850821d-5*x809
      x811 = x371*x83
      x812 = x811*(x115 - x370*(x110 + x810) + 3.34438567425411d-5*x809)
     &
      x813 = x139*x144
      x814 = 6.68877134850821d-5*x813
      x815 = x121*x402
      x816 = x815*(x147 - x401*(x141 + x814) + 3.34438567425411d-5*x813)
     &
      x817 = x812 + x816
      x818 = x176*x39
      x819 = x817*x818
      x820 = x175*x39
      x821 = gamma_cc*x153
      x822 = x83*(-x119 + x98)
      x823 = x121*(x133 - x151)
      x824 = x822 + x823
      x825 = x179*x824
      x826 = x820*x821*x825
      x827 = 0.022020833726518d0*x825
      x828 = 0.240374928384568d0*x827
      x829 = x163/x171**5
      x830 = x154*x161*x829
      x831 = rho_c**(-7.66666666666667d0)*x830
      x832 = 0.120187464192284d0*x168
      x833 = x164*x419*x421
      x834 = x762*x79 + x82
      x835 = x173*x821*x834
      x836 = x760*x807
      x837 = 0.125d0*x836
      x838 = x185*x837 - 2.0d0*x835*(-x160 + x832*x833)
      x839 = 0.13125d0*x199
      x840 = 0.0105164031168249d0*x168
      x841 = x827*x840
      x842 = 0.0875d0*x160
      x843 = x518*x837 - x835*(x833*x840 - x842)
      x844 = x205*x61
      x845 = 0.000244140625d0*x213
      x846 = x351/x350
      x847 = x846**5.56d0
      x848 = -x338 + 1.0d0
      x849 = x72*x846**2.78d0
      x850 = x75/x848
      x851 = -x849*x850 + 1.0d0
      x852 = 1d0/(x68*x847/x848**2 + x851)
      x853 = x849*x851
      x854 = x852*x853
      x855 = 64.0d0*x806*x852*(-x341*x854*(x327 + x690*x804) - x67*x847*
     &x850 + x853)/ (x75*x854 + 1.0d0)
      x856 = 0.001953125d0*rho_c
      x857 = x415*x829
      x858 = x762*x78
      x859 = -5.6d0*x208 + 0.0553050111671826d0*x51
      x860 = -2.8d0*x81 + 0.35d0*x858 + x859
      x861 = x170*x829
      x862 = 0.00375585825600888d0*x860*x861
      x863 = x817*x820
      x864 = x185*x834
      x865 = 64.0d0*x206
      x866 = x865*(x176*x180*x824 + x864)
      x867 = x705 - 32768.0d0
      x868 = x168*x857
      x869 = x834*x861
      x870 = rho_c*x154
      x871 = x214*x870
      x872 = x420*x829
      x873 = x827*x872
      x874 = 0.0210328062336497d0*x168
      x875 = x167*x868
      x876 = rho_s*x875
      x877 = -0.00751171651201775d0*x860*x876
      x878 = 0.0115661327577615d0*x69
      x879 = x108*x878
      x880 = x811*(x113*x879 - x370*(0.000133775426970164d0*x809 + x879)
     & + x810)
      x881 = x140*x878
      x882 = x815*(x145*x881 - x401*(-0.000133775426970164d0*x139*x143*(
     &-gamma_cc - gamma_ss + x99) + x881) + x814)
      x883 = x820*(x880 - x882)
      x884 = 0.480749856769136d0*x168
      x885 = x834*x876
      x886 = x0*x174
      x887 = x181*x886
      x888 = x173*x198
      x889 = -4.2d0*x199
      x890 = x173*(0.35d0*x209 + 2.8d0*x42 + x859)
      x891 = gamma_cc*x0/tau_c
      x892 = 9.5367431640625d-7*x0*x204*x213*x61
      x893 = x121*x198 + x177*x40
      x894 = 0.5d0*x893
      x895 = x176*x894
      x896 = x199 - x895
      x897 = x176*x893
      x898 = 1.4d0*x897
      x899 = x202 - x898
      x900 = x206*x899 + 32768.0d0
      x901 = x210 - x896*x900
      x902 = x214*x901
      x903 = x404*x513
      x904 = 0.0625d0*x209
      x905 = x158*x6
      x906 = x41*x905
      x907 = 0.0098758948512826d0*x51
      x908 = 1.0d0*x208
      x909 = x907 - x908
      x910 = x904 + x906 + x909
      x911 = x176*(x376 + x903 + x910)
      x912 = x410*x893
      x913 = 1.4d0*x427
      x914 = x206*(x704 + x898)
      x915 = -x914 + 32768.0d0
      x916 = x176*(x376 + x404*x740 - x904 - x906 - x907 + x908)
      x917 = x744 + x895
      x918 = 0.021875d0*x711
      x919 = x121*x834 - x40*x822
      x920 = x175*x821
      x921 = x919*x920
      x922 = 0.0625d0*x55*x760*x805
      x923 = x818*(x121*x922 + 0.25d0*x812)
      x924 = x167*x919
      x925 = x831*x924
      x926 = 0.00525820155841243d0*x168
      x927 = 64.0d0*x844
      x928 = x812*x820
      x929 = x868*x924
      x930 = -x864
      x931 = 0.5d0*x176
      x932 = x206*(x919*x931 + x930)
      x933 = x914 - 32768.0d0
      x934 = 0.00375585825600888d0*x869
      x935 = 3.0517578125d-5*x213*x870
      x936 = x820*x880
      x937 = x872*x924
      x938 = -0.00751171651201775d0*x885
      x939 = 0.00375585825600888d0*x168
      x940 = 1.4d0*x886
      x941 = x178*x40 + x198*x83
      x942 = 0.5d0*x941
      x943 = x176*x942
      x944 = x199 - x943
      x945 = x176*x941
      x946 = 1.4d0*x945
      x947 = x202 - x946
      x948 = x206*x947 + 32768.0d0
      x949 = x210 - x944*x948
      x950 = x214*x949
      x951 = x373*x513
      x952 = x176*(x407 + x910 + x951)
      x953 = x410*x941
      x954 = x206*(x704 + x946)
      x955 = -x954 + 32768.0d0
      x956 = x176*(x373*x740 + x708 + x910)
      x957 = x744 + x943
      x958 = -x40*x823 + x83*x834
      x959 = x920*x958
      x960 = x818*(0.25d0*x816 + x83*x922)
      x961 = x167*x958
      x962 = x831*x961
      x963 = x816*x820
      x964 = x875*x958
      x965 = x206*(x930 + x931*x958)
      x966 = x954 - 32768.0d0
      x967 = x820*x882
      x968 = x872*x961
      x969 = x174*x910
      x970 = 2.0d0*x969
      x971 = x183*x970
      x972 = x199 - x971
      x973 = 5.6d0*x969
      x974 = x183*x973
      x975 = x202 - x974
      x976 = x206*x975 + 32768.0d0
      x977 = x210 - x972*x976
      x978 = x214*x977
      x979 = x176*(x197 - 2.0d0*x208 + 1.0d0*x42 + x52 + x903 + x951)
      x980 = x425*x969
      x981 = x424*x514*x910
      x982 = -x206*(x704 + x974) + 32768.0d0
      x983 = 0.0253304990246769d0*x0*x50*x716 + 0.125d0*x196*x723 + 1.28
     &244070062102d0*x41*x434*x716 + 2.0d0*x41*x720 - 0.0395035794051304
     &d0*x49*x720 + 0.0625d0*x739
      x984 = x174*x183
      x985 = x710*x910
      x986 = x744 + x971
      x987 = -x80*x905 + 0.0625d0*x858 + x909
      x988 = x160*x174*x821*x987
      x989 = x167*x419*x830*x987
      x990 = x183*x969
      PBE_value = x79*log(x76*x77/(x77 + x68/x60**2) + 1.0d0) + x82
      PBE_beta_value = x120*x40
      PBE_alpha_value = x152*x40
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x215 + x352)
      d1E(1) = -x215 + x702*(-x211*x527 - x213*(-x200*x206*(-x203*x517 -
     & 2.8d0*x408 + 8.4d0*x411 - 2.8d0*x428 + x519) - x207*(-x408 + 3.0d
     &0*x411 - x428 + x516) + x526) + x701) + x703
      d1E(2) = -rho_c*(-x214*(x206*x746*(0.04375d0*x709 + 0.04375d0*x712
     & + x743) - x706 *(x709 + x712 + x742) + x764) + x803)
      d1E(3) = x856*(-x845*(-x706*(-x168*x828*x831 - 0.25d0*x819 + 2.0d0
     &*x826 + x838) + x746*x844*(gamma_cc*(-0.0109375d0*x819 + 0.0875d0*
     &x826 - x831* x841 + x843) - 0.13125d0*x182 + x839) + x808) + x855)
     &
      d1E(4) = x871*(-x866*(x841*x857 + x862 + 0.0109375d0*x863) + x867*
     &(x828*x868 + 0.25d0*x863 + 0.240374928384568d0*x869))
      d1E(5) = x871*(-x866*(-x873*x874 + x877 + 0.0109375d0*x883) + x867
     &*(-x873*x884 + 0.25d0*x883 - 0.480749856769136d0*x885))
      d1E(6) = -x892*(x706*(-x887 + x888) + x745*x891*(4.2d0*x182 - x183
     &*(-2.8d0*x887 + x890) + x889))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x352 + x902)
      d1E(1) = x702*(-x213*(-x206*x896*(-x517*x899 + x519 - x893*x913 - 
     &2.8d0*x911 + 4.2d0*x912) + x526 - x900*(-x427*x894 + x516 - x911 +
     & 1.5d0*x912 )) - x527*x901 + x701) + x703 - x902
      d1E(2) = -rho_c*(-x214*(x764 + x865*x917*(x743 + x893*x918 + 0.043
     &75d0*x916) - x915*(x711*x894 + x742 + x916)) + x803)
      d1E(3) = x856*(-x845*(x808 - x915*(x832*x925 + x838 - 1.0d0*x921 -
     & x923) + x917* x927*(gamma_cc*(x843 - 0.04375d0*x921 - 0.04375d0*x
     &923 + x925* x926) + x839 - 0.065625d0*x897)) + x855)
      d1E(4) = x935*(x932*(x862 + 0.0109375d0*x928 - 0.00525820155841243
     &d0*x929) + x933 *(0.00390625d0*x928 - 0.00187792912800444d0*x929 +
     & x934))
      d1E(5) = x935*(x932*(x840*x937 + x877 + 0.0109375d0*x936) + x933*(
     &0.00390625d0* x936 + x937*x939 + x938))
      d1E(6) = -x892*(x891*x917*(-x183*(x890 - x893*x940) + x889 + 2.1d0
     &*x897) + x915*( -x886*x894 + x888))
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x352 + x950)
      d1E(1) = x702*(-x213*(-x206*x944*(-x517*x947 + x519 - x913*x941 - 
     &2.8d0*x952 + 4.2d0*x953) + x526 - x948*(-x427*x942 + x516 - x952 +
     & 1.5d0*x953 )) - x527*x949 + x701) + x703 - x950
      d1E(2) = -rho_c*(-x214*(x764 + x865*x957*(x743 + x918*x941 + 0.043
     &75d0*x956) - x955*(x711*x942 + x742 + x956)) + x803)
      d1E(3) = x856*(-x845*(x808 + x927*x957*(gamma_cc*(x843 + x926*x962
     & - 0.04375d0* x959 - 0.04375d0*x960) + x839 - 0.065625d0*x945) - x
     &955*(x832* x962 + x838 - 1.0d0*x959 - x960)) + x855)
      d1E(4) = x935*(x965*(x862 + 0.0109375d0*x963 - 0.00525820155841243
     &d0*x964) + x966 *(x934 + 0.00390625d0*x963 - 0.00187792912800444d0
     &*x964))
      d1E(5) = x935*(x965*(x840*x968 + x877 - 0.0109375d0*x967) + x966*(
     &x938 + x939* x968 - 0.00390625d0*x967))
      d1E(6) = -x892*(x891*x957*(-x183*(x890 - x940*x941) + x889 + 2.1d0
     &*x945) + x955*( -x886*x942 + x888))
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x352 + x978)
      d1E(1) = x702*(-x213*(-x206*x972*(-x517*x975 + x519 - 2.8d0*x979 +
     & 16.8d0*x980 - 5.6d0*x981) + x526 - x976*(x516 - x979 + 6.0d0*x980
     & - 2.0d0*x981 )) - x527*x977 + x701) + x703 - x978
      d1E(2) = -rho_c*(-x214*(x764 + x865*x986*(x154*x172*x842*x985 + x4
     &09*x842*x983 + x743) - x982*(2.0d0*x514*x985 + x742 + 2.0d0*x983*x
     &984)) + x803)
      d1E(3) = x856*(-x845*(x808 + x927*x986*(gamma_cc*(-0.00546875d0*x8
     &36*x984 + x843 + x874*x989 - 0.175d0*x988) + x839 - 0.2625d0*x990)
     & - x982*(-x837 *x984 + x838 + x884*x989 - 4.0d0*x988)) + x855)
      d1E(4) = 0
      d1E(5) = 0
      d1E(6) = -x892*(x891*x986*(-x183*(x890 - x973) + x889 + 8.4d0*x990
     &) + x982*(x888 - x970))
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, g
     &amma_ss, gamma_cs, tau_c, mu, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: January 19, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, gamma_ss, gamma_cs, 
     &tau_c, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
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
     &524, x1525, x1526, x1527, x1528, x1529, x1530, x1531, x1532, x1533
     &, x1534, x1535, x1536, x1537, x1538, x1539, x1540, x1541, x1542, x
     &1543, x1544, x1545, x1546, x1547, x1548, x1549, x1550, x1551, x155
     &2, x1553, x1554, x1555, x1556, x1557, x1558, x1559, x1560, x1561, 
     &x1562, x1563, x1564, x1565, x1566, x1567, x1568, x1569, x1570, x15
     &71, x1572, x1573, x1574, x1575, x1576, x1577, x1578, x1579, x1580,
     & x1581, x1582, x1583, x1584, x1585, x1586, x1587, x1588, x1589, x1
     &590, x1591, x1592, x1593, x1594, x1595, x1596, x1597, x1598, x1599
     &, x1600, x1601, x1602, x1603, x1604, x1605, x1606, x1607, x1608, x
     &1609, x1610, x1611, x1612, x1613, x1614, x1615, x1616, x1617, x161
     &8, x1619, x1620, x1621, x1622, x1623, x1624, x1625, x1626, x1627, 
     &x1628, x1629, x1630, x1631, x1632, x1633, x1634, x1635, x1636, x16
     &37, x1638, x1639, x1640, x1641, x1642, x1643, x1644, x1645, x1646,
     & x1647, x1648, x1649, x1650, x1651, x1652, x1653, x1654, x1655, x1
     &656, x1657, x1658, x1659, x1660, x1661, x1662, x1663, x1664, x1665
     &, x1666, x1667, x1668, x1669, x1670, x1671, x1672, x1673, x1674, x
     &1675, x1676, x1677, x1678, x1679, x1680, x1681, x1682, x1683, x168
     &4, x1685, x1686, x1687, x1688, x1689, x1690, x1691, x1692, x1693, 
     &x1694, x1695, x1696, x1697, x1698, x1699, x1700, x1701, x1702, x17
     &03, x1704, x1705, x1706, x1707, x1708, x1709, x1710, x1711, x1712,
     & x1713, x1714, x1715, x1716, x1717, x1718, x1719, x1720, x1721, x1
     &722, x1723, x1724, x1725, x1726, x1727, x1728, x1729, x1730, x1731
     &, x1732, x1733, x1734, x1735, x1736, x1737, x1738, x1739, x1740, x
     &1741, x1742, x1743, x1744, x1745, x1746, x1747, x1748, x1749, x175
     &0, x1751, x1752, x1753, x1754, x1755, x1756, x1757, x1758, x1759, 
     &x1760, x1761, x1762, x1763, x1764, x1765, x1766, x1767, x1768, x17
     &69, x1770, x1771, x1772, x1773, x1774, x1775, x1776, x1777, x1778,
     & x1779, x1780, x1781, x1782, x1783, x1784, x1785, x1786, x1787, x1
     &788, x1789, x1790, x1791, x1792, x1793, x1794, x1795, x1796, x1797
     &, x1798, x1799, x1800, x1801, x1802, x1803, x1804, x1805, x1806, x
     &1807, x1808, x1809, x1810, x1811, x1812, x1813, x1814, x1815, x181
     &6, x1817, x1818, x1819, x1820, x1821, x1822, x1823, x1824, x1825, 
     &x1826, x1827, x1828, x1829, x1830, x1831, x1832, x1833, x1834, x18
     &35, x1836, x1837, x1838, x1839, x1840, x1841, x1842, x1843, x1844,
     & x1845, x1846, x1847, x1848, x1849, x1850, x1851, x1852, x1853, x1
     &854, x1855, x1856, x1857, x1858, x1859, x1860, x1861, x1862, x1863
     &, x1864, x1865, x1866, x1867, x1868, x1869, x1870, x1871, x1872, x
     &1873, x1874, x1875, x1876, x1877, x1878, x1879, x1880, x1881, x188
     &2, x1883, x1884, x1885, x1886, x1887, x1888, x1889, x1890, x1891, 
     &x1892, x1893, x1894, x1895, x1896, x1897, x1898, x1899, x1900, x19
     &01, x1902, x1903, x1904, x1905, x1906, x1907, x1908, x1909, x1910,
     & x1911, x1912, x1913, x1914, x1915, x1916, x1917, x1918, x1919, x1
     &920, x1921, x1922, x1923, x1924, x1925, x1926, x1927, x1928, x1929
     &, x1930, x1931, x1932, x1933, x1934, x1935, x1936, x1937, x1938, x
     &1939, x1940, x1941, x1942, x1943, x1944, x1945, x1946, x1947, x194
     &8, x1949, x1950, x1951, x1952, x1953, x1954, x1955, x1956, x1957, 
     &x1958, x1959, x1960, x1961, x1962, x1963, x1964, x1965, x1966, x19
     &67, x1968, x1969, x1970, x1971, x1972, x1973, x1974, x1975, x1976,
     & x1977, x1978, x1979, x1980, x1981, x1982, x1983, x1984, x1985, x1
     &986, x1987, x1988, x1989, x1990, x1991, x1992, x1993, x1994, x1995
     &, x1996, x1997, x1998, x1999, x2000, x2001, x2002, x2003, x2004, x
     &2005, x2006, x2007, x2008, x2009, x2010, x2011, x2012, x2013, x201
     &4, x2015, x2016, x2017, x2018, x2019, x2020, x2021, x2022, x2023, 
     &x2024, x2025, x2026, x2027, x2028, x2029, x2030, x2031, x2032, x20
     &33, x2034, x2035, x2036, x2037, x2038, x2039, x2040, x2041, x2042,
     & x2043, x2044, x2045, x2046, x2047, x2048, x2049, x2050, x2051, x2
     &052, x2053, x2054, x2055, x2056, x2057, x2058, x2059, x2060, x2061
     &, x2062, x2063, x2064, x2065, x2066, x2067, x2068, x2069, x2070, x
     &2071, x2072, x2073, x2074, x2075, x2076, x2077, x2078, x2079, x208
     &0, x2081, x2082, x2083, x2084, x2085, x2086, x2087, x2088, x2089, 
     &x2090, x2091, x2092, x2093, x2094, x2095, x2096, x2097, x2098, x20
     &99, x2100, x2101, x2102, x2103, x2104, x2105, x2106, x2107, x2108,
     & x2109, x2110, x2111, x2112, x2113, x2114, x2115, x2116, x2117, x2
     &118, x2119, x2120, x2121, x2122, x2123, x2124, x2125, x2126, x2127
     &, x2128, x2129, x2130, x2131, x2132, x2133, x2134, x2135, x2136, x
     &2137, x2138, x2139, x2140, x2141, x2142, x2143, x2144, x2145, x214
     &6, x2147, x2148, x2149, x2150, x2151, x2152, x2153, x2154, x2155, 
     &x2156, x2157, x2158, x2159, x2160, x2161, x2162, x2163, x2164, x21
     &65, x2166, x2167, x2168, x2169, x2170, x2171, x2172, x2173, x2174,
     & x2175, x2176, x2177, x2178, x2179, x2180, x2181, x2182, x2183, x2
     &184, x2185, x2186, x2187, x2188, x2189, x2190, x2191, x2192, x2193
     &, x2194, x2195, x2196, x2197, x2198, x2199, x2200, x2201, x2202, x
     &2203, x2204, x2205, x2206, x2207, x2208, x2209, x2210, x2211, x221
     &2, x2213, x2214, x2215, x2216, x2217, x2218, x2219, x2220, x2221, 
     &x2222, x2223, x2224, x2225, x2226, x2227, x2228, x2229, x2230, x22
     &31, x2232, x2233, x2234, x2235, x2236, x2237, x2238, x2239, x2240,
     & x2241, x2242, x2243, x2244, x2245, x2246, x2247, x2248, x2249, x2
     &250, x2251, x2252, x2253, x2254, x2255, x2256, x2257, x2258, x2259
     &, x2260, x2261, x2262, x2263, x2264, x2265, x2266, x2267, x2268, x
     &2269, x2270, x2271, x2272, x2273, x2274, x2275, x2276, x2277, x227
     &8, x2279, x2280, x2281, x2282, x2283, x2284, x2285, x2286, x2287, 
     &x2288, x2289, x2290, x2291, x2292, x2293, x2294, x2295, x2296, x22
     &97, x2298, x2299, x2300, x2301, x2302, x2303, x2304, x2305, x2306,
     & x2307, x2308, x2309, x2310, x2311, x2312, x2313, x2314, x2315, x2
     &316, x2317, x2318, x2319, x2320, x2321, x2322, x2323, x2324, x2325
     &, x2326, x2327, x2328, x2329, x2330, x2331, x2332, x2333, x2334, x
     &2335, x2336, x2337, x2338, x2339, x2340, x2341, x2342, x2343, x234
     &4, x2345, x2346, x2347, x2348, x2349, x2350, x2351, x2352, x2353, 
     &x2354, x2355, x2356, x2357, x2358, x2359, x2360, x2361, x2362, x23
     &63, x2364, x2365, x2366, x2367, x2368, x2369, x2370, x2371, x2372,
     & x2373, x2374, x2375, x2376, x2377, x2378, x2379, x2380, x2381, x2
     &382, x2383, x2384, x2385, x2386, x2387, x2388, x2389, x2390, x2391
     &, x2392, x2393, x2394, x2395, x2396, x2397, x2398, x2399, x2400, x
     &2401, x2402, x2403, x2404, x2405, x2406, x2407, x2408, x2409, x241
     &0, x2411, x2412, x2413, x2414, x2415, x2416, x2417, x2418, x2419, 
     &x2420, x2421, x2422, x2423, x2424, x2425, x2426, x2427, x2428, x24
     &29, x2430, x2431, x2432, x2433, x2434, x2435, x2436, x2437, x2438,
     & x2439, x2440, x2441, x2442, x2443, x2444, x2445, x2446, x2447, x2
     &448, x2449, x2450, x2451, x2452, x2453, x2454, x2455, x2456, x2457
     &, x2458, x2459, x2460, x2461, x2462, x2463, x2464, x2465, x2466, x
     &2467, x2468, x2469, x2470, x2471, x2472, x2473, x2474, x2475, x247
     &6, x2477, x2478, x2479, x2480, x2481, x2482, x2483, x2484, x2485, 
     &x2486, x2487, x2488, x2489, x2490, x2491, x2492, x2493, x2494, x24
     &95, x2496, x2497, x2498, x2499, x2500, x2501, x2502, x2503, x2504,
     & x2505, x2506, x2507, x2508, x2509, x2510, x2511, x2512, x2513, x2
     &514, x2515, x2516, x2517, x2518, x2519, x2520, x2521, x2522, x2523
     &, x2524, x2525, x2526, x2527, x2528, x2529, x2530, x2531, x2532, x
     &2533, x2534, x2535, x2536, x2537, x2538, x2539, x2540, x2541, x254
     &2, x2543, x2544, x2545, x2546, x2547, x2548, x2549, x2550, x2551, 
     &x2552, x2553, x2554, x2555, x2556, x2557, x2558, x2559, x2560, x25
     &61, x2562, x2563, x2564, x2565, x2566, x2567, x2568, x2569, x2570,
     & x2571, x2572, x2573, x2574, x2575, x2576, x2577, x2578, x2579, x2
     &580, x2581, x2582, x2583, x2584, x2585, x2586, x2587, x2588, x2589
     &, x2590, x2591, x2592, x2593, x2594, x2595, x2596, x2597, x2598, x
     &2599, x2600, x2601, x2602, x2603, x2604, x2605, x2606, x2607, x260
     &8, x2609, x2610, x2611, x2612, x2613, x2614, x2615, x2616, x2617, 
     &x2618, x2619, x2620, x2621, x2622, x2623, x2624, x2625, x2626, x26
     &27, x2628, x2629, x2630, x2631, x2632, x2633, x2634, x2635, x2636,
     & x2637, x2638, x2639, x2640, x2641, x2642, x2643, x2644, x2645, x2
     &646, x2647, x2648, x2649, x2650, x2651, x2652, x2653, x2654, x2655
     &, x2656, x2657, x2658, x2659, x2660, x2661, x2662, x2663, x2664, x
     &2665, x2666, x2667, x2668, x2669, x2670, x2671, x2672, x2673, x267
     &4, x2675, x2676, x2677, x2678, x2679, x2680, x2681, x2682, x2683, 
     &x2684, x2685, x2686, x2687, x2688, x2689, x2690, x2691, x2692, x26
     &93, x2694, x2695, x2696, x2697, x2698, x2699, x2700, x2701, x2702,
     & x2703, x2704, x2705, x2706, x2707, x2708, x2709, x2710, x2711, x2
     &712, x2713, x2714, x2715, x2716, x2717, x2718, x2719, x2720, x2721
     &, x2722, x2723, x2724, x2725, x2726, x2727, x2728, x2729, x2730, x
     &2731, x2732, x2733, x2734, x2735, x2736, x2737, x2738, x2739, x274
     &0, x2741, x2742, x2743, x2744, x2745, x2746, x2747, x2748, x2749, 
     &x2750, x2751, x2752, x2753, x2754, x2755, x2756, x2757, x2758, x27
     &59, x2760, x2761, x2762, x2763, x2764, x2765, x2766, x2767, x2768,
     & x2769, x2770, x2771, x2772, x2773, x2774, x2775, x2776, x2777, x2
     &778, x2779, x2780, x2781, x2782, x2783, x2784, x2785, x2786, x2787
     &, x2788, x2789, x2790, x2791, x2792, x2793, x2794, x2795, x2796, x
     &2797, x2798, x2799, x2800, x2801, x2802, x2803, x2804, x2805, x280
     &6, x2807, x2808, x2809, x2810, x2811, x2812, x2813, x2814, x2815, 
     &x2816, x2817, x2818, x2819, x2820, x2821, x2822, x2823, x2824, x28
     &25, x2826, x2827, x2828, x2829, x2830, x2831, x2832, x2833, x2834,
     & x2835, x2836, x2837, x2838, x2839, x2840, x2841, x2842, x2843, x2
     &844, x2845, x2846, x2847, x2848, x2849, x2850, x2851, x2852, x2853
     &, x2854, x2855, x2856, x2857, x2858, x2859, x2860, x2861, x2862, x
     &2863, x2864, x2865, x2866, x2867, x2868, x2869, x2870, x2871, x287
     &2, x2873, x2874, x2875, x2876, x2877, x2878, x2879, x2880, x2881, 
     &x2882, x2883, x2884, x2885, x2886, x2887, x2888, x2889, x2890, x28
     &91, x2892, x2893, x2894, x2895, x2896, x2897, x2898, x2899, x2900,
     & x2901, x2902, x2903, x2904, x2905, x2906, x2907, x2908, x2909, x2
     &910, x2911, x2912, x2913, x2914, x2915, x2916, x2917, x2918, x2919
     &, x2920, x2921, x2922, x2923, x2924, x2925, x2926, x2927, x2928, x
     &2929, x2930, x2931, x2932, x2933, x2934, x2935, x2936, x2937, x293
     &8, x2939, x2940, x2941, x2942, x2943, x2944, x2945, x2946, x2947, 
     &x2948, x2949, x2950, x2951, x2952, x2953, x2954, x2955, x2956, x29
     &57, x2958, x2959, x2960, x2961, x2962, x2963, x2964, x2965, x2966,
     & x2967, x2968, x2969, x2970, x2971, x2972, x2973, x2974, x2975, x2
     &976, x2977, x2978, x2979, x2980, x2981, x2982, x2983, x2984, x2985
     &, x2986, x2987, x2988, x2989, x2990, x2991, x2992, x2993, x2994, x
     &2995, x2996, x2997, x2998, x2999, x3000, x3001, x3002, x3003, x300
     &4, x3005, x3006, x3007, x3008, x3009, x3010, x3011, x3012, x3013, 
     &x3014, x3015, x3016, x3017, x3018, x3019, x3020, x3021, x3022, x30
     &23, x3024, x3025, x3026, x3027, x3028, x3029, x3030, x3031, x3032,
     & x3033, x3034, x3035, x3036, x3037, x3038, x3039, x3040, x3041, x3
     &042, x3043, x3044, x3045, x3046, x3047, x3048, x3049, x3050, x3051
     &, x3052, x3053, x3054, x3055, x3056, x3057, x3058, x3059, x3060, x
     &3061, x3062, x3063, x3064, x3065, x3066, x3067, x3068, x3069, x307
     &0, x3071, x3072, x3073, x3074, x3075, x3076, x3077, x3078, x3079, 
     &x3080, x3081, x3082, x3083, x3084, x3085, x3086, x3087, x3088, x30
     &89, x3090, x3091, x3092, x3093, x3094, x3095, x3096, x3097, x3098,
     & x3099, x3100, x3101, x3102, x3103, x3104, x3105, x3106, x3107, x3
     &108, x3109, x3110, x3111, x3112, x3113, x3114, x3115, x3116, x3117
     &, x3118, x3119, x3120, x3121, x3122, x3123, x3124, x3125, x3126, x
     &3127, x3128, x3129, x3130, x3131, x3132, x3133, x3134, x3135, x313
     &6, x3137, x3138, x3139, x3140, x3141, x3142, x3143, x3144, x3145, 
     &x3146, x3147, x3148, x3149, x3150, x3151, x3152, x3153, x3154, x31
     &55, x3156, x3157, x3158, x3159, x3160, x3161, x3162, x3163, x3164,
     & x3165, x3166, x3167, x3168, x3169, x3170, x3171, x3172, x3173, x3
     &174, x3175, x3176, x3177, x3178, x3179, x3180, x3181, x3182, x3183
     &, x3184, x3185, x3186, x3187, x3188, x3189, x3190, x3191, x3192, x
     &3193, x3194, x3195, x3196, x3197, x3198, x3199, x3200, x3201, x320
     &2, x3203, x3204, x3205, x3206, x3207, x3208, x3209, x3210, x3211, 
     &x3212, x3213, x3214, x3215, x3216, x3217, x3218, x3219, x3220, x32
     &21, x3222, x3223, x3224, x3225, x3226, x3227, x3228, x3229, x3230,
     & x3231, x3232, x3233, x3234, x3235, x3236, x3237, x3238, x3239, x3
     &240, x3241, x3242, x3243, x3244, x3245, x3246, x3247, x3248, x3249
     &, x3250, x3251, x3252, x3253, x3254, x3255, x3256, x3257, x3258, x
     &3259, x3260, x3261, x3262, x3263, x3264, x3265, x3266, x3267, x326
     &8, x3269, x3270, x3271, x3272, x3273, x3274, x3275, x3276, x3277, 
     &x3278, x3279, x3280, x3281, x3282, x3283, x3284, x3285, x3286, x32
     &87, x3288, x3289, x3290, x3291, x3292, x3293, x3294, x3295, x3296,
     & x3297, x3298, x3299, x3300, x3301, x3302, x3303, x3304, x3305, x3
     &306, x3307, x3308, x3309, x3310, x3311, x3312, x3313, x3314, x3315
     &, x3316, x3317, x3318, x3319, x3320, x3321, x3322, x3323, x3324, x
     &3325, x3326, x3327, x3328, x3329, x3330, x3331, x3332, x3333, x333
     &4, x3335, x3336, x3337, x3338, x3339, x3340, x3341, x3342, x3343, 
     &x3344, x3345, x3346, x3347, x3348, x3349, x3350, x3351, x3352, x33
     &53, x3354, x3355, x3356, x3357, x3358, x3359, x3360, x3361, x3362,
     & x3363, x3364, x3365, x3366, x3367, x3368, x3369, x3370, x3371, x3
     &372, x3373, x3374, x3375, x3376, x3377, x3378, x3379, x3380, x3381
     &, x3382, x3383, x3384, x3385, x3386, x3387, x3388, x3389, x3390, x
     &3391, x3392, x3393, x3394, x3395, x3396, x3397, x3398, x3399, x340
     &0, x3401, x3402, x3403, x3404, x3405, x3406, x3407, x3408, x3409, 
     &x3410, x3411, x3412, x3413, x3414, x3415, x3416, x3417, x3418, x34
     &19, x3420, x3421, x3422, x3423, x3424, x3425, x3426, x3427, x3428,
     & x3429, x3430, x3431, x3432, x3433, x3434, x3435, x3436, x3437, x3
     &438, x3439, x3440, x3441, x3442, x3443, x3444, x3445, x3446, x3447
     &, x3448, x3449, x3450, x3451, x3452, x3453, x3454, x3455, x3456, x
     &3457, x3458, x3459, x3460, x3461, x3462, x3463, x3464, x3465, x346
     &6, x3467, x3468, x3469, x3470, x3471, x3472, x3473, x3474, x3475, 
     &x3476, x3477, x3478, x3479, x3480, x3481, x3482, x3483, x3484, x34
     &85, x3486, x3487, x3488, x3489, x3490, x3491, x3492, x3493, x3494,
     & x3495, x3496, x3497, x3498, x3499, x3500, x3501, x3502, x3503, x3
     &504, x3505, x3506, x3507, x3508, x3509, x3510, x3511, x3512, x3513
     &, x3514, x3515, x3516, x3517, x3518, x3519, x3520, x3521, x3522, x
     &3523, x3524, x3525, x3526, x3527, x3528, x3529, x3530, x3531, x353
     &2, x3533, x3534, x3535, x3536, x3537, x3538, x3539, x3540, x3541, 
     &x3542, x3543, x3544, x3545, x3546, x3547, x3548, x3549, x3550, x35
     &51, x3552, x3553, x3554, x3555, x3556, x3557, x3558, x3559, x3560,
     & x3561, x3562, x3563, x3564, x3565, x3566, x3567, x3568, x3569, x3
     &570, x3571, x3572, x3573, x3574, x3575, x3576, x3577, x3578, x3579
     &, x3580, x3581, x3582, x3583, x3584, x3585, x3586, x3587, x3588, x
     &3589, x3590, x3591, x3592, x3593, x3594, x3595, x3596, x3597, x359
     &8, x3599, x3600, x3601, x3602, x3603, x3604, x3605, x3606, x3607, 
     &x3608, x3609, x3610, x3611, x3612, x3613, x3614, x3615, x3616, x36
     &17, x3618, x3619, x3620, x3621, x3622, x3623, x3624, x3625, x3626,
     & x3627, x3628, x3629, x3630, x3631, x3632, x3633, x3634, x3635, x3
     &636, x3637, x3638, x3639, x3640, x3641, x3642, x3643, x3644, x3645
     &, x3646, x3647, x3648, x3649, x3650, x3651, x3652, x3653, x3654, x
     &3655, x3656, x3657, x3658, x3659, x3660, x3661, x3662, x3663, x366
     &4, x3665, x3666, x3667, x3668, x3669, x3670, x3671, x3672, x3673, 
     &x3674, x3675, x3676, x3677, x3678, x3679, x3680, x3681, x3682, x36
     &83, x3684, x3685, x3686, x3687, x3688, x3689, x3690, x3691, x3692,
     & x3693, x3694, x3695, x3696, x3697, x3698, x3699, x3700, x3701, x3
     &702, x3703, x3704, x3705, x3706, x3707, x3708, x3709, x3710, x3711
     &, x3712, x3713, x3714, x3715, x3716, x3717, x3718, x3719, x3720, x
     &3721, x3722, x3723, x3724, x3725, x3726, x3727, x3728, x3729, x373
     &0, x3731, x3732, x3733, x3734, x3735, x3736, x3737, x3738, x3739, 
     &x3740, x3741, x3742, x3743, x3744, x3745, x3746, x3747, x3748, x37
     &49, x3750, x3751, x3752, x3753, x3754, x3755, x3756, x3757, x3758,
     & x3759, x3760, x3761, x3762, x3763, x3764, x3765, x3766, x3767, x3
     &768, x3769, x3770, x3771, x3772, x3773, x3774, x3775, x3776, x3777
     &, x3778, x3779, x3780, x3781, x3782, x3783, x3784, x3785, x3786, x
     &3787, x3788, x3789, x3790, x3791, x3792, x3793, x3794, x3795, x379
     &6, x3797, x3798, x3799, x3800, x3801, x3802, x3803, x3804, x3805, 
     &x3806, x3807, x3808, x3809, x3810, x3811, x3812, x3813, x3814, x38
     &15, x3816, x3817, x3818, x3819, x3820, x3821, x3822, x3823, x3824,
     & x3825, x3826, x3827, x3828, x3829, x3830, x3831, x3832, x3833, x3
     &834, x3835, x3836, x3837, x3838, x3839, x3840, x3841, x3842, x3843
     &, x3844, x3845, x3846, x3847, x3848, x3849, x3850, x3851, x3852, x
     &3853, x3854, x3855, x3856, x3857, x3858, x3859, x3860, x3861, x386
     &2, x3863, x3864, x3865, x3866, x3867, x3868, x3869, x3870, x3871, 
     &x3872, x3873, x3874, x3875, x3876, x3877, x3878, x3879, x3880, x38
     &81, x3882, x3883, x3884, x3885, x3886, x3887, x3888, x3889, x3890,
     & x3891, x3892, x3893, x3894, x3895, x3896, x3897, x3898, x3899, x3
     &900, x3901, x3902, x3903, x3904, x3905, x3906, x3907, x3908, x3909
     &, x3910, x3911, x3912, x3913, x3914, x3915, x3916, x3917, x3918, x
     &3919, x3920, x3921, x3922, x3923, x3924, x3925, x3926, x3927, x392
     &8, x3929, x3930, x3931, x3932, x3933, x3934, x3935, x3936, x3937, 
     &x3938, x3939, x3940, x3941, x3942, x3943, x3944, x3945, x3946, x39
     &47, x3948, x3949, x3950, x3951, x3952, x3953, x3954, x3955, x3956,
     & x3957, x3958, x3959, x3960, x3961, x3962, x3963, x3964, x3965, x3
     &966, x3967, x3968, x3969, x3970, x3971, x3972, x3973, x3974, x3975
     &, x3976, x3977, x3978, x3979, x3980, x3981, x3982, x3983, x3984, x
     &3985, x3986, x3987, x3988, x3989, x3990, x3991, x3992, x3993, x399
     &4, x3995, x3996, x3997, x3998, x3999, x4000, x4001, x4002, x4003, 
     &x4004, x4005, x4006, x4007, x4008, x4009, x4010, x4011, x4012, x40
     &13, x4014, x4015, x4016, x4017, x4018, x4019, x4020, x4021, x4022,
     & x4023, x4024, x4025, x4026, x4027, x4028, x4029, x4030, x4031, x4
     &032, x4033, x4034, x4035, x4036, x4037, x4038, x4039, x4040, x4041
     &, x4042, x4043, x4044, x4045, x4046, x4047, x4048, x4049, x4050, x
     &4051, x4052, x4053, x4054, x4055, x4056, x4057, x4058, x4059, x406
     &0, x4061, x4062, x4063, x4064, x4065, x4066, x4067, x4068, x4069, 
     &x4070, x4071, x4072, x4073, x4074, x4075, x4076, x4077, x4078, x40
     &79, x4080, x4081, x4082, x4083, x4084, x4085, x4086, x4087, x4088,
     & x4089, x4090, x4091, x4092, x4093, x4094, x4095, x4096, x4097, x4
     &098, x4099, x4100, x4101, x4102, x4103, x4104, x4105, x4106, x4107
     &, x4108, x4109, x4110, x4111, x4112, x4113, x4114, x4115, x4116, x
     &4117, x4118, x4119, x4120, x4121, x4122, x4123, x4124, x4125, x412
     &6, x4127, x4128, x4129, x4130, x4131, x4132, x4133, x4134, x4135, 
     &x4136, x4137, x4138, x4139, x4140, x4141, x4142, x4143, x4144, x41
     &45, x4146, x4147, x4148, x4149, x4150, x4151, x4152, x4153, x4154,
     & x4155, x4156, x4157, x4158, x4159, x4160, x4161, x4162, x4163, x4
     &164, x4165, x4166, x4167, x4168, x4169, x4170, x4171, x4172, x4173
     &, x4174, x4175, x4176, x4177, x4178, x4179, x4180, x4181, x4182, x
     &4183, x4184, x4185, x4186, x4187, x4188, x4189, x4190, x4191, x419
     &2, x4193, x4194, x4195, x4196, x4197, x4198, x4199, x4200, x4201, 
     &x4202, x4203, x4204, x4205, x4206, x4207, x4208, x4209, x4210, x42
     &11, x4212, x4213, x4214, x4215, x4216, x4217, x4218, x4219, x4220,
     & x4221, x4222, x4223, x4224, x4225, x4226, x4227, x4228, x4229, x4
     &230, x4231, x4232, x4233, x4234, x4235, x4236, x4237, x4238, x4239
     &, x4240, x4241, x4242, x4243, x4244, x4245, x4246
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 0.693147180559945d0
      x1 = -x0 + 1.0d0
      x2 = 1d0/x1
      x3 = 1d0/rho_c
      x4 = rho_s*x3
      x5 = x4 + 1.0d0
      x6 = -x4 + 1.0d0
      x7 = x5**0.666666666666667d0 + x6**0.666666666666667d0
      x8 = x7**2
      x9 = 1d0/x8
      x10 = x2*x9
      x11 = gamma_cc*x10
      x12 = 1.46459188756152d0
      x13 = 1d0/x12
      x14 = x3**0.333333333333333d0
      x15 = x13*x14
      x16 = 0.186690969707574d0*x15 + 1.0d0
      x17 = 2.14502939711103d0
      x18 = 1d0/x17
      x19 = x3**0.666666666666667d0
      x20 = x18*x19
      x21 = 0.564189583547756d0
      x22 = x21*sqrt(x3)
      x23 = 2.91521471421918d0*x22
      x24 = 0.826307487110758d0
      x25 = x24*x3**0.166666666666667d0
      x26 = 5.63098414909788d0*x15 + 0.516066464547863d0*x20 + x23 + 13.
     &4579137143944d0*x25
      x27 = 9.86960440108936d0
      x28 = x2*x27
      x29 = 1.0d0 + x28/x26
      x30 = log(x29)
      x31 = x16*x30
      x32 = 0.194159335344114d0*x15 + 1.0d0
      x33 = 1.41872281647967d0*x22
      x34 = 3.25955091942229d0*x15 + 0.406913004517529d0*x20 + 7.2401019
     &3431683d0* x25 + x33
      x35 = x28/x34
      x36 = log(0.5d0*x35 + 1.0d0)
      x37 = x32*x36
      x38 = 2.0d0*x37
      x39 = -x31 + x38
      x40 = 1d0/x27
      x41 = x1*x40
      x42 = x5**1.33333333333333d0
      x43 = x6**1.33333333333333d0
      x44 = 1.92366105093154d0*x42 + 1.92366105093154d0*x43 - 3.84732210
     &186307d0
      x45 = rho_c**(-4)
      x46 = rho_s**4
      x47 = x45*x46
      x48 = x44*x47
      x49 = x41*x48
      x50 = x39*x49
      x51 = x38*x41
      x52 = -x47 + 1.0d0
      x53 = x44*x52
      x54 = 0.101077332976288d0*x15 + 1.0d0
      x55 = 0.76232752193529d0*x22
      x56 = 3.29180480994506d0*x15 + 0.410025070949612d0*x20 + 9.8721297
     &2256927d0* x25 + x55
      x57 = 1.0d0 + 29.6085746432167d0/x56
      x58 = log(x57)
      x59 = x54*x58
      x60 = 0.0197517897025652d0*x59
      x61 = x53*x60
      x62 = -x51 + x61
      x63 = x50 + x62
      x64 = x7**3
      x65 = 1d0/x64
      x66 = x28*x65
      x67 = 8.0d0*x66
      x68 = exp(-x63*x67)
      x69 = -x68 + 1.0d0
      x70 = rho_c**(-4.66666666666667d0)
      x71 = 208.945363813799d0
      x72 = x70*x71
      x73 = 0.000133775426970164d0*x72
      x74 = gamma_cc**2
      x75 = x1**(-2)
      x76 = x7**(-4)
      x77 = x75*x76
      x78 = x74*x77
      x79 = x73*x78
      x80 = 1d0/x69
      x81 = 14.4549425392770d0
      x82 = rho_c**(-2.33333333333333d0)
      x83 = x81*x82
      x84 = 0.0115661327577615d0*x83
      x85 = x11*x84
      x86 = -x80*x85 + 1.0d0
      x87 = 1d0/(x86 + x79/x69**2)
      x88 = x84*x86
      x89 = x87*x88
      x90 = x11*x89 + 1.0d0
      x91 = log(x90)
      x92 = x41*x64
      x93 = 0.125d0*x92
      x94 = x91*x93
      x95 = x31 - x38
      x96 = x49*x95
      x97 = x62 - x96
      x98 = x94 + x97
      x99 = rho_c + rho_s
      x100 = 1d0/x99
      x101 = x100**0.333333333333333d0
      x102 = x101*x13
      x103 = 0.235215882559859d0*x102 + 1.0d0
      x104 = x100**0.666666666666667d0
      x105 = 0.819204448708788d0*x18
      x106 = sqrt(x100)
      x107 = 4.12273618607837d0*x21
      x108 = x106*x107
      x109 = x100**0.166666666666667d0
      x110 = 15.10599739383d0*x24
      x111 = 7.09459546107279d0*x102 + x104*x105 + x108 + x109*x110
      x112 = 1.0d0 + x28/x111
      x113 = log(x112)
      x114 = x103*x113
      x115 = 1.0d0*x114
      x116 = 2.0d0*gamma_cs
      x117 = gamma_cc + gamma_ss
      x118 = x116 + x117
      x119 = x118*x2
      x120 = x99**(-4.66666666666667d0)
      x121 = 3.34438567425411d-5*x71
      x122 = x121*x75
      x123 = exp(2.0000000000000004d0*x114)
      x124 = x123 - 1.0d0
      x125 = x124**(-2)
      x126 = x118**2
      x127 = x125*x126
      x128 = x120*x122*x127
      x129 = x99**(-2.33333333333333d0)
      x130 = 0.00578306637888076d0*x81
      x131 = x129*x130
      x132 = 1d0/x124
      x133 = x119*x132
      x134 = x131*x133 + 1.0d0
      x135 = x128 + x134
      x136 = 1d0/x135
      x137 = x131*x134
      x138 = x136*x137
      x139 = x119*x138 + 1.0d0
      x140 = log(x139)
      x141 = 0.5d0*x140
      x142 = -x115 + x141
      x143 = rho_c - rho_s
      x144 = 1d0/x143
      x145 = x144**0.333333333333333d0
      x146 = x13*x145
      x147 = 0.235215882559859d0*x146 + 1.0d0
      x148 = x144**0.666666666666667d0
      x149 = sqrt(x144)
      x150 = x107*x149
      x151 = x144**0.166666666666667d0
      x152 = x105*x148 + x110*x151 + 7.09459546107279d0*x146 + x150
      x153 = 1.0d0 + x28/x152
      x154 = log(x153)
      x155 = x147*x154
      x156 = 1.0d0*x155
      x157 = -x116 + x117
      x158 = x157*x2
      x159 = x143**(-4.66666666666667d0)
      x160 = exp(2.0000000000000004d0*x155)
      x161 = x160 - 1.0d0
      x162 = x161**(-2)
      x163 = x157**2
      x164 = x162*x163
      x165 = x122*x159*x164
      x166 = x143**(-2.33333333333333d0)
      x167 = x130*x166
      x168 = 1d0/x161
      x169 = x168*x2
      x170 = x157*x169
      x171 = x167*x170 + 1.0d0
      x172 = x165 + x171
      x173 = 1d0/x172
      x174 = x167*x171
      x175 = x173*x174
      x176 = x158*x175 + 1.0d0
      x177 = log(x176)
      x178 = 0.5d0*x177
      x179 = -x156 + x178
      x180 = rho_c**(-3)
      x181 = rho_c**(-6)
      x182 = rho_s**6*x181
      x183 = 0.5d0*x47
      x184 = rho_c**2
      x185 = 1d0/x184
      x186 = rho_s**2
      x187 = x185*x186
      x188 = 2.26d0*x182 + x183 + 0.87d0*x187 + 0.53d0
      x189 = gamma_cc*x186
      x190 = rho_c*rho_s
      x191 = gamma_ss*x184 - x116*x190 + x189
      x192 = 4.60115111447049d0
      x193 = 1d0/x192
      x194 = 1d0/x43 + 1d0/x42
      x195 = x194*x70
      x196 = x193*x195
      x197 = 0.060093732096142d0*x191*x196 + 1.0d0
      x198 = x197**(-4)
      x199 = x188*x198
      x200 = x199 + 1.0d0
      x201 = x180*x200
      x202 = tau_c**(-2)
      x203 = x202*x74
      x204 = x201*x203
      x205 = x142*x99
      x206 = x143*x179
      x207 = x205 + x206
      x208 = -0.5d0*x0 + 0.5d0
      x209 = x208*x40
      x210 = x207*x209
      x211 = x204*x210
      x212 = x185*x203
      x213 = x199*x212
      x214 = x213 + 64.0d0
      x215 = x51 - x61
      x216 = x215 - x50
      x217 = exp(x216*x67)
      x218 = x217 - 1.0d0
      x219 = x218**(-2)
      x220 = 1d0/x218
      x221 = x220*x85 + 1.0d0
      x222 = x219*x79 + x221
      x223 = 1d0/x222
      x224 = x221*x84
      x225 = x223*x224
      x226 = x11*x225 + 1.0d0
      x227 = log(x226)
      x228 = x227*x93
      x229 = x228 + x63
      x230 = x214*x229
      x231 = -x211 + x230
      x232 = 2.8d0*x211
      x233 = 2.8d0*x230
      x234 = -x232 + x233
      x235 = gamma_cc**3
      x236 = tau_c**(-3)
      x237 = x180*x236
      x238 = x235*x237
      x239 = x234*x238 + 32768.0d0
      x240 = x37*x41
      x241 = x227*x92
      x242 = x53*x59
      x243 = -4194304.0d0*x240 + 262144.0d0*x241 + 41422.505278314d0*x24
     &2 + 2097152.0d0*x50
      x244 = -x231*x239 + x243
      x245 = mu*x18
      x246 = 1.00537384780842d0*x245
      x247 = exp(-rho_c**(-0.33333333333333331d0)*x246)
      x248 = 4.76837158203125d-7*x247
      x249 = x244*x248
      x250 = mu**2
      x251 = x20*x250
      x252 = 0.508616435555896d0*x251 + 1.0d0
      x253 = x252**(-4)
      x254 = x187 - 1.0d0
      x255 = exp(-0.6836107611867116d0*x15)
      x256 = x3**1.33333333333333d0
      x257 = x193*x256
      x258 = 3.14159265358979d0
      x259 = x3**1.0d0
      x260 = x259/x258
      x261 = 0.009578475d0*x260
      x262 = -x261
      x263 = 0.0188171922990732d0*x15 + 0.0676322201645715d0*x20 + 0.001
     &26674656487366d0*x257 + x262 + 1.0d0
      x264 = x255*x263
      x265 = x254*x264
      x266 = x3**1.66666666666667d0
      x267 = 0.0837224526465878d0*x266
      x268 = 1.41421356237310d0
      x269 = mu**5*x268
      x270 = x267*x269
      x271 = 0.0892278228691848d0*x270
      x272 = x5**0.666666666666667d0
      x273 = x6**0.666666666666667d0
      x274 = x272 + x273
      x275 = x274**3
      x276 = x274**2
      x277 = 1d0/x276
      x278 = x250*x277
      x279 = x15*x278
      x280 = 1d0/x274
      x281 = mu*x280
      x282 = 11.1447260721495d0*x25*x281 + 1.0d0
      x283 = 12.5327170711751d0*x279 + x282
      x284 = 1d0/x283
      x285 = mu**3
      x286 = 1d0/x275
      x287 = x22*x286
      x288 = 27.0733719592201d0*x279 + x282 + 27.1408204624105d0*x285*x2
     &87
      x289 = x284*x288
      x290 = log(x289)
      x291 = x275*x290
      x292 = x40*(x0 - 1.0d0)
      x293 = 0.179587122125167d0*x259
      x294 = x264*x293
      x295 = 0.175432441092201d0*x294
      x296 = x254*x295
      x297 = exp(-0.28165369188898165d0*x15)
      x298 = x297*(-4.49737346725955d0*x15 + 0.825481812223657d0*x20)
      x299 = x258/x259
      x300 = x298*x299
      x301 = -x187 + 1.0d0
      x302 = 12.0d0*x301
      x303 = x300*x302
      x304 = x5**2
      x305 = 1d0/x5
      x306 = x305**0.333333333333333d0
      x307 = x15*x306
      x308 = 0.0259335011650457d0*x307
      x309 = -x308 + 1.0d0
      x310 = x305**0.666666666666667d0
      x311 = 0.0524148278841779d0*x20
      x312 = 0.494402081358784d0*x307 + x310*x311 + 1.0d0
      x313 = 1d0/x312
      x314 = 1d0/x310
      x315 = x313*x314
      x316 = x309*x315
      x317 = x304*x316
      x318 = 1d0/x19
      x319 = 2.28942848510666d0
      x320 = x192*x319
      x321 = x318*x320
      x322 = 0.817704266774463d0*x321
      x323 = x6**2
      x324 = 1d0/x6
      x325 = x324**0.333333333333333d0
      x326 = x15*x325
      x327 = 0.0259335011650457d0*x326
      x328 = -x327 + 1.0d0
      x329 = x324**0.666666666666667d0
      x330 = x311*x329 + 0.494402081358784d0*x326 + 1.0d0
      x331 = 1d0/x330
      x332 = 1d0/x329
      x333 = x331*x332
      x334 = x328*x333
      x335 = x323*x334
      x336 = x267*(x303 + x317*x322 + x322*x335)
      x337 = 0.00111534778586481d0*x336
      x338 = x268*x285
      x339 = 4.0d0*x338
      x340 = 0.0472353356922751d0
      x341 = mu**8*x3**2.66666666666667d0*x340
      x342 = 0.267683468607554d0*x341
      x343 = x264 - 1.0d0
      x344 = 0.148394183600699d0
      x345 = x266*x344
      x346 = x343*x345
      x347 = x301*x346
      x348 = 0.0167302167879722d0*x347
      x349 = 0.101321183642338d0*x3**2.0d0
      x350 = x349*x63
      x351 = 0.526297323276602d0*x350
      x352 = mu**6
      x353 = 4.0d0*x352
      x354 = x260*x343
      x355 = 0.13157433081915d0*x301*x354
      x356 = -x355
      x357 = x320*(x5**2.66666666666667d0 + x6**2.66666666666667d0)
      x358 = 0.817704266774463d0*x357
      x359 = 0.558025705063192d0*x20
      x360 = 0.352521395009435d0*x15
      x361 = exp(-0.49698248213959023d0*x15)
      x362 = x361*(x359 - x360)
      x363 = x301*x362
      x364 = x17*x363
      x365 = 10.9027235569928d0*x364
      x366 = 0.817704266774463d0*x320
      x367 = x317*x366
      x368 = x335*x366
      x369 = x358 - x365 - x367 - x368
      x370 = x3**1.0d0*x344
      x371 = 0.00139418473233101d0*x370
      x372 = 1.55214407110551d0*x257
      x373 = mu**4
      x374 = 4.0d0*x373
      x375 = x265*x271 + x291*x292 + x339*(x296 - x337) + x342*x63 + x35
     &3*(-x348 + x351) + x374*(x356 + x369*x371 + x372*x63)
      x376 = x253*x375
      x377 = 0.25d0*x376
      x378 = -x377 + x63
      x379 = exp(-x378*x67)
      x380 = x379 - 1.0d0
      x381 = x380**(-2)
      x382 = 1d0/x97
      x383 = x264*x301
      x384 = x295*x301
      x385 = -x358 + x365 + x367 + x368
      x386 = -x271*x383 - x291*x41 - x339*(x337 + x384) + x342*x97 - x35
     &3*(x348 - x351) + x374*(x356 - x371*x385 + x372*x97)
      x387 = x253*x386
      x388 = 0.25d0*x387
      x389 = -x388 + x97
      x390 = x382*x389
      x391 = x390**5.56d0
      x392 = x381*x391
      x393 = x390**2.78d0
      x394 = x393*x84
      x395 = 1d0/x380
      x396 = x11*x395
      x397 = x394*x396 + 1.0d0
      x398 = x392*x79 + x397
      x399 = 1d0/x398
      x400 = x394*x397
      x401 = x399*x400
      x402 = x11*x401 + 1.0d0
      x403 = log(x402)
      x404 = 0.125d0*x403*x92
      x405 = x215 + x96
      x406 = x388 + x405
      x407 = -x404 + x406
      x408 = x99**(-3.33333333333333d0)
      x409 = 0.0115661327577615d0*x408
      x410 = x99**(-3.33333333333333d0)
      x411 = 0.00192768879296025d0*x410
      x412 = x102*x113
      x413 = 0.546136299139192d0*x18
      x414 = 2.06136809303918d0*x21
      x415 = 2.517666232305d0*x24
      x416 = 2.36486515369093d0*x102 + x104*x413 + x106*x414 + x109*x415
     &
      x417 = 1d0/x112
      x418 = x417/x111**2
      x419 = x416*x418
      x420 = 0.470431765119717d0*x102 + 2.0d0
      x421 = x28*x420
      x422 = -0.156810588373239d0*x412 + x419*x421
      x423 = x123*x132
      x424 = x422*x423
      x425 = x408*x424
      x426 = x409 + x411 + 0.00578306637888076d0*x425
      x427 = x133*x71
      x428 = x426*x427
      x429 = 0.00578306637888076d0*x129
      x430 = x428*x429
      x431 = x409*x81
      x432 = x134*x431
      x433 = x411*x81
      x434 = x134*x433
      x435 = 0.000133775426970164d0*x71
      x436 = x99**(-5.66666666666667d0)
      x437 = x133*x436
      x438 = x435*x437
      x439 = 2.22959044950274d-5*x71
      x440 = x133*x99**(-5.66666666666667d0)
      x441 = x439*x440
      x442 = 6.68877134850821d-5*x71
      x443 = x119*x123*x125*x442
      x444 = x422*x436
      x445 = x130*x425 + x431 + x433
      x446 = x438 + x441 + x443*x444 + x445
      x447 = x133*x138
      x448 = x446*x447
      x449 = -x432 - x434 + x448
      x450 = -x430 + x449
      x451 = 0.5d0*x450
      x452 = x136/x139
      x453 = x40*x452
      x454 = x118*x453
      x455 = x1/x81
      x456 = 0.0784052941866196d0*x455
      x457 = x101*x113
      x458 = x100*x456*x457
      x459 = 1.0d0*x103
      x460 = x100*x419*x459
      x461 = x458 - x460
      x462 = x451*x454 + x461
      x463 = 0.5d0*x99
      x464 = x462*x463
      x465 = 0.25d0*x41
      x466 = x140*x465
      x467 = 0.5d0*x41
      x468 = x114*x467
      x469 = x466 - x468
      x470 = x464 + x469
      x471 = x173/x176
      x472 = x40*x471
      x473 = x157*x472
      x474 = x146*x154
      x475 = 0.156810588373239d0*x474
      x476 = 2.36486515369093d0*x146 + x148*x413 + x149*x414 + x151*x415
     &
      x477 = 1d0/x153
      x478 = x477/x152**2
      x479 = x476*x478
      x480 = 0.470431765119717d0*x146 + 2.0d0
      x481 = x28*x480
      x482 = x479*x481
      x483 = -x475 + x482
      x484 = x160*x168
      x485 = x483*x484
      x486 = x143**(-3.33333333333333d0)
      x487 = 0.00578306637888076d0*x486
      x488 = 0.0115661327577615d0*x486
      x489 = x143**(-3.33333333333333d0)
      x490 = 0.00192768879296025d0*x489
      x491 = x488 + x490
      x492 = x485*x487 + x491
      x493 = x170*x71
      x494 = x492*x493
      x495 = 0.00578306637888076d0*x166
      x496 = x494*x495
      x497 = x488*x81
      x498 = x171*x497
      x499 = x490*x81
      x500 = x171*x499
      x501 = x158*x160*x162
      x502 = x442*x501
      x503 = x143**(-5.66666666666667d0)
      x504 = x483*x503
      x505 = x170*x503
      x506 = x435*x505
      x507 = x143**(-5.66666666666667d0)
      x508 = x170*x507
      x509 = x439*x508
      x510 = x506 + x509
      x511 = x130*x486
      x512 = x497 + x499
      x513 = x485*x511 + x512
      x514 = x502*x504 + x510 + x513
      x515 = x170*x175
      x516 = x514*x515
      x517 = -x498 - x500 + x516
      x518 = -x496 + x517
      x519 = 0.5d0*x518
      x520 = x145*x154
      x521 = x144*x456*x520
      x522 = 1.0d0*x147
      x523 = x144*x479*x522
      x524 = x521 - x523
      x525 = x473*x519 + x524
      x526 = 0.5d0*x143
      x527 = x177*x465
      x528 = x155*x467
      x529 = x527 - x528
      x530 = x525*x526 + x529
      x531 = x470 + x530
      x532 = x204*x531
      x533 = x200*x45
      x534 = x203*x533
      x535 = 3.0d0*x534
      x536 = 2.0d0*x187 + 13.56d0*x47 + 1.74d0
      x537 = x186*x536
      x538 = x180*x537
      x539 = gamma_ss*rho_c
      x540 = gamma_cs*rho_s
      x541 = x539 - x540
      x542 = 0.480749856769136d0*x195
      x543 = x541*x542
      x544 = rho_c**(-5.66666666666667d0)
      x545 = x194*x544
      x546 = 1.12174966579465d0*x545
      x547 = x191*x546
      x548 = rho_c**(-6.66666666666667d0)
      x549 = x191*x548
      x550 = x5**(-2.33333333333333d0)
      x551 = x6**(-2.33333333333333d0)
      x552 = x550 - x551
      x553 = 0.320499904512757d0*rho_s
      x554 = x552*x553
      x555 = x549*x554
      x556 = -x543 + x547 - x555
      x557 = 1d0/x197
      x558 = x188*x193
      x559 = x557*x558
      x560 = x556*x559
      x561 = -x538 + x560
      x562 = x180*x198
      x563 = x203*x562
      x564 = x561*x563
      x565 = x210*x564
      x566 = rho_c**(-5)
      x567 = x46*x566
      x568 = x44*x567
      x569 = x39*x41
      x570 = x568*x569
      x571 = 2.0d0*x3
      x572 = 1.08651697314076d0*x15 + 0.27127533634502d0*x20 + 0.7093614
     &08239834d0* x22 + 1.20668365571947d0*x25
      x573 = x35 + 2.0d0
      x574 = 1d0/x573
      x575 = x574/x34**2
      x576 = x572*x575
      x577 = x32*x576
      x578 = 0.0790071588102608d0*x59
      x579 = x568*x578
      x580 = 0.129439556896076d0*x14
      x581 = x36*x455
      x582 = x580*x581
      x583 = x3*x582
      x584 = x6**0.333333333333333d0
      x585 = x5**0.333333333333333d0
      x586 = x584 - x585
      x587 = 0.0506609980493537d0*x59
      x588 = x52*x587
      x589 = x586*x588
      x590 = rho_s*x185
      x591 = x589*x590
      x592 = x15*x58
      x593 = 0.000665486074881265d0*x592
      x594 = x53*x593
      x595 = -x3*x594
      x596 = 1.09726826998169d0*x15 + 0.273350047299742d0*x20 + 0.381163
     &760967645d0* x22 + 1.64535495376155d0*x25
      x597 = 1d0/x57
      x598 = x597/x56**2
      x599 = x596*x598
      x600 = 0.58482233974552d0*x599
      x601 = x53*x54
      x602 = x600*x601
      x603 = x3*x602
      x604 = -x571*x577 + x579 + x583 + x591 + x595 + x603
      x605 = 1.87699471636596d0*x15 + 0.344044309698576d0*x20 + 1.457607
     &35710959d0* x22 + 2.24298561906574d0*x25
      x606 = 1d0/x29
      x607 = x26**(-2)
      x608 = x16*x606*x607
      x609 = x605*x608
      x610 = 0.062230323235858d0*x14
      x611 = x30*x455
      x612 = x610*x611
      x613 = 2.0d0*x32
      x614 = x576*x613
      x615 = -x582 + x614
      x616 = -x609 + x612 + x615
      x617 = x44*x616
      x618 = x567*x617
      x619 = rho_c**(-3.33333333333333d0)
      x620 = 0.0269876431014436d0*x619
      x621 = x221*x620
      x622 = 1d0/x7
      x623 = rho_s*x622
      x624 = rho_c**(-4.33333333333333d0)
      x625 = x221*x624
      x626 = x5**(-0.333333333333333d0)
      x627 = x6**(-0.333333333333333d0)
      x628 = x626 - x627
      x629 = 0.015421510343682d0*x628
      x630 = x625*x629
      x631 = x623*x630
      x632 = x620*x81
      x633 = x624*x81
      x634 = x622*x633
      x635 = rho_s*x634
      x636 = x629*x635
      x637 = 2.0d0*x216
      x638 = 1.0d0*x627
      x639 = x626 - x638
      x640 = x4*x622
      x641 = x639*x640
      x642 = x47 - 1.0d0
      x643 = x44*x642
      x644 = x54*x643
      x645 = x4*x642
      x646 = x586*x645
      x647 = 4.0d0*x50
      x648 = x47*x617
      x649 = x48*x578
      x650 = x566*x569
      x651 = rho_s**5
      x652 = x586*x651
      x653 = 2.56488140124205d0*x652
      x654 = x650*x653
      x655 = x615 + x647 - x648 - x649 - x654
      x656 = x587*x646 - x593*x643 + x600*x644 + x655
      x657 = x637*x641 + x656
      x658 = x217*x220
      x659 = x657*x658
      x660 = 142.664564503142d0
      x661 = 0.0925290620620922d0*x660
      x662 = x2*x65
      x663 = x661*x662
      x664 = x619*x663
      x665 = x659*x664
      x666 = -x632 + x636 - x665
      x667 = x666*x82
      x668 = x11*x220
      x669 = 0.0115661327577615d0*x668
      x670 = x667*x669
      x671 = 0.0115661327577615d0*x221
      x672 = x671*x82
      x673 = x668*x71
      x674 = x544*x673
      x675 = 2062.20808228389d0
      x676 = gamma_cc*x217
      x677 = x675*x676
      x678 = 0.00214040683152263d0*x677
      x679 = x219*x657
      x680 = x7**(-5)
      x681 = x680*x75
      x682 = x544*x681
      x683 = x679*x682
      x684 = gamma_cc*rho_s
      x685 = x628*x684
      x686 = x220*x685
      x687 = x548*x71
      x688 = x662*x687
      x689 = 0.000356734471920438d0*x688
      x690 = x632 - x636 + x665
      x691 = 0.000624285325860766d0*x674 + x678*x683 - x686*x689 + x690
      x692 = x223*x691
      x693 = x668*x692
      x694 = x672*x693
      x695 = -x621 + x631 + x670 + x694
      x696 = x223/x226
      x697 = 1.46459188756152d0
      x698 = gamma_cc*x697
      x699 = x696*x698
      x700 = x699*x7
      x701 = x695*x700
      x702 = x181*x569
      x703 = -1.0d0*x626 + x627
      x704 = x590*x703
      x705 = 0.25d0*x704
      x706 = x41*x8
      x707 = x227*x706
      x708 = x618 + x653*x702 + 0.125d0*x701 + x705*x707
      x709 = -4.0d0*x570 + x604 + x708
      x710 = x188*x571
      x711 = x561 - x710
      x712 = x198*x229
      x713 = x212*x712
      x714 = x214*x709 + x711*x713
      x715 = x210*x535 - x532 - x565 + x714
      x716 = 3.0d0*x234
      x717 = 8.4d0*x210
      x718 = 0.04375d0*x213 + 2.8d0
      x719 = x709*x718
      x720 = 0.04375d0*x538
      x721 = x3*(0.19775d0*x182 + 0.076125d0*x187 + 0.04375d0*x47 + 0.04
     &6375d0)
      x722 = 0.04375d0*x560
      x723 = -x720 - x721 + x722
      x724 = 64.0d0*x212
      x725 = x712*x723*x724 + 64.0d0*x719
      x726 = -x3*x716 - 2.8d0*x532 + x534*x717 - 2.8d0*x565 + x725
      x727 = x238*x726
      x728 = x592*x643
      x729 = 1395.62545290939d0*x3
      x730 = 5378946.16037756d0*x652
      x731 = x185*x59
      x732 = x642*x731
      x733 = 106243.813381198d0*rho_s
      x734 = 524288.0d0*x704
      x735 = 1226461.339442d0*x3
      x736 = x3*x577
      x737 = 165690.021113256d0*x59
      x738 = x14*x581
      x739 = x3*x738
      x740 = x568*x737 - 4194304.0d0*x736 + 271454.42562372d0*x739
      x741 = -8388608.0d0*x570 - x586*x732*x733 - x599*x644*x735 + 20971
     &52.0d0*x618 + 262144.0d0*x701 + x702*x730 + x707*x734 + x728*x729 
     &+ x740
      x742 = x247*(-x231*x727 - x239*x715 + x741)
      x743 = rho_c**(-1.33333333333333d0)
      x744 = x245*x743
      x745 = 0.33512461593614d0*x744
      x746 = x244*x247
      x747 = x187*x264
      x748 = 0.17845564573837d0*x270
      x749 = -x747*x748
      x750 = 0.00168899541983155d0*x257
      x751 = 0.045088146776381d0*x20
      x752 = 0.00627239743302441d0*x15
      x753 = x262 + x750 + x751 + x752
      x754 = x255*x301
      x755 = x753*x754
      x756 = x269*x383
      x757 = 0.0571643564037363d0*x3**2.0d0
      x758 = 0.0203323666368788d0*x757
      x759 = 0.148713038115308d0*x270
      x760 = x276*x41
      x761 = x5**(-0.333333333333333d0)
      x762 = x6**(-0.333333333333333d0)
      x763 = x761 - x762
      x764 = 2.0d0*x290
      x765 = x763*x764
      x766 = x760*x765
      x767 = 0.713822582953479d0*x341
      x768 = x15*x281
      x769 = 4.17757235705837d0*x768
      x770 = mu*x15
      x771 = x277*x770
      x772 = x4*x763
      x773 = x771*x772
      x774 = 16.7102894282335d0*x773
      x775 = 1.85745434535825d0*x25
      x776 = 7.42981738143299d0*x25
      x777 = x280*x4
      x778 = x763*x777
      x779 = x776*x778
      x780 = -x775 + x779
      x781 = -x769 + x774 + x780
      x782 = 9.02445731974005d0*x768
      x783 = x22*x278
      x784 = 13.5704102312052d0*x783
      x785 = 36.0978292789602d0*x773
      x786 = x250*x287
      x787 = x772*x786
      x788 = 54.281640924821d0*x787
      x789 = x780 - x782 - x784 + x785 + x788
      x790 = -x289*x781 + x789
      x791 = 1d0/x288
      x792 = mu*x791
      x793 = x760*x792
      x794 = x582 - x614
      x795 = x609 - x612 + x794
      x796 = x48*x795
      x797 = 4.0d0*x96
      x798 = x41*x95
      x799 = x566*x798
      x800 = x4*x589
      x801 = -x594 + x602 + x649 + x794 + x800
      x802 = -x653*x799 - x796 + x797 + x801
      x803 = x297*(1.49912448908652d0*x15 - 0.550321208149104d0*x20)
      x804 = x299*x302
      x805 = x803*x804
      x806 = x298*x301
      x807 = 2.14502939711103d0*x3**(-0.666666666666667d0)
      x808 = x806*x807
      x809 = 1.12661476755593d0*x808
      x810 = 24.0d0*x300
      x811 = x187*x810
      x812 = 1d0/x306
      x813 = x313*x812
      x814 = x304*x813
      x815 = x305*x4
      x816 = -x815 + 1.0d0
      x817 = 1d0/x14
      x818 = x258*x319*x817
      x819 = 0.00706864485168613d0*x818
      x820 = x816*x819
      x821 = x814*x820
      x822 = 1d0/x325
      x823 = x331*x822
      x824 = x323*x823
      x825 = x324*x4
      x826 = x825 + 1.0d0
      x827 = x819*x826
      x828 = x824*x827
      x829 = 0.545136177849642d0*x321
      x830 = x317*x829
      x831 = x335*x829
      x832 = x312**(-2)
      x833 = x314*x832
      x834 = x309*x833
      x835 = x322*x834
      x836 = 0.034943218589452d0*x20
      x837 = x310*x836
      x838 = 0.164800693786261d0*x307
      x839 = x837 + x838
      x840 = -x815*x837 - x815*x838 + x839
      x841 = x304*x840
      x842 = x835*x841
      x843 = x330**(-2)
      x844 = x332*x843
      x845 = x328*x844
      x846 = x322*x845
      x847 = x329*x836
      x848 = 0.164800693786261d0*x326
      x849 = x847 + x848
      x850 = x825*x847 + x825*x848 + x849
      x851 = x323*x850
      x852 = x846*x851
      x853 = 2.18054471139857d0*x321
      x854 = x334*x6
      x855 = x853*x854
      x856 = x4*x855
      x857 = x316*x5
      x858 = x853*x857
      x859 = x4*x858
      x860 = x821 + x828 + x830 + x831 + x842 + x852 + x856 - x859
      x861 = x303 + x805 + x809 + x811 + x860
      x862 = 0.00111534778586481d0*x267
      x863 = 0.175432441092201d0*x293
      x864 = 0.122619224952946d0*x256
      x865 = 0.0399758348639607d0*x864
      x866 = 0.350864882184401d0*x294
      x867 = -x187*x866 + 0.00185891297644135d0*x336
      x868 = 0.0334604335759443d0*x346
      x869 = x187*x868
      x870 = 0.0278836946466203d0*x347
      x871 = x15*x263
      x872 = 0.000844497709915773d0*x257
      x873 = 0.0047892375d0*x260
      x874 = 0.0225440733881905d0*x20
      x875 = 0.0031361987165122d0*x15
      x876 = -x872 + x873 - x874 - x875
      x877 = 0.113935126864452d0*x871 + x876
      x878 = x754*x877
      x879 = 0.0334604335759443d0*x345
      x880 = x878*x879
      x881 = 1.0525946465532d0*x350
      x882 = -x647
      x883 = x648 + x654 + x801 + x882
      x884 = 0.526297323276602d0*x349
      x885 = x883*x884
      x886 = 0.263148661638301d0*x354
      x887 = x187*x886
      x888 = x260*x878
      x889 = 0.263148661638301d0*x888
      x890 = x369*x370
      x891 = 0.00232364122055169d0*x890
      x892 = x257*x63
      x893 = 2.06952542814068d0*x892
      x894 = x318*x357
      x895 = 0.545136177849642d0*x894
      x896 = x318*x364
      x897 = 7.26848237132856d0*x896
      x898 = x17*x318
      x899 = x361*x898
      x900 = 0.117507131669812d0*x15 - 0.372017136708795d0*x20
      x901 = x301*x900
      x902 = x899*x901
      x903 = 10.9027235569928d0*x902
      x904 = x12*x817
      x905 = x363*x904
      x906 = 1.80615420514536d0*x905
      x907 = x362*x898
      x908 = 21.8054471139857d0*x907
      x909 = x187*x908
      x910 = x6**1.66666666666667d0
      x911 = x5**1.66666666666667d0
      x912 = x321*(x910 - x911)
      x913 = x4*x912
      x914 = 2.18054471139857d0*x913
      x915 = -x821 - x828 - x830 - x831 - x842 - x852 - x856 + x859
      x916 = x895 - x897 - x903 - x906 - x909 + x914 + x915
      x917 = 0.00139418473233101d0*x345
      x918 = x916*x917
      x919 = x372*x883
      x920 = x271*x755 + x339*(-x383*x865 + x384 + x755*x863 - x861*x862
     & + x867) + x342*x802 - x353*(x869 - x870 + x880 + x881 - x885) - x
     &374*(x356 + x887 + x889 + x891 + x893 - x918 - x919) + x383*x759 +
     & x4*x766 + x749 - x756*x758 - x767*x97 - x790*x793
      x921 = x253*x920
      x922 = x3*x921
      x923 = x44*x795
      x924 = x567*x923
      x925 = 8388608.0d0*x798
      x926 = x52*x586
      x927 = x731*x926
      x928 = x53*x592
      x929 = x181*x798
      x930 = x251/x252**5
      x931 = x386*x930
      x932 = x3*x931
      x933 = 1d0/x402
      x934 = 262144.0d0*x399
      x935 = x933*x934
      x936 = x698*x7
      x937 = x397*x620
      x938 = -x626
      x939 = x627 + x938
      x940 = 0.015421510343682d0*x939
      x941 = x623*x624
      x942 = x397*x940*x941
      x943 = 8.0d0*x240
      x944 = 0.0790071588102608d0*x242
      x945 = -x943 + x944
      x946 = -x387 - x797 + x945
      x947 = 1d0/x946
      x948 = x382*x946
      x949 = 4.0d0*x796
      x950 = 1.35631049481572d0*x930
      x951 = x386*x950
      x952 = 16.0d0*x96
      x953 = 10.2595256049682d0*x652
      x954 = x799*x953
      x955 = 8.0d0*x577
      x956 = 0.517758227584304d0*x738
      x957 = 0.316028635241043d0*x59
      x958 = x48*x957
      x959 = 0.00266194429952506d0*x592
      x960 = x53*x959
      x961 = 2.33928935898208d0*x599
      x962 = x601*x961
      x963 = 0.202643992197415d0*x59
      x964 = x926*x963
      x965 = x4*x964
      x966 = -x955 + x956 + x958 - x960 + x962 + x965
      x967 = -x921 - x949 - x951 + x952 - x954 + x966
      x968 = -x883*x948 + x967
      x969 = x947*x968
      x970 = 0.032153849066577d0*x619
      x971 = x397*x969*x970
      x972 = 1d0/x63
      x973 = x378*x972
      x974 = x973**2.78d0
      x975 = x619*x81
      x976 = 0.107950572405774d0*x975
      x977 = x628*x635
      x978 = 0.0616860413747281d0*x977
      x979 = 0.128615396266308d0*x975
      x980 = -x376 + x647 + x945
      x981 = 1d0/x980
      x982 = x594 - x602 + x655 - x800
      x983 = x972*x980
      x984 = x255*(x261 - x750 - x751 - x752)
      x985 = x254*x984
      x986 = x265*x269
      x987 = -x761 + x762
      x988 = x276*x292
      x989 = x4*x764
      x990 = x988*x989
      x991 = x775 - x779
      x992 = x769 - x774 + x991
      x993 = x289*x992 + x789
      x994 = x792*x988
      x995 = -x811
      x996 = -x303 - x805 - x809 + x915 + x995
      x997 = -x265*x759 + x271*x985 + x339*(x265*x865 - x296 + x862*x996
     & + x863*x985 + x867) + x342*x883 + x353*(-x869 + x870 - x880 - x88
     &1 + x885) + x374*(x355 - x887 - x889 - x891 - x893 + x918 + x919) 
     &- x63*x767 + x749 + x758*x986 + x987*x990 + x993*x994
      x998 = x253*x997
      x999 = 4.0d0*x648
      x1000 = x375*x950
      x1001 = 16.0d0*x50
      x1002 = x650*x953
      x1003 = -x1000 - x1001 + x1002 + x966 - x998 + x999
      x1004 = x1003 + x982*x983
      x1005 = x1004*x981
      x1006 = x1005*x979
      x1007 = x376 + x882 + x943 - x944
      x1008 = 2.0d0*x641
      x1009 = x955 - x956 - x958
      x1010 = x1000 + x1001 - x1002 + x998 - x999
      x1011 = x1009 + x1010 - x643*x959 + x644*x961 + x646*x963
      x1012 = x1007*x1008 + x1011
      x1013 = x216 + x377
      x1014 = exp(x1013*x67)
      x1015 = x1014 - 1.0d0
      x1016 = 1d0/x1015
      x1017 = x1014*x1016
      x1018 = x1017*x663
      x1019 = x1012*x1018*x619
      x1020 = x974*(-x1006 + x1019 + x976 - x978)
      x1021 = x1016*x11
      x1022 = 0.00289153318944038d0*x82
      x1023 = x1021*x1022
      x1024 = x1020*x1023
      x1025 = x391*x396
      x1026 = x544*x71
      x1027 = 0.00249714130344307d0*x1026
      x1028 = gamma_cc*x675
      x1029 = x1028*x682
      x1030 = 0.00214040683152263d0*x1029
      x1031 = x622*x980
      x1032 = x638 + x938
      x1033 = 2.0d0*x1032
      x1034 = x1009 + x960 - x962 - x965
      x1035 = x1034 + x921 + x949 + x951 - x952 + x954
      x1036 = x1031*x1033*x4 + x1035
      x1037 = x1036*x379*x392
      x1038 = 0.00142693788768175d0*x939
      x1039 = x391*x395
      x1040 = x684*x688
      x1041 = x1039*x1040
      x1042 = x1025*x1026
      x1043 = x1042*x969
      x1044 = x393*x975
      x1045 = x1044*x969
      x1046 = 0.0616860413747281d0*x939
      x1047 = x379*x395
      x1048 = x1036*x1047
      x1049 = -0.128615396266308d0*x1045 + x1046*x393*x635 + x1048*x393*
     &x664 + x393* x976
      x1050 = x1025*x1027 + x1030*x1037 + x1038*x1041 - 0.00297516549581
     &645d0*x1043 + x1049
      x1051 = x1050*x399
      x1052 = x1051*x396
      x1053 = x1022*x1052*x397
      x1054 = x393*(-x1024 + x1053 - x937 - x942 + x971)
      x1055 = x403*x706
      x1056 = 524288.0d0*x185
      x1057 = rho_s*x1056*x703
      x1058 = x599*x601
      x1059 = x1054*x935*x936 + x1055*x1057 + x1058*x735 + x568*x925 - x
     &729*x928 - x730*x929 + x733*x927 + x740 - 524288.0d0*x922 - 209715
     &2.0d0*x924 - 711097.316705946d0*x932
      x1060 = 4.76837158203125d-7*rho_c
      x1061 = x389 + x404
      x1062 = -x233
      x1063 = x1062 + x232
      x1064 = x1063*x238
      x1065 = -x1064 + 32768.0d0
      x1066 = x475 - x482
      x1067 = x1066*x484
      x1068 = x1067*x511
      x1069 = x1066*x503
      x1070 = x1069*x502
      x1071 = -x506 - x509
      x1072 = x1068 + x1070 + x1071 - x497 - x499
      x1073 = -x1067*x487 + x491
      x1074 = x493*x495
      x1075 = x498 + x500
      x1076 = x1073*x1074 + x1075
      x1077 = x1072*x515 + x1076
      x1078 = 0.5d0*x473
      x1079 = x1077*x1078
      x1080 = -x521 + x523
      x1081 = -x527 + x528
      x1082 = x1081 + x526*(x1079 + x1080)
      x1083 = x1082 + x470
      x1084 = x1083*x204
      x1085 = rho_s*x536
      x1086 = x1085*x185
      x1087 = -gamma_cs*rho_c + x684
      x1088 = -x550 + x551
      x1089 = x1088*x191
      x1090 = x1089*x544
      x1091 = x1087*x542 + 0.320499904512757d0*x1090
      x1092 = x1091*x559
      x1093 = x1086 - x1092
      x1094 = x1093*x563
      x1095 = x1094*x210
      x1096 = rho_s**3
      x1097 = x1096*x44
      x1098 = x1097*x45
      x1099 = 4.0d0*x1098
      x1100 = x619*x940
      x1101 = x1100*x221
      x1102 = x81*x940
      x1103 = x10*x661
      x1104 = x622*x703
      x1105 = 2.0d0*x1104
      x1106 = -x584 + x585
      x1107 = x1106*x642
      x1108 = x1097*x578
      x1109 = x1108*x180
      x1110 = 4.0d0*x180
      x1111 = x1097*x1110
      x1112 = x1111*x569
      x1113 = 2.56488140124205d0*x47
      x1114 = x1106*x569
      x1115 = x1113*x1114
      x1116 = x1109 - x1112 - x1115
      x1117 = x1107*x587 + x1116
      x1118 = x1105*x216 + x1117
      x1119 = x1118*x658
      x1120 = x1103*x1119
      x1121 = x1102 - x1120
      x1122 = x1121*x544
      x1123 = x940*x975
      x1124 = x9*x939
      x1125 = gamma_cc*x220
      x1126 = x1124*x1125*x2
      x1127 = x1026*x1126
      x1128 = 0.000356734471920438d0*x1127
      x1129 = x1120*x619
      x1130 = x219*x678
      x1131 = x1130*x77
      x1132 = x1118*x544
      x1133 = x1131*x1132
      x1134 = -x1123 - x1128 + x1129 + x1133
      x1135 = x1134*x223
      x1136 = x668*x672
      x1137 = x1101 + x1122*x669 + x1135*x1136
      x1138 = 0.125d0*x1137
      x1139 = x1138*x699
      x1140 = x1106*x588
      x1141 = x1108*x45
      x1142 = -x1141
      x1143 = 0.25d0*x3
      x1144 = x1143*x639
      x1145 = x1144*x707
      x1146 = 2.56488140124205d0*x567
      x1147 = x1114*x1146
      x1148 = x1139 + x1140*x3 + x1142 + x1145 + x1147
      x1149 = x1099*x569 + x1148
      x1150 = x1149*x214
      x1151 = x1093*x713
      x1152 = -x1150 - x1151
      x1153 = x1084 + x1095 + x1152
      x1154 = 0.04375d0*x1084
      x1155 = 0.04375d0*x1095
      x1156 = x1149*x718
      x1157 = 0.04375d0*x1151
      x1158 = -x1156 - x1157
      x1159 = x1154 + x1155 + x1158
      x1160 = -x230
      x1161 = x1160 + x211
      x1162 = x1161*x238
      x1163 = 64.0d0*x1162
      x1164 = x3*x59
      x1165 = x1164*x926
      x1166 = x68 - 1.0d0
      x1167 = x1166**(-2)
      x1168 = x1167*x79
      x1169 = 1d0/x1166
      x1170 = x1169*x85 + 1.0d0
      x1171 = x1168 + x1170
      x1172 = 1d0/x1171
      x1173 = x1170*x84
      x1174 = x1172*x1173
      x1175 = x11*x1174 + 1.0d0
      x1176 = 1d0/x1175
      x1177 = x1172*x1176
      x1178 = 262144.0d0*x1177
      x1179 = x11*x1169
      x1180 = 0.0115661327577615d0*x1179
      x1181 = x586*x798
      x1182 = x1113*x1181
      x1183 = x1111*x798
      x1184 = x1109 - x1182 + x1183 + x589
      x1185 = -x1105*x63 + x1184
      x1186 = x1169*x68
      x1187 = x1103*x1186
      x1188 = x1102 - x1185*x1187
      x1189 = x1188*x544
      x1190 = x1123 - x1129
      x1191 = x1128 - x1133 + x1190
      x1192 = x669*x82
      x1193 = x1170*x1172
      x1194 = x1192*x1193
      x1195 = x1100*x1170 + x1180*x1189 - x1191*x1194
      x1196 = x1195*x698
      x1197 = log(x1175)
      x1198 = x1197*x706
      x1199 = x1032*x3
      x1200 = 524288.0d0*x1199
      x1201 = x1181*x567
      x1202 = x1098*x737 + x1098*x925 + 106243.813381198d0*x1165 - x1178
     &*x1196 + x1198 *x1200 - 5378946.16037756d0*x1201
      x1203 = -x1065*x1153 + x1159*x1163 + x1202
      x1204 = x1203*x248
      x1205 = x264*x4
      x1206 = x1205*x748
      x1207 = 36.0978292789602d0*x768 + x776 + 54.281640924821d0*x783
      x1208 = x1207*x763
      x1209 = 16.7102894282335d0*x768 + x776
      x1210 = x1209*x289
      x1211 = x1208 + x1210*x987
      x1212 = x274*x41
      x1213 = x1212*x792
      x1214 = -x1109
      x1215 = x1182 - x1183 + x1214 - x589
      x1216 = x4*x868
      x1217 = x4*x866
      x1218 = x4*x810
      x1219 = -x1218
      x1220 = x5*x813
      x1221 = x1220*x819
      x1222 = x6*x823
      x1223 = x1222*x819
      x1224 = x834*x839
      x1225 = x1224*x5
      x1226 = x1225*x322
      x1227 = x845*x849
      x1228 = x1227*x6
      x1229 = x1228*x322
      x1230 = x1221 - x1223 + x1226 - x1229 - x855 + x858
      x1231 = x1219 + x1230
      x1232 = x4*x886
      x1233 = x321*(-x910 + x911)
      x1234 = 2.18054471139857d0*x1233
      x1235 = x4*x908
      x1236 = x1230 - x1234 - x1235
      x1237 = x1206 + x1211*x1213 + x1215*x342 + x339*(x1217 - x1231*x86
     &2) + x353*( x1215*x884 + x1216) + x374*(x1215*x372 + x1232 - x1236
     &*x917) - x766
      x1238 = x1237*x253
      x1239 = 0.25d0*x3
      x1240 = x1238*x1239
      x1241 = x3*x589
      x1242 = x1099*x798
      x1243 = x393*x397
      x1244 = x1243*x622
      x1245 = x1112 + x1115 + x1140 + x1214
      x1246 = x1097*x957
      x1247 = x1246*x180
      x1248 = -x1247
      x1249 = x47*x798
      x1250 = 10.2595256049682d0*x586
      x1251 = x1249*x1250
      x1252 = x1097*x180
      x1253 = 16.0d0*x1252
      x1254 = x1253*x798
      x1255 = -x1238 + x1248 + x1251 - x1254 - x964
      x1256 = -x1245*x948 + x1255
      x1257 = x1256*x947
      x1258 = x1243*x1257
      x1259 = x622*x81
      x1260 = x1046*x1259
      x1261 = 0.128615396266308d0*x81
      x1262 = x1247 + x1254
      x1263 = x1238 - x1251 + x1262 + x964
      x1264 = -x1105*x980 + x1263
      x1265 = x1047*x1264
      x1266 = x1257*x1261 + x1260 - x1265*x663
      x1267 = x1025*x1266
      x1268 = 0.00289153318944038d0*x544
      x1269 = x973**5.56d0
      x1270 = x1016*x1269
      x1271 = gamma_cc*x1270
      x1272 = x1026*x662
      x1273 = x1271*x1272
      x1274 = x1038*x1273
      x1275 = x1116 - x1140
      x1276 = x274*x292
      x1277 = x1276*x792
      x1278 = x1207 - x1210
      x1279 = x1278*x987
      x1280 = -x1221 + x1223 - x1226 + x1229 + x855 - x858
      x1281 = x1218 + x1280
      x1282 = x1234 + x1235 + x1280
      x1283 = x1206 + x1245*x342 + x1277*x1279 + x339*(x1217 + x1281*x86
     &2) + x353*( x1216 + x1245*x884) + x374*(x1232 + x1245*x372 + x1282
     &*x917) + x765*x988
      x1284 = x1283*x253
      x1285 = x1106*x52
      x1286 = x1285*x963
      x1287 = x1253*x569
      x1288 = 10.2595256049682d0*x1114
      x1289 = x1288*x47
      x1290 = x1248 + x1287 + x1289
      x1291 = -x1284 + x1286 + x1290
      x1292 = x1275*x983 + x1291
      x1293 = x1292*x981
      x1294 = x1021*x1269
      x1295 = 0.00297516549581645d0*x1026*x1294
      x1296 = x1293*x1295
      x1297 = x1247 - x1287 - x1289
      x1298 = x1284 + x1297
      x1299 = x1107*x963 + x1298
      x1300 = x1007*x1105 + x1299
      x1301 = x1015**(-2)
      x1302 = x1269*x1301
      x1303 = x1014*x1302
      x1304 = x1030*x1303
      x1305 = x1300*x1304
      x1306 = x622*x974
      x1307 = x1046*x1306*x975
      x1308 = x1293*x974
      x1309 = x1308*x979
      x1310 = x1018*x1300
      x1311 = x1310*x619*x974
      x1312 = x1307 + x1309 - x1311
      x1313 = x1274 + x1296 - x1305 + x1312
      x1314 = x1243*x399
      x1315 = x1313*x1314
      x1316 = -x1023*x1315 + x1100*x1244 + x1258*x970 + x1267*x1268
      x1317 = x399*x933
      x1318 = x1316*x1317
      x1319 = 0.125d0*x936
      x1320 = x1318*x1319
      x1321 = x1032*x1055*x1239
      x1322 = x1146*x1181
      x1323 = x1141 + x1240 + x1241 + x1242 - x1320 + x1321 - x1322
      x1324 = 0.000267550853940328d0*x72
      x1325 = x1324*x668
      x1326 = x1325 + x84
      x1327 = x1326*x668
      x1328 = x1173 - x1174*x1327 + x1179*x73
      x1329 = x1328*x40
      x1330 = x1329*x7
      x1331 = x1178*x1330
      x1332 = x120*x133
      x1333 = x121*x1332
      x1334 = x120*x442
      x1335 = x133*x1334
      x1336 = x131 + x1335
      x1337 = x1336*x136
      x1338 = x133*x1337
      x1339 = x1338*x137
      x1340 = x1333 - x1339 + x137
      x1341 = x1340*x452
      x1342 = x1341*x99
      x1343 = x159*x170
      x1344 = x121*x1343
      x1345 = x159*x442
      x1346 = x1345*x170
      x1347 = x1346 + x167
      x1348 = x1347*x173
      x1349 = x1348*x170
      x1350 = x1349*x174
      x1351 = x1344 - x1350 + x174
      x1352 = x1351*x471
      x1353 = x1352*x143
      x1354 = x1342 + x1353
      x1355 = x1354*x40
      x1356 = x1355*x204
      x1357 = 0.25d0*x1356
      x1358 = x201*x202
      x1359 = gamma_cc*x1358
      x1360 = x99*(x115 - x141)
      x1361 = x143*(x156 - x178)
      x1362 = x1360 + x1361
      x1363 = x1362*x209
      x1364 = x1359*x1363
      x1365 = 2.0d0*x1364
      x1366 = rho_c**(-7.66666666666667d0)
      x1367 = x1362*x208
      x1368 = x1366*x1367
      x1369 = x188*x194
      x1370 = x197**(-5)
      x1371 = 0.0220208337265180d0
      x1372 = x1370*x1371
      x1373 = x1369*x1372
      x1374 = x186*x203
      x1375 = x1373*x1374
      x1376 = 0.240374928384568d0*x1368*x1375
      x1377 = x1197*x93
      x1378 = x1377 + x97
      x1379 = x1378*x199
      x1380 = gamma_cc*x202
      x1381 = x193*x557
      x1382 = x194*x548
      x1383 = x1382*x189
      x1384 = x1381*x1383
      x1385 = 0.120187464192284d0*x1384
      x1386 = 2.0d0*x1385 - 2.0d0*x185
      x1387 = x1380*x1386
      x1388 = x1379*x1387
      x1389 = 0.125d0*x40
      x1390 = x1389*x7
      x1391 = x1177*x1328
      x1392 = x1390*x1391
      x1393 = x1392*x214
      x1394 = -x1388 + x1393
      x1395 = -x1357 + x1365 - x1376 + x1394
      x1396 = 0.13125d0*x230
      x1397 = 0.0109375d0*x1356
      x1398 = 0.0875d0*x1364
      x1399 = 0.0105164031168249d0*x1366
      x1400 = x1367*x1375*x1399
      x1401 = 0.0875d0*x185
      x1402 = 0.0105164031168249d0*x1384
      x1403 = x1380*(-x1401 + x1402)
      x1404 = x1379*x1403
      x1405 = x1392*x718
      x1406 = -x1404 + x1405
      x1407 = gamma_cc*(-x1397 + x1398 - x1400 + x1406)
      x1408 = x1396 + x1407 - 0.13125d0*x211
      x1409 = x237*x74
      x1410 = 64.0d0*x1409
      x1411 = x1408*x1410
      x1412 = x247*(-x1065*x1395 + x1161*x1411 + x1331)
      x1413 = x406/x405
      x1414 = x1413**5.56d0
      x1415 = -x379 + 1.0d0
      x1416 = 1d0/x1415
      x1417 = x1413**2.78d0*x84
      x1418 = -x11*x1416*x1417 + 1.0d0
      x1419 = 1d0/(x1414*x79/x1415**2 + x1418)
      x1420 = x1417*x1418
      x1421 = x1419*x1420
      x1422 = x11*x1421 + 1.0d0
      x1423 = 1d0/x1422
      x1424 = x40*x7
      x1425 = x1423*x1424
      x1426 = x11*x73
      x1427 = x1419*(-x1414*x1416*x1426 + x1420 - x1421*x396*(x1025*x132
     &4 + x394))
      x1428 = x1425*x1427
      x1429 = 64.0d0*x1428
      x1430 = 0.001953125d0*rho_c
      x1431 = x1372*x188
      x1432 = x1431*x545
      x1433 = x1367*x1432
      x1434 = x1197*x92
      x1435 = -5.6d0*x240 + 0.0553050111671826d0*x242
      x1436 = x1435 - 2.8d0*x96
      x1437 = 0.35d0*x1434 + x1436
      x1438 = x195*x558
      x1439 = x1370*x1438
      x1440 = 0.00375585825600888d0*x1439
      x1441 = x1437*x1440
      x1442 = x201*x40
      x1443 = x1354*x1442
      x1444 = 0.0109375d0*x1443
      x1445 = 0.0105164031168249d0*x1433 + x1441 + x1444
      x1446 = x1363*x204
      x1447 = x1378*x214
      x1448 = x1446 + x1447
      x1449 = x1448*x238
      x1450 = x1445*x1449
      x1451 = x1064 - 32768.0d0
      x1452 = 0.240374928384568d0*x1438
      x1453 = x1370*x1452
      x1454 = 0.25d0*x1443
      x1455 = x1378*x1453 + 0.240374928384568d0*x1433 + x1454
      x1456 = x1451*x1455
      x1457 = -64.0d0*x1450 + x1456
      x1458 = x203*x248
      x1459 = rho_c*x1458
      x1460 = rho_s*x1370
      x1461 = x1369*x1460
      x1462 = x1371*x1461
      x1463 = x1462*x548
      x1464 = x1367*x1463
      x1465 = x545*x558
      x1466 = x1460*x1465
      x1467 = 0.00751171651201775d0*x1466
      x1468 = -x1437*x1467
      x1469 = 0.0115661327577615d0*x81
      x1470 = x129*x1469
      x1471 = x134*x1470
      x1472 = x120*x435
      x1473 = x133*x1472
      x1474 = x1470 + x1473
      x1475 = x136*x1474
      x1476 = x137*x1475
      x1477 = -x133*x1476 + x1335 + x1471
      x1478 = x1477*x452
      x1479 = x1478*x99
      x1480 = x1469*x166
      x1481 = x159*x435
      x1482 = -gamma_cc - gamma_ss + x116
      x1483 = x1482*x169
      x1484 = x1480 - x1481*x1483
      x1485 = x1484*x515
      x1486 = x1480*x171
      x1487 = x1346 + x1486
      x1488 = -x1485 + x1487
      x1489 = x1488*x471
      x1490 = x143*x1489
      x1491 = x1479 - x1490
      x1492 = x1442*x1491
      x1493 = -0.0210328062336497d0*x1464 + x1468 + 0.0109375d0*x1492
      x1494 = x1449*x1493
      x1495 = 0.480749856769136d0*x1465
      x1496 = x1460*x1495
      x1497 = -x1378*x1496 - 0.480749856769136d0*x1464 + 0.25d0*x1492
      x1498 = x1451*x1497
      x1499 = -64.0d0*x1494 + x1498
      x1500 = x200*x3
      x1501 = x1500*x209
      x1502 = x1501*x207
      x1503 = x199*x229
      x1504 = -x1502 + x1503
      x1505 = 4.2d0*x211
      x1506 = 4.2d0*x230
      x1507 = -x1506
      x1508 = x1435 + 0.35d0*x241 + 2.8d0*x50
      x1509 = x1508*x199
      x1510 = x212*(-2.8d0*x1502 + x1509)
      x1511 = x1505 + x1507 - x1510
      x1512 = 1d0/tau_c
      x1513 = gamma_cc*x1512
      x1514 = x1513*x3
      x1515 = x1511*x1514
      x1516 = x1065*x1504 + x1161*x1515
      x1517 = 9.5367431640625d-7*x247
      x1518 = x236*x74
      x1519 = x1518*x3
      x1520 = x1517*x1519
      x1521 = x99**(-2)
      x1522 = 0.104540392248826d0*x455
      x1523 = 1d0/x697
      x1524 = x101*x419
      x1525 = 0.0784052941866196d0*x1521
      x1526 = x104*x18
      x1527 = 3.09205213955878d0*x21
      x1528 = x109*x24
      x1529 = x418*(3.15315353825457d0*x102 + x106*x1527 + 0.91022716523
     &1986d0*x1526 + 2.9372772710225d0*x1528)
      x1530 = x1521*x459
      x1531 = x416**2/(x111**4*x112**2)
      x1532 = x416*x417*(4.72973030738186d0*x102 + x108 + 1.092272598278
     &38d0*x1526 + 5.03533246461d0*x1528)/x111**3
      x1533 = x432 + x434 - x448
      x1534 = x1533 + x430
      x1535 = x135**(-2)
      x1536 = x1535/x139**2
      x1537 = x2*x40
      x1538 = x99**(-4.33333333333333d0)
      x1539 = 0.0449794051690726d0*x1538
      x1540 = x1539*x81
      x1541 = 0.0231322655155231d0*x408
      x1542 = 0.00385537758592051d0*x410
      x1543 = x1541*x81
      x1544 = x134*x1543
      x1545 = x136*x446
      x1546 = x133*x1545
      x1547 = x1542*x81
      x1548 = x134*x1547
      x1549 = 0.0115661327577615d0*x125
      x1550 = x1545*x75
      x1551 = x426*x71
      x1552 = 0.000267550853940329d0*x71
      x1553 = 4.45918089900548d-5*x71
      x1554 = x123*x125
      x1555 = x119*x1554*x444
      x1556 = x1552*x437 + x1553*x440 + x1555*x435
      x1557 = x1543 + x1547 + x424*x431
      x1558 = x127*x75
      x1559 = x137*x1535
      x1560 = x1558*x1559*x446
      x1561 = x99**(-4.33333333333333d0)
      x1562 = x1561*x424
      x1563 = 0.00385537758592051d0*x1538*x424
      x1564 = exp(4.0000000000000009d0*x114)
      x1565 = x422**2
      x1566 = x1561*x1565
      x1567 = x1564*x1566
      x1568 = 0.00578306637888076d0*x423
      x1569 = 6.73880859569814d0*x2
      x1570 = 0.313621176746478d0*x1569
      x1571 = 97.4090910340024d0*x75
      x1572 = -x1524*x1570 - x1529*x421 - x1531*x1571*x420 + x1532*x421 
     &+ 0.209080784497652d0*x412
      x1573 = x1561*x1572
      x1574 = x427*x429
      x1575 = x99**(-6.66666666666667d0)
      x1576 = x99**(-6.66666666666667d0)
      x1577 = 0.0231322655155231d0*x81
      x1578 = x130*x423
      x1579 = 0.000535101707880657d0*x71
      x1580 = x119*x1575
      x1581 = x1554*x422
      x1582 = 8.91836179801095d-5*x71
      x1583 = 0.000200663140455246d0*x71
      x1584 = x1575*x443
      x1585 = 0.5d0*x454
      x1586 = x462 + x463*(x102*x1525*x419 + x126*x1534*x1536*x1537*x451
     & - x1521*x1522 *x457 + x1523*x1524*x1525 + x1529*x1530 + x1530*x15
     &31*x28 - x1530 *x1532 + x1585*(-x126*x129*x1549*x1550*x1551 + x134
     &*x1540 + x1541 *x428 + x1542*x428 - x1544*x1546 - x1546*x1548 + x1
     &560*(x1556 + x1557) + x1574*(x1539 + x1549*x1567 + 0.0231322655155
     &231d0*x1562 + x1563 - x1566*x1568 - x1568*x1573) + x447*(-x119*x15
     &76*x1581* x1582 - x125*x1469*x1567 - x1540 - x1562*x1577 - x1563*x
     &81 + x1565*x1584 + x1566*x1578 + x1572*x1584 + x1573*x1578 - 0.000
     &668877134850821d0*x1575*x427 - 0.000215527076785265d0*x1576* x427 
     &- x1579*x1580*x1581 - x1564*x1565*x1580*x1583/x124**3)))
      x1587 = x143**(-4.33333333333333d0)
      x1588 = 0.0449794051690726d0*x1587
      x1589 = x1588*x81
      x1590 = x1589*x171
      x1591 = 0.0231322655155231d0*x486
      x1592 = 0.00385537758592051d0*x489
      x1593 = x1591*x81
      x1594 = x1593*x171
      x1595 = x170*x173
      x1596 = x1595*x514
      x1597 = x1592*x81
      x1598 = x1597*x171
      x1599 = x514*x71
      x1600 = 0.0115661327577615d0*x162
      x1601 = x173*x75
      x1602 = x1600*x1601*x163*x166
      x1603 = x435*x501
      x1604 = x1552*x505
      x1605 = x1553*x508
      x1606 = x1604 + x1605
      x1607 = x1603*x504 + x1606
      x1608 = x1593 + x1597
      x1609 = x1608 + x485*x497
      x1610 = x164*x75
      x1611 = x172**(-2)
      x1612 = x1611*x174
      x1613 = x1610*x1612
      x1614 = x1613*x514
      x1615 = x143**(-4.33333333333333d0)
      x1616 = x483**2
      x1617 = x1615*x1616
      x1618 = exp(4.0000000000000009d0*x155)
      x1619 = x1600*x1618
      x1620 = x1615*x485
      x1621 = x1587*x485
      x1622 = 0.00385537758592051d0*x1621
      x1623 = 0.00578306637888076d0*x484
      x1624 = 0.209080784497652d0*x474
      x1625 = x145*x479
      x1626 = x148*x18
      x1627 = x151*x24
      x1628 = x478*(3.15315353825457d0*x146 + x149*x1527 + 0.91022716523
     &1986d0*x1626 + 2.9372772710225d0*x1627)
      x1629 = x1628*x481
      x1630 = x476**2/(x152**4*x153**2)
      x1631 = x1571*x1630*x480
      x1632 = x476*x477*(4.72973030738186d0*x146 + x150 + 1.092272598278
     &38d0*x1626 + 5.03533246461d0*x1627)/x152**3
      x1633 = x1632*x481
      x1634 = -x1570*x1625 + x1624 - x1629 - x1631 + x1633
      x1635 = x1615*x1623
      x1636 = x1588 - x1634*x1635
      x1637 = x130*x484
      x1638 = x1469*x1618*x162
      x1639 = x143**(-6.66666666666667d0)
      x1640 = x1639*x502
      x1641 = x483*x501
      x1642 = x1639*x1641
      x1643 = x143**(-6.66666666666667d0)
      x1644 = x1641*x1643
      x1645 = x158*x1583*x1618*x1639/x161**3
      x1646 = x1615*x1637
      x1647 = 0.000215527076785265d0*x1643*x493
      x1648 = 0.000668877134850821d0*x1639*x493
      x1649 = -x1589 + x1634*x1640 + x1634*x1646 - x1647 - x1648
      x1650 = x1075 - x516
      x1651 = x1650 + x496
      x1652 = x1611/x176**2
      x1653 = x1537*x163*x1652
      x1654 = x1653*x519
      x1655 = x143**(-2)
      x1656 = x1655*x522
      x1657 = x1628*x1656
      x1658 = x1522*x1655*x520
      x1659 = 0.0784052941866196d0*x1655
      x1660 = x1523*x1625*x1659
      x1661 = x146*x1659*x479
      x1662 = x1632*x1656
      x1663 = x1630*x1656*x28
      x1664 = x1657 - x1658 + x1660 + x1661 - x1662 + x1663
      x1665 = x525 + x526*(x1078*(x1074*(x1617*x1619 - x1617*x1623 + 0.0
     &231322655155231d0*x1620 + x1622 + x1636) + x1590 + x1591*x494 + x1
     &592*x494 - x1594*x1596 - x1596*x1598 - x1599*x1602*x492 + x1614*(x
     &1607 + x1609) + x515*(-x1577*x1620 - x1579*x1642 - x1582* x1644 + 
     &x1616*x1640 - x1616*x1645 + x1617*x1637 - x1617*x1638 - x1622*x81 
     &+ x1649)) + x1651*x1654 + x1664)
      x1666 = x204*(x1586 + x1665)
      x1667 = x531*x534
      x1668 = x198*x561
      x1669 = x1668*x203
      x1670 = x1110*x1669
      x1671 = x200*x203
      x1672 = x1671*x566
      x1673 = x186*x45
      x1674 = x1673*(10.0d0*x187 + 94.92d0*x47 + 5.22d0)
      x1675 = x1381*x556
      x1676 = 0.60093732096142d0*x195
      x1677 = rho_s*x552
      x1678 = x197**(-2)
      x1679 = x188*x340
      x1680 = x1678*x1679
      x1681 = x1680*x556
      x1682 = x1681*(-x1676*x541 - 0.400624880640947d0*x1677*x549 + 1.40
     &218708224331d0 *x191*x545)
      x1683 = x1366*x191
      x1684 = x5**(-3.33333333333333d0)
      x1685 = 3.11111111111111d0*x4
      x1686 = x1684*x1685
      x1687 = x6**(-3.33333333333333d0)
      x1688 = x1685*x1687
      x1689 = x559*(-gamma_ss*x542 - 0.240374928384568d0*rho_s*x1683*(x1
     &686 + x1688 - 2.66666666666667d0*x550 + 2.66666666666667d0*x551) +
     & 2.99133244211907d0*x1677*x1683 - 1.28199961805103d0*x1677*x541* x
     &548 - 6.35658143950302d0*x194*x549 + 4.48699866317861d0*x541* x545
     &)
      x1690 = x1674 - 2.0d0*x1675*x538 + x1682 + x1689
      x1691 = 2.0d0*x198
      x1692 = x1691*x180
      x1693 = x1692*x203
      x1694 = x203*x45
      x1695 = x1694*x210
      x1696 = x1668*x1695
      x1697 = x185*x738
      x1698 = x181*x46
      x1699 = x1698*x44
      x1700 = x1699*x592
      x1701 = rho_c**(-7)
      x1702 = x59*x652
      x1703 = x1701*x1702
      x1704 = x1699*x59
      x1705 = rho_s*x180
      x1706 = x1705*x592
      x1707 = x1706*x926
      x1708 = x185*x928
      x1709 = 0.0395035794051304d0*x59
      x1710 = x5**(-0.666666666666667d0)
      x1711 = 0.854960467080683d0*x4
      x1712 = x1710*x1711
      x1713 = x6**(-0.666666666666667d0)
      x1714 = x1711*x1713
      x1715 = x1712 + x1714
      x1716 = x1715 - 5.1297628024841d0*x584 + 5.1297628024841d0*x585
      x1717 = x1705*x1716
      x1718 = x1717*x52
      x1719 = x1523*x576
      x1720 = x14*x1719
      x1721 = 0.258879113792152d0*x185
      x1722 = x15*x576
      x1723 = x575*(1.44868929752102d0*x15 + 0.452125560575033d0*x20 + 1
     &.06404211235975d0*x22 + 1.40779759833938d0*x25)
      x1724 = x1723*x32
      x1725 = 4.0d0*x185
      x1726 = x28*x572**2/(x34**4*x573**2)
      x1727 = x1726*x32
      x1728 = x54*x599
      x1729 = x1699*x1728
      x1730 = x15*x599
      x1731 = x1730*x53
      x1732 = x1731*x185
      x1733 = x1705*x1728
      x1734 = x1733*x926
      x1735 = x1701*x569
      x1736 = x1735*x652
      x1737 = x598*(1.46302435997558d0*x15 + 0.455583412166236d0*x20 + 0
     &.571745641451468d0*x22 + 1.91958077938847d0*x25)
      x1738 = x1737*x601
      x1739 = x1738*x185
      x1740 = x596**2/(x56**4*x57**2)
      x1741 = x1740*x601
      x1742 = x1741*x185
      x1743 = x1699*x569
      x1744 = x572*x574*(2.17303394628153d0*x15 + 0.542550672690039d0*x2
     &0 + 2.41336731143894d0*x25 + x33)/x34**3
      x1745 = x1744*x32
      x1746 = x1716*x651
      x1747 = x1735*x1746
      x1748 = x596*x597*(2.19453653996337d0*x15 + 0.546700094599483d0*x2
     &0 + 3.29070990752309d0*x25 + x55)/x56**3
      x1749 = x1748*x601
      x1750 = x1749*x185
      x1751 = x616*x652
      x1752 = x1701*x1751
      x1753 = x1698*x617
      x1754 = x41*x7
      x1755 = x1754*x703
      x1756 = x1755*x939
      x1757 = x1673*x1756
      x1758 = x1757*x227
      x1759 = 2.0d0*x627
      x1760 = 2.0d0*x626
      x1761 = x5**(-1.33333333333333d0)
      x1762 = 0.333333333333333d0*x4
      x1763 = x1761*x1762
      x1764 = x6**(-1.33333333333333d0)
      x1765 = x1762*x1764
      x1766 = -x1763 - x1765
      x1767 = x1705*(-x1759 + x1760 + x1766)
      x1768 = x227*x467*x8
      x1769 = x608*(2.50265962182128d0*x15 + 0.573407182830959d0*x20 + 2
     &.18641103566438d0*x22 + 2.6168165555767d0*x25)
      x1770 = 0.0829737643144773d0*x14*x611
      x1771 = x605*x606
      x1772 = x16*x1771*(3.75398943273192d0*x15 + 0.688088619397151d0*x2
     &0 + x23 + 4.48597123813148d0*x25)/x26**3
      x1773 = x1771*x607
      x1774 = x1523*x1773*x610
      x1775 = 0.062230323235858d0*x15*x1773
      x1776 = x16*x28*x605**2/(x26**4*x29**2)
      x1777 = x1723*x613
      x1778 = 0.172586075861435d0*x738
      x1779 = x1744*x613
      x1780 = x1719*x580
      x1781 = 0.129439556896076d0*x1722
      x1782 = x1726*x613
      x1783 = -x1777 + x1778 + x1779 - x1780 - x1781 - x1782
      x1784 = x1769 - x1770 - x1772 + x1774 + x1775 + x1776 + x1783
      x1785 = x1699*x1784
      x1786 = x185*x684*x697*x703
      x1787 = x1786*x695*x696
      x1788 = x222**(-2)
      x1789 = x1788/x226**2
      x1790 = 21.170591578193d0*x74
      x1791 = x1789*x1790
      x1792 = x1791*x695
      x1793 = x2*x622
      x1794 = x1792*x1793*(x621 - x631 - x670 - x694)
      x1795 = 0.0899588103381452d0*x624
      x1796 = rho_c**(-5.33333333333333d0)
      x1797 = x1796*x628
      x1798 = 0.0719670482705162d0*x1797
      x1799 = x628**2
      x1800 = 0.0308430206873641d0*x221
      x1801 = rho_c**(-6.33333333333333d0)*x186
      x1802 = x1801*x9
      x1803 = 2.66666666666667d0*x627
      x1804 = 2.66666666666667d0*x626
      x1805 = 0.444444444444444d0*x4
      x1806 = x1761*x1805
      x1807 = x1764*x1805
      x1808 = x1806 + x1807
      x1809 = x1796*(x1803 - x1804 + x1808)
      x1810 = x1809*x623
      x1811 = x666*x668
      x1812 = 0.0539752862028871d0*x619
      x1813 = x662*x686
      x1814 = 0.0308430206873641d0*x1813
      x1815 = x221*x223
      x1816 = x1815*x668
      x1817 = 0.0231322655155231d0*x219
      x1818 = x667*x78
      x1819 = 0.000713468943840876d0*x688
      x1820 = x660*x662
      x1821 = x1820*x659
      x1822 = 0.185058124124184d0*x619
      x1823 = 0.00428081366304526d0*x677
      x1824 = x219*x78
      x1825 = x1788*x672
      x1826 = x1824*x1825*x691
      x1827 = x1795*x81
      x1828 = rho_s*x1259
      x1829 = x1798*x1828
      x1830 = x81*x9
      x1831 = x1799*x1801*x1830
      x1832 = 0.0308430206873641d0*x1831
      x1833 = x1809*x1828
      x1834 = 0.0115661327577615d0*x1833
      x1835 = x1821*x624
      x1836 = 0.431802289623097d0*x1835
      x1837 = 16.0d0*x66
      x1838 = exp(x1837*x216)
      x1839 = x7**(-6)
      x1840 = 1408.04281369971d0*x75
      x1841 = x1839*x1840
      x1842 = x1841*x624
      x1843 = x657**2
      x1844 = x1843*x219
      x1845 = 1.48046499299348d0*x1838*x1842*x1844
      x1846 = 0.740232496496738d0*x1842*x1843*x658
      x1847 = x2*x76
      x1848 = x1847*x660
      x1849 = rho_s*x1848
      x1850 = 0.246744165498913d0*x1797*x1849*x659
      x1851 = x1763 + x1765
      x1852 = x1851 - 2.0d0*x626 + 2.0d0*x627
      x1853 = x1852*x640
      x1854 = 0.000887314766508353d0*x592
      x1855 = x1716*x60
      x1856 = 5.33333333333333d0*x216
      x1857 = x1856*x639
      x1858 = x187*x628*x9
      x1859 = 17.3157558993759d0*x1740
      x1860 = 0.58482233974552d0*x644
      x1861 = 0.0394081882442864d0*x1730
      x1862 = 0.00341378571249706d0*x592
      x1863 = 2.99999988448829d0*x1728
      x1864 = x1784*x48
      x1865 = 8.0d0*x648
      x1866 = 5.1297628024841d0*x566
      x1867 = x1751*x1866
      x1868 = x48*x59
      x1869 = 0.395035794051304d0*x1868
      x1870 = x1702*x566
      x1871 = 0.40528798439483d0*x1870
      x1872 = x1746*x650
      x1873 = 20.0d0*x50
      x1874 = x48*x592
      x1875 = 0.00532388859905012d0*x1874
      x1876 = 20.5190512099364d0*x652
      x1877 = x1876*x650
      x1878 = x1728*x48
      x1879 = 4.67857871796416d0*x1878
      x1880 = x1783 - x1864 + x1865 - x1867 + x1869 - x1871 - x1872 - x1
     &873 + x1875 + x1877 - x1879
      x1881 = -x1737*x1860 + x1748*x1860 + x1853*x637 + x1854*x643 + x18
     &55*x645 + x1857*x1858 - x1859*x644 - x1861*x643 - x1862*x646 + x18
     &63*x646 + x1880 + 4.0d0*x641*x656
      x1882 = x1881*x624*x658*x663
      x1883 = x548*x673
      x1884 = x1366*x71
      x1885 = x1813*x1884
      x1886 = rho_c**(-8.66666666666667d0)
      x1887 = x1799*x1847*x1886*x189*x71
      x1888 = x220*x684
      x1889 = 0.888888888888889d0*x4
      x1890 = x1761*x1889
      x1891 = x1764*x1889
      x1892 = x1884*x662
      x1893 = x1892*(x1890 + x1891 - 5.33333333333333d0*x626 + 5.3333333
     &3333333d0*x627 )
      x1894 = x548*x681
      x1895 = x1894*x677
      x1896 = x1895*x679
      x1897 = 20353.1779648712d0*x548/x1**3
      x1898 = x1897/x7**8
      x1899 = gamma_cc*x1898
      x1900 = 0.0513697639565431d0*x1899
      x1901 = x1838/x218**3
      x1902 = 0.017123254652181d0*x676
      x1903 = x1839*x675*x75
      x1904 = x1366*x217*x685
      x1905 = x1136*x223
      x1906 = x700*(x1192*(x1827 - x1829 + x1832 + x1834 + x1836 + x1845
     & - x1846 - x1850 - x1882) + x1795*x221 - x1798*x221*x623 + x1799*x
     &1800*x1802 + x1810*x671 - x1811*x1812 - x1812*x1816*x691 + x1814*x
     &624*x666 + x1814*x625*x692 + x1817*x1818*x692 + x1826*(-x1819*x686
     & + x1821* x1822 + x1823*x683 + 0.00124857065172153d0*x674 + 0.0539
     &752862028871d0*x975 - 0.0308430206873641d0*x977) + x1905*( x1130*x
     &1881*x1894 - x1827 + x1829 - x1832 - x1834 - x1836 - x1843 *x1900*
     &x1901 + x1844*x1898*x1902 - x1845 + x1846 + x1850 + x1882 - 0.0035
     &3761684654434d0*x1883 + 0.00332952173792409d0*x1885 - 0.0011891149
     &0640146d0*x1887*x220 - 0.000133775426970164d0*x1888* x1893 - 0.019
     &9771304275445d0*x1896 + 0.011415503101454d0*x1903* x1904*x679))
      x1907 = -0.345172151722869d0*x1697 - 0.0106477771981002d0*x1700 + 
     &0.810575968789659d0*x1703 - 0.790071588102608d0*x1704 - 0.00682757
     &142499412d0*x1707 + 0.00177462953301671d0*x1708 + x1709 *x1718 + x
     &1720*x1721 + x1721*x1722 + x1724*x1725 + x1725*x1727 - x1725*x1745
     & + 9.35715743592833d0*x1729 - 0.0788163764885729d0* x1732 + 5.9999
     &9976897657d0*x1734 - 41.0381024198728d0*x1736 - 1.16964467949104d0
     &*x1739 - 34.6315117987517d0*x1742 + 40.0d0* x1743 + 2.0d0*x1747 + 
     &1.16964467949104d0*x1750 + 10.2595256049682d0*x1752 - 16.0d0*x1753
     & + 0.666666666666667d0* x1758 + x1767*x1768 + 2.0d0*x1785 + x1787 
     &+ 0.25d0*x1794 + 0.25d0 *x1906
      x1908 = x185*x188
      x1909 = x45*x537
      x1910 = x3*x560
      x1911 = 2.0d0*x185
      x1912 = x1911*x203
      x1913 = x198*x212
      x1914 = x1913*x709
      x1915 = x1907*x214 + x1912*x712*(x1690 + 6.0d0*x1908 + 4.0d0*x1909
     & - 4.0d0*x1910 ) + 4.0d0*x1914*x711
      x1916 = 0.13125d0*x534
      x1917 = x713*x723 + x719
      x1918 = 96.0d0*x3
      x1919 = 1.4d0*x564
      x1920 = x1690*x563
      x1921 = 0.0875d0*x538
      x1922 = 8.0d0*x1907*x718 + 32.0d0*x1914*x723 + 16.0d0*x713*(0.0437
     &5d0*x1674 - x1675*x1921 + 0.04375d0*x1682 + 0.04375d0*x1689 + x185
     &*(0.59325d0 *x182 + 0.228375d0*x187 + 0.13125d0*x47 + 0.139125d0) 
     &+ 0.175d0* x1909 - 0.175d0*x1910)
      x1923 = x231*x238
      x1924 = x235*x236
      x1925 = x1110*x1924
      x1926 = x1724*x185
      x1927 = x1745*x185
      x1928 = 542908.851247439d0*x185
      x1929 = x1727*x185
      x1930 = x1740*x644
      x1931 = x1737*x644
      x1932 = 2452922.678884d0*x185
      x1933 = x586*x642
      x1934 = x1730*x643
      x1935 = x1748*x644
      x1936 = -723878.468329919d0*x1697 - 22330.0072465503d0*x1700 + 169
     &9901.01409917d0*x1703 - 1656900.21113256d0*x1704 + 14318.455069069
     &3d0*x1706*x1933 - 82845.0105566281d0*x1717*x59* x642 + x1720*x1928
     & + x1722*x1928 + 19623381.431072d0*x1729 - 12582911.5155088d0*x173
     &3*x1933 - 86063138.566041d0*x1736 + 83886080.0d0*x1743 + 4194304.0
     &d0*x1747 + 21515784.6415103d0*x1752 - 33554432.0d0*x1753 + 1398101
     &.33333333d0*x1758 + 1048576.0d0* x1767*x707 + 4194304.0d0*x1785 + 
     &2097152.0d0*x1787 + 524288.0d0* x1794 + 72627544.2317758d0*x185*x1
     &930 + 165289.921585764d0*x185* x1934 - 3721.66787442505d0*x185*x72
     &8 + 524288.0d0*x1906 + 8388608.0d0*x1926 - 8388608.0d0*x1927 + 838
     &8608.0d0*x1929 + x1931 *x1932 - x1932*x1935
      x1937 = 3.0d0*x247
      x1938 = 2.68099692748912d0*rho_c**(-2.33333333333333d0)*x245
      x1939 = 4.02149539123368d0*x744
      x1940 = 0.673851049238073d0*rho_c**(-2.66666666666667d0)*x193*x250
     &
      x1941 = 0.211916079314314d0*x187*x346
      x1942 = 0.074356519057654d0*x347
      x1943 = x187*x255
      x1944 = x1943*x877
      x1945 = 0.133841734303777d0*x1944*x345
      x1946 = 0.111534778586481d0*x345*x878
      x1947 = 0.455740507457808d0*x15
      x1948 = 0.00418159828868294d0*x15 + 0.0259624262672375d0*x20*x263 
     &+ 0.0375734556469842d0*x20 + 0.00197049465647014d0*x257 + x262 - 0
     &.151913502485936d0*x871
      x1949 = x754*(x1947*x876 + x1948)
      x1950 = x1949*x879
      x1951 = 3.15778393965961d0*x350
      x1952 = 2.10518929310641d0*x349*x883
      x1953 = x1854*x53
      x1954 = x4*x52
      x1955 = x1855*x1954
      x1956 = x1859*x601
      x1957 = 0.58482233974552d0*x1738
      x1958 = 0.58482233974552d0*x1749
      x1959 = x1861*x53
      x1960 = x4*x926
      x1961 = x1862*x1960
      x1962 = x1863*x1960
      x1963 = x1777 - x1778 - x1779 + x1780 + x1781 + x1782
      x1964 = -x1869 + x1871 - x1875 + x1879 + x1953 + x1955 - x1956 - x
     &1957 + x1958 - x1959 - x1961 + x1962 + x1963
      x1965 = x1864 - x1865 + x1867 + x1872 + x1873 - x1877 + x1964
      x1966 = x1965*x884
      x1967 = 12.0d0*x352
      x1968 = 0.00464728244110338d0*x345
      x1969 = -x1769 + x1770 + x1772 - x1774 - x1775 - x1776 + x1963
      x1970 = x1969*x48
      x1971 = x652*x795
      x1972 = x1746*x799
      x1973 = -x1866*x1971 + x1876*x799 + x1964 - x1970 - x1972 + 8.0d0*
     &x796 - 20.0d0* x96
      x1974 = 4.13905085628136d0*x257
      x1975 = 0.299209d0*x363
      x1976 = 0.181712059283214d0*x894
      x1977 = 14.5369647426571d0*x902
      x1978 = 10.9027235569928d0*x301*x899*(-0.156676175559749d0*x15 + 0
     &.620028561181324d0*x20)
      x1979 = 3.61230841029072d0*x904
      x1980 = x1979*x361
      x1981 = x1980*x901
      x1982 = 2.22044604925031d-16*x905
      x1983 = 2.42282745710952d0*x896
      x1984 = x899*x900
      x1985 = 43.6108942279714d0*x187*x1984
      x1986 = 7.22461682058144d0*x187*x362*x904
      x1987 = 36.3424118566428d0*x187*x907
      x1988 = 5.0d0*x4
      x1989 = x1988*x272
      x1990 = x1988*x273
      x1991 = x1989 + x1990
      x1992 = 0.726848237132856d0*x321
      x1993 = x1992*x4*(x1991 - 6.0d0*x910 + 6.0d0*x911)
      x1994 = 2.90739294853142d0*x913
      x1995 = 0.00942485980224818d0*x818
      x1996 = x1995*x816
      x1997 = x1996*x814
      x1998 = x1995*x826
      x1999 = x1998*x824
      x2000 = 1d0/x304
      x2001 = 0.0115260005177981d0*x187
      x2002 = 0.817704266774463d0*x818
      x2003 = x2002*x814*(-x2000*x2001 + 0.0230520010355962d0*x815 - 0.0
     &115260005177981d0)
      x2004 = 0.181712059283214d0*x321
      x2005 = x2004*x317
      x2006 = x2004*x335
      x2007 = 1d0/x323
      x2008 = x2002*x824*(x2001*x2007 + 0.0230520010355962d0*x825 + 0.01
     &15260005177981d0)
      x2009 = 0.0141372897033723d0*x818
      x2010 = x812*x832
      x2011 = x2010*x840
      x2012 = x2009*x2011*x304*x816
      x2013 = x822*x843
      x2014 = x2009*x2013
      x2015 = x2014*x826*x851
      x2016 = 3.63424118566428d0*x321
      x2017 = x2016*x316
      x2018 = x187*x2017
      x2019 = x2016*x334
      x2020 = x187*x2019
      x2021 = 1.09027235569928d0*x321
      x2022 = x2021*x834*x841
      x2023 = x2021*x845*x851
      x2024 = 0.116477395298173d0*x20
      x2025 = 0.0582386976490866d0*x20
      x2026 = x2025*x310
      x2027 = x187*x2000
      x2028 = 0.219734258381682d0*x307
      x2029 = -x2026 - x2028
      x2030 = x304*x835*(x2024*x310*x815 - x2026*x2027 - x2027*x2028 + x
     &2029 + 0.439468516763364d0*x307*x815)
      x2031 = x2025*x329
      x2032 = x187*x2007
      x2033 = 0.219734258381682d0*x326
      x2034 = x2031 + x2033
      x2035 = x323*x846*(x2024*x329*x825 + x2031*x2032 + x2032*x2033 + x
     &2034 + 0.439468516763364d0*x326*x825)
      x2036 = 0.0376994392089927d0*x4*x818
      x2037 = x1222*x2036*x826
      x2038 = 1.45369647426571d0*x321*x4
      x2039 = x2038*x857
      x2040 = 0.0698864371789039d0*x20
      x2041 = x2040*x310
      x2042 = 0.329601387572523d0*x307
      x2043 = x2041 + x2042
      x2044 = x312**(-3)
      x2045 = x2044*x309*x314*x322
      x2046 = x2045*x841*(-x2041*x815 - x2042*x815 + x2043)
      x2047 = x2040*x329
      x2048 = 0.329601387572523d0*x326
      x2049 = x2047 + x2048
      x2050 = x330**(-3)
      x2051 = x2050*x322*x328*x332
      x2052 = x2051*x851*(x2047*x825 + x2048*x825 + x2049)
      x2053 = x1220*x816
      x2054 = x2036*x2053
      x2055 = x2038*x854
      x2056 = x6*x845
      x2057 = 4.36108942279713d0*x321
      x2058 = x2057*x4
      x2059 = x2056*x2058*x850
      x2060 = x5*x834
      x2061 = x2058*x2060*x840
      x2062 = x1997 + x1999 + x2003 - x2005 - x2006 - x2008 + x2012 + x2
     &015 + x2018 + x2020 + x2022 + x2023 + x2030 - x2035 + x2037 + x203
     &9 + x2046 + x2052 - x2054 - x2055 + x2059 - x2061
      x2063 = 0.263148661638301d0*x260
      x2064 = 1.3157433081915d0*x187*x354 - 1.0525946465532d0*x1944*x260
     & - x301*x886 + 0.526297323276602d0*x888
      x2065 = 12.0d0*x373
      x2066 = 0.0037178259528827d0*x267
      x2067 = x297*x804*(-1.99883265211535d0*x15 + 0.917202013581841d0*x
     &20)
      x2068 = x301*x803
      x2069 = 24.0d0*x299
      x2070 = x2068*x2069
      x2071 = 0.154912426780983d0*x3**(-0.333333333333333d0)*x806
      x2072 = 2.25322953511185d0*x807
      x2073 = x2068*x2072
      x2074 = 0.751076511703951d0*x808
      x2075 = 48.0d0*x187*x299*x803
      x2076 = 4.50645907022371d0*x187*x298*x807
      x2077 = 0.133252782879869d0*x864
      x2078 = 0.0799516697279215d0*x864
      x2079 = 0.00910930363347549d0*x267
      x2080 = 0.00836319657736588d0*x15 + 0.0751469112939684d0*x20 + 0.0
     &0394098931294027d0*x257 - 0.01915695d0*x260
      x2081 = x2080*x863
      x2082 = 0.350864882184401d0*x293
      x2083 = 0.701729764368802d0*x293
      x2084 = x1943*x753
      x2085 = 1.75432441092201d0*x187*x294 - 0.00495710127051027d0*x336 
     &- 0.159903339455843d0*x747*x864
      x2086 = 12.0d0*x338
      x2087 = 0.803050405822664d0*x341
      x2088 = 7.85204841248826d0*x341
      x2089 = 4.28293549772087d0*x341
      x2090 = 5.57009647607783d0*x768
      x2091 = x763**2
      x2092 = x187*x2091
      x2093 = x286*x770
      x2094 = 33.420578856467d0*x2092*x2093
      x2095 = 6.0d0*x761
      x2096 = 6.0d0*x762
      x2097 = x5**(-1.33333333333333d0)
      x2098 = x2097*x4
      x2099 = x6**(-1.33333333333333d0)
      x2100 = x2099*x4
      x2101 = x2098 + x2100
      x2102 = -x2095 + x2096 + x2101
      x2103 = x2102*x4
      x2104 = x2103*x771
      x2105 = 5.57009647607783d0*x2104
      x2106 = 11.1401929521557d0*x773
      x2107 = 2.16703006958462d0*x25
      x2108 = 2.476605793811d0*x25
      x2109 = x2102*x2108*x777
      x2110 = 9.90642317524398d0*x25
      x2111 = x2091*x2110
      x2112 = x187*x2111*x277
      x2113 = x2108*x778
      x2114 = x2107 + x2109 + x2112 - x2113
      x2115 = x281*x781
      x2116 = 2.0d0*x284
      x2117 = x2116*x789
      x2118 = 3.71490869071649d0*x25
      x2119 = 8.35514471411675d0*x768
      x2120 = 14.859634762866d0*x25
      x2121 = x2120*x778
      x2122 = 33.420578856467d0*x773
      x2123 = x288/x283**2
      x2124 = x2115*x2123
      x2125 = 12.0326097596534d0*x768
      x2126 = 72.1956585579204d0*x2091
      x2127 = x187*x2093*x2126 + 144.751042466189d0*x2092*x22*x250/x274*
     &*4 + 18.0938803082737d0*x2103*x786 + 12.0326097596534d0*x2104 + x2
     &114 + x2125 - 24.0652195193068d0*x773 + 20.3556153468079d0*x783 - 
     &54.281640924821d0*x787
      x2128 = 8.0d0*x290
      x2129 = x2091*x2128
      x2130 = x1212*x2129
      x2131 = 3.0d0*x250
      x2132 = x2131/x288**2
      x2133 = x2132*x790
      x2134 = x269*x757
      x2135 = 0.121994199821273d0*x2134
      x2136 = 0.284653132916304d0*x757
      x2137 = 0.267683468607554d0*x2080*x270
      x2138 = 0.892278228691848d0*x270
      x2139 = 0.000542507213303895d0*x3**2.33333333333333d0
      x2140 = 1.18970430492246d0*x270
      x2141 = x2131*x791
      x2142 = x2141*x284
      x2143 = x2142*x781
      x2144 = 1.07073387443022d0*x270
      x2145 = x1213*x790
      x2146 = -0.243988399642546d0*x2134*x747 + 3.39065726902902d0*x270*
     &x747
      x2147 = x253*(x1212*x2133*x789 - x1212*x2143*x790 - x187*x2130 - x
     &1967*(-x1941 + x1942 + x1945 - x1946 + x1950 - x1951 + x1952 - x19
     &66) + x1973* x2087 + x2065*(x1968*(x860 - x895 + x897 + x903 + x90
     &6 + x909 - x914) + x1973*x372 - x1974*x802 - x2063*x754*(-x1947*(x
     &872 - x873 + x874 + x875) + x1948) + x2064 + 4.82889266566159d0*x2
     &57*x97 - 0.00619637658813783d0*x370*x385 - x917*(x1975 + x1976 + x
     &1977 + x1978 + x1981 - x1982 - x1983 + x1985 + x1986 - x1987 - x19
     &93 - x1994 + x2062)) + x2084*x2144 + x2086*(x2066*x861 + x2077*x38
     &3 + x2078*x755 - x2079*x383 - x2081*x754 - x2082*x755 + x2083*x208
     &4 + x2085 - x301*x866 - x862*(x2062 + x2067 + x2070 + x2071 + x207
     &3 + x2074 + x2075 + x2076 + x995)) + x2088*x97 - x2089*x802 + x210
     &2* x760*x989 + x2135*x755 + x2136*x756 - x2137*x754 - x2138*x755 -
     & x2139*x756 - x2140*x383 + 12.0d0*x2145*x772 + x2146 - 3.0d0*x793*
     & (-x2115*x2117 + x2124*(-x2118 - x2119 + x2121 + x2122) + x2127 - 
     &x289*(x2090 + x2094 + x2105 - x2106 + x2114)))
      x2148 = x1701*x798
      x2149 = x1796*x939
      x2150 = x2149*x623
      x2151 = x1243*x2150
      x2152 = x939**2
      x2153 = 0.0308430206873641d0*x2152
      x2154 = x1243*x624
      x2155 = x2154*x969
      x2156 = x968**2
      x2157 = x946**(-2)
      x2158 = x2154*x2157
      x2159 = 0.0893877004050842d0*x2158
      x2160 = 0.032153849066577d0*x382
      x2161 = x2155*x2160
      x2162 = 0.032153849066577d0*x2158
      x2163 = x2162*x968
      x2164 = 0.0134938215507218d0*x619
      x2165 = x1020*x393
      x2166 = x1021*x2165
      x2167 = x1016*x662
      x2168 = x2165*x2167
      x2169 = x684*x939
      x2170 = 0.00771075517184102d0*x2169
      x2171 = x2166*x619
      x2172 = 8.34d0*x948
      x2173 = x382*x967
      x2174 = 0.258879113792152d0*x738
      x2175 = 0.00133097214976253d0*x928
      x2176 = 0.158014317620522d0*x1868
      x2177 = 0.101321996098707d0*x1960*x59
      x2178 = 4.0d0*x577
      x2179 = 1.16964467949104d0*x1058
      x2180 = x1866*x569*x652
      x2181 = 8.0d0*x50
      x2182 = 2.0d0*x648
      x2183 = x946/x97**2
      x2184 = 8.34d0*x2183
      x2185 = x1971*x566
      x2186 = x257*x373/x252**6
      x2187 = x2186*x386
      x2188 = x920*x930
      x2189 = x652*x799
      x2190 = x1716*x1954
      x2191 = x2190*x59
      x2192 = x1960*x592
      x2193 = x1728*x1960
      x2194 = 4.3181036180531d0*x1720 + 4.3181036180531d0*x1722 + 66.72d
     &0*x1724 + 66.72d0*x1727 - 1.3146571598294d0*x1731 - 19.50967325391
     &06d0* x1738 - 577.653616803179d0*x1741 - 66.72d0*x1745 + 19.509673
     &2539106d0*x1749 - 13.1783940895515d0*x1868 + 13.5204071594115d0*x1
     &870 - 0.177604923664312d0*x1874 + 156.077386031284d0*x1878 + 0.658
     &919704477575d0*x2191 - 0.113883891368902d0*x2192 + 100.07999614652
     &9d0*x2193 - 5.75747149073746d0*x738 + 0.0296008206107187d0*x928
      x2195 = -x1965*x2172 - 33.36d0*x1970 - 33.36d0*x1972 - 2.78d0*x214
     &7 - 16.68d0* x2173*x883 + x2184*x883*(x2174 - x2175 + x2176 + x217
     &7 - x2178 + x2179 + x2180 - x2181 + x2182) - 171.128887090869d0*x2
     &185 - 19.1776023007703d0*x2187 - 22.6232590535263d0*x2188 + 684.51
     &5548363478d0*x2189 + x2194 + 266.88d0*x796 + 18.8527158779386d0*x9
     &31 - 667.2d0*x96
      x2196 = 0.00385537758592051d0*x2154
      x2197 = x2196*x947
      x2198 = x1052*x1243
      x2199 = x1051*x395
      x2200 = x2199*x662
      x2201 = x1243*x969
      x2202 = x1052*x619
      x2203 = 0.00144576659472019d0*x78
      x2204 = x634*x939
      x2205 = 0.123372082749456d0*rho_s
      x2206 = x1047*x1820
      x2207 = x1036*x2206
      x2208 = 0.00428081366304526d0*x1029
      x2209 = 0.0028538757753635d0*x939
      x2210 = x1037*x2208 + x1041*x2209 + 0.00499428260688613d0*x1042 - 
     &0.0059503309916329d0*x1043
      x2211 = 0.000722883297360095d0*x78
      x2212 = x2211*x82
      x2213 = x398**(-2)
      x2214 = x1243*x2213
      x2215 = x1050*x2214
      x2216 = x2212*x2215
      x2217 = 4.31802289623097d0*x633
      x2218 = 3.45441831698478d0*x1828
      x2219 = 1.48046499299348d0*x2152
      x2220 = x1830*x2219
      x2221 = 0.555174372372553d0*x1833
      x2222 = 7.20246219091326d0*x633
      x2223 = 4.29060961944404d0*x2157
      x2224 = x2223*x633
      x2225 = x1828*x2149
      x2226 = 5.18162747547716d0*x624
      x2227 = 4.44139497898043d0*x381*exp(-x1837*x378)
      x2228 = x1036**2*x1842
      x2229 = 2.22069748949021d0*x1047
      x2230 = x1849*x2149
      x2231 = 1.5433847551957d0*x382
      x2232 = x633*x969
      x2233 = 1.5433847551957d0*x2157
      x2234 = x633*x968
      x2235 = x2207*x624
      x2236 = x1032*x640
      x2237 = 6.0d0*x1031
      x2238 = 16.0d0*x1124
      x2239 = x1032*x2238*x980
      x2240 = 0.237021476430782d0*x59
      x2241 = -1.55327468275291d0*x1720 - 1.55327468275291d0*x1722 - 24.
     &0d0*x1724 - 24.0d0*x1727 + 24.0d0*x1745 + 4.74042952861565d0*x1868
     & - 4.86345581273796d0*x1870 + 0.0638866631886015d0*x1874 - 56.1429
     &4461557d0*x1878 + 2.07103291033722d0*x738
      x2242 = 0.370116248248369d0*x2206
      x2243 = x2242*x624
      x2244 = 0.185058124124184d0*x947
      x2245 = 0.000240961099120032d0*x1025
      x2246 = x2245*x82
      x2247 = x1797*x974
      x2248 = 1.48046499299348d0*x974
      x2249 = x1294*x687
      x2250 = x1270*x685
      x2251 = x1892*x2250
      x2252 = x1005*x974
      x2253 = x980**(-2)
      x2254 = x1004**2*x2253
      x2255 = x633*x974
      x2256 = 4.29060961944404d0*x2255
      x2257 = x1797*x1828
      x2258 = x1017*x974
      x2259 = x1820*x2258
      x2260 = x1012*x1303
      x2261 = x1028*x1894
      x2262 = x2260*x2261
      x2263 = x1012**2
      x2264 = x1842*x2263
      x2265 = exp(x1013*x1837)
      x2266 = 4.44139497898043d0*x1301*x2265*x974
      x2267 = 2.22069748949021d0*x2258
      x2268 = 0.154109291869629d0*x1269*x1899*x2265/x1015**3
      x2269 = x1303*x1900
      x2270 = x1012*x1017
      x2271 = x1903*x2260
      x2272 = x1366*x685
      x2273 = 1.5433847551957d0*x2255
      x2274 = x1005*x972
      x2275 = x2273*x2274
      x2276 = 0.198503041880874d0*x2249
      x2277 = 0.0357019859497974d0*x2249
      x2278 = x2274*x2277
      x2279 = x1004*x2253
      x2280 = x2279*(x1010 + x1034)
      x2281 = x1820*x624
      x2282 = 6.17353902078279d0*x2281
      x2283 = -x1997 - x1999 - x2003 + x2005 + x2006 + x2008 - x2012 - x
     &2015 - x2018 - x2020 - x2022 - x2023 - x2030 + x2035 - x2037 - x20
     &39 - x2046 - x2052 + x2054 + x2055 - x2059 + x2061
      x2284 = x254*x255
      x2285 = x187*x984
      x2286 = x281*x992
      x2287 = x2123*x2286
      x2288 = -x2098 - x2100
      x2289 = x2132*x993
      x2290 = x987**2
      x2291 = x1276*x2128
      x2292 = x1277*x993
      x2293 = x253*(x1276*x2143*x993 + x1276*x2289*(x782 + x784 - x785 -
     & x788 + x991) + x187*x2290*x2291 + x1965*x2087 + x1967*(x1941 - x1
     &942 - x1945 + x1946 - x1950 + x1951 - x1952 + x1966) + x2065*(-x19
     &49*x2063 + x1965*x372 - x1968*x916 - x1974*x883 + x2064 + 0.006196
     &37658813783d0*x890 + 4.82889266566159d0*x892 + x917*( -x1975 - x19
     &76 - x1977 - x1978 - x1981 + x1982 + x1983 - x1985 - x1986 + x1987
     & + x1993 + x1994 + x2283)) + x2086*(-x2066*x996 - x2077*x265 + x20
     &78*x985 + x2079*x265 + x2081*x2284 - x2082*x985 - x2083*x2285 + x2
     &085 + x254*x866 + x862*(-x2067 - x2070 - x2071 - x2073 - x2074 - x
     &2075 - x2076 + x2283 + x811)) + x2088*x63 - x2089*x883 + x2135*x98
     &5 - x2136*x986 + x2137*x2284 - x2138*x985 + x2139*x986 + x2140*x26
     &5 - x2144*x2285 + x2146 + 12.0d0*x2292*x4* x987 + x990*(x2095 - x2
     &096 + x2288) + 3.0d0*x994*(x2117*x2286 + x2127 + x2287*(x2118 + x2
     &119 - x2121 - x2122) + x289*(-x2090 - x2094 - x2105 + x2106 - x210
     &7 - x2109 - x2112 + x2113)))
      x2294 = x1751*x566
      x2295 = 6.0d0*x1007
      x2296 = x2186*x375
      x2297 = x930*x997
      x2298 = x375*x930
      x2299 = x1007*x639
      x2300 = x650*x652
      x2301 = 12.0d0*x1011*x641 + x1716*x2240*x645 + 35.9999986138594d0*
     &x1728*x646 + x1853*x2295 + 16.0d0*x1858*x2299 - 12.0d0*x1864 - 12.
     &0d0*x1872 - 207.78907079251d0*x1930 - 7.01786807694624d0*x1931 - 0
     &.472898258931437d0*x1934 + 7.01786807694624d0*x1935 + x2241 + x229
     &3 - 61.5571536298092d0*x2294 + 6.89841809380227d0*x2296 + 8.137862
     &96889434d0*x2297 - 6.78155247407862d0*x2298 + 246.228614519237d0*x
     &2300 - 240.0d0*x50 - 0.0409654285499647d0* x592*x646 + 96.0d0*x648
     & + 0.0106477771981002d0*x728
      x2302 = 0.370116248248369d0*x2258*x2281
      x2303 = x1303*x2261
      x2304 = 0.00856162732609051d0*x2303
      x2305 = x983*(x1880 - x1953 - x1955 + x1956 + x1957 - x1958 + x195
     &9 + x1961 - x1962)
      x2306 = x1003*x972
      x2307 = 16.68d0*x982
      x2308 = x980/x63**2
      x2309 = x2308*(-x2174 + x2175 - x2176 - x2177 + x2178 - x2179 - x2
     &180 + x2181 - x2182)
      x2310 = 0.185058124124184d0*x2255*x981
      x2311 = 0.00214040683152263d0*x2249*x981
      x2312 = x1021*x1314
      x2313 = 0.000240961099120032d0*x2312*x82
      x2314 = x1317*x936
      x2315 = 1628726.55374232d0*x185
      x2316 = x84*x974
      x2317 = x1021*x2316 + 1.0d0
      x2318 = x1302*x79 + x2317
      x2319 = 1d0/x2318
      x2320 = x2317*x974
      x2321 = x2320*x970
      x2322 = x2319*x2320
      x2323 = x1023*x2322
      x2324 = 1572864.0d0*x1054
      x2325 = x402**(-2)
      x2326 = x2325*x399
      x2327 = x1790*x1793*x2326
      x2328 = 6291456.0d0*x1054*x1317*x1786 - 3145728.0d0*x1055*x1705*(x
     &1759 - x1760 + x1851) - 2171635.40498976d0*x1697 + 100663296.0d0*x
     &1698*x923 - 12582912.0d0*x1699*x1969 - 251658240.0d0*x1699*x798 - 
     &66990.021739651d0*x1700 - 64547353.9245308d0*x1701*x1971 + 5099703
     &.04229752d0*x1703 - 4970700.63339768d0*x1704 - 42955.3652072078d0*
     &x1707 + 11165.0036232752d0*x1708 + 248535.031669884d0*x1718*x59 + 
     &x1720*x2315 + x1722*x2315 + 58870144.2932159d0*x1729 - 495869.7647
     &57291d0*x1732 + 37748734.5465263d0*x1734 - 7358768.03665198d0*x173
     &9 - 217882632.695327d0*x1742 - 12582912.0d0*x1746*x2148 + 7358768.
     &03665198d0*x1750 + 4194304.0d0*x1757*x403 - 1048576.0d0* x185*x214
     &7 - 7233515.65112681d0*x185*x2187 - 8533167.80047135d0* x185*x2188
     & + 7110973.16705946d0*x185*x931 + 25165824.0d0*x1926 - 25165824.0d
     &0*x1927 + 25165824.0d0*x1929 + 258189415.698123d0* x2148*x652 + 15
     &72864.0d0*x2314*(-x1003*x2163 - x1016*x2165*x2199* x2203*x82 + x12
     &43*x1795 + x1243*x1802*x2153 + 0.0115661327577615d0*x1243*x1810 - 
     &0.0857435975108721d0*x2151* x969 + 0.0719670482705162d0*x2151 - x2
     &154*x2170*x2200 - 0.150051295644026d0*x2155 + x2156*x2159 + x2161*
     &x802 + x2164* x2166 - x2164*x2198 + x2168*x2170*x624 - 0.016076924
     &5332885d0* x2171*x969 + x2195*x2197 + 0.0160769245332885d0*x2201*x
     &2202 + x2216*x381*(0.215901144811548d0*x1044 - 0.257230792532616d0
     &*x1045 + x1822*x2207*x393 + x2204*x2205*x393 + x2210) + x2246*(-x1
     &003* x2233*x2234 + 2.96092998598695d0*x1048*x2230 + x1801*x2220 + 
     &x2149*x2218 + x2156*x2224 + x2195*x2244*x633 + x2207*x2226 + x2217
     & + x2221 - x2222*x969 - 4.11569268052186d0*x2225*x969 + x2227*x222
     &8 - x2228*x2229 + x2231*x2232*x802 - 6.17353902078279d0 *x2235*x96
     &9 - x2243*(-12.0d0*x1035*x2236 + 0.472898258931437d0* x1731 + 7.01
     &786807694624d0*x1738 + 207.78907079251d0*x1741 - 7.01786807694624d
     &0*x1749 - x1852*x2237*x4 - x187*x2239 + 12.0d0* x1970 + 12.0d0*x19
     &72 + x2147 + 61.5571536298092d0*x2185 + 6.89841809380227d0*x2187 +
     & 8.13786296889434d0*x2188 - 246.228614519237d0*x2189 - x2190*x2240
     & + 0.0409654285499647d0* x2192 - 35.9999986138594d0*x2193 + x2241 
     &- 96.0d0*x796 - 0.0106477771981002d0*x928 - 6.78155247407862d0*x93
     &1 + 240.0d0*x96 )) - x2313*(-0.333218535531443d0*x1005*x2249 + 0.1
     &90410591732253d0*x1005*x2251 - 0.285615887598379d0*x1005*x2262 + x
     &1012*x2226*x2259 + 0.0570775155072701d0*x1270*x1887 + 0.0064212204
     &9456788d0*x1270*x1893*x684 + x1831*x2248 - 2.96092998598695d0*x184
     &9*x2247*x2270 + x2217*x974 - x2218*x2247 + x2221*x974 - x2222*x225
     &2 + 0.169805608634128d0*x2249 - 0.159817043420356d0*x2251 + 4.1156
     &9268052186d0*x2252*x2257 - x2252*x2270*x2282 + x2254*x2256 + x2254
     &*x2276 + 0.239725565130534d0*x2262 + x2263*x2268 - x2263*x2269 + x
     &2264* x2266 - x2264*x2267 - 0.136986037217448d0*x2271*x2272 + x227
     &3* x2280 + x2275*x883 + x2277*x2280 + x2278*x883 - x2301*x2302 - x
     &2301*x2304 + x2310*(33.36d0*x1864 + 33.36d0*x1872 + x2194 - 2.78d0
     &*x2293 + 171.128887090869d0*x2294 - 19.1776023007703d0* x2296 - 22
     &.6232590535263d0*x2297 + 18.8527158779386d0*x2298 - 684.5155483634
     &78d0*x2300 + 8.34d0*x2305 + x2306*x2307 + 8.34d0* x2309*x982 + 667
     &.2d0*x50 - 266.88d0*x648) + x2311*( 8.63620723610619d0*x1720 + 8.6
     &3620723610619d0*x1722 + 133.44d0* x1724 + 133.44d0*x1727 - 2.62931
     &431965879d0*x1731 - 39.0193465078211d0*x1738 - 1155.30723360636d0*
     &x1741 - 133.44d0* x1745 + 39.0193465078211d0*x1749 + 66.72d0*x1864
     & - 26.356788179103d0*x1868 + 27.040814318823d0*x1870 + 66.72d0*x18
     &72 - 0.355209847328624d0*x1874 + 312.154772062569d0*x1878 + 1.3178
     &3940895515d0*x2191 - 0.227767782737804d0*x2192 + 200.159992293059d
     &0*x2193 - 5.56d0*x2293 + 342.257774181739d0* x2294 - 38.3552046015
     &406d0*x2296 - 45.2465181070525d0*x2297 + 37.7054317558771d0*x2298 
     &- 1369.03109672696d0*x2300 + 16.68d0* x2305 + 33.36d0*x2306*x982 +
     & x2307*x2309 + 1334.4d0*x50 - 533.76d0*x648 - 11.5149429814749d0*x
     &738 + 0.0592016412214373d0* x928))) - x2319*x2324*x2327*(x1005*x23
     &21 + x1023*x1269*(x1006 - x1019 - x976 + x978) - x2320*x620 + x232
     &0*x629*x941 + x2323*( -x1005*x1295 - x1006*x974 + x1012*x1304 + x1
     &019*x974 + x1027* x1294 - 0.00142693788768175d0*x2250*x688 + x974*
     &x976 - x974*x978 ))
      x2329 = 7.94728597005208d-8*rho_c
      x2330 = 3.19599739013806d-7*x744
      x2331 = x568*x798
      x2332 = 5.1297628024841d0*x652
      x2333 = x1058*x3
      x2334 = 0.158014317620522d0*x59
      x2335 = rho_s*x927
      x2336 = x2174*x3 - x2175*x3 + x2334*x568 + 0.101321996098707d0*x23
     &35
      x2337 = 0.25d0*x1054*x2314 + 0.5d0*x1055*x704 - x2178*x3 + 8.0d0*x
     &2331 - x2332* x929 + 1.16964467949104d0*x2333 + x2336 - 0.5d0*x922
     & - 2.0d0*x924 - 0.678155247407862d0*x932
      x2338 = x1073*x493
      x2339 = x1072*x1595
      x2340 = x1601*x164
      x2341 = x2340*x495
      x2342 = x1073*x2341
      x2343 = x1072*x71
      x2344 = x2341*x492
      x2345 = x1067*x497
      x2346 = x1069*x1603
      x2347 = -x1593 - x1597 - x1604 - x1605 + x2345 + x2346
      x2348 = 0.00192768879296025d0*x1621
      x2349 = x1067*x1615
      x2350 = x1067*x1587
      x2351 = 0.00192768879296025d0*x2350
      x2352 = x1615*x1619
      x2353 = x1066*x483
      x2354 = x1066*x1620
      x2355 = 0.313621176746478d0*x1569*x1625 - x1624 + x1629 + x1631 - 
     &x1633
      x2356 = x1066*x501
      x2357 = x1639*x2356
      x2358 = x1643*x2356
      x2359 = x1615*x1638
      x2360 = 0.25d0*x473
      x2361 = x1077*x2360 + x1080 - x2360*x518 + x526*(-x1077*x1654 + x1
     &078*(x1074*( -x1588 - 0.0115661327577615d0*x1620 - x1635*x2355 - x
     &2348 + 0.0115661327577615d0*x2349 + x2351 + x2352*x2353 - 0.005783
     &06637888076d0*x2354) - x1590 + x1596*x498 + x1596*x500 + x1599*x23
     &42 + x1614*x2347 - x2338*x488 - x2338*x490 - x2339*x498 - x2339*x5
     &00 - x2343*x2344 - x488*x494 - x490*x494 + x515*(x130* x2354 + x14
     &69*x1620 - x1469*x2349 + x1552*x1642 - x1552*x2357 + x1553*x1644 -
     & x1553*x2358 + x1589 + x1640*x2353 + x1640*x2355 - x1645*x2353 + x
     &1646*x2355 + x1647 + x1648 + x2348*x81 - x2351*x81 - x2353*x2359))
     & - x1657 + x1658 - x1660 - x1661 + x1662 - x1663)
      x2362 = x204*(x1586 + x2361)
      x2363 = x1083*x534
      x2364 = x1094*x531
      x2365 = x1083*x564
      x2366 = x1705*(8.0d0*x187 + 81.36d0*x47 + 3.48d0)
      x2367 = x1091*x1381
      x2368 = x1086*x1675
      x2369 = x1087*x1676 + 0.400624880640947d0*x1090
      x2370 = x1681*x2369
      x2371 = x1088*x544
      x2372 = 0.640999809025515d0*rho_s
      x2373 = x2372*x552
      x2374 = x1089*x548
      x2375 = x559*(gamma_cs*x542 - x1087*x2373*x548 + 2.2434993315893d0
     &*x1087*x545 - 0.640999809025515d0*x2371*x541 + 1.49566622105954d0*
     &x2374 - 0.240374928384568d0*x549*(-x1686 - x1688 + 1.3333333333333
     &3d0* x550 - 1.33333333333333d0*x551))
      x2376 = -x2366 + x2367*x538 + x2368 - x2370 + x2375
      x2377 = x2376*x563
      x2378 = x210*x2377
      x2379 = x1093*x198
      x2380 = x1695*x2379
      x2381 = x1097*x566
      x2382 = x1698*x963
      x2383 = 0.00170689285624853d0*x592
      x2384 = x1285*x2383
      x2385 = 2.56488140124205d0*x585
      x2386 = 2.56488140124205d0*x584
      x2387 = -x1712 - x1714 - x2385 + x2386
      x2388 = x2387*x60
      x2389 = x2388*x52
      x2390 = x54*x961
      x2391 = 1.49999994224414d0*x1728
      x2392 = x1285*x2391
      x2393 = x1250*x569
      x2394 = x1097*x650
      x2395 = x2387*x569
      x2396 = x1698*x2395
      x2397 = x1106*x616
      x2398 = 2.56488140124205d0*x2397
      x2399 = x1096*x617
      x2400 = x2399*x566
      x2401 = x227*x628
      x2402 = x1705*x1755*x2401
      x2403 = x185*(x1851 + x703)
      x2404 = x2403*x707
      x2405 = 0.25d0*x696
      x2406 = x1104*x1137*x185*x684*x697
      x2407 = x695*x699
      x2408 = x10*x1792
      x2409 = 0.0359835241352581d0*x939
      x2410 = 1.33333333333333d0*x626
      x2411 = 1.33333333333333d0*x627
      x2412 = -x1806 - x1807 + x2410 - x2411
      x2413 = 0.0115661327577615d0*x2412
      x2414 = x1797*x939
      x2415 = x1121*x548
      x2416 = rho_c**(-7.66666666666667d0)
      x2417 = x2416*x662
      x2418 = x1888*x662
      x2419 = 0.0115661327577615d0*x219
      x2420 = x939*x975
      x2421 = 0.0308430206873641d0*x2420
      x2422 = 0.000713468943840876d0*x1127
      x2423 = x10*x660
      x2424 = x1119*x1822*x2423
      x2425 = x219*x77
      x2426 = x1132*x1823*x2425
      x2427 = -x2421 - x2422 + x2424 + x2426
      x2428 = x2409*x633
      x2429 = x2413*x633
      x2430 = 0.0308430206873641d0*x1828*x2414
      x2431 = x1119*x624
      x2432 = 0.215901144811548d0*x2423*x2431
      x2433 = x1119*x1820
      x2434 = x1797*x2205*x2433
      x2435 = 0.123372082749456d0*x1835*x939
      x2436 = x1840*x680
      x2437 = 1.48046499299348d0*x2436
      x2438 = x1838*x2437
      x2439 = x1118*x679
      x2440 = x2438*x2439*x624
      x2441 = 0.740232496496738d0*x2436
      x2442 = x2431*x2441*x657
      x2443 = x622*(x1766 + x639)
      x2444 = x1110*x2399
      x2445 = x2398*x47
      x2446 = x47*x963
      x2447 = x2446*x586
      x2448 = x1106*x2446
      x2449 = x2395*x47
      x2450 = x1252*x959
      x2451 = x2393*x47
      x2452 = x1252*x2390
      x2453 = x1290 - x2444 - x2445 + x2447 - x2448 - x2449 - x2450 - x2
     &451 + x2452
      x2454 = x1008*x1117 + x1105*x656 - x1107*x2383 + x1107*x2391 + x11
     &24*x1857*x4 + x2388*x642 + x2443*x637 + x2453
      x2455 = x1103*x658
      x2456 = x2454*x2455*x624
      x2457 = -x1803 + x1804 - x1890 - x1891
      x2458 = x1118*x219
      x2459 = x1897/x7**7
      x2460 = 0.0513697639565431d0*gamma_cc*x1901*x2459
      x2461 = x1902*x2459
      x2462 = x1131*x548
      x2463 = x699*(x1100*x1811 + x1101*x693 + x1121*x1888*x2417*x629 + 
     &x1122*x2419* x692*x78 + x1135*x1818*x2419 + x1135*x2418*x630 - x11
     &35*x621*x668 + x1192*(-x2428 + x2429 + x2430 + x2432 - x2434 - x24
     &35 + x2440 - x2442 - x2456) + x1800*x2414*x623 + x1826*x2427 + x19
     &05*( -0.00998856521377226d0*x1118*x2425*x548*x677 - x1118*x2460*x6
     &57 + 0.00166476086896204d0*x1126*x687 - 0.000133775426970164d0*x18
     &83* x2457 - 0.00118911490640146d0*x1885*x939 + 0.00570775155072701
     &d0* x1896*x939 + 0.00570775155072701d0*x1904*x2458*x675*x681 + x24
     &28 - x2429 - x2430 - x2432 + x2434 + x2435 + x2439*x2461 - x2440 +
     & x2442 + x2454*x2462 + x2456) - x2409*x625 + x2413*x625 - 0.026987
     &6431014436d0*x2415*x668)
      x2464 = x1106*x2382 - x1138*x2408 + x1144*x2407 + x1246*x566 - x12
     &88*x1698 + x1698*x2393 + x1698*x2398 - x185*x2384 + x185*x2389 + x
     &185*x2392 - x2381*x2390 + x2381*x959 - x2382*x586 - 16.0d0*x2394 +
     & x2396 + 4.0d0*x2400 + 0.333333333333333d0*x2402 + 0.25d0*x2404 + 
     &x2405* x2406 + 0.125d0*x2463
      x2465 = x1093*x1914
      x2466 = x1149*x1913
      x2467 = x1085*x180
      x2468 = x214*x2464 + x2465 + x2466*x711 + x713*(x2367*x710 + x2376
     & - 2.0d0*x2467 )
      x2469 = x1150 + x1151
      x2470 = x1156 + x1157
      x2471 = -x1154 - x1155 + x2470
      x2472 = 3.0d0*x3
      x2473 = x2464*x718 + 0.04375d0*x2465 + x2466*x723 + x713*(0.0875d0
     &*x1092*x3 - 0.04375d0*x2366 + x2367*x720 + 0.04375d0*x2368 - 0.043
     &75d0*x2370 + 0.04375d0*x2375 - 0.0875d0*x2467)
      x2474 = 64.0d0*x238
      x2475 = x1097*x59
      x2476 = x2475*x566
      x2477 = x1698*x59
      x2478 = 424975.253524793d0*x2477
      x2479 = x2381*x592
      x2480 = 21515784.6415103d0*x1698
      x2481 = x569*x586
      x2482 = x185*x592
      x2483 = x3*x639
      x2484 = x1728*x2381
      x2485 = x1728*x185
      x2486 = 262144.0d0*x1137
      x2487 = x1106*x2478 + 3579.61376726732d0*x1107*x2482 - 3145727.878
     &87719d0*x1107* x2485 - x1114*x2480 + 5378946.16037756d0*x1698*x239
     &7 - 41422.505278314d0*x2387*x732 - 33554432.0d0*x2394 + 2097152.0d
     &0* x2396 + 8388608.0d0*x2400 + 699050.666666667d0*x2402 + 524288.0
     &d0 *x2404 + 524288.0d0*x2406*x696 + 524288.0d0*x2407*x2483 - x2408
     &* x2486 + 262144.0d0*x2463 + 662760.084453024d0*x2476 - x2478*x586
     & + 5582.50181163758d0*x2479 + x2480*x2481 - 4905845.35776799d0* x2
     &484
      x2488 = x246*x247*x743
      x2489 = 3.0d0*x910
      x2490 = 3.0d0*x911
      x2491 = 21.8054471139857d0*x1984*x4
      x2492 = x4*(-x359 + x360)
      x2493 = 0.726848237132856d0*x321
      x2494 = x327 - 1.0d0
      x2495 = x2494*x6
      x2496 = x308 - 1.0d0
      x2497 = x2496*x5
      x2498 = x2497*x833
      x2499 = x840*x853
      x2500 = x2031*x825 + x2033*x825 + x2034
      x2501 = x2495*x844
      x2502 = x2016*x4
      x2503 = x850*x853
      x2504 = x2026*x815
      x2505 = x2028*x815
      x2506 = x2026 + x2028
      x2507 = x4*x853
      x2508 = 0.00471242990112409d0*x818
      x2509 = x1222*x2508
      x2510 = x1220*x2508
      x2511 = x1222*x1998
      x2512 = 0.0188497196044964d0*x818
      x2513 = x2512*x4
      x2514 = x2513*x813
      x2515 = x2513*x823
      x2516 = x2013*x6
      x2517 = x2516*x819*x850
      x2518 = x2011*x5*x819
      x2519 = x2516*x827*x849
      x2520 = x2010*x839
      x2521 = x2520*x5*x820
      x2522 = x2509 - x2510 + x2511 + x2514 + x2515 + x2517 - x2518 + x2
     &519 - x2521
      x2523 = 0.00232364122055169d0*x345
      x2524 = x1096*x923
      x2525 = x586*x795
      x2526 = x1715 + x2385 - x2386
      x2527 = x2526*x52
      x2528 = x47*x59
      x2529 = x2528*x586
      x2530 = x1249*x2526
      x2531 = x1249*x586
      x2532 = x2450 - x2452
      x2533 = -x1110*x2524 + x1113*x2525 + x1262 + x2383*x926 - x2391*x9
     &26 - x2527*x60 - 0.40528798439483d0*x2529 + x2530 - 20.51905120993
     &64d0*x2531 + x2532
      x2534 = 2.06952542814068d0*x257
      x2535 = x255*x4
      x2536 = x2535*x877
      x2537 = 0.526297323276602d0*x2536*x260 - 0.789445984914903d0*x354*
     &x4
      x2538 = 1.0525946465532d0*x349
      x2539 = 0.0669208671518886d0*x2536*x345 - 0.122688256445129d0*x346
     &*x4
      x2540 = 0.00185891297644135d0*x267
      x2541 = x2072*x298*x4
      x2542 = x2069*x4*x803
      x2543 = x1220*x1995*(x815 - 1.0d0)
      x2544 = x2053*x2512
      x2545 = x2017*x4
      x2546 = x2493*x857
      x2547 = x2019*x4
      x2548 = x2493*x854
      x2549 = x1224*x2507
      x2550 = x1225*x829
      x2551 = x1227*x2507
      x2552 = x1228*x829
      x2553 = x5*x835*(x2029 + x2504 + x2505)
      x2554 = x2060*x2499
      x2555 = x2056*x2503
      x2556 = x2500*x6*x846
      x2557 = x2043*x2045
      x2558 = x2557*x5*x840
      x2559 = x2049*x2051
      x2560 = x2559*x6*x850
      x2561 = x2535*x753
      x2562 = x1205*x2078 - 1.0525946465532d0*x294*x4
      x2563 = 2.14146774886044d0*x341
      x2564 = 3.0d0*x761
      x2565 = 3.0d0*x762
      x2566 = x2101 - x2564 + x2565
      x2567 = x2566*x764
      x2568 = 1.2383028969055d0*x25
      x2569 = 27.1408204624105d0*x783
      x2570 = x2111*x280
      x2571 = x2126*x771
      x2572 = 144.751042466189d0*x786
      x2573 = x2091*x2572
      x2574 = 12.0326097596534d0*x768
      x2575 = 18.0938803082737d0*x783
      x2576 = x284*x987
      x2577 = x1209*x2576
      x2578 = -x2577*x281*x789
      x2579 = x2120 + 33.420578856467d0*x768
      x2580 = x2579*x987
      x2581 = x2288 + x2564 - x2565
      x2582 = x2581*x768
      x2583 = 5.57009647607783d0*x2582
      x2584 = x2090*x987
      x2585 = x2122*x987
      x2586 = x2108*x2581
      x2587 = x2568*x987
      x2588 = x2110*x778*x987
      x2589 = x2586 - x2587 + x2588
      x2590 = 3.0d0*x1213
      x2591 = 6.0d0*x763
      x2592 = x1207*x987
      x2593 = x2592*x41
      x2594 = x41*x763
      x2595 = x1211*x2594
      x2596 = 6.0d0*x4*x792
      x2597 = x1209*x2142
      x2598 = 0.535366937215109d0*x270
      x2599 = x1205*x2135 - 1.96301210312207d0*x1205*x270
      x2600 = x253*(-x1215*x2563 + x1967*(-x1215*x2538 + x2533*x884 + x2
     &539) + x2065*( -x1215*x2534 + x1236*x2523 + x2533*x372 + x2537 + x
     &917*(-x1220* x1996 - x1980*x2492 - x1992*(x1991 - x2489 + x2490) +
     & x2043*x2044 *x2497*x314*x322*x840 - x2049*x2050*x2495*x322*x332*x
     &850 + x2491 + 29.0739294853142d0*x2492*x899 + x2493*x2495*x333 - x
     &2493*x2497* x315 - x2494*x2502*x333 - x2494*x2507*x844*x849 - x249
     &6*x2502* x315 - x2496*x2507*x833*x839 + x2498*x2499 - x2498*x322*(
     &-x2504 - x2505 + x2506) + x2498*x829*x839 + x2500*x2501*x322 - x25
     &01*x2503 - x2501*x829*x849 + x2522 - 1.45369647426571d0*x912)) + x
     &2086*( x1231*x2540 - x2082*x2561 + x2562 - x862*(x1218 - x2509 + x
     &2510 - x2511 - x2514 - x2515 - x2517 + x2518 - x2519 + x2521 - x25
     &41 - x2542 + x2543 + x2544 - x2545 - x2546 - x2547 + x2548 - x2549
     & + x2550 - x2551 - x2552 + x2553 + x2554 - x2555 + x2556 + x2558 -
     & x2560)) + x2087*x2533 + x2130*x4 + x2133*x2593 - x2145*x2591 - x2
     &561*x2598 - x2567*x760 - x2590*(x1208*x2115*x284 - x2108*x2566 + x
     &2124*x2580 + x2125*x763 - x2566*x2574 - x2566*x2575 + x2568* x763 
     &+ x2569*x763 - x2570*x4 - x2571*x4 - x2573*x4 + x2578 - x289 *(x25
     &83 - x2584 + x2585 + x2589)) + x2594*x2597*x790 - x2595* x2596 + x
     &2599)
      x2601 = x1244*x624
      x2602 = rho_s*x1796
      x2603 = x2153*x9
      x2604 = x1808 - x2410 + x2411
      x2605 = x1258*x624
      x2606 = x1267*x548
      x2607 = x1039*x1266
      x2608 = x1256*x968
      x2609 = gamma_cc*x939
      x2610 = 0.00385537758592051d0*x2609*x619
      x2611 = x1297 - x2384 + x2389 + x2392 + x2444 + x2445 - x2447 + x2
     &448 + x2449 + x2451 + x2532
      x2612 = x1255*x382
      x2613 = 3.0d0*x883
      x2614 = 48.0d0*x180
      x2615 = x2524*x2614
      x2616 = x1252*x2334
      x2617 = x1285*x59
      x2618 = 0.101321996098707d0*x2617
      x2619 = x1114*x47
      x2620 = 5.1297628024841d0*x2619
      x2621 = x1252*x569
      x2622 = 8.0d0*x2621
      x2623 = -x2616 + x2618 + x2620 + x2622
      x2624 = 30.7785768149046d0*x47
      x2625 = x2525*x2624
      x2626 = x2240*x2527
      x2627 = 0.0204827142749824d0*x592
      x2628 = x2627*x926
      x2629 = 4.06893148444717d0*x930
      x2630 = x1237*x2629
      x2631 = 4.86345581273796d0*x2529
      x2632 = 12.0d0*x2530
      x2633 = 192.0d0*x1252*x798
      x2634 = 17.9999993069297d0*x1728
      x2635 = x2634*x926
      x2636 = 246.228614519237d0*x2531
      x2637 = 3.79234362289252d0*x180*x2475
      x2638 = 0.0319433315943007d0*x1252*x592
      x2639 = 28.071472307785d0*x1252*x1728
      x2640 = x2637 + x2638 - x2639
      x2641 = x947*(-3.0d0*x1245*x2173 + x2183*x2613*x2623 - x2600 - 3.0
     &d0*x2611*x948 - x2612*x2613 - x2615 + x2625 - x2626 + x2628 - x263
     &0 - x2631 + x2632 + x2633 - x2635 - x2636 + x2640)
      x2642 = x1315*x619
      x2643 = x1313*x399
      x2644 = x1021*x2643*x619
      x2645 = x1301*x2212
      x2646 = x1051*x392
      x2647 = x1300*x1303
      x2648 = 0.0059503309916329d0*x1026*x1293*x1294 + x1273*x2209 - x22
     &08*x2647
      x2649 = -x1300*x1822*x2259 + 0.123372082749456d0*x1306*x2420 + 0.2
     &57230792532616d0*x1308*x975 + x2648
      x2650 = x1016*x395
      x2651 = 1.72720915849239d0*x2204
      x2652 = 0.555174372372553d0*x634
      x2653 = 3.60123109545663d0*x633
      x2654 = 2.59081373773858d0*x2281
      x2655 = 2.05784634026093d0*x2204
      x2656 = x1848*x939
      x2657 = x2656*x624
      x2658 = x1215*x2231
      x2659 = x1291*x2233
      x2660 = x1036*x1842
      x2661 = 3.0867695103914d0*x2281
      x2662 = 6.0d0*x1104
      x2663 = -x2637 - x2638 + x2639
      x2664 = 0.514461585065233d0*x633
      x2665 = rho_s*x2248
      x2666 = x1271*x688
      x2667 = x2666*x939
      x2668 = 0.0570775155072701d0*x1847
      x2669 = x1300*x2258
      x2670 = x1300*x2303
      x2671 = x1017*x1300
      x2672 = x1903*x2647
      x2673 = x1293*x2249
      x2674 = x2609*x548
      x2675 = -x1286 + x1298
      x2676 = x2273*x2675
      x2677 = x2277*x2675
      x2678 = x1012*x1300
      x2679 = x1842*x2678
      x2680 = x2253*x2256
      x2681 = x2522 - x2543 - x2544 + x2545 + x2546 + x2547 - x2548 + x2
     &549 - x2550 + x2551 + x2552 - x2553 - x2554 + x2555 - x2556 - x255
     &8 + x2560
      x2682 = x4*x984
      x2683 = x772*x987
      x2684 = 3.0d0*x1277
      x2685 = x1278*x292
      x2686 = x253*(-x1245*x2563 + x1967*(-x1245*x2538 + x2539 + x2611*x
     &884) + x2065*( -x1245*x2534 - x1282*x2523 + x2537 + x2611*x372 + x
     &917*( 1.45369647426571d0*x1233 + x1979*x362*x4 + x1992*(-x1989 - x
     &1990 + x2489 - x2490) + x2491 + x2681 - 29.0739294853142d0*x4*x907
     &)) + x2086*(-x1281*x2540 + x2082*x2682 + x2562 + x862*(x1219 + x25
     &41 + x2542 + x2681)) + x2087*x2611 + x2141*x2577*x292*x993 - x2289
     &* x2592*x292 + x2290*x2596*x2685 + x2291*x2683 + x2292*x2591 + x25
     &67*x988 + x2598*x2682 + x2599 + x2684*(x1207*x2286*x2576 - x2125*x
     &987 - x2287*x2580 - x2569*x987 + x2572*x2683 + x2575*x2581 + x2578
     & + 12.0326097596534d0*x2582 + x2589 + x289*(-x2583 + x2584 - x2585
     & - x2586 + x2587 - x2588) + 72.1956585579204d0*x773*x987))
      x2687 = x2399*x2614
      x2688 = x2240*x2387
      x2689 = x2397*x2624
      x2690 = x1283*x2629
      x2691 = 2.43172790636898d0*x2529
      x2692 = 2.43172790636898d0*x1106*x2528
      x2693 = 12.0d0*x2449
      x2694 = 192.0d0*x2621
      x2695 = 123.114307259618d0*x2619
      x2696 = 123.114307259618d0*x2481*x47
      x2697 = x1011*x2662 - x1107*x2627 + x1107*x2634 + 6.0d0*x1299*x641
     & + x2238*x2299 *x4 + x2295*x2443 + x2663 + x2686 - x2687 + x2688*x
     &642 - x2689 + x2690 + x2691 - x2692 - x2693 + x2694 + x2695 - x269
     &6
      x2698 = x1291*x972
      x2699 = 3.0d0*x982
      x2700 = x2308*(x2616 - x2618 - x2620 - x2622)
      x2701 = x981*(3.0d0*x1275*x2306 - x1285*x2627 + x1285*x2634 + x264
     &0 - x2686 + x2687 + x2688*x52 + x2689 - x2690 - x2691 + x2692 + x2
     &693 - x2694 - x2695 + x2696 + x2698*x2699 + x2699*x2700 + 3.0d0*x9
     &83*(x2384 - x2389 - x2392 + x2453))
      x2702 = x1698*x798
      x2703 = x1199*x1317
      x2704 = x2327*(x1100*x2320*x622 + x1268*x1294*(x1260 + x1261*x1293
     & - x1310) + x1293*x2321 + x2323*(-x1274 - x1296 + x1305 - x1307 - 
     &x1309 + x1311))
      x2705 = -786432.0d0*x1054*x2319*x2704 + 1572864.0d0*x1055*x2403 - 
     &x1056*x2600 + 100663296.0d0*x1097*x799 - 129094707.849062d0*x1181*
     &x1698 - 2133291.95011784d0*x1237*x185*x930 + 1572864.0d0*x1318*x17
     &86 + 16136838.4811327d0*x1698*x2525 - 2097152.0d0*x1705*x1756*x403
     & + 786432.0d0*x2314*(0.00674691077536089d0*x1021*x2642 + x1215*x21
     &61 + x1243*x2200*x2610 - x1243*x2602*x2603 - 0.00803846226664426d0
     &* x1257*x2171 - 0.042871798755436d0*x1258*x2150 + 0.00803846226664
     &426d0*x1258*x2202 + x1266*x2211*x2646*x544 - x1291*x2163 + 0.01071
     &7949688859d0*x2154*x2641 + x2159*x2608 + x2165*x2643*x2645 + x2167
     &*x2169*x2196*x2643 - x2168*x2610 - 0.00385537758592051d0*x2169*x24
     &17*x2607 - 0.00803846226664426d0* x2201*x2644 - x2216*x2649*x2650 
     &+ x2246*(-1.48046499299348d0* x1048*x2657 - 2.05784634026093d0*x12
     &57*x2225 - 3.0867695103914d0* x1257*x2235 - x1257*x2653 + x1264*x2
     &227*x2660 + 1.48046499299348d0*x1265*x2230 + x1265*x2654 - 2.22069
     &748949021d0 *x1265*x2660 - x1265*x2661*x969 - x2220*x2602 + x2224*
     &x2608 + x2232*x2658 - x2234*x2659 - x2243*(x1035*x2662 - 6.0d0*x12
     &63* x2236 + x2237*(x1032 + x1851) + x2239*x4 + x2600 + x2615 - x26
     &25 + x2626 - x2628 + x2630 + x2631 - x2632 - x2633 + x2635 + x2636
     & + x2663) - x2604*x2652 + x2641*x2664 - x2651 + x2655*x969) - x231
     &3* (x1004*x1292*x2680 + 0.0952052958661265d0*x1005*x2667 - 0.14280
     &794379919d0*x1005*x2670 + x1245*x2275 + x1245*x2278 + x1292*x2276*
     &x2279 + 0.0952052958661265d0*x1293*x2251 - 0.14280794379919d0*x129
     &3*x2262 - x1300*x1900*x2260 + 2.05784634026093d0*x1308*x2257 - x13
     &08*x2270*x2661 - x1308*x2653 - x1797*x1848*x2665*x2671 + x1830*x24
     &14*x2665 + x1884*x2250*x2668 *x939 - x2248*x2270*x2657 + 0.0119006
     &619832658d0*x2249*x2701 + x2252*x2655 - x2252*x2661*x2671 + x2266*
     &x2679 - x2267*x2679 + x2268*x2678 - 0.0684930186087241d0*x2271*x26
     &74 - 0.0684930186087241d0*x2272*x2672 + x2279*x2676 + x2279*x2677 
     &- x2302*x2697 - x2304*x2697 + x2412*x2652*x974 + 0.006421220494567
     &88d0*x2457*x2666 - x2651*x974 + x2654*x2669 + x2664*x2701*x974 - 0
     &.0799085217101781d0*x2667 + 0.119862782565267d0*x2670 - 0.16660926
     &7765721d0*x2673) - x2409* x2601 - 0.0115661327577615d0*x2601*x2604
     & + 0.042871798755436d0* x2601*x939*x969 - 0.0750256478220131d0*x26
     &05 + 0.00803846226664426d0*x2606*x969 - 0.00674691077536089d0*x260
     &6) - x2324*x2703*x698 + 1988280.25335907d0*x2476 - 2549851.5211487
     &6d0* x2477*x586 + 16747.5054349127d0*x2479 + 10738.8413018019d0*x2
     &482* x926 - 14717536.073304d0*x2484 - 9437183.63663157d0*x2485*x92
     &6 - 25165824.0d0*x2524*x566 + 6291456.0d0*x2526*x2702 - 124267.515
     &834942d0*x2527*x731
      x2706 = 1.58945719401042d-7*rho_c
      x2707 = x1142 - x1241 - x1242 + x1322
      x2708 = -x1240 + x1320 - x1321 + x2707
      x2709 = x1072*x1613
      x2710 = x1066**2
      x2711 = 0.00385537758592051d0*x2350
      x2712 = -x1079 + x524 + x526*(-0.5d0*x1077**2*x1653 + x1078*(x1073
     &*x1602*x2343 + x1074*(-x1635*x2710 + x1636 - 0.0231322655155231d0*
     &x2349 + x2352* x2710 - x2711) + x1590 + x1591*x2338 + x1592*x2338 
     &+ x1594*x2339 + x1598*x2339 + x2347*x2709 + x515*(x1577*x2349 + x1
     &579*x2357 + x1582*x2358 + x1640*x2710 - x1645*x2710 + x1646*x2710 
     &+ x1649 - x2359*x2710 + x2711*x81)) + x1664)
      x2713 = x204*(x1586 + x2712)
      x2714 = x1093*x1693
      x2715 = x559*(gamma_cc*x542 + 1.28199961805103d0*x1087*x2371 + 0.7
     &47833110529767d0*x549*(x1684 + x1687))
      x2716 = x1091*x1680
      x2717 = x2369*x2716
      x2718 = -2.0d0*x1086*x2367 + x185*(6.0d0*x187 + 67.8d0*x47 + 1.74d
     &0) - x2715 + x2717
      x2719 = x2718*x563
      x2720 = x210*x2719
      x2721 = x1096*x59
      x2722 = 0.40528798439483d0*x2721
      x2723 = x2722*x566
      x2724 = x2240*x44
      x2725 = x1710 + x1713
      x2726 = x2725*x59
      x2727 = 0.0168869993497846d0*x2726
      x2728 = x2727*x52
      x2729 = x185*x2728
      x2730 = 0.854960467080683d0*x2725
      x2731 = x2730*x569
      x2732 = x1096*x1106
      x2733 = x2732*x650
      x2734 = 12.0d0*x44
      x2735 = x2734*x569
      x2736 = x1761 + x1764
      x2737 = x185*x2736
      x2738 = x2737*x707
      x2739 = x1754*x185
      x2740 = 0.333333333333333d0*x2739
      x2741 = x2401*x639
      x2742 = x1137**2*x1791*x662
      x2743 = x1137*x622*x699
      x2744 = 0.00514050344789401d0*x2736
      x2745 = x2153*x622
      x2746 = x1125*x662
      x2747 = 0.0308430206873641d0*x939
      x2748 = x2746*x2747
      x2749 = x681*x74
      x2750 = rho_c**(-6.66666666666667d0)
      x2751 = x1118**2
      x2752 = x219*x2751
      x2753 = x2438*x2752
      x2754 = x187*x2724
      x2755 = 4.0d0*x1104
      x2756 = x2736*x622
      x2757 = 0.666666666666667d0*x2756
      x2758 = x1124*x703
      x2759 = x1106*x180
      x2760 = x2722*x2759
      x2761 = x187*x2735
      x2762 = 20.5190512099364d0*x1096
      x2763 = x2759*x569
      x2764 = x2762*x2763
      x2765 = x2731*x47
      x2766 = x2760 - x2761 - x2764 - x2765
      x2767 = x1117*x2755 + x1856*x2758 + x216*x2757 + x2727*x642 + x275
     &4 + x2766
      x2768 = x2455*x2767
      x2769 = x2441*x2751*x658
      x2770 = 0.246744165498913d0*x939
      x2771 = x1259*x2153 + x2744*x81
      x2772 = x2744*x633
      x2773 = x2153*x634
      x2774 = 0.000118911490640146d0*x1883*x2736
      x2775 = x2152*x687
      x2776 = 0.00118911490640146d0*x2746*x2775
      x2777 = x2753*x624
      x2778 = x2769*x624
      x2779 = x2460*x2751
      x2780 = x2461*x2752
      x2781 = x1820*x2770
      x2782 = x2431*x2781
      x2783 = 0.011415503101454d0*x1895*x2458*x939
      x2784 = x2768*x624
      x2785 = x2462*x2767
      x2786 = x699*(x1122*x1135*x1817*x2749 + x1134*x1825*x219*x2427*x27
     &49 + x1135* x1800*x2746*x619*x939 + x1905*(-x2772 - x2773 - x2774 
     &- x2776 - x2777 + x2778 - x2779 + x2780 + x2782 + x2783 + x2784 + 
     &x2785) + x2415*x2748 + x2744*x625 + x2745*x625 + x2750*x669*(-x243
     &3*x2770 + x2753 - x2768 - x2769 + x2771))
      x2787 = -x1106*x2723 - x1673*x2724 + x1673*x2735 + x1698*x2731 + 0
     &.5d0*x2483* x2743 + x2729 + 20.5190512099364d0*x2733 - 0.083333333
     &3333333d0* x2738 + x2740*x2741 - 0.125d0*x2742 + 0.125d0*x2786
      x2788 = x1149*x2379
      x2789 = -x1912*x2788 - x214*x2787 - x2718*x713
      x2790 = 0.0875d0*x1094
      x2791 = x1401*x203
      x2792 = -x2787*x718 - x2788*x2791 - x713*(-x1085*x1401*x2367 + x18
     &5*(0.2625d0* x187 + 2.96625d0*x47 + 0.076125d0) - 0.04375d0*x2715 
     &+ 0.04375d0* x2717)
      x2793 = 128.0d0*x238
      x2794 = x44*x59
      x2795 = x1673*x2794
      x2796 = x1170*x624
      x2797 = gamma_cc*x1169
      x2798 = x548*x662
      x2799 = x1191*x1193
      x2800 = x1191*x2749
      x2801 = x1169*x220
      x2802 = x1172*x2801
      x2803 = x1189*x2802
      x2804 = x1185**2
      x2805 = x180*x586
      x2806 = x2722*x2805
      x2807 = x2734*x798
      x2808 = x187*x2807
      x2809 = x1249*x2730
      x2810 = x2805*x798
      x2811 = x2762*x2810
      x2812 = -x2728 + x2754
      x2813 = x2422 - x2426
      x2814 = x1171**(-2)
      x2815 = x1170*x2814
      x2816 = x2726*x52
      x2817 = x185*x2816
      x2818 = x44*x798
      x2819 = x2702*x2725
      x2820 = x586*x799
      x2821 = x1032*x939
      x2822 = x1177*x1196
      x2823 = x1172*x223/x1175**2
      x2824 = -43031569.2830205d0*x1096*x2820 - x1178*x698*(x1180*x2750*
     &(x1167*x2437* x2804*exp(-x1837*x63) - x1185*x1186*x2781 - x1186*x2
     &441*x2804 - x1187*(x1184*x2755 - x2757*x63 - 5.33333333333333d0*x2
     &758*x63 - x2806 + x2808 + x2809 - x2811 + x2812) + x2771) + x1188*
     &x2747* x2797*x2798 - x1194*(x2772 + x2773 + x2774 + x2776 + x2777 
     &- x2778 + x2779 - x2780 - x2782 - x2783 - x2784 - x2785) + x2419* 
     &x2800*x2815*x82*(x2421 - x2424 + x2813) + x2744*x2796 + x2745* x27
     &96 - x2748*x2799*x619 - 0.0231322655155231d0*x2800*x2803) + x1195*
     &x1790*x2486*x2823*x662 - 699050.666666667d0*x1197*x2739* x2821 + 1
     &74762.666666667d0*x1198*x2737 + 1048576.0d0*x1199*x2822* x622 + 25
     &165824.0d0*x1673*x2818 - 849950.507049586d0*x2721*x566* x586 + 497
     &070.063339768d0*x2795 - 35414.6044603994d0*x2817 + 1792982.0534591
     &9d0*x2819
      x2825 = x2598*x264
      x2826 = x2097 + x2099
      x2827 = x2826*x764
      x2828 = 12.0d0*x792
      x2829 = x2728 - x2754
      x2830 = x2806 - x2808 - x2809 + x2811 + x2829
      x2831 = x2108*x2826
      x2832 = mu*x1209*x277
      x2833 = x2290*x771
      x2834 = x2290*x280
      x2835 = x2110*x2834 + x2831
      x2836 = x2290*x2832
      x2837 = x2123*x2579*x2836 + x2574*x2826 + x2575*x2826 - x289*(5.57
     &009647607783d0 *x2826*x768 + 33.420578856467d0*x2833 + x2835)
      x2838 = 0.0282745794067445d0*x818
      x2839 = x2838*x813
      x2840 = x2838*x823
      x2841 = x2009*x2520
      x2842 = x2014*x849
      x2843 = x1224*x2057
      x2844 = x1227*x2057
      x2845 = x2506*x835
      x2846 = x2034*x846
      x2847 = x2557*x839
      x2848 = x2559*x849
      x2849 = x2017 + x2019 + x2839 + x2840 + x2841 + x2842 + x2843 + x2
     &844 - x2845 - x2846 + x2847 + x2848
      x2850 = x2016*x274
      x2851 = x253*(-x1211*x2132*x2593*x280 + x1967*(x2830*x884 + x868) 
     &+ x2065*(x2830 *x372 + x886 - x917*(x2849 - x2850 - x908)) + x2086
     &*(-x862*(x2849 - x810) + x866) + x2087*x2830 - x2130 - x2590*(x120
     &8*x2116*x2832* x987 + x2570 + x2571 + x2573 + x2831 + x2837) - x25
     &95*x2597*x280 + x2595*x2828 + x2825 + x2827*x760)
      x2852 = x1256**2
      x2853 = x1257*x939
      x2854 = 0.00771075517184102d0*x2609
      x2855 = -x2760 + x2761 + x2764 + x2765 + x2829
      x2856 = x2721*x2805
      x2857 = x1096*x2810
      x2858 = x187*x2818
      x2859 = x1249*x2725
      x2860 = x187*x2794
      x2861 = 0.563350298308813d0*x2816 - 7.9070364537309d0*x2860
      x2862 = x1245*x2184*x2623 - 16.68d0*x1245*x2612 - x2172*x2855 - 2.
     &78d0*x2851 + 13.5204071594115d0*x2856 + 684.515548363478d0*x2857 -
     & 400.32d0* x2858 - 28.5214811818116d0*x2859 + x2861
      x2863 = x2607*x544
      x2864 = x1016*x2643
      x2865 = x1313*x2214
      x2866 = x1264**2*x1841
      x2867 = 2.96092998598695d0*x2656
      x2868 = 2.84425771716939d0*x2860
      x2869 = 16.0d0*x2758
      x2870 = 12.0d0*x1104
      x2871 = x1292**2
      x2872 = x1245*x972
      x2873 = x1300**2
      x2874 = x1842*x2873
      x2875 = x1292*x2253
      x2876 = x2721*x2759
      x2877 = x2725*x569
      x2878 = x2877*x47
      x2879 = x1096*x2763
      x2880 = x44*x569
      x2881 = x187*x2880
      x2882 = x2685*x2834
      x2883 = -x2017 - x2019 - x2839 - x2840 - x2841 - x2842 - x2843 - x
     &2844 + x2845 + x2846 - x2847 - x2848
      x2884 = x253*(-x1207*x2132*x2882 + x1276*x2129 + x1279*x2828*x292*
     &x763 + x1967*( x2855*x884 + x868) + x2065*(x2855*x372 + x886 + x91
     &7*(x2850 + x2883 + x908)) + x2086*(x862*(x2883 + x810) + x866) + x
     &2087*x2855 + x2597*x2882 + x2684*(-x1207*x2116*x2836 + x2290*x2572
     & + 72.1956585579204d0*x2833 + x2835 + x2837) + x2825 - x2827*x988)
     &
      x2885 = 2.0d0*x1007*x2756 + x1007*x2869 + x1299*x2870 + 0.20264399
     &2197415d0* x2726*x642 + x2868 + 4.86345581273796d0*x2876 - 10.2595
     &256049682d0*x2878 - 246.228614519237d0*x2879 - 144.0d0* x2881 + x2
     &884
      x2886 = x983*(x2766 + x2812)
      x2887 = x1275*x2698
      x2888 = x1275*x2700
      x2889 = x1316*x2319
      x2890 = 0.0833333333333333d0*x1055*x2737 + 0.5d0*x1199*x1318*x698 
     &- x1317*x1319* (x1215*x2160*x2605 - x1256*x1291*x2162 + 0.01607692
     &45332885d0* x1257*x2606 - 0.0160769245332885d0*x1258*x2644 + x2154
     &*x2603 + x2159*x2852 - x2167*x2642*x2854 + x2197*x2862 - x2203*x28
     &63*x2864 + x2245*x2750*(-x1256*x2659*x81 - 6.17353902078279d0*x125
     &7*x1265* x1820 + x1257*x2658*x81 + 0.246744165498913d0*x1259*x2736
     & + 4.11569268052186d0*x1259*x2853 - x1265*x2867 + x2220 + x2223* x
     &2852*x81 + x2227*x2866 - x2229*x2866 - x2242*(-2.0d0*x1031*x2736 +
     & x1263*x2870 - 0.202643992197415d0*x2816 + x2851 - 4.8634558127379
     &6d0*x2856 - 246.228614519237d0*x2857 + 144.0d0* x2858 + 10.2595256
     &049682d0*x2859 + x2868 - x2869*x980) + x2244* x2862*x81) - x2313*(
     &x1271*x2668*x2775 + 0.190410591732253d0*x1293 *x2667 - 0.285615887
     &598379d0*x1293*x2670 + 4.11569268052186d0* x1308*x2204 - x1308*x22
     &82*x2671 + 1.5433847551957d0*x1308*x2872* x633 + x2219*x2255*x9 + 
     &x2253*x2276*x2871 + 0.246744165498913d0* x2255*x2756 + x2266*x2874
     & - x2267*x2874 + x2268*x2873 - x2269* x2873 - x2302*x2885 - x2304*
     &x2885 + x2310*(x2861 - 13.5204071594115d0*x2876 + 28.5214811818116
     &d0*x2878 + 684.515548363478d0*x2879 + 400.32d0*x2881 - 2.78d0*x288
     &4 + 8.34d0 *x2886 + 16.68d0*x2887 + 8.34d0*x2888) + x2311*( 1.1267
     &0059661763d0*x2816 - 15.8140729074618d0*x2860 - 27.040814318823d0*
     &x2876 + 57.0429623636231d0*x2878 + 1369.03109672696d0*x2879 + 800.
     &64d0*x2881 - 5.56d0*x2884 + 16.68d0*x2886 + 33.36d0*x2887 + 16.68d
     &0*x2888) + 0.00570775155072701d0*x2666*x2736 - x2669*x2867*x624 - 
     &0.136986037217448d0*x2672*x2674 + 0.0357019859497974d0*x2673* x287
     &2 + x2676*x2875 + x2677*x2875 + x2680*x2871) + x2601*x2744 + 0.085
     &7435975108721d0*x2601*x2853 + x2607*x2798*x2854 + x2645* x2649*x28
     &65) + x1673*x2807 + 0.0833333333333333d0*x185*x2851 + 0.125d0*x270
     &4*x2889 - x2723*x586 - x2729 - x2740*x2821*x403 - x2762*x2820 + 0.
     &237021476430782d0*x2795 + 0.854960467080682d0* x2819
      x2891 = 262144.0d0*x1428
      x2892 = x1534*x1585
      x2893 = -x463*(x2892 - x458 + x460)
      x2894 = x1078*x1651
      x2895 = -x526*(x1080 + x2894)
      x2896 = x469 + x529
      x2897 = x2893 + x2895 + x2896
      x2898 = x1359*x2897
      x2899 = x1355*x533
      x2900 = 0.75d0*x2899
      x2901 = x437*x442
      x2902 = 1.11479522475137d-5*x71
      x2903 = x2902*x440
      x2904 = x119*x1534
      x2905 = x1341*x2904
      x2906 = 0.0115661327577615d0*x129*x428
      x2907 = x128*x1545
      x2908 = x447*(x1556 + x445)
      x2909 = x1551*x1558*x429
      x2910 = x1337*x2909
      x2911 = x1338*x432
      x2912 = x1338*x434
      x2913 = x1474*x1560
      x2914 = x452*x99
      x2915 = x1352*x158
      x2916 = x1651*x2915
      x2917 = 0.0115661327577615d0*x166
      x2918 = x2917*x494
      x2919 = x165*x173
      x2920 = x2919*x514
      x2921 = x515*(x1607 + x513)
      x2922 = x1348*x1610*x495*x71
      x2923 = x2922*x492
      x2924 = x1481*x170
      x2925 = x1480 + x2924
      x2926 = x1614*x2925
      x2927 = x442*x505
      x2928 = x2902*x508
      x2929 = x1349*x498
      x2930 = x1349*x500
      x2931 = x2927 + x2928 - x2929 - x2930
      x2932 = x143*x471
      x2933 = x1341 + x1352
      x2934 = x204*x40
      x2935 = x2934*(-x2914*(x1533 + x2901 + x2903 - x2905 + x2906 - x29
     &07 - x2908 - x2910 - x2911 - x2912 + x2913) - x2932*(x1650 - x2916
     & + x2918 - x2920 - x2921 - x2923 + x2926 + x2931) + x2933)
      x2936 = x1380*x533
      x2937 = x543 - x547 + x555
      x2938 = x2937*x559
      x2939 = x2938 + x538
      x2940 = 0.25d0*x2939
      x2941 = x1363*x2939
      x2942 = x1380*x1692
      x2943 = x2374*x553 + x542*(-x539 + x540) + x547
      x2944 = x2943*x557
      x2945 = 0.30046866048071d0*x1679
      x2946 = x195*x2945
      x2947 = x1088*x553
      x2948 = x558*(x2947*x548 + x546)
      x2949 = x193*x537
      x2950 = 0.240374928384568d0*x194
      x2951 = x1366*x2950
      x2952 = x2949*x2951
      x2953 = x2948 + x2952
      x2954 = -x2944*x2946 + x2953
      x2955 = x1370*x180
      x2956 = x1363*x2955
      x2957 = x1374*x2954*x2956
      x2958 = x2951*x558
      x2959 = x1370*x203
      x2960 = x186*x2959
      x2961 = x2897*x2960
      x2962 = 0.721124785153704d0*x1373
      x2963 = x1367*x1374*x1886
      x2964 = 1d0/x90
      x2965 = x87*(-x1179*x89*(x1179*x1324 + x84) - x1426*x80 + x88)
      x2966 = x1170*x620
      x2967 = x2796*x623
      x2968 = x2967*x940
      x2969 = x1192*x690
      x2970 = x1026*x1179
      x2971 = x1033*x63*x640 + x656
      x2972 = x1167*x2971*x68
      x2973 = x1169*x939
      x2974 = x2973*x684
      x2975 = x1102*x941
      x2976 = x1186*x2971*x664 + x2975 + x632
      x2977 = x1030*x2972 + 0.000624285325860766d0*x2970 + x2974*x689 + 
     &x2976
      x2978 = x1180*x1193*x2977*x82
      x2979 = x2966 + x2968 + x2969 - x2978
      x2980 = x2964*x2979
      x2981 = x1793*x698
      x2982 = 0.125d0*x2981
      x2983 = x1328*x2814
      x2984 = 0.000312142662930383d0*x2970
      x2985 = x1170*x632
      x2986 = 0.000178367235960219d0*x1040
      x2987 = x2973*x2986
      x2988 = x1102*x2967
      x2989 = x1172*x1327
      x2990 = x2985*x2989
      x2991 = x1172*x1326
      x2992 = x1102*x2418*x2796*x2991
      x2993 = 0.0231322655155231d0*x83
      x2994 = x2993*x668*x690
      x2995 = x1824*x2991*x690*x84
      x2996 = x1168*x1172*x2977
      x2997 = x1174*x1179
      x2998 = x2977*x2997
      x2999 = x2997*(x1819*x2974 + x2208*x2972 + 0.00124857065172153d0*x
     &2970 + x2976)
      x3000 = 0.000535101707880657d0*x72
      x3001 = x2993 + x3000*x668
      x3002 = x1173*x2801*x2814*x2977*x3001*x78
      x3003 = x1172*(x2984 + x2985 + x2987 + x2988 - x2990 - x2992 + x29
     &94 - x2995 - x2996 - x2998 - x2999 + x3002)
      x3004 = -x1390*x3003 + x2965*x40*x705 + x2980*x2982*x2983
      x3005 = x2964*x3004
      x3006 = x1908*x193
      x3007 = 2.0d0*x3006
      x3008 = 0.30046866048071d0*x1383
      x3009 = x1680*x2943
      x3010 = x537*x566
      x3011 = gamma_cc*x187*x2948*x557
      x3012 = gamma_cc*x1381*x536
      x3013 = x2950*x3012
      x3014 = rho_c**(-9.66666666666667d0)*x46
      x3015 = x189*x194
      x3016 = x1366*x3015
      x3017 = x3016*x559
      x3018 = -x1110*x188 - 2.0d0*x3010 + x3011 + x3013*x3014 + 0.480749
     &856769136d0* x3017
      x3019 = x198*x98
      x3020 = x1380*x3019
      x3021 = 0.125d0*x936
      x3022 = 4.0d0*x2331 + x604 - x653*x929 - x924
      x3023 = -x1172*x2980*x3021 + x3022 + x705*x706*x91
      x3024 = x199*x3023
      x3025 = x2939 + x710
      x3026 = x1390*x2965
      x3027 = x2964*x3026
      x3028 = x1913*x3027
      x3029 = -x1387*x3024 + x214*x3005 + x3020*(x2944*x3007 - x3008*x30
     &09 + x3018) - x3025*x3028
      x3030 = 2.8d0*x1447
      x3031 = 2.8d0*x1446 + x3030
      x3032 = 9.0d0*x3
      x3033 = x3031*x3032
      x3034 = 0.0328125d0*x2899
      x3035 = x2939*x563
      x3036 = 0.0109375d0*x1355
      x3037 = 0.2625d0*x1363
      x3038 = x1380*x562
      x3039 = 0.0875d0*x3038
      x3040 = x194*x558
      x3041 = x1399*x3040
      x3042 = 0.0315492093504746d0*x1373
      x3043 = x2943*x559
      x3044 = 0.0131455038960311d0*x1383
      x3045 = 0.3955d0*x182 + 0.15225d0*x187 + 0.0875d0*x47 + 0.09275d0
      x3046 = 0.0105164031168249d0*x194*x3012
      x3047 = -x180*x3045 - 0.0875d0*x3010 + 0.04375d0*x3011 + x3014*x30
     &46 + 0.0210328062336497d0*x3017
      x3048 = x720 + x721
      x3049 = 0.04375d0*x2938 + x3048
      x3050 = -x1403*x3024 + x3005*x718 + x3020*(x1401*x3043 - x3009*x30
     &44 + x3047) - x3028*x3049
      x3051 = 64.0d0*gamma_cc
      x3052 = x1404 - x1405
      x3053 = gamma_cc*(x1397 - x1398 + x1400 + x3052)
      x3054 = 192.0d0*x3
      x3055 = x463*(-x2892 + x461)
      x3056 = x526*(-x2894 + x524)
      x3057 = x2896 + x3055 + x3056
      x3058 = x204*x3057
      x3059 = 8.4d0*x3058
      x3060 = x1363*x534
      x3061 = 25.2d0*x3060
      x3062 = x1363*x3035
      x3063 = 8.4d0*x3062
      x3064 = -x2966 - x2968 - x2969 + x2978
      x3065 = x1177*x3021*x3064 + x1198*x705 + x3022
      x3066 = x3065*x718
      x3067 = 192.0d0*x3066
      x3068 = x3048 - x722
      x3069 = x1378*x198
      x3070 = x212*x3069
      x3071 = 192.0d0*x3068*x3070
      x3072 = x3067 - x3071
      x3073 = x1388 - x1393
      x3074 = x1357 - x1365 + x1376 + x3073
      x3075 = 8.4d0*x3060
      x3076 = 2.8d0*x3062
      x3077 = x3068*x3069*x724
      x3078 = 2.8d0*x3058
      x3079 = 64.0d0*x3066
      x3080 = x238*(x2472*x3031 + x3075 + x3076 + x3077 + x3078 - x3079)
     &
      x3081 = 0.13125d0*x1447
      x3082 = 0.13125d0*x1446 - x3053 + x3081
      x3083 = x1363*x535
      x3084 = x214*x3065
      x3085 = x538 + x710
      x3086 = x3085 - x560
      x3087 = x3070*x3086
      x3088 = -x3084 + x3087
      x3089 = x3058 + x3062 + x3083 + x3088
      x3090 = x1424*x2964
      x3091 = x2964*x2965
      x3092 = -262144.0d0*x2979*x2981*x2983/x90**2 + 262144.0d0*x3003*x3
     &090 - x3091* x40*x734
      x3093 = x1243*x632
      x3094 = 0.032153849066577d0*x975
      x3095 = x2201*x3094
      x3096 = x1269*x1324
      x3097 = x1021*x3096
      x3098 = x2316 + x3097
      x3099 = x1021*x3098
      x3100 = x3099*x399
      x3101 = x3098*x399
      x3102 = 0.00289153318944038d0*x83
      x3103 = x3102*x78
      x3104 = x1301*x3103
      x3105 = x72*x78
      x3106 = 3.34438567425411d-5*x3105
      x3107 = x1294*x3000 + x2993*x974
      x3108 = x1039*x1426 - x3099*x401 + x400
      x3109 = x1057*x1423*x1427*x40 - x1425*x934*(x1039*x2986*x939 + 0.0
     &00312142662930383d0*x1042 - 0.000371895686977057d0*x1043 - x1102*x
     &2154*x2167*x3101*x684 + x1243*x2975 - x1314*x3102*x396*( x1049 + x
     &2210) - x2165*x3101*x3104 + 0.00578306637888076d0*x2166* x83 - x21
     &98*x3102 + x2215*x2650*x3103*x3107 - x2646*x3106 - x3093 *x3100 + 
     &x3093 + x3095*x3100 - x3095) + 262144.0d0*x2213*x2981* x3108*x393*
     &(x1024 - x1053 + x937 + x942 - x971)/x1422**2
      x3110 = x3055 + x469
      x3111 = -x1068 + x512
      x3112 = -x1070 + x3111 + x510
      x3113 = -x3112*x515
      x3114 = x1076 + x3113
      x3115 = x1081 + x526*(x1078*x3114 + x1080)
      x3116 = x3110 + x3115
      x3117 = x1359*x3116
      x3118 = -x2906
      x3119 = -x2901 - x2903 + x2907 + x2908 + x2910 + x2911 + x2912 - x
     &2913 + x3118 + x449
      x3120 = x99*(x2905 + x3119)
      x3121 = x3120*x452
      x3122 = x2338*x2917
      x3123 = x1613*x3112
      x3124 = x143*(-x1073*x2922 + x1075 - x2915*x3114 - x2919*x3112 + x
     &2925*x3123 + x2931 + x3113 + x3122 - x515*(x1606 - x2346 + x3111))
     &
      x3125 = x3124*x471
      x3126 = x1341 - x1352 + x3121 + x3125
      x3127 = x2934*x3126
      x3128 = x1093*x1363
      x3129 = 0.240374928384568d0*x2949
      x3130 = x2371*x553 + x542
      x3131 = x1091*x557
      x3132 = rho_s*x3131
      x3133 = x1382*x3129 - x2946*x3132 + x3130*x558
      x3134 = x1705*x3133
      x3135 = x2959*x3134
      x3136 = x1363*x3135
      x3137 = x1370*x3116
      x3138 = x1374*x3137
      x3139 = x1093*x1913
      x3140 = x1032*x1143
      x3141 = x668*x73
      x3142 = x1327*x225
      x3143 = x224 + x3141 - x3142
      x3144 = x1176*x223
      x3145 = x3143*x3144
      x3146 = x10*x1196
      x3147 = x3145*x3146
      x3148 = x2*x2797
      x3149 = 0.00192768879296025d0*x1170*x2420
      x3150 = 0.00144576659472019d0*x81
      x3151 = 0.00144576659472019d0*x83
      x3152 = x3151*x668
      x3153 = x40*(2.22959044950274d-5*x1026*x1124*x3148 + 0.00289153318
     &944038d0*x1189 *x1830*x3148 + x1191*x1824*x2815*x3001*x3151 - 1.67
     &219283712705d-5*x1191*x2802*x3105 - x1193*x3152*(x1190 + x2813) - 
     &x1326*x2803*x3150*x78 - x2799*x3152 - x2989*x3149 + x3149)
      x3154 = x1177*(-x1329*x3140 - 0.125d0*x3147 + x3153)
      x3155 = x1085*x45
      x3156 = x1096*x1886
      x3157 = x3006*x3130*x557*x684
      x3158 = x3007*x557
      x3159 = x1380*x3069
      x3160 = -x1198*x3140 + x2707 + 0.125d0*x2822
      x3161 = x199*x3160
      x3162 = -x1387*x3161 + x214*x3154 - x3159*(x1091*x3158 - x2716*x30
     &08 + x3013* x3156 - 2.0d0*x3155 + x3157)
      x3163 = x1392*x3139 + x3162
      x3164 = 0.00546875d0*x1177*x1330
      x3165 = -x1403*x3161 + x3154*x718 - x3159*(x1092*x1401 - x2716*x30
     &44 + x3046* x3156 - 0.0875d0*x3155 + 0.04375d0*x3157)
      x3166 = x3139*x3164 + x3165
      x3167 = x1159*x1161
      x3168 = 192.0d0*x1409
      x3169 = 262144.0d0*x2823*x3143
      x3170 = x1177*x1200*x1329 - 2097152.0d0*x1177*x3153 + x3146*x3169
      x3171 = x1258*x3094
      x3172 = gamma_cc*x2
      x3173 = x2320*x84
      x3174 = x2319*x3098
      x3175 = -x1021*x3173*x3174 + x1270*x1426 + x3173
      x3176 = -x1424*x935*(-gamma_cc*x1123*x1314*x2167*x3098 - x1016*x12
     &68*x2607*x3101 *x78*x81 - x1021*x1315*x3102 + 0.000178367235960219
     &d0*x1039*x1272 *x2609 - x1039*x2864*x3106 + 0.000371895686977057d0
     &*x1042*x1257 + x1123*x1244 + 0.00578306637888076d0*x1830*x2863*x31
     &72 - x2312* x3102*(x1312 + x2648) + x2865*x3104*x3107 - x3100*x317
     &1 + x3171) + x2325*x2889*x2981*x3175*x934 + 524288.0d0*x2703*x3108
     &*x40
      x3177 = x1358*x210
      x3178 = x1355*x1359
      x3179 = -x1333 + x1339 - x137
      x3180 = x1536*x99
      x3181 = x1340*x3180
      x3182 = x3179*x3181
      x3183 = -x1344 + x1350 - x174
      x3184 = x143*x1652
      x3185 = x1351*x3184
      x3186 = x3183*x3185
      x3187 = 3020.29322777679d0
      x3188 = 3.86816087015788d-7*x3187
      x3189 = x133*x134*x136*x99**(-7.0d0)
      x3190 = x3188*x3189
      x3191 = x1335*x1337
      x3192 = x1337*x1471
      x3193 = x133*x1336*x1559
      x3194 = x1474*x3193
      x3195 = x132*(x1334 - x3190 - x3191 - x3192 + x3194)
      x3196 = x2914*x3195
      x3197 = x171*x173
      x3198 = x143**(-7.0d0)*x170*x3197
      x3199 = x3188*x3198
      x3200 = x1346*x1348
      x3201 = x1348*x1486
      x3202 = x1612*x170
      x3203 = x1347*x3202
      x3204 = x2925*x3203
      x3205 = x168*(x1345 - x3199 - x3200 - x3201 + x3204)
      x3206 = x2932*x3205
      x3207 = x2*(x3182 + x3186 + x3196 + x3206)
      x3208 = x2934*x3207
      x3209 = 0.120187464192284d0*x1354
      x3210 = x1366*x1375
      x3211 = x207*x208
      x3212 = x194**2
      x3213 = x197**(-6)
      x3214 = 0.00478594012208448d0*x188*x3213
      x3215 = x3212*x3214
      x3216 = 0.0722251327448578d0*x3215
      x3217 = rho_c**(-12.3333333333333d0)*x203
      x3218 = x3217*x46
      x3219 = x1366*x189
      x3220 = x1373*x3219
      x3221 = x202*x3211*x3220
      x3222 = 0.0361125663724289d0*x3212
      x3223 = x3222*x46
      x3224 = rho_c**(-11.3333333333333d0)
      x3225 = x1678*x3224*x340*x74
      x3226 = x1503*x202
      x3227 = x1389*x1793
      x3228 = x214*x3227
      x3229 = rho_c**(-7.0d0)
      x3230 = x1816*x3187*x3229
      x3231 = x1326*x223
      x3232 = x1326*x1815
      x3233 = x1327*x1788*x3001
      x3234 = x220*(x1324 - x1325*x3231 + x224*x3233 - x2993*x3232 - 3.0
     &945286961263d-6*x3230)
      x3235 = x3234*x696
      x3236 = x3143*(-x224 - x3141 + x3142)
      x3237 = x1789*x3236
      x3238 = x3143*x696
      x3239 = x1380*x1424
      x3240 = x199*x3239
      x3241 = -2.0d0*x3226*(-0.480749856769136d0*x1384 + x185 + x3223*x3
     &225) - x3228* x3235 - x3228*x3237 - 0.5d0*x3238*x3240*(-x1385 + x1
     &85)
      x3242 = 0.2625d0*x230
      x3243 = 0.00525820155841243d0*x1354
      x3244 = 0.00315984955758753d0*x3212
      x3245 = x3244*x46
      x3246 = x3227*x718
      x3247 = -x2405*x3143*x3240*(x1401 - x1402) - x3226*(-0.04206561246
     &72994d0*x1384 + x1401 + x3225*x3245) - x3235*x3246 - x3237*x3246
      x3248 = x237*x3051
      x3249 = 128.0d0*x1409
      x3250 = 0.547944148869793d0*x72
      x3251 = x1537*x622
      x3252 = x1329*x1793*x3169 + 128.0d0*x220*x3144*x3251*(-23.68743988
     &78956d0*x221* x3233*x83 + 0.00633759476966667d0*x3230 + x3231*x325
     &0*x668 + 47.3748797757912d0*x3232*x83 - x3250)
      x3253 = 0.0078125d0*x247
      x3254 = 2048.0d0*x2319*x3251
      x3255 = x1016*x3254*x933*(3.0945286961263d-6*x1021*x2317*x2319*x31
     &87*x3229*x973 **8.34d0 + x2322*x2993*x3098 - x3096 + x3097*x3174 -
     & x3099*x3107* x3173/x2318**2) + x2326*x3108*x3175*x3254
      x3256 = 6.103515625d-5*rho_c
      x3257 = x1432*x3211
      x3258 = x1440*x1508
      x3259 = x1444 - 0.0105164031168249d0*x3257 + x3258
      x3260 = x3084 - x3087
      x3261 = x2474*(-x3058 - x3062 - x3083 + x3260)
      x3262 = -x3077 + x3079
      x3263 = x1063*x2472 - x3075 - x3076 - x3078 + x3262
      x3264 = x1370*x229
      x3265 = x238*(x1452*x3264 + x1454 - 0.240374928384568d0*x3257)
      x3266 = x99*(x119*x3179*x450*x452 + x3119)
      x3267 = x158*x471
      x3268 = -x2918 + x2920 + x2921 + x2923 - x2926 - x2927 - x2928 + x
     &2929 + x2930 + x517
      x3269 = x143*(x3183*x3267*x518 + x3268)
      x3270 = x1442*(x2933 + x3266*x452 + x3269*x471)
      x3271 = x1355*x562
      x3272 = x3271*x561
      x3273 = x3129*x545
      x3274 = rho_c**(-3.66666666666667d0)
      x3275 = x194*x3274
      x3276 = 0.640999809025515d0*x3275
      x3277 = x3276 - x554*x70
      x3278 = x556*x557
      x3279 = rho_c**(-2.66666666666667d0)*x194*x2945
      x3280 = x3273 + x3277*x558 - x3278*x3279
      x3281 = x2955*x3280
      x3282 = x210*x3281
      x3283 = 0.240374928384568d0*x1465
      x3284 = x1370*x531
      x3285 = x3211*x548
      x3286 = x3006*x3277
      x3287 = x1495 + x2952
      x3288 = x1453*x709 - x3264*(-x2946*x3278 + x3286 + x3287)
      x3289 = 0.0105164031168249d0*x1465
      x3290 = x194*x2949
      x3291 = 0.0210328062336497d0*x1465
      x3292 = x1679*x195
      x3293 = 0.0131455038960311d0*x3292
      x3294 = x229*(x1399*x3290 - x3278*x3293 + 0.04375d0*x3286 + x3291)
     &
      x3295 = 0.0105164031168249d0*x1438
      x3296 = x3295*x709
      x3297 = x1370*x3294 - x1370*x3296
      x3298 = 0.015625d0*rho_c
      x3299 = 3.0517578125d-5*x247
      x3300 = x203*x3299
      x3301 = x204*x3116
      x3302 = x1094*x1363
      x3303 = x1093*x3070
      x3304 = -x214*x3160 - x3303
      x3305 = x238*(x3301 - x3302 + x3304)
      x3306 = -x3160*x718 - 0.04375d0*x3303
      x3307 = x238*(0.04375d0*x3301 - 0.04375d0*x3302 + x3306)
      x3308 = x1442*x3126
      x3309 = x1093*x3271
      x3310 = 0.240374928384568d0*x1085*x196 + 0.320499904512757d0*x1088
     &*x3274*x558 - x3131*x3279
      x3311 = x2956*x3310
      x3312 = x1085*x1382*x193
      x3313 = x2371*x558
      x3314 = x1378*(-x3131*x3293 + 0.0105164031168249d0*x3312 + 0.01402
     &18708224331d0* x3313)
      x3315 = x3160*x3295
      x3316 = x1370*x3314 + x1370*x3315
      x3317 = x1378*(-0.0046948228200111d0*x3131*x3292 + 0.0037558582560
     &0888d0*x3312 + 0.00500781100801184d0*x3313)
      x3318 = x1370*x3317
      x3319 = gamma_cc**4*x237
      x3320 = x1395*x3319
      x3321 = rho_c**(-10.3333333333333d0)
      x3322 = x189*x3321
      x3323 = x3211*x3322
      x3324 = x193*x70
      x3325 = rho_c**(-9.33333333333333d0)
      x3326 = 0.0011285176991384d0*x3325
      x3327 = x3015*x340*x557
      x3328 = -0.00751171651201775d0*x3324 + x3326*x3327
      x3329 = x1369*x3264
      x3330 = x3328*x3329
      x3331 = x1354*x1432
      x3332 = gamma_cc*x3331
      x3333 = x1354*x1373
      x3334 = x3219*x3333
      x3335 = x3238*x7
      x3336 = gamma_cc*x1431*x195*x3335
      x3337 = gamma_cc*x1442
      x3338 = 0.00390625d0*x3337
      x3339 = 0.0109375d0*x1442
      x3340 = x3214*x3244
      x3341 = x3325*x3327
      x3342 = -0.0210328062336497d0*x3324 + 0.00315984955758753d0*x3341
      x3343 = 0.0112675747680266d0*x1439*x1508 - x3329*x3342 + 0.0013145
     &5038960311d0* x3336
      x3344 = rho_c*x1380
      x3345 = x3299*x3344
      x3346 = gamma_cc**5/tau_c**5
      x3347 = x180*x3346
      x3348 = x1455*x3347
      x3349 = rho_c**(-8.33333333333333d0)
      x3350 = x1367*x3349
      x3351 = x1679*x3213
      x3352 = x3212*x3351
      x3353 = rho_c**(-7.33333333333333d0)*x3352
      x3354 = 0.0011285176991384d0*x3353
      x3355 = x1340**2*x3180
      x3356 = x1351**2*x3184
      x3357 = x3191 + x3192
      x3358 = x132*x2914
      x3359 = x3358*(-x1334 + x3190 - x3194 + x3357)
      x3360 = x3200 + x3201
      x3361 = x168*x2932
      x3362 = x3361*(-x1345 + x3199 - x3204 + x3360)
      x3363 = x2*(x3355 + x3356 + x3359 + x3362)
      x3364 = 0.5d0*x1449
      x3365 = 0.25d0*x1442
      x3366 = 0.0078125d0*x1451
      x3367 = x203*x247
      x3368 = x3256*x3367
      x3369 = x1462*x3285
      x3370 = -x1467*x1508
      x3371 = -x1346 + x1485 - x1486
      x3372 = x2932*x3371
      x3373 = x1479 + x3372
      x3374 = x3339*x3373 + 0.0210328062336497d0*x3369 + x3370
      x3375 = x238*(-rho_s*x1495*x3264 + x3365*x3373 + 0.480749856769136
     &d0*x3369)
      x3376 = x1491*x40
      x3377 = x3376*x533
      x3378 = x175*x2925
      x3379 = x1487 - x170*x3378
      x3380 = x3267*x3379
      x3381 = x1067*x488
      x3382 = x1345*x2340
      x3383 = x1482*x1552*x160*x162*x2
      x3384 = -x1483*x1579*x503 - x1483*x1582*x507
      x3385 = x173*x2925
      x3386 = x170*x3385
      x3387 = x1577*x166
      x3388 = x1552*x159
      x3389 = x170*x3388 + x3387
      x3390 = x143*(x1074*(x1591 + x1592 - x3381) - x1486*x1595*x3112 + 
     &x1594 + x1598 - x1651*x3380 - x2344*x2925*x71 + x2918 - x3112*x338
     &2 + x3123* x3389 - x3386*x498 - x3386*x500 + x510 - x515*(x1609 - 
     &x3383*x504 + x3384))
      x3391 = x133*x1475
      x3392 = x129*x1577
      x3393 = x120*x1552
      x3394 = x133*x3393 + x3392
      x3395 = x99*(x127*x1334*x1550 + x1471*x1546 + x1475*x2909 + x1478*
     &x2904 - x1544 - x1548 - x1560*x3394 - x1574*(x1541 + x1542 + x409*
     &x424) + x3118 + x3391*x432 + x3391*x434 - x438 - x441 + x447*(x155
     &2*x1555 + x1557 + x1579*x437 + x1582*x440))
      x3396 = x1478 + x3395*x452
      x3397 = -x1489 + x3390*x471 + x3396
      x3398 = x3376*x562
      x3399 = 0.480749856769136d0*x2949
      x3400 = 1.76274947482017d0*x195
      x3401 = x2373*x544
      x3402 = x2937*x557
      x3403 = 0.60093732096142d0*x1679
      x3404 = x3275*x3403
      x3405 = x1382*x3399 + x3402*x3404 + x558*(x3400 - x3401)
      x3406 = x1705*x3405
      x3407 = x1370*x3406
      x3408 = x1363*x3407
      x3409 = x1382*x558
      x3410 = x1460*x3409
      x3411 = x3057*x3410
      x3412 = x1368*x1462
      x3413 = x1886*x3290
      x3414 = x3006*(-x3400 + x3401)
      x3415 = x3278*x545
      x3416 = x1378*(-x3403*x3415 + 0.961499713538272d0*x3409 + 0.480749
     &856769136d0* x3413 - x3414)
      x3417 = x1495*x3065
      x3418 = -x1460*x3416 + x1460*x3417
      x3419 = 0.0109375d0*x3398
      x3420 = 0.0420656124672994d0*x3409
      x3421 = 0.0262910077920621d0*x1679
      x3422 = x1378*(0.0210328062336497d0*x3413 - 0.04375d0*x3414 - x341
     &5*x3421 + x3420)
      x3423 = x3065*x3291
      x3424 = -x1460*x3422 + x1460*x3423
      x3425 = x1484*x173
      x3426 = x170*x3425
      x3427 = -x1483*x3388 + x3387
      x3428 = x143*(x1071 - x1072*x3382 + x1074*(-x1591 - x1592 + x3381)
     & + x1484*x2342 *x71 - x1486*x2339 - x1594 - x1598 + x2709*x3427 + 
     &x3114*x3380 - x3122 + x3426*x498 + x3426*x500 + x515*(x1069*x3383 
     &+ x1608 - x2345 + x3384))
      x3429 = x1489 + x3396 + x3428*x471
      x3430 = x1088*x2372*x70 + 0.480749856769136d0*x3275
      x3431 = -x3132*x3404 + x3399*x545 + x3430*x558
      x3432 = x2956*x3431
      x3433 = x3116*x3410
      x3434 = x1366*x3290
      x3435 = x3006*x3430
      x3436 = x3132*x545
      x3437 = x1378*(-x3421*x3436 + 0.0210328062336497d0*x3434 + 0.04375
     &d0*x3435)
      x3438 = -x1370*x3437 - x1460*x3160*x3291
      x3439 = 0.25d0*x1093
      x3440 = x1378*(-x3403*x3436 + 0.480749856769136d0*x3434 + x3435)
      x3441 = -x1370*x3440 - x1496*x3160
      x3442 = 0.015625d0*x1451
      x3443 = rho_c*x3300
      x3444 = x1367*x3215
      x3445 = gamma_cc*x1096
      x3446 = x3224*x3444*x3445
      x3447 = x193*x544
      x3448 = x3321*x3327
      x3449 = x1378*(-0.0150234330240355d0*x3447 + 0.00225703539827681d0
     &*x3448)
      x3450 = x1461*x3449
      x3451 = x548*x684
      x3452 = x1491*x3220
      x3453 = x1391*x1432*x684*x7
      x3454 = x1477*x3181
      x3455 = x3185*x3379
      x3456 = 7.73632174031576d-7*x3187
      x3457 = x132*(x1333*x1475 - x1472 + x1476 + x3189*x3456 - x3193*x3
     &394 + x3357)
      x3458 = x2914*x3457
      x3459 = x168*(x1344*x3385 - x1481 + x3198*x3456 - x3203*x3389 + x3
     &360 + x3378)
      x3460 = x2932*x3459
      x3461 = x2*(-x3454 + x3455 - x3458 + x3460)
      x3462 = 0.0328125d0*x1442
      x3463 = 0.021875d0*x1492
      x3464 = x3339*x3461
      x3465 = 0.0420656124672994d0*x1464
      x3466 = x1373*x3451
      x3467 = x1378*(-0.0420656124672994d0*x3447 + 0.00631969911517506d0
     &*x3448)
      x3468 = -x1461*x3467 + 0.0225351495360533d0*x1466*x1508 + 0.002629
     &10077920621d0* x3453
      x3469 = x1497*x3347
      x3470 = rho_s*x3325*x3444
      x3471 = x3349*x3352
      x3472 = rho_s*x3471
      x3473 = 0.00225703539827681d0*x3472
      x3474 = x1437*x3473
      x3475 = x1432*x1491
      x3476 = x1378*x3472
      x3477 = x186*x3321*x3444
      x3478 = x186*x3352
      x3479 = x3325*x3478
      x3480 = -0.00451407079655361d0*x1437*x3479
      x3481 = x1463*x1491
      x3482 = x1477**2*x3180
      x3483 = x1488*x3184*x3379
      x3484 = 1.54726434806315d-6*x3187
      x3485 = x3358*(-x133*x1474*x1559*x3394 + x134*x1475*x3392 + x1473*
     &x1475 + x3189* x3484 - x3393)
      x3486 = x3361*(x1484*x3197*x3387 - x1484*x3202*x3427 + x2924*x3425
     & + x3198*x3484 - x3388)
      x3487 = x2*(x3482 + x3483 + x3485 + x3486)
      x3488 = x1378*x3479
      x3489 = 0.67024923187228d0*rho_c**(-3.33333333333333d0)*x245
      x3490 = x1500*x2897
      x3491 = 3.0d0*x185
      x3492 = x198*x3
      x3493 = x2939*x3492
      x3494 = x1363*x3493
      x3495 = x3019*(-x3043 + x3085) - x3024
      x3496 = x1362*x1501
      x3497 = x1379 + x3496
      x3498 = -x1505 + x1506 + x1510
      x3499 = x1513*x571
      x3500 = x1437*x199
      x3501 = x203*(2.8d0*x3496 + x3500)
      x3502 = 6.0d0*x180
      x3503 = x185*x200
      x3504 = x1921 + x3*x3045
      x3505 = x3019*(-0.0875d0*x3043 + x3504) - 0.0875d0*x3024
      x3506 = -x3067 + x3071
      x3507 = x1518*x248
      x3508 = x1093*x3492
      x3509 = x1149*x199
      x3510 = x712*(-x1086 + x1092)
      x3511 = x3509 - x3510
      x3512 = -x1083*x1500 - x210*x3508 + x3511
      x3513 = x1504*x238
      x3514 = 1.5d0*x1514
      x3515 = 6.103515625d-5*x1519*x247
      x3516 = x1367*x1442
      x3517 = x1367*x3220
      x3518 = x199*x98
      x3519 = gamma_cc*x185*x199
      x3520 = x1386*x3518 - x3027*x3519
      x3521 = 64.0d0*x566
      x3522 = x1924*x3521
      x3523 = x1512*x180*x74
      x3524 = -12.6d0*x1447
      x3525 = 0.021875d0*x1443
      x3526 = 0.175d0*x185
      x3527 = x2965*x3090
      x3528 = 0.0109375d0*x3527
      x3529 = x3518*(0.0210328062336497d0*x1384 - x3526) - x3519*x3528
      x3530 = 32.0d0*x203
      x3531 = x1513*x180
      x3532 = rho_c*x1517
      x3533 = gamma_cc*x236*x3532
      x3534 = x1504*x3346*x3521
      x3535 = 0.00751171651201775d0*x1439
      x3536 = x1437*x3535
      x3537 = 32.0d0*x1162
      x3538 = 96.0d0*x1162
      x3539 = x1518*x3532
      x3540 = 0.0150234330240355d0*x1466
      x3541 = -x1437*x3540
      x3542 = 16.8d0*x230
      x3543 = 3.0d0*x3500
      x3544 = x1239*x1513
      x3545 = 3.814697265625d-6*x247*x3*x74/tau_c**4
      x3546 = x143*x229 + x205*x41
      x3547 = 0.5d0*x3546
      x3548 = x204*x3547
      x3549 = x230 - x3548
      x3550 = x204*x3546
      x3551 = 1.4d0*x3550
      x3552 = x233 - x3551
      x3553 = x238*x3552 + 32768.0d0
      x3554 = x243 - x3549*x3553
      x3555 = x248*x3554
      x3556 = x526*x709
      x3557 = 0.0625d0*x241
      x3558 = x183*x2880
      x3559 = 0.0098758948512826d0*x242
      x3560 = 1.0d0*x240
      x3561 = x3559 - x3560
      x3562 = x3557 + x3558 + x3561
      x3563 = x3556 + x3562 + x470
      x3564 = x204*x3563
      x3565 = x3546*x534
      x3566 = -x3547*x564 - x3564 + 1.5d0*x3565 + x714
      x3567 = -x1919*x3546 - x2472*x3552 - 2.8d0*x3564 + 4.2d0*x3565 + x
     &725
      x3568 = x238*x3567
      x3569 = x247*(-x3549*x3568 - x3553*x3566 + x741)
      x3570 = x247*x3554
      x3571 = x1062 + x3551
      x3572 = x238*x3571
      x3573 = -x3572 + 32768.0d0
      x3574 = x1149*x526
      x3575 = -x3557 - x3558 - x3559 + x3560 + x3574 + x470
      x3576 = x204*x3575
      x3577 = x1094*x3547
      x3578 = x1152 + x3576 + x3577
      x3579 = 0.04375d0*x3576
      x3580 = 0.021875d0*x3546
      x3581 = x1094*x3580
      x3582 = x1158 + x3579 + x3581
      x3583 = x1160 + x3548
      x3584 = x238*x3583
      x3585 = 64.0d0*x3584
      x3586 = x1202 - x3573*x3578 + x3582*x3585
      x3587 = x248*x3586
      x3588 = -x1360*x41
      x3589 = x1378*x143 + x3588
      x3590 = x1359*x3589
      x3591 = 0.25d0*x1342
      x3592 = x143*x7
      x3593 = 0.0625d0*x1177
      x3594 = x1328*x3593
      x3595 = x3591 + x3592*x3594
      x3596 = x2934*x3595
      x3597 = x1370*x3589
      x3598 = x1374*x3597
      x3599 = x1366*x3040
      x3600 = 0.120187464192284d0*x3599
      x3601 = x3598*x3600
      x3602 = x1394 - 1.0d0*x3590 - x3596 + x3601
      x3603 = 0.04375d0*x3590
      x3604 = 0.04375d0*x3596
      x3605 = 0.00525820155841243d0*x3599
      x3606 = x3598*x3605
      x3607 = gamma_cc*(x1406 - x3603 - x3604 + x3606)
      x3608 = x1396 - 0.065625d0*x3550 + x3607
      x3609 = x1410*x3608
      x3610 = x247*(x1331 - x3573*x3602 + x3583*x3609)
      x3611 = x1342*x3339
      x3612 = x1465*x3597
      x3613 = x1441 + x3611 - 0.00525820155841243d0*x3612
      x3614 = -x1447
      x3615 = x204*x3589
      x3616 = x3614 + 0.5d0*x3615
      x3617 = x238*x3616
      x3618 = x3613*x3617
      x3619 = x3572 - 32768.0d0
      x3620 = 0.00390625d0*x1442
      x3621 = x1342*x3620
      x3622 = x1378*x1440
      x3623 = -0.00187792912800444d0*x3612 + x3621 + x3622
      x3624 = x3619*x3623
      x3625 = x3618 + x3624
      x3626 = x1479*x3339
      x3627 = x3410*x3589
      x3628 = x1468 + x3626 + 0.0105164031168249d0*x3627
      x3629 = x3617*x3628
      x3630 = x1479*x3620
      x3631 = -x1378*x1467
      x3632 = 0.00375585825600888d0*x3627 + x3630 + x3631
      x3633 = x3619*x3632
      x3634 = x3629 + x3633
      x3635 = x1500*x3546
      x3636 = x1503 - 0.5d0*x3635
      x3637 = 2.1d0*x3550
      x3638 = x212*(x1509 - 1.4d0*x3635)
      x3639 = x1507 + x3637 - x3638
      x3640 = x1514*x3639
      x3641 = x3573*x3636 + x3583*x3640
      x3642 = 0.25d0*x1907
      x3643 = x143*x3642
      x3644 = -4.0d0*x570 + x579 + x583 + x591 + x595 + x603 + x708 - 2.
     &0d0*x736
      x3645 = x204*(x1586 + x3643 + x3644)
      x3646 = x3563*x534
      x3647 = x1672*x3546
      x3648 = x1694*x3546
      x3649 = x1668*x3648
      x3650 = x1920*x3546
      x3651 = 8.0d0*x238
      x3652 = x2464*x526
      x3653 = 0.5d0*x618
      x3654 = 1.0d0*x736
      x3655 = x1709*x568
      x3656 = 0.064719778448038d0*x739
      x3657 = 2.0d0*x570
      x3658 = 0.000332743037440633d0*x3*x928
      x3659 = 0.0625d0*x699
      x3660 = x3659*x695*x7
      x3661 = 1.28244070062102d0*x652*x702
      x3662 = 0.0253304990246769d0*x2335
      x3663 = 0.125d0*x707
      x3664 = x3663*x704
      x3665 = 0.29241116987276d0*x2333
      x3666 = -x1098*x1709
      x3667 = 0.0253304990246769d0*x2617*x3
      x3668 = x1098*x569
      x3669 = 1.28244070062102d0*x1114*x567 + x1137*x3659 + x2483*x3663 
     &+ x3666 + x3667 + 2.0d0*x3668
      x3670 = x204*(x1586 + x3652 - x3653 + x3654 - x3655 - x3656 + x365
     &7 + x3658 - x3660 - x3661 - x3662 - x3664 - x3665 + x3669)
      x3671 = x3575*x534
      x3672 = x1094*x3563
      x3673 = x3575*x564
      x3674 = x2379*x3648
      x3675 = x2470 - x3579 - x3581
      x3676 = 4.0d0*x3668
      x3677 = x204*(0.0253304990246769d0*x1107*x1164 - x1139 + x1141 - x
     &1145 - x1147 + x1586 + x2787*x526 - x3667 - x3676)
      x3678 = x1177*(x1329*x705 + x1390*(-x2984 - x2985 - x2987 - x2988 
     &+ x2990 + x2992 - x2994 + x2995 + x2996 + x2998 + x2999 - x3002) -
     & x2982* x3064*x3145)
      x3679 = x1330*x3593
      x3680 = 8.36096418563526d-6*x71
      x3681 = x129*x134*x3150
      x3682 = 0.25d0*x40
      x3683 = x3121*x3682 + x453*(x1332*x3680 - x1338*x3681 + x3681)
      x3684 = x204*(x3678*x526 + x3679 + x3683)
      x3685 = x3065*x526
      x3686 = -x183*x2818 + x3561
      x3687 = 0.0625d0*x1434 + x3686
      x3688 = x3110 + x3685 + x3687
      x3689 = x1359*x3688
      x3690 = x2936*x3589
      x3691 = x3595*x40
      x3692 = x2939*x3589
      x3693 = 1.0d0*x3038
      x3694 = x3035*x3691
      x3695 = 0.5d0*x203
      x3696 = x2955*x3589
      x3697 = x2946*x3402
      x3698 = x186*(x2953 + x3697)
      x3699 = x3696*x3698
      x3700 = x1886*x3040
      x3701 = x3598*x3700
      x3702 = x2960*x3688
      x3703 = x1680*x2937
      x3704 = x199*x3065
      x3705 = x1392*x1913
      x3706 = x1387*x3704 - x214*x3678 + x3086*x3705 - x3159*(-x2937*x31
     &58 + x3008* x3703 + x3018)
      x3707 = 0.04375d0*x3038
      x3708 = 0.021875d0*x203
      x3709 = x1403*x3704 + x3068*x3705 - x3159*(-x1401*x2938 + x3044*x3
     &703 + x3047) - x3678*x718
      x3710 = -x3030
      x3711 = 1.4d0*x3615 + x3710
      x3712 = x3032*x3711
      x3713 = gamma_cc*(x3052 + x3603 + x3604 - x3606)
      x3714 = x3589*x534
      x3715 = 12.6d0*x3714
      x3716 = x204*x3688
      x3717 = 8.4d0*x3716
      x3718 = x3692*x563
      x3719 = 4.2d0*x3718
      x3720 = x3073 + x3590 + x3596 - x3601
      x3721 = x3262 + 4.2d0*x3714 - 2.8d0*x3716 + 1.4d0*x3718
      x3722 = x238*(x2472*x3711 + x3721)
      x3723 = -x3081
      x3724 = 0.065625d0*x3615 + x3713 + x3723
      x3725 = 1.5d0*x3714
      x3726 = 0.5d0*x3718
      x3727 = x3088 + x3716 - x3725 - x3726
      x3728 = x204*(x3154*x526 - x3679 + x3683)
      x3729 = x3562 - x3574 - x464 - x466 + x468
      x3730 = x1359*x3729
      x3731 = x1093*x3589
      x3732 = x1094*x3691
      x3733 = x3134*x3597
      x3734 = x2960*x3729
      x3735 = x3582*x3583
      x3736 = x1358*x3546
      x3737 = x143*x696
      x3738 = 0.0625d0*x7
      x3739 = x3143*x3738
      x3740 = x3591 + x3737*x3739
      x3741 = x3740*x40
      x3742 = x1110*x1380*x200
      x3743 = 0.0625d0*x622
      x3744 = x1789*x3743
      x3745 = x143*x3234
      x3746 = x3743*x696
      x3747 = x2*x2934
      x3748 = x3747*(x143*x3236*x3744 + 0.25d0*x3182 + 0.25d0*x3196 + x3
     &745*x3746)
      x3749 = x3351*x3546
      x3750 = x3217*x3223
      x3751 = x1370*x3546
      x3752 = x202*x3016*x558
      x3753 = x3751*x3752
      x3754 = x3210*x3740
      x3755 = 0.175d0*x1359
      x3756 = 0.00157992477879376d0*x3352
      x3757 = x3218*x3756
      x3758 = 21.447975419913d0*x744
      x3759 = x1465*x3751
      x3760 = x3258 + x3611 - 0.00525820155841243d0*x3759
      x3761 = x2474*(x3260 - x3716 + x3725 + x3726)
      x3762 = 0.00375585825600888d0*x1438*x3264
      x3763 = x3621 - 0.00187792912800444d0*x3759 + x3762
      x3764 = x2474*(x2472*x3571 + x3721)
      x3765 = x3409*x3751
      x3766 = x201*x453
      x3767 = x3766*(x1340 + x3266)
      x3768 = x1370*x3563
      x3769 = x40*x533
      x3770 = x1342*x3769
      x3771 = x40*x562
      x3772 = x3771*x561
      x3773 = 0.021875d0*x3281
      x3774 = 0.0109375d0*x1342
      x3775 = x204*x3729
      x3776 = x1094*x3589
      x3777 = x238*(x3304 - x3775 + 0.5d0*x3776)
      x3778 = x238*(x3306 - 0.04375d0*x3775 + 0.021875d0*x3776)
      x3779 = 64.0d0*x3623
      x3780 = x3310*x3696
      x3781 = x1370*x3289
      x3782 = x3766*(x1340 + x3120)
      x3783 = x1093*x3771
      x3784 = 0.00375585825600888d0*x1370
      x3785 = x1465*x3784
      x3786 = 0.00390625d0*x3783
      x3787 = x1438*x3160
      x3788 = x3318 + x3784*x3787
      x3789 = 64.0d0*x3319
      x3790 = x3602*x3789
      x3791 = x238*x3763
      x3792 = 4096.0d0*x3608
      x3793 = x1342*x1442
      x3794 = gamma_cc*x1432
      x3795 = x3740*x3794
      x3796 = x3222*x3322
      x3797 = 0.25d0*rho_c
      x3798 = 0.25d0*rho_s
      x3799 = x3172*x3766*(x3797 + x3798)
      x3800 = x3799*(x1341*x3179 + x3195)
      x3801 = x1342*x3220
      x3802 = -64.0d0*x3330 + 0.030046866048071d0*x3336
      x3803 = x3322*x3756
      x3804 = x248*x3344
      x3805 = x3347*x3613
      x3806 = x2*x3355
      x3807 = x1342*x1432
      x3808 = x91*x92
      x3809 = x1436 + 0.35d0*x3808
      x3810 = x3354*x3809
      x3811 = x143*x98 + x3588
      x3812 = x3471*x3811
      x3813 = x2*x3359
      x3814 = x3354*x98
      x3815 = x3410*x3546
      x3816 = x3370 + x3626 + 0.0105164031168249d0*x3815
      x3817 = rho_s*x1465
      x3818 = -0.00751171651201775d0*x3264*x3817
      x3819 = x3630 + 0.00375585825600888d0*x3815 + x3818
      x3820 = x3406*x3597
      x3821 = x3766*(x1477 + x3395)
      x3822 = 0.25d0*x3821
      x3823 = x1460*x3599
      x3824 = x3589*x3823
      x3825 = x1479*x3769
      x3826 = x3410*x3688
      x3827 = x1479*x3771
      x3828 = 0.0109375d0*x3821
      x3829 = 0.0109375d0*x2939
      x3830 = x3431*x3696
      x3831 = x3410*x3729
      x3832 = 0.0109375d0*x3783
      x3833 = x1430*x3367
      x3834 = x238*x3819
      x3835 = x1442*x1479
      x3836 = x3445*x3589
      x3837 = x3224*x3352
      x3838 = 0.0722251327448578d0*x3837
      x3839 = x3466*x3595
      x3840 = x3799*(x1341*x1477 + x3457)
      x3841 = x1479*x3220
      x3842 = -64.0d0*x3450 + 0.060093732096142d0*x3453
      x3843 = 0.0210328062336497d0*x3627
      x3844 = x3244*x3351
      x3845 = x3224*x3844
      x3846 = x3347*x3628
      x3847 = x1432*x1479
      x3848 = x1342*x1463
      x3849 = x2*x3339
      x3850 = rho_s*x3811
      x3851 = x3325*x3844
      x3852 = x2*x3620
      x3853 = x3326*x3352
      x3854 = x2*x3482
      x3855 = x1463*x1479
      x3856 = x3321*x3478
      x3857 = x3589*x3856
      x3858 = -0.00451407079655361d0*x3488
      x3859 = x3023*x526
      x3860 = x3686 + 0.0625d0*x3808
      x3861 = x1500*(x2893 + x3859 + x3860 + x469)
      x3862 = x3503*x3811
      x3863 = x3493*x3811
      x3864 = x1500*x3589
      x3865 = x1379 - 0.5d0*x3864
      x3866 = x1506 - x3637 + x3638
      x3867 = x203*(x3500 - 1.4d0*x3864)
      x3868 = x1500*x3575
      x3869 = 0.04375d0*x3508
      x3870 = 0.0875d0*x3509 - 0.0875d0*x3510
      x3871 = x238*x3636
      x3872 = x201*x3811
      x3873 = x3337*(0.0625d0*x3091*x3592 + x3591)
      x3874 = x1370*x558
      x3875 = x3811*x3874
      x3876 = 0.120187464192284d0*x3016
      x3877 = x1399*x3015
      x3878 = 2.0d0*x3346*x566
      x3879 = x3636*x3878
      x3880 = x1378*x3535
      x3881 = 3.0d0*x3584
      x3882 = 2.0d0*x3639
      x3883 = rho_c*x3299
      x3884 = x1518*x3883
      x3885 = -x1378*x3540
      x3886 = x206*x41 + x229*x99
      x3887 = 0.5d0*x3886
      x3888 = x204*x3887
      x3889 = x230 - x3888
      x3890 = x204*x3886
      x3891 = 1.4d0*x3890
      x3892 = x233 - x3891
      x3893 = x238*x3892 + 32768.0d0
      x3894 = x243 - x3889*x3893
      x3895 = x248*x3894
      x3896 = x463*x709
      x3897 = x3562 + x3896 + x530
      x3898 = x204*x3897
      x3899 = x3886*x534
      x3900 = -x3887*x564 - x3898 + 1.5d0*x3899 + x714
      x3901 = -x1919*x3886 - x2472*x3892 - 2.8d0*x3898 + 4.2d0*x3899 + x
     &725
      x3902 = x238*x3901
      x3903 = x247*(-x3889*x3902 - x3893*x3900 + x741)
      x3904 = x247*x3894
      x3905 = x1062 + x3891
      x3906 = x238*x3905
      x3907 = -x3906 + 32768.0d0
      x3908 = x1082 + x1149*x463 + x3562
      x3909 = x204*x3908
      x3910 = x1094*x3887
      x3911 = x1152 + x3909 + x3910
      x3912 = 0.04375d0*x3909
      x3913 = 0.021875d0*x3886
      x3914 = x1094*x3913
      x3915 = x1158 + x3912 + x3914
      x3916 = x1160 + x3888
      x3917 = x238*x3916
      x3918 = 64.0d0*x3917
      x3919 = x1202 - x3907*x3911 + x3915*x3918
      x3920 = x248*x3919
      x3921 = x1361*x41
      x3922 = x1378*x99
      x3923 = -x3921 + x3922
      x3924 = x1359*x3923
      x3925 = 0.25d0*x1353
      x3926 = x3594*x7*x99 + x3925
      x3927 = x2934*x3926
      x3928 = x2960*x3923
      x3929 = x3600*x3928
      x3930 = x1394 - 1.0d0*x3924 - x3927 + x3929
      x3931 = 0.04375d0*x3924
      x3932 = 0.04375d0*x3927
      x3933 = x3605*x3928
      x3934 = gamma_cc*(x1406 - x3931 - x3932 + x3933)
      x3935 = x1396 - 0.065625d0*x3890 + x3934
      x3936 = x1410*x3935
      x3937 = x247*(x1331 - x3907*x3930 + x3916*x3936)
      x3938 = x1353*x3339
      x3939 = x1370*x1465
      x3940 = x3923*x3939
      x3941 = x1441 + x3938 - 0.00525820155841243d0*x3940
      x3942 = x204*x3923
      x3943 = x3614 + 0.5d0*x3942
      x3944 = x238*x3943
      x3945 = x3941*x3944
      x3946 = x3906 - 32768.0d0
      x3947 = x1353*x3620
      x3948 = x3622 - 0.00187792912800444d0*x3940 + x3947
      x3949 = x3946*x3948
      x3950 = x3945 + x3949
      x3951 = x3410*x3923
      x3952 = x1468 - x1490*x3339 + 0.0105164031168249d0*x3951
      x3953 = x3944*x3952
      x3954 = -x1490*x3620 + x3631 + 0.00375585825600888d0*x3951
      x3955 = x3946*x3954
      x3956 = x3953 + x3955
      x3957 = x1500*x3886
      x3958 = x1503 - 0.5d0*x3957
      x3959 = 2.1d0*x3890
      x3960 = x212*(x1509 - 1.4d0*x3957)
      x3961 = x1507 + x3959 - x3960
      x3962 = x1514*x3961
      x3963 = x3907*x3958 + x3916*x3962
      x3964 = x3642*x99
      x3965 = x204*(x1665 + x3644 + x3964)
      x3966 = x3897*x534
      x3967 = x1672*x3886
      x3968 = x1694*x3886
      x3969 = x1668*x3968
      x3970 = x1920*x3886
      x3971 = x2464*x463
      x3972 = x204*(x2361 + x3653 - x3654 + x3655 + x3656 - x3657 - x365
     &8 + x3660 + x3661 + x3662 + x3664 + x3665 + x3669 + x3971)
      x3973 = x3908*x534
      x3974 = x1094*x3897
      x3975 = x3908*x564
      x3976 = x2379*x3968
      x3977 = x2470 - x3912 - x3914
      x3978 = x1148 + x3676
      x3979 = x204*(x2712 + x2787*x463 + x3978)
      x3980 = x166*x171*x3150
      x3981 = x472*(x1343*x3680 - x1349*x3980 + x3980)
      x3982 = x204*(x2932*x3682*(x2916 + x3268) + x3678*x463 + x3679 + x
     &3981)
      x3983 = x3065*x463
      x3984 = x3056 + x3687 + x3983 + x529
      x3985 = x1359*x3984
      x3986 = x2936*x3923
      x3987 = x3926*x40
      x3988 = x2939*x3923
      x3989 = x3035*x3987
      x3990 = x2955*x3923
      x3991 = x3698*x3990
      x3992 = x3700*x3928
      x3993 = x2960*x3984
      x3994 = x3710 + 1.4d0*x3942
      x3995 = x3032*x3994
      x3996 = gamma_cc*(x3052 + x3931 + x3932 - x3933)
      x3997 = x3923*x534
      x3998 = 12.6d0*x3997
      x3999 = x204*x3984
      x4000 = 8.4d0*x3999
      x4001 = x3035*x3923
      x4002 = 4.2d0*x4001
      x4003 = x3073 + x3924 + x3927 - x3929
      x4004 = x3262 + 4.2d0*x3997 - 2.8d0*x3999 + 1.4d0*x4001
      x4005 = x238*(x2472*x3994 + x4004)
      x4006 = x3723 + 0.065625d0*x3942 + x3996
      x4007 = 1.5d0*x3997
      x4008 = 0.5d0*x3923
      x4009 = x3035*x4008
      x4010 = x3088 + x3999 - x4007 - x4009
      x4011 = x204*(x3125*x3682 + x3154*x463 + x3679 - x3981)
      x4012 = x3115 + x3160*x463 + x3687
      x4013 = x1359*x4012
      x4014 = x1093*x3923
      x4015 = x1094*x3987
      x4016 = x2960*x4012
      x4017 = 0.021875d0*x3923
      x4018 = x3915*x3916
      x4019 = x1358*x3886
      x4020 = x3739*x696*x99 + x3925
      x4021 = x40*x4020
      x4022 = x3236*x99
      x4023 = x3234*x99
      x4024 = x3747*(0.25d0*x3186 + 0.25d0*x3206 + x3744*x4022 + x3746*x
     &4023)
      x4025 = x3351*x3886
      x4026 = x1370*x3886
      x4027 = x3752*x4026
      x4028 = x3210*x4020
      x4029 = x3886*x3939
      x4030 = x3258 + x3938 - 0.00525820155841243d0*x4029
      x4031 = x2474*(x3260 - x3999 + x4007 + x4009)
      x4032 = x3762 + x3947 - 0.00187792912800444d0*x4029
      x4033 = x2474*(x2472*x3905 + x4004)
      x4034 = x3409*x4026
      x4035 = x201*x472
      x4036 = x4035*(x1351 + x3269)
      x4037 = x1353*x3769
      x4038 = x204*x4012
      x4039 = x238*(x1094*x4008 + x3304 + x4038)
      x4040 = x238*(x1094*x4017 + x3306 + 0.04375d0*x4038)
      x4041 = 64.0d0*x3948
      x4042 = x3310*x3990
      x4043 = x4035*(x3124 + x3183)
      x4044 = x3789*x3930
      x4045 = x238*x4032
      x4046 = 4096.0d0*x3935
      x4047 = x1353*x1442
      x4048 = x3794*x4020
      x4049 = x3172*x4035*(x3797 - x3798)
      x4050 = x4049*(x1352*x3183 + x3205)
      x4051 = x1353*x3220
      x4052 = x3347*x3941
      x4053 = x2*x3356
      x4054 = x1353*x1432
      x4055 = x3921 - x3922
      x4056 = x3471*x4055
      x4057 = x3410*x3886
      x4058 = x3339*x3372 + x3370 + 0.0105164031168249d0*x4057
      x4059 = x3372*x3620 + x3818 + 0.00375585825600888d0*x4057
      x4060 = x4035*(x3371 + x3390)
      x4061 = x3823*x3923
      x4062 = x1490*x3769
      x4063 = x3410*x3984
      x4064 = x1490*x3771
      x4065 = x3431*x3990
      x4066 = x4035*(x1488 + x3428)
      x4067 = x3410*x4012
      x4068 = x238*x4059
      x4069 = x1442*x1490
      x4070 = x3445*x3923
      x4071 = x3466*x3926
      x4072 = x4049*(x1352*x3379 + x3459)
      x4073 = x1490*x3220
      x4074 = 0.021875d0*x4069
      x4075 = 0.0210328062336497d0*x3951
      x4076 = x3347*x3952
      x4077 = x1353*x1463
      x4078 = x1432*x1490
      x4079 = rho_s*x3923
      x4080 = x1463*x1490
      x4081 = x3856*x3923
      x4082 = x3023*x463
      x4083 = x1500*(x2895 + x3860 + x4082 + x529)
      x4084 = x3503*x4055
      x4085 = x3493*x4055
      x4086 = x1500*x3923
      x4087 = x1379 - 0.5d0*x4086
      x4088 = x1506 - x3959 + x3960
      x4089 = x203*(x3500 - 1.4d0*x4086)
      x4090 = x1500*x3908
      x4091 = x238*x3958
      x4092 = x201*x4055
      x4093 = x3337*(x3091*x3738*x99 + x3925)
      x4094 = x3874*x4055
      x4095 = x3878*x3958
      x4096 = 3.0d0*x3917
      x4097 = 2.0d0*x3961
      x4098 = x200*x3562
      x4099 = x1912*x4098
      x4100 = x230 - x4099
      x4101 = 5.6d0*x4098
      x4102 = x212*x4101
      x4103 = x233 - x4102
      x4104 = x238*x4103 + 32768.0d0
      x4105 = x243 - x4100*x4104
      x4106 = x248*x4105
      x4107 = -2.0d0*x240 + x61
      x4108 = x228 + x3556 + x3896 + x4107 + 1.0d0*x50
      x4109 = x204*x4108
      x4110 = x204*x3562
      x4111 = x1668*x1912
      x4112 = -x3562*x4111 - x4109 + 6.0d0*x4110 + x714
      x4113 = x1913*x3562
      x4114 = -x2472*x4103 - 2.8d0*x4109 + 16.8d0*x4110 - 5.6d0*x4113*x5
     &61 + x725
      x4115 = x238*x4114
      x4116 = x247*(-x4100*x4115 - x4104*x4112 + x741)
      x4117 = x247*x4105
      x4118 = x238*(x1062 + x4102)
      x4119 = -x4118 + 32768.0d0
      x4120 = x200*x3669
      x4121 = x1912*x4120
      x4122 = x1093*x3562
      x4123 = x1912*x198
      x4124 = x4122*x4123
      x4125 = x1152 + x4121 + x4124
      x4126 = x2791*x4120
      x4127 = x1401*x3562
      x4128 = x203*x4127
      x4129 = x2379*x4128
      x4130 = x1158 + x4126 + x4129
      x4131 = x1160 + x4099
      x4132 = x238*x4131
      x4133 = 64.0d0*x4132
      x4134 = x1202 - x4119*x4125 + x4130*x4133
      x4135 = x248*x4134
      x4136 = x200*x3687
      x4137 = x1380*x4136
      x4138 = x1725*x4137
      x4139 = x3409*x3687
      x4140 = x2960*x4139
      x4141 = 0.480749856769136d0*x4140
      x4142 = x200*x212
      x4143 = x1392*x4142
      x4144 = x1394 - x4138 + x4141 - x4143
      x4145 = x212*x4098
      x4146 = x3526*x4137
      x4147 = 0.0210328062336497d0*x4140
      x4148 = x3164*x4142
      x4149 = gamma_cc*(x1406 - x4146 + x4147 - x4148)
      x4150 = x1396 - 0.2625d0*x4145 + x4149
      x4151 = x1410*x4150
      x4152 = x247*(x1331 - x4119*x4144 + x4131*x4151)
      x4153 = x1503 - 2.0d0*x4098
      x4154 = 8.4d0*x4145
      x4155 = x212*(x1509 - x4101)
      x4156 = x1507 + x4154 - x4155
      x4157 = x1514*x4156
      x4158 = x4119*x4153 + x4131*x4157
      x4159 = x1694*x4098
      x4160 = x204*(x1768*x704 + x2179*x3 + x2332*x702 + x2336 + x3643 +
     & x3964 - 8.0d0 *x570 + 2.0d0*x618 + 0.25d0*x701 - 4.0d0*x736)
      x4161 = x4108*x534
      x4162 = x1690*x4113
      x4163 = x3562*x564
      x4164 = x204*(x3652 + x3971 + x3978)
      x4165 = x204*x3669
      x4166 = x1094*x4108
      x4167 = x3562*x4123
      x4168 = x2470 - x4126 - x4129
      x4169 = x198*x4128
      x4170 = x1671*(x1144*x2743 + 6.0d0*x1673*x2880 + 0.427480233540341
     &d0*x1698*x2877 - x2732*x566*x963 + 10.2595256049682d0*x2733 - 0.04
     &16666666666667d0*x2738 + 0.166666666666667d0*x2739*x2741 - 0.0625d
     &0*x2742 + 0.0625d0*x2786 - 0.118510738215391d0*x2795 + 0.008443499
     &67489228d0*x2817)
      x4171 = x2379*x3526
      x4172 = x4118 - 32768.0d0
      x4173 = x4107 - 1.0d0*x96
      x4174 = x3859 + x4082 + x4173 + x94
      x4175 = x1359*x4174
      x4176 = x1359*x3860
      x4177 = x204*x2964*(x3004*x463 + x3004*x526 + x3026)
      x4178 = x1380*x1725
      x4179 = x2939*x3860
      x4180 = x198*x4179
      x4181 = x2959*x3860
      x4182 = x187*x2954*x4181
      x4183 = x204*x3527
      x4184 = x2960*x4174
      x4185 = x186*x2416*x3040*x4181
      x4186 = x1912*x4136 + x3614
      x4187 = x204*x3860
      x4188 = 50.4d0*x4187
      x4189 = x1913*x4179
      x4190 = 16.8d0*x4189
      x4191 = x212*x3019
      x4192 = 192.0d0*x3049*x4191
      x4193 = 5.6d0*x4136
      x4194 = x212*x4193 + x3710
      x4195 = x3032*x4194
      x4196 = gamma_cc*(x3052 + x4146 - x4147 + x4148)
      x4197 = x3023*x718
      x4198 = 192.0d0*x4197
      x4199 = x204*x4174
      x4200 = 8.4d0*x4199
      x4201 = x3073 + x4138 - x4141 + x4143
      x4202 = x238*(x2472*x4194 - x3019*x3049*x724 + 16.8d0*x4187 + 5.6d
     &0*x4189 + 64.0d0*x4197 - 2.8d0*x4199)
      x4203 = x212*x4136
      x4204 = x3723 + x4196 + 0.2625d0*x4203
      x4205 = x214*x3023 - x3025*x4191 + x4123*x4179 + 6.0d0*x4187 - x41
     &99
      x4206 = -2.0d0*x1098*x798 - 0.0253304990246769d0*x1165 - 0.125d0*x
     &1198*x1199 + 1.28244070062102d0*x1201 + 0.0625d0*x2822 + x3666
      x4207 = x1380*x4206
      x4208 = x1177*x1671*(-x1199*x1328*x1389 - 0.0625d0*x3147 + 0.5d0*x
     &3153)
      x4209 = x1911*x3687
      x4210 = x1460*x203*x3133
      x4211 = 0.480749856769136d0*x3409
      x4212 = x1401*x3687
      x4213 = 0.0210328062336497d0*x4206
      x4214 = x4130*x4131
      x4215 = x202*x4098
      x4216 = x3238*x3239*x3503
      x4217 = x203*x3562*x3837*x46
      x4218 = x1370*x189*x202*x3409*x3562
      x4219 = x3747*(x3236*x3737 + x3745 + x4022*x696 + x4023)
      x4220 = x1375*x3335*x548
      x4221 = x2947*x70 + x3276
      x4222 = x1377 + x3685 + x3983 + x4173
      x4223 = x238*x4186
      x4224 = 64.0d0*x4223
      x4225 = x1370*x1459
      x4226 = x3310*x3687
      x4227 = x3324*x3860
      x4228 = x3341*x3860
      x4229 = x1369*x1370
      x4230 = x3299*x3346*x3352
      x4231 = x3409*x4222
      x4232 = x3817*x4206
      x4233 = rho_s*x3160
      x4234 = x3447*x3687
      x4235 = x3448*x3687
      x4236 = x1500*x3687
      x4237 = x2939*x3687
      x4238 = x1500*x4222
      x4239 = x1379 - 2.0d0*x4136
      x4240 = x1506 - x4154 + x4155
      x4241 = x203*(x3500 - x4193)
      x4242 = 0.175d0*x198
      x4243 = x238*x4153
      x4244 = x3503*x3860
      x4245 = gamma_cc*x3503
      x4246 = x1370*x189*x3860
      PBE_value = x98
      PBE_beta_value = x142*x41
      PBE_alpha_value = x179*x41
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x249 + x407)
      d1E(1) = x1060*(x1059 - x742 - x745*x746) + x1061 - x249
      d1E(2) = -rho_c*(-x1204 + x1323)
      d1E(3) = x1430*(-0.000244140625d0*x1412 + x1429)
      d1E(4) = x1457*x1459
      d1E(5) = x1459*x1499
      d1E(6) = -x1516*x1520
      d2E(1) = x2329*(-x1937*(-8.0d0*x1923*(-0.7d0*x1666 + 4.2d0*x1667 -
     & x1672*x717 + 4.2d0*x1696 + x185*x716 - x1918*(x1916*x210 + x1917 
     &- 0.04375d0* x532 - 0.04375d0*x565) - x1919*x531 - 0.7d0*x1920*x21
     &0 + x1922) - x1925*x715*x726 + x1936 - x239*(-2.0d0*x1666 + 12.0d0
     &*x1667 - x1670*x531 - 24.0d0*x1672*x210 - x1690*x1693*x210 + 12.0d
     &0*x1696 + x1915)) + x1938*x746 - x1939*x742 - x1940*x746 + x2328) 
     &- x2330 *x746 + x2337 - 9.5367431640625d-7*x742
      d2E(2) = x1204 + x2706*(x1203*x2488 - x1937*(-64.0d0*x1923*(-0.043
     &75d0*x2362 + 0.13125d0*x2363 - 0.04375d0*x2364 - 0.04375d0*x2365 -
     & 0.04375d0* x2378 + 0.13125d0*x2380 - x2471*x2472 + x2473) - x239*
     &(-x2362 + 3.0d0*x2363 - x2364 - x2365 - x2378 + 3.0d0*x2380 + x246
     &8) - x2471*x2474*x715 + x2487 - x727*(-x1084 - x1095 + x2469)) + x
     &2705 ) + x2708
      d2E(3) = -rho_c*(-x248*(-x1065*(x1083*x2714 + x2713 + x2720 + x278
     &9) + x1153* x1159*x2793 + x1163*(x1083*x2790 + 0.04375d0*x2713 + 0
     &.04375d0* x2720 + x2792) + x2824) + x2890)
      d2E(4) = 4.76837158203125d-7*rho_c*(-x1412*x745 + x247*(x1409*x144
     &8*(-x3033 + x3051*(x203*x3034 - 0.0875d0*x2898 - 0.0109375d0*x2935
     & - x2936* x3037 - x2941*x3039 + 0.04375d0*x2957 + x2961*x3041 + x2
     &963*x3042 + x3035*x3036 + x3050) + x3053*x3054 - x3059 - x3061 - x
     &3063 + x3072) - x1410*x3082*x3089 - x1451*(x1355*x2940*x563 - 6.0d
     &0* x1363*x2936 + x203*x2900 - 2.0d0*x2898 - 0.25d0*x2935 - x2941* 
     &x2942 + x2957 + x2958*x2961 + x2962*x2963 + x3029) + x3074*x3080 +
     & x3092) + x3109) - 4.76837158203125d-7*x1412 + 4.76837158203125d-7
     &*x2891
      d2E(5) = -x1060*(-x247*(x1065*(-0.25d0*x1094*x1355 + x2942*x3128 +
     & x2958*x3138 - 2.0d0*x3117 - 0.25d0*x3127 - x3136 + x3163) - x1153
     &*x1411 - x1159 *x1395*x2474 - x1163*(-x1094*x3036 + x3039*x3128 + 
     &x3041*x3138 - 0.0875d0*x3117 - 0.0109375d0*x3127 - 0.04375d0*x3136
     & + x3166) + x3167*x3168 + x3170) + x3176)
      d2E(6) = -x3256*(-x3253*(-x1065*(2.0d0*x3177 + 1.0d0*x3178 + 0.25d
     &0*x3208 - x3209 *x3210 + x3211*x3216*x3218 - 0.961499713538272d0*x
     &3221 + x3241) - x1161*x3248*(6.0d0*x1407 - 0.2625d0*x211 + x3242 -
     & x74*(0.0875d0* x3177 + 0.04375d0*x3178 + 0.0109375d0*x3208 - x321
     &0*x3243 + x3211 *x3214*x3217*x3245 - 0.0420656124672994d0*x3221 + 
     &x3247)) + x1395 *x1408*x3249 + x3252) + x3255)
      d2E(7) = x3300*(-x1450 + 0.015625d0*x1456 + x3298*(-x1065*(-x2900 
     &+ x2962*x3285 + 0.25d0*x3270 + 0.25d0*x3272 + x3282 - x3283*x3284 
     &+ x3288) - x1163*(x2472*x3259 + x3034 - x3042*x3285 - 0.0109375d0*
     &x3270 - 0.0109375d0*x3272 - 0.04375d0*x3282 + x3284*x3289 + x3297)
     & + x1457*x745 - x3259*x3261 - x3263*x3265))
      d2E(8) = x1459*(64.0d0*x1445*x3305 - 64.0d0*x1449*(-x3137*x3289 + 
     &0.0109375d0* x3308 + 0.0109375d0*x3309 + 0.04375d0*x3311 + x3316) 
     &+ x1451*( x1453*x3160 - x3137*x3283 + 0.25d0*x3308 + 0.25d0*x3309 
     &+ x3311 + 64.0d0*x3318) + 64.0d0*x1455*x3307)
      d2E(9) = -x3345*(x1065*(0.0078125d0*x1443 + x3207*x3338 + 0.001128
     &5176991384d0* x3215*x3323 - 0.00751171651201775d0*x3257 - x3330 - 
     &0.000938964564002219d0*x3332 - 0.000938964564002219d0*x3334 + 0.00
     &046948228200111d0*x3336) - x1162*(gamma_cc*x3207*x3339 + 0.0546875
     &d0*x1443 - 0.0525820155841243d0*x3257 + x3323*x3340 - 0.0026291007
     &7920621d0*x3332 - 0.00262910077920621d0*x3334 + x3343 ) + x1408*x3
     &265 + x3259*x3320)
      d2E(10) = x3368*(x1445*x3348 + x3364*(x1437*x3354 + 0.005258201558
     &41243d0*x3331 + x3339*x3363 + x3340*x3350) - x3366*(0.072225132744
     &8578d0*x1378* x3353 + x1432*x3209 + x3216*x3350 + x3363*x3365))
      d2E(11) = x3300*(-x1494 + 0.015625d0*x1498 + x3298*(x1065*(x2940*x
     &3398 - x3365* x3397 + 0.75d0*x3377 - x3408 - 0.480749856769136d0*x
     &3411 - 1.44224957030741d0*x3412 + x3418) - x1163*(x2472*x3374 + x2
     &939* x3419 - x3339*x3397 + 0.0328125d0*x3377 - 0.04375d0*x3408 - 0
     &.0210328062336497d0*x3411 - 0.0630984187009491d0*x3412 + x3424) + 
     &x1499*x745 - x3261*x3374 - x3263*x3375))
      d2E(12) = x3443*(-x1449*(x1093*x3419 + x3339*x3429 - 0.04375d0*x34
     &32 + 0.0210328062336497d0*x3433 + x3438) + x1493*x3305 + x1497*x33
     &07 + x3442*(x3365*x3429 + x3398*x3439 - x3432 + 0.480749856769136d
     &0* x3433 + x3441))
      d2E(13) = -x3345*(-x1065*(0.0150234330240355d0*x1464 - 0.0078125d0
     &*x1492 - 0.00187792912800444d0*x3333*x3451 - x3338*x3461 - 0.00225
     &703539827681d0*x3446 - x3450 + 0.000938964564002219d0* x3452 + 0.0
     &00938964564002219d0*x3453) + x1162*(-gamma_cc*x3464 - x3243*x3466 
     &- 0.0630984187009491d0*x3369 - x3373*x3462 - 0.00631969911517506d0
     &*x3446 + 0.00262910077920621d0*x3452 - x3463 + x3465 + x3468) + x1
     &408*x3375 + x3320*x3374)
      d2E(14) = x3443*(x1445*x3469 - x1449*(x1463*x3243 + x3464 + 0.0063
     &1969911517506d0* x3470 + x3474 - 0.00262910077920621d0*x3475) + x1
     &493*x3348 + x3442*(x1463*x3209 + x3365*x3461 + 0.144450265489716d0
     &*x3470 - 0.060093732096142d0*x3475 + 0.144450265489716d0*x3476))
      d2E(15) = x3368*(x1493*x3469 - x3364*(-x3339*x3487 - 0.01263939823
     &03501d0*x3477 + x3480 + 0.0105164031168249d0*x3481) + x3366*(-x336
     &5*x3487 - 0.288900530979431d0*x3477 + 0.240374928384568d0*x3481 - 
     &0.288900530979431d0*x3488))
      d2E(16) = x3507*(rho_c*(-x1516*x3489 + x185*(x1448*x1514*(x3033 + 
     &x3059 + x3061 + x3063 + x3501*x3502 + x3506 + x724*(x3037*x3503 + 
     &0.0875d0*x3490 + 0.0875d0*x3494 + x3505)) - 2.0d0*x1451*(x1363*x20
     &0*x3491 + x3490 + x3494 + x3495) + 2.0d0*x3080*x3497 + x3089*x3498
     &*x3499)) - x1516*x1911)
      d2E(17) = -x3515*(0.015625d0*x1065*x3512 + 0.015625d0*x1153*x1515 
     &- x1159*x3513 - 0.04375d0*x1162*x3512 + x3167*x3514)
      d2E(18) = x3533*(x1448*x3531*(-12.6d0*x1446 + 96.0d0*x3053 - x3491
     &*x3501 + x3524 + x3530*(gamma_cc*x3525 - 0.175d0*x3516 + 0.0210328
     &062336497d0* x3517 + x3529)) - x1451*(gamma_cc*x1454 - 2.0d0*x3516
     & + 0.240374928384568d0*x3517 + x3520) + x3074*x3498*x3523 - x3082*
     & x3497*x3522)
      d2E(19) = -x3539*(-x1065*x1455 + x1511*x3265 - x3259*x3534 + x3259
     &*x3538 + x3537*( 0.0210328062336497d0*x1433 + x3525 + x3536))
      d2E(20) = -x3539*(-x1065*x1497 + x1511*x3375 - x3374*x3534 + x3374
     &*x3538 + x3537*( x3463 - x3465 + x3541))
      d2E(21) = -x3545*(-0.75d0*x1065*x3497 + x1161*x3544*(6.0d0*x1510 -
     & 16.8d0*x211 + x212*(8.4d0*x3496 + x3543) + x3542) + x1511*x3513)
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x3555 + x407)
      d1E(1) = x1060*(x1059 - x3569 - x3570*x745) + x1061 - x3555
      d1E(2) = -rho_c*(x1323 - x3587)
      d1E(3) = x1430*(x1429 - 0.000244140625d0*x3610)
      d1E(4) = x3443*x3625
      d1E(5) = x3443*x3634
      d1E(6) = -x1520*x3641
      d2E(1) = x2329*(-x1937*(-x1925*x3566*x3567 + x1936 - x3549*x3651*(
     &-x1918*(x1917 - 0.04375d0*x3564 + 0.065625d0*x3565 - x3580*x564) -
     & x1919*x3563 + x1922 + x3491*x3552 - 0.7d0*x3645 + 4.2d0*x3646 - 4
     &.2d0*x3647 + 2.1d0*x3649 - 0.35d0*x3650) - x3553*(-x1670*x3563 + x
     &1915 - 2.0d0 *x3645 + 12.0d0*x3646 - 12.0d0*x3647 + 6.0d0*x3649 - 
     &1.0d0*x3650 )) + x1938*x3570 - x1939*x3569 - x1940*x3570 + x2328) 
     &- x2330* x3570 + x2337 - 9.5367431640625d-7*x3569
      d2E(2) = x2706*(-x1937*(-x2474*x3549*(-x2377*x3580 - x2472*x3675 +
     & x2473 - 0.04375d0*x3670 + 0.13125d0*x3671 - 0.04375d0*x3672 - 0.0
     &4375d0* x3673 + 0.065625d0*x3674) - x2474*x3566*x3675 + x2487 - x3
     &553*( -x2377*x3547 + x2468 - x3670 + 3.0d0*x3671 - x3672 - x3673 +
     & 1.5d0*x3674) - x3568*(x2469 - x3576 - x3577)) + x2488*x3586 + x27
     &05) + x2708 + x3587
      d2E(3) = -rho_c*(-x248*(x2793*x3578*x3582 + x2824 - x3573*(x2714*x
     &3575 + x2719* x3547 + x2789 + x3677) + x3585*(x2719*x3580 + x2790*
     &x3575 + x2792 + 0.04375d0*x3677)) + x2890)
      d2E(4) = 4.76837158203125d-7*rho_c*(x247*(-x1409*x3616*(-x3051*(-x
     &1916*x3691 - x3041*x3702 + 0.04375d0*x3684 + 0.0875d0*x3689 - 0.13
     &125d0*x3690 - x3692*x3707 - 0.04375d0*x3694 + x3699*x3708 + 0.0157
     &746046752373d0*x3701 + x3709) + x3054*x3713 + x3072 + x3712 + x371
     &5 - x3717 + x3719) + x1410*x3724*x3727 + x3092 + x3619*( -x2958*x3
     &702 + x3684 + 2.0d0*x3689 - 3.0d0*x3690 - x3691*x535 - x3692*x3693
     & - x3694 + x3695*x3699 + 0.360562392576852d0*x3701 + x3706) - x372
     &0*x3722) + x3109 - x3610*x745) + 4.76837158203125d-7 *x2891 - 4.76
     &837158203125d-7*x3610
      d2E(5) = -x1060*(-x247*(-x2474*x3582*x3602 + x3168*x3735 + x3170 +
     & x3573*(-x2958* x3734 + x3163 - x3693*x3731 + x3695*x3733 - x3728 
     &+ 2.0d0*x3730 - x3732) - x3578*x3609 - x3585*(-x3041*x3734 + x3166
     & - x3707*x3731 + x3708*x3733 - 0.04375d0*x3728 + 0.0875d0*x3730 - 
     &0.04375d0* x3732)) + x3176)
      d2E(6) = -x3256*(-x3253*(-x3248*x3583*(x3242 - 0.13125d0*x3550 + 6
     &.0d0*x3607 - x74*(x3247 + x3546*x3757 + 0.04375d0*x3736 + x3741*x3
     &755 + 0.04375d0*x3748 - 0.0210328062336497d0*x3753 - 0.02103280623
     &36497d0*x3754)) + x3249*x3602*x3608 + x3252 - x3573* (x3241 + x373
     &6 + x3741*x3742 + x3748 + x3749*x3750 - 0.480749856769136d0*x3753 
     &- 0.480749856769136d0*x3754)) + x3255)
      d2E(7) = x1458*(rho_c*(-x3573*(x3281*x3547 - x3283*x3768 + x3288 +
     & x3591*x3772 + 0.360562392576852d0*x3765 + 0.25d0*x3767 - 0.75d0*x
     &3770) - x3585* (x2472*x3760 + x3289*x3768 + x3297 - x3546*x3773 - 
     &0.0157746046752373d0*x3765 - 0.0109375d0*x3767 + 0.0328125d0* x377
     &0 - x3772*x3774) + x3625*x3758 - x3760*x3761 - x3763*x3764) + 64.0
     &d0*x3618 + 64.0d0*x3624)
      d2E(8) = x3443*(x3613*x3777 + x3617*(x3316 + x3729*x3781 + x3774*x
     &3783 - 0.021875d0*x3780 + 0.0109375d0*x3782) + x3619*(x1342*x3786 
     &+ x3729*x3785 - 0.0078125d0*x3780 + 0.00390625d0*x3782 + x3788) + 
     &x3778*x3779)
      d2E(9) = -x3804*(x3573*(-x3283*x3751 + x3749*x3796 + 0.5d0*x3793 -
     & 0.240374928384568d0*x3795 + x3800 - 0.060093732096142d0*x3801 + x
     &3802) - x3585*(x3343 + x3546*x3803 - 0.0262910077920621d0*x3759 + 
     &0.0546875d0*x3793 - 0.0105164031168249d0*x3795 + 0.04375d0* x3800 
     &- 0.00262910077920621d0*x3801) + x3760*x3790 + x3791*x3792)
      d2E(10) = x3443*(-x3617*(x3339*x3806 + x3339*x3813 + 0.00525820155
     &841243d0*x3807 + x3810 - 0.00157992477879376d0*x3812) - x3619*(x36
     &20*x3806 + x3620 *x3813 + 0.00187792912800444d0*x3807 - 0.00056425
     &8849569202d0* x3812 + x3814) + 128.0d0*x3623*x3805)
      d2E(11) = x1458*(rho_c*(x3573*(x2940*x3827 + x3418 + 0.5d0*x3820 -
     & x3822 + 0.721124785153704d0*x3824 + 0.75d0*x3825 - 0.480749856769
     &136d0* x3826) - x3585*(x2472*x3816 + x3424 + 0.021875d0*x3820 + 0.
     &0315492093504746d0*x3824 + 0.0328125d0*x3825 - 0.0210328062336497d
     &0*x3826 + x3827*x3829 - x3828) + x3634*x3758 - x3761*x3816 - x3764
     &*x3819) + 64.0d0*x3629 + 64.0d0*x3633)
      d2E(12) = x3833*(0.015625d0*x3617*(x1479*x3832 + x3438 + x3828 + 0
     &.021875d0*x3830 - 0.0210328062336497d0*x3831) + 0.000244140625d0*x
     &3619*(x3439* x3827 + x3441 + x3822 + 0.5d0*x3830 - 0.4807498567691
     &36d0*x3831) + 0.015625d0*x3628*x3777 + x3632*x3778)
      d2E(13) = -x3804*(-x3573*(-0.480749856769136d0*x3627 - 0.5d0*x3835
     & + x3836*x3838 - 0.480749856769136d0*x3839 + x3840 + 0.06009373209
     &6142d0*x3841 + x3842) + x3585*(x3468 - 0.0315492093504746d0*x3815 
     &- 0.0546875d0* x3835 + x3836*x3845 - 0.0210328062336497d0*x3839 + 
     &0.04375d0* x3840 + 0.00262910077920621d0*x3841 - x3843) + x3790*x3
     &816 + x3792*x3834)
      d2E(14) = x3443*(-x3617*(x3454*x3849 + x3458*x3849 - x3473*x3809 +
     & 0.00262910077920621d0*x3847 - 0.00525820155841243d0*x3848 + x3850
     & *x3851) - x3619*(x3454*x3852 + x3458*x3852 - x3473*x98 + 0.000938
     &964564002219d0*x3847 - 0.00187792912800444d0*x3848 + x3850*x3853) 
     &+ 64.0d0*x3632*x3805 + x3779*x3846)
      d2E(15) = x3443*(x3617*(-x3339*x3854 + x3480 - x3485*x3849 + 0.010
     &5164031168249d0* x3855 + 0.00631969911517506d0*x3857) + x3619*(-x3
     &485*x3852 - x3620*x3854 + 0.00375585825600888d0*x3855 + 0.00225703
     &539827681d0 *x3857 + x3858) + 128.0d0*x3632*x3846)
      d2E(16) = x3507*(rho_c*(x185*(-x1514*x3616*(x3502*x3867 + x3506 - 
     &x3712 - x3715 + x3717 - x3719 + x724*(x3505 + 0.0875d0*x3861 - 0.1
     &3125d0*x3862 - 0.04375d0*x3863)) + x3499*x3727*x3866 - 2.0d0*x3619
     &*(x3495 + x3861 - 1.5d0*x3862 - 0.5d0*x3863) - 2.0d0*x3722*x3865) 
     &- x3489* x3641) - x1911*x3641)
      d2E(17) = -x3515*(x3514*x3735 + 0.015625d0*x3573*(-x3508*x3547 + x
     &3511 - x3868) + 0.015625d0*x3578*x3640 - x3582*x3871 - 0.5d0*x3584
     &*(-x3546*x3869 - 0.0875d0*x3868 + x3870))
      d2E(18) = x3533*(x3522*x3724*x3865 + x3523*x3720*x3866 - x3531*x36
     &16*(-x3491*x3867 + x3524 + x3530*(x3529 + 0.0875d0*x3872 + 0.0875d
     &0*x3873 - x3875* x3877) + 6.3d0*x3615 + 96.0d0*x3713) - x3619*(x35
     &20 + x3872 + x3873 - x3875*x3876))
      d2E(19) = -x3884*(-x3573*(-0.00375585825600888d0*x3612 + 0.0078125
     &d0*x3793 + x3880 ) + x3584*(-x3289*x3597 + x3536 + 0.021875d0*x379
     &3) - x3760*x3879 + x3760*x3881 + x3791*x3882)
      d2E(20) = -x3884*(-x3573*(0.00751171651201775d0*x3627 + 0.0078125d
     &0*x3835 + x3885 ) + x3584*(x3541 + 0.021875d0*x3835 + x3843) - x38
     &16*x3879 + x3816*x3881 + x3834*x3882)
      d2E(21) = -x3545*(x3544*x3583*(x212*(x3543 - 4.2d0*x3864) + x3542 
     &- 8.4d0*x3550 + 6.0d0*x3638) - 0.75d0*x3573*x3865 + x3639*x3871)
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = -rho_c*(x3895 + x407)
      d1E(1) = x1060*(x1059 - x3903 - x3904*x745) + x1061 - x3895
      d1E(2) = -rho_c*(x1323 - x3920)
      d1E(3) = x1430*(x1429 - 0.000244140625d0*x3937)
      d1E(4) = x3443*x3950
      d1E(5) = x3443*x3956
      d1E(6) = -x1520*x3963
      d2E(1) = x2329*(-x1937*(-x1925*x3900*x3901 + x1936 - x3651*x3889*(
     &-x1918*(x1917 - 0.04375d0*x3898 + 0.065625d0*x3899 - x3913*x564) -
     & x1919*x3897 + x1922 + x3491*x3892 - 0.7d0*x3965 + 4.2d0*x3966 - 4
     &.2d0*x3967 + 2.1d0*x3969 - 0.35d0*x3970) - x3893*(-x1670*x3897 + x
     &1915 - 2.0d0 *x3965 + 12.0d0*x3966 - 12.0d0*x3967 + 6.0d0*x3969 - 
     &1.0d0*x3970 )) + x1938*x3904 - x1939*x3903 - x1940*x3904 + x2328) 
     &- x2330* x3904 + x2337 - 9.5367431640625d-7*x3903
      d2E(2) = x2706*(-x1937*(-x2474*x3889*(-x2377*x3913 - x2472*x3977 +
     & x2473 - 0.04375d0*x3972 + 0.13125d0*x3973 - 0.04375d0*x3974 - 0.0
     &4375d0* x3975 + 0.065625d0*x3976) - x2474*x3900*x3977 + x2487 - x3
     &893*( -x2377*x3887 + x2468 - x3972 + 3.0d0*x3973 - x3974 - x3975 +
     & 1.5d0*x3976) - x3902*(x2469 - x3909 - x3910)) + x2488*x3919 + x27
     &05) + x2708 + x3920
      d2E(3) = -rho_c*(-x248*(x2793*x3911*x3915 + x2824 - x3907*(x2714*x
     &3908 + x2719* x3887 + x2789 + x3979) + x3918*(x2719*x3913 + x2790*
     &x3908 + x2792 + 0.04375d0*x3979)) + x2890)
      d2E(4) = 4.76837158203125d-7*rho_c*(x247*(-x1409*x3943*(-x3051*(-x
     &1916*x3987 - x3041*x3993 - x3707*x3988 + x3708*x3991 + x3709 + 0.0
     &4375d0*x3982 + 0.0875d0*x3985 - 0.13125d0*x3986 - 0.04375d0*x3989 
     &+ 0.0157746046752373d0*x3992) + x3054*x3996 + x3072 + x3995 + x399
     &8 - x4000 + x4002) + x1410*x4006*x4010 + x3092 + x3946*(-x2958* x3
     &993 - x3693*x3988 + x3695*x3991 + x3706 + x3982 + 2.0d0*x3985 - 3.
     &0d0*x3986 - x3987*x535 - x3989 + 0.360562392576852d0*x3992) - x400
     &3*x4005) + x3109 - x3937*x745) + 4.76837158203125d-7*x2891 - 4.768
     &37158203125d-7*x3937
      d2E(5) = -x1060*(-x247*(-x2474*x3915*x3930 + x3168*x4018 + x3170 +
     & x3907*(x2958* x4016 + x3135*x4008 + x3163 - x3693*x4014 - x4011 -
     & 2.0d0*x4013 - x4015) - x3911*x3936 - x3918*(x3041*x4016 + x3135*x
     &4017 + x3166 - x3707*x4014 - 0.04375d0*x4011 - 0.0875d0*x4013 - 0.
     &04375d0*x4015 )) + x3176)
      d2E(6) = -x3256*(-x3253*(-x3248*x3916*(x3242 - 0.13125d0*x3890 + 6
     &.0d0*x3934 - x74*(x3247 + x3755*x4021 + x3757*x3886 + 0.04375d0*x4
     &019 + 0.04375d0*x4024 - 0.0210328062336497d0*x4027 - 0.02103280623
     &36497d0*x4028)) + x3249*x3930*x3935 + x3252 - x3907* (x3241 + x374
     &2*x4021 + x3750*x4025 + x4019 + x4024 - 0.480749856769136d0*x4027 
     &- 0.480749856769136d0*x4028)) + x3255)
      d2E(7) = x1458*(rho_c*(x3758*x3950 - x3907*(-x1370*x3283*x3897 + x
     &3281*x3887 + x3288 + x3772*x3925 + 0.360562392576852d0*x4034 + 0.2
     &5d0*x4036 - 0.75d0*x4037) - x3918*(-0.0109375d0*x1353*x3772 + x247
     &2*x4030 + x3297 - x3773*x3886 + x3781*x3897 - 0.0157746046752373d0
     &*x4034 - 0.0109375d0*x4036 + 0.0328125d0*x4037) - x4030*x4031 - x4
     &032* x4033) + 64.0d0*x3945 + 64.0d0*x3949)
      d2E(8) = x3443*(x3941*x4039 + x3944*(x1353*x3832 + x3316 - x3781*x
     &4012 - 0.021875d0*x4042 + 0.0109375d0*x4043) + x3946*(x1353*x3786 
     &- x3785*x4012 + x3788 - 0.0078125d0*x4042 + 0.00390625d0*x4043) + 
     &x4040*x4041)
      d2E(9) = -x3804*(x3907*(-x3283*x4026 + x3796*x4025 + x3802 + 0.5d0
     &*x4047 - 0.240374928384568d0*x4048 + x4050 - 0.060093732096142d0*x
     &4051) - x3918*(x3343 + x3803*x3886 - 0.0262910077920621d0*x4029 + 
     &0.0546875d0*x4047 - 0.0105164031168249d0*x4048 + 0.04375d0*x4050 -
     & 0.00262910077920621d0*x4051) + x4030*x4044 + x4045*x4046)
      d2E(10) = x3443*(-x3944*(x3339*x4053 + x3362*x3849 + x3810 + 0.005
     &25820155841243d0 *x4054 + 0.00157992477879376d0*x4056) - x3946*(x3
     &362*x3852 + x3620*x4053 + x3814 + 0.00187792912800444d0*x4054 + 0.
     &000564258849569202d0*x4056) + 128.0d0*x3948*x4052)
      d2E(11) = x1458*(rho_c*(x3758*x3956 + x3907*(-x2940*x4064 + x3407*
     &x4008 + x3418 - 0.25d0*x4060 + 0.721124785153704d0*x4061 - 0.75d0*
     &x4062 - 0.480749856769136d0*x4063) - x3918*(x2472*x4058 + x3407*x4
     &017 + x3424 - x3829*x4064 - 0.0109375d0*x4060 + 0.0315492093504746
     &d0* x4061 - 0.0328125d0*x4062 - 0.0210328062336497d0*x4063) - x403
     &1* x4058 - x4033*x4059) + 64.0d0*x3953 + 64.0d0*x3955)
      d2E(12) = x3833*(0.015625d0*x3944*(-x1490*x3832 + x3438 + 0.021875
     &d0*x4065 + 0.0109375d0*x4066 + 0.0210328062336497d0*x4067) + 0.000
     &244140625d0*x3946*(-x3439*x4064 + x3441 + 0.5d0*x4065 + 0.25d0*x40
     &66 + 0.480749856769136d0*x4067) + 0.015625d0*x3952* x4039 + x3954*
     &x4040)
      d2E(13) = -x3804*(-x3907*(x3838*x4070 + x3842 - 0.480749856769136d
     &0*x3951 + 0.5d0* x4069 - 0.480749856769136d0*x4071 - x4072 - 0.060
     &093732096142d0* x4073) + x3918*(-x3372*x3462 + x3468 + x3845*x4070
     & - 0.0315492093504746d0*x4057 - 0.0210328062336497d0*x4071 - 0.043
     &75d0*x4072 - 0.00262910077920621d0*x4073 + x4074 - x4075) + x4044*
     &x4058 + x4046*x4068)
      d2E(14) = x3443*(x3944*(x3455*x3849 + x3460*x3849 + x3474 - x3851*
     &x4079 + 0.00525820155841243d0*x4077 + 0.00262910077920621d0*x4078)
     & + x3946*(x3455*x3852 + x3460*x3852 + 0.00225703539827681d0*x3476 
     &- x3853*x4079 + 0.00187792912800444d0*x4077 + 0.000938964564002219
     &d0*x4078) + 64.0d0*x3954*x4052 + x4041*x4076)
      d2E(15) = x3443*(x3944*(x3480 - x3483*x3849 - x3486*x3849 - 0.0105
     &164031168249d0* x4080 + 0.00631969911517506d0*x4081) + x3946*(-x34
     &83*x3852 - x3486*x3852 + x3858 - 0.00375585825600888d0*x4080 + 0.0
     &0225703539827681d0*x4081) + 128.0d0*x3954*x4076)
      d2E(16) = x3507*(rho_c*(x185*(-x1514*x3943*(x3502*x4089 + x3506 - 
     &x3995 - x3998 + x4000 - x4002 + x724*(x3505 + 0.0875d0*x4083 + 0.1
     &3125d0*x4084 + 0.04375d0*x4085)) + x3499*x4010*x4088 - 2.0d0*x3946
     &*(x3495 + x4083 + 1.5d0*x4084 + 0.5d0*x4085) - 2.0d0*x4005*x4087) 
     &- x3489* x3963) - x1911*x3963)
      d2E(17) = -x3515*(x3514*x4018 + 0.015625d0*x3907*(-x3508*x3887 + x
     &3511 - x4090) + 0.015625d0*x3911*x3962 - x3915*x4091 - 0.5d0*x3917
     &*(-x3869*x3886 + x3870 - 0.0875d0*x4090))
      d2E(18) = x3533*(x3522*x4006*x4087 + x3523*x4003*x4088 - x3531*x39
     &43*(-x3491*x4089 + x3524 + x3530*(x3529 + x3877*x4094 - 0.0875d0*x
     &4092 + 0.0875d0* x4093) + 6.3d0*x3942 + 96.0d0*x3996) - x3946*(x35
     &20 + x3876*x4094 - 1.0d0*x4092 + x4093))
      d2E(19) = -x3884*(-x3907*(-x3785*x3923 + x3880 + 0.0078125d0*x4047
     &) + x3917*(x3536 - x3781*x3923 + 0.021875d0*x4047) - x4030*x4095 +
     & x4030*x4096 + x4045*x4097)
      d2E(20) = -x3884*(-x3907*(x3885 + 0.00751171651201775d0*x3951 - 0.
     &0078125d0*x4069 ) + x3917*(x3541 - x4074 + x4075) - x4058*x4095 + 
     &x4058*x4096 + x4068*x4097)
      d2E(21) = -x3545*(x3544*x3916*(x212*(x3543 - 4.2d0*x4086) + x3542 
     &- 8.4d0*x3890 + 6.0d0*x3960) - 0.75d0*x3907*x4087 + x3961*x4091)
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = -rho_c*(x407 + x4106)
      d1E(1) = x1060*(x1059 - x4116 - x4117*x745) + x1061 - x4106
      d1E(2) = -rho_c*(x1323 - x4135)
      d1E(3) = x1430*(x1429 - 0.000244140625d0*x4152)
      d1E(4) = 0
      d1E(5) = 0
      d1E(6) = -x1520*x4158
      d2E(1) = x2329*(-x1937*(-x1925*x4112*x4114 + x1936 - x3651*x4100*(
     &-x1918*(-x1669* x4127 + x1917 - 0.04375d0*x4109 + 0.2625d0*x4110) 
     &- x1919*x4108 + x1922 + x3491*x4103 - 16.8d0*x4159 - 0.7d0*x4160 +
     & 4.2d0*x4161 - 1.4d0*x4162 + 8.4d0*x4163) - x4104*(-x1670*x4108 + 
     &x1915 - 48.0d0 *x4159 - 2.0d0*x4160 + 12.0d0*x4161 - 4.0d0*x4162 +
     & 24.0d0*x4163 )) + x1938*x4117 - x1939*x4116 - x1940*x4117 + x2328
     &) - x2330* x4117 + x2337 - 9.5367431640625d-7*x4116
      d2E(2) = x2706*(-x1937*(-x2474*x4100*(-x1668*x2791*x3669 - x2376*x
     &4169 - x2472* x4168 + x2473 + 0.2625d0*x4122*x563 - 0.04375d0*x416
     &4 + 0.2625d0* x4165 - 0.04375d0*x4166) - x2474*x4112*x4168 + x2487
     & - x4104*( x198*x203*x3502*x4122 - x2376*x4167 + x2468 - x3669*x41
     &11 - x4164 + 6.0d0*x4165 - x4166) - x4115*(x2469 - x4121 - x4124))
     & + x2488* x4134 + x2705) + x2708 + x4135
      d2E(3) = -rho_c*(-x248*(x2793*x4125*x4130 + x2824 - x4119*(x1911*x
     &4170 + x2718* x4167 + x2789 + 4.0d0*x3139*x3669) + x4133*(x1401*x4
     &170 + x203* x3669*x4171 + x2718*x4169 + x2792)) + x2890)
      d2E(4) = 4.76837158203125d-7*rho_c*(x247*(-x1409*x4186*(x3051*(x13
     &80*x3526*x4180 + 0.00546875d0*x1913*x2939*x3527 + x3041*x4184 + x3
     &050 - 0.0875d0 *x4175 + 0.525d0*x4176 - 0.04375d0*x4177 - 0.0875d0
     &*x4182 + 0.01640625d0*x4183 - 0.0630984187009491d0*x4185) + x3054*
     &x4196 + x4188 + x4190 - x4192 + x4195 + x4198 - x4200) - x1410*x42
     &04* x4205 + x3092 - x4172*(x2939*x3028 + x2958*x4184 + x3029 - 2.0
     &d0* x4175 + 12.0d0*x4176 - x4177 + x4178*x4180 - 2.0d0*x4182 + 0.3
     &75d0*x4183 - 1.44224957030741d0*x4185) - x4201*x4202) + x3109 - x4
     &152*x745) + 4.76837158203125d-7*x2891 - 4.76837158203125d-7* x4152
     &
      d2E(5) = -x1060*(-x247*(-x2474*x4130*x4144 + x3168*x4214 + x3170 +
     & x4119*(-x1911* x4208 - x2379*x3687*x4178 + x2960*x4206*x4211 + x3
     &162 - 4.0d0* x3503*x4207 + x4209*x4210) - x4125*x4151 - x4133*(-x1
     &380*x3687* x4171 - x1401*x4208 - x200*x3526*x4207 + x2960*x3409*x4
     &213 + x3165 + x4210*x4212)) + x3176)
      d2E(6) = -x3256*(-x3253*(-x3248*x4131*(x3242 - 0.525d0*x4145 + 6.0
     &d0*x4149 - x74* (x3247 + x3526*x4215 + 0.021875d0*x4216 + 0.006319
     &69911517506d0* x4217 - 0.0841312249345988d0*x4218 + 0.002734375d0*
     &x4219*x622* x696 - 0.00262910077920621d0*x4220)) + x3249*x4144*x41
     &50 + x3252 - x4119*(x1725*x4215 + x3241 + x3746*x4219 + 0.5d0*x421
     &6 + 0.144450265489716d0*x4217 - 1.92299942707654d0*x4218 - 0.06009
     &3732096142d0*x4220)) + x3255)
      d2E(7) = x4225*(x4172*(-x1378*(x3006*x4221 + x3287 + x3697) + x145
     &2*x3065 + 1.44224957030741d0*x1465*x3687 - x3283*x4222 + x4209*(x3
     &273 + x3279*x3402 + x4221*x558)) + x4224*(0.0630984187009491d0*x14
     &65* x3562 + x3280*x4127 - x3289*x4108 - x3294 + x3296))
      d2E(8) = -x2959*x3532*(-x4119*(x1452*x4206 + x185*x4226 - 32.0d0*x
     &3317 - 0.120187464192284d0*x3787) + 32.0d0*x4132*(x1401*x4226 + x1
     &438* x4213 - x3314 - x3315))
      d2E(9) = x3345*x4229*(-x4172*(x3328*x98 + 0.0150234330240355d0*x42
     &27 - 0.00225703539827681d0*x4228) - x4223*(x3342*x98 + 0.042065612
     &4672994d0*x4227 - 0.00631969911517506d0*x4228))
      d2E(10) = x3325*x4186*x4230*(-5.42101086242752d-20*x1434 + 8.67361
     &737988404d-19* x240 - 1.35525271560688d-20*x242 + 4.33680868994202
     &d-19*x96)
      d2E(11) = -x1370*x1458*x190*(-x4119*(x3405*x4209 - x3416 + x3417 +
     & 2.88449914061482d0*x4139 - 0.480749856769136d0*x4231) + x4133*( x
     &3405*x4212 - x3422 + x3423 + 0.126196837401898d0*x4139 - 0.0210328
     &062336497d0*x4231))
      d2E(12) = x4225*(x4172*(-x1495*x4233 + x3431*x4209 - x3440 + 0.961
     &499713538272d0* x4232) + x4224*(-x3291*x4233 + x3431*x4212 - x3437
     & + 0.0420656124672994d0*x4232))
      d2E(13) = -x202*x3883*x4229*x684*(-x4119*(-x3449 - 0.0300468660480
     &71d0*x4234 + 0.00451407079655361d0*x4235) + x4132*(-x3467 - 0.0841
     &312249345988d0*x4234 + 0.0126393982303501d0*x4235))
      d2E(14) = -rho_s*x3321*x4131*x4230*(-1.0842021724855d-19*x1434 + 1
     &.73472347597681d-18*x240 - 2.71050543121376d-20*x242 + 8.673617379
     &88404d-19*x96)
      d2E(15) = x3224*x3299*x3346*x3478*x4186*(-2.16840434497101d-19*x14
     &34 + 3.46944695195361d-18*x240 - 5.42101086242752d-20*x242 + 1.734
     &72347597681d-18*x96)
      d2E(16) = x3507*(rho_c*(x185*(-x1514*x4186*(x3502*x4241 - x4188 - 
     &x4190 + x4192 - x4195 - x4198 + x4200 - x724*(-x3069*(0.0875d0*x29
     &38 + x3504) + 0.0875d0*x3704 + 0.525d0*x4236 + x4237*x4242 - 0.087
     &5d0*x4238)) - x3499*x4205*x4240 + 2.0d0*x4172*(x1691*x4237 - x3025
     &*x3069 + x3704 + 6.0d0*x4236 - x4238) - 2.0d0*x4202*x4239) - x3489
     &*x4158) - x1911*x4158)
      d2E(17) = -x3515*(x3514*x4214 + 0.015625d0*x4119*(-x1691*x4122 + x
     &3511 - 2.0d0* x4120) + 0.015625d0*x4125*x4157 - x4130*x4243 - 0.5d
     &0*x4132*( x3870 - 0.175d0*x4120 - x4122*x4242))
      d2E(18) = x3533*(x3522*x4204*x4239 + x3523*x4201*x4240 - x3531*x41
     &86*(-x3491*x4241 + x3524 + x3530*(-x3420*x4246 + x3528*x4245 + x35
     &29 + 0.35d0* x4244) + 96.0d0*x4196 + 25.2d0*x4203) - x4172*(x3027*
     &x4245 + x3520 - x4211*x4246 + 4.0d0*x4244))
      d2E(19) = 0
      d2E(20) = 0
      d2E(21) = -x3545*(x3544*x4131*(x212*(x3543 - 16.8d0*x4136) + x3542
     & - 33.6d0*x4145 + 6.0d0*x4155) - 0.75d0*x4119*x4239 + x4156*x4243)
     &
      end if
      end subroutine
