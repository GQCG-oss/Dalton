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
C   Generated: November 27, 2018
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
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
      pure subroutine ESRX_TPSS_ERF_case_1(rho_a, gamma_aa, tau_a, mu, E
     &a)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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


C*****************************************************************************
      pure subroutine ESRC_PW92_ERF(rho_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = (0.132568899905202d0*x1 + 1.0d0)*log(1.0d0 + 16.0819842214574
     &07d0/( 2.2255694211506878d0*x1 + 5.9825504357710795d0*x2 + 0.80042
     &863499936356d0*x3 + 0.18970043257475591d0*x4))
      x6 = mu**2
      x7 = x0**1.66666666666667d0
      x8 = x0**1.0d0
      x9 = x0**1.33333333333333d0
      x10 = (0.0128480790170175d0*x1 + 0.0315297404574782d0*x4 - 0.00304
     &892328706429d0*x8 + 0.000275310793616358d0*x9 + 1.0d0)*exp (-0.466
     &75853320810856d0*x1)
      x11 = x1*x6
      x12 = 4.60448529760779d0*mu*x2 + 1.0d0
      x13 = mu**3
      x14 = x10 - 1.0d0
      x15 = 1d0/(0.337569861992026d0*x1 + 0.0244354823084343d0*x4 + 1.0d
     &0)
      E = -rho_c*(0.0621813817393098d0*x5 - (0.000196557189785076d0*mu**
     &8*x0** 2.66666666666667d0*x5 + mu**6*(0.00331582639355101d0*x0**2.
     &0d0*x5 + 0.00248266686171384d0*x14*x7) + 0.00264117541143784d0*mu*
     &*5*x10 *x7 - mu**4*(0.00020688890514282d0*x0**1.0d0*(-x15*( -0.305
     &045320493918d0*x1 + 17.2274040011292d0) - ( -5.62907590342551d0*x1
     & + 6.084d0*x4)*exp(-0.33933171852197186d0* x1) + 17.2274040011292d
     &0) - 0.0418814102677522d0*x14*x8 - 0.0209761558789675d0*x5*x9) + x
     &13*(0.0445553741830535d0*x10*x8 + 0.000132058770571892d0*x7*(17.22
     &74040011292d0*x15*( -0.0177069812999059d0*x1 + 1.0d0)/x4 + 37.6991
     &118430775d0*( -3.07073492995203d0*x1 + 0.384834731559127d0*x4)*exp
     &( -0.19230865217881402d0*x1)/x8)) + 0.0621813817393098d0*log(( 4.6
     &213167280845804d0*x11 + x12 + 1.9140710242289749d0*x13*x3)/( 2.139
     &2848713714905d0*x11 + x12)))/(0.237113969739022d0*x4*x6 + 1.0d0)**
     &4)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_PW92_ERF(rho_c, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, mu
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67, x68, x69, x70, x71, x72, x73
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = x0**0.666666666666667d0
      x6 = 2.22556942115069d0*x1 + 5.98255043577108d0*x3 + 0.80042863499
     &9364d0*x4 + 0.189700432574756d0*x5
      x7 = 1.0d0 + 16.0819842214574d0/x6
      x8 = log(x7)
      x9 = x2*x8
      x10 = 0.0621813817393098d0*x9
      x11 = mu**2
      x12 = x11*x5
      x13 = 0.237113969739022d0*x12 + 1.0d0
      x14 = x13**(-4)
      x15 = mu**5
      x16 = x0**1.66666666666667d0
      x17 = exp(-0.46675853320810856d0*x1)
      x18 = x0**1.0d0
      x19 = 0.00304892328706429d0*x18
      x20 = x0**1.33333333333333d0
      x21 = 0.0128480790170175d0*x1 - x19 + 0.000275310793616358d0*x20 +
     & 0.0315297404574782d0*x5 + 1.0d0
      x22 = x17*x21
      x23 = x16*x22
      x24 = 0.00264117541143784d0*x15*x23
      x25 = x0**2.66666666666667d0
      x26 = x25*x9
      x27 = 0.000196557189785076d0*mu**8*x26
      x28 = x1*x11
      x29 = 4.60448529760779d0*mu*x3 + 1.0d0
      x30 = mu**3
      x31 = 4.62131672808458d0*x28 + x29 + 1.91407102422897d0*x30*x4
      x32 = x31/(2.13928487137149d0*x28 + x29)
      x33 = -0.0621813817393098d0
      x34 = x33*log(x32)
      x35 = x22 - 1.0d0
      x36 = x16*x35
      x37 = x0**2.0d0
      x38 = x37*x9
      x39 = mu**6*(0.00248266686171384d0*x36 + 0.00331582639355101d0*x38
     &)
      x40 = 0.0445553741830535d0*x18
      x41 = x22*x40
      x42 = exp(-0.19230865217881402d0*x1)
      x43 = x42*(-3.07073492995203d0*x1 + 0.384834731559127d0*x5)
      x44 = 37.6991118430775d0/x18
      x45 = x43*x44
      x46 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x5 + 1.0d0
      x47 = 1d0/x46
      x48 = -0.0177069812999059d0*x1 + 1.0d0
      x49 = x47*x48
      x50 = 1d0/x5
      x51 = 17.2274040011292d0*x50
      x52 = x16*(x45 + x49*x51)
      x53 = x30*(x41 + 0.000132058770571892d0*x52)
      x54 = mu**4
      x55 = 0.0418814102677522d0*x18*x35
      x56 = x20*x9
      x57 = exp(-0.33933171852197186d0*x1)
      x58 = x0**1.0d0*(-x47*(-0.305045320493918d0*x1 + 17.2274040011292d
     &0) - x57*( -5.62907590342551d0*x1 + 6.084d0*x5) + 17.2274040011292
     &d0)
      x59 = x54*(-x55 - 0.0209761558789675d0*x56 + 0.00020688890514282d0
     &*x58)
      x60 = x14*(x24 + x27 - x34 + x39 + x53 - x59)
      x61 = (0.741856473716896d0*x1 + 0.99709173929518d0*x3 + 0.40021431
     &7499682d0*x4 + 0.126466955049837d0*x5)/(x6**2*x7)
      x62 = mu**7
      x63 = x17*(-0.00428269300567249d0*x1 + x19 - 0.000367081058155144d
     &0*x20 - 0.0210198269716521d0*x5)
      x64 = x2*x61
      x65 = 0.767414216267966d0*x3
      x66 = mu*x1
      x67 = x17*(0.0777930888680181d0*x1*x21 - 0.00214134650283625d0*x1 
     &+ 0.00152446164353214d0*x18 - 0.000183540529077572d0*x20 - 0.01050
     &99134858261d0*x5)
      x68 = 1d0/x1
      x69 = 0.101681773497972d0*x47*x68
      x70 = 11.4849360007528d0*x50
      x71 = x49*x70
      x72 = x48*x51*(0.112523287330675d0*x1 + 0.0162903215389562d0*x5)/x
     &46**2
      x73 = x57*(-0.240695990468967d0*x1 + 0.26014827853397d0*x5)
      E = -rho_c*(x10 - x60)
      d1E(1) = -mu*x14*(-0.000410930386996015d0*x0**2.0d0*x22*x54 + 8.68
     &579013942185d-6 *x0**3.0d0*x62*x8 + x11*(-0.000132058770571892d0*x
     &16*( 2.41662179562687d0*x0**(-0.666666666666667d0)*x43 + x42*x44*(
     & 1.02357830998401d0*x1 - 0.256556487706084d0*x5) + x45 + x69 + x71
     & + x72) - 0.00693220036674015d0*x20*x22 - x40*x63 + x41 + 0.000220
     &097950953153d0*x52) + x15*(0.00014652515242323d0*x0** 2.3333333333
     &3333d0*x8 - 0.00496533372342768d0*x16*x67 + 0.0041377781028564d0*x
     &36 - 0.0533250677421794d0*x37*x64 + 0.00663165278710203d0*x38) - 0
     &.00264117541143784d0*x16*x54*x63 + 0.00440195901906306d0*x23*x54 -
     & 0.00316102962473761d0*x25*x62*x64 + 0.000524152506093537d0*x26*x6
     &2 + x30*(0.000926928636371586d0* x16*x8 + 0.00020688890514282d0*x1
     &6*(-23.3866625383245d0*x50*x57*( 0.0802319968229891d0*x1 - 0.17343
     &2185689313d0*x5) - 15.591108358883d0*x50*x73 - 2.64527879654103d0*
     &x68*x73 - x69 + x70 - x71 - x72) - 0.0837628205355044d0*x18*x67 - 
     &0.337338207872387d0*x20*x64 + x55 + 0.0279682078386234d0*x56 - 0.0
     &003448148419047d0*x58) + x33*(-0.957035512114487d0*x11*x4 + x32*(x
     &65 + 0.71309495712383d0*x66) - x65 - 1.54043890936153d0*x66 )/x31)
     & + 0.00274777245725523d0*x1*x8 - x10 - 0.632303919304059d0* x12*(-
     &x24 - x27 + x34 - x39 - x53 + x59)/x13**5 + x60 - 16.0819842214574
     &d0*x61*(0.0082433173717657d0*x1 + 0.0621813817393098d0)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_PW92_ERF(rho_c, mu, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x171, x172, x173
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = 0.800428634999364d0*x4
      x6 = x0**0.666666666666667d0
      x7 = 2.22556942115069d0*x1 + 5.98255043577108d0*x3 + x5 + 0.189700
     &432574756d0 *x6
      x8 = 1.0d0 + 16.0819842214574d0/x7
      x9 = log(x8)
      x10 = x2*x9
      x11 = 0.0621813817393098d0*x10
      x12 = mu**2
      x13 = x12*x6
      x14 = 0.237113969739022d0*x13 + 1.0d0
      x15 = x14**(-4)
      x16 = mu**5
      x17 = x0**1.66666666666667d0
      x18 = exp(-0.46675853320810856d0*x1)
      x19 = x0**1.0d0
      x20 = 0.00304892328706429d0*x19
      x21 = -x20
      x22 = x0**1.33333333333333d0
      x23 = 0.0128480790170175d0*x1 + x21 + 0.000275310793616358d0*x22 +
     & 0.0315297404574782d0*x6 + 1.0d0
      x24 = x18*x23
      x25 = x17*x24
      x26 = 0.00264117541143784d0*x16*x25
      x27 = x0**2.66666666666667d0
      x28 = x10*x27
      x29 = 0.000196557189785076d0*mu**8*x28
      x30 = x1*x12
      x31 = 4.60448529760779d0*mu*x3 + 1.0d0
      x32 = 2.13928487137149d0*x30 + x31
      x33 = 1d0/x32
      x34 = mu**3
      x35 = 4.62131672808458d0*x30 + x31 + 1.91407102422897d0*x34*x4
      x36 = x33*x35
      x37 = 0.693147180559945d0
      x38 = 0.202642367284676d0*x37 - 0.202642367284676d0
      x39 = x38*log(x36)
      x40 = x24 - 1.0d0
      x41 = x17*x40
      x42 = x0**2.0d0
      x43 = x10*x42
      x44 = mu**6*(0.00248266686171384d0*x41 + 0.00331582639355101d0*x43
     &)
      x45 = 0.0445553741830535d0*x19
      x46 = x24*x45
      x47 = exp(-0.19230865217881402d0*x1)
      x48 = x47*(-3.07073492995203d0*x1 + 0.384834731559127d0*x6)
      x49 = 1d0/x19
      x50 = 37.6991118430775d0*x49
      x51 = x48*x50
      x52 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x6 + 1.0d0
      x53 = 1d0/x52
      x54 = -0.0177069812999059d0*x1 + 1.0d0
      x55 = x53*x54
      x56 = 1d0/x6
      x57 = 17.2274040011292d0*x56
      x58 = x17*(x51 + x55*x57)
      x59 = x34*(x46 + 0.000132058770571892d0*x58)
      x60 = mu**4
      x61 = x19*x40
      x62 = 0.0418814102677522d0*x61
      x63 = x10*x22
      x64 = x0**1.0d0
      x65 = exp(-0.33933171852197186d0*x1)
      x66 = x65*(-5.62907590342551d0*x1 + 6.084d0*x6)
      x67 = x53*(-0.305045320493918d0*x1 + 17.2274040011292d0)
      x68 = x64*(-x66 - x67 + 17.2274040011292d0)
      x69 = x60*(-x62 - 0.0209761558789675d0*x63 + 0.00020688890514282d0
     &*x68)
      x70 = x26 + x29 - x39 + x44 + x59 - x69
      x71 = x15*x70
      x72 = x1*x9
      x73 = 0.741856473716896d0*x1 + 0.99709173929518d0*x3 + 0.400214317
     &499682d0*x4 + 0.126466955049837d0*x6
      x74 = 1d0/x8
      x75 = x74/x7**2
      x76 = x73*x75
      x77 = 0.0082433173717657d0*x1 + 0.0621813817393098d0
      x78 = 16.0819842214574d0*x77
      x79 = -x26 - x29 + x39 - x44 - x59 + x69
      x80 = x14**(-5)
      x81 = x13*x80
      x82 = x79*x81
      x83 = 0.00428269300567249d0*x1
      x84 = 0.0210198269716521d0*x6
      x85 = 0.000367081058155144d0*x22
      x86 = x18*(x20 - x83 - x84 - x85)
      x87 = x45*x86
      x88 = x22*x24
      x89 = 0.00693220036674015d0*x88
      x90 = 0.000220097950953153d0*x58
      x91 = 1.02357830998401d0*x1 - 0.256556487706084d0*x6
      x92 = x47*x50
      x93 = x0**(-0.666666666666667d0)
      x94 = x48*x93
      x95 = 1d0/x1
      x96 = x53*x95
      x97 = 0.101681773497972d0*x96
      x98 = 11.4849360007528d0*x56
      x99 = x55*x98
      x100 = x52**(-2)
      x101 = 0.112523287330675d0*x1 + 0.0162903215389562d0*x6
      x102 = x100*x101
      x103 = x54*x57
      x104 = x102*x103
      x105 = x104 + x97 + x99
      x106 = x105 + x51 + x91*x92 + 2.41662179562687d0*x94
      x107 = 0.000132058770571892d0*x17
      x108 = x106*x107
      x109 = x65*(-0.240695990468967d0*x1 + 0.26014827853397d0*x6)
      x110 = x109*x56
      x111 = 15.591108358883d0*x110
      x112 = 2.64527879654103d0*x109*x95
      x113 = x65*(0.0802319968229891d0*x1 - 0.173432185689313d0*x6)
      x114 = x113*x56
      x115 = 23.3866625383245d0*x114
      x116 = -x104 - x111 - x112 - x115 - x97 + x98 - x99
      x117 = 0.00020688890514282d0*x17
      x118 = x17*x9
      x119 = x1*x23
      x120 = -0.00214134650283625d0*x1 + 0.00152446164353214d0*x19 - 0.0
     &00183540529077572d0*x22 - 0.0105099134858261d0*x6
      x121 = x18*(0.0777930888680181d0*x119 + x120)
      x122 = x121*x19
      x123 = x2*x76
      x124 = x123*x22
      x125 = 0.000926928636371586d0*x118 - 0.0837628205355044d0*x122 - 0
     &.337338207872387d0*x124 + x62 + 0.0279682078386234d0*x63
      x126 = 0.767414216267966d0*x3
      x127 = mu*x1
      x128 = 1.54043890936153d0*x127
      x129 = x12*x4
      x130 = 0.957035512114487d0*x129
      x131 = x126 + 0.71309495712383d0*x127
      x132 = -x126 - x128 - x130 + x131*x36
      x133 = 1d0/x35
      x134 = x133*x38
      x135 = x132*x134
      x136 = x60*x86
      x137 = x136*x17
      x138 = x0**2.33333333333333d0
      x139 = x138*x9
      x140 = x121*x17
      x141 = 0.0533250677421794d0*x42
      x142 = mu**7
      x143 = x0**3.0d0*x142
      x144 = x143*x9
      x145 = x142*x28
      x146 = x25*x60
      x147 = x0**2.0d0
      x148 = x24*x60
      x149 = x147*x148
      x150 = x142*x27
      x151 = x123*x150
      x152 = 8.68579013942185d-6*x144 + 0.000524152506093537d0*x145 + 0.
     &00440195901906306d0*x146 - 0.000410930386996015d0*x149 - 0.0031610
     &2962473761d0*x151 + x16*(-x123*x141 + 0.00014652515242323d0*x139 -
     & 0.00496533372342768d0*x140 + 0.0041377781028564d0*x41 + 0.0066316
     &5278710203d0*x43)
      x153 = x12*(-x108 + x46 - x87 - x89 + x90) + x135 - 0.002641175411
     &43784d0*x137 + x152 + x34*(x116*x117 + x125 - 0.0003448148419047d0
     &*x68)
      x154 = mu*x15
      x155 = x75*(0.989141964955861d0*x1 + 1.16327369584438d0*x3 + 0.600
     &321476249523d0*x4 + 0.210778258416395d0*x6)
      x156 = x73**2/(x7**4*x8**2)
      x157 = x73*x74*(1.48371294743379d0*x1 + 1.99418347859036d0*x3 + x5
     & + 0.252933910099675d0*x6)/x7**3
      x158 = 0.00264117541143784d0*x17*x60
      x159 = x18*(0.00571025734089666d0*x1 - 0.00609784657412857d0*x19 +
     & 0.000856522469028669d0*x22 + 0.0350330449527535d0*x6)
      x160 = 0.00316102962473761d0*x150*x2
      x161 = x156*x2
      x162 = x126 + x128 + x130
      x163 = mu*x131
      x164 = x163*x33
      x165 = 0.89531658564596d0*x3
      x166 = 0.0891107483661069d0*x19
      x167 = x47*x91
      x168 = 3.82831200025092d0*x56
      x169 = x100*x103*(0.150031049774234d0*x1 + 0.0271505358982603d0*x6
     &) - x101*x103 *(0.225046574661351d0*x1 + 0.0325806430779124d0*x6)/
     &x52**3 - 22.9698720015055d0*x102*x54*x56 - 0.203363546995945d0*x10
     &2*x95 + x168*x55 + 1.38777878078145d-17*x96
      x170 = x141*x2
      x171 = x18*(0.311172355472072d0*x1*x120 + 0.00285512867044833d0*x1
     & - 0.103724118490691d0*x119 + x21 + 0.000428261234514334d0*x22 + 0
     &.0121035293512547d0*x23*x6 + 0.0175165224763768d0*x6)
      x172 = x2*x22
      x173 = 0.337338207872387d0*x172
      E = -rho_c*(x11 - x71)
      d1E(1) = -x11 - x153*x154 + x71 + 0.00274777245725523d0*x72 - x76*
     &x78 - 0.632303919304059d0*x82
      d2E(1) = x0*(0.0883792666034679d0*x1*x76 - 1.26460783860812d0*x153
     &*x34*x6*x80 - 2.0d0*x154*(-x12*(x108 - x46 + x87 + x89 - x90) - 0.
     &202642367284676d0*x132*x133*(-x37 + 1.0d0) + x152 + x158*x18*( x21
     & + x83 + x84 + x85) + x34*(-x117*(x105 + x111 + x112 + x115 - x98)
     & + x125 + 0.0003448148419047d0*x64*(x66 + x67 - 17.2274040011292d0
     &))) - x154*(mu*x132*x162*x38/x35**2 - 6.39350882283002d-5*x0**2.33
     &333333333333d0*x148 + x12*( 0.000440195901906306d0*x106*x17 + x107
     &*(-0.154912426780983d0*x0** (-0.333333333333333d0)*x48 - 75.398223
     &686155d0*x167*x49 - 4.83324359125373d0*x167*x93 + x169 - x92*(-1.3
     &6477107997868d0*x1 + 0.427594146176807d0*x6) - 1.61108119708458d0*
     &x94) - x159*x45 - x166*x24 + x166*x86 - 0.0138644007334803d0*x22*x
     &86 - 0.00107855455836145d0*x25 - 0.000586927869208408d0*x58 + 0.02
     &31073345558005d0*x88) + x134*(2.05391854581537d0*x127 + 1.43555326
     &817173d0*x129 - 2.0d0*x162*x164 + x163*x35*( 1.42618991424766d0*x1
     &27 + 1.53482843253593d0*x3)/x32**2 + x165 - x36*(0.950793276165107
     &d0*x127 + x165)) - x135*x164 - 0.000821860773992031d0*x136*x147 + 
     &0.00880391803812612d0*x137 + 0.000279369479946145d0*x143*x76 - 5.7
     &9052675961456d-5*x144 - 0.00192189252234297d0*x145 - 0.01173855738
     &41682d0*x146 + 0.00191767513931474d0*x149 + 0.0508356285485896d0*x
     &150*x161 + 0.0168588246652672d0*x151 + x155*x160 - x157*x160 - x15
     &8*x159 + x16*(0.213300270968718d0*x123*x42 + 0.00471283037863406d0
     &*x138* x76 - 0.000781467479590561d0*x139 + 0.0165511124114256d0*x1
     &40 + x155*x170 - x157*x170 + 0.857572898037877d0*x161*x42 - 0.0049
     &6533372342768d0*x17*x171 - 0.0110340749409504d0*x41 - 0.0198949583
     &613061d0*x43) + x34*(-0.0006896296838094d0*x116*x17 + x117*(5.1970
     &3611962767d0*x110 - 5.29055759308205d0*x113*x95 - 31.182216717766d
     &0*x114 - x168 + x169 - 23.3866625383245d0*x56*x65 *(-0.10697599576
     &3985d0*x1 + 0.289053642815522d0*x6) - x65*( -0.0720184066122292d0*
     &x1 + 0.0778387062718705d0*x6)) - 0.00370771454548634d0*x118 + 0.16
     &7525641071009d0*x122 + 0.899568554326366d0*x124 + x155*x173 + 5.42
     &506773629845d0*x156* x172 - x157*x173 + 0.0298137034090898d0*x17*x
     &76 - 0.0837628205355044d0*x171*x19 - 0.0837628205355044d0*x61 - 0.
     &0652591516234546d0*x63 + 0.0009195062450792d0*x68)) + x155*x78 + 2
     &58.630216499205d0*x156*x77 - x157*x78 + 1.26460783860812d0*x70 *x8
     &1 + 0.00183184830483682d0*x72 - 32.1639684429148d0*x76*x77 + 1.053
     &83986550676d0*x82 - 0.499760307959092d0*x22*x60*x79/x14**6)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = 0.0621813817393098d0*(0.132568899905202d0*x1 + 1.0d0)*log(1.0
     &d0 + 16.081984221457407d0/(2.2255694211506878d0*x1 + 5.98255043577
     &10795d0*x2 + 0.80042863499936356d0*x3 + 0.18970043257475591d0*x4))
     &
      x6 = rho_s**4/rho_c**4
      x7 = rho_s*x0
      x8 = x7 + 1.0d0
      x9 = -x7 + 1.0d0
      x10 = 1.92366105093154d0*x8**1.33333333333333d0 + 1.92366105093154
     &d0*x9** 1.33333333333333d0 - 3.84732210186307d0
      x11 = 0.0197517897025652d0*x10*(0.0690139921125583d0*x1 + 1.0d0)*(
     &-x6 + 1.0d0) *log(1.0d0 + 29.608574643216677d0/(2.2475918635776164
     &d0*x1 + 8.1574147034876408d0*x2 + 0.43009724712766434d0*x3 + 0.191
     &15125951273382d0*x4))
      x12 = 0.0310906908696548d0*(0.127469618870009d0*x1 + 1.0d0)*log(1.
     &0d0 + 32.163968442914921d0/(3.8447462374472119d0*x1 + 11.120374863
     &094685d0*x2 + 1.6447337755676092d0*x3 + 0.24058712912881922d0*x4))
     &
      x13 = -x5
      x14 = x10*x6
      x15 = x14*(x12 + x13)
      x16 = mu**2
      x17 = x0**1.66666666666667d0
      x18 = 1.0d0 - rho_s**2/rho_c**2
      x19 = x0**1.0d0
      x20 = x0**1.33333333333333d0
      x21 = (0.0128480790170175d0*x1 - 0.00304892328706429d0*x19 + 0.000
     &275310793616358d0*x20 + 0.0315297404574782d0*x4 + 1.0d0)*exp (-0.4
     &6675853320810856d0*x1)
      x22 = x18*x21
      x23 = x8**0.666666666666667d0 + x9**0.666666666666667d0
      x24 = x23**3
      x25 = x1*x16/x23**2
      x26 = 9.20897059521559d0*mu*x2/x23 + 1.0d0
      x27 = mu**3
      x28 = x11 + x13
      x29 = -x15 + x28
      x30 = 1d0/x4
      x31 = 1d0/x8
      x32 = x31**0.666666666666667d0
      x33 = x1*x31**0.333333333333333d0
      x34 = 0.0244354823084343d0*x4
      x35 = 8.61370200056458d0*x8**2*(-0.0177069812999059d0*x33 + 1.0d0)
     &/(x32*(x32* x34 + 0.337569861992026d0*x33 + 1.0d0))
      x36 = 1d0/x9
      x37 = x36**0.666666666666667d0
      x38 = x1*x36**0.333333333333333d0
      x39 = 8.61370200056458d0*x9**2*(-0.0177069812999059d0*x38 + 1.0d0)
     &/(x37*(x34* x37 + 0.337569861992026d0*x38 + 1.0d0))
      x40 = x18*(x21 - 1.0d0)
      E = -rho_c*(-x11 + x15 + x5 + (0.00316102962473761d0*mu**8*x0** 2.
     &66666666666667d0*x29 - mu**6*(-0.0533250677421794d0*x0**2.0d0*( x1
     &4*(-x12 + x5) + x28) + 0.00248266686171384d0*x17*x40) - 0.00264117
     &541143784d0*mu**5*x17*x22 + mu**4*( -0.00020688890514282d0*x0**1.0
     &d0*(23.3866625383245d0*x18*( -0.240695990468967d0*x1 + 0.260148278
     &53397d0*x4)*exp( -0.33933171852197186d0*x1) + x35 + x39 - 8.613702
     &00056458d0*x8** 2.66666666666667d0 - 8.61370200056458d0*x9**2.6666
     &6666666667d0) - 0.0418814102677522d0*x19*x40 + 0.337338207872387d0
     &*x20*x29) - 0.00777267271741372d0*x24*log((18.485266912338322d0*x2
     &5 + x26 + 15.312568193831799d0*x27*x3/x24)/(8.5571394854859619d0*x
     &25 + x26 )) - x27*(0.000132058770571892d0*x17*(37.6991118430775d0*
     &x18*( -3.07073492995203d0*x1 + 0.384834731559127d0*x4)*exp( -0.192
     &30865217881402d0*x1)/x19 + x30*x35 + x30*x39) + 0.0445553741830535
     &d0*x19*x22))/(0.237113969739022d0*x16*x4 + 1.0d0)**4)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x193
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = 2.22556942115069d0*x1 + 5.98255043577108d0*x2 + 0.80042863499
     &9364d0*x3 + 0.189700432574756d0*x4
      x6 = 1.0d0 + 16.0819842214574d0/x5
      x7 = log(x6)
      x8 = 0.0621813817393098d0*x7*(0.132568899905202d0*x1 + 1.0d0)
      x9 = rho_s*x0
      x10 = x9 + 1.0d0
      x11 = -x9 + 1.0d0
      x12 = 1.92366105093154d0*x10**1.33333333333333d0 + 1.9236610509315
     &4d0*x11** 1.33333333333333d0 - 3.84732210186307d0
      x13 = rho_s**4/rho_c**4
      x14 = -x13 + 1.0d0
      x15 = 2.24759186357762d0*x1 + 8.15741470348764d0*x2 + 0.4300972471
     &27664d0*x3 + 0.191151259512734d0*x4
      x16 = 1.0d0 + 29.6085746432167d0/x15
      x17 = log(x16)
      x18 = x14*x17
      x19 = x12*x18
      x20 = 0.0197517897025652d0*x19*(0.0690139921125583d0*x1 + 1.0d0)
      x21 = 3.84474623744721d0*x1 + 11.1203748630947d0*x2 + 1.6447337755
     &6761d0*x3 + 0.240587129128819d0*x4
      x22 = 1.0d0 + 32.1639684429149d0/x21
      x23 = log(x22)
      x24 = 0.0310906908696548d0*x23*(0.127469618870009d0*x1 + 1.0d0)
      x25 = -x8
      x26 = x24 + x25
      x27 = x12*x13
      x28 = x26*x27
      x29 = mu**2
      x30 = x29*x4
      x31 = 0.237113969739022d0*x30 + 1.0d0
      x32 = x31**(-4)
      x33 = mu**5
      x34 = rho_s**2/rho_c**2
      x35 = -x34 + 1.0d0
      x36 = exp(-0.46675853320810856d0*x1)
      x37 = x0**1.0d0
      x38 = -0.00304892328706429d0*x37
      x39 = x0**1.33333333333333d0
      x40 = 0.0128480790170175d0*x1 + x38 + 0.000275310793616358d0*x39 +
     & 0.0315297404574782d0*x4 + 1.0d0
      x41 = x36*x40
      x42 = x35*x41
      x43 = x33*x42
      x44 = x0**1.66666666666667d0
      x45 = 0.00264117541143784d0*x44
      x46 = x10**0.666666666666667d0 + x11**0.666666666666667d0
      x47 = x46**3
      x48 = 0.306852819440055d0
      x49 = x46**2
      x50 = 1d0/x49
      x51 = x29*x50
      x52 = x1*x51
      x53 = 1d0/x46
      x54 = mu*x53
      x55 = 9.20897059521559d0*x2*x54 + 1.0d0
      x56 = mu**3
      x57 = x3/x47
      x58 = 18.4852669123383d0*x52 + x55 + 15.3125681938318d0*x56*x57
      x59 = x58/(8.55713948548596d0*x52 + x55)
      x60 = x48*log(x59)
      x61 = x20 + x25
      x62 = -x28 + x61
      x63 = mu**8
      x64 = x0**2.66666666666667d0
      x65 = 0.00316102962473761d0*x63*x64
      x66 = 0.0445553741830535d0*x37
      x67 = x42*x66
      x68 = -3.07073492995203d0*x1 + 0.384834731559127d0*x4
      x69 = exp(-0.19230865217881402d0*x1)
      x70 = x35*x68*x69
      x71 = 1d0/x37
      x72 = 37.6991118430775d0*x71
      x73 = x70*x72
      x74 = 1d0/x4
      x75 = x10**2
      x76 = 1d0/x10
      x77 = x76**0.333333333333333d0
      x78 = x1*x77
      x79 = x76**0.666666666666667d0
      x80 = 0.0244354823084343d0*x4
      x81 = 0.337569861992026d0*x78 + x79*x80 + 1.0d0
      x82 = 1d0/x81
      x83 = x75*x82
      x84 = (-0.0177069812999059d0*x78 + 1.0d0)/x79
      x85 = 8.61370200056458d0*x84
      x86 = x83*x85
      x87 = x11**2
      x88 = 1d0/x11
      x89 = x88**0.333333333333333d0
      x90 = x1*x89
      x91 = x88**0.666666666666667d0
      x92 = x80*x91 + 0.337569861992026d0*x90 + 1.0d0
      x93 = 1d0/x92
      x94 = x87*x93
      x95 = (-0.0177069812999059d0*x90 + 1.0d0)/x91
      x96 = 8.61370200056458d0*x95
      x97 = x94*x96
      x98 = x44*(x73 + x74*x86 + x74*x97)
      x99 = mu**6
      x100 = x41 - 1.0d0
      x101 = x100*x35
      x102 = x101*x44
      x103 = x0**2.0d0
      x104 = -x24 + x8
      x105 = x104*x27
      x106 = x105 + x61
      x107 = x103*x106
      x108 = mu**4
      x109 = -0.0418814102677522d0*x101*x37
      x110 = x0**1.0d0
      x111 = x10**2.66666666666667d0
      x112 = 8.61370200056458d0*x111
      x113 = x11**2.66666666666667d0
      x114 = 8.61370200056458d0*x113
      x115 = exp(-0.33933171852197186d0*x1)
      x116 = x115*(-0.240695990468967d0*x1 + 0.26014827853397d0*x4)
      x117 = x116*x35
      x118 = 23.3866625383245d0*x117
      x119 = 0.337338207872387d0*x39
      x120 = x108*(x109 - 0.00020688890514282d0*x110*(-x112 - x114 + x11
     &8 + x86 + x97 ) + x119*x62) - x43*x45 - 0.0253302959105844d0*x47*x
     &60 - x56*(x67 + 0.000132058770571892d0*x98) + x62*x65 - x99*( 0.00
     &248266686171384d0*x102 - 0.0533250677421794d0*x107)
      x121 = x120*x32
      x122 = x34*x41
      x123 = 0.00528235082287567d0*x33*x44
      x124 = x35*x36
      x125 = x124*(0.00428269300567249d0*x1 + x38 + 0.000367081058155144
     &d0*x39 + 0.0210198269716521d0*x4)
      x126 = x10**(-0.333333333333333d0)
      x127 = x11**(-0.333333333333333d0)
      x128 = x126 - x127
      x129 = 0.0506605918211689d0*x128*x49*x60
      x130 = x1*x54
      x131 = x128*x9
      x132 = mu*x1*x131*x50
      x133 = 6.13931373014372d0*x2
      x134 = x131*x133*x53 - 1.53482843253593d0*x2
      x135 = x3*x51
      x136 = 0.0253302959105844d0*mu*x48/x58
      x137 = -2.56488140124205d0*x10**0.333333333333333d0 + 2.5648814012
     &4205d0*x11** 0.333333333333333d0
      x138 = rho_s**5*x137/rho_c**5
      x139 = 0.00274777245725523d0*x1*x7
      x140 = 16.0819842214574d0*(0.0082433173717657d0*x1 + 0.06218138173
     &93098d0)*( 0.741856473716896d0*x1 + 0.99709173929518d0*x2 + 0.4002
     &14317499682d0*x3 + 0.126466955049837d0*x4)/(x5**2*x6)
      x141 = 0.00136314985874174d0*x1 + 0.0197517897025652d0
      x142 = x141*x17
      x143 = x137*x141*x18
      x144 = -0.000454383286247248d0*x1*x19 + 29.6085746432167d0*x12*x14
     &*x141*( 0.749197287859205d0*x1 + 1.35956911724794d0*x2 + 0.2150486
     &23563832d0*x3 + 0.127434173008489d0*x4)/(x15**2*x16) + x139 - x140
     & + 4.0d0*x142*x27 + x143*x9 + x27*( 0.00132103950518672d0*x1*x23 -
     & x139 + x140 - 32.1639684429149d0*( 0.00396311851556016d0*x1 + 0.0
     &310906908696548d0)*( 1.28158207914907d0*x1 + 1.85339581051578d0*x2
     & + 0.822366887783805d0*x3 + 0.160391419419213d0*x4)/(x21**2*x22))
      x145 = -x138*x26 + x144
      x146 = 0.0891107483661069d0*x37
      x147 = 75.398223686155d0*x68*x69*x71
      x148 = x76*x9
      x149 = 1d0/x1
      x150 = 0.0508408867489862d0*x149
      x151 = x150/x77
      x152 = x151*x83*(-x148 + 1.0d0)
      x153 = 5.74246800037638d0*x74
      x154 = x153*x83*x84
      x155 = x88*x9
      x156 = x150/x89
      x157 = x156*x94*(x155 + 1.0d0)
      x158 = 22.9698720015055d0*x74
      x159 = x158*x9
      x160 = x10*x82
      x161 = x160*x84
      x162 = x159*x161
      x163 = x153*x94*x95
      x164 = x11*x93
      x165 = x164*x95
      x166 = x159*x165
      x167 = 0.112523287330675d0*x78
      x168 = 0.0162903215389562d0*x4
      x169 = x168*x79
      x170 = x167 + x169
      x171 = x74*x85/x81**2
      x172 = x171*x75*(-x148*x167 - x148*x169 + x170)
      x173 = 0.112523287330675d0*x90
      x174 = x168*x91
      x175 = x173 + x174
      x176 = x74*x96/x92**2
      x177 = x176*x87*(x155*x173 + x155*x174 + x175)
      x178 = 0.000132058770571892d0*x44
      x179 = x100*x34
      x180 = 0.00496533372342768d0*x44
      x181 = x124*(0.0777930888680181d0*x1*x40 - 0.00214134650283625d0*x
     &1 + 0.00152446164353214d0*x37 - 0.000183540529077572d0*x39 - 0.010
     &5099134858261d0*x4)
      x182 = x104*x138 - 4.0d0*x105 + x144
      x183 = 0.0533250677421794d0*x103
      x184 = 0.0837628205355044d0*x37
      x185 = x11**1.66666666666667d0
      x186 = x10**1.66666666666667d0
      x187 = 46.773325076649d0*x116*x74
      x188 = 0.00020688890514282d0*x44
      x189 = x41*x9
      x190 = 4.0d0*rho_s**3*x12/rho_c**3
      x191 = x13*x137*x26 - x142*x190 - x143 - x190*x26
      x192 = x100*x9
      x193 = x10*x170*x171 - x11*x175*x176 + x151*x160 - x156*x164 + x15
     &8*x161 - x158 *x165
      E = -rho_c*(x121 - x20 + x28 + x8)
      d1E(1) = -0.632303919304059d0*x120*x30/x31**5 - x121 + x145 + 3.0d
     &0*x28 - x32*( -0.000410930386996015d0*x0**2.0d0*x43 - x108*(0.4497
     &84277163183d0 *x106*x39 + x109 + 0.0003448148419047d0*x110*(x112 +
     & x114 - x118 - x86 - x97) - x119*x182 + x179*x184 + x181*x184 - x1
     &88*( -23.3866625383245d0*x115*x35*x74*(0.0802319968229891d0*x1 - 0
     &.173432185689313d0*x4) - 2.64527879654103d0*x117*x149 - 15.5911083
     &58883d0*x117*x74 - x152 + x153*(x111 + x113) - x154 - x157 + x159*
     &(x185 - x186) + x162 - x163 - x166 - x172 - x177 - x187*x34)) - x1
     &22*x123 + x125*x33*x45 + x129*x9 - x136*x49*( -6.16175563744611d0*
     &x130 + 30.6251363876636d0*x131*x29*x57 + 24.6470225497844d0*x132 +
     & x134 - 7.6562840969159d0*x135 - x59*( -2.85237982849532d0*x130 + 
     &11.4095193139813d0*x132 + x134)) + 0.00440195901906306d0*x43*x44 +
     & x56*(-x122*x146 + x125*x66 - x178 *(2.41662179562687d0*x0**(-0.66
     &6666666666667d0)*x70 + x147*x34 + x152 + x154 + x157 - x162 + x163
     & + x166 + x172 + x177 + x35*x69* x72*(1.02357830998401d0*x1 - 0.25
     &6556487706084d0*x4) + x73) - 0.00693220036674015d0*x39*x42 + x67 +
     & 0.000220097950953153d0*x98 ) - 0.00842941233263362d0*x62*x63*x64 
     &+ x65*(x145 + 4.0d0*x28) - x99*(-0.0041377781028564d0*x102 + 0.106
     &650135484359d0*x107 + x179 *x180 + x180*x181 - x182*x183)) + x61
      d1E(2) = x191 - x32*(x108*(x119*x191 + x184*x192 - x188*(-x158*(-x
     &185 + x186) - x187*x9 + x193)) + x123*x189 - x129 + x136*x46*(x128
     &*( 24.6470225497844d0*x130 + x133 + 30.6251363876636d0*x135) + x59
     &*( -x126 + x127)*(11.4095193139813d0*x130 + x133)) + x191*x65 + x5
     &6* (x146*x189 - x178*(-x147*x9 + x193)) + x99*(x180*x192 + x183*x1
     &91 ))
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_PW92_ERF(rho_c, rho_s, mu, E, d1E, d2E
     &)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x578, x579, x580, x581, x582, x583
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = 0.800428634999364d0*x3
      x5 = x0**0.666666666666667d0
      x6 = 2.22556942115069d0*x1 + 5.98255043577108d0*x2 + x4 + 0.189700
     &432574756d0 *x5
      x7 = 1.0d0 + 16.0819842214574d0/x6
      x8 = log(x7)
      x9 = 0.0621813817393098d0*x8*(0.132568899905202d0*x1 + 1.0d0)
      x10 = rho_s*x0
      x11 = x10 + 1.0d0
      x12 = -x10 + 1.0d0
      x13 = 1.92366105093154d0*x11**1.33333333333333d0 + 1.9236610509315
     &4d0*x12** 1.33333333333333d0 - 3.84732210186307d0
      x14 = rho_s**4/rho_c**4
      x15 = -x14 + 1.0d0
      x16 = 0.430097247127664d0*x3
      x17 = 2.24759186357762d0*x1 + x16 + 8.15741470348764d0*x2 + 0.1911
     &51259512734d0*x5
      x18 = 1.0d0 + 29.6085746432167d0/x17
      x19 = log(x18)
      x20 = x15*x19
      x21 = x13*x20
      x22 = 0.0197517897025652d0*x21*(0.0690139921125583d0*x1 + 1.0d0)
      x23 = 1.64473377556761d0*x3
      x24 = 3.84474623744721d0*x1 + 11.1203748630947d0*x2 + x23 + 0.2405
     &87129128819d0*x5
      x25 = 1.0d0 + 32.1639684429149d0/x24
      x26 = log(x25)
      x27 = 0.0310906908696548d0*x26*(0.127469618870009d0*x1 + 1.0d0)
      x28 = -x9
      x29 = x27 + x28
      x30 = x13*x14
      x31 = x29*x30
      x32 = mu**2
      x33 = x32*x5
      x34 = 0.237113969739022d0*x33 + 1.0d0
      x35 = x34**(-4)
      x36 = rho_s**2/rho_c**2
      x37 = -x36 + 1.0d0
      x38 = x0**1.66666666666667d0
      x39 = exp(-0.46675853320810856d0*x1)
      x40 = x0**1.0d0
      x41 = -0.00304892328706429d0*x40
      x42 = x0**1.33333333333333d0
      x43 = 0.0128480790170175d0*x1 + x41 + 0.000275310793616358d0*x42 +
     & 0.0315297404574782d0*x5 + 1.0d0
      x44 = x39*x43
      x45 = x38*x44
      x46 = x37*x45
      x47 = mu**5
      x48 = 0.00264117541143784d0*x47
      x49 = x46*x48
      x50 = x11**0.666666666666667d0
      x51 = x12**0.666666666666667d0
      x52 = x50 + x51
      x53 = x52**3
      x54 = 0.306852819440055d0
      x55 = x52**2
      x56 = 1d0/x55
      x57 = x32*x56
      x58 = x1*x57
      x59 = 1d0/x52
      x60 = mu*x59
      x61 = 9.20897059521559d0*x2*x60 + 1.0d0
      x62 = 8.55713948548596d0*x58 + x61
      x63 = 1d0/x62
      x64 = mu**3
      x65 = 1d0/x53
      x66 = x3*x65
      x67 = 18.4852669123383d0*x58 + x61 + 15.3125681938318d0*x64*x66
      x68 = x63*x67
      x69 = x54*log(x68)
      x70 = 0.0253302959105844d0*x53*x69
      x71 = x22 + x28
      x72 = -x31 + x71
      x73 = mu**8*x0**2.66666666666667d0
      x74 = 0.00316102962473761d0*x73
      x75 = x72*x74
      x76 = x40*x44
      x77 = 0.0445553741830535d0*x37
      x78 = x76*x77
      x79 = -3.07073492995203d0*x1 + 0.384834731559127d0*x5
      x80 = exp(-0.19230865217881402d0*x1)
      x81 = x80/x40
      x82 = x79*x81
      x83 = 37.6991118430775d0*x37
      x84 = x82*x83
      x85 = 1d0/x5
      x86 = x11**2
      x87 = 8.61370200056458d0*x86
      x88 = 1d0/x11
      x89 = x88**0.333333333333333d0
      x90 = x1*x89
      x91 = 0.0177069812999059d0*x90
      x92 = -x91 + 1.0d0
      x93 = x88**0.666666666666667d0
      x94 = 0.0244354823084343d0*x5
      x95 = 0.337569861992026d0*x90 + x93*x94 + 1.0d0
      x96 = 1d0/x95
      x97 = 1d0/x93
      x98 = x96*x97
      x99 = x92*x98
      x100 = x87*x99
      x101 = x12**2
      x102 = 8.61370200056458d0*x101
      x103 = 1d0/x12
      x104 = x103**0.333333333333333d0
      x105 = x1*x104
      x106 = 0.0177069812999059d0*x105
      x107 = -x106 + 1.0d0
      x108 = x103**0.666666666666667d0
      x109 = 0.337569861992026d0*x105 + x108*x94 + 1.0d0
      x110 = 1d0/x109
      x111 = 1d0/x108
      x112 = x110*x111
      x113 = x107*x112
      x114 = x102*x113
      x115 = x38*(x100*x85 + x114*x85 + x84)
      x116 = x64*(0.000132058770571892d0*x115 + x78)
      x117 = mu**6
      x118 = x44 - 1.0d0
      x119 = x118*x38
      x120 = x119*x37
      x121 = x0**2.0d0
      x122 = -x27 + x9
      x123 = x122*x30
      x124 = x123 + x71
      x125 = x121*x124
      x126 = x117*(0.00248266686171384d0*x120 - 0.0533250677421794d0*x12
     &5)
      x127 = mu**4
      x128 = x118*x40
      x129 = -0.0418814102677522d0*x128*x37
      x130 = x0**1.0d0
      x131 = x11**2.66666666666667d0
      x132 = 8.61370200056458d0*x131
      x133 = x12**2.66666666666667d0
      x134 = 8.61370200056458d0*x133
      x135 = 0.240695990468967d0*x1
      x136 = 0.26014827853397d0*x5
      x137 = exp(-0.33933171852197186d0*x1)
      x138 = x137*(-x135 + x136)
      x139 = x138*x37
      x140 = 23.3866625383245d0*x139
      x141 = x130*(x100 + x114 - x132 - x134 + x140)
      x142 = 0.337338207872387d0*x42
      x143 = x127*(x129 - 0.00020688890514282d0*x141 + x142*x72)
      x144 = -x116 - x126 + x143 - x49 - x70 + x75
      x145 = x144*x35
      x146 = x45*x47
      x147 = 0.00528235082287567d0*x146
      x148 = x147*x36
      x149 = 0.00428269300567249d0*x1 + x41 + 0.000367081058155144d0*x42
     & + 0.0210198269716521d0*x5
      x150 = x38*x39
      x151 = x149*x150
      x152 = x37*x48
      x153 = x151*x152
      x154 = x46*x47
      x155 = 0.00440195901906306d0*x154
      x156 = x0**2.0d0
      x157 = x37*x47
      x158 = x156*x157
      x159 = x158*x44
      x160 = 0.000410930386996015d0*x159
      x161 = x55*x69
      x162 = x11**(-0.333333333333333d0)
      x163 = x12**(-0.333333333333333d0)
      x164 = x162 - x163
      x165 = 0.0506605918211689d0*x164
      x166 = x161*x165
      x167 = x10*x166
      x168 = x72*x73
      x169 = 0.00842941233263362d0*x168
      x170 = x1*x60
      x171 = mu*x1
      x172 = x171*x56
      x173 = x10*x164
      x174 = x172*x173
      x175 = 6.13931373014372d0*x2
      x176 = x10*x59
      x177 = x164*x176
      x178 = x175*x177 - 1.53482843253593d0*x2
      x179 = -2.85237982849532d0*x170 + 11.4095193139813d0*x174 + x178
      x180 = x3*x57
      x181 = x32*x66
      x182 = x173*x181
      x183 = -6.16175563744611d0*x170 + 24.6470225497844d0*x174 + x178 -
     & 7.6562840969159d0*x180 + 30.6251363876636d0*x182
      x184 = -x179*x68 + x183
      x185 = 0.0253302959105844d0*x54
      x186 = x184*x185
      x187 = 1d0/x67
      x188 = mu*x187
      x189 = x188*x55
      x190 = x186*x189
      x191 = rho_s**5/rho_c**5
      x192 = x12**0.333333333333333d0
      x193 = x11**0.333333333333333d0
      x194 = x192 - x193
      x195 = 2.56488140124205d0*x194
      x196 = x191*x195
      x197 = x1*x8
      x198 = 0.00274777245725523d0*x197
      x199 = x1*x26
      x200 = 0.00132103950518672d0*x199
      x201 = 1d0/x7
      x202 = x6**(-2)
      x203 = 0.741856473716896d0*x1 + 0.99709173929518d0*x2 + 0.40021431
     &7499682d0*x3 + 0.126466955049837d0*x5
      x204 = x201*x202*x203
      x205 = 0.0082433173717657d0*x1 + 0.0621813817393098d0
      x206 = 16.0819842214574d0*x205
      x207 = x204*x206
      x208 = 1d0/x25
      x209 = x24**(-2)
      x210 = 1.28158207914907d0*x1 + 1.85339581051578d0*x2 + 0.822366887
     &783805d0*x3 + 0.160391419419213d0*x5
      x211 = x208*x209*x210
      x212 = 0.00396311851556016d0*x1 + 0.0310906908696548d0
      x213 = 32.1639684429149d0*x212
      x214 = x211*x213
      x215 = -x198 + x200 + x207 - x214
      x216 = x215*x30
      x217 = x1*x21
      x218 = 0.00136314985874174d0*x1 + 0.0197517897025652d0
      x219 = x19*x218
      x220 = x219*x30
      x221 = x20*x218
      x222 = x195*x221
      x223 = 1d0/x18
      x224 = x17**(-2)
      x225 = 0.749197287859205d0*x1 + 1.35956911724794d0*x2 + 0.21504862
     &3563832d0*x3 + 0.127434173008489d0*x5
      x226 = x223*x224*x225
      x227 = x13*x15*x218
      x228 = 29.6085746432167d0*x227
      x229 = x198 - x207
      x230 = x10*x222 + x216 - 0.000454383286247248d0*x217 + 4.0d0*x220 
     &+ x226*x228 + x229
      x231 = -x196*x29 + x230
      x232 = x231 + 4.0d0*x31
      x233 = x232*x74
      x234 = 0.0891107483661069d0*x76
      x235 = x39*x40
      x236 = x235*x77
      x237 = x37*x42
      x238 = x237*x44
      x239 = 1.02357830998401d0*x1 - 0.256556487706084d0*x5
      x240 = x81*x83
      x241 = x0**(-0.666666666666667d0)
      x242 = x37*x79*x80
      x243 = x241*x242
      x244 = 75.398223686155d0*x82
      x245 = x244*x36
      x246 = 5.74246800037638d0*x85
      x247 = x86*x99
      x248 = x246*x247
      x249 = x101*x113
      x250 = x246*x249
      x251 = 1d0/x1
      x252 = 0.0508408867489862d0*x251
      x253 = 1d0/x89
      x254 = x253*x96
      x255 = x10*x88
      x256 = -x255 + 1.0d0
      x257 = x256*x86
      x258 = x254*x257
      x259 = x252*x258
      x260 = 1d0/x104
      x261 = x110*x260
      x262 = x10*x103
      x263 = x262 + 1.0d0
      x264 = x101*x263
      x265 = x261*x264
      x266 = x252*x265
      x267 = 0.112523287330675d0*x90
      x268 = 0.0162903215389562d0*x5
      x269 = x268*x93
      x270 = x267 + x269
      x271 = -x255*x267 - x255*x269 + x270
      x272 = x95**(-2)
      x273 = x272*x97
      x274 = x273*x92
      x275 = x271*x274
      x276 = x275*x85
      x277 = x276*x87
      x278 = 0.112523287330675d0*x105
      x279 = x108*x268
      x280 = x278 + x279
      x281 = x262*x278 + x262*x279 + x280
      x282 = x109**(-2)
      x283 = x111*x282
      x284 = x107*x283
      x285 = x281*x284*x85
      x286 = x102*x285
      x287 = 22.9698720015055d0*x85
      x288 = x10*x287
      x289 = x113*x12
      x290 = x288*x289
      x291 = x11*x99
      x292 = x288*x291
      x293 = x248 + x250 + x259 + x266 + x277 + x286 + x290 - x292
      x294 = x239*x240 + 2.41662179562687d0*x243 + x245 + x293 + x84
      x295 = 0.000132058770571892d0*x38
      x296 = x64*(0.000220097950953153d0*x115 + x149*x236 - x234*x36 - 0
     &.00693220036674015d0*x238 - x294*x295 + x78)
      x297 = 0.00496533372342768d0*x119
      x298 = x1*x43
      x299 = 0.00214134650283625d0*x1
      x300 = 0.0105099134858261d0*x5
      x301 = 0.00152446164353214d0*x40
      x302 = 0.000183540529077572d0*x42
      x303 = -x299 - x300 + x301 - x302
      x304 = 0.0777930888680181d0*x298 + x303
      x305 = x304*x37
      x306 = 0.00496533372342768d0*x150
      x307 = x122*x196 - 4.0d0*x123 + x230
      x308 = 0.0533250677421794d0*x121
      x309 = x117*(-0.0041377781028564d0*x120 + 0.106650135484359d0*x125
     & + x297*x36 + x305*x306 - x307*x308)
      x310 = 0.0837628205355044d0*x128
      x311 = x235*x305
      x312 = 0.449784277163183d0*x42
      x313 = x85*(x131 + x133)
      x314 = 5.74246800037638d0*x313
      x315 = x12**1.66666666666667d0
      x316 = x11**1.66666666666667d0
      x317 = x315 - x316
      x318 = x288*x317
      x319 = x138*x85
      x320 = 46.773325076649d0*x319
      x321 = x320*x36
      x322 = x139*x85
      x323 = 15.591108358883d0*x322
      x324 = 2.64527879654103d0*x139*x251
      x325 = x137*x85
      x326 = 0.0802319968229891d0*x1 - 0.173432185689313d0*x5
      x327 = x326*x37
      x328 = x325*x327
      x329 = 23.3866625383245d0*x328
      x330 = 0.00020688890514282d0*x38
      x331 = x127*(x124*x312 + x129 + 0.0003448148419047d0*x130*(-x100 -
     & x114 + x132 + x134 - x140) - x142*x307 + x310*x36 + 0.08376282053
     &55044d0*x311 - x330*(-x248 - x250 - x259 - x266 - x277 - x286 - x2
     &90 + x292 + x314 + x318 - x321 - x323 - x324 - x329))
      x332 = -x148 + x153 + x155 - x160 + x167 - x169 - x190 + x233 + x2
     &96 - x309 - x331
      x333 = x33/x34**5
      x334 = x144*x333
      x335 = rho_s**3/rho_c**3
      x336 = 4.0d0*x335
      x337 = x13*x336
      x338 = -x222
      x339 = x14*x29
      x340 = x13*x29
      x341 = x195*x339 - x219*x337 - x336*x340 + x338
      x342 = 24.6470225497844d0*x170 + x175 + 30.6251363876636d0*x180
      x343 = x164*x342
      x344 = -x162 + x163
      x345 = 11.4095193139813d0*x170 + x175
      x346 = x344*x345
      x347 = x343 + x346*x68
      x348 = x188*x347
      x349 = x185*x52
      x350 = x10*x244
      x351 = x11*x254
      x352 = x12*x261
      x353 = 8.61370200056458d0*x85
      x354 = x270*x274
      x355 = x11*x354
      x356 = x280*x284
      x357 = x12*x356
      x358 = x252*x351 - x252*x352 - x287*x289 + x287*x291 + x353*x355 -
     & x353*x357
      x359 = -x350 + x358
      x360 = -x10*x320 - x287*(-x315 + x316) + x358
      x361 = x10*x147 + x117*(x10*x297 + x308*x341) + x127*(x10*x310 + x
     &142*x341 - x330*x360) - x166 + x341*x74 + x348*x349 + x64*(x10*x23
     &4 - x295* x359)
      x362 = -x35*x361
      x363 = x36*x47
      x364 = x156*x44
      x365 = x149*x39
      x366 = 0.00571025734089666d0*x1 - 0.00609784657412857d0*x40 + 0.00
     &0856522469028669d0*x42 + 0.0350330449527535d0*x5
      x367 = x164**2
      x368 = 0.0675474557615585d0*x367*x52*x69
      x369 = 0.0168868639403896d0*x161
      x370 = x11**(-1.33333333333333d0)
      x371 = x10*x370
      x372 = x12**(-1.33333333333333d0)
      x373 = x10*x372
      x374 = x371 + x373
      x375 = -6.0d0*x162 + 6.0d0*x163 + x374
      x376 = x10*x375
      x377 = 0.101321183642338d0*x164
      x378 = x184*x188*x52*x54
      x379 = x186*x32
      x380 = x187*x63
      x381 = x379*x380
      x382 = x67**(-2)
      x383 = x379*x382
      x384 = 3.80317310466043d0*x170
      x385 = x172*x376
      x386 = x36*x367
      x387 = x171*x65
      x388 = 2.04643791004791d0*x2
      x389 = 8.18575164019163d0*x2
      x390 = x176*x375*x388 - x177*x388 + x386*x389*x56
      x391 = 8.21567418326147d0*x170
      x392 = x179*x60
      x393 = 2.0d0*x63
      x394 = 12.2786274602874d0*x2
      x395 = 22.8190386279626d0*x172
      x396 = x173*x395
      x397 = x67/x62**2
      x398 = x392*x397
      x399 = 49.2940450995689d0*x367
      x400 = x191*x194
      x401 = 20.5190512099364d0*x400
      x402 = 0.00366369660967365d0*x197
      x403 = 0.0883792666034679d0*x1*x204
      x404 = x201*x206
      x405 = x202*x404*(0.989141964955861d0*x1 + 1.16327369584438d0*x2 +
     & 0.600321476249523d0*x3 + 0.210778258416395d0*x5)
      x406 = x208*x213
      x407 = 258.630216499205d0*x203**2*x205/(x6**4*x7**2)
      x408 = x203*x404*(1.48371294743379d0*x1 + 1.99418347859036d0*x2 + 
     &x4 + 0.252933910099675d0*x5)/x6**3
      x409 = x30*(0.0849797459133392d0*x1*x211 - 0.00176138600691563d0*x
     &199 + x209* x406*(1.70877610553209d0*x1 + 2.16229511226841d0*x2 + 
     &1.23355033167571d0*x3 + 0.267319032365355d0*x5) + 1034.52086599683
     &d0*x210**2*x212/(x24**4*x25**2) - x210*x406*( 2.56316415829814d0*x
     &1 + 3.70679162103156d0*x2 + x23 + 0.320782838838426d0*x5)/x24**3 +
     & x402 - x403 - x405 - x407 + x408 )
      x410 = x11**(-0.666666666666667d0)
      x411 = 0.854960467080683d0*x10
      x412 = x12**(-0.666666666666667d0)
      x413 = x410*x411 + x411*x412
      x414 = -5.1297628024841d0*x192 + 5.1297628024841d0*x193 + x413
      x415 = x191*x414
      x416 = -8.0d0*x216
      x417 = 5.1297628024841d0*x215*x400
      x418 = x219*x401
      x419 = x1*x19
      x420 = -0.00363506628997799d0*x30*x419
      x421 = -x408
      x422 = x10*x221
      x423 = x414*x422
      x424 = x1*x194*x20
      x425 = -0.00233087847986162d0*x10*x424
      x426 = x15*x226
      x427 = -0.0269072828949635d0*x1*x13*x426
      x428 = x223*x228
      x429 = -x224*x428*(0.998929717145607d0*x1 + 1.5861639701226d0*x2 +
     & 0.322572935345748d0*x3 + 0.212390288347482d0*x5)
      x430 = -876.667692402934d0*x225**2*x227/(x17**4*x18**2)
      x431 = x225*x428*(1.49839457571841d0*x1 + x16 + 2.71913823449588d0
     &*x2 + 0.254868346016978d0*x5)/x17**3
      x432 = x218*x226
      x433 = 236.868597145733d0*x30*x432
      x434 = x194*x218*x426
      x435 = 151.884964839347d0*x10*x434
      x436 = -x29*x415 + x403 + x405 + x407 + x409 + x416 + x417 + x418 
     &+ x420 + x421 + x423 + x425 + x427 + x429 + x430 + x431 + x433 + x
     &435
      x437 = 0.000605844381662998d0*x217 - 20.0d0*x220 - x402
      x438 = x29*x401 - 20.0d0*x31 + x436 + x437
      x439 = 0.00316102962473761d0*x73
      x440 = x235*x36
      x441 = x42*x44
      x442 = x235*x37
      x443 = 0.0891107483661069d0*x149
      x444 = x239*x37
      x445 = 75.398223686155d0*x81
      x446 = x241*x80
      x447 = 4.83324359125373d0*x446
      x448 = 1d0/x86
      x449 = 0.00786976946662483d0*x36
      x450 = x251*x254
      x451 = 0.0677878489986483d0*x251
      x452 = 1.91415600012546d0*x85
      x453 = 1d0/x101
      x454 = x251*x261
      x455 = 0.0543010717965206d0*x5
      x456 = 0.150031049774234d0*x90
      x457 = x36*x448
      x458 = 0.0271505358982603d0*x5
      x459 = x458*x93
      x460 = -x456 - x459
      x461 = x85*x87
      x462 = 0.101681773497972d0*x251
      x463 = x253*x272
      x464 = x271*x463
      x465 = x260*x282
      x466 = x462*x465
      x467 = 38.2831200025092d0*x85
      x468 = x467*x99
      x469 = x113*x467
      x470 = 0.150031049774234d0*x105
      x471 = x36*x453
      x472 = x108*x458
      x473 = x470 + x472
      x474 = x102*x85
      x475 = 0.225046574661351d0*x90
      x476 = 0.0325806430779124d0*x5
      x477 = x476*x93
      x478 = x475 + x477
      x479 = x95**(-3)
      x480 = x479*x92*x97
      x481 = 0.225046574661351d0*x105
      x482 = x108*x476
      x483 = x481 + x482
      x484 = x109**(-3)
      x485 = x107*x111*x484
      x486 = x263*x352
      x487 = 0.271151395994593d0*x10
      x488 = x10*x85
      x489 = 15.3132480010037d0*x488
      x490 = x251*x256*x351
      x491 = x12*x281
      x492 = 45.9397440030111d0*x85
      x493 = x10*x492
      x494 = x11*x275
      x495 = 11.4849360007528d0*x101*x285 - x102*x454*(0.015739538933249
     &7d0*x262 + x449*x453 + 0.00786976946662483d0) - x247*x452 - x249*x
     &452 + x251 *x486*x487 + x257*x462*x464 + x258*x451 + x264*x281*x46
     &6 + x265* x451 + x271*x461*x480*(-x255*x475 - x255*x477 + x478) + 
     &x274*x461 *(x255*x455*x93 + 0.300062099548467d0*x255*x90 - x456*x4
     &57 - x457 *x459 + x460) + 11.4849360007528d0*x276*x86 + x281*x474*
     &x485*( x262*x481 + x262*x482 + x483) - x284*x474*(0.30006209954846
     &7d0* x105*x262 + x108*x262*x455 + x470*x471 + x471*x472 + x473) + 
     &x284 *x491*x493 - x289*x489 + x291*x489 + x36*x468 + x36*x469 + x4
     &50* x87*(0.0157395389332497d0*x255 - x448*x449 - 0.007869769466624
     &83d0) - x487*x490 - x493*x494
      x496 = x150*x304
      x497 = 0.311172355472072d0*x1
      x498 = 0.00285512867044833d0*x1 - 0.103724118490691d0*x298 + x41 +
     & 0.000428261234514334d0*x42 + 0.0121035293512547d0*x43*x5 + 0.0175
     &165224763768d0*x5
      x499 = 5.29055759308205d0*x137*x251
      x500 = 5.0d0*x10
      x501 = x50*x500 + x500*x51
      x502 = 7.65662400050184d0*x85
      x503 = x317*x85
      x504 = x325*x326
      x505 = 1.26460783860812d0*x333
      x506 = x10*x47
      x507 = 3.0d0*x162
      x508 = 3.0d0*x163
      x509 = x374 - x507 + x508
      x510 = x348*x54
      x511 = x342*x344
      x512 = x164*x345
      x513 = x335*x340
      x514 = x194*x339
      x515 = x13*x335
      x516 = x219*x515
      x517 = -2.56488140124205d0*x192 + 2.56488140124205d0*x193 + x413
      x518 = 20.5190512099364d0*x194
      x519 = -x14*x195*x215 - x14*x219*x518 + x215*x337 - x221*x517 + x3
     &39*x517 + 0.00181753314498899d0*x419*x515 + 0.00116543923993081d0*
     &x424 - 118.434298572867d0*x432*x515 - 75.9424824196734d0*x434
      x520 = 16.0d0*x513 - 20.5190512099364d0*x514 + 16.0d0*x516 + x519
      x521 = 1.02321895502395d0*x2
      x522 = x389*x59
      x523 = x367*x522
      x524 = x172*x399
      x525 = x181*x367
      x526 = 10.2083787958879d0*x180
      x527 = 22.8190386279626d0*x170 + x394
      x528 = -x371 - x373 + x507 - x508
      x529 = x188*x349
      x530 = x10*x235
      x531 = 0.135575697997297d0*x10
      x532 = x450*x531
      x533 = 0.0338939244993242d0*x251
      x534 = x351*x533
      x535 = x454*x531
      x536 = x352*x533
      x537 = x351*x451
      x538 = 38.2831200025092d0*x488
      x539 = 7.65662400050185d0*x85
      x540 = x451*x486
      x541 = x270*x463
      x542 = x11*x252
      x543 = x256*x541*x542
      x544 = x252*x465
      x545 = x12*x263*x280*x544
      x546 = x464*x542
      x547 = x491*x544
      x548 = x255*x456
      x549 = x255*x459
      x550 = x274*x353
      x551 = x287*x491
      x552 = x262*x470 + x262*x472 + x473
      x553 = x284*x353
      x554 = x353*x478
      x555 = x480*x554
      x556 = x353*x483
      x557 = x485*x556
      x558 = x10*(x135 - x136)
      x559 = x91 - 1.0d0
      x560 = x11*x559
      x561 = x106 - 1.0d0
      x562 = x12*x561
      x563 = x273*x560
      x564 = x283*x561
      x565 = x283*x562
      x566 = x456 + x459
      x567 = x335*x518
      x568 = x219*x567
      x569 = 0.854960467080683d0*x410 + 0.854960467080683d0*x412
      x570 = x221*x569
      x571 = 12.0d0*x36
      x572 = x13*x219*x571
      x573 = x339*x569
      x574 = x29*x567
      x575 = x340*x571
      x576 = x370 + x372
      x577 = x185*x32*x347*x59
      x578 = x568 + x570 - x572 - x573 + x574 - x575
      x579 = x388*x576
      x580 = mu*x56
      x581 = x344**2
      x582 = 0.203363546995945d0*x251
      x583 = x254*x582 + x261*x582 + x270*x555 + x280*x466 + x280*x557 +
     & x354*x492 + x356*x492 + x462*x541 + x468 + x469 - x473*x553 - x55
     &0*x566
      E = -rho_c*(x145 - x22 + x31 + x9)
      d1E(1) = -x145 + x231 + 3.0d0*x31 - x332*x35 - 0.632303919304059d0
     &*x334 + x71
      d1E(2) = x341 + x362
      d2E(1) = x0*(-0.499760307959092d0*x127*x144*x42/x34**6 + 5.1297628
     &024841d0*x194* x422 + 0.00183184830483682d0*x197 - 32.163968442914
     &8d0*x204*x205 - 0.000302922190831499d0*x217 - 12.0d0*x220 + 59.217
     &1492864334d0* x226*x227 + 15.3892884074523d0*x29*x400 - 2.0d0*x30*
     &(-x200 + x214 + x229) - 12.0d0*x31 - x332*x505 + 1.05383986550676d
     &0*x334 + 2.0d0*x35*(x148 - x153 - x155 + x160 - x167 + x169 + x190
     & - x233 - x296 + x309 + x331) - x35*(-6.39350882283002d-5*x0** 2.3
     &3333333333333d0*x157*x44 + x10*x377*x378 - x117*( -0.0314471135817
     &087d0*x119*x36 + 0.0110340749409504d0*x120 + 0.213300270968718d0*x
     &121*x307 - 0.319950406453076d0*x125 - 0.0165511124114256d0*x150*x3
     &05 + x306*x37*(x303*x497 + x498) - x308*(-x122*x401 + x122*x415 + 
     &20.0d0*x123 + x403 + x405 + x407 + x409 + x416 + x417 + x418 + x42
     &0 + x421 + x423 + x425 + x427 + x429 + x430 + x431 + x433 + x435 +
     & x437) + 0.0198613348937107d0* x36*x496) + x127*(0.418814102677522
     &d0*x128*x36 - 0.0009195062450792d0*x141 + x142*x438 - 0.8995685543
     &26366d0*x232* x42 - 0.335051282142018d0*x304*x440 - x310*x37 + 0.1
     &67525641071009d0*x311 - x330*(-x10*x502*(-6.0d0*x315 + 6.0d0* x316
     & + x501) - 30.6264960020074d0*x10*x503 + 10.5811151861641d0* x138*
     &x251*x36 + 0.299209d0*x139 + 1.91415600012546d0*x313 - 77.95554179
     &44151d0*x319*x36 - 5.19703611962767d0*x322 + 23.3866625383245d0*x3
     &25*x37*(-0.106975995763985d0*x1 + 0.289053642815522d0*x5) + x327*x
     &499 + 31.182216717766d0*x328 + 93.5466501532981d0*x36*x504 + x495)
     & + 0.0006896296838094d0*x38*( x293 - x314 - x318 + x321 + x323 + x
     &324 + x329) + 1.04949664671409d0*x42*x72 - 0.0837628205355044d0*x4
     &42*(-x497*( x299 + x300 - x301 + x302) + x498)) + 0.03345488854487
     &92d0*x146* x36 - x150*x152*x366 - 0.00880391803812612d0*x151*x157 
     &+ 0.0105647016457513d0*x151*x363 - 0.0117385573841682d0*x154 + 0.0
     &00821860773992031d0*x158*x365 + 0.00191767513931474d0*x159 + 0.030
     &9078452196566d0*x168 - x179*x381*x52 + x183*x383*x52 - x185* x189*
     &(-16.4313483665229d0*x174 + 11.4844261453738d0*x180 + 10.208378795
     &8879d0*x181*x376 - 30.6251363876636d0*x182 - x183* x392*x393 + 1.7
     &9063317129192d0*x2 + 81.6670303671029d0*x3*x32* x386/x52**4 + x36*
     &x387*x399 + 8.21567418326147d0*x385 + x390 + x391 + x398*(-5.70475
     &965699064d0*x170 + x177*x394 - 3.06965686507186d0*x2 + x396) - x68
     &*(-7.60634620932085d0*x174 + 1.79063317129192d0*x2 + x384 + 3.8031
     &7310466043d0*x385 + 22.8190386279626d0*x386*x387 + x390)) - 0.0168
     &588246652672d0*x232 *x73 - x36*x368 - 0.00164372154798406d0*x363*x
     &364 + x369*x376 + x438*x439 + x64*(-0.000586927869208408d0*x115 + 
     &0.178221496732214d0*x149*x440 - x234*x37 - x236*x366 + 0.013864400
     &7334803d0*x237*x365 + 0.0231073345558005d0*x238 + 0.00044019590190
     &6306d0*x294*x38 - x295*(0.154912426780983d0*x0**( -0.3333333333333
     &33d0)*x242 + 150.79644737231d0*x239*x36*x81 + x240*(-1.36477107997
     &868d0*x1 + 0.427594146176807d0*x5) + 1.61108119708458d0*x243 - x24
     &5 + 9.66648718250747d0*x36*x446*x79 + x444*x445 + x444*x447 + x495
     &) - 0.0277288014669606d0*x36*x441 + 0.445553741830535d0*x36*x76 - 
     &x442*x443 - 0.00107855455836145d0* x46)) + x436 + x505*(x116 + x12
     &6 - x143 + x49 + x70 - x75))
      d2E(2) = x0*(-0.632303919304059d0*x333*x361 + x338 - x35*(-0.01936
     &86196838775d0* x10*x146 - x10*x165*x510 + x10*x368 + x117*(-0.0182
     &062236525682d0 *x10*x119 + 0.00993066744685536d0*x10*x496 - 0.1066
     &50135484359d0* x121*x341 + x308*x520) + x127*(-0.251288461606513d0
     &*x10*x128 + x142*x520 + 0.167525641071009d0*x304*x530 - x312*x341 
     &+ x330*( 46.773325076649d0*x10*x504 - x111*x484*x491*x556*x561 - x
     &112*x538 *x561 + x112*x539*x562 + x246*x270*x563 - x246*x280*x565 
     &- x256* x537 - x270*x273*x288*x559 + x271*x287*x563 + x271*x479*x5
     &54*x560 *x97 - x280*x288*x564 + 62.364433435532d0*x325*x558 + x353
     &*x552* x565 - x353*x563*(-x548 - x549 + x566) - x499*x558 - x502*(
     &-3.0d0 *x315 + 3.0d0*x316 + x501) - 15.3132480010037d0*x503 + x532
     & - x534 + x535 + x536 - x538*x559*x98 - x539*x560*x98 + x540 - x54
     &3 + x545 - x546 + x547 - x551*x564) + 0.0003448148419047d0*x360*x3
     &8 ) - 0.00528235082287567d0*x151*x506 - x165*x378 - 0.008429412332
     &63362d0*x341*x73 + 0.000821860773992031d0*x364*x506 - x369*x509 + 
     &x381*x512 + x383*x511 + x439*x520 - x529*(-x10*x523 - x10*x524 - 8
     &1.6670303671029d0*x10*x525 + 15.3125681938318d0* x164*x180 + x164*
     &x391 + x164*x521 - x183*x346*x60*x63 + x343*x392 *x63 + x344*x398*
     &x527 - x388*x509 - x391*x509 - x509*x526 - x68*( x173*x344*x522 - 
     &x344*x384 + x344*x396 - x344*x521 + x384*x528 + x388*x528)) + x64*
     &(0.0138644007334803d0*x10*x441 - 0.267332245098321d0*x10*x76 - x29
     &5*(-x10*x239*x445 - x10*x447*x79 + x11*x271*x555 + x11*x550*(x460 
     &+ x548 + x549) - x113*x538 + x12 *x552*x553 + x246*x355 - x246*x35
     &7 - x284*x551 + x287*x494 - x288 *x354 - x288*x356 + x289*x539 - x
     &291*x539 + x350 + 0.135575697997297d0*x490 - x491*x557 - x532 + x5
     &34 - x535 - x536 + x537*(x255 - 1.0d0) - x538*x99 - x540 + x543 - 
     &x545 + x546 - x547) + 0.000220097950953153d0*x359*x38 - x443*x530)
     &) + x362 + 12.0d0*x513 - 17.9541698086943d0*x514 + 12.0d0*x516 + x
     &519)
      d2E(3) = -x0*(x35*(x117*(x297 + x308*x578) + x127*(x142*x578 + x31
     &0 - x330*(-x320 - x467*x52 + x583)) + x147 - x368 + x369*x576 + x3
     &77*x510 - x380* x512*x577 - x382*x511*x577 + x439*x578 - x529*(8.2
     &1567418326148d0 *x170*x576 + x343*x346*x393*x580 + x345*x397*x527*
     &x580*x581 + x523 + x524 + 81.6670303671029d0*x525 + x526*x576 + x5
     &79 - x68*( x384*x576 + x395*x581 + x522*x581 + x579)) + x64*(x234 
     &- x295*( -x244 + x583))) - x568 - x570 + x572 + x573 - x574 + x575
     &)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = x0**0.666666666666667d0
      x6 = log(1.0d0 + 16.081984221457407d0/(2.2255694211506878d0*x1 + 5
     &.9825504357710795d0*x3 + 0.80042863499936356d0*x4 + 0.189700432574
     &75591d0*x5))
      x7 = x2*x6
      x8 = 0.0621813817393098d0*x7
      x9 = mu**2
      x10 = (0.237113969739022d0*x5*x9 + 1.0d0)**(-4)
      x11 = x0**1.66666666666667d0
      x12 = x0**1.0d0
      x13 = x0**1.33333333333333d0
      x14 = (0.0128480790170175d0*x1 - 0.00304892328706429d0*x12 + 0.000
     &275310793616358d0*x13 + 0.0315297404574782d0*x5 + 1.0d0)*exp (-0.4
     &6675853320810856d0*x1)
      x15 = mu**5*x11*x14
      x16 = 0.00264117541143784d0*x15
      x17 = mu**8*x0**2.66666666666667d0*x7
      x18 = 0.000196557189785076d0*x17
      x19 = x1*x9
      x20 = 4.60448529760779d0*mu*x3 + 1.0d0
      x21 = mu**3
      x22 = -0.306852819440055d0*log((4.6213167280845804d0*x19 + x20 + 1
     &.9140710242289749d0*x21*x4)/(2.1392848713714905d0*x19 + x20))
      x23 = 0.202642367284676d0*x22
      x24 = x14 - 1.0d0
      x25 = mu**6*(0.00331582639355101d0*x0**2.0d0*x7 + 0.00248266686171
     &384d0*x11* x24)
      x26 = 1d0/(0.337569861992026d0*x1 + 0.0244354823084343d0*x5 + 1.0d
     &0)
      x27 = x21*(0.000132058770571892d0*x11*(17.2274040011292d0*x26*( -0
     &.0177069812999059d0*x1 + 1.0d0)/x5 + 37.6991118430775d0*( -3.07073
     &492995203d0*x1 + 0.384834731559127d0*x5)*exp( -0.19230865217881402
     &d0*x1)/x12) + 0.0445553741830535d0*x12*x14)
      x28 = mu**4*(0.00020688890514282d0*x0**1.0d0*(-x26*(-0.30504532049
     &3918d0*x1 + 17.2274040011292d0) - (-5.62907590342551d0*x1 + 6.084d
     &0*x5)*exp( -0.33933171852197186d0*x1) + 17.2274040011292d0) - 0.04
     &18814102677522d0*x12*x24 - 0.0209761558789675d0*x13*x7)
      x29 = exp(-x10*(0.08495068258568908d0*x15 + 0.0063220592494752145d
     &0*x17 - 6.5177827065418574d0*x22 + 32.163968442914815d0*x25 + 32.1
     &63968442914815d0*x27 - 32.163968442914815d0*x28) + 2.0d0*x7) - 1.0
     &d0
      x30 = (x10*(-x16 - x18 + x23 - x25 - x27 + x28) + x8)/(x2*x6)
      x31 = 307.49539002419d0*gamma_cc*rho_c**(-2.33333333333333d0)*x30*
     &*2.78d0
      x32 = 1.0d0 + x31/x29
      E = -rho_c*(-x10*(x16 + x18 - x23 + x25 + x27 - x28) + x8 - 0.0310
     &906908696549d0*log(x31*x32/(94553.414886128725d0*gamma_cc** 2*rho_
     &c**(-4.666666666666667d0)*x30**5.5599999999999996d0/x29**2 + x32) 
     &+ 1.0d0))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x138, x139, x140, x141
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = x0**0.666666666666667d0
      x6 = 2.22556942115069d0*x1 + 5.98255043577108d0*x3 + 0.80042863499
     &9364d0*x4 + 0.189700432574756d0*x5
      x7 = 1.0d0 + 16.0819842214574d0/x6
      x8 = log(x7)
      x9 = x2*x8
      x10 = 2.0d0*x9
      x11 = mu**2
      x12 = x11*x5
      x13 = 0.237113969739022d0*x12 + 1.0d0
      x14 = x13**(-4)
      x15 = mu**4
      x16 = x0**1.0d0
      x17 = exp(-0.46675853320810856d0*x1)
      x18 = 0.00304892328706429d0*x16
      x19 = x0**1.33333333333333d0
      x20 = 0.0128480790170175d0*x1 - x18 + 0.000275310793616358d0*x19 +
     & 0.0315297404574782d0*x5 + 1.0d0
      x21 = x17*x20
      x22 = x21 - 1.0d0
      x23 = 0.0418814102677522d0*x16*x22
      x24 = x19*x9
      x25 = 0.0209761558789675d0*x24
      x26 = exp(-0.33933171852197186d0*x1)
      x27 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x5 + 1.0d0
      x28 = 1d0/x27
      x29 = x0**1.0d0*(-x26*(-5.62907590342551d0*x1 + 6.084d0*x5) - x28*
     &( -0.305045320493918d0*x1 + 17.2274040011292d0) + 17.2274040011292
     &d0)
      x30 = 0.00020688890514282d0*x29
      x31 = x15*(-x23 - x25 + x30)
      x32 = 32.1639684429148d0*x31
      x33 = 0.693147180559945d0
      x34 = x33 - 1.0d0
      x35 = x1*x11
      x36 = 4.60448529760779d0*mu*x3 + 1.0d0
      x37 = mu**3
      x38 = 4.62131672808458d0*x35 + x36 + 1.91407102422897d0*x37*x4
      x39 = x38/(2.13928487137149d0*x35 + x36)
      x40 = log(x39)
      x41 = x34*x40
      x42 = 6.51778270654186d0*x41
      x43 = 0.0445553741830535d0*x16
      x44 = x21*x43
      x45 = x0**1.66666666666667d0
      x46 = exp(-0.19230865217881402d0*x1)
      x47 = x46*(-3.07073492995203d0*x1 + 0.384834731559127d0*x5)
      x48 = 37.6991118430775d0/x16
      x49 = x47*x48
      x50 = -0.0177069812999059d0*x1 + 1.0d0
      x51 = x28*x50
      x52 = 1d0/x5
      x53 = 17.2274040011292d0*x52
      x54 = x45*(x49 + x51*x53)
      x55 = x37*(x44 + 0.000132058770571892d0*x54)
      x56 = 32.1639684429148d0*x55
      x57 = x22*x45
      x58 = x0**2.0d0
      x59 = x58*x9
      x60 = mu**6*(0.00248266686171384d0*x57 + 0.00331582639355101d0*x59
     &)
      x61 = 32.1639684429148d0*x60
      x62 = mu**5
      x63 = x21*x45
      x64 = x62*x63
      x65 = 0.0849506825856891d0*x64
      x66 = x0**2.66666666666667d0
      x67 = x66*x9
      x68 = mu**8*x67
      x69 = 0.00632205924947521d0*x68
      x70 = x56 + x61 + x65 + x69
      x71 = exp(x10 - x14*(-x32 - x42 + x70)) - 1.0d0
      x72 = rho_c**(-4.66666666666667d0)
      x73 = 0.0621813817393098d0*x9
      x74 = 0.00264117541143784d0*x64
      x75 = 0.000196557189785076d0*x68
      x76 = 0.202642367284676d0*x41
      x77 = x31 - x55 - x60 - x74 - x75 + x76
      x78 = x14*x77 + x73
      x79 = 1d0/x2
      x80 = 1d0/x8
      x81 = x79*x80
      x82 = x78*x81
      x83 = x82**5.56d0
      x84 = x72*x83
      x85 = 94553.4148861287d0*gamma_cc**2
      x86 = x84*x85
      x87 = 1d0/x71
      x88 = x82**2.78d0
      x89 = 307.49539002419d0*rho_c**(-2.33333333333333d0)
      x90 = x88*x89
      x91 = gamma_cc*x90
      x92 = x87*x91 + 1.0d0
      x93 = 0.0310906908696549d0*log(x91*x92/(x92 + x86/x71**2) + 1.0d0)
     &
      x94 = x14*(-x31 + x55 + x60 + x74 + x75 - x76)
      x95 = x73 - x94
      x96 = x1*x8
      x97 = (0.741856473716896d0*x1 + 0.99709173929518d0*x3 + 0.40021431
     &7499682d0*x4 + 0.126466955049837d0*x5)/(x6**2*x7)
      x98 = x97*(0.0082433173717657d0*x1 + 0.0621813817393098d0)
      x99 = x12*x77/x13**5
      x100 = mu**7
      x101 = x17*(-0.00428269300567249d0*x1 + x18 - 0.000367081058155144
     &d0*x19 - 0.0210198269716521d0*x5)
      x102 = x2*x97
      x103 = 0.767414216267966d0*x3
      x104 = mu*x1
      x105 = x17*(0.0777930888680181d0*x1*x20 - 0.00214134650283625d0*x1
     & + 0.00152446164353214d0*x16 - 0.000183540529077572d0*x19 - 0.0105
     &099134858261d0*x5)
      x106 = 1d0/x1
      x107 = 0.101681773497972d0*x106*x28
      x108 = 11.4849360007528d0*x52
      x109 = x108*x51
      x110 = x50*x53*(0.112523287330675d0*x1 + 0.0162903215389562d0*x5)/
     &x27**2
      x111 = x26*(-0.240695990468967d0*x1 + 0.26014827853397d0*x5)
      x112 = mu*x14*(-0.000410930386996015d0*x0**2.0d0*x15*x21 + 8.68579
     &013942185d-6* x0**3.0d0*x100*x8 - 0.00316102962473761d0*x100*x102*
     &x66 + 0.000524152506093537d0*x100*x67 - 0.00264117541143784d0*x101
     &*x15* x45 + x11*(-x101*x43 - 0.00693220036674015d0*x19*x21 + x44 -
     & 0.000132058770571892d0*x45*(2.41662179562687d0*x0**( -0.666666666
     &666667d0)*x47 + x107 + x109 + x110 + x46*x48*( 1.02357830998401d0*
     &x1 - 0.256556487706084d0*x5) + x49) + 0.000220097950953153d0*x54) 
     &+ 0.00440195901906306d0*x15*x63 + 0.202642367284676d0*x34*(-x103 -
     & 1.54043890936153d0*x104 - 0.957035512114487d0*x11*x4 + x39*(x103 
     &+ 0.71309495712383d0*x104 ))/x38 + x37*(-0.337338207872387d0*x102*
     &x19 - 0.0837628205355044d0*x105*x16 + x23 + 0.0279682078386234d0*x
     &24 - 0.0003448148419047d0*x29 + 0.000926928636371586d0*x45*x8 + 0.
     &00020688890514282d0*x45*(-2.64527879654103d0*x106*x111 - x107 + x1
     &08 - x109 - x110 - 15.591108358883d0*x111*x52 - 23.3866625383245d0
     &*x26*x52*(0.0802319968229891d0*x1 - 0.173432185689313d0*x5))) + x6
     &2*(0.00014652515242323d0*x0** 2.33333333333333d0*x8 - 0.0533250677
     &421794d0*x102*x58 - 0.00496533372342768d0*x105*x45 + 0.00413777810
     &28564d0*x57 + 0.00663165278710203d0*x59))
      x113 = exp(x10 + x14*(x32 + x42 - x56 - x61 - x65 - x69))
      x114 = x113 - 1.0d0
      x115 = x114**(-2)
      x116 = 1d0/x114
      x117 = x116*x91
      x118 = x117 + 1.0d0
      x119 = 1d0/(x115*x86 + x118)
      x120 = x118*x119
      x121 = rho_c**(-3.33333333333333d0)*x88
      x122 = 717.489243389776d0*x121
      x123 = 1d0/x78
      x124 = x1*x78*x79
      x125 = x78*x80*x97
      x126 = 307.49539002419d0*x121*x123*(2.78d0*x112 + 0.12284718057882
     &d0*x124 - 44.7079161356516d0*x125 - 0.00763880743116955d0*x96 + 44
     &.7079161356516d0*x98 + 1.75780489566528d0*x99)
      x127 = 307.49539002419d0*x116
      x128 = x113*(32.1639684429148d0*x112 - 0.0883792666034679d0*x96 + 
     &16.0819842214574d0*x97*(0.265137799810404d0*x1 + 2.0d0) + 20.33740
     &33068271d0*x99)
      x129 = x127*x128
      x130 = gamma_cc*rho_c**(-5.66666666666667d0)*x83
      x131 = x116*x130
      x132 = exp(x10 - x14*(32.163968442914815d0*x15*(x23 + x25 - x30) +
     & 6.5177827065418574d0*x40*(-x33 + 1.0d0) + x70)) - 1.0d0
      x133 = x81*x95
      x134 = x133**5.56d0*x72
      x135 = x133**2.78d0*x89
      x136 = gamma_cc/x132
      x137 = x135*x136 + 1.0d0
      x138 = 1d0/(x137 + x134*x85/x132**2)
      x139 = x135*x137
      x140 = x138*x139
      x141 = gamma_cc*x87
      E = -rho_c*(-x93 + x95)
      d1E(1) = -rho_c*(-0.0310906908696549d0*gamma_cc*x119*(x117*x120*( 
     &189106.829772257d0*x115*x128*x130 + x121*x129 + x122 - 94553.41488
     &61287d0*x123*x131*(5.56d0*x112 + 0.245694361157641d0* x124 - 89.41
     &58322713032d0*x125 - 0.0152776148623391d0*x96 + 89.4158322713032d0
     &*x98 + 3.51560979133057d0*x99) - x126 + 441249.269468601d0*x131) -
     & x118*x122 + x118*x126 + x127*x130*( x123*(854.837184267248d0*x112
     & + 37.7749417054564d0*x124 - 13747.478109301d0*x125 - 2.3488980703
     &6716d0*x96 + 13747.478109301d0*x98 + 540.516901979027d0*x99) - x12
     &9 - 717.489243389776d0))/(x120*x91 + 1.0d0) + x0*x112 - 0.00274777
     &245725523d0*x0*x96 + 16.0819842214574d0*x0*x98 + 0.632303919304059
     &d0*x0*x99) - x73 + x93 + x94
      d1E(3) = 0.0310906908696549d0*rho_c*x138*(94553.4148861287d0*x134*
     &x136 + x139 - x140*x141*(189106.829772257d0*x141*x84 + x90))/(gamm
     &a_cc*x140 + 1.0d0)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_PBE_GWS_ERF(rho_c, gamma_cc, mu, E, d1E, d2
     &E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x435, x436, x437, x438
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = 0.800428634999364d0*x4
      x6 = x0**0.666666666666667d0
      x7 = 2.22556942115069d0*x1 + 5.98255043577108d0*x3 + x5 + 0.189700
     &432574756d0 *x6
      x8 = 1.0d0 + 16.0819842214574d0/x7
      x9 = log(x8)
      x10 = x2*x9
      x11 = 2.0d0*x10
      x12 = mu**2
      x13 = x12*x6
      x14 = 0.237113969739022d0*x13 + 1.0d0
      x15 = x14**(-4)
      x16 = mu**4
      x17 = x0**1.0d0
      x18 = exp(-0.46675853320810856d0*x1)
      x19 = 0.00304892328706429d0*x17
      x20 = -x19
      x21 = x0**1.33333333333333d0
      x22 = 0.0128480790170175d0*x1 + x20 + 0.000275310793616358d0*x21 +
     & 0.0315297404574782d0*x6 + 1.0d0
      x23 = x18*x22
      x24 = x23 - 1.0d0
      x25 = x17*x24
      x26 = 0.0418814102677522d0*x25
      x27 = x10*x21
      x28 = 0.0209761558789675d0*x27
      x29 = x0**1.0d0
      x30 = exp(-0.33933171852197186d0*x1)
      x31 = x30*(-5.62907590342551d0*x1 + 6.084d0*x6)
      x32 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x6 + 1.0d0
      x33 = 1d0/x32
      x34 = x33*(-0.305045320493918d0*x1 + 17.2274040011292d0)
      x35 = x29*(-x31 - x34 + 17.2274040011292d0)
      x36 = 0.00020688890514282d0*x35
      x37 = x16*(-x26 - x28 + x36)
      x38 = 32.1639684429148d0*x37
      x39 = 0.693147180559945d0
      x40 = x39 - 1.0d0
      x41 = x1*x12
      x42 = 4.60448529760779d0*mu*x3 + 1.0d0
      x43 = 2.13928487137149d0*x41 + x42
      x44 = 1d0/x43
      x45 = mu**3
      x46 = 1.91407102422897d0*x4*x45 + 4.62131672808458d0*x41 + x42
      x47 = x44*x46
      x48 = log(x47)
      x49 = x40*x48
      x50 = 6.51778270654186d0*x49
      x51 = 0.0445553741830535d0*x17
      x52 = x23*x51
      x53 = x0**1.66666666666667d0
      x54 = exp(-0.19230865217881402d0*x1)
      x55 = x54*(-3.07073492995203d0*x1 + 0.384834731559127d0*x6)
      x56 = 1d0/x17
      x57 = 37.6991118430775d0*x56
      x58 = x55*x57
      x59 = -0.0177069812999059d0*x1 + 1.0d0
      x60 = x33*x59
      x61 = 1d0/x6
      x62 = 17.2274040011292d0*x61
      x63 = x53*(x58 + x60*x62)
      x64 = x45*(x52 + 0.000132058770571892d0*x63)
      x65 = 32.1639684429148d0*x64
      x66 = x24*x53
      x67 = x0**2.0d0
      x68 = x10*x67
      x69 = mu**6*(0.00248266686171384d0*x66 + 0.00331582639355101d0*x68
     &)
      x70 = 32.1639684429148d0*x69
      x71 = mu**5
      x72 = x23*x53
      x73 = x71*x72
      x74 = 0.0849506825856891d0*x73
      x75 = x0**2.66666666666667d0
      x76 = x10*x75
      x77 = mu**8*x76
      x78 = 0.00632205924947521d0*x77
      x79 = x65 + x70 + x74 + x78
      x80 = exp(x11 - x15*(-x38 - x50 + x79))
      x81 = x80 - 1.0d0
      x82 = 1d0/x81
      x83 = 0.0621813817393098d0*x10
      x84 = 0.00264117541143784d0*x73
      x85 = 0.000196557189785076d0*x77
      x86 = 0.202642367284676d0*x49
      x87 = x37 - x64 - x69 - x84 - x85 + x86
      x88 = x15*x87 + x83
      x89 = 1d0/x2
      x90 = 1d0/x9
      x91 = x89*x90
      x92 = x88*x91
      x93 = x92**2.78d0
      x94 = rho_c**(-2.33333333333333d0)
      x95 = x93*x94
      x96 = 307.49539002419d0*x95
      x97 = gamma_cc*x96
      x98 = x82*x97 + 1.0d0
      x99 = x96*x98
      x100 = x81**(-2)
      x101 = rho_c**(-4.66666666666667d0)
      x102 = x92**5.56d0
      x103 = x101*x102
      x104 = gamma_cc**2
      x105 = 94553.4148861287d0*x104
      x106 = x103*x105
      x107 = x100*x106 + x98
      x108 = 1d0/x107
      x109 = gamma_cc*x108
      x110 = x109*x99 + 1.0d0
      x111 = 0.0310906908696549d0*log(x110)
      x112 = -x37 + x64 + x69 + x84 + x85 - x86
      x113 = x112*x15
      x114 = -x113 + x83
      x115 = x1*x9
      x116 = 0.00274777245725523d0*x115
      x117 = 0.0082433173717657d0*x1 + 0.0621813817393098d0
      x118 = 0.741856473716896d0*x1 + 0.99709173929518d0*x3 + 0.40021431
     &7499682d0*x4 + 0.126466955049837d0*x6
      x119 = 1d0/x8
      x120 = x119/x7**2
      x121 = x118*x120
      x122 = x117*x121
      x123 = 16.0819842214574d0*x122
      x124 = x14**(-5)
      x125 = x124*x13
      x126 = x125*x87
      x127 = 0.632303919304059d0*x126
      x128 = 0.00428269300567249d0*x1
      x129 = 0.0210198269716521d0*x6
      x130 = 0.000367081058155144d0*x21
      x131 = -x128 - x129 - x130 + x19
      x132 = x18*x51
      x133 = x131*x132
      x134 = x21*x23
      x135 = 0.00693220036674015d0*x134
      x136 = 0.000220097950953153d0*x63
      x137 = 1.02357830998401d0*x1 - 0.256556487706084d0*x6
      x138 = x54*x57
      x139 = x0**(-0.666666666666667d0)
      x140 = x139*x55
      x141 = 1d0/x1
      x142 = x141*x33
      x143 = 0.101681773497972d0*x142
      x144 = 11.4849360007528d0*x61
      x145 = x144*x60
      x146 = x32**(-2)
      x147 = 0.112523287330675d0*x1 + 0.0162903215389562d0*x6
      x148 = x146*x147
      x149 = x148*x59
      x150 = x149*x62
      x151 = x143 + x145 + x150
      x152 = x137*x138 + 2.41662179562687d0*x140 + x151 + x58
      x153 = 0.000132058770571892d0*x53
      x154 = x152*x153
      x155 = x30*(-0.240695990468967d0*x1 + 0.26014827853397d0*x6)
      x156 = x155*x61
      x157 = 15.591108358883d0*x156
      x158 = 2.64527879654103d0*x141*x155
      x159 = x30*(0.0802319968229891d0*x1 - 0.173432185689313d0*x6)
      x160 = x159*x61
      x161 = 23.3866625383245d0*x160
      x162 = -x143 + x144 - x145 - x150 - x157 - x158 - x161
      x163 = 0.00020688890514282d0*x53
      x164 = x53*x9
      x165 = x1*x22
      x166 = 0.00214134650283625d0*x1
      x167 = 0.0105099134858261d0*x6
      x168 = 0.00152446164353214d0*x17
      x169 = 0.000183540529077572d0*x21
      x170 = -x166 - x167 + x168 - x169
      x171 = x18*(0.0777930888680181d0*x165 + x170)
      x172 = x17*x171
      x173 = x121*x2
      x174 = x173*x21
      x175 = 0.000926928636371586d0*x164 - 0.0837628205355044d0*x172 - 0
     &.337338207872387d0*x174 + x26 + 0.0279682078386234d0*x27
      x176 = 0.767414216267966d0*x3
      x177 = mu*x1
      x178 = 1.54043890936153d0*x177
      x179 = x12*x4
      x180 = 0.957035512114487d0*x179
      x181 = x176 + 0.71309495712383d0*x177
      x182 = -x176 - x178 - x180 + x181*x47
      x183 = 0.202642367284676d0/x46
      x184 = x183*x40
      x185 = x182*x184
      x186 = x131*x18
      x187 = x16*x53
      x188 = 0.00264117541143784d0*x187
      x189 = x0**2.33333333333333d0
      x190 = x189*x9
      x191 = 0.00496533372342768d0*x53
      x192 = 0.0533250677421794d0*x67
      x193 = mu**7
      x194 = x0**3.0d0*x193
      x195 = x194*x9
      x196 = x193*x76
      x197 = x16*x72
      x198 = x0**2.0d0*x16
      x199 = x198*x23
      x200 = x193*x75
      x201 = x173*x200
      x202 = 8.68579013942185d-6*x195 + 0.000524152506093537d0*x196 + 0.
     &00440195901906306d0*x197 - 0.000410930386996015d0*x199 - 0.0031610
     &2962473761d0*x201 + x71*(-x171*x191 - x173*x192 + 0.00014652515242
     &323d0*x190 + 0.0041377781028564d0*x66 + 0.00663165278710203d0*x68)
     &
      x203 = x12*(-x133 - x135 + x136 - x154 + x52) + x185 - x186*x188 +
     & x202 + x45*( x162*x163 + x175 - 0.0003448148419047d0*x35)
      x204 = mu*x15
      x205 = x203*x204
      x206 = exp(x11 + x15*(x38 + x50 - x65 - x70 - x74 - x78))
      x207 = x206 - 1.0d0
      x208 = x207**(-2)
      x209 = 1d0/x207
      x210 = x209*x97 + 1.0d0
      x211 = x106*x208 + x210
      x212 = 1d0/x211
      x213 = x210*x96
      x214 = x212*x213
      x215 = rho_c**(-3.33333333333333d0)*x93
      x216 = 717.489243389776d0*x215
      x217 = 1d0/x88
      x218 = 0.00763880743116955d0*x115
      x219 = 44.7079161356516d0*x122
      x220 = 1.75780489566528d0*x126
      x221 = x1*x89
      x222 = x221*x88
      x223 = 0.12284718057882d0*x222
      x224 = x88*x90
      x225 = x121*x224
      x226 = 44.7079161356516d0*x225
      x227 = 2.78d0*x205
      x228 = -x218 + x219 + x220 + x223 - x226 + x227
      x229 = x217*x228
      x230 = 307.49539002419d0*x215
      x231 = x229*x230
      x232 = 307.49539002419d0*x209
      x233 = 0.265137799810404d0*x1 + 2.0d0
      x234 = 16.0819842214574d0*x233
      x235 = -0.0883792666034679d0*x115 + x121*x234
      x236 = 20.3374033068271d0*x126 + 32.1639684429148d0*x205 + x235
      x237 = x206*x236
      x238 = x232*x237
      x239 = x217*(-2.34889807036716d0*x115 + 13747.478109301d0*x122 + 5
     &40.516901979027d0*x126 + 854.837184267248d0*x205 + 37.774941705456
     &4d0*x222 - 13747.478109301d0*x225)
      x240 = -x238 + x239 - 717.489243389776d0
      x241 = gamma_cc*x240
      x242 = rho_c**(-5.66666666666667d0)*x102
      x243 = x232*x242
      x244 = gamma_cc*x209
      x245 = x242*x244
      x246 = 189106.829772257d0*gamma_cc
      x247 = x208*x242
      x248 = x237*x247
      x249 = 0.0152776148623391d0*x115
      x250 = 89.4158322713032d0*x122
      x251 = 3.51560979133057d0*x126
      x252 = 0.245694361157641d0*x222
      x253 = 89.4158322713032d0*x225
      x254 = 5.56d0*x205
      x255 = -x249 + x250 + x251 + x252 - x253 + x254
      x256 = x217*x255
      x257 = x245*x256
      x258 = x215*x238 + x216 - x231 + 441249.269468601d0*x245 + x246*x2
     &48 - 94553.4148861287d0*x257
      x259 = x214*x244
      x260 = -x210*x216 + x210*x231 + x241*x243 + x258*x259
      x261 = 0.0310906908696549d0*x212
      x262 = gamma_cc*x261
      x263 = -x39 + 1.0d0
      x264 = exp(x11 - x15*(32.163968442914815d0*x16*(x26 + x28 - x36) +
     & 6.5177827065418574d0*x263*x48 + x79)) - 1.0d0
      x265 = x114*x91
      x266 = x101*x265**5.56d0
      x267 = 307.49539002419d0*x265**2.78d0*x94
      x268 = gamma_cc/x264
      x269 = x267*x268 + 1.0d0
      x270 = 1d0/(x105*x266/x264**2 + x269)
      x271 = x267*x269
      x272 = gamma_cc*x82
      x273 = 189106.829772257d0*x103
      x274 = x270*x271
      x275 = x270*(94553.4148861287d0*x266*x268 + x271 - x272*x274*(x272
     &*x273 + x96))
      x276 = 1d0/(gamma_cc*x274 + 1.0d0)
      x277 = 0.0310906908696549d0*x276
      x278 = x112*x125
      x279 = x29*(x31 + x34 - 17.2274040011292d0)
      x280 = -x144 + x151 + x157 + x158 + x161
      x281 = x183*x263
      x282 = x182*x281
      x283 = x18*(x128 + x129 + x130 + x20)
      x284 = -x12*(x133 + x135 - x136 + x154 - x52) + x188*x283 + x202 -
     & x282 + x45*( -x163*x280 + x175 + 0.0003448148419047d0*x279)
      x285 = x204*x284
      x286 = x238 - x239 + 717.489243389776d0
      x287 = gamma_cc*x243*x286
      x288 = x216*x98
      x289 = x217*(x218 - x219 - x220 - x223 + x226 - x227)
      x290 = x230*x289
      x291 = x290*x98
      x292 = x242*x272
      x293 = x217*x292*(x249 - x250 - x251 - x252 + x253 - x254)
      x294 = x80*(x235 - 20.3374033068271d0*x278 + 32.1639684429148d0*x2
     &85)
      x295 = x100*x242*x294
      x296 = x216 + x230*x294*x82 + x290
      x297 = x246*x295 + 441249.269468601d0*x292 + 94553.4148861287d0*x2
     &93 + x296
      x298 = x108*x99
      x299 = x272*x298
      x300 = x297*x299
      x301 = -x288 - x291 + x300
      x302 = -x287 + x301
      x303 = 1d0/x110
      x304 = x108*x303
      x305 = rho_c**(-2)
      x306 = x1*x121
      x307 = x120*(0.989141964955861d0*x1 + 1.16327369584438d0*x3 + 0.60
     &0321476249523d0*x4 + 0.210778258416395d0*x6)
      x308 = x117*x307
      x309 = 16.0819842214574d0*x305
      x310 = x118**2/(x7**4*x8**2)
      x311 = x117*x310
      x312 = x118*x119*(1.48371294743379d0*x1 + 1.99418347859036d0*x3 + 
     &x5 + 0.252933910099675d0*x6)/x7**3
      x313 = x117*x312
      x314 = x16*x21/x14**6
      x315 = x124*x45*x6
      x316 = 0.0006896296838094d0*x53
      x317 = 3.82831200025092d0*x61
      x318 = x30*(-0.0720184066122292d0*x1 + 0.0778387062718705d0*x6)
      x319 = 1.38777878078145d-17*x142
      x320 = 5.19703611962767d0*x156
      x321 = 23.3866625383245d0*x30*x61*(-0.106975995763985d0*x1 + 0.289
     &053642815522d0*x6)
      x322 = 31.182216717766d0*x160
      x323 = 5.29055759308205d0*x141*x159
      x324 = x317*x60
      x325 = 0.203363546995945d0*x141*x148
      x326 = 22.9698720015055d0*x149*x61
      x327 = x59*x62
      x328 = x146*x327*(0.150031049774234d0*x1 + 0.0271505358982603d0*x6
     &)
      x329 = x147*x327*(0.225046574661351d0*x1 + 0.0325806430779124d0*x6
     &)/x32**3
      x330 = 0.311172355472072d0*x1
      x331 = 0.00285512867044833d0*x1 - 0.103724118490691d0*x165 + x20 +
     & 0.000428261234514334d0*x21 + 0.0121035293512547d0*x22*x6 + 0.0175
     &165224763768d0*x6
      x332 = x17*x18
      x333 = 0.0837628205355044d0*x332
      x334 = x2*x21
      x335 = 0.337338207872387d0*x334
      x336 = 0.0298137034090898d0*x121*x53 - 0.00370771454548634d0*x164 
     &+ 0.167525641071009d0*x172 + 0.899568554326366d0*x174 - 0.08376282
     &05355044d0*x25 - 0.0652591516234546d0*x27 + x307*x335 + 5.42506773
     &629845d0*x310*x334 - x312*x335
      x337 = 0.0891107483661069d0*x131*x332
      x338 = 0.0138644007334803d0*x186*x21
      x339 = 0.00571025734089666d0*x1 - 0.00609784657412857d0*x17 + 0.00
     &0856522469028669d0*x21 + 0.0350330449527535d0*x6
      x340 = x132*x339
      x341 = 0.0891107483661069d0*x17*x23
      x342 = 0.0231073345558005d0*x134
      x343 = 0.00107855455836145d0*x72
      x344 = 0.000586927869208408d0*x63
      x345 = 0.000440195901906306d0*x152*x53
      x346 = x137*x54
      x347 = x319 + x324 - x325 - x326 + x328 - x329
      x348 = x153*(-0.154912426780983d0*x0**(-0.333333333333333d0)*x55 -
     & x138*( -1.36477107997868d0*x1 + 0.427594146176807d0*x6) - 4.83324
     &359125373d0*x139*x346 - 1.61108119708458d0*x140 - 75.398223686155d
     &0*x346*x56 + x347)
      x349 = 0.000781467479590561d0*x190
      x350 = 0.0110340749409504d0*x66
      x351 = 0.0198949583613061d0*x68
      x352 = 0.0165511124114256d0*x171*x53
      x353 = 0.00471283037863406d0*x121*x189
      x354 = 0.213300270968718d0*x173*x67
      x355 = x192*x2
      x356 = x307*x355
      x357 = x2*x310
      x358 = 0.857572898037877d0*x357*x67
      x359 = x312*x355
      x360 = x170*x330 + x331
      x361 = x18*x191*x360
      x362 = 0.89531658564596d0*x3
      x363 = x176 + x178 + x180
      x364 = mu*x181
      x365 = x364*x44
      x366 = 2.05391854581537d0*x177 + 1.43555326817173d0*x179 + x362 - 
     &2.0d0*x363* x365 + x364*x46*(1.42618991424766d0*x177 + 1.534828432
     &53593d0*x3) /x43**2 - x47*(0.950793276165107d0*x177 + x362)
      x367 = 0.000821860773992031d0*x198
      x368 = 0.00880391803812612d0*x187
      x369 = 0.202642367284676d0*mu*x182*x363/x46**2
      x370 = 0.00316102962473761d0*x2*x200
      x371 = -6.39350882283002d-5*x0**2.33333333333333d0*x16*x23 + 0.000
     &279369479946145d0*x121*x194 - x18*x188*x339 - 5.79052675961456d-5*
     &x195 - 0.00192189252234297d0*x196 - 0.0117385573841682d0*x197 + 0.
     &00191767513931474d0*x199 + 0.0508356285485896d0*x200*x357 + 0.0168
     &588246652672d0*x201 + x307 *x370 - x312*x370
      x372 = x104*x108
      x373 = rho_c**(-4.33333333333333d0)*x93
      x374 = 2391.63081129926d0*x373
      x375 = x229*x373
      x376 = 1434.97848677955d0*x375
      x377 = x88**(-2)
      x378 = 307.49539002419d0*x377
      x379 = x228**2*x378
      x380 = x373*x379
      x381 = 13.5881085271426d0*x221*x229
      x382 = x373*x381
      x383 = x121*x90
      x384 = 4945.13601053991d0*x383
      x385 = x375*x384
      x386 = -x116 + x123 + x127 + x205
      x387 = x228*x378*x386
      x388 = x373*x387
      x389 = 1434.97848677955d0*x209
      x390 = rho_c**(-6.66666666666667d0)*x102
      x391 = x244*x390
      x392 = 614.99078004838d0*x229
      x393 = x314*x87
      x394 = x6*x88/x2**2
      x395 = x224*x307
      x396 = x310*x88/x9**2
      x397 = x224*x310
      x398 = x306*x92
      x399 = x224*x312
      x400 = x203*x315
      x401 = x116 - x123 - x127 - x205
      x402 = x221*x401
      x403 = x383*x401
      x404 = x204*(x12*(x337 - x338 - x340 - x341 + x342 - x343 - x344 +
     & x345 + x348 ) + x184*x366 - x185*x365 - x186*x367 + x186*x368 + x
     &369*x40 + x371 + x45*(-x162*x316 + x163*(-x317 - x318 + x320 - x32
     &1 - x322 - x323 + x347) - x333*x360 + x336 + 0.0009195062450792d0*
     &x35) + x71*(-x349 - x350 - x351 + x352 + x353 + x354 + x356 + x358
     & - x359 - x361))
      x405 = 307.49539002419d0*x217*x373*(0.0101850765748927d0*x115 - 2.
     &9296748261088d0*x126 - 0.16379624077176d0*x222 - 0.245694361157641
     &d0*x306 - 44.7079161356516d0*x308 - 718.99200186779d0*x311 + 44.70
     &79161356516d0*x313 + 1.38933365612628d0*x393 + 0.0108571437238599d
     &0*x394 + 44.7079161356516d0*x395 + 1437.98400373558d0*x396 + 718.9
     &9200186779d0*x397 - 3.95125283943823d0*x398 - 44.7079161356516d0*x
     &399 + 3.51560979133057d0*x400 - 0.245694361157641d0*x402 + 89.4158
     &322713032d0*x403 + 2.78d0*x404)
      x406 = 1434.97848677955d0*x215
      x407 = x210*x212
      x408 = x244*x407
      x409 = x258*x408
      x410 = x215*x392
      x411 = x104*x258
      x412 = x209*x237
      x413 = 614.99078004838d0*x412
      x414 = 378213.659544515d0*gamma_cc
      x415 = x213/x211**2
      x416 = x237*x389
      x417 = exp(4.0d0*x10 + x15*(64.327936885829629d0*x37 + 13.03556541
     &3083715d0*x49 - 64.327936885829629d0*x64 - 64.327936885829629d0*x6
     &9 - 0.16990136517137816d0*x73 - 0.012644118498950429d0*x77))
      x418 = x236**2
      x419 = x208*x418
      x420 = 614.99078004838d0*x417*x419
      x421 = x206*x232
      x422 = x418*x421
      x423 = 0.117839022137957d0*x115 - 33.8956721780452d0*x126 - 258.63
     &0216499205d0* x233*x310 - x234*x307 + x234*x312 - 2.8426279420419d
     &0*x306 + 16.0742747742176d0*x393 + 40.6748066136542d0*x400 + 32.16
     &39684429148d0*x404
      x424 = x421*x423
      x425 = 13747.478109301d0*x117
      x426 = gamma_cc*x390
      x427 = x208*x237*x426
      x428 = x206*x246*x390
      x429 = x256*x391
      x430 = 94553.4148861287d0*x391
      x431 = x377*x430
      x432 = 94553.4148861287d0*x103
      x433 = x244*x273 + x96
      x434 = x244*x433
      x435 = x272*x432 - x298*x434 + x99
      x436 = 614.99078004838d0*x95
      x437 = 378213.659544515d0*x103*x244 + x436
      x438 = 189106.829772257d0*x103
      E = -rho_c*(-x111 + x114)
      d1E(1) = -rho_c*(-x0*x116 + x0*x123 + x0*x127 + x0*x205 - x260*x26
     &2/(gamma_cc* x214 + 1.0d0)) + x111 + x113 - x83
      d1E(3) = rho_c*x275*x277
      d2E(1) = 0.0621813817393098d0*gamma_cc*x302*x304 + rho_c*(0.499760
     &307959092d0* x112*x305*x314 - 0.00366369660967365d0*x115*x305 - x2
     &04*x305*( -x12*(-x337 + x338 + x340 + x341 - x342 + x343 + x344 - 
     &x345 - x348) - x263*x369 - x281*x366 + x282*x365 + x283*x367 - x28
     &3*x368 + x371 + x45*(-x163*(x317 + x318 - x319 - x320 + x321 + x32
     &2 + x323 - x324 + x325 + x326 - x328 + x329) - 0.0009195062450792d
     &0* x279 + x280*x316 - x333*(-x330*(x166 + x167 - x168 + x169) + x3
     &31 ) + x336) - x71*(x349 + x350 + x351 - x352 - x353 - x354 - x356
     & - x358 + x359 + x361)) - x262*x303*(-gamma_cc*rho_c**( -6.6666666
     &6666667d0)*x102*x232*(-x217*(-10.9615243283801d0*x115 + 64154.8978
     &434045d0*x122 + 2522.41220923546d0*x126 + 3989.24019324716d0*x205 
     &+ 176.28306129213d0*x222 - 64154.8978434045d0*x225) + x217*(3.1318
     &6409382288d0*x115 - 900.861503298378d0*x126 - 50.3665889406086d0*x
     &222 - 75.5498834109129d0*x306 - x307*x425 - 221086.726038609d0*x31
     &1 + x312*x425 + 427.213694464283d0*x393 + 3.33852164391699d0*x394 
     &+ 13747.478109301d0*x395 + 442173.452077218d0*x396 + 221086.726038
     &609d0*x397 - 1214.99203294725d0*x398 - 13747.478109301d0*x399 + 10
     &81.03380395805d0*x400 - 75.5498834109129d0*x402 + 27494.9562186019
     &d0*x403 + 854.837184267248d0*x404) + x229*x384 + x379 - x381 - x38
     &7 - x392* x412 + x416 + x420 - x422 - x424 + 2391.63081129926d0) -
     & x208* x411*x415*(x215*x413 + 882498.538937201d0*x245 + x248*x414 
     &- 189106.829772257d0*x257 + x406 - x410) - x210*x374 + x210*x376 -
     & x210*x380 + x210*x382 - x210*x385 + x210*x388 - x210*x405 - 614.9
     &9078004838d0*x212*x240*x247*x411 - x240*x391*x392 + x241* x389*x39
     &0 - x259*(x208*x423*x428 - x217*x430*( 0.0203701531497855d0*x115 -
     & 5.85934965221761d0*x126 - 0.327592481543521d0*x222 - 0.4913887223
     &15281d0*x306 - 89.4158322713032d0*x308 - 1437.98400373558d0*x311 +
     & 89.4158322713032d0*x313 + 2.77866731225255d0*x393 + 0.02171428744
     &77198d0*x394 + 89.4158322713032d0*x395 + 2875.96800747116d0*x396 +
     & 1437.98400373558d0*x397 - 7.90250567887647d0*x398 - 89.4158322713
     &032d0*x399 + 7.03121958266113d0*x400 - 0.491388722315281d0*x402 + 
     &178.831664542606d0*x403 + 5.56d0*x404) + 4178.28073124474d0*x221* 
     &x429 - x255**2*x431 + x255*x386*x431 + 378213.659544515d0*x256* x4
     &27 - x373*x416 - x373*x420 + x373*x422 + x373*x424 - x374 + x375*x
     &413 + x376 - x380 + x382 - 1520606.52628364d0*x383*x429 - x385 + x
     &388 - 2500412.52698874d0*x391 - x405 + x419*x428 - 1764997.0778744
     &d0*x427 + 882498.538937201d0*x429 - 567320.489316772d0*x417*x418*x
     &426/x207**3) + x406*x409 - x409* x410) - 1.05383986550676d0*x278*x
     &305 - 1.26460783860812d0*x284* x305*x315 + 0.0883792666034679d0*x3
     &05*x306 + 258.630216499205d0* x305*x311 + x308*x309 - x309*x313 - 
     &x260*x261*x302*x372/x110**2) + 0.00549554491451047d0*x0*x115 - 32.
     &1639684429148d0*x0*x122 + 1.26460783860812d0*x0*x278 - 2.0d0*x0*x2
     &85
      d2E(4) = x277*(rho_c*x108*(x100*x297*x372*x432 - x104*x209*x297*x4
     &37*x82*x99/x107 **2 + x108*x288*x434 + x109*x215*x232*x289*x433*x9
     &8 + x109*x276* x435*(x287 + x288 + x291 - x300) - 614.99078004838d
     &0*x245*x286 + 307.49539002419d0*x247*x286*x372*x433 - 94553.414886
     &1287d0*x289* x292 - 220624.6347343d0*x292 + x299*(882498.538937201
     &d0*x292 + 189106.829772257d0*x293 + x295*x414 + x296) + x301) + x2
     &75)
      d2E(6) = -rho_c*x261*x303*(x209*(58149478.3770584d0*rho_c**(-7.0d0
     &)*x408*x92** 8.34d0 + x212*x434*x438 + x407*x433*x436 - x415*x434*
     &x437 - x438 ) + x304*x435*(x213 - x214*x434 + x244*x432))
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu, 
     &E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s, gamma_cc, mu
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, 
     &x66, x67
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = 0.0621813817393098d0*(0.132568899905202d0*x1 + 1.0d0)*log(1.0
     &d0 + 16.081984221457407d0/(2.2255694211506878d0*x1 + 5.98255043577
     &10795d0*x2 + 0.80042863499936356d0*x3 + 0.18970043257475591d0*x4))
     &
      x6 = rho_s**4/rho_c**4
      x7 = rho_s*x0
      x8 = x7 + 1.0d0
      x9 = -x7 + 1.0d0
      x10 = 1.92366105093154d0*x8**1.33333333333333d0 + 1.92366105093154
     &d0*x9** 1.33333333333333d0 - 3.84732210186307d0
      x11 = 0.0197517897025652d0*x10*(0.0690139921125583d0*x1 + 1.0d0)*(
     &-x6 + 1.0d0) *log(1.0d0 + 29.608574643216677d0/(2.2475918635776164
     &d0*x1 + 8.1574147034876408d0*x2 + 0.43009724712766434d0*x3 + 0.191
     &15125951273382d0*x4))
      x12 = 0.0310906908696548d0*(0.127469618870009d0*x1 + 1.0d0)*log(1.
     &0d0 + 32.163968442914921d0/(3.8447462374472119d0*x1 + 11.120374863
     &094685d0*x2 + 1.6447337755676092d0*x3 + 0.24058712912881922d0*x4))
     &
      x13 = -x5
      x14 = x10*x6
      x15 = x14*(x12 + x13)
      x16 = mu**2
      x17 = (0.237113969739022d0*x16*x4 + 1.0d0)**(-4)
      x18 = rho_s**2/rho_c**2
      x19 = -x18 + 1.0d0
      x20 = x0**1.0d0
      x21 = x0**1.33333333333333d0
      x22 = (0.0128480790170175d0*x1 - 0.00304892328706429d0*x20 + 0.000
     &275310793616358d0*x21 + 0.0315297404574782d0*x4 + 1.0d0)*exp (-0.4
     &6675853320810856d0*x1)
      x23 = x19*x22
      x24 = x0**1.66666666666667d0
      x25 = 0.00264117541143784d0*mu**5*x24
      x26 = 0.693147180559945d0
      x27 = x8**0.666666666666667d0 + x9**0.666666666666667d0
      x28 = x27**3
      x29 = x1*x16/x27**2
      x30 = 9.20897059521559d0*mu*x2/x27 + 1.0d0
      x31 = mu**3
      x32 = 0.0253302959105844d0*x28*log((18.485266912338322d0*x29 + x30
     & + 15.312568193831799d0*x3*x31/x28)/(8.5571394854859619d0*x29 + x3
     &0 ))
      x33 = x11 + x13
      x34 = -x15 + x33
      x35 = 0.00316102962473761d0*mu**8*x0**2.66666666666667d0
      x36 = 0.0445553741830535d0*x20
      x37 = 1d0/x4
      x38 = 1d0/x8
      x39 = x38**0.666666666666667d0
      x40 = x1*x38**0.333333333333333d0
      x41 = 0.0244354823084343d0*x4
      x42 = 8.61370200056458d0*x8**2*(-0.0177069812999059d0*x40 + 1.0d0)
     &/(x39*(x39* x41 + 0.337569861992026d0*x40 + 1.0d0))
      x43 = 1d0/x9
      x44 = x43**0.666666666666667d0
      x45 = x1*x43**0.333333333333333d0
      x46 = 8.61370200056458d0*x9**2*(-0.0177069812999059d0*x45 + 1.0d0)
     &/(x44*(x41* x44 + 0.337569861992026d0*x45 + 1.0d0))
      x47 = 0.000132058770571892d0*x24*(37.6991118430775d0*x19*(-3.07073
     &492995203d0* x1 + 0.384834731559127d0*x4)*exp(-0.19230865217881402
     &d0*x1)/x20 + x37*x42 + x37*x46)
      x48 = mu**6
      x49 = x19*(x22 - 1.0d0)
      x50 = 0.00248266686171384d0*x24*x49
      x51 = x14*(-x12 + x5) + x33
      x52 = 0.0533250677421794d0*x0**2.0d0*x51
      x53 = mu**4
      x54 = -0.0418814102677522d0*x20*x49
      x55 = 8.61370200056458d0*x8**2.66666666666667d0
      x56 = 8.61370200056458d0*x9**2.66666666666667d0
      x57 = 23.3866625383245d0*x19*(-0.240695990468967d0*x1 + 0.26014827
     &853397d0*x4) *exp(-0.33933171852197186d0*x1)
      x58 = 0.00020688890514282d0*x0**1.0d0
      x59 = 0.337338207872387d0*x21
      x60 = x17*(-x23*x25 - x31*(x23*x36 + x47) - x32*(-x26 + 1.0d0) + x
     &34*x35 - x48 *(x50 - x52) + x53*(x34*x59 + x54 - x58*(x42 + x46 - 
     &x55 - x56 + x57)))
      x61 = x8**0.666666666666667d0 + x9**0.666666666666667d0
      x62 = x61**3
      x63 = x22*(x18 - 1.0d0)
      x64 = -1.0d0 + exp(-257.31174754331852d0*(-x17*(x25*x63 + x31*(x36
     &*x63 - x47) + x32*(x26 - 1.0d0) + x35*x51 + x48*(-x50 + x52) + x53
     &*(x51*x59 + x54 + x58*(-x42 - x46 + x55 + x56 - x57))) + x51)/x62)
     &
      x65 = (x34 - x60)/x34
      x66 = 0.544846825002868d0*gamma_cc*rho_c**(-2.33333333333333d0)*x6
     &5**2.78d0/ x61**2
      x67 = 1.0d0 + x66/x64
      E = -rho_c*(-x11 + x15 + x5 + x60 - 0.00388633635870686d0*x62*log(
     &x66*x67/( 0.29685806271570547d0*gamma_cc**2*rho_c**(-4.66666666666
     &6667d0)* x65**5.5599999999999996d0/(x61**4*x64**2) + x67) + 1.0d0)
     &)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu
     &, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x424, x425, x426, x427, x428, x429, x430, x431, x432, x433
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = 3.84474623744721d0*x1 + 11.1203748630947d0*x2 + 1.64473377556
     &761d0*x3 + 0.240587129128819d0*x4
      x6 = 1.0d0 + 32.1639684429149d0/x5
      x7 = log(x6)
      x8 = 0.0310906908696548d0*x7*(0.127469618870009d0*x1 + 1.0d0)
      x9 = 2.22556942115069d0*x1 + 5.98255043577108d0*x2 + 0.80042863499
     &9364d0*x3 + 0.189700432574756d0*x4
      x10 = 1.0d0 + 16.0819842214574d0/x9
      x11 = log(x10)
      x12 = 0.0621813817393098d0*x11*(0.132568899905202d0*x1 + 1.0d0)
      x13 = -x12
      x14 = x13 + x8
      x15 = rho_s*x0
      x16 = x15 + 1.0d0
      x17 = -x15 + 1.0d0
      x18 = 1.92366105093154d0*x16**1.33333333333333d0 + 1.9236610509315
     &4d0*x17** 1.33333333333333d0 - 3.84732210186307d0
      x19 = rho_c**(-4)
      x20 = rho_s**4
      x21 = x19*x20
      x22 = x18*x21
      x23 = x14*x22
      x24 = -x21 + 1.0d0
      x25 = 2.24759186357762d0*x1 + 8.15741470348764d0*x2 + 0.4300972471
     &27664d0*x3 + 0.191151259512734d0*x4
      x26 = 1.0d0 + 29.6085746432167d0/x25
      x27 = log(x26)
      x28 = x24*x27
      x29 = x18*x28
      x30 = 0.0197517897025652d0*x29*(0.0690139921125583d0*x1 + 1.0d0)
      x31 = x13 + x30
      x32 = -x23 + x31
      x33 = mu**2
      x34 = x33*x4
      x35 = 0.237113969739022d0*x34 + 1.0d0
      x36 = x35**(-4)
      x37 = mu**5
      x38 = rho_c**(-2)
      x39 = rho_s**2*x38
      x40 = -x39 + 1.0d0
      x41 = exp(-0.46675853320810856d0*x1)
      x42 = x0**1.0d0
      x43 = 0.00304892328706429d0*x42
      x44 = -x43
      x45 = x0**1.33333333333333d0
      x46 = 0.0128480790170175d0*x1 + 0.0315297404574782d0*x4 + x44 + 0.
     &000275310793616358d0*x45 + 1.0d0
      x47 = x41*x46
      x48 = x40*x47
      x49 = x37*x48
      x50 = x0**1.66666666666667d0
      x51 = 0.00264117541143784d0*x50
      x52 = x49*x51
      x53 = 0.693147180559945d0
      x54 = -x53 + 1.0d0
      x55 = x16**0.666666666666667d0 + x17**0.666666666666667d0
      x56 = x55**2
      x57 = 1d0/x56
      x58 = x33*x57
      x59 = x1*x58
      x60 = 1d0/x55
      x61 = mu*x60
      x62 = 9.20897059521559d0*x2*x61 + 1.0d0
      x63 = mu**3
      x64 = x55**3
      x65 = x3/x64
      x66 = 18.4852669123383d0*x59 + x62 + 15.3125681938318d0*x63*x65
      x67 = x66/(8.55713948548596d0*x59 + x62)
      x68 = log(x67)
      x69 = x54*x68
      x70 = 0.0253302959105844d0*x64
      x71 = x69*x70
      x72 = mu**8*x0**2.66666666666667d0
      x73 = 0.00316102962473761d0*x72
      x74 = x32*x73
      x75 = x40*x42
      x76 = 0.0445553741830535d0*x47
      x77 = x75*x76
      x78 = -3.07073492995203d0*x1 + 0.384834731559127d0*x4
      x79 = exp(-0.19230865217881402d0*x1)
      x80 = x40*x78*x79
      x81 = 1d0/x42
      x82 = 37.6991118430775d0*x81
      x83 = x80*x82
      x84 = 1d0/x4
      x85 = x16**2
      x86 = 1d0/x16
      x87 = x86**0.333333333333333d0
      x88 = x1*x87
      x89 = x86**0.666666666666667d0
      x90 = 0.0244354823084343d0*x4
      x91 = 0.337569861992026d0*x88 + x89*x90 + 1.0d0
      x92 = 1d0/x91
      x93 = x85*x92
      x94 = (-0.0177069812999059d0*x88 + 1.0d0)/x89
      x95 = 8.61370200056458d0*x94
      x96 = x93*x95
      x97 = x17**2
      x98 = 1d0/x17
      x99 = x98**0.333333333333333d0
      x100 = x1*x99
      x101 = x98**0.666666666666667d0
      x102 = 0.337569861992026d0*x100 + x101*x90 + 1.0d0
      x103 = 1d0/x102
      x104 = x103*x97
      x105 = (-0.0177069812999059d0*x100 + 1.0d0)/x101
      x106 = 8.61370200056458d0*x105
      x107 = x104*x106
      x108 = x50*(x107*x84 + x83 + x84*x96)
      x109 = 0.000132058770571892d0*x108
      x110 = x63*(x109 + x77)
      x111 = mu**6
      x112 = x47 - 1.0d0
      x113 = x112*x40
      x114 = x113*x50
      x115 = 0.00248266686171384d0*x114
      x116 = x0**2.0d0
      x117 = x12 - x8
      x118 = x117*x22
      x119 = x118 + x31
      x120 = x116*x119
      x121 = 0.0533250677421794d0*x120
      x122 = x111*(x115 - x121)
      x123 = mu**4
      x124 = 0.0418814102677522d0*x113*x42
      x125 = -x124
      x126 = x16**2.66666666666667d0
      x127 = 8.61370200056458d0*x126
      x128 = x17**2.66666666666667d0
      x129 = 8.61370200056458d0*x128
      x130 = exp(-0.33933171852197186d0*x1)
      x131 = x130*(-0.240695990468967d0*x1 + 0.26014827853397d0*x4)
      x132 = x131*x40
      x133 = 23.3866625383245d0*x132
      x134 = x0**1.0d0
      x135 = 0.00020688890514282d0*x134
      x136 = 0.337338207872387d0*x45
      x137 = x123*(x125 - x135*(x107 - x127 - x129 + x133 + x96) + x136*
     &x32)
      x138 = -x110 - x122 + x137 - x52 - x71 + x74
      x139 = x138*x36
      x140 = -x139 + x32
      x141 = x140/x32
      x142 = x141**2.78d0
      x143 = 0.544846825002868d0*rho_c**(-2.33333333333333d0)
      x144 = x142*x143
      x145 = x16**0.666666666666667d0 + x17**0.666666666666667d0
      x146 = x145**2
      x147 = gamma_cc/x146
      x148 = x39 - 1.0d0
      x149 = x148*x47
      x150 = x37*x51
      x151 = x53 - 1.0d0
      x152 = x151*x68
      x153 = x148*x42
      x154 = x153*x76
      x155 = -x107 + x127 + x129 - x133 - x96
      x156 = x111*(-x115 + x121) + x119*x73 + x123*(x119*x136 + x125 + x
     &135*x155) + x149*x150 + x152*x70 + x63*(-x109 + x154)
      x157 = x156*x36
      x158 = x119 - x157
      x159 = x145**3
      x160 = 1d0/x159
      x161 = 257.311747543319d0*x160
      x162 = exp(-x158*x161)
      x163 = x162 - 1.0d0
      x164 = 1d0/x163
      x165 = x147*x164
      x166 = x144*x165 + 1.0d0
      x167 = x163**(-2)
      x168 = rho_c**(-4.66666666666667d0)
      x169 = x141**5.56d0
      x170 = x168*x169
      x171 = 0.296858062715705d0*gamma_cc**2/x145**4
      x172 = 1d0/(x166 + x167*x170*x171)
      x173 = x166*x172
      x174 = x144*x147*x173
      x175 = x174 + 1.0d0
      x176 = log(x175)
      x177 = 0.00388633635870686d0*x159
      x178 = x12 - x30
      x179 = x178 + x23
      x180 = x139 + x179
      x181 = x1*x11
      x182 = 0.00274777245725523d0*x181
      x183 = rho_c**(-5)
      x184 = x183*x20
      x185 = x18*x184
      x186 = 0.00136314985874174d0*x1 + 0.0197517897025652d0
      x187 = x186*x27
      x188 = 4.0d0*x187
      x189 = x1*x29
      x190 = 0.000454383286247248d0*x189
      x191 = x17**0.333333333333333d0
      x192 = x16**0.333333333333333d0
      x193 = x191 - x192
      x194 = 2.56488140124205d0*x193
      x195 = x186*x28
      x196 = x194*x195
      x197 = rho_s*x38
      x198 = (0.0082433173717657d0*x1 + 0.0621813817393098d0)*(0.7418564
     &73716896d0*x1 + 0.99709173929518d0*x2 + 0.400214317499682d0*x3 + 0
     &.126466955049837d0*x4)/(x10*x9**2)
      x199 = 16.0819842214574d0*x198
      x200 = rho_s**5
      x201 = x194*x200
      x202 = 4.0d0*x14
      x203 = x18*x186*x24*(0.749197287859205d0*x1 + 1.35956911724794d0*x
     &2 + 0.215048623563832d0*x3 + 0.127434173008489d0*x4)/(x25**2*x26)
      x204 = 29.6085746432167d0*x203
      x205 = 0.00132103950518672d0*x1*x7
      x206 = 32.1639684429149d0*(0.00396311851556016d0*x1 + 0.0310906908
     &696548d0)*( 1.28158207914907d0*x1 + 1.85339581051578d0*x2 + 0.8223
     &66887783805d0*x3 + 0.160391419419213d0*x4)/(x5**2*x6)
      x207 = x182 - x199
      x208 = x110 + x122 - x137 + x52 + x71 - x74
      x209 = x34/x35**5
      x210 = 0.632303919304059d0*x209
      x211 = x39*x47
      x212 = 0.00528235082287567d0*x37*x50
      x213 = x211*x212
      x214 = 0.00428269300567249d0*x1
      x215 = 0.0210198269716521d0*x4
      x216 = 0.000367081058155144d0*x45
      x217 = x214 + x215 + x216 + x44
      x218 = x150*x41
      x219 = x217*x218*x40
      x220 = 0.00440195901906306d0*x50
      x221 = x220*x49
      x222 = 0.000410930386996015d0*x0**2.0d0
      x223 = x222*x49
      x224 = x16**(-0.333333333333333d0)
      x225 = x17**(-0.333333333333333d0)
      x226 = x224 - x225
      x227 = 0.0506605918211689d0*x56
      x228 = x226*x227
      x229 = x228*x69
      x230 = x15*x229
      x231 = 0.00842941233263362d0*x72
      x232 = x231*x32
      x233 = x1*x61
      x234 = 2.85237982849532d0*x233
      x235 = x15*x226
      x236 = mu*x1*x235*x57
      x237 = 11.4095193139813d0*x236
      x238 = 1.53482843253593d0*x2
      x239 = 6.13931373014372d0*x2
      x240 = x235*x239*x60
      x241 = -x238 + x240
      x242 = x3*x58
      x243 = -6.16175563744611d0*x233 + 30.6251363876636d0*x235*x33*x65 
     &+ 24.6470225497844d0*x236 + x241 - 7.6562840969159d0*x242
      x244 = 0.0253302959105844d0*mu/x66
      x245 = x244*x56
      x246 = x245*x54*(x243 - x67*(-x234 + x237 + x241))
      x247 = x183*x201
      x248 = -x182 + x199
      x249 = x22*(x205 - x206 + x248)
      x250 = x188*x22
      x251 = x15*x196
      x252 = -x190 + x204 + x207 + x249 + x250 + x251
      x253 = -x14*x247 + 4.0d0*x23 + x252
      x254 = x253*x73
      x255 = 75.398223686155d0*x78*x79*x81
      x256 = x255*x39
      x257 = 2.41662179562687d0*x0**(-0.666666666666667d0)*x80
      x258 = x40*x79*x82*(1.02357830998401d0*x1 - 0.256556487706084d0*x4
     &)
      x259 = x15*x86
      x260 = 1d0/x1
      x261 = 0.0508408867489862d0*x260
      x262 = x261/x87
      x263 = x262*x93*(-x259 + 1.0d0)
      x264 = 5.74246800037638d0*x84
      x265 = x264*x93*x94
      x266 = x15*x98
      x267 = x261/x99
      x268 = x104*x267*(x266 + 1.0d0)
      x269 = 22.9698720015055d0*x84
      x270 = x15*x269
      x271 = x16*x92
      x272 = x271*x94
      x273 = x270*x272
      x274 = x104*x105*x264
      x275 = x103*x17
      x276 = x105*x275
      x277 = x270*x276
      x278 = 0.112523287330675d0*x88
      x279 = 0.0162903215389562d0*x4
      x280 = x279*x89
      x281 = x278 + x280
      x282 = x84*x95/x91**2
      x283 = x282*x85*(-x259*x278 - x259*x280 + x281)
      x284 = 0.112523287330675d0*x100
      x285 = x101*x279
      x286 = x284 + x285
      x287 = x106*x84/x102**2
      x288 = x287*x97*(x266*x284 + x266*x285 + x286)
      x289 = 0.000132058770571892d0*x50
      x290 = 0.0445553741830535d0*x41
      x291 = 0.00693220036674015d0*x45
      x292 = 0.0891107483661069d0*x42
      x293 = 0.000220097950953153d0*x108 - x211*x292
      x294 = x63*(x217*x290*x75 - x289*(x256 + x257 + x258 + x263 + x265
     & + x268 - x273 + x274 + x277 + x283 + x288 + x83) - x291*x48 + x29
     &3 + x77)
      x295 = x112*x39
      x296 = 0.00496533372342768d0*x50
      x297 = x295*x296
      x298 = 0.0041377781028564d0*x114
      x299 = x41*(0.0777930888680181d0*x1*x46 - 0.00214134650283625d0*x1
     & - 0.0105099134858261d0*x4 + 0.00152446164353214d0*x42 - 0.0001835
     &40529077572d0*x45)
      x300 = x296*x299*x40
      x301 = 0.106650135484359d0*x120
      x302 = 4.0d0*x118
      x303 = x117*x247
      x304 = x252 - x302 + x303
      x305 = 0.0533250677421794d0*x116
      x306 = x304*x305
      x307 = x111*(x297 - x298 + x300 + x301 - x306)
      x308 = 0.0837628205355044d0*x42
      x309 = x295*x308
      x310 = 0.0837628205355044d0*x299*x75
      x311 = 0.0003448148419047d0*x134*x155
      x312 = 0.449784277163183d0*x119*x45
      x313 = x136*x304
      x314 = x17**1.66666666666667d0
      x315 = x16**1.66666666666667d0
      x316 = 46.773325076649d0*x131*x84
      x317 = -x263 - x265 - x268 + x273 - x274 - x277 - x283 - x288
      x318 = 0.00020688890514282d0*x50
      x319 = x318*(-23.3866625383245d0*x130*x40*x84*(0.0802319968229891d
     &0*x1 - 0.173432185689313d0*x4) - 2.64527879654103d0*x132*x260 - 15
     &.591108358883d0*x132*x84 + x264*(x126 + x128) + x270*(x314 - x315)
     & - x316*x39 + x317)
      x320 = x123*(x125 + x309 + x310 + x311 + x312 - x313 - x319)
      x321 = -x162 + 1.0d0
      x322 = x180/x179
      x323 = x168*x322**5.56d0
      x324 = x143*x322**2.78d0
      x325 = x147/x321
      x326 = -x324*x325 + 1.0d0
      x327 = 1d0/(x171*x323/x321**2 + x326)
      x328 = x324*x326
      x329 = x327*x328
      x330 = x147*x329 + 1.0d0
      x331 = log(x330)
      x332 = x17**(-0.333333333333333d0)
      x333 = x16**(-0.333333333333333d0)
      x334 = x332 - 1.0d0*x333
      x335 = 0.00777267271741373d0*x146
      x336 = 1d0/x330
      x337 = rho_c**(-3.33333333333333d0)
      x338 = x166*x337
      x339 = -x333
      x340 = x332 + x339
      x341 = 1d0/x145
      x342 = 0.726462433337157d0*x341
      x343 = rho_c**(-4.33333333333333d0)*rho_s*x342
      x344 = x340*x343
      x345 = 1d0/x140
      x346 = -x213
      x347 = x345*(-x138*x210 - x141*x304 + x253 - x36*(x219 + x221 - x2
     &23 + x230 - x232 - x246 + x254 + x294 - x307 - x320 + x346))
      x348 = 1.51467417350797d0*x338
      x349 = x158/x119
      x350 = x349**2.78d0
      x351 = 1.27130925834002d0*x337
      x352 = -x214 - x215 - x216 + x43
      x353 = x149*x37
      x354 = -x224 + x225
      x355 = x36*(x111*(-x297 + x298 - x300 - x301 + x306) - x119*x231 +
     & x123*(x124 - x309 - x310 - x311 - x312 + x313 + x319) + x148*x218
     &*x352 + x15* x152*x227*x354 + x151*x245*(x243 + x67*(x234 - x237 +
     & x238 - x240 )) - x220*x353 + x222*x353 + x304*x73 + x346 + x63*(x
     &149*x291 + x153*x290*x352 - x154 + x289*(-x256 - x257 - x258 + x31
     &7 - x83) + x293))
      x356 = 1d0/x158
      x357 = 1.51467417350797d0*x337
      x358 = x356*x357
      x359 = 1.0d0*x332
      x360 = -x118 + x157 + x178
      x361 = 514.623495086637d0*x360
      x362 = x15*x341
      x363 = 659.974115594947d0*x117
      x364 = 1029.24699017327d0*x187
      x365 = 659.974115594947d0*x195
      x366 = 1029.24699017327d0*x118 - x15*x193*x365 + 162.699226454617d
     &0*x156*x209 - 0.707034132827743d0*x181 - x183*x193*x200*x363 + 0.1
     &16918157438755d0*x189 + 4138.08346398728d0*x198 - 7618.63408371287
     &d0*x203 - x22*x364 - 257.311747543319d0*x249 + 257.311747543319d0*
     &x355
      x367 = 0.544846825002867d0*x337
      x368 = exp(x161*x360)
      x369 = x368 - 1.0d0
      x370 = 1d0/x369
      x371 = x160*x370
      x372 = x368*x371
      x373 = x147*x370
      x374 = rho_c**(-5.66666666666667d0)
      x375 = x169*x374
      x376 = x165*x375
      x377 = x160*x164
      x378 = 0.791621500575215d0*gamma_cc*x340
      x379 = x339 + x359
      x380 = 514.623495086637d0*x158*x362*x379 + x366
      x381 = x162*x377
      x382 = 0.593716125431411d0*gamma_cc/x145**5
      x383 = 0.00388633635870686d0*x145
      x384 = gamma_cc*x172*x383
      x385 = rho_s**3*x18
      x386 = x14*x194
      x387 = x202*x385
      x388 = x15*x47
      x389 = x212*x388
      x390 = rho_c**(-3)
      x391 = x385*x390
      x392 = x188*x391
      x393 = -x392
      x394 = -x196 + x21*x386 - x387*x390 + x393
      x395 = 24.6470225497844d0*x233 + x239 + 30.6251363876636d0*x242
      x396 = x67*(11.4095193139813d0*x233 + x239)
      x397 = x244*x55
      x398 = x112*x15
      x399 = x296*x398
      x400 = x292*x388
      x401 = x15*x255
      x402 = x262*x271
      x403 = x267*x275
      x404 = x269*x272
      x405 = x269*x276
      x406 = x16*x281*x282
      x407 = x17*x286*x287
      x408 = x402 - x403 + x404 - x405 + x406 - x407
      x409 = x308*x398
      x410 = x269*(-x314 + x315)
      x411 = x15*x316
      x412 = x36*(x111*(x305*x394 + x399) + x123*(x136*x394 - x318*(x408
     & - x410 - x411) + x409) - x229 + x389 + x394*x73 + x397*x54*(x226*
     &x395 + x354*x396) + x63*(-x289*(-x401 + x408) + x400))
      x413 = x340*x342
      x414 = -x191 + x192
      x415 = 2.56488140124205d0*x414
      x416 = x195*x415
      x417 = x117*x21*x415
      x418 = x117*x391
      x419 = 4.0d0*x418
      x420 = x393 + x416 + x417 + x419
      x421 = x141*x420
      x422 = x187*x391
      x423 = 3.88495961657227d0*x193
      x424 = x14*x21
      x425 = x14*x391
      x426 = 359.584801465982d0*x193
      x427 = x334*x341
      x428 = x337*x350
      x429 = x349**5.56d0*x374
      x430 = -x402 + x403 - x404 + x405 - x406 + x407
      x431 = x36*(x111*(x305*x420 + x399) + x123*(x136*x420 + x318*(x410
     & + x411 + x430) + x409) + x151*x354*x397*(x395 - x396) + x152*x228
     & + x389 + x420*x73 + x63*(x289*(x401 + x430) + x400))
      x432 = x349*(x392 - x416 - x417 - x419) + x420 - x431
      x433 = -x21*x363*x414 + x361*x427 + x364*x391 - x365*x414 - 1029.2
     &4699017327d0* x418 + 257.311747543319d0*x431
      E = -rho_c*(-x176*x177 + x180)
      d1E(1) = rho_c*(x0*x182 - x0*x190 - x0*x199 + x0*x204 + x0*x208*x2
     &10 + x0*x36*( x213 - x219 - x221 + x223 - x230 + x232 + x246 - x25
     &4 - x294 + x307 + x320) - x142*x336*x384*(-x143*x165*x173*(rho_c**
     &( -6.66666666666667d0)*rho_s*x169*x377*x378 + x142*x344 - x142*x34
     &7 *x357 + x142*x351 + x142*x367*x380*x381 + x162*x167*x375*x380* x
     &382 - 1.65053082869932d0*x347*x376 + 1.38533762600663d0*x376) + x1
     &43*x350*x373*(-x343*(-x332 + x333) + x351 - x358*(-x156*x210 + x30
     &4 + x349*(x190 - x204 + x248 - x249 - x250 - x251 + x302 - x303) -
     & x355) + x367*x372*(x361*x362*(x333 - x359) + x366)) + x166*x344 +
     & 1.27130925834002d0*x338 - x347*x348) + x185*x188 + x185*x202 - x1
     &85*(-x205 + x206 + x207) + x196*x197 + x197*x331* x334*x335 - x14*
     &x201/rho_c**6) + x177*x331 + x208*x36 + x32
      d1E(2) = -rho_c*(x0*x176*x335*x379 + x0*x196 + x0*x412 - x184*x386
     & + x188*x19* x385 + x19*x387 - x384*(x142*x338*x413 + x142*x345*x3
     &48*(x394 - x412 - x421) - x174*x370*(x350*x358*x432 + 1.6505308286
     &9932d0* x356*x373*x429*x432 - x368*x382*x429*x433/x369**2 + x371*x
     &378* x429 - 0.544846825002867d0*x372*x428*x433 + x413*x428) + 0.54
     &4846825002868d0*x376*(x341*(0.726462433337157d0*x332 - 0.726462433
     &337157d0*x333) + x345*(-x195*x423 - 1.51467417350797d0 *x412 - 1.5
     &1467417350797d0*x421 - 6.05869669403189d0*x422 + x423* x424 - 6.05
     &869669403189d0*x425) - x381*(-280.390977369833d0*x158* x427 + x195
     &*x426 + 140.195488684916d0*x412 + 560.781954739666d0* x422 - x424*
     &x426 + 560.781954739666d0*x425)))/x175)
      d1E(3) = rho_c*x327*x336*x383*(-x165*x329*(x144 + 0.59371612543141
     &1d0*x165*x170) - 0.296858062715705d0*x323*x325 + x328)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_PBE_GWS_ERF(rho_c, rho_s, gamma_cc, mu
     &, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &05, x1506, x1507, x1508, x1509, x1510, x1511, x1512
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
      x11 = 0.296858062715705d0*x10
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
      x39 = x36**2
      x40 = 1d0/x39
      x41 = x13*x40
      x42 = x21*x41
      x43 = x0**0.166666666666667d0
      x44 = 1d0/x36
      x45 = mu*x44
      x46 = 9.20897059521559d0*x43*x45 + 1.0d0
      x47 = 8.55713948548596d0*x42 + x46
      x48 = 1d0/x47
      x49 = mu**3
      x50 = sqrt(x0)
      x51 = 1d0/x37
      x52 = x50*x51
      x53 = 18.4852669123383d0*x42 + x46 + 15.3125681938318d0*x49*x52
      x54 = x48*x53
      x55 = log(x54)
      x56 = 0.693147180559945d0
      x57 = x56 - 1.0d0
      x58 = x55*x57
      x59 = 1.64473377556761d0*x50
      x60 = 0.240587129128819d0*x12 + 3.84474623744721d0*x21 + 11.120374
     &8630947d0* x43 + x59
      x61 = 1.0d0 + 32.1639684429149d0/x60
      x62 = log(x61)
      x63 = 0.0310906908696548d0*x62*(0.127469618870009d0*x21 + 1.0d0)
      x64 = 0.800428634999364d0*x50
      x65 = 0.189700432574756d0*x12 + 2.22556942115069d0*x21 + 5.9825504
     &3577108d0* x43 + x64
      x66 = 1.0d0 + 16.0819842214574d0/x65
      x67 = log(x66)
      x68 = 0.0621813817393098d0*x67*(0.132568899905202d0*x21 + 1.0d0)
      x69 = -x63 + x68
      x70 = 1.92366105093154d0*x2**1.33333333333333d0 + 1.92366105093154
     &d0*x3** 1.33333333333333d0 - 3.84732210186307d0
      x71 = rho_c**(-4)
      x72 = rho_s**4
      x73 = x71*x72
      x74 = x70*x73
      x75 = x69*x74
      x76 = -x68
      x77 = -x73 + 1.0d0
      x78 = 0.430097247127664d0*x50
      x79 = 0.191151259512734d0*x12 + 2.24759186357762d0*x21 + 8.1574147
     &0348764d0* x43 + x78
      x80 = 1.0d0 + 29.6085746432167d0/x79
      x81 = log(x80)
      x82 = x77*x81
      x83 = x70*x82
      x84 = 0.0197517897025652d0*x83*(0.0690139921125583d0*x21 + 1.0d0)
      x85 = x76 + x84
      x86 = x75 + x85
      x87 = mu**8*x0**2.66666666666667d0
      x88 = 0.00316102962473761d0*x87
      x89 = x23*x28
      x90 = 0.0445553741830535d0*x89
      x91 = x20*x90
      x92 = 0.384834731559127d0*x12 - 3.07073492995203d0*x21
      x93 = exp(-0.19230865217881402d0*x21)
      x94 = x93/x23
      x95 = x92*x94
      x96 = -x19 + 1.0d0
      x97 = 37.6991118430775d0*x96
      x98 = x95*x97
      x99 = 1d0/x12
      x100 = x2**2
      x101 = 8.61370200056458d0*x100
      x102 = 1d0/x2
      x103 = x102**0.333333333333333d0
      x104 = x103*x21
      x105 = 0.0177069812999059d0*x104
      x106 = -x105 + 1.0d0
      x107 = x102**0.666666666666667d0
      x108 = 0.0244354823084343d0*x12
      x109 = 0.337569861992026d0*x104 + x107*x108 + 1.0d0
      x110 = 1d0/x109
      x111 = 1d0/x107
      x112 = x110*x111
      x113 = x106*x112
      x114 = x101*x113
      x115 = x3**2
      x116 = 8.61370200056458d0*x115
      x117 = 1d0/x3
      x118 = x117**0.333333333333333d0
      x119 = x118*x21
      x120 = 0.0177069812999059d0*x119
      x121 = -x120 + 1.0d0
      x122 = x117**0.666666666666667d0
      x123 = x108*x122 + 0.337569861992026d0*x119 + 1.0d0
      x124 = 1d0/x123
      x125 = 1d0/x122
      x126 = x124*x125
      x127 = x121*x126
      x128 = x116*x127
      x129 = x30*(x114*x99 + x128*x99 + x98)
      x130 = 0.000132058770571892d0*x129
      x131 = mu**6
      x132 = x28 - 1.0d0
      x133 = x132*x30
      x134 = x133*x96
      x135 = 0.00248266686171384d0*x134
      x136 = x0**2.0d0
      x137 = x136*x86
      x138 = 0.0533250677421794d0*x137
      x139 = mu**4
      x140 = x132*x23
      x141 = 0.0418814102677522d0*x140*x96
      x142 = -x141
      x143 = x2**2.66666666666667d0
      x144 = 8.61370200056458d0*x143
      x145 = x3**2.66666666666667d0
      x146 = 8.61370200056458d0*x145
      x147 = 0.240695990468967d0*x21
      x148 = 0.26014827853397d0*x12
      x149 = exp(-0.33933171852197186d0*x21)
      x150 = x149*(-x147 + x148)
      x151 = x150*x96
      x152 = 23.3866625383245d0*x151
      x153 = -x114 - x128 + x144 + x146 - x152
      x154 = x0**1.0d0
      x155 = 0.00020688890514282d0*x154
      x156 = 0.337338207872387d0*x26
      x157 = x131*(-x135 + x138) + x139*(x142 + x153*x155 + x156*x86) + 
     &x29*x33 + x38 *x58 + x49*(-x130 + x91) + x86*x88
      x158 = x157*x16
      x159 = -x158 + x86
      x160 = x4**3
      x161 = 1d0/x160
      x162 = 257.311747543319d0*x161
      x163 = exp(-x159*x162)
      x164 = x163 - 1.0d0
      x165 = x164**(-2)
      x166 = rho_c**(-4.66666666666667d0)
      x167 = x63 + x76
      x168 = x167*x74
      x169 = -x168 + x85
      x170 = 1d0/x169
      x171 = x28*x96
      x172 = x171*x33
      x173 = -x56 + 1.0d0
      x174 = x173*x55
      x175 = x174*x38
      x176 = x169*x88
      x177 = x90*x96
      x178 = x49*(x130 + x177)
      x179 = x131*(x135 - x138)
      x180 = x114 + x128 - x144 - x146 + x152
      x181 = x139*(x142 - x155*x180 + x156*x169)
      x182 = -x172 - x175 + x176 - x178 - x179 + x181
      x183 = x16*x182
      x184 = x169 - x183
      x185 = x170*x184
      x186 = x185**5.56d0
      x187 = x166*x186
      x188 = x165*x187
      x189 = x185**2.78d0
      x190 = rho_c**(-2.33333333333333d0)
      x191 = 0.544846825002868d0*x190
      x192 = x189*x191
      x193 = 1d0/x164
      x194 = x193*x7
      x195 = x192*x194 + 1.0d0
      x196 = x11*x188 + x195
      x197 = 1d0/x196
      x198 = x192*x195
      x199 = x197*x198
      x200 = x199*x7 + 1.0d0
      x201 = log(x200)
      x202 = 0.00388633635870686d0*x160
      x203 = x68 - x84
      x204 = x168 + x203
      x205 = x183 + x204
      x206 = x21*x67
      x207 = 0.00274777245725523d0*x206
      x208 = rho_c**(-5)
      x209 = x70*x72
      x210 = x208*x209
      x211 = 0.00136314985874174d0*x21 + 0.0197517897025652d0
      x212 = x211*x81
      x213 = 4.0d0*x212
      x214 = x21*x83
      x215 = 0.000454383286247248d0*x214
      x216 = x212*x77
      x217 = x3**0.333333333333333d0
      x218 = x2**0.333333333333333d0
      x219 = x217 - x218
      x220 = 2.56488140124205d0*x219
      x221 = x216*x220
      x222 = rho_s*x17
      x223 = 0.0082433173717657d0*x21 + 0.0621813817393098d0
      x224 = 1d0/x66
      x225 = x65**(-2)
      x226 = 0.126466955049837d0*x12 + 0.741856473716896d0*x21 + 0.99709
     &173929518d0* x43 + 0.400214317499682d0*x50
      x227 = x224*x225*x226
      x228 = x223*x227
      x229 = 16.0819842214574d0*x228
      x230 = rho_c**(-6)
      x231 = rho_s**5
      x232 = x167*x231
      x233 = x230*x232
      x234 = 4.0d0*x167
      x235 = x211*x70*x77
      x236 = 1d0/x80
      x237 = x79**(-2)
      x238 = 0.127434173008489d0*x12 + 0.749197287859205d0*x21 + 1.35956
     &911724794d0* x43 + 0.215048623563832d0*x50
      x239 = x236*x237*x238
      x240 = x235*x239
      x241 = 29.6085746432167d0*x240
      x242 = x21*x62
      x243 = 0.00132103950518672d0*x242
      x244 = 1d0/x61
      x245 = x60**(-2)
      x246 = 0.160391419419213d0*x12 + 1.28158207914907d0*x21 + 1.853395
     &81051578d0* x43 + 0.822366887783805d0*x50
      x247 = x244*x245*x246
      x248 = 0.00396311851556016d0*x21 + 0.0310906908696548d0
      x249 = 32.1639684429149d0*x248
      x250 = x247*x249
      x251 = x207 - x229
      x252 = x210*(-x243 + x250 + x251)
      x253 = x14/x15**5
      x254 = 0.632303919304059d0*x253
      x255 = x172 + x175 - x176 + x178 + x179 - x181
      x256 = x0*x255
      x257 = x28*x32
      x258 = 0.00528235082287567d0*x257
      x259 = x19*x258
      x260 = 0.00428269300567249d0*x21
      x261 = 0.0210198269716521d0*x12
      x262 = 0.000367081058155144d0*x26
      x263 = x25 + x260 + x261 + x262
      x264 = x263*x96
      x265 = x22*x30
      x266 = x265*x31
      x267 = 0.00264117541143784d0*x266
      x268 = x264*x267
      x269 = 0.00440195901906306d0*x32
      x270 = x171*x269
      x271 = x0**2.0d0*x31
      x272 = 0.000410930386996015d0*x271
      x273 = x171*x272
      x274 = x174*x39
      x275 = x2**(-0.333333333333333d0)
      x276 = x3**(-0.333333333333333d0)
      x277 = x275 - x276
      x278 = 0.0506605918211689d0*x277
      x279 = x274*x278
      x280 = x1*x279
      x281 = 0.00842941233263362d0*x87
      x282 = x169*x281
      x283 = x21*x45
      x284 = 2.85237982849532d0*x283
      x285 = mu*x21
      x286 = x285*x40
      x287 = x1*x277
      x288 = x286*x287
      x289 = 11.4095193139813d0*x288
      x290 = 1.53482843253593d0*x43
      x291 = 6.13931373014372d0*x43
      x292 = x1*x44
      x293 = x277*x292
      x294 = x291*x293
      x295 = -x290 + x294
      x296 = -x284 + x289 + x295
      x297 = 6.16175563744611d0*x283
      x298 = x41*x50
      x299 = 7.6562840969159d0*x298
      x300 = 24.6470225497844d0*x288
      x301 = x13*x52
      x302 = x287*x301
      x303 = 30.6251363876636d0*x302
      x304 = x295 - x297 - x299 + x300 + x303
      x305 = -x296*x54 + x304
      x306 = 0.0253302959105844d0*x173
      x307 = x305*x306
      x308 = 1d0/x53
      x309 = mu*x308
      x310 = x309*x39
      x311 = x307*x310
      x312 = x208*x231
      x313 = x220*x312
      x314 = -x207 + x229
      x315 = x243 - x250 + x314
      x316 = x315*x74
      x317 = x213*x74
      x318 = x1*x221
      x319 = -x215 + x241 + x251 + x316 + x317 + x318
      x320 = -x167*x313 + 4.0d0*x168 + x319
      x321 = x320*x88
      x322 = -0.256556487706084d0*x12 + 1.02357830998401d0*x21
      x323 = x94*x97
      x324 = x322*x323
      x325 = x0**(-0.666666666666667d0)
      x326 = x92*x93*x96
      x327 = x325*x326
      x328 = 2.41662179562687d0*x327
      x329 = 75.398223686155d0*x95
      x330 = x19*x329
      x331 = 5.74246800037638d0*x99
      x332 = x100*x113
      x333 = x331*x332
      x334 = x115*x127
      x335 = x331*x334
      x336 = 1d0/x21
      x337 = 0.0508408867489862d0*x336
      x338 = 1d0/x103
      x339 = x110*x338
      x340 = x1*x102
      x341 = -x340 + 1.0d0
      x342 = x100*x341
      x343 = x339*x342
      x344 = x337*x343
      x345 = 1d0/x118
      x346 = x124*x345
      x347 = x1*x117
      x348 = x347 + 1.0d0
      x349 = x115*x348
      x350 = x346*x349
      x351 = x337*x350
      x352 = 0.112523287330675d0*x104
      x353 = 0.0162903215389562d0*x12
      x354 = x107*x353
      x355 = x352 + x354
      x356 = -x340*x352 - x340*x354 + x355
      x357 = x109**(-2)
      x358 = x111*x357
      x359 = x106*x358
      x360 = x356*x359
      x361 = x360*x99
      x362 = x101*x361
      x363 = 0.112523287330675d0*x119
      x364 = x122*x353
      x365 = x363 + x364
      x366 = x347*x363 + x347*x364 + x365
      x367 = x123**(-2)
      x368 = x125*x367
      x369 = x121*x368
      x370 = x366*x369*x99
      x371 = x116*x370
      x372 = 22.9698720015055d0*x99
      x373 = x127*x3
      x374 = x372*x373
      x375 = x1*x374
      x376 = x113*x2
      x377 = x372*x376
      x378 = x1*x377
      x379 = x333 + x335 + x344 + x351 + x362 + x371 + x375 - x378
      x380 = x324 + x328 + x330 + x379 + x98
      x381 = 0.000132058770571892d0*x30
      x382 = x22*x23
      x383 = 0.0445553741830535d0*x382
      x384 = 0.00693220036674015d0*x26
      x385 = 0.0891107483661069d0*x89
      x386 = 0.000220097950953153d0*x129 - x19*x385
      x387 = x49*(-x171*x384 + x177 + x264*x383 - x380*x381 + x386)
      x388 = 0.00496533372342768d0*x133
      x389 = x19*x388
      x390 = 0.0041377781028564d0*x134
      x391 = x21*x27
      x392 = 0.00214134650283625d0*x21
      x393 = 0.0105099134858261d0*x12
      x394 = 0.00152446164353214d0*x23
      x395 = 0.000183540529077572d0*x26
      x396 = -x392 - x393 + x394 - x395
      x397 = 0.0777930888680181d0*x391 + x396
      x398 = x397*x96
      x399 = 0.00496533372342768d0*x265
      x400 = x398*x399
      x401 = 0.106650135484359d0*x137
      x402 = 4.0d0*x75
      x403 = x313*x69
      x404 = x319 - x402 + x403
      x405 = 0.0533250677421794d0*x136
      x406 = x404*x405
      x407 = x131*(x389 - x390 + x400 + x401 - x406)
      x408 = 0.0837628205355044d0*x140
      x409 = x19*x408
      x410 = x382*x398
      x411 = 0.0837628205355044d0*x410
      x412 = x153*x154
      x413 = 0.0003448148419047d0*x412
      x414 = x26*x86
      x415 = 0.449784277163183d0*x414
      x416 = x156*x404
      x417 = x99*(x143 + x145)
      x418 = 5.74246800037638d0*x417
      x419 = x3**1.66666666666667d0
      x420 = x2**1.66666666666667d0
      x421 = x99*(x419 - x420)
      x422 = x1*x421
      x423 = 22.9698720015055d0*x422
      x424 = 2.64527879654103d0*x151*x336
      x425 = x151*x99
      x426 = 15.591108358883d0*x425
      x427 = x149*x99
      x428 = -0.173432185689313d0*x12 + 0.0802319968229891d0*x21
      x429 = x428*x96
      x430 = x427*x429
      x431 = 23.3866625383245d0*x430
      x432 = x150*x99
      x433 = 46.773325076649d0*x432
      x434 = x19*x433
      x435 = -x333 - x335 - x344 - x351 - x362 - x371 - x375 + x378
      x436 = x418 + x423 - x424 - x426 - x431 - x434 + x435
      x437 = 0.00020688890514282d0*x30
      x438 = x436*x437
      x439 = x139*(x142 + x409 + x411 + x413 + x415 - x416 - x438)
      x440 = x0*x16*(x259 - x268 - x270 + x273 - x280 + x282 + x311 - x3
     &21 - x387 + x407 + x439)
      x441 = 0.00777267271741373d0*rho_s
      x442 = x3**(-0.333333333333333d0)
      x443 = x2**(-0.333333333333333d0)
      x444 = x442 - 1.0d0*x443
      x445 = x17*x444
      x446 = x441*x445
      x447 = -x163 + 1.0d0
      x448 = x205/x204
      x449 = x166*x448**5.56d0
      x450 = x191*x448**2.78d0
      x451 = x7/x447
      x452 = -x450*x451 + 1.0d0
      x453 = 1d0/(x11*x449/x447**2 + x452)
      x454 = x450*x452
      x455 = x453*x454
      x456 = x455*x7 + 1.0d0
      x457 = log(x456)
      x458 = x457*x5
      x459 = 0.00388633635870686d0*x4
      x460 = x197*x459
      x461 = gamma_cc*x460
      x462 = 1d0/x456
      x463 = rho_c**(-3.33333333333333d0)
      x464 = 1.27130925834002d0*x463
      x465 = x195*x464
      x466 = rho_c**(-4.33333333333333d0)
      x467 = x195*x466
      x468 = 0.726462433337157d0*rho_s
      x469 = 1d0/x4
      x470 = -x443
      x471 = x442 + x470
      x472 = x469*x471
      x473 = x468*x472
      x474 = x467*x473
      x475 = 1.51467417350797d0*x463
      x476 = 1d0/x184
      x477 = -x259
      x478 = x268 + x270 - x273 + x280 - x282 - x311 + x321 + x387 - x40
     &7 - x439 + x477
      x479 = -x16*x478 - x182*x254 + x320
      x480 = -x185*x404 + x479
      x481 = x476*x480
      x482 = x195*x481
      x483 = x475*x482
      x484 = 1.27130925834002d0*x463
      x485 = -x442 + x443
      x486 = x466*x469
      x487 = x468*x485*x486
      x488 = 1d0/x159
      x489 = x215 - x241 + x314 - x316 - x317 - x318 + x402 - x403
      x490 = 1d0/x86
      x491 = x159*x490
      x492 = x24 - x260 - x261 - x262
      x493 = x20*x492
      x494 = -x275 + x276
      x495 = x39*x58
      x496 = x1*x495
      x497 = x290 - x294
      x498 = x284 - x289 + x497
      x499 = x304 + x498*x54
      x500 = 0.0253302959105844d0*x57
      x501 = x310*x500
      x502 = -x330
      x503 = -x324 - x328 + x435 + x502 - x98
      x504 = x20*x383
      x505 = x131*(-x389 + x390 - x400 - x401 + x406) + x139*(x141 - x40
     &9 - x411 - x413 - x415 + x416 + x438) + x267*x493 - x269*x29 + x27
     &2*x29 - x281*x86 + x404*x88 + x477 + x49*(x29*x384 + x381*x503 + x
     &386 + x492*x504 - x91) + 0.0506605918211689d0*x494*x496 + x499*x50
     &1
      x506 = x16*x505
      x507 = x157*x254
      x508 = x404 - x506 - x507
      x509 = x489*x491 + x508
      x510 = x488*x509
      x511 = x475*x510
      x512 = x158 + x203 - x75
      x513 = exp(x162*x512)
      x514 = x513 - 1.0d0
      x515 = 1d0/x514
      x516 = x161*x515
      x517 = 0.544846825002867d0*x516
      x518 = 1.0d0*x442
      x519 = x443 - x518
      x520 = 514.623495086637d0*x469
      x521 = x512*x520
      x522 = x1*x521
      x523 = x157*x253
      x524 = x219*x312
      x525 = x524*x69
      x526 = 1029.24699017327d0*x212
      x527 = 659.974115594947d0*x216*x219
      x528 = -x1*x527 - 0.707034132827743d0*x206 + 0.116918157438755d0*x
     &214 + 4138.08346398728d0*x228 - 7618.63408371287d0*x240 - 257.3117
     &47543319d0*x316 + 257.311747543319d0*x506 + 162.699226454617d0*x52
     &3 - 659.974115594947d0*x525 - x526*x74 + 1029.24699017327d0*x75
      x529 = x519*x522 + x528
      x530 = x513*x529
      x531 = x517*x530
      x532 = x463*x531
      x533 = x484 - x487 - x511 + x532
      x534 = x491**2.78d0
      x535 = x191*x534
      x536 = x515*x7
      x537 = x535*x536
      x538 = x533*x537
      x539 = rho_c**(-5.66666666666667d0)
      x540 = 1.38533762600663d0*x539
      x541 = x186*x194
      x542 = x4**(-5)
      x543 = 0.593716125431411d0*x542
      x544 = x159*x520
      x545 = x470 + x518
      x546 = x1*x545
      x547 = x528 + x544*x546
      x548 = x186*x539
      x549 = gamma_cc*x548
      x550 = x163*x165*x547*x549
      x551 = 0.791621500575215d0*rho_s
      x552 = rho_c**(-6.66666666666667d0)
      x553 = x161*x193*x471
      x554 = gamma_cc*x553
      x555 = x186*x554
      x556 = x552*x555
      x557 = x539*x541
      x558 = x481*x557
      x559 = x189*x463
      x560 = x481*x559
      x561 = x471*x486
      x562 = x189*x561
      x563 = x163*x193
      x564 = x161*x563
      x565 = x547*x564
      x566 = x559*x565
      x567 = x189*x484 + x468*x562 - 1.51467417350797d0*x560 + 0.5448468
     &25002867d0* x566
      x568 = x540*x541 + x543*x550 + x551*x556 - 1.65053082869932d0*x558
     & + x567
      x569 = x197*x568
      x570 = x194*x569
      x571 = x191*x195*x570
      x572 = x189*x462*(x465 + x474 - x483 + x538 - x571)
      x573 = rho_s**3
      x574 = x573*x70
      x575 = x213*x574*x71
      x576 = x0*x221
      x577 = x167*x220
      x578 = x208*x577*x72
      x579 = x234*x574
      x580 = x579*x71
      x581 = x1*x258
      x582 = rho_c**(-3)
      x583 = x574*x582
      x584 = x213*x583
      x585 = -x584
      x586 = x577*x73
      x587 = x579*x582
      x588 = -x221 + x585 + x586 - x587
      x589 = x309*x36
      x590 = 24.6470225497844d0*x283 + x291 + 30.6251363876636d0*x298
      x591 = x277*x590
      x592 = 11.4095193139813d0*x283 + x291
      x593 = x54*x592
      x594 = x494*x593 + x591
      x595 = x306*x594
      x596 = x1*x388
      x597 = x1*x385
      x598 = x1*x329
      x599 = -x598
      x600 = x2*x339
      x601 = x337*x600
      x602 = x3*x346
      x603 = x337*x602
      x604 = 8.61370200056458d0*x99
      x605 = x355*x359
      x606 = x2*x605
      x607 = x604*x606
      x608 = x365*x369
      x609 = x3*x608
      x610 = x604*x609
      x611 = -x374 + x377 + x601 - x603 + x607 - x610
      x612 = x599 + x611
      x613 = x1*x408
      x614 = x99*(-x419 + x420)
      x615 = 22.9698720015055d0*x614
      x616 = x1*x433
      x617 = x611 - x615 - x616
      x618 = x131*(x405*x588 + x596) + x139*(x156*x588 - x437*x617 + x61
     &3) - x279 + x49*(-x381*x612 + x597) + x581 + x588*x88 + x589*x595
      x619 = x16*x618
      x620 = x0*x619
      x621 = x0*x545
      x622 = 0.00777267271741373d0*x621
      x623 = x201*x5
      x624 = x622*x623
      x625 = x469*(0.726462433337157d0*x442 - 0.726462433337157d0*x443)
      x626 = x212*x583
      x627 = -6.05869669403189d0*x626
      x628 = 3.88495961657227d0*x219
      x629 = x167*x73
      x630 = x167*x583
      x631 = -x217 + x218
      x632 = 2.56488140124205d0*x631
      x633 = x216*x632
      x634 = x69*x73
      x635 = x632*x634
      x636 = x583*x69
      x637 = 4.0d0*x636
      x638 = x585 + x633 + x635 + x637
      x639 = x185*x638
      x640 = 359.584801465982d0*x219
      x641 = x444*x469
      x642 = x159*x641
      x643 = x476*(-x216*x628 - 1.51467417350797d0*x619 + x627 + x628*x6
     &29 - 6.05869669403189d0*x630 - 1.51467417350797d0*x639) - x564*(x2
     &16* x640 + 140.195488684916d0*x619 + 560.781954739666d0*x626 - x62
     &9* x640 + 560.781954739666d0*x630 - 280.390977369833d0*x642) + x62
     &5
      x644 = x557*x643
      x645 = 0.726462433337157d0*x472
      x646 = x195*x559
      x647 = x588 - x619
      x648 = -x639 + x647
      x649 = x476*x648
      x650 = 1.51467417350797d0*x649
      x651 = x646*x650
      x652 = gamma_cc*x516
      x653 = x471*x652
      x654 = x491**5.56d0
      x655 = x539*x654
      x656 = x653*x655
      x657 = 0.791621500575215d0*x656
      x658 = x584 - x633 - x635 - x637
      x659 = x491*x658
      x660 = x500*x589
      x661 = x590 - x593
      x662 = x494*x661
      x663 = x374 - x377 - x601 + x603 - x607 + x610
      x664 = x598 + x663
      x665 = x615 + x616 + x663
      x666 = x131*(x405*x638 + x596) + x139*(x156*x638 + x437*x665 + x61
     &3) + x278* x495 + x49*(x381*x664 + x597) + x581 + x638*x88 + x660*
     &x662
      x667 = x16*x666
      x668 = x638 - x667
      x669 = x659 + x668
      x670 = x488*x669
      x671 = x536*x654
      x672 = x539*x671
      x673 = 1.65053082869932d0*x672
      x674 = x670*x673
      x675 = x526*x583
      x676 = 659.974115594947d0*x631
      x677 = 514.623495086637d0*x641
      x678 = -x216*x676 + x512*x677 - x634*x676 - 1029.24699017327d0*x63
     &6 + 257.311747543319d0*x667 + x675
      x679 = x513*x678
      x680 = gamma_cc*x679
      x681 = x514**(-2)
      x682 = x543*x681
      x683 = x655*x682
      x684 = x680*x683
      x685 = x463*x534
      x686 = x645*x685
      x687 = x488*x534
      x688 = x669*x687
      x689 = x475*x688
      x690 = x517*x679
      x691 = x685*x690
      x692 = x686 + x689 - x691
      x693 = x657 + x674 - x684 + x692
      x694 = x199*x536
      x695 = x645*x646 + x651 - x693*x694
      x696 = 0.544846825002868d0*x644 + x695
      x697 = 1d0/x200
      x698 = x461*x697
      x699 = x696*x698
      x700 = x187*x194
      x701 = x453*(-x194*x455*(x192 + 0.593716125431411d0*x700) - 0.2968
     &58062715705d0 *x449*x451 + x454)
      x702 = x459*x701
      x703 = 0.00549554491451047d0*x206
      x704 = 8.0d0*x210
      x705 = 0.000908766572494496d0*x214
      x706 = 5.1297628024841d0*x219
      x707 = x216*x706
      x708 = 32.1639684429148d0*x228
      x709 = 59.2171492864334d0*x240
      x710 = 1.26460783860812d0*x253
      x711 = 0.0155453454348275d0*rho_s*x445
      x712 = gamma_cc*x572
      x713 = 0.00366369660967365d0*x206
      x714 = x209*x230
      x715 = x21*x81
      x716 = 0.00363506628997799d0*x715
      x717 = rho_c**(-7)
      x718 = x219*x231*x717
      x719 = 20.5190512099364d0*x212
      x720 = 20.0d0*x714
      x721 = x21*x82
      x722 = x219*x721
      x723 = rho_s*x582
      x724 = 0.000605844381662998d0*x214
      x725 = x21*x227
      x726 = 0.0883792666034679d0*x725
      x727 = x223*x224
      x728 = 16.0819842214574d0*x727
      x729 = x225*(0.210778258416395d0*x12 + 0.989141964955861d0*x21 + 1
     &.16327369584438d0*x43 + 0.600321476249523d0*x50)
      x730 = x728*x729
      x731 = x223*x226**2/(x65**4*x66**2)
      x732 = 258.630216499205d0*x731
      x733 = x2**(-0.666666666666667d0)
      x734 = 0.854960467080683d0*x1
      x735 = x733*x734
      x736 = x3**(-0.666666666666667d0)
      x737 = x734*x736
      x738 = x735 + x737
      x739 = -5.1297628024841d0*x217 + 5.1297628024841d0*x218 + x738
      x740 = x232*x717
      x741 = 20.5190512099364d0*x219
      x742 = x211*x239
      x743 = x226*(0.252933910099675d0*x12 + 1.48371294743379d0*x21 + 1.
     &99418347859036d0*x43 + x64)/x65**3
      x744 = x728*x743
      x745 = x239*x77
      x746 = x21*x70*x745
      x747 = 0.0269072828949635d0*x746
      x748 = x211*x745
      x749 = x219*x748
      x750 = x235*x236
      x751 = 29.6085746432167d0*x750
      x752 = x237*(0.212390288347482d0*x12 + 0.998929717145607d0*x21 + 1
     &.5861639701226d0*x43 + 0.322572935345748d0*x50)
      x753 = x751*x752
      x754 = x235*x238**2/(x79**4*x80**2)
      x755 = 876.667692402934d0*x754
      x756 = x238*(0.254868346016978d0*x12 + 1.49839457571841d0*x21 + 2.
     &71913823449588d0*x43 + x78)/x79**3
      x757 = x751*x756
      x758 = 5.1297628024841d0*x315
      x759 = x244*x249
      x760 = x713 - x726 - x730 - x732 + x744
      x761 = 0.0849797459133392d0*x21*x247 - 0.00176138600691563d0*x242 
     &+ x245*x759*( 0.267319032365355d0*x12 + 1.70877610553209d0*x21 + 2
     &.16229511226841d0*x43 + 1.23355033167571d0*x50) + 1034.52086599683
     &d0*x246**2*x248/(x60**4*x61**2) - x246*x759*( 0.320782838838426d0*
     &x12 + 2.56316415829814d0*x21 + 3.70679162103156d0*x43 + x59)/x60**
     &3 + x760
      x762 = x182*x253
      x763 = x139*x26/x15**6
      x764 = x18*x71
      x765 = x444*x471
      x766 = 0.0103635636232183d0*x201*x4
      x767 = x765*x766
      x768 = x2**(-1.33333333333333d0)
      x769 = 0.333333333333333d0*x1
      x770 = x768*x769
      x771 = x3**(-1.33333333333333d0)
      x772 = x769*x771
      x773 = x770 + x772
      x774 = 0.000440195901906306d0*x30
      x775 = x323*(0.427594146176807d0*x12 - 1.36477107997868d0*x21)
      x776 = x322*x96
      x777 = 75.398223686155d0*x94
      x778 = x776*x777
      x779 = 0.154912426780983d0*x0**(-0.333333333333333d0)*x326
      x780 = x325*x93
      x781 = 4.83324359125373d0*x780
      x782 = x776*x781
      x783 = 1.61108119708458d0*x327
      x784 = 150.79644737231d0*x19*x322*x94
      x785 = 9.66648718250747d0*x19*x780*x92
      x786 = 1d0/x100
      x787 = 0.00786976946662483d0*x19
      x788 = x336*x339
      x789 = x101*x788*(0.0157395389332497d0*x340 - x786*x787 - 0.007869
     &76946662483d0 )
      x790 = 0.0677878489986483d0*x336
      x791 = x343*x790
      x792 = x350*x790
      x793 = 1.91415600012546d0*x99
      x794 = x332*x793
      x795 = x334*x793
      x796 = 1d0/x115
      x797 = x336*x346
      x798 = x116*x797*(0.0157395389332497d0*x347 + x787*x796 + 0.007869
     &76946662483d0 )
      x799 = 11.4849360007528d0*x100*x361
      x800 = 11.4849360007528d0*x115*x370
      x801 = 0.0543010717965206d0*x12
      x802 = 0.150031049774234d0*x104
      x803 = x19*x786
      x804 = 0.0271505358982603d0*x12
      x805 = x107*x804
      x806 = -x802 - x805
      x807 = x101*x99
      x808 = x359*x807*(0.300062099548467d0*x104*x340 + x107*x340*x801 -
     & x802*x803 - x803*x805 + x806)
      x809 = 0.101681773497972d0*x336
      x810 = x338*x357
      x811 = x356*x810
      x812 = x342*x809*x811
      x813 = x345*x367
      x814 = x809*x813
      x815 = x349*x366*x814
      x816 = 38.2831200025092d0*x99
      x817 = x113*x816
      x818 = x127*x816
      x819 = 0.150031049774234d0*x119
      x820 = x19*x796
      x821 = x122*x804
      x822 = x819 + x821
      x823 = x116*x99
      x824 = x369*x823*(0.300062099548467d0*x119*x347 + x122*x347*x801 +
     & x819*x820 + x820*x821 + x822)
      x825 = 0.225046574661351d0*x104
      x826 = 0.0325806430779124d0*x12
      x827 = x107*x826
      x828 = x825 + x827
      x829 = x109**(-3)
      x830 = x106*x111*x829
      x831 = x356*x807*x830*(-x340*x825 - x340*x827 + x828)
      x832 = 0.225046574661351d0*x119
      x833 = x122*x826
      x834 = x832 + x833
      x835 = x123**(-3)
      x836 = x121*x125*x835
      x837 = x366*x823*x836*(x347*x832 + x347*x833 + x834)
      x838 = x348*x602
      x839 = 0.271151395994593d0*x1
      x840 = x336*x838*x839
      x841 = 15.3132480010037d0*x1*x99
      x842 = x376*x841
      x843 = x336*x341*x600
      x844 = x839*x843
      x845 = x373*x841
      x846 = x3*x366
      x847 = 45.9397440030111d0*x99
      x848 = x1*x847
      x849 = x369*x846*x848
      x850 = x2*x360
      x851 = x848*x850
      x852 = x19*x817 + x19*x818 + x789 + x791 + x792 - x794 - x795 - x7
     &98 + x799 + x800 + x808 + x812 + x815 - x824 + x831 + x837 + x840 
     &+ x842 - x844 - x845 + x849 - x851
      x853 = 0.0231073345558005d0*x26
      x854 = x22*x264
      x855 = 0.0138644007334803d0*x26
      x856 = 0.00107855455836145d0*x30
      x857 = 0.0350330449527535d0*x12 + 0.00571025734089666d0*x21 - 0.00
     &609784657412857d0*x23 + 0.000856522469028669d0*x26
      x858 = x857*x96
      x859 = 0.0891107483661069d0*x382
      x860 = x19*x382
      x861 = 0.178221496732214d0*x860
      x862 = x19*x28
      x863 = -0.000586927869208408d0*x129 + 0.445553741830535d0*x19*x89 
     &- 0.0277288014669606d0*x26*x862
      x864 = x312*x739
      x865 = x167*x864
      x866 = x167*x524
      x867 = x1*x739
      x868 = x216*x867
      x869 = x1*x722
      x870 = x1*x219
      x871 = x748*x870
      x872 = x74*x761
      x873 = 8.0d0*x316
      x874 = x524*x758
      x875 = x212*x74
      x876 = 20.0d0*x875
      x877 = x524*x719
      x878 = x716*x74
      x879 = x74*x742
      x880 = 236.868597145733d0*x879
      x881 = -x713 + x724 + x726 + x730 + x732 - x744 - x747 - x753 - x7
     &55 + x757 + x872 - x873 + x874 - x876 + x877 - x878 + x880
      x882 = x868 - 0.00233087847986162d0*x869 + 151.884964839347d0*x871
     & + x881
      x883 = -20.0d0*x168 - x865 + 20.5190512099364d0*x866 + x882
      x884 = 0.0006896296838094d0*x30
      x885 = 0.899568554326366d0*x26
      x886 = 1.91415600012546d0*x417
      x887 = 0.299209d0*x151
      x888 = 5.29055759308205d0*x336
      x889 = x149*x888
      x890 = x429*x889
      x891 = 31.182216717766d0*x430
      x892 = 23.3866625383245d0*x427*x96*(0.289053642815522d0*x12 - 0.10
     &6975995763985d0*x21)
      x893 = 5.19703611962767d0*x425
      x894 = 5.0d0*x1
      x895 = x34*x894
      x896 = x35*x894
      x897 = x895 + x896
      x898 = 7.65662400050184d0*x99
      x899 = x1*x898*(-6.0d0*x419 + 6.0d0*x420 + x897)
      x900 = 30.6264960020074d0*x422
      x901 = 10.5811151861641d0*x150*x19*x336
      x902 = x427*x428
      x903 = 93.5466501532981d0*x19*x902
      x904 = 77.9555417944151d0*x19*x432
      x905 = 0.311172355472072d0*x21
      x906 = 0.0121035293512547d0*x12*x27 + 0.0175165224763768d0*x12 + 0
     &.00285512867044833d0*x21 + x25 + 0.000428261234514334d0*x26 - 0.10
     &3724118490691d0*x391
      x907 = 0.0837628205355044d0*x382*x96
      x908 = 0.418814102677522d0*x140*x19 - 0.335051282142018d0*x397*x86
     &0 - x408*x96 + 0.167525641071009d0*x410
      x909 = 0.0314471135817087d0*x133*x19
      x910 = 0.0110340749409504d0*x134
      x911 = x265*x397
      x912 = 0.0198613348937107d0*x19*x911
      x913 = 0.0165511124114256d0*x265*x398
      x914 = x396*x905 + x906
      x915 = x399*x914*x96
      x916 = 0.319950406453076d0*x137
      x917 = 0.213300270968718d0*x136*x404
      x918 = x69*x864
      x919 = 20.0d0*x75
      x920 = 20.5190512099364d0*x525
      x921 = x918 + x919 - x920
      x922 = x882 + x921
      x923 = x405*x922
      x924 = 0.0309078452196566d0*x87
      x925 = 0.00316102962473761d0*x87
      x926 = 0.0168588246652672d0*x87
      x927 = 0.000821860773992031d0*x271
      x928 = 0.00191767513931474d0*x271
      x929 = 1.79063317129192d0*x43
      x930 = 3.80317310466043d0*x283
      x931 = 6.0d0*x275
      x932 = 6.0d0*x276
      x933 = x2**(-1.33333333333333d0)
      x934 = x1*x933
      x935 = x3**(-1.33333333333333d0)
      x936 = x1*x935
      x937 = x934 + x936
      x938 = -x931 + x932 + x937
      x939 = x1*x938
      x940 = x286*x939
      x941 = 3.80317310466043d0*x940
      x942 = x277**2
      x943 = x19*x942
      x944 = x285*x51
      x945 = 22.8190386279626d0*x943*x944
      x946 = 7.60634620932085d0*x288
      x947 = 2.04643791004791d0*x43
      x948 = x292*x938*x947
      x949 = 8.18575164019163d0*x43
      x950 = x942*x949
      x951 = x19*x40*x950
      x952 = x293*x947
      x953 = x948 + x951 - x952
      x954 = x296*x45
      x955 = 2.0d0*x48
      x956 = x304*x955
      x957 = 3.06965686507186d0*x43
      x958 = 5.70475965699064d0*x283
      x959 = 12.2786274602874d0*x43
      x960 = x293*x959
      x961 = 22.8190386279626d0*x288
      x962 = x53/x47**2
      x963 = x954*x962
      x964 = 8.21567418326147d0*x283
      x965 = 49.2940450995689d0*x942
      x966 = 81.6670303671029d0*x13*x50*x943/x36**4 + x19*x944*x965 - 16
     &.4313483665229d0*x288 + 11.4844261453738d0*x298 + 10.2083787958879
     &d0*x301*x939 - 30.6251363876636d0*x302 + 1.79063317129192d0*x43 + 
     &8.21567418326147d0*x940 + x953 + x964
      x967 = 6.39350882283002d-5*x0**2.33333333333333d0*x31
      x968 = 0.00880391803812612d0*x266
      x969 = 0.0117385573841682d0*x32
      x970 = x13/x53**2
      x971 = x36*x970
      x972 = x263*x266
      x973 = 0.0105647016457513d0*x19
      x974 = 0.0675474557615585d0*x36
      x975 = x942*x974
      x976 = x174*x975
      x977 = x13*x308
      x978 = x48*x977
      x979 = x307*x978
      x980 = x296*x36
      x981 = x173*x305
      x982 = 0.101321183642338d0*x589
      x983 = 0.0334548885448792d0*x19*x257 - 0.00164372154798406d0*x271*
     &x862
      x984 = x16*(-x131*(-x909 + x910 + x912 - x913 + x915 - x916 + x917
     & - x923) + x139*(-0.0009195062450792d0*x154*x180 + x156*x883 + 1.0
     &4949664671409d0*x169*x26 - x320*x885 - x437*(x852 + x886 + x887 + 
     &x890 + x891 + x892 - x893 - x899 - x900 + x901 + x903 - x904) + x8
     &84*(x379 - x418 - x423 + x424 + x426 + x431 + x434) - x907*(-x905*
     &(x392 + x393 - x394 + x395) + x906) + x908) + x169* x924 + x171*x9
     &28 - x171*x967 - x171*x969 - x19*x976 - x264*x968 - x267*x858 + 0.
     &0168868639403896d0*x274*x939 + x287*x981*x982 + x304*x307*x971 - x
     &306*x310*(-x54*(x929 + x930 + x941 + x945 - x946 + x953) - x954*x9
     &56 + x963*(-x957 - x958 + x960 + x961) + x966) - x320*x926 + x49*(
     &x171*x853 - x171*x856 + x263*x861 - x264 *x859 + x380*x774 - x381*
     &(x502 + x775 + x778 + x779 + x782 + x783 + x784 + x785 + x852) - x
     &383*x858 - x385*x96 + x854*x855 + x863) + x854*x927 + x883*x925 + 
     &x972*x973 - x979*x980 + x983)
      x985 = x197*x697
      x986 = x189*(-x465 - x474 + x483 - x538 + x571)
      x987 = gamma_cc*x985*x986
      x988 = x537 + 1.0d0
      x989 = x534*x988
      x990 = x552*x654
      x991 = x652*x990
      x992 = gamma_cc*x530
      x993 = x535*x988
      x994 = x166*x654
      x995 = x11*x681*x994 + x988
      x996 = 1d0/x995
      x997 = x536*x996
      x998 = x993*x997
      x999 = 0.00388633635870686d0*x469
      x1000 = x996*x999
      x1001 = x1000*x197*x8/x200**2
      x1002 = x189*x467
      x1003 = rho_c**(-5.33333333333333d0)*rho_s
      x1004 = x189*x195
      x1005 = x1003*x1004
      x1006 = x471**2
      x1007 = 1.45292486667431d0*x6
      x1008 = x1006*x1007
      x1009 = rho_c**(-6.33333333333333d0)*x18
      x1010 = x1008*x1009
      x1011 = 2.66666666666667d0*x442
      x1012 = 2.66666666666667d0*x443
      x1013 = 0.444444444444444d0*x1
      x1014 = x1013*x768
      x1015 = x1013*x771
      x1016 = x1014 + x1015
      x1017 = 0.544846825002868d0*x469
      x1018 = x1003*x1017*(x1011 - x1012 + x1016)
      x1019 = x1002*x481
      x1020 = x480**2
      x1021 = x184**(-2)
      x1022 = 4.21079420235216d0*x1021
      x1023 = x1002*x1022
      x1024 = 4.03913112935459d0*x472
      x1025 = x1003*x1024
      x1026 = x170*x320
      x1027 = 1.51467417350797d0*x1019
      x1028 = 1.51467417350797d0*x1021
      x1029 = x1028*x508
      x1030 = x1002*x480
      x1031 = x533*x536
      x1032 = x534*x559
      x1033 = x1031*x1032
      x1034 = x466*x534
      x1035 = 1.45292486667431d0*rho_s
      x1036 = x533*x653
      x1037 = 3.02934834701594d0*x560
      x1038 = x73 - 1.0d0
      x1039 = x1038*x212
      x1040 = x1039*x867
      x1041 = x1038*x870
      x1042 = 0.00116543923993081d0*x1041*x715
      x1043 = 0.00116543923993081d0*x721
      x1044 = x1043*x219
      x1045 = x1*x1044
      x1046 = 75.9424824196734d0*x748
      x1047 = x1046*x870
      x1048 = 75.9424824196734d0*x1041*x742
      x1049 = 2.78d0*x185
      x1050 = x170*x479
      x1051 = 5.56d0*x404
      x1052 = 8.0d0*x875
      x1053 = x1*x707
      x1054 = 5.1297628024841d0*x525
      x1055 = 8.0d0*x75
      x1056 = 2.0d0*x316
      x1057 = 2.78d0*x404
      x1058 = x184/x169**2
      x1059 = 1.38933365612628d0*x763
      x1060 = 3.51560979133057d0*x253
      x1061 = x727*x729
      x1062 = x315*x524
      x1063 = x212*x524
      x1064 = x715*x74
      x1065 = x727*x743
      x1066 = x750*x752
      x1067 = x750*x756
      x1068 = 44.7079161356516d0*x1061 + 14.2607405909058d0*x1062 + 57.0
     &429623636231d0 *x1063 - 0.0101054842861388d0*x1064 - 44.7079161356
     &516d0*x1065 - 82.3118375081424d0*x1066 + 82.3118375081424d0*x1067 
     &- 0.0101850765748927d0*x206 + 0.00168424738102313d0*x214 - 22.24d0
     &* x316 + 0.245694361157641d0*x725 + 718.99200186779d0*x731 - 0.074
     &8022464479985d0*x746 - 2437.13618488016d0*x754 + 2.78d0*x868 - 0.0
     &064798421740153d0*x869 + 422.240202253384d0*x871 + 2.78d0* x872 - 
     &55.6d0*x875 + 658.494700065139d0*x879
      x1069 = -x1049*(-x1040 + x1042 - x1045 + x1047 - x1048 + x881 + x9
     &21) - x1050* x1051 + x1057*x1058*(x1052 + x1053 + x1054 - x1055 + 
     &x1056 + x703 - x705 - x708 + x709) - x1059*x182 - x1060*x478 + x10
     &68 - 55.6d0* x168 + 2.9296748261088d0*x762 - 2.78d0*x865 + 57.0429
     &623636231d0* x866 - 2.78d0*x984
      x1070 = 0.544846825002868d0*x1002*x476
      x1071 = x554*x569
      x1072 = x482*x559
      x1073 = 1.08969365000574d0*x190*x534
      x1074 = x1073*x189
      x1075 = x10*x569
      x1076 = x193*x515
      x1077 = 1.58324300115043d0*rho_s*x556 + 1.18743225086282d0*x542*x5
     &50 + 2.77067525201325d0*x557 - 3.30106165739864d0*x558
      x1078 = x196**(-2)
      x1079 = x10*x1078*x198
      x1080 = x1079*x568
      x1081 = 4.23769752780008d0*x466
      x1082 = 3.39015802224007d0*x1003
      x1083 = x466*x476
      x1084 = x1083*x480
      x1085 = x1022*x466
      x1086 = x466*x564
      x1087 = 514.623495086637d0*x161
      x1088 = 1.08969365000573d0*x165*exp(-x1087*x159)
      x1089 = x4**(-6)
      x1090 = x1089*x466
      x1091 = x1090*x547**2
      x1092 = 0.544846825002867d0*x563
      x1093 = x471*x9
      x1094 = x1093*x563
      x1095 = 1.45292486667431d0*x1094
      x1096 = 1.51467417350797d0*x1084
      x1097 = x466*x480
      x1098 = 1029.24699017327d0*x469
      x1099 = 2.0d0*x442 - 2.0d0*x443 + x773
      x1100 = x471*x545
      x1101 = x19*x6
      x1102 = 1372.32932023103d0*x159
      x1103 = 38.2831200025092d0*x99
      x1104 = x1103*x113
      x1105 = x1103*x127
      x1106 = -x1104*x19 - x1105*x19 - x789 - x791 - x792 + x794 + x795 
     &+ x798 - x799 - x800 - x808 - x812 - x815 + x824 - x831 - x837 - x
     &840 - x842 + x844 + x845 - x849 + x851
      x1107 = x22*x493
      x1108 = x45*x498
      x1109 = x1108*x962
      x1110 = -x934 - x936
      x1111 = x494**2
      x1112 = x58*x974
      x1113 = x499*x500
      x1114 = x266*x492
      x1115 = x499*x57
      x1116 = x16*(x1*x1115*x494*x982 + x1107*x927 + x1111*x1112*x19 + x
     &1113*x971*( x297 + x299 - x300 - x303 + x497) + x1113*x978*x980 - 
     &x1114*x973 + x131*(x909 - x910 - x912 + x913 - x915 + x916 - x917 
     &+ x923) + x139*(x156*x922 - x404*x885 + 0.0009195062450792d0*x412 
     &+ 1.04949664671409d0*x414 - x436*x884 + x437*(x1106 - x886 - x887 
     &- x890 - x891 - x892 + x893 + x899 + x900 - x901 - x903 + x904) - 
     &x907*x914 + x908) + x20*x267*x857 - x29*x928 + x29*x967 + x29* x96
     &9 - x404*x926 + x49*(x1107*x855 + x20*x385 - x29*x853 + x29* x856 
     &+ x381*(x1106 + x330 - x775 - x778 - x779 - x782 - x783 - x784 - x
     &785) - x492*x861 - x493*x859 - x503*x774 + x504*x857 + x863) - x49
     &3*x968 + 0.0168868639403896d0*x496*(x1110 + x931 - x932) + x501*(x
     &1108*x956 + x1109*(x957 + x958 - x960 - x961) + x54*(-x929 - x930 
     &- x941 - x945 + x946 - x948 - x951 + x952) + x966) + x86*x924 + x9
     &22*x925 + x983)
      x1117 = x157*x763
      x1118 = x253*x505
      x1119 = -4138.08346398728d0*x1061 - 1319.94823118989d0*x1062 - 527
     &9.79292475958d0*x1063 + 0.935345259510043d0*x1064 + 4138.083463987
     &28d0*x1065 + 7618.63408371287d0*x1066 - 7618.63408371287d0*x1067 +
     & 257.311747543319d0*x1116 + 128.594198193741d0*x1117 + 325.3984529
     &09234d0*x1118 + 0.942712177103657d0*x206 - 0.155890876585007d0*x21
     &4 + 2058.49398034655d0*x316 - 271.165377424361d0*x523 + 5279.79292
     &475958d0*x525 - 22.7410235363352d0*x725 - 66548.5929749173d0*x731 
     &+ 6.92355998334549d0*x746 - 5146.23495086637d0*x75 + 225576.895946
     &967d0*x754 - 257.311747543319d0*x868 + 0.599762414964306d0*x869 - 
     &39081.7857283678d0*x871 - 257.311747543319d0*x872 + 5146.234950866
     &37d0*x875 - 60949.072669703d0*x879 - 257.311747543319d0*x918
      x1120 = 0.544846825002867d0*x1086
      x1121 = 0.544846825002868d0*x1083
      x1122 = x191*x541
      x1123 = x469*x485
      x1124 = x1123*x534
      x1125 = x485**2
      x1126 = x1007*x534
      x1127 = x536*x990
      x1128 = rho_c**(-7.66666666666667d0)*rho_s*x654
      x1129 = x1128*x652
      x1130 = x1129*x485
      x1131 = x515*x9
      x1132 = 2.63873833525072d0*x1131
      x1133 = 0.888888888888889d0*x1
      x1134 = x1133*x768
      x1135 = x1133*x771
      x1136 = x1034*x510
      x1137 = x159**(-2)
      x1138 = x1137*x509**2
      x1139 = 4.21079420235216d0*x1034
      x1140 = x1034*x516
      x1141 = gamma_cc*x990
      x1142 = x542*x681
      x1143 = x1141*x1142
      x1144 = x1143*x530
      x1145 = x529**2
      x1146 = exp(x1087*x512)
      x1147 = x1145*x1146
      x1148 = x1034*x1089
      x1149 = 1.08969365000573d0*x1148*x681
      x1150 = x1145*x513
      x1151 = 0.544846825002867d0*x1148*x515
      x1152 = x1141/x4**8
      x1153 = 1.78114837629423d0*x1152/x514**3
      x1154 = 0.593716125431411d0*x1152*x681
      x1155 = x1003*x485
      x1156 = x1131*x1155*x534
      x1157 = x1128*x485
      x1158 = x1089*x681
      x1159 = 3.16648600230086d0*x1158
      x1160 = 1.51467417350797d0*x490
      x1161 = x1136*x1160
      x1162 = x1127*x510
      x1163 = 9.17695140756823d0*x1127
      x1164 = 1.65053082869932d0*x490
      x1165 = x1162*x1164
      x1166 = 1.51467417350797d0*x1034
      x1167 = x489 + x506 + x507
      x1168 = x1137*x509
      x1169 = x1167*x1168
      x1170 = 1.65053082869932d0*x1127
      x1171 = x516*x530
      x1172 = x1*x519
      x1173 = 1372.32932023103d0*x512
      x1174 = x513*(x1098*x1167*x1172 + x1099*x522 + x1101*x1173*x485*x5
     &19 + x1119)
      x1175 = x1034*x517
      x1176 = x1141*x682
      x1177 = x491*(x1040 - x1042 + x1045 - x1047 + x1048 - x724 + x747 
     &+ x753 + x755 - x757 + x760 - x872 + x873 - x874 + x876 - x877 + x
     &878 - x880 - x918 - x919 + x920)
      x1178 = x490*x508
      x1179 = 5.56d0*x489
      x1180 = 2.78d0*x489
      x1181 = x159/x86**2
      x1182 = x1181*(-x1052 - x1053 - x1054 + x1055 - x1056 - x703 + x70
     &5 + x708 - x709)
      x1183 = 0.544846825002868d0*x466*x687
      x1184 = 11.12d0*x489
      x1185 = 0.296858062715705d0*x1127*x488
      x1186 = x208*x574
      x1187 = 0.00181753314498899d0*x715
      x1188 = x230*x72
      x1189 = x212*x741
      x1190 = 16.0d0*x1186
      x1191 = 2.56488140124205d0*x218
      x1192 = 2.56488140124205d0*x217
      x1193 = x1191 - x1192 + x738
      x1194 = x1193*x216
      x1195 = x167*x741
      x1196 = x1046*x219
      x1197 = x1193*x167
      x1198 = x220*x315
      x1199 = 4.0d0*x315
      x1200 = x17*x623
      x1201 = 0.000220097950953153d0*x30
      x1202 = x1*x781*x92
      x1203 = x1*x322*x777
      x1204 = 0.135575697997297d0*x1
      x1205 = x1204*x788
      x1206 = 0.0338939244993242d0*x336
      x1207 = x1206*x600
      x1208 = x1204*x797
      x1209 = x1206*x602
      x1210 = x600*x790
      x1211 = x1210*(x340 - 1.0d0)
      x1212 = 0.135575697997297d0*x843
      x1213 = 7.65662400050185d0*x99
      x1214 = x1213*x376
      x1215 = x790*x838
      x1216 = x1213*x373
      x1217 = x355*x810
      x1218 = x2*x337
      x1219 = x1217*x1218*x341
      x1220 = x1*x372
      x1221 = x1220*x605
      x1222 = x331*x606
      x1223 = x337*x813
      x1224 = x1223*x3*x348*x365
      x1225 = x1220*x608
      x1226 = x331*x609
      x1227 = x1218*x811
      x1228 = x1223*x846
      x1229 = x340*x802
      x1230 = x340*x805
      x1231 = x359*x604
      x1232 = x1231*x2*(x1229 + x1230 + x806)
      x1233 = x372*x850
      x1234 = x372*x846
      x1235 = x1234*x369
      x1236 = x347*x819 + x347*x821 + x822
      x1237 = x369*x604
      x1238 = x1236*x1237*x3
      x1239 = x604*x828
      x1240 = x1239*x830
      x1241 = x1240*x2*x356
      x1242 = x604*x834
      x1243 = x1242*x836
      x1244 = x1243*x846
      x1245 = x1*x859
      x1246 = x1*x28
      x1247 = -0.267332245098321d0*x1*x89 + x1246*x855
      x1248 = 3.0d0*x419
      x1249 = 3.0d0*x420
      x1250 = x1*(x147 - x148)
      x1251 = 46.773325076649d0*x1*x902
      x1252 = x120 - 1.0d0
      x1253 = x1252*x3
      x1254 = x105 - 1.0d0
      x1255 = x1254*x2
      x1256 = x1255*x358
      x1257 = x1253*x368
      x1258 = x1252*x368
      x1259 = x802 + x805
      x1260 = x1*x1103
      x1261 = x1205 - x1207 + x1208 + x1209 + x1215 - x1219 + x1224 - x1
     &227 + x1228
      x1262 = x1197*x73
      x1263 = x212*x73
      x1264 = x1199*x583
      x1265 = 16.0d0*x626
      x1266 = x1187*x583
      x1267 = x583*x742
      x1268 = 118.434298572867d0*x1267
      x1269 = x1264 + x1265 + x1266 - x1268
      x1270 = x1044 - x1194 - x1196 - x1198*x73 + x1262 - x1263*x741 + x
     &1269 - x629* x741 + 16.0d0*x630
      x1271 = 0.0003448148419047d0*x30
      x1272 = 0.449784277163183d0*x26
      x1273 = -0.251288461606513d0*x1*x140 + 0.167525641071009d0*x1*x382
     &*x397
      x1274 = 0.106650135484359d0*x136
      x1275 = -0.0182062236525682d0*x1*x133 + 0.00993066744685536d0*x1*x
     &911
      x1276 = 3.0d0*x275
      x1277 = 3.0d0*x276
      x1278 = -x1276 + x1277 + x937
      x1279 = 0.0168868639403896d0*x1278
      x1280 = 1.02321895502395d0*x43
      x1281 = 15.3125681938318d0*x298
      x1282 = x44*x950
      x1283 = x286*x965
      x1284 = 10.2083787958879d0*x298
      x1285 = x48*x494
      x1286 = x1285*x592
      x1287 = -x1286*x304*x45
      x1288 = 22.8190386279626d0*x283 + x959
      x1289 = x1288*x494
      x1290 = x1110 + x1276 - x1277
      x1291 = x1290*x930
      x1292 = x494*x930
      x1293 = x494*x961
      x1294 = x1290*x947
      x1295 = x1280*x494
      x1296 = x293*x494*x949
      x1297 = x1294 - x1295 + x1296
      x1298 = x306*x589
      x1299 = x590*x970
      x1300 = x1299*x494
      x1301 = x278*x589
      x1302 = 0.00528235082287567d0*x1
      x1303 = x277*x592
      x1304 = x173*x594
      x1305 = x1*x309
      x1306 = -0.0193686196838775d0*x1*x257 + x1246*x927
      x1307 = x16*(x1*x976 + x1270*x925 - x1279*x274 - x1298*(-x1*x1282 
     &- x1*x1283 - 81.6670303671029d0*x1*x301*x942 - x1278*x1284 - x1278
     &*x947 - x1278*x964 + x1280*x277 + x1281*x277 + x1287 + x1289*x963 
     &+ x277* x964 + x48*x591*x954 - x54*(x1291 - x1292 + x1293 + x1297)
     &) + x1300*x307 - x1301*x981 - x1302*x972 + x1303*x979 - x1304*x130
     &5* x278 + x1306 + x131*(x1270*x405 - x1274*x588 + x1275) + x139*( 
     &x1270*x156 + x1271*x617 - x1272*x588 + x1273 + x437*(x111*x1239* x
     &1255*x356*x829 - x112*x1213*x1255 - x112*x1254*x1260 - x1210* x341
     & + x1213*x1253*x126 - x1220*x1254*x355*x358 - x1220*x1258* x365 - 
     &x1234*x1258 + x1236*x1257*x604 - x1242*x125*x1252*x835* x846 + 62.
     &364433435532d0*x1250*x427 - x1250*x889 + x1251 - x1252* x126*x1260
     & + x1256*x331*x355 + x1256*x356*x372 - x1256*x604*( -x1229 - x1230
     & + x1259) - x1257*x331*x365 + x1261 - 15.3132480010037d0*x421 - x8
     &98*(-x1248 + x1249 + x897))) - x281* x588 + x49*(x1201*x612 - x124
     &5*x263 + x1247 - x381*(-x1*x1104 - x1*x1105 - x1202 - x1203 - x120
     &5 + x1207 - x1208 - x1209 + x1211 + x1212 - x1214 - x1215 + x1216 
     &+ x1219 - x1221 + x1222 - x1224 - x1225 - x1226 + x1227 - x1228 + 
     &x1232 + x1233 - x1235 + x1238 + x1241 - x1244 + x598)))
      x1308 = gamma_cc*x696
      x1309 = x1308*x985
      x1310 = 3.88495961657227d0*x631
      x1311 = x1001*(x463*x645*x989 + x475*x670*x989 + 0.544846825002868
     &d0*x672*(x488* (x1310*x216 + x1310*x634 + x627 + 6.05869669403189d
     &0*x636 + 1.51467417350797d0*x659 - 1.51467417350797d0*x667) + x625
     & - x690 ) + x998*(-x657 - x674 + x684 - x686 - x689 + x691))
      x1312 = x1003*x1008
      x1313 = 1.33333333333333d0*x442
      x1314 = 1.33333333333333d0*x443
      x1315 = x1016 + x1313 - x1314
      x1316 = 2.0195655646773d0*x472
      x1317 = x1316*x649
      x1318 = x552*x643
      x1319 = x1318*x541
      x1320 = x170*x588
      x1321 = x1028*x668
      x1322 = x480*x648
      x1323 = -x1191 + x1192 - x735 - x737
      x1324 = x1039*x1323
      x1325 = x1323*x634
      x1326 = 16.0d0*x636
      x1327 = x315*x73
      x1328 = x1327*x632
      x1329 = 10.2595256049682d0*x634
      x1330 = x1329*x219
      x1331 = x1043*x631
      x1332 = x1329*x631
      x1333 = 10.2595256049682d0*x1263
      x1334 = x1333*x631
      x1335 = x1333*x219
      x1336 = x1046*x631
      x1337 = x1269 + x1325 - x1326 + x1328 + x1330 - x1331 - x1332 + x1
     &334 - x1335 + x1336
      x1338 = x170*x647
      x1339 = x1327*x219
      x1340 = x219*x629
      x1341 = x253*x618
      x1342 = 57.0429623636231d0*x219
      x1343 = -x1263*x1342
      x1344 = x315*x583
      x1345 = x583*x715
      x1346 = -329.247350032569d0*x1267 + 11.12d0*x1344 + 0.005052742143
     &0694d0*x1345 + 44.48d0*x626
      x1347 = -x1049*(-x1324 + x1337) - 2.78d0*x1050*x638 + x1051*x1058*
     &x638 - x1057* x1338 - 2.78d0*x1194 + 2.78d0*x1262 - 2.78d0*x1307 -
     & 7.13037029545289d0*x1339 - 57.0429623636231d0*x1340 - 1.757804895
     &66528d0*x1341 + x1343 + x1346 + 44.48d0*x630 + 0.00323992108700765
     &d0*x722 - 211.120101126692d0*x749
      x1348 = x197*x693
      x1349 = x189*x465
      x1350 = x1349*x536
      x1351 = x1348*x653
      x1352 = x1002*x468
      x1353 = 1.51467417350797d0*x1072
      x1354 = x1348*x536
      x1355 = x10*x1348
      x1356 = x192*x533*x534*x681
      x1357 = 0.726462433337157d0*x646
      x1358 = x548*x643
      x1359 = 0.544846825002868d0*x1358
      x1360 = x516*x679
      x1361 = 3.30106165739864d0*x670
      x1362 = -1.18743225086282d0*x1142*x655*x680 + x1361*x672 + 1.58324
     &300115043d0* x656
      x1363 = -1.08969365000573d0*x1360*x685 + x1362 + 1.45292486667431d
     &0*x472*x685 + 3.02934834701594d0*x670*x685
      x1364 = x1076*x1080
      x1365 = 1.69507901112003d0*x561
      x1366 = 0.544846825002868d0*x486
      x1367 = x1083*x648
      x1368 = -659.974115594947d0*x1340 + x527 + 257.311747543319d0*x619
     & + 1029.24699017327d0*x630 - 514.623495086637d0*x642 + x675
      x1369 = 0.726462433337157d0*x1094
      x1370 = x1090*x1368*x547
      x1371 = x1368*x564
      x1372 = x1263*x219
      x1373 = x221 + x584 - x586 + x587 + x619
      x1374 = x471*x6
      x1375 = 30474.5363348515d0*x1267 - 1029.24699017327d0*x1344 - 0.46
     &7672629755021d0*x1345 - 4116.9879606931d0*x626
      x1376 = x653*x990
      x1377 = x1093*x515
      x1378 = x466*x688
      x1379 = x1143*x679
      x1380 = 1.58324300115043d0*x1158
      x1381 = x1127*x670
      x1382 = x1034*x1377
      x1383 = x1141*x471
      x1384 = x658 + x667
      x1385 = x1166*x1384
      x1386 = x1170*x1384
      x1387 = x1137*x1139
      x1388 = x1146*x1149
      x1389 = x529*x678
      x1390 = x529*x679
      x1391 = x1146*x1153
      x1392 = x1*x817 + x1*x818 - x1211 - x1212 + x1214 - x1216 + x1221 
     &- x1222 + x1225 + x1226 - x1232 - x1233 + x1235 - x1238 - x1241 + 
     &x1244 + x1261
      x1393 = x1323*x216
      x1394 = x1337 + x1393
      x1395 = x16*(0.0506605918211689d0*x1111*x1305*x57*x661 + x1112*x28
     &7*x494 + x1113 *x1286*x977 - x1113*x1300 + x1114*x1302 + x1115*x13
     &01 + x1279* x495 + x1306 + x131*(-x1274*x638 + x1275 + x1394*x405)
     & + x139*( -x1271*x665 - x1272*x638 + x1273 + x1394*x156 + x437*(x1
     &*x150* x888 - 62.364433435532d0*x1*x432 + x1251 + x1392 + 15.31324
     &80010037d0*x614 + x898*(x1248 - x1249 - x895 - x896))) + x1394*x92
     &5 - x281*x638 + x49*(-x1201*x664 + x1245*x492 + x1247 + x381*(x120
     &2 + x1203 + x1392 + x599)) + x660*(x1108*x1285*x590 - x1109*x1289 
     &- x1281*x494 + x1284*x1290 + x1287 + x1290*x964 + x1297 + 49.29404
     &50995689d0*x288*x494 + 81.6670303671029d0*x302* x494 - x494*x964 +
     & x54*(-x1291 + x1292 - x1293 - x1294 + x1295 - x1296)))
      x1396 = x253*x666
      x1397 = x631*x721
      x1398 = 2639.89646237979d0*x634
      x1399 = x1263*x631
      x1400 = x631*x748
      x1401 = x1167*x677 + x1172*x1173*x1374 + x1172*x1384*x520 - 257.31
     &1747543319d0* x1325 - x1327*x676 + 2639.89646237979d0*x1372 + x137
     &5 - 257.311747543319d0*x1393 + 257.311747543319d0*x1395 + 162.6992
     &26454617d0*x1396 + 0.299881207482153d0*x1397 - x1398*x219 + x1398*
     &x631 - 2639.89646237979d0*x1399 - 19540.8928641839d0* x1400 + x521
     &*(x519 - x770 - x772) + 4116.9879606931d0*x636
      x1402 = x1175*x513
      x1403 = x1176*x513
      x1404 = x491*(-x1264 - x1265 - x1266 + x1268 + x1324 - x1325 + x13
     &26 - x1328 - x1330 + x1331 + x1332 - x1334 + x1335 - x1336)
      x1405 = x1178*x658
      x1406 = x490*x668
      x1407 = x1327*x631
      x1408 = 28.5214811818116d0*x634
      x1409 = x1181*x658
      x1410 = 57.0429623636231d0*x634
      x1411 = x70*x764
      x1412 = 12.0d0*x212
      x1413 = x733 + x736
      x1414 = 0.854960467080683d0*x1413
      x1415 = x1414*x216
      x1416 = x1195*x573
      x1417 = 12.0d0*x167
      x1418 = x933 + x935
      x1419 = 0.0168868639403896d0*x1418
      x1420 = 0.101321183642338d0*x277*x309
      x1421 = x44*x595
      x1422 = x19*x70
      x1423 = x1417*x1422
      x1424 = x1416*x582
      x1425 = x1414*x629
      x1426 = x573*x582
      x1427 = x1189*x1426
      x1428 = x1412*x1422
      x1429 = x1415 - x1428
      x1430 = -x1423 + x1424 - x1425 + x1427 + x1429
      x1431 = x1418*x947
      x1432 = 81.6670303671029d0*x301
      x1433 = mu*x40*x592
      x1434 = x1111*x286
      x1435 = x1111*x44
      x1436 = x1431 + x1435*x949
      x1437 = x1111*x1433
      x1438 = x1284*x1418 + x1288*x1437*x962 + 8.21567418326148d0*x1418*
     &x283 - x54*( x1418*x930 + 22.8190386279626d0*x1434 + x1436)
      x1439 = 0.203363546995945d0*x336
      x1440 = x1439*x339
      x1441 = x1439*x346
      x1442 = x1217*x809
      x1443 = x365*x814
      x1444 = x605*x847
      x1445 = x608*x847
      x1446 = x1231*x1259
      x1447 = x1237*x822
      x1448 = x1240*x355
      x1449 = x1243*x365
      x1450 = x1440 + x1441 + x1442 + x1443 + x1444 + x1445 - x1446 - x1
     &447 + x1448 + x1449 + x817 + x818
      x1451 = x36*x816
      x1452 = x16*(-x1298*(x1282 + x1283 + x1431 + x1432*x942 + x1433*x4
     &94*x591*x955 + x1438) - x1300*x1421 - x1303*x1421*x978 + x1304*x14
     &20 + x131*( x1430*x405 + x388) + x139*(x1430*x156 + x408 - x437*(x
     &1450 - x1451 - x433)) + x1419*x274 + x1430*x925 + x258 + x49*(-x38
     &1*( x1450 - x329) + x385) - x976)
      x1453 = x768 + x771
      x1454 = x1453*x469
      x1455 = x1022*x648**2
      x1456 = x1024*x649
      x1457 = x1320*x650
      x1458 = x1321*x648
      x1459 = -x1415 + x1423 - x1424 + x1425 - x1427 + x1428
      x1460 = x1058*x638**2
      x1461 = x1338*x638
      x1462 = x1426*x167
      x1463 = x1422*x167
      x1464 = 2.3767900984843d0*x1413
      x1465 = x1426*x212
      x1466 = x1422*x212
      x1467 = x1464*x216 - 33.36d0*x1466
      x1468 = 3.02934834701594d0*x649
      x1469 = x1076*x1355
      x1470 = x1079*x681*x693
      x1471 = x1089*x1368**2
      x1472 = x1426*x631
      x1473 = x1472*x212
      x1474 = 119.861600488661d0*x1413
      x1475 = x1472*x69
      x1476 = x1422*x69
      x1477 = x6*x765
      x1478 = x1435*x500*x661
      x1479 = 12.0d0*x1476
      x1480 = 20.5190512099364d0*x1475
      x1481 = x1414*x634
      x1482 = x1472*x719
      x1483 = x1429 + x1479 + x1480 + x1481 - x1482
      x1484 = -x1104 - x1105 - x1440 - x1441 - x1442 - x1443 - x1444 - x
     &1445 + x1446 + x1447 - x1448 - x1449
      x1485 = x16*(-x1299*x1478 + x131*(x1483*x405 + x388) + x139*(x1483
     &*x156 + x408 + x437*(x1451 + x1484 + x433)) - x1419*x495 + x1420*x
     &57*x662 + x1478*x592*x978 + x1483*x925 + x258 + x49*(x381*(x1484 +
     & x329) + x385) + x58*x975 + x660*(x1111*x1432 + 49.2940450995689d0
     &*x1434 + x1436 - x1437*x590*x955 + x1438))
      x1486 = 31.0796769325781d0*x219
      x1487 = 1.29498653885742d0*x1413
      x1488 = x669**2
      x1489 = x678**2
      x1490 = x1489*x513
      x1491 = x1137*x669
      x1492 = 219.991371864982d0*x1413
      x1493 = x1173*x1477 + 1029.24699017327d0*x1384*x641 + 171.54116502
     &8879d0*x1454* x512 + 3087.74097051982d0*x1466 + 5279.79292475958d0
     &*x1473 - 5279.79292475958d0*x1475 - 3087.74097051982d0*x1476 + 257
     &.311747543319d0*x1485 - x1492*x216 - x1492*x634
      x1494 = x491*(x1039*x1414 + x1428 - x1479 - x1480 - x1481 + x1482)
     &
      x1495 = x1406*x658
      x1496 = x1181*x658**2
      x1497 = 4.75358019696859d0*x1413
      x1498 = x536*x994
      x1499 = 0.593716125431411d0*x1498 + x535
      x1500 = -x1499*x694 + x198 + 0.296858062715705d0*x700
      x1501 = 0.395810750287607d0*x553
      x1502 = x1499*x197
      x1503 = x1502*x653
      x1504 = x1502*x536
      x1505 = x10*x1502
      x1506 = x194*x199
      x1507 = x1073 + 1.18743225086282d0*x1498
      x1508 = x1499*x996
      x1509 = x1508*x536
      x1510 = 0.296858062715705d0*x1498 - x1509*x993 + x993
      x1511 = x1000*x697
      x1512 = 0.593716125431411d0*x994
      E = -rho_c*(-x201*x202 + x205)
      d1E(1) = rho_c*(x0*x207 - x0*x215 - x0*x229 + x0*x241 + x210*x213 
     &+ x210*x234 - x220*x233 + x221*x222 - x252 + x254*x256 + x440 + x4
     &46*x458 - x461*x572) + x16*x255 + x169 + x202*x457
      d1E(2) = -rho_c*(x575 + x576 - x578 + x580 + x620 + x624 - x699)
      d1E(3) = rho_c*x462*x702
      d2E(1) = -rho_c*(x1001*x986*(x191*x671*(-x484 + x487 + x511 - x532
     &) - x464*x989 + x487*x989 + x511*x989 + x998*(x484*x534 - x485*x55
     &1*x991 - x487* x534 - x510*x673 - x511*x534 + x531*x685 + x540*x67
     &1 + x683*x992 )) + x167*x720 + 0.499760307959092d0*x17*x182*x763 +
     & x17*x478* x710 + x17*x713 - x17*x724 - x17*x726 - x17*x730 - x17*
     &x732 + x17 *x744 + x17*x747 + x17*x753 + x17*x755 - x17*x757 - 1.0
     &5383986550676d0*x17*x762 + x17*x984 + x212*x720 - x216*x723* x739 
     &+ 8.0d0*x315*x714 + x441*x582*x623*(2.0d0*x442 - 2.0d0*x443 + x773
     &) - x698*(-x1002*x1035*x1071 + 4.23769752780008d0*x1002 + x1004*x1
     &010 + x1004*x1018 + 3.39015802224006d0*x1005*x472 - 7.068479476370
     &53d0*x1019 + x1020*x1023 - x1025*x189*x482 + x1026* x1027 - x1029*
     &x1030 - x1031*x1037*x534 + 2.54261851668005d0*x1033 + x1034*x1035*
     &x1036*x189 + x1069*x1070 + 3.02934834701594d0*x1072 *x570 - x1074*
     &x1075*x1076*x533 + x1080*x165*(x1035*x562 - x1037 + x1077 + 2.5426
     &1851668005d0*x559 + 1.08969365000573d0*x566) + x1122*(x1003*x1095*
     &x547 + x1010 + x1018 + x1020*x1085 - x1025* x481 + x1026*x1096 - x
     &1029*x1097 + x1069*x1121 + x1081 + x1082* x472 - 3.02934834701594d
     &0*x1084*x565 - 7.06847947637054d0*x1084 + 2.54261851668005d0*x1086
     &*x547 + x1088*x1091 - x1091*x1092 - x1120 *(-x1*x1099*x544 + x1098
     &*x508*x546 - x1100*x1101*x1102 + x1119)) - 2.54261851668005d0*x570
     &*x646 - x694*(gamma_cc*rho_c**( -8.66666666666667d0)*x1125*x1132*x
     &18*x654 + 4.03913112935459d0* x1003*x1124*x510 + x1009*x1125*x1126
     & + x1018*x534 + x1081*x534 - x1082*x1124 + 7.85024654737088d0*x112
     &7 + 0.296858062715705d0* x1129*(x1134 + x1135 + 5.33333333333333d0
     &*x442 - 5.33333333333333d0*x443) + 8.80283108639638d0*x1130*x510 -
     & 7.388467338702d0*x1130 - 3.02934834701594d0*x1136*x1171 - 7.06847
     &947637054d0*x1136 + x1138*x1139 + x1138*x1163 + 2.54261851668005d0
     &*x1140*x530 - 6.60212331479729d0*x1144*x510 + 5.5413505040265d0*x1
     &144 + x1147*x1149 + x1147*x1153 - x1150*x1151 - x1150*x1154 - 1.45
     &292486667431d0*x1156*x530 - x1157*x1159*x992 + x1161*x404 - 15.404
     &9544011937d0*x1162 + x1165*x404 + x1166* x1169 + x1169*x1170 - x11
     &74*x1175 - x1174*x1176 + x1183*(-x1059* x157 - x1060*x505 + x1068 
     &- 2.78d0*x1116 + 2.78d0*x1177 + x1178* x1179 + x1180*x1182 + 2.929
     &6748261088d0*x523 - 57.0429623636231d0 *x525 + 55.6d0*x75 + 2.78d0
     &*x918) + x1185*(89.4158322713032d0* x1061 + 28.5214811818116d0*x10
     &62 + 114.085924727246d0*x1063 - 0.0202109685722776d0*x1064 - 89.41
     &58322713032d0*x1065 - 164.623675016285d0*x1066 + 164.623675016285d
     &0*x1067 - 5.56d0* x1116 - 2.77866731225255d0*x1117 - 7.03121958266
     &113d0*x1118 + 5.56d0*x1177 + x1178*x1184 + x1179*x1182 - 0.0203701
     &531497855d0* x206 + 0.00336849476204627d0*x214 - 44.48d0*x316 + 5.
     &85934965221761d0*x523 - 114.085924727246d0*x525 + 0.49138872231528
     &1d0*x725 + 1437.98400373558d0*x731 - 0.149604492895997d0*x746 + 11
     &1.2d0*x75 - 4874.27236976031d0*x754 + 5.56d0*x868 - 0.012959684348
     &0306d0*x869 + 844.480404506768d0* x871 + 5.56d0*x872 - 111.2d0*x87
     &5 + 1316.98940013028d0*x879 + 5.56d0*x918))) - x711*x987 + x714*x7
     &16 - 236.868597145733d0*x714* x742 - x714*x761 - x718*x719 - x718*
     &x758 + 0.00233087847986162d0* x722*x723 - 151.884964839347d0*x723*
     &x749 + x739*x740 - x740*x741 - x764*x767) + x0*x703 - x0*x705 - x0
     &*x708 + x0*x709 + x167*x704 - 0.00777267271741373d0*x197*x4*x712 +
     & x212*x704 + x222*x707 - x233*x706 - 2.0d0*x252 + x256*x710 + 2.0d
     &0*x440 + x458*x711
      d2E(2) = rho_c*(x1044*x17 + x1186*x1187 + x1186*x1199 - 118.434298
     &572867d0*x1186* x742 - x1188*x1189 - x1188*x1195 + x1188*x1197 - x
     &1188*x1198 + x1190*x167 + x1190*x212 - x1194*x17 - x1196*x17 + 0.0
     &0777267271741373d0*x1200*(x444 + x773) - x1307*x17 + x1309* x446 -
     & x1311*x986 - x17*x254*x618 - x622*x987 + x698*(-rho_c**( -7.66666
     &666666667d0)*x468*x555*x643 - x1002*x1017*x1315 - 1.69507901112003
     &d0*x1002*x472 - 3.53423973818527d0*x1002*x649 - x1004*x1312 - x100
     &5*x1317 + x1019*x1316 + x1023*x1322 + x1027* x1320 - x1030*x1321 -
     & 0.726462433337157d0*x1032*x1036 - x1033* x650 + x1070*x1347 + x10
     &71*x1357 + x1075*x1359*x165 + x1122*( -x1003*x1317 + x1003*x1368*x
     &1369 + x1084*x1316 - 1.51467417350797d0*x1084*x1371 + x1085*x1322 
     &+ 1.27130925834002d0 *x1086*x1368 + x1088*x1370 - x1092*x1370 + x1
     &096*x1320 - x1097* x1321 - x1120*(x1102*x1374*x546 + 257.311747543
     &319d0*x1194 - 257.311747543319d0*x1262 + 257.311747543319d0*x1307 
     &+ 659.974115594947d0*x1339 + 5279.79292475958d0*x1340 + 162.699226
     &454617d0*x1341 + 5279.79292475958d0*x1372 - x1373*x520* x546 + x13
     &75 - x508*x677 + x544*(x545 + x773) - 4116.9879606931d0 *x630 - 0.
     &299881207482153d0*x722 + 19540.8928641839d0*x749) + x1121*x1347 - 
     &x1312 - x1315*x1366 - x1365 - 1.51467417350797d0* x1367*x565 - 3.5
     &3423973818527d0*x1367 - x1369*x466*x547) + 1.51467417350797d0*x131
     &9*x481 - 1.27130925834002d0*x1319 + x1348* x1350 + x1351*x1352 - x
     &1353*x1354 + x1355*x1356 - x1363*x1364 + x570*x651 - x694*(2.63873
     &833525072d0*gamma_cc*x1157*x1377 + 2.0195655646773d0*x1003*x1123*x
     &688 + x1126*x1155*x471 + 9.17695140756823d0*x1127*x1168*x669 + 4.4
     &0141554319819d0*x1130* x670 - 1.51467417350797d0*x1136*x1360 + 1.2
     &7130925834002d0*x1140* x679 - x1144*x1361 - x1151*x1390 - x1154*x1
     &390 - 0.726462433337157d0*x1156*x679 - x1157*x1380*x680 + x1161*x6
     &38 + x1165*x638 + x1168*x1385 + x1168*x1386 - 1.51467417350797d0*x
     &1171 *x1378 + x1183*(x1179*x1409 + x1180*x1406 + 2.78d0*x1325 + x1
     &346 - 28.5214811818116d0*x1372 + 2.78d0*x1393 - 2.78d0*x1395 - 1.7
     &5780489566528d0*x1396 - 0.00323992108700765d0*x1397 + 28.521481181
     &8116d0*x1399 + 211.120101126692d0*x1400 + 2.78d0* x1404 + 2.78d0*x
     &1405 + 7.13037029545289d0*x1407 + x1408*x219 - x1408*x631 - 44.48d
     &0*x636) + x1185*(-x1060*x666 + x1179*x1406 + x1184*x1409 - 658.494
     &700065139d0*x1267 + 5.56d0*x1325 + x1343 + 22.24d0*x1344 + 0.01010
     &54842861388d0*x1345 + 5.56d0*x1393 - 5.56d0*x1395 - 0.006479842174
     &0153d0*x1397 + 57.0429623636231d0* x1399 + 422.240202253384d0*x140
     &0 + 5.56d0*x1404 + 5.56d0*x1405 + 14.2607405909058d0*x1407 + x1410
     &*x219 - x1410*x631 + 88.96d0*x626 - 88.96d0*x636) - x1365*x534 + x
     &1366*x534*(-x1014 - x1015 - x1313 + x1314) + 4.40141554319819d0*x1
     &376*x510 - 3.694233669351d0*x1376 - 3.53423973818527d0*x1378 - 3.3
     &0106165739864d0*x1379*x510 + 2.77067525201325d0*x1379 - x1380*x138
     &3*x530 - 7.70247720059684d0* x1381 - 0.726462433337157d0*x1382*x53
     &0 + x1387*x509*x669 + x1388* x1389 + x1389*x1391 - x1401*x1402 - x
     &1401*x1403 + 2.0195655646773d0*x510*x534*x561 + 0.296858062715705d
     &0*x991*( -x1011 + x1012 - x1134 - x1135))) - x723*x767) - x575 - x
     &576 + x578 - x580 - x620 - x624 + x699
      d2E(3) = -rho_c*(-x1100*x17*x766 + x1188*x1414*x167 - x1189*x208*x
     &573 + 0.00259089090580458d0*x1200*x1453 + 0.0155453454348275d0*x13
     &09* x621 + x1311*x696 + x1411*x1412 + x1411*x1417 - x1415*x17 - x1
     &416 *x208 + x1452*x17 - x698*(0.544846825002868d0*rho_c**( -6.6666
     &6666666667d0)*x541*(x1008 + x1088*x1471 - x1092*x1471 - x1095*x136
     &8 - 3.02934834701594d0*x1371*x649 + x1455 + x1456 + x1457 - x1458 
     &+ x469*(0.242154144445719d0*x768 + 0.242154144445719d0*x771) + x47
     &6*(-1.51467417350797d0*x1452 + 1.51467417350797d0*x1459*x185 + 3.0
     &2934834701594d0*x1460 - 3.02934834701594d0*x1461 + x1462*x1486 - 1
     &8.1760900820957d0*x1463 + x1465*x1486 - 18.1760900820957d0*x1466 +
     & x1487*x216 - x1487* x629) - x564*(560.781954739666d0*x1373*x641 -
     & 93.4636591232777d0* x1454*x159 + 1682.345864219d0*x1466 + 2876.67
     &841172786d0*x1473 - x1474*x216 - x1474*x634 - 2876.67841172786d0*x
     &1475 - 1682.345864219d0*x1476 - 747.709272986221d0*x1477*x159 + 14
     &0.195488684916d0*x1485)) + x1002*x1008 + 0.242154144445719d0* x100
     &2*x1454 + x1002*x1455 + x1002*x1456 + x1002*x1457 - x1002* x1458 +
     & x1070*(x1049*x1459 + x1342*x1462 + x1342*x1465 - 2.78d0* x1452 + 
     &5.56d0*x1460 - 5.56d0*x1461 - 33.36d0*x1463 - x1464*x629 + x1467) 
     &+ 1.45292486667431d0*x1318*x555 + x1319*x1468 - 1.45292486667431d0
     &*x1351*x646 - x1354*x1468*x646 - 1.08969365000574d0*x1358*x1469 + 
     &x1363*x1470 - x694*(x1006*x1132* x1141 + x1008*x1034 + x1137*x1163
     &*x1488 - x1151*x1490 - x1154* x1490 - x1159*x1383*x679 + x1160*x13
     &78*x638 + x1164*x1381*x638 + x1183*(x1464*x634 + x1467 - 57.042962
     &3636231d0*x1473 + 57.0429623636231d0*x1475 + 33.36d0*x1476 - 2.78d
     &0*x1485 + 2.78d0* x1494 + 5.56d0*x1495 + 5.56d0*x1496) + x1185*(-6
     &6.72d0*x1466 - 114.085924727246d0*x1473 + 114.085924727246d0*x1475
     & + 66.72d0* x1476 - 5.56d0*x1485 + 5.56d0*x1494 + 11.12d0*x1495 + 
     &11.12d0* x1496 + x1497*x216 + x1497*x634) - 3.02934834701594d0*x13
     &60*x1378 + 8.80283108639638d0*x1376*x670 - 6.60212331479729d0*x137
     &9*x670 - 1.45292486667431d0*x1382*x679 + x1385*x1491 + x1386*x1491
     & + x1387 *x1488 + x1388*x1489 + x1391*x1489 - x1402*x1493 - x1403*
     &x1493 + 0.242154144445719d0*x1453*x486*x534 + 0.263873833525072d0*
     &x1453* x991 + 4.03913112935459d0*x561*x688)))
      d2E(4) = x462*(rho_c*(x1078*x1500*x712*x999 + x446*x701 - x460*(ga
     &mma_cc*rho_s* x1501*x186*x552 + x1002*x473 + x1031*x1074 - 0.29685
     &8062715705d0* x1075*x188 + x1349 - x1350*x1502 - x1352*x1503 + x13
     &53*x1504 - x1353 - x1356*x1505 + x1364*x1507 - x1506*x568 - x1506*
     &(x1077 + x567) + 0.692668813003313d0*x557 - 0.825265414349661d0*x5
     &58)) + x702)
      d2E(5) = -rho_c*x985*(x1308*x1510*x1511 + x1500*x622 - x459*(-x107
     &6*x1359*x1505 - x1357*x1503 - 0.296858062715705d0*x1469*x187 + x14
     &70*x1507 + x1501*x549 - x1504*x651 + 0.825265414349661d0*x557*x649
     & + 1.08969365000574d0*x644 - x694*(x1362 + x692) + x695))
      d2E(6) = -rho_c*x1511*(x1500*x1510*x985 + x515*(0.323484345894309d
     &0*rho_c**( -7.0d0)*x491**8.34d0*x988*x997 + x1073*x1508*x988 - x14
     &99*x1507* x536*x993/x995**2 + x1509*x1512 - x1512))
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, gamma_cc, tau_c, mu
      real*8, intent(out) :: E
      real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52,
     & x53, x54, x55
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = x0**0.666666666666667d0
      x6 = log(1.0d0 + 16.081984221457407d0/(2.2255694211506878d0*x1 + 5
     &.9825504357710795d0*x3 + 0.80042863499936356d0*x4 + 0.189700432574
     &75591d0*x5))
      x7 = x2*x6
      x8 = 0.0621813817393098d0*x7
      x9 = -x8
      x10 = 2.0d0*x7
      x11 = exp(x10) - 1.0d0
      x12 = gamma_cc**2
      x13 = rho_c**(-4.66666666666667d0)*x12
      x14 = gamma_cc*rho_c**(-2.33333333333333d0)
      x15 = 0.136211706250717d0*x14
      x16 = 1.0d0 + x15/x11
      x17 = 0.0310906908696549d0*log(x15*x16/(x16 + 0.018553628919731596
     &d0*x13/x11** 2) + 1.0d0)
      x18 = x17 + x9
      x19 = (0.160601656036401d0*x1 + 1.0d0)*log(1.0d0 + 32.163968442914
     &921d0/( 4.8440767160638538d0*x1 + 12.482198746797325d0*x3 + 2.3260
     &048119008192d0*x4 + 0.3819082618690966d0*x5))
      x20 = exp(1.999999999999994d0*x19) - 1.0d0
      x21 = 0.272423412501434d0*x14
      x22 = 1.0d0 + x21/x20
      x23 = -0.0310906908696548d0*x19 + 0.0155453454348274d0*log(x21*x22
     &/( 0.074214515678926424d0*x13/x20**2 + x22) + 1.0d0)
      x24 = exp(-0.46869933305458689d0*mu*rho_c**(-0.33333333333333331d0
     &))
      x25 = x12/tau_c**2
      x26 = x25/rho_c**2
      x27 = x23*x26
      x28 = x18*(0.00828125d0*x26 + 1.0d0)
      x29 = 2.8d0*x28
      x30 = rho_c**(-3)
      x31 = 0.001953125d0*gamma_cc**3*x30/tau_c**3
      x32 = -x17 + x8
      x33 = mu**2
      x34 = x0**1.66666666666667d0
      x35 = x0**1.0d0
      x36 = x0**1.33333333333333d0
      x37 = (0.0128480790170175d0*x1 - 0.00304892328706429d0*x35 + 0.000
     &275310793616358d0*x36 + 0.0315297404574782d0*x5 + 1.0d0)*exp (-0.4
     &6675853320810856d0*x1)
      x38 = x1*x33
      x39 = 4.60448529760779d0*mu*x3 + 1.0d0
      x40 = mu**3
      x41 = 0.5d0*x37 - 0.5d0
      x42 = 1d0/x5
      x43 = 17.2274040011292d0*x42
      x44 = -x43*(-0.0177069812999059d0*x1 + 1.0d0)/(0.337569861992026d0
     &*x1 + 0.0244354823084343d0*x5 + 1.0d0)
      x45 = (-0.000196557189785076d0*mu**8*x0**2.66666666666667d0*x7 + m
     &u**6*( -0.00331582639355101d0*x0**2.0d0*x7 - 0.00496533372342768d0
     &*x34* x41) - 0.00264117541143784d0*mu**5*x34*x37 + mu**4*( 0.00020
     &688890514282d0*x34*(-23.3866625383245d0*x42*( -0.240695990468967d0
     &*x1 + 0.26014827853397d0*x5)*exp( -0.33933171852197186d0*x1) + x43
     & + x44) - 0.0837628205355044d0* x35*x41 - 0.0209761558789675d0*x36
     &*x7) + x40*( 0.000132058770571892d0*x34*(x44 - 37.6991118430775d0*
     &( -3.07073492995203d0*x1 + 0.384834731559127d0*x5)*exp( -0.1923086
     &5217881402d0*x1)/x35) - 0.0445553741830535d0*x35*x37) - 0.06218138
     &17393098d0*log((4.6213167280845804d0*x38 + x39 + 1.914071024228974
     &9d0*x4*x40)/(2.1392848713714905d0*x38 + x39)))/( 0.237113969739022
     &d0*x33*x5 + 1.0d0)**4
      x46 = -x45 + x9
      x47 = -x46/(x2*x6)
      x48 = 0.0310906908696549d0*exp(x10 + 32.163968442914815d0*x45) - 0
     &.0310906908696549d0
      x49 = x14*x47**2.78d0
      x50 = 1.0d0 + 9.56024411508605d0*x49/x48
      x51 = x46 + 0.0310906908696549d0*log(307.49539002418987d0*x49*x50/
     &( 91.39826754003758d0*x13*x47**5.5599999999999996d0/x48**2 + x50) 
     &+ 1.0d0)
      x52 = 0.5d0*rho_c
      x53 = x25*x30*(x18*x52 + x23*x52)
      x54 = rho_c*(x24*(x32 + (x28 - 0.02390625d0*x53)*(x31*(x29 - 0.066
     &9375d0*x53) + 1.0d0)) + x51)
      x55 = x18*x26
      PBE_value = x18
      PBE_beta_value = x23
      PBE_alpha_value = x23
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*(x24*(x32 + (-0.02390625d0*x27 + x28)*(x31*(-0.0669375d0
     &*x27 + x29 ) + 1.0d0)) + x51)
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = x54
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = x54
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*(x24*(x32 + (x28 - 0.02390625d0*x55)*(x31*(x29 - 0.06693
     &75d0*x55) + 1.0d0)) + x51)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E,
     & d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x281, x282, x283
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = x0**0.666666666666667d0
      x6 = 2.22556942115069d0*x1 + 5.98255043577108d0*x3 + 0.80042863499
     &9364d0*x4 + 0.189700432574756d0*x5
      x7 = 1.0d0 + 16.0819842214574d0/x6
      x8 = log(x7)
      x9 = x2*x8
      x10 = 0.0621813817393098d0*x9
      x11 = -x10
      x12 = 2.0d0*x9
      x13 = exp(x12)
      x14 = x13 - 1.0d0
      x15 = 1d0/x14
      x16 = rho_c**(-2.33333333333333d0)
      x17 = 0.136211706250717d0*x16
      x18 = x15*x17
      x19 = gamma_cc*x18 + 1.0d0
      x20 = x14**(-2)
      x21 = gamma_cc**2
      x22 = rho_c**(-4.66666666666667d0)
      x23 = x21*x22
      x24 = x19 + 0.0185536289197316d0*x20*x23
      x25 = 1d0/x24
      x26 = x17*x25
      x27 = x19*x26
      x28 = gamma_cc*x27 + 1.0d0
      x29 = log(x28)
      x30 = 0.0310906908696549d0*x29
      x31 = x11 + x30
      x32 = 4.84407671606385d0*x1 + 12.4821987467973d0*x3 + 2.3260048119
     &0082d0*x4 + 0.381908261869097d0*x5
      x33 = 1.0d0 + 32.1639684429149d0/x32
      x34 = log(x33)
      x35 = x34*(0.160601656036401d0*x1 + 1.0d0)
      x36 = exp(1.999999999999994d0*x35)
      x37 = x36 - 1.0d0
      x38 = 1d0/x37
      x39 = 0.272423412501434d0*x16
      x40 = x38*x39
      x41 = gamma_cc*x40 + 1.0d0
      x42 = x37**(-2)
      x43 = 0.0742145156789264d0*x23*x42 + x41
      x44 = 1d0/x43
      x45 = x39*x44
      x46 = x41*x45
      x47 = gamma_cc*x46 + 1.0d0
      x48 = log(x47)
      x49 = -0.0310906908696548d0*x35 + 0.0155453454348274d0*x48
      x50 = exp(-0.46869933305458689d0*mu*rho_c**(-0.33333333333333331d0
     &))
      x51 = rho_c**(-2)
      x52 = tau_c**(-2)
      x53 = x21*x52
      x54 = x51*x53
      x55 = x49*x54
      x56 = 0.00828125d0*x54 + 1.0d0
      x57 = x31*x56
      x58 = -0.02390625d0*x55 + x57
      x59 = 2.8d0*x57
      x60 = -0.0669375d0*x55 + x59
      x61 = gamma_cc**3
      x62 = tau_c**(-3)
      x63 = x61*x62
      x64 = x60*x63
      x65 = rho_c**(-3)
      x66 = 0.001953125d0*x65
      x67 = x64*x66 + 1.0d0
      x68 = x10 - x30
      x69 = x50*(x58*x67 + x68)
      x70 = mu**2
      x71 = x5*x70
      x72 = 0.237113969739022d0*x71 + 1.0d0
      x73 = x72**(-4)
      x74 = mu**5
      x75 = exp(-0.46675853320810856d0*x1)
      x76 = x0**1.0d0
      x77 = 0.00304892328706429d0*x76
      x78 = x0**1.33333333333333d0
      x79 = x75*(0.0128480790170175d0*x1 + 0.0315297404574782d0*x5 - x77
     & + 0.000275310793616358d0*x78 + 1.0d0)
      x80 = x74*x79
      x81 = x0**1.66666666666667d0
      x82 = 0.00264117541143784d0*x81
      x83 = mu**8
      x84 = x0**2.66666666666667d0*x83
      x85 = x84*x9
      x86 = -0.0621813817393098d0
      x87 = x1*x70
      x88 = mu*x3
      x89 = 4.60448529760779d0*x88 + 1.0d0
      x90 = 2.13928487137149d0*x87 + x89
      x91 = 1d0/x90
      x92 = mu**3
      x93 = x4*x92
      x94 = 4.62131672808458d0*x87 + x89 + 1.91407102422897d0*x93
      x95 = mu**6
      x96 = 0.5d0*x79 - 0.5d0
      x97 = 0.00496533372342768d0*x81
      x98 = x0**2.0d0
      x99 = x9*x98
      x100 = 0.0445553741830535d0*x76
      x101 = x100*x79
      x102 = exp(-0.19230865217881402d0*x1)
      x103 = x102*(-3.07073492995203d0*x1 + 0.384834731559127d0*x5)
      x104 = 37.6991118430775d0/x76
      x105 = x103*x104
      x106 = 1d0/x5
      x107 = 17.2274040011292d0*x106
      x108 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x5 + 1.0d0
      x109 = 1d0/x108
      x110 = -0.0177069812999059d0*x1 + 1.0d0
      x111 = x109*x110
      x112 = -x107*x111
      x113 = -x105 + x112
      x114 = 0.000132058770571892d0*x81
      x115 = mu**4
      x116 = 0.0837628205355044d0*x76
      x117 = x116*x96
      x118 = exp(-0.33933171852197186d0*x1)
      x119 = x118*(-0.240695990468967d0*x1 + 0.26014827853397d0*x5)
      x120 = 23.3866625383245d0*x106
      x121 = x107 + x112 - x119*x120
      x122 = 0.00020688890514282d0*x81
      x123 = x115*(-x117 + x121*x122 - 0.0209761558789675d0*x78*x9) - x8
     &0*x82 - 0.000196557189785076d0*x85 + x86*log(x91*x94) + x92*(-x101
     & + x113 *x114) + x95*(-x96*x97 - 0.00331582639355101d0*x99)
      x124 = x123*x73
      x125 = x11 - x124
      x126 = 1d0/x2
      x127 = 1d0/x8
      x128 = x126*x127
      x129 = -x125*x128
      x130 = x129**5.56d0
      x131 = exp(x12 + 32.163968442914815d0*x124)
      x132 = 0.0310906908696549d0*x131 - 0.0310906908696549d0
      x133 = x132**(-2)
      x134 = x130*x133
      x135 = 91.3982675400376d0*x134*x23
      x136 = x129**2.78d0
      x137 = 1d0/x132
      x138 = x136*x137
      x139 = 9.56024411508605d0*x138*x16
      x140 = gamma_cc*x139
      x141 = x140 + 1.0d0
      x142 = x135 + x141
      x143 = 1d0/x142
      x144 = x141*x143
      x145 = x136*x16
      x146 = 307.49539002419d0*x145
      x147 = x144*x146
      x148 = gamma_cc*x147
      x149 = x148 + 1.0d0
      x150 = x125 + 0.0310906908696549d0*log(x149)
      x151 = x150 + x69
      x152 = x53*x65
      x153 = 0.0478125d0*x49
      x154 = x0*x1
      x155 = x154*x34
      x156 = x0*x3
      x157 = x0*x4
      x158 = x0*x5
      x159 = 32.1639684429149d0*(1.61469223868795d0*x154 + 2.08036645779
     &955d0*x156 + 1.16300240595041d0*x157 + 0.254605507912731d0*x158)/(
     &x32**2*x33)
      x160 = 1d0/x47
      x161 = gamma_cc*rho_c**(-3.33333333333333d0)
      x162 = 0.635654629170013d0*x161
      x163 = x162*x38
      x164 = x36*(-0.107067770690933d0*x155 + x159*(0.321203312072801d0*
     &x1 + 1.99999999999999d0))
      x165 = gamma_cc*x39
      x166 = x164*x165*x42
      x167 = rho_c**(-5.66666666666667d0)*x21
      x168 = x165*x41/x43**2
      x169 = 0.00166440548032745d0*x155 + x159*(-0.00499321644098237d0*x
     &1 - 0.0310906908696548d0) + 0.0155453454348274d0*x160*(gamma_cc*x4
     &5*( -x163 - x166) - x162*x41*x44 + x168*(x163 + 0.148429031357853d
     &0* x164*x23/x37**3 + x166 + 0.346334406501657d0*x167*x42))
      x170 = x169*x54
      x171 = 0.317827314585006d0*x161
      x172 = x15*x171
      x173 = x154*x8
      x174 = (0.741856473716896d0*x154 + 0.99709173929518d0*x156 + 0.400
     &214317499682d0*x157 + 0.126466955049837d0*x158)/(x6**2*x7)
      x175 = 16.0819842214574d0*x174
      x176 = -0.0883792666034679d0*x173 + x175*(0.265137799810404d0*x1 +
     & 2.0d0)
      x177 = x13*x176
      x178 = gamma_cc*x17
      x179 = x177*x178*x20
      x180 = x178*x19/x24**2
      x181 = 1d0/x28
      x182 = 0.0310906908696549d0*x181
      x183 = x182*(gamma_cc*x26*(-x172 - x179) - x171*x19*x25 + x180*( 0
     &.0865836016254141d0*x167*x20 + x172 + x179 + 0.0371072578394632d0*
     &x177*x23/x14**3))
      x184 = 0.00274777245725523d0*x173
      x185 = 0.0082433173717657d0*x1
      x186 = x175*(-x185 - 0.0621813817393098d0) + x184
      x187 = x183 + x186
      x188 = x187*x56
      x189 = x152*x31
      x190 = x188 - 0.0165625d0*x189
      x191 = rho_c**(-4)
      x192 = 0.005859375d0*x191
      x193 = 0.133875d0*x49
      x194 = 0.0231875d0*x54 + 2.8d0
      x195 = x187*x194
      x196 = -0.046375d0*x189 + x195
      x197 = x63*x66
      x198 = x175*(x185 + 0.0621813817393098d0) - x183 - x184
      x199 = 0.156233111018196d0*mu*rho_c**(-1.33333333333333d0)
      x200 = x125*x127*x154/x2**2
      x201 = x125*x126*x174/x8**2
      x202 = x0*x8
      x203 = x0*x81
      x204 = x0*x78
      x205 = x75*(x0*x77 - 0.00428269300567249d0*x154 - 0.02101982697165
     &21d0*x158 - 0.000367081058155144d0*x204)
      x206 = x174*x2
      x207 = 0.767414216267966d0*x0*x88
      x208 = x0*x87
      x209 = 0.00152446164353214d0*x76
      x210 = 0.155586177736036d0*x154*x75*(0.00642403950850874d0*x1 - x2
     &09 + 0.0157648702287391d0*x5 + 0.000137655396808179d0*x78 + 0.5d0)
     & + x75*(x0*x209 - 0.00214134650283625d0*x154 - 0.0105099134858261d
     &0* x158 - 0.000183540529077572d0*x204)
      x211 = x0/x1
      x212 = x0*x106
      x213 = 11.4849360007528d0*x212
      x214 = -x107*x110*(0.112523287330675d0*x154 + 0.0162903215389562d0
     &*x158)/x108** 2 - 0.101681773497972d0*x109*x211 - x111*x213
      x215 = x73*(0.000524152506093537d0*x0*x85 + 8.68579013942185d-6*x0
     &**3.0d0*x202* x83 - 0.000410930386996015d0*x0**3.0d0*x80 + x115*(x
     &0*x117 - x116 *x210 - 0.0003448148419047d0*x121*x203 + x122*(-x118
     &*x120*( 0.0802319968229891d0*x154 - 0.173432185689313d0*x158) - 2.
     &64527879654103d0*x119*x211 - 15.591108358883d0*x119*x212 + x213 + 
     &x214) + 0.000926928636371586d0*x203*x8 + 0.0279682078386234d0* x20
     &4*x9 - 0.337338207872387d0*x206*x78) + 0.00440195901906306d0* x203
     &*x80 - x205*x74*x82 - 0.00316102962473761d0*x206*x84 + x86* x90*(x
     &91*(-0.957035512114487d0*x0*x93 - x207 - 1.54043890936153d0 *x208)
     & + x94*(x207 + 0.71309495712383d0*x208)/x90**2)/x94 + x92*( x0*x10
     &1 - x100*x205 - 0.000220097950953153d0*x113*x203 + x114*( -2.41662
     &179562687d0*x0**0.333333333333333d0*x103 - x0*x105 - x102 *x104*(1
     &.02357830998401d0*x154 - 0.256556487706084d0*x158) + x214 ) - 0.00
     &693220036674015d0*x204*x79) + x95*(0.00663165278710203d0* x0*x99 +
     & 0.00014652515242323d0*x0**2.33333333333333d0*x202 + 0.00827555620
     &57128d0*x203*x96 - 0.0533250677421794d0*x206*x98 - x210*x97))
      x216 = x0*x123*x71/x72**5
      x217 = x186 - x215 - 0.632303919304059d0*x216
      x218 = x128*x217
      x219 = x9/x125
      x220 = x219*(-0.12284718057882d0*x200 + 44.7079161356516d0*x201 - 
     &2.78d0*x218)
      x221 = 22.3072362685341d0*x138*x161
      x222 = x131*(x176 + 32.1639684429148d0*x215 + 20.3374033068271d0*x
     &216)
      x223 = 0.297234594420578d0*gamma_cc*x133*x145*x222
      x224 = x140*x220
      x225 = gamma_cc*x146
      x226 = x141*x225/x142**2
      x227 = 0.0310906908696549d0/x149
      x228 = x217 + x227*(-717.489243389776d0*x136*x144*x161 + x143*x225
     &*(-x221 - x223 - x224) - x148*x220 + x226*(5.68327056421864d0*x130
     &*x222*x23 /x132**3 + 426.525248520175d0*x134*x167 + x135*x219*( -0
     &.245694361157641d0*x200 + 89.4158322713032d0*x201 - 5.56d0*x218 ) 
     &+ x221 + x223 + x224))
      x229 = gamma_cc*x22
      x230 = 0.0185536289197316d0*x15*x229*x25 + x180*(-x18 - 0.03710725
     &78394632d0* x20*x229) + x27
      x231 = x182*x230
      x232 = -x231
      x233 = gamma_cc*x52
      x234 = x233*x51
      x235 = x160*(x168*(-0.148429031357853d0*x229*x42 - x40) + 0.074214
     &5156789264d0* x229*x38*x44 + x46)
      x236 = x235*x54
      x237 = x231*x56
      x238 = 0.0165625d0*x51
      x239 = x238*x31
      x240 = x233*x239 + x237
      x241 = x21*x62
      x242 = 0.005859375d0*x65
      x243 = x242*x60
      x244 = x194*x231
      x245 = x234*x31
      x246 = x244 + 0.046375d0*x245
      x247 = x227*(2939.73099289485d0*x130*x137*x143*x229 + x147 + x226*
     &( -182.796535080075d0*x134*x229 - x139))
      x248 = -x239*x241
      x249 = x241*x51
      x250 = -x238*x241*(0.0870539344350337d0*x29 - 0.174107868870067d0*
     &x9)
      x251 = x61/tau_c**4
      x252 = rho_c*x50
      x253 = 0.5d0*rho_c
      x254 = x253*x31 + x253*x49
      x255 = x152*x254
      x256 = -0.02390625d0*x255 + x57
      x257 = -0.0669375d0*x255 + x59
      x258 = x197*x257 + 1.0d0
      x259 = x50*(x256*x258 + x68)
      x260 = x150 + x259
      x261 = rho_c*x260
      x262 = x191*x254*x53
      x263 = x152*(x169*x253 + x187*x253 + 0.0155453454348275d0*x29 - 0.
     &0155453454348274d0*x35 + 0.00777267271741372d0*x48 - 0.03109069086
     &96549d0*x9)
      x264 = x192*x63
      x265 = rho_c*(x199*x259 + x228 + x50*(x198 + x256*(x197*(x196 + 0.
     &2008125d0* x262 - 0.0669375d0*x263) - x257*x264) + x258*(x190 + 0.
     &07171875d0 *x262 - 0.02390625d0*x263))) + x260
      x266 = x181*x230
      x267 = x152*(0.00777267271741372d0*rho_c*x235 + 0.0155453454348275
     &d0*rho_c*x266 )
      x268 = x254*x65
      x269 = x233*x268
      x270 = x242*x257
      x271 = rho_c*(x247 + x50*(x232 + x256*(x197*(x246 - 0.0669375d0*x2
     &67 - 0.133875d0*x269) + x241*x270) + x258*(x240 - 0.02390625d0*x26
     &7 - 0.0478125d0*x269)))
      x272 = x241*x268
      x273 = x252*(x256*(x197*(x250 + 0.133875d0*x272) - x251*x270) + x2
     &58*(x248 + 0.0478125d0*x272))
      x274 = x31*x54
      x275 = -0.02390625d0*x274 + x57
      x276 = -0.0669375d0*x274 + x59
      x277 = x197*x276 + 1.0d0
      x278 = x50*(x275*x277 + x68)
      x279 = x150 + x278
      x280 = x187*x54
      x281 = x266*x54
      x282 = x242*x276
      x283 = x249*x31
      PBE_value = x31
      PBE_beta_value = x49
      PBE_alpha_value = x49
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*x151
      d1E(1) = rho_c*(x199*x69 + x228 + x50*(x198 + x58*(-x192*x64 + x19
     &7*(x152*x193 - 0.0669375d0*x170 + x196)) + x67*(x152*x153 - 0.0239
     &0625d0*x170 + x190))) + x151
      d1E(3) = rho_c*(x247 + x50*(x232 + x58*(x197*(-x193*x234 - 0.00104
     &056656004376d0* x236 + x246) + x241*x243) + x67*(-x153*x234 - 0.00
     &0371630914301344d0*x236 + x240)))
      d1E(6) = x252*(x58*(x197*(x193*x249 + x250) - x243*x251) + x67*(x1
     &53*x249 + x248 ))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = x261
      d1E(1) = x265
      d1E(3) = x271
      d1E(6) = x273
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = x261
      d1E(1) = x265
      d1E(3) = x271
      d1E(6) = x273
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*x279
      d1E(1) = rho_c*(x199*x278 + x228 + x50*(x198 + x275*(x197*(0.0875d
     &0*x189 + x195 - 0.0669375d0*x280) - x264*x276) + x277*(x188 + 0.03
     &125d0*x189 - 0.02390625d0*x280))) + x279
      d1E(3) = rho_c*(x247 + x50*(x232 + x275*(x197*(x244 - 0.0875d0*x24
     &5 - 0.00208113312008752d0*x281) + x241*x282) + x277*(x237 - 0.0312
     &5d0 *x245 - 0.000743261828602687d0*x281)))
      d1E(6) = x252*(x275*(x197*(x250 + 0.133875d0*x283) - x251*x282) + 
     &0.03125d0*x277* x283)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_TPSS_GWS_ERF(rho_c, gamma_cc, tau_c, mu, E,
     & d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x809, x810, x811, x812, x813, x814, x815
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = 0.132568899905202d0*x1 + 1.0d0
      x3 = x0**0.166666666666667d0
      x4 = sqrt(x0)
      x5 = 0.800428634999364d0*x4
      x6 = x0**0.666666666666667d0
      x7 = 2.22556942115069d0*x1 + 5.98255043577108d0*x3 + x5 + 0.189700
     &432574756d0 *x6
      x8 = 1.0d0 + 16.0819842214574d0/x7
      x9 = log(x8)
      x10 = x2*x9
      x11 = 0.0621813817393098d0*x10
      x12 = -x11
      x13 = 2.0d0*x10
      x14 = exp(x13)
      x15 = x14 - 1.0d0
      x16 = x15**(-2)
      x17 = gamma_cc**2
      x18 = rho_c**(-4.66666666666667d0)
      x19 = x17*x18
      x20 = rho_c**(-2.33333333333333d0)
      x21 = 1d0/x15
      x22 = x20*x21
      x23 = 0.136211706250717d0*x22
      x24 = gamma_cc*x23 + 1.0d0
      x25 = 0.0185536289197316d0*x16*x19 + x24
      x26 = 1d0/x25
      x27 = x24*x26
      x28 = 0.136211706250717d0*x20
      x29 = x27*x28
      x30 = gamma_cc*x29 + 1.0d0
      x31 = log(x30)
      x32 = 0.0310906908696549d0*x31
      x33 = x12 + x32
      x34 = 2.32600481190082d0*x4
      x35 = 4.84407671606385d0*x1 + 12.4821987467973d0*x3 + x34 + 0.3819
     &08261869097d0*x6
      x36 = 1.0d0 + 32.1639684429149d0/x35
      x37 = log(x36)
      x38 = x37*(0.160601656036401d0*x1 + 1.0d0)
      x39 = exp(1.999999999999994d0*x38)
      x40 = x39 - 1.0d0
      x41 = x40**(-2)
      x42 = 1d0/x40
      x43 = x20*x42
      x44 = 0.272423412501434d0*x43
      x45 = gamma_cc*x44 + 1.0d0
      x46 = 0.0742145156789264d0*x19*x41 + x45
      x47 = 1d0/x46
      x48 = x45*x47
      x49 = 0.272423412501434d0*x20
      x50 = x48*x49
      x51 = gamma_cc*x50 + 1.0d0
      x52 = log(x51)
      x53 = -0.0310906908696548d0*x38 + 0.0155453454348274d0*x52
      x54 = exp(-0.46869933305458689d0*mu*rho_c**(-0.33333333333333331d0
     &))
      x55 = rho_c**(-2)
      x56 = tau_c**(-2)
      x57 = x55*x56
      x58 = x17*x57
      x59 = x53*x58
      x60 = 0.00828125d0*x58 + 1.0d0
      x61 = x33*x60
      x62 = -0.02390625d0*x59 + x61
      x63 = rho_c**(-3)
      x64 = 2.8d0*x61
      x65 = -0.0669375d0*x59 + x64
      x66 = x63*x65
      x67 = gamma_cc**3
      x68 = tau_c**(-3)
      x69 = x67*x68
      x70 = 0.001953125d0*x69
      x71 = x66*x70 + 1.0d0
      x72 = x11 - x32
      x73 = x54*(x62*x71 + x72)
      x74 = mu**2
      x75 = x6*x74
      x76 = 0.237113969739022d0*x75 + 1.0d0
      x77 = x76**(-4)
      x78 = mu**5
      x79 = exp(-0.46675853320810856d0*x1)
      x80 = x0**1.0d0
      x81 = 0.00304892328706429d0*x80
      x82 = x0**1.33333333333333d0
      x83 = x79*(0.0128480790170175d0*x1 + 0.0315297404574782d0*x6 - x81
     & + 0.000275310793616358d0*x82 + 1.0d0)
      x84 = x78*x83
      x85 = x0**1.66666666666667d0
      x86 = 0.00264117541143784d0*x85
      x87 = mu**8
      x88 = x0**2.66666666666667d0*x87
      x89 = x10*x88
      x90 = -0.0621813817393098d0
      x91 = x1*x74
      x92 = mu*x3
      x93 = 4.60448529760779d0*x92 + 1.0d0
      x94 = 2.13928487137149d0*x91 + x93
      x95 = 1d0/x94
      x96 = mu**3
      x97 = x4*x96
      x98 = 4.62131672808458d0*x91 + x93 + 1.91407102422897d0*x97
      x99 = mu**6
      x100 = 0.5d0*x83 - 0.5d0
      x101 = 0.00496533372342768d0*x85
      x102 = x0**2.0d0
      x103 = x10*x102
      x104 = 0.0445553741830535d0*x80
      x105 = x104*x83
      x106 = exp(-0.19230865217881402d0*x1)
      x107 = x106*(-3.07073492995203d0*x1 + 0.384834731559127d0*x6)
      x108 = 1d0/x80
      x109 = 37.6991118430775d0*x108
      x110 = x107*x109
      x111 = 1d0/x6
      x112 = 17.2274040011292d0*x111
      x113 = 0.337569861992026d0*x1 + 0.0244354823084343d0*x6 + 1.0d0
      x114 = 1d0/x113
      x115 = -0.0177069812999059d0*x1 + 1.0d0
      x116 = x114*x115
      x117 = -x112*x116
      x118 = -x110 + x117
      x119 = 0.000132058770571892d0*x85
      x120 = mu**4
      x121 = 0.0837628205355044d0*x80
      x122 = x100*x121
      x123 = exp(-0.33933171852197186d0*x1)
      x124 = x123*(-0.240695990468967d0*x1 + 0.26014827853397d0*x6)
      x125 = 23.3866625383245d0*x111
      x126 = x112 + x117 - x124*x125
      x127 = 0.00020688890514282d0*x85
      x128 = x120*(-0.0209761558789675d0*x10*x82 - x122 + x126*x127) - x
     &84*x86 - 0.000196557189785076d0*x89 + x90*log(x95*x98) + x96*(-x10
     &5 + x118 *x119) + x99*(-x100*x101 - 0.00331582639355101d0*x103)
      x129 = x128*x77
      x130 = x12 - x129
      x131 = 1d0/x2
      x132 = 1d0/x9
      x133 = x131*x132
      x134 = -x130*x133
      x135 = x134**2.78d0
      x136 = exp(32.163968442914815d0*x129 + x13)
      x137 = 0.0310906908696549d0*x136 - 0.0310906908696549d0
      x138 = 1d0/x137
      x139 = x135*x138
      x140 = x139*x20
      x141 = 9.56024411508605d0*x140
      x142 = gamma_cc*x141
      x143 = x142 + 1.0d0
      x144 = x134**5.56d0
      x145 = x137**(-2)
      x146 = x144*x145
      x147 = x146*x19
      x148 = 91.3982675400376d0*x147
      x149 = x143 + x148
      x150 = 1d0/x149
      x151 = x135*x150
      x152 = x143*x151
      x153 = 307.49539002419d0*x20
      x154 = x152*x153
      x155 = gamma_cc*x154
      x156 = x155 + 1.0d0
      x157 = x130 + 0.0310906908696549d0*log(x156)
      x158 = x157 + x73
      x159 = x56*x63
      x160 = x159*x17
      x161 = x160*x53
      x162 = -0.00499321644098237d0*x1 - 0.0310906908696548d0
      x163 = x0*x3
      x164 = x0*x1
      x165 = x0*x4
      x166 = x0*x6
      x167 = 2.08036645779955d0*x163 + 1.61469223868795d0*x164 + 1.16300
     &240595041d0* x165 + 0.254605507912731d0*x166
      x168 = 1d0/x36
      x169 = x35**(-2)
      x170 = 32.1639684429149d0*x168*x169
      x171 = x167*x170
      x172 = rho_c**(-3.33333333333333d0)
      x173 = 0.635654629170013d0*x172
      x174 = x173*x48
      x175 = gamma_cc*x174
      x176 = x172*x42
      x177 = 0.635654629170013d0*x176
      x178 = gamma_cc*x177
      x179 = x164*x37
      x180 = 0.321203312072801d0*x1 + 1.99999999999999d0
      x181 = x171*x180 - 0.107067770690933d0*x179
      x182 = x181*x39
      x183 = x182*x41
      x184 = x183*x49
      x185 = gamma_cc*x184
      x186 = -x178 - x185
      x187 = x47*x49
      x188 = x186*x187
      x189 = gamma_cc*x188
      x190 = rho_c**(-5.66666666666667d0)
      x191 = x17*x190
      x192 = x191*x41
      x193 = x40**(-3)
      x194 = x182*x193
      x195 = 0.148429031357853d0*x19
      x196 = x178 + x185 + 0.346334406501657d0*x192 + x194*x195
      x197 = x46**(-2)
      x198 = x197*x45
      x199 = x198*x49
      x200 = x196*x199
      x201 = gamma_cc*x200
      x202 = -x175 + x189 + x201
      x203 = 1d0/x51
      x204 = 0.0155453454348274d0*x203
      x205 = 0.00166440548032745d0*x179 + x202*x204
      x206 = x162*x171 + x205
      x207 = x206*x58
      x208 = 0.317827314585006d0*x172
      x209 = x208*x27
      x210 = gamma_cc*x209
      x211 = x172*x21
      x212 = 0.317827314585006d0*x211
      x213 = gamma_cc*x212
      x214 = x164*x9
      x215 = 0.265137799810404d0*x1 + 2.0d0
      x216 = 0.99709173929518d0*x163 + 0.741856473716896d0*x164 + 0.4002
     &14317499682d0 *x165 + 0.126466955049837d0*x166
      x217 = 1d0/x8
      x218 = x217/x7**2
      x219 = 16.0819842214574d0*x218
      x220 = x216*x219
      x221 = -0.0883792666034679d0*x214 + x215*x220
      x222 = x14*x221
      x223 = x16*x222
      x224 = x223*x28
      x225 = gamma_cc*x224
      x226 = -x213 - x225
      x227 = x26*x28
      x228 = x226*x227
      x229 = gamma_cc*x228
      x230 = x16*x191
      x231 = x15**(-3)
      x232 = x222*x231
      x233 = 0.0371072578394632d0*x19
      x234 = x213 + x225 + 0.0865836016254141d0*x230 + x232*x233
      x235 = x25**(-2)
      x236 = x235*x24
      x237 = x236*x28
      x238 = x234*x237
      x239 = gamma_cc*x238
      x240 = -x210 + x229 + x239
      x241 = 1d0/x30
      x242 = 0.0310906908696549d0*x241
      x243 = x240*x242
      x244 = 0.00274777245725523d0*x214
      x245 = 0.0082433173717657d0*x1
      x246 = -x245 - 0.0621813817393098d0
      x247 = x220*x246
      x248 = x244 + x247
      x249 = x243 + x248
      x250 = x249*x60
      x251 = x160*x33
      x252 = x250 - 0.0165625d0*x251
      x253 = 0.0478125d0*x161 - 0.02390625d0*x207 + x252
      x254 = rho_c**(-4)
      x255 = x254*x69
      x256 = 0.005859375d0*x255
      x257 = 0.0231875d0*x58 + 2.8d0
      x258 = x249*x257
      x259 = -0.046375d0*x251 + x258
      x260 = 0.133875d0*x161 - 0.0669375d0*x207 + x259
      x261 = x63*x70
      x262 = -x256*x65 + x260*x261
      x263 = -x244
      x264 = x245 + 0.0621813817393098d0
      x265 = x220*x264 - x243 + x263
      x266 = x54*(x253*x71 + x262*x62 + x265)
      x267 = mu*rho_c**(-1.33333333333333d0)
      x268 = x267*x73
      x269 = 717.489243389776d0*x172
      x270 = x152*x269
      x271 = gamma_cc*x270
      x272 = x2**2
      x273 = 1d0/x272
      x274 = x130*x132
      x275 = x273*x274
      x276 = x164*x275
      x277 = x130*x131
      x278 = x9**2
      x279 = 1d0/x278
      x280 = x216*x218
      x281 = x279*x280
      x282 = x277*x281
      x283 = x0**3.0d0
      x284 = x283*x87*x9
      x285 = x0*x85
      x286 = x0**2.0d0
      x287 = x286*x84
      x288 = x0*x82
      x289 = x0*x81 - 0.00428269300567249d0*x164 - 0.0210198269716521d0*
     &x166 - 0.000367081058155144d0*x288
      x290 = x289*x79
      x291 = x78*x86
      x292 = x2*x280
      x293 = 0.00316102962473761d0*x88
      x294 = 1d0/x98
      x295 = x0*x92
      x296 = 0.767414216267966d0*x295
      x297 = -x296
      x298 = x0*x91
      x299 = 1.54043890936153d0*x298
      x300 = 0.957035512114487d0*x0*x97
      x301 = x297 - x299 - x300
      x302 = 0.71309495712383d0*x298
      x303 = x296 + x302
      x304 = x94**(-2)
      x305 = x304*x98
      x306 = x301*x95 + x303*x305
      x307 = x294*x306*x90
      x308 = x0**2.33333333333333d0
      x309 = x308*x9
      x310 = 0.00152446164353214d0*x80
      x311 = x79*(0.00642403950850874d0*x1 - x310 + 0.0157648702287391d0
     &*x6 + 0.000137655396808179d0*x82 + 0.5d0)
      x312 = x79*(x0*x310 - 0.00214134650283625d0*x164 - 0.0105099134858
     &261d0*x166 - 0.000183540529077572d0*x288)
      x313 = 0.155586177736036d0*x164*x311 + x312
      x314 = 0.0533250677421794d0*x102
      x315 = x104*x79
      x316 = 1.02357830998401d0*x164 - 0.256556487706084d0*x166
      x317 = x106*x109
      x318 = x0**(-0.666666666666667d0)
      x319 = x0*x318
      x320 = 1d0/x1
      x321 = x114*x320
      x322 = x0*x111
      x323 = 11.4849360007528d0*x322
      x324 = x113**(-2)
      x325 = 0.112523287330675d0*x164 + 0.0162903215389562d0*x166
      x326 = x324*x325
      x327 = x112*x115
      x328 = -0.101681773497972d0*x0*x321 - x116*x323 - x326*x327
      x329 = -x0*x110 - 2.41662179562687d0*x107*x319 - x316*x317 + x328
      x330 = 0.337338207872387d0*x82
      x331 = 0.0802319968229891d0*x164 - 0.173432185689313d0*x166
      x332 = x123*x125
      x333 = x0*x320
      x334 = -15.591108358883d0*x124*x322 - 2.64527879654103d0*x124*x333
     & + x323 + x328 - x331*x332
      x335 = 8.68579013942185d-6*x0*x284 - 0.000410930386996015d0*x0*x28
     &7 + 0.000524152506093537d0*x0*x89 + x120*(x0*x122 + 0.027968207838
     &6234d0*x10*x288 - x121*x313 - 0.0003448148419047d0* x126*x285 + x1
     &27*x334 + 0.000926928636371586d0*x285*x9 - x292* x330) + 0.0044019
     &5901906306d0*x285*x84 - x290*x291 - x292*x293 + x307*x94 + x96*(x0
     &*x105 - 0.000220097950953153d0*x118*x285 + x119 *x329 - 0.00693220
     &036674015d0*x288*x83 - x289*x315) + x99*( 0.00663165278710203d0*x0
     &*x103 + 0.00014652515242323d0*x0*x309 + 0.0082755562057128d0*x100*
     &x285 - x101*x313 - x292*x314)
      x336 = x335*x77
      x337 = x76**(-5)
      x338 = x0*x337*x75
      x339 = x128*x338
      x340 = 0.632303919304059d0*x339
      x341 = x248 - x336 - x340
      x342 = x133*x341
      x343 = -0.12284718057882d0*x276 + 44.7079161356516d0*x282 - 2.78d0
     &*x342
      x344 = 1d0/x130
      x345 = x10*x344
      x346 = x343*x345
      x347 = x155*x346
      x348 = x139*x172
      x349 = 22.3072362685341d0*x348
      x350 = gamma_cc*x349
      x351 = x221 + 32.1639684429148d0*x336 + 20.3374033068271d0*x339
      x352 = x136*x351
      x353 = x145*x352
      x354 = x135*x20
      x355 = 0.297234594420578d0*x354
      x356 = x353*x355
      x357 = gamma_cc*x356
      x358 = x142*x346
      x359 = -x350 - x357 - x358
      x360 = x151*x153
      x361 = x359*x360
      x362 = gamma_cc*x361
      x363 = x135*x153
      x364 = x146*x191
      x365 = x137**(-3)
      x366 = x352*x365
      x367 = x144*x19
      x368 = 5.68327056421864d0*x367
      x369 = -0.245694361157641d0*x276 + 89.4158322713032d0*x282 - 5.56d
     &0*x342
      x370 = x148*x345
      x371 = x350 + x357 + x358 + 426.525248520175d0*x364 + x366*x368 + 
     &x369*x370
      x372 = x149**(-2)
      x373 = x143*x372
      x374 = x371*x373
      x375 = x363*x374
      x376 = gamma_cc*x375
      x377 = -x271 - x347 + x362 + x376
      x378 = 1d0/x156
      x379 = 0.0310906908696549d0*x378
      x380 = x341 + x377*x379
      x381 = x21*x26
      x382 = gamma_cc*x18
      x383 = 0.0185536289197316d0*x382
      x384 = x381*x383
      x385 = x16*x382
      x386 = -x23 - 0.0371072578394632d0*x385
      x387 = gamma_cc*x237
      x388 = x386*x387
      x389 = x29 + x384 + x388
      x390 = x242*x389
      x391 = -x390
      x392 = x53*x57
      x393 = 0.0478125d0*x392
      x394 = x42*x47
      x395 = 0.0742145156789264d0*x382
      x396 = x394*x395
      x397 = x382*x41
      x398 = -0.148429031357853d0*x397 - x44
      x399 = gamma_cc*x199
      x400 = x398*x399
      x401 = x396 + x400 + x50
      x402 = x203*x401
      x403 = 0.000371630914301344d0*x58
      x404 = x390*x60
      x405 = x33*x57
      x406 = 0.0165625d0*x405
      x407 = gamma_cc*x406 + x404
      x408 = -gamma_cc*x393 - x402*x403 + x407
      x409 = x17*x68
      x410 = 0.005859375d0*x66
      x411 = 0.133875d0*x392
      x412 = 0.00104056656004376d0*x58
      x413 = x257*x390
      x414 = 0.046375d0*x405
      x415 = gamma_cc*x414 + x413
      x416 = -gamma_cc*x411 - x402*x412 + x415
      x417 = x261*x416 + x409*x410
      x418 = x54*(x391 + x408*x71 + x417*x62)
      x419 = x138*x144
      x420 = x150*x419
      x421 = 2939.73099289485d0*x382
      x422 = x420*x421
      x423 = x146*x382
      x424 = 182.796535080075d0*x423
      x425 = -x141 - x424
      x426 = x373*x425
      x427 = gamma_cc*x363
      x428 = x426*x427
      x429 = x154 + x422 + x428
      x430 = x379*x429
      x431 = x418 + x430
      x432 = x409*x55
      x433 = 0.0165625d0*x432
      x434 = -x33*x433
      x435 = x432*x53
      x436 = x434 + 0.0478125d0*x435
      x437 = -0.174107868870067d0*x10 + 0.0870539344350337d0*x31
      x438 = -x433*x437
      x439 = 0.133875d0*x435 + x438
      x440 = tau_c**(-4)
      x441 = x440*x67
      x442 = x261*x439 - x410*x441
      x443 = x54*(x436*x71 + x442*x62)
      x444 = x160*x206
      x445 = x254*x56
      x446 = x17*x445
      x447 = x446*x53
      x448 = x1*x55
      x449 = x37*x448
      x450 = x167*x168
      x451 = x169*x450
      x452 = x164*x451
      x453 = x3*x55
      x454 = x4*x55
      x455 = x55*x6
      x456 = x170*(-2.15292298491727d0*x448 - 2.42709420076615d0*x453 - 
     &1.74450360892561d0*x454 - 0.424342513187885d0*x455)
      x457 = 1034.52086599683d0*x167**2/(x35**4*x36**2)
      x458 = 32.1639684429149d0*x450*(x0*x34 + 4.16073291559911d0*x163 +
     & 3.2293844773759d0*x164 + 0.509211015825462d0*x166)/x35**3
      x459 = x51**(-2)
      x460 = 0.0155453454348274d0*x202*x459
      x461 = gamma_cc*rho_c**(-4.33333333333333d0)
      x462 = 2.11884876390004d0*x461
      x463 = 1.27130925834003d0*gamma_cc
      x464 = x172*x463
      x465 = 0.544846825002868d0*x20
      x466 = gamma_cc*x465
      x467 = x196*x197
      x468 = 0.296858062715706d0*x194
      x469 = gamma_cc*x49
      x470 = x45*x469/x46**3
      x471 = x196*x470
      x472 = x42*x462
      x473 = x183*x464
      x474 = exp(3.999999999999988d0*x38)
      x475 = x181**2
      x476 = gamma_cc*x475
      x477 = x193*x465*x474*x476
      x478 = x39*x41*x49
      x479 = x476*x478
      x480 = x180*x456 - x180*x457 + x180*x458 + 0.142757027587911d0*x44
     &9 - 6.88744879551287d0*x452
      x481 = gamma_cc*x478*x480
      x482 = gamma_cc*x187
      x483 = rho_c**(-6.66666666666667d0)*x17
      x484 = x193*x195*x39
      x485 = x162*x456 - x162*x457 + x162*x458 + x204*(-x186*x464*x47 + 
     &x186*x466* x467 - x196*x198*x464 + x399*(-0.445287094073559d0*x19*
     &x474*x475/ x40**4 - 1.38533762600663d0*x191*x194 - 1.9625616368427
     &2d0*x41* x483 - x472 - x473 + x475*x484 - x477 + x479 + x480*x484 
     &+ x481) + x462*x48 + x471*(x176*x463 + x183*x466 + x19*x468 + 0.69
     &2668813003313d0*x192) + x482*(x472 + x473 + x477 - x479 - x481)) -
     & 0.00221920730710327d0*x449 + 0.107067770690934d0*x452 + x460*(x17
     &5 - x189 - x201)
      x486 = x485*x58
      x487 = 1.05942438195002d0*x461
      x488 = 0.635654629170012d0*gamma_cc
      x489 = x172*x488
      x490 = 0.272423412501434d0*x20
      x491 = gamma_cc*x490
      x492 = x234*x235
      x493 = 0.0742145156789264d0*x232
      x494 = gamma_cc*x28
      x495 = x24*x494/x25**3
      x496 = x234*x495
      x497 = x21*x487
      x498 = x223*x489
      x499 = 4.0d0*x10
      x500 = exp(x499)
      x501 = x221**2
      x502 = gamma_cc*x501
      x503 = x231*x490*x500*x502
      x504 = x14*x16*x28
      x505 = x502*x504
      x506 = x448*x9
      x507 = x164*x280
      x508 = -0.989141964955861d0*x448 - 1.16327369584438d0*x453 - 0.600
     &321476249523d0*x454 - 0.210778258416395d0*x455
      x509 = x219*x508
      x510 = x216**2/(x7**4*x8**2)
      x511 = 258.630216499205d0*x510
      x512 = x216*x217*(x0*x5 + 1.99418347859036d0*x163 + 1.483712947433
     &79d0*x164 + 0.252933910099675d0*x166)/x7**3
      x513 = 16.0819842214574d0*x512
      x514 = x215*x509 - x215*x511 + x215*x513 + 0.117839022137957d0*x50
     &6 - 2.8426279420419d0*x507
      x515 = gamma_cc*x504*x514
      x516 = gamma_cc*x227
      x517 = x14*x231*x233
      x518 = x242*(-x226*x26*x489 + x226*x491*x492 - x234*x236*x489 + x2
     &7*x487 + x387 *(-0.49064040921068d0*x16*x483 - 0.346334406501656d0
     &*x191*x232 - x497 - x498 + x501*x517 - x503 + x505 + x514*x517 + x
     &515 - 0.11132177351839d0*x19*x500*x501/x15**4) + x496*(x19*x493 + 
     &x211* x488 + x223*x491 + 0.173167203250828d0*x230) + x516*(x497 + 
     &x498 + x503 - x505 - x515))
      x519 = x30**(-2)
      x520 = 0.0310906908696549d0*x519
      x521 = x240*x520
      x522 = x521*(x210 - x229 - x239)
      x523 = 0.00366369660967365d0*x506
      x524 = 0.0883792666034679d0*x507
      x525 = x246*x509 - x246*x511 + x246*x513 - x523 + x524
      x526 = x518 + x522 + x525
      x527 = x526*x60
      x528 = 0.0496875d0*x17
      x529 = x33*x528
      x530 = x160*x249
      x531 = x445*x529 + x527 - 0.033125d0*x530
      x532 = rho_c**(-5)
      x533 = 0.0234375d0*x532*x69
      x534 = 0.01171875d0*x255
      x535 = x257*x526
      x536 = x33*x446
      x537 = -0.09275d0*x530 + x535 + 0.139125d0*x536
      x538 = x264*x509 - x264*x511 + x264*x513 - x518 - x522 + x523 - x5
     &24
      x539 = 0.0244087849784238d0*rho_c**(-2.66666666666667d0)*x74
      x540 = 0.312466222036391d0*x267
      x541 = 0.208310814690927d0*mu*rho_c**(-2.33333333333333d0)
      x542 = gamma_cc*x172
      x543 = x343*x344
      x544 = x1*x9
      x545 = x543*x544
      x546 = 1434.97848677955d0*x542
      x547 = x272*x278
      x548 = x343**2*x547
      x549 = x130**(-2)
      x550 = x155*x549
      x551 = gamma_cc*x292*x543
      x552 = x10*(-x247 + x263 + x336 + x340)
      x553 = x343*x552
      x554 = x151*x359
      x555 = gamma_cc*x346
      x556 = x274*x455/x2**3
      x557 = x275*x448
      x558 = x218*x508
      x559 = x277*x279
      x560 = x558*x559
      x561 = x277*x510/x9**3
      x562 = x510*x559
      x563 = x130*x273*x279*x507
      x564 = x512*x559
      x565 = x132*x164*x273*x341
      x566 = x131*x281*x341
      x567 = x55*x85
      x568 = x55*x80
      x569 = x55*x82
      x570 = 0.00571025734089666d0*x448 + 0.0350330449527535d0*x455 - 0.
     &00609784657412857d0*x568 + 0.000856522469028669d0*x569
      x571 = x290*x78
      x572 = x0*x280
      x573 = x2*x293
      x574 = x2*x510
      x575 = x0*x292
      x576 = x90*x94
      x577 = 0.89531658564596d0*x55*x92
      x578 = x55*x91
      x579 = x0*x80
      x580 = x106*x316
      x581 = x107*x55
      x582 = x111*x55
      x583 = 3.82831200025092d0*x582
      x584 = -22.9698720015055d0*x115*x322*x326 + x116*x583 + 1.38777878
     &078145d-17* x321*x55 - x324*x327*(-0.150031049774234d0*x448 - 0.02
     &71505358982603d0*x455) - 0.203363546995945d0*x326*x333 - x325 *x32
     &7*(0.225046574661351d0*x164 + 0.0325806430779124d0*x166)/x113 **3
      x585 = x2*x314
      x586 = 0.311172355472072d0*x164*x312 - 0.207448236981382d0*x311*x4
     &48 + 0.0242070587025094d0*x311*x455 + x79*(0.00285512867044833d0*x
     &448 + 0.0175165224763768d0*x455 - x55*x81 + 0.000428261234514334d0
     &* x569)
      x587 = x2*x330
      x588 = x123*x331
      x589 = x77*(-0.000821860773992031d0*x0*x286*x571 - 6.3935088228300
     &2d-5*x0** 2.33333333333333d0*x55*x84 + x120*(-0.0652591516234546d0
     &*x10*x569 - 0.167525641071009d0*x100*x568 - x121*x586 + 0.00091950
     &62450792d0*x126*x567 + x127*(-0.299209d0*x124*x55 + 5.197036119627
     &67d0*x124*x582 - 31.182216717766d0*x322*x588 - x332 *(-0.106975995
     &763985d0*x448 + 0.289053642815522d0*x455) - 5.29055759308205d0*x33
     &3*x588 - x583 + x584) + 0.0298137034090898d0*x280*x285 - 0.0006896
     &296838094d0*x285*x334 + 0.899568554326366d0*x288*x292 + 0.16752564
     &1071009d0*x313*x579 - x512*x587 - x558*x587 - 0.00370771454548634d
     &0*x567*x9 + 5.42506773629845d0*x574*x82) + 0.000279369479946145d0*
     &x283*x572* x87 - 5.79052675961456d-5*x284*x55 + 0.0088039180381261
     &2d0*x285* x571 + 0.00191767513931474d0*x287*x55 - x291*x570*x79 + 
     &x294*x576 *(2.0d0*x301*x303*x304 + x303*x98*(1.53482843253593d0*x2
     &95 + 1.42618991424766d0*x298)/x94**3 + x305*(-x577 - 0.95079327616
     &5107d0*x578) + x95*(1.43555326817173d0*x55*x97 + x577 + 2.05391854
     &581537d0*x578)) + x306*x576*(x296 + x299 + x300) /x98**2 + x307*(x
     &297 - x302) - x512*x573 - 0.00192189252234297d0* x55*x89 - x558*x5
     &73 - 0.0117385573841682d0*x567*x84 + 0.0508356285485896d0*x574*x88
     & + 0.0168588246652672d0*x575*x88 + x96*(0.000586927869208408d0*x11
     &8*x567 + x119*( -0.154912426780983d0*x0**(-0.333333333333333d0)*x5
     &81 - 75.398223686155d0*x0*x108*x580 - x317*(-1.36477107997868d0*x4
     &48 + 0.427594146176807d0*x455) - 1.61108119708458d0*x318*x581 - 4.
     &83324359125373d0*x319*x580 + x584) - 0.000440195901906306d0* x285*
     &x329 - 0.0138644007334803d0*x288*x290 + 0.0891107483661069d0 *x290
     &*x579 - x315*x570 - 0.00107855455836145d0*x567*x83 - 0.08911074836
     &61069d0*x568*x83 + 0.0231073345558005d0*x569*x83) + x99*(-0.022068
     &1498819008d0*x100*x567 - x101*x586 + 0.857572898037877d0*x102*x574
     & + 0.213300270968718d0*x102*x575 - 0.0198949583613061d0*x103*x55 +
     & 0.0165511124114256d0*x285*x313 + 0.00471283037863406d0*x308*x572 
     &- 0.000781467479590561d0*x309*x55 - x512*x585 - x558*x585))
      x590 = x128*x337*x55*x75
      x591 = x335*x338
      x592 = x120*x128*x569/x76**6
      x593 = x525 - x589 + 1.05383986550676d0*x590 - 1.26460783860812d0*
     &x591 - 0.499760307959092d0*x592
      x594 = x133*x593
      x595 = x345*(-0.0108571437238599d0*x556 + 0.16379624077176d0*x557 
     &+ 44.7079161356516d0*x560 - 1437.98400373558d0*x561 - 718.99200186
     &779d0*x562 + 3.95125283943823d0*x563 + 44.7079161356516d0*x564 - 0
     &.245694361157641d0*x565 + 89.4158322713032d0*x566 - 2.78d0*x594)
      x596 = 614.99078004838d0*x354
      x597 = gamma_cc*x596
      x598 = x359*x372
      x599 = gamma_cc*x348
      x600 = 44.6144725370683d0*x599
      x601 = 853.050497040351d0*x364
      x602 = 0.594469188841156d0*gamma_cc*x353*x354
      x603 = 11.3665411284373d0*x366*x367
      x604 = 19.1204882301721d0*x140
      x605 = x345*x369
      x606 = x143/x149**3
      x607 = x371*x427*x606
      x608 = 74.3574542284471d0*x139*x461
      x609 = 1.3870947739627d0*x135*x353*x542
      x610 = exp(64.327936885829629d0*x129 + x499)
      x611 = x351**2
      x612 = gamma_cc*x611
      x613 = 0.0184824577817949d0*x354*x365*x610*x612
      x614 = x136*x145*x355
      x615 = x612*x614
      x616 = 0.422463681720712d0*x545*x599
      x617 = x346*x600
      x618 = x142*x549
      x619 = x548*x618
      x620 = 153.747695012095d0*x140*x551
      x621 = x553*x618
      x622 = x346*x602
      x623 = x514 + 32.1639684429148d0*x589 - 33.8956721780452d0*x590 + 
     &40.6748066136542d0*x591 + 16.0742747742176d0*x592
      x624 = gamma_cc*x614*x623
      x625 = x142*x595
      x626 = gamma_cc*x360
      x627 = x144*x366
      x628 = x136*x365*x368
      x629 = x344*x369
      x630 = x148*x549
      x631 = x373*x427
      x632 = 0.0310906908696549d0/x156**2
      x633 = x377*x632
      x634 = x379*(-x135*x374*x546 - 4945.13601053991d0*x152*x20*x551 + 
     &x152*x346* x546 + 2391.63081129926d0*x152*x461 + 13.5881085271426d
     &0*x152* x542*x545 - x155*x595 - 614.99078004838d0*x20*x554*x555 - 
     &x346* x374*x597 + x371*x597*x598 - x546*x554 + x548*x550 - x550*x5
     &53 + x607*(182.796535080075d0*x147*x605 + x555*x604 + x600 + x601 
     &+ x602 + x603) + x626*(x608 + x609 + x613 - x615 + x616 + x617 + x
     &619 - x620 - x621 + x622 - x624 - x625) + x631*( -2416.97640828099
     &d0*x146*x483 + 1469.86549644743d0*x147*x292*x629 - 53.043858599374
     &d0*x191*x627 - 4.03885592700803d0*x364*x544*x629 - x369**2*x547*x6
     &30 + x369*x552*x630 + x370*( -0.0217142874477198d0*x556 + 0.327592
     &481543521d0*x557 + 89.4158322713032d0*x560 - 2875.96800747116d0*x5
     &61 - 1437.98400373558d0*x562 + 7.90250567887647d0*x563 + 89.415832
     &2713032d0*x564 - 0.491388722315281d0*x565 + 178.831664542606d0*x56
     &6 - 5.56d0*x594) - x601*x605 - x603*x605 - x608 - x609 + x611*x628
     & - x613 + x615 - x616 - x617 - x619 + x620 + x621 - x622 + x623*x6
     &28 + x624 + x625 - 0.530090424722193d0* x367*x610*x611/x137**4)) +
     & x593 + x633*(x271 + x347 - x362 - x376 )
      x635 = x246*x280
      x636 = 0.00549554491451047d0*x214 - 2.0d0*x336 - 1.26460783860812d
     &0*x339 + 0.0621813817393098d0*x377*x378 + 32.1639684429148d0*x635
      x637 = gamma_cc*x159
      x638 = 0.095625d0*x53
      x639 = -x396 - x400 - x50
      x640 = gamma_cc*x190
      x641 = x198*x398
      x642 = x197*x398
      x643 = -0.296858062715706d0*x397 - 0.544846825002868d0*x43
      x644 = x204*(-gamma_cc*x173*x641 - x174 + x186*x469*x642 + x188 + 
     &x200 - 0.173167203250828d0*x394*x640 + x395*x42*x467 + x399*(x177 
     &+ x184 + x382*x468 + 0.692668813003313d0*x41*x640) + x471*x643 + x
     &482*( -x177 - x184)) + x460*x639
      x645 = x58*x644
      x646 = gamma_cc*x57
      x647 = x206*x646
      x648 = x160*x402
      x649 = -x29 - x384 - x388
      x650 = x521*x649
      x651 = x236*x386
      x652 = x235*x386
      x653 = -0.272423412501434d0*x22 - 0.0742145156789264d0*x385
      x654 = x242*(-gamma_cc*x208*x651 - x209 + x21*x383*x492 + x226*x49
     &4*x652 + x228 + x238 - 0.043291800812707d0*x381*x640 + x387*( 0.17
     &3167203250828d0*x16*x640 + x212 + x224 + x382*x493) + x496* x653 +
     & x516*(-x212 - x224))
      x655 = x650 + x654
      x656 = x60*x655
      x657 = x249*x646
      x658 = x33*x637
      x659 = x241*x389
      x660 = x160*x659
      x661 = x656 + 0.0165625d0*x657 - 0.033125d0*x658 - 0.0005149395675
     &28659d0*x660
      x662 = x254*x409
      x663 = 0.017578125d0*x65
      x664 = x409*x63
      x665 = 0.005859375d0*x260
      x666 = 0.26775d0*x53
      x667 = x257*x655
      x668 = 0.046375d0*x657 - 0.09275d0*x658 - 0.00144183078908025d0*x6
     &60 + x667
      x669 = -x650 - x654
      x670 = 0.156233111018196d0*x267
      x671 = x141*x346
      x672 = x425*x427
      x673 = x372*x419
      x674 = -365.59307016015d0*x423 - x604
      x675 = -x154 - x422 - x428
      x676 = x379*(-gamma_cc*x135*x269*x426 - x154*x346 - x270 - x346*x4
     &22 - x346* x428 + x361 + x371*x421*x673 + x375 - 6859.37231675466d
     &0*x420* x640 + x598*x672 + x607*x674 + x626*(-x349 - x356 - x671) 
     &+ x631* (853.050497040351d0*x146*x640 + x349 + x356 + 11.366541128
     &4373d0* x382*x627 + x424*x605 + x671)) + x633*x675
      x677 = gamma_cc*rho_c**(-7.0d0)
      x678 = 0.148429031357853d0*x18*x394 - 0.0404355432367886d0*x198*x4
     &1*x677 + 0.148429031357853d0*x382*x42*x642 + x398*x470*x643 + x465
     &*x641
      x679 = x203*x678
      x680 = x402*x646
      x681 = x401*x459*x639
      x682 = x646*x659
      x683 = -0.00505444290459857d0*x16*x236*x677 + 0.0371072578394632d0
     &*x18*x381 + 0.0371072578394632d0*x21*x382*x652 + x386*x495*x653 + 
     &x490*x651
      x684 = x242*x683
      x685 = x389*x649
      x686 = x520*x685
      x687 = x60*x684 + x60*x686
      x688 = x406 + 0.00102987913505732d0*x682 + x687
      x689 = gamma_cc*x68
      x690 = 0.01171875d0*x689
      x691 = 0.01171875d0*x664
      x692 = x257*x684 + x257*x686
      x693 = x414 + 0.00288366157816049d0*x682 + x692
      x694 = -x684 - x686
      x695 = x379*(-56209.0918495182d0*x134**8.34d0*x145*x373*x677 + 587
     &9.4619857897d0*x18*x420 + 5879.4619857897d0*x382*x425*x673 + x426*
     &x596 + x606*x672*x674) + x429*x632*x675
      x696 = x206*x432
      x697 = x33*x664
      x698 = -x249*x433 + 0.033125d0*x697
      x699 = x254*x441
      x700 = x249*x432
      x701 = 0.09275d0*x697 - 0.046375d0*x700
      x702 = x441*x63
      x703 = x55*x689
      x704 = x402*x432
      x705 = x33*x703
      x706 = x432*x659
      x707 = -0.033125d0*x705 - 0.000514939567528659d0*x706
      x708 = 0.005859375d0*x664
      x709 = x17*x440
      x710 = 0.017578125d0*x709
      x711 = -0.09275d0*x705 - 0.00144183078908025d0*x706
      x712 = 0.005859375d0*x702
      x713 = rho_c*x54
      x714 = x440*x55
      x715 = x529*x714
      x716 = 0.1434375d0*x709
      x717 = x53*x55
      x718 = 0.01171875d0*x702
      x719 = x437*x528*x714
      x720 = 0.401625d0*x709
      x721 = 0.0234375d0*x67/tau_c**5
      x722 = 0.5d0*rho_c
      x723 = x33*x722 + x53*x722
      x724 = x159*x723
      x725 = x17*x724
      x726 = x61 - 0.02390625d0*x725
      x727 = x64 - 0.0669375d0*x725
      x728 = x261*x727 + 1.0d0
      x729 = x54*(x72 + x726*x728)
      x730 = x157 + x729
      x731 = rho_c*x730
      x732 = x446*x723
      x733 = -0.0310906908696549d0*x10 + x206*x722 + x249*x722 + 0.01554
     &53454348275d0 *x31 - 0.0155453454348274d0*x38 + 0.0077726727174137
     &2d0*x52
      x734 = x160*x733
      x735 = x252 + 0.07171875d0*x732 - 0.02390625d0*x734
      x736 = x259 + 0.2008125d0*x732 - 0.0669375d0*x734
      x737 = -x256*x727 + x261*x736
      x738 = x54*(x265 + x726*x737 + x728*x735)
      x739 = rho_c*(x380 + x670*x729 + x738) + x730
      x740 = 0.0155453454348275d0*x241
      x741 = x389*x740
      x742 = 0.00777267271741372d0*x203
      x743 = x401*x742
      x744 = rho_c*x741 + rho_c*x743
      x745 = x160*x744
      x746 = 0.0478125d0*x724
      x747 = -gamma_cc*x746 + x407 - 0.02390625d0*x745
      x748 = 0.133875d0*x724
      x749 = -gamma_cc*x748 + x415 - 0.0669375d0*x745
      x750 = x261*x749 + x708*x727
      x751 = x54*(x391 + x726*x750 + x728*x747)
      x752 = x430 + x751
      x753 = rho_c*x752
      x754 = x664*x723
      x755 = x434 + 0.0478125d0*x754
      x756 = x438 + 0.133875d0*x754
      x757 = x261*x756 - x712*x727
      x758 = x54*(x726*x757 + x728*x755)
      x759 = rho_c*x758
      x760 = x446*x733
      x761 = x17*x532*x56*x723
      x762 = x160*(16.0819842214575d0*x162*x451 + x171*(-0.0024966082204
     &9118d0*x1 - 0.0155453454348274d0) + x205 + x220*(-0.00412165868588
     &285d0*x1 - 0.0310906908696549d0) + x243 + x244 + x485*x722 + x526*
     &x722 + 8.0409921107287d0*x635)
      x763 = rho_c*(x539*x729 + x54*(x538 + x726*(x261*(x537 + 0.401625d
     &0*x760 - 0.80325d0*x761 - 0.0669375d0*x762) + x533*x727 - x534*x73
     &6) + x728*(x531 + 0.1434375d0*x760 - 0.286875d0*x761 - 0.02390625d
     &0* x762) + 2.0d0*x735*x737) + x540*x738 - x541*x729 + x634) + x540
     &* x729 + x636 + 2.0d0*x738
      x764 = x446*x744
      x765 = 0.1434375d0*x723
      x766 = gamma_cc*x445
      x767 = x160*(x644*x722 + x655*x722 + x741 + x743)
      x768 = x637*x733
      x769 = 0.017578125d0*x727
      x770 = 0.401625d0*x723
      x771 = rho_c*(x54*(x669 + x726*(-x256*x749 + x261*(x668 + 0.200812
     &5d0*x764 + x766*x770 - 0.0669375d0*x767 - 0.133875d0*x768) - x662*
     &x769 + x708*x736) + x728*(x661 + 0.07171875d0*x764 + x765*x766 - 0
     &.02390625d0*x767 - 0.0478125d0*x768) + x735*x750 + x737*x747) + x6
     &70*x751 + x676) + x752
      x772 = x519*x685
      x773 = x160*(rho_c*x678*x742 + 0.00777267271741372d0*rho_c*x681 + 
     &rho_c*x683* x740 + 0.0155453454348275d0*rho_c*x772)
      x774 = x637*x744
      x775 = x63*x727
      x776 = rho_c*(x54*(x694 + x726*(x261*(x693 - x748 - 0.0669375d0*x7
     &73 - 0.26775d0*x774) + x690*x775 + x691*x749) + x728*(x688 - x746 
     &- 0.02390625d0*x773 - 0.095625d0*x774) + 2.0d0*x747*x750) + x695)
      x777 = x664*x733
      x778 = rho_c*(x54*(x726*(-x256*x756 + x261*(-x662*x770 + x701 + 0.
     &133875d0*x777 ) + x699*x769 - x712*x736) + x728*(-x662*x765 + x698
     & + 0.0478125d0*x777) + x735*x757 + x737*x755) + x670*x758) + x758
      x779 = x664*x744
      x780 = x63*x723
      x781 = x689*x780
      x782 = x713*(x726*(x261*(x711 + 0.133875d0*x779 + 0.26775d0*x781) 
     &+ x708*x756 - x710*x775 - x712*x749) + x728*(x707 + 0.0478125d0*x7
     &79 + 0.095625d0*x781) + x747*x757 + x750*x755)
      x783 = x713*(x726*(x261*(x719 - x720*x780) - x718*x756 + x721*x775
     &) + x728*( x715 - x716*x780) + 2.0d0*x755*x757)
      x784 = x33*x58
      x785 = x61 - 0.02390625d0*x784
      x786 = x64 - 0.0669375d0*x784
      x787 = x261*x786 + 1.0d0
      x788 = x54*(x72 + x785*x787)
      x789 = x157 + x788
      x790 = x249*x58
      x791 = x250 + 0.03125d0*x251 - 0.02390625d0*x790
      x792 = 0.0875d0*x251 + x258 - 0.0669375d0*x790
      x793 = -x256*x786 + x261*x792
      x794 = x54*(x265 + x785*x793 + x787*x791)
      x795 = 0.03125d0*x405
      x796 = x58*x659
      x797 = -gamma_cc*x795 + x404 - 0.000743261828602687d0*x796
      x798 = 0.0875d0*x405
      x799 = -gamma_cc*x798 + x413 - 0.00208113312008752d0*x796
      x800 = x261*x799 + x708*x786
      x801 = x54*(x391 + x785*x800 + x787*x797)
      x802 = x430 + x801
      x803 = x33*x432
      x804 = 0.03125d0*x803
      x805 = x438 + 0.133875d0*x803
      x806 = x261*x805 - x712*x786
      x807 = x54*(x785*x806 + x787*x804)
      x808 = x526*x58
      x809 = x58*x655
      x810 = 0.017578125d0*x786
      x811 = 0.000743261828602687d0*x58
      x812 = x241*x683
      x813 = x63*x786
      x814 = 0.00208113312008752d0*x58
      x815 = x33*x55
      PBE_value = x33
      PBE_beta_value = x53
      PBE_alpha_value = x53
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*x158
      d1E(1) = rho_c*(x266 + 0.156233111018196d0*x268 + x380) + x158
      d1E(3) = rho_c*x431
      d1E(6) = rho_c*x443
      d2E(1) = rho_c*(x266*x540 + x539*x73 + x54*(2.0d0*x253*x262 + x538
     & + x62*(-x260* x534 + x261*(0.26775d0*x444 - 0.401625d0*x447 - 0.0
     &669375d0*x486 + x537) + x533*x65) + x71*(0.095625d0*x444 - 0.14343
     &75d0*x447 - 0.02390625d0*x486 + x531)) - x541*x73 + x634) + 2.0d0*
     &x266 + 0.312466222036391d0*x268 + x636
      d2E(4) = rho_c*(x418*x670 + x54*(x253*x417 + x262*x408 + x62*(-x25
     &6*x416 + x261*( x637*x666 - 0.0669375d0*x645 - 0.133875d0*x647 + 0
     &.00208113312008752d0*x648 + x668) - x662*x663 + x664*x665) + x669 
     &+ x71*(x637*x638 - 0.02390625d0*x645 - 0.0478125d0*x647 + 0.000743
     &261828602687d0*x648 + x661)) + x676) + x431
      d2E(6) = rho_c*(x54*(2.0d0*x408*x417 + x62*(x261*(-x411 - x412*x67
     &9 - x412*x681 - 0.00416226624017505d0*x680 + x693) + x416*x691 + x
     &66*x690) + x694 + x71*(-x393 - x403*x679 - x403*x681 - 0.001486523
     &65720537d0*x680 + x688)) + x695)
      d2E(16) = rho_c*(x443*x670 + x54*(x253*x442 + x262*x436 + x62*(-x2
     &56*x439 + x261*( -x664*x666 + 0.133875d0*x696 + x701) + x663*x699 
     &- x665*x702) + x71*(-x638*x664 + 0.0478125d0*x696 + x698))) + x443
     &
      d2E(18) = x713*(x408*x442 + x417*x436 + x62*(x261*(x666*x703 + 0.0
     &0208113312008752d0*x704 + x711) - x416*x712 + x439*x708 - x66* x71
     &0) + x71*(x638*x703 + 0.000743261828602687d0*x704 + x707))
      d2E(21) = x713*(2.0d0*x436*x442 + x62*(x261*(-x717*x720 + x719) - 
     &x439*x718 + x66* x721) + x71*(x715 - x716*x717))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = x731
      d1E(1) = x739
      d1E(3) = x753
      d1E(6) = x759
      d2E(1) = x763
      d2E(4) = x771
      d2E(6) = x776
      d2E(16) = x778
      d2E(18) = x782
      d2E(21) = x783
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = x731
      d1E(1) = x739
      d1E(3) = x753
      d1E(6) = x759
      d2E(1) = x763
      d2E(4) = x771
      d2E(6) = x776
      d2E(16) = x778
      d2E(18) = x782
      d2E(21) = x783
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*x789
      d1E(1) = rho_c*(x380 + x670*x788 + x794) + x789
      d1E(3) = rho_c*x802
      d1E(6) = rho_c*x807
      d2E(1) = rho_c*(x539*x788 + x54*(x538 + x785*(x261*(0.175d0*x530 +
     & x535 - 0.2625d0*x536 - 0.0669375d0*x808) + x533*x786 - x534*x792)
     & + x787 *(x527 + 0.0625d0*x530 - 0.09375d0*x536 - 0.02390625d0*x80
     &8) + 2.0d0*x791*x793) + x540*x794 - x541*x788 + x634) + x540*x788 
     &+ x636 + 2.0d0*x794
      d2E(4) = rho_c*(x54*(x669 + x785*(-x256*x799 + x261*(-0.0875d0*x65
     &7 + 0.175d0* x658 + 0.0027204354510948d0*x660 + x667 - 0.0669375d0
     &*x809) - x662*x810 + x708*x792) + x787*(x656 - 0.03125d0*x657 + 0.
     &0625d0* x658 + 0.000971584089676716d0*x660 - 0.02390625d0*x809) + 
     &x791* x800 + x793*x797) + x670*x801 + x676) + x802
      d2E(6) = rho_c*(x54*(x694 + x785*(x261*(-0.00544087090218961d0*x68
     &2 + x692 - x772 *x814 - x798 - x812*x814) + x690*x813 + x691*x799)
     & + x787*( -0.00194316817935343d0*x682 + x687 - x772*x811 - x795 - 
     &x811*x812 ) + 2.0d0*x797*x800) + x695)
      d2E(16) = rho_c*(x54*(x785*(-x256*x805 + x261*(-0.175d0*x697 + 0.0
     &875d0*x700) + x699*x810 - x712*x792) + x787*(-0.0625d0*x697 + 0.03
     &125d0*x700) + x791*x806 + x793*x804) + x670*x807) + x807
      d2E(18) = x713*(x785*(x261*(0.175d0*x705 + 0.0027204354510948d0*x7
     &06) + x708*x805 - x710*x813 - x712*x799) + x787*(0.0625d0*x705 + 0
     &.000971584089676716d0*x706) + x797*x806 + x800*x804)
      d2E(21) = x713*(-0.09375d0*x709*x787*x815 + x785*(x261*(x719 - x72
     &0*x815) - x718* x805 + x721*x813) + 0.0625d0*x803*x806)
      end if
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, gam
     &ma_ss, gamma_cs, tau_c, mu, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x105, x106, x107, x108, x109, x110, x111, x112, x113
      E = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = 0.0621813817393098d0*(0.132568899905202d0*x1 + 1.0d0)*log(1.0
     &d0 + 16.081984221457407d0/(2.2255694211506878d0*x1 + 5.98255043577
     &10795d0*x2 + 0.80042863499936356d0*x3 + 0.18970043257475591d0*x4))
     &
      x6 = rho_s**4/rho_c**4
      x7 = rho_s*x0
      x8 = x7 + 1.0d0
      x9 = x8**1.33333333333333d0
      x10 = -x7 + 1.0d0
      x11 = x10**1.33333333333333d0
      x12 = 1.92366105093154d0*x11 + 1.92366105093154d0*x9 - 3.847322101
     &86307d0
      x13 = 0.0197517897025652d0*x12*(0.0690139921125583d0*x1 + 1.0d0)*(
     &-x6 + 1.0d0) *log(1.0d0 + 29.608574643216677d0/(2.2475918635776164
     &d0*x1 + 8.1574147034876408d0*x2 + 0.43009724712766434d0*x3 + 0.191
     &15125951273382d0*x4))
      x14 = x12*x6*(x5 - 0.0310906908696548d0*(0.127469618870009d0*x1 + 
     &1.0d0)*log( 1.0d0 + 32.163968442914921d0/(3.8447462374472119d0*x1 
     &+ 11.120374863094685d0*x2 + 1.6447337755676092d0*x3 + 0.2405871291
     &2881922d0*x4)))
      x15 = -x13 - x14 + x5
      x16 = x8**0.666666666666667d0
      x17 = x10**0.666666666666667d0
      x18 = (0.5d0*x16 + 0.5d0*x17)**3
      x19 = 32.1639684429148d0/x18
      x20 = exp(x15*x19) - 1.0d0
      x21 = rho_c**(-4.66666666666667d0)
      x22 = gamma_cc**2
      x23 = x16 + x17
      x24 = x21*x22/x23**4
      x25 = gamma_cc*rho_c**(-2.33333333333333d0)/x23**2
      x26 = 0.544846825002868d0*x25
      x27 = 1.0d0 + x26/x20
      x28 = 0.0310906908696549d0*x18
      x29 = x28*log(x26*x27/(x27 + 0.29685806271570553d0*x24/x20**2) + 1
     &.0d0)
      x30 = x13 + x14 - x5
      x31 = x29 + x30
      x32 = 6.28318530717959d0*rho_c
      x33 = 6.28318530717959d0*rho_s
      x34 = 1d0/(x32 + x33)
      x35 = x34**0.333333333333333d0
      x36 = (0.296353441706766d0*x35 + 1.0d0)*log(1.0d0 + 32.16396844291
     &4921d0/( 16.955908776434466d0*x34**0.16666666666666666d0 + 5.83042
     &94284383547d0*sqrt(x34) + 1.3004060036573588d0*x34** 0.66666666666
     &666663d0 + 8.9386301618942241d0*x35))
      x37 = 0.5d0*rho_c
      x38 = 0.5d0*rho_s
      x39 = x37 + x38
      x40 = 14.804406601634d0*rho_c
      x41 = 14.804406601634d0*rho_s
      x42 = x40 + x41
      x43 = 0.5d0*gamma_cs
      x44 = 0.25d0*gamma_cc + 0.25d0*gamma_ss
      x45 = x43 + x44
      x46 = exp(1.999999999999994d0*x36) - 1.0d0
      x47 = 0.668920897233289d0*x42**(-0.333333333333333d0)*x45/x39**2
      x48 = 1.0d0 + x47/x46
      x49 = -0.0310906908696548d0*x36 + 0.0155453454348274d0*log(x47*x48
     &/(x48 + 0.44745516675538899d0*x42**(-0.66666666666666663d0)*x45**2
     &/(x39** 4*x46**2)) + 1.0d0)
      x50 = 1d0/(x32 - x33)
      x51 = x50**0.333333333333333d0
      x52 = (0.296353441706766d0*x51 + 1.0d0)*log(1.0d0 + 32.16396844291
     &4921d0/( 16.955908776434466d0*x50**0.16666666666666666d0 + 5.83042
     &94284383547d0*sqrt(x50) + 1.3004060036573588d0*x50** 0.66666666666
     &666663d0 + 8.9386301618942241d0*x51))
      x53 = x37 - x38
      x54 = x40 - x41
      x55 = -x43 + x44
      x56 = exp(1.999999999999994d0*x52) - 1.0d0
      x57 = 0.668920897233289d0*x54**(-0.333333333333333d0)*x55/x53**2
      x58 = 1.0d0 + x57/x56
      x59 = -0.0310906908696548d0*x52 + 0.0155453454348274d0*log(x57*x58
     &/(x58 + 0.44745516675538899d0*x54**(-0.66666666666666663d0)*x55**2
     &/(x53** 4*x56**2)) + 1.0d0)
      x60 = exp(-0.46869933305458689d0*mu*rho_c**(-0.33333333333333331d0
     &))
      x61 = x22/tau_c**2
      x62 = 0.015625d0*x61
      x63 = x53*x59
      x64 = x39*x49
      x65 = rho_c**(-3)
      x66 = rho_c**2
      x67 = 1d0/x66
      x68 = rho_s**2
      x69 = x67*x68
      x70 = (0.5d0*x6 + 0.87d0*x69 + 0.53d0 + 2.26d0*rho_s**6/rho_c**6)/
     &( 0.0130605864926168d0*x21*(1d0/x9 + 1d0/x11)*(gamma_cc*x68 - 2.0d
     &0 *gamma_cs*rho_c*rho_s + gamma_ss*x66) + 1.0d0)**4
      x71 = x65*(x70 + 1.0d0)
      x72 = x71*(x63 + x64)
      x73 = x31*(x62*x67*x70 + 1.0d0)
      x74 = 0.04375d0*x61
      x75 = 2.8d0*x73
      x76 = 0.001953125d0*gamma_cc**3*x65/tau_c**3
      x77 = x15 - x29
      x78 = mu**2
      x79 = x0**1.66666666666667d0
      x80 = x0**1.0d0
      x81 = x0**1.33333333333333d0
      x82 = (0.0128480790170175d0*x1 + 0.0315297404574782d0*x4 - 0.00304
     &892328706429d0*x80 + 0.000275310793616358d0*x81 + 1.0d0)* exp(-0.4
     &6675853320810856d0*x1)
      x83 = x82*(x69 - 1.0d0)
      x84 = 0.5d0*x10**0.666666666666667d0 + 0.5d0*x8**0.666666666666667
     &d0
      x85 = x84**3
      x86 = x1*x78/x84**2
      x87 = 4.60448529760779d0*mu*x2/x84 + 1.0d0
      x88 = mu**3
      x89 = -x69 + 1.0d0
      x90 = 1d0/x8
      x91 = x90**0.666666666666667d0
      x92 = 0.5d0*x7
      x93 = x1*x90**0.333333333333333d0
      x94 = 0.0244354823084343d0*x4
      x95 = 1d0/x4
      x96 = 34.4548080022583d0*x95
      x97 = 1d0/x10
      x98 = x97**0.666666666666667d0
      x99 = x1*x97**0.333333333333333d0
      x100 = -x96*(-x92 + 0.5d0)**2*(-0.0177069812999059d0*x99 + 1.0d0)/
     &(x98*(x94*x98 + 0.337569861992026d0*x99 + 1.0d0)) - x96*(x92 + 0.5
     &d0)**2*( -0.0177069812999059d0*x93 + 1.0d0)/(x91*(x91*x94 + 0.3375
     &69861992026d0*x93 + 1.0d0))
      x101 = (-3.0d0*x69 + 3.0d0)*(0.5d0*x82 - 0.5d0)
      x102 = (0.00316102962473761d0*mu**8*x0**2.66666666666667d0*x30 + m
     &u**6*( 0.0533250677421794d0*x0**2.0d0*x30 - 0.00165511124114256d0*
     &x101* x79) + 0.00264117541143784d0*mu**5*x79*x83 + mu**4*( -0.0279
     &209401785015d0*x101*x80 + 0.337338207872387d0*x30*x81 + 0.00020688
     &890514282d0*x79*(x100 - 23.3866625383245d0*x89*x95*( -0.2406959904
     &68967d0*x1 + 0.26014827853397d0*x4)*exp( -0.33933171852197186d0*x1
     &) + 17.2274040011292d0*x95*(0.5d0*x10** 2.66666666666667d0 + 0.5d0
     &*x8**2.66666666666667d0))) - 0.0621813817393098d0*x85*log((1.91407
     &10242289749d0*x3*x88/x85 + 4.6213167280845804d0*x86 + x87)/(2.1392
     &848713714905d0*x86 + x87 )) + x88*(0.000132058770571892d0*x79*(x10
     &0 - 37.6991118430775d0* x89*(-3.07073492995203d0*x1 + 0.3848347315
     &59127d0*x4)*exp( -0.19230865217881402d0*x1)/x80) + 0.0445553741830
     &535d0*x80*x83))/ (0.237113969739022d0*x4*x78 + 1.0d0)**4
      x103 = -x102 + x30
      x104 = x103/x30
      x105 = 0.0310906908696549d0*exp(x19*(x102 + x15)) - 0.031090690869
     &6549d0
      x106 = x104**2.78d0*x25
      x107 = 1.0d0 + 0.0169396642074771d0*x106/x105
      x108 = x103 + x28*log(0.54484682500286752d0*x106*x107/(0.000286952
     &2234620814d0* x104**5.5599999999999996d0*x24/x105**2 + x107) + 1.0
     &d0)
      x109 = x31*x53
      x110 = x71*(x109 + x64)
      x111 = x31*x39
      x112 = x71*(x111 + x63)
      x113 = x71*(x109 + x111)
      PBE_value = x31
      PBE_beta_value = x49
      PBE_alpha_value = x59
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*(x108 + x60*(x77 + (-x62*x72 + x73)*(x76*(-x72*x74 + x75
     &) + 1.0d0 )))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*(x108 + x60*(x77 + (-x110*x62 + x73)*(x76*(-x110*x74 + x
     &75) + 1.0d0)))
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*(x108 + x60*(x77 + (-x112*x62 + x73)*(x76*(-x112*x74 + x
     &75) + 1.0d0)))
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*(x108 + x60*(x77 + (-x113*x62 + x73)*(x76*(-x113*x74 + x
     &75) + 1.0d0)))
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, g
     &amma_ss, gamma_cs, tau_c, mu, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &, x754, x755, x756, x757, x758, x759, x760
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = x0**0.333333333333333d0
      x2 = x0**0.166666666666667d0
      x3 = sqrt(x0)
      x4 = x0**0.666666666666667d0
      x5 = 2.22556942115069d0*x1 + 5.98255043577108d0*x2 + 0.80042863499
     &9364d0*x3 + 0.189700432574756d0*x4
      x6 = 1.0d0 + 16.0819842214574d0/x5
      x7 = log(x6)
      x8 = x7*(0.132568899905202d0*x1 + 1.0d0)
      x9 = 0.0621813817393098d0*x8
      x10 = rho_s*x0
      x11 = x10 + 1.0d0
      x12 = x11**1.33333333333333d0
      x13 = -x10 + 1.0d0
      x14 = x13**1.33333333333333d0
      x15 = 1.92366105093154d0*x12 + 1.92366105093154d0*x14 - 3.84732210
     &186307d0
      x16 = rho_c**(-4)
      x17 = rho_s**4
      x18 = x16*x17
      x19 = -x18 + 1.0d0
      x20 = 2.24759186357762d0*x1 + 8.15741470348764d0*x2 + 0.4300972471
     &27664d0*x3 + 0.191151259512734d0*x4
      x21 = 1.0d0 + 29.6085746432167d0/x20
      x22 = log(x21)
      x23 = x19*x22
      x24 = x15*x23
      x25 = x24*(0.0690139921125583d0*x1 + 1.0d0)
      x26 = 0.0197517897025652d0*x25
      x27 = 3.84474623744721d0*x1 + 11.1203748630947d0*x2 + 1.6447337755
     &6761d0*x3 + 0.240587129128819d0*x4
      x28 = 1.0d0 + 32.1639684429149d0/x27
      x29 = log(x28)
      x30 = -0.0310906908696548d0*x29*(0.127469618870009d0*x1 + 1.0d0) +
     & x9
      x31 = x18*x30
      x32 = x15*x31
      x33 = -x26 - x32 + x9
      x34 = x11**0.666666666666667d0
      x35 = x13**0.666666666666667d0
      x36 = 0.5d0*x34 + 0.5d0*x35
      x37 = x36**3
      x38 = 32.1639684429148d0/x37
      x39 = exp(x33*x38)
      x40 = x39 - 1.0d0
      x41 = 1d0/x40
      x42 = x34 + x35
      x43 = x42**(-2)
      x44 = x41*x43
      x45 = rho_c**(-2.33333333333333d0)
      x46 = 0.544846825002868d0*x45
      x47 = x44*x46
      x48 = gamma_cc*x47 + 1.0d0
      x49 = x40**(-2)
      x50 = gamma_cc**2
      x51 = x42**(-4)
      x52 = x50*x51
      x53 = x49*x52
      x54 = rho_c**(-4.66666666666667d0)
      x55 = 0.296858062715706d0*x54
      x56 = x48 + x53*x55
      x57 = 1d0/x56
      x58 = x43*x46
      x59 = x57*x58
      x60 = x48*x59
      x61 = gamma_cc*x60 + 1.0d0
      x62 = log(x61)
      x63 = 0.0310906908696549d0*x37
      x64 = x62*x63
      x65 = x26 - x9
      x66 = x32 + x65
      x67 = x64 + x66
      x68 = 6.28318530717959d0*rho_c
      x69 = 6.28318530717959d0*rho_s
      x70 = 1d0/(x68 + x69)
      x71 = x70**0.333333333333333d0
      x72 = x70**0.166666666666667d0
      x73 = sqrt(x70)
      x74 = x70**0.666666666666667d0
      x75 = 8.93863016189422d0*x71 + 16.9559087764345d0*x72 + 5.83042942
     &843835d0*x73 + 1.30040600365736d0*x74
      x76 = 1.0d0 + 32.1639684429149d0/x75
      x77 = log(x76)
      x78 = x77*(0.296353441706766d0*x71 + 1.0d0)
      x79 = 0.5d0*gamma_cs
      x80 = 0.25d0*gamma_cc
      x81 = 0.25d0*gamma_ss
      x82 = x80 + x81
      x83 = x79 + x82
      x84 = 0.668920897233289d0*x83
      x85 = exp(1.999999999999994d0*x78)
      x86 = x85 - 1.0d0
      x87 = x86**(-2)
      x88 = x83**2
      x89 = x87*x88
      x90 = 0.5d0*rho_c
      x91 = 0.5d0*rho_s
      x92 = x90 + x91
      x93 = x92**(-4)
      x94 = 14.804406601634d0*rho_c
      x95 = 14.804406601634d0*rho_s
      x96 = x94 + x95
      x97 = x96**(-0.666666666666667d0)
      x98 = x93*x97
      x99 = 0.447455166755389d0*x98
      x100 = x92**(-2)
      x101 = 1d0/x86
      x102 = x96**(-0.333333333333333d0)
      x103 = x101*x102
      x104 = x100*x103
      x105 = x104*x84 + 1.0d0
      x106 = x105 + x89*x99
      x107 = 1d0/x106
      x108 = x105*x107
      x109 = x100*x102
      x110 = x108*x109
      x111 = x110*x84 + 1.0d0
      x112 = log(x111)
      x113 = 0.0155453454348274d0*x112 - 0.0310906908696548d0*x78
      x114 = 1d0/(x68 - x69)
      x115 = x114**0.333333333333333d0
      x116 = x114**0.166666666666667d0
      x117 = sqrt(x114)
      x118 = x114**0.666666666666667d0
      x119 = 8.93863016189422d0*x115 + 16.9559087764345d0*x116 + 5.83042
     &942843835d0* x117 + 1.30040600365736d0*x118
      x120 = 1.0d0 + 32.1639684429149d0/x119
      x121 = log(x120)
      x122 = x121*(0.296353441706766d0*x115 + 1.0d0)
      x123 = -x79 + x82
      x124 = 0.668920897233289d0*x123
      x125 = exp(1.999999999999994d0*x122)
      x126 = x125 - 1.0d0
      x127 = x126**(-2)
      x128 = x123**2
      x129 = x127*x128
      x130 = x90 - x91
      x131 = x130**(-4)
      x132 = x94 - x95
      x133 = x132**(-0.666666666666667d0)
      x134 = x131*x133
      x135 = 0.447455166755389d0*x134
      x136 = x130**(-2)
      x137 = 1d0/x126
      x138 = x132**(-0.333333333333333d0)
      x139 = x137*x138
      x140 = x136*x139
      x141 = x124*x140 + 1.0d0
      x142 = x129*x135 + x141
      x143 = 1d0/x142
      x144 = x141*x143
      x145 = x136*x138
      x146 = x144*x145
      x147 = x124*x146 + 1.0d0
      x148 = log(x147)
      x149 = -0.0310906908696548d0*x122 + 0.0155453454348274d0*x148
      x150 = exp(-0.46869933305458689d0*mu*rho_c**(-0.33333333333333331d
     &0))
      x151 = rho_c**(-6)
      x152 = rho_s**6
      x153 = 0.5d0*x18
      x154 = rho_c**2
      x155 = 1d0/x154
      x156 = rho_s**2
      x157 = x155*x156
      x158 = 2.26d0*x151*x152 + x153 + 0.87d0*x157 + 0.53d0
      x159 = 1d0/x14 + 1d0/x12
      x160 = 2.0d0*gamma_cs
      x161 = rho_s*x160
      x162 = gamma_cc*x156 + gamma_ss*x154 - rho_c*x161
      x163 = x159*x162
      x164 = 0.0130605864926168d0*x163*x54 + 1.0d0
      x165 = x164**(-4)
      x166 = x158*x165
      x167 = x166 + 1.0d0
      x168 = x130*x149
      x169 = x113*x92
      x170 = x168 + x169
      x171 = x167*x170
      x172 = rho_c**(-3)
      x173 = tau_c**(-2)
      x174 = x173*x50
      x175 = 0.015625d0*x174
      x176 = x172*x175
      x177 = x155*x166
      x178 = x175*x177 + 1.0d0
      x179 = x178*x67
      x180 = -x171*x176 + x179
      x181 = x172*x174
      x182 = 0.04375d0*x181
      x183 = 2.8d0*x179
      x184 = -x171*x182 + x183
      x185 = gamma_cc**3
      x186 = tau_c**(-3)
      x187 = x185*x186
      x188 = x184*x187
      x189 = 0.001953125d0*x172
      x190 = x188*x189 + 1.0d0
      x191 = x33 - x64
      x192 = x150*(x180*x190 + x191)
      x193 = mu**2
      x194 = x193*x4
      x195 = 0.237113969739022d0*x194 + 1.0d0
      x196 = x195**(-4)
      x197 = mu**5
      x198 = x157 - 1.0d0
      x199 = exp(-0.46675853320810856d0*x1)
      x200 = x0**1.0d0
      x201 = 0.00304892328706429d0*x200
      x202 = x0**1.33333333333333d0
      x203 = x199*(0.0128480790170175d0*x1 - x201 + 0.000275310793616358
     &d0*x202 + 0.0315297404574782d0*x4 + 1.0d0)
      x204 = x198*x203
      x205 = x197*x204
      x206 = x0**1.66666666666667d0
      x207 = 0.00264117541143784d0*x206
      x208 = 0.5d0*x11**0.666666666666667d0 + 0.5d0*x13**0.6666666666666
     &67d0
      x209 = x208**3
      x210 = -0.0621813817393098d0
      x211 = x208**2
      x212 = 1d0/x211
      x213 = x1*x193
      x214 = x212*x213
      x215 = mu/x208
      x216 = 4.60448529760779d0*x2
      x217 = x215*x216 + 1.0d0
      x218 = 2.13928487137149d0*x214 + x217
      x219 = 1d0/x218
      x220 = mu**3
      x221 = 1d0/x209
      x222 = x220*x221
      x223 = 1.91407102422897d0*x3
      x224 = 4.62131672808458d0*x214 + x217 + x222*x223
      x225 = x210*log(x219*x224)
      x226 = mu**8
      x227 = x0**2.66666666666667d0
      x228 = 0.00316102962473761d0*x226*x227
      x229 = 0.0445553741830535d0*x200
      x230 = x204*x229
      x231 = 1d0/x200
      x232 = -3.07073492995203d0*x1 + 0.384834731559127d0*x4
      x233 = exp(-0.19230865217881402d0*x1)
      x234 = x231*x232*x233
      x235 = -x157 + 1.0d0
      x236 = 37.6991118430775d0*x235
      x237 = x234*x236
      x238 = 1d0/x4
      x239 = 34.4548080022583d0*x238
      x240 = 1d0/x11
      x241 = x240**0.666666666666667d0
      x242 = 1d0/x241
      x243 = x240**0.333333333333333d0
      x244 = x1*x243
      x245 = -0.0177069812999059d0*x244 + 1.0d0
      x246 = 0.5d0*x10
      x247 = x246 + 0.5d0
      x248 = x247**2
      x249 = 0.0244354823084343d0*x4
      x250 = x241*x249 + 0.337569861992026d0*x244 + 1.0d0
      x251 = 1d0/x250
      x252 = x248*x251
      x253 = x242*x245*x252
      x254 = 1d0/x13
      x255 = x254**0.666666666666667d0
      x256 = 1d0/x255
      x257 = x254**0.333333333333333d0
      x258 = x1*x257
      x259 = x256*(-0.0177069812999059d0*x258 + 1.0d0)
      x260 = -x246 + 0.5d0
      x261 = x260**2
      x262 = x249*x255 + 0.337569861992026d0*x258 + 1.0d0
      x263 = 1d0/x262
      x264 = x261*x263
      x265 = x239*x264
      x266 = -x239*x253 - x259*x265
      x267 = -x237 + x266
      x268 = 0.000132058770571892d0*x206
      x269 = mu**6
      x270 = -3.0d0*x157 + 3.0d0
      x271 = 0.5d0*x203 - 0.5d0
      x272 = x270*x271
      x273 = 0.00165511124114256d0*x206
      x274 = x0**2.0d0
      x275 = 0.0533250677421794d0*x274
      x276 = mu**4
      x277 = 0.0279209401785015d0*x200
      x278 = x272*x277
      x279 = 0.5d0*x11**2.66666666666667d0 + 0.5d0*x13**2.66666666666667
     &d0
      x280 = 17.2274040011292d0*x238
      x281 = -0.240695990468967d0*x1 + 0.26014827853397d0*x4
      x282 = exp(-0.33933171852197186d0*x1)
      x283 = x235*x281*x282
      x284 = 23.3866625383245d0*x238
      x285 = x266 + x279*x280 - x283*x284
      x286 = 0.00020688890514282d0*x206
      x287 = 0.337338207872387d0*x202
      x288 = x205*x207 + x209*x225 + x220*(x230 + x267*x268) + x228*x66 
     &+ x269*(-x272 *x273 + x275*x66) + x276*(-x278 + x285*x286 + x287*x
     &66)
      x289 = x196*x288
      x290 = x289 + x33
      x291 = exp(x290*x38)
      x292 = 0.0310906908696549d0*x291 - 0.0310906908696549d0
      x293 = x292**(-2)
      x294 = x293*x52
      x295 = 1d0/x66
      x296 = -x289 + x66
      x297 = x295*x296
      x298 = x297**5.56d0
      x299 = x298*x54
      x300 = 0.000286952223462081d0*x299
      x301 = x294*x300
      x302 = 1d0/x292
      x303 = x297**2.78d0
      x304 = x303*x43
      x305 = x302*x304
      x306 = 0.0169396642074771d0*x45
      x307 = x305*x306
      x308 = gamma_cc*x307
      x309 = x308 + 1.0d0
      x310 = x301 + x309
      x311 = 1d0/x310
      x312 = x309*x311
      x313 = x304*x45
      x314 = 0.544846825002868d0*x313
      x315 = x312*x314
      x316 = gamma_cc*x315
      x317 = x316 + 1.0d0
      x318 = log(x317)
      x319 = x296 + x318*x63
      x320 = x192 + x319
      x321 = x70*x71
      x322 = x321*x77
      x323 = 32.1639684429149d0*(18.7210232331754d0*x321 + 17.7561861489
     &951d0*x70*x72 + 18.3168342596557d0*x70*x73 + 5.44712793036536d0*x7
     &0*x74)/(x75** 2*x76)
      x324 = x84/x92**3
      x325 = 3.30099231565716d0*x100*x96**(-1.33333333333333d0)
      x326 = x103*x324
      x327 = x101*x83
      x328 = x325*x327
      x329 = x109*x84
      x330 = x85*(-1.24136239377603d0*x322 + x323*(0.592706883413531d0*x
     &71 + 1.99999999999999d0))
      x331 = x329*x330*x87
      x332 = x105*x329/x106**2
      x333 = 1d0/x111
      x334 = 0.0155453454348274d0*x333
      x335 = 0.00777267271741372d0*x112 - 0.0155453454348274d0*x78 + x92
     &*( 0.0192974072210527d0*x322 + x323*(-0.00921383324426333d0*x71 - 
     &0.0310906908696548d0) + x334*(-x102*x108*x324 + x107*x329*(-x326 -
     & x328 - x331) - x108*x325*x83 + x332*(x326 + x328 + 0.894910333510
     &778d0*x330*x88*x98/x86**3 + x331 + 4.41620548309916d0*x89*x93*x96*
     &*(-1.66666666666667d0) + 0.894910333510778d0*x89*x97/x92**5)))
      x336 = 0.00777267271741372d0*x148
      x337 = x114*x115
      x338 = x121*x337
      x339 = 0.0192974072210527d0*x338
      x340 = -0.00921383324426333d0*x115 - 0.0310906908696548d0
      x341 = 17.7561861489951d0*x114*x116
      x342 = 18.7210232331754d0*x337
      x343 = 18.3168342596557d0*x114*x117
      x344 = 5.44712793036536d0*x114*x118
      x345 = 32.1639684429149d0/(x119**2*x120)
      x346 = x345*(x341 + x342 + x343 + x344)
      x347 = x124/x130**3
      x348 = x138*x144*x347
      x349 = 3.30099231565716d0*x132**(-1.33333333333333d0)*x136
      x350 = x123*x144*x349
      x351 = 1.24136239377603d0*x338
      x352 = 0.592706883413531d0*x115 + 1.99999999999999d0
      x353 = x125*(x346*x352 - x351)
      x354 = x124*x145
      x355 = x127*x354
      x356 = x353*x355
      x357 = x139*x347
      x358 = x123*x137
      x359 = x349*x358
      x360 = -x357 - x359
      x361 = x143*x354
      x362 = 0.894910333510778d0*x129*x133/x130**5
      x363 = 4.41620548309916d0*x129*x131*x132**(-1.66666666666667d0)
      x364 = 0.894910333510778d0*x128*x134/x126**3
      x365 = x357 + x359
      x366 = x141*x354/x142**2
      x367 = 1d0/x147
      x368 = 0.0155453454348274d0*x367
      x369 = 0.0155453454348274d0*x122
      x370 = x130*(x339 + x340*x346 + x368*(-x348 - x350 + x361*(-x356 +
     & x360) + x366 *(x353*x364 + x356 + x362 + x363 + x365))) + x336 - 
     &x369
      x371 = x167*(x335 + x370)
      x372 = x17/rho_c**5
      x373 = x156*x172
      x374 = x165*(-2.0d0*x372 - 1.74d0*x373 - 13.56d0*x152/rho_c**7)
      x375 = 0.0522423459704672d0*x54
      x376 = x159*x375
      x377 = rho_c**(-5.66666666666667d0)
      x378 = x11**(-2.33333333333333d0)
      x379 = rho_s*x155
      x380 = 1.33333333333333d0*x379
      x381 = x13**(-2.33333333333333d0)
      x382 = x162*x375
      x383 = x158/x164**5
      x384 = x383*(0.243797614528847d0*x163*x377 - x376*(2.0d0*gamma_ss*
     &rho_c - x161 ) - x382*(x378*x380 - x380*x381))
      x385 = x374 + x384
      x386 = x170*x385
      x387 = x16*x174
      x388 = x171*x387
      x389 = gamma_cc*rho_c**(-3.33333333333333d0)
      x390 = 1.27130925834002d0*x389
      x391 = x48*x57
      x392 = x42**(-3)
      x393 = x11**(-0.333333333333333d0)
      x394 = x379*x393
      x395 = x13**(-0.333333333333333d0)
      x396 = x379*x395
      x397 = x392*(1.33333333333333d0*x394 - 1.33333333333333d0*x396)
      x398 = gamma_cc*x397
      x399 = x391*x46
      x400 = x390*x44
      x401 = gamma_cc*x41
      x402 = x401*x46
      x403 = x397*x402
      x404 = 32.1639684429148d0/x36**4
      x405 = x404*(x394 - 1.0d0*x396)
      x406 = 0.00136314985874174d0*x1
      x407 = -x406 - 0.0197517897025652d0
      x408 = x11**0.333333333333333d0
      x409 = 2.56488140124205d0*x379
      x410 = x13**0.333333333333333d0
      x411 = -x408*x409 + x409*x410
      x412 = x23*x411
      x413 = x15*x407
      x414 = 4.0d0*x22
      x415 = x372*x414
      x416 = x0*x2
      x417 = x0*x1
      x418 = x0*x3
      x419 = x0*x4
      x420 = 29.6085746432167d0*x19*(1.35956911724794d0*x416 + 0.7491972
     &87859205d0* x417 + 0.215048623563832d0*x418 + 0.127434173008489d0*
     &x419)/(x20 **2*x21)
      x421 = 0.00274777245725523d0*x417*x7
      x422 = -x421
      x423 = 0.0082433173717657d0*x1
      x424 = 16.0819842214574d0*(0.99709173929518d0*x416 + 0.74185647371
     &6896d0*x417 + 0.400214317499682d0*x418 + 0.126466955049837d0*x419)
     &/(x5**2*x6)
      x425 = x422 + x424*(x423 + 0.0621813817393098d0)
      x426 = x31*x411
      x427 = x15*x18*(0.00132103950518672d0*x29*x417 + x425 + 32.1639684
     &429149d0*( -0.00396311851556016d0*x1 - 0.0310906908696548d0)*( 1.8
     &5339581051578d0*x416 + 1.28158207914907d0*x417 + 0.822366887783805
     &d0*x418 + 0.160391419419213d0*x419)/(x27**2*x28 ))
      x428 = x15*x30
      x429 = 4.0d0*x428
      x430 = x372*x429
      x431 = 0.000454383286247248d0*x24*x417
      x432 = -x426 - x427 + x430 + x431
      x433 = x407*x412 + x413*x415 + x413*x420 + x425 + x432
      x434 = x39*(x33*x405 + x38*x433)
      x435 = gamma_cc*x49
      x436 = x435*x58
      x437 = x434*x436
      x438 = gamma_cc*x59
      x439 = x50/x42**5
      x440 = x439*(2.66666666666667d0*x394 - 2.66666666666667d0*x396)
      x441 = x49*x55
      x442 = 0.593716125431411d0*x54
      x443 = x442*x52/x40**3
      x444 = gamma_cc*x48*x58/x56**2
      x445 = x63/x61
      x446 = x445*(-x390*x391*x43 + x398*x399 + x438*(-x400 + x403 - x43
     &7) + x444*( 1.38533762600663d0*x377*x53 + x400 - x403 + x434*x443 
     &+ x437 - x440*x441))
      x447 = 0.0310906908696549d0*x36**2
      x448 = x447*(-1.0d0*x394 + x396)
      x449 = x448*x62
      x450 = x406 + 0.0197517897025652d0
      x451 = x15*x450
      x452 = x415*x451
      x453 = x412*x450
      x454 = x424*(-x423 - 0.0621813817393098d0)
      x455 = x420*x451
      x456 = x421 + x426 + x427 - x430 - x431 + x452 + x453 + x454 + x45
     &5
      x457 = x446 + x449 + x456
      x458 = x155*x374
      x459 = x166*x181
      x460 = x155*x384
      x461 = x178*x457 + x67*(x175*x458 + x175*x460 - 0.03125d0*x459)
      x462 = 0.005859375d0*x16
      x463 = 0.04375d0*x174
      x464 = x177*x463 + 2.8d0
      x465 = x457*x464 + x67*(x458*x463 - 0.0875d0*x459 + x460*x463)
      x466 = x187*x189
      x467 = x433 - x446 - x449
      x468 = 0.156233111018196d0*mu*rho_c**(-1.33333333333333d0)
      x469 = x303*x398
      x470 = 0.544846825002868d0*x312*x45
      x471 = x296/x66**2
      x472 = x471*(x422 + x432 - x452 - x453 - x454 - x455)
      x473 = x203*x373
      x474 = 0.00528235082287567d0*x197*x206
      x475 = x0*x206
      x476 = x0*x202
      x477 = x198*x199*(x0*x201 - 0.00428269300567249d0*x417 - 0.0210198
     &269716521d0* x419 - 0.000367081058155144d0*x476)
      x478 = x13**(-0.333333333333333d0)
      x479 = x379*x478
      x480 = x11**(-0.333333333333333d0)
      x481 = x379*x480
      x482 = x211*x225
      x483 = x0*x66
      x484 = 0.767414216267966d0*x215*x416
      x485 = x193*x212*x417
      x486 = mu*x212*x216
      x487 = x486*(-0.333333333333333d0*x479 + 0.333333333333333d0*x481)
     &
      x488 = x220*x223/x208**4
      x489 = x213*x221
      x490 = x489*(-0.666666666666667d0*x479 + 0.666666666666667d0*x481)
     &
      x491 = x224/x218**2
      x492 = x209*x210*x218/x224
      x493 = x271*x373
      x494 = 0.00993066744685536d0*x206
      x495 = 0.00152446164353214d0*x200
      x496 = x270*(0.155586177736036d0*x199*x417*(0.00642403950850874d0*
     &x1 + 0.000137655396808179d0*x202 + 0.0157648702287391d0*x4 - x495 
     &+ 0.5d0) + x199*(x0*x495 - 0.00214134650283625d0*x417 - 0.01050991
     &34858261d0*x419 - 0.000183540529077572d0*x476))
      x497 = 0.0891107483661069d0*x200
      x498 = 75.398223686155d0*x234
      x499 = 0.00590232709996862d0*x0
      x500 = 0.00590232709996862d0*x379
      x501 = x240*x244
      x502 = x239*x242
      x503 = x254*x258
      x504 = x0*x238
      x505 = 22.9698720015055d0*x504
      x506 = x253*x505
      x507 = x259*x264
      x508 = x505*x507
      x509 = 0.112523287330675d0*x0
      x510 = x244*x509
      x511 = 0.0162903215389562d0*x419
      x512 = x241*x511
      x513 = 0.112523287330675d0*x379
      x514 = x240*x379
      x515 = 0.0162903215389562d0*x4
      x516 = x245*x502
      x517 = x248*x516/x250**2
      x518 = x258*x509
      x519 = x255*x511
      x520 = x254*x379
      x521 = x239*x259
      x522 = x261*x521/x262**2
      x523 = x247*x251
      x524 = x260*x263
      x525 = 22.9698720015055d0*x238
      x526 = -x252*x502*(x244*x499 - x500*x501) + x253*x514*x525 - x256*
     &x265*(x258* x499 + x500*x503) + x379*x516*x523 - x379*x521*x524 - 
     &x506 - x507 *x520*x525 - x508 - x517*(-x241*x514*x515 - x501*x513 
     &+ x510 + x512) - x522*(x255*x515*x520 + x503*x513 + x518 + x519)
      x527 = 0.167525641071009d0*x200
      x528 = x11**1.66666666666667d0
      x529 = 1.33333333333333d0*x379
      x530 = x13**1.66666666666667d0
      x531 = x0/x1
      x532 = 46.773325076649d0*x238*x281*x282
      x533 = x196*(0.000410930386996015d0*x0**3.0d0*x205 + x197*x207*x47
     &7 - 0.00440195901906306d0*x205*x475 + x220*(-x0*x230 + 0.006932200
     &36674015d0*x204*x476 + x229*x477 - 0.000220097950953153d0*x267*x47
     &5 + x268*(-2.41662179562687d0*x0** 0.333333333333333d0*x232*x233*x
     &235 - x0*x237 - x231*x233*x236*( 1.02357830998401d0*x417 - 0.25655
     &6487706084d0*x419) - x373*x498 + x526) - x473*x497) - 0.0084294123
     &3263362d0*x226*x227*x483 + x228* x456 + x269*(0.0027585187352376d0
     &*x272*x475 - x273*x496 - 0.106650135484359d0*x274*x483 + x275*x456
     & - x493*x494) + x276*(x0 *x278 - x277*x496 - 0.0003448148419047d0*
     &x285*x475 + x286*(-x235* x282*x284*(0.0802319968229891d0*x417 - 0.
     &173432185689313d0*x419) + 11.4849360007528d0*x279*x504 + x280*(-x5
     &28*x529 + x529*x530) - 15.591108358883d0*x283*x504 - 2.64527879654
     &103d0*x283*x531 - x373 *x532 + x526) + x287*x456 - 0.4497842771631
     &83d0*x476*x66 - x493* x527) - x473*x474 + x482*(x479 - x481) + x49
     &2*(x219*( -0.957035512114487d0*x222*x418 - x484 - 1.54043890936153
     &d0*x485 + x487 + x488*(-x479 + x481) + 4.62131672808458d0*x490) + 
     &x491*( x484 + 0.71309495712383d0*x485 - x487 - 2.13928487137149d0*
     &x490 )))
      x534 = 0.632303919304059d0*x0*x194*x288/x195**5
      x535 = x456 - x533 - x534
      x536 = x295*x535
      x537 = 2.78d0*x472 + 2.78d0*x536
      x538 = x66/x296
      x539 = x316*x538
      x540 = 0.0395258831507799d0*x305*x389
      x541 = x302*x306
      x542 = x469*x541
      x543 = x291*(x290*x405 + x38*(x433 + x533 + x534))
      x544 = gamma_cc*x293
      x545 = 0.000526665863310429d0*x313*x544
      x546 = x543*x545
      x547 = x308*x538
      x548 = x537*x547
      x549 = gamma_cc*x311
      x550 = x314*x549
      x551 = x293*x300
      x552 = 1.78430857480394d-5*x299*x52/x292**3
      x553 = x301*x538
      x554 = gamma_cc*x309*x314/x310**2
      x555 = x63/x317
      x556 = x318*x448 + x535 + x555*(-1.27130925834002d0*x304*x312*x389
     & + x469*x470 + x537*x539 + x550*(-x540 + x542 - x546 + x548) + x55
     &4*( 0.00133911037615638d0*x294*x298*x377 - x440*x551 + x540 - x542
     & + x543*x552 + x546 - x548 - x553*(5.56d0*x472 + 5.56d0*x536)))
      x557 = x345*(-x341 - x342 - x343 - x344)
      x558 = x125*(x351 + x352*x557)
      x559 = x355*x558
      x560 = x130*(-x339 + x340*x557 + x368*(x348 + x350 + x361*(x365 - 
     &x559) + x366* (x360 - x362 - x363 + x364*x558 + x559))) - x336 + x
     &369
      x561 = x167*(x335 + x560)
      x562 = rho_s**3*x16
      x563 = x165*(13.56d0*rho_s**5*x151 + 1.74d0*x379 + 2.0d0*x562)
      x564 = 1.33333333333333d0*x0
      x565 = x383*(-x376*(2.0d0*gamma_cc*rho_s - rho_c*x160) - x382*(-x3
     &78*x564 + x381*x564))
      x566 = x563 + x565
      x567 = x170*x566
      x568 = x0*x393
      x569 = x0*x395
      x570 = x392*(-1.33333333333333d0*x568 + 1.33333333333333d0*x569)
      x571 = gamma_cc*x570
      x572 = x402*x570
      x573 = x404*(-1.0d0*x568 + x569)
      x574 = 2.56488140124205d0*x0
      x575 = x408*x574 - x410*x574
      x576 = x23*x575
      x577 = x414*x562
      x578 = x31*x575
      x579 = x429*x562
      x580 = -x578 - x579
      x581 = x407*x576 - x413*x577 + x580
      x582 = x39*(x33*x573 + x38*x581)
      x583 = x436*x582
      x584 = x439*(-2.66666666666667d0*x568 + 2.66666666666667d0*x569)
      x585 = x445*(x399*x571 + x438*(x572 - x583) + x444*(-x441*x584 + x
     &443*x582 - x572 + x583))
      x586 = x447*(x568 - 1.0d0*x569)
      x587 = x586*x62
      x588 = x451*x577
      x589 = x450*x576
      x590 = x578 + x579 - x588 + x589
      x591 = x585 + x587 + x590
      x592 = x155*x563
      x593 = x155*x565
      x594 = x178*x591 + x67*(x175*x592 + x175*x593)
      x595 = x464*x591 + x67*(x463*x592 + x463*x593)
      x596 = x180*x466
      x597 = x581 - x585 - x587
      x598 = x303*x571
      x599 = x471*(x580 + x588 - x589)
      x600 = x203*x379
      x601 = x0*x480
      x602 = x0*x478
      x603 = x271*x379
      x604 = x486*(-0.333333333333333d0*x601 + 0.333333333333333d0*x602)
     &
      x605 = x489*(-0.666666666666667d0*x601 + 0.666666666666667d0*x602)
     &
      x606 = 34.4548080022583d0*x504
      x607 = 0.203363546995945d0*x531
      x608 = -x240*x506 - x240*x252*x607/x243 - x242*x245*x523*x606 + x2
     &54*x508 + x254*x264*x607/x257 + x259*x524*x606 - x517*(x240*x510 +
     & x240* x512) - x522*(-x254*x518 - x254*x519)
      x609 = 1.33333333333333d0*x0
      x610 = x196*(x220*(x268*(x379*x498 + x608) + x497*x600) + x228*x59
     &0 + x269*( x275*x590 + x494*x603) + x276*(x286*(x280*(x528*x609 - 
     &x530*x609 ) + x379*x532 + x608) + x287*x590 + x527*x603) + x474*x6
     &00 + x482 *(x601 - x602) + x492*(x219*(x488*(-x601 + x602) + x604 
     &+ 4.62131672808458d0*x605) + x491*(-x604 - 2.13928487137149d0*x605
     & )))
      x611 = x590 - x610
      x612 = x295*x611
      x613 = 2.78d0*x599 + 2.78d0*x612
      x614 = x541*x598
      x615 = x291*(x290*x573 + x38*(x581 + x610))
      x616 = x545*x615
      x617 = x547*x613
      x618 = x318*x586 + x555*(x470*x598 + x539*x613 + x550*(x614 - x616
     & + x617) + x554*(-x551*x584 + x552*x615 - x553*(5.56d0*x599 + 5.56
     &d0*x612) - x614 + x616 - x617)) + x611
      x619 = x445*(0.296858062715705d0*x401*x51*x54*x57 + x444*(-x435*x4
     &42*x51 - x47 ) + x60)
      x620 = -x619
      x621 = x134*x143*x358
      x622 = 0.25d0*gamma_cs
      x623 = 0.125d0*gamma_cc + 0.125d0*gamma_ss
      x624 = x127*x135
      x625 = 0.167230224308322d0*x146 + x366*(-0.167230224308322d0*x140 
     &- x624*(-x622 + x623)) + 0.111863791688847d0*x621
      x626 = x130*x368
      x627 = x625*x626
      x628 = x107*x327*x98
      x629 = x87*x99
      x630 = 0.167230224308322d0*x110 + x332*(-0.167230224308322d0*x104 
     &- x629*(x622 + x623)) + 0.111863791688847d0*x628
      x631 = x334*x92
      x632 = x630*x631
      x633 = x167*(x627 + x632)
      x634 = -x176*x633
      x635 = gamma_cc*x173
      x636 = x171*x172
      x637 = x635*x636
      x638 = x159*x174*x383
      x639 = 0.00081628665578855d0*x638
      x640 = rho_c**(-7.66666666666667d0)*x156
      x641 = x170*x640
      x642 = 0.03125d0*x177
      x643 = rho_c**(-6.66666666666667d0)
      x644 = x156*x643
      x645 = x178*x619 + x67*(x635*x642 - x639*x644)
      x646 = x186*x50
      x647 = 0.005859375d0*x172
      x648 = x184*x647
      x649 = -x182*x633
      x650 = 0.00228560263620794d0*x638
      x651 = 0.0875d0*x635
      x652 = x464*x619 + x67*(x177*x651 - x644*x650)
      x653 = x299*x51
      x654 = x555*(0.00922952226005862d0*x302*x549*x653 + x315 + x554*(-
     &x307 - 0.000573904446924163d0*x544*x653))
      x655 = x170*x377
      x656 = x37*x62
      x657 = 0.0553050111671826d0*x25 + 2.8d0*x32 + 0.0870539344350337d0
     &*x656 - 0.174107868870067d0*x8
      x658 = x54*x639
      x659 = -x657*x658
      x660 = -x658*x67
      x661 = rho_c*x150
      x662 = x170*x643
      x663 = rho_s*x638
      x664 = 0.00457120527241588d0*x663
      x665 = 0.0016325733115771d0*x663
      x666 = x377*x665
      x667 = x657*x666
      x668 = -0.334460448616645d0*x146 + x366*(0.334460448616645d0*x140 
     &- x624*(x79 - x80 - x81)) - 0.223727583377694d0*x621
      x669 = 0.334460448616645d0*x110 + x332*(-0.334460448616645d0*x104 
     &- x629*x83) + 0.223727583377694d0*x628
      x670 = x167*(x626*x668 + x631*x669)
      x671 = x666*x67
      x672 = x636*x646
      x673 = x642*x646
      x674 = -x67*x673
      x675 = -x657*x673
      x676 = x185/tau_c**4
      x677 = x130*x67
      x678 = x169 + x677
      x679 = x167*x678
      x680 = -x176*x679 + x179
      x681 = -x182*x679 + x183
      x682 = x466*x681 + 1.0d0
      x683 = x150*(x191 + x680*x682)
      x684 = x319 + x683
      x685 = x130*x457
      x686 = 0.0155453454348275d0*x656
      x687 = 0.0310906908696549d0*x8
      x688 = x153*x428
      x689 = 0.0098758948512826d0*x25
      x690 = x686 - x687 + x688 + x689
      x691 = x167*(x335 + x685 + x690)
      x692 = x385*x678
      x693 = x387*x679
      x694 = x187*x462
      x695 = x130*x591
      x696 = x167*(x335 - x686 + x687 - x688 - x689 + x695)
      x697 = x566*x678
      x698 = x466*x680
      x699 = x172*x679
      x700 = 0.03125d0*x699
      x701 = x130*x619
      x702 = x167*(x632 + x701)
      x703 = x640*x678
      x704 = x647*x681
      x705 = x333*x92
      x706 = x630*x705
      x707 = x167*x181
      x708 = 0.000680108862773701d0*x707
      x709 = x377*x678
      x710 = 0.000242896022419179d0*x707
      x711 = x669*x705
      x712 = x643*x678
      x713 = 0.0875d0*x646
      x714 = x67*x92
      x715 = x168 + x714
      x716 = x167*x715
      x717 = -x176*x716 + x179
      x718 = -x182*x716 + x183
      x719 = x466*x718 + 1.0d0
      x720 = x150*(x191 + x717*x719)
      x721 = x319 + x720
      x722 = x457*x92
      x723 = x167*(x370 + x690 + x722)
      x724 = x385*x715
      x725 = x387*x716
      x726 = x591*x92
      x727 = x167*(x560 + x690 + x726)
      x728 = x566*x715
      x729 = x466*x717
      x730 = x172*x716
      x731 = 0.03125d0*x730
      x732 = x619*x92
      x733 = x167*(x627 + x732)
      x734 = x640*x715
      x735 = x647*x718
      x736 = x130*x367
      x737 = x625*x736
      x738 = x377*x715
      x739 = x668*x736
      x740 = x643*x715
      x741 = x677 + x714
      x742 = x167*x741
      x743 = -x176*x742 + x179
      x744 = -x182*x742 + x183
      x745 = x466*x744 + 1.0d0
      x746 = x150*(x191 + x743*x745)
      x747 = x319 + x746
      x748 = x167*(1.0d0*x32 + x64 + x65 + x685 + x722)
      x749 = x385*x741
      x750 = x387*x742
      x751 = x167*(x695 + x726)
      x752 = x566*x741
      x753 = x466*x743
      x754 = x172*x742
      x755 = 0.03125d0*x754
      x756 = x167*(x701 + x732)
      x757 = x640*x741
      x758 = x647*x744
      x759 = x377*x741
      x760 = x643*x741
      PBE_value = x67
      PBE_beta_value = x113
      PBE_alpha_value = x149
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*x320
      d1E(1) = rho_c*(x150*(x180*(-x188*x462 + x466*(-x182*x371 - x182*x
     &386 + 0.13125d0 *x388 + x465)) + x190*(-x176*x371 - x176*x386 + 0.
     &046875d0*x388 + x461) + x467) + x192*x468 + x556) + x320
      d1E(2) = rho_c*(x150*(x190*(-x176*x561 - x176*x567 + x594) + x596*
     &(-x182*x561 - x182*x567 + x595) + x597) + x618)
      d1E(3) = rho_c*(x150*(x180*(x466*(-0.0875d0*x637 + x641*x650 + x64
     &9 + x652) + x646*x648) + x190*(x634 - 0.03125d0*x637 + x639*x641 +
     & x645) + x620) + x654)
      d1E(4) = x661*(x190*(x634 + x639*x655 + x660) + x596*(x649 + x650*
     &x655 + x659))
      d1E(5) = x661*(x190*(-x176*x670 - x662*x665 + x671) + x596*(-x182*
     &x670 - x662* x664 + x667))
      d1E(6) = x661*(x180*(x466*(0.0875d0*x672 + x675) - x648*x676) + x1
     &90*(0.03125d0* x672 + x674))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*x684
      d1E(1) = rho_c*(x150*(x467 + x680*(x466*(-x182*x691 - x182*x692 + 
     &x465 + 0.13125d0*x693) - x681*x694) + x682*(-x176*x691 - x176*x692
     & + x461 + 0.046875d0*x693)) + x468*x683 + x556) + x684
      d1E(2) = rho_c*(x150*(x597 + x682*(-x176*x696 - x176*x697 + x594) 
     &+ x698*(-x182* x696 - x182*x697 + x595)) + x618)
      d1E(3) = rho_c*(x150*(x620 + x680*(x466*(-x182*x702 + x650*x703 - 
     &x651*x699 + x652) + x646*x704) + x682*(-x176*x702 - x635*x700 + x6
     &39*x703 + x645)) + x654)
      d1E(4) = x661*(x682*(x639*x709 + x660 - x706*x710) + x698*(x650*x7
     &09 + x659 - x706*x708))
      d1E(5) = x661*(x682*(-x665*x712 + x671 - x710*x711) + x698*(-x664*
     &x712 + x667 - x708*x711))
      d1E(6) = x661*(x680*(x466*(x675 + x699*x713) - x676*x704) + x682*(
     &x646*x700 + x674))
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*x721
      d1E(1) = rho_c*(x150*(x467 + x717*(x466*(-x182*x723 - x182*x724 + 
     &x465 + 0.13125d0*x725) - x694*x718) + x719*(-x176*x723 - x176*x724
     & + x461 + 0.046875d0*x725)) + x468*x720 + x556) + x721
      d1E(2) = rho_c*(x150*(x597 + x719*(-x176*x727 - x176*x728 + x594) 
     &+ x729*(-x182* x727 - x182*x728 + x595)) + x618)
      d1E(3) = rho_c*(x150*(x620 + x717*(x466*(-x182*x733 + x650*x734 - 
     &x651*x730 + x652) + x646*x735) + x719*(-x176*x733 - x635*x731 + x6
     &39*x734 + x645)) + x654)
      d1E(4) = x661*(x719*(x639*x738 + x660 - x710*x737) + x729*(x650*x7
     &38 + x659 - x708*x737))
      d1E(5) = x661*(x719*(-x665*x740 + x671 - x710*x739) + x729*(-x664*
     &x740 + x667 - x708*x739))
      d1E(6) = x661*(x717*(x466*(x675 + x713*x730) - x676*x735) + x719*(
     &x646*x731 + x674))
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*x747
      d1E(1) = rho_c*(x150*(x467 + x743*(x466*(-x182*x748 - x182*x749 + 
     &x465 + 0.13125d0*x750) - x694*x744) + x745*(-x176*x748 - x176*x749
     & + x461 + 0.046875d0*x750)) + x468*x746 + x556) + x747
      d1E(2) = rho_c*(x150*(x597 + x745*(-x176*x751 - x176*x752 + x594) 
     &+ x753*(-x182* x751 - x182*x752 + x595)) + x618)
      d1E(3) = rho_c*(x150*(x620 + x743*(x466*(-x182*x756 + x650*x757 - 
     &x651*x754 + x652) + x646*x758) + x745*(-x176*x756 - x635*x755 + x6
     &39*x757 + x645)) + x654)
      d1E(4) = x661*(x745*(x639*x759 + x660) + x753*(x650*x759 + x659))
      d1E(5) = x661*(x745*(-x665*x760 + x671) + x753*(-x664*x760 + x667)
     &)
      d1E(6) = x661*(x743*(x466*(x675 + x713*x754) - x676*x758) + x745*(
     &x646*x755 + x674))
      end if
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_TPSS_GWS_ERF(rho_c, rho_s, gamma_cc, g
     &amma_ss, gamma_cs, tau_c, mu, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: November 27, 2018
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
     &57, x2958, x2959, x2960, x2961, x2962, x2963, x2964
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 1d0/rho_c
      x1 = rho_s*x0
      x2 = x1 + 1.0d0
      x3 = x2**0.666666666666667d0
      x4 = -x1 + 1.0d0
      x5 = x4**0.666666666666667d0
      x6 = x3 + x5
      x7 = x6**(-2)
      x8 = x0**0.333333333333333d0
      x9 = x0**0.166666666666667d0
      x10 = sqrt(x0)
      x11 = 0.800428634999364d0*x10
      x12 = x0**0.666666666666667d0
      x13 = x11 + 0.189700432574756d0*x12 + 2.22556942115069d0*x8 + 5.98
     &255043577108d0*x9
      x14 = 1.0d0 + 16.0819842214574d0/x13
      x15 = log(x14)
      x16 = x15*(0.132568899905202d0*x8 + 1.0d0)
      x17 = 0.0621813817393098d0*x16
      x18 = x2**1.33333333333333d0
      x19 = x4**1.33333333333333d0
      x20 = 1.92366105093154d0*x18 + 1.92366105093154d0*x19 - 3.84732210
     &186307d0
      x21 = rho_c**(-4)
      x22 = rho_s**4
      x23 = x21*x22
      x24 = -x23 + 1.0d0
      x25 = 0.430097247127664d0*x10
      x26 = 0.191151259512734d0*x12 + x25 + 2.24759186357762d0*x8 + 8.15
     &741470348764d0*x9
      x27 = 1.0d0 + 29.6085746432167d0/x26
      x28 = log(x27)
      x29 = x24*x28
      x30 = x20*x29
      x31 = x30*(0.0690139921125583d0*x8 + 1.0d0)
      x32 = 0.0197517897025652d0*x31
      x33 = 1.64473377556761d0*x10
      x34 = 0.240587129128819d0*x12 + x33 + 3.84474623744721d0*x8 + 11.1
     &203748630947d0*x9
      x35 = 1.0d0 + 32.1639684429149d0/x34
      x36 = log(x35)
      x37 = x17 - 0.0310906908696548d0*x36*(0.127469618870009d0*x8 + 1.0
     &d0)
      x38 = x23*x37
      x39 = x20*x38
      x40 = x17 - x32 - x39
      x41 = 0.5d0*x3 + 0.5d0*x5
      x42 = x41**3
      x43 = 1d0/x42
      x44 = 32.1639684429148d0*x43
      x45 = exp(x40*x44)
      x46 = x45 - 1.0d0
      x47 = x46**(-2)
      x48 = x6**(-4)
      x49 = x47*x48
      x50 = rho_c**(-4.66666666666667d0)
      x51 = gamma_cc**2
      x52 = x50*x51
      x53 = 0.296858062715706d0*x52
      x54 = 1d0/x46
      x55 = x54*x7
      x56 = rho_c**(-2.33333333333333d0)
      x57 = 0.544846825002868d0*x56
      x58 = x55*x57
      x59 = gamma_cc*x58 + 1.0d0
      x60 = x49*x53 + x59
      x61 = 1d0/x60
      x62 = x61*x7
      x63 = x57*x59
      x64 = x62*x63
      x65 = gamma_cc*x64 + 1.0d0
      x66 = log(x65)
      x67 = 0.0310906908696549d0*x42
      x68 = x66*x67
      x69 = -x17 + x32
      x70 = x39 + x69
      x71 = x68 + x70
      x72 = 6.28318530717959d0*rho_c
      x73 = 6.28318530717959d0*rho_s
      x74 = x72 + x73
      x75 = 1d0/x74
      x76 = x75**0.333333333333333d0
      x77 = x75**0.166666666666667d0
      x78 = sqrt(x75)
      x79 = x75**0.666666666666667d0
      x80 = 8.93863016189422d0*x76 + 16.9559087764345d0*x77 + 5.83042942
     &843835d0*x78 + 1.30040600365736d0*x79
      x81 = 1.0d0 + 32.1639684429149d0/x80
      x82 = log(x81)
      x83 = x82*(0.296353441706766d0*x76 + 1.0d0)
      x84 = 0.5d0*gamma_cs
      x85 = 0.25d0*gamma_cc
      x86 = 0.25d0*gamma_ss
      x87 = x85 + x86
      x88 = x84 + x87
      x89 = exp(1.999999999999994d0*x83)
      x90 = x89 - 1.0d0
      x91 = 1d0/x90
      x92 = 0.5d0*rho_c
      x93 = 0.5d0*rho_s
      x94 = x92 + x93
      x95 = x94**(-2)
      x96 = x91*x95
      x97 = 14.804406601634d0*rho_c
      x98 = 14.804406601634d0*rho_s
      x99 = x97 + x98
      x100 = x99**(-0.333333333333333d0)
      x101 = 0.668920897233289d0*x100
      x102 = x101*x96
      x103 = x102*x88 + 1.0d0
      x104 = x90**(-2)
      x105 = x88**2
      x106 = x104*x105
      x107 = x94**(-4)
      x108 = x99**(-0.666666666666667d0)
      x109 = x107*x108
      x110 = 0.447455166755389d0*x109
      x111 = x103 + x106*x110
      x112 = 1d0/x111
      x113 = x112*x95
      x114 = x103*x113
      x115 = x101*x88
      x116 = x114*x115 + 1.0d0
      x117 = log(x116)
      x118 = 0.0155453454348274d0*x117 - 0.0310906908696548d0*x83
      x119 = x72 - x73
      x120 = 1d0/x119
      x121 = x120**0.333333333333333d0
      x122 = x120**0.166666666666667d0
      x123 = sqrt(x120)
      x124 = x120**0.666666666666667d0
      x125 = 8.93863016189422d0*x121 + 16.9559087764345d0*x122 + 5.83042
     &942843835d0* x123 + 1.30040600365736d0*x124
      x126 = 1.0d0 + 32.1639684429149d0/x125
      x127 = log(x126)
      x128 = x127*(0.296353441706766d0*x121 + 1.0d0)
      x129 = x92 - x93
      x130 = x129**(-2)
      x131 = exp(1.999999999999994d0*x128)
      x132 = x131 - 1.0d0
      x133 = x132**(-2)
      x134 = -x84 + x87
      x135 = x134**2
      x136 = x133*x135
      x137 = x129**(-4)
      x138 = x97 - x98
      x139 = x138**(-0.666666666666667d0)
      x140 = x137*x139
      x141 = 0.447455166755389d0*x140
      x142 = 1d0/x132
      x143 = x130*x142
      x144 = x138**(-0.333333333333333d0)
      x145 = 0.668920897233289d0*x144
      x146 = x143*x145
      x147 = x134*x146 + 1.0d0
      x148 = x136*x141 + x147
      x149 = 1d0/x148
      x150 = x147*x149
      x151 = x130*x150
      x152 = x134*x145
      x153 = x151*x152 + 1.0d0
      x154 = log(x153)
      x155 = -0.0310906908696548d0*x128 + 0.0155453454348274d0*x154
      x156 = exp(-0.46869933305458689d0*mu*rho_c**(-0.33333333333333331d
     &0))
      x157 = x129*x155
      x158 = x118*x94
      x159 = x157 + x158
      x160 = rho_c**(-3)
      x161 = rho_c**(-6)
      x162 = rho_s**6
      x163 = 0.5d0*x23
      x164 = rho_c**2
      x165 = 1d0/x164
      x166 = rho_s**2
      x167 = x165*x166
      x168 = 2.26d0*x161*x162 + x163 + 0.87d0*x167 + 0.53d0
      x169 = 1d0/x19 + 1d0/x18
      x170 = gamma_cc*x166
      x171 = 2.0d0*gamma_cs
      x172 = rho_s*x171
      x173 = gamma_ss*x164 - rho_c*x172 + x170
      x174 = x169*x173
      x175 = 0.0130605864926168d0*x174*x50 + 1.0d0
      x176 = x175**(-4)
      x177 = x168*x176
      x178 = x177 + 1.0d0
      x179 = x160*x178
      x180 = x159*x179
      x181 = tau_c**(-2)
      x182 = x181*x51
      x183 = 0.015625d0*x182
      x184 = x165*x177
      x185 = x183*x184 + 1.0d0
      x186 = x185*x71
      x187 = -x180*x183 + x186
      x188 = 0.04375d0*x182
      x189 = 2.8d0*x186
      x190 = -x180*x188 + x189
      x191 = x160*x190
      x192 = gamma_cc**3
      x193 = tau_c**(-3)
      x194 = x192*x193
      x195 = 0.001953125d0*x194
      x196 = x191*x195 + 1.0d0
      x197 = x40 - x68
      x198 = x156*(x187*x196 + x197)
      x199 = x48*x52
      x200 = 1d0/x70
      x201 = mu**2
      x202 = x12*x201
      x203 = 0.237113969739022d0*x202 + 1.0d0
      x204 = x203**(-4)
      x205 = mu**5
      x206 = x167 - 1.0d0
      x207 = exp(-0.46675853320810856d0*x8)
      x208 = x0**1.0d0
      x209 = 0.00304892328706429d0*x208
      x210 = x0**1.33333333333333d0
      x211 = x207*(0.0315297404574782d0*x12 - x209 + 0.00027531079361635
     &8d0*x210 + 0.0128480790170175d0*x8 + 1.0d0)
      x212 = x206*x211
      x213 = x205*x212
      x214 = x0**1.66666666666667d0
      x215 = 0.00264117541143784d0*x214
      x216 = x2**0.666666666666667d0
      x217 = x4**0.666666666666667d0
      x218 = 0.5d0*x216 + 0.5d0*x217
      x219 = x218**3
      x220 = -0.0621813817393098d0
      x221 = x218**2
      x222 = 1d0/x221
      x223 = x201*x8
      x224 = x222*x223
      x225 = mu/x218
      x226 = 4.60448529760779d0*x9
      x227 = x225*x226 + 1.0d0
      x228 = 2.13928487137149d0*x224 + x227
      x229 = 1d0/x228
      x230 = mu**3
      x231 = 1d0/x219
      x232 = x230*x231
      x233 = 1.91407102422897d0*x10
      x234 = 4.62131672808458d0*x224 + x227 + x232*x233
      x235 = x220*log(x229*x234)
      x236 = mu**8*x0**2.66666666666667d0
      x237 = 0.00316102962473761d0*x236
      x238 = 0.0445553741830535d0*x208
      x239 = x212*x238
      x240 = 0.384834731559127d0*x12 - 3.07073492995203d0*x8
      x241 = exp(-0.19230865217881402d0*x8)
      x242 = x241/x208
      x243 = x240*x242
      x244 = -x167 + 1.0d0
      x245 = 37.6991118430775d0*x244
      x246 = x243*x245
      x247 = 1d0/x12
      x248 = 34.4548080022583d0*x247
      x249 = 0.5d0*x1
      x250 = x249 + 0.5d0
      x251 = x250**2
      x252 = 1d0/x2
      x253 = x252**0.333333333333333d0
      x254 = x253*x8
      x255 = x252**0.666666666666667d0
      x256 = 0.0244354823084343d0*x12
      x257 = 0.337569861992026d0*x254 + x255*x256 + 1.0d0
      x258 = 1d0/x257
      x259 = x251*x258
      x260 = 1d0/x255
      x261 = x260*(-0.0177069812999059d0*x254 + 1.0d0)
      x262 = x259*x261
      x263 = 1d0/x4
      x264 = x263**0.666666666666667d0
      x265 = 1d0/x264
      x266 = x263**0.333333333333333d0
      x267 = x266*x8
      x268 = x265*(-0.0177069812999059d0*x267 + 1.0d0)
      x269 = -x249 + 0.5d0
      x270 = x269**2
      x271 = x256*x264 + 0.337569861992026d0*x267 + 1.0d0
      x272 = 1d0/x271
      x273 = x270*x272
      x274 = x248*x273
      x275 = -x248*x262 - x268*x274
      x276 = -x246 + x275
      x277 = 0.000132058770571892d0*x214
      x278 = mu**6
      x279 = -3.0d0*x167 + 3.0d0
      x280 = 0.5d0*x211 - 0.5d0
      x281 = x279*x280
      x282 = 0.00165511124114256d0*x214
      x283 = x0**2.0d0
      x284 = 0.0533250677421794d0*x283
      x285 = mu**4
      x286 = 0.0279209401785015d0*x208
      x287 = x281*x286
      x288 = 0.5d0*x2**2.66666666666667d0 + 0.5d0*x4**2.66666666666667d0
     &
      x289 = 17.2274040011292d0*x247
      x290 = exp(-0.33933171852197186d0*x8)
      x291 = x290*(0.26014827853397d0*x12 - 0.240695990468967d0*x8)
      x292 = x244*x291
      x293 = 23.3866625383245d0*x247
      x294 = x275 + x288*x289 - x292*x293
      x295 = 0.00020688890514282d0*x214
      x296 = 0.337338207872387d0*x210
      x297 = x213*x215 + x219*x235 + x230*(x239 + x276*x277) + x237*x70 
     &+ x278*(-x281 *x282 + x284*x70) + x285*(-x287 + x294*x295 + x296*x
     &70)
      x298 = x204*x297
      x299 = -x298 + x70
      x300 = x200*x299
      x301 = x300**5.56d0
      x302 = x298 + x40
      x303 = exp(x302*x44)
      x304 = 0.0310906908696549d0*x303 - 0.0310906908696549d0
      x305 = x304**(-2)
      x306 = x301*x305
      x307 = 0.000286952223462081d0*x306
      x308 = x199*x307
      x309 = 1d0/x304
      x310 = x300**2.78d0
      x311 = x310*x7
      x312 = x309*x311
      x313 = x312*x56
      x314 = 0.0169396642074771d0*x313
      x315 = gamma_cc*x314
      x316 = x315 + 1.0d0
      x317 = x308 + x316
      x318 = 1d0/x317
      x319 = x311*x318
      x320 = x316*x56
      x321 = 0.544846825002868d0*x320
      x322 = x319*x321
      x323 = gamma_cc*x322
      x324 = x323 + 1.0d0
      x325 = log(x324)
      x326 = x299 + x325*x67
      x327 = x198 + x326
      x328 = -0.00921383324426333d0*x76 - 0.0310906908696548d0
      x329 = x75*x77
      x330 = x75*x76
      x331 = x75*x78
      x332 = x75*x79
      x333 = 17.7561861489951d0*x329 + 18.7210232331754d0*x330 + 18.3168
     &342596557d0* x331 + 5.44712793036536d0*x332
      x334 = 1d0/x81
      x335 = x80**(-2)
      x336 = 32.1639684429149d0*x334*x335
      x337 = x333*x336
      x338 = x94**(-3)
      x339 = x103*x112
      x340 = x338*x339
      x341 = x115*x340
      x342 = x99**(-1.33333333333333d0)
      x343 = x114*x342
      x344 = 3.30099231565716d0*x88
      x345 = x343*x344
      x346 = x338*x91
      x347 = x346*x88
      x348 = x101*x347
      x349 = x342*x96
      x350 = x349*x88
      x351 = 3.30099231565716d0*x350
      x352 = x330*x82
      x353 = 0.592706883413531d0*x76 + 1.99999999999999d0
      x354 = x337*x353 - 1.24136239377603d0*x352
      x355 = x354*x89
      x356 = x104*x355
      x357 = x115*x95
      x358 = x356*x357
      x359 = -x348 - x351 - x358
      x360 = x113*x359
      x361 = x115*x360
      x362 = x111**(-2)
      x363 = x94**(-5)
      x364 = x106*x108
      x365 = x363*x364
      x366 = x99**(-1.66666666666667d0)
      x367 = x106*x107
      x368 = x366*x367
      x369 = x105*x109
      x370 = x90**(-3)
      x371 = x355*x370
      x372 = 0.894910333510778d0*x371
      x373 = x348 + x351 + x358 + 0.894910333510778d0*x365 + 4.416205483
     &09916d0*x368 + x369*x372
      x374 = x362*x373
      x375 = x103*x374
      x376 = x357*x375
      x377 = -x341 - x345 + x361 + x376
      x378 = 1d0/x116
      x379 = 0.0155453454348274d0*x378
      x380 = 0.0192974072210527d0*x352 + x377*x379
      x381 = 0.00777267271741372d0*x117 - 0.0155453454348274d0*x83 + x94
     &*(x328*x337 + x380)
      x382 = 0.00777267271741372d0*x154
      x383 = x120*x122
      x384 = 17.7561861489951d0*x383
      x385 = x120*x121
      x386 = 18.7210232331754d0*x385
      x387 = x120*x123
      x388 = 18.3168342596557d0*x387
      x389 = x120*x124
      x390 = 5.44712793036536d0*x389
      x391 = x384 + x386 + x388 + x390
      x392 = -0.00921383324426333d0*x121 - 0.0310906908696548d0
      x393 = 1d0/x126
      x394 = x393/x125**2
      x395 = 32.1639684429149d0*x394
      x396 = x392*x395
      x397 = x130*x152
      x398 = x127*x385
      x399 = 1.24136239377603d0*x398
      x400 = 0.592706883413531d0*x121 + 1.99999999999999d0
      x401 = x391*x400
      x402 = x395*x401 - x399
      x403 = x131*x133
      x404 = x402*x403
      x405 = x397*x404
      x406 = x129**(-3)
      x407 = x142*x406
      x408 = x134*x407
      x409 = x145*x408
      x410 = x138**(-1.33333333333333d0)
      x411 = x143*x410
      x412 = x134*x411
      x413 = 3.30099231565716d0*x412
      x414 = -x409 - x413
      x415 = -x405 + x414
      x416 = x130*x149
      x417 = x152*x416
      x418 = x415*x417
      x419 = x129**(-5)
      x420 = x136*x139
      x421 = x419*x420
      x422 = 0.894910333510778d0*x421
      x423 = x138**(-1.66666666666667d0)
      x424 = x136*x137
      x425 = x423*x424
      x426 = 4.41620548309916d0*x425
      x427 = x132**(-3)
      x428 = x402*x427
      x429 = x131*x428
      x430 = x135*x140
      x431 = 0.894910333510778d0*x430
      x432 = x429*x431
      x433 = x409 + x413
      x434 = x405 + x422 + x426 + x432 + x433
      x435 = x148**(-2)
      x436 = x147*x435
      x437 = x130*x436
      x438 = x152*x437
      x439 = x434*x438
      x440 = x150*x406
      x441 = x152*x440
      x442 = x151*x410
      x443 = 3.30099231565716d0*x134
      x444 = x442*x443
      x445 = -x441 - x444
      x446 = x418 + x439 + x445
      x447 = 1d0/x153
      x448 = 0.0155453454348274d0*x447
      x449 = 0.0192974072210527d0*x398
      x450 = x446*x448 + x449
      x451 = 0.0155453454348274d0*x128
      x452 = x129*(x391*x396 + x450) + x382 - x451
      x453 = x381 + x452
      x454 = x179*x453
      x455 = rho_c**(-7)
      x456 = rho_c**(-5)
      x457 = x22*x456
      x458 = 2.0d0*x457
      x459 = x160*x166
      x460 = -13.56d0*x162*x455 - x458 - 1.74d0*x459
      x461 = x176*x460
      x462 = 0.0522423459704672d0*x50
      x463 = 2.0d0*gamma_ss*rho_c - x172
      x464 = x169*x463
      x465 = rho_c**(-5.66666666666667d0)
      x466 = x174*x465
      x467 = x2**(-2.33333333333333d0)
      x468 = 1.33333333333333d0*x165
      x469 = x467*x468
      x470 = x4**(-2.33333333333333d0)
      x471 = x468*x470
      x472 = rho_s*x469 - rho_s*x471
      x473 = 0.0522423459704672d0*x472
      x474 = x173*x50
      x475 = -x462*x464 + 0.243797614528847d0*x466 - x473*x474
      x476 = x175**(-5)
      x477 = x168*x476
      x478 = x475*x477
      x479 = x461 + x478
      x480 = x159*x160
      x481 = x479*x480
      x482 = x159*x178
      x483 = x182*x21
      x484 = 0.046875d0*x483
      x485 = x2**0.333333333333333d0
      x486 = 2.56488140124205d0*x165
      x487 = x485*x486
      x488 = x4**0.333333333333333d0
      x489 = x486*x488
      x490 = -rho_s*x487 + rho_s*x489
      x491 = x38*x490
      x492 = x0*x8
      x493 = 1d0/x35
      x494 = x34**(-2)
      x495 = x0*x9
      x496 = x0*x10
      x497 = x0*x12
      x498 = 1.28158207914907d0*x492 + 1.85339581051578d0*x495 + 0.82236
     &6887783805d0* x496 + 0.160391419419213d0*x497
      x499 = x493*x494*x498
      x500 = -0.00396311851556016d0*x8 - 0.0310906908696548d0
      x501 = 32.1639684429149d0*x500
      x502 = x15*x492
      x503 = 0.00274777245725523d0*x502
      x504 = -x503
      x505 = 0.0082433173717657d0*x8
      x506 = x505 + 0.0621813817393098d0
      x507 = 0.741856473716896d0*x492 + 0.99709173929518d0*x495 + 0.4002
     &14317499682d0 *x496 + 0.126466955049837d0*x497
      x508 = 1d0/x14
      x509 = x13**(-2)
      x510 = 16.0819842214574d0*x508*x509
      x511 = x507*x510
      x512 = x504 + x506*x511
      x513 = 0.00132103950518672d0*x36*x492 + x499*x501 + x512
      x514 = x20*x23
      x515 = x513*x514
      x516 = 0.00136314985874174d0*x8
      x517 = x516 + 0.0197517897025652d0
      x518 = x29*x490
      x519 = x517*x518
      x520 = x20*x37
      x521 = 4.0d0*x457
      x522 = x520*x521
      x523 = -x505 - 0.0621813817393098d0
      x524 = x511*x523
      x525 = x28*x517
      x526 = x20*x521
      x527 = x525*x526
      x528 = 1d0/x27
      x529 = x26**(-2)
      x530 = 0.749197287859205d0*x492 + 1.35956911724794d0*x495 + 0.2150
     &48623563832d0 *x496 + 0.127434173008489d0*x497
      x531 = x528*x529*x530
      x532 = 29.6085746432167d0*x531
      x533 = x20*x24
      x534 = x517*x533
      x535 = x532*x534
      x536 = 0.000454383286247248d0*x492
      x537 = x30*x536
      x538 = x503 - x537
      x539 = x491 + x515 + x519 - x522 + x524 + x527 + x535 + x538
      x540 = x59*x62
      x541 = rho_c**(-3.33333333333333d0)
      x542 = 1.27130925834002d0*x541
      x543 = x540*x542
      x544 = gamma_cc*x543
      x545 = x2**(-0.333333333333333d0)
      x546 = x165*x545
      x547 = rho_s*x546
      x548 = x4**(-0.333333333333333d0)
      x549 = x165*x548
      x550 = rho_s*x549
      x551 = 1.33333333333333d0*x547 - 1.33333333333333d0*x550
      x552 = x6**(-3)
      x553 = x551*x552
      x554 = x61*x63
      x555 = x553*x554
      x556 = gamma_cc*x555
      x557 = x541*x55
      x558 = 1.27130925834002d0*x557
      x559 = gamma_cc*x558
      x560 = x54*x552
      x561 = x560*x57
      x562 = x551*x561
      x563 = gamma_cc*x562
      x564 = x57*x7
      x565 = x547 - 1.0d0*x550
      x566 = x41**(-4)
      x567 = 32.1639684429148d0*x566
      x568 = x565*x567
      x569 = -x516 - 0.0197517897025652d0
      x570 = x28*x569
      x571 = x533*x569
      x572 = -x491 - x515 + x522 + x537
      x573 = x512 + x518*x569 + x526*x570 + x532*x571 + x572
      x574 = x40*x568 + x44*x573
      x575 = x45*x574
      x576 = x47*x575
      x577 = x564*x576
      x578 = gamma_cc*x577
      x579 = -x559 + x563 - x578
      x580 = x57*x62
      x581 = x579*x580
      x582 = gamma_cc*x581
      x583 = x49*x51
      x584 = x465*x583
      x585 = 2.66666666666667d0*x547 - 2.66666666666667d0*x550
      x586 = x6**(-5)
      x587 = x47*x586
      x588 = x53*x587
      x589 = x46**(-3)
      x590 = x575*x589
      x591 = 0.593716125431411d0*x199
      x592 = x590*x591
      x593 = x559 - x563 + x578 + 1.38533762600663d0*x584 - x585*x588 + 
     &x592
      x594 = x60**(-2)
      x595 = x594*x7
      x596 = x595*x63
      x597 = x593*x596
      x598 = gamma_cc*x597
      x599 = -x544 + x556 + x582 + x598
      x600 = 1d0/x65
      x601 = x600*x67
      x602 = x599*x601
      x603 = -1.0d0*x547 + x550
      x604 = x41**2
      x605 = 0.0310906908696549d0*x604
      x606 = x603*x605
      x607 = x606*x66
      x608 = x602 + x607
      x609 = x539 + x608
      x610 = x165*x461
      x611 = 0.03125d0*x181
      x612 = x160*x51
      x613 = x177*x612
      x614 = x165*x478
      x615 = x183*x610 + x183*x614 - x611*x613
      x616 = x185*x609 + x615*x71
      x617 = -x183*x454 - x183*x481 + x482*x484 + x616
      x618 = x194*x21
      x619 = 0.005859375d0*x618
      x620 = 0.13125d0*x483
      x621 = x184*x188 + 2.8d0
      x622 = 0.0875d0*x181
      x623 = x188*x610 + x188*x614 - x613*x622
      x624 = x609*x621 + x623*x71
      x625 = -x188*x454 - x188*x481 + x482*x620 + x624
      x626 = x160*x195
      x627 = -x190*x619 + x625*x626
      x628 = x573 - x602 - x607
      x629 = x156*(x187*x627 + x196*x617 + x628)
      x630 = mu*rho_c**(-1.33333333333333d0)
      x631 = x198*x630
      x632 = x316*x319
      x633 = 1.27130925834002d0*x541
      x634 = x632*x633
      x635 = gamma_cc*x634
      x636 = x310*x318
      x637 = x321*x636
      x638 = x553*x637
      x639 = gamma_cc*x638
      x640 = x504 - x519 - x524 - x527 - x535 + x572
      x641 = x70**2
      x642 = 1d0/x641
      x643 = x299*x642
      x644 = x640*x643
      x645 = x214*x459
      x646 = 0.00528235082287567d0*x205*x211
      x647 = x0*x214
      x648 = x0**2.0d0
      x649 = x0*x648
      x650 = x0*x210
      x651 = x207*(x0*x209 - 0.00428269300567249d0*x492 - 0.021019826971
     &6521d0*x497 - 0.000367081058155144d0*x650)
      x652 = x206*x651
      x653 = x205*x215
      x654 = x4**(-0.333333333333333d0)
      x655 = x165*x654
      x656 = rho_s*x655
      x657 = x2**(-0.333333333333333d0)
      x658 = x165*x657
      x659 = rho_s*x658
      x660 = x656 - x659
      x661 = x221*x235
      x662 = x236*x70
      x663 = 0.00842941233263362d0*x0
      x664 = 0.666666666666667d0*x659
      x665 = 0.666666666666667d0*x656
      x666 = x664 - x665
      x667 = x223*x231
      x668 = x666*x667
      x669 = 4.62131672808458d0*x668
      x670 = x218**(-4)
      x671 = -x656 + x659
      x672 = x670*x671
      x673 = x230*x233
      x674 = x672*x673
      x675 = x201*x222
      x676 = x492*x675
      x677 = 1.54043890936153d0*x676
      x678 = 0.957035512114487d0*x496
      x679 = x232*x678
      x680 = -0.333333333333333d0*x656 + 0.333333333333333d0*x659
      x681 = mu*x222
      x682 = x226*x681
      x683 = x680*x682
      x684 = x225*x495
      x685 = 0.767414216267966d0*x684
      x686 = x683 - x685
      x687 = x669 + x674 - x677 - x679 + x686
      x688 = 0.71309495712383d0*x676
      x689 = 2.13928487137149d0*x668
      x690 = -x683 + x685
      x691 = x688 - x689 + x690
      x692 = x228**(-2)
      x693 = x234*x692
      x694 = x229*x687 + x691*x693
      x695 = 1d0/x234
      x696 = x219*x220*x695
      x697 = x228*x696
      x698 = 0.00993066744685536d0*x280
      x699 = 0.00152446164353214d0*x208
      x700 = x207*(0.0157648702287391d0*x12 + 0.000137655396808179d0*x21
     &0 - x699 + 0.00642403950850874d0*x8 + 0.5d0)
      x701 = x207*(x0*x699 - 0.00214134650283625d0*x492 - 0.010509913485
     &8261d0*x497 - 0.000183540529077572d0*x650)
      x702 = 0.155586177736036d0*x492*x700 + x701
      x703 = x279*x702
      x704 = x283*x70
      x705 = 0.106650135484359d0*x0
      x706 = x208*x459
      x707 = 0.0891107483661069d0*x211
      x708 = x206*x238
      x709 = 0.000220097950953153d0*x647
      x710 = 1.02357830998401d0*x492 - 0.256556487706084d0*x497
      x711 = x242*x245
      x712 = 75.398223686155d0*x243
      x713 = x0*x244
      x714 = x0**(-0.666666666666667d0)
      x715 = x240*x241
      x716 = x714*x715
      x717 = 0.00590232709996862d0*x0
      x718 = 0.00590232709996862d0*x165
      x719 = rho_s*x252
      x720 = x254*x719
      x721 = x260*(x254*x717 - x718*x720)
      x722 = x248*x259
      x723 = rho_s*x263
      x724 = x267*x723
      x725 = x267*x717 + x718*x724
      x726 = x265*x274
      x727 = x0*x247
      x728 = 22.9698720015055d0*x727
      x729 = x262*x728
      x730 = x268*x273
      x731 = x728*x730
      x732 = 0.112523287330675d0*x0
      x733 = x254*x732
      x734 = 0.0162903215389562d0*x497
      x735 = x255*x734
      x736 = 0.112523287330675d0*x165
      x737 = x12*x165
      x738 = 0.0162903215389562d0*x737
      x739 = x255*x719
      x740 = -x720*x736 + x733 + x735 - x738*x739
      x741 = x251*x740
      x742 = x257**(-2)
      x743 = x261*x742
      x744 = x248*x743
      x745 = x267*x732
      x746 = x264*x734
      x747 = x264*x723
      x748 = x724*x736 + x738*x747 + x745 + x746
      x749 = x270*x748
      x750 = x271**(-2)
      x751 = x268*x750
      x752 = x248*x751
      x753 = x258*x261
      x754 = x250*x753
      x755 = x165*x247
      x756 = rho_s*x755
      x757 = 34.4548080022583d0*x756
      x758 = x268*x272
      x759 = x269*x758
      x760 = 22.9698720015055d0*x755
      x761 = x719*x760
      x762 = x723*x760
      x763 = x262*x761 - x721*x722 - x725*x726 - x729 - x730*x762 - x731
     & - x741*x744 - x749*x752 + x754*x757 - x757*x759
      x764 = -x0*x246 - x459*x712 - x710*x711 - 2.41662179562687d0*x713*
     &x716 + x763
      x765 = 0.167525641071009d0*x280
      x766 = 0.0003448148419047d0*x647
      x767 = 0.449784277163183d0*x650
      x768 = x2**1.66666666666667d0
      x769 = 1.33333333333333d0*x165
      x770 = x768*x769
      x771 = x4**1.66666666666667d0
      x772 = x769*x771
      x773 = -rho_s*x770 + rho_s*x772
      x774 = 11.4849360007528d0*x727
      x775 = x290*(0.0802319968229891d0*x492 - 0.173432185689313d0*x497)
     &
      x776 = x244*x293
      x777 = 1d0/x8
      x778 = x0*x777
      x779 = 46.773325076649d0*x291
      x780 = x247*x459
      x781 = x288*x774 + x289*x773 - 15.591108358883d0*x292*x727 - 2.645
     &27879654103d0 *x292*x778 + x763 - x775*x776 - x779*x780
      x782 = -0.00440195901906306d0*x213*x647 + 0.000410930386996015d0*x
     &213*x649 + x230*(-x0*x239 + 0.00693220036674015d0*x212*x650 - x276
     &*x709 + x277*x764 + x651*x708 - x706*x707) + x237*x539 + x278*( 0.
     &0027585187352376d0*x281*x647 - x282*x703 + x284*x539 - x645* x698 
     &- x704*x705) + x285*(x0*x287 - x286*x703 - x294*x766 + x295* x781 
     &+ x296*x539 - x70*x767 - x706*x765) - x645*x646 + x652*x653 + x660
     &*x661 - x662*x663 + x694*x697
      x783 = x204*x782
      x784 = x203**(-5)
      x785 = x0*x202*x784
      x786 = x297*x785
      x787 = 0.632303919304059d0*x786
      x788 = x539 - x783 - x787
      x789 = x200*x788
      x790 = 2.78d0*x644 + 2.78d0*x789
      x791 = 1d0/x299
      x792 = x70*x791
      x793 = x790*x792
      x794 = x323*x793
      x795 = x312*x541
      x796 = 0.0395258831507799d0*x795
      x797 = gamma_cc*x796
      x798 = x309*x310
      x799 = x553*x798
      x800 = 0.0169396642074771d0*x56
      x801 = x799*x800
      x802 = gamma_cc*x801
      x803 = x783 + x787
      x804 = x573 + x803
      x805 = x302*x568 + x44*x804
      x806 = x303*x805
      x807 = gamma_cc*x56
      x808 = x305*x311
      x809 = 0.000526665863310429d0*x808
      x810 = x807*x809
      x811 = x806*x810
      x812 = x315*x792
      x813 = x790*x812
      x814 = -x797 + x802 - x811 + x813
      x815 = x319*x814
      x816 = 0.544846825002868d0*x56
      x817 = x815*x816
      x818 = gamma_cc*x817
      x819 = x465*x51
      x820 = x306*x48
      x821 = x819*x820
      x822 = 0.00133911037615638d0*x821
      x823 = x585*x586
      x824 = x307*x52
      x825 = x823*x824
      x826 = x304**(-3)
      x827 = x806*x826
      x828 = x199*x301
      x829 = 1.78430857480394d-5*x828
      x830 = x827*x829
      x831 = 5.56d0*x644 + 5.56d0*x789
      x832 = x308*x792
      x833 = x797 - x802 + x811 - x813 + x822 - x825 + x830 - x831*x832
      x834 = x317**(-2)
      x835 = x311*x834
      x836 = x321*x835
      x837 = x833*x836
      x838 = gamma_cc*x837
      x839 = -x635 + x639 + x794 + x818 + x838
      x840 = 1d0/x324
      x841 = x67*x840
      x842 = x325*x606 + x788 + x839*x841
      x843 = -x449
      x844 = -x384 - x386 - x388 - x390
      x845 = x395*x400
      x846 = x399 + x844*x845
      x847 = x397*x403
      x848 = x846*x847
      x849 = x433 - x848
      x850 = x417*x849
      x851 = x131*x427*x431
      x852 = x414 - x422 - x426 + x846*x851 + x848
      x853 = x438*x852
      x854 = x441 + x444
      x855 = x850 + x853 + x854
      x856 = x448*x855
      x857 = x129*(x396*x844 + x843 + x856) - x382 + x451
      x858 = x381 + x857
      x859 = x179*x858
      x860 = 1.74d0*x165
      x861 = rho_s**5
      x862 = rho_s**3
      x863 = x21*x862
      x864 = 2.0d0*x863
      x865 = rho_s*x860 + 13.56d0*x161*x861 + x864
      x866 = x176*x865
      x867 = 2.0d0*rho_s
      x868 = gamma_cc*x867 - rho_c*x171
      x869 = x169*x868
      x870 = 1.33333333333333d0*x0
      x871 = -x467*x870 + x470*x870
      x872 = x173*x462
      x873 = -x462*x869 - x871*x872
      x874 = x477*x873
      x875 = x866 + x874
      x876 = x480*x875
      x877 = 2.56488140124205d0*x0
      x878 = x485*x877 - x488*x877
      x879 = x38*x878
      x880 = 4.0d0*x863
      x881 = x520*x880
      x882 = x29*x878
      x883 = x517*x882
      x884 = x20*x880
      x885 = x525*x884
      x886 = x883 - x885
      x887 = x879 + x881 + x886
      x888 = x0*x545
      x889 = x0*x548
      x890 = -1.33333333333333d0*x888 + 1.33333333333333d0*x889
      x891 = x552*x890
      x892 = x554*x891
      x893 = gamma_cc*x892
      x894 = x561*x890
      x895 = gamma_cc*x894
      x896 = -1.0d0*x888 + x889
      x897 = x567*x896
      x898 = -x879 - x881
      x899 = x569*x882 - x570*x884 + x898
      x900 = x40*x897 + x44*x899
      x901 = x45*x900
      x902 = x47*x564
      x903 = x901*x902
      x904 = gamma_cc*x903
      x905 = x895 - x904
      x906 = x580*x905
      x907 = gamma_cc*x906
      x908 = -2.66666666666667d0*x888 + 2.66666666666667d0*x889
      x909 = x589*x591
      x910 = -x588*x908 - x895 + x901*x909 + x904
      x911 = x596*x910
      x912 = gamma_cc*x911
      x913 = x893 + x907 + x912
      x914 = x601*x913
      x915 = x888 - 1.0d0*x889
      x916 = x605*x915
      x917 = x66*x916
      x918 = x914 + x917
      x919 = x887 + x918
      x920 = x165*x866
      x921 = x165*x874
      x922 = x183*x920 + x183*x921
      x923 = x185*x919 + x71*x922
      x924 = -x183*x859 - x183*x876 + x923
      x925 = x188*x920 + x188*x921
      x926 = x621*x919 + x71*x925
      x927 = -x188*x859 - x188*x876 + x926
      x928 = x187*x626
      x929 = -x914 - x917
      x930 = x899 + x929
      x931 = x156*(x196*x924 + x927*x928 + x930)
      x932 = x637*x891
      x933 = gamma_cc*x932
      x934 = -x883 + x885 + x898
      x935 = x643*x934
      x936 = x165*x214
      x937 = x646*x936
      x938 = x0*x657
      x939 = x0*x654
      x940 = x938 - x939
      x941 = x698*x936
      x942 = -0.333333333333333d0*x938 + 0.333333333333333d0*x939
      x943 = x682*x942
      x944 = -x938 + x939
      x945 = x670*x673
      x946 = x944*x945
      x947 = 0.666666666666667d0*x938
      x948 = 0.666666666666667d0*x939
      x949 = -x947 + x948
      x950 = x667*x949
      x951 = 4.62131672808458d0*x950
      x952 = x943 + x946 + x951
      x953 = -x943
      x954 = 2.13928487137149d0*x950
      x955 = x953 - x954
      x956 = x229*x952 + x693*x955
      x957 = x165*x208
      x958 = x707*x957
      x959 = x165*x712
      x960 = 34.4548080022583d0*x727
      x961 = 0.203363546995945d0*x778
      x962 = 1d0/x266
      x963 = x273*x962
      x964 = x263*x963
      x965 = x252*x733 + x252*x735
      x966 = x251*x744
      x967 = -x263*x745 - x263*x746
      x968 = x270*x752
      x969 = 1d0/x253
      x970 = x259*x969
      x971 = x252*x970
      x972 = -x252*x729 + x263*x731 - x754*x960 + x759*x960 + x961*x964 
     &- x961*x971 - x965*x966 - x967*x968
      x973 = rho_s*x959 + x972
      x974 = x765*x957
      x975 = 1.33333333333333d0*x0
      x976 = x768*x975 - x771*x975
      x977 = x755*x779
      x978 = rho_s*x977 + x289*x976 + x972
      x979 = rho_s*x937 + x230*(rho_s*x958 + x277*x973) + x237*x887 + x2
     &78*(rho_s* x941 + x284*x887) + x285*(rho_s*x974 + x295*x978 + x296
     &*x887) + x661*x940 + x697*x956
      x980 = x204*x979
      x981 = x887 - x980
      x982 = x200*x981
      x983 = 2.78d0*x935 + 2.78d0*x982
      x984 = x792*x983
      x985 = x323*x984
      x986 = x798*x800
      x987 = x891*x986
      x988 = gamma_cc*x987
      x989 = x899 + x980
      x990 = x302*x897 + x44*x989
      x991 = x303*x990
      x992 = x810*x991
      x993 = x812*x983
      x994 = x988 - x992 + x993
      x995 = x319*x816
      x996 = x994*x995
      x997 = gamma_cc*x996
      x998 = x586*x824
      x999 = x908*x998
      x1000 = x826*x829
      x1001 = x1000*x991
      x1002 = 5.56d0*x935 + 5.56d0*x982
      x1003 = x1001 - x1002*x832 - x988 + x992 - x993 - x999
      x1004 = x1003*x836
      x1005 = gamma_cc*x1004
      x1006 = x1005 + x933 + x985 + x997
      x1007 = x1006*x841 + x325*x916 + x981
      x1008 = x1007 + x931
      x1009 = x48*x50
      x1010 = x54*x61
      x1011 = x1009*x1010
      x1012 = 0.296858062715705d0*gamma_cc
      x1013 = x1011*x1012
      x1014 = gamma_cc*x49
      x1015 = x1014*x50
      x1016 = -0.593716125431411d0*x1015 - x58
      x1017 = gamma_cc*x596
      x1018 = x1016*x1017
      x1019 = x1013 + x1018 + x64
      x1020 = x1019*x601
      x1021 = -x1020
      x1022 = x134*x142
      x1023 = x140*x149
      x1024 = x1022*x1023
      x1025 = 0.111863791688847d0*x1024
      x1026 = 0.167230224308322d0*x144
      x1027 = x1026*x130
      x1028 = x1027*x150
      x1029 = 0.25d0*gamma_cs
      x1030 = 0.125d0*gamma_cc + 0.125d0*gamma_ss
      x1031 = -x1029 + x1030
      x1032 = x133*x141
      x1033 = -x1026*x143 - x1031*x1032
      x1034 = x1033*x436
      x1035 = x1034*x397
      x1036 = x1025 + x1028 + x1035
      x1037 = x129*x448
      x1038 = x1036*x1037
      x1039 = x112*x91
      x1040 = x1039*x109
      x1041 = 0.111863791688847d0*x1040
      x1042 = x1041*x88
      x1043 = 0.167230224308322d0*x100
      x1044 = x1043*x114
      x1045 = x1029 + x1030
      x1046 = x104*x110
      x1047 = -x1043*x96 - x1045*x1046
      x1048 = x103*x362
      x1049 = x1047*x1048
      x1050 = x1049*x357
      x1051 = x1042 + x1044 + x1050
      x1052 = x379*x94
      x1053 = x1051*x1052
      x1054 = x1038 + x1053
      x1055 = x1054*x179
      x1056 = -x1055*x183
      x1057 = x180*x611
      x1058 = 0.00081628665578855d0*x477
      x1059 = rho_c**(-7.66666666666667d0)
      x1060 = x166*x169
      x1061 = x1060*x182
      x1062 = x1059*x1061
      x1063 = x1062*x159
      x1064 = x184*x611
      x1065 = rho_c**(-6.66666666666667d0)
      x1066 = x1061*x1065
      x1067 = gamma_cc*x1064 - x1058*x1066
      x1068 = x1020*x185 + x1067*x71
      x1069 = -gamma_cc*x1057 + x1056 + x1058*x1063 + x1068
      x1070 = 0.005859375d0*x191
      x1071 = x193*x51
      x1072 = -x1055*x188
      x1073 = x180*x622
      x1074 = 0.00228560263620794d0*x477
      x1075 = x184*x622
      x1076 = gamma_cc*x1075 - x1066*x1074
      x1077 = x1020*x621 + x1076*x71
      x1078 = -gamma_cc*x1073 + x1063*x1074 + x1072 + x1077
      x1079 = x1070*x1071 + x1078*x626
      x1080 = x156*(x1021 + x1069*x196 + x1079*x187)
      x1081 = gamma_cc*x1009
      x1082 = x301*x309
      x1083 = x1082*x318
      x1084 = 0.00922952226005862d0*x1083
      x1085 = x1081*x1084
      x1086 = 0.000573904446924163d0*x306
      x1087 = x1081*x1086
      x1088 = -x1087 - x314
      x1089 = gamma_cc*x836
      x1090 = x1088*x1089
      x1091 = x1085 + x1090 + x322
      x1092 = x1091*x841
      x1093 = x1080 + x1092
      x1094 = x182*x465
      x1095 = x1094*x169
      x1096 = x1095*x159
      x1097 = x42*x66
      x1098 = 0.0870539344350337d0*x1097 - 0.174107868870067d0*x16 + 0.0
     &553050111671826d0*x31 + 2.8d0*x39
      x1099 = x169*x50
      x1100 = x1099*x182
      x1101 = x1058*x1100
      x1102 = -x1098*x1101
      x1103 = x1072 + x1074*x1096 + x1102
      x1104 = -x1101*x71
      x1105 = x1056 + x1058*x1096 + x1104
      x1106 = x156*(x1103*x928 + x1105*x196)
      x1107 = x1065*x159
      x1108 = 0.00457120527241588d0*x477
      x1109 = x169*x182
      x1110 = rho_s*x1109
      x1111 = x1108*x1110
      x1112 = 0.0016325733115771d0*x477
      x1113 = x1112*x169
      x1114 = x1094*x1113
      x1115 = rho_s*x1114
      x1116 = x1098*x1115
      x1117 = 0.223727583377694d0*x1024
      x1118 = 0.334460448616645d0*x144
      x1119 = x1118*x151
      x1120 = x1118*x143
      x1121 = x84 - x85 - x86
      x1122 = -x1032*x1121 + x1120
      x1123 = x1122*x437
      x1124 = x1123*x152
      x1125 = -x1117 - x1119 + x1124
      x1126 = 0.223727583377694d0*x1040
      x1127 = x1126*x88
      x1128 = 0.334460448616645d0*x100
      x1129 = x1128*x114
      x1130 = -x1128*x96
      x1131 = -x1046*x88 + x1130
      x1132 = x1048*x1131
      x1133 = x1132*x95
      x1134 = x1133*x115
      x1135 = x1127 + x1129 + x1134
      x1136 = x1037*x1125 + x1052*x1135
      x1137 = x1136*x179
      x1138 = -x1107*x1111 + x1116 - x1137*x188
      x1139 = rho_s*x1113*x182
      x1140 = x1115*x71
      x1141 = -x1107*x1139 - x1137*x183 + x1140
      x1142 = x156*(x1138*x928 + x1141*x196)
      x1143 = x1071*x160
      x1144 = x1143*x482
      x1145 = 0.03125d0*x1071
      x1146 = x184*x71
      x1147 = -x1145*x1146
      x1148 = 0.03125d0*x1144 + x1147
      x1149 = x1145*x184
      x1150 = -x1098*x1149
      x1151 = 0.0875d0*x1144 + x1150
      x1152 = tau_c**(-4)
      x1153 = x1152*x192
      x1154 = -x1070*x1153 + x1151*x626
      x1155 = x156*(x1148*x196 + x1154*x187)
      x1156 = x182*x456
      x1157 = x1156*x482
      x1158 = x611*x612
      x1159 = x453*x479
      x1160 = x74**(-2)
      x1161 = x1160*x76
      x1162 = x1161*x82
      x1163 = x336*(-130.15964257673d0*x1160*x77 - 172.632095841468d0*x1
     &160*x78 - 57.042190297332d0*x1160*x79 - 156.836877485407d0*x1161)
      x1164 = 1034.52086599683d0*x333**2/(x80**4*x81**2)
      x1165 = x333*x334
      x1166 = x1165*x335
      x1167 = x1166*x330
      x1168 = 32.1639684429149d0*x1165*(35.5123722979901d0*x329 + 37.442
     &0464663507d0* x330 + 36.6336685193113d0*x331 + 10.8942558607307d0*
     &x332)/x80**3
      x1169 = x116**(-2)
      x1170 = 0.0155453454348274d0*x1169
      x1171 = x1170*x377
      x1172 = x107*x88
      x1173 = 1.00338134584993d0*x100
      x1174 = 6.60198463131432d0*x342
      x1175 = x1174*x88
      x1176 = 65.1589765731441d0*x88*x99**(-2.33333333333333d0)
      x1177 = 1.33784179446658d0*x100
      x1178 = x1177*x88
      x1179 = x1178*x338
      x1180 = x1175*x95
      x1181 = x1178*x95
      x1182 = x103*x357/x111**3
      x1183 = x1182*x373
      x1184 = x1172*x91
      x1185 = x1173*x1184
      x1186 = x1174*x347
      x1187 = x1176*x96
      x1188 = x1179*x356
      x1189 = x1180*x356
      x1190 = exp(3.999999999999988d0*x83)
      x1191 = x354**2
      x1192 = x1191*x370
      x1193 = x1181*x1190*x1192
      x1194 = x104*x357*x89
      x1195 = x1191*x1194
      x1196 = 10.3996132712785d0*x1162 + x1163*x353 - x1164*x353 - 79.85
     &42817192673d0* x1167 + x1168*x353
      x1197 = x1194*x1196
      x1198 = x113*x115
      x1199 = x94**(-6)
      x1200 = 17.6648219323966d0*x366
      x1201 = x108*x363
      x1202 = x105*x371
      x1203 = 0.894910333510778d0*x369*x89
      x1204 = x1048*x357
      x1205 = 16.0819842214575d0*x1166*x328 + x337*(-0.00460691662213167
     &d0*x76 - 0.0155453454348274d0) + x380 + x94*(-0.16166558069064d0*x
     &1162 + x1163*x328 - x1164*x328 + 1.24136239377604d0*x1167 + x1168*
     &x328 + x1171*(x341 + x345 - x361 - x376) + x379*(-x112*x1179*x359 
     &+ x114 *x1176 + x1172*x1173*x339 + x1175*x340 - x1175*x360 - x1179
     &*x375 - x1180*x375 + x1181*x359*x374 + x1183*(x1177*x347 + x1181*x
     &356 + 6.60198463131432d0*x350 + 1.78982066702156d0*x365 + 8.832410
     &96619832d0*x368 + 1.78982066702156d0*x369*x371) + x1198*( x1185 + 
     &x1186 + x1187 + x1188 + x1189 + x1193 - x1195 - x1197) + x1204*(-x
     &106*x1200*x363 - x107*x1200*x1202 - x1185 - x1186 - x1187 - x1188 
     &- x1189 - 2.68473100053233d0*x1190*x1191*x369/x90** 4 + x1192*x120
     &3 - x1193 + x1195 + x1196*x1203*x370 + x1197 - 2.23727583377694d0*
     &x1199*x364 - 3.57964133404311d0*x1201*x1202 - 108.965502680276d0*x
     &367*x99**(-2.66666666666667d0))))
      x1206 = x402**2
      x1207 = exp(3.999999999999988d0*x128)
      x1208 = x1206*x1207
      x1209 = 1.33784179446658d0*x144
      x1210 = x1209*x134
      x1211 = x1210*x130
      x1212 = x1211*x427
      x1213 = x1208*x1212
      x1214 = x1210*x406
      x1215 = x1214*x404
      x1216 = 6.60198463131432d0*x410
      x1217 = x1216*x134
      x1218 = x130*x404
      x1219 = x1217*x1218
      x1220 = x1206*x847
      x1221 = 1034.52086599683d0/(x125**4*x126**2)
      x1222 = x1221*x391**2
      x1223 = 35.5123722979901d0*x383
      x1224 = 37.4420464663507d0*x385
      x1225 = 36.6336685193113d0*x387
      x1226 = 10.8942558607307d0*x389
      x1227 = 32.1639684429149d0*x393/x125**3
      x1228 = x1227*(x1223 + x1224 + x1225 + x1226)
      x1229 = x391*x394
      x1230 = x1229*x385
      x1231 = x119**(-2)
      x1232 = x121*x1231
      x1233 = x1232*x127
      x1234 = 10.3996132712785d0*x1233
      x1235 = 130.15964257673d0*x122*x1231
      x1236 = 156.836877485407d0*x1232
      x1237 = 172.632095841468d0*x123*x1231
      x1238 = 57.042190297332d0*x1231*x124
      x1239 = -x1235 - x1236 - x1237 - x1238
      x1240 = x1234 + x1239*x845
      x1241 = -x1222*x400 + x1228*x401 - 79.8542817192673d0*x1230 + x124
     &0
      x1242 = x1241*x847
      x1243 = x1022*x137
      x1244 = 1.00338134584993d0*x144
      x1245 = x1243*x1244
      x1246 = 65.1589765731441d0*x134*x138**(-2.33333333333333d0)
      x1247 = x1246*x143
      x1248 = x1216*x408
      x1249 = x1245 + x1247 + x1248
      x1250 = x149*x406
      x1251 = x1250*x415
      x1252 = x415*x416
      x1253 = x139*x419
      x1254 = x135*x429
      x1255 = x1253*x1254
      x1256 = 17.6648219323966d0*x423
      x1257 = x1256*x137
      x1258 = 2.68473100053233d0*x430/x132**4
      x1259 = x129**(-6)
      x1260 = 2.23727583377694d0*x1259*x420
      x1261 = 108.965502680276d0*x138**(-2.66666666666667d0)*x424
      x1262 = x1256*x136*x419
      x1263 = -x1245 - x1247 - x1248
      x1264 = -x1260 - x1261 - x1262 + x1263
      x1265 = x434*x435
      x1266 = x434*x436
      x1267 = x434*x437
      x1268 = x1209*x408
      x1269 = 6.60198463131432d0*x412
      x1270 = 1.78982066702156d0*x421
      x1271 = 8.83241096619832d0*x425
      x1272 = 1.78982066702156d0*x430
      x1273 = x147*x397/x148**3
      x1274 = x1273*x434
      x1275 = x1244*x134*x137*x150
      x1276 = x1246*x151
      x1277 = x1217*x440
      x1278 = x1275 + x1276 + x1277
      x1279 = x153**(-2)
      x1280 = 0.0155453454348274d0*x1279
      x1281 = x1280*x446
      x1282 = x391*x392
      x1283 = 0.16166558069064d0*x1233
      x1284 = x1239*x396 - x1283
      x1285 = 16.0819842214575d0*x392
      x1286 = x1229*x1285
      x1287 = 0.00460691662213167d0*x121
      x1288 = x395*(-x1287 - 0.0155453454348274d0)
      x1289 = x1286 + x1288*x391 + x129*(-x1222*x392 + x1228*x1282 + 1.2
     &4136239377604d0*x1230 + x1281*(-x418 - x439 + x854) + x1284 + x448
     &*(-x1210*x1251 + x1211*x1265*x415 - x1214*x1266 - x1217*x1252 - x1
     &217*x1267 + x1274*(x1211*x404 + x1268 + x1269 + x1270 + x1271 + x1
     &272*x429) + x1278 + x417*(x1213 + x1215 + x1219 - x1220 - x1242 + 
     &x1249) + x438*(x1206*x851 - x1208*x1258 - x1213 - x1215 - x1219 + 
     &x1220 + x1241*x851 + x1242 - x1254*x1257 - 3.57964133404311d0*x125
     &5 + x1264))) + x450
      x1290 = x179*(x1205 + x1289)
      x1291 = x161*x22
      x1292 = x166*x21
      x1293 = x176*(10.0d0*x1291 + 5.22d0*x1292 + 94.92d0*x162/rho_c**8)
     &
      x1294 = x460*x476
      x1295 = x1294*x475
      x1296 = 0.065302932463084d0*x50
      x1297 = 0.065302932463084d0*x474
      x1298 = x168/x175**6
      x1299 = x1298*x475
      x1300 = x1299*(-x1296*x464 - x1297*x472 + 0.304747018161059d0*x466
     &)
      x1301 = 0.104484691940934d0*x50
      x1302 = x1301*x169
      x1303 = 0.487595229057694d0*x465
      x1304 = rho_s*x160
      x1305 = 2.66666666666667d0*x1304
      x1306 = x2**(-3.33333333333333d0)
      x1307 = 3.11111111111111d0*x1292
      x1308 = x4**(-3.33333333333333d0)
      x1309 = x477*(-gamma_ss*x1302 - 1.38151981566347d0*x1065*x174 - x1
     &301*x463*x472 + x1303*x173*x472 + x1303*x464 - x872*(-x1305*x467 +
     & x1305*x470 + x1306*x1307 + x1307*x1308))
      x1310 = x1293 + 2.0d0*x1295 + x1300 + x1309
      x1311 = x1310*x480
      x1312 = 0.09375d0*x483
      x1313 = x178*x453
      x1314 = x159*x479
      x1315 = gamma_cc*rho_c**(-4.33333333333333d0)
      x1316 = 4.23769752780008d0*x1315
      x1317 = 2.54261851668005d0*gamma_cc
      x1318 = x1317*x541
      x1319 = x553*x61
      x1320 = x48*x551
      x1321 = x1320*(2.0d0*x547 - 2.0d0*x550)
      x1322 = gamma_cc*x1321
      x1323 = x2**(-1.33333333333333d0)
      x1324 = 0.444444444444444d0*x1292
      x1325 = x4**(-1.33333333333333d0)
      x1326 = -x1305*x545 + x1305*x548 + x1323*x1324 + x1324*x1325
      x1327 = gamma_cc*x552
      x1328 = x1326*x1327
      x1329 = 1.08969365000574d0*x56
      x1330 = gamma_cc*x1329
      x1331 = x1330*x579
      x1332 = x59*x595
      x1333 = x593*x594
      x1334 = x1330*x553
      x1335 = x551*x560
      x1336 = x585*x587
      x1337 = 0.593716125431411d0*x52
      x1338 = x576*x7
      x1339 = 1.18743225086282d0*x199
      x1340 = gamma_cc*x63
      x1341 = x1340*x7/x60**3
      x1342 = x1341*x593
      x1343 = x1316*x55
      x1344 = x1318*x1335
      x1345 = gamma_cc*x57
      x1346 = x1345*x54
      x1347 = x1321*x1346
      x1348 = gamma_cc*x561
      x1349 = x1326*x1348
      x1350 = x1318*x1338
      x1351 = x574**2
      x1352 = 64.3279368858296d0*x43
      x1353 = exp(x1352*x40)
      x1354 = x1351*x1353
      x1355 = x1330*x7
      x1356 = x1355*x589
      x1357 = x1354*x1356
      x1358 = x1351*x45
      x1359 = gamma_cc*x902
      x1360 = x1358*x1359
      x1361 = x1334*x576
      x1362 = 32.1639684429148d0/x41**5
      x1363 = x1362*x565
      x1364 = x1363*x551
      x1365 = x160*x867
      x1366 = x1304*x548
      x1367 = 0.333333333333333d0*x1292
      x1368 = x1323*x1367
      x1369 = x1325*x1367
      x1370 = x567*(-x1365*x545 + 2.0d0*x1366 + x1368 + x1369)
      x1371 = 64.3279368858296d0*x566
      x1372 = x1371*x565
      x1373 = 5.1297628024841d0*x1304
      x1374 = x2**(-0.666666666666667d0)
      x1375 = 0.854960467080683d0*x1292
      x1376 = x4**(-0.666666666666667d0)
      x1377 = x1373*x485 - x1373*x488 + x1374*x1375 + x1375*x1376
      x1378 = x1377*x29
      x1379 = 20.0d0*x1291
      x1380 = x1379*x20
      x1381 = 8.0d0*x457
      x1382 = x1381*x490
      x1383 = 0.000908766572494496d0*x492
      x1384 = x1383*x518
      x1385 = 29.6085746432167d0*x528
      x1386 = x1385*x571
      x1387 = x165*x9
      x1388 = x165*x8
      x1389 = x10*x165
      x1390 = x529*(-1.5861639701226d0*x1387 - 0.998929717145607d0*x1388
     & - 0.322572935345748d0*x1389 - 0.212390288347482d0*x737)
      x1391 = x24*x490
      x1392 = 59.2171492864334d0*x531
      x1393 = x1391*x1392
      x1394 = 876.667692402934d0*x530**2/(x26**4*x27**2)
      x1395 = x530*(x0*x25 + 1.49839457571841d0*x492 + 2.71913823449588d
     &0*x495 + 0.254868346016978d0*x497)/x26**3
      x1396 = x20*x531
      x1397 = x1396*x569
      x1398 = 236.868597145733d0*x457
      x1399 = x510*(-1.16327369584438d0*x1387 - 0.989141964955861d0*x138
     &8 - 0.600321476249523d0*x1389 - 0.210778258416395d0*x737)
      x1400 = 258.630216499205d0*x507**2/(x13**4*x14**2)
      x1401 = x507*x508
      x1402 = 16.0819842214574d0*x1401*(x0*x11 + 1.48371294743379d0*x492
     & + 1.99418347859036d0*x495 + 0.252933910099675d0*x497)/x13**3
      x1403 = 0.00366369660967365d0*x1388*x15
      x1404 = x1401*x509
      x1405 = 0.0883792666034679d0*x1404*x492
      x1406 = x1403 - x1405
      x1407 = x1399*x506 - x1400*x506 + x1402*x506 + x1406
      x1408 = x1377*x38
      x1409 = x493*x501
      x1410 = x514*(-0.00176138600691563d0*x1388*x36 + x1407 + x1409*x49
     &4*( -2.16229511226841d0*x1387 - 1.70877610553209d0*x1388 - 1.23355
     &033167571d0*x1389 - 0.267319032365355d0*x737) + x1409*x498 *(x0*x3
     &3 + 2.56316415829814d0*x492 + 3.70679162103156d0*x495 + 0.32078283
     &8838426d0*x497)/x34**3 + 0.0849797459133392d0*x492*x499 - 1034.520
     &86599683d0*x498**2*x500/(x34**4*x35**2))
      x1411 = x1379*x520
      x1412 = x23*x513
      x1413 = 2.0d0*x1412*x490
      x1414 = x37*x490
      x1415 = x1381*x1414
      x1416 = x20*x513
      x1417 = x1381*x1416
      x1418 = x20*x28*x8
      x1419 = 0.00363506628997799d0*x1291*x1418
      x1420 = 0.000605844381662998d0*x1388*x30
      x1421 = 0.0269072828949635d0*x492*x531*x533
      x1422 = -x1408 - x1410 - x1411 - x1413 + x1415 + x1417 + x1419 - x
     &1420 + x1421
      x1423 = x1378*x569 - x1380*x570 + x1382*x570 + x1384 + x1386*x1390
     & + x1386*x1395 + x1393*x569 - x1394*x571 + x1397*x1398 + x1407 + x
     &1422
      x1424 = x45*(x1364*x40 + x1370*x40 + x1372*x573 + x1423*x44)
      x1425 = x1359*x1424
      x1426 = gamma_cc*x580
      x1427 = 2.77067525201325d0*x465
      x1428 = x1427*x51
      x1429 = x6**(-6)
      x1430 = x1429*x47
      x1431 = x1430*x53
      x1432 = x585*(3.33333333333333d0*x547 - 3.33333333333333d0*x550)
      x1433 = x1304*x545
      x1434 = 0.888888888888889d0*x1292
      x1435 = x1323*x1434 + x1325*x1434 + 5.33333333333333d0*x1366 - 5.3
     &3333333333333d0*x1433
      x1436 = x48*x819
      x1437 = 1.78114837629423d0*x199/x46**4
      x1438 = 1.18743225086282d0*x52
      x1439 = x601*(x1017*(-7.85024654737088d0*x1065*x583 + x1336*x1428 
     &- x1343 + x1344 - x1347 - x1349 - x1350 - x1354*x1437 - x1357 + x1
     &358*x909 + x1360 + x1361 + x1424*x909 + x1425 - x1431*x1432 - x143
     &5*x588 - 5.5413505040265d0*x1436*x590 + x1438*x590*x823) + x1316*x
     &540 - x1318*x1319*x59 - x1318*x1332*x593 - x1318*x579*x62 + x1319*
     &x1331 + x1322*x554 + x1328*x554 + x1331*x593*x595 + x1333*x1334*x5
     &9 + x1342*(x1317*x557 - x1330*x1335 + x1330*x1338 - x1336*x1337 + 
     &x1339*x590 + 2.77067525201325d0*x584) + x1426*(x1343 - x1344 + x13
     &47 + x1349 + x1350 + x1357 - x1360 - x1361 - x1425))
      x1440 = x605*(-x1365*x548 - x1368 - x1369 + 2.0d0*x1433)
      x1441 = x1440*x66
      x1442 = x67/x65**2
      x1443 = x1442*x599
      x1444 = x1443*(x544 - x556 - x582 - x598)
      x1445 = 0.0310906908696549d0*x41
      x1446 = x1445*x603
      x1447 = x1446*(-0.666666666666667d0*x547 + 0.666666666666667d0*x55
     &0)
      x1448 = x1447*x66
      x1449 = x603*x604
      x1450 = 0.0621813817393098d0*x1449
      x1451 = x599*x600
      x1452 = x1450*x1451
      x1453 = x1380*x525
      x1454 = x1382*x525
      x1455 = x1399*x523
      x1456 = x1400*x523
      x1457 = x1396*x517
      x1458 = x1398*x1457
      x1459 = x1385*x534
      x1460 = x1390*x1459
      x1461 = x1394*x534
      x1462 = x1402*x523
      x1463 = x1395*x1459
      x1464 = x1378*x517 - x1384 + x1393*x517 - x1403 + x1405 + x1408 + 
     &x1410 + x1411 + x1413 - x1415 - x1417 - x1419 + x1420 - x1421 - x1
     &453 + x1454 + x1455 - x1456 + x1458 + x1460 - x1461 + x1462 + x146
     &3
      x1465 = x1439 + x1441 + x1444 + x1448 + x1452 + x1464
      x1466 = 0.0625d0*x160
      x1467 = x182*x461
      x1468 = x1293*x165
      x1469 = x165*x51
      x1470 = x1295*x1469
      x1471 = x182*x478
      x1472 = x165*x183
      x1473 = x1309*x165
      x1474 = 2.0d0*x609
      x1475 = x1465*x185 + x1474*x615 + x71*(x1300*x1472 + x1312*x177 - 
     &x1466*x1467 - x1466*x1471 + x1468*x183 + x1470*x611 + x1473*x183)
      x1476 = 0.0234375d0*x194*x456
      x1477 = 0.01171875d0*x618
      x1478 = 0.2625d0*x483
      x1479 = x612*x622
      x1480 = 0.175d0*x160
      x1481 = x165*x188
      x1482 = x1465*x621 + x1474*x623 + x71*(x1300*x1481 - x1467*x1480 +
     & x1468*x188 + x1470*x622 - x1471*x1480 + x1473*x188 + x1478*x177)
      x1483 = x1423 - x1439 - x1441 - x1444 - x1448 - x1452
      x1484 = 0.0244087849784238d0*rho_c**(-2.66666666666667d0)*x201
      x1485 = 0.312466222036391d0*x630
      x1486 = 0.208310814690927d0*mu*rho_c**(-2.33333333333333d0)
      x1487 = x553*x636
      x1488 = gamma_cc*x541
      x1489 = 2.54261851668005d0*x1488
      x1490 = x1489*x316
      x1491 = x790**2
      x1492 = x299**(-2)
      x1493 = x1492*x641
      x1494 = x1493*x323
      x1495 = 1.08969365000574d0*x1487
      x1496 = gamma_cc*x320
      x1497 = x539*x791
      x1498 = x323*x790
      x1499 = x1492*x70
      x1500 = x1499*(x640 + x803)
      x1501 = x807*x814
      x1502 = 2.78d0*x640
      x1503 = 0.00549554491451047d0*x502
      x1504 = x20*x525
      x1505 = x1381*x1504
      x1506 = x1383*x30
      x1507 = 2.0d0*x519
      x1508 = x1404*x523
      x1509 = 32.1639684429148d0*x1508
      x1510 = x1381*x520
      x1511 = 2.0d0*x491
      x1512 = x1392*x534
      x1513 = 2.0d0*x515
      x1514 = x299/x70**3
      x1515 = x1514*(-x1503 - x1505 + x1506 - x1507 - x1509 + x1510 - x1
     &511 - x1512 - x1513)
      x1516 = x23 - 1.0d0
      x1517 = x1516*x525
      x1518 = x1516*x490
      x1519 = x517*x532
      x1520 = x643*(x1377*x1517 - x1391*x1519 + x1406 + x1422 + x1453 - 
     &x1454 - x1455 + x1456 - x1458 - x1460 + x1461 - x1462 - x1463 + x1
     &518*x1519 - 0.000454383286247248d0*x1518*x28*x492 + x518*x536)
      x1521 = x642*x788
      x1522 = 5.56d0*x640
      x1523 = x1292*x211
      x1524 = x1523*x205
      x1525 = x165*x213
      x1526 = x165*x210
      x1527 = x207*(0.00571025734089666d0*x1388 + 0.000856522469028669d0
     &*x1526 + 0.0350330449527535d0*x737 - 0.00609784657412857d0*x957)
      x1528 = x205*x652
      x1529 = x218*x235
      x1530 = x1529*x660
      x1531 = x1365*x657
      x1532 = x1365*x654
      x1533 = x2**(-1.33333333333333d0)
      x1534 = 0.333333333333333d0*x1292
      x1535 = x1533*x1534
      x1536 = x4**(-1.33333333333333d0)
      x1537 = x1534*x1536
      x1538 = x220*x228
      x1539 = x1538*x694
      x1540 = x221*x695
      x1541 = x1539*x1540
      x1542 = x694*x696
      x1543 = x0*x539
      x1544 = x219/x234**2
      x1545 = x1539*x1544
      x1546 = x687*x692
      x1547 = x680*x681
      x1548 = 9.20897059521559d0*x9
      x1549 = x234/x228**3
      x1550 = x1549*x691
      x1551 = 0.89531658564596d0*x1387*x225
      x1552 = x1388*x675
      x1553 = 1.53482843253593d0*x1547*x495
      x1554 = x231*x666
      x1555 = x201*x492
      x1556 = x1554*x1555
      x1557 = mu*x226
      x1558 = x1557*x680
      x1559 = x1554*x1558
      x1560 = x223*x666
      x1561 = x1560*x672
      x1562 = 0.666666666666667d0*x1304
      x1563 = 0.111111111111111d0*x1292
      x1564 = x682*(x1533*x1563 + x1536*x1563 + x1562*x654 - x1562*x657)
     &
      x1565 = 1.33333333333333d0*x1304
      x1566 = 0.222222222222222d0*x1292
      x1567 = x667*(x1533*x1566 + x1536*x1566 + x1565*x654 - x1565*x657)
     &
      x1568 = x673/x218**5
      x1569 = x1568*x671
      x1570 = x1292*x280
      x1571 = x279*(-0.207448236981382d0*x1388*x700 + x207*(0.0028551286
     &7044833d0* x1388 + 0.000428261234514334d0*x1526 - x165*x209 + 0.01
     &75165224763768d0*x737) + 0.311172355472072d0*x492*x701 + 0.0242070
     &587025094d0*x700*x737)
      x1572 = 0.0891107483661069d0*x957
      x1573 = x0*x208
      x1574 = x242*x710
      x1575 = 75.398223686155d0*x1574
      x1576 = x165*x244
      x1577 = 0.00786976946662483d0*x165
      x1578 = x1577*x254
      x1579 = 0.0157395389332497d0*x1304
      x1580 = x252*x254
      x1581 = 0.00786976946662483d0*x1292
      x1582 = x2**(-2)
      x1583 = x1582*x254
      x1584 = x260*x722
      x1585 = x1577*x267
      x1586 = x263*x267
      x1587 = x4**(-2)
      x1588 = x1587*x267
      x1589 = 7.65662400050185d0*x755
      x1590 = x1589*x262
      x1591 = x1589*x730
      x1592 = 45.9397440030111d0*x727
      x1593 = x259*x721
      x1594 = x265*x725
      x1595 = x1594*x273
      x1596 = x1292*x289
      x1597 = 0.150031049774234d0*x165
      x1598 = x1597*x254
      x1599 = 0.0271505358982603d0*x737
      x1600 = x1599*x255
      x1601 = 0.300062099548467d0*x1304
      x1602 = x12*x255
      x1603 = x1304*x252
      x1604 = 0.150031049774234d0*x1292
      x1605 = 0.0271505358982603d0*x1292
      x1606 = x1582*x1602
      x1607 = x1597*x267
      x1608 = x1599*x264
      x1609 = x12*x264
      x1610 = x1304*x263
      x1611 = x1587*x1609
      x1612 = 68.9096160045166d0*x247
      x1613 = x721*x742
      x1614 = x1594*x750
      x1615 = x247*x759
      x1616 = 22.9698720015055d0*x1304
      x1617 = x250*x258
      x1618 = x1617*x721
      x1619 = 68.9096160045166d0*x756
      x1620 = x247*x754
      x1621 = x741*x743
      x1622 = x749*x751
      x1623 = 0.225046574661351d0*x0
      x1624 = x1623*x254
      x1625 = 0.0325806430779124d0*x497
      x1626 = x1625*x255
      x1627 = 0.225046574661351d0*x165
      x1628 = 0.0325806430779124d0*x737
      x1629 = x261/x257**3
      x1630 = x1629*x248*x741
      x1631 = x1623*x267
      x1632 = x1625*x264
      x1633 = x268/x271**3
      x1634 = x1633*x248*x749
      x1635 = x269*x272
      x1636 = x1594*x1635
      x1637 = 45.9397440030111d0*x755
      x1638 = x1637*x719
      x1639 = 7.65662400050185d0*x1292
      x1640 = x247*x262
      x1641 = x1582*x1640
      x1642 = x247*x730
      x1643 = x1587*x1642
      x1644 = x250*x743
      x1645 = x1644*x740
      x1646 = x1637*x723
      x1647 = 45.9397440030111d0*x1292
      x1648 = x269*x751
      x1649 = x1648*x748
      x1650 = -x1584*(-x1578 + x1579*x1580 - x1581*x1583) + x1590 + x159
     &1 - x1592* x1593 - x1592*x1595 - x1592*x1621 - x1592*x1622 + x1593
     &*x1638 - x1595*x1646 - x1596*x753 - x1596*x758 - 15.3132480010037d
     &0*x1603* x1640 + 15.3132480010037d0*x1610*x1642 - x1612*x1613*x741
     & - x1612 *x1614*x749 + x1615*x1616 - x1615*x1647*x263 - x1616*x162
     &0 + x1618*x1619 - x1619*x1636 + x1619*x1645 - x1619*x1649 - x1620*
     & x1647*x252 + x1621*x1638 - x1622*x1646 - x1630*(x1624 + x1626 - x
     &1627*x720 - x1628*x739) - x1634*(x1627*x724 + x1628*x747 + x1631 +
     & x1632) + x1639*x1641 + x1639*x1643 - x726*(-x1579*x1586 - x1581 *
     &x1588 - x1585) - x966*(x1580*x1601 - x1583*x1604 - x1598 - x1600 +
     & 0.0543010717965206d0*x1602*x1603 - x1605*x1606) - x968*(-x1586* x
     &1601 - x1588*x1604 - x1605*x1611 - x1607 - x1608 - 0.0543010717965
     &206d0*x1609*x1610)
      x1651 = 2.66666666666667d0*x1304
      x1652 = 2.22222222222222d0*x1292
      x1653 = x1292*x291
      x1654 = 5.29055759308205d0*x777
      x1655 = x204*(6.39350882283002d-5*x0**2.33333333333333d0*x1525 + x
     &1464*x237 + 0.0334548885448792d0*x1524*x214 - 0.00164372154798406d
     &0*x1524* x648 - 0.00191767513931474d0*x1525*x648 + x1527*x206*x653
     & - 0.00880391803812612d0*x1528*x647 + 0.000821860773992031d0*x1528
     &* x649 + x1530*(-x664 + x665) + 2.0d0*x1541*x660 + x1542*(x686 - x
     &688 + x689) - 0.0168588246652672d0*x1543*x236 + x1545*(-x669 - x67
     &4 + x677 + x679 + x690) + 0.0309078452196566d0*x165*x662 - 0.01056
     &47016457513d0*x205*x645*x651 + 0.0117385573841682d0*x213* x936 + x
     &230*(0.445553741830535d0*x1523*x208 - 0.0277288014669606d0*x1523*x
     &210 - 0.0231073345558005d0*x1526*x212 + x1527*x708 + x1572*x212 - 
     &0.0891107483661069d0*x1573*x652 + 0.00107855455836145d0*x212*x936 
     &+ 0.000586927869208408d0*x276* x936 + x277*(-0.154912426780983d0*x
     &0**(-0.333333333333333d0)* x1576*x715 + x1292*x712 - 9.66648718250
     &747d0*x1292*x716 - 150.79644737231d0*x1574*x459 - x1575*x713 - 1.6
     &1108119708458d0* x1576*x716 + x1650 - 4.83324359125373d0*x241*x710
     &*x713*x714 - x711*(-1.36477107997868d0*x1388 + 0.427594146176807d0
     &*x737)) - 0.000440195901906306d0*x647*x764 + 0.0138644007334803d0*
     &x650*x652 - 0.178221496732214d0*x651*x706) + x278*(x1464*x284 - 0.
     &213300270968718d0*x1543*x283 + 0.0628942271634173d0*x1570*x214 - x
     &1571*x282 + 0.319950406453076d0*x165*x704 - 0.0073560499606336d0*x
     &281*x936 - 0.0198613348937107d0*x645*x702 + 0.0055170374704752d0*x
     &647*x703) + x285*(x1464*x296 + 1.04949664671409d0*x1526*x70 + 0.83
     &7628205355044d0*x1570*x208 - x1571*x286 + 0.055841880357003d0*x157
     &3*x703 - 0.055841880357003d0 *x281*x957 + 0.0009195062450792d0*x29
     &4*x936 + x295*(-0.299209d0* x165*x292 + x1650 + 77.9555417944151d0
     &*x1653*x247 - 10.5811151861641d0*x1653*x777 - x1654*x713*x775 - 31
     &.182216717766d0*x244*x727*x775 - 3.82831200025092d0*x288*x755 + x2
     &89*(x1651*x768 - x1651*x771 + x1652*x216 + x1652*x217) - x290* x77
     &6*(-0.106975995763985d0*x1388 + 0.289053642815522d0*x737) + 5.1970
     &3611962767d0*x292*x755 + x728*x773 - 93.5466501532981d0* x775*x780
     &) - 0.899568554326366d0*x539*x650 - 0.0006896296838094d0 *x647*x78
     &1 - 0.335051282142018d0*x702*x706) + x661*(x1531 - x1532 - x1535 -
     & x1537) + x697*(2.0d0*x1546*x691 + x1550*(-x1547*x1548 - 4.2785697
     &4274298d0*x668 + 1.42618991424766d0*x676 + 1.53482843253593d0*x684
     &) + x229*(1.43555326817173d0*x1389*x232 + x1551 + 2.05391854581537
     &d0*x1552 - x1553 - 3.08087781872305d0* x1556 + x1559 + 4.621316728
     &08458d0*x1561 + x1564 + 4.62131672808458d0*x1567 + x1569*(-1.33333
     &333333333d0*x656 + 1.33333333333333d0*x659) - 1.91407102422897d0*x
     &230*x496*x672 + x945*(-x1531 + x1532 + x1535 + x1537)) + x693*(-x1
     &551 - 0.950793276165107d0*x1552 + x1553 + 1.42618991424766d0*x1556
     & - x1559 - 2.13928487137149d0*x1561 - x1564 - 2.13928487137149d0* 
     &x1567)))
      x1656 = 1.05383986550676d0*x165*x202*x297*x784
      x1657 = 1.26460783860812d0*x782*x785
      x1658 = 0.499760307959092d0*x1526*x285*x297/x203**6
      x1659 = x1464 - x1655 + x1656 - x1657 - x1658
      x1660 = x1659*x200
      x1661 = x1502*x1515 + 2.78d0*x1520 + x1521*x1522 + 2.78d0*x1660
      x1662 = x323*x792
      x1663 = x833*x835
      x1664 = x833*x834
      x1665 = x310*x553
      x1666 = 1.08969365000574d0*x1496
      x1667 = 1.08969365000574d0*x835
      x1668 = x1667*x320
      x1669 = gamma_cc*x1668
      x1670 = 0.0790517663015599d0*gamma_cc
      x1671 = x1670*x795
      x1672 = 0.00267822075231276d0*x821
      x1673 = 0.0338793284149542d0*x807
      x1674 = x1673*x799
      x1675 = x1086*x52
      x1676 = x1675*x823
      x1677 = x806*x808
      x1678 = 0.00105333172662086d0*x807
      x1679 = x1677*x1678
      x1680 = 3.56861714960788d-5*x828
      x1681 = x1680*x827
      x1682 = 0.0338793284149542d0*x313
      x1683 = gamma_cc*x1682
      x1684 = x792*x831
      x1685 = x1086*x199
      x1686 = gamma_cc*x321
      x1687 = x1686*x311/x317**3
      x1688 = x1687*x833
      x1689 = 0.131752943835933d0*x1315*x312
      x1690 = x1670*x541*x799
      x1691 = x1322*x986
      x1692 = x1328*x986
      x1693 = 0.002457774028782d0*x1488*x1677
      x1694 = x805**2
      x1695 = exp(x1352*x302)
      x1696 = x1694*x1695
      x1697 = 3.27488110955689d-5*x311*x807*x826
      x1698 = x1696*x1697
      x1699 = x303*x810
      x1700 = x1694*x1699
      x1701 = x1665*x305
      x1702 = x1678*x1701*x806
      x1703 = x1671*x793
      x1704 = x1493*x315
      x1705 = x1491*x1704
      x1706 = x1674*x793
      x1707 = x315*x790
      x1708 = x1497*x1707
      x1709 = x1500*x1707
      x1710 = x1679*x793
      x1711 = x1364*x302 + x1370*x302 + x1372*x804 + x44*(x1423 + x1655 
     &- x1656 + x1657 + x1658)
      x1712 = x1699*x1711
      x1713 = x1661*x812
      x1714 = gamma_cc*x995
      x1715 = x306*x819
      x1716 = x1429*x824
      x1717 = x1436*x301
      x1718 = 1.66426158945912d-6*x828/x304**4
      x1719 = x1000*x303
      x1720 = x301*x827
      x1721 = x52*x823
      x1722 = x1493*x308
      x1723 = x308*x831
      x1724 = x67/x324**2
      x1725 = x1724*x839
      x1726 = x839*x840
      x1727 = x1440*x325 + x1447*x325 + x1450*x1726 + x1659 + x1725*(x63
     &5 - x639 - x794 - x818 - x838) + x841*(x1089*(-0.00758829213155282
     &d0*x1065* x51*x820 - x1432*x1716 - x1435*x998 - x1497*x1723 - x150
     &0*x1723 + x1672*x1684 - x1676*x1684 + x1681*x1684 - x1689 + x1690 
     &- x1691 - x1692 - x1693 + x1694*x1719 - x1696*x1718 - x1698 + x170
     &0 + x1702 + x1703 - x1705 - x1706 - x1708 - x1709 + x1710 + x1711*
     &x1719 + x1712 - x1713 + 0.00267822075231276d0*x1715*x823 - 0.00016
     &6535466981701d0*x1717*x827 + 3.56861714960788d-5*x1720* x1721 - x1
     &722*x831**2 - x832*(x1515*x1522 + 5.56d0*x1520 + 11.12d0*x1521*x64
     &0 + 5.56d0*x1660)) + 4.23769752780008d0*x1315* x632 + x1322*x637 +
     & x1328*x637 - x1487*x1490 - x1489*x632*x793 - x1489*x815 - x1490*x
     &1663 + x1491*x1494 + x1495*x1496*x793 + x1495 *x1501 + x1497*x1498
     & + x1498*x1500 + x1501*x1667*x833 + x1661* x1662 + x1664*x1665*x16
     &66 + x1669*x793*x833 + x1688*(x1671 + x1672 - x1674 - x1676 + x167
     &9 + x1681 - x1683*x793 - x1684*x1685 ) + x1714*(x1689 - x1690 + x1
     &691 + x1692 + x1693 + x1698 - x1700 - x1702 - x1703 + x1705 + x170
     &6 + x1708 + x1709 - x1710 - x1712 + x1713) + 1.08969365000574d0*x7
     &93*x807*x815)
      x1728 = 0.0621813817393098d0*x42
      x1729 = x1503 - x1506 + x1509 + x1512
      x1730 = x1450*x325 + x1505 + x1507 - x1510 + x1511 + x1513 + x1726
     &*x1728 + x1729 - 2.0d0*x783 - 1.26460783860812d0*x786
      x1731 = x1235 + x1236 + x1237 + x1238
      x1732 = x394*x844
      x1733 = x1732*x385
      x1734 = x1221*x844
      x1735 = x1227*(-x1223 - x1224 - x1225 - x1226)
      x1736 = x445 - x850 - x853
      x1737 = x1250*x849
      x1738 = x410*x443
      x1739 = x416*x849
      x1740 = x152*x406
      x1741 = x436*x852
      x1742 = x437*x852
      x1743 = x435*x852
      x1744 = x397*x415
      x1745 = x397*x849
      x1746 = x403*x846
      x1747 = x131*x427*x846
      x1748 = x1211*x1746 - x1268 - x1269 - x1270 - x1271 + x1272*x1747
      x1749 = x1740*x1746
      x1750 = x130*x1746
      x1751 = x1738*x1750
      x1752 = x1740*x404
      x1753 = 39.9271408596337d0*x1230 - x1234 + x1731*x845 - 39.9271408
     &596337d0*x1733 - x1734*x401 + x1735*x401
      x1754 = x1753*x847
      x1755 = x1218*x1738
      x1756 = x1207*x846
      x1757 = x1211*x1756*x428
      x1758 = x405*x846
      x1759 = x137*x423
      x1760 = 8.83241096619832d0*x1759
      x1761 = x135*x1747
      x1762 = x1253*x1761
      x1763 = 0.00777267271741372d0*x447
      x1764 = -x1286 + x1288*x844 + x129*(-0.620681196888019d0*x1230 + x
     &1281*x1736 - x1282*x1734 + x1282*x1735 + x1283 + x1731*x396 + 0.62
     &0681196888019d0*x1733 + x448*(x1251*x152 + x1252*x1738 + x1265*x17
     &45 + x1266*x1740 + x1267*x1738 + x1274*x1748 - x1275 - x1276 - x12
     &77 - x152*x1737 - x1738*x1739 - x1738*x1742 - x1740* x1741 + x1743
     &*x1744 + x417*(x1263 + x1749 + x1751 - x1752 - x1754 - x1755 + x17
     &57 - x1758) + x438*(x1249 + x1254*x1760 + 1.78982066702156d0*x1255
     & - x1258*x1756*x402 + x1260 + x1261 + x1262 - x1749 - x1751 + x175
     &2 + x1753*x851 + x1754 + x1755 - x1757 + x1758 - x1760*x1761 - 1.7
     &8982066702156d0*x1762 + x432* x846))) - x1763*x446 + x1763*x855 + 
     &x843
      x1765 = x179*(x1205 + x1764)
      x1766 = x160*x183
      x1767 = x453*x875
      x1768 = x479*x858
      x1769 = x456*x862
      x1770 = x176*(-3.48d0*x1304 - 8.0d0*x1769 - 81.36d0*x455*x861)
      x1771 = x1294*x873
      x1772 = x476*x865
      x1773 = x1772*x475
      x1774 = -x1296*x869 - x1297*x871
      x1775 = x1299*x1774
      x1776 = 0.243797614528847d0*x465
      x1777 = x462*x871
      x1778 = x473*x50
      x1779 = 3.11111111111111d0*x1304
      x1780 = x477*(gamma_cs*x1302 + x173*x1776*x871 + x1776*x869 - x177
     &7*x463 - x1778 *x868 - x872*(-x1306*x1779 - x1308*x1779 + x469 - x
     &471))
      x1781 = x1770 + x1771 + x1773 + x1775 + x1780
      x1782 = x1781*x480
      x1783 = x178*x858
      x1784 = x159*x875
      x1785 = gamma_cc*x542
      x1786 = x61*x891
      x1787 = 0.444444444444444d0*x1304
      x1788 = -x1323*x1787 - x1325*x1787 + 1.33333333333333d0*x546 - 1.3
     &3333333333333d0*x549
      x1789 = x1327*x1788
      x1790 = -2.0d0*x888 + 2.0d0*x889
      x1791 = gamma_cc*x1790
      x1792 = x1320*x1791
      x1793 = x594*x910
      x1794 = x1340*x553
      x1795 = x1340*x891
      x1796 = gamma_cc*x564
      x1797 = x1796*x579
      x1798 = x1796*x905
      x1799 = x560*x890
      x1800 = x587*x908
      x1801 = x47*x901
      x1802 = x589*x901
      x1803 = -x1330*x1799 - x1337*x1800 + x1339*x1802 + x1355*x1801
      x1804 = x1785*x1799
      x1805 = x1348*x1788
      x1806 = x1346*x1790
      x1807 = x1320*x1806
      x1808 = x1785*x1801*x7
      x1809 = x1345*x1801*x553
      x1810 = x1345*x576*x891
      x1811 = x1353*x574*x900
      x1812 = x1356*x1811
      x1813 = x578*x900
      x1814 = 0.333333333333333d0*x1304
      x1815 = x1323*x1814
      x1816 = x1325*x1814
      x1817 = x567*(-x1815 - x1816 + x546 - 1.0d0*x549)
      x1818 = x1363*x890
      x1819 = 0.854960467080683d0*x1304
      x1820 = -x1374*x1819 - x1376*x1819 - x487 + x489
      x1821 = x1820*x29
      x1822 = x490*x880
      x1823 = x521*x878
      x1824 = 16.0d0*x1769
      x1825 = x20*x570
      x1826 = x24*x878
      x1827 = 118.434298572867d0*x863
      x1828 = x1412*x878
      x1829 = x1820*x38
      x1830 = x1414*x880
      x1831 = x1416*x880
      x1832 = x37*x878
      x1833 = x1832*x521
      x1834 = x1824*x520
      x1835 = x536*x882
      x1836 = 0.00181753314498899d0*x1418*x1769
      x1837 = -x1828 - x1829 - x1830 - x1831 + x1833 + x1834 + x1835 - x
     &1836
      x1838 = -x1397*x1827 + x1821*x569 - x1822*x570 + x1823*x570 + x182
     &4*x1825 + x1826*x532*x569 + x1837
      x1839 = x45*(x1817*x40 + x1818*x40 + x1838*x44 + x568*x899 + x573*
     &x897)
      x1840 = x1359*x1839
      x1841 = 0.888888888888889d0*x1304
      x1842 = -x1323*x1841 - x1325*x1841 + 2.66666666666667d0*x546 - 2.6
     &6666666666667d0*x549
      x1843 = -3.33333333333333d0*x888 + 3.33333333333333d0*x889
      x1844 = x1431*x1843
      x1845 = x586*x908
      x1846 = x601*(x1017*(x1337*x1802*x823 + x1337*x1845*x590 - x1428*x
     &1802*x48 - x1437*x1811 + 1.38533762600663d0*x1800*x819 + x1804 - x
     &1805 - x1807 - x1808 + x1809 + x1810 - x1812 + x1813 + x1839*x909 
     &+ x1840 - x1842*x588 - x1844*x585 + x592*x900) + x1319*x1345*x905 
     &- x1332*x1785*x910 + x1333*x1795 + x1333*x1798 + x1342*x1803 + x13
     &45*x1786*x579 + x1426*(-x1804 + x1805 + x1807 + x1808 - x1809 - x1
     &810 + x1812 - x1813 - x1840) - x1785*x1786*x59 - x1785*x62* x905 +
     & x1789*x554 + x1792*x554 + x1793*x1794 + x1793*x1797)
      x1847 = x605*(x1815 + x1816 - 1.0d0*x546 + x549)
      x1848 = x1847*x66
      x1849 = x1451*x916
      x1850 = x600*x913
      x1851 = x1850*x606
      x1852 = -x893 - x907 - x912
      x1853 = x1443*x1852
      x1854 = 0.666666666666667d0*x888 - 0.666666666666667d0*x889
      x1855 = x1446*x1854
      x1856 = x1855*x66
      x1857 = x1504*x1824
      x1858 = x1823*x525
      x1859 = x1822*x525
      x1860 = x1457*x1827
      x1861 = x1519*x1826
      x1862 = x1821*x517 + x1828 + x1829 + x1830 + x1831 - x1833 - x1834
     & - x1835 + x1836 + x1857 + x1858 - x1859 - x1860 + x1861
      x1863 = x1846 + x1848 + x1849 + x1851 + x1853 + x1856 + x1862
      x1864 = x165*x1770
      x1865 = x185*x1863 + x609*x922 + x615*x919 + x71*(-x1158*x866 - x1
     &158*x874 + x1472*x1771 + x1472*x1773 + x1472*x1775 + x1472*x1780 +
     & x183* x1864)
      x1866 = x160*x188
      x1867 = x1863*x621 + x609*x925 + x623*x919 + x71*(-x1479*x866 - x1
     &479*x874 + x1481*x1771 + x1481*x1773 + x1481*x1775 + x1481*x1780 +
     & x1864* x188)
      x1868 = x626*x927
      x1869 = x1838 - x1846 - x1848 - x1849 - x1851 - x1853 - x1856
      x1870 = 0.156233111018196d0*x630
      x1871 = x636*x891
      x1872 = gamma_cc*x1871
      x1873 = x316*x633
      x1874 = x791*x887
      x1875 = x1499*(x934 + x980)
      x1876 = x790*x983
      x1877 = gamma_cc*x994
      x1878 = x1877*x816
      x1879 = gamma_cc*x814*x816
      x1880 = x1003*x835
      x1881 = gamma_cc*x1873
      x1882 = x1003*x834
      x1883 = x1665*x1686
      x1884 = 8.0d0*x863
      x1885 = x1514*(x1504*x1884 - x1884*x520 - 2.0d0*x879 - 2.0d0*x883)
     &
      x1886 = x643*(x1517*x1820 + x1837 - x1857 - x1858 + x1859 + x1860 
     &- x1861)
      x1887 = x642*x981
      x1888 = x1521*x934
      x1889 = x1304*x211
      x1890 = x1889*x205
      x1891 = rho_s*x651
      x1892 = x947 - x948
      x1893 = 0.333333333333333d0*x1304
      x1894 = x1533*x1893
      x1895 = x1536*x1893
      x1896 = x1538*x956
      x1897 = x1540*x1896
      x1898 = x943 + x954
      x1899 = -x946 - x951 + x953
      x1900 = x1304*x280
      x1901 = rho_s*x702
      x1902 = x692*x952
      x1903 = x681*x942
      x1904 = -x1548*x1903 - 4.27856974274298d0*x950
      x1905 = 0.767414216267966d0*x1903*x495
      x1906 = x231*x949
      x1907 = x1555*x1906
      x1908 = 0.111111111111111d0*x1304
      x1909 = x682*(-x1533*x1908 - x1536*x1908 - 0.333333333333333d0*x65
     &5 + 0.333333333333333d0*x658)
      x1910 = x1558*x1906
      x1911 = 0.222222222222222d0*x1304
      x1912 = x667*(-x1533*x1911 - x1536*x1911 - 0.666666666666667d0*x65
     &5 + 0.666666666666667d0*x658)
      x1913 = x670*x944
      x1914 = x1560*x1913
      x1915 = -1.33333333333333d0*x938 + 1.33333333333333d0*x939
      x1916 = 0.00786976946662483d0*x1304
      x1917 = 11.4849360007528d0*x755
      x1918 = x1304*x289
      x1919 = x165*x777
      x1920 = 0.135575697997297d0*x1919
      x1921 = 0.150031049774234d0*x1304
      x1922 = 0.0271505358982603d0*x1304
      x1923 = x251*x965
      x1924 = x1923*x248
      x1925 = x270*x967
      x1926 = x1925*x248
      x1927 = x263*x728
      x1928 = x1617*x969
      x1929 = 0.203363546995945d0*x777
      x1930 = x1635*x962
      x1931 = x263*x750*x962
      x1932 = x1582*x970
      x1933 = 0.135575697997297d0*x1304*x777
      x1934 = x1587*x963
      x1935 = x252*x728
      x1936 = x1923*x743
      x1937 = x1925*x751
      x1938 = x1624*x252 + x1626*x252
      x1939 = -x1631*x263 - x1632*x263
      x1940 = x252*x742*x969
      x1941 = x1644*x965
      x1942 = 7.65662400050185d0*x1304
      x1943 = x1648*x967
      x1944 = -x1584*(-x1578*x252 + x1583*x1916) + x1590*x252 - x1591*x2
     &63 - x1593* x1935 + x1595*x1927 + 45.9397440030111d0*x1603*x1620 +
     & x1603* x1928*x1929 + 45.9397440030111d0*x1610*x1615 + x1610*x1929
     &*x1930 - x1613*x1924 - x1614*x1926 - x1618*x960 - x1621*x1935 + x1
     &622* x1927 - x1630*x1938 - x1634*x1939 + x1636*x960 - x1641*x1942 
     &- x1643*x1942 - x1645*x960 + x1649*x960 + x1917*x754 - x1917*x759 
     &+ x1918*x753 + x1918*x758 + x1920*x964 - x1920*x971 + x1931*x749* 
     &x961 + x1932*x1933 + x1933*x1934 - x1936*x728 + x1936*x761 - x1937
     &*x728 - x1937*x762 - x1940*x741*x961 + x1941*x757 - x1943* x757 - 
     &x726*(x1585*x263 + x1588*x1916) - x966*(x1583*x1921 - x1598*x252 -
     & x1600*x252 + x1606*x1922) - x968*(x1588*x1921 + x1607*x263 + x160
     &8*x263 + x1611*x1922)
      x1945 = 2.22222222222222d0*x1304
      x1946 = x1304*x291
      x1947 = x204*(x1530*x1892 + x1541*x940 + x1542*x1898 + x1545*x1899
     & + x1862*x237 - 0.0193686196838775d0*x1890*x214 + 0.00082186077399
     &2031d0*x1890* x648 + 0.00528235082287567d0*x1891*x205*x936 + x1897
     &*x660 + x230* (x1572*x1891 - 0.267332245098321d0*x1889*x208 + 0.01
     &38644007334803d0*x1889*x210 + x277*(rho_s*x1575*x165 - x1304* x712
     & + 4.83324359125373d0*x1304*x716 + x1944) - x709*x973) - x236 *x66
     &3*x887 + x278*(x1862*x284 - 0.0364124473051363d0*x1900*x214 + 0.00
     &993066744685536d0*x1901*x936 - x283*x705*x887) + x285*(x1862* x296
     & - 0.502576923213027d0*x1900*x208 + 0.167525641071009d0*x1901 *x95
     &7 + x295*(x1654*x1946 + x1944 - 62.364433435532d0*x1946*x247 + x28
     &9*(-x1945*x216 - x1945*x217 - x770 + x772) + 46.773325076649d0*x75
     &6*x775 + x774*x976) - x766*x978 - x767*x887 ) + x661*(x1894 + x189
     &5 + x655 - x658) + x697*(x1546*x955 + x1550 *x1904 + x1902*x691 + 
     &x229*(x1569*x1915 - x1905 - 1.54043890936153d0*x1907 + x1909 + x19
     &10 + 4.62131672808458d0* x1912 - x1913*x230*x678 + 4.6213167280845
     &8d0*x1914 + x945*(-x1894 - x1895 - x655 + x658)) + x693*(x1905 + 0
     &.71309495712383d0*x1907 - x1909 - x1910 - 2.13928487137149d0*x1912
     & - 2.13928487137149d0* x1914)))
      x1948 = 0.632303919304059d0*x785*x979
      x1949 = x1862 - x1947 - x1948
      x1950 = x1949*x200
      x1951 = x1502*x1885 + x1502*x1887 + 2.78d0*x1886 + 2.78d0*x1888 + 
     &2.78d0*x1950
      x1952 = x310*x891
      x1953 = x1686*x1952
      x1954 = x798*x891
      x1955 = x1673*x1954
      x1956 = x1675*x1845
      x1957 = x808*x991
      x1958 = x1678*x1957
      x1959 = x826*x991
      x1960 = x1680*x1959
      x1961 = x1002*x792
      x1962 = -x1683*x984 - x1685*x1961 - x1955 - x1956 + x1958 + x1960
      x1963 = 0.0395258831507799d0*x1488*x1954
      x1964 = x1789*x986
      x1965 = x1792*x986
      x1966 = 0.001228887014391d0*x1488*x1957
      x1967 = 0.000526665863310429d0*x807
      x1968 = x1701*x1967*x991
      x1969 = x797*x984
      x1970 = x802*x984
      x1971 = x1952*x305
      x1972 = x1967*x1971*x806
      x1973 = x793*x988
      x1974 = x1707*x1874
      x1975 = x1707*x1875
      x1976 = x1695*x805*x990
      x1977 = x1697*x1976
      x1978 = x811*x990
      x1979 = x1704*x1876
      x1980 = x811*x984
      x1981 = x793*x992
      x1982 = x1817*x302 + x1818*x302 + x44*(x1838 + x1947 + x1948) + x5
     &68*x989 + x804 *x897
      x1983 = x1699*x1982
      x1984 = x1951*x812
      x1985 = x1716*x1843
      x1986 = x1959*x301
      x1987 = x1845*x52
      x1988 = x1006*x840
      x1989 = -x1005 - x933 - x985 - x997
      x1990 = x1725*x1989 + x1726*x916 + x1847*x325 + x1855*x325 + x1949
     & + x1988*x606 + x841*(x1005*x793 + x1089*(x1001*x1684 - x1002*x172
     &2*x831 - x1684*x999 + 0.00133911037615638d0*x1715*x1845 - 8.326773
     &34908506d-5*x1717*x1959 - x1718*x1976 + x1719*x1982 + 1.7843085748
     &0394d-5*x1720*x1987 + 1.78430857480394d-5*x1721*x1986 - x1723*x187
     &4 - x1723*x1875 - x1842*x998 + x1961*x822 - x1961* x825 + x1961*x8
     &30 + x1963 - x1964 - x1965 - x1966 + x1968 + x1969 - x1970 + x1972
     & - x1973 - x1974 - x1975 - x1977 + x1978 - x1979 + x1980 + x1981 +
     & x1983 - x1984 - x1985*x585 + x830*x990 - x832*( x1522*x1885 + x15
     &22*x1887 + 5.56d0*x1886 + 5.56d0*x1888 + 5.56d0* x1950)) + x1487*x
     &1878 + x1494*x1876 + x1498*x1874 + x1498*x1875 + x1662*x1951 + x16
     &63*x1878 + x1664*x1953 + x1688*x1962 + x1714*( -x1963 + x1964 + x1
     &965 + x1966 - x1968 - x1969 + x1970 - x1972 + x1973 + x1974 + x197
     &5 + x1977 - x1978 + x1979 - x1980 - x1981 - x1983 + x1984) + x1789
     &*x637 + x1792*x637 + x1871*x1879 - x1872* x1873 - x1877*x319*x633 
     &+ x1879*x1880 - x1880*x1881 + x1882*x1883 - x635*x984 + x639*x984 
     &+ x793*x933 + x793*x997 + x818*x984 + x838*x984)
      x1991 = x858*x875
      x1992 = x846**2
      x1993 = x1207*x1992
      x1994 = x1212*x1993
      x1995 = x1992*x847
      x1996 = x1221*x844**2
      x1997 = x1735*x844
      x1998 = x1240 + 79.8542817192673d0*x1733 - x1996*x400 + x1997*x400
     &
      x1999 = x1998*x847
      x2000 = x1214*x1746
      x2001 = x1217*x1750
      x2002 = x1273*x852
      x2003 = x1280*x855
      x2004 = -x1285*x1732 + x129*(x1284 - 1.24136239377604d0*x1733 + x1
     &736*x2003 - x1996*x392 + x1997*x392 + x448*(x1210*x1737 + x1211*x1
     &743*x849 + x1214*x1741 + x1217*x1739 + x1217*x1742 + x1278 + x1748
     &*x2002 + x417*(x1249 + x1994 - x1995 - x1999 - x2000 - x2001) + x4
     &38*( x1257*x1761 - x1258*x1993 + x1264 + 3.57964133404311d0*x1762 
     &+ x1992*x851 - x1994 + x1995 + x1998*x851 + x1999 + x2000 + x2001 
     &))) + x395*x844*(x1287 + 0.0155453454348274d0) + x449 - x856
      x2005 = x179*(x1205 + x2004)
      x2006 = x176*(67.8d0*x1291 + 6.0d0*x1292 + x860)
      x2007 = x1772*x873
      x2008 = 3.11111111111111d0*x165
      x2009 = x477*(-gamma_cc*x1302 - x1301*x868*x871 - x872*(x1306*x200
     &8 + x1308* x2008))
      x2010 = x1298*x873
      x2011 = x1774*x2010
      x2012 = x2006 + 2.0d0*x2007 + x2009 + x2011
      x2013 = x2012*x480
      x2014 = 0.444444444444444d0*x165
      x2015 = x1323*x2014 + x1325*x2014
      x2016 = x1327*x2015
      x2017 = x48*x890
      x2018 = x1791*x2017
      x2019 = x1330*x905
      x2020 = x1330*x891
      x2021 = x1341*x910
      x2022 = x1348*x2015
      x2023 = x1806*x2017
      x2024 = x900**2
      x2025 = x1353*x2024
      x2026 = x1356*x2025
      x2027 = x1359*x45
      x2028 = x2024*x2027
      x2029 = x1801*x2020
      x2030 = 0.333333333333333d0*x165
      x2031 = x1323*x2030
      x2032 = x1325*x2030
      x2033 = x567*(x2031 + x2032)
      x2034 = x1362*x890*x896
      x2035 = x1371*x896
      x2036 = 0.854960467080683d0*x165
      x2037 = x1374*x2036 + x1376*x2036
      x2038 = x2037*x29
      x2039 = 12.0d0*x1292
      x2040 = x1884*x878
      x2041 = x2037*x38
      x2042 = x2039*x520
      x2043 = x1832*x1884
      x2044 = -x2041 - x2042 - x2043
      x2045 = -x1825*x2039 + x2038*x569 - x2040*x570 + x2044
      x2046 = x2033*x40 + x2034*x40 + x2035*x899 + x2045*x44
      x2047 = x2027*x2046
      x2048 = 0.888888888888889d0*x165
      x2049 = x1323*x2048 + x1325*x2048
      x2050 = x45*x909
      x2051 = x601*(x1017*(-x1437*x2025 + x1438*x1802*x1845 - x1844*x908
     & - x2022 - x2023 + x2024*x2050 - x2026 + x2028 + x2029 + x2046*x20
     &50 + x2047 - x2049*x588) + x1426*(x2022 + x2023 + x2026 - x2028 - 
     &x2029 - x2047) + x1786*x2019 + x1793*x2020*x59 + x1803*x2021 + x20
     &16*x554 + x2018*x554 + x2019*x595*x910)
      x2052 = x605*(-x2031 - x2032)
      x2053 = x2052*x66
      x2054 = x1442*x913
      x2055 = x1852*x2054
      x2056 = x1445*x1854*x915
      x2057 = x2056*x66
      x2058 = x604*x915
      x2059 = 0.0621813817393098d0*x2058
      x2060 = x1850*x2059
      x2061 = x1504*x2039
      x2062 = x2040*x525
      x2063 = x2038*x517 + x2041 + x2042 + x2043 - x2061 - x2062
      x2064 = x2051 + x2053 + x2055 + x2057 + x2060 + x2063
      x2065 = x165*x2006
      x2066 = x1469*x2007
      x2067 = 2.0d0*x919
      x2068 = x185*x2064 + x2067*x922 + x71*(x1472*x2009 + x1472*x2011 +
     & x183*x2065 + x2066*x611)
      x2069 = 0.00390625d0*x160*x194
      x2070 = x2064*x621 + x2067*x925 + x71*(x1481*x2009 + x1481*x2011 +
     & x188*x2065 + x2066*x622)
      x2071 = x2045 - x2051 - x2053 - x2055 - x2057 - x2060
      x2072 = x983**2
      x2073 = x323*x983
      x2074 = x643*(x1517*x2037 + x2044 + x2061 + x2062)
      x2075 = x1885*x934
      x2076 = x1887*x934
      x2077 = 0.333333333333333d0*x165
      x2078 = x1533*x2077
      x2079 = x1536*x2077
      x2080 = 0.111111111111111d0*x165
      x2081 = x682*(x1533*x2080 + x1536*x2080)
      x2082 = 0.222222222222222d0*x165
      x2083 = x667*(x1533*x2082 + x1536*x2082)
      x2084 = x1557*x1906*x942
      x2085 = x1913*x223*x949
      x2086 = x165*x289
      x2087 = 2.77555756156289d-17*x1919
      x2088 = 0.40672709399189d0*x1919
      x2089 = 68.9096160045166d0*x727
      x2090 = 0.40672709399189d0*x778
      x2091 = x1582*x1590 + x1587*x1591 - x1592*x1936*x252 + x1592*x1937
     &*x263 - x1629* x1924*x1938 - x1633*x1926*x1939 - x1637*x252*x754 -
     & x1637*x263* x759 - x1923*x1940*x2090 + x1925*x1931*x2090 - x1928*
     &x2088*x252 - x1930*x2088*x263 + x1932*x2087 + x1934*x2087 - x1941*
     &x2089 + x1943*x2089 - x2086*x753 - x2086*x758 - x966*(-x1582*x1598
     & - x1582*x1600) - x968*(-x1587*x1607 - x1587*x1608)
      x2092 = 2.22222222222222d0*x165
      x2093 = x204*(x1529*x1892*x940 + x1544*x1896*x1899 + 2.0d0*x1897*x
     &940 + x1898* x696*x956 + x2063*x237 + x230*(x277*(x2091 + x959) + 
     &x958) + x278 *(x2063*x284 + x941) + x285*(x2063*x296 + x295*(x2091
     & + x289*( x2092*x216 + x2092*x217) + x977) + x974) + x661*(-x2078 
     &- x2079) + x697*(x1549*x1904*x955 + 2.0d0*x1902*x955 + x229*(x1568
     &*x1915* x944 + x2081 + 4.62131672808458d0*x2083 + x2084 + 4.621316
     &72808458d0*x2085 + x945*(x2078 + x2079)) + x693*(-x2081 - 2.139284
     &87137149d0*x2083 - x2084 - 2.13928487137149d0*x2085)) + x937)
      x2094 = x2063 - x2093
      x2095 = x200*x2094
      x2096 = 2.78d0*x2074 + 2.78d0*x2075 + 5.56d0*x2076 + 2.78d0*x2095
      x2097 = x807*x994
      x2098 = 1.08969365000574d0*x2097
      x2099 = x1003*x1687
      x2100 = x2016*x986
      x2101 = x2018*x986
      x2102 = x990**2
      x2103 = x1695*x2102
      x2104 = x1697*x2103
      x2105 = x1699*x2102
      x2106 = x1678*x1971*x991
      x2107 = x1704*x2072
      x2108 = x1955*x984
      x2109 = x315*x983
      x2110 = x1874*x2109
      x2111 = x1875*x2109
      x2112 = x1958*x984
      x2113 = x2033*x302 + x2034*x302 + x2035*x989 + x44*(x2045 + x2093)
     &
      x2114 = x1699*x2113
      x2115 = x2096*x812
      x2116 = x1002*x308
      x2117 = x1006*x1724
      x2118 = x1988*x2059 + x1989*x2117 + x2052*x325 + x2056*x325 + x209
     &4 + x841*( x1003*x1667*x2097 + x1003*x1669*x984 + x1089*(-x1002**2
     &*x1722 - x1718*x2103 + x1719*x2102 + x1719*x2113 - x1874*x2116 - x
     &1875* x2116 - x1956*x1961 + x1960*x1961 - x1985*x908 + 3.568617149
     &60788d-5*x1986*x1987 - x2049*x998 - x2100 - x2101 - x2104 + x2105 
     &+ x2106 - x2107 - x2108 - x2110 - x2111 + x2112 + x2114 - x2115 - 
     &x832*(5.56d0*x2074 + 5.56d0*x2075 + 11.12d0*x2076 + 5.56d0*x2095))
     & + x1494*x2072 + x1662*x2096 + x1666*x1882*x1952 + x1714*(x2100 + 
     &x2101 + x2104 - x2105 - x2106 + x2107 + x2108 + x2110 + x2111 - x2
     &112 - x2114 + x2115) + x1871*x2098 + 1.08969365000574d0*x1872*x320
     &*x984 + x1874*x2073 + x1875*x2073 + x1962*x2099 + x2016*x637 + x20
     &18*x637 + x2098*x319*x984)
      x2119 = gamma_cc*x611
      x2120 = x1060*x462
      x2121 = x477*(x1060*x1776 - x166*x1778)
      x2122 = x1060*x1296
      x2123 = -x1294*x2120 - x1299*x2122 + x2121
      x2124 = x2123*x480
      x2125 = gamma_cc*x181
      x2126 = x21*x482
      x2127 = x2125*x2126
      x2128 = x1062*x453
      x2129 = 0.00244885996736565d0*x159
      x2130 = rho_c**(-8.66666666666667d0)*x1061*x477
      x2131 = 0.00777267271741372d0*x378
      x2132 = -x1042 - x1044 - x1050
      x2133 = x1201*x88
      x2134 = x1039*x2133
      x2135 = x112*x1184*x366
      x2136 = x115*x338
      x2137 = x342*x344
      x2138 = x1043*x346
      x2139 = 0.825248078914289d0*x349
      x2140 = x1043*x95
      x2141 = x2140*x356
      x2142 = x1047*x362
      x2143 = x357*x359
      x2144 = x104*x1045
      x2145 = 0.894910333510778d0*x2144
      x2146 = 4.41620548309916d0*x366
      x2147 = x109*x88
      x2148 = x2147*x91
      x2149 = x2148*x374
      x2150 = -x109*x2145 + x1130
      x2151 = x1051*x2131 + x94*(x1171*x2132 + x379*(-x1043*x340 + x1043
     &*x360 - x1049* x2136 - x1049*x2137*x95 + x1183*x2150 + x1198*(-x21
     &38 - x2139 - x2141) + x1204*(x1045*x109*x372 + x107*x2144*x2146 + 
     &x1201*x2145 + x2138 + x2139 + x2141) - 0.111863791688847d0*x2134 -
     & 0.552025685387395d0*x2135 + x2140*x375 + x2142*x2143 + 0.11186379
     &1688847d0*x2149 - 0.825248078914289d0*x343))
      x2152 = x1036*x1763
      x2153 = -x1025 - x1028 - x1035
      x2154 = 0.111863791688847d0*x1022
      x2155 = x1253*x149
      x2156 = x2154*x2155
      x2157 = x1243*x149*x423
      x2158 = 0.552025685387395d0*x2157
      x2159 = x1026*x440
      x2160 = 0.825248078914289d0*x442
      x2161 = x1034*x1740
      x2162 = x1034*x130*x1738
      x2163 = x1027*x403
      x2164 = x2163*x402
      x2165 = x1026*x407
      x2166 = 0.825248078914289d0*x411
      x2167 = -x2165 - x2166
      x2168 = x1027*x149
      x2169 = x1033*x435
      x2170 = x1031*x133
      x2171 = 0.894910333510778d0*x2170
      x2172 = x1253*x2171
      x2173 = 4.41620548309916d0*x1759
      x2174 = x2170*x2173
      x2175 = 0.894910333510778d0*x140
      x2176 = x1031*x2175
      x2177 = x2165 + x2166
      x2178 = x140*x2154
      x2179 = x1027*x436
      x2180 = -x1120 - x140*x2171
      x2181 = x129*(x1281*x2153 + x448*(x1265*x2178 + x1274*x2180 + x174
     &4*x2169 - x2156 - x2158 - x2159 - x2160 - x2161 - x2162 + x2168*x4
     &15 + x2179*x434 + x417*(-x2164 + x2167) + x438*(x2164 + x2172 + x2
     &174 + x2176*x429 + x2177))) + x2152
      x2182 = x179*(x2151 + x2181)
      x2183 = x1054*x479
      x2184 = x1054*x178
      x2185 = -x1766*x2183 - x183*x2182 + x2184*x484
      x2186 = x1019*x600
      x2187 = x2186*x606
      x2188 = -x1013 - x1018 - x64
      x2189 = x1443*x2188
      x2190 = gamma_cc*x465
      x2191 = x2190*x48
      x2192 = x551*x586
      x2193 = x1010*x1012*x50
      x2194 = x1016*x1332
      x2195 = x1016*x594
      x2196 = x1009*x1012*x54
      x2197 = gamma_cc*x50
      x2198 = 0.593716125431411d0*x2197
      x2199 = 1.18743225086282d0*x1081
      x2200 = -1.18743225086282d0*x1015 - x1329*x55
      x2201 = x601*(-0.692668813003313d0*x1010*x2191 + x1017*(x1014*x142
     &7 - x1336* x2198 + x2199*x590 + x558 - x562 + x577) + x1333*x2196 
     &+ x1342* x2200 + x1426*(-x558 + x562 - x577) - x1785*x2194 + x1794
     &*x2195 + x1797*x2195 + x2192*x2193 - x543 + x555 + x581 + x597)
      x2202 = x2187 + x2189 + x2201
      x2203 = x177*x2125
      x2204 = x1066*x1294
      x2205 = x1062*x1112
      x2206 = x1066*x1299
      x2207 = x1020*x615 + x1067*x609 + x185*x2202 + x71*(-x1466*x2203 +
     & x1472*x2121 + x2119*x610 + x2119*x614 - 0.00081628665578855d0*x22
     &04 + x2205 - 0.00102035831973569d0*x2206)
      x2208 = 0.017578125d0*x21
      x2209 = x190*x2208
      x2210 = 0.005859375d0*x1143
      x2211 = gamma_cc*x622
      x2212 = 0.00685680790862382d0*x2130
      x2213 = -x1866*x2183 - x188*x2182 + x2184*x620
      x2214 = x1062*x1108
      x2215 = x1020*x623 + x1076*x609 + x2202*x621 + x71*(-x1480*x2203 +
     & x1481*x2121 - 0.00228560263620794d0*x2204 - 0.00285700329525993d0
     &*x2206 + x2211 *x610 + x2211*x614 + x2214)
      x2216 = -x2187 - x2189 - x2201
      x2217 = x1084*x2197
      x2218 = x1088*x835
      x2219 = x1088*x834
      x2220 = x56*x809
      x2221 = x2220*x806
      x2222 = x314*x792
      x2223 = x2222*x790
      x2224 = x1081*x1082
      x2225 = 0.00922952226005862d0*x2224
      x2226 = x1086*x2197
      x2227 = 3.56861714960788d-5*x1081
      x2228 = -0.00114780889384833d0*x1081*x306 - x1682
      x2229 = x1091*x840
      x2230 = -x1085 - x1090 - x322
      x2231 = x1725*x2230 + x2229*x606 + x841*(-0.0215355519401368d0*x10
     &83*x2191 + x1085*x793 + x1089*(-x1087*x1684 + x1720*x2227 + 0.0026
     &7822075231276d0*x2190*x820 + x2221 - x2223 - x2226*x823 + x796 - x
     &801) + x1090*x793 + x1664*x2225 + x1688*x2228 + x1714*( -x2221 + x
     &2223 - x796 + x801) + x1879*x2218 - x1881*x2218 + x1883 *x2219 + x
     &2192*x2217 + x322*x793 - x634 + x638 + x817 + x837)
      x2232 = x477*(-rho_s*x1302 - x166*x1777)
      x2233 = -x1772*x2120 - x2010*x2122 + x2232
      x2234 = x2233*x480
      x2235 = x1062*x858
      x2236 = x2163*x846
      x2237 = x129*(x2003*x2153 + x448*(x1743*x2178 + x1745*x2169 + x200
     &2*x2180 + x2156 + x2158 + x2159 + x2160 + x2161 + x2162 + x2168*x8
     &49 + x2179*x852 + x417*(x2177 - x2236) + x438*(x1747*x2176 + x2167
     & - x2172 - x2174 + x2236))) - x2152
      x2238 = x179*(x2151 + x2237)
      x2239 = x1054*x875
      x2240 = -x1766*x2239 - x183*x2238
      x2241 = x2186*x916
      x2242 = x2054*x2188
      x2243 = x586*x890
      x2244 = x601*(x1017*(-x1800*x2198 + x1802*x2199 - x894 + x903) + x
     &1426*(x894 - x903) + x1793*x2196 + x1795*x2195 + x1798*x2195 + x20
     &21*x2200 + x2193*x2243 + x892 + x906 + x911)
      x2245 = x2241 + x2242 + x2244
      x2246 = x1066*x1772
      x2247 = x1066*x2010
      x2248 = x1020*x922 + x1067*x919 + x185*x2245 + x71*(x1472*x2232 + 
     &x2119*x920 + x2119*x921 - 0.00081628665578855d0*x2246 - 0.00102035
     &831973569d0* x2247)
      x2249 = -x1866*x2239 - x188*x2238
      x2250 = x1020*x925 + x1076*x919 + x2245*x621 + x71*(x1481*x2232 + 
     &x2211*x920 + x2211*x921 - 0.00228560263620794d0*x2246 - 0.00285700
     &329525993d0* x2247)
      x2251 = x187*x927
      x2252 = -x2241 - x2242 - x2244
      x2253 = x2220*x991
      x2254 = x2222*x983
      x2255 = x2117*x2230 + x2229*x916 + x841*(x1004 + x1085*x984 + x108
     &9*(-x1087* x1961 - x1845*x2226 + x1986*x2227 + x2253 - x2254 - x98
     &7) + x1090 *x984 + x1714*(-x2253 + x2254 + x987) + x1878*x2218 + x
     &1882*x2225 + x1953*x2219 + x2099*x2228 + x2217*x2243 + x322*x984 +
     & x932 + x996)
      x2256 = x1170*x94
      x2257 = x1051*x2256
      x2258 = x104*x1048*x1199*x88*1d0/x99
      x2259 = x1047*x1182
      x2260 = x1128*x95
      x2261 = x1049*x2260 + 0.223727583377694d0*x2142*x2148
      x2262 = 0.0559318958444236d0*x1040 + x2150*x2259 - 0.0374140139522
     &107d0*x2258 + x2261
      x2263 = x1052*x2262 + x2132*x2257
      x2264 = x1280*x129
      x2265 = x1036*x2264
      x2266 = x1023*x142
      x2267 = x1259*x133*x134*1d0/x138*x436
      x2268 = 0.223727583377694d0*x1022
      x2269 = x140*x2268
      x2270 = x2169*x2269
      x2271 = x1118*x130
      x2272 = x1034*x2271
      x2273 = x1033*x1273
      x2274 = x2180*x2273 + 0.0559318958444236d0*x2266 - 0.0374140139522
     &107d0*x2267 + x2270 + x2272
      x2275 = x1037*x2274 + x2153*x2265
      x2276 = x179*(x2263 + x2275)
      x2277 = -x183*x2276
      x2278 = x1055*x2125
      x2279 = 0.0032651466231542d0*x477
      x2280 = x169*x181
      x2281 = x170*x2280
      x2282 = x1059*x159
      x2283 = x2281*x2282
      x2284 = x1298*x169**2*x182
      x2285 = 5.33059123534764d-5*x2284
      x2286 = x159*x2285
      x2287 = rho_c**(-12.3333333333333d0)*x22
      x2288 = x1065*x477
      x2289 = 0.0032651466231542d0*x2288
      x2290 = rho_c**(-11.3333333333333d0)
      x2291 = x22*x2290
      x2292 = gamma_cc*rho_c**(-7.0d0)
      x2293 = x601*(0.593716125431411d0*x1011 + x1016*x1341*x2200 + 0.59
     &3716125431411d0*x1081*x2195*x54 + x1329*x2194 - 0.323484345894309d
     &0*x1430*x2292*x59*x594)
      x2294 = x1728*x2186
      x2295 = x1019*x1442*x2188
      x2296 = x1067*x2294 + x185*x2293 + x185*x2295 + x71*(x1064 - x2281
     &*x2289 + x2285 *x2291)
      x2297 = gamma_cc*x193
      x2298 = 0.01171875d0*x2297
      x2299 = 0.01171875d0*x1143
      x2300 = -x188*x2276
      x2301 = 0.00914241054483176d0*x477
      x2302 = 0.000149256554589734d0*x2284
      x2303 = x159*x2302
      x2304 = 0.00914241054483176d0*x2288
      x2305 = x1076*x2294 + x2293*x621 + x2295*x621 + x71*(x1075 - x2281
     &*x2304 + x2291 *x2302)
      x2306 = -x2293 - x2295
      x2307 = x1091*x1724*x2230 + x841*(0.0184590445201172d0*x1009*x1083
     & + x1088*x1668 + x1088*x1687*x2228 - 0.000312690015761657d0*x1429*
     &x2292*x300** 8.34d0*x305*x316*x834 + 0.0184590445201172d0*x2219*x2
     &224)
      x2308 = x1103*x626
      x2309 = rho_c**(-2.66666666666667d0)
      x2310 = 0.0522423459704672d0*x2309
      x2311 = x1294*x169
      x2312 = rho_c**(-3.66666666666667d0)
      x2313 = x169*x2312
      x2314 = x477*(-x2309*x473 + 0.139312922587913d0*x2313)
      x2315 = x1299*x169
      x2316 = 0.065302932463084d0*x2309
      x2317 = -x2310*x2311 + x2314 - x2315*x2316
      x2318 = x2317*x480
      x2319 = x1095*x453
      x2320 = x1109*x2288
      x2321 = x1100*x1294
      x2322 = x1100*x1299
      x2323 = -x1101*x609 + x71*(x1114 + x1472*x2314 - 0.000816286655788
     &55d0*x2321 - 0.00102035831973569d0*x2322)
      x2324 = 0.00685680790862382d0*x2320
      x2325 = x1095*x1108
      x2326 = x1074*x1100
      x2327 = -x2326*x609 + x71*(x1481*x2314 - 0.00228560263620794d0*x23
     &21 - 0.00285700329525993d0*x2322 + x2325)
      x2328 = x169*x1772
      x2329 = x477*x871
      x2330 = x169*x2010
      x2331 = -x2310*x2328 - x2310*x2329 - x2316*x2330
      x2332 = x2331*x480
      x2333 = x1095*x858
      x2334 = x1100*x1772
      x2335 = x181*x2329*x52
      x2336 = x1100*x2010
      x2337 = -x2326*x919 + x71*(-0.00228560263620794d0*x2334 - 0.002285
     &60263620794d0* x2335 - 0.00285700329525993d0*x2336)
      x2338 = -x1101*x919 + x71*(-0.00081628665578855d0*x2334 - 0.000816
     &28665578855d0* x2335 - 0.00102035831973569d0*x2336)
      x2339 = rho_c*x156
      x2340 = x1054*x1095
      x2341 = x159*x2190
      x2342 = x1113*x181
      x2343 = x1054*x1062
      x2344 = rho_c**(-10.3333333333333d0)
      x2345 = x166*x2344
      x2346 = x1099*x2125
      x2347 = rho_c**(-9.33333333333333d0)
      x2348 = x166*x2347
      x2349 = x2186*x42
      x2350 = x2349*x477
      x2351 = x1100*x2350
      x2352 = -2.53789160761462d-5*x2351 + x71*(-x1112*x2346 + x2285*x23
     &48)
      x2353 = x1108*x2280
      x2354 = -7.10609650132094d-5*x2351 + x71*(-x1108*x2346 + x2302*x23
     &48)
      x2355 = rho_c**(-8.33333333333333d0)
      x2356 = x159*x2355
      x2357 = rho_c**(-7.33333333333333d0)*x2285
      x2358 = x1098*x2357
      x2359 = x2357*x71
      x2360 = x1138*x626
      x2361 = -x1127 - x1129 - x1134
      x2362 = x1128*x346
      x2363 = 1.65049615782858d0*x349
      x2364 = x2260*x356
      x2365 = x1131*x362
      x2366 = 0.894910333510778d0*x104
      x2367 = -x102 - x2147*x2366
      x2368 = x1135*x2131 + x94*(x1171*x2361 + x379*(-x1128*x340 + x1128
     &*x360 - x1132* x2136 - x1133*x2137 + x1183*x2367 + x1198*(-x2362 -
     & x2363 - x2364 ) + x1204*(x104*x1172*x2146 + x2133*x2366 + x2147*x
     &372 + x2362 + x2363 + x2364) - 0.223727583377694d0*x2134 - 1.10405
     &137077479d0* x2135 + x2143*x2365 + 0.223727583377694d0*x2149 + x22
     &60*x375 - 1.65049615782858d0*x343))
      x2369 = x1125*x1763
      x2370 = x1117 + x1119 - x1124
      x2371 = x2155*x2268
      x2372 = 1.10405137077479d0*x2157
      x2373 = x1118*x440
      x2374 = 1.65049615782858d0*x442
      x2375 = x1122*x1740*x436
      x2376 = x1123*x1738
      x2377 = x2271*x403
      x2378 = x2377*x402
      x2379 = x1118*x407
      x2380 = 1.65049615782858d0*x411
      x2381 = x2379 + x2380
      x2382 = x1118*x416
      x2383 = x1122*x435
      x2384 = x1121*x133
      x2385 = 0.894910333510778d0*x2384
      x2386 = x1253*x2385
      x2387 = x2173*x2384
      x2388 = x1121*x2175
      x2389 = -x2379 - x2380
      x2390 = x2271*x436
      x2391 = -x140*x2385 + x146
      x2392 = x129*(x1281*x2370 + x448*(-x1265*x2269 + x1274*x2391 + x17
     &44*x2383 + x2371 + x2372 + x2373 + x2374 - x2375 - x2376 - x2382*x
     &415 - x2390*x434 + x417*(x2378 + x2381) + x438*(-x2378 + x2386 + x
     &2387 + x2388*x429 + x2389))) + x2369
      x2393 = x179*(x2368 + x2392)
      x2394 = x1136*x479
      x2395 = 0.104484691940934d0*x2313
      x2396 = rho_s*x2395
      x2397 = 0.104484691940934d0*rho_s*x2312
      x2398 = x477*(-0.38311053711676d0*rho_s*x1099 + x2397*x472)
      x2399 = 0.130605864926168d0*rho_s*x2313
      x2400 = x1294*x2396 + x1299*x2399 + x2398
      x2401 = x2400*x480
      x2402 = x1136*x178
      x2403 = x1110*x477
      x2404 = x2282*x2403
      x2405 = x1065*x453
      x2406 = rho_s*x1094
      x2407 = x2311*x2406
      x2408 = x1110*x2289
      x2409 = x2315*x2406
      x2410 = x1115*x609 + x71*(x1472*x2398 + 0.0016325733115771d0*x2407
     & - x2408 + 0.00204071663947138d0*x2409)
      x2411 = x1110*x2304
      x2412 = rho_s*x2325
      x2413 = x2412*x609 + x71*(x1481*x2398 + 0.00457120527241588d0*x240
     &7 + 0.00571400659051985d0*x2409 - x2411)
      x2414 = x2377*x846
      x2415 = x129*(x2003*x2370 + x448*(-x1743*x2269 + x1745*x2383 + x20
     &02*x2391 - x2371 - x2372 - x2373 - x2374 + x2375 + x2376 - x2382*x
     &849 - x2390*x852 + x417*(x2389 + x2414) + x438*(x1747*x2388 + x238
     &1 - x2386 - x2387 - x2414))) - x2369
      x2416 = x179*(x2368 + x2415)
      x2417 = x1136*x875
      x2418 = x477*(x2395 + x2397*x871)
      x2419 = x1772*x2396 + x2010*x2399 + x2418
      x2420 = x2419*x480
      x2421 = x1065*x858
      x2422 = x2328*x2406
      x2423 = x2330*x2406
      x2424 = x2412*x919 + x71*(x1481*x2418 + 0.00457120527241588d0*x242
     &2 + 0.00571400659051985d0*x2423)
      x2425 = x1115*x919 + x71*(x1472*x2418 + 0.0016325733115771d0*x2422
     & + 0.00204071663947138d0*x2423)
      x2426 = 0.000106611824706953d0*x159
      x2427 = x2284*x862
      x2428 = x2290*x2427
      x2429 = x1062*x1136
      x2430 = rho_s*x159
      x2431 = x169*x2430
      x2432 = x2125*x2431
      x2433 = x2148*x2365
      x2434 = x1041 + x1132*x2140 - 0.0748280279044215d0*x2258 + x2259*x
     &2367 + x2261 + 0.111863791688847d0*x2433
      x2435 = x1052*x2434 + x2257*x2361
      x2436 = x1122*x2179 + x2178*x2383 - 0.111863791688847d0*x2266 + 0.
     &0748280279044215d0*x2267 - x2270 - x2272 + x2273*x2391
      x2437 = x1037*x2436 + x2265*x2370
      x2438 = x179*(x2435 + x2437)
      x2439 = x1054*x1065
      x2440 = -x1139*x2439 - x183*x2438
      x2441 = rho_s*x2279
      x2442 = x2190*x2280
      x2443 = x2344*x2427
      x2444 = rho_s*x1095*x2350
      x2445 = 5.07578321522924d-5*x2444 + x71*(x2441*x2442 - 0.000106611
     &824706953d0* x2443)
      x2446 = 0.000298513109179468d0*x2428
      x2447 = -x1111*x2439 - x188*x2438
      x2448 = 0.000142121930026419d0*x2444 + x71*(rho_s*x2301*x2442 - 0.
     &000298513109179468d0*x2443)
      x2449 = rho_s*x2284
      x2450 = 0.000106611824706953d0*x2355*x2449
      x2451 = -x1098*x2450
      x2452 = x1095*x1136
      x2453 = -x2450*x71
      x2454 = x2347*x2449
      x2455 = x159*x2345
      x2456 = 0.000597026218358936d0*x2284
      x2457 = 0.000213223649413906d0*x2284
      x2458 = x2348*x2457
      x2459 = x1098*x2458
      x2460 = x1125*x2370
      x2461 = x1135*x2361
      x2462 = -0.447455166755389d0*x1022*x140*x2383 + x1122*x1273*x2391 
     &- x1123*x145 + 0.223727583377694d0*x2266 - 0.149656055808843d0*x22
     &67
      x2463 = x101*x1133 + x1126 + x1131*x1182*x2367 - 0.149656055808843
     &d0*x2258 + 0.447455166755389d0*x2433
      x2464 = x179*(x1037*x2462 + x1052*x2463 + x2256*x2461 + x2264*x246
     &0)
      x2465 = x2458*x71
      x2466 = x1071*x2126
      x2467 = 0.03125d0*x1143
      x2468 = x1071*x610
      x2469 = x1143*x177
      x2470 = -x1149*x609 + x71*(-x1145*x614 - 0.03125d0*x2468 + 0.0625d
     &0*x2469)
      x2471 = 0.0875d0*x1143
      x2472 = 0.0875d0*x1071
      x2473 = x184*x2472
      x2474 = -x2473*x609 + x71*(-0.0875d0*x2468 + 0.175d0*x2469 - x2472
     &*x614)
      x2475 = x1153*x160
      x2476 = 0.005859375d0*x2475
      x2477 = x1071*x920
      x2478 = -x2473*x919 + x71*(-x2472*x921 - 0.0875d0*x2477)
      x2479 = -x1149*x919 + x71*(-x1145*x921 - 0.03125d0*x2477)
      x2480 = x2184*x2467
      x2481 = 0.0625d0*x2297
      x2482 = x1071*x1112
      x2483 = x1060*x2282
      x2484 = x1060*x1065
      x2485 = x1071*x184*x2349
      x2486 = -0.000971584089676716d0*x2485 + x71*(-x184*x2481 + x2482*x
     &2484)
      x2487 = x1152*x51
      x2488 = 0.017578125d0*x2487
      x2489 = x2184*x2471
      x2490 = 0.175d0*x2297
      x2491 = x1071*x1108
      x2492 = -0.0027204354510948d0*x2485 + x71*(-x184*x2490 + x2484*x24
     &91)
      x2493 = x187*x2476
      x2494 = x159*x465
      x2495 = x169*x2491
      x2496 = x1099*x2482
      x2497 = x1098*x2496
      x2498 = x1071*x1113
      x2499 = x2496*x71
      x2500 = x1071*x2431
      x2501 = x1071*x169
      x2502 = x2441*x2501*x465
      x2503 = -x1098*x2502
      x2504 = -x2502*x71
      x2505 = x180*x2487
      x2506 = 0.09375d0*x2487
      x2507 = x1146*x2506
      x2508 = 0.01171875d0*x2475
      x2509 = x1098*x184*x2506
      x2510 = 0.0234375d0*x192/tau_c**5
      x2511 = x129*x71
      x2512 = x158 + x2511
      x2513 = x179*x2512
      x2514 = -x183*x2513 + x186
      x2515 = -x188*x2513 + x189
      x2516 = x2515*x626 + 1.0d0
      x2517 = x156*(x197 + x2514*x2516)
      x2518 = x2517 + x326
      x2519 = x129*x609
      x2520 = 0.0155453454348275d0*x1097
      x2521 = 0.0310906908696549d0*x16
      x2522 = x163*x520
      x2523 = 0.0098758948512826d0*x31
      x2524 = x2520 - x2521 + x2522 + x2523
      x2525 = x2519 + x2524 + x381
      x2526 = x179*x2525
      x2527 = x2512*x479
      x2528 = x178*x2512
      x2529 = -x1766*x2527 - x183*x2526 + x2528*x484 + x616
      x2530 = -x1866*x2527 - x188*x2526 + x2528*x620 + x624
      x2531 = -x2515*x619 + x2530*x626
      x2532 = x156*(x2514*x2531 + x2516*x2529 + x628)
      x2533 = x129*x919
      x2534 = -x2520 + x2521 - x2522 - x2523 + x2533 + x381
      x2535 = x179*x2534
      x2536 = x2512*x875
      x2537 = -x1766*x2536 - x183*x2535 + x923
      x2538 = -x1866*x2536 - x188*x2535 + x926
      x2539 = x2514*x626
      x2540 = x156*(x2516*x2537 + x2538*x2539 + x930)
      x2541 = x1007 + x2540
      x2542 = x2513*x611
      x2543 = x1020*x129
      x2544 = x1053 + x2543
      x2545 = x179*x2544
      x2546 = x1062*x2512
      x2547 = -gamma_cc*x2542 + x1058*x2546 + x1068 - x183*x2545
      x2548 = x2513*x622
      x2549 = -gamma_cc*x2548 + x1074*x2546 + x1077 - x188*x2545
      x2550 = x2210*x2515 + x2549*x626
      x2551 = x156*(x1021 + x2514*x2550 + x2516*x2547)
      x2552 = x1092 + x2551
      x2553 = x378*x94
      x2554 = x1051*x2553
      x2555 = x179*x182
      x2556 = 0.000680108862773701d0*x2555
      x2557 = x1095*x2512
      x2558 = x1074*x2557 + x1102 - x2554*x2556
      x2559 = 0.000242896022419179d0*x2555
      x2560 = x1058*x2557 + x1104 - x2554*x2559
      x2561 = x156*(x2516*x2560 + x2539*x2558)
      x2562 = x1135*x2553
      x2563 = x1065*x2512
      x2564 = -x1111*x2563 + x1116 - x2556*x2562
      x2565 = -x1139*x2563 + x1140 - x2559*x2562
      x2566 = x156*(x2516*x2565 + x2539*x2564)
      x2567 = x1147 + x2467*x2528
      x2568 = x1150 + x2471*x2528
      x2569 = -x2476*x2515 + x2568*x626
      x2570 = x156*(x2514*x2569 + x2516*x2567)
      x2571 = x1156*x2528
      x2572 = x2525*x479
      x2573 = x129*x1465
      x2574 = 0.000681574929370872d0*x8
      x2575 = x2574 + 0.0098758948512826d0
      x2576 = 4.0d0*x520
      x2577 = x2575*x28
      x2578 = 0.5d0*x519
      x2579 = 8.0409921107287d0*x1508
      x2580 = x1504*x458
      x2581 = 14.8042873216083d0*x531*x534
      x2582 = x2578 + x2579 + x2580 + x2581
      x2583 = x2575*x518 + x2575*x532*x533 - x2576*x457 + x2577*x526 + x
     &2582 + 1.0d0* x491 + x511*(-0.00412165868588285d0*x8 - 0.031090690
     &8696549d0) + 1.0d0*x515 + x538 + x608
      x2584 = x179*(x1205 + x2573 + x2583)
      x2585 = x1310*x2512
      x2586 = x178*x2525
      x2587 = x129*x1863
      x2588 = 0.00137388622862762d0*x502
      x2589 = 0.0155453454348275d0*x42
      x2590 = x2589*x600
      x2591 = x2590*x599
      x2592 = 0.0155453454348275d0*x66
      x2593 = x1449*x2592
      x2594 = x458*x520
      x2595 = x1414*x163
      x2596 = x1416*x163
      x2597 = 0.000227191643123624d0*x30*x492
      x2598 = x2575*x882 - x2577*x884
      x2599 = x163*x1832 + x1850*x2589 + x2058*x2592 + x2598 + x520*x864
     &
      x2600 = x179*(x1205 - x2578 - x2579 - x2580 - x2581 + x2587 - x258
     &8 - x2591 - x2593 + x2594 - x2595 - x2596 + x2597 + x2599)
      x2601 = x2525*x875
      x2602 = x2534*x479
      x2603 = x1781*x2512
      x2604 = x178*x2534
      x2605 = x2538*x626
      x2606 = x2534*x875
      x2607 = x129*x2064
      x2608 = -x2574 - 0.0098758948512826d0
      x2609 = 0.5d0*x883
      x2610 = 1.0d0*x879
      x2611 = x2576*x863
      x2612 = x1504*x864
      x2613 = x179*(x1205 + x2607 - x2608*x28*x884 + x2608*x882 - x2609 
     &- x2610 - x2611 + x2612 + x929)
      x2614 = x2012*x2512
      x2615 = x160*x2527
      x2616 = x129*x2202
      x2617 = x1019*x2590
      x2618 = x179*(x2151 + x2616 + x2617)
      x2619 = x2544*x479
      x2620 = x2123*x2512
      x2621 = x21*x2528
      x2622 = x2125*x2621
      x2623 = x178*x2544
      x2624 = x1062*x2525
      x2625 = 0.00244885996736565d0*x2130
      x2626 = x2208*x2515
      x2627 = x160*x2536
      x2628 = x129*x2245
      x2629 = x179*(x2151 - x2617 + x2628)
      x2630 = x2544*x875
      x2631 = x2233*x2512
      x2632 = x1062*x2534
      x2633 = x2514*x2538
      x2634 = x2125*x2545
      x2635 = x129*x2293 + x129*x2295
      x2636 = x179*(x2263 + x2635)
      x2637 = x1059*x2512
      x2638 = x2281*x2637
      x2639 = x2287*x2512
      x2640 = x160*x2515
      x2641 = x2558*x626
      x2642 = x179*x2151
      x2643 = -x183*x2642
      x2644 = x2317*x2512
      x2645 = x178*x483
      x2646 = x2554*x2645
      x2647 = x1095*x2525
      x2648 = 0.00244885996736565d0*x2320
      x2649 = x160*x182
      x2650 = x2554*x2649
      x2651 = x2650*x479
      x2652 = -x188*x2642
      x2653 = x2331*x2512
      x2654 = x1095*x2534
      x2655 = x2650*x875
      x2656 = x179*x2263
      x2657 = x1095*x2544
      x2658 = x2190*x2512
      x2659 = x179*x2125
      x2660 = 0.000485792044838358d0*x2554
      x2661 = x2345*x2512
      x2662 = x1062*x477
      x2663 = x2554*x2662
      x2664 = 0.0013602177255474d0*x2659
      x2665 = x1095*x477
      x2666 = 7.10609650132094d-5*x2665
      x2667 = x2355*x2512
      x2668 = x1169*x94
      x2669 = x1051*x2668
      x2670 = x2132*x2669
      x2671 = x2262*x2553
      x2672 = 2.53789160761462d-5*x2554
      x2673 = x2564*x626
      x2674 = x179*x2368
      x2675 = -x183*x2674
      x2676 = x2400*x2512
      x2677 = x2562*x2645
      x2678 = x2562*x2649
      x2679 = x2678*x479
      x2680 = x2403*x2637
      x2681 = x1065*x2525
      x2682 = -x188*x2674
      x2683 = x2419*x2512
      x2684 = x2678*x875
      x2685 = x1065*x2534
      x2686 = x179*x2435
      x2687 = 0.000485792044838358d0*x2562
      x2688 = 0.000106611824706953d0*x2512
      x2689 = x1065*x2544
      x2690 = rho_s*x2512
      x2691 = x2289*x2690
      x2692 = x169*x2125
      x2693 = x2562*x2662
      x2694 = x2304*x2690
      x2695 = x2562*x2665
      x2696 = rho_s*x2320
      x2697 = 7.10609650132094d-5*x2696
      x2698 = 0.000298513109179468d0*x2454
      x2699 = x2361*x2669
      x2700 = x2434*x2553
      x2701 = x2562*x2696
      x2702 = x2461*x2668
      x2703 = x2463*x2553
      x2704 = x1071*x2621
      x2705 = x1060*x2637
      x2706 = x2476*x2514
      x2707 = x1143*x178
      x2708 = 0.0013602177255474d0*x2707
      x2709 = x2512*x465
      x2710 = x2487*x2513
      x2711 = x71*x94
      x2712 = x157 + x2711
      x2713 = x179*x2712
      x2714 = -x183*x2713 + x186
      x2715 = -x188*x2713 + x189
      x2716 = x2715*x626 + 1.0d0
      x2717 = x156*(x197 + x2714*x2716)
      x2718 = x2717 + x326
      x2719 = x609*x94
      x2720 = x2524 + x2719 + x452
      x2721 = x179*x2720
      x2722 = x2712*x479
      x2723 = x178*x2712
      x2724 = -x1766*x2722 - x183*x2721 + x2723*x484 + x616
      x2725 = -x1866*x2722 - x188*x2721 + x2723*x620 + x624
      x2726 = -x2715*x619 + x2725*x626
      x2727 = x156*(x2714*x2726 + x2716*x2724 + x628)
      x2728 = x919*x94
      x2729 = x2524 + x2728 + x857
      x2730 = x179*x2729
      x2731 = x2712*x875
      x2732 = -x1766*x2731 - x183*x2730 + x923
      x2733 = -x1866*x2731 - x188*x2730 + x926
      x2734 = x2714*x626
      x2735 = x156*(x2716*x2732 + x2733*x2734 + x930)
      x2736 = x1007 + x2735
      x2737 = x2713*x611
      x2738 = x1020*x94
      x2739 = x1038 + x2738
      x2740 = x179*x2739
      x2741 = x1062*x2712
      x2742 = -gamma_cc*x2737 + x1058*x2741 + x1068 - x183*x2740
      x2743 = x2713*x622
      x2744 = -gamma_cc*x2743 + x1074*x2741 + x1077 - x188*x2740
      x2745 = x2210*x2715 + x2744*x626
      x2746 = x156*(x1021 + x2714*x2745 + x2716*x2742)
      x2747 = x1092 + x2746
      x2748 = x129*x447
      x2749 = x1036*x2748
      x2750 = x1095*x2712
      x2751 = x1074*x2750 + x1102 - x2556*x2749
      x2752 = x1058*x2750 + x1104 - x2559*x2749
      x2753 = x156*(x2716*x2752 + x2734*x2751)
      x2754 = x1125*x2748
      x2755 = x1065*x2712
      x2756 = -x1111*x2755 + x1116 - x2556*x2754
      x2757 = -x1139*x2755 + x1140 - x2559*x2754
      x2758 = x156*(x2716*x2757 + x2734*x2756)
      x2759 = x1147 + x2467*x2723
      x2760 = x1150 + x2471*x2723
      x2761 = -x2476*x2715 + x2760*x626
      x2762 = x156*(x2714*x2761 + x2716*x2759)
      x2763 = x1156*x2723
      x2764 = x2720*x479
      x2765 = x1465*x94
      x2766 = x179*(x1289 + x2583 + x2765)
      x2767 = x1310*x2712
      x2768 = x178*x2720
      x2769 = x1863*x94
      x2770 = x179*(x1764 + x2582 + x2588 + x2591 + x2593 - x2594 + x259
     &5 + x2596 - x2597 + x2599 + x2769)
      x2771 = x2720*x875
      x2772 = x2729*x479
      x2773 = x1781*x2712
      x2774 = x178*x2729
      x2775 = x2733*x626
      x2776 = x2729*x875
      x2777 = x2064*x94
      x2778 = x2610 + x2611 + x918
      x2779 = x179*(x2004 + x2598 + x2609 - x2612 + x2777 + x2778)
      x2780 = x2012*x2712
      x2781 = x160*x2722
      x2782 = x2202*x94
      x2783 = x179*(x2181 + x2617 + x2782)
      x2784 = x2739*x479
      x2785 = x2123*x2712
      x2786 = x21*x2723
      x2787 = x2125*x2786
      x2788 = x178*x2739
      x2789 = x1062*x2720
      x2790 = x2208*x2715
      x2791 = x160*x2731
      x2792 = x2245*x94
      x2793 = x179*(x2237 + x2617 + x2792)
      x2794 = x2739*x875
      x2795 = x2233*x2712
      x2796 = x1062*x2729
      x2797 = x2714*x2733
      x2798 = x2125*x2740
      x2799 = x2293*x94 + x2295*x94
      x2800 = x179*(x2275 + x2799)
      x2801 = x1059*x2712
      x2802 = x2281*x2801
      x2803 = x2287*x2712
      x2804 = x160*x2715
      x2805 = x2751*x626
      x2806 = x179*x2181
      x2807 = x2317*x2712
      x2808 = x2645*x2749
      x2809 = x1095*x2720
      x2810 = x2649*x2749
      x2811 = x2810*x479
      x2812 = x179*x2237
      x2813 = x2331*x2712
      x2814 = x1095*x2729
      x2815 = x2810*x875
      x2816 = x179*x2275
      x2817 = x1095*x2739
      x2818 = x2190*x2712
      x2819 = 0.000485792044838358d0*x2659
      x2820 = x2345*x2712
      x2821 = x2662*x2749
      x2822 = x2355*x2712
      x2823 = x1279*x129
      x2824 = x1036*x2823
      x2825 = x2153*x2824
      x2826 = x2274*x2748
      x2827 = 2.53789160761462d-5*x2749
      x2828 = x2756*x626
      x2829 = x179*x2392
      x2830 = x2400*x2712
      x2831 = x2645*x2754
      x2832 = x2649*x2754
      x2833 = x2832*x479
      x2834 = x2403*x2801
      x2835 = x1065*x2720
      x2836 = x179*x2415
      x2837 = x2419*x2712
      x2838 = x2832*x875
      x2839 = x1065*x2729
      x2840 = x179*x2437
      x2841 = 0.000106611824706953d0*x2712
      x2842 = x1065*x2739
      x2843 = rho_s*x2712
      x2844 = x2692*x2843
      x2845 = x2662*x2754
      x2846 = x2665*x2754
      x2847 = x2370*x2824
      x2848 = x2436*x2748
      x2849 = x2696*x2754
      x2850 = x2460*x2823
      x2851 = x2462*x2748
      x2852 = x1071*x2786
      x2853 = x1060*x2801
      x2854 = x2476*x2714
      x2855 = x2712*x465
      x2856 = 0.000485792044838358d0*x2707
      x2857 = x2501*x2843
      x2858 = x2487*x2713
      x2859 = x2511 + x2711
      x2860 = x179*x2859
      x2861 = -x183*x2860 + x186
      x2862 = -x188*x2860 + x189
      x2863 = x2862*x626 + 1.0d0
      x2864 = x156*(x197 + x2861*x2863)
      x2865 = x2864 + x326
      x2866 = x2519 + x2719 + 1.0d0*x39 + x68 + x69
      x2867 = x179*x2866
      x2868 = x2859*x479
      x2869 = x178*x2859
      x2870 = -x1766*x2868 - x183*x2867 + x2869*x484 + x616
      x2871 = -x1866*x2868 - x188*x2867 + x2869*x620 + x624
      x2872 = -x2862*x619 + x2871*x626
      x2873 = x156*(x2861*x2872 + x2863*x2870 + x628)
      x2874 = x2533 + x2728
      x2875 = x179*x2874
      x2876 = x2859*x875
      x2877 = -x1766*x2876 - x183*x2875 + x923
      x2878 = -x1866*x2876 - x188*x2875 + x926
      x2879 = x2861*x626
      x2880 = x156*(x2863*x2877 + x2878*x2879 + x930)
      x2881 = x1007 + x2880
      x2882 = x2860*x611
      x2883 = x2543 + x2738
      x2884 = x179*x2883
      x2885 = x1062*x2859
      x2886 = -gamma_cc*x2882 + x1058*x2885 + x1068 - x183*x2884
      x2887 = x2860*x622
      x2888 = -gamma_cc*x2887 + x1074*x2885 + x1077 - x188*x2884
      x2889 = x2210*x2862 + x2888*x626
      x2890 = x156*(x1021 + x2861*x2889 + x2863*x2886)
      x2891 = x1092 + x2890
      x2892 = x1095*x2859
      x2893 = x1074*x2892 + x1102
      x2894 = x1058*x2892 + x1104
      x2895 = x156*(x2863*x2894 + x2879*x2893)
      x2896 = x1065*x2859
      x2897 = -x1111*x2896 + x1116
      x2898 = -x1139*x2896 + x1140
      x2899 = x156*(x2863*x2898 + x2879*x2897)
      x2900 = x1147 + x2467*x2869
      x2901 = x1150 + x2471*x2869
      x2902 = -x2476*x2862 + x2901*x626
      x2903 = x156*(x2861*x2902 + x2863*x2900)
      x2904 = x1156*x2869
      x2905 = x2866*x479
      x2906 = 8.0d0*x457
      x2907 = x179*(x1450*x66 + x1451*x1728 + x1504*x2906 + x1729 + x257
     &3 + x2765 - x2906*x520 + 2.0d0*x491 + 2.0d0*x515 + 2.0d0*x519)
      x2908 = x1310*x2859
      x2909 = x178*x2866
      x2910 = x179*(x2587 + x2769 + x2778 + x886)
      x2911 = x2866*x875
      x2912 = x2874*x479
      x2913 = x1781*x2859
      x2914 = x178*x2874
      x2915 = x2878*x626
      x2916 = x2874*x875
      x2917 = x179*(x2607 + x2777)
      x2918 = x2012*x2859
      x2919 = x160*x2868
      x2920 = x179*(x1020 + x2616 + x2782)
      x2921 = x2883*x479
      x2922 = x2123*x2859
      x2923 = x21*x2869
      x2924 = x2125*x2923
      x2925 = x178*x2883
      x2926 = x1062*x2866
      x2927 = x2208*x2862
      x2928 = x160*x2876
      x2929 = x179*(x2628 + x2792)
      x2930 = x2883*x875
      x2931 = x2233*x2859
      x2932 = x1062*x2874
      x2933 = x2861*x2878
      x2934 = x2125*x2884
      x2935 = x179*(x2635 + x2799)
      x2936 = x1059*x2859
      x2937 = x2281*x2936
      x2938 = x2287*x2859
      x2939 = x160*x2862
      x2940 = x2893*x626
      x2941 = x2317*x2859
      x2942 = x1095*x2866
      x2943 = x2331*x2859
      x2944 = x1095*x2874
      x2945 = x1095*x2883
      x2946 = x2190*x2859
      x2947 = x2345*x2859
      x2948 = x2355*x2859
      x2949 = x2897*x626
      x2950 = x2400*x2859
      x2951 = x2403*x2936
      x2952 = x1065*x2866
      x2953 = x2419*x2859
      x2954 = x1065*x2874
      x2955 = 0.000106611824706953d0*x2859
      x2956 = x1065*x2883
      x2957 = rho_s*x2859
      x2958 = x2692*x2957
      x2959 = x1071*x2923
      x2960 = x1060*x2936
      x2961 = x2476*x2861
      x2962 = x2859*x465
      x2963 = x2501*x2957
      x2964 = x2487*x2860
      PBE_value = x71
      PBE_beta_value = x118
      PBE_alpha_value = x155
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*x327
      d1E(1) = rho_c*(x629 + 0.156233111018196d0*x631 + x842) + x327
      d1E(2) = rho_c*x1008
      d1E(3) = rho_c*x1093
      d1E(4) = rho_c*x1106
      d1E(5) = rho_c*x1142
      d1E(6) = rho_c*x1155
      d2E(1) = rho_c*(x1484*x198 + x1485*x629 - x1486*x198 + x156*(x1483
     & + x187*(x1476* x190 - x1477*x625 + x626*(-0.525d0*x1157 - x1159*x
     &1479 - x1290* x188 - x1311*x188 + x1313*x1478 + x1314*x1478 + x148
     &2)) + x196*( -0.1875d0*x1157 - x1158*x1159 - x1290*x183 - x1311*x1
     &83 + x1312* x1313 + x1312*x1314 + x1475) + 2.0d0*x617*x627) + x172
     &7) + x1730 + 2.0d0*x629 + 0.312466222036391d0*x631
      d2E(2) = rho_c*(x156*(x1868*x617 + x1869 + x187*(-x619*x927 + x626
     &*(-x1765*x188 - x1767*x1866 - x1768*x1866 - x1782*x188 + x1783*x62
     &0 + x1784*x620 + x1867)) + x196*(-x1765*x183 - x1766*x1767 - x1766
     &*x1768 - x1782 *x183 + x1783*x484 + x1784*x484 + x1865) + x627*x92
     &4) + x1870* x931 + x1990) + x1008
      d2E(3) = rho_c*(x156*(x196*(-x1158*x1991 - x183*x2005 - x183*x2013
     & + x2068) + x2069*x924*x927 + x2071 + x928*(-x1479*x1991 - x188*x2
     &005 - x188* x2013 + x2070)) + x2118)
      d2E(4) = rho_c*(x1080*x1870 + x156*(x1069*x627 + x1079*x617 + x187
     &*(-x1071*x2209 - x1078*x619 + x2210*x625 + x626*(x1074*x2128 - x15
     &9*x2212 - x188 *x2124 + 0.2625d0*x2127 - x2211*x454 - x2211*x481 +
     & x2213 + x2215 )) + x196*(x1058*x2128 - x183*x2124 - x2119*x454 - 
     &x2119*x481 + 0.09375d0*x2127 - x2129*x2130 + x2185 + x2207) + x221
     &6) + x2231) + x1093
      d2E(5) = rho_c*(x156*(x1069*x1868 + x1079*x924 + x196*(x1058*x2235
     & - x183*x2234 - x2119*x859 - x2119*x876 + x2240 + x2248) + x2210*x
     &2251 + x2252 + x928*(x1074*x2235 - x188*x2234 - x2211*x859 - x2211
     &*x876 + x2249 + x2250)) + x2255)
      d2E(6) = rho_c*(x156*(2.0d0*x1069*x1079 + x187*(x1078*x2299 + x191
     &*x2298 + x626*( x1054*x2214 - x1073 - 0.175d0*x2278 + x2283*x2301 
     &- x2287*x2303 + x2300 + x2305)) + x196*(x1054*x2205 - x1057 + x227
     &7 - 0.0625d0* x2278 + x2279*x2283 - x2286*x2287 + x2296) + x2306) 
     &+ x2307)
      d2E(7) = rho_c*(x1106*x1870 + x156*(x1105*x627 + x187*(-x1103*x619
     & + x626*(x1074* x2319 - x159*x2324 - x188*x2318 + x2213 + x2327)) 
     &+ x196*(x1058* x2319 - x183*x2318 - x2129*x2320 + x2185 + x2323) +
     & x2308*x617)) + x1106
      d2E(8) = x2339*(x1105*x1868 + x196*(x1058*x2333 - x183*x2332 + x22
     &40 + x2338) + x2308*x924 + x928*(x1074*x2333 - x188*x2332 + x2249 
     &+ x2337))
      d2E(9) = x2339*(x1069*x2308 + x1079*x1105 + x187*(x1103*x2210 + x6
     &26*(-x1055* x2211 + x1074*x2340 + x1074*x2343 + x2300 - x2303*x234
     &5 + x2341* x2353 + x2354)) + x196*(-x1055*x2119 + x1058*x2340 + x1
     &058*x2343 + x2277 - x2286*x2345 + x2341*x2342 + x2352))
      d2E(10) = x2339*(x1103*x1105*x2069 + x196*(x1054*x1114 + x2277 - x
     &2285*x2356 + x2359) + x928*(x1054*x2325 + x2300 - x2302*x2356 + x2
     &358))
      d2E(11) = rho_c*(x1142*x1870 + x156*(x1141*x627 + x187*(-x1138*x61
     &9 + x626*(-x1111 *x2405 - x1866*x2394 - x188*x2393 - x188*x2401 + 
     &x2402*x620 + 0.0137136158172476d0*x2404 + x2413)) + x196*(-x1139*x
     &2405 - x1766 *x2394 - x183*x2393 - x183*x2401 + x2402*x484 + 0.004
     &8977199347313d0*x2404 + x2410) + x2360*x617)) + x1142
      d2E(12) = x2339*(x1141*x1868 + x196*(-x1139*x2421 - x1766*x2417 - 
     &x183*x2416 - x183*x2420 + x2425) + x2360*x924 + x928*(-x1111*x2421
     & - x1866* x2417 - x188*x2416 - x188*x2420 + x2424))
      d2E(13) = x2339*(x1069*x2360 + x1079*x1141 + x187*(x1138*x2210 + x
     &626*(x1074*x2429 - x1137*x2211 + x159*x2446 - x2304*x2432 + x2447 
     &+ x2448)) + x196 *(x1058*x2429 - x1137*x2119 - x2289*x2432 + x2426
     &*x2428 + x2440 + x2445))
      d2E(14) = x2339*(x1105*x2360 + x1141*x2308 + x196*(x1058*x2452 + x
     &2426*x2454 + x2440 + x2453) + x928*(x1074*x2452 + 0.00029851310917
     &9468d0*x2284 *x2347*x2430 + x2447 + x2451))
      d2E(15) = x2339*(x1138*x1141*x2069 + x196*(-x1136*x2408 - x183*x24
     &64 - x2455*x2457 + x2465) + x928*(-x1136*x2411 - x188*x2464 - x245
     &5*x2456 + x2459 ))
      d2E(16) = rho_c*(x1155*x1870 + x156*(x1148*x627 + x1154*x617 + x18
     &7*(-x1151*x619 + x1153*x2209 - x2476*x625 + x626*(x1313*x2471 + x1
     &314*x2471 - 0.2625d0*x2466 + x2474)) + x196*(x1313*x2467 + x1314*x
     &2467 - 0.09375d0*x2466 + x2470))) + x1155
      d2E(17) = x2339*(x1148*x1868 + x1154*x924 + x196*(x1783*x2467 + x1
     &784*x2467 + x2479) - x2251*x2476 + x928*(x1783*x2471 + x1784*x2471
     & + x2478))
      d2E(18) = x2339*(x1069*x1154 + x1079*x1148 + x187*(-x1078*x2476 + 
     &x1151*x2210 - x191*x2488 + x626*(x180*x2490 - x2483*x2491 + x2489 
     &+ x2492)) + x196*(x180*x2481 + x2480 - x2482*x2483 + x2486))
      d2E(19) = x2339*(-x1103*x2493 + x1105*x1154 + x1148*x2308 + x196*(
     &x2480 - x2494* x2498 + x2499) + x928*(x2489 - x2494*x2495 + x2497)
     &)
      d2E(20) = x2339*(-x1138*x2493 + x1141*x1154 + x1148*x2360 + x196*(
     &x2289*x2500 + x2402*x2467 + x2504) + x928*(x2304*x2500 + x2402*x24
     &71 + x2503))
      d2E(21) = x2339*(2.0d0*x1148*x1154 + x187*(-x1151*x2508 + x191*x25
     &10 + x626*( -0.2625d0*x2505 + x2509)) + x196*(-0.09375d0*x2505 + x
     &2507))
      end if
      if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*x2518
      d1E(1) = rho_c*(x1870*x2517 + x2532 + x842) + x2518
      d1E(2) = rho_c*x2541
      d1E(3) = rho_c*x2552
      d1E(4) = rho_c*x2561
      d1E(5) = rho_c*x2566
      d1E(6) = rho_c*x2570
      d2E(1) = rho_c*(x1484*x2517 + x1485*x2532 - x1486*x2517 + x156*(x1
     &483 + x2514*( x1476*x2515 - x1477*x2530 + x626*(x1478*x2527 + x147
     &8*x2586 - x1479*x2572 + x1482 - x1866*x2585 - x188*x2584 - 0.525d0
     &*x2571)) + x2516*(-x1158*x2572 + x1312*x2527 + x1312*x2586 + x1475
     & - x1766 *x2585 - x183*x2584 - 0.1875d0*x2571) + 2.0d0*x2529*x2531
     &) + x1727) + x1485*x2517 + x1730 + 2.0d0*x2532
      d2E(2) = rho_c*(x156*(x1869 + x2514*(-x2538*x619 + x626*(-x1866*x2
     &601 - x1866* x2602 - x1866*x2603 + x1867 - x188*x2600 + x2536*x620
     & + x2604* x620)) + x2516*(-x1766*x2601 - x1766*x2602 - x1766*x2603
     & - x183* x2600 + x1865 + x2536*x484 + x2604*x484) + x2529*x2605 + 
     &x2531* x2537) + x1870*x2540 + x1990) + x2541
      d2E(3) = rho_c*(x156*(x2069*x2537*x2538 + x2071 + x2516*(-x1158*x2
     &606 - x1766* x2614 - x183*x2613 + x2068) + x2539*(-x1479*x2606 - x
     &1866*x2614 - x188*x2613 + x2070)) + x2118)
      d2E(4) = rho_c*(x156*(x2216 + x2514*(-x1071*x2626 + x2210*x2530 - 
     &x2549*x619 + x626*(x1074*x2624 - x1866*x2619 - x1866*x2620 - x188*
     &x2618 - x2211*x2526 - x2211*x2615 - x2212*x2512 + x2215 + 0.2625d0
     &*x2622 + x2623*x620)) + x2516*(x1058*x2624 - x1766*x2619 - x1766*x
     &2620 - x183*x2618 - x2119*x2526 - x2119*x2615 + x2207 - x2512*x262
     &5 + 0.09375d0*x2622 + x2623*x484) + x2529*x2550 + x2531*x2547) + x
     &1870*x2551 + x2231) + x2552
      d2E(5) = rho_c*(x156*(x2210*x2633 + x2252 + x2516*(x1058*x2632 - x
     &1766*x2630 - x1766*x2631 - x183*x2629 - x2119*x2535 - x2119*x2627 
     &+ x2248) + x2537*x2550 + x2539*(x1074*x2632 - x1866*x2630 - x1866*
     &x2631 - x188*x2629 - x2211*x2535 - x2211*x2627 + x2250) + x2547*x2
     &605) + x2255)
      d2E(6) = rho_c*(x156*(x2306 + x2514*(x2298*x2640 + x2299*x2549 + x
     &626*(-x188* x2636 + x2214*x2544 + x2301*x2638 - x2302*x2639 + x230
     &5 - x2548 - 0.175d0*x2634)) + x2516*(-x183*x2636 + x2205*x2544 + x
     &2279*x2638 - x2285*x2639 + x2296 - x2542 - 0.0625d0*x2634) + 2.0d0
     &*x2547* x2550) + x2307)
      d2E(7) = rho_c*(x156*(x2514*(-x2558*x619 + x626*(x1074*x2647 - x18
     &66*x2644 - x2324*x2512 + x2327 + 0.0020403265883211d0*x2646 - 0.00
     &0680108862773701d0*x2651 + x2652)) + x2516*(x1058*x2647 - x1766*x2
     &644 + x2323 - x2512*x2648 + x2643 + 0.000728688067257537d0*x2646 -
     & 0.000242896022419179d0*x2651) + x2529*x2641 + x2531*x2560) + x187
     &0*x2561) + x2561
      d2E(8) = x2339*(x2516*(x1058*x2654 - x1766*x2653 + x2338 + x2643 -
     & 0.000242896022419179d0*x2655) + x2537*x2641 + x2539*(x1074*x2654 
     &- x1866*x2653 + x2337 + x2652 - 0.000680108862773701d0*x2655) + x2
     &560*x2605)
      d2E(9) = x2339*(x2514*(x2210*x2558 + x626*(x1074*x2657 - x188*x265
     &6 - x2302*x2661 + x2353*x2658 + x2354 - x2554*x2664 + 3.5530482506
     &6047d-5*x2663 )) + x2516*(x1058*x2657 - x183*x2656 - x2285*x2661 +
     & x2342*x2658 + x2352 - x2659*x2660 + 1.26894580380731d-5*x2663) + 
     &x2547*x2641 + x2550*x2560)
      d2E(10) = x2339*(x2069*x2558*x2560 + x2516*(-x2285*x2667 + x2359 -
     & x2559*x2670 - x2559*x2671 + x2665*x2672) + x2539*(-x2302*x2667 + 
     &x2358 + x2554* x2666 - x2556*x2670 - x2556*x2671))
      d2E(11) = rho_c*(x156*(x2514*(-x2564*x619 + x626*(-x1111*x2681 - x
     &1866*x2676 + x2413 + 0.0020403265883211d0*x2677 - 0.00068010886277
     &3701d0*x2679 + 0.0137136158172476d0*x2680 + x2682)) + x2516*(-x113
     &9*x2681 - x1766*x2676 + x2410 + x2675 + 0.000728688067257537d0*x26
     &77 - 0.000242896022419179d0*x2679 + 0.0048977199347313d0*x2680) + 
     &x2529*x2673 + x2531*x2565) + x1870*x2566) + x2566
      d2E(12) = x2339*(x2516*(-x1139*x2685 - x1766*x2683 + x2425 + x2675
     & - 0.000242896022419179d0*x2684) + x2537*x2673 + x2539*(-x1111*x26
     &85 - x1866*x2683 + x2424 + x2682 - 0.000680108862773701d0*x2684) +
     & x2565*x2605)
      d2E(13) = x2339*(x2514*(x2210*x2564 + x626*(-x1111*x2689 - x188*x2
     &686 + x2446* x2512 + x2448 - x2562*x2664 - x2692*x2694 + 3.5530482
     &5066047d-5* x2693)) + x2516*(-x1139*x2689 - x183*x2686 + x2428*x26
     &88 + x2445 - x2659*x2687 - x2691*x2692 + 1.26894580380731d-5*x2693
     &) + x2547* x2673 + x2550*x2565)
      d2E(14) = x2339*(x2516*(x2453 + x2454*x2688 - x2559*x2699 - x2559*
     &x2700 - x2672* x2696 + 1.26894580380731d-5*x2695) + x2539*(x2451 +
     & x2512*x2698 - x2554*x2697 - x2556*x2699 - x2556*x2700 + 3.5530482
     &5066047d-5* x2695) + x2560*x2673 + x2565*x2641)
      d2E(15) = x2339*(x2069*x2564*x2565 + x2516*(-x2457*x2661 + x2465 -
     & x2559*x2702 - x2559*x2703 - 5.07578321522924d-5*x2701) + x2539*(-
     &x2456*x2661 + x2459 - x2556*x2702 - x2556*x2703 - 0.00014212193002
     &6419d0*x2701 ))
      d2E(16) = rho_c*(x156*(x2514*(x1153*x2626 - x2476*x2530 - x2568*x6
     &19 + x626*(x2471 *x2527 + x2471*x2586 + x2474 - 0.2625d0*x2704)) +
     & x2516*(x2467* x2527 + x2467*x2586 + x2470 - 0.09375d0*x2704) + x2
     &529*x2569 + x2531*x2567) + x1870*x2570) + x2570
      d2E(17) = x2339*(-x2476*x2633 + x2516*(x2467*x2536 + x2467*x2604 +
     & x2479) + x2537* x2569 + x2539*(x2471*x2536 + x2471*x2604 + x2478)
     & + x2567*x2605)
      d2E(18) = x2339*(x2514*(x2210*x2568 - x2476*x2549 - x2488*x2640 + 
     &x626*(x2471* x2623 + x2490*x2513 - x2491*x2705 + x2492)) + x2516*(
     &x2467*x2623 + x2481*x2513 - x2482*x2705 + x2486) + x2547*x2569 + x
     &2550*x2567)
      d2E(19) = x2339*(x2516*(-x2498*x2709 + x2499 + x2660*x2707) + x253
     &9*(-x2495*x2709 + x2497 + x2554*x2708) - x2558*x2706 + x2560*x2569
     & + x2567*x2641)
      d2E(20) = x2339*(x2516*(x2501*x2691 + x2504 + x2687*x2707) + x2539
     &*(x2501*x2694 + x2503 + x2562*x2708) - x2564*x2706 + x2565*x2569 +
     & x2567*x2673)
      d2E(21) = x2339*(x2514*(-x2508*x2568 + x2510*x2640 + x626*(x2509 -
     & 0.2625d0*x2710 )) + x2516*(x2507 - 0.09375d0*x2710) + 2.0d0*x2567
     &*x2569)
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PB
     &E_value)) then
      E = rho_c*x2718
      d1E(1) = rho_c*(x1870*x2717 + x2727 + x842) + x2718
      d1E(2) = rho_c*x2736
      d1E(3) = rho_c*x2747
      d1E(4) = rho_c*x2753
      d1E(5) = rho_c*x2758
      d1E(6) = rho_c*x2762
      d2E(1) = rho_c*(x1484*x2717 + x1485*x2727 - x1486*x2717 + x156*(x1
     &483 + x2714*( x1476*x2715 - x1477*x2725 + x626*(x1478*x2722 + x147
     &8*x2768 - x1479*x2764 + x1482 - x1866*x2767 - x188*x2766 - 0.525d0
     &*x2763)) + x2716*(-x1158*x2764 + x1312*x2722 + x1312*x2768 + x1475
     & - x1766 *x2767 - x183*x2766 - 0.1875d0*x2763) + 2.0d0*x2724*x2726
     &) + x1727) + x1485*x2717 + x1730 + 2.0d0*x2727
      d2E(2) = rho_c*(x156*(x1869 + x2714*(-x2733*x619 + x626*(-x1866*x2
     &771 - x1866* x2772 - x1866*x2773 + x1867 - x188*x2770 + x2731*x620
     & + x2774* x620)) + x2716*(-x1766*x2771 - x1766*x2772 - x1766*x2773
     & - x183* x2770 + x1865 + x2731*x484 + x2774*x484) + x2724*x2775 + 
     &x2726* x2732) + x1870*x2735 + x1990) + x2736
      d2E(3) = rho_c*(x156*(x2069*x2732*x2733 + x2071 + x2716*(-x1158*x2
     &776 - x1766* x2780 - x183*x2779 + x2068) + x2734*(-x1479*x2776 - x
     &1866*x2780 - x188*x2779 + x2070)) + x2118)
      d2E(4) = rho_c*(x156*(x2216 + x2714*(-x1071*x2790 + x2210*x2725 - 
     &x2744*x619 + x626*(x1074*x2789 - x1866*x2784 - x1866*x2785 - x188*
     &x2783 - x2211*x2721 - x2211*x2781 - x2212*x2712 + x2215 + 0.2625d0
     &*x2787 + x2788*x620)) + x2716*(x1058*x2789 - x1766*x2784 - x1766*x
     &2785 - x183*x2783 - x2119*x2721 - x2119*x2781 + x2207 - x2625*x271
     &2 + 0.09375d0*x2787 + x2788*x484) + x2724*x2745 + x2726*x2742) + x
     &1870*x2746 + x2231) + x2747
      d2E(5) = rho_c*(x156*(x2210*x2797 + x2252 + x2716*(x1058*x2796 - x
     &1766*x2794 - x1766*x2795 - x183*x2793 - x2119*x2730 - x2119*x2791 
     &+ x2248) + x2732*x2745 + x2734*(x1074*x2796 - x1866*x2794 - x1866*
     &x2795 - x188*x2793 - x2211*x2730 - x2211*x2791 + x2250) + x2742*x2
     &775) + x2255)
      d2E(6) = rho_c*(x156*(x2306 + x2714*(x2298*x2804 + x2299*x2744 + x
     &626*(-x188* x2800 + x2214*x2739 + x2301*x2802 - x2302*x2803 + x230
     &5 - x2743 - 0.175d0*x2798)) + x2716*(-x183*x2800 + x2205*x2739 + x
     &2279*x2802 - x2285*x2803 + x2296 - x2737 - 0.0625d0*x2798) + 2.0d0
     &*x2742* x2745) + x2307)
      d2E(7) = rho_c*(x156*(x2714*(-x2751*x619 + x626*(x1074*x2809 - x18
     &66*x2807 - x188 *x2806 - x2324*x2712 + x2327 + 0.0020403265883211d
     &0*x2808 - 0.000680108862773701d0*x2811)) + x2716*(x1058*x2809 - x1
     &766*x2807 - x183*x2806 + x2323 - x2648*x2712 + 0.00072868806725753
     &7d0*x2808 - 0.000242896022419179d0*x2811) + x2724*x2805 + x2726*x2
     &752) + x1870*x2753) + x2753
      d2E(8) = x2339*(x2716*(x1058*x2814 - x1766*x2813 - x183*x2812 + x2
     &338 - 0.000242896022419179d0*x2815) + x2732*x2805 + x2734*(x1074*x
     &2814 - x1866*x2813 - x188*x2812 + x2337 - 0.000680108862773701d0*x
     &2815 ) + x2752*x2775)
      d2E(9) = x2339*(x2714*(x2210*x2751 + x626*(x1074*x2817 - x188*x281
     &6 - x2302*x2820 + x2353*x2818 + x2354 - x2664*x2749 + 3.5530482506
     &6047d-5*x2821 )) + x2716*(x1058*x2817 - x183*x2816 - x2285*x2820 +
     & x2342*x2818 + x2352 - x2749*x2819 + 1.26894580380731d-5*x2821) + 
     &x2742*x2805 + x2745*x2752)
      d2E(10) = x2339*(x2069*x2751*x2752 + x2716*(-x2285*x2822 + x2359 -
     & x2559*x2825 - x2559*x2826 + x2665*x2827) + x2734*(-x2302*x2822 + 
     &x2358 - x2556* x2825 - x2556*x2826 + x2666*x2749))
      d2E(11) = rho_c*(x156*(x2714*(-x2756*x619 + x626*(-x1111*x2835 - x
     &1866*x2830 - x188*x2829 + x2413 + 0.0020403265883211d0*x2831 - 0.0
     &00680108862773701d0*x2833 + 0.0137136158172476d0*x2834)) + x2716*(
     &-x1139*x2835 - x1766*x2830 - x183*x2829 + x2410 + 0.00072868806725
     &7537d0*x2831 - 0.000242896022419179d0*x2833 + 0.0048977199347313d0
     &*x2834) + x2724*x2828 + x2726*x2757) + x1870* x2758) + x2758
      d2E(12) = x2339*(x2716*(-x1139*x2839 - x1766*x2837 - x183*x2836 + 
     &x2425 - 0.000242896022419179d0*x2838) + x2732*x2828 + x2734*(-x111
     &1*x2839 - x1866*x2837 - x188*x2836 + x2424 - 0.000680108862773701d
     &0*x2838 ) + x2757*x2775)
      d2E(13) = x2339*(x2714*(x2210*x2756 + x626*(-x1111*x2842 - x188*x2
     &840 - x2304* x2844 + x2446*x2712 + x2448 - x2664*x2754 + 3.5530482
     &5066047d-5* x2845)) + x2716*(-x1139*x2842 - x183*x2840 - x2289*x28
     &44 + x2428* x2841 + x2445 - x2754*x2819 + 1.26894580380731d-5*x284
     &5) + x2742* x2828 + x2745*x2757)
      d2E(14) = x2339*(x2716*(x2453 + x2454*x2841 - x2559*x2847 - x2559*
     &x2848 - x2696* x2827 + 1.26894580380731d-5*x2846) + x2734*(x2451 -
     & x2556*x2847 - x2556*x2848 - x2697*x2749 + x2698*x2712 + 3.5530482
     &5066047d-5* x2846) + x2752*x2828 + x2757*x2805)
      d2E(15) = x2339*(x2069*x2756*x2757 + x2716*(-x2457*x2820 + x2465 -
     & x2559*x2850 - x2559*x2851 - 5.07578321522924d-5*x2849) + x2734*(-
     &x2456*x2820 + x2459 - x2556*x2850 - x2556*x2851 - 0.00014212193002
     &6419d0*x2849 ))
      d2E(16) = rho_c*(x156*(x2714*(x1153*x2790 - x2476*x2725 - x2760*x6
     &19 + x626*(x2471 *x2722 + x2471*x2768 + x2474 - 0.2625d0*x2852)) +
     & x2716*(x2467* x2722 + x2467*x2768 + x2470 - 0.09375d0*x2852) + x2
     &724*x2761 + x2726*x2759) + x1870*x2762) + x2762
      d2E(17) = x2339*(-x2476*x2797 + x2716*(x2467*x2731 + x2467*x2774 +
     & x2479) + x2732* x2761 + x2734*(x2471*x2731 + x2471*x2774 + x2478)
     & + x2759*x2775)
      d2E(18) = x2339*(x2714*(x2210*x2760 - x2476*x2744 - x2488*x2804 + 
     &x626*(x2471* x2788 + x2490*x2713 - x2491*x2853 + x2492)) + x2716*(
     &x2467*x2788 + x2481*x2713 - x2482*x2853 + x2486) + x2742*x2761 + x
     &2745*x2759)
      d2E(19) = x2339*(x2716*(-x2498*x2855 + x2499 + x2749*x2856) + x273
     &4*(-x2495*x2855 + x2497 + x2708*x2749) - x2751*x2854 + x2752*x2761
     & + x2759*x2805)
      d2E(20) = x2339*(x2716*(x2289*x2857 + x2504 + x2754*x2856) + x2734
     &*(x2304*x2857 + x2503 + x2708*x2754) - x2756*x2854 + x2757*x2761 +
     & x2759*x2828)
      d2E(21) = x2339*(x2714*(-x2508*x2760 + x2510*x2804 + x626*(x2509 -
     & 0.2625d0*x2858 )) + x2716*(x2507 - 0.09375d0*x2858) + 2.0d0*x2759
     &*x2761)
      end if
      if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PB
     &E_value)) then
      E = rho_c*x2865
      d1E(1) = rho_c*(x1870*x2864 + x2873 + x842) + x2865
      d1E(2) = rho_c*x2881
      d1E(3) = rho_c*x2891
      d1E(4) = rho_c*x2895
      d1E(5) = rho_c*x2899
      d1E(6) = rho_c*x2903
      d2E(1) = rho_c*(x1484*x2864 + x1485*x2873 - x1486*x2864 + x156*(x1
     &483 + x2861*( x1476*x2862 - x1477*x2871 + x626*(x1478*x2868 + x147
     &8*x2909 - x1479*x2905 + x1482 - x1866*x2908 - x188*x2907 - 0.525d0
     &*x2904)) + x2863*(-x1158*x2905 + x1312*x2868 + x1312*x2909 + x1475
     & - x1766 *x2908 - x183*x2907 - 0.1875d0*x2904) + 2.0d0*x2870*x2872
     &) + x1727) + x1485*x2864 + x1730 + 2.0d0*x2873
      d2E(2) = rho_c*(x156*(x1869 + x2861*(-x2878*x619 + x626*(-x1866*x2
     &911 - x1866* x2912 - x1866*x2913 + x1867 - x188*x2910 + x2876*x620
     & + x2914* x620)) + x2863*(-x1766*x2911 - x1766*x2912 - x1766*x2913
     & - x183* x2910 + x1865 + x2876*x484 + x2914*x484) + x2870*x2915 + 
     &x2872* x2877) + x1870*x2880 + x1990) + x2881
      d2E(3) = rho_c*(x156*(x2069*x2877*x2878 + x2071 + x2863*(-x1158*x2
     &916 - x1766* x2918 - x183*x2917 + x2068) + x2879*(-x1479*x2916 - x
     &1866*x2918 - x188*x2917 + x2070)) + x2118)
      d2E(4) = rho_c*(x156*(x2216 + x2861*(-x1071*x2927 + x2210*x2871 - 
     &x2888*x619 + x626*(x1074*x2926 - x1866*x2921 - x1866*x2922 - x188*
     &x2920 - x2211*x2867 - x2211*x2919 - x2212*x2859 + x2215 + 0.2625d0
     &*x2924 + x2925*x620)) + x2863*(x1058*x2926 - x1766*x2921 - x1766*x
     &2922 - x183*x2920 - x2119*x2867 - x2119*x2919 + x2207 - x2625*x285
     &9 + 0.09375d0*x2924 + x2925*x484) + x2870*x2889 + x2872*x2886) + x
     &1870*x2890 + x2231) + x2891
      d2E(5) = rho_c*(x156*(x2210*x2933 + x2252 + x2863*(x1058*x2932 - x
     &1766*x2930 - x1766*x2931 - x183*x2929 - x2119*x2875 - x2119*x2928 
     &+ x2248) + x2877*x2889 + x2879*(x1074*x2932 - x1866*x2930 - x1866*
     &x2931 - x188*x2929 - x2211*x2875 - x2211*x2928 + x2250) + x2886*x2
     &915) + x2255)
      d2E(6) = rho_c*(x156*(x2306 + x2861*(x2298*x2939 + x2299*x2888 + x
     &626*(-x188* x2935 + x2214*x2883 + x2301*x2937 - x2302*x2938 + x230
     &5 - x2887 - 0.175d0*x2934)) + x2863*(-x183*x2935 + x2205*x2883 + x
     &2279*x2937 - x2285*x2938 + x2296 - x2882 - 0.0625d0*x2934) + 2.0d0
     &*x2886* x2889) + x2307)
      d2E(7) = rho_c*(x156*(x2861*(-x2893*x619 + x626*(x1074*x2942 - x18
     &66*x2941 - x2324*x2859 + x2327)) + x2863*(x1058*x2942 - x1766*x294
     &1 + x2323 - x2648*x2859) + x2870*x2940 + x2872*x2894) + x1870*x289
     &5) + x2895
      d2E(8) = x2339*(x2863*(x1058*x2944 - x1766*x2943 + x2338) + x2877*
     &x2940 + x2879*( x1074*x2944 - x1866*x2943 + x2337) + x2894*x2915)
      d2E(9) = x2339*(x2861*(x2210*x2893 + x626*(x1074*x2945 - x2302*x29
     &47 + x2353* x2946 + x2354)) + x2863*(x1058*x2945 - x2285*x2947 + x
     &2342*x2946 + x2352) + x2886*x2940 + x2889*x2894)
      d2E(10) = x2339*(x2069*x2893*x2894 + x2863*(-x2285*x2948 + x2359) 
     &+ x2879*(-x2302* x2948 + x2358))
      d2E(11) = rho_c*(x156*(x2861*(-x2897*x619 + x626*(-x1111*x2952 - x
     &1866*x2950 + x2413 + 0.0137136158172476d0*x2951)) + x2863*(-x1139*
     &x2952 - x1766*x2950 + x2410 + 0.0048977199347313d0*x2951) + x2870*
     &x2949 + x2872*x2898) + x1870*x2899) + x2899
      d2E(12) = x2339*(x2863*(-x1139*x2954 - x1766*x2953 + x2425) + x287
     &7*x2949 + x2879* (-x1111*x2954 - x1866*x2953 + x2424) + x2898*x291
     &5)
      d2E(13) = x2339*(x2861*(x2210*x2897 + x626*(-x1111*x2956 - x2304*x
     &2958 + x2446* x2859 + x2448)) + x2863*(-x1139*x2956 - x2289*x2958 
     &+ x2428*x2955 + x2445) + x2886*x2949 + x2889*x2898)
      d2E(14) = x2339*(x2863*(x2453 + x2454*x2955) + x2879*(x2451 + x269
     &8*x2859) + x2894 *x2949 + x2898*x2940)
      d2E(15) = x2339*(x2069*x2897*x2898 + x2863*(-x2457*x2947 + x2465) 
     &+ x2879*(-x2456* x2947 + x2459))
      d2E(16) = rho_c*(x156*(x2861*(x1153*x2927 - x2476*x2871 - x2901*x6
     &19 + x626*(x2471 *x2868 + x2471*x2909 + x2474 - 0.2625d0*x2959)) +
     & x2863*(x2467* x2868 + x2467*x2909 + x2470 - 0.09375d0*x2959) + x2
     &870*x2902 + x2872*x2900) + x1870*x2903) + x2903
      d2E(17) = x2339*(-x2476*x2933 + x2863*(x2467*x2876 + x2467*x2914 +
     & x2479) + x2877* x2902 + x2879*(x2471*x2876 + x2471*x2914 + x2478)
     & + x2900*x2915)
      d2E(18) = x2339*(x2861*(x2210*x2901 - x2476*x2888 - x2488*x2939 + 
     &x626*(x2471* x2925 + x2490*x2860 - x2491*x2960 + x2492)) + x2863*(
     &x2467*x2925 + x2481*x2860 - x2482*x2960 + x2486) + x2886*x2902 + x
     &2889*x2900)
      d2E(19) = x2339*(x2863*(-x2498*x2962 + x2499) + x2879*(-x2495*x296
     &2 + x2497) - x2893*x2961 + x2894*x2902 + x2900*x2940)
      d2E(20) = x2339*(x2863*(x2289*x2963 + x2504) + x2879*(x2304*x2963 
     &+ x2503) - x2897 *x2961 + x2898*x2902 + x2900*x2949)
      d2E(21) = x2339*(x2861*(-x2508*x2901 + x2510*x2939 + x626*(x2509 -
     & 0.2625d0*x2964 )) + x2863*(x2507 - 0.09375d0*x2964) + 2.0d0*x2900
     &*x2902)
      end if
      end subroutine
