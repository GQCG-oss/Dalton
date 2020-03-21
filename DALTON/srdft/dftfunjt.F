C File: srdft/dftfunjt.F


! hjaaj April 2018: search for TODO and fix them all !!!

C*****************************************************************************
      subroutine DESRC_MGGA(ESRC_MGGA, rho, grad2, tau, laplace, mu,
     & norder, E, d1E, d2E)
C*****************************************************************************
C   Finite difference code for spin-unpolarized meta-GGA correlation functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(1) : total density
C            grad2(1) : grad_c . grad_c
C            tau(1) : total kinetic energy density
C            laplace(1) : laplacian of total density
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRC_MGGA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRC_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: September 06, 2018
C*****************************************************************************
      implicit none
      external ESRC_MGGA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(2), grad2(3), tau(2), laplace(2)
      real*8, intent(out) :: E, d1E(9), d2E(45)
      
      real*8 :: rho_c, rho_cp, rho_cm, gamma_cc, gamma_ccp, gamma_ccm,
     & tau_c, tau_cp, tau_cm, laplace_c, laplace_cp, laplace_cm, Ep_1,
     & Em_1, Epp_3_1, Emm_3_1, Epm_3_1, Emp_3_1, Ep_3, Em_3, Epp_6_1,
     & Emm_6_1, Epm_6_1, Emp_6_1, Epp_6_3, Emm_6_3, Epm_6_3, Emp_6_3,
     & Ep_6, Em_6, Epp_8_1, Emm_8_1, Epm_8_1, Emp_8_1, Epp_8_3, Emm_8_3,
     & Epm_8_3, Emp_8_3, Epp_8_6, Emm_8_6, Epm_8_6, Emp_8_6, Ep_8,
     & Em_8, hr_c, hg_cc, ht_c, hl_c
      
      rho_c = rho(1)
      gamma_cc = grad2(1)
      tau_c = tau(1)
      laplace_c = laplace(1)
      ! Do checks of input
      if (gamma_cc .lt. 0.0d0) then
         call quit('DESRC_MGGA ERROR: gamma_cc negative')
      end if
      
      E = 0.0d0
      if (norder .ge. 1) then
         d1E(:) = 0.0d0
      end if
      if (norder .ge. 2) then
         d2E(:) = 0.0d0
      end if
      if (rho_c .lt. 1.d-10) then
         E = 0.0d0
         return
      end if
      
      hr_c = 1.d-4*rho_c
      if (gamma_cc .lt. 1.d-10) then
         hg_cc = 0.0d0
      else
         hg_cc = 1.d-4*gamma_cc
      end if
      if (tau_c .lt. 1.d-10) then
         ht_c = 0.0d0
      else
         ht_c = 1.d-4*tau_c
      end if
      if (laplace_c .lt. 1.d-10) then
         hl_c = 0.0d0
      else
         hl_c = 1.d-4*laplace_c
      end if
      
      ! Variation of variables
      rho_cp = rho_c + hr_c
      rho_cm = rho_c - hr_c
      gamma_ccp = gamma_cc + hg_cc
      gamma_ccm = gamma_cc - hg_cc
      tau_cp = tau_c + ht_c
      tau_cm = tau_c - ht_c
      laplace_cp = laplace_c + hl_c
      laplace_cm = laplace_c - hl_c
      
      call ESRC_MGGA(rho_c,gamma_cc,tau_c,laplace_c,mu,E)
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_c .gt. 0.0d0) then
            call ESRC_MGGA(rho_cp,gamma_cc,tau_c,laplace_c,mu,Ep_1)
            call ESRC_MGGA(rho_cm,gamma_cc,tau_c,laplace_c,mu,Em_1)
         end if
         if (hg_cc .gt. 0.0d0) then
            call ESRC_MGGA(rho_c,gamma_ccp,tau_c,laplace_c,mu,Ep_3)
            call ESRC_MGGA(rho_c,gamma_ccm,tau_c,laplace_c,mu,Em_3)
         end if
         if (ht_c .gt. 0.0d0) then
            call ESRC_MGGA(rho_c,gamma_cc,tau_cp,laplace_c,mu,Ep_6)
            call ESRC_MGGA(rho_c,gamma_cc,tau_cm,laplace_c,mu,Em_6)
         end if
         if (hl_c .gt. 0.0d0) then
            call ESRC_MGGA(rho_c,gamma_cc,tau_c,laplace_cp,mu,Ep_8)
            call ESRC_MGGA(rho_c,gamma_cc,tau_c,laplace_cm,mu,Em_8)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hg_cc .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_MGGA(rho_cp,gamma_ccp,tau_c,laplace_c,mu,
     &                     Epp_3_1)
               call ESRC_MGGA(rho_cm,gamma_ccp,tau_c,laplace_c,mu,
     &                     Epm_3_1)
               call ESRC_MGGA(rho_cp,gamma_ccm,tau_c,laplace_c,mu,
     &                     Emp_3_1)
               call ESRC_MGGA(rho_cm,gamma_ccm,tau_c,laplace_c,mu,
     &                     Emm_3_1)
            end if
         end if
         if (ht_c .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_MGGA(rho_cp,gamma_cc,tau_cp,laplace_c,mu,
     &                     Epp_6_1)
               call ESRC_MGGA(rho_cm,gamma_cc,tau_cp,laplace_c,mu,
     &                     Epm_6_1)
               call ESRC_MGGA(rho_cp,gamma_cc,tau_cm,laplace_c,mu,
     &                     Emp_6_1)
               call ESRC_MGGA(rho_cm,gamma_cc,tau_cm,laplace_c,mu,
     &                     Emm_6_1)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_MGGA(rho_c,gamma_ccp,tau_cp,laplace_c,mu,
     &                     Epp_6_3)
               call ESRC_MGGA(rho_c,gamma_ccm,tau_cp,laplace_c,mu,
     &                     Epm_6_3)
               call ESRC_MGGA(rho_c,gamma_ccp,tau_cm,laplace_c,mu,
     &                     Emp_6_3)
               call ESRC_MGGA(rho_c,gamma_ccm,tau_cm,laplace_c,mu,
     &                     Emm_6_3)
            end if
         end if
         if (hl_c .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_MGGA(rho_cp,gamma_cc,tau_c,laplace_cp,mu,
     &                     Epp_8_1)
               call ESRC_MGGA(rho_cm,gamma_cc,tau_c,laplace_cp,mu,
     &                     Epm_8_1)
               call ESRC_MGGA(rho_cp,gamma_cc,tau_c,laplace_cm,mu,
     &                     Emp_8_1)
               call ESRC_MGGA(rho_cm,gamma_cc,tau_c,laplace_cm,mu,
     &                     Emm_8_1)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_MGGA(rho_c,gamma_ccp,tau_c,laplace_cp,mu,
     &                     Epp_8_3)
               call ESRC_MGGA(rho_c,gamma_ccm,tau_c,laplace_cp,mu,
     &                     Epm_8_3)
               call ESRC_MGGA(rho_c,gamma_ccp,tau_c,laplace_cm,mu,
     &                     Emp_8_3)
               call ESRC_MGGA(rho_c,gamma_ccm,tau_c,laplace_cm,mu,
     &                     Emm_8_3)
            end if
            if (ht_c .gt. 0.0d0) then
               call ESRC_MGGA(rho_c,gamma_cc,tau_cp,laplace_cp,mu,
     &                     Epp_8_6)
               call ESRC_MGGA(rho_c,gamma_cc,tau_cm,laplace_cp,mu,
     &                     Epm_8_6)
               call ESRC_MGGA(rho_c,gamma_cc,tau_cp,laplace_cm,mu,
     &                     Emp_8_6)
               call ESRC_MGGA(rho_c,gamma_cc,tau_cm,laplace_cm,mu,
     &                     Emm_8_6)
            end if
         end if
      end if
      ! Calculate first derivatives
      if (norder .ge. 1) then
         ! d1E / drho_c
         if (hr_c .gt. 0.0d0) then
            d1E(1)=(Ep_1-Em_1)/(2.0d0*hr_c)
         end if
         ! d1E / dgamma_cc
         if (hg_cc .gt. 0.0d0) then
            d1E(3)=(Ep_3-Em_3)/(2.0d0*hg_cc)
         end if
         ! d1E / dtau_c
         if (ht_c .gt. 0.0d0) then
            d1E(6)=(Ep_6-Em_6)/(2.0d0*ht_c)
         end if
         ! d1E / dlaplace_c
         if (hl_c .gt. 0.0d0) then
            d1E(8)=(Ep_8-Em_8)/(2.0d0*hl_c)
         end if
      end if
      ! Calculate second derivatives
      if (norder .ge. 2) then
         if (hr_c .gt. 0.0d0) then
            ! d2E / (drho_c**2)
            d2E(1)=(Ep_1-2.0d0*E+Em_1)/(hr_c**2)
         end if
         if (hg_cc .gt. 0.0d0) then
            ! d2E / (dgamma_cc * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(4)=(Epp_3_1-Epm_3_1-Emp_3_1+Emm_3_1)/(4.0d0*hg_cc*
     &                     hr_c)
            end if
            ! d2E / (dgamma_cc**2)
            d2E(6)=(Ep_3-2.0d0*E+Em_3)/(hg_cc**2)
         end if
         if (ht_c .gt. 0.0d0) then
            ! d2E / (dtau_c * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(16)=(Epp_6_1-Epm_6_1-Emp_6_1+Emm_6_1)/(4.0d0*ht_c*
     &                     hr_c)
            end if
            ! d2E / (dtau_c * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(18)=(Epp_6_3-Epm_6_3-Emp_6_3+Emm_6_3)/(4.0d0*ht_c*
     &                     hg_cc)
            end if
            ! d2E / (dtau_c**2)
            d2E(21)=(Ep_6-2.0d0*E+Em_6)/(ht_c**2)
         end if
         if (hl_c .gt. 0.0d0) then
            ! d2E / (dlaplace_c * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(29)=(Epp_8_1-Epm_8_1-Emp_8_1+Emm_8_1)/(4.0d0*hl_c*
     &                     hr_c)
            end if
            ! d2E / (dlaplace_c * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(31)=(Epp_8_3-Epm_8_3-Emp_8_3+Emm_8_3)/(4.0d0*hl_c*
     &                     hg_cc)
            end if
            ! d2E / (dlaplace_c * dtau_c)
            if (ht_c .gt. 0.0d0) then
               d2E(34)=(Epp_8_6-Epm_8_6-Emp_8_6+Emm_8_6)/(4.0d0*hl_c*
     &                     ht_c)
            end if
            ! d2E / (dlaplace_c**2)
            d2E(36)=(Ep_8-2.0d0*E+Em_8)/(hl_c**2)
         end if
      end if
      
      end subroutine


C*****************************************************************************
      subroutine DESRC_SPINMGGA(ESRC_SPINMGGA, rho, grad2, tau, laplace,
     & mu, norder, E, d1E, d2E)
C*****************************************************************************
C   Finite difference code for spin-unpolarized meta-GGA correlation functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(1) : total density
C            rho(2) : total density
C            grad2(1) : grad_c . grad_c
C            grad2(2) : grad_s . grad_s
C            grad2(3) : grad_c . grad_s
C            tau(1) : total kinetic energy density
C            tau(2) : spin kinetic energy density
C            laplace(1) : laplacian of total density
C            laplace(2) : laplacian of spin density
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRC_MGGA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRC_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: September 06, 2018
C*****************************************************************************
      implicit none
      external ESRC_SPINMGGA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(2), grad2(3), tau(2), laplace(2)
      real*8, intent(out) :: E, d1E(9), d2E(45)
      
      real*8 :: rho_c, rho_cp, rho_cm, rho_s, rho_sp, rho_sm, gamma_cc,
     & gamma_ccp, gamma_ccm, gamma_ss, gamma_ssp, gamma_ssm, gamma_cs,
     & gamma_csp, gamma_csm, tau_c, tau_cp, tau_cm, tau_s, tau_sp,
     & tau_sm, laplace_c, laplace_cp, laplace_cm, laplace_s, laplace_sp,
     & laplace_sm, Ep_1, Em_1, Epp_2_1, Emm_2_1, Epm_2_1, Emp_2_1,
     & Ep_2, Em_2, Epp_3_1, Emm_3_1, Epm_3_1, Emp_3_1, Epp_3_2, Emm_3_2,
     & Epm_3_2, Emp_3_2, Ep_3, Em_3, Epp_4_1, Emm_4_1, Epm_4_1, Emp_4_1,
     & Epp_4_2, Emm_4_2, Epm_4_2, Emp_4_2, Epp_4_3, Emm_4_3, Epm_4_3,
     & Emp_4_3, Ep_4, Em_4, Epp_5_1, Emm_5_1, Epm_5_1, Emp_5_1, Epp_5_2,
     & Emm_5_2, Epm_5_2, Emp_5_2, Epp_5_3, Emm_5_3, Epm_5_3, Emp_5_3,
     & Epp_5_4, Emm_5_4, Epm_5_4, Emp_5_4, Ep_5, Em_5, Epp_6_1, Emm_6_1,
     & Epm_6_1, Emp_6_1, Epp_6_2, Emm_6_2, Epm_6_2, Emp_6_2, Epp_6_3,
     & Emm_6_3, Epm_6_3, Emp_6_3, Epp_6_4, Emm_6_4, Epm_6_4, Emp_6_4,
     & Epp_6_5, Emm_6_5, Epm_6_5, Emp_6_5, Ep_6, Em_6, Epp_7_1, Emm_7_1,
     & Epm_7_1, Emp_7_1, Epp_7_2, Emm_7_2, Epm_7_2, Emp_7_2, Epp_7_3,
     & Emm_7_3, Epm_7_3, Emp_7_3, Epp_7_4, Emm_7_4, Epm_7_4, Emp_7_4,
     & Epp_7_5, Emm_7_5, Epm_7_5, Emp_7_5, Epp_7_6, Emm_7_6, Epm_7_6,
     & Emp_7_6, Ep_7, Em_7, Epp_8_1, Emm_8_1, Epm_8_1, Emp_8_1, Epp_8_2,
     & Emm_8_2, Epm_8_2, Emp_8_2, Epp_8_3, Emm_8_3, Epm_8_3, Emp_8_3,
     & Epp_8_4, Emm_8_4, Epm_8_4, Emp_8_4, Epp_8_5, Emm_8_5, Epm_8_5,
     & Emp_8_5, Epp_8_6, Emm_8_6, Epm_8_6, Emp_8_6, Epp_8_7, Emm_8_7,
     & Epm_8_7, Emp_8_7, Ep_8, Em_8, Epp_9_1, Emm_9_1, Epm_9_1, Emp_9_1,
     & Epp_9_2, Emm_9_2, Epm_9_2, Emp_9_2, Epp_9_3, Emm_9_3, Epm_9_3,
     & Emp_9_3, Epp_9_4, Emm_9_4, Epm_9_4, Emp_9_4, Epp_9_5, Emm_9_5,
     & Epm_9_5, Emp_9_5, Epp_9_6, Emm_9_6, Epm_9_6, Emp_9_6, Epp_9_7,
     & Emm_9_7, Epm_9_7, Emp_9_7, Epp_9_8, Emm_9_8, Epm_9_8, Emp_9_8,
     & Ep_9, Em_9, hr_c, hr_s, hg_cc, hg_ss, hg_cs, ht_c, ht_s, hl_c,
     & hl_s
      
      rho_c = rho(1)
      rho_s = rho(2)
      gamma_cc = grad2(1)
      gamma_ss = grad2(2)
      gamma_cs = grad2(3)
      tau_c = tau(1)
      tau_s = tau(2)
      laplace_c = laplace(1)
      laplace_s = laplace(2)
      ! Do checks of input
      if (dabs(rho_s) .gt. rho_c) then
         call quit('DESRC_SPINMGGA ERROR: |rho_s| > rho_c')
      end if
      if (gamma_cc .lt. 0.0d0) then
         call quit('DESRC_SPINMGGA ERROR: gamma_cc negative')
      end if
      if (gamma_ss .lt. 0.0d0) then
         call quit('DESRC_SPINMGGA ERROR: gamma_ss negative')
      end if
      
      E = 0.0d0
      if (norder .ge. 1) then
         d1E(:) = 0.0d0
      end if
      if (norder .ge. 2) then
         d2E(:) = 0.0d0
      end if
      if (rho_c .lt. 1.d-10) then
         E = 0.0d0
         return
      end if
      
      hr_c = 1.d-4*rho_c
      hr_s = 1.d-4*rho_c
      if (gamma_cc .lt. 1.d-10) then
         hg_cc = 0.0d0
      else
         hg_cc = 1.d-4*gamma_cc
      end if
      if (gamma_ss .lt. 1.d-10) then
         hg_ss = 0.0d0
      else
         hg_ss = 1.d-4*gamma_ss
      end if
      ! abc(), gamma_cs can be negative
      if (dabs(gamma_cs) .lt. 1.d-10) then
         hg_cs = 0.0d0
      else
         hg_cs = 1.d-4*dabs(gamma_cs)
      end if
      if (tau_c .lt. 1.d-10) then
         ht_c = 0.0d0
      else
         ht_c = 1.d-4*tau_c
      end if
      if (tau_s .lt. 1.d-10) then
         ht_s = 0.0d0
      else
         ht_s = 1.d-4*tau_s
      end if
      if (laplace_c .lt. 1.d-10) then
         hl_c = 0.0d0
      else
         hl_c = 1.d-4*laplace_c
      end if
      if (laplace_s .lt. 1.d-10) then
         hl_s = 0.0d0
      else
         hl_s = 1.d-4*laplace_s
      end if
      
      ! Variation of variables
      rho_cp = rho_c + hr_c
      rho_cm = rho_c - hr_c
      rho_sp = rho_s + hr_s
      rho_sm = rho_s - hr_s
      gamma_ccp = gamma_cc + hg_cc
      gamma_ccm = gamma_cc - hg_cc
      gamma_ssp = gamma_ss + hg_ss
      gamma_ssm = gamma_ss - hg_ss
      gamma_csp = gamma_cs + hg_cs
      gamma_csm = gamma_cs - hg_cs
      tau_cp = tau_c + ht_c
      tau_cm = tau_c - ht_c
      tau_sp = tau_s + ht_s
      tau_sm = tau_s - ht_s
      laplace_cp = laplace_c + hl_c
      laplace_cm = laplace_c - hl_c
      laplace_sp = laplace_s + hl_s
      laplace_sm = laplace_s - hl_s
      
      call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,tau_c,
     &tau_s,laplace_c,laplace_s,mu,E)
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_c .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Ep_1)
            call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Em_1)
         end if
         if (hr_s .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Ep_2)
            call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Em_2)
         end if
         if (hg_cc .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Ep_3)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Em_3)
         end if
         if (hg_ss .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Ep_4)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Em_4)
         end if
         if (hg_cs .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_csp,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Ep_5)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_csm,
     &              tau_c,tau_s,laplace_c,laplace_s,mu,Em_5)
         end if
         if (ht_c .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_cp,tau_s,laplace_c,laplace_s,mu,Ep_6)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_cm,tau_s,laplace_c,laplace_s,mu,Em_6)
         end if
         if (ht_s .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_sp,laplace_c,laplace_s,mu,Ep_7)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_sm,laplace_c,laplace_s,mu,Em_7)
         end if
         if (hl_c .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_cp,laplace_s,mu,Ep_8)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_cm,laplace_s,mu,Em_8)
         end if
         if (hl_s .gt. 0.0d0) then
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_sp,mu,Ep_9)
            call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              tau_c,tau_s,laplace_c,laplace_sm,mu,Em_9)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hr_s .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_2_1)
               call ESRC_SPINMGGA(rho_cm,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_2_1)
               call ESRC_SPINMGGA(rho_cp,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_2_1)
               call ESRC_SPINMGGA(rho_cm,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_2_1)
            end if
         end if
         if (hg_cc .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_3_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_3_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_3_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_3_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_3_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_3_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_3_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_3_2)
            end if
         end if
         if (hg_ss .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_4_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_4_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_4_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_4_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_4_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_4_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_4_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_4_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_4_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_4_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_4_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_4_3)
            end if
         end if
         if (hg_cs .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_5_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_5_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_5_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_5_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_5_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_5_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_5_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_5_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_5_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_5_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_5_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_5_3)
            end if
            if (hg_ss .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_5_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_csp,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_5_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_5_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_csm,tau_c,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_5_4)
            end if
         end if
         if (ht_c .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_6_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_6_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_6_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_6_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_6_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_6_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_6_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_6_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_6_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_6_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_6_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_6_3)
            end if
            if (hg_ss .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epp_6_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_cp,tau_s,laplace_c,laplace_s,
     &                     mu,Epm_6_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emp_6_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_cm,tau_s,laplace_c,laplace_s,
     &                     mu,Emm_6_4)
            end if
            if (hg_cs .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_cp,tau_s,laplace_c,laplace_s,mu,Epp_6_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_cp,tau_s,laplace_c,laplace_s,mu,Epm_6_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_cm,tau_s,laplace_c,laplace_s,mu,Emp_6_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_cm,tau_s,laplace_c,laplace_s,mu,Emm_6_5)
            end if
         end if
         if (ht_s .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epp_7_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epm_7_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emp_7_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emm_7_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epp_7_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epm_7_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emp_7_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emm_7_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epp_7_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epm_7_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emp_7_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emm_7_3)
            end if
            if (hg_ss .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epp_7_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_sp,laplace_c,laplace_s,
     &                     mu,Epm_7_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emp_7_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_sm,laplace_c,laplace_s,
     &                     mu,Emm_7_4)
            end if
            if (hg_cs .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_c,tau_sp,laplace_c,laplace_s,mu,Epp_7_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_c,tau_sp,laplace_c,laplace_s,mu,Epm_7_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_c,tau_sm,laplace_c,laplace_s,mu,Emp_7_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_c,tau_sm,laplace_c,laplace_s,mu,Emm_7_5)
            end if
            if (ht_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cp,tau_sp,laplace_c,laplace_s,mu,Epp_7_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cm,tau_sp,laplace_c,laplace_s,mu,Epm_7_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cp,tau_sm,laplace_c,laplace_s,mu,Emp_7_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cm,tau_sm,laplace_c,laplace_s,mu,Emm_7_6)
            end if
         end if
         if (hl_c .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epp_8_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epm_8_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emp_8_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emm_8_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epp_8_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epm_8_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emp_8_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emm_8_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epp_8_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epm_8_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emp_8_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emm_8_3)
            end if
            if (hg_ss .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epp_8_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_cp,laplace_s,
     &                     mu,Epm_8_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emp_8_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_cm,laplace_s,
     &                     mu,Emm_8_4)
            end if
            if (hg_cs .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_c,tau_s,laplace_cp,laplace_s,mu,Epp_8_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_c,tau_s,laplace_cp,laplace_s,mu,Epm_8_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_c,tau_s,laplace_cm,laplace_s,mu,Emp_8_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_c,tau_s,laplace_cm,laplace_s,mu,Emm_8_5)
            end if
            if (ht_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cp,tau_s,laplace_cp,laplace_s,mu,Epp_8_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cm,tau_s,laplace_cp,laplace_s,mu,Epm_8_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cp,tau_s,laplace_cm,laplace_s,mu,Emp_8_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cm,tau_s,laplace_cm,laplace_s,mu,Emm_8_6)
            end if
            if (ht_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sp,laplace_cp,laplace_s,mu,Epp_8_7)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sm,laplace_cp,laplace_s,mu,Epm_8_7)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sp,laplace_cm,laplace_s,mu,Emp_8_7)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sm,laplace_cm,laplace_s,mu,Emm_8_7)
            end if
         end if
         if (hl_s .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epp_9_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epm_9_1)
               call ESRC_SPINMGGA(rho_cp,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emp_9_1)
               call ESRC_SPINMGGA(rho_cm,rho_s,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emm_9_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epp_9_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epm_9_2)
               call ESRC_SPINMGGA(rho_c,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emp_9_2)
               call ESRC_SPINMGGA(rho_c,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emm_9_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epp_9_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epm_9_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emp_9_3)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emm_9_3)
            end if
            if (hg_ss .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epp_9_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sp,
     &                     mu,Epm_9_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emp_9_4)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,tau_c,tau_s,laplace_c,laplace_sm,
     &                     mu,Emm_9_4)
            end if
            if (hg_cs .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_c,tau_s,laplace_c,laplace_sp,mu,Epp_9_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_c,tau_s,laplace_c,laplace_sp,mu,Epm_9_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,tau_c,tau_s,laplace_c,laplace_sm,mu,Emp_9_5)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,tau_c,tau_s,laplace_c,laplace_sm,mu,Emm_9_5)
            end if
            if (ht_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cp,tau_s,laplace_c,laplace_sp,mu,Epp_9_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cm,tau_s,laplace_c,laplace_sp,mu,Epm_9_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cp,tau_s,laplace_c,laplace_sm,mu,Emp_9_6)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_cm,tau_s,laplace_c,laplace_sm,mu,Emm_9_6)
            end if
            if (ht_s .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sp,laplace_c,laplace_sp,mu,Epp_9_7)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sm,laplace_c,laplace_sp,mu,Epm_9_7)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sp,laplace_c,laplace_sm,mu,Emp_9_7)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_sm,laplace_c,laplace_sm,mu,Emm_9_7)
            end if
            if (hl_c .gt. 0.0d0) then
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_s,laplace_cp,laplace_sp,mu,Epp_9_8)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_s,laplace_cm,laplace_sp,mu,Epm_9_8)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_s,laplace_cp,laplace_sm,mu,Emp_9_8)
               call ESRC_SPINMGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         ,tau_c,tau_s,laplace_cm,laplace_sm,mu,Emm_9_8)
            end if
         end if
      end if
      ! Calculate first derivatives
      if (norder .ge. 1) then
         ! d1E / drho_c
         if (hr_c .gt. 0.0d0) then
            d1E(1)=(Ep_1-Em_1)/(2.0d0*hr_c)
         end if
         ! d1E / drho_s
         if (hr_s .gt. 0.0d0) then
            d1E(2)=(Ep_2-Em_2)/(2.0d0*hr_s)
         end if
         ! d1E / dgamma_cc
         if (hg_cc .gt. 0.0d0) then
            d1E(3)=(Ep_3-Em_3)/(2.0d0*hg_cc)
         end if
         ! d1E / dgamma_ss
         if (hg_ss .gt. 0.0d0) then
            d1E(4)=(Ep_4-Em_4)/(2.0d0*hg_ss)
         end if
         ! d1E / dgamma_cs
         if (hg_cs .gt. 0.0d0) then
            d1E(5)=(Ep_5-Em_5)/(2.0d0*hg_cs)
         end if
         ! d1E / dtau_c
         if (ht_c .gt. 0.0d0) then
            d1E(6)=(Ep_6-Em_6)/(2.0d0*ht_c)
         end if
         ! d1E / dtau_s
         if (ht_s .gt. 0.0d0) then
            d1E(7)=(Ep_7-Em_7)/(2.0d0*ht_s)
         end if
         ! d1E / dlaplace_c
         if (hl_c .gt. 0.0d0) then
            d1E(8)=(Ep_8-Em_8)/(2.0d0*hl_c)
         end if
         ! d1E / dlaplace_s
         if (hl_s .gt. 0.0d0) then
            d1E(9)=(Ep_9-Em_9)/(2.0d0*hl_s)
         end if
      end if
      ! Calculate second derivatives
      if (norder .ge. 2) then
         if (hr_c .gt. 0.0d0) then
            ! d2E / (drho_c**2)
            d2E(1)=(Ep_1-2.0d0*E+Em_1)/(hr_c**2)
         end if
         if (hr_s .gt. 0.0d0) then
            ! d2E / (drho_s * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(2)=(Epp_2_1-Epm_2_1-Emp_2_1+Emm_2_1)/(4.0d0*hr_s*
     &                     hr_c)
            end if
            ! d2E / (drho_s**2)
            d2E(3)=(Ep_2-2.0d0*E+Em_2)/(hr_s**2)
         end if
         if (hg_cc .gt. 0.0d0) then
            ! d2E / (dgamma_cc * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(4)=(Epp_3_1-Epm_3_1-Emp_3_1+Emm_3_1)/(4.0d0*hg_cc*
     &                     hr_c)
            end if
            ! d2E / (dgamma_cc * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(5)=(Epp_3_2-Epm_3_2-Emp_3_2+Emm_3_2)/(4.0d0*hg_cc*
     &                     hr_s)
            end if
            ! d2E / (dgamma_cc**2)
            d2E(6)=(Ep_3-2.0d0*E+Em_3)/(hg_cc**2)
         end if
         if (hg_ss .gt. 0.0d0) then
            ! d2E / (dgamma_ss * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(7)=(Epp_4_1-Epm_4_1-Emp_4_1+Emm_4_1)/(4.0d0*hg_ss*
     &                     hr_c)
            end if
            ! d2E / (dgamma_ss * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(8)=(Epp_4_2-Epm_4_2-Emp_4_2+Emm_4_2)/(4.0d0*hg_ss*
     &                     hr_s)
            end if
            ! d2E / (dgamma_ss * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(9)=(Epp_4_3-Epm_4_3-Emp_4_3+Emm_4_3)/(4.0d0*hg_ss*
     &                     hg_cc)
            end if
            ! d2E / (dgamma_ss**2)
            d2E(10)=(Ep_4-2.0d0*E+Em_4)/(hg_ss**2)
         end if
         if (hg_cs .gt. 0.0d0) then
            ! d2E / (dgamma_cs * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(11)=(Epp_5_1-Epm_5_1-Emp_5_1+Emm_5_1)/(4.0d0*hg_cs*
     &                     hr_c)
            end if
            ! d2E / (dgamma_cs * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(12)=(Epp_5_2-Epm_5_2-Emp_5_2+Emm_5_2)/(4.0d0*hg_cs*
     &                     hr_s)
            end if
            ! d2E / (dgamma_cs * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(13)=(Epp_5_3-Epm_5_3-Emp_5_3+Emm_5_3)/(4.0d0*hg_cs*
     &                     hg_cc)
            end if
            ! d2E / (dgamma_cs * dgamma_ss)
            if (hg_ss .gt. 0.0d0) then
               d2E(14)=(Epp_5_4-Epm_5_4-Emp_5_4+Emm_5_4)/(4.0d0*hg_cs*
     &                     hg_ss)
            end if
            ! d2E / (dgamma_cs**2)
            d2E(15)=(Ep_5-2.0d0*E+Em_5)/(hg_cs**2)
         end if
         if (ht_c .gt. 0.0d0) then
            ! d2E / (dtau_c * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(16)=(Epp_6_1-Epm_6_1-Emp_6_1+Emm_6_1)/(4.0d0*ht_c*
     &                     hr_c)
            end if
            ! d2E / (dtau_c * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(17)=(Epp_6_2-Epm_6_2-Emp_6_2+Emm_6_2)/(4.0d0*ht_c*
     &                     hr_s)
            end if
            ! d2E / (dtau_c * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(18)=(Epp_6_3-Epm_6_3-Emp_6_3+Emm_6_3)/(4.0d0*ht_c*
     &                     hg_cc)
            end if
            ! d2E / (dtau_c * dgamma_ss)
            if (hg_ss .gt. 0.0d0) then
               d2E(19)=(Epp_6_4-Epm_6_4-Emp_6_4+Emm_6_4)/(4.0d0*ht_c*
     &                     hg_ss)
            end if
            ! d2E / (dtau_c * dgamma_cs)
            if (hg_cs .gt. 0.0d0) then
               d2E(20)=(Epp_6_5-Epm_6_5-Emp_6_5+Emm_6_5)/(4.0d0*ht_c*
     &                     hg_cs)
            end if
            ! d2E / (dtau_c**2)
            d2E(21)=(Ep_6-2.0d0*E+Em_6)/(ht_c**2)
         end if
         if (ht_s .gt. 0.0d0) then
            ! d2E / (dtau_s * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(22)=(Epp_7_1-Epm_7_1-Emp_7_1+Emm_7_1)/(4.0d0*ht_s*
     &                     hr_c)
            end if
            ! d2E / (dtau_s * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(23)=(Epp_7_2-Epm_7_2-Emp_7_2+Emm_7_2)/(4.0d0*ht_s*
     &                     hr_s)
            end if
            ! d2E / (dtau_s * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(24)=(Epp_7_3-Epm_7_3-Emp_7_3+Emm_7_3)/(4.0d0*ht_s*
     &                     hg_cc)
            end if
            ! d2E / (dtau_s * dgamma_ss)
            if (hg_ss .gt. 0.0d0) then
               d2E(25)=(Epp_7_4-Epm_7_4-Emp_7_4+Emm_7_4)/(4.0d0*ht_s*
     &                     hg_ss)
            end if
            ! d2E / (dtau_s * dgamma_cs)
            if (hg_cs .gt. 0.0d0) then
               d2E(26)=(Epp_7_5-Epm_7_5-Emp_7_5+Emm_7_5)/(4.0d0*ht_s*
     &                     hg_cs)
            end if
            ! d2E / (dtau_s * dtau_c)
            if (ht_c .gt. 0.0d0) then
               d2E(27)=(Epp_7_6-Epm_7_6-Emp_7_6+Emm_7_6)/(4.0d0*ht_s*
     &                     ht_c)
            end if
            ! d2E / (dtau_s**2)
            d2E(28)=(Ep_7-2.0d0*E+Em_7)/(ht_s**2)
         end if
         if (hl_c .gt. 0.0d0) then
            ! d2E / (dlaplace_c * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(29)=(Epp_8_1-Epm_8_1-Emp_8_1+Emm_8_1)/(4.0d0*hl_c*
     &                     hr_c)
            end if
            ! d2E / (dlaplace_c * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(30)=(Epp_8_2-Epm_8_2-Emp_8_2+Emm_8_2)/(4.0d0*hl_c*
     &                     hr_s)
            end if
            ! d2E / (dlaplace_c * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(31)=(Epp_8_3-Epm_8_3-Emp_8_3+Emm_8_3)/(4.0d0*hl_c*
     &                     hg_cc)
            end if
            ! d2E / (dlaplace_c * dgamma_ss)
            if (hg_ss .gt. 0.0d0) then
               d2E(32)=(Epp_8_4-Epm_8_4-Emp_8_4+Emm_8_4)/(4.0d0*hl_c*
     &                     hg_ss)
            end if
            ! d2E / (dlaplace_c * dgamma_cs)
            if (hg_cs .gt. 0.0d0) then
               d2E(33)=(Epp_8_5-Epm_8_5-Emp_8_5+Emm_8_5)/(4.0d0*hl_c*
     &                     hg_cs)
            end if
            ! d2E / (dlaplace_c * dtau_c)
            if (ht_c .gt. 0.0d0) then
               d2E(34)=(Epp_8_6-Epm_8_6-Emp_8_6+Emm_8_6)/(4.0d0*hl_c*
     &                     ht_c)
            end if
            ! d2E / (dlaplace_c * dtau_s)
            if (ht_s .gt. 0.0d0) then
               d2E(35)=(Epp_8_7-Epm_8_7-Emp_8_7+Emm_8_7)/(4.0d0*hl_c*
     &                     ht_s)
            end if
            ! d2E / (dlaplace_c**2)
            d2E(36)=(Ep_8-2.0d0*E+Em_8)/(hl_c**2)
         end if
         if (hl_s .gt. 0.0d0) then
            ! d2E / (dlaplace_s * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(37)=(Epp_9_1-Epm_9_1-Emp_9_1+Emm_9_1)/(4.0d0*hl_s*
     &                     hr_c)
            end if
            ! d2E / (dlaplace_s * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(38)=(Epp_9_2-Epm_9_2-Emp_9_2+Emm_9_2)/(4.0d0*hl_s*
     &                     hr_s)
            end if
            ! d2E / (dlaplace_s * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(39)=(Epp_9_3-Epm_9_3-Emp_9_3+Emm_9_3)/(4.0d0*hl_s*
     &                     hg_cc)
            end if
            ! d2E / (dlaplace_s * dgamma_ss)
            if (hg_ss .gt. 0.0d0) then
               d2E(40)=(Epp_9_4-Epm_9_4-Emp_9_4+Emm_9_4)/(4.0d0*hl_s*
     &                     hg_ss)
            end if
            ! d2E / (dlaplace_s * dgamma_cs)
            if (hg_cs .gt. 0.0d0) then
               d2E(41)=(Epp_9_5-Epm_9_5-Emp_9_5+Emm_9_5)/(4.0d0*hl_s*
     &                     hg_cs)
            end if
            ! d2E / (dlaplace_s * dtau_c)
            if (ht_c .gt. 0.0d0) then
               d2E(42)=(Epp_9_6-Epm_9_6-Emp_9_6+Emm_9_6)/(4.0d0*hl_s*
     &                     ht_c)
            end if
            ! d2E / (dlaplace_s * dtau_s)
            if (ht_s .gt. 0.0d0) then
               d2E(43)=(Epp_9_7-Epm_9_7-Emp_9_7+Emm_9_7)/(4.0d0*hl_s*
     &                     ht_s)
            end if
            ! d2E / (dlaplace_s * dlaplace_c)
            if (hl_c .gt. 0.0d0) then
               d2E(44)=(Epp_9_8-Epm_9_8-Emp_9_8+Emm_9_8)/(4.0d0*hl_s*
     &                     hl_c)
            end if
            ! d2E / (dlaplace_s**2)
            d2E(45)=(Ep_9-2.0d0*E+Em_9)/(hl_s**2)
         end if
      end if
      
      end subroutine


C*****************************************************************************
      subroutine DESRX_BOTHMGGA(ESRX_BOTHMGGA, rho, grad2, tau, laplace,
     & mu, norder, E, d1E, d2E)
C*****************************************************************************
C   Finite difference code for spin-polarized and spin-unpolarized meta-GGA exchange functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(3) : alpha density
C            rho(4) : beta density
C            grad2(4) : grad_alpha . grad_alpha
C            grad2(5) : grad_beta . grad_beta
C            tau(3) : alpha kinetic energy density
C            tau(4) : beta kinetic energy density
C            laplace(3) : laplacian of alpha density
C            laplace(4) : laplacian of beta density
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRX_BOTHGGA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRX_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: September 06, 2018
C*****************************************************************************
      implicit none
      external ESRX_BOTHMGGA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(4), grad2(6), tau(4), laplace(4)
      real*8, intent(out) :: E, d1E(9), d2E(45)
      
      logical :: alpha_not_equal_beta
      real*8 :: rho_a, rho_ap, rho_am, rho_b, rho_bp, rho_bm, gamma_aa,
     & gamma_aap, gamma_aam, gamma_bb, gamma_bbp, gamma_bbm, tau_a,
     & tau_ap, tau_am, tau_b, tau_bp, tau_bm, laplace_a, laplace_ap,
     & laplace_am, laplace_b, laplace_bp, laplace_bm, Ep_1, Em_1,
     & Epp_2_1, Emm_2_1, Epm_2_1, Emp_2_1, Ep_2, Em_2, Epp_3_1, Emm_3_1,
     & Epm_3_1, Emp_3_1, Epp_3_2, Emm_3_2, Epm_3_2, Emp_3_2, Ep_3,
     & Em_3, Epp_4_1, Emm_4_1, Epm_4_1, Emp_4_1, Epp_4_2, Emm_4_2,
     & Epm_4_2, Emp_4_2, Epp_4_3, Emm_4_3, Epm_4_3, Emp_4_3, Ep_4,
     & Em_4, Epp_6_1, Emm_6_1, Epm_6_1, Emp_6_1, Epp_6_2, Emm_6_2,
     & Epm_6_2, Emp_6_2, Epp_6_3, Emm_6_3, Epm_6_3, Emp_6_3, Epp_6_4,
     & Emm_6_4, Epm_6_4, Emp_6_4, Ep_6, Em_6, Epp_7_1, Emm_7_1, Epm_7_1,
     & Emp_7_1, Epp_7_2, Emm_7_2, Epm_7_2, Emp_7_2, Epp_7_3, Emm_7_3,
     & Epm_7_3, Emp_7_3, Epp_7_4, Emm_7_4, Epm_7_4, Emp_7_4, Epp_7_6,
     & Emm_7_6, Epm_7_6, Emp_7_6, Ep_7, Em_7, Epp_8_1, Emm_8_1, Epm_8_1,
     & Emp_8_1, Epp_8_2, Emm_8_2, Epm_8_2, Emp_8_2, Epp_8_3, Emm_8_3,
     & Epm_8_3, Emp_8_3, Epp_8_4, Emm_8_4, Epm_8_4, Emp_8_4, Epp_8_6,
     & Emm_8_6, Epm_8_6, Emp_8_6, Epp_8_7, Emm_8_7, Epm_8_7, Emp_8_7,
     & Ep_8, Em_8, Epp_9_1, Emm_9_1, Epm_9_1, Emp_9_1, Epp_9_2, Emm_9_2,
     & Epm_9_2, Emp_9_2, Epp_9_3, Emm_9_3, Epm_9_3, Emp_9_3, Epp_9_4,
     & Emm_9_4, Epm_9_4, Emp_9_4, Epp_9_6, Emm_9_6, Epm_9_6, Emp_9_6,
     & Epp_9_7, Emm_9_7, Epm_9_7, Emp_9_7, Epp_9_8, Emm_9_8, Epm_9_8,
     & Emp_9_8, Ep_9, Em_9, hr_a, hr_b, hg_aa, hg_bb, hg_ab, ht_a,
     & ht_b, hl_a, hl_b, Ea, Eb, d1Eab(9), d2Eab(45)
      
      rho_a = rho(3)
      rho_b = rho(4)
      gamma_aa = grad2(4)
      gamma_bb = grad2(5)
      tau_a = tau(3)
      tau_b = tau(4)
      laplace_a = laplace(3)
      laplace_b = laplace(4)
      
      ! Check if alpha is different from beta
      ! If alpha equal beta, then alot of derivatives are the same
      if (dabs(laplace_a - laplace_b) .le. 1.d-14 .and. dabs(tau_a - tau
     &   _b) .le. 1.d-14 .and. dabs(gamma_aa - gamma_bb) .le. 1.d-14 .an
     &   d. dabs(rho_a - rho_b) .le. 1.d-14) then
         alpha_not_equal_beta = .false.
      else
         alpha_not_equal_beta = .true.
      end if
      
      E = 0.0d0
      Ea = 0.0d0
      if (norder .ge. 1) then
         d1Eab(:) = 0.0d0
         d1E(:) = 0.0d0
      end if
      if (norder .ge. 2) then
         d2Eab(:) = 0.0d0
         d2E(:) = 0.0d0
      end if
      
      if (rho_a .lt. 1.d-10) then
         hr_a = 0.0d0
      else
         hr_a = 1.d-4*rho_a
      end if
      if (gamma_aa .lt. 1.d-10) then
         hg_aa = 0.0d0
      else
         hg_aa = 1.d-4*gamma_aa
      end if
      if (tau_a .lt. 1.d-10) then
         ht_a = 0.0d0
      else
         ht_a = 1.d-4*tau_a
      end if
      if (laplace_a .lt. 1.d-10) then
         hl_a = 0.0d0
      else
         hl_a = 1.d-4*laplace_a
      end if
      ! Variation of variables
      rho_ap = rho_a + hr_a
      rho_am = rho_a - hr_a
      gamma_aap = gamma_aa + hg_aa
      gamma_aam = gamma_aa - hg_aa
      tau_ap = tau_a + ht_a
      tau_am = tau_a - ht_a
      laplace_ap = laplace_a + hl_a
      laplace_am = laplace_a - hl_a
      
      if (alpha_not_equal_beta) then
         Eb = 0.0d0
         if (rho_b .lt. 1.d-10) then
            hr_b = 0.0d0
         else
            hr_b = 1.d-4*rho_b
         end if
         if (gamma_bb .lt. 1.d-10) then
            hg_bb = 0.0d0
         else
            hg_bb = 1.d-4*gamma_bb
         end if
         if (tau_b .lt. 1.d-10) then
            ht_b = 0.0d0
         else
            ht_b = 1.d-4*tau_b
         end if
         if (laplace_b .lt. 1.d-10) then
            hl_b = 0.0d0
         else
            hl_b = 1.d-4*laplace_b
         end if
         rho_bp = rho_b + hr_b
         rho_bm = rho_b - hr_b
         gamma_bbp = gamma_bb + hg_bb
         gamma_bbm = gamma_bb - hg_bb
         tau_bp = tau_b + ht_b
         tau_bm = tau_b - ht_b
         laplace_bp = laplace_b + hl_b
         laplace_bm = laplace_b - hl_b
      end if
      
      call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_a,laplace_a,mu,Ea)
      if (alpha_not_equal_beta) then
         call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_b,laplace_b,mu,Eb)
         E = Ea + Eb
      else
         E = 2.0d0*Ea
      end if
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_a .gt. 0.0d0) then
            call ESRX_BOTHMGGA(rho_ap,gamma_aa,tau_a,laplace_a,mu,
     &              Ep_1)
            call ESRX_BOTHMGGA(rho_am,gamma_aa,tau_a,laplace_a,mu,
     &              Em_1)
         end if
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHMGGA(rho_bp,gamma_bb,tau_b,laplace_b,mu,
     &              Ep_2)
            call ESRX_BOTHMGGA(rho_bm,gamma_bb,tau_b,laplace_b,mu,
     &              Em_2)
         end if
         if (hg_aa .gt. 0.0d0) then
            call ESRX_BOTHMGGA(rho_a,gamma_aap,tau_a,laplace_a,mu,
     &              Ep_3)
            call ESRX_BOTHMGGA(rho_a,gamma_aam,tau_a,laplace_a,mu,
     &              Em_3)
         end if
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHMGGA(rho_b,gamma_bbp,tau_b,laplace_b,mu,
     &              Ep_4)
            call ESRX_BOTHMGGA(rho_b,gamma_bbm,tau_b,laplace_b,mu,
     &              Em_4)
         end if
         if (ht_a .gt. 0.0d0) then
            call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_ap,laplace_a,mu,
     &              Ep_6)
            call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_am,laplace_a,mu,
     &              Em_6)
         end if
         if (ht_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_bp,laplace_b,mu,
     &              Ep_7)
            call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_bm,laplace_b,mu,
     &              Em_7)
         end if
         if (hl_a .gt. 0.0d0) then
            call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_a,laplace_ap,mu,
     &              Ep_8)
            call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_a,laplace_am,mu,
     &              Em_8)
         end if
         if (hl_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_b,laplace_bp,mu,
     &              Ep_9)
            call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_b,laplace_bm,mu,
     &              Em_9)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
         end if
         if (hg_aa .gt. 0.0d0) then
            if (hr_a .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_ap,gamma_aap,tau_a,laplace_a,
     &                     mu,Epp_3_1)
               call ESRX_BOTHMGGA(rho_am,gamma_aap,tau_a,laplace_a,
     &                     mu,Epm_3_1)
               call ESRX_BOTHMGGA(rho_ap,gamma_aam,tau_a,laplace_a,
     &                     mu,Emp_3_1)
               call ESRX_BOTHMGGA(rho_am,gamma_aam,tau_a,laplace_a,
     &                     mu,Emm_3_1)
            end if
         end if
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            if (hr_b .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_bp,gamma_bbp,tau_b,laplace_b,
     &                     mu,Epp_4_2)
               call ESRX_BOTHMGGA(rho_bm,gamma_bbp,tau_b,laplace_b,
     &                     mu,Epm_4_2)
               call ESRX_BOTHMGGA(rho_bp,gamma_bbm,tau_b,laplace_b,
     &                     mu,Emp_4_2)
               call ESRX_BOTHMGGA(rho_bm,gamma_bbm,tau_b,laplace_b,
     &                     mu,Emm_4_2)
            end if
         end if
         if (ht_a .gt. 0.0d0) then
            if (hr_a .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_ap,gamma_aa,tau_ap,laplace_a,
     &                     mu,Epp_6_1)
               call ESRX_BOTHMGGA(rho_am,gamma_aa,tau_ap,laplace_a,
     &                     mu,Epm_6_1)
               call ESRX_BOTHMGGA(rho_ap,gamma_aa,tau_am,laplace_a,
     &                     mu,Emp_6_1)
               call ESRX_BOTHMGGA(rho_am,gamma_aa,tau_am,laplace_a,
     &                     mu,Emm_6_1)
            end if
            if (hg_aa .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_a,gamma_aap,tau_ap,laplace_a,
     &                     mu,Epp_6_3)
               call ESRX_BOTHMGGA(rho_a,gamma_aam,tau_ap,laplace_a,
     &                     mu,Epm_6_3)
               call ESRX_BOTHMGGA(rho_a,gamma_aap,tau_am,laplace_a,
     &                     mu,Emp_6_3)
               call ESRX_BOTHMGGA(rho_a,gamma_aam,tau_am,laplace_a,
     &                     mu,Emm_6_3)
            end if
         end if
         if (ht_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            if (hr_b .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_bp,gamma_bb,tau_bp,laplace_b,
     &                     mu,Epp_7_2)
               call ESRX_BOTHMGGA(rho_bm,gamma_bb,tau_bp,laplace_b,
     &                     mu,Epm_7_2)
               call ESRX_BOTHMGGA(rho_bp,gamma_bb,tau_bm,laplace_b,
     &                     mu,Emp_7_2)
               call ESRX_BOTHMGGA(rho_bm,gamma_bb,tau_bm,laplace_b,
     &                     mu,Emm_7_2)
            end if
            if (hg_bb .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_b,gamma_bbp,tau_bp,laplace_b,
     &                     mu,Epp_7_4)
               call ESRX_BOTHMGGA(rho_b,gamma_bbm,tau_bp,laplace_b,
     &                     mu,Epm_7_4)
               call ESRX_BOTHMGGA(rho_b,gamma_bbp,tau_bm,laplace_b,
     &                     mu,Emp_7_4)
               call ESRX_BOTHMGGA(rho_b,gamma_bbm,tau_bm,laplace_b,
     &                     mu,Emm_7_4)
            end if
         end if
         if (hl_a .gt. 0.0d0) then
            if (hr_a .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_ap,gamma_aa,tau_a,laplace_ap,
     &                     mu,Epp_8_1)
               call ESRX_BOTHMGGA(rho_am,gamma_aa,tau_a,laplace_ap,
     &                     mu,Epm_8_1)
               call ESRX_BOTHMGGA(rho_ap,gamma_aa,tau_a,laplace_am,
     &                     mu,Emp_8_1)
               call ESRX_BOTHMGGA(rho_am,gamma_aa,tau_a,laplace_am,
     &                     mu,Emm_8_1)
            end if
            if (hg_aa .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_a,gamma_aap,tau_a,laplace_ap,
     &                     mu,Epp_8_3)
               call ESRX_BOTHMGGA(rho_a,gamma_aam,tau_a,laplace_ap,
     &                     mu,Epm_8_3)
               call ESRX_BOTHMGGA(rho_a,gamma_aap,tau_a,laplace_am,
     &                     mu,Emp_8_3)
               call ESRX_BOTHMGGA(rho_a,gamma_aam,tau_a,laplace_am,
     &                     mu,Emm_8_3)
            end if
            if (ht_a .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_ap,laplace_ap,
     &                     mu,Epp_8_6)
               call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_am,laplace_ap,
     &                     mu,Epm_8_6)
               call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_ap,laplace_am,
     &                     mu,Emp_8_6)
               call ESRX_BOTHMGGA(rho_a,gamma_aa,tau_am,laplace_am,
     &                     mu,Emm_8_6)
            end if
         end if
         if (hl_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            if (hr_b .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_bp,gamma_bb,tau_b,laplace_bp,
     &                     mu,Epp_9_2)
               call ESRX_BOTHMGGA(rho_bm,gamma_bb,tau_b,laplace_bp,
     &                     mu,Epm_9_2)
               call ESRX_BOTHMGGA(rho_bp,gamma_bb,tau_b,laplace_bm,
     &                     mu,Emp_9_2)
               call ESRX_BOTHMGGA(rho_bm,gamma_bb,tau_b,laplace_bm,
     &                     mu,Emm_9_2)
            end if
            if (hg_bb .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_b,gamma_bbp,tau_b,laplace_bp,
     &                     mu,Epp_9_4)
               call ESRX_BOTHMGGA(rho_b,gamma_bbm,tau_b,laplace_bp,
     &                     mu,Epm_9_4)
               call ESRX_BOTHMGGA(rho_b,gamma_bbp,tau_b,laplace_bm,
     &                     mu,Emp_9_4)
               call ESRX_BOTHMGGA(rho_b,gamma_bbm,tau_b,laplace_bm,
     &                     mu,Emm_9_4)
            end if
            if (ht_b .gt. 0.0d0) then
               call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_bp,laplace_bp,
     &                     mu,Epp_9_7)
               call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_bm,laplace_bp,
     &                     mu,Epm_9_7)
               call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_bp,laplace_bm,
     &                     mu,Emp_9_7)
               call ESRX_BOTHMGGA(rho_b,gamma_bb,tau_bm,laplace_bm,
     &                     mu,Emm_9_7)
            end if
         end if
      end if
      ! Calculate first derivatives
      if (norder .ge. 1) then
         ! d1E / drho_a
         if (hr_a .gt. 0.0d0) then
            d1Eab(1)=(Ep_1-Em_1)/(2.0d0*hr_a)
         end if
         ! d1E / drho_b
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            d1Eab(2)=(Ep_2-Em_2)/(2.0d0*hr_b)
         end if
         ! d1E / dgamma_aa
         if (hg_aa .gt. 0.0d0) then
            d1Eab(3)=(Ep_3-Em_3)/(2.0d0*hg_aa)
         end if
         ! d1E / dgamma_bb
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            d1Eab(4)=(Ep_4-Em_4)/(2.0d0*hg_bb)
         end if
         ! d1E / dtau_a
         if (ht_a .gt. 0.0d0) then
            d1Eab(6)=(Ep_6-Em_6)/(2.0d0*ht_a)
         end if
         ! d1E / dtau_b
         if (ht_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            d1Eab(7)=(Ep_7-Em_7)/(2.0d0*ht_b)
         end if
         ! d1E / dlaplace_a
         if (hl_a .gt. 0.0d0) then
            d1Eab(8)=(Ep_8-Em_8)/(2.0d0*hl_a)
         end if
         ! d1E / dlaplace_b
         if (hl_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            d1Eab(9)=(Ep_9-Em_9)/(2.0d0*hl_b)
         end if
      end if
      ! Calculate second derivatives
      if (norder .ge. 2) then
         if (hr_a .gt. 0.0d0) then
            ! d2E / (drho_a**2)
            d2Eab(1)=(Ep_1-2.0d0*Ea+Em_1)/(hr_a**2)
         end if
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            ! d2E / (drho_b**2)
            d2Eab(3)=(Ep_2-2.0d0*Eb+Em_2)/(hr_b**2)
         end if
         if (hg_aa .gt. 0.0d0) then
            ! d2E / (dgamma_aa * drho_a)
            if (hr_a .gt. 0.0d0) then
               d2Eab(4)=(Epp_3_1-Epm_3_1-Emp_3_1+Emm_3_1)/(4.0d0*
     &                     hg_aa*hr_a)
            end if
            ! d2E / (dgamma_aa**2)
            d2Eab(6)=(Ep_3-2.0d0*Ea+Em_3)/(hg_aa**2)
         end if
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            ! d2E / (dgamma_bb * drho_b)
            if (hr_b .gt. 0.0d0) then
               d2Eab(8)=(Epp_4_2-Epm_4_2-Emp_4_2+Emm_4_2)/(4.0d0*
     &                     hg_bb*hr_b)
            end if
            ! d2E / (dgamma_bb**2)
            d2Eab(10)=(Ep_4-2.0d0*Eb+Em_4)/(hg_bb**2)
         end if
         if (ht_a .gt. 0.0d0) then
            ! d2E / (dtau_a * drho_a)
            if (hr_a .gt. 0.0d0) then
               d2Eab(16)=(Epp_6_1-Epm_6_1-Emp_6_1+Emm_6_1)/(4.0d0*
     &                     ht_a*hr_a)
            end if
            ! d2E / (dtau_a * dgamma_aa)
            if (hg_aa .gt. 0.0d0) then
               d2Eab(18)=(Epp_6_3-Epm_6_3-Emp_6_3+Emm_6_3)/(4.0d0*
     &                     ht_a*hg_aa)
            end if
            ! d2E / (dtau_a**2)
            d2Eab(21)=(Ep_6-2.0d0*Ea+Em_6)/(ht_a**2)
         end if
         if (ht_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            ! d2E / (dtau_b * drho_b)
            if (hr_b .gt. 0.0d0) then
               d2Eab(23)=(Epp_7_2-Epm_7_2-Emp_7_2+Emm_7_2)/(4.0d0*
     &                     ht_b*hr_b)
            end if
            ! d2E / (dtau_b * dgamma_bb)
            if (hg_bb .gt. 0.0d0) then
               d2Eab(25)=(Epp_7_4-Epm_7_4-Emp_7_4+Emm_7_4)/(4.0d0*
     &                     ht_b*hg_bb)
            end if
            ! d2E / (dtau_b**2)
            d2Eab(28)=(Ep_7-2.0d0*Eb+Em_7)/(ht_b**2)
         end if
         if (hl_a .gt. 0.0d0) then
            ! d2E / (dlaplace_a * drho_a)
            if (hr_a .gt. 0.0d0) then
               d2Eab(29)=(Epp_8_1-Epm_8_1-Emp_8_1+Emm_8_1)/(4.0d0*
     &                     hl_a*hr_a)
            end if
            ! d2E / (dlaplace_a * dgamma_aa)
            if (hg_aa .gt. 0.0d0) then
               d2Eab(31)=(Epp_8_3-Epm_8_3-Emp_8_3+Emm_8_3)/(4.0d0*
     &                     hl_a*hg_aa)
            end if
            ! d2E / (dlaplace_a * dtau_a)
            if (ht_a .gt. 0.0d0) then
               d2Eab(34)=(Epp_8_6-Epm_8_6-Emp_8_6+Emm_8_6)/(4.0d0*
     &                     hl_a*ht_a)
            end if
            ! d2E / (dlaplace_a**2)
            d2Eab(36)=(Ep_8-2.0d0*Ea+Em_8)/(hl_a**2)
         end if
         if (hl_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            ! d2E / (dlaplace_b * drho_b)
            if (hr_b .gt. 0.0d0) then
               d2Eab(38)=(Epp_9_2-Epm_9_2-Emp_9_2+Emm_9_2)/(4.0d0*
     &                     hl_b*hr_b)
            end if
            ! d2E / (dlaplace_b * dgamma_bb)
            if (hg_bb .gt. 0.0d0) then
               d2Eab(40)=(Epp_9_4-Epm_9_4-Emp_9_4+Emm_9_4)/(4.0d0*
     &                     hl_b*hg_bb)
            end if
            ! d2E / (dlaplace_b * dtau_b)
            if (ht_b .gt. 0.0d0) then
               d2Eab(43)=(Epp_9_7-Epm_9_7-Emp_9_7+Emm_9_7)/(4.0d0*
     &                     hl_b*ht_b)
            end if
            ! d2E / (dlaplace_b**2)
            d2Eab(45)=(Ep_9-2.0d0*Eb+Em_9)/(hl_b**2)
         end if
      end if
      
      ! If alpha is equal to beta, assing derivatives
      ! that was not calculated above
      if (.not. alpha_not_equal_beta) then
         if (norder .ge. 1) then
            d1Eab(2)=d1Eab(1)
            d1Eab(4)=d1Eab(3)
            d1Eab(7)=d1Eab(6)
            d1Eab(9)=d1Eab(8)
         end if
         if (norder .ge. 2) then
            d2Eab(3)=d2Eab(1)
            d2Eab(8)=d2Eab(4)
            d2Eab(10)=d2Eab(6)
            d2Eab(23)=d2Eab(16)
            d2Eab(25)=d2Eab(18)
            d2Eab(28)=d2Eab(21)
            d2Eab(38)=d2Eab(29)
            d2Eab(40)=d2Eab(31)
            d2Eab(43)=d2Eab(34)
            d2Eab(45)=d2Eab(36)
         end if
      end if
      
      ! Transform from derivatess of alpha and beta
      ! to derivatives of total and charge density.
      ! The following transformations are only valid for
      ! exchnage functional since there are no
      ! cross-derivatives between alpha and beta.
      if (norder .ge. 1) then
         d1E(1) = 0.5d0*(d1Eab(1) + d1Eab(2))
         d1E(2) = 0.5d0*(d1Eab(1) - d1Eab(2))
         if (hg_aa .gt. 0.0d0 .or. hg_bb .gt. 0.0d0) then
            d1E(3) = 0.25d0*(d1Eab(3) + d1Eab(4))
            d1E(4) = d1E(3)
            d1E(5) = 0.5d0*(d1Eab(3) - d1Eab(4))
         end if
         if (ht_a .gt. 0.0d0 .or. ht_b .gt. 0.0d0) then
            d1E(6) = 0.5d0*(d1Eab(6) + d1Eab(7))
            d1E(7) = 0.5d0*(d1Eab(6) - d1Eab(7))
         end if
         if (hl_a .gt. 0.0d0 .or. hl_b .gt. 0.0d0) then
            d1E(8) = 0.5d0*(d1Eab(8) + d1Eab(9))
            d1E(9) = 0.5d0*(d1Eab(8) - d1Eab(9))
         end if
      end if
      
      if (norder .ge. 2) then
         d2E(1) = 0.25d0*(d2Eab(1) + d2Eab(3))
         d2E(2) = 0.25d0*(d2Eab(1) - d2Eab(3))
         d2E(3) = d2E(1)
         if (hg_aa .gt. 0.0d0 .or. hg_bb .gt. 0.0d0) then
            d2E(4) = 0.125d0*(d2Eab(4) + d2Eab(8))
            d2E(5) = 0.125d0*(d2Eab(4) - d2Eab(8))
            d2E(6) = 0.0625d0*(d2Eab(6) + d2Eab(10))
            d2E(7) = d2E(4)
            d2E(8) = d2E(5)
            d2E(9) = d2E(6)
            d2E(10) = d2E(6)
            d2E(11) = 2.0d0*d2E(5)
            d2E(12) = 2.0d0*d2E(4)
            d2E(13) = 0.125d0*(d2Eab(6) - d2Eab(10))
            d2E(14) = d2E(13)
            d2E(15) = 4.0d0*d2E(6)
         end if
         if (ht_a .gt. 0.0d0 .or. ht_b .gt. 0.0d0) then
            d2E(16) = 0.25d0*(d2Eab(16) + d2Eab(23))
            d2E(17) = 0.25d0*(d2Eab(16) - d2Eab(23))
            d2E(18) = 0.125d0*(d2Eab(18) + d2Eab(25))
            d2E(19) = d2E(18)
            d2E(20) = 2.0d0*d2E(24)
            d2E(21) = 0.25d0*(d2Eab(21) + d2Eab(28))
            d2E(22) = 0.25d0*(d2Eab(16) - d2Eab(23))
            d2E(23) = d2E(16)
            d2E(24) = 0.125d0*(d2Eab(18) - d2Eab(25))
            d2E(25) = d2E(24)
            d2E(26) = 2.0d0*d2E(18)
            d2E(27) = 0.25d0*(d2Eab(21) - d2Eab(28))
            d2E(28) = d2E(21)
         end if
         if (hl_a .gt. 0.0d0 .or. hl_b .gt. 0.0d0) then
            d2E(29) = 0.25d0*(d2Eab(29) + d2Eab(38))
            d2E(30) = 0.25d0*(d2Eab(29) - d2Eab(38))
            d2E(31) = 0.125d0*(d2Eab(31) + d2Eab(40))
            d2E(32) = d2E(31)
            d2E(33) = 2.0d0*d2E(39)
            d2E(34) = 0.25d0*(d2Eab(34) + d2Eab(43))
            d2E(35) = 0.25d0*(d2Eab(34) - d2Eab(43))
            d2E(36) = 0.25d0*(d2Eab(36) + d2Eab(45))
            d2E(37) = 0.25d0*(d2Eab(29) - d2Eab(38))
            d2E(38) = d2E(29)
            d2E(39) = 0.125d0*(d2Eab(31) - d2Eab(40))
            d2E(40) = d2E(39)
            d2E(41) = 2.0d0*d2E(31)
            d2E(42) = 0.25d0*(d2Eab(34) - d2Eab(43))
            d2E(43) = d2E(34)
            d2E(44) = 0.25d0*(d2Eab(36) - d2Eab(45))
            d2E(45) = d2E(36)
         end if
      end if
      
      end subroutine


C*****************************************************************************
      SUBROUTINE DESR(ESR,rhoc,grdcc,mu,norder,E,d1E,d2E)
C*****************************************************************************
C     Compute first-order and/or second-order numerical derivatives 
C     of a generic closed-shell short-range GGA functional ESR(rhoc,grdcc,mu,E)
C
C     Input: ESR    : functional routine
C            rhoc   : charge density = alpha density + beta density
C            grdcc  : grad(rhoc) . grad(rhoc)
C            mu     : Interaction parameter
C            norder : order of derivatives
C
C     Ouput: E                         : energy
C            d1E(1) =  d1Edrhoc        : first derivative wrt rhoc
C            d1E(3) =  d1Edgrdcc       : first derivative wrt grdcc
C
C            d2E(01) = d2Edrhocdrhoc   : second derivative wrt rhoc and rhoc
C            d2E(04) = d2Edrhocdgrdcc  : second derivative wrt rhoc and grdcc
C            d2E(06) = d2Edgrdccdgrdcc : second derivative wrt grdcc and grdcc
C
C     Author: J. Toulouse
C     Date  : 16-02-05
C     Revision: Introduce d1E and d2E arrays (17-02-16), E. D. Hedegaard
C*****************************************************************************
      implicit none
#include "priunit.h"
      external ESR

      integer, intent(in) :: norder
      real*8, intent(in) :: rhoc, grdcc, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)

      real*8 :: hr, hg
      real*8 :: rhocm, rhocp, rhocpp, grdccm, grdccp

      real*8 :: Em0, Ep0, Epp0
      real*8 :: E0m, E0p, E0pp
      real*8 :: Emm, Epp, Epm, Emp

C Check
      if(norder .lt. 0 .or. norder .gt. 2) then
       print *,'DESR error: norder=',norder
       call quit('DESR error: norder must be 0, 1, or 2')
      endif

C Energy -------------------------------------------------------------
      call ESR(rhoc,grdcc,mu,E)

C First-order derivatives --------------------------------------------
      if (norder .le. 0) go to 9000 ! return
      d1E(:) = 0.0d0

C     (optimal?) numerical step
      hr = 1.d-4*rhoc
      hg = 1.d-4*grdcc

C     variations of variables
      rhocm = rhoc-hr
      rhocp = rhoc+hr
      grdccm = grdcc-hg
      grdccp = grdcc+hg

C     derivative wrt rhoc ---------------------
!      if(rhocm .ge. 0.d0) then 
!      central 2nd-order approximation
       call ESR(rhocm,grdcc,mu,Em0)
       call ESR(rhocp,grdcc,mu,Ep0)
       d1E(1) = (Ep0 - Em0)/(2.d0*hr) ! d1Edrhoc
!      else 
!       stop 'rhocm < 0'
!!      forward 1st-order approximation
!       call ESR(rhocp,grdcc,mu,Ep0)
!       d1Edrhoc = (Ep0 - E)/hr
!      endif

C     derivative wrt grdcc ---------------------
!      if(grdccm .ge. 0.d0) then
!      central 2nd-order approximation
       call ESR(rhoc,grdccm,mu,E0m)
       call ESR(rhoc,grdccp,mu,E0p)
       d1E(3) = (E0p - E0m)/(2.d0*hg) ! d1Edgrdcc
!      else
!       stop 'grdccm < 0'
!!      forward 1st-order approximation
!       call ESR(rhoc,grdccp,mu,E0p)
!       d1Edgrdcc = (E0p - E)/hg
!      endif

C Second-order derivatives --------------------------------------------
      if (norder .le. 1) go to 9000 ! return
      if (norder .ge. 2) d2E(:) = 0.0d0

C     derivative wrt rhoc and rhoc -----------------
!      if(rhocm .ge. 0.d0) then 
!      central 2nd-order approximation
       
       d2E(1) = (Ep0 - 2.d0*E + Em0)/(hr*hr)
!      else 
!       stop 'rhocm < 0'
!!      forward 1st-order approximation
!       rhocpp = rhoc + 2.d0*hr
!       call ESR(rhocp,grdcc,mu,Ep0)
!       call ESR(rhocpp,grdcc,mu,Epp0)
!       d2E(1) = (Epp0 -2.d0*Ep0 + E)/hr2
!      endif

C     derivative wrt grdcc and grdcc -----------------------
!      if(grdccm .ge. 0.d0) then
!      central 2nd-order approximation
       d2E(6) = (E0p - 2.d0*E + E0m)/(hg*hg)
!      else
!       stop 'grdccm < 0'
!!      forward 1st-order approximation
!       grdccpp = grdcc + 2.d0*hg
!       call ESR(rhoc,grdccp,mu,E0p)
!       call ESR(rhoc,grdccpp,mu,E0pp)
!       d2E(2) = (E0pp -2.d0*E0p + E)/hg2
!      endif

C     derivative wrt rhoc and grdcc -----------------------
!      if(rhocm.ge. 0.d0 .and. grdccm .ge. 0.d0) then
       call ESR(rhocp,grdccp,mu,Epp)
       call ESR(rhocm,grdccp,mu,Emp)
       call ESR(rhocp,grdccm,mu,Epm)
       call ESR(rhocm,grdccm,mu,Emm)
!      central 2nd-order approximation
       d2E(4) = (Epp - Emp - Epm + Emm)/(4.d0*hr*hg)
!      else
!       stop 'rhocm < 0 or grdccm < 0'
!       call ESR(rhocp,grdccp,mu,Epp)
!!      forward 1st-order approximation
!       d2E(3) = (Epp - E0p - Ep0 + E)/hrg
!      endif

 9000 return
      end


C*****************************************************************************
      SUBROUTINE DELA(ESR,rhoc,grdcc,VLAMBDA,norder,d1E,d2E)
C*****************************************************************************
C Linear complement exchange PBE functional, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      implicit none

      external ESR

      real*8, intent(in) :: rhoc, grdcc, VLAMBDA
      integer, intent(in) :: norder
      real*8, intent(out) :: d1E(9), d2E(45)
      real*8 :: E

      call DESR(ESR,rhoc,grdcc,0.d0,norder,E,d1E,d2E)

      d1E(1) = (1.d0- VLAMBDA) * d1E(1)
      d1E(3) = (1.d0- VLAMBDA) * d1E(3)

      d2E(1) = (1.d0- VLAMBDA) * d2E(1)
      d2E(4) = (1.d0- VLAMBDA) * d2E(4)
      d2E(6) = (1.d0- VLAMBDA) * d2E(6)

      return
      end


C*****************************************************************************
      SUBROUTINE DELANSC(ESR,rhoc,grdcc,VLAMBDA,norder,d1E,d2E)
C*****************************************************************************
C Linear complement (NON-Scaled) correlation PBE GWS functional, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      implicit none

      external ESR
      real*8, intent(in) :: rhoc, grdcc, VLAMBDA
      integer, intent(in) :: norder
      real*8, intent(out) :: d1E(9), d2E(45)
      real*8 :: E

      call DESR(ESR,rhoc,grdcc,0.d0,norder,E,d1E,d2E)

      d1E(1) = (1.d0- VLAMBDA**2) * d1E(1)
      d1E(3) = (1.d0- VLAMBDA**2) * d1E(3)

      d2E(1) = (1.d0- VLAMBDA**2) * d2E(1)
      d2E(4) = (1.d0- VLAMBDA**2) * d2E(4)
      d2E(6) = (1.d0- VLAMBDA**2) * d2E(6)

      return
      end


C*****************************************************************************
      SUBROUTINE DELASC(ESR,rhoc,grdcc,VLAMBDA,norder,d1E,d2E)
C*****************************************************************************
C Linear complement Scaled correlation PBE functional, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      implicit none

      external ESR
      integer, intent(in) :: norder
      real*8, intent(in) :: rhoc, grdcc, VLAMBDA
      real*8, intent(out) :: d1E(9), d2E(45)
      real*8 :: rhocscaled, grdccscaled
      real*8 :: Ecoul, d1Ecoul(9), d2Ecoul(45)
      real*8 :: Esca,  d1Esca(9),  d2Esca(45)

! TODO HJAAJ 
      call quit('hjaaj: need to define new d1E and d2E addresses')

!     d1Ecoul(1) = d1Edrhoc_coul
!     d1Ecoul(3) = d1Edgrdcc_coul
!     d2Ecoul(1) = d2Edrhocdrhoc_coul
!     d2Ecoul(4) = d2Edrhocdgrdcc_coul
!     d2Ecoul(6) = d2Edgrdccdgrdcc_coul 
      d1Ecoul(:) = 0.0d0
      d2Ecoul(:) = 0.0d0

      call  DESR(ESR,rhoc,grdcc,0.d0,norder,Ecoul,d1Ecoul,d2Ecoul)

      rhocscaled = rhoc/VLAMBDA**3
      grdccscaled  = grdcc/VLAMBDA**8

!     d1Esca(1) = d1Edrhoc_sca
!     d1Esca(3) = d1Edgrdcc_sca 
!     d2Esca(1) = d2Edrhocdrhoc_sca 
!     d2Esca(4) = d2Edrhocdgrdcc_sca
!     d2Esca(6) = d2Edgrdccdgrdcc_sca 
      d1Esca(:) = 0.0d0
      d2Esca(:) = 0.0d0

      call  DESR(ESR,rhocscaled,grdccscaled,0.d0,norder,
     &   Esca,d1Esca,d2Esca)

      d1E(1) = d1E(1) + (d1Ecoul(1) - VLAMBDA**2 * d1Esca(1))
      d1E(3) = d1E(3) + (d1Ecoul(3) - (d1Esca(3)/VLAMBDA**3))

      d2E(1) = d2E(1) + (d2Ecoul(1) - (d2Esca(1)/VLAMBDA))
      d2E(4) = d2E(4) + (d2Ecoul(4) - (d2Esca(4)/VLAMBDA**2))
      d2E(6) = d2E(6) + (d2Ecoul(6) - (d2Esca(6)/VLAMBDA**3))

      return
      end


C*****************************************************************************
      SUBROUTINE DESRSPIN(ESRSPIN,rhoc,rhos,mu,dospin,norder,E,d1E,d2E)
C*****************************************************************************
C     Compute first-order and/or second-order numerical derivatives 
C     of a generic spin-dependent short-range LDA functional ESRSPIN(rhoc,rhos,mu,dospin,E)
C
C     Input: ESRSPIN: functional routine
C            rhoc   : total charge density
C            rhos   : spin density
C            dospin : use spin density
C            mu     : Interaction parameter
C            norder : order of derivatives
C
C     Ouput: d1E(1) = d1Edrhoc      : first derivative wrt rhoc
C            d1E(2) = d1Edrhos      : first derivative wrt rhos
C            d1E(3-5) not calculated, see DESRX_SPINGGA and
C                     DESRC_SPINGGA routines
C
C            d2E(1) = d2Edrhocdrhoc : second derivative wrt rhoc and rhoc
C            d2E(2) = d2Edrhocdrhos : second derivative wrt rhoc and rhos
C            d2E(3) = d2Edrhosdrhos : second derivative wrt rhos and rhos
C            d2E(4-15) not calculated, see DESRX_SPINGGA and
C                   DESRC_SPINGGA routines
C              
C     Created: 17-08-09, J. Toulouse
C     Revision: Introduce d1E and d2E arrays (17-02-16), E. D. Hedegaard
C*****************************************************************************
      implicit none

      external ESRSPIN
      integer, intent(in) :: norder
      logical, intent(in) :: dospin
      real*8, intent(in) :: rhoc, rhos, mu
      real*8, intent(out) :: E, d1E(9), d2E(45)

      real*8 :: hrc, hrs
      real*8 :: rhocm, rhocp, rhosp, rhosm
      real*8 :: arhos

      real*8 :: Em0, Ep0, Epp0
      real*8 :: E0m, E0p, E0pp
      real*8 :: Emm, Epp, Epm, Emp

      logical :: pure_spin

C     check
      if(norder .le. 0 .or. norder .ge. 3) then
       print *,'DESRSPIN: norder=',norder
       call quit('DESRSPIN: norder must be 1 or 2')
      endif

      if (norder .ge. 1) d1E(:) = 0.0d0
      if (norder .ge. 2) d2E(:) = 0.0d0

      arhos = abs(rhos)
      if (arhos .gt. rhoc) then
         call quit('DESRSPIN error: |rhos| > rhoc')
      end if
      pure_spin = (rhoc - arhos .le. 1.d-14 ) ! pure alpha or pure beta

C     energy
      call ESRSPIN(rhoc,rhos,mu,dospin,E)

C     (optimal?) numerical step
      hrc = 1.d-4*rhoc
      hrs = hrc
!      arhos_min = 1.d-10*max(1.d0,rhoc) ! for numerical accuracy
!      if (arhos .gt. 1.d-10*rhoc) then
!      !   not closed shell
!         hrs = 1.d-4*rhos
!      else
!         hrs = 0.0d0
!      end if

C     variations of variables
      rhocm = rhoc-hrc
      rhocp = rhoc+hrc
      rhosm = rhos-hrs
      rhosp = rhos+hrs

C First-order derivatives --------------------------------------------

C     derivative wrt rhoc ---------------------
!      central 2nd-order approximation

       if (pure_spin) then
          call ESRSPIN(rhocm,rhosm,mu,dospin,Em0)
          call ESRSPIN(rhocp,rhosp,mu,dospin,Ep0)
          d1E(1) = (Ep0 - Em0)/(2.d0*hrc)
!         if (rhos .gt. 0.0d0) then
!            d1E(2) = d1E(1)
!         else
!            d1E(2) = -d1E(1)
!         end if
       else
          call ESRSPIN(rhocm,rhos,mu,dospin,Em0)
          call ESRSPIN(rhocp,rhos,mu,dospin,Ep0)
          d1E(1) = (Ep0 - Em0)/(2.d0*hrc)
          d1E(2) = 0.0d0
       end if

C     derivative wrt rhos ---------------------
!      central 2nd-order approximation
       if (dospin .and. .not. pure_spin) then
         if (hrs .gt. 0.0d0) then
           call ESRSPIN(rhoc,rhosm,mu,dospin,E0m)
           call ESRSPIN(rhoc,rhosp,mu,dospin,E0p)
           d1E(2) = (E0p - E0m)/(2.d0*hrs)
         else
           d1E(2) = 0.0d0
         end if
       endif

C Second-order derivatives --------------------------------------------
      if(norder .ge. 2) then

C       derivative wrt rhoc and rhoc -----------------
!       central 2nd-order approximation
        d2E(1) = (Ep0 - 2.d0*E + Em0)/(hrc*hrc)
        if (dospin) then
          if (pure_spin) then
            d2E(3) = 0.0d0
            d2E(2) = 0.0d0
          else if (hrs .gt. 0.0d0) then
C       derivative wrt rhos and rhos -----------------------
!       central 2nd-order approximation
            d2E(3) = (E0p - 2.d0*E + E0m)/(hrs*hrs)
       
C         derivative wrt rhoc and rhos -----------------------
          call ESRSPIN(rhocp,rhosp,mu,dospin,Epp)
          call ESRSPIN(rhocm,rhosp,mu,dospin,Emp)
          call ESRSPIN(rhocp,rhosm,mu,dospin,Epm)
          call ESRSPIN(rhocm,rhosm,mu,dospin,Emm)
!         central 2nd-order approximation
            d2E(2) = (Epp - Emp - Epm + Emm)/(4.d0*hrc*hrs)
          else
            d2E(3) = 0.0d0
            d2E(2) = 0.0d0
          end if
        endif

      endif

      return
      end


C*****************************************************************************
      SUBROUTINE DESRC_SPINGGA(ESRC_SPINGGA,rho,grd2,mu,norder,
     &   E,d1E,d2E)
C*****************************************************************************
C     Compute first-order and/or second-order numerical derivatives 
C     of a generic spin-dependent short-range GGA functional
C     ESRC_SPINGGA(rhoc,rhos,grdcc,grdss,grdcs,mu,E)
C
C NOTE: (April 2018, hjaaj)
C     This generic routine can be called for any LDA_SPIN and GGA_SPIN;
C     derivatives for grdcc, grdss, and grdcs are only calculated
C     for those which are non-zero.
C     For example for PBE_SPIN (which does not depend on grdss and grdcs) the
C     grdss = grd2(2) and grdcs = GRD2(3) should be zero on input to avoid unnecssary
C     derivatives wrt grdss and grdcs, although the true grdss, grdcs are non-zero.
C     Can also be used for LDA by setting GRD2(1:3) elements to zero.
C     Derivatives wrt rhoc and rhos are always calculated.
C
C     This routine performes the derivatives with respect to charge and spin variables,
C     normally best suited for correlation functionals.
C     The routine DESRX_SPINGGA performs the derivatives with respect to alpha and beta
C     variables, always best suited for exchange functionals because these never mix
C     alpha and beta variables, thus all cross derivatives between alpha and beta are zero.
C
C Input:     ESRC_SPINGGA : spin-density GGA functional routine
C            rho(1)     : rhoc; total density
C            rho(2)     : rhos; spin density
C            grd2(1)    : grdcc = grad(rhoc) . grad(rhoc)
C            grd2(2)    : grdss = grad(rhos) . grad(rhos)
C            grd2(3)    : grdcs = grad(rhoc) . grad(rhos)
C            mu         : Interaction parameter
C            norder     : order of derivatives
C
C Output:    E      = correlation energy density
C
C            d1E(1) = d1Edrhoc      : first derivative wrt rhoc
C            d1E(2) = d1Edrhos      : first derivative wrt rhos
C            d1E(3) = d1Edgrdcc     : first derivative wrt grdcc
C            d1E(4) = d1Edgrdss     : first derivative wrt grdss
C            d1E(5) = d1Edgrdcs     : first derivative wrt grdcs
C
C            d2E(1) = d2Edrhocdrhoc : second derivative wrt rhoc and rhoc
C            d2E(2) = d2Edrhocdrhos : second derivative wrt rhoc and rhos
C            d2E(3) = d2Edrhosdrhos : second derivative wrt rhos and rhos
C            d2E(4) = d2Edrhocdgrdcc : second derivative wrt rhoc and grdcc
C            d2E(5) = d2Edrhosdgrdcc : second derivative wrt rhos and grdcc
C            d2E(6) = d2Edgrdccdgrdcc : second derivative wrt grdcc and grdcc
C            d2E(7) = d2Edrhocdgrdss : second derivative wrt rhoc and grdss
C            d2E(8) = d2Edrhosdgrdss : second derivative wrt rhos and grdss
C            d2E(9) = d2Edgrdccdgrdss : second derivative wrt grdcc and grdss
C            d2E(10)= d2Edgrdssdgrdss : second derivative wrt grdss and grdss
C            d2E(11)= d2Edrhocdgrdcs:
C            d2E(12)= d2Edrhosdgrdcs:
C            d2E(13)= d2Edgrdccdgrdcs:
C            d2E(14)= d2Edgrdssdgrdcs:
C            d2E(15)= d2Edgrdcsdgrdcs:
C
C     Created: 29-01-16, H. J. Aa. Jensen + Erik D. Hedegaard)
C     (based on DESRSPIN routine by J. Toulouse)
C*****************************************************************************
      implicit none

      external ESRC_SPINGGA
      integer, intent(in) :: norder
      real*8, intent(in) :: rho(2), grd2(3), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)

      real*8 :: rhoc, rhos, grdcc, grdss, grdcs
      real*8 :: arhos, agrdcs
      real*8 :: hrc, hrs, hgc, hgs, hgcs
      real*8 :: rhocm, rhocp, rhosp, rhosm
      real*8 :: grdccm, grdccp, grdssp, grdssm, grdcsm, grdcsp

      real*8 :: Em000, Ep000, E0m00, E0p00
      real*8 :: E00m0, E00p0, E000m, E000p
      real*8 :: E0000m,E0000p

      real*8 :: Emm00, Epp00, Epm00, Emp00
      real*8 :: Em0m0, Ep0p0, Ep0m0, Em0p0
      real*8 :: Em00m, Ep00p, Ep00m, Em00p
      real*8 :: E0mm0, E0pp0, E0pm0, E0mp0
      real*8 :: E0m0m, E0p0p, E0p0m, E0m0p
      real*8 :: E00mm, E00pp, E00pm, E00mp
      real*8 :: E_test

      logical :: pure_spin

C     check

      if(norder .lt. 0 .or. norder .ge. 3) then
       print *,'DESRC_SPINGGA: norder=',norder
       call quit('DESRC_SPINGGA: norder must be 0, 1, or 2')
      endif

      if (norder .ge. 1) d1E(:) = 0.0d0
      if (norder .ge. 2) d2E(:) = 0.0d0
      
      rhoc  = rho(1)
      rhos  = rho(2)
      grdcc = grd2(1)
      grdss = grd2(2)
      grdcs = grd2(3)

      arhos = abs(rhos)
      agrdcs = abs(grdcs)

      if (arhos .gt. rhoc) then
         call quit('DESRC_SPINGGA error: |rhos| > rhoc')
      end if
      if (grdcc .lt. 0.0d0 .or. grdss .lt. 0.0d0) then
         call quit('DESRC_SPINGGA error: grdcc or grdss negative')
      end if

      if (rhoc .lt. 1.d-10) then ! no significant density in this point
         E = 0.0d0
         return
      end if


C     energy
      call ESRC_SPINGGA(rhoc,rhos,grdcc,grdss,grdcs,mu,E)

C     (optimal?) numerical step
C     (step length of hx=1.d-4*X should give a relative accuracy of approx. 8 digits
C      because we use central 2nd-order approximation)
C     NB! The zeroing of grdcs, grdss, grdcc depending on type means these will be skipped below.

      hrc = 1.d-4*rhoc
      hrs = hrc

      if (grdcc .lt. 1.d-20) then
         hgc = 0.0d0
      else
         hgc = 1.d-4*grdcc ! grdcc .ge. 0 so hgc .ge. 0
      end if
      if (grdss .lt. 1.d-20) then
         hgs = 0.0d0
      else
         hgs = 1.d-4*grdss
      end if
      if (agrdcs .lt. 1.d-20) then
         hgcs = 0.0d0
      else
         hgcs = 1.d-4*agrdcs
      end if
      
!      arhos_min = 1.d-10*max(1.d0,rhoc) ! for numerical accuracy
!      if (arhos .gt. 1.d-10*rhoc) then
!      !   not closed shell
!         hrs = 1.d-4*rhos
!      else
!         hrs = 0.0d0
!      end if
!      hg = 1.d-4*grd2

C     variations of variables
      rhocm  = rhoc-hrc
      rhocp  = rhoc+hrc
      rhosm  = rhos-hrs
      rhosp  = rhos+hrs
      grdccm = grdcc-hgc
      grdccp = grdcc+hgc
      grdssm = grdss-hgs
      grdssp = grdss+hgs
      grdcsm = grdcs-hgcs
      grdcsp = grdcs+hgcs
C
C     Special treatment of pure spin (because rhosm or rhosp will correspond
C     to negative rhoa or negative rhob)

      pure_spin = ( rhosp .gt. rhoc .or. -rhosm .gt. rhoc ) ! pure alpha or pure beta needs special handling
      if (pure_spin) then

       call quit('DESRC_SPINGGA: pure_spin not implemented yet')

         if (norder .ge. 1) then
           d1E(:) = 0.0d0
           call ESRC_SPINGGA(rhocm,rhosm,grdcc,grdss,grdcs,mu,Em000)
           call ESRC_SPINGGA(rhocp,rhosp,grdcc,grdss,grdcs,mu,Ep000)
           d1E(1) = (Ep000 - Em000)/(2.d0*hrc) ! d1Edrhoc
!         if (rhos .gt. 0.0d0) then
!            d1Edrhos = d1Edrhoc
!         else
!            d1Edrhos = -d1Edrhoc
!         end if
         end if

         if (norder .ge. 2) then
            d2E(:) = 0.0d0
         end if

         return
      end if

! (All derivatives in the following are calculated with central 2nd-order approximation)

C --- First-order derivatives --------------------------------------------
      if (norder .lt. 1) return

      d1E(:) = 0.0d0

C     derivatives wrt rhoc, rhos, grdcc, grdss ---------------------

      call ESRC_SPINGGA(rhocm,rhos,grdcc,grdss,grdcs,mu,Em000)
      call ESRC_SPINGGA(rhocp,rhos,grdcc,grdss,grdcs,mu,Ep000)
      d1E(1) = (Ep000 - Em000)/(2.d0*hrc) ! d1Edrhoc

      if (hrs .gt. 0.0d0) then
         call ESRC_SPINGGA(rhoc,rhosm,grdcc,grdss,grdcs,mu,E0m00)
         call ESRC_SPINGGA(rhoc,rhosp,grdcc,grdss,grdcs,mu,E0p00)
         d1E(2) = (E0p00 - E0m00)/(2.d0*hrs) ! d1Edrhos
      end if

      if (hgc .gt. 0.0d0) then
         call ESRC_SPINGGA(rhoc,rhos,grdccm,grdss,grdcs,mu,E00m0)
         call ESRC_SPINGGA(rhoc,rhos,grdccp,grdss,grdcs,mu,E00p0)
         d1E(3) = (E00p0 - E00m0)/(2.d0*hgc) ! d1Edgrdcc
      end if

      if (hgs .gt. 0.0d0) then
         call ESRC_SPINGGA(rhoc,rhos,grdcc,grdssm,grdcs,mu,E000m)
         call ESRC_SPINGGA(rhoc,rhos,grdcc,grdssp,grdcs,mu,E000p)
         d1E(4) = (E000p - E000m)/(2.d0*hgs) ! d1Edgrdss
      end if

      if (hgcs .gt. 0.0d0) then
         call ESRC_SPINGGA(rhoc,rhos,grdcc,grdss,grdcsm,mu,E0000m)
         call ESRC_SPINGGA(rhoc,rhos,grdcc,grdss,grdcsp,mu,E0000p)
         d1E(5) = (E0000p - E0000m)/(2.d0*hgcs) ! d1Edgrdcs
      end if

C --- Second-order derivatives --------------------------------------------
      if (norder .lt. 2) return

      d2E(:) = 0.0d0

C     derivative wrt rhoc and rhoc -----------------
      d2E(1) = (Ep000 - 2.d0*E + Em000)/(hrc*hrc) ! d2Edrhocdrhoc

      if (hrs .gt. 0.0d0) then

C        derivative wrt rhoc and rhos -----------------------
         call ESRC_SPINGGA(rhocp,rhosp,grdcc,grdss,grdcs,mu,Epp00)
         call ESRC_SPINGGA(rhocm,rhosp,grdcc,grdss,grdcs,mu,Emp00)
         call ESRC_SPINGGA(rhocp,rhosm,grdcc,grdss,grdcs,mu,Epm00)
         call ESRC_SPINGGA(rhocm,rhosm,grdcc,grdss,grdcs,mu,Emm00)
         d2E(2) = (Epp00 - Emp00 - Epm00 + Emm00)/(4.d0*hrc*hrs) ! d2Edrhocdrhos

C        derivative wrt rhos and rhos -----------------------
         d2E(3) = (E0p00 - 2.d0*E + E0m00)/(hrs*hrs) ! d2Edrhosdrhos

      end if

      if (hgc .gt. 0.0d0) then

C        derivative wrt rhoc and grdcc -----------------------
         call ESRC_SPINGGA(rhocp,rhos,grdccp,grdss,grdcs,mu,Ep0p0)
         call ESRC_SPINGGA(rhocm,rhos,grdccp,grdss,grdcs,mu,Em0p0)
         call ESRC_SPINGGA(rhocp,rhos,grdccm,grdss,grdcs,mu,Ep0m0)
         call ESRC_SPINGGA(rhocm,rhos,grdccm,grdss,grdcs,mu,Em0m0)
         d2E(4) = (Ep0p0 - Em0p0 - Ep0m0 + Em0m0)/(4.d0*hrc*hgc) ! d2Edrhocdgrdcc

         if (hrs .gt. 0.0d00) then
C        derivative wrt rhos and grdcc -----------------------
         call ESRC_SPINGGA(rhoc,rhosp,grdccp,grdss,grdcs,mu,E0pp0)
         call ESRC_SPINGGA(rhoc,rhosm,grdccp,grdss,grdcs,mu,E0mp0)
         call ESRC_SPINGGA(rhoc,rhosp,grdccm,grdss,grdcs,mu,E0pm0)
         call ESRC_SPINGGA(rhoc,rhosm,grdccm,grdss,grdcs,mu,E0mm0)
         d2E(5) = (E0pp0 - E0mp0 - E0pm0 + E0mm0)/(4.d0*hrs*hgc) ! d2Edrhosdgrdcc
         end if

C        derivative wrt grdcc and grdcc -----------------------
         d2E(6) = (E00p0 - 2.d0*E + E00m0)/(hgc*hgc) ! d2Edgrdcgrdcc

       end if ! hgc .gt. 0

       if (hgs .gt. 0.0d0) then
C        derivative wrt rhoc and grdss -----------------------
         call ESRC_SPINGGA(rhocp,rhos,grdcc,grdssp,grdcs,mu,Ep00p)
         call ESRC_SPINGGA(rhocm,rhos,grdcc,grdssp,grdcs,mu,Em00p)
         call ESRC_SPINGGA(rhocp,rhos,grdcc,grdssm,grdcs,mu,Ep00m)
         call ESRC_SPINGGA(rhocm,rhos,grdcc,grdssm,grdcs,mu,Em00m)
         d2E(7) = (Ep00p - Em00p - Ep00m + Em00m)/(4.d0*hrc*hgs) ! d2Edrhocdgrdss

       if (hrs .gt. 0.0d00) then
C        derivative wrt rhos and grdss -----------------------
         call ESRC_SPINGGA(rhoc,rhosp,grdcc,grdssp,grdcs,mu,E0p0p)
         call ESRC_SPINGGA(rhoc,rhosm,grdcc,grdssp,grdcs,mu,E0m0p)
         call ESRC_SPINGGA(rhoc,rhosp,grdcc,grdssm,grdcs,mu,E0p0m)
         call ESRC_SPINGGA(rhoc,rhosm,grdcc,grdssm,grdcs,mu,E0m0m)
         d2E(8) = (E0p0p - E0m0p - E0p0m + E0m0m)/(4.d0*hrs*hgs) ! d2Edrhosdgrdss
       end if

       if (hgc .gt. 0.0d00) then
C        derivative wrt grdcc and grdss -----------------------
         call ESRC_SPINGGA(rhoc,rhos,grdccp,grdssp,grdcs,mu,E00pp)
         call ESRC_SPINGGA(rhoc,rhos,grdccm,grdssp,grdcs,mu,E00mp)
         call ESRC_SPINGGA(rhoc,rhos,grdccp,grdssm,grdcs,mu,E00pm)
         call ESRC_SPINGGA(rhoc,rhos,grdccm,grdssm,grdcs,mu,E00mm)
         d2E(9) = (E00pp - E00mp - E00pm + E00mm)/(4.d0*hgc*hgs) ! d2Edgrdccdgrdss
       end if

C        derivative wrt grdss and grdss -----------------------
         d2E(10) = (E000p - 2.d0*E + E000m)/(hgs*hgs) ! d2Edgrdssgrdss

      end if ! hgs .gt. 0
       
      if (hgcs .gt. 0) then
       ! TODO: d2E(11:15) involving grdcs and hgcs
          call quit(
     &    'DESRC_SPINGGA error: hgcs .gt. 0 not immplemented yet')
      end if ! hgcs .gt. 0

      return
      end ! SUBROUTINE DESRC_SPINGGA


C*****************************************************************************
      SUBROUTINE DESRX_SPINGGA(ESRX_GGA,rho,grd2,mu,norder,E,d1E,d2E)
C*****************************************************************************
C     Compute first-order and/or second-order numerical derivatives
C     of a generic spin-dependent short-range GGA exchange functional
C     using ESRX_GGA(rhoa,grdaa,mu,Ea) and ESRX_GGA(rhob,grdbb,mu,Eb)
C
C     Exchange derivatives will be done based on alpha and beta values
C     because there are no mixed alpha-beta terms for exchange.
C     The routine DESRC_SPINGGA performes the derivatives with respect to charge and spin variables
C     for correlation functionals.
C
C Input: ESRX_GGA  : GGA exchange functional routine for ESRX_GGA(rhoa,grdaa,mu,Ea) and ESRX_GGA(rhob,grdbb,mu,Eb)
C        rho(1:4)  : total density, spin density, rho_a, rho_b (rho(3:4) used)
C        grd2(1:6) : grdcc, grdss, grdcs, grdaa, grdbbb, grdab (grd2(4:5) used)
C        mu        : Interaction parameter
C        norder    : order of derivatives
C
C Output:    d1E(1) = d1Edrhoc        : first derivative wrt rhoc
C            d1E(2) = d1Edrhos        : first derivative wrt rhos
C            d1E(3) = d1Edgrdcc       : first derivative wrt grdcc
C            d1E(4) = d1Edgrdss       : first derivative wrt grdss
C            d1E(5) = d1Edgrdcs       : first derivative wrt grdcs
C
C            d2E(1) = d2Edrhocdrhoc   : second derivative wrt rhoc and rhoc
C            d2E(2) = d2Edrhocdrhos   : second derivative wrt rhoc and rhos
C            d2E(3) = d2Edrhosdrhos   : second derivative wrt rhos and rhos
C            d2E(4) = d2Edrhocdgrdcc  : second derivative wrt rhoc and grdcc
C            d2E(5) = d2Edrhosdgrdcc  : second derivative wrt rhos and grdcc
C            d2E(6) = d2Edgrdccdgrdcc : second derivative wrt grdcc and grdcc
C            d2E(7) = d2Edrhocdgrdss  : second derivative wrt rhoc and grdss
C            d2E(8) = d2Edrhosdgrdss  : second derivative wrt rhos and grdss
C            d2E(9) = d2Edgrdccdgrdss : second derivative wrt grdcc and grdss
C            d2E(10)= d2Edgrdssdgrdss : second derivative wrt grdss and grdss
C            d2E(11)= d2Edrhocdgrdcs
C            d2E(12)= d2Edrhosdgrdcs
C            d2E(13)= d2Edgrdccdgrdcs
C            d2E(14)= d2Edgrdssdgrdcs
C            d2E(15)= d2Edgrdcsdgrdcs
C
C     Created: Apr-18, H. J. Aa. Jensen (based on previous DESRSPINGGA routine, now removed)
C*****************************************************************************
      implicit none

      external ESRX_GGA
      integer, intent(in) :: norder
      real*8, intent(in) :: rho(4), grd2(6), mu
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: Ea, Eb, d1W(4), d2W(6) ! for derivative wrt rhoa, rhob, grdaa, grdbb

      real*8 :: rhoa,  rhoam,  rhoap,  hra
      real*8 :: grdaa, grdaam, grdaap, hga

      real*8 :: rhob,  rhobm,  rhobp,  hrb
      real*8 :: grdbb, grdbbp, grdbbm, hgb

      real*8 :: Em0, Ep0, E0m, E0p
      real*8 :: Emm, Epp, Epm, Emp

      real*8, parameter :: D4I = 0.25d0, D8I = 0.125d0, D16I = 0.0625d0

      logical :: Eb_eq_Ea

C     check
      if(norder .lt. 0 .or. norder .ge. 3) then
        print *,'Fatal error DESRX_SPINGGA: norder=',norder
        call quit('ERROR DESRX_SPINGGA: norder must be 0, 1, or 2')
      endif

      if (norder .ge. 1) d1E(:) = 0.0d0
      if (norder .ge. 2) d2E(:) = 0.0d0

      if (rho(1) .lt. 1.d-10) then ! no significant density in this point
         E = 0.0d0
         return
      end if

      rhoa  = rho(3)
      rhob  = rho(4)
      grdaa = grd2(4)
      grdbb = grd2(5)
      if (rhoa .lt. 0.0d0 .or. rhob .lt. 0.0d0) then
         call quit('DESRX_SPINGGA error: rhoa or rhob negative')
      end if
      if (grdaa .lt. 0.0d0 .or. grdbb .lt. 0.0d0) then
         call quit('DESRX_SPINGGA error: grdaa or grdbb negative')
      end if
      Eb_eq_Ea  = (abs(rhoa-rhob)   .le. 1.d-14
     &       .and. abs(grdaa-grdbb) .le. 1.d-14)

C     energy
      call ESRX_GGA(rhoa,grdaa,mu,Ea)
      if (Eb_eq_Ea) then
         E = 2.0d0 * Ea
      else
         call ESRX_GGA(rhob,grdbb,mu,Eb)
         E = Ea + Eb
      end if

C     (optimal?) numerical step
C     (step length of hx=1.d-4*X should give a relative accuracy of approx. 16 digits
C      because we use central 2nd-order approximation=)

      hra = 1.d-4*rhoa
      hrb = 1.d-4*rhob
      hga = 1.d-4*grdaa
      hgb = 1.d-4*grdbb
      if (hra .lt. 1.d-12) hra = 0.d0
      if (hrb .lt. 1.d-12) hrb = 0.d0
      if (hga .lt. 1.d-24) hga = 0.d0
      if (hgb .lt. 1.d-24) hgb = 0.d0

C     variations of variables for numerical differentiation
      rhoam  = rhoa-hra
      rhoap  = rhoa+hra
      rhobm  = rhob-hrb
      rhobp  = rhob+hrb
      grdaam = grdaa-hga
      grdaap = grdaa+hga
      grdbbm = grdbb-hgb
      grdbbp = grdbb+hgb

! (All derivatives in the following calculated with central 2nd-order approximation)

C First-order derivatives --------------------------------------------
      if (norder .lt. 1) return

      d1W(:) = 0.0d0
      d2W(:) = 0.0d0

C     derivatives wrt rhoa, rhob, grdaa, grdbb ---------------------

      if (hra .gt. 0.0d0) then
         call ESRX_GGA(rhoam,grdaa,mu,Em0)
         call ESRX_GGA(rhoap,grdaa,mu,Ep0)
         d1W(1) = (Ep0 - Em0)/(2.d0*hra) ! d1Edrhoa
         d2W(1) = (Ep0 - 2.d0*Ea + Em0)/(hra*hra) ! d2Edrhoadrhoa
      end if

      if (hga .gt. 0.0d0) then
         call ESRX_GGA(rhoa,grdaam,mu,E0m)
         call ESRX_GGA(rhoa,grdaap,mu,E0p)
         d1W(3) = (E0p - E0m)/(2.d0*hga) ! d1Edgrdaa
         d2W(3) = (E0p - 2.d0*Ea + E0m)/(hga*hga) ! d2Edgrdaagrdaa
      end if

      IF (Eb_eq_Ea) THEN
         d1W(2) = d1W(1)
         d2W(2) = d2W(1)
         d1W(4) = d1W(3)
         d2W(4) = d2W(3)
      ELSE
       if (hrb .gt. 0.0d0) then
         call ESRX_GGA(rhobm,grdbb,mu,Em0)
         call ESRX_GGA(rhobp,grdbb,mu,Ep0)
         d1W(2) = (Ep0 - Em0)/(2.d0*hrb) ! d1Edrhob
         d2W(2) = (Ep0 - 2.d0*Eb + Em0)/(hrb*hrb) ! d2Edrhobdrhob
       end if

       if (hgb .gt. 0.0d0) then
         call ESRX_GGA(rhob,grdbbm,mu,E0m)
         call ESRX_GGA(rhob,grdbbp,mu,E0p)
         d1W(4) = (E0p - E0m)/(2.d0*hgb) ! d1Edgrdbb
         d2W(4) = (E0p - 2.d0*Eb + E0m)/(hgb*hgb) ! d2Edgrdbbdgrdbb
       end if

      END IF ! Eb_eq_Ea

       d1E(1) = 0.50d0*(d1W(1) + d1W(2)) ! d1Edrhoc
       d1E(2) = 0.50d0*(d1W(1) - d1W(2)) ! d1Edrhos
       d1E(3) = 0.25d0*(d1W(3) + d1W(4)) ! d1Edgrdcc
       d1E(4) = d1E(3)                   ! d1Edgrdss
       d1E(5) = 0.50d0*(d1W(3) - d1W(4)) ! d1Edgrdcs


C Second-order derivatives --------------------------------------------
      if (norder .lt. 2) return

C      cross-derivative wrt rhoa and grdaa -----------------------
       IF (hga .gt. 0.0d0) then
       if (hra .gt. 0.0d0) then
         call ESRX_GGA(rhoap,grdaap,mu,Epp)
         call ESRX_GGA(rhoam,grdaap,mu,Emp)
         call ESRX_GGA(rhoap,grdaam,mu,Epm)
         call ESRX_GGA(rhoam,grdaam,mu,Emm)
         d2W(5) = (Epp - Emp - Epm + Emm)/(4.d0*hra*hga) ! d2Edrhoadgrdaa
       end if
       END IF

C      cross-derivative wrt rhob and grdbb -----------------------
       IF (Eb_eq_Ea) THEN
          d2W(6) = d2W(5)
       ELSE
        IF (hgb .gt. 0.0d0 .and. hrb .gt. 0.0d0) then
          call ESRX_GGA(rhobp,grdbbp,mu,Epp)
          call ESRX_GGA(rhobm,grdbbp,mu,Emp)
          call ESRX_GGA(rhobp,grdbbm,mu,Epm)
          call ESRX_GGA(rhobm,grdbbm,mu,Emm)
          d2W(6) = (Epp - Emp - Epm + Emm)/(4.d0*hrb*hgb) ! d2Edrhobdgrdbb
        END IF
       END IF ! Eb_eq_Ea

! The following transformations from d2E(a,b) to d2E(c,s) are only valid for exchange energies
! for which there never are cross-derivatives terms between a and b.

       d2E( 1)   = D4I  * (d2W(1) + d2W(2)) !      d2Edrc2   = (d2Edra2 + d2Edrb2) / 4
       d2E( 2)   = D4I  * (d2W(1) - d2W(2)) !      d2Edrcdrs = (d2Edra2 - d2Edrb2) / 4
       d2E( 3)   = d2E(1)                   !      d2Edrs2   = d2Edrc2

       IF (hga .gt. 0.0d0 .or. hgb .gt. 0.0d0) THEN

       d2E( 4)   = D8I  * (d2W(5) + d2W(6)) !      d2Edrcdgcc = (d2Edradgaa + d2Edrbdgbb) / 8
       d2E( 5)   = D8I  * (d2W(5) - d2W(6)) !      d2Edrsdgcc = (d2Edradgaa - d2Edrbdgbb) / 8
       d2E( 6)   = D16I * (d2W(3) + d2W(4)) !      d2Edgcc2   = (d2Edgaa2 + d2Edgbb2) / 16
       d2E( 7)   = d2E(4)                   !      d2Edrcdgss = d2Edrcdgcc
       d2E( 8)   = d2E(5)                   !      d2Edrsdgss = d2Edrsdgcc
       d2E( 9)   = d2E(6)                   !      d2Edgccdgss= d2Edgcc2
       d2E(10)   = d2E(6)                   !      d2Edgss2   = d2Edgcc2
       d2E(11)   = 2.0d0 * d2E(5)           !      d2Edrcgcs  = (d2Edradgaa - d2Edrbdgbb) / 4 = 2 d2Edrsdgcc
       d2E(12)   = 2.0d0 * d2E(4)           !      d2Edrsgcs  = (d2Edradgaa + d2Edrbdgbb) / 4 = 2 d2Edrcdgcc
       d2E(13)   = D8I  * (d2W(3) - d2W(4)) !      d2Edgccdgcs= (d2Edgaa2 - d2Edgbb2) / 8
       d2E(14)   = d2E(13)                  !      d2Edgssdgcs= d2EDgccdgcs
       d2E(15)   = 4.0d0 * d2E(6)           !      d2Edgcs2   = (d2Edgaa2 + d2Edgbb2) / 4 = 4 d2Edgc2

       END IF

      return
      end ! SUBROUTINE DESRX_SPINGGA


C*****************************************************************************
      SUBROUTINE DELSPIN(ESRSPIN,rhoc,rhos,VLAMBDA,dospin,
     &                   norder,d1E,d2E)
C*****************************************************************************
C linear complement spin-dependent exchange LDA functional, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      implicit none

      external ESRSPIN
      integer, intent(in) :: norder
      logical, intent(in) :: dospin
      real*8, intent(in) :: rhoc, rhos, VLAMBDA
      real*8, intent(out) :: d1E(9), d2E(45)
      real*8 :: E

      call DESRSPIN(ESRSPIN,rhoc,rhos,0.d0,dospin,norder,E,d1E,d2E)

      d1E(1) = (1.d0- VLAMBDA) * d1E(1)

      if (dospin) d1E(2) = (1.d0- VLAMBDA) * d1E(2)

      d2E(1) = (1.d0- VLAMBDA) * d2E(1)

      if (dospin) d2E(2) = (1.d0 - VLAMBDA) * d2E(2)
      if (dospin) d2E(3) = (1.d0 - VLAMBDA) * d2E(3)

      return
      end


C*****************************************************************************
      SUBROUTINE DELSPINnSC(ESRSPIN,rhoc,rhos,VLAMBDA,dospin,
     &                      norder,d1E,d2E)
C*****************************************************************************
C linear complement non-scaled spin-dependent correlation LDA functional, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      implicit none

      external ESRSPIN
      integer, intent(in) :: norder
      logical, intent(in) :: dospin
      real*8, intent(in) :: rhoc, rhos, VLAMBDA
      real*8, intent(out) :: d1E(9), d2E(45)
      real*8 :: E

      call DESRSPIN(ESRSPIN,rhoc,rhos,0.d0,dospin,norder,E,d1E,d2E)

      d1E(1) = (1.d0- VLAMBDA**2)* d1E(1)
      if (dospin) d1E(2) = (1.d0- VLAMBDA**2)* d1E(2)

      d2E(1) = (1.d0- VLAMBDA**2)* d2E(1)
      if (dospin) d2E(2) = (1.d0- VLAMBDA**2)* d2E(2)
      if (dospin) d2E(3) = (1.d0- VLAMBDA**2)* d2E(3)

      return
      end


C*****************************************************************************
      subroutine ESRX_PBETCSERF(rhoa,grdaa,mu,e)
C*****************************************************************************
C     Short-range PBE exchange functional from 
C       Toulouse, Colonna, Savin JCP 122, 014110 (2005)
C
C     Input: rhoa   : alpha or beta density
C            grdaa  : grad(rhoa) . grad(rhoa)
C            mu     : Interaction parameter
C
C     Ouput: e      : PBETCS alpha or beta exchange energy
C
C     Author: J. Toulouse
C     Date  : 19-02-05; revised 22-Apr-18 hjaaj
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8 :: exerflda, exerfpbe
      real*8 :: rho, rho13, grd2, t1, t2, t3, t4

      real*8 :: XKSR(3)
      logical :: ERFEXP(0:2)

! function 
      real*8 :: berft

C     LDA energy
      rho = 2.0d0*rhoa ! VXSRLDA assumes rhoa .eq. rhob and rho = rhoc = rhoa + rhob
      rho13 = rho**(1.d0/3.d0)
      ERFEXP(0:2) = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      exerflda = XKSR(1)

      grd2 = 4.0d0*grdaa
!     Code generated by Mathematica
      exerfpbe = (-5.521381337364587693627982d1*exerflda*(1.63580000000
     & 000014281909d0*rho**4 + 4.272901476924512165211922d-2*grd2*berft
     &  (1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(4.d0/
     &  3.d0) + 1.d0*exerflda*rho**(8.d0/3.d0)))/(-9.0318755916609930523
     & 0091d1*rho**4 + 1.442249570307408301772512d0*grd2*exerflda*berft
     &  (1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)) - 5.5213813
     &  37364587693627982d1*exerflda*rho**(8.d0/3.d0))

      ! exerfpbe is total exchange energy for rhoa .eq. rhob;
      ! thus divide by 2 to get alpha (or beta) exchange energy
      e = 0.5d0*exerfpbe

      return
      end


C*****************************************************************************
      subroutine ESRX_PBEHSEERF(rhoa,grdaa,mu,e)
C*****************************************************************************
C     Short-range PBE exchange functional from 
C       Heyd and Scuseria, JCP 120, 7274 (2004)
C
C     Input: rhoa   : alpha or beta density
C            grdaa  : grad(rhoa) . grad(rhoa)
C            mu     : Interaction parameter
C
C     Ouput: e      : PBEHSE alpha or beta exchange energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05; revised 22-Apr-18 hjaaj
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8 :: rho13, kF, s
      real*8 :: rho, grd, exLDA, FxHSE

C     LDA exchange
      rho = 2.0d0*rhoa
      rho13 = rho**(1.d0/3.d0)
      CALL EDRC(exLDA,rho,rho13)

C     Gradient
      kF = (3.d0*(pi**2)*rho)**(1.d0/3.d0)
      grd = 2.0d0*sqrt(grdaa)
      s = grd/(2.d0*kF*rho)

C     Enhancement factor (over LDA) of short-range PBE from HSE
      call wpbe_analytical_erfc_approx(rho,s,mu,FxHSE)

C     Energy
      e = 0.5d0*exLDA*FxHSE
      ! exLDA*FxHSE is total exchange energy for rhoa .eq. rhob;
      ! thus divide by 2 to get alpha (or beta) exchange energy
       
      return
      end


C*****************************************************************************
      subroutine wpbe_analytical_erfc_approx(rho,s,omega,Fx_wpbe)
C*****************************************************************************
C    Short-range exchange PBE enhancement factor
C    from Heyd and Scuseria, JCP 120, 7274 (2004)
C
C    The Erfc function is approximated instead of the PBE exchange hole
C    to perform the integration over r12
C    
C    'wPBE Enhancement Factor (erfc approx.,analytical, no gradients)'
C    routine from J. Heyd
C
C*****************************************************************************
      Implicit None
#include "pi.h"
      Real*8 rho,s,omega,Fx_wpbe

      Real*8 f12,f13,f14,f18,f23,f43,f32,f72,f34,f94,f1516,f98
      Real*8 pi_23
      Real*8 Three_13
    
      Real*8 ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8
      Real*8 eb1
      Real*8 A,B,C,D,E
      Real*8 Ha1,Ha2,Ha3,Ha4,Ha5
      Real*8 Fc1,Fc2
      Real*8 EGa1,EGa2,EGa3
      Real*8 EGscut,wcutoff,expfcutoff

      Real*8 xkf, xkfrho
      Real*8 w,w2,w3,w4,w5,w6,w7,w8
      Real*8 A2,A3,A4,A12,A32,A52,A72
      Real*8 X
      Real*8 s2,s3,s4,s5,s6

      Real*8 H,xF,C_0,H_0,DH_0,DH_072
      Real*8 aDH,aDH2,aDHw2,sraDH,sraDH2,sraDH3,sraDH5
      Real*8 G_a,G_b,EG,E_0

      Real*8 Hsbw,Hsbw2,Hsbw3,Hsbw12,Hsbw32,Hsbw52
      Real*8 DHsbw,DHsbw2,DHsbw3,DHsbw4
      Real*8 DHsbw12,DHsbw32,DHsbw52,DHsbw72
      Real*8 HsbwA94,HsbwA9412
      Real*8 HsbwA942,HsbwA943,HsbwA945
      Real*8 experf,expei
      Real*8 expei1,expei2,expei3,expei4

      Real*8 np1,np2
      Real*8 factor1,t1,t2f1,t3f1,t4f1,t9,t10f1
      Real*8 t5f1,t6f1,t7f1,t8f1,t9f1

      Real*8 atb3H,atb4H

      Real*8 term1,term2,term3,term4
c     Real*8 term1p,p0,p1,p2,p3,p4,p5,p6,p7,p8

      Real*8 ax,um,uk,ul
      Real*8 gc1,gc2

      Real*8 ei,fexpei

      Real*8 Zero,One,Two,Three,Four,Six,Eight,Nine,Ten
      Real*8 Fifteen,Sixteen
      Real*8 r64,r36,r81,r256,r384,r864,r1944,r4374
      Real*8 r27,r48,r120,r128,r144,r288,r324,r512,r729
      Real*8 r32,r243,r2187,r6561,r40

      Save Zero,One,Two,Three,Four,Six,Eight,Nine,Ten
      Data Zero,One,Two,Three,Four,Six,Eight,Nine,Ten
     &  / 0D0,1D0,2D0,3D0,4D0,6D0,8D0,9D0,10D0 /
      Save Fifteen,Sixteen
      Data Fifteen,Sixteen / 1.5D1, 1.6D1 /
      Save r36,r64,r81,r256,r384,r864,r1944,r4374
      Data r36,r64,r81,r256,r384,r864,r1944,r4374
     &  / 3.6D1,6.4D1,8.1D1,2.56D2,3.84D2,8.64D2,1.944D3,4.374D3 /
      Save r27,r48,r120,r128,r144,r288,r324,r512,r729
      Data r27,r48,r120,r128,r144,r288,r324,r512,r729
     &  / 2.7D1,4.8D1,1.2D2,1.28D2,1.44D2,2.88D2,3.24D2,5.12D2,7.29D2 /
      Save r32,r243,r2187,r6561,r40
      Data r32,r243,r2187,r6561,r40
     &  / 3.2D1,2.43D2,2.187D3,6.561D3,4.0d1 /

c     General constants

      f12    = 0.5d0
      f13    = One/Three
      f14    = 0.25d0
      f18    = 0.125d0

      f23    = Two * f13
      f43    = Two * f23

      f32    = 1.5d0
      f72    = 3.5d0
      f34    = 0.75d0
      f94    = 2.25d0
      f98    = 1.125d0
      f1516  = Fifteen / Sixteen

      pi_23  = pi2**f13

      Three_13 = Three**f13

c     Constants from fit

      ea1 = -1.128223946706117d0
      ea2 = 1.452736265762971d0
      ea3 = -1.243162299390327d0
      ea4 = 0.971824836115601d0
      ea5 = -0.568861079687373d0
      ea6 = 0.246880514820192d0
      ea7 = -0.065032363850763d0
      ea8 = 0.008401793031216d0

      eb1 = 1.455915450052607d0

c     Constants for PBE hole

      A      =  1.0161144d0
      B      = -3.7170836d-1
      C      = -7.7215461d-2
      D      =  5.7786348d-1
      E      = -5.1955731d-2
      X      = - Eight/Nine

c     Constants for fit of H(s) (PBE)

      Ha1    = 9.79681d-3
      Ha2    = 4.10834d-2
      Ha3    = 1.87440d-1
      Ha4    = 1.20824d-3
      Ha5    = 3.47188d-2

c     Constants for F(H) (PBE)

      Fc1    = 6.4753871d0
      Fc2    = 4.7965830d-1

c     Constants for polynomial expansion for EG for small s

      EGa1   = -2.628417880d-2
      EGa2   = -7.117647788d-2
      EGa3   =  8.534541323d-2

c     Constants for large x expansion of exp(x)*ei(-x)

      expei1 = 4.03640D0
      expei2 = 1.15198D0
      expei3 = 5.03627D0
      expei4 = 4.19160D0

c     Cutoff criterion below which to use polynomial expansion

      EGscut     = 8.0d-2
      wcutoff    = 1.4D1
      expfcutoff = 7.0D2

c     Calculate prelim variables

      xkf    = (Three*pi2*rho) ** f13
      xkfrho = xkf * rho

      A2 = A*A
      A3 = A2*A
      A4 = A3*A
      A12 = Sqrt(A)
      A32 = A12*A
      A52 = A32*A
      A72 = A52*A

      w      = omega / xkf
      w2    = w * w
      w3    = w2 * w
      w4    = w2 * w2
      w5    = w3 * w2
      w6    = w5 * w
      w7    = w6 * w
      w8    = w7 * w

      X      = - Eight/Nine

      s2     = s*s
      s3     = s2*s
      s4     = s2*s2
      s5     = s4*s
      s6     = s5*s

c     Calculate wPBE enhancement factor

      H      = (Ha1*s2 + Ha2*s4) /
     &         (One + Ha3*s4 + Ha4*s4*s + Ha5*s4*s2)

      xF     = Fc1*H + Fc2

      C_0    = C * (One + s2*xF)
      H_0    = s2 * H
      DH_0   = D + H_0
      DH_072 = DH_0**f72

      aDH    = w2 + D + H_0
      aDH2   = aDH*aDH
      aDHw2  = aDH - w2
      sraDH  = sqrt(aDH)
      sraDH2 = sraDH*sraDH
      sraDH3 = sraDH2*sraDH
      sraDH5 = sraDH2*sraDH3

      if(s .gt. EGscut) then

        G_a    = sqrtpi * (Fifteen*E + Six*C_0*DH_0 +
     &                   Four*B*(DH_0**2) + Eight*A*(DH_0**3))
     &                * (One / (Sixteen * DH_072))
     &                 - f34*pi*sqrt(A) * exp(f94*H_0/A) *
     &                   (One - erf(f32*s*sqrt(H/A)))

        G_b    = (f1516 * sqrtpi * s2) / DH_072

        EG     = - (f34*pi + G_a) / G_b

      else

        EG = EGa1 + EGa2*s2 + EGa3*s4

      endif

      E_0    = E + s2*EG

c     Change exponent of Gaussian if we're using the simple approx.

      if(w .gt. wcutoff) then

        eb1 = 2.0d0

      endif

c     Calculate helper variables (should be moved later on...)

      Hsbw = s2*H + eb1*w2
      Hsbw2 = Hsbw*Hsbw
      Hsbw3 = Hsbw2*Hsbw
      Hsbw12 = Sqrt(Hsbw)
      Hsbw32 = Hsbw12*Hsbw
      Hsbw52 = Hsbw32*Hsbw
      
      DHsbw = D + s2*H + eb1*w2
      DHsbw2 = DHsbw*DHsbw
      DHsbw3 = DHsbw2*DHsbw
      DHsbw4 = DHsbw3*DHsbw
      DHsbw12 = Sqrt(DHsbw)
      DHsbw32 = DHsbw12*DHsbw
      DHsbw52 = DHsbw32*DHsbw
      DHsbw72 = DHsbw52*DHsbw
      
      HsbwA94   = f94 * Hsbw / A
      HsbwA942  = HsbwA94*HsbwA94
      HsbwA943  = HsbwA942*HsbwA94
      HsbwA945  = HsbwA943*HsbwA942
      HsbwA9412 = Sqrt(HsbwA94)

c    Calculate the terms needed in any case

      term2  = f12 * (B * (w - sraDH)) / 
     &               (sraDH * aDHw2)

      term3  = f14 * (C_0 * (Three*aDH*w - w3 - Two*sraDH3)) /
     &               (sraDH3 * aDHw2**2)

      term4  = f18 * (E_0 * (Ten*aDH*w3 - Three*w5 -
     &               Fifteen*aDH2*w + Eight*sraDH5)) /
     &               (sraDH5 * aDHw2**3) 


      t10f1 = (f12)*A*Log(Hsbw / DHsbw)

c     Calculate exp(x)*f(x) depending on size of x

      if(HsbwA94 .lt. expfcutoff) then

        experf = Exp(HsbwA94)*erfc(HsbwA9412)
        expei  = Exp(HsbwA94)*Ei(-HsbwA94)

      else

        experf = One/(sqrtpi*HsbwA9412) 
     &           - One/(Two*Sqrt(pi*HsbwA943))
     &           + Three/(Four*Sqrt(pi*HsbwA945))

        expei  = - (One/HsbwA94) *
     &             (HsbwA942 + expei1*HsbwA94 + expei2) /
     &             (HsbwA942 + expei3*HsbwA94 + expei4)

      endif

      if (w .eq. Zero) then

c       Fall back to original expression for the PBE hole

        t1 = -f12*A*expei

        if(s .gt. 0.0D0) then

          t10f1 = (f12)*A*(Log(Hsbw) - Log(DHsbw))

          term1 = t1 + t10f1

          Fx_wpbe = X * (term1 - term2 - term3 + term4)

        else

c         term1 = f12*A* -8.23742D-1

          Fx_wpbe = 1.0d0

        endif


      elseif(w .gt. wcutoff) then

c       Use simple Gaussian approximation for large w

        term1 = -f12*A*(expei+log(DHsbw)-log(Hsbw))

c       Fx_wpbe = X * (term1 - term2 - term3 + term4)
        Fx_wpbe = term1

      else

c       For everything else, use the full blown expression

        np1 = -f32*ea1*A12*w + r27*ea3*w3/(Eight*A12)
     &        - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)

        np2 = -A + f94*ea2*w2 - r81*ea4*w4/(Sixteen*A)
     &        + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)


        t1 = f12*(np1*pi*experf + np2*expei)

        factor1 = One / (r128*DHsbw4*Hsbw3*A2)
        
        t2f1 = (f12)*ea1*sqrtpi*A*w / DHsbw12
        
        t3f1 = (f12)*ea2*A*w2 / DHsbw
        
        t4f1 = sqrtpi*(-f98*ea3 / Hsbw12 
     &           + f14*ea3*A / DHsbw32)*w3
        
        t5f1 = (One/r128) * (-r144*(One/Hsbw) 
     &        + r64*(One/DHsbw2)*A)*ea4*w4
        
        t6f1 = (Three*ea5*sqrtpi*(Three*DHsbw52*(Nine*Hsbw-Two*A)
     &           + Four*Hsbw32*A2) * w5)
     &           / (r32*DHsbw52*Hsbw32*A)
        
        t7f1 = (ea6*((r32*A)/DHsbw3
     &         + (-r36 + (r81*s2*H)/A)/Hsbw2) * w6) / r32
        
        t8f1 = (-Three*ea7*sqrtpi*(-r40*Hsbw52*A3
     &          +Nine*DHsbw72*(r27*Hsbw2-Six*Hsbw*A+Four*A2)) * w7)
     &         / (r128 * DHsbw72*Hsbw52*A2)
        
        t9 = (r324*ea6*eb1*DHsbw4*Hsbw*A 
     &        + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2
     &        + r324*Hsbw*A - r288*A2)))*w8
      
c       The final value of term1 for 0 < omega < wcutoff is:

        term1 = t1 + t2f1 + t3f1 + t4f1 + t5f1 + t6f1 + t7f1 + t8f1
     &             + factor1*t9 + t10f1

        Fx_wpbe = X * (term1 - term2 - term3 + term4)
c       Fx_wpbe = term1

      endif

      end


C*****************************************************************************
      subroutine ECPBE(rho,grdcc,e)
C*****************************************************************************
C     PBE correlation functional
C
C     Input: rho    : density
C            grdcc  : grad(rho) . grad(rho)
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rho, grdcc
      real*8, intent(out) :: e

      real*8 :: rho13, eLDA, epsLDA
      real*8 :: kF, ks, t2
      real*8 :: A, At2
      real*8 :: frac, H

C     PBE constants
      real*8, parameter :: gam = 0.031090690869654897d0
      real*8, parameter :: bet = 0.066725d0

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      call EVWN(eLDA,rho,rho13)
      epsLDA = eLDA/rho

      A = exp(-epsLDA/gam) - 1.d0

      if (abs(A) .lt. 1.0d-14) then
         H = 0.0d0
      else
         A = (bet/gam)/A
      
         kF = (3.d0*pi2*rho)**(1.d0/3.d0)
         ks = sqrt(4.d0*kF/pi) 
         t2  = grdcc / (2.d0*ks*rho)**2
         At2 = A*t2

         frac = 1.d0 + At2
         frac = frac/(frac + At2*At2)
         H = gam*log(1.d0+(bet/gam)*t2*frac)
      end if

      e = rho*(epsLDA + H)

      return
      end


C*****************************************************************************
      subroutine ESRC_PBERIERF(rho,grdcc,mu,e)
C*****************************************************************************
C     rational interpolation for short-range erf correlation functional
C     Toulouse, Colonna, Savin PRA 70, 062505 (2004)
C
C     Input: rho    : density
C            grdcc  : grad(rho) . grad(rho)
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rho, grdcc, mu
      real*8, intent(out) :: e

      real*8 :: mu2, rs
      real*8 :: grs12, expgrs12, g0
      real*8 :: e_pbe, eps_pbe, eps_LDA_fac, mu_ref
      real*8 :: d1, d2, denom, p6, q7, q8, d2od1, mu6

C     General constants
      real*8 :: mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srcpbe_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4

C     Constants for g0
      real*8, parameter :: d = 32.D0/(3.D0*pi)
      real*8, parameter :: a = 3.2581D0
      real*8, parameter :: bet = 163.44D0
      real*8, parameter :: gam = 4.7125D0
      real*8, parameter :: sqrt2 = sqrt(2.0d0)
      
!     Variables
      rs = (3.d0/(4.d0*pi*rho))**(1.d0/3.d0)
      mu2 =mu*mu

!     PBE correlation energy
      call ECPBE(rho,grdcc,e_pbe)
      if (abs(e_pbe) .lt. 1.0D-15) then
         e = 0.0d0
         return
      end if
      eps_pbe = e_pbe/rho

!     g0 from Burke, Perdew, Ernzerhof
      grs12 = sqrt(gam+rs)
      expgrs12 = exp(-a*grs12)
      g0 = d*(grs12**3+bet)*expgrs12

C    modified d1 coefficient from Gori-Giorgi and Savin (submitted)
C    added by Emmanuel Fromager 13.12.05    
C      d2 = 2.d0*eps_pbe/(pi*rho*(g0 - 0.5d0))
C      d1 = -2.d0*(d2**2)*sqrtpi*rho*g0/(3.d0*eps_pbe)
Cefr   d1 = -2.d0*(d2**2)*sqrt2*sqrtpi*rho*g0/(3.d0*eps_pbe)
Chjaaj 4-oct-06: this d1 code below gave NaN if eps_pbe .eq. 0.0d0
C                thus the new code
      d2 = 2.d0/(pi*rho*(g0 - 0.5d0))
      d1 = -2.d0*eps_pbe*(d2**2)*sqrt2*sqrtpi*rho*g0/3.d0
      d2 = d2*eps_pbe
      if (mu_0 .eq. 0.0d0) then
         e = e_pbe
      else if (mu_0 .lt. 0.0d0) then
C        using e = e_pbe * ( 1 / (1 + d1*mu + d2*mu^2) )
C        based on the asymptotic Taylor expansion to 3. order
C        eps_LDA_asymptot3 = 1/d2 *  mu^(-2) - d1/(d2^2) * mu^(-3)
         denom = 1.d0 + d1*mu + d2*mu2
         e = e_pbe/denom
      else
         denom = 1.d0 + d1*mu + d2*mu2
         e = e_pbe
     &     * ( 1.0d0 + ((erf(mu_0*mu))**5) * (1.0D0/denom - 1.0D0) )
      end if

      return
      end


C*****************************************************************************
      subroutine ESRC_LYPRIERF(rho,grdcc,mu,e)
C*****************************************************************************
C     rational interpolation for short-range erf correlation functional
C     Toulouse, Colonna, Savin PRA 70, 062505 (2004)
C
C     Input: rho    : rhoc, charge density
C            grdcc  : grad(rho) . grad(rho)
C            mu     : Interaction parameter
C
C     Ouput: e      : LYPRI energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rho, grdcc, mu
      real*8, intent(out) :: e

      real*8 :: mu2, rs
      real*8 :: rhogrd
      real*8 :: grs12, expgrs12, g0, rho13
      real*8 :: e_lyp, eps_lyp, eps_LDA_fac, mu_ref
      real*8 :: d1, d2, denom, p6, q7, q8, d2od1, mu6

C     General constants
      real*8 :: mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srclyp_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4

C     Constants for g0
      real*8, parameter :: d = 32.D0/(3.D0*pi)
      real*8, parameter :: a = 3.2581D0
      real*8, parameter :: bet = 163.44D0
      real*8, parameter :: gam = 4.7125D0
      real*8, parameter :: sqrt2 = sqrt(2.0d0)
      
!     set code for do rational interpolation
      mu_0 = -1.0D0

!     Variables
      rs = (3.d0/(4.d0*pi*rho))**(1.d0/3.d0)
      mu2 =mu*mu

!     LYP correlation energy
      RHO13 = RHO**(1.0D0/3.0D0)
      rhogrd = sqrt(grdcc)
      CALL ELYP(e_lyp,rho,RHO13,rhogrd)
      if (abs(e_lyp) .lt. 1.0D-15) then
         e = 0.0d0
         return
      end if
      eps_lyp = e_lyp/rho

!     g0 from Burke, Perdew, Ernzerhof
      grs12 = sqrt(gam+rs)
      expgrs12 = exp(-a*grs12)
      g0 = d*(grs12**3+bet)*expgrs12


C    modified d1 coefficient from Gori-Giorgi and Savin (submitted)
C    added by Emmanuel Fromager 13.12.05    
C      d2 = 2.d0*eps_lyp/(pi*rho*(g0 - 0.5d0))
C      d1 = -2.d0*(d2**2)*sqrtpi*rho*g0/(3.d0*eps_lyp)
Cefr   d1 = -2.d0*(d2**2)*sqrt2*sqrtpi*rho*g0/(3.d0*eps_lyp)
Chjaaj 4-oct-06: this d1 code below gave NaN if eps_lyp .eq. 0.0d0
C                thus the new code
      d2 = 2.d0/(pi*rho*(g0 - 0.5d0))
      d1 = -2.d0*eps_lyp*(d2**2)*sqrt2*sqrtpi*rho*g0/3.d0
      d2 = d2*eps_lyp
      if (mu_0 .eq. 0.0d0) then
         e = e_lyp
      else if (mu_0 .lt. 0.0d0) then
C        using e = e_lyp * ( 1 / (1 + d1*mu + d2*mu^2) )
C        based on the asymptotic Taylor expansion to 3. order
C        eps_LDA_asymptot3 = 1/d2 *  mu^(-2) - d1/(d2^2) * mu^(-3)
         denom = 1.d0 + d1*mu + d2*mu2
         e = e_lyp/denom
      else
C        using e = e_lyp * ( (1 + p6*mu^6) / (1 + q7*mu^7 + q8*mu^8) )
C        based on 
C        1) the asymptotic Taylor expansion to 3. order
C           eps_LDA_asymptot3 = 1/d2 *  mu^(-2) - d1/(d2^2) * mu^(-3)
C        2) the condition that eps(mu) = eps(0) + O(mu^6)
C        2) the condition that eps(mu_0) = eps_LDA_asymptot3(mu_0)
C        4) mu_0 is determined by the inflection point
C           d2(eps_LDA_asymptot3)/dmu^2 : mu_0 = 2 d1/d2
C        (e = eps * rho)
         d2od1 = d2 / d1
C ----- first attempt. fixed mu_0, does not work!!!!
C        p6 = d2od1/d1 * mu_0_m6 - mu_0_m7 / d1 - d2od1**2 * mu_0_m4
c --- second attmpt
C        using mu_ref = alfa * d1/d2:
C        p6 = alfa**(-4) * d2od1**6 * ( alfa**(-3)*(alfa-1.0d0)*d2od1/d1
C        - 1.0d0 )
C        with alfa = 2 (the inflection point):
!        p6 = 0.0625D0 * d2od1**6 * (0.125D0*d2od1/d1 - 1.0D0 )
!        q7 = d1 * p6
!        q8 = d2 * p6
!        mu6 = mu2**3
!        e = e_lyp * (1.0D0 + p6 * mu6)
!    &             / (1.0D0 + mu6 * (q7 * mu + q8 * mu2))
c --- third attmpt
c        mu_ref = alfa * d1/d2 = 2.D0 * d1/d2 ! the inflection point
c        eps_LDA = e_lyp/d2 * ( mu_ref**(-2) - d1od2 * mu_ref**(-3))
c                = e_lyp/d2 * ( 0.25d0 * d2od1**2 - 0.125d0 * d2od1**2)
         mu_ref  = 2.0d0 / d2od1
         if (mu .ge. mu_ref) then
            e = e_lyp/(d2*mu2) * (1.0d0 - d1 / (d2*mu))
         else
c           e_lyp * (1 + p6 * mu_ref**6) =  eps_LDA(mu_ref)
            eps_LDA_fac = 0.125d0/d2 * d2od1**2
            p6 = (eps_LDA_fac - 1.0d0) * d2od1**6 * 0.015625d0
            e = e_lyp * (1.0d0 + p6 * mu2**3)
         end if
      end if

      return
      end


C*****************************************************************************
      subroutine ESRX_PBERIERF(rhoa,grdaa,mu,e)
C*****************************************************************************
C     rational interpolation for short-range erf exchange functional
C     using the two first terms of both developments at small and large mu
C     (eqs. 17 and 33 in Ref : Toulouse, Colonna, Savin PRA 70, 062505 (2004)) 
C
C     Input: rhoa   : alpha or beta density
C            grdaa  : grad(rhoa) . grad(rhoa)
C            mu     : Interaction parameter
C
C     Ouput: e      : PBERI  alpha or beta exchange energy
C
C     Author: Emmanuel Fromager (Southern Denmark University) 
C     Date  : 31-01-2006
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8 :: mu2, mu6
      real*8 :: rho, epbe
      real*8 :: d, d2, denom1, denom2
      
C     General constants
      real*8, parameter :: three53 = 3.d0**(5.d0/3.d0) 
      real*8, parameter :: pi13 = pi**(1.d0/3.d0)
      real*8, parameter :: pi43 = pi**(4.d0/3.d0) 

!     Variables
      real*8, parameter :: zero = 0.d0 
      mu2  = mu*mu
      mu6  = mu2*mu2*mu2
      rho  = 2.0d0*rhoa ! as if rhoa .eq. rhob

      d = -5.d0*(rho**(4.d0/3.d0))/(pi13*three53)      
      d2 = 20.d0/(pi43*three53*(rho**(2.d0/3.d0)))
      denom1 = 1.d0 + d2*mu2
      denom2 = 1.d0 + mu6

C Get standard PBE exchange functional calling HSEPBE with mu=0
      call ESRX_PBEHSEERF(rhoa,grdaa,zero,epbe)

      epbe = 2.0d0*epbe ! total as if rhoa .eq. rhob

      e = (d/denom1) + ((epbe + (rho*mu/sqrtpi) - d )/denom2)  
      e = 0.5d0 * e ! exchange energy for rho_sigma
C      eps_pbe = epbe/rho
      return
      end


C*****************************************************************************
      subroutine ESRC_PBELOERF(rho,grdcc,mu,e)
C*****************************************************************************
C     Our new local energy formula for sr correlation functional.
C     With PBE functional.
C
C     E_sr_c = max( E_c_pbe, E_c_srlda(mu) )
C
C     Input: rho    : density
C            grdcc  : grad(rho) . grad(rho)
C            mu     : Interaction parameter
C
C     Ouput: e      : energy
C
C     Author: M. N. Pedersen
C     Date  : Nov. 2011
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rho, grdcc, mu
      real*8, intent(out) :: e

      real*8 :: e_pbe, e_lda,p,z,Ha,Hb
      real*8 :: rho13
      real*8 :: ZKSR(3)
      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.

      rho13 = rho**(1.0d0/3.0d0)

!     PBE correlation energy (mu=0)
      call ECPBE(rho,grdcc,e_pbe)
      if (abs(e_pbe) .lt. 1.0D-15) then
         e = 0.0d0
         return
      end if

!    short-range LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      e_lda = ZKSR(1)

!      e = max(e_lda,e_pbe)
!     Instead of max(e_lda,e_pbe) we use a smoothing function HS
      z = (e_pbe - e_lda)/(e_pbe**2 + e_lda**2)

      call HS(z,1.5d0,Ha)
      call HS(-z,1.5d0,Hb)
      e = e_pbe*Ha + e_lda*Hb

      return
      end


C*****************************************************************************
      subroutine ESRC_PBEWIERF(rho,grdcc,mu,e)
C*****************************************************************************
C     weighted interpolation for for short-range erf correlation functional
C     Toulouse, Colonna, Savin PRA 70, 062505 (2004)
C
C     Input: rho    : density
C            grdcc  : grad(rho) . grad(rho)
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rho, grdcc, mu
      real*8, intent(out) :: e

      real*8 :: rho13, rs
      real*8 :: elda, epbe, eldamu
      real*8 :: ZKSR(3)
      real*8, parameter :: mu0 = 0.d0
      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.

      rho13 = rho**(1.d0/3.d0)
      rs = (3.d0/(4.d0*pi*rho))**(1.d0/3.d0)

!     LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu0,MULOCAL,0,ERFEXP)
      elda = ZKSR(1)

!     PBE energy
      call ECPBE(rho,grdcc,epbe)

!     short-range LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      eldamu = ZKSR(1)

      e = (epbe - elda)*erfc(mu*rs) + eldamu

      return
      end


C*****************************************************************************
      subroutine ESRC_PBETCSERF(rho,grdcc,mu,e)
C*****************************************************************************
C     Short-range PBE correlation functional
C       Toulouse, Colonna, Savin JCP 122, 014110 (2005)
C
C     Input: rho    : density
C            grdcc  : grad(rho) . grad(rho)
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rho, grdcc, mu
      real*8, intent(out) :: e

      real*8 :: rho13
      real*8 :: exerflda, ecerflda, ecerfpbe, t1
      real*8 :: XKSR(3), ZKSR(3)

      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.
! function 
      real*8  :: berft

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     Code generated by Mathematica
C     Energy
       t1 = 9.869604401089358618834491d0*ecerflda
     &    + 3.068528194400546905827679d-1*rho*
     & log(1.d0 - (8.401605835889726292406297d-1*grdcc*exerflda
     &  *berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &  rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*grdcc*
     & exerflda*berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )*rho**(-1.1d1/3.d0))/(-1.d0 +
     &   exp((-3.216396844291482112072009d1*ecerflda)/rho))))/
     & (1.d0 + (7.058698062165630644545894d-1*grdcc**2*exerflda**2*
     &  berft(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))**2*rho**(-2.2d1/3.d0))/
     &  (-1.d0 + exp((-3.216396844291482112072009d1*ecerflda)/rho))**2
     &  - (8.401605835889726292406297d-1*grdcc*exerflda*
     &  berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))
     &  *rho**(-1.1d1/3.d0))/(-1.d0 +
     &   exp((-3.216396844291482112072009d1*ecerflda)/rho))))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      e = ecerfpbe

      return
      end


C*****************************************************************************
      subroutine ESRC_PBETCSJERF(rho,grdcc,mu,e)
C*****************************************************************************
C     Short-range PBE correlation functional
C       Toulouse, Colonna, Savin JCP 122, 014110 (2005)
C       with modification by Jensen
C
C     Input: rho    : density
C            grdcc  : grad(rho).grad(rho)
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rho, grdcc, mu
      real*8, intent(out) :: e

      real*8 :: rho13, t1
      real*8 :: exerflda, ecerflda, ecerfpbe, epbe

      real*8 :: XKSR(3), ZKSR(3)

      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.
! function 
      real*8 :: berft

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     Code generated by Mathematica
C     Energy
      t1 = 9.869604401089358618834491d0*ecerflda
     &   + 3.068528194400546905827679d-1*rho*
     &    log(1.d0 - (8.401605835889726292406297d-1*grdcc*exerflda*
     &      berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &      rho**(-1.1d1/3.d0)*
     &      (1.d0 - (8.401605835889726292406297d-1*grdcc*exerflda*
     &      berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &      rho**(-1.1d1/3.d0))/
     &      (-1.d0 + exp((-3.216396844291482112072009d1*ecerflda)/rho))
     &      ) )/
     &      (1.d0 + (7.058698062165630644545894d-1*grdcc**2*exerflda**2*
     &      berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2
     &      *rho**(-2.2d1/3.d0))/
     &   (-1.d0 + exp((-3.216396844291482112072009d1*ecerflda)/rho))**2
     & - (8.401605835889726292406297d-1*grdcc*exerflda*
     &    berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))
     &          *rho**(-1.1d1/3.d0))/
     &   (-1.d0 + exp((-3.216396844291482112072009d1*ecerflda)/rho)) ))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      call ECPBE(rho,grdcc,epbe)

! modification of ESRC_PBETCS:
      e = max(ecerfpbe,epbe)

      return
      end


C*****************************************************************************
      real*8 pure function berft(a)
C*****************************************************************************
C  Second-order exchange gradient expansion coefficient for erf 
C  interaction
C  a = mu/(2*kF)
C
C  Author : J. Toulouse
C  Date   : 10-03-04
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: a

      if(a .lt. 0.02d0) then
!      expansion for small mu to avoid numerical problems
       berft = 7.d0/81.d0 + 56.d0*sqrt(pi)*a/243.d0 
     &        + (-128.d0/81.d0 + 448.d0*pi/729.d0)*a*a

      else

!      Code generated by Mathematica
       berft =(1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0
     &   + exp(2.5d-1/a**2)) - 2.d0*a**2*(1.1d1 + 7.d0*exp(2.5d-1/a**2
     &  ))) )/
     & (a**2*(3.2d1*a**4*(-1.d0 + exp(2.5d-1/a**2))
     &  - 3.d0*exp(2.5d-1/a**2)
     &  + 1.417963080724412821838534d1*a*erf(5.d-1/a)*exp(2.5d-1/a**2)
     &  - 8.d0*a**2*(-2.d0 + 3.d0*exp(2.5d-1/a**2)) ) )

      endif

      return
      end


C*****************************************************************************
      real*8 pure function dberftda(a)
C*****************************************************************************
!  Derivative of second-order exchange gradient
!  expansion coefficient for erf interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: a
      real*8 :: t1, t2, t3, t4, t5

      if(a .lt. 0.02d0) then
!      expansion for small mu to avoid numerical problems
       dberftda = 56.d0*sqrt(pi)/243.d0
     &        + 2.d0*(-128.d0/81.d0 + 448.d0*pi/729.d0)*a

      else

!      Code generated by Mathematica
       t1 = (1.851851851851851851851852d-2*(5.76d2*a**3*(-1.d0 + exp(2.
     &  5d-1/a**2)) + (7.d0*exp(2.5d-1/a**2))/a - 7.2d1*a*exp(2.5d-1/a
     &  **2) - 4.d0*a*(1.1d1 + 7.d0*exp(2.5d-1/a**2))))/(a**2*(3.2d1*a*
     &  *4*(-1.d0 + exp(2.5d-1/a**2)) - 3.d0*exp(2.5d-1/a**2) + 1.4179
     &  63080724412821838534d1*a*erf(5.d-1/a)*exp(2.5d-1/a**2) - 8.d0
     &  *a**2*(-2.d0 + 3.d0*exp(2.5d-1/a**2))))
       t2 = -1.851851851851851851851852d-2/a**2
       t3 = -8.d0/a + 1.28d2*a**3*(-1.d0 + exp(2.5d-1/a**2)) + (1.5d0*d
     &  exp(2.5d-1/a**2))/a**3 + (1.2d1*exp(2.5d-1/a**2))/a - 1.6d1*a*d
     &  exp(2.5d-1/a**2) + 1.417963080724412821838534d1*erf(5.d-1/a)*d
     &  exp(2.5d-1/a**2) - (7.08981540362206410919267d0*erf(5.d-1/a)*d
     &  exp(2.5d-1/a**2))/a**2 - 1.6d1*a*(-2.d0 + 3.d0*exp(2.5d-1/a**2)
     &  )
       t4 = (-1.d0 + 1.44d2*a**4*(-1.d0 + exp(2.5d-1/a**2)) - 2.d0*a**2
     &  *(1.1d1 + 7.d0*exp(2.5d-1/a**2)))/(3.2d1*a**4*(-1.d0 + exp(2.5
     &  d-1/a**2)) - 3.d0*exp(2.5d-1/a**2) + 1.417963080724412821838534
     &  d1*a*erf(5.d-1/a)*exp(2.5d-1/a**2) - 8.d0*a**2*(-2.d0 + 3.d0*
     &  exp(2.5d-1/a**2)))**2
       t5 = (-3.703703703703703703703704d-2*(-1.d0 + 1.44d2*a**4*(-1.d0 
     &  + exp(2.5d-1/a**2)) - 2.d0*a**2*(1.1d1 + 7.d0*exp(2.5d-1/a**2)
     &  )))/(a**3*(3.2d1*a**4*(-1.d0 + exp(2.5d-1/a**2)) - 3.d0*exp(2.
     &  5d-1/a**2) + 1.417963080724412821838534d1*a*erf(5.d-1/a)*exp(
     &  2.5d-1/a**2) - 8.d0*a**2*(-2.d0 + 3.d0*exp(2.5d-1/a**2))))
       dberftda = t1 + t2*t3*t4 + t5

      endif

      return
      end


C*****************************************************************************
      subroutine ESRX_PBEGWSERF(rhoa,grdaa,mu,e)
C*****************************************************************************
C     Short-range PBE exchange functional from 
C       Goll, Werner, Stoll, PCCP 7, 3917 (2005)
C
C     Input: rhoa  : alpha or beta density
C            grdaa : grad(rhoa) . grad(rhoa)
C            mu    : Interaction parameter
C
C     Ouput: e     : PBEGWS exchange energy
C
C     Author: J. Toulouse
C     Date  : 11-08-09
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8 :: rho, rho13, grd2
      real*8 :: exerflda, exerfpbe
      real*8 :: sq, fx

      real*8 :: XKSR(3)
! function 
      real*8 :: berf
! parameters
      real*8, parameter :: kappa=0.804d0
      logical, parameter :: ERFEXP(0:2) = .false.

C     LDA energy
      rho = 2.0d0*rhoa 
      rho13 = rho**(1.d0/3.d0)
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      exerflda = XKSR(1)

      grd2 = 4.0d0*grdaa ! as if grdaa .eq. grdbb
      sq=grd2*2.6121172985233599567768d-2*rho**(-8.d0/3.d0)
      fx=1.d0+kappa
     &    -kappa/(1.d0+berf(1.616204596739954813d-1
     &            *mu*rho**(-1.d0/3.d0))*sq/kappa)
      exerfpbe=exerflda*fx

      e = 0.5d0*exerfpbe ! only sigma exchange energy

      return
      end


C*****************************************************************************
      subroutine ELAX_PBEGWS(rho,grd2,VLAMBDA,e)
C*****************************************************************************
C linear complement exchange PBE energy, Kamal SHARKAS 19/05/2011
C*****************************************************************************
      real*8, intent(in) :: rho(4), grd2(6), VLAMBDA
      real*8, intent(out) :: e

      call ESRX_PBEGWSERF(rho(3),grd2(4),0.d0,e)
      e = (1.d0 - VLAMBDA) * e

      return
      end


C*****************************************************************************
      subroutine ESRC_PBEGWSERF(rhoc,grdcc,mu,e)
C*****************************************************************************
C     Modified version of Short-range PBE correlation functional from
C     Goll, Werner, Stoll, PCCP 7, 3917 (2005)
C     Based on the LDA parameterization of Perdew and Wang (PW92);
C     Perdew, Wang, PRB 45, 13244 (1992), rather than VWN
C
C     Input: rhoc   : density
C            grdcc  : grad(rhoc) . grad(rhoc)
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: E. Hedegaard (based on routine by J. Toulouse)
C     Date  : 21-03-16
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoc, grdcc, mu
      real*8, intent(out) :: e

      real*8, parameter :: tol=1d-12
      real*8 :: rho13

      real*8 :: ecerflda,eclda
      real*8 :: ecerfpbe
      real*8 :: beta, Aa, Ab, Ac, tq, arglog
      real*8, parameter :: alpha = 2.78d0
      real*8, parameter :: gamma_symbol = 3.1091d-2
 
!     short-range LDA energy (order = 0)
      rho13 = rhoc**(1.d0/3.d0)
      call ESRC_LDAERFSPIN(rhoc,0.0d0,mu,.false.,e)
      ecerflda =  e

!     full-range LDA energy
      call ESRC_LDAERFSPIN(rhoc,0.0d0,0.0d0,.false.,e)
      eclda = e

      if ((ecerflda/eclda).le.0.d0) then
         beta=0.d0
      else
         beta=6.6725d-2*(ecerflda/eclda)**alpha
      endif
      tq=grdcc*6.346820607d-2*rhoc**(-7.d0/3.d0)
      Ab=exp(-ecerflda/(rhoc*gamma_symbol))-1.d0
      if (abs(Ab).le.abs(beta*tol)) then
         ecerfpbe=ecerflda
       else
         Aa=beta/(gamma_symbol*Ab)
         if (Aa.lt.tol) Aa=tol
         Ac=1.d0+Aa*tq+Aa**2*tq**2
         arglog=1.d0+beta*(1.d0-1.d0/Ac)/(gamma_symbol*Aa)
         ecerfpbe=ecerflda+rhoc*gamma_symbol*log(arglog)
       endif

      e = ecerfpbe

      return
      end


C*****************************************************************************
      subroutine ELAC_PBEGWS(rhoc,grdcc,VLAMBDA,e)
C*****************************************************************************
C linear complement non-scaled correlation PBE functional energy, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      real*8, intent(in) :: rhoc, grdcc, VLAMBDA
      real*8, intent(out) :: e

      call ESRC_PBEGWSERF(rhoc,grdcc,0.d0,e)
      e = (1.d0 - VLAMBDA**2) * e

      return
      end


C*****************************************************************************
      subroutine ELASC_PBEGWS(rhoc,grdcc,VLAMBDA,e)
C*****************************************************************************
C linear complement Scaled correlation PBE functional energy, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      real*8, intent(in) :: rhoc, grdcc, VLAMBDA
      real*8, intent(out) :: e
      real*8 :: ecoul, escaled,rhocscaled,grdscaled

      ecoul = 0.d0

      call ESRC_PBEGWSERF(rhoc,grdcc,0.d0,ecoul)

      rhocscaled = rhoc/VLAMBDA**3
      grdscaled  = grdcc/VLAMBDA**4

      escaled = 0.d0

      call ESRC_PBEGWSERF(rhocscaled,grdscaled,0.d0,escaled)

      e = ecoul - (VLAMBDA**5)*escaled

      return
      end


C*****************************************************************************
      real*8 pure function berf(a)
C*****************************************************************************
C  Second-order exchange gradient expansion coefficient for erf
C  interaction
C  a = mu/(2*kF)
C 
C  Modified version from Goll, Werner, Stoll, PCCP 7, 3917 (2005)
C
C  Author : J. Toulouse
C  Date   : 11-08-09
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: a
      real*8 :: fak, fak2
! function
      ! eta=19.0d0
      fak=2.540118935556d0*exp(-19.0d0*a*a)

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-7.d0+72.d0*a*a)
     &        /(27.d0*(-3.d0-24.d0*a*a+32.d0*a**4+8.d0*sqrt(pi)*a))

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)
     &        - 23.d0/(358400.d0*a**6)

      else

!     Code generated by Mathematica
      fak2 = exp(0.25d0/a**2)
      berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0
     &  + fak2) - 2.d0*a**2*(1.1d1 + 7.d0*fak2 ))) /
     &  (a**2*(3.2d1*a**4*(-1.d0 + fak2) - 3.d0*fak2
     &   + 1.417963080724412821838534d1*a*erf(5.d-1/a)*fak2
     &   - 8.d0*a**2*(-2.d0 + 3.d0*fak2)))

      end if

      berf=berf*fak

      return
      end


C*****************************************************************************
      subroutine ESRC_MULOC_GGA(rho,grd,mu,e)
C*****************************************************************************
C     Modified version of local mu formula C short-range correlation functional from
C     J. Toulouse, F. Colonna and A. Savin, J. Chem. Phys., 122, 014110 (2005)
C
C     Our modification: factor alpha below
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: M.N.P
C     Date  : 20-09-11
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rho, grd, mu
      real*8, intent(out) :: e
      real*8 :: mu_new, amul
      real*8 :: alpha, x
      real*8, parameter :: a1 = 0.559463d0
      real*8, parameter :: b1 = -0.134086d0
      real*8, parameter :: c1 = 0.0296811d0
      real*8 :: rho13, ZKSR(3)

      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.
!     local mu parameter 
      rho13 = rho**(1.0d0/3.0d0)
!     Multiplication factor for mu parameter
      x = grd/rho
      alpha = a1 + b1*exp(-c1*x)
      amul = x*alpha
      mu_new = max(amul,mu)
!     short-range LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu_new,MULOCAL,0,ERFEXP)
      e = ZKSR(1)
      return
      end


C*****************************************************************************
      subroutine ESRC_MULOD_GGA(rho,grd,mu,e)
C*****************************************************************************
C     Modified version of SRCMULO_A functional
C
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: M.N.P
C     Date  : 20-09-11
C*****************************************************************************

      implicit none

      real*8, intent(in) :: rho, grd, mu
      real*8, intent(out) :: e
      real*8 :: amul, mu_new
      real*8 :: rs, x, a
      real*8, parameter :: rsfac = 0.62035049089940008660d0
      real*8, parameter :: b = 1.741930d0
      real*8, parameter :: c = -0.528057d0
      real*8, parameter :: d = 0.0720226d0
      real*8 :: rho13

      real*8 :: ZKSR(3)
      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.

      rho13 = rho**(1.0d0/3.0d0)
      rs   = rsfac/rho13
!     Multiplying factor for mu parameter
      x = grd/rho
      a = b + c*exp(-d*x)
!     local mu parameter 
      amul = a/rs
      mu_new = max(amul,mu)
!     short-range LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu_new,MULOCAL,0,ERFEXP)
      e = ZKSR(1)
      return
      end


C*****************************************************************************
      subroutine ESRC_MULOE_GGA(rhoc,grdcc,mu,e)
C*****************************************************************************
C     Modified version of SRCMULO_B correlation functional 
C
C
C
C     Input: rhoc   : density
C            grdcc  : grad(rhoc) . grad(rhoc)
C            mu     : Interaction parameter
C
C     Ouput: e      : energy
C
C     Author: M.N.P
C     Date  : 20-09-11
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoc, grdcc, mu
      real*8, intent(out) :: e
      real*8 :: amul, mu_new, x, rs, alpha
      real*8, parameter :: rsfac = 0.62035049089940008660d0
      real*8, parameter :: a = 0.878435d0
      real*8, parameter :: b = 0.5676020d0
      real*8, parameter :: bfac = 0.814516076948D0
      real*8 :: rho13

      real*8 :: ZKSR(3)
      logical, parameter :: ERFEXP(0:2) = .false.
      logical, parameter :: MULOCAL(0:2) = .false.

      rho13 = rhoc**(1.0d0/3.0d0)
      rs   = rsfac/rho13
      x = sqrt(grdcc)/rhoc
      alpha = a*x**b
!     local mu parameter 
      amul = (alpha*bfac)/sqrt(rs)
      mu_new = max(amul,mu)
!     short-range LDA energy
      call VCSRLDA(ZKSR,rhoc,rho13,mu_new,MULOCAL,0,ERFEXP)
      e = ZKSR(1)
      return
      end


C*****************************************************************************
      subroutine ESRX_LDAERF(rhoa,grdaa,mu,e)
C*****************************************************************************
C     Short-range spin-dependent LDA exchange functional
C       obtained by spin-scaling Ex[na,nb] = 1/2 (Ex[na]+Ex[nb])
C       where Ex[n] is from Toulouse, Savin, Flad, IJQC 100, 1047 (2004)
C
C     Input: rhoa   : alpha or beta density
C            grdaa  : gradient squared, zero for LDA
C            mu     : Interaction parameter
C
C     Ouput: e      : LDA exchange energy
C
C     Created: 17-08-09, J. Toulouse, subroutine adapted from Molpro
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8 :: rho, akf, a, a2, a3, ea

      real*8, parameter :: ckf = 3.093667726280136d0
      real*8, parameter :: f13 = 0.333333333333333d0

C       Density and kF
        rho = 2*rhoa ! spin scaling, rho is twice the real rhoa
        akf = ckf*(rho**f13)
        a = mu/(2.d0*akf)
        a2 = a*a
        a3 = a2*a

C       Test on the value of a

C       Limit for small a (expansion not so important as for large a)
        if (a.lt.1.d-9) then
          e = -3.d0/8.d0*rho*(24.d0*rho/pi)**f13

C       Intermediate values of a
        elseif (a.le.100) then
          e = - (rho*(24.d0*rho/pi)**f13)
     &       *(3.d0/8.d0-a*(sqrtpi*erf(0.5d0/a)+
     &       (2.d0*a-4.d0*a3)*exp(-0.25d0/a2)-3.d0*a+4.d0*a3))

C       Expansion for large a
        elseif (a.lt.1.d+9) then
          e = - (rho*(24.d0*rho/pi)**f13)*1.d0/(96.d0*a2)

C       Limit for large a
        else
          e = 0.d0
        end if

C       divide by two because of spin scaling
        e = 0.5d0 * e

      return
      end


C*****************************************************************************
      subroutine ESRX_LDAERFSPIN(rhoc,rhos,mu,dospin,e)
C*****************************************************************************
C     Short-range spin-dependent LDA exchange functional
C       obtained by spin-scaling Ex[na,nb] = 1/2 (Ex[na]+Ex[nb])
C       where Ex[n] is from Toulouse, Savin, Flad, IJQC 100, 1047 (2004)
C
C     subroutine adapted from Molpro
C
C     Input: rhoc   : total density
C            rhos   : spin density
C            mu     : Interation parameter
C            dospin : use spin density
C
C     Ouput: e      : energy
C
C     Created: 17-08-09, J. Toulouse
C              22-05-18, Re-introduced for LDA exchange as temporary
C              solution until DESRX_GGA works for both LDA and GGA.              
C*****************************************************************************
      implicit none
!pi, sqrtpi
#include "pi.h"
      logical, intent(in) :: dospin
      real*8, intent(in) :: rhoc, rhos, mu
      real*8, intent(out) :: e

      real*8 :: rho, akf, a, a2, a3
      real*8 :: rhoa, rhob
      real*8 :: ea, eb
      real*8, parameter :: ckf = 3.093667726280136d0
      real*8, parameter :: f13 = 0.333333333333333d0

C     Spin-unpolarized case
      if (.not.dospin) then
C       Density and kF
        rho = rhoc
        akf = ckf*(rho**f13)
        a = mu/(2.d0*akf)
        a2 = a*a
        a3 = a2*a

C       Test on the value of a

C       Limit for small a (expansion not so important as for large a)
        if (a.lt.1.d-9) then
          e = -3.d0/8.d0*rho*(24.d0*rho/pi)**f13

C       Intermediate values of a
        elseif (a.le.100) then
          e = - (rho*(24.d0*rho/pi)**f13)
     &       *(3.d0/8.d0-a*(sqrtpi*erf(0.5d0/a)+
     &       (2.d0*a-4.d0*a3)*exp(-0.25d0/a2)-3.d0*a+4.d0*a3))

C       Expansion for large a
        elseif (a.lt.1.d+9) then
          e = - (rho*(24.d0*rho/pi)**f13)*1.d0/(96.d0*a2)

C       Limit for large a
        else
          e = 0.d0
        end if
      
C     Spin-polarized case
      else

C       Density and kF for spin alpha
!       rhoa=max(0.5d0*(rhoc+rhos),0d0)
!       rhoa=2.d0*rhoa ! spin-scaling
        rhoa=max((rhoc+rhos),0d0) ! rhoa is twice real rho_a
        akf = ckf*(rhoa**f13)
        a = mu/(2.d0*akf)
        a2 = a*a
        a3 = a2*a

C       Test on the value of a

C       Limit for small a (expansion not so important as for large a)
        if (a.lt.1.d-9) then
          ea = -3.d0/8.d0*rhoa*(24.d0*rhoa/pi)**f13

C       Intermediate values of a
        elseif (a.le.100) then
          ea = - (rhoa*(24.d0*rhoa/pi)**f13)
     &       *(3.d0/8.d0-a*(sqrtpi*erf(0.5d0/a)+
     &       (2.d0*a-4.d0*a3)*exp(-0.25d0/a2)-3.d0*a+4.d0*a3))

C       Expansion for large a
        elseif (a.lt.1.d+9) then
          ea = - (rhoa*(24.d0*rhoa/pi)**f13)*1.d0/(96.d0*a2)

C       Limit for large a
        else
          ea = 0.d0
        end if

C       Density and kF for spin beta
!       rhob = max(0.5d0*(rhoc-rhos),0d0)
!       rhob = 2.d0*rhob ! spin-scaling
        rhob=max((rhoc-rhos),0d0) ! rhob is twice real rho_b
C
        akf = ckf*(rhob**f13)
        a = mu/(2.d0*akf)
        a2 = a*a
        a3 = a2*a

C       Test on the value of a

C       Limit for small a (expansion not so important as for large a)
        if (a.lt.1.d-9) then
          eb = -3.d0/8.d0*rhob*(24.d0*rhob/pi)**f13

C       Intermediate values of a
        elseif (a.le.100) then
          eb = - (rhob*(24.d0*rhob/pi)**f13)
     &       *(3.d0/8.d0-a*(sqrtpi*erf(0.5d0/a)+
     &       (2.d0*a-4.d0*a3)*exp(-0.25d0/a2)-3.d0*a+4.d0*a3))

C       Expansion for large a
        elseif (a.lt.1.d+9) then
          eb = - (rhob*(24.d0*rhob/pi)**f13)*1.d0/(96.d0*a2)

C       Limit for large a
        else
          eb = 0.d0
        end if

C       divide by two because of spin scaling
        e = 0.5d0 * (ea + eb)

      endif

      return
      end


C*****************************************************************************
      subroutine ELX_LDALAMSPIN(rhoc,rhos,VLAMBDA,e)
C*****************************************************************************
C linear complement spin-dependent exchange LDA functional energy, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      real*8, intent(in) :: rhoc, rhos, VLAMBDA
      real*8, intent(out) :: e
      real*8 :: rhoa, rhob, ea, eb

      rhoa=0.5d0*max((rhoc+rhos),0d0)
      call ESRX_LDAERF(rhoa,0.0d0,0.d0,ea)
      rhob=0.5d0*max((rhoc-rhos),0d0)
      call ESRX_LDAERF(rhob,0.0d0,0.d0,eb)

      e = (1.d0 - VLAMBDA) * (ea + eb)

      return
      end


C*****************************************************************************
      subroutine ESRC_LDAERFSPIN(rhoc,rhos,mu,dospin,e)
C*****************************************************************************
C     Short-range spin-dependent LDA correlation functional from 
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     subroutine adapted from Molpro
C
C     Input: rhoc   : total density
C            rhos   : spin density
C            mu     : Interaction parameter
C            dospin : use spin density
C
C     Ouput: e      : energy
C
C     Created: 18-08-09, J. Toulouse
C*****************************************************************************
      implicit none
#include "pi.h"
      logical, intent(in) :: dospin
      real*8, intent(in) :: rhoc, rhos, mu
      real*8, intent(out) :: e

      real*8 :: rs
      real*8 :: ec, ecd, ecz, eclr
      real*8 :: rhoa, rhob, z
      real*8 :: ea, eb
      real*8, parameter :: f13 = 0.333333333333333d0
      real*8, parameter :: rsfac = 0.620350490899400d0

C     Spin-unpolarized case

      if (.not.dospin) then

        rs = rsfac/(rhoc**f13)
        call ecPW92(rs,0.0d0,ec,ecd,ecz)
        call ecorrlr(rs,0.0d0,mu,eclr)
        e=(ec-eclr)*rhoc

C     Spin-polarized case
      else

        rs=rsfac/(rhoc**f13)
        rhoa=max((rhoc+rhos)*.5d0,1.0d-15)
        rhob=max((rhoc-rhos)*.5d0,1.0d-15)
        z=(rhoa-rhob)/(rhoa+rhob)
        if (abs(z) > 1.0d0) 
     &   call quit('in ESRC_LDAERFSPIN numerical problem for |z| = 1.0')

        call ecPW92(rs,z,ec,ecd,ecz)
        call ecorrlr(rs,z,mu,eclr)
        e=(ec-eclr)*rhoc

      endif

      return
      end


C*****************************************************************************
      subroutine ELC_LDALAMSPIN(rhoc,rhos,VLAMBDA,dospin,e)
C*****************************************************************************
C linear complement spin-dependent non-scaled correlation LDA functional, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      real*8, intent(in) :: rhoc, rhos, VLAMBDA
      logical, intent(in) :: dospin
      real*8, intent(out) :: e

      call ESRC_LDAERFSPIN(rhoc,rhos,0.d0,dospin,e)
      e = (1.d0 - VLAMBDA**2) * e

      return
       end


C*****************************************************************************
      subroutine ESRC_MD_LDAERF(rhoc,rhos,mu,dospin,e)
C*****************************************************************************
C     Short-range spin-dependent LDA correlation functional with multideterminant reference
C       for OEP calculations from Section V of 
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     Input: rhoc   : total density
C            rhos   : spin density
C            mu     : Interaction parameter
C            dospin : use spin density
C
C     Ouput: e      : energy
C
C     Created: 26-08-11, J. Toulouse
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoc, rhos, mu
      logical, intent(in) :: dospin
      real*8, intent(out) :: e

      real*8 :: e1

      call ESRC_LDAERFSPIN(rhoc,rhos,mu,dospin,e1)
      call DELTA_LRSR_LDAERF(rhoc,rhos,mu,dospin,e)
      e = e1 + e

      end


C*****************************************************************************
      subroutine DELTA_LRSR_LDAERF(rhoc,rhos,mu,dospin,e)
C*****************************************************************************
C     LDA approximation to term Delta_(LR-SR) from  Eq. 42 of
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     Input: rhoc   : total density
C            rhos   : spin density
C            mu     : Interaction parameter
C            dospin : use spin density
C
C     Ouput: e      : energy
C
C     Warning: not tested for z != 0
C
C     Created: 26-08-11, J. Toulouse
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rhoc, rhos, mu
      logical, intent(in) :: dospin
      real*8, intent(out) :: e

      real*8 :: rs, rs2, rs3

      real*8 :: rhoa, rhob, z, z2, onepz, onemz, zp, zm, phi8
      real*8 :: g0f, g0s
      real*8 :: bd2, bd3
      real*8 :: c45, c4, c5
      real*8 :: bc2, bc4, bc3t, bc5t, d0
      real*8 :: delta2, delta3, delta4, delta5, delta6
      real*8 :: delta

      real*8, parameter :: f13 = 0.333333333333333d0
      real*8, parameter :: f83 = 2.6666666666666665d0
      real*8, parameter :: rsfac = 0.620350490899400d0
      real*8, parameter :: alpha2 = 0.2715053589826032d0

      rs = rsfac/(rhoc**f13)
      rs2 = rs*rs
      rs3 = rs2*rs

C     Spin-unpolarized case
      if (.not.dospin) then
       z = 0.d0

C     Spin-polarized case
      else
        rhoa=max((rhoc+rhos)*.5d0,1.0d-15)
        rhob=max((rhoc-rhos)*.5d0,1.0d-15)
        z=(rhoa-rhob)/(rhoa+rhob)
      endif

      z2=z*z
 
      bd2=exp(-0.547d0*rs)*(-0.388d0*rs+0.676*rs2)/rs2
      bd3=exp(-0.31d0*rs)*(-4.95d0*rs+rs2)/rs3

      onepz=1.d0+z
      onemz=1.d0-z
      phi8=0.5d0*(onepz**f83+onemz**f83)

      zp=onepz/2.d0
      zm=onemz/2.d0
      c45=(zp**2)*g0s(rs*zp**(-f13))+(zm**2)*g0s(rs*zm**(-f13))
      c4=c45+(1.d0-z2)*bd2-phi8/(5.d0*alpha2*rs2)
      c5=c45+(1.d0-z2)*bd3
 
      bc2=-3.d0*(1-z2)*(g0f(rs)-0.5d0)/(8.d0*rs3)
      bc4=-9.d0*c4/(64.d0*rs3)
      bc3t=-(1-z2)*g0f(rs)*(2.d0*sqrt(2.d0)-1)/(2.d0*sqrt(pi)*rs3)
      bc5t = -3.d0*c5*(3.d0-sqrt(2.d0))/(20.d0*sqrt(2.d0*pi)*rs3)

      d0=(0.70605d0+0.12927d0*z2)*rs
      delta2=0.073867d0*(rs**(1.5d0))
      delta3=4*(d0**6)*bc3t+(d0**8)*bc5t;
      delta4=4*(d0**6)*bc2+(d0**8)*bc4;
      delta5=(d0**8)*bc3t;
      delta6=(d0**8)*bc2;

      delta=(delta2*(mu**2)+delta3*(mu**3)+delta4*(mu**4)+delta5*(mu**5)
     &      +delta6*(mu**6))/((1+(d0**2)*(mu**2))**4)

C     multiply by rhoc to get energy density
      e=delta*rhoc

      end


C*****************************************************************************
      subroutine ecorrlr(rs,z,mu,eclr)
C*****************************************************************************
C     Hartree atomic units used
C     for given density parameter rs, spin polarization z
C     and cutoff parameter mu
C     gives the correlation energy of the LR gas
C     => eclr
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rs, z, mu
      real*8, intent(out) :: eclr
      real*8 :: ec,ecd,ecz
      real*8 :: phi
      real*8 :: g0f,dpol,d2anti,d3anti,Qrpa
      real*8 :: coe2,coe3,coe4,coe5
      real*8 :: a1,a2,a3,a4,b0
C     parameters from the fit
      real*8, parameter :: adib = 0.784949d0
      real*8, parameter :: q1a  = -0.388d0
      real*8, parameter :: q2a  = 0.676d0
      real*8, parameter :: q3a  = 0.547d0
      real*8, parameter :: t1a  = -4.95d0
      real*8, parameter :: t2a  = 1.d0
      real*8, parameter :: t3a  = 0.31d0
      
      real*8, parameter :: alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      real*8, parameter :: cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0

      b0=adib*rs

      d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
      d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     &        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     &        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     &        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     &        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     &        +(1.-z)**(8.d0/3.d0))/rs**2)

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     &        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     &        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     &        d3anti)
      end if

      call ecPW92(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a3=b0**8*coe3
      a4=b0**6*(b0**2*coe2+4.d0*ec)

      eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+
     &     a4*mu**6+b0**8*mu**8*ec)/((1.d0+b0**2*mu**2)**4)

      return
      end


C*****************************************************************************
      subroutine vcorrlr(rs,z,mu,vclrup,vclrdown)
C*****************************************************************************
C     Hartree atomic units used
C     for given density parameter rs, spin polarization z
C     and cutoff mu it gives the correlation LSD potential for LR interaction
C     => vclrup (spin-up electrons), vclrdown (spin-down electrons)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rs,z,mu
      real*8, intent(out) :: vclrup, vclrdown
      real*8 :: eclr,eclrrs,eclrz
      real*8 :: ec,ecd,ecz
      real*8 :: phi
      real*8 :: g0f,dpol,d2anti,d3anti,Qrpa
      real*8 :: g0d,dpold,d2antid,d3antid,Qrpad,x
      real*8 :: coe2,coe3,coe4,coe5
      real*8 :: coe2rs,coe3rs,coe4rs,coe5rs
      real*8 :: coe2z,coe3z,coe4z,coe5z
      real*8 :: a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs,
     &     b0rs,a1z,a2z,a3z,a4z,a5z,b0z

      real*8, parameter :: alpha = (4.d0/9.d0/pi)**(1.d0/3.d0)
      real*8, parameter :: cf = 1.d0/alpha

C     parameters from the fit
      real*8, parameter :: adib = 0.784949d0
      real*8, parameter :: q1a = -0.388d0
      real*8, parameter :: q2a = 0.676d0
      real*8, parameter :: q3a = 0.547d0
      real*8, parameter :: t1a = -4.95d0
      real*8, parameter :: t2a = 1.d0
      real*8, parameter :: t3a = 0.31d0

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
      b0=adib*rs

      d2anti=(q1a+q2a*rs)*exp(-q3a*rs)/rs
      d3anti=(t1a+t2a*rs)*exp(-t3a*rs)/rs**2

      d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/
     &    rs**2)*exp(-q3a*rs)
      d3antid=-((rs*t2a*(1 + rs*t3a) + t1a*(2 + rs*t3a))/
     &    rs**3)*exp(-rs*t3a)

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)
      coe2rs=-3.d0/8.d0/rs**3*(1.d0-z**2)*g0d(rs)+
     &     9.d0/8.d0/rs**4*(1.d0-z**2)*(g0f(rs)-0.5d0)
      coe2z=-3.d0/8.d0/rs**3*(-2.d0*z)*(g0f(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**3)
      coe3rs=-(1.d0-z**2)*g0d(rs)/(sqrt(2.d0*pi)*rs**3)+
     &    3.d0*(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**4)
      coe3z=2.d0*z*g0f(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     &        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(dpold(rs)
     &        +2.d0*cf**2*2**(5.d0/3.d0)/5.d0/rs**3)
        coe4z=-9.d0/64.d0/rs**3*(dpol(rs)-rs/6.d0*dpold(rs)-2.d0*d2anti
     &       -4.d0/15.d0/rs**2*cf**2*2.d0**(5.d0/3.d0))*z
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)
        coe5rs=-3.d0/rs*coe5-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpold(rs)
        coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(dpol(rs)-rs/6.d0*
     &       dpold(rs)-2.d0*d3anti)*z

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     &        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     &        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     &        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     &        +(1.-z)**(8.d0/3.d0))/rs**2)
         coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(
     &        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     &        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     &        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     &        d2antid+cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)
     &        +(1.d0-z)**(8.d0/3.d0))/rs**3)
         coe4z=-9.d0/64.d0/rs**3*(1.d0/2.d0*(1.d0+z)*
     &        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     &        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     &        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     &        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     &        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d2anti-
     &        4.d0/15.d0/rs**2*cf**2*((1.d0+z)**(5.d0/3.d0)-
     &        (1.d0-z)**(5.d0/3.d0)))

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     &        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     &        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     &        d3anti)
         coe5rs=-3.d0/rs*coe5-9.d0/(40.d0*sqrt(2.d0*pi)*rs**3)*(
     &        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     &        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     &        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     &        d3antid)
         coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(1.d0/2.d0*(1.d0+z)*
     &        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     &        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     &        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     &        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     &        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d3anti)

      end if

      call ecPW92(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a1rs=24.d0*adib*b0**5*coe3+4.d0*b0**6*coe3rs+8.d0*adib*b0**7*
     &     coe5+b0**8*coe5rs
      a1z=4.d0*b0**6*coe3z+b0**8*coe5z

      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a2rs=24.d0*adib*b0**5*coe2+4.d0*b0**6*coe2rs+8.d0*adib*b0**7*
     &     coe4+b0**8*coe4rs+24.d0*adib*b0**3*ec+6.d0*b0**4*ecd
      a2z=4.d0*b0**6*coe2z+b0**8*coe4z+6.d0*b0**4*ecz

      a3=b0**8*coe3
      a3rs=8.d0*adib*b0**7*coe3+b0**8*coe3rs
      a3z=b0**8*coe3z

      a4=b0**6*(b0**2*coe2+4.d0*ec)
      a4rs=8.d0*adib*b0**7*coe2+b0**8*coe2rs+24.d0*adib*b0**5*ec+
     &     4.d0*b0**6*ecd
      a4z=b0**6*(b0**2*coe2z+4.d0*ecz)

      a5=b0**8*ec
      a5rs=8.d0*adib*b0**7*ec+b0**8*ecd
      a5z=b0**8*ecz

      x=mu*sqrt(rs)/phi

      eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+
     &     a4*mu**6+a5*mu**8)/((1.d0+b0**2*mu**2)**4)

      eclrrs=-4.d0/(1.d0+b0**2*mu**2)*2.d0*adib*b0*mu**2*eclr+
     &     1.d0/((1.d0+b0**2*mu**2)**4)*(phi**2*mu/(2.d0*sqrt(rs))
     &     *Qrpad(x)+
     &     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)


      if(z.eq.1.d0) then
         vclrup=eclr-rs/3.d0*eclrrs
         vclrdown=0.d0
      elseif(z.eq.-1.d0) then
         vclrup=0.d0
         vclrdown=eclr-rs/3.d0*eclrrs
      else

         eclrz=(phi**2*((1.d0+z)**(-1.d0/3.d0)-(1.d0-z)**(-1.d0/3.d0))
     &        *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1.d0+z)**(-1.d0/3.d0)
     &        -(1.d0-z)**(-1.d0/3.d0))/3.d0+
     &        a1z*mu**3+a2z*mu**4+a3z*mu**5+
     &        a4z*mu**6+a5z*mu**8)/((1.d0+b0**2*mu**2)**4)

         vclrup=eclr-rs/3.d0*eclrrs-(z-1.d0)*eclrz
         vclrdown=eclr-rs/3.d0*eclrrs-(z+1.d0)*eclrz
      end if
      return
      end


C*****************************************************************************
      real*8 pure function g0f(x)
C*****************************************************************************
C     on-top pair-distribution function
C     Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
C     x -> rs
C*****************************************************************************
      implicit none
      real*8, intent(in) :: x
      real*8, parameter :: C0f = 0.0819306d0
      real*8, parameter :: D0f = 0.752411d0
      real*8, parameter :: E0f = -0.0127713d0
      real*8, parameter :: F0f = 0.00185898d0
      g0f=(1.d0-(0.7317d0-D0f)*x+C0f*x**2+E0f*x**3+
     &     F0f*x**4)*exp(-abs(D0f)*x)/2.d0
      return
      end


C*****************************************************************************
      real*8 pure function g0d(rs)
C*****************************************************************************
C     derivative of on-top pair-distribution function
C     Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rs
      real*8, parameter :: Cg0 = 0.0819306d0
      real*8, parameter :: Fg0 = 0.752411d0
      real*8, parameter :: Dg0 = -0.0127713d0
      real*8, parameter :: Eg0 = 0.00185898d0
      real*8, parameter :: Bg0 =0.7317d0-Fg0
      g0d=(-Bg0+2*Cg0*rs+3*Dg0*rs**2+4*Eg0*rs**3)/2.d0*exp(-Fg0*rs)
     &   - (Fg0*(1 - Bg0*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4))/
     &   2.d0*exp(-Fg0*rs)
      return
      end


C*****************************************************************************
      real*8 pure function dpol(rs)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rs
      real*8, parameter :: cf = (9.d0*pi/4.d0)**(1.d0/3.d0)
      real*8, parameter :: p2p = 0.04d0
      real*8, parameter :: p3p = 0.4319d0
      dpol=2.d0**(5.d0/3.d0)/5.d0*cf**2/rs**2*(1.d0+(p3p-0.454555d0)*rs)
     &     /(1.d0+p3p*rs+p2p*rs**2)
      return
      end


C*****************************************************************************
      real*8 pure function dpold(rs)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: rs
      real*8, parameter :: cf = (9.d0*pi/4.d0)**(1.d0/3.d0)
      real*8, parameter :: p2p = 0.04d0
      real*8, parameter :: p3p = 0.4319d0
      dpold=2.d0**(5.d0/3.d0)/5.d0*cf**2*
     & (-2. + (0.454555 - 4.*p3p)*rs +
     &    (-4.*p2p +
     &       (0.90911 - 2.*p3p)*p3p)*rs**2
     &      + p2p*(1.363665 - 3.*p3p)*
     &     rs**3)/
     &  (rs**3*(1. + p3p*rs + p2p*rs**2)**2)
      return
      end


C*****************************************************************************
      real*8 pure function Qrpa(x)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: x
      real*8, parameter :: Acoul = 2.d0*(log(2.d0)-1.d0)/pi2
      real*8, parameter :: a2 = 5.84605d0
      real*8, parameter :: c2 = 3.91744d0
      real*8, parameter :: d2 = 3.44851d0
      real*8, parameter :: b2 = d2-3.d0/(2.d0*pi*Acoul)*
     &                           (4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpa=Acoul*log((1.d0+a2*x+b2*x**2+c2*x**3)/(1.d0+a2*x+d2*x**2))
      return
      end


C*****************************************************************************
      real*8 pure function Qrpad(x)
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: x
      real*8, parameter :: Acoul = 2.d0*(log(2.d0)-1.d0)/pi2
      real*8, parameter :: a2 = 5.84605d0
      real*8, parameter :: c2 = 3.91744d0
      real*8, parameter :: d2 = 3.44851d0
      real*8, parameter :: b2 = d2-3.d0/(2.d0*pi*Acoul) * 
     &                           (4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpad=Acoul*((x*(b2*(2.d0 + a2*x) +
     &      c2*x*(3.d0 + 2.d0*a2*x) +
     &      d2*(-2.d0 - a2*x + c2*x**3)))/
     &  ((1.d0 + a2*x + d2*x**2)*
     &    (1.d0 + a2*x + b2*x**2 + c2*x**3)))
      return
      end


C*****************************************************************************
      subroutine ecPW92(x,y,ec,ecd,ecz)
C*****************************************************************************
C     correlation energy and its derivative w.r.t. rs and z at mu=infinity
C     Perdew & Wang PRB 45, 13244 (1992)
C     in Hartree; ec=ec(rs,zeta)
C     x -> rs; y -> zeta
C     ecd is d/drs ec
C     ecz is d/dz ec
C*****************************************************************************
      implicit none
#include "pi.h"
      real*8, intent(in) :: x, y
      real*8, intent(out) :: ec, ecd, ecz
      real*8 :: f02,ff,ec0,ec0d,ec1,ec1d,
     &          aaa,G,Gd,alfac,alfacd

      f02=4.d0/(9.d0*(2.d0**(1.d0/3.d0)-1.d0))

      ff=((1.d0+y)**(4.d0/3.d0)+(1.d0-y)**(4.d0/3.d0)-
     &     2.d0)/(2.d0**(4.d0/3.d0)-2.d0)

      aaa=(1.d0-log(2.d0))/pi**2
      call GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0,
     &     1.6382d0,0.49294d0,G,Gd)
      ec0 =G
      ec0d=Gd

      aaa=aaa*0.5d0
      call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0,
     &     3.3662d0,0.62517d0,G,Gd)
      ec1 =G
      ec1d=Gd
      call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,
     &     0.88026d0,0.49671d0,G,Gd)
      alfac =-G
      alfacd=-Gd

      ec =ec0 +alfac*ff/f02*(1.d0-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1.d0-y**4)+(ec1d-ec0d)*
     &     ff*y**4
      ecz=alfac*(-4.d0*y**3)*ff/f02+alfac*(1.d0-y**4)/f02*
     &     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     &     (2.d0**(4.d0/3.d0)-2.d0)+(ec1-ec0)*(4.d0*y**3*ff+
     &     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     &     (2.d0**(4.d0/3.d0)-2.d0)*y**4)

      return
      end


C*****************************************************************************
      pure subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd)
C*****************************************************************************
C     Gd is d/drs G
C*****************************************************************************
      implicit none
      real*8, intent(in) :: Ac,alfa1,beta1,beta2,beta3,beta4,x
      real*8, intent(out) :: G, Gd
      G=-2.d0*Ac*(1.d0+alfa1*x)*log(1.d0+1.d0/(2.d0*
     &     Ac*(beta1*x**0.5d0+
     &     beta2*x+beta3*x**1.5d0+beta4*x**2)))
      Gd=(1.d0+alfa1*x)*(beta2+beta1/(2.d0*sqrt(x))+3.d0*beta3*
     &     sqrt(x)/2.d0+2.d0*beta4*x)/((beta1*sqrt(x)+beta2*x+
     &     beta3*x**(3.d0/2.d0)+beta4*x**2)**2*(1.d0+1.d0/
     &     (2.d0*Ac*(beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     &     beta4*x**2))))-2.d0*Ac*alfa1*log(1.d0+1.d0/(2.d0*Ac*
     &     (beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     &     beta4*x**2)))
      return
      end


C*****************************************************************************
      real*8 pure function g0s(rs)
C*****************************************************************************
C     g"(0,rs,z=1) from Eq. 32 of
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     Created: 26-08-11, J. Toulouse
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rs
      real*8, parameter :: f53 = 1.6666666666666667d0
      real*8, parameter :: alpha2 = 0.2715053589826032d0
      real*8 :: rs2
      rs2 = rs*rs
      g0s=(2.d0**f53)*(1.d0-0.02267d0*rs)/((5.d0*alpha2*rs2)
     &    *(1.d0+0.4319d0*rs+0.04d0*rs2))
      end


C*****************************************************************************
      pure subroutine HS(z,k,H)
C*****************************************************************************
C     Smoothening formula for the SRCMULO and SRCPBELO functionals 
C     Cancels the discontinuety arising from the max() function.
C
C     INPUT:   
C             z   :   max(a,b)
C             k   :   
C     OUTPUT:
C             H  
C*****************************************************************************
      real*8, intent(in) :: z, k
      real*8, intent(out) :: H
! dftcom.h : heaviside_pvalue
#include "dftcom.h"

      H = 1/(1+k*exp(-heaviside_pvalue*z))

      return
      end
C ... end of dftfunjt.F
