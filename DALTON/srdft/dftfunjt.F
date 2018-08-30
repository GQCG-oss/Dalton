C File: srdft/dftfunjt.F


! hjaaj April 2018: search for TODO and fix them all !!!

C*****************************************************************************
      subroutine DESRC_LDA(ESRC_LDA, rho, mu, norder, E, d1E, d2E)
C*****************************************************************************
C   Finite difference code for spin-unpolarized LDA correlation functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(1) : total density
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRC_LDA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRC_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: August 30, 2018
C*****************************************************************************
      implicit none
      external ESRC_LDA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(2)
      real*8, intent(out) :: E, d1E(5), d2E(15)
      
      real*8 :: rho_c, rho_cp, rho_cm, Ep_1, Em_1, hr_c
      ! dummy_rho_s is introduced to make code
      ! compatible with existing code
      real*8, parameter :: dummy_rho_s = 0.0d0
      
      rho_c = rho(1)
      
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
      
      ! Variation of variables
      rho_cp = rho_c + hr_c
      rho_cm = rho_c - hr_c
      
      call ESRC_LDA(rho_c,dummy_rho_s,mu,E)
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_c .gt. 0.0d0) then
            call ESRC_LDA(rho_cp,dummy_rho_s,mu,Ep_1)
            call ESRC_LDA(rho_cm,dummy_rho_s,mu,Em_1)
         end if
      end if
      ! Calculate first derivatives
      if (norder .ge. 1) then
         ! d1E / drho_c
         if (hr_c .gt. 0.0d0) then
            d1E(1)=(Ep_1-Em_1)/(2.0d0*hr_c)
         end if
      end if
      ! Calculate second derivatives
      if (norder .ge. 2) then
         if (hr_c .gt. 0.0d0) then
            ! d2E / (drho_c**2)
            d2E(1)=(Ep_1-2.0d0*E+Em_1)/(hr_c**2)
         end if
      end if
      
      end subroutine

C*****************************************************************************
      subroutine DESRC_SPINLDA(ESRC_SPINLDA, rho, mu, norder, E, d1E,
     & d2E)
C*****************************************************************************
C   Finite difference code for spin-polarized LDA correlation functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(1) : total density
C            rho(2) : spin density
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRC_SPINLDA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRC_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: August 30, 2018
C*****************************************************************************
      implicit none
      external ESRC_SPINLDA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(2)
      real*8, intent(out) :: E, d1E(5), d2E(15)
      
      real*8 :: rho_c, rho_cp, rho_cm, rho_s, rho_sp, rho_sm, Ep_1,
     & Em_1, Epp_2_1, Emm_2_1, Ep_2, Em_2, hr_c, hr_s
      
      rho_c = rho(1)
      rho_s = rho(2)
      ! Do checks of input
      if (dabs(rho_s) .gt. rho_c) then
         call quit('DESRC_SPINLDA ERROR: |rho_s| > rho_c')
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
      
      ! Variation of variables
      rho_cp = rho_c + hr_c
      rho_cm = rho_c - hr_c
      rho_sp = rho_s + hr_s
      rho_sm = rho_s - hr_s
      
      call ESRC_SPINLDA(rho_c,rho_s,mu,E)
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_c .gt. 0.0d0) then
            call ESRC_SPINLDA(rho_cp,rho_s,mu,Ep_1)
            call ESRC_SPINLDA(rho_cm,rho_s,mu,Em_1)
         end if
         if (hr_s .gt. 0.0d0) then
            call ESRC_SPINLDA(rho_c,rho_sp,mu,Ep_2)
            call ESRC_SPINLDA(rho_c,rho_sm,mu,Em_2)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hr_s .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINLDA(rho_cp,rho_sp,mu,Epp_2_1)
               call ESRC_SPINLDA(rho_cm,rho_sm,mu,Emm_2_1)
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
               d2E(2)=(Epp_2_1-Ep_2-Ep_1+2.0d0*E-Em_2-Em_1+Emm_2_1)/
     &                     (2.0d0*hr_s*hr_c)
            end if
            ! d2E / (drho_s**2)
            d2E(3)=(Ep_2-2.0d0*E+Em_2)/(hr_s**2)
         end if
      end if
      
      end subroutine

C*****************************************************************************
      subroutine DESRC_GGA(ESRC_GGA, rho, grad2, mu, norder, E, d1E,
     & d2E)
C*****************************************************************************
C   Finite difference code for spin-unpolarized GGA correlation functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(1) : total density
C            grad2(1) : grad_c . grad_c
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRC_GGA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRC_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: August 30, 2018
C*****************************************************************************
      implicit none
      external ESRC_GGA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(2), grad2(3)
      real*8, intent(out) :: E, d1E(5), d2E(15)
      
      real*8 :: rho_c, rho_cp, rho_cm, gamma_cc, gamma_ccp, gamma_ccm,
     & Ep_1, Em_1, Epp_3_1, Emm_3_1, Ep_3, Em_3, hr_c, hg_cc
      
      rho_c = rho(1)
      gamma_cc = grad2(1)
      ! Do checks of input
      if (gamma_cc .lt. 0.0d0) then
         call quit('DESRC_GGA ERROR: gamma_cc negative')
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
      
      ! Variation of variables
      rho_cp = rho_c + hr_c
      rho_cm = rho_c - hr_c
      gamma_ccp = gamma_cc + hg_cc
      gamma_ccm = gamma_cc - hg_cc
      
      call ESRC_GGA(rho_c,gamma_cc,mu,E)
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_c .gt. 0.0d0) then
            call ESRC_GGA(rho_cp,gamma_cc,mu,Ep_1)
            call ESRC_GGA(rho_cm,gamma_cc,mu,Em_1)
         end if
         if (hg_cc .gt. 0.0d0) then
            call ESRC_GGA(rho_c,gamma_ccp,mu,Ep_3)
            call ESRC_GGA(rho_c,gamma_ccm,mu,Em_3)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hg_cc .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_GGA(rho_cp,gamma_ccp,mu,Epp_3_1)
               call ESRC_GGA(rho_cm,gamma_ccm,mu,Emm_3_1)
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
               d2E(4)=(Epp_3_1-Ep_3-Ep_1+2.0d0*E-Em_3-Em_1+Emm_3_1)/
     &                     (2.0d0*hg_cc*hr_c)
            end if
            ! d2E / (dgamma_cc**2)
            d2E(6)=(Ep_3-2.0d0*E+Em_3)/(hg_cc**2)
         end if
      end if
      
      end subroutine

C*****************************************************************************
      subroutine DESRC_SPINGGA(ESRC_SPINGGA, rho, grad2, mu, norder,
     & E, d1E, d2E)
C*****************************************************************************
C   Finite difference code for spin-polarized GGA correlation functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(1) : total density
C            rho(2) : spin density
C            grad2(1) : grad_c . grad_c
C            grad2(2) : grad_s . grad_s
C            grad2(3) : grad_c . grad_s
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRC_SPINGGA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRC_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: August 30, 2018
C*****************************************************************************
      implicit none
      external ESRC_SPINGGA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(2), grad2(3)
      real*8, intent(out) :: E, d1E(5), d2E(15)
      
      real*8 :: rho_c, rho_cp, rho_cm, rho_s, rho_sp, rho_sm, gamma_cc,
     & gamma_ccp, gamma_ccm, gamma_ss, gamma_ssp, gamma_ssm, gamma_cs,
     & gamma_csp, gamma_csm, Ep_1, Em_1, Epp_2_1, Emm_2_1, Ep_2, Em_2,
     & Epp_3_1, Emm_3_1, Epp_3_2, Emm_3_2, Ep_3, Em_3, Epp_4_1, Emm_4_1,
     & Epp_4_2, Emm_4_2, Epp_4_3, Emm_4_3, Ep_4, Em_4, Epp_5_1, Emm_5_1,
     & Epp_5_2, Emm_5_2, Epp_5_3, Emm_5_3, Epp_5_4, Emm_5_4, Ep_5,
     & Em_5, hr_c, hr_s, hg_cc, hg_ss, hg_cs
      
      rho_c = rho(1)
      rho_s = rho(2)
      gamma_cc = grad2(1)
      gamma_ss = grad2(2)
      gamma_cs = grad2(3)
      ! Do checks of input
      if (dabs(rho_s) .gt. rho_c) then
         call quit('DESRC_SPINGGA ERROR: |rho_s| > rho_c')
      end if
      if (gamma_cc .lt. 0.0d0) then
         call quit('DESRC_SPINGGA ERROR: gamma_cc negative')
      end if
      if (gamma_ss .lt. 0.0d0) then
         call quit('DESRC_SPINGGA ERROR: gamma_ss negative')
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
      
      call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_cs,mu,
     &E)
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_c .gt. 0.0d0) then
            call ESRC_SPINGGA(rho_cp,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              mu,Ep_1)
            call ESRC_SPINGGA(rho_cm,rho_s,gamma_cc,gamma_ss,gamma_cs,
     &              mu,Em_1)
         end if
         if (hr_s .gt. 0.0d0) then
            call ESRC_SPINGGA(rho_c,rho_sp,gamma_cc,gamma_ss,gamma_cs,
     &              mu,Ep_2)
            call ESRC_SPINGGA(rho_c,rho_sm,gamma_cc,gamma_ss,gamma_cs,
     &              mu,Em_2)
         end if
         if (hg_cc .gt. 0.0d0) then
            call ESRC_SPINGGA(rho_c,rho_s,gamma_ccp,gamma_ss,gamma_cs,
     &              mu,Ep_3)
            call ESRC_SPINGGA(rho_c,rho_s,gamma_ccm,gamma_ss,gamma_cs,
     &              mu,Em_3)
         end if
         if (hg_ss .gt. 0.0d0) then
            call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ssp,gamma_cs,
     &              mu,Ep_4)
            call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ssm,gamma_cs,
     &              mu,Em_4)
         end if
         if (hg_cs .gt. 0.0d0) then
            call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_csp,
     &              mu,Ep_5)
            call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ss,gamma_csm,
     &              mu,Em_5)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hr_s .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_cp,rho_sp,gamma_cc,gamma_ss,
     &                     gamma_cs,mu,Epp_2_1)
               call ESRC_SPINGGA(rho_cm,rho_sm,gamma_cc,gamma_ss,
     &                     gamma_cs,mu,Emm_2_1)
            end if
         end if
         if (hg_cc .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_cp,rho_s,gamma_ccp,gamma_ss,
     &                     gamma_cs,mu,Epp_3_1)
               call ESRC_SPINGGA(rho_cm,rho_s,gamma_ccm,gamma_ss,
     &                     gamma_cs,mu,Emm_3_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_c,rho_sp,gamma_ccp,gamma_ss,
     &                     gamma_cs,mu,Epp_3_2)
               call ESRC_SPINGGA(rho_c,rho_sm,gamma_ccm,gamma_ss,
     &                     gamma_cs,mu,Emm_3_2)
            end if
         end if
         if (hg_ss .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_cp,rho_s,gamma_cc,gamma_ssp,
     &                     gamma_cs,mu,Epp_4_1)
               call ESRC_SPINGGA(rho_cm,rho_s,gamma_cc,gamma_ssm,
     &                     gamma_cs,mu,Emm_4_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_c,rho_sp,gamma_cc,gamma_ssp,
     &                     gamma_cs,mu,Epp_4_2)
               call ESRC_SPINGGA(rho_c,rho_sm,gamma_cc,gamma_ssm,
     &                     gamma_cs,mu,Emm_4_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_c,rho_s,gamma_ccp,gamma_ssp,
     &                     gamma_cs,mu,Epp_4_3)
               call ESRC_SPINGGA(rho_c,rho_s,gamma_ccm,gamma_ssm,
     &                     gamma_cs,mu,Emm_4_3)
            end if
         end if
         if (hg_cs .gt. 0.0d0) then
            if (hr_c .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_cp,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         p,mu,Epp_5_1)
               call ESRC_SPINGGA(rho_cm,rho_s,gamma_cc,gamma_ss,gamma_cs
     &         m,mu,Emm_5_1)
            end if
            if (hr_s .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_c,rho_sp,gamma_cc,gamma_ss,gamma_cs
     &         p,mu,Epp_5_2)
               call ESRC_SPINGGA(rho_c,rho_sm,gamma_cc,gamma_ss,gamma_cs
     &         m,mu,Emm_5_2)
            end if
            if (hg_cc .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_c,rho_s,gamma_ccp,gamma_ss,gamma_cs
     &         p,mu,Epp_5_3)
               call ESRC_SPINGGA(rho_c,rho_s,gamma_ccm,gamma_ss,gamma_cs
     &         m,mu,Emm_5_3)
            end if
            if (hg_ss .gt. 0.0d0) then
               call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ssp,gamma_cs
     &         p,mu,Epp_5_4)
               call ESRC_SPINGGA(rho_c,rho_s,gamma_cc,gamma_ssm,gamma_cs
     &         m,mu,Emm_5_4)
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
               d2E(2)=(Epp_2_1-Ep_2-Ep_1+2.0d0*E-Em_2-Em_1+Emm_2_1)/
     &                     (2.0d0*hr_s*hr_c)
            end if
            ! d2E / (drho_s**2)
            d2E(3)=(Ep_2-2.0d0*E+Em_2)/(hr_s**2)
         end if
         if (hg_cc .gt. 0.0d0) then
            ! d2E / (dgamma_cc * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(4)=(Epp_3_1-Ep_3-Ep_1+2.0d0*E-Em_3-Em_1+Emm_3_1)/
     &                     (2.0d0*hg_cc*hr_c)
            end if
            ! d2E / (dgamma_cc * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(5)=(Epp_3_2-Ep_3-Ep_2+2.0d0*E-Em_3-Em_2+Emm_3_2)/
     &                     (2.0d0*hg_cc*hr_s)
            end if
            ! d2E / (dgamma_cc**2)
            d2E(6)=(Ep_3-2.0d0*E+Em_3)/(hg_cc**2)
         end if
         if (hg_ss .gt. 0.0d0) then
            ! d2E / (dgamma_ss * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(7)=(Epp_4_1-Ep_4-Ep_1+2.0d0*E-Em_4-Em_1+Emm_4_1)/
     &                     (2.0d0*hg_ss*hr_c)
            end if
            ! d2E / (dgamma_ss * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(8)=(Epp_4_2-Ep_4-Ep_2+2.0d0*E-Em_4-Em_2+Emm_4_2)/
     &                     (2.0d0*hg_ss*hr_s)
            end if
            ! d2E / (dgamma_ss * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(9)=(Epp_4_3-Ep_4-Ep_3+2.0d0*E-Em_4-Em_3+Emm_4_3)/
     &                     (2.0d0*hg_ss*hg_cc)
            end if
            ! d2E / (dgamma_ss**2)
            d2E(10)=(Ep_4-2.0d0*E+Em_4)/(hg_ss**2)
         end if
         if (hg_cs .gt. 0.0d0) then
            ! d2E / (dgamma_cs * drho_c)
            if (hr_c .gt. 0.0d0) then
               d2E(11)=(Epp_5_1-Ep_5-Ep_1+2.0d0*E-Em_5-Em_1+Emm_5_1)/
     &                     (2.0d0*hg_cs*hr_c)
            end if
            ! d2E / (dgamma_cs * drho_s)
            if (hr_s .gt. 0.0d0) then
               d2E(12)=(Epp_5_2-Ep_5-Ep_2+2.0d0*E-Em_5-Em_2+Emm_5_2)/
     &                     (2.0d0*hg_cs*hr_s)
            end if
            ! d2E / (dgamma_cs * dgamma_cc)
            if (hg_cc .gt. 0.0d0) then
               d2E(13)=(Epp_5_3-Ep_5-Ep_3+2.0d0*E-Em_5-Em_3+Emm_5_3)/
     &                     (2.0d0*hg_cs*hg_cc)
            end if
            ! d2E / (dgamma_cs * dgamma_ss)
            if (hg_ss .gt. 0.0d0) then
               d2E(14)=(Epp_5_4-Ep_5-Ep_4+2.0d0*E-Em_5-Em_4+Emm_5_4)/
     &                     (2.0d0*hg_cs*hg_ss)
            end if
            ! d2E / (dgamma_cs**2)
            d2E(15)=(Ep_5-2.0d0*E+Em_5)/(hg_cs**2)
         end if
      end if
      
      end subroutine

C*****************************************************************************
      subroutine DESRX_BOTHLDA(ESRX_BOTHLDA, rho, mu, norder, E, d1E,
     & d2E)
C*****************************************************************************
C   Finite difference code for spin-polarized and spin-unpolarized LDA exchange functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(3) : alpha density
C            rho(4) : beta density
C            mu     : range seperation parameter
C            norder : order of derivative (up to two)
C            ESRX_BOTHLDA : Functional routine
C   Output : E      : energy
C            d1E    : first derivative of energy
C            d2E    : second derivative of energy
C   See notes in SRDFTLTRJT() for ordering of derivatives.
C   
C   Written by, Erik R. Kjellgren
C   Based DESRX_SPINGGA by H. J. Aa. Jensen and Erik D. Hedegaard
C   Generated: August 30, 2018
C*****************************************************************************
      implicit none
      external ESRX_BOTHLDA
#include "priunit.h"
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(4)
      real*8, intent(out) :: E, d1E(5), d2E(15)
      
      logical :: alpha_not_equal_beta
      real*8 :: rho_a, rho_ap, rho_am, rho_b, rho_bp, rho_bm, Ep_1,
     & Em_1, Epp_2_1, Emm_2_1, Ep_2, Em_2, hr_a, hr_b, Ea, Eb, d1Eab(5),
     & d2Eab(15)
      integer, save :: call_counter = 0
      
      call_counter = call_counter +1
      if (call_counter .eq. 1) then
         write(lupri,*)'DESRX_BOTHLDA: WARNING!'
         write(lupri,*)'DESRX_BOTHLDA: routine not covered by ctest'
      end if
      
      rho_a = rho(3)
      rho_b = rho(4)
      
      ! Check if alpha is different from beta
      ! If alpha equal beta, then alot of derivatives are the same
      if (dabs(rho_a - rho_b) .le. 1.d-14) then
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
      ! Variation of variables
      rho_ap = rho_a + hr_a
      rho_am = rho_a - hr_a
      
      if (alpha_not_equal_beta) then
         Eb = 0.0d0
         if (rho_b .lt. 1.d-10) then
            hr_b = 0.0d0
         else
            hr_b = 1.d-4*rho_b
         end if
         rho_bp = rho_b + hr_b
         rho_bm = rho_b - hr_b
      end if
      
      call ESRX_BOTHLDA(rho_a,mu,Ea)
      if (alpha_not_equal_beta) then
         call ESRX_BOTHLDA(rho_b,mu,Eb)
         E = Ea + Eb
      else
         E = 2.0d0*Ea
      end if
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_a .gt. 0.0d0) then
            call ESRX_BOTHLDA(rho_ap,mu,Ep_1)
            call ESRX_BOTHLDA(rho_am,mu,Em_1)
         end if
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHLDA(rho_bp,mu,Ep_2)
            call ESRX_BOTHLDA(rho_bm,mu,Em_2)
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
      end if
      
      ! If alpha is equal to beta, assing derivatives
      ! that was not calculated above
      if (.not. alpha_not_equal_beta) then
         if (norder .ge. 1) then
            d1Eab(2)=d1Eab(1)
         end if
         if (norder .ge. 2) then
            d2Eab(3)=d2Eab(1)
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
      end if
      
      if (norder .ge. 2) then
         d2E(1) = 0.25d0*(d2Eab(1) + d2Eab(3))
         d2E(2) = 0.25d0*(d2Eab(1) - d2Eab(3))
         d2E(3) = d2E(1)
      end if
      
      end subroutine

C*****************************************************************************
      subroutine DESRX_BOTHGGA(ESRX_BOTHGGA, rho, grad2, mu, norder,
     & E, d1E, d2E)
C*****************************************************************************
C   Finite difference code for spin-polarized and spin-unpolarized GGA exchange functiontionals.
C   Return up to second order derivative.
C   
C   Input  : rho(3) : alpha density
C            rho(4) : beta density
C            grad2(4) : grad_alpha . grad_alpha
C            grad2(5) : grad_beta . grad_beta
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
C   Generated: August 30, 2018
C*****************************************************************************
      implicit none
      external ESRX_BOTHGGA
      
      integer, intent(in) :: norder
      real*8, intent(in) :: mu, rho(4), grad2(6)
      real*8, intent(out) :: E, d1E(5), d2E(15)
      
      logical :: alpha_not_equal_beta
      real*8 :: rho_a, rho_ap, rho_am, rho_b, rho_bp, rho_bm, gamma_aa,
     & gamma_aap, gamma_aam, gamma_bb, gamma_bbp, gamma_bbm, Ep_1,
     & Em_1, Epp_2_1, Emm_2_1, Ep_2, Em_2, Epp_3_1, Emm_3_1, Epp_3_2,
     & Emm_3_2, Ep_3, Em_3, Epp_4_1, Emm_4_1, Epp_4_2, Emm_4_2, Epp_4_3,
     & Emm_4_3, Ep_4, Em_4, hr_a, hr_b, hg_aa, hg_bb, hg_ab, Ea, Eb,
     & d1Eab(5), d2Eab(15)
      
      rho_a = rho(3)
      rho_b = rho(4)
      gamma_aa = grad2(4)
      gamma_bb = grad2(5)
      
      ! Check if alpha is different from beta
      ! If alpha equal beta, then alot of derivatives are the same
      if (dabs(gamma_aa - gamma_bb) .le. 1.d-14 .and. dabs(rho_a - rho_b
     &   ) .le. 1.d-14) then
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
      ! Variation of variables
      rho_ap = rho_a + hr_a
      rho_am = rho_a - hr_a
      gamma_aap = gamma_aa + hg_aa
      gamma_aam = gamma_aa - hg_aa
      
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
         rho_bp = rho_b + hr_b
         rho_bm = rho_b - hr_b
         gamma_bbp = gamma_bb + hg_bb
         gamma_bbm = gamma_bb - hg_bb
      end if
      
      call ESRX_BOTHGGA(rho_a,gamma_aa,mu,Ea)
      if (alpha_not_equal_beta) then
         call ESRX_BOTHGGA(rho_b,gamma_bb,mu,Eb)
         E = Ea + Eb
      else
         E = 2.0d0*Ea
      end if
      ! Calculate energies for first derivatives
      if (norder .ge. 1) then
         if (hr_a .gt. 0.0d0) then
            call ESRX_BOTHGGA(rho_ap,gamma_aa,mu,Ep_1)
            call ESRX_BOTHGGA(rho_am,gamma_aa,mu,Em_1)
         end if
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHGGA(rho_bp,gamma_bb,mu,Ep_2)
            call ESRX_BOTHGGA(rho_bm,gamma_bb,mu,Em_2)
         end if
         if (hg_aa .gt. 0.0d0) then
            call ESRX_BOTHGGA(rho_a,gamma_aap,mu,Ep_3)
            call ESRX_BOTHGGA(rho_a,gamma_aam,mu,Em_3)
         end if
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            call ESRX_BOTHGGA(rho_b,gamma_bbp,mu,Ep_4)
            call ESRX_BOTHGGA(rho_b,gamma_bbm,mu,Em_4)
         end if
      end if
      ! Calculate energies for second derivatives
      if (norder .ge. 2) then
         if (hr_b .gt. 0.0d0 .and. alpha_not_equal_beta) then
         end if
         if (hg_aa .gt. 0.0d0) then
            if (hr_a .gt. 0.0d0) then
               call ESRX_BOTHGGA(rho_ap,gamma_aap,mu,Epp_3_1)
               call ESRX_BOTHGGA(rho_am,gamma_aam,mu,Emm_3_1)
            end if
         end if
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            if (hr_b .gt. 0.0d0) then
               call ESRX_BOTHGGA(rho_bp,gamma_bbp,mu,Epp_4_2)
               call ESRX_BOTHGGA(rho_bm,gamma_bbm,mu,Emm_4_2)
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
               d2Eab(4)=(Epp_3_1-Ep_3-Ep_1+2.0d0*Ea-Em_3-Em_1+Emm_3_1)/
     &                     (2.0d0*hg_aa*hr_a)
            end if
            ! d2E / (dgamma_aa**2)
            d2Eab(6)=(Ep_3-2.0d0*Ea+Em_3)/(hg_aa**2)
         end if
         if (hg_bb .gt. 0.0d0 .and. alpha_not_equal_beta) then
            ! d2E / (dgamma_bb * drho_b)
            if (hr_b .gt. 0.0d0) then
               d2Eab(8)=(Epp_4_2-Ep_4-Ep_2+2.0d0*Eb-Em_4-Em_2+Emm_4_2)/
     &                     (2.0d0*hg_bb*hr_b)
            end if
            ! d2E / (dgamma_bb**2)
            d2Eab(10)=(Ep_4-2.0d0*Eb+Em_4)/(hg_bb**2)
         end if
      end if
      
      ! If alpha is equal to beta, assing derivatives
      ! that was not calculated above
      if (.not. alpha_not_equal_beta) then
         if (norder .ge. 1) then
            d1Eab(2)=d1Eab(1)
            d1Eab(4)=d1Eab(3)
         end if
         if (norder .ge. 2) then
            d2Eab(3)=d2Eab(1)
            d2Eab(8)=d2Eab(4)
            d2Eab(10)=d2Eab(6)
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
      end if
      
      end subroutine

C****************************************************************************
      SUBROUTINE DELA(ESR,rhoc,grdcc,VLAMBDA,norder,d1E,d2E)
C Linear complement exchange PBE functional, Kamal SHARKAS 20/05/2011
C*****************************************************************************

      implicit none

      external ESR

      real*8   rhoc, grdcc, VLAMBDA
      integer  norder
      real*8   E, d1E(5), d2E(15) 

      call   DESRC_GGA(ESR,rhoc,grdcc,0.d0,norder,E,d1E,d2E)

      d1E(1) = (1.d0- VLAMBDA) * d1E(1)
      d1E(3) = (1.d0- VLAMBDA) * d1E(3)

      d2E(1) = (1.d0- VLAMBDA) * d2E(1)
      d2E(4) = (1.d0- VLAMBDA) * d2E(4)
      d2E(6) = (1.d0- VLAMBDA) * d2E(6)

      return
      end

C****************************************************************************
      SUBROUTINE DELANSC(ESR,rhoc,grdcc,VLAMBDA,norder,d1E,d2E)

C Linear complement (NON-Scaled) correlation PBE GWS functional, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      implicit none

      external ESR
      real*8   rhoc, grdcc, VLAMBDA
      integer  norder
      real*8   E, d1E(5), d2E(15)

      call   DESRC_GGA(ESR,rhoc,grdcc,0.d0,norder,E,d1E,d2E)

      d1E(1) = (1.d0- VLAMBDA**2) * d1E(1)
      d1E(3) = (1.d0- VLAMBDA**2) * d1E(3)

      d2E(1) = (1.d0- VLAMBDA**2) * d2E(1)
      d2E(4) = (1.d0- VLAMBDA**2) * d2E(4)
      d2E(6) = (1.d0- VLAMBDA**2) * d2E(6)

      return
      end

C****************************************************************************
      SUBROUTINE DELASC(ESR,rhoc,grdcc,VLAMBDA,norder,d1E,d2E)

C Linear complement Scaled correlation PBE functional, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      implicit none

      external ESR
      integer norder

      real*8   rhoc, rhocscaled, grdcc, grdccscaled, VLAMBDA
      real*8   d1E(5), d2E(15)
      real*8   Ecoul, d1Ecoul(5), d2Ecoul(15)
      real*8   Esca,  d1Esca(5),  d2Esca(15)

! TODO HJAAJ 
      call quit('hjaaj: need to define new d1E and d2E addresses')

!     d1Ecoul(1) = d1Edrhoc_coul
!     d1Ecoul(3) = d1Edgrdcc_coul
!     d2Ecoul(1) = d2Edrhocdrhoc_coul
!     d2Ecoul(4) = d2Edrhocdgrdcc_coul
!     d2Ecoul(6) = d2Edgrdccdgrdcc_coul 
      d1Ecoul(:) = 0.0d0
      d2Ecoul(:) = 0.0d0

      call  DESRC_GGA(ESR,rhoc,grdcc,0.d0,norder,Ecoul,d1Ecoul,d2Ecoul)

      rhocscaled = rhoc/VLAMBDA**3
      grdccscaled  = grdcc/VLAMBDA**8

!     d1Esca(1) = d1Edrhoc_sca
!     d1Esca(3) = d1Edgrdcc_sca 
!     d2Esca(1) = d2Edrhocdrhoc_sca 
!     d2Esca(4) = d2Edrhocdgrdcc_sca
!     d2Esca(6) = d2Edgrdccdgrdcc_sca 
      d1Esca(:) = 0.0d0
      d2Esca(:) = 0.0d0

      call  DESRC_GGA(ESR,rhocscaled,grdccscaled,0.d0,norder,
     &   Esca,d1Esca,d2Esca)

      d1E(1) = d1E(1) + (d1Ecoul(1) - VLAMBDA**2 * d1Esca(1))
      d1E(3) = d1E(3) + (d1Ecoul(3) - (d1Esca(3)/VLAMBDA**3))

      d2E(1) = d2E(1) + (d2Ecoul(1) - (d2Esca(1)/VLAMBDA))
      d2E(4) = d2E(4) + (d2Ecoul(4) - (d2Esca(4)/VLAMBDA**2))
      d2E(6) = d2E(6) + (d2Ecoul(6) - (d2Esca(6)/VLAMBDA**3))

      return
      end

C****************************************************************************
      SUBROUTINE DELSPIN(ESRSPIN,rhoc,rhos,VLAMBDA,dospin,
     &                   norder,d1E,d2E)

C linear complement spin-dependent exchange LDA functional, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      implicit none

      external ESRSPIN
      integer norder
      logical dospin
      real*8   rhoc, rhos, VLAMBDA, rho(2)
      real*8   E, d1E(5), d2E(15)
      
      rho(1) = rhoc
      rho(2) = rhos
      if (dospin) then
         call DESRC_SPINLDA(ESRSPIN,rho,0.d0,norder,E,d1E,d2E)
      else
         call DESRC_LDA(ESRSPIN,rho,0.d0,norder,E,d1E,d2E)
      end if

      d1E(1) = (1.d0- VLAMBDA) * d1E(1)

      if (dospin) d1E(2) = (1.d0- VLAMBDA) * d1E(2)

      d2E(1) = (1.d0- VLAMBDA) * d2E(1)

      if (dospin) d2E(2) = (1.d0 - VLAMBDA) * d2E(2)
      if (dospin) d2E(3) = (1.d0 - VLAMBDA) * d2E(3)

      return
      end

C****************************************************************************
      SUBROUTINE DELSPINnSC(ESRSPIN,rhoc,rhos,VLAMBDA,dospin,
     &                      norder,d1E,d2E)
C linear complement non-scaled spin-dependent correlation LDA functional, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      implicit none

      external ESRSPIN
      integer norder
      logical dospin
      real*8   rhoc, rhos, VLAMBDA, rho(2)
      real*8   E, d1E(5), d2E(15)

      rho(1) = rhoc
      rho(2) = rhos
      if (dospin) then
         call DESRC_SPINLDA(ESRSPIN,rho,0.d0,norder,E,d1E,d2E)
      else
         call DESRC_LDA(ESRSPIN,rho,0.d0,norder,E,d1E,d2E)
      end if

      d1E(1) = (1.d0- VLAMBDA**2)* d1E(1)
      if (dospin) d1E(2) = (1.d0- VLAMBDA**2)* d1E(2)

      d2E(1) = (1.d0- VLAMBDA**2)* d2E(1)
      if (dospin) d2E(2) = (1.d0- VLAMBDA**2)* d2E(2)
      if (dospin) d2E(3) = (1.d0- VLAMBDA**2)* d2E(3)

      return
      end

C****************************************************************************
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

      real*8  rhoa, grdaa, mu, e

      real*8  exerflda, exerfpbe
      real*8  rho,rho13,grd2,t1,t2,t3,t4

      real*8  XKSR(3)
      logical ERFEXP(0:2)

! function 
      real*8  berft

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

C****************************************************************************
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

      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8   rho13,kF,s
      real*8   rho, grd, exLDA, FxHSE

C     constants
      real*8, parameter :: pi = 3.141592653589793d0

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
C------------------------------------------------------------------
      Subroutine HSEFx(Rho,s,Omega,FxHSE)
C------------------------------------------------------------------
C    Short-range exchange PBE enhancement factor
C    from Heyd and Scuseria, JCP 120, 7274 (2004)
C
C    The Erfc function is approximated instead of the PBE exchange hole
C    to perform the integration over r12
C
C    routine taken from PhD thesis of J. Heyd
C    modified by J. Toulouse:
C    gexp  -> dexp
C    gerf  -> derf
C    gerfc -> derfc
C    gei(x)-> ei
C
C    Should be the same as wpbe_analytical_erfc_approx
C------------------------------------------------------------------
      Implicit Real*8(A-H,O-Z)
      Real*8 rho,s,omega,Fxhse
C
      Real*8 A,B,C,D,E
      Real*8 Ha1,Ha2,Ha3,Ha4,Ha5
      Real*8 Fc1,Fc2
      Real*8 ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8,eb1,wcut
      Real*8 EGscut,EGa1,EGa2,EGa3
      Real*8 expcut,exei1,exei2,exei3,exei4
      Real*8 smax,strans,sconst
C
      Real*8 Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
      Real*8 FIfteen,Sixteen
      Real*8 r12,r64,r36,r81,r256,r384
      Real*8 r18,r20,r25,r27,r128,r144,r288,r324,r729
      Real*8 r30,r32,r75,r243,r2187,r6561,r40,r105,r54,r135
      Real*8 r1215,r15309
      Real*8 f12,f14,f32,f34,f94,f98,f1516
C
      Save A,B,C,D,E
      Save Ha1,Ha2,Ha3,Ha4,Ha5
      Save Fc1,Fc2
      Save ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8,wcut
      Save EGscut,EGa1,EGa2,EGa3
      Save expcut,exei1,exei2,exei3,exei4
      Save smax,strans,sconst
      Save Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
      Save Fifteen,Sixteen
      Save r36,r64,r81,r256,r384
      Save r27,r128,r144,r288,r324,r729
      Save r18,r20,r32,r243,r2187,r6561,r40
      Save r12,r25,r30,r54,r75,r105,r135,r1215,r15309
      Save f12,f14,f32,f34,f94,f98,f1516
C
C         Constants for PBE hole
C

      Data A,B,C,D,E
     $     / 1.0161144D0,-3.7170836D-1,-7.7215461D-2,
     $        5.7786348D-1,-5.1955731D-2 /
C
C         Constants for fit of H(s) (PBE hole)
C
      Data Ha1,Ha2,Ha3,Ha4,Ha5
     $     / 9.79681d-3,4.10834d-2,1.87440d-1,1.20824d-3,3.47188d-2 /
C
C         Constants for F(H) (PBE hole)
C
      Data Fc1,Fc2 / 6.4753871d0,4.7965830d-1 /
C
C         Constants for expansion of erfc(x) (eb1 set later
C         depending on wcut)
C
      Data ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8,wcut
     $     / -1.128223946706117d0,1.452736265762971d0,
     $        -1.243162299390327d0,0.971824836115601d0,
     $        -0.568861079687373d0,0.246880514820192d0,
     $        -0.065032363850763d0,0.008401793031216d0,1.4D1 /
C
C         Constants for polynomial expansion of EG for small s
C
      Data EGscut,EGa1,EGa2,EGa3
     $     / 8.0d-2,-2.628417880d-2,-7.117647788d-2,8.534541323d-2 /
C
C         Constants for large x expansion of exp(x)*ei(-x)
C
      Data expcut,exei1,exei2,exei3,exei4
     $     / 7.0D2,4.03640D0,1.15198D0,5.03627D0,4.19160D0 /
C
C         Constants for enforcement of local Lieb-Oxford bound
C
      Data smax,strans,sconst
     $     / 8.572844D0,8.3D0,1.879622316D1 /
C
C         Whole numbers used during evaluation
C
      Data Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
     $     / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
      Data FIfteen,Sixteen / 1.5D1, 1.6D1 /
      Data r36,r64,r81,r256,r384
     $     / 3.6D1,6.4D1,8.1D1,2.56D2,3.84D2 /
      Data r27,r128,r144,r288,r324,r729
     $         / 2.7D1,1.28D2,1.44D2,2.88D2,3.24D2,7.29D2 /
      Data r18,r20,r32,r243,r2187,r6561,r40
     $         / 1.8d1,2.0d1,3.2D1,2.43D2,2.187D3,6.561D3,4.0d1 /
      Data r12,r25,r30,r54,r75,r105,r135,r1215,r15309
     $         / 1.2D1,2.5d1,3.0d1,5.4D1,7.5d1,1.05D2,1.35D2,1.215D3,
     $            1.5309D4 /
C
C         Fractions used during evaluation
C
      Data f12,f14,f32,f34,f94,f98,f1516
     $         / 0.5D0,0.25D0,1.5D0,0.75D0,2.25D0,1.125D0,0.9375D0 /
C
C         General constants
C
          f13      = One/Three
          pi       = ACos(-One)
          pi2      = pi*pi
          srpi = sqrt(pi)
          X        = -Eight/Nine
C
C         Cutoff criterion to enforce local Lieb-Oxford bound
C         This ensures that the enhancement factor does not exceed the
C         one of the original PBE functional (Fx(max)=1.804).
C
          sreal               = s
          If(s.gt.strans) then
               s = smax-(sconst/sreal**2)
               stra1s = Two*sconst/sreal**3
          else
               stra1s = One
          endIf
C
C    Calculate prelim variables
C
      xkf      = (Three*pi2*rho) ** f13
      A2      = A*A
      A3      = A2*A
      A4      = A3*A
      A12     = Sqrt(A)
      A32     = A12*A
      A52     = A32*A
      w       = omega / xkf
      w2      = w * w
      w3      = w2 * w
      w4      = w2 * w2
      w5      = w3 * w2
      w6      = w5 * w
      w7      = w6 * w
      w8      = w7 * w
      s2      = s*s
      s3      = s2*s
      s4      = s2*s2
      s5      = s4*s
      s6      = s5*s
C
C    Calculate H(s) and F(H(s)) for the PBE hole
C
      Hnum    = Ha1*s2 + Ha2*s4
      Hden    = One + Ha3*s4 + Ha4*s5 + Ha5*s6
      H        = Hnum/Hden
      Hnu1s = Two*Ha1*s + Four*Ha2*s3
      Hde1s = Four*Ha3*s3 + Five*Ha4*s4 + Six*Ha5*s5
      H1s      = (Hden*Hnu1s - Hnum*Hde1s) / (Hden*Hden)
      F        = Fc1*H + Fc2
C
C    Set exponent of Gaussian in approx depending on which approx
C    we use for erfc(x)
C
      If(w .lt. wcut) then
          eb1 = 1.455915450052607d0
      else
          eb1 = 2.0d0
      endIf
C
C    Calculate helper variables
C
      Hsbw = s2*H + eb1*w2
      Hsbw2     = Hsbw*Hsbw
      Hsbw3     = Hsbw2*Hsbw
      Hsbw4     = Hsbw3*Hsbw
      Hsbw6     = Hsbw3*Hsbw3
      Hsbw12 = Sqrt(Hsbw)
      Hsbw32 = Hsbw12*Hsbw
      Hsbw52 = Hsbw32*Hsbw
      Hsbw72 = Hsbw52*Hsbw
      DHsb      = D + s2*H + eb1*w2
      DHsb2     = DHsb*DHsb
      DHsb3     = DHsb2*DHsb
      DHsb4     = DHsb3*DHsb
      DHsb5     = DHsb4*DHsb
      DHsb12 = Sqrt(DHsb)
      DHsb32 = DHsb12*DHsb
      DHsb52 = DHsb32*DHsb
      DHsb72 = DHsb52*DHsb
      DHsb92 = DHsb72*DHsb
      HA94     = f94 * Hsbw / A
      HA942    = HA94*HA94
      HA943    = HA942*HA94
      HA945    = HA943*HA942
      HA9412 = Sqrt(HA94)
      DHs      = D + s2*H
      DHs2     = DHs*DHs
      DHs3     = DHs2*DHs
      DHs4     = DHs3*DHs
      DHs72    = DHs3*sqrt(DHs)
      DHs92    = DHs72*DHs
      DHsw     = DHs + w2
      DHsw2    = DHsw*DHsw
      DHsw52 = sqrt(DHsw)*DHsw2
      DHsw72 = DHsw52*DHsw
C
C         Calculate EG(s), using expansion for small s If necessary
C
          If(s .gt. EGscut) then
           G_a     = srpi * (FIfteen*E + Six*C*(One+F*s2)*DHs +
     $                           Four*B*(DHs2) + Eight*A*(DHs3))
     $                    * (One / (Sixteen * DHs72))
     $                    - f34*pi*sqrt(A) * dexp(f94*H*s2/A) *
     $                          (One - derf(f32*s*sqrt(H/A)))
           G_b     = (f1516 * srpi * s2) / DHs72
           EG      = - (f34*pi + G_a) / G_b
          else
           EG      = EGa1 + EGa2*s2 + EGa3*s4
          endIf
C
C    Calculate the terms needed in any case
C
      tm2      = (DHs2*B + DHs*C + Two*E + DHs*s2*C*F + Two*s2*EG) /
     $               (Two*DHs3)
      tm3      = - w    * (Four*DHsw2*B + Six*DHsw*C + FIfteen*E
     $                         + Six*DHsw*s2*C*F + FIfteen*s2*EG) /
     $                         (Eight*DHs*DHsw52)
      tm4      = - w3 * (DHsw*C + Five*E + DHsw*s2*C*F + Five*s2*EG) /
     $                         (Two*DHs2*DHsw52)
      tm5      = - w5 * (E + s2*EG) / (DHs3*DHsw52)
C
C         Calculate t10 unless that would generate a div. by zero
C
          If((s.gt.0.0d0).or.(w.gt.0.0d0)) then
           t10       = (f12)*A*Log(Hsbw / DHsb)
          endIf
C
C         Calculate exp(x)*f(x) depending on size of x
C
          If(HA94 .lt. expcut) then
           exer      = pi*dexp(HA94)*derfc(HA9412)
           exHA94 = dexp(HA94)
           eiHA94 = ei(-HA94)
           exei      = exHA94*eiHA94
          else
           exer      = pi*(One/(srpi*HA9412)
     $                      - One/(Two*Sqrt(pi*HA943))
     $                      + Three/(Four*Sqrt(pi*HA945)))
           exei      = - (One/HA94) *
     $                      (HA942 + exei1*HA94 + exei2) /
     $                         (HA942 + exei3*HA94 + exei4)
          endIf
C
          If (w .eq. Zero) then
C
C          Fall back to original expression for the PBE hole
C
           t1            = -f12*A*exei
C
           If(s .gt. 0.0D0) then
                  tm1        = t1 + t10
                  Fxhse      = X * (tm1 + tm2)
           else
                  Fxhse      = 1.0d0
           endIf
C
          else If(w .gt. wcut) then
C
C          Use simple Gaussian approximation for large w
C
           tm1             = -f12*A*(exei+log(DHsb)-log(Hsbw))
           Fxhse           = X * (tm1 + tm2 + tm3 + tm4 + tm5)
          else
C
C          For everything else, use the full blown expression
C
C          First, we calculate the polynomials for the first term
C
           pn1             = -f32*ea1*A12*w + r27*ea3*w3/(Eight*A12)
     $                 - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)
          pn2     = -A + f94*ea2*w2 - r81*ea4*w4/(Sixteen*A)
     $                 + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)
C
C         The first term is
C
          t1      = f12*(pn1*exer + pn2*exei)
C
C         The factors for the main polynomoal in w
C
          f2      = (f12)*ea1*srpi*A / DHsb12
          f3      = (f12)*ea2*A / DHsb
          f4      =    ea3*srpi*(-f98 / Hsbw12
     $                  + f14*A / DHsb32)
          f5      = ea4*(One/r128) * (-r144*(One/Hsbw)
     $                  + r64*(One/DHsb2)*A)
          f6      = ea5*(Three*srpi*(Three*DHsb52*(Nine*Hsbw-Two*A)
     $                 + Four*Hsbw32*A2))
     $                 / (r32*DHsb52*Hsbw32*A)
          f7      = ea6*(((r32*A)/DHsb3
     $                  + (-r36 + (r81*s2*H)/A)/Hsbw2)) / r32
          f8      = ea7*(-Three*srpi*(-r40*Hsbw52*A3
     $                  +Nine*DHsb72*(r27*Hsbw2-Six*Hsbw*A+Four*A2)))
     $                  / (r128 * DHsb72*Hsbw52*A2)
       f9      = (r324*ea6*eb1*DHsb4*Hsbw*A
     $             + ea8*(r384*Hsbw3*A3 + DHsb4*(-r729*Hsbw2
     $             + r324*Hsbw*A - r288*A2))) / (r128*DHsb4*Hsbw3*A2)
C
          t2t9    = f2*w     + f3*w2 + f4*w3 + f5*w4 + f6*w5
     $                 + f7*w6 + f8*w7 + f9*w8
C
C     The final value of term1 for 0 < omega < wcut is:
C
      tm1     = t1 + t2t9 + t10
C
C     The final value for the enhancement factor is
C
      Fxhse = X * (tm1 + tm2 + tm3 + tm4 + tm5)
      endIf
C
      end

C------------------------------------------------------------------
      subroutine wpbe_analytical(rho,s,omega,Fx_wpbe)
C------------------------------------------------------------------
C    Short-range exchange PBE enhancement factor
C    from Heyd, Scuseria, Ernzerhof, JCP 118, 8207 (2003)
C
C    The PBE exchange hole is approximated to perform the integration
C    over r12
C    
c    'wPBE Enhancement Factor (analytical, no gradients)'
C    routine from J. Heyd
C
C    Gives sligthly different results compared to wpbe_analytical_erfc_approx
c--------------------------------------------------------------------
      Implicit None

      Real*8 rho,s,omega,Fx_wpbe

      Real*8 f12,f13,f14,f18,f23,f43,f32,f72,f34,f94,f1516
      Real*8 pi,pi2,pi_23,srpi
      Real*8 Three_13
    
      Real*8 fa1,fa2,fa3,fa4,fa5
      Real*8 fb1,fb2,fb3,fb4,fb5
      Real*8 A,B,C,D,E
      Real*8 Ha1,Ha2,Ha3,Ha4,Ha5
      Real*8 Fc1,Fc2
      Real*8 EGa1,EGa2,EGa3
      Real*8 EGscut

      Real*8 xkf, xkfrho
      Real*8 wk,wk2,wk3,wk5
      Real*8 X
      Real*8 s2,s4

      Real*8 H,xF,C_0,H_0,DH_0,DH_072
      Real*8 aDH,aDH2,aDHwk2,sraDH,sraDH2,sraDH3,sraDH5
      Real*8 G_a,G_b,EG,E_0

      Real*8 ab1H,srab1H,ab2H,srab2H,b3H,srb3H,srb3H3
      Real*8 b4H,srb4H,srb4H3,ab5H,srab5H,srb5H3

      Real*8 atb3H,atb4H

      Real*8 term1,term2,term3,term4,term5,term6,term7,term8

      Real*8 ax,um,uk,ul
      Real*8 gc1,gc2

      Real*8 derf

      Real*8 One,Two,Three,Four,Six,Eight,Nine,Ten
      Real*8 Fifteen,Sixteen

      Save One,Two,Three,Four,Six,Eight,Nine,Ten
      Data One,Two,Three,Four,Six,Eight,Nine,Ten
     $  / 1D0,2D0,3D0,4D0,6D0,8D0,9D0,10D0 /
      Save Fifteen,Sixteen
      Data Fifteen,Sixteen / 1.5D1, 1.6D1 /

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
      f1516  = Fifteen / Sixteen

      pi     = ACos(-One)
      pi2    = pi*pi
      pi_23  = pi2**f13
      srpi   = sqrt(pi)

      Three_13 = Three**f13

c     Constants from fit

      fa1    = 2.054837628d-4
      fa2    = 1.094652397d-1
      fa3    = 6.407879991d-2
      fa4    = 8.181734694d-3
      fa5    = 1.106664933d-4

      fb1    = 6.601358926d-3
      fb2    = 2.599311408d-1
      fb3    = 5.203522409d-1
      fb4    = 1.185510426d-1
      fb5    = 4.600377701d-2

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

c     Cutoff criterion below which to use polynomial expansion

      EGscut =  8.0d-2

c     Constants for LDA and PBE (from paper)

      ax     = -0.930525736348d0
      um     = 0.2195149727645171d0
      uk     = 0.8040d0
      ul     = um/uk

c     Constants for PBE (from Gaussian99)

      gc1    = 0.00336116d0
      gc2    = 0.00449267d0

c     Calculate the wPBE enhancement factor

      xkf    = (Three*pi2*rho) ** f13
      xkfrho = xkf * rho

      wk     = omega / xkf
      wk2    = wk * wk
      wk3    = wk2 * wk
      wk5    = wk3 * wk2

      X      = - Eight/Nine

      s2     = s*s
      s4     = s2*s2

c     Calculate wPBE enhancement factor

      H      = (Ha1*s2 + Ha2*s4) /
     1         (One + Ha3*s4 + Ha4*s4*s + Ha5*s4*s2)

      xF     = Fc1*H + Fc2

      C_0    = C * (One + s2*xF)
      H_0    = s2 * H
      DH_0   = D + H_0
      DH_072 = DH_0**f72

      aDH    = wk2 + D + H_0
      aDH2   = aDH*aDH
      aDHwk2 = aDH - wk2
      sraDH  = sqrt(aDH)
      sraDH2 = sraDH*sraDH
      sraDH3 = sraDH2*sraDH
      sraDH5 = sraDH2*sraDH3

      if(s .gt. EGscut) then

        G_a    = srpi * (Fifteen*E + Six*C_0*DH_0 +
     1                   Four*B*(DH_0**2) + Eight*A*(DH_0**3))
     2                * (One / (Sixteen * DH_072))
     3                 - f34*pi*sqrt(A) * exp(f94*H_0/A) *
     4                   (One - derf(f32*s*sqrt(H/A)))

        G_b    = (f1516 * srpi * s2) / DH_072

        EG     = - (f34*pi + G_a) / G_b

      else

        EG = EGa1 + EGa2*s2 + EGa3*s4

      endif

      E_0    = E + s2*EG

      ab1H   = wk2 + fb1 + H_0
      srab1H = sqrt(ab1H)
      ab2H   = wk2 + fb2 + H_0
      srab2H = sqrt(ab2H)
      b3H    = fb3 + H_0
      srb3H  = sqrt(b3H)
      srb3H3 = srb3H*srb3H*srb3H
      b4H    = fb4 + H_0
      srb4H  = sqrt(b4H)
      srb4H3 = srb4H*srb4H*srb4H
      ab5H   = wk2 + fb5 + H_0
      srab5H = sqrt(ab5H)
      srb5H3 = srab5H*srab5H*srab5H

      term1  = f12 * (fa1 * (wk - srab1H)) /
     1               (srab1H * (wk2 - ab1H))

      term2  = f12 * (fa2 * (wk - srab2H)) /
     1               (srab2H * (ab2H - wk2))

      atb3H  = atan(srb3H / wk)

      term3  = f12 * (fa3 * (atb3H*wk2 + atb3H*b3H - wk*srb3H)) /
     1               (srpi*srb3H3*(wk2 + b3H))

      atb4H  = atan(srb4H / wk)

      term4  = f12 * (fa4 * (atb4H*wk2 + atb4H*b4H - wk*srb4H)) /
     1               (srpi*srb4H3*(wk2 + b4H))

      term5  = f14 * (fa5 * (wk3 - Three*ab5H*wk + Two*srb5H3)) / 
     1               (srb5H3 * (wk2 - ab5H)**2)

      term6  = f12 * (B * (wk - sraDH)) / 
     1               (sraDH * aDHwk2)

      term7  = f14 * (C_0 * (Three*aDH*wk - wk3 - Two*sraDH3)) /
     1               (sraDH3 * aDHwk2**2)

      term8  = f18 * (E_0 * (Ten*aDH*wk3 - Three*wk5 -
     1               Fifteen*aDH2*wk + Eight*sraDH5)) /
     1               (sraDH5 * aDHwk2**3) 

      Fx_wpbe = X * (-term1 + term2 - term3 - term4 - term5 -
     1                term6 - term7 + term8)

      end

C------------------------------------------------------------------
      subroutine wpbe_analytical_erfc_approx(rho,s,omega,Fx_wpbe)
C------------------------------------------------------------------
C    Short-range exchange PBE enhancement factor
C    from Heyd and Scuseria, JCP 120, 7274 (2004)
C
C    The Erfc function is approximated instead of the PBE exchange hole
C    to perform the integration over r12
C    
c    'wPBE Enhancement Factor (erfc approx.,analytical, no gradients)'
C    routine from J. Heyd
C
C    It is the preferred routine. Should be the same as HSEFx
c--------------------------------------------------------------------
      Implicit None

      Real*8 rho,s,omega,Fx_wpbe

      Real*8 f12,f13,f14,f18,f23,f43,f32,f72,f34,f94,f1516,f98
      Real*8 pi,pi2,pi_23,srpi
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

      Real*8 derf,derfc,ei,fexpei

      Real*8 Zero,One,Two,Three,Four,Six,Eight,Nine,Ten
      Real*8 Fifteen,Sixteen
      Real*8 r64,r36,r81,r256,r384,r864,r1944,r4374
      Real*8 r27,r48,r120,r128,r144,r288,r324,r512,r729
      Real*8 r32,r243,r2187,r6561,r40

      Save Zero,One,Two,Three,Four,Six,Eight,Nine,Ten
      Data Zero,One,Two,Three,Four,Six,Eight,Nine,Ten
     $  / 0D0,1D0,2D0,3D0,4D0,6D0,8D0,9D0,10D0 /
      Save Fifteen,Sixteen
      Data Fifteen,Sixteen / 1.5D1, 1.6D1 /
      Save r36,r64,r81,r256,r384,r864,r1944,r4374
      Data r36,r64,r81,r256,r384,r864,r1944,r4374
     $  / 3.6D1,6.4D1,8.1D1,2.56D2,3.84D2,8.64D2,1.944D3,4.374D3 /
      Save r27,r48,r120,r128,r144,r288,r324,r512,r729
      Data r27,r48,r120,r128,r144,r288,r324,r512,r729
     $  / 2.7D1,4.8D1,1.2D2,1.28D2,1.44D2,2.88D2,3.24D2,5.12D2,7.29D2 /
      Save r32,r243,r2187,r6561,r40
      Data r32,r243,r2187,r6561,r40
     $  / 3.2D1,2.43D2,2.187D3,6.561D3,4.0d1 /

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

      pi     = ACos(-One)
      pi2    = pi*pi
      pi_23  = pi2**f13
      srpi   = sqrt(pi)

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
     1         (One + Ha3*s4 + Ha4*s4*s + Ha5*s4*s2)

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

        G_a    = srpi * (Fifteen*E + Six*C_0*DH_0 +
     1                   Four*B*(DH_0**2) + Eight*A*(DH_0**3))
     2                * (One / (Sixteen * DH_072))
     3                 - f34*pi*sqrt(A) * exp(f94*H_0/A) *
     4                   (One - derf(f32*s*sqrt(H/A)))

        G_b    = (f1516 * srpi * s2) / DH_072

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
     1               (sraDH * aDHw2)

      term3  = f14 * (C_0 * (Three*aDH*w - w3 - Two*sraDH3)) /
     1               (sraDH3 * aDHw2**2)

      term4  = f18 * (E_0 * (Ten*aDH*w3 - Three*w5 -
     1               Fifteen*aDH2*w + Eight*sraDH5)) /
     2               (sraDH5 * aDHw2**3) 


      t10f1 = (f12)*A*Log(Hsbw / DHsbw)

c     Calculate exp(x)*f(x) depending on size of x

      if(HsbwA94 .lt. expfcutoff) then

        experf = Exp(HsbwA94)*dErfc(HsbwA9412)
        expei  = Exp(HsbwA94)*Ei(-HsbwA94)

      else

c       print *,rho,s," LARGE HsbwA94"

        experf = One/(srpi*HsbwA9412) 
     1           - One/(Two*Sqrt(pi*HsbwA943))
     2           + Three/(Four*Sqrt(pi*HsbwA945))

        expei  = - (One/HsbwA94) *
     1             (HsbwA942 + expei1*HsbwA94 + expei2) /
     2             (HsbwA942 + expei3*HsbwA94 + expei4)

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

c       print *,rho,s," LARGE w"

        term1 = -f12*A*(expei+log(DHsbw)-log(Hsbw))

c       Fx_wpbe = X * (term1 - term2 - term3 + term4)
        Fx_wpbe = term1

      else

c       For everything else, use the full blown expression

        np1 = -f32*ea1*A12*w + r27*ea3*w3/(Eight*A12)
     $        - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)

        np2 = -A + f94*ea2*w2 - r81*ea4*w4/(Sixteen*A)
     $        + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)


        t1 = f12*(np1*pi*experf + np2*expei)

        factor1 = One / (r128*DHsbw4*Hsbw3*A2)
        
        t2f1 = (f12)*ea1*srpi*A*w / DHsbw12
        
        t3f1 = (f12)*ea2*A*w2 / DHsbw
        
        t4f1 = srpi*(-f98*ea3 / Hsbw12 
     $           + f14*ea3*A / DHsbw32)*w3
        
        t5f1 = (One/r128) * (-r144*(One/Hsbw) 
     $        + r64*(One/DHsbw2)*A)*ea4*w4
        
        t6f1 = (Three*ea5*srpi*(Three*DHsbw52*(Nine*Hsbw-Two*A)
     $           + Four*Hsbw32*A2) * w5)
     $           / (r32*DHsbw52*Hsbw32*A)
        
        t7f1 = (ea6*((r32*A)/DHsbw3
     $         + (-r36 + (r81*s2*H)/A)/Hsbw2) * w6) / r32
        
        t8f1 = (-Three*ea7*srpi*(-r40*Hsbw52*A3
     $          +Nine*DHsbw72*(r27*Hsbw2-Six*Hsbw*A+Four*A2)) * w7)
     $         / (r128 * DHsbw72*Hsbw52*A2)
        
        t9 = (r324*ea6*eb1*DHsbw4*Hsbw*A 
     $        + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2
     $        + r324*Hsbw*A - r288*A2)))*w8
      
c       The final value of term1 for 0 < omega < wcutoff is:

        term1 = t1 + t2f1 + t3f1 + t4f1 + t5f1 + t6f1 + t7f1 + t8f1
     $             + factor1*t9 + t10f1

        Fx_wpbe = X * (term1 - term2 - term3 + term4)
c       Fx_wpbe = term1

      endif

      end

C****************************************************************************
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

      real*8, intent(in) :: rho, grdcc
      real*8, intent(out) :: e

      real*8   rho13, eLDA, epsLDA
      real*8   kF, ks, t2
      real*8   A, At2
      real*8   frac, H

C     General constants
      real*8, parameter :: pi = 3.141592653589793d0
      real*8, parameter :: pi2 = 9.869604401089358d0
C     PBE constants
      real*8, parameter :: gam = 0.031090690869654897d0
      real*8, parameter :: bet = 0.066725d0

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      call EVWN(eLDA,rho,rho13)
      epsLDA = eLDA/rho

      A = dexp(-epsLDA/gam) - 1.d0

      if (abs(A) .lt. 1.0d-14) then
         H = 0.0d0
      else
         A = (bet/gam)/A
      
         kF = (3.d0*pi2*rho)**(1.d0/3.d0)
         ks = dsqrt(4.d0*kF/pi) 
         t2  = grdcc / (2.d0*ks*rho)**2
         At2 = A*t2

         frac = 1.d0 + At2
         frac = frac/(frac + At2*At2)
         H = gam*dlog(1.d0+(bet/gam)*t2*frac)
      end if

      e = rho*(epsLDA + H)

      return
      end

C****************************************************************************
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

      real*8   rho, grdcc, mu, e

      real*8   mu2, rs, sqrt2
      real*8   d, a, bet, gam
      real*8   grs12, expgrs12, g0
      real*8   e_pbe, eps_pbe, eps_LDA_fac, mu_ref
      real*8   d1, d2, denom, p6, q7, q8, d2od1, mu6

C     General constants
#include "pi.h"
      real*8   mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srcpbe_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4

C     Constants for g0
      parameter ( d = 32.D0/(3.D0*pi) )
      parameter ( a = 3.2581D0 )
      parameter ( bet = 163.44D0 )
      parameter ( gam = 4.7125D0 )
      sqrt2 = sqrt(2.0d0)
      
!      print *,'ESRC_PBERIERF: mu=',mu
!      print *,'ESRC_PBERIERF: rho=',rho,' grdcc=',grdcc
     

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
!      print *,'ESRC_PBERIERF: e_pbe=',e_pbe

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
     &     * ( 1.0d0 + ((DERF(mu_0*mu))**5) * (1.0D0/denom - 1.0D0) )
      end if

!      print *,'ESRC_PBERIERF: e=',e
      return
      end

C****************************************************************************
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

      real*8   rho, grdcc, mu, e

      real*8   mu2, rs, sqrt2
      real*8   d, a, bet, gam, rhogrd
      real*8   grs12, expgrs12, g0, rho13
      real*8   e_lyp, eps_lyp, eps_LDA_fac, mu_ref
      real*8   d1, d2, denom, p6, q7, q8, d2od1, mu6

C     General constants
#include "pi.h"
      real*8   mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srclyp_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4

C     Constants for g0
      parameter ( d = 32.D0/(3.D0*pi) )
      parameter ( a = 3.2581D0 )
      parameter ( bet = 163.44D0 )
      parameter ( gam = 4.7125D0 )
      sqrt2 = sqrt(2.0d0)
      
!      print *,'ESRC_LYPRIERF: mu=',mu
!      print *,'ESRC_LYPRIERF: rho=',rho,' grdcc=',grdcc

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
!      print *,'ESRC_LYPRIERF: e_lyp=',e_lyp

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

!      print *,'ESRC_LYPRIERF: e=',e
      return
      end

C
C  begin Manu
C****************************************************************************
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

      real*8   rhoa, grdaa, mu, e

      real*8   mu2, mu6, pi, sqrtpi 
      real*8   three53, pi13, pi43, zero 
      real*8   rho, epbe
      real*8   d, d2, denom1, denom2
      
!      print *,'ESRX_PBERIERF: mu=',mu
!      print *,'ESRX_PBERIERF: rho=',rho,' grd=',grd
     
C     General constants
      pi = 3.141592653589793d0
      sqrtpi = dsqrt(pi)
      three53 = 3.d0**(5.d0/3.d0) 
      pi13 = pi**(1.d0/3.d0)
      pi43 = pi**(4.d0/3.d0) 

!     Variables
      zero = 0.d0 
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
!      print *,'ESRX_PBERIERF: epbe=',epbe
!      print *,'ESRX_PBERIERF: e=',e
      return
      end
C
C end Manu
C
C****************************************************************************
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

      real*8   rho, grdcc, mu, e

      real*8   e_pbe, e_lda,p,z,Ha,Hb
      real*8   rho13
      logical  ERFEXP(0:2), MULOCAL(0:2)
      real*8   ZKSR(3)

      rho13 = rho**(1.0d0/3.0d0)

!     PBE correlation energy (mu=0)
      call ECPBE(rho,grdcc,e_pbe)
      if (abs(e_pbe) .lt. 1.0D-15) then
         e = 0.0d0
         return
      end if

!    short-range LDA energy
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
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

C****************************************************************************
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

      real*8   rho, grdcc, mu, e

      real*8   derfc
      logical  ERFEXP(0:2), MULOCAL(0:2)
      real*8   mu0,rho13, rs, pi
      real*8   elda,epbe,eldamu
      real*8   ZKSR(3)
C     
!      print *,'ESRC_PBEWIERF: mu=',mu
!      print *,'ESRC_PBEWIERF: rho=',rho,' grdcc=',grdcc

      pi = 3.141592653589793d0
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
      rho13 = rho**(1.d0/3.d0)
      rs = (3.d0/(4.d0*pi*rho))**(1.d0/3.d0)

!     LDA energy
      mu0 = 0.d0
      call VCSRLDA(ZKSR,rho,rho13,mu0,MULOCAL,0,ERFEXP)
      elda = ZKSR(1)
!      print *,'ESRC_PBEWIERF: elda=',elda

!     PBE energy
      call ECPBE(rho,grdcc,epbe)
!      print *,'ESRC_PBEWIERF: epbe=',epbe

!     short-range LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      eldamu = ZKSR(1)
!      print *,'ESRC_PBEWIERF: eldamu=',eldamu

      e = (epbe - elda)*derfc(mu*rs) + eldamu

!      print *,'ESRC_PBEWIERF: e=',e
      return
      end
C****************************************************************************
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

      real*8   rho, grdcc, mu, e

      real*8   rho13
      real*8   exerflda, ecerflda, ecerfpbe, t1

      logical  ERFEXP(0:2), MULOCAL(0:2)
      real*8   XKSR(3), ZKSR(3)

! function 
      real*8   berft

!      print *,'ESRC_PBETCSERF: mu=',mu
!      print *,'ESRC_PBETCSERF: rho=',rho,' grdcc=',grdcc

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     Code generated by Mathematica
C     Energy
       t1 = 9.869604401089358618834491d0*ecerflda
     &    + 3.068528194400546905827679d-1*rho*
     & dlog(1.d0 - (8.401605835889726292406297d-1*grdcc*exerflda
     &  *berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &  rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*grdcc*
     & exerflda*berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )*rho**(-1.1d1/3.d0))/(-1.d0 +
     &   dexp((-3.216396844291482112072009d1*ecerflda)/rho))))/
     & (1.d0 + (7.058698062165630644545894d-1*grdcc**2*exerflda**2*
     &  berft(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))**2*rho**(-2.2d1/3.d0))/
     &  (-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))**2
     &  - (8.401605835889726292406297d-1*grdcc*exerflda*
     &  berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))
     &  *rho**(-1.1d1/3.d0))/(-1.d0 +
     &   dexp((-3.216396844291482112072009d1*ecerflda)/rho))))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      e = ecerfpbe
!      print *,'ESRC_PBETCSERF: e=',e

      return
      end
C****************************************************************************
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

      real*8   rho, grdcc, mu, e

      real*8   rho13, t1
      real*8   exerflda, ecerflda, ecerfpbe, epbe

      logical  ERFEXP(0:2), MULOCAL(0:2)
      real*8   XKSR(3), ZKSR(3)

! function 
      real*8   berft

!      print *,'ESRC_PBETCSERF: mu=',mu
!      print *,'ESRC_PBETCSERF: rho=',rho,' grdcc=',grdcc

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     Code generated by Mathematica
C     Energy
      t1 = 9.869604401089358618834491d0*ecerflda
     &   + 3.068528194400546905827679d-1*rho*
     &    dlog(1.d0 - (8.401605835889726292406297d-1*grdcc*exerflda*
     &      berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &      rho**(-1.1d1/3.d0)*
     &      (1.d0 - (8.401605835889726292406297d-1*grdcc*exerflda*
     &      berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &      rho**(-1.1d1/3.d0))/
     &      (-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))
     &      ) )/
     &      (1.d0 + (7.058698062165630644545894d-1*grdcc**2*exerflda**2*
     &      berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2
     &      *rho**(-2.2d1/3.d0))/
     &   (-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))**2
     & - (8.401605835889726292406297d-1*grdcc*exerflda*
     &    berft(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))
     &          *rho**(-1.1d1/3.d0))/
     &   (-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho)) ))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      call ECPBE(rho,grdcc,epbe)

! modification of ESRC_PBETCS:
      e = max(ecerfpbe,epbe)
!      print *,'ESRC_PBETCSERF: e=',e

      return
      end
!-------------------------------------------
      real*8 function berft(a)
!-------------------------------------------
!  Second-order exchange gradient expansion coefficient for erf 
!  interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
 
      real*8   a
      real*8   derf

      real*8   pi

C     constants
      pi = 3.141592653589793d0

      if(a .lt. 0.02d0) then
!      expansion for small mu to avoid numerical problems
       berft = 7.d0/81.d0 + 56.d0*dsqrt(pi)*a/243.d0 
     >        + (-128.d0/81.d0 + 448.d0*pi/729.d0)*a*a

      else

!      Code generated by Mathematica
       berft =(1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0
     &   + dexp(2.5d-1/a**2)) - 2.d0*a**2*(1.1d1 + 7.d0*dexp(2.5d-1/a**2
     &  ))) )/
     & (a**2*(3.2d1*a**4*(-1.d0 + dexp(2.5d-1/a**2))
     &  - 3.d0*dexp(2.5d-1/a**2)
     &  + 1.417963080724412821838534d1*a*derf(5.d-1/a)*dexp(2.5d-1/a**2)
     &  - 8.d0*a**2*(-2.d0 + 3.d0*dexp(2.5d-1/a**2)) ) )

      endif

      return
      end

!-------------------------------------------
      function dberftda(a)
!-------------------------------------------
!  Derivative of second-order exchange gradient
!  expansion coefficient for erf interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
 
      real*8   a
      real*8   dberftda, derf
      real*8   t1,t2,t3,t4,t5
      real*8   pi

C     constants
      pi = 3.141592653589793d0

      if(a .lt. 0.02d0) then
!      expansion for small mu to avoid numerical problems
       dberftda = 56.d0*dsqrt(pi)/243.d0
     >        + 2.d0*(-128.d0/81.d0 + 448.d0*pi/729.d0)*a

      else

!      Code generated by Mathematica
       t1 = (1.851851851851851851851852d-2*(5.76d2*a**3*(-1.d0 + dexp(2.
     &  5d-1/a**2)) + (7.d0*dexp(2.5d-1/a**2))/a - 7.2d1*a*dexp(2.5d-1/a
     &  **2) - 4.d0*a*(1.1d1 + 7.d0*dexp(2.5d-1/a**2))))/(a**2*(3.2d1*a*
     &  *4*(-1.d0 + dexp(2.5d-1/a**2)) - 3.d0*dexp(2.5d-1/a**2) + 1.4179
     &  63080724412821838534d1*a*derf(5.d-1/a)*dexp(2.5d-1/a**2) - 8.d0
     &  *a**2*(-2.d0 + 3.d0*dexp(2.5d-1/a**2))))
       t2 = -1.851851851851851851851852d-2/a**2
       t3 = -8.d0/a + 1.28d2*a**3*(-1.d0 + dexp(2.5d-1/a**2)) + (1.5d0*d
     &  exp(2.5d-1/a**2))/a**3 + (1.2d1*dexp(2.5d-1/a**2))/a - 1.6d1*a*d
     &  exp(2.5d-1/a**2) + 1.417963080724412821838534d1*derf(5.d-1/a)*d
     &  exp(2.5d-1/a**2) - (7.08981540362206410919267d0*derf(5.d-1/a)*d
     &  exp(2.5d-1/a**2))/a**2 - 1.6d1*a*(-2.d0 + 3.d0*dexp(2.5d-1/a**2)
     &  )
       t4 = (-1.d0 + 1.44d2*a**4*(-1.d0 + dexp(2.5d-1/a**2)) - 2.d0*a**2
     &  *(1.1d1 + 7.d0*dexp(2.5d-1/a**2)))/(3.2d1*a**4*(-1.d0 + dexp(2.5
     &  d-1/a**2)) - 3.d0*dexp(2.5d-1/a**2) + 1.417963080724412821838534
     &  d1*a*derf(5.d-1/a)*dexp(2.5d-1/a**2) - 8.d0*a**2*(-2.d0 + 3.d0*
     &  dexp(2.5d-1/a**2)))**2
       t5 = (-3.703703703703703703703704d-2*(-1.d0 + 1.44d2*a**4*(-1.d0 
     &  + dexp(2.5d-1/a**2)) - 2.d0*a**2*(1.1d1 + 7.d0*dexp(2.5d-1/a**2)
     &  )))/(a**3*(3.2d1*a**4*(-1.d0 + dexp(2.5d-1/a**2)) - 3.d0*dexp(2.
     &  5d-1/a**2) + 1.417963080724412821838534d1*a*derf(5.d-1/a)*dexp(
     &  2.5d-1/a**2) - 8.d0*a**2*(-2.d0 + 3.d0*dexp(2.5d-1/a**2))))
       dberftda = t1 + t2*t3*t4 + t5

      endif

      return
      end

C****************************************************************************
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

      real*8   rho, rho13, grd2
      real*8   exerflda, exerfpbe
      real*8   t1,t2,t3,t4
      real*8   kappa,sq,fx

      real*8   XKSR(3)
      logical  ERFEXP(0:2)

! function 
      real*8   berf

C     LDA energy
      rho = 2.0d0*rhoa ! as if rhoa .eq. rhob
      rho13 = rho**(1.d0/3.d0)
      ERFEXP(0:2) = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      exerflda = XKSR(1)

      grd2 = 4.0d0*grdaa ! as if grdaa .eq. grdbb
      kappa=0.804d0
      sq=grd2*2.6121172985233599567768d-2*rho**(-8.d0/3.d0)
      fx=1.d0+kappa
     &    -kappa/(1.d0+berf(1.616204596739954813d-1
     &            *mu*rho**(-1.d0/3.d0))*sq/kappa)
      exerfpbe=exerflda*fx

      e = 0.5d0*exerfpbe ! only sigma exchange energy

      return
      end
C****************************************************************************
      subroutine ELAX_PBEGWS(rho,grd2,VLAMBDA,e)
C linear complement exchange PBE energy, Kamal SHARKAS 19/05/2011 beg
C*****************************************************************************
      real*8  rho(4), grd2(6), VLAMBDA
      real*8  e

      call ESRX_PBEGWSERF(rho(3),grd2(4),0.d0,e)
!     call ESRX_PBEGWSERF(rhoa  ,grdaa  ,0.d0,e)
      e = (1.d0 - VLAMBDA) * e

      return
       end
C****************************************************************************
      subroutine ESRC_PBEGWSERF(rhoc,grdcc,mu,e)
C*****************************************************************************
C     Short-range PBE correlation functional from
C     Goll, Werner, Stoll, PCCP 7, 3917 (2005)
C     Based on the VWN LDA parameterization
C
C     Input: rhoc   : density
C            grdcc  : grad(rhoc) . grad(rhoc)
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 11-08-09
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoc, grdcc, mu
      real*8, intent(out) :: e

      real*8, parameter :: tol=1d-12
      real*8  rho13

      real*8  ecerflda,eclda
      real*8  ecerfpbe
      real*8  t1
      real*8  alpha, beta, gamma_symbol, Aa,Ab, Ac, tq, arglog
 
      logical ERFEXP(0:2), MULOCAL(0:2)
      real*8  XKSR(3), ZKSR(3)

!     short-range LDA energy (order = 0)
      rho13 = rhoc**(1.d0/3.d0)
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
      call VCSRLDA(ZKSR,rhoc,rho13,mu,MULOCAL,0,ERFEXP)
      ecerflda = ZKSR(1)

!     full-range LDA energy
      call VCSRLDA(ZKSR,rhoc,rho13,0.d0,MULOCAL,0,ERFEXP)
      eclda = ZKSR(1)

      alpha=2.78d0
      gamma_symbol=3.1091d-2

      if ((ecerflda/eclda).le.0.d0) then
         beta=0.d0
      else
         beta=6.6725d-2*(ecerflda/eclda)**alpha
      endif
      tq=grdcc*6.346820607d-2*rhoc**(-7.d0/3.d0)
      Ab=dexp(-ecerflda/(rhoc*gamma_symbol))-1.d0
      if (dabs(Ab).le.dabs(beta*tol)) then
         ecerfpbe=ecerflda
       else
         Aa=beta/(gamma_symbol*Ab)
         if (Aa.lt.tol) Aa=tol
         Ac=1.d0+Aa*tq+Aa**2*tq**2
         arglog=1.d0+beta*(1.d0-1.d0/Ac)/(gamma_symbol*Aa)
         ecerfpbe=ecerflda+rhoc*gamma_symbol*dlog(arglog)
       endif

      e = ecerfpbe

      return
      end
C****************************************************************************
      subroutine ESRC_PBEGWSERF_PW92(rhoc,grdcc,mu,e)
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

      real*8   rhoc, grdcc, mu
      real*8   e

      real*8   tol
      parameter(tol=1d-12)
      real*8   rho13

      real*8   ecerflda,eclda
      real*8   ecerfpbe
      real*8   t1
      real*8   alpha, beta, gamma_symbol, Aa,Ab, Ac, tq, arglog
 
!     short-range LDA energy (order = 0)
      rho13 = rhoc**(1.d0/3.d0)
      call ESRC_LDAERFSPIN(rhoc,0.0d0,mu,e)
      ecerflda =  e

!     full-range LDA energy
      call ESRC_LDAERFSPIN(rhoc,0.0d0,0.0d0,e)
      eclda = e

      alpha=2.78d0
      gamma_symbol=3.1091d-2

      if ((ecerflda/eclda).le.0.d0) then
         beta=0.d0
      else
         beta=6.6725d-2*(ecerflda/eclda)**alpha
      endif
      tq=grdcc*6.346820607d-2*rhoc**(-7.d0/3.d0)
      Ab=dexp(-ecerflda/(rhoc*gamma_symbol))-1.d0
      if (dabs(Ab).le.dabs(beta*tol)) then
         ecerfpbe=ecerflda
       else
         Aa=beta/(gamma_symbol*Ab)
         if (Aa.lt.tol) Aa=tol
         Ac=1.d0+Aa*tq+Aa**2*tq**2
         arglog=1.d0+beta*(1.d0-1.d0/Ac)/(gamma_symbol*Aa)
         ecerfpbe=ecerflda+rhoc*gamma_symbol*dlog(arglog)
       endif

      e = ecerfpbe

      return
      end

C****************************************************************************
      subroutine ESRC_PBEGWSERFSPIN(rhoc,rhos,grdcc,grdss,grdcs,mu,e)
C*****************************************************************************
C     Short-range spin-dependent PBE correlation functional from
C       Goll et al. Chem. Phys. 329 (2006), 276
C
C       This functional depends only on rhoc, rhos, grdcc
C       and _not_ on grdss and grdcs
C
C     Input: rho*   : densities
C            grd*   : grad rho_i . grad rho_j
C            mu     : Interaction parameter
C
C     Ouput: e      : PBEGWS correlation energy density
C
C     Author: Erik Hedegaard (after Molpro code by J. Toulouse + A. Savin)
C     Date  : 11-02-16; input parameters revised apr.2018 hjaaj
C     Date  : 4th-July-18; LDA was not doing spin, changed that. E.R.K  
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoc,rhos,grdcc,grdss,grdcs,mu
      real*8, intent(out) :: e
      real*8   rho13, zeta

      real*8, parameter :: tol=1d-12

      real*8   ecerflda, eclda

      real*8   ecerfpbe
      real*8   t1
      real*8   alpha, beta, gamma_symbol
      real*8   phi, phi2, phi3, phi4
      real*8   Aa, Ab, Ac, tq, arglog

      rho13 = rhoc**(1.d0/3.d0)
C     LDA energy ----
!     In ESRC_PBEGWS VCSRLDA is used for the same purpose, but
!     as VCSRLDA does not include rhos ESRX_LDAERFSPIN is used instead
      call ESRC_LDAERFSPIN(rhoc,rhos,mu,ecerflda)

      ! full range LDA energy (mu = 0.0)
      call ESRC_LDAERFSPIN(rhoc,rhos,0.0d0,eclda)

!     Spin polarisation
!     Define zeta
      zeta = rhos / rhoc ! zeta = (rhoa-rhob)/(rhoa+rhob)

      alpha=2.78d0
      gamma_symbol=3.1091d-2

      if ((ecerflda/eclda).le.0d0) then
!       ! edh define beta and check alpha!
        beta=0d0
      else
        beta=6.6725d-2*(ecerflda/eclda)**alpha
      endif
      phi=((1d0+zeta)**(2d0/3d0)+(1d0-zeta)**(2d0/3d0))/2d0
      phi2=phi*phi
      phi3=phi2*phi
      phi4=phi3*phi
      tq=grdcc*6.346820607d-2*rhoc**(-7d0/3d0)/phi2
      Ab=dexp(-ecerflda/(rhoc*gamma_symbol*phi3))-1d0
      if(dabs(Ab).le.dabs(beta*tol)) then
        ecerfpbe=ecerflda
      else
        Aa=beta/(gamma_symbol*Ab)
        Ac=1d0+Aa*tq+Aa**2*tq**2
        if (Aa.lt.tol) Aa=tol
        arglog=1d0+beta*(1d0-1d0/Ac)/(gamma_symbol*Aa)
        ecerfpbe=ecerflda+rhoc*phi3*gamma_symbol*dlog(arglog)
      endif

       e = ecerfpbe

      return
      end
C****************************************************************************
      subroutine ELAC_PBEGWS(rhoc,grdcc,VLAMBDA,e)
C linear complement non-scaled correlation PBE functional energy, Kamal SHARKAS 20/05/2011 beg
C*****************************************************************************
      real*8   rhoc, grdcc, VLAMBDA
      real*8   e

      call ESRC_PBEGWSERF(rhoc,grdcc,0.d0,e)
      e = (1.d0 - VLAMBDA**2) * e

      return
       end

C****************************************************************************
      subroutine ELASC_PBEGWS(rhoc,grdcc,VLAMBDA,e)
C linear complement Scaled correlation PBE functional energy, Kamal SHARKAS 20/05/2011
C*****************************************************************************
      real*8   rhoc, grdcc, VLAMBDA
      real*8   e, ecoul, escaled,rhocscaled,grdscaled

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

!-------------------------------------------
      real*8 function berf(a)
!-------------------------------------------
!  Second-order exchange gradient expansion coefficient for erf
!  interaction
!  a = mu/(2*kF)
! 
!  Modified version from Goll, Werner, Stoll, PCCP 7, 3917 (2005)
!
!  Author : J. Toulouse
!  Date   : 11-08-09
!-------------------------------------------
      implicit none

      real*8, intent(in) :: a
      real*8, parameter :: pi = 3.141592653589793d0
      real*8  fak, fak2

! function
      real*8   derf

      ! eta=19.0d0
      fak=2.540118935556d0*dexp(-19.0d0*a*a)

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-7.d0+72.d0*a*a)
     >        /(27.d0*(-3.d0-24.d0*a*a+32.d0*a**4+8.d0*dsqrt(pi)*a))

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)
     >        - 23.d0/(358400.d0*a**6)

      else


!     Code generated by Mathematica
      fak2 = dexp(0.25d0/a**2)
      berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0
     &  + fak2) - 2.d0*a**2*(1.1d1 + 7.d0*fak2 ))) /
     &  (a**2*(3.2d1*a**4*(-1.d0 + fak2) - 3.d0*fak2
     &   + 1.417963080724412821838534d1*a*derf(5.d-1/a)*fak2
     &   - 8.d0*a**2*(-2.d0 + 3.d0*fak2)))

      end if

      berf=berf*fak

      return
      end

C****************************************************************************
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

      real*8   rho, grd, mu, mu_new, amul
      real*8   e
      real*8   alpha, x, a1, b1, c1
      parameter(a1 = 0.559463d0, b1 = -0.134086d0, c1 = 0.0296811d0)
      real*8   rho13

      logical ERFEXP(0:2), MULOCAL(0:2)
      real*8   ZKSR(3)
!     local mu parameter 
      rho13 = rho**(1.0d0/3.0d0)
!     Multiplication factor for mu parameter
      x = grd/rho
      alpha = a1 + b1*exp(-c1*x)
      amul = x*alpha
      mu_new = max(amul,mu)
!     short-range LDA energy
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
      call VCSRLDA(ZKSR,rho,rho13,mu_new,MULOCAL,0,ERFEXP)
      e = ZKSR(1)
      return
      end
C****************************************************************************
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

      real*8   rho, grd, mu, amul, mu_new
      real*8   e
      real*8   rsfac, rs, x, a, b, c, d
      parameter(rsfac = 0.62035049089940008660d0)
      parameter(b = 1.741930d0, c = -0.528057d0, d = 0.0720226d0 )
      real*8   rho13

      logical ERFEXP(0:2), MULOCAL(0:2)
      real*8   ZKSR(3)

      rho13 = rho**(1.0d0/3.0d0)
      rs   = rsfac/rho13
!     Multiplying factor for mu parameter
      x = grd/rho
      a = b + c*exp(-d*x)
!     local mu parameter 
      amul = a/rs
      mu_new = max(amul,mu)
!     short-range LDA energy
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
!      print *,'mu = 0.1' , mu
      call VCSRLDA(ZKSR,rho,rho13,mu_new,MULOCAL,0,ERFEXP)
      e = ZKSR(1)
      return
      end

C****************************************************************************
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

      real*8   rhoc, grdcc, mu, e
      real*8   amul, mu_new, x, rs, alpha
      real*8   rsfac, a, b, bfac
      parameter(rsfac = 0.62035049089940008660d0, a = 0.878435d0)
      parameter(b = 0.5676020d0, bfac = 0.814516076948D0)
      real*8   rho13

      logical ERFEXP(0:2), MULOCAL(0:2)
      real*8   ZKSR(3)

      rho13 = rhoc**(1.0d0/3.0d0)
      rs   = rsfac/rho13
      x = sqrt(grdcc)/rhoc
      alpha = a*x**b
!     local mu parameter 
      amul = (alpha*bfac)/sqrt(rs)
      mu_new = max(amul,mu)
!     short-range LDA energy
      ERFEXP(0:2) = .false.
      MULOCAL(0:2) = .false.
!      print *,'mu = 0.1' , mu
      call VCSRLDA(ZKSR,rhoc,rho13,mu_new,MULOCAL,0,ERFEXP)
      e = ZKSR(1)
      return
      end
C****************************************************************************
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

      real*8, intent(in) :: rhoa, grdaa, mu
      real*8, intent(out) :: e

      real*8   derf

      real*8   rho, akf, a, a2, a3, ea

      real*8, parameter :: pi = 3.141592653589793d0
      real*8, parameter :: pisqrt = 1.772453850905516d0
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
     >       *(3.d0/8.d0-a*(pisqrt*derf(0.5d0/a)+
     >       (2.d0*a-4.d0*a3)*dexp(-0.25d0/a2)-3.d0*a+4.d0*a3))

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
C****************************************************************************
      subroutine ESRX_LDAERFSPIN(rhoc,rhos,mu,e)
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

      real*8, intent(in) :: rhoc, rhos, mu
      real*8, intent(out) :: e

      real*8   derf
      real*8   rho, akf, a, a2, a3
      real*8   rhoa, rhob
      real*8   ea, eb
      real*8, parameter :: pi = 3.141592653589793d0
      real*8, parameter :: pisqrt = 1.772453850905516d0
      real*8, parameter :: ckf = 3.093667726280136d0
      real*8, parameter :: f13 = 0.333333333333333d0

C     Spin-unpolarized case
      ! 1.0d-10 is basically zero
      if (dabs(rhos) .lt. 1.0d-10) then
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
     >       *(3.d0/8.d0-a*(pisqrt*derf(0.5d0/a)+
     >       (2.d0*a-4.d0*a3)*dexp(-0.25d0/a2)-3.d0*a+4.d0*a3))

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
     >       *(3.d0/8.d0-a*(pisqrt*derf(0.5d0/a)+
     >       (2.d0*a-4.d0*a3)*dexp(-0.25d0/a2)-3.d0*a+4.d0*a3))

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
     >       *(3.d0/8.d0-a*(pisqrt*derf(0.5d0/a)+
     >       (2.d0*a-4.d0*a3)*dexp(-0.25d0/a2)-3.d0*a+4.d0*a3))

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
C****************************************************************************
      subroutine ELX_LDALAMSPIN (rhoc,rhos,VLAMBDA,e)
C linear complement spin-dependent exchange LDA functional energy, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      real*8   rhoc, rhos, VLAMBDA, e

      real*8   rhoa, rhob, ea, eb

      rhoa=0.5d0*max((rhoc+rhos),0d0)
      call ESRX_LDAERF(rhoa,0.0d0,0.d0,ea)
      rhob=0.5d0*max((rhoc-rhos),0d0)
      call ESRX_LDAERF(rhob,0.0d0,0.d0,eb)

      e = (1.d0 - VLAMBDA) * (ea + eb)

      return
       end
C*****************************************************************************

C****************************************************************************
      subroutine ESRC_LDAERFSPIN(rhoc,rhos,mu,e)
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

      real*8, intent(in) :: rhoc, rhos, mu
      real*8, intent(out) :: e

      real*8   rs
      real*8   ec,ecd,ecz,eclr
      real*8   rhoa, rhob, z
      real*8   ea, eb
      real*8, parameter :: f13 = 0.333333333333333d0
      real*8, parameter :: pi = 3.141592653589793d0
      real*8, parameter :: rsfac = 0.620350490899400d0

C     Spin-unpolarized case
      ! this should be the condition instead
      ! of giving it as an input
      if (dabs(rhos) .lt. 1.0d-10) then

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
C****************************************************************************
      subroutine ESRC_LDAERFSPIN_VWN5(rhoc,rhos,mu,e)
C*****************************************************************************
C     subroutine is modified version of ESRC_LDAERFSPIN
C
C     Input: rhoc   : total density
C            rhos   : spin density
C            mu     : Interaction parameter
C            dospin : use spin density
C
C     Ouput: e      : energy
C
C     Created: 2nd July 2018, Erik Kjellgren
C*****************************************************************************
      implicit none

      real*8, intent(in) :: rhoc, rhos, mu
      real*8, intent(out) :: e

      real*8   rs
      real*8   ec,eclr
      real*8   rhoa, rhob, z
      real*8   ea, eb
      real*8   ecVWN5
      real*8, parameter :: f13 = 0.333333333333333d0
      real*8, parameter :: pi = 3.141592653589793d0
      real*8, parameter :: rsfac = 0.620350490899400d0

C     Spin-unpolarized case

      ! 1.0d-10 is basically zero
      if (dabs(rhos) .lt. 1.0d-10) then

        rs = rsfac/(rhoc**f13)
        ec = ecVWN5(rs,0.0d0)
        call ecorrlrVWN5(rs,0.0d0,mu,eclr)
        e=(ec-eclr)*rhoc

C     Spin-polarized case
      else

        rs=rsfac/(rhoc**f13)
        rhoa=max((rhoc+rhos)*.5d0,1.0d-15)
        rhob=max((rhoc-rhos)*.5d0,1.0d-15)
        z=(rhoa-rhob)/(rhoa+rhob)
        if (abs(z) > 1.0d0) 
     &   call quit('in ESRC_LDAERFSPIN numerical problem for |z| = 1.0')

        ec = ecVWN5(rs,z)
        call ecorrlrVWN5(rs,z,mu,eclr)
        e=(ec-eclr)*rhoc

      endif

      return
      end
C****************************************************************************
      subroutine ELC_LDALAMSPIN (rhoc,rhos,VLAMBDA,dospin,e)
C linear complement spin-dependent non-scaled correlation LDA functional, Kamal SHARKAS 23/04/2011
C*****************************************************************************
      real*8   rhoc, rhos, VLAMBDA
      real*8   e
      logical dospin
      if (dospin) then
         call ESRC_LDAERFSPIN(rhoc,rhos,0.d0,e)
      else
         call ESRC_LDAERFSPIN(rhoc,0.0d0,0.d0,e)
      end if 
      e = (1.d0 - VLAMBDA**2) * e

      return
       end
C*****************************************************************************
C****************************************************************************
      subroutine ESRC_MD_LDAERF (rhoc,rhos,mu,dospin,e)
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

      real*8   rhoc, rhos, mu
      logical dospin
      real*8   e

      real*8   e1

      call ESRC_LDAERFSPIN (rhoc,rhos,mu,dospin,e1)
      call DELTA_LRSR_LDAERF (rhoc,rhos,mu,dospin,e)
      e = e1 + e

      end

C****************************************************************************
      subroutine DELTA_LRSR_LDAERF (rhoc,rhos,mu,dospin,e)
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

      real*8   rhoc, rhos, mu
      logical dospin
      real*8   e

      real*8   f13, f83, pi, rsfac, alpha2
      real*8   rs, rs2, rs3

      real*8   rhoa, rhob, z, z2, onepz, onemz, zp, zm, phi8
      real*8   g0f, g0s
      real*8   bd2, bd3
      real*8   c45, c4, c5
      real*8   bc2, bc4, bc3t, bc5t, d0
      real*8   delta2,delta3,delta4,delta5,delta6
      real*8   delta

      parameter(f13 = 0.333333333333333d0)
      parameter(f83 = 2.6666666666666665d0)
      parameter(pi = 3.141592653589793d0)
      parameter(rsfac = 0.620350490899400d0)
      parameter(alpha2 = 0.2715053589826032d0)

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
 
      bd2=dexp(-0.547d0*rs)*(-0.388d0*rs+0.676*rs2)/rs2
      bd3=dexp(-0.31d0*rs)*(-4.95d0*rs+rs2)/rs3

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
      bc3t=-(1-z2)*g0f(rs)*(2.d0*dsqrt(2.d0)-1)/(2.d0*dsqrt(pi)*rs3)
      bc5t = -3.d0*c5*(3.d0-dsqrt(2.d0))/(20.d0*dsqrt(2.d0*pi)*rs3)

      d0=(0.70605d0+0.12927d0*z2)*rs
      delta2=0.073867d0*(rs**(1.5d0))
      delta3=4*(d0**6)*bc3t+(d0**8)*bc5t;
      delta4=4*(d0**6)*bc2+(d0**8)*bc4;
      delta5=(d0**8)*bc3t;
      delta6=(d0**8)*bc2;

      delta=(delta2*(mu**2)+delta3*(mu**3)+delta4*(mu**4)+delta5*(mu**5)
     >      +delta6*(mu**6))/((1+(d0**2)*(mu**2))**4)

C     multiply by rhoc to get energy density
      e=delta*rhoc

      end

C*****************************************************************************
      subroutine ecorrlr(rs,z,mu,eclr)
C*****************************************************************************
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff parameter mu
ccc gives the correlation energy of the LR gas
ccc  => eclr
      implicit none
      real*8, intent(in) :: rs, z, mu
      real*8, intent(out) :: eclr
      real*8   ec,ecd,ecz
      real*8   alpha,cf,phi
      real*8   g0f,dpol,d2anti,d3anti,Qrpa
      real*8   coe2,coe3,coe4,coe5
      real*8   a1,a2,a3,a4,b0
      real*8, parameter :: pi = 3.141592653589793d0
cc parameters from the fit
      real*8, parameter :: adib = 0.784949d0
      real*8, parameter :: q1a  = -0.388d0
      real*8, parameter :: q2a  = 0.676d0
      real*8, parameter :: q3a  = 0.547d0
      real*8, parameter :: t1a  = -4.95d0
      real*8, parameter :: t2a  = 1.d0
      real*8, parameter :: t3a  = 0.31d0
      
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0

      b0=adib*rs

      d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
      d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
      end if

      call ecPW92(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a3=b0**8*coe3
      a4=b0**6*(b0**2*coe2+4.d0*ec)

      eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+b0**8*mu**8*ec)/((1.d0+b0**2*mu**2)**4)

!      Write(6,*)'Qrpa=',Qrpa(mu*sqrt(rs)/phi)
!      Write(6,*)'sqrt(rs)',sqrt(rs)
!      Write(6,*)'phi=',phi
!      Write(6,*)'mu=',mu
!      Write(6,*)'mu*sqrt(rs)/phi=',mu*sqrt(rs)/phi
      return
      end

C*****************************************************************************
      subroutine ecorrlrVWN5(rs,z,mu,eclr)
C
C     Modified version of ecorrlr for VWN5 functional
C
C     Created 2nd July 2018, Erik Kjellgren
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rs, z, mu
      real*8, intent(out) :: eclr
      real*8   ec
      real*8   alpha,cf,phi
      real*8   g0f,dpol,d2anti,d3anti,Qrpa
      real*8   coe2,coe3,coe4,coe5
      real*8   a1,a2,a3,a4,b0
      real*8   ecVWN5
      real*8, parameter :: pi = 3.141592653589793d0
cc parameters from the fit
      real*8, parameter :: adib = 0.784949d0
      real*8, parameter :: q1a  = -0.388d0
      real*8, parameter :: q2a  = 0.676d0
      real*8, parameter :: q3a  = 0.547d0
      real*8, parameter :: t1a  = -4.95d0
      real*8, parameter :: t2a  = 1.d0
      real*8, parameter :: t3a  = 0.31d0

      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0

      b0=adib*rs

      d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
      d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
      end if

      ec = ecVWN5(rs,z)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a3=b0**8*coe3
      a4=b0**6*(b0**2*coe2+4.d0*ec)

      eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+b0**8*mu**8*ec)/((1.d0+b0**2*mu**2)**4)

      return
      end

C*****************************************************************************
      subroutine vcorrlr(rs,z,mu,vclrup,vclrdown)
C*****************************************************************************
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the correlation LSD potential for LR interaction
ccc  => vclrup (spin-up electrons), vclrdown (spin-down electrons)
      implicit none
      real*8   rs,z,mu,eclr,eclrrs,eclrz,vclrup,vclrdown
      real*8   ec,ecd,ecz
      real*8   alpha,cf,phi
      real*8   g0f,dpol,d2anti,d3anti,Qrpa
      real*8   g0d,dpold,d2antid,d3antid,Qrpad,x
      real*8   coe2,coe3,coe4,coe5
      real*8   coe2rs,coe3rs,coe4rs,coe5rs
      real*8   coe2z,coe3z,coe4z,coe5z
      real*8   a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs,
     $     b0rs,a1z,a2z,a3z,a4z,a5z,b0z
      real*8   q1a,q2a,q3a,t1a,t2a,t3a,adib

#include "pi.h"
!     pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
cc parameters from the fit
      adib   = 0.784949d0
      q1a    = -0.388d0
      q2a    = 0.676d0
      q3a    = 0.547d0
      t1a    = -4.95d0
      t2a    = 1.d0
      t3a    = 0.31d0

      b0=adib*rs

      d2anti=(q1a+q2a*rs)*exp(-q3a*rs)/rs
      d3anti=(t1a+t2a*rs)*exp(-t3a*rs)/rs**2

      d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/
     -    rs**2)*exp(-q3a*rs)
      d3antid=-((rs*t2a*(1 + rs*t3a) + t1a*(2 + rs*t3a))/
     -    rs**3)*exp(-rs*t3a)

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0f(rs)-0.5d0)
      coe2rs=-3.d0/8.d0/rs**3*(1.d0-z**2)*g0d(rs)+
     $     9.d0/8.d0/rs**4*(1.d0-z**2)*(g0f(rs)-0.5d0)
      coe2z=-3.d0/8.d0/rs**3*(-2.d0*z)*(g0f(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**3)
      coe3rs=-(1.d0-z**2)*g0d(rs)/(sqrt(2.d0*pi)*rs**3)+
     $    3.d0*(1.d0-z**2)*g0f(rs)/(sqrt(2.d0*pi)*rs**4)
      coe3z=2.d0*z*g0f(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(dpold(rs)
     $        +2.d0*cf**2*2**(5.d0/3.d0)/5.d0/rs**3)
        coe4z=-9.d0/64.d0/rs**3*(dpol(rs)-rs/6.d0*dpold(rs)-2.d0*d2anti
     $       -4.d0/15.d0/rs**2*cf**2*2.d0**(5.d0/3.d0))*z
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)
        coe5rs=-3.d0/rs*coe5-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpold(rs)
        coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(dpol(rs)-rs/6.d0*
     $       dpold(rs)-2.d0*d3anti)*z

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)
         coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(
     $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d2antid+cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.d0-z)**(8.d0/3.d0))/rs**3)
         coe4z=-9.d0/64.d0/rs**3*(1.d0/2.d0*(1.d0+z)*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     $        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     $        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     $        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     $        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d2anti-
     $        4.d0/15.d0/rs**2*cf**2*((1.d0+z)**(5.d0/3.d0)-
     $        (1.d0-z)**(5.d0/3.d0)))

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
         coe5rs=-3.d0/rs*coe5-9.d0/(40.d0*sqrt(2.d0*pi)*rs**3)*(
     $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3antid)
         coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(1.d0/2.d0*(1.d0+z)*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     $        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     $        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     $        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     $        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d3anti)

      end if

      call ecPW92(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a1rs=24.d0*adib*b0**5*coe3+4.d0*b0**6*coe3rs+8.d0*adib*b0**7*
     $     coe5+b0**8*coe5rs
      a1z=4.d0*b0**6*coe3z+b0**8*coe5z

      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a2rs=24.d0*adib*b0**5*coe2+4.d0*b0**6*coe2rs+8.d0*adib*b0**7*
     $     coe4+b0**8*coe4rs+24.d0*adib*b0**3*ec+6.d0*b0**4*ecd
      a2z=4.d0*b0**6*coe2z+b0**8*coe4z+6.d0*b0**4*ecz

      a3=b0**8*coe3
      a3rs=8.d0*adib*b0**7*coe3+b0**8*coe3rs
      a3z=b0**8*coe3z

      a4=b0**6*(b0**2*coe2+4.d0*ec)
      a4rs=8.d0*adib*b0**7*coe2+b0**8*coe2rs+24.d0*adib*b0**5*ec+
     $     4.d0*b0**6*ecd
      a4z=b0**6*(b0**2*coe2z+4.d0*ecz)

      a5=b0**8*ec
      a5rs=8.d0*adib*b0**7*ec+b0**8*ecd
      a5z=b0**8*ecz

      x=mu*sqrt(rs)/phi

      eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+a5*mu**8)/((1.d0+b0**2*mu**2)**4)

      eclrrs=-4.d0/(1.d0+b0**2*mu**2)*2.d0*adib*b0*mu**2*eclr+
     $     1.d0/((1.d0+b0**2*mu**2)**4)*(phi**2*mu/(2.d0*sqrt(rs))
     $     *Qrpad(x)+
     $     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)


      if(z.eq.1.d0) then
         vclrup=eclr-rs/3.d0*eclrrs
         vclrdown=0.d0
      elseif(z.eq.-1.d0) then
         vclrup=0.d0
         vclrdown=eclr-rs/3.d0*eclrrs
      else

         eclrz=(phi**2*((1.d0+z)**(-1.d0/3.d0)-(1.d0-z)**(-1.d0/3.d0))
     $        *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1.d0+z)**(-1.d0/3.d0)
     $        -(1.d0-z)**(-1.d0/3.d0))/3.d0+
     $        a1z*mu**3+a2z*mu**4+a3z*mu**5+
     $        a4z*mu**6+a5z*mu**8)/((1.d0+b0**2*mu**2)**4)

         vclrup=eclr-rs/3.d0*eclrrs-(z-1.d0)*eclrz
         vclrdown=eclr-rs/3.d0*eclrrs-(z+1.d0)*eclrz
      end if
      return
      end

C*****************************************************************************
      real*8 pure function g0f(x)
C*****************************************************************************
ccc on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
ccc x -> rs
      implicit none
      real*8, intent(in) :: x
      real*8, parameter :: C0f = 0.0819306d0
      real*8, parameter :: D0f = 0.752411d0
      real*8, parameter :: E0f = -0.0127713d0
      real*8, parameter :: F0f = 0.00185898d0
      g0f=(1.d0-(0.7317d0-D0f)*x+C0f*x**2+E0f*x**3+
     $     F0f*x**4)*exp(-abs(D0f)*x)/2.d0
      return
      end

C*****************************************************************************
      real*8   function g0d(rs)
C*****************************************************************************
ccc derivative of on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
      implicit none
      real*8   Bg0,Cg0,Dg0,Eg0,Fg0,rs
      Cg0             = 0.0819306d0
      Fg0             = 0.752411d0
      Dg0             = -0.0127713d0
      Eg0             = 0.00185898d0
      Bg0             =0.7317d0-Fg0
      g0d=(-Bg0+2*Cg0*rs+3*Dg0*rs**2+4*Eg0*rs**3)/2.d0*exp(-Fg0*rs)
     -   - (Fg0*(1 - Bg0*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4))/
     -   2.d0*exp(-Fg0*rs)
      return
      end

C*****************************************************************************
      real*8   function dpol(rs)
C*****************************************************************************
      implicit none
      real*8   cf,rs,p2p,p3p
#include "pi.h"
!     pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0
      p3p    = 0.4319d0
      dpol=2.d0**(5.d0/3.d0)/5.d0*cf**2/rs**2*(1.d0+(p3p-0.454555d0)*rs)
     $     /(1.d0+p3p*rs+p2p*rs**2)
      return
      end

      real*8   function dpold(rs)
      implicit none
      real*8   cf,rs,p2p,p3p
#include "pi.h"
!     pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0
      p3p    = 0.4319d0
      dpold=2.d0**(5.d0/3.d0)/5.d0*cf**2*
     - (-2. + (0.454555 - 4.*p3p)*rs +
     -    (-4.*p2p +
     -       (0.90911 - 2.*p3p)*p3p)*rs**2
     -      + p2p*(1.363665 - 3.*p3p)*
     -     rs**3)/
     -  (rs**3*(1. + p3p*rs + p2p*rs**2)**2)
      return
      end

C*****************************************************************************
      real*8   function Qrpa(x)
C*****************************************************************************
      implicit none
      real*8   a2,b2,c2,d2,x,Acoul
#include "pi.h"
!     pi=dacos(-1.d0)
      Acoul=2.d0*(log(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0
      c2              = 3.91744d0
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpa=Acoul*log((1.d0+a2*x+b2*x**2+c2*x**3)/(1.d0+a2*x+d2*x**2))
      return
      end

C*****************************************************************************
      real*8   function Qrpad(x)
C*****************************************************************************
      implicit none
      real*8   a2,b2,c2,d2,x,Acoul
#include "pi.h"
!     pi=dacos(-1.d0)
      Acoul=2.d0*(log(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0
      c2              = 3.91744d0
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpad=Acoul*((x*(b2*(2.d0 + a2*x) +
     -      c2*x*(3.d0 + 2.d0*a2*x) +
     -      d2*(-2.d0 - a2*x + c2*x**3)))/
     -  ((1.d0 + a2*x + d2*x**2)*
     -    (1.d0 + a2*x + b2*x**2 + c2*x**3)))
      return
      end

C*****************************************************************************
      subroutine ecPW92(x,y,ec,ecd,ecz)
C*****************************************************************************
c correlation energy and its derivative w.r.t. rs and z at mu=infinity
c Perdew & Wang PRB 45, 13244 (1992)
c in Hartree; ec=ec(rs,zeta)
c x -> rs; y -> zeta
ccc ecd is d/drs ec
ccc ecz is d/dz ec
      implicit none
      real*8, intent(in) :: x, y
      real*8, intent(out) :: ec, ecd, ecz
      real*8   f02,ff,ec0,ec0d,ec1,ec1d,
     $     aaa,G,Gd,alfac,alfacd
      real*8, parameter :: pi = 3.141592653589793d0

      f02=4.d0/(9.d0*(2.d0**(1.d0/3.d0)-1.d0))

      ff=((1.d0+y)**(4.d0/3.d0)+(1.d0-y)**(4.d0/3.d0)-
     $     2.d0)/(2.d0**(4.d0/3.d0)-2.d0)

      aaa=(1.d0-log(2.d0))/pi**2
      call  GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0,
     $     1.6382d0,0.49294d0,G,Gd)
      ec0 =G
      ec0d=Gd

      aaa=aaa*0.5d0
      call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0,
     $     3.3662d0,0.62517d0,G,Gd)
      ec1 =G
      ec1d=Gd
      call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,
     $     0.88026d0,0.49671d0,G,Gd)
      alfac =-G
      alfacd=-Gd

      ec =ec0 +alfac*ff/f02*(1.d0-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1.d0-y**4)+(ec1d-ec0d)*
     $     ff*y**4
      ecz=alfac*(-4.d0*y**3)*ff/f02+alfac*(1.d0-y**4)/f02*
     $     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     $     (2.d0**(4.d0/3.d0)-2.d0)+(ec1-ec0)*(4.d0*y**3*ff+
     $     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     $     (2.d0**(4.d0/3.d0)-2.d0)*y**4)

      return
      end

C*****************************************************************************
      pure subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd)
C*****************************************************************************
ccc Gd is d/drs G
      implicit none
      real*8, intent(in) :: Ac,alfa1,beta1,beta2,beta3,beta4,x
      real*8, intent(out) :: G, Gd
      G=-2.d0*Ac*(1.d0+alfa1*x)*dlog(1.d0+1.d0/(2.d0*
     $     Ac*(beta1*x**0.5d0+
     $     beta2*x+beta3*x**1.5d0+beta4*x**2)))
      Gd=(1.d0+alfa1*x)*(beta2+beta1/(2.d0*sqrt(x))+3.d0*beta3*
     $     sqrt(x)/2.d0+2.d0*beta4*x)/((beta1*sqrt(x)+beta2*x+
     $     beta3*x**(3.d0/2.d0)+beta4*x**2)**2*(1.d0+1.d0/
     $     (2.d0*Ac*(beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     $     beta4*x**2))))-2.d0*Ac*alfa1*dlog(1.d0+1.d0/(2.d0*Ac*
     $     (beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     $     beta4*x**2)))
      return
      end
C*****************************************************************************
      real*8   function g0s(rs)
C*****************************************************************************
C     g"(0,rs,z=1) from Eq. 32 of
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     Created: 26-08-11, J. Toulouse
C*****************************************************************************
      implicit none
      real*8   rs
      real*8   rs2, f53, alpha2
      parameter(f53 = 1.6666666666666667d0)
      parameter(alpha2 = 0.2715053589826032d0)
      rs2=rs*rs
      g0s=(2.d0**f53)*(1.d0-0.02267d0*rs)/((5.d0*alpha2*rs2)
     >    *(1.d0+0.4319d0*rs+0.04d0*rs2))
      end
C*****************************************************************************
      subroutine HS(z,k,H)
C
C     Smoothening formula for the SRCMULO and SRCPBELO functionals 
C     Cancels the discontinuety arising from the max() function.
C
C     INPUT:   
C             z   :   max(a,b)
C             k   :   
C     OUTPUT:
C             H  
C
C*****************************************************************************
      real*8   z,h,k
! dftcom.h : heaviside_pvalue
#include "dftcom.h"

      H = 1/(1+k*exp(-heaviside_pvalue*z))

      return
      end

C*****************************************************************************
      real*8 pure function ecVWN5(rs, zeta)
C
C     Heavily based on the implementation in libxc
C     
C     Input  :
C              zeta : rhoc / rhos
C              rs   : Wieger-Seitz radius
C     Output :
C              ecVWN5 : Correlation energy
C
C     Erik Kjellgren, 2nd July 2018
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rs, zeta
      real*8 :: t1, t3, t4, t6, t7, t8, t9, t12, t13, t14, t15, t17
      real*8 :: t18, t19, t20, t22, t23, t25, t26, t29, t31, t32, t34
      real*8 :: t37, t38, t39, t41, t42, t44, t45, t46, t47, t49, t50
      real*8 :: t51, t54, t55, t59, t60, t62, t64, t67, t69, t70, t72
      real*8 :: t74, t75, t76

      t1 = sqrt(rs)
      t3 = rs + 0.372744d1 * t1 + 0.129352d2
      t4 = 0.1d1 / t3
      t6 = log(rs * t4)
      t7 = 0.310907d-1 * t6
      t8 = 0.2d1 * t1
      t9 = t8 + 0.372744d1
      t12 = atan(0.61519908197590802322d1 / t9)
      t13 = 0.38783294878113014393d-1 * t12
      t14 = t1 + 0.10498d0
      t15 = t14 * t14
      t17 = log(t15 * t4)
      t18 = 0.96902277115443742139d-3 * t17
      t19 = 0.31415926535897932385d1 * 0.31415926535897932385d1
      t20 = 0.1d1 / t19
      t22 = rs + 0.113107d1 * t1 + 0.130045d2
      t23 = 0.1d1 / t22
      t25 = log(rs * t23)
      t26 = t8 + 0.113107d1
      t29 = atan(0.71231089178181179908d1 / t26)
      t31 = t1 + 0.47584d-2
      t32 = t31 * t31
      t34 = log(t32 * t23)
      t37 = t20 * (t25 + 0.31770800474394146400d0 * t29 + 0.41403379428
     &206274608d-3 * t34)
      t38 = 0.1d1 + zeta
      t39 = (t38)**(1.0d0/3.0d0)
      t41 = 0.1d1 - zeta
      t42 = (t41)**(1.0d0/3.0d0)
      t44 = t39 * t38 + t42 * t41 - 0.2d1
      t45 = t37 * t44
      t46 = 2.0d0**(1.0d0/3.0d0)
      t47 = t46 - 0.1d1
      t49 = 0.1d1 / t47 / 0.2d1
      t50 = zeta * zeta
      t51 = t50 * t50
      t54 = 0.9d1 * t47
      t55 = t49 * (-t51 + 0.1d1) * t54
      t59 = rs + 0.706042d1 * t1 + 0.180578d2
      t60 = 0.1d1 / t59
      t62 = log(rs * t60)
      t64 = t8 + 0.706042d1
      t67 = atan(0.47309269095601128300d1 / t64)
      t69 = t1 + 0.32500d0
      t70 = t69 * t69
      t72 = log(t70 * t60)
      t74 = 0.1554535d-1 * t62 + 0.52491393169780936218d-1 * t67 + 0.22
     &478670955426118383d-2 * t72 - t7 - t13 - t18
      t75 = t74 * t44
      t76 = t49 * t51
      ecVWN5 = t7 + t13 + t18 - t45 * t55 / 0.24d2 + t75 * t76
      end function

C ... end of dftfunjt.F
