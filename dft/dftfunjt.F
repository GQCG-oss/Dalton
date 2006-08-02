C****************************************************************************
      SUBROUTINE DESR(ESR,rho,grd,mu,norder,
     >                d1Edrho,d1Edgrd,d1Edgrda2,
     >                d2Edrhodrho,d2Edgrddgrd,d2Edrhodgrd)
C*****************************************************************************
C     Compute first-order and/or second-order numerical derivatives 
C     of a generic short-range functional ESR(rho,grd,mu,E)
C
C     Input: ESR    : functional routine
C            rho    : density
C            grd    : norm of gradient of density
C            mu     : Interation parameter
C            norder : order of derivatives
C
C     Ouput: d1Edrho       : first derivative wrt to rho
C            d1Edgrd       : first derivative wrt to grd  (for 2nd-order optimisation)
C            d1Edgrda2     : first derivative wrt to grda2 (for HF optimisation)
C            (grda2 is the square of gradient of spin-alpha density)
C            d2Edrhodrho   : second derivative wrt to rho and rho
C            d2Edgrddgrd   : second derivative wrt to grd and grd
C            d2Edrhodgrd   : second derivative wrt to rho and grd
C
C     Author: J. Toulouse
C     Date  : 16-02-05
C*****************************************************************************
      implicit none

      external ESR
      integer norder
      double precision rho, grd, mu
      double precision E, d1Edrho,d1Edgrd,d1Edgrda2
      double precision d2Edrhodrho,d2Edgrddgrd,d2Edrhodgrd

      double precision hr, hg, hr2, hg2, hrg
      double precision rhom, rhop, rhopp, grdm, grdp, grdpp

      double precision Em0, Ep0, Epp0
      double precision E0m, E0p, E0pp
      double precision Emm, Epp, Epm, Emp
      double precision d1Edgrd2, dsdgrd, dgrddgrd2

C     check
      if(norder .le. 0 .or. norder .ge. 3) then
       print *,'DESR: norder=',norder
       call quit('DESR: norder must be 1 or 2')
      endif

C     energy
      call ESR(rho,grd,mu,E)
!      print *,'D1SR: E=',E

C     (optimal?) numerical step
      hr = 1.d-4*rho
      hg = 1.d-4*grd

C     variations of variables
      rhom = rho-hr
      rhop = rho+hr
      grdm = grd-hg
      grdp = grd+hg

C First-order derivatives --------------------------------------------

C     derivative wrt rho ---------------------
!      if(rhom .ge. 0.d0) then 
!      central 2nd-order approximation
       call ESR(rhom,grd,mu,Em0)
       call ESR(rhop,grd,mu,Ep0)
       d1Edrho = (Ep0 - Em0)/(2.d0*hr)
!      else 
!       stop 'rhom < 0'
!!      forward 1st-order approximation
!       call ESR(rhop,grd,mu,Ep0)
!       d1Edrho = (Ep0 - E)/hr
!      endif
!      print *,'D1SR: d1Edrho=',d1Edrho

C     derivative wrt grd ---------------------
!      if(grdm .ge. 0.d0) then
!      central 2nd-order approximation
       call ESR(rho,grdm,mu,E0m)
       call ESR(rho,grdp,mu,E0p)
       d1Edgrd = (E0p - E0m)/(2.d0*hg)
!      else
!       stop 'grdm < 0'
!!      forward 1st-order approximation
!       call ESR(rho,grdp,mu,E0p)
!       d1Edgrd = (E0p - E)/hg
!      endif
!      print *,'D1SR: d1Edgrd=',d1Edgrd

C     derivative wrt grd2
      dgrddgrd2 = 1.d0/(2.d0*grd)
      d1Edgrd2 = d1Edgrd * dgrddgrd2

C     derivative wrt grda2
      d1Edgrda2 = 4.d0*d1Edgrd2

C Second-order derivatives --------------------------------------------
      if(norder .ge. 2) then

      hr2= hr*hr
      hg2= hg*hg
      hrg= hr*hg

C     derivative wrt rho and rho -----------------
!      if(rhom .ge. 0.d0) then 
!      central 2nd-order approximation
       d2Edrhodrho = (Ep0 - 2.d0*E + Em0)/hr2
!      else 
!       stop 'rhom < 0'
!!      forward 1st-order approximation
!       rhopp = rho + 2.d0*hr
!       call ESR(rhop,grd,mu,Ep0)
!       call ESR(rhopp,grd,mu,Epp0)
!       d2Edrhodrho = (Epp0 -2.d0*Ep0 + E)/hr2
!      endif

C     derivative wrt grd and grd -----------------------
!      if(grdm .ge. 0.d0) then
!      central 2nd-order approximation
       d2Edgrddgrd = (E0p - 2.d0*E + E0m)/hg2
!      else
!       stop 'grdm < 0'
!!      forward 1st-order approximation
!       grdpp = grd + 2.d0*hg
!       call ESR(rho,grdp,mu,E0p)
!       call ESR(rho,grd,mu,E0pp)
!       d2Edgrddgrd = (E0pp -2.d0*E0p + E)/hg2
!      endif

C     derivative wrt rho and grd -----------------------
!      if(rhom.ge. 0.d0 .and. grdm .ge. 0.d0) then
       call ESR(rhop,grdp,mu,Epp)
       call ESR(rhom,grdp,mu,Emp)
       call ESR(rhop,grdm,mu,Epm)
       call ESR(rhom,grdm,mu,Emm)
!      central 2nd-order approximation
       d2Edrhodgrd = (Epp - Emp - Epm + Emm)/(4.d0*hrg)
!      else
!       stop 'rhom < 0 or grdm < 0'
!       call ESR(rhop,grdp,mu,Epp)
!!      forward 1st-order approximation
!       d2Edrhodgrd = (Epp - E0p - Ep0 + E)/hrg
!      endif
     
      endif

      return
      end

C****************************************************************************
      subroutine ESRXPBETCSERF(rho,grd,mu,e)
C*****************************************************************************
C     Short-range PBE exchange functional from 
C       Toulouse, Colonna, Savin JCP 122, 014110 (2005)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interation parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 19-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision rho13
      double precision drho2

      double precision exerflda
      double precision exerfpbe
      double precision t1,t2,t3,t4

      double precision XKSR(3)
      logical ERFEXP

! function 
      double precision berf

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      ERFEXP = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      exerflda = XKSR(1)

      drho2 = grd*grd

!      Code generated by Mathematica
       exerfpbe = (-5.521381337364587693627982d1*exerflda*(1.63580000000
     &  000014281909d0*rho**4 + 4.272901476924512165211922d-2*drho2*berf
     &  (1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(4.d0/
     &  3.d0) + 1.d0*exerflda*rho**(8.d0/3.d0)))/(-9.0318755916609930523
     &  0091d1*rho**4 + 1.442249570307408301772512d0*drho2*exerflda*berf
     &  (1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)) - 5.5213813
     &  37364587693627982d1*exerflda*rho**(8.d0/3.d0))

        e = exerfpbe

      return
      end
C****************************************************************************
      subroutine ESRXPBEHSEERF(rho,grd,mu,e)
C*****************************************************************************
C     Short-range PBE exchange functional from 
C       Heyd and Scuseria, JCP 120, 7274 (2004)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interation parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision pi,rho13,kF,s
      double precision exLDA
      double precision FxHSE

C     constants
      pi = 3.141592653589793d0

C     LDA exchange
      rho13 = rho**(1.d0/3.d0)
      CALL EDRC(exLDA,rho,rho13)

C     Gradient
      kF = (3.d0*(pi**2)*rho)**(1.d0/3.d0)
      s = grd/(2.d0*kF*rho)

C     Enhancement factor (over LDA) of short-range PBE from HSE
      call wpbe_analytical_erfc_approx(rho,s,mu,FxHSE)

C     Energy
      e = exLDA*FxHSE
       
      return
      end
C****************************************************************************
      SUBROUTINE VSRXPBEERF(e,d1edrho,d1edgrd,d1edgrda2,d2edrhodrho,
     >                      d2edgrddgrd,d2edrhodgrd,rho,grd,mu,
     >                      norder,FxSR)
C*****************************************************************************
C     Short-range PBE exchange functional from 
C       Heyd and Scuseria, JCP 120, 7274 (2004)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interation parameter
C            norder : Order of derivative
C            FxSR   : routine for evaluation of Fx(short-range PBE)
C
C     Ouput: e             : energy
C            d1edrho       : first derivative wrt to rho
C            d1edgrd       : first derivative wrt to grd   (for 2nd-order optimisation)
C            d1edgrda2     : first derivative wrt to grda2 (for HF optimisation)
C            (grda2 is the square of gradient of spin-alpha density)
C            d2edrhodrho   : second derivative wrt to rho and rho
C            d2edgrddgrd   : second derivative wrt to grd and grd
C            d2edrhodgrd   : second derivative wrt to rho and grd
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C
C     Derivatives wrong. Not used.
C*****************************************************************************
      implicit none

      integer norder

      double precision rho, grd, mu
      double precision e, d1edrho, d1edgrd, d1edgrda2
      double precision d2edrhodrho,d2edgrddgrd,d2edrhodgrd

      double precision pi,rho13,kF,s
      double precision exLDA
      double precision FxHSE

      integer i
      double precision rho1, rho2, s1, s2
      double precision FxHSE1, FxHSE2

      double precision h, h2
      double precision rhom, rhop, rhopp, sm, sp, spp

      double precision d1exLDAdrho
      double precision FxHSEm0, FxHSEp0, FxHSEpp0
      double precision FxHSE0m, FxHSE0p, FxHSE0pp
      double precision FxHSEmm, FxHSEmp, FxHSEpm, FxHSEpp
      double precision d1FxHSEdrho
      double precision d1FxHSEds, d1FxHSEdgrd, d1FxHSEdgrd2
      double precision d1FxHSEdgrda2
      double precision dsdrho, dsdgrd, dgrddgrd2

      double precision d2exLDAdrhodrho
      double precision d2FxHSEdrhodrho
      double precision d2FxHSEdrhods, d2FxHSEdgrddgrd
      double precision d2FxHSEdsds, d2FxHSEdrhodgrd
      double precision dsdgrd2

      external FxSR

      double precision rhom13, rhop13, rhopp13
      double precision DF1000,DF0010,DF2000,DF1010,DF0020

!      print *,'VSRPBEXERF: mu=',mu
!      print *,'VSRPBEXERF: rho=',rho
!      print *,'VSRPBEXERF: grd=',grd

      call quit('routine is wrong. Use ESRXPBEERF instead')
C     constants
      pi = 3.141592653589793d0

!      print *,'VSRPBEXERF: mu=',mu
!      print *,'VSRPBEXERF: rho=',rho,' grd=',grd

C     LDA exchange
      rho13 = rho**(1.d0/3.d0)
      CALL EDRC(exLDA,rho,rho13)

C     Gradient
      kF = (3.d0*(pi**2)*rho)**(1.d0/3.d0)
      s = grd/(2.d0*kF*rho)

C     Enhancement factor (over LDA) of short-range PBE from HSE
      call FxSR(rho,s,mu,FxHSE)

C     Energy
      e = exLDA*FxHSE
!debug LDA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      FxHSE = exLDA
!      e = FxHSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      print *,'VSRXPBEERF: E=',E
       
C     First-order (numerical) derivatives ---------------------------------
      if(norder .ge. 1) then

!     LDA
      CALL VDRC(d1exLDAdrho,rho13)

      dsdgrd = 1.d0/(2.d0*kF*rho)

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C debug: test values for h
      if(.false.) then
      rho = 1.003d0
      s   = 1.000d0
      mu  = 1.000d0

      print *,'rho=',rho
      print *,'s=',s
      print *,'mu=',mu

      call FxSR(rho,s,mu,FxHSE)
      print *,'FxHSE=',FxHSE

      do i = 0, 20
      h = 10.d0**(-i)

      rho1 = rho-h
      rho2 = rho+h
      if(rho1 .ge. 0.d0) then
       call FxSR(rho1,s,mu,FxHSE1)
       call FxSR(rho2,s,mu,FxHSE2)
       d1FxHSEdrho = (FxHSE2 - FxHSE1)/(2.d0*h)
      else
       call FxSR(rho2,s,mu,FxHSE2)
       d1FxHSEdrho = (FxHSE2 - FxHSE)/h
      endif

C     derivative wrt s
      s1 = s-h
      s2 = s+h
      if(s1 .ge. 0.d0) then
       call FxSR(rho,s1,mu,FxHSE1)
       call FxSR(rho,s2,mu,FxHSE2)
       d1FxHSEds = (FxHSE2 - FxHSE1)/(2.d0*h)
      else
       call FxSR(rho,s2,mu,FxHSE2)
       d1FxHSEds = (FxHSE2 - FxHSE)/h
      endif

C     Print
      write(*,'(A,D10.5,A,F25.20,A,F25.20)') 'h=',h,
     >               ' d1FxHSEdrho=',d1FxHSEdrho,' d1FxHSEds=',d1FxHSEds
      enddo
      stop
      endif
C end debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C     numerical step
      h = 1.d-4

C     variations of variables
      rhom = rho-h
      rhop = rho+h
      sm = s-h
      sp = s+h
!debug LDA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      rhom13 = rhom**(1.d0/3.d0)
!      rhop13 = rhop**(1.d0/3.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C     derivative wrt rho ---------------------
      call FxSR(rhop,s,mu,FxHSEp0)
      if(rhom .ge. 0.d0) then 
!      central 2nd-order approximation
       call FxSR(rhom,s,mu,FxHSEm0)
       d1FxHSEdrho = (FxHSEp0 - FxHSEm0)/(2.d0*h)
      else 
!      forward 1st-order approximation
       d1FxHSEdrho = (FxHSEp0 - FxHSE)/h
      endif

      d1edrho = d1exLDAdrho*FxHSE + exLDA*d1FxHSEdrho
!      print *,'VSRXPBEERF: d1edrho=',d1edrho

C     derivative wrt grd ---------------------

C     derivative wrt s
      call FxSR(rho,sp,mu,FxHSE0p)
      if(sm .ge. 0.d0) then
!      central 2nd-order approximation
       call FxSR(rho,sm,mu,FxHSE0m)
       d1FxHSEds = (FxHSE0p - FxHSE0m)/(2.d0*h)
      else
!      forward 1st-order approximation
       d1FxHSEds = (FxHSE0p - FxHSE)/h
      endif

C     derivative wrt grd
      d1FxHSEdgrd = d1FxHSEds * dsdgrd 
      d1edgrd   = exLDA*d1FxHSEdgrd
!      print *,'VSRXPBERF: d1Edgrd=',d1Edgrd
      
C     derivative wrt grd2
      dgrddgrd2 = 1.d0/(2.d0*grd)
      d1FxHSEdgrd2 = d1FxHSEdgrd * dgrddgrd2

C     derivative wrt grda2
      d1FxHSEdgrda2 = 4.d0*d1FxHSEdgrd2
      d1edgrda2 = exLDA*d1FxHSEdgrda2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dsdrho = -grd*2.d0/(3.d0**(4.d0/3.d0)*
     >  pi**(2.d0/3.d0)*rho**(7.d0/3.d0))
      d1FxHSEdrho = d1FxHSEdrho + d1FxHSEds * dsdrho
      d1edrho = d1exLDAdrho*FxHSE + exLDA*d1FxHSEdrho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

!debug LDA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if(rhom .ge. 0.d0) then 
!!      central 2nd-order approximation
!       CALL EDRC(FxHSEm0,rhom,rhom13)
!       CALL EDRC(FxHSEp0,rhop,rhop13)
!       d1FxHSEdrho = (FxHSEp0 - FxHSEm0)/(2.d0*h)
!      else 
!!      forward 1st-order approximation
!       CALL EDRC(FxHSEp0,rhop,rhop13)
!       d1FxHSEdrho = (FxHSEp0 - FxHSE)/h
!      endif
!      d1edrho   = d1FxHSEdrho
!      d1edgrd   = 0.d0
!      d1edgrda2 = 0.d0
!!      print *,'VSRXPBERF: d1edrho=',d1edrho
!!      print *,'VSRXPBERF: d1exLDAdrho=',d1exLDAdrho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif

C     Second-order (numerical) derivatives ------------------------------
      if(norder .ge. 2) then
      
C     LDA
      call V1DRC(d1exLDAdrho,d2exLDAdrhodrho,rho,rho13)

      h2 = h*h

C     derivative wrt rho and rho -----------------
      if(rhom .ge. 0.d0) then 
!      central 2nd-order approximation
       d2FxHSEdrhodrho = (FxHSEp0 - 2.d0*FxHSE + FxHSEm0)/h2
      else 
!      forward 1st-order approximation
       rhopp = rho + 2.d0*h
       call FxSR(rhopp,s,mu,FxHSEpp0)
       d2FxHSEdrhodrho = (FxHSEpp0 -2.d0*FxHSEp0 + FxHSE)/h2
      endif

      d2edrhodrho = d2exLDAdrhodrho*FxHSE + 2.d0*d1exLDAdrho*d1FxHSEdrho
     >            + exLDA*d2FxHSEdrhodrho

C     derivative wrt grd and grd -----------------------

C     derivative wrt s and s
      if(sm .ge. 0.d0) then
!      central 2nd-order approximation
       d2FxHSEdsds = (FxHSE0p - 2.d0*FxHSE + FxHSE0m)/h2
      else
!      forward 1st-order approximation
       spp = s + 2.d0*h
       call FxSR(rho,spp,mu,FxHSE0pp)
       d2FxHSEdsds = (FxHSE0pp -2.d0*FxHSE0p + FxHSE)/h2
      endif

C     derivative wrt grd and grd
      dsdgrd2 = dsdgrd*dsdgrd
      d2FxHSEdgrddgrd = d2FxHSEdsds * dsdgrd2
      
      d2edgrddgrd = exLDA*d2FxHSEdgrddgrd

C     derivative wrt rho and grd -----------------------

C     derivative wrt rho and s
      if(rhom.ge. 0.d0 .and. sm .ge. 0.d0) then
       call FxSR(rhop,sp,mu,FxHSEpp)
       call FxSR(rhom,sp,mu,FxHSEmp)
       call FxSR(rhop,sm,mu,FxHSEpm)
       call FxSR(rhom,sm,mu,FxHSEmm)
!      central 2nd-order approximation
       d2FxHSEdrhods = (FxHSEpp - FxHSEmp - FxHSEpm + FxHSEmm)/(4.d0*h2)
      else
       call FxSR(rhop,sp,mu,FxHSEpp)
!      forward 1st-order approximation
       d2FxHSEdrhods = (FxHSEpp - FxHSE0p - FxHSEp0 + FxHSE)/h2
      endif

      d2FxHSEdrhodgrd = d2FxHSEdrhods * dsdgrd 

      d2edrhodgrd = d1exLDAdrho*d1FxHSEdgrd + exLDA*d2FxHSEdrhodgrd


      Call V1BCK(DF1000,DF0010,DF2000,DF1010,DF0020,
     &                 RHO,GRD)
      DF1000 = DF1000 + d1exLDAdrho
      DF2000 = DF2000 + d2exLDAdrhodrho

      print *, 'rho=',rho,'grd=',grd
      print *, 'DF1000=',DF1000,d1Edrho,DF1000-d1Edrho
      print *, 'DF0010=',DF0010,d1Edgrd,DF0010-d1Edgrd
      print *, 'DF2000=',DF2000,d2Edrhodrho,DF2000-d2Edrhodrho
      print *, 'DF0020=',DF0020,d2Edgrddgrd,DF0020-d2Edgrddgrd
      print *, 'DF1010=',DF1010,d2Edrhodgrd,DF1010-d2Edrhodgrd

!debug LDA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if(rhom .ge. 0.d0) then 
!!      central 2nd-order approximation
!       d2FxHSEdrhodrho = (FxHSEp0 - 2.d0*FxHSE + FxHSEm0)/h2
!      else 
!!      forward 1st-order approximation
!       rhopp = rho + 2.d0*h
!       rhopp13 = rhopp**(1.d0/3.d0)
!       CALL EDRC(FxHSEpp0,rhopp,rhopp13)
!       d2FxHSEdrhodrho = (FxHSEpp0 -2.d0*FxHSEp0 + FxHSE)/h2
!      endif
!      d2edrhodrho   = d2FxHSEdrhodrho
!      d2edgrddgrd   = 0.d0
!      d2edrhodgrd   = 0.d0
!!      print *,'VSRXPBERF: d2edrhodrho=',d2edrhodrho
!!      print *,'VSRXPBERF: d2exLDAdrhodrho=',d2exLDAdrhodrho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif

      RETURN
      END
C------------------------------------------------------------------
      Subroutine FxB(Rho,s,Omega,Fx)
C     for debug
C------------------------------------------------------------------
      implicit none
      
      double precision Rho,s,Omega,Fx
      double precision EXLDA, BECKI, RHO13, kF, RHOGRD, pi

C     constants
      pi = 3.141592653589793d0

      RHO13 = RHO**(1.d0/3.d0)
      kF = (3.d0*(pi**2)*rho)**(1.d0/3.d0)
      RHOGRD = s*(2.d0*kF*RHO)

      call EDRC(EXLDA,RHO,RHO13) 
      call EBCK(BECKI,RHO,RHO13,RHOGRD)
   
      Fx = 1.d0 + BECKI/EXLDA

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
      subroutine ECPBE(rho,grd,e)
C*****************************************************************************
C     PBE correlation functional
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision pi, pi2
      double precision rho13, eLDA, epsLDA
      double precision gam, bet
      double precision kF, ks, t, t2, t4
      double precision A, A2
      double precision frac, H

C     General constants
      pi = 3.141592653589793d0
      pi2 = 9.869604401089358d0

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      call EVWN(eLDA,rho,rho13)
      epsLDA = eLDA/rho

      gam = 0.031090690869654897d0
      bet = 0.066725

      kF = (3.d0*pi2*rho)**(1.d0/3.d0)
      ks = dsqrt(4.d0*kF/pi) 
      t = dabs(grd)/(2.d0*ks*rho)
      t2 = t*t
      t4 = t2*t2

      A = (bet/gam)/(dexp(-epsLDA/gam) - 1.d0)
      A2 = A*A

      frac = (1.d0 + A*t2)/(1.d0 + A*t2 + A2*t4)
      H = gam*dlog(1.d0+(bet/gam)*t2*frac)

      e = rho*(epsLDA + H)

      return
      end

C****************************************************************************
      subroutine ECPBE_old(rho,grd,e)
C*****************************************************************************
C     PBE correlation functional (old version)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision rho13
      double precision exLDA
      double precision ecLDA

      double precision drho2
      double precision t1, bPBE

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      call EDRC(exLDA,rho,rho13)
      call EVWN(ecLDA,rho,rho13)

!     square of density gradient 
      drho2 = grd*grd

      bPBE = 0.21951d0

!     Code generated by Mathematica
C     Energy
       t1 = 9.869604401089358618834491d0*ecLDA + 3.068528194400546905
     & 827679d-1*rho*dlog(1.d0 - (8.401605835889726292406297d-1*drho2*
     &  exLDA*bPBE*
     &  rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*drho2*
     &  exLDA*bPBE
     &  *rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009
     &  d1*ecLDA)/rho))))/(1.d0 + (7.058698062165630644545894d-1*drho
     &  2**2*exLDA**2*bPBE
     &  **2*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.2163968442914
     &  82112072009d1*ecLDA)/rho))**2 - (8.401605835889726292406297d-
     &  1*drho2*exLDA*bPBE
     &  *rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.2163968442914821
     &  12072009d1*ecLDA)/rho))))
       e = 1.013211836423377714438795d-1*t1

      return
      end

C****************************************************************************
      subroutine ESRCRIERF(rho,grd,mu,e)
C*****************************************************************************
C     rational interpolation for short-range erf correlation functional
C     Toulouse, Colonna, Savin PRA 70, 062505 (2004)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision mu2, rs, pi, sqrtpi, sqrt2
      double precision d, a, bet, gam
      double precision grs32, grs12, expgrs12, g0
      double precision epbe, epspbe
      double precision d1, d2, denom
      
!      print *,'ESRCRIERF: mu=',mu
!      print *,'ESRCRIERF: rho=',rho,' grd=',grd
     
C     General constants
      pi = 3.141592653589793d0
      sqrtpi = dsqrt(pi)
      sqrt2 = sqrt(2.0d0)

C     Constants for g0
      d = 32.D0/(3.D0*pi)
      a = 3.2581D0
      bet = 163.44D0
      gam = 4.7125D0

!     Variables
      rs = (3.d0/(4.d0*pi*rho))**(1.d0/3.d0)
      mu2 =mu*mu

!     PBE correlation energy
      call ECPBE(rho,grd,epbe)
      epspbe = epbe/rho
!      print *,'ESRCRIERF: epbe=',epbe

!     g0 from Burke, Perdew, Ernzerhof
      grs32 = (gam+rs)**(3.d0/2.d0)
      grs12 = dsqrt(gam+rs)
      expgrs12 = dexp(-a*grs12)
      g0 = d*(grs32+bet)*expgrs12

      d2 = 2.d0*epspbe/(pi*rho*(g0 - 0.5d0))

C    modified d1 coefficient from Gori-Giorgi and Savin (submitted)
C    added by Emmanuel Fromager 13.12.05    
      d1 = -2.d0*(d2**2)*sqrt2*sqrtpi*rho*g0/(3.d0*epspbe)
C      d1 = -2.d0*(d2**2)*sqrtpi*rho*g0/(3.d0*epspbe)
      denom = 1.d0 + d1*mu + d2*mu2
      e = epbe/denom

!      print *,'ESRCRIERF: e=',e
      return
      end
C
C  begin Manu
C****************************************************************************
      subroutine ESRXRIERF(rho,grd,mu,e)
C*****************************************************************************
C     rational interpolation for short-range erf exchange functional
C     using the two first terms of both developments at small and large mu
C     (eqs. 17 and 33 in Ref : Toulouse, Colonna, Savin PRA 70, 062505 (2004)) 
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: Emmanuel Fromager (Southern Denmark University) 
C     Date  : 31-01-2006
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision mu2, mu6, pi, sqrtpi 
      double precision three53, pi13, pi43, zero 
      double precision epbe
      double precision d, d2, denom1, denom2
      
!      print *,'ESRCRIERF: mu=',mu
!      print *,'ESRCRIERF: rho=',rho,' grd=',grd
     
C     General constants
      pi = 3.141592653589793d0
      sqrtpi = dsqrt(pi)
      three53 = 3.d0**(5.d0/3.d0) 
      pi13 = pi**(1.d0/3.d0)
      pi43 = pi**(4.d0/3.d0) 

!     Variables
      zero = 0.d0 
      mu2 =mu*mu
      mu6 =mu2*mu2*mu2

      d = -5.d0*(rho**(4.d0/3.d0))/(pi13*three53)      
      d2 = 20.d0/(pi43*three53*(rho**(2.d0/3.d0)))
      denom1 = 1.d0 + d2*mu2
      denom2 = 1.d0 + mu6

C Get standard PBE exchange functional calling HSEPBE with mu=0
      call ESRXPBEHSEERF(rho,grd,zero,epbe)

      e = (d/denom1) + ((epbe + (rho*mu/sqrtpi) - d )/denom2)  
C      epspbe = epbe/rho
!      print *,'ESRCRIERF: epbe=',epbe
!      print *,'ESRCRIERF: e=',e
      return
      end
C
C end Manu
C
C****************************************************************************
      subroutine ESRCWIERF(rho,grd,mu,e)
C*****************************************************************************
C     weighted interpolation for for short-range erf correlation functional
C     Toulouse, Colonna, Savin PRA 70, 062505 (2004)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e, derfc

      logical ERFEXP, MULOCAL
      double precision mu0,rho13, rs, pi
      double precision elda,epbe,eldamu
      double precision ZKSR(3)
C     
!      print *,'ESRCWIERF: mu=',mu
!      print *,'ESRCWIERF: rho=',rho,' grd=',grd

      pi = 3.141592653589793d0
      ERFEXP = .false.
      MULOCAL = .false.
      rho13 = rho**(1.d0/3.d0)
      rs = (3.d0/(4.d0*pi*rho))**(1.d0/3.d0)

!     LDA energy
      mu0 = 0.d0
      call VCSRLDA(ZKSR,rho,rho13,mu0,MULOCAL,0,ERFEXP)
      elda = ZKSR(1)
!      print *,'ESRCWIERF: elda=',elda

!     PBE energy
      call ECPBE(rho,grd,epbe)
!      print *,'ESRCWIERF: epbe=',epbe

!     short-range LDA energy
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      eldamu = ZKSR(1)
!      print *,'ESRCWIERF: eldamu=',eldamu

      e = (epbe - elda)*derfc(mu*rs) + eldamu

!      print *,'ESRCWIERF: e=',e
      return
      end
C****************************************************************************
      subroutine ESRCPBETCSERF(rho,grd,mu,e)
C*****************************************************************************
C     Short-range PBE correlation functional
C       Toulouse, Colonna, Savin JCP 122, 014110 (2005)
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision rho13

      double precision drho2
      double precision exerflda
      double precision ecerflda
      double precision ecerfpbe
      double precision t1

      logical ERFEXP, MULOCAL
      double precision XKSR(3), ZKSR(3)

! function 
      double precision berf

!      print *,'ESRCPBETCSERF: mu=',mu
!      print *,'ESRCPBETCSERF: rho=',rho,' grd=',grd

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      ERFEXP = .false.
      MULOCAL = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     square of density gradient 
      drho2 = grd*grd

!     Code generated by Mathematica
C     Energy
       t1 = 9.869604401089358618834491d0*ecerflda + 3.068528194400546905
     & 827679d-1*rho*dlog(1.d0 - (8.401605835889726292406297d-1*drho2*ex
     &  erflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &  rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*drho2*
     &  exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009
     &  d1*ecerflda)/rho))))/(1.d0 + (7.058698062165630644545894d-1*drho
     &  2**2*exerflda**2*berf(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))**2*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.2163968442914
     &  82112072009d1*ecerflda)/rho))**2 - (8.401605835889726292406297d-
     &  1*drho2*exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.2163968442914821
     &  12072009d1*ecerflda)/rho))))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      e = ecerfpbe
!      print *,'ESRCPBETCSERF: e=',e

      return
      end
C****************************************************************************
      subroutine ESRCPBETCSJERF(rho,grd,mu,e)
C*****************************************************************************
C     Short-range PBE correlation functional
C       Toulouse, Colonna, Savin JCP 122, 014110 (2005)
C       with modification by Jensen
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interaction parameter
C
C     Ouput: e             : energy
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      double precision rho, grd, mu
      double precision e

      double precision rho13

      double precision drho2
      double precision exerflda
      double precision ecerflda
      double precision ecerfpbe
      double precision t1

      logical ERFEXP, MULOCAL
      double precision XKSR(3), ZKSR(3)

      double precision epbe

! function 
      double precision berf

!      print *,'ESRCPBETCSERF: mu=',mu
!      print *,'ESRCPBETCSERF: rho=',rho,' grd=',grd

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      ERFEXP = .false.
      MULOCAL = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     square of density gradient 
      drho2 = grd*grd

!     Code generated by Mathematica
C     Energy
       t1 = 9.869604401089358618834491d0*ecerflda + 3.068528194400546905
     & 827679d-1*rho*dlog(1.d0 - (8.401605835889726292406297d-1*drho2*ex
     &  erflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &  rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*drho2*
     &  exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009
     &  d1*ecerflda)/rho))))/(1.d0 + (7.058698062165630644545894d-1*drho
     &  2**2*exerflda**2*berf(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))**2*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.2163968442914
     &  82112072009d1*ecerflda)/rho))**2 - (8.401605835889726292406297d-
     &  1*drho2*exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.2163968442914821
     &  12072009d1*ecerflda)/rho))))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      call ECPBE(rho,grd,epbe)

! modification of ESRCPBETCS:
      e = max(ecerfpbe,epbe)
!      print *,'ESRCPBETCSERF: e=',e

      return
      end
C****************************************************************************
      SUBROUTINE VSRCPBEERF(e,d1edrho,d1edgrd,d1edgrda2,d2edrhodrho,
     >                      d2edgrddgrd,d2edrhodgrd,rho,grd,mu,
     >                      norder)
C*****************************************************************************
C     Short-range PBE correlation functional from 
C
C     Input: rho    : density
C            grd    : norm of gradient of density
C            mu     : Interation parameter
C            norder : Order of derivative
C
C     Ouput: e             : energy
C            d1edrho       : first derivative wrt to rho
C            d1edgrd       : first derivative wrt to grd   (for 2nd-order optimisation)
C            d1edgrda2     : first derivative wrt to grda2 (for HF optimisation)
C            (grda2 is the square of gradient of spin-alpha density)
C            d2edrhodrho   : second derivative wrt to rho and rho
C            d2edgrddgrd   : second derivative wrt to grd and grd
C            d2edrhodgrd   : second derivative wrt to rho and grd
C
C     Author: J. Toulouse
C     Date  : 12-02-05
C*****************************************************************************
      implicit none

      integer norder
      double precision rho, grd, mu
      double precision e, d1edrho, d1edgrd, d1edgrda2
      double precision d2edrhodrho, d2edgrddgrd, d2edrhodgrd

      double precision rho13
      double precision exLDA, d1exLDAdrho
      double precision ecLDA, d1ecLDAdrho

      double precision drho2
      double precision exerflda,dexerfldadrho
      double precision ecerflda,decerfldadrho
      double precision ecerfpbe,decerfpbedrho,decerfpbeddrho2 
      double precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      double precision t11,t12,t13,t14,t15,t16,t17,t18,t19

      logical ERFEXP, MULOCAL
      double precision XKSR(3), ZKSR(3)

! function 
      double precision berf
      double precision dberfda

!      print *,'VSRPBECERF: mu=',mu
!      print *,'VSRPBECERF: rho=',rho,' grd=',grd

      if(rho .lt. 1.d-10) then
       return
      endif

C     LDA energy
      rho13 = rho**(1.d0/3.d0)
      ERFEXP = .false.
      MULOCAL = .false.
      call VXSRLDA(XKSR,rho,rho13,mu,0,ERFEXP)
      call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,0,ERFEXP)
      exerflda = XKSR(1)
      ecerflda = ZKSR(1)

!     square of density gradient 
      drho2 = grd*grd

!     Code generated by Mathematica
C     Energy
       t1 = 9.869604401089358618834491d0*ecerflda + 3.068528194400546905
     & 827679d-1*rho*dlog(1.d0 - (8.401605835889726292406297d-1*drho2*ex
     &  erflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*
     &  rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*drho2*
     &  exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009
     &  d1*ecerflda)/rho))))/(1.d0 + (7.058698062165630644545894d-1*drho
     &  2**2*exerflda**2*berf(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))**2*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.2163968442914
     &  82112072009d1*ecerflda)/rho))**2 - (8.401605835889726292406297d-
     &  1*drho2*exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.
     &  d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.2163968442914821
     &  12072009d1*ecerflda)/rho))))
       ecerfpbe = 1.013211836423377714438795d-1*t1

      e = ecerfpbe
!      print *,'VSRCPBERF: e=',e

C Derivative
      if (norder .ge. 1) then
            
!      LDA energy density derivative
       call VXSRLDA(XKSR,rho,rho13,mu,1,ERFEXP)
       call VCSRLDA(ZKSR,rho,rho13,mu,MULOCAL,1,ERFEXP)
       dexerfldadrho = XKSR(2)
       decerfldadrho = ZKSR(2)

!      Code generated by Mathematica
C      Derivative wrt rho
       t1 = ecerflda/rho
       t2 = -(ecerflda/rho**2)
       t3 = decerfldadrho/rho
       t4 = 1/(1.d0 - (8.401605835889726292406297d-1*drho2*exerflda*berf
     &  (1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d
     &  1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*drho2*exerflda*be
     &  rf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.
     &  1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda
     &  )/rho))))/(1.d0 + (7.058698062165630644545894d-1*drho2**2*exerfl
     &  da**2*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**
     &  2*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009
     &  d1*ecerflda)/rho))**2 - (8.401605835889726292406297d-1*drho2*exe
     &  rflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*r
     &  ho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*
     &  ecerflda)/rho))))
       t5 = (4.526237990654068686941552d-2*drho2*exerflda*mu*dberfda(1.6
     &  16204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*(1.d0 - (8.4016
     &  05835889726292406297d-1*drho2*exerflda*berf(1.616204596739954813
     &  316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dex
     &  p((-3.216396844291482112072009d1*ecerflda)/rho))))/(rho**5*(1.d0
     &   + (7.058698062165630644545894d-1*drho2**2*exerflda**2*berf(1.61
     &  6204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*rho**(-2.2d1/
     &  3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rh
     &  o))**2 - (8.401605835889726292406297d-1*drho2*exerflda*berf(1.61
     &  6204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d
     &  0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))
     &  ))
       t6 = (3.080588806492899640548976d0*drho2*exerflda*berf(1.61620459
     &  6739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.4d1/3.d0)*(1.
     &  d0 - (8.401605835889726292406297d-1*drho2*exerflda*berf(1.616204
     &  596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d0))/
     &  (-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))))/(
     &  1.d0 + (7.058698062165630644545894d-1*drho2**2*exerflda**2*berf(
     &  1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*rho**(-2.
     &  2d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda
     &  )/rho))**2 - (8.401605835889726292406297d-1*drho2*exerflda*berf(
     &  1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1
     &  /3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/r
     &  ho)))
       t7 = (-8.401605835889726292406297d-1*dexerfldadrho*drho2*berf(1.6
     &  16204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.
     &  d0)*(1.d0 - (8.401605835889726292406297d-1*drho2*exerflda*berf(1
     &  .616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/
     &  3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rh
     &  o))))/(1.d0 + (7.058698062165630644545894d-1*drho2**2*exerflda**
     &  2*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*rh
     &  o**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*e
     &  cerflda)/rho))**2 - (8.401605835889726292406297d-1*drho2*exerfld
     &  a*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**
     &  (-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecer
     &  flda)/rho)))
       t8 = -8.401605835889726292406297d-1*drho2*exerflda
       t9 = berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho
     &  **(-1.1d1/3.d0)
       t10 = 1/(1.d0 + (7.058698062165630644545894d-1*drho2**2*exerflda*
     &  *2*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*r
     &  ho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*
     &  ecerflda)/rho))**2 - (8.401605835889726292406297d-1*drho2*exerfl
     &  da*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho*
     &  *(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ece
     &  rflda)/rho)))
       t11 = (4.526237990654068686941552d-2*drho2*exerflda*mu*dberfda(1.
     &  616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)))/(rho**5*(-1.d
     &  0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))) + (3.08
     &  0588806492899640548976d0*drho2*exerflda*berf(1.61620459673995481
     &  3316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.4d1/3.d0))/(-1.d0 + de
     &  xp((-3.216396844291482112072009d1*ecerflda)/rho)) - (8.401605835
     &  889726292406297d-1*dexerfldadrho*drho2*berf(1.616204596739954813
     &  316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dex
     &  p((-3.216396844291482112072009d1*ecerflda)/rho)) + (8.4016058358
     &  89726292406297d-1*drho2*exerflda*((3.216396844291482112072009d1*
     &  ecerflda)/rho**2 - (3.216396844291482112072009d1*decerfldadrho)/
     &  rho)*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*de
     &  xp((-3.216396844291482112072009d1*ecerflda)/rho)*rho**(-1.1d1/3.
     &  d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho)
     &  )**2
       t12 = 8.401605835889726292406297d-1*drho2*exerflda*berf(1.6162045
     &  96739954813316614d-1*mu*rho**(-1.d0/3.d0))
       t13 = rho**(-1.1d1/3.d0)*(1.d0 - (8.401605835889726292406297d-1*d
     &  rho2*exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/
     &  3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.2163968442914821120
     &  72009d1*ecerflda)/rho)))
       t14 = (1.d0 + (7.058698062165630644545894d-1*drho2**2*exerflda**2
     &  *berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*rho
     &  **(-2.2d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ec
     &  erflda)/rho))**2 - (8.401605835889726292406297d-1*drho2*exerflda
     &  *berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(
     &  -1.1d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerf
     &  lda)/rho)))**(-2)
       t15 = (4.526237990654068686941552d-2*drho2*exerflda*mu*dberfda(1.
     &  616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)))/(rho**5*(-1.d
     &  0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))) - (7.60
     &  5533503381002378534421d-2*drho2**2*exerflda**2*mu*berf(1.6162045
     &  96739954813316614d-1*mu*rho**(-1.d0/3.d0))*dberfda(1.61620459673
     &  9954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-2.6d1/3.d0))/(-1.d
     &  0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho))**2
       t16 = (-5.176378578921462472666989d0*drho2**2*exerflda**2*berf(1.
     &  616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*rho**(-2.5d
     &  1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/
     &  rho))**2 + (1.411739612433126128909179d0*dexerfldadrho*drho2**2*
     &  exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )**2*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112072
     &  009d1*ecerflda)/rho))**2
       t17 = (-1.411739612433126128909179d0*drho2**2*exerflda**2*((3.216
     &  396844291482112072009d1*ecerflda)/rho**2 - (3.216396844291482112
     &  072009d1*decerfldadrho)/rho)*berf(1.616204596739954813316614d-1*
     &  mu*rho**(-1.d0/3.d0))**2*dexp((-3.216396844291482112072009d1*ece
     &  rflda)/rho)*rho**(-2.2d1/3.d0))/(-1.d0 + dexp((-3.21639684429148
     &  2112072009d1*ecerflda)/rho))**3 + (3.080588806492899640548976d0*
     &  drho2*exerflda*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0
     &  /3.d0))*rho**(-1.4d1/3.d0))/(-1.d0 + dexp((-3.216396844291482112
     &  072009d1*ecerflda)/rho))
       t18 = (-8.401605835889726292406297d-1*dexerfldadrho*drho2*berf(1.
     &  616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3
     &  .d0))/(-1.d0 + dexp((-3.216396844291482112072009d1*ecerflda)/rho
     &  )) + (8.401605835889726292406297d-1*drho2*exerflda*((3.216396844
     &  291482112072009d1*ecerflda)/rho**2 - (3.216396844291482112072009
     &  d1*decerfldadrho)/rho)*berf(1.616204596739954813316614d-1*mu*rho
     &  **(-1.d0/3.d0))*dexp((-3.216396844291482112072009d1*ecerflda)/rh
     &  o)*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.21639684429148211207200
     &  9d1*ecerflda)/rho))**2
       t19 = 3.109069086965489503494086d-2*log(1.d0 - (8.401605835889726
     &  292406297d-1*drho2*exerflda*berf(1.616204596739954813316614d-1*m
     &  u*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d0)*(1.d0 - (8.4016058358897
     &  26292406297d-1*drho2*exerflda*berf(1.616204596739954813316614d-1
     &  *mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dexp((-3.216
     &  396844291482112072009d1*ecerflda)/rho))))/(1.d0 + (7.05869806216
     &  5630644545894d-1*drho2**2*exerflda**2*berf(1.6162045967399548133
     &  16614d-1*mu*rho**(-1.d0/3.d0))**2*rho**(-2.2d1/3.d0))/(-1.d0 + d
     &  exp((-3.216396844291482112072009d1*ecerflda)/rho))**2 - (8.40160
     &  5835889726292406297d-1*drho2*exerflda*berf(1.6162045967399548133
     &  16614d-1*mu*rho**(-1.d0/3.d0))*rho**(-1.1d1/3.d0))/(-1.d0 + dexp
     &  ((-3.216396844291482112072009d1*ecerflda)/rho))))
       decerfpbedrho = t1 + rho*(t2 + t3 + 3.109069086965489503494086d-2
     &  *t4*(t5 + t6 + t7 + t8*t9*t10*t11 + t12*t13*t14*(t15 + t16 + t17
     &   + t18))) + t19

C       Derivative wrt drho2
       t1 = -2.608295326454957778144342d0*exerflda*berf(1.61620459673995
     &  4813316614d-1*mu*rho**(-1.d0/3.d0))*(-1.d0 + dexp((3.21639684429
     &  1482112072009d1*ecerflda)/rho))**3
       t2 = (3.093667726280135930968945d0*drho2*exerflda*berf(1.61620459
     &  6739954813316614d-1*mu*rho**(-1.d0/3.d0))*dexp((3.21639684429148
     &  2112072009d1*ecerflda)/rho) + 1.841116916640328143496607d0*(-1.d
     &  0 + dexp((3.216396844291482112072009d1*ecerflda)/rho))*rho**(1.1
     &  d1/3.d0))/(9.570780000627306053141656d0*drho2**2*exerflda**2*ber
     &  f(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0))**2*dexp((3
     &  .216396844291482112072009d1*ecerflda)/rho) + 1.13916079706371570
     &  5810133d1*drho2*exerflda*berf(1.616204596739954813316614d-1*mu*r
     &  ho**(-1.d0/3.d0))*(-1.d0 + dexp((3.216396844291482112072009d1*ec
     &  erflda)/rho))*rho**(1.1d1/3.d0) + 1.355884600295675603989991d1*(
     &  -1.d0 + dexp((3.216396844291482112072009d1*ecerflda)/rho))**2*rh
     &  o**(2.2d1/3.d0))
       t3 = rho**(2.5d1/3.d0)/(9.570780000627306053141656d0*drho2**2*exe
     &  rflda**2*berf(1.616204596739954813316614d-1*mu*rho**(-1.d0/3.d0)
     &  )**2*dexp((6.432793688582964224144018d1*ecerflda)/rho) + 1.13916
     &  0797063715705810133d1*drho2*exerflda*berf(1.61620459673995481331
     &  6614d-1*mu*rho**(-1.d0/3.d0))*(-1.d0 + dexp((3.21639684429148211
     &  2072009d1*ecerflda)/rho))*dexp((3.216396844291482112072009d1*ece
     &  rflda)/rho)*rho**(1.1d1/3.d0) + 1.355884600295675603989991d1*(-1
     &  .d0 + dexp((3.216396844291482112072009d1*ecerflda)/rho))**2*rho*
     &  *(2.2d1/3.d0))
       decerfpbeddrho2 = t1*t2*t3

C      derivatives
       d1edrho  =  decerfpbedrho
       d1edgrda2  =  4.d0*decerfpbeddrho2

C     d1edgrd to do !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C     End of derivative
      end if

C     Second derivative
      if(norder .ge. 2) then
       CALL QUIT ('Second derivative of VSRCPBEERF not yet implemented')
      endif

      RETURN
      END
!-------------------------------------------
      function berf(a)
!-------------------------------------------
!  Second-order exchange gradient expansion coefficient for erf 
!  interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
 
      double precision a
      double precision berf, derf

      double precision pi

C     constants
      pi = 3.141592653589793d0

      if(a .lt. 0.02d0) then
!      expansion for small mu to avoid numerical problems
       berf = 7.d0/81.d0 + 56.d0*dsqrt(pi)*a/243.d0 
     >        + (-128.d0/81.d0 + 448.d0*pi/729.d0)*a*a

      else

!      Code generated by Mathematica
       berf = (1.851851851851851851851852d-2*(-1.d0 + 1.44d2*a**4*(-1.d0
     &   + dexp(2.5d-1/a**2)) - 2.d0*a**2*(1.1d1 + 7.d0*dexp(2.5d-1/a**2
     &  ))))/(a**2*(3.2d1*a**4*(-1.d0 + dexp(2.5d-1/a**2)) - 3.d0*dexp(2
     &  .5d-1/a**2) + 1.417963080724412821838534d1*a*derf(5.d-1/a)*dexp
     &  (2.5d-1/a**2) - 8.d0*a**2*(-2.d0 + 3.d0*dexp(2.5d-1/a**2))))

      endif

      return
      end

!-------------------------------------------
      function dberfda(a)
!-------------------------------------------
!  Derivative of second-order exchange gradient
!  expansion coefficient for erf interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none
 
      double precision a
      double precision dberfda, derf
      double precision t1,t2,t3,t4,t5
      double precision pi

C     constants
      pi = 3.141592653589793d0

      if(a .lt. 0.02d0) then
!      expansion for small mu to avoid numerical problems
       dberfda = 56.d0*dsqrt(pi)/243.d0
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
       dberfda = t1 + t2*t3*t4 + t5

      endif

      return
      end

