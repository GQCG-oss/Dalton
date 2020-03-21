C SOURCES:
C    Simone Paziani, Saverio Moroni, Paola Gori-Giorgi, and Giovanni B. Bachelet.
C    Local-spin-densityfunctional for multideterminant density functional theory.
C    Physical Review B, 73(15), apr 2006.



C*****************************************************************************
      pure subroutine ESRX_LDA_ERF_case_1(rho_a, mu, Ea)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
C   Generated: March 21, 2019
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
