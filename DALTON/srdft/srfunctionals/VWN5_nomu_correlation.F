C SOURCES:
C    S. H. Vosko, L. Wilk, and M. Nusair.
C    Accurate spin-dependent electron liquid correlation energies for local spin density calculations: a critical analysis.
C    Canadian Journal of Physics, 58(8):1200-1211, aug 1980.



C*****************************************************************************
      pure subroutine ESRC_VWN5(rho_c, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2
      E = 0.0d0
      x0 = 0.826307487110758d0*rho_c**(-0.166666666666667d0)
      x1 = 0.6203504908994d0*rho_c**(-0.333333333333333d0)
      x2 = 1d0/(3.5529372610885d0*x0 + x1 + 12.9352d0)
      E = rho_c*(0.0310907d0*log(x1*x2) + 0.000969022771154437d0*log(x2*
     &( 0.95318429299693652d0*x0 + 0.10498d0)**2) + 0.038783294878113d0*
     & atan(6.1519908197590798d0/(1.906368585993873d0*x0 + 3.72744000000
     &00001d0)))
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_VWN5(rho_c, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 0.826307487110758d0
      x1 = rho_c**(-0.166666666666667d0)*x0
      x2 = 1.90636858599387d0*x1 + 3.72744d0
      x3 = 0.682784063255296d0
      x4 = rho_c**0.333333333333333d0
      x5 = 0.90856029641607d0/x4
      x6 = x3*x5
      x7 = 1d0/(3.5529372610885d0*x1 + x6 + 12.9352d0)
      x8 = 0.953184292996937d0*x1 + 0.10498d0
      x9 = 0.0310907d0*log(x6*x7) + 0.000969022771154437d0*log(x7*x8**2)
     & + 0.038783294878113d0*atan(6.1519908197590798d0/x2)
      x10 = x2**(-2)
      x11 = rho_c**(-1.16666666666667d0)*x0
      x12 = 0.30285343213869d0*rho_c**(-1.33333333333333d0)
      x13 = x7*(0.592156210181417d0*x11 + x12*x3)
      E = rho_c*x9
      d1E(1) = -rho_c*(-0.0758081683534927d0*x10*x11/(37.8469910464d0*x1
     &0 + 1.0d0) + 0.0342197431724027d0*x4*(x12 - x13*x5) + (0.000307885
     &761673592d0* x11 - 0.000969022771154437d0*x13*x8)/x8) + x9
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_VWN5(rho_c, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 0.826307487110758d0
      x1 = rho_c**(-0.166666666666667d0)*x0
      x2 = 1.90636858599387d0*x1 + 3.72744d0
      x3 = 0.682784063255296d0
      x4 = rho_c**0.333333333333333d0
      x5 = 0.90856029641607d0/x4
      x6 = x3*x5
      x7 = 3.5529372610885d0*x1 + x6 + 12.9352d0
      x8 = 1d0/x7
      x9 = 0.953184292996937d0*x1 + 0.10498d0
      x10 = x9**2
      x11 = x10*x8
      x12 = 0.000969022771154437d0*log(x11) + 0.0310907d0*log(x6*x8) + 0
     &.038783294878113d0*atan(6.1519908197590798d0/x2)
      x13 = rho_c**(-1.16666666666667d0)*x0
      x14 = x2**(-2)
      x15 = 37.8469910464d0*x14 + 1.0d0
      x16 = 1d0/x15
      x17 = x14*x16
      x18 = x13*x17
      x19 = rho_c**(-1.33333333333333d0)
      x20 = 0.30285343213869d0*x19
      x21 = 0.592156210181417d0*x13 + x20*x3
      x22 = x21*x8
      x23 = x20 - x22*x5
      x24 = 0.0342197431724027d0*x4
      x25 = x23*x24
      x26 = 1d0/x9
      x27 = 0.000307885761673592d0*x13
      x28 = x22*x9
      x29 = 0.000615771523347183d0*x13
      x30 = rho_c**(-2.33333333333333d0)*x3
      x31 = rho_c**(-2.16666666666667d0)*x0
      x32 = 1d0/x10
      x33 = 0.317728097665645d0*x13 - x28
      x34 = 0.40380457618492d0*rho_c**(-2.33333333333333d0)
      x35 = 0.60570686427738d0*x19
      x36 = x3*x34 + 0.69084891187832d0*x31
      x37 = x21*(1.18431242036283d0*x13 + x3*x35)/x7**2
      E = rho_c*x12
      d1E(1) = -rho_c*(-0.0758081683534927d0*x18 + x25 + x26*(x27 - 0.00
     &0969022771154437d0*x28)) + x12
      d2E(1) = rho_c*(-0.0114065810574676d0*rho_c**(-0.666666666666667d0
     &)*x23 + 0.0481727702369445d0*x16*x30/x2**3 - 0.0884428630790749d0*
     &x17*x31 + x22*x25 + 0.000969022771154437d0*x22*x26*x33 + x24*(-x22
     &*x35 + x34 - x36*x5*x8 + x37*x5) - x27*x32*x33 + x32*( 0.000969022
     &771154437d0*x10*x37 - 0.000969022771154437d0*x11*x36 - x28*x29 + 4
     &.89119786774443d-5*x30 + 0.000359200055285857d0*x31*x9 ) - 1.82319
     &440383792d0*x30/(x15**2*x2**5)) + 0.151616336706985d0* x18 - 0.068
     &4394863448054d0*x23*x4 - x26*(-0.00193804554230887d0* x28 + x29)
      end subroutine


C*****************************************************************************
      pure subroutine ESRC_SPIN_VWN5(rho_c, rho_s, E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s
      real*8, intent(out) :: E
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
      E = 0.0d0
      x0 = rho_s/rho_c
      x1 = 1.92366105093154d0*(-x0 + 1.0d0)**1.33333333333333d0 + 1.9236
     &6105093154d0*(x0 + 1.0d0)**1.33333333333333d0 - 3.84732210186307d0
     &
      x2 = 0.6203504908994d0*rho_c**(-0.333333333333333d0)
      x3 = 0.826307487110758d0*rho_c**(-0.166666666666667d0)
      x4 = 1d0/(x2 + 6.72988144596143d0*x3 + 18.0578d0)
      x5 = 0.953184292996937d0*x3
      x6 = 1.90636858599387d0*x3
      x7 = 1d0/(x2 + 3.5529372610885d0*x3 + 12.9352d0)
      x8 = 0.0310907d0*log(x2*x7) + 0.000969022771154437d0*log(x7*(x5 + 
     &0.10498d0) **2) + 0.038783294878113d0*atan(6.1519908197590798d0/(x
     &6 + 3.7274400000000001d0))
      x9 = rho_s**4/rho_c**4
      x10 = 1d0/(x2 + 1.07811815828004d0*x3 + 13.0045d0)
      E = rho_c*(-x1*x9*(x8 - 0.01554535d0*log(x2*x4) - 0.00224786709554
     &261d0*log( x4*(x5 + 0.32500000000000001d0)**2) - 0.052491393169780
     &9d0*atan( 4.7309269095601136d0/(x6 + 7.0604199999999997d0))) - 0.0
     &0987581566084038d0*x1*(-x9 + 1.0d0)*(log(x10*x2) + 0.0004140337942
     &82063d0*log(x10*(x5 + 0.0047584000000000003d0)**2 ) + 0.3177080047
     &43941d0*atan(7.1231089178181177d0/(x6 + 1.13107d0 ))) + x8)
      end subroutine


C*****************************************************************************
      pure subroutine D1ESRC_SPIN_VWN5(rho_c, rho_s, E, d1E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s
      real*8, intent(out) :: E, d1E(9)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40, x41, x42, x43, x44, x45, x46, x47, x48
      E = 0.0d0
      d1E(:) = 0.0d0
      x0 = 0.682784063255296d0
      x1 = rho_c**0.333333333333333d0
      x2 = 0.90856029641607d0/x1
      x3 = x0*x2
      x4 = 0.826307487110758d0
      x5 = rho_c**(-0.166666666666667d0)*x4
      x6 = 1d0/(x3 + 6.72988144596143d0*x5 + 18.0578d0)
      x7 = 0.953184292996937d0*x5
      x8 = x7 + 0.325d0
      x9 = 1.90636858599387d0*x5
      x10 = x9 + 7.06042d0
      x11 = 1d0/(x3 + 3.5529372610885d0*x5 + 12.9352d0)
      x12 = x9 + 3.72744d0
      x13 = x7 + 0.10498d0
      x14 = 0.000969022771154437d0*log(x11*x13**2) + 0.0310907d0*log(x11
     &*x3) + 0.038783294878113d0*atan(6.1519908197590798d0/x12)
      x15 = x14 - 0.01554535d0*log(x3*x6) - 0.00224786709554261d0*log(x6
     &*x8**2) - 0.0524913931697809d0*atan(4.7309269095601136d0/x10)
      x16 = rho_s/rho_c
      x17 = x16 + 1.0d0
      x18 = -x16 + 1.0d0
      x19 = 1.92366105093154d0*x17**1.33333333333333d0 + 1.9236610509315
     &4d0*x18** 1.33333333333333d0 - 3.84732210186307d0
      x20 = rho_s**4
      x21 = x20/rho_c**4
      x22 = x19*x21
      x23 = 0.101321183642338d0
      x24 = -x21 + 1.0d0
      x25 = x9 + 1.13107d0
      x26 = 1d0/(x3 + 1.07811815828004d0*x5 + 13.0045d0)
      x27 = x7 + 0.0047584d0
      x28 = 0.000414033794282063d0*log(x26*x27**2) + log(x26*x3) + 0.317
     &708004743941d0*atan(7.1231089178181177d0/x25)
      x29 = x23*x24*x28
      x30 = 0.0974703937105774d0*x19
      x31 = x14 - x15*x22 - x29*x30
      x32 = rho_c**(-1.16666666666667d0)*x4
      x33 = 0.30285343213869d0*rho_c**(-1.33333333333333d0)
      x34 = x0*x33
      x35 = x6*(1.12164690766024d0*x32 + x34)
      x36 = x10**(-2)
      x37 = x11*(0.592156210181417d0*x32 + x34)
      x38 = x12**(-2)
      x39 = 0.0342197431724027d0*x1*(-x2*x37 + x33) - 0.0758081683534927
     &d0*x32*x38/( 37.8469910464d0*x38 + 1.0d0) + (-0.000969022771154437
     &d0*x13*x37 + 0.000307885761673592d0*x32)/x13
      x40 = 4.0d0*x15*x19
      x41 = x20/rho_c**5
      x42 = -x17**0.333333333333333d0 + x18**0.333333333333333d0
      x43 = 2.56488140124205d0*x15*x42
      x44 = x25**(-2)
      x45 = x26*(0.179686359713341d0*x32 + x34)
      x46 = 0.38988157484231d0*x19*x23*x28
      x47 = 0.25d0*x29*x42
      x48 = rho_s**3/rho_c**3
      E = rho_c*x31
      d1E(1) = -rho_c*(-x22*(-0.0171098715862014d0*x1*(-x2*x35 + x33) + 
     &0.0789023540332771d0*x32*x36/(22.3816694236d0*x36 + 1.0d0) + x39 -
     & (0.000714210536071954d0*x32 - 0.00224786709554261d0*x35*x8)/x8 ) 
     &+ x23*x24*x30*(-1.10064241629821d0*x1*(-x2*x45 + x33) + 0.71904051
     &9881222d0*x32*x44/(50.7386806551d0*x44 + 1.0d0) - ( -0.00041403379
     &4282063d0*x27*x45 + 0.000131550169826529d0*x32)/x27 ) + x39 - x40*
     &x41 + x41*x46 + rho_s*x47/rho_c**2 + rho_s**5*x43/ rho_c**6) + x31
     &
      d1E(2) = x21*x43 - x40*x48 + x46*x48 + x47
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_SPIN_VWN5(rho_c, rho_s, E, d1E, d2E)
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c, rho_s
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
     &, x171, x172, x173, x174, x175
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 0.682784063255296d0
      x1 = rho_c**0.333333333333333d0
      x2 = 0.90856029641607d0/x1
      x3 = x0*x2
      x4 = 0.826307487110758d0
      x5 = rho_c**(-0.166666666666667d0)*x4
      x6 = x3 + 6.72988144596143d0*x5 + 18.0578d0
      x7 = 1d0/x6
      x8 = 0.01554535d0*log(x3*x7)
      x9 = 0.953184292996937d0*x5
      x10 = x9 + 0.325d0
      x11 = x10**2
      x12 = x11*x7
      x13 = 0.00224786709554261d0*log(x12)
      x14 = 1.90636858599387d0*x5
      x15 = x14 + 7.06042d0
      x16 = 0.0524913931697809d0*atan(4.7309269095601136d0/x15)
      x17 = x3 + 3.5529372610885d0*x5 + 12.9352d0
      x18 = 1d0/x17
      x19 = 0.0310907d0*log(x18*x3)
      x20 = x14 + 3.72744d0
      x21 = 0.038783294878113d0*atan(6.1519908197590798d0/x20)
      x22 = x9 + 0.10498d0
      x23 = x22**2
      x24 = x18*x23
      x25 = 0.000969022771154437d0*log(x24)
      x26 = x19 + x21 + x25
      x27 = -x13 - x16 + x26 - x8
      x28 = 1d0/rho_c
      x29 = rho_s*x28
      x30 = x29 + 1.0d0
      x31 = -x29 + 1.0d0
      x32 = 1.92366105093154d0*x30**1.33333333333333d0 + 1.9236610509315
     &4d0*x31** 1.33333333333333d0 - 3.84732210186307d0
      x33 = rho_c**(-4)
      x34 = rho_s**4
      x35 = x33*x34
      x36 = x32*x35
      x37 = x14 + 1.13107d0
      x38 = x3 + 1.07811815828004d0*x5 + 13.0045d0
      x39 = 1d0/x38
      x40 = x9 + 0.0047584d0
      x41 = x40**2
      x42 = x39*x41
      x43 = 0.000414033794282063d0*log(x42) + log(x3*x39) + 0.3177080047
     &43941d0*atan (7.1231089178181177d0/x37)
      x44 = 0.101321183642338d0
      x45 = x44*(-x35 + 1.0d0)
      x46 = x43*x45
      x47 = 0.0974703937105774d0*x32
      x48 = x26 - x27*x36 - x46*x47
      x49 = rho_c**(-1.16666666666667d0)*x4
      x50 = x20**(-2)
      x51 = 37.8469910464d0*x50 + 1.0d0
      x52 = 1d0/x51
      x53 = x50*x52
      x54 = x49*x53
      x55 = -0.0758081683534927d0*x54
      x56 = 1d0/x10
      x57 = 0.000714210536071954d0*x49
      x58 = rho_c**(-1.33333333333333d0)
      x59 = 0.30285343213869d0*x58
      x60 = x0*x59
      x61 = 1.12164690766024d0*x49 + x60
      x62 = x61*x7
      x63 = x10*x62
      x64 = 0.00224786709554261d0*x63
      x65 = x2*x62
      x66 = 0.0171098715862014d0*x1
      x67 = x15**(-2)
      x68 = 22.3816694236d0*x67 + 1.0d0
      x69 = 1d0/x68
      x70 = x67*x69
      x71 = 0.0789023540332771d0*x49*x70 + x55
      x72 = 1d0/x22
      x73 = 0.000307885761673592d0*x49
      x74 = 0.592156210181417d0*x49 + x60
      x75 = x18*x74
      x76 = x22*x75
      x77 = 0.000969022771154437d0*x76
      x78 = x2*x75
      x79 = x59 - x78
      x80 = 0.0342197431724027d0*x1
      x81 = x72*(x73 - x77) + x79*x80
      x82 = -x56*(x57 - x64) - x66*(x59 - x65) + x71 + x81
      x83 = x36*x82
      x84 = x34/rho_c**5
      x85 = x27*x84
      x86 = 4.0d0*x32
      x87 = x31**0.333333333333333d0
      x88 = x30**0.333333333333333d0
      x89 = x87 - x88
      x90 = 2.56488140124205d0*x89
      x91 = x27*x90
      x92 = rho_s**5
      x93 = rho_c**(-6)
      x94 = x92*x93
      x95 = x37**(-2)
      x96 = 50.7386806551d0*x95 + 1.0d0
      x97 = 1d0/x96
      x98 = x95*x97
      x99 = 0.719040519881222d0*x49*x98
      x100 = 0.179686359713341d0*x49 + x60
      x101 = x100*x39
      x102 = x101*x2
      x103 = 1.10064241629821d0*x1
      x104 = 1d0/x40
      x105 = 0.000131550169826529d0*x49
      x106 = x101*x40
      x107 = 0.000414033794282063d0*x106
      x108 = -x103*(-x102 + x59) - x104*(x105 - x107) + x99
      x109 = x45*x47
      x110 = x32*x84
      x111 = x43*x44
      x112 = 0.38988157484231d0*x111
      x113 = 0.25d0*x89
      x114 = x113*x46
      x115 = rho_c**(-2)
      x116 = rho_s*x115
      x117 = rho_c**(-3)
      x118 = rho_s**3
      x119 = x118*x32
      x120 = x117*x119
      x121 = x118*x27
      x122 = x117*x86
      x123 = 0.000615771523347183d0*x49
      x124 = 0.77976314968462d0*x110
      x125 = 0.5d0*x116*x89
      x126 = 5.1297628024841d0*x89*x94
      x127 = 8.0d0*x32
      x128 = x108*x45
      x129 = 1d0/x11
      x130 = rho_c**(-2.33333333333333d0)*x0
      x131 = rho_c**(-2.16666666666667d0)*x4
      x132 = 0.40380457618492d0*rho_c**(-2.33333333333333d0)
      x133 = x0*x132
      x134 = 1.30858805893694d0*x131 + x133
      x135 = 0.60570686427738d0*x58
      x136 = x0*x135
      x137 = x61*(x136 + 2.24329381532048d0*x49)/x6**2
      x138 = rho_c**(-0.666666666666667d0)
      x139 = -x59
      x140 = x139 + x65
      x141 = -0.317728097665645d0*x49
      x142 = x141 + x63
      x143 = x140*x66
      x144 = 1d0/x23
      x145 = 0.69084891187832d0*x131 + x133
      x146 = x74*(x136 + 1.18431242036283d0*x49)/x17**2
      x147 = x139 + x78
      x148 = x147*x80
      x149 = x141 + x76
      x150 = -0.0481727702369445d0*x130*x52/x20**3 + 1.82319440383792d0*
     &x130/(x20**5* x51**2) + 0.0884428630790749d0*x131*x53 - 0.01140658
     &10574676d0* x138*x147 - x144*x149*x73 - x144*(-x123*x76 + 4.891197
     &86774443d-5 *x130 + 0.000359200055285857d0*x131*x22 - 0.0009690227
     &71154437d0* x145*x24 + 0.000969022771154437d0*x146*x23) + x148*x75
     & + 0.000969022771154437d0*x149*x72*x75 - x80*(x132 - x135*x75 - x1
     &45 *x18*x2 + x146*x2)
      x151 = x30**(-0.666666666666667d0)
      x152 = 0.854960467080683d0*x29
      x153 = x31**(-0.666666666666667d0)
      x154 = x151*x152 + x152*x153
      x155 = x154 - 5.1297628024841d0*x87 + 5.1297628024841d0*x88
      x156 = x13 + x16 - x19 - x21 - x25 + x8
      x157 = x92/rho_c**7
      x158 = x156*x157
      x159 = x32*x34*x93
      x160 = x143 - x148 + x56*(-x57 + x64) + x71 - x72*(-x73 + x77)
      x161 = 20.5190512099364d0*x89
      x162 = x102 + x139
      x163 = 1d0/x41
      x164 = x106 + x141
      x165 = x103*x162
      x166 = 0.209634086332231d0*x131 + x133
      x167 = x100*(x136 + 0.359372719426682d0*x49)/x38**2
      x168 = 2.0d0*x111*x89
      x169 = x104*(-x105 + x107) + x165 + x99
      x170 = 0.0974703937105774d0*x46
      x171 = 1.16964472452693d0*x111
      x172 = x154 - 2.56488140124205d0*x87 + 2.56488140124205d0*x88
      x173 = 12.0d0*x32
      x174 = x151 + x153
      x175 = rho_s**2*x115
      E = rho_c*x48
      d1E(1) = -rho_c*(x108*x109 + x110*x112 + x114*x116 + x55 + x81 - x
     &83 - x85*x86 + x91*x94) + x48
      d1E(2) = x112*x120 + x114 - x121*x122 + x35*x91
      d2E(1) = -rho_c*(rho_s*x117*x155*x170 + x109*(-0.00041403379428206
     &3d0*x101*x104* x164 - x101*x165 + x103*(-x101*x135 + x132 - x166*x
     &2*x39 + x167* x2) + x105*x163*x164 + 0.456918753052755d0*x130*x97/
     &x37**3 - 23.1834546964702d0*x130/(x37**5*x96**2) - 0.8388806065280
     &93d0* x131*x98 + 0.366880805432736d0*x138*x162 + x163*( -0.0002631
     &00339653058d0*x106*x49 + 2.08985926032878d-5*x130 + 0.000153475198
     &130951d0*x131*x40 - 0.000414033794282063d0*x166*x42 + 0.0004140337
     &94282063d0*x167*x41)) - 1.94940787421155d0*x111* x159 + x124*x169*
     &x44 + x125*x169*x45 - x126*x160 + x127*x160*x84 + x150 - x155*x158
     & - 20.0d0*x156*x159 + x157*x168 + x158*x161 - x36*(x129*x142*x57 +
     & x129*(0.000833245625417279d0*x10*x131 + 0.00224786709554261d0*x11
     &*x137 - 0.00224786709554261d0*x12*x134 + 0.000113462377479451d0*x1
     &30 - 0.00142842107214391d0*x49*x63) + 0.0501389896966688d0*x130*x6
     &9/x15**3 - 1.12219429262413d0*x130/( x15**5*x68**2) - 0.0920527463
     &721566d0*x131*x70 + 0.00570329052873379d0*x138*x140 - 0.0022478670
     &9554261d0*x142*x56* x62 - x143*x62 + x150 + x66*(x132 - x134*x2*x7
     & - x135*x62 + x137* x2))) - 0.0684394863448054d0*x1*x79 - x111*x12
     &4 - x125*x46 - x126 *x27 + x127*x85 - 0.194940787421155d0*x128*x32
     & + 0.151616336706985d0*x54 - x72*(x123 - 0.00193804554230887d0*x76
     &) + 2.0d0*x83
      d2E(2) = 0.38988157484231d0*x108*x120*x44 + x113*x128 + x114*x28 +
     & x118*x122*x82 - x119*x171*x33 + x121*x173*x33 + x168*x84 + x170*x
     &172*x28 + x172 *x85 - x35*x82*x90 - 17.9541698086943d0*x85*x89
      d2E(3) = -x28*(x117*x118*x168 - x117*x121*x161 - x171*x175*x32 + x
     &173*x175*x27 + 0.854960467080683d0*x174*x27*x35 + 0.08333333333333
     &33d0*x174*x46)
      end subroutine


C*****************************************************************************
      pure subroutine D2ESRC_VWN5_singletref_triplet(rho_c, E, d1E, d2E)
     &
C*****************************************************************************
C   Implemented by E.R. Kjellgren.
C
C   Subroutine generated using Sympy 1.3
C   Generated: March 25, 2019
C*****************************************************************************
      implicit none
      real*8, intent(in) :: rho_c
      real*8, intent(out) :: E, d1E(9), d2E(45)
      real*8 :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x
     &13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x2
     &6, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39
     &, x40
      E = 0.0d0
      d1E(:) = 0.0d0
      d2E(:) = 0.0d0
      x0 = 0.826307487110758d0
      x1 = rho_c**(-0.166666666666667d0)*x0
      x2 = 1.90636858599387d0*x1
      x3 = x2 + 3.72744d0
      x4 = 0.682784063255296d0
      x5 = rho_c**0.333333333333333d0
      x6 = 0.90856029641607d0/x5
      x7 = x4*x6
      x8 = 3.5529372610885d0*x1 + x7 + 12.9352d0
      x9 = 1d0/x8
      x10 = 0.953184292996937d0*x1
      x11 = x10 + 0.10498d0
      x12 = x11**2
      x13 = x12*x9
      x14 = 0.000969022771154437d0*log(x13) + 0.0310907d0*log(x7*x9) + 0
     &.038783294878113d0*atan(6.1519908197590798d0/x3)
      x15 = rho_c**(-1.16666666666667d0)*x0
      x16 = x3**(-2)
      x17 = 37.8469910464d0*x16 + 1.0d0
      x18 = 1d0/x17
      x19 = x16*x18
      x20 = x15*x19
      x21 = rho_c**(-1.33333333333333d0)
      x22 = 0.30285343213869d0*x21
      x23 = 0.592156210181417d0*x15 + x22*x4
      x24 = x23*x9
      x25 = x22 - x24*x6
      x26 = 0.0342197431724027d0*x5
      x27 = x25*x26
      x28 = 1d0/x11
      x29 = 0.000307885761673592d0*x15
      x30 = x11*x24
      x31 = 0.000615771523347183d0*x15
      x32 = rho_c**(-2.33333333333333d0)*x4
      x33 = rho_c**(-2.16666666666667d0)*x0
      x34 = 1d0/x12
      x35 = 0.317728097665645d0*x15 - x30
      x36 = 0.40380457618492d0*rho_c**(-2.33333333333333d0)
      x37 = 0.60570686427738d0*x21
      x38 = 0.69084891187832d0*x33 + x36*x4
      x39 = x23*(1.18431242036283d0*x15 + x37*x4)/x8**2
      x40 = 1d0/(1.07811815828004d0*x1 + x7 + 13.0045d0)
      E = rho_c*x14
      d1E(1) = -rho_c*(-0.0758081683534927d0*x20 + x27 + x28*(x29 - 0.00
     &0969022771154437d0*x30)) + x14
      d2E(1) = rho_c*(-0.0114065810574676d0*rho_c**(-0.666666666666667d0
     &)*x25 + 0.0481727702369445d0*x18*x32/x3**3 - 0.0884428630790749d0*
     &x19*x33 + x24*x27 + 0.000969022771154437d0*x24*x28*x35 + x26*(-x24
     &*x37 + x36 - x38*x6*x9 + x39*x6) - x29*x34*x35 + x34*( 0.000359200
     &055285857d0*x11*x33 + 0.000969022771154437d0*x12*x39 - 0.000969022
     &771154437d0*x13*x38 - x30*x31 + 4.89119786774443d-5* x32) - 1.8231
     &9440383792d0*x32/(x17**2*x3**5)) + 0.151616336706985d0*x20 - 0.068
     &4394863448054d0*x25*x5 - x28*( -0.00193804554230887d0*x30 + x31)
      d2E(3) = -0.0168868639403896d0*(log(x40*x7) + 0.000414033794282063
     &d0*log(x40*(x10 + 0.0047584000000000003d0)**2) + 0.317708004743941
     &d0*atan( 7.1231089178181177d0/(x2 + 1.13107d0)))/rho_c
      end subroutine
