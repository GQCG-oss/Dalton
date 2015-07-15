!============================!
!  Module fundamental data   !
!============================!
!
! Common block codata.h and DISOTP function which returns
! properties of nuclei both rewritten by Vladimir Rybkin in 2010
! for linear scaling
!  codata.h revised 1999/04/29 and combined with units.h 
!
!
!  From
!  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
!                     Peter J. Mohr and Barry N. Taylor
!  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
!
! -- Fundamental constants, atomic units, Avogadro's constant
Module Fundamental
Use precision
Implicit none
REAL(realk) CVEL, ALPHAC, ALPHA2, PI, SQRTPI,R2PI52
REAL(realk) bohr_to_angstrom,XFAMU,ECHARGE,HBAR,XFMOL
REAL(realk) XTJ,XTKAYS,XTHZ,XTEV,XKJMOL,XKCMOL,XTKJML,XTKCML,XTNM,XAJOUL,bohr_to_angstromM10,XPRTMAS
REAL(realk) XFSEC,XTKMML,TESLA,AUTK,DEBYE,PMASS,EMASS,CCM,hartree_to_eV,hartree_to_cm1
!
PARAMETER ( bohr_to_angstrom = 0.5291772083E0_realk, XFAMU  = 1822.88848E0_realk, &   
     &    ECHARGE = 1.602176462E-19_realk, HBAR = 1.054571596E-34_realk)
PARAMETER ( XFMOL = 6.02214199E23_realk )
PARAMETER ( bohr_to_angstromM10 = bohr_to_angstrom*1.0E-10_realk)
!
PARAMETER ( PMASS = 1.007276470E0_realk, EMASS = 9.10938188E-31_realk)
PARAMETER (hartree_to_eV = 27.211399)
PARAMETER (hartree_to_cm1 = 219474.63)
!
!-- Pi number
      PARAMETER (PI     = 3.14159265358979323846E00_realk, &
     &           SQRTPI = 1.77245385090551602730E00_realk, &
     &           R2PI52 = 5.91496717279561287782E00_realk)
!     R2PI52 = sqrt(2 * sqrt(PI^5) ) -- used in calc. of 2-el. integrals
!
!     Fine structure constant
!
      PARAMETER (CCM = 299792458.0E0_realk)
      PARAMETER (CVEL = CCM*bohr_to_angstromM10*EMASS/(HBAR))
      PARAMETER (ALPHAC = 1.0E0_realk/CVEL, ALPHA2 = ALPHAC*ALPHAC)
!
! -- conversion from hartree (au) to:
#ifdef SYS_REAL
PARAMETER ( XTJ  = 4.35974380425140E-18_realk, &
     &    XTHZ =   6.57968391802650E+15_realk, &  
#else
PARAMETER ( XTJ  = HBAR**2/(bohr_to_angstromM10*bohr_to_angstromM10*EMASS), &
     &    XTHZ =  HBAR/(2.0E0_realk*PI*bohr_to_angstromM10*bohr_to_angstromM10*EMASS), &  
#endif
     &    XTKAYS = 1.0E-2_realk*XTHZ/CCM, &
     &    XTEV = XTJ/ECHARGE, & 
     &    XKJMOL = XTJ*XFMOL*1E-3_realk, XKCMOL = XKJMOL/4.184E0_realk, &
     &             XTKJML = XKJMOL, XTKCML = XKCMOL, &  
     &    XTNM = 1E7_realk/XTKAYS, XAJOUL = 1.0E18_realk*XTJ)
!
! -- other
      PARAMETER ( XFSEC = HBAR/XTJ)
      PARAMETER ( XTKMML = 974.864E0_realk)
      PARAMETER ( TESLA=(bohr_to_angstrom*bohr_to_angstrom*ECHARGE/HBAR)*1E-20_realk )
      PARAMETER ( AUTK = 3.1577465E5_realk )
      PARAMETER ( DEBYE = ECHARGE*bohr_to_angstrom*CCM*1E11_realk )
      PARAMETER ( XPRTMAS = 1836.1526675E0_realk )   
!
Interface Isotopes
    Module procedure Isotopes
End interface
Interface CovRad
    Module procedure CovRad
End Interface
! 
Contains
!
Function Isotopes(Atom,NIsotop,TYPE,LUPRI)
Use Precision
Implicit none
!
!     NOTE: Isotopes are sorted according to abundance,
!     i.e. Isotops(ATOM,1,TYPE,LUPRI) will return the most abundant
!     isotope etc.
!
!     Extended to Z=54 on march 10, 1994, K.Ruud
!     Extended to Z=86 feb. 1996, S. Kirpekar
!
!     Proton mass and electron charge:
!        1986 CODATA Recommended Values
!
!     Nuclear masses:
!        A. H. Wapstra and K. Bos,
!        Atomic Data and Nuclear Tables 19 (1977) 177
!
!     Abundancies:
!        Handbook of Chemistry and Physics, 73rd Edition
!        
!     Nuclear moments and spins:
!
!        P. Raghavan,
!        Atomic Data and Nuclear Data Tables 42 (1989) 189
!
!     Quadrupole moments:
!
!        P.Pykkoe and J.Li
!        Report HUKI 1-92
!        Updated  to P.Pyykkoe, Mol.Phys. (2001) by K.Ruud, Aug.2001
!
!     Nuclear masses, Abundancies, nuclear moments, spins 
!     and quadrupole moments for Z= 55 to Z = 86:
!
!       I. Mills, T. Cvitas, K. Homann, N. Kallay, and K. Kuchitsu
!       Quantities, Units and Symbols in Physical Chemistry
!       (IUPAC, Blackwell Scientific, Oxford, 1988)
!
!
      INTEGER :: MAXISO, MAXCHR
      INTEGER :: NIsotop,LUPri,Atom
      INTEGER :: I,J,K
      REAL(realk) :: Isotopes, D0, D4, DMP,THRESH, SPIN
      REAL(realk) :: DATNUC
      PARAMETER (D0 = 0.0E0_realk, D4= 4.0E0_realk, MAXISO = 6, MAXCHR = 86)
      PARAMETER (DMP = PMASS*XFAMU*EMASS)
      PARAMETER (THRESH = 1.0E-10_realk)
!
      CHARACTER*(*) TYPE
      DIMENSION DATNUC(5,MAXISO,0:MAXCHR)
!
!     H - Ne
!     ======
!
!
!
!     Dummy:
!
      DATA ((DATNUC(I,J,0),I=1,5),J=1,MAXISO) / &
       0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk/ 
!
!     H:
!
      DATA ((DATNUC(I,J,1),I=1,5),J=1,MAXISO) / &
 &     1.007825E0_realk, 99.985000E0_realk,   .500000E0_realk,  2.792847E0_realk,   .000000E0_realk, &
 &     2.014102E0_realk,   .015000E0_realk,  1.000000E0_realk,   .857438E0_realk,   .002860E0_realk, &
 &     3.016049E0_realk,   .000000E0_realk,   .500000E0_realk,  2.978962E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     He:
!
      DATA ((DATNUC(I,J,2),I=1,5),J=1,MAXISO) / &
 &     4.002603E0_realk, 99.999870E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &     3.016029E0_realk,   .000130E0_realk,   .500000E0_realk, -2.127625E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/

!     Li:
!
      DATA ((DATNUC(I,J,3),I=1,5),J=1,MAXISO) / &
 &     7.016005E0_realk, 92.500000E0_realk,  1.500000E0_realk,  3.256427E0_realk,  -.040100E0_realk, &
 &     6.015123E0_realk,  7.500000E0_realk,  1.000000E0_realk,   .822047E0_realk,  -.000808E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/

!     Be:
!
      DATA ((DATNUC(I,J,4),I=1,5),J=1,MAXISO) / &
 &     9.012183E0_realk,100.000000E0_realk,  1.500000E0_realk, -1.177800E0_realk,   .052880E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     B:
!
      DATA ((DATNUC(I,J,5),I=1,5),J=1,MAXISO) / &
 &    11.009305E0_realk, 80.100000E0_realk,  1.500000E0_realk,  2.688649E0_realk,   .040590E0_realk, &
 &    10.012938E0_realk, 19.900000E0_realk,  3.000000E0_realk,  1.800645E0_realk,   .084590E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     C:
!
      DATA ((DATNUC(I,J,6),I=1,5),J=1,MAXISO) / &
 &    12.000000E0_realk, 98.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    13.003355E0_realk,  1.100000E0_realk,   .500000E0_realk,   .702412E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
! 
!     N:
!
      DATA ((DATNUC(I,J,7),I=1,5),J=1,MAXISO) / &
 &    14.003074E0_realk, 99.630000E0_realk,  1.000000E0_realk,   .403761E0_realk,   .020440E0_realk, &
 &    15.000109E0_realk,   .370000E0_realk,   .500000E0_realk,  -.283189E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    12.000000E0_realk, 98.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
! 
!     O:
!
      DATA ((DATNUC(I,J,8),I=1,5),J=1,MAXISO) / &
 &    15.994915E0_realk, 99.760000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    17.999159E0_realk,   .200000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    16.999131E0_realk,   .040000E0_realk,  2.500000E0_realk, -1.893790E0_realk,  -.025580E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    14.003074E0_realk, 99.630000E0_realk,  1.000000E0_realk,   .403761E0_realk,   .020200E0_realk/
!
!     F:
!
      DATA ((DATNUC(I,J,9),I=1,5),J=1,MAXISO) / &
 &    18.998403E0_realk,100.000000E0_realk,   .500000E0_realk,  2.628868E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    15.994915E0_realk, 99.760000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Ne:
!
      DATA ((DATNUC(I,J,10),I=1,5),J=1,MAXISO) / &
 &    19.992439E0_realk, 90.480000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, & 
 &    21.991384E0_realk,  9.250000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    20.993845E0_realk,   .270000E0_realk,  1.500000E0_realk,  -.661797E0_realk,   .101550E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/ 
!
!    Na - Ar
!    =======
!
!
!     Na:
!
      DATA ((DATNUC(I,J,11),I=1,5),J=1,MAXISO) / &
 &    22.989770E0_realk,100.000000E0_realk,  1.500000E0_realk,  2.217656E0_realk,   .104000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Mg:
!
      DATA ((DATNUC(I,J,12),I=1,5),J=1,MAXISO) / &
 &    23.985045E0_realk, 78.990000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    25.982595E0_realk, 11.010000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    24.985839E0_realk, 10.000000E0_realk,  2.500000E0_realk,  -.855450E0_realk,   .199400E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Al:
!
      DATA ((DATNUC(I,J,13),I=1,5),J=1,MAXISO) / &
 &    26.981541E0_realk,100.000000E0_realk,  2.500000E0_realk,  3.641507E0_realk,   .146600E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Si:
!
      DATA ((DATNUC(I,J,14),I=1,5),J=1,MAXISO) / &
 &    27.976928E0_realk, 92.230000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    28.976496E0_realk,  4.670000E0_realk,   .500000E0_realk,  -.555290E0_realk,   .000000E0_realk, &
 &    29.973772E0_realk,  3.100000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     P:
!
      DATA ((DATNUC(I,J,15),I=1,5),J=1,MAXISO) / &
 &    30.973763E0_realk,100.000000E0_realk,   .500000E0_realk,  1.131600E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     S:
!
      DATA ((DATNUC(I,J,16),I=1,5),J=1,MAXISO) / &
 &    31.972072E0_realk, 95.020000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    33.967868E0_realk,  4.210000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    32.971459E0_realk,   .750000E0_realk,  1.500000E0_realk,   .643821E0_realk,  -.067800E0_realk, &
 &    35.967079E0_realk,   .020000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Cl:
!
      DATA ((DATNUC(I,J,17),I=1,5),J=1,MAXISO) / &
 &    34.968853E0_realk, 75.770000E0_realk,  1.500000E0_realk,   .821874E0_realk,  -.081650E0_realk, &
 &    36.965903E0_realk, 24.230000E0_realk,  1.500000E0_realk,   .684124E0_realk,  -.064350E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Ar:
!
      DATA ((DATNUC(I,J,18),I=1,5),J=1,MAXISO) / &
 &    39.962383E0_realk, 99.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    35.967546E0_realk,   .337000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    37.962732E0_realk,   .063000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/ 
!
!     K - Ca
!     ======
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=19,20) / &
!
!     K:
!
      38.963708E0_realk, 93.258100E0_realk,  1.500000E0_realk,   .391507E0_realk,   .058500E0_realk, &
 &    40.961825E0_realk,  6.730200E0_realk,  1.500000E0_realk,   .214893E0_realk,   .071100E0_realk, &
 &    39.963999E0_realk,   .011700E0_realk,  4.000000E0_realk, -1.298100E0_realk,  -.073000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Ca:
!
 &    39.962591E0_realk, 96.941000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    43.955485E0_realk,  2.086000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    41.958622E0_realk,   .647000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    47.952532E0_realk,   .187000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    42.958770E0_realk,   .135000E0_realk,  3.500000E0_realk, -1.317643E0_realk,  -.040800E0_realk, &
 &    45.953689E0_realk,   .004000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/ 
!
!     Sc - Zn
!     =======
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=21,25) / &
!
!     Sc:
!
      44.955914E0_realk,100.000000E0_realk,  3.500000E0_realk,  4.756487E0_realk,  -.220000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Ti:
!
 &    47.947947E0_realk, 73.800000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    45.952633E0_realk,  8.000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    46.951765E0_realk,  7.300000E0_realk,  2.500000E0_realk,  -.788480E0_realk,   .302000E0_realk, &
 &    48.947871E0_realk,  5.500000E0_realk,  3.500000E0_realk, -1.104170E0_realk,   .247000E0_realk, &
 &    49.944786E0_realk,  5.400000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     V:
!
 &    50.943963E0_realk, 99.750000E0_realk,  3.500000E0_realk,  5.148706E0_realk,  -.052000E0_realk, &
 &    49.947161E0_realk,   .250000E0_realk,  6.000000E0_realk,  3.345689E0_realk,   .210000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     Cr:
!
 &    51.940510E0_realk, 83.790000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    52.940651E0_realk,  9.500000E0_realk,  1.500000E0_realk,  -.474540E0_realk,  -.150000E0_realk, &
 &    49.946046E0_realk,  4.345000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    53.938882E0_realk,  2.365000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Mn:
!
 &    54.938046E0_realk,100.000000E0_realk,  2.500000E0_realk,  3.468719E0_realk,   .330000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/


      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=26,30) / &
!
!     Fe:
!
 &    55.934939E0_realk, 91.720000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    53.939612E0_realk,  5.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    56.935396E0_realk,  2.100000E0_realk,   .500000E0_realk,   .090623E0_realk,   .000000E0_realk, &
 &    57.933278E0_realk,   .280000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Co:
!
 &    58.933198E0_realk,100.000000E0_realk,  3.500000E0_realk,  4.627000E0_realk,   .420000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Ni:
!
 &    57.935347E0_realk, 68.077000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    59.930789E0_realk, 26.223000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    61.928346E0_realk,  3.634000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    60.931059E0_realk,  1.140000E0_realk,  1.500000E0_realk,  -.750020E0_realk,   .162000E0_realk, &
 &    63.927968E0_realk,  0.926000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Cu:
!
 &    62.929599E0_realk, 69.170000E0_realk,  1.500000E0_realk,  2.227206E0_realk,  -.220000E0_realk, &
 &    64.927792E0_realk, 30.830000E0_realk,  1.500000E0_realk,  2.381610E0_realk,  -.204000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Zn:
!
 &    63.929145E0_realk, 48.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    65.926035E0_realk, 27.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    67.924846E0_realk, 18.800000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    66.927129E0_realk,  4.100000E0_realk,  2.500000E0_realk,   .875479E0_realk,   .150000E0_realk, &
 &    69.925325E0_realk,   .600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Ga - Kr
!     =======
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=31,36) / &
!
!     Ga:
!
      68.925581E0_realk, 60.108000E0_realk,  1.500000E0_realk,  2.016589E0_realk,   .171000E0_realk, &
 &    70.924701E0_realk, 39.892000E0_realk,  1.500000E0_realk,  2.562266E0_realk,   .107000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Ge:
!
 &    73.921179E0_realk, 35.940000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    71.922080E0_realk, 27.660000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    69.924250E0_realk, 21.240000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    72.923464E0_realk,  7.720000E0_realk,  4.500000E0_realk,  -.879468E0_realk,  -.196000E0_realk, &
 &    75.921403E0_realk,  7.440000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     As:
!
 &    74.921596E0_realk,100.000000E0_realk,  1.500000E0_realk,  1.439475E0_realk,   .314000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     Se:
!
 &    79.916521E0_realk, 49.610000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    77.917304E0_realk, 23.770000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    75.919207E0_realk,  9.360000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    81.916709E0_realk,  8.740000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    76.919908E0_realk,  7.630000E0_realk,   .500000E0_realk,   .535042E0_realk,   .000000E0_realk, &
 &    73.922477E0_realk,  0.890000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Br:
!
 &    78.918336E0_realk, 50.690000E0_realk,  1.500000E0_realk,  2.106400E0_realk,   .313000E0_realk, &
 &    80.916290E0_realk, 49.310000E0_realk,  1.500000E0_realk,  2.270562E0_realk,   .261500E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Kr:
!
 &    83.911506E0_realk, 57.000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    85.910614E0_realk, 17.300000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    81.913483E0_realk, 11.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    82.914134E0_realk, 11.500000E0_realk,  4.500000E0_realk,  -.970669E0_realk,   .259000E0_realk, &
 &    79.916375E0_realk,  2.250000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    77.920397E0_realk,  0.350000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Rb:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=37,40) / &
!
      84.911800E0_realk, 72.170000E0_realk,  2.500000E0_realk,  1.353352E0_realk,   .276000E0_realk, & 
 &    86.909184E0_realk, 27.830000E0_realk,  1.500000E0_realk,  2.751818E0_realk,   .133500E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Sr:
!
 &    87.905625E0_realk, 82.580000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    85.909273E0_realk,  9.860000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    86.908890E0_realk,  7.000000E0_realk,  4.500000E0_realk, -1.093603E0_realk,   .335000E0_realk, &
 &    83.913428E0_realk,   .560000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     Y:
!
 &    88.905856E0_realk,100.000000E0_realk,   .500000E0_realk,  -.137415E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     Zr:
!
 &    89.904708E0_realk, 51.450000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    93.906319E0_realk, 17.380000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    91.905039E0_realk, 17.150000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    90.905644E0_realk, 11.220000E0_realk,  2.500000E0_realk, -1.303620E0_realk,  -.176000E0_realk, &
 &    95.908272E0_realk,  2.800000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!
!     Nb:
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=41,45) / &
!
 &    92.906378E0_realk,100.000000E0_realk,  4.500000E0_realk,  6.170500E0_realk,  -.320000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Mo:
!
 &    97.905405E0_realk, 24.130000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    95.904676E0_realk, 16.680000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    94.905838E0_realk, 15.920000E0_realk,  2.500000E0_realk,  -.914200E0_realk,  -.022000E0_realk, &
 &    93.905086E0_realk, 14.840000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    99.907473E0_realk,  9.630000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    96.906018E0_realk,  9.550000E0_realk,  2.500000E0_realk,  -.933500E0_realk,  0.255000E0_realk, &
!
!     Tc:
!
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Ru:
!
 &   101.904348E0_realk, 31.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   103.905422E0_realk, 18.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   100.905581E0_realk, 17.100000E0_realk,  2.500000E0_realk,  -.718800E0_realk,   .457000E0_realk, &
 &    98.905937E0_realk, 12.700000E0_realk,  2.500000E0_realk,  -.641300E0_realk,   .079000E0_realk, &
 &    99.904218E0_realk, 12.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &    95.907596E0_realk,  5.540000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Rh:
!
 &   102.905503E0_realk,100.000000E0_realk,   .500000E0_realk,  -.088400E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/ 
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=46,50) / &
!
!     Pd:
!
     105.903475E0_realk, 27.330000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   107.903894E0_realk, 26.460000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   104.905075E0_realk, 22.330000E0_realk,  2.500000E0_realk,  -.642000E0_realk,   .660000E0_realk, &
 &   109.905169E0_realk, 11.720000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   103.904026E0_realk, 11.140000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   101.905609E0_realk,  1.020000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     Ag:
!
 &   106.905095E0_realk, 51.839000E0_realk,   .500000E0_realk,  -.113570E0_realk,   .000000E0_realk, &
 &   108.904754E0_realk, 48.161000E0_realk,   .500000E0_realk,   .130563E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Cd:
!
 &   113.903361E0_realk, 28.730000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   111.902761E0_realk, 24.130000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   110.904182E0_realk, 12.800000E0_realk,   .500000E0_realk,  -.594886E0_realk,   .000000E0_realk, &
 &   109.903007E0_realk, 12.490000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   112.904401E0_realk, 12.220000E0_realk,   .500000E0_realk,  -.622301E0_realk,   .000000E0_realk, &
 &   115.904758E0_realk,  7.490000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     In:
!
 &   114.903875E0_realk, 95.700000E0_realk,  4.500000E0_realk,  5.540800E0_realk,   .810000E0_realk, &
 &   112.904056E0_realk,  4.300000E0_realk,  4.500000E0_realk,  5.528900E0_realk,   .799000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Sn:
!
 &   119.902199E0_realk, 32.590000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   117.901607E0_realk, 24.220000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   115.901744E0_realk, 14.530000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   118.903310E0_realk,  8.580000E0_realk,   .500000E0_realk, -1.047280E0_realk,   .000000E0_realk, &
 &   116.902954E0_realk,  7.680000E0_realk,   .500000E0_realk, -1.001040E0_realk,   .000000E0_realk, &
 &   123.905271E0_realk,  5.790000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
! 
!     Sb:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=51,54) / &
 &   120.903824E0_realk, 57.360000E0_realk,  2.500000E0_realk,  3.363400E0_realk,  -.360000E0_realk, &
 &   122.904222E0_realk, 42.640000E0_realk,  3.500000E0_realk,  2.549800E0_realk,  -.490000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     Te:
!
 &   129.906229E0_realk, 33.870000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   127.904464E0_realk, 31.700000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   125.903310E0_realk, 18.930000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   124.904435E0_realk,  7.120000E0_realk,   .500000E0_realk,  -.888505E0_realk,   .000000E0_realk, &
 &   123.902825E0_realk,  4.790000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   121.903055E0_realk,  2.590000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
! 
!     I:
!
 &   126.904477E0_realk,100.000000E0_realk,  2.500000E0_realk,  2.813273E0_realk,  -.710000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &      .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
!
!     Xe:
!
 &   131.904148E0_realk, 26.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   128.904780E0_realk, 26.400000E0_realk,   .500000E0_realk,  -.777976E0_realk,   .000000E0_realk, &
 &   130.905076E0_realk, 21.200000E0_realk,  1.500000E0_realk,   .691862E0_realk,  -.114000E0_realk, &
 &   133.905395E0_realk, 10.400000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   135.907219E0_realk,  8.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk, &
 &   129.903510E0_realk,  4.100000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/ 
!
!
!
!
!    Cs - Rn*
!    =======
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=55,60) / &
!
!     Cs:
!
     132.905429E0_realk,100.000000E0_realk,  3.50000E0_realk,   2.582025E0_realk,  -0.00343E0_realk, &
 &     0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk, &
 &     0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk, &
 &     0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk, &
 &     0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk, &
 &     0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk, &
!
!     Ba:
!   
 &   137.905232E0_realk, 71.70000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &   136.905812E0_realk, 11.23000E0_realk,   1.50000E0_realk,  0.937365E0_realk,   0.245000E0_realk, &
 &   135.904553E0_realk,  7.85400E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &   134.905665E0_realk,  6.59200E0_realk,   1.50000E0_realk,  0.837943E0_realk,   0.160000E0_realk, &
 &   133.904486E0_realk,  2.41700E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &   131.905042E0_realk,  0.10100E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
!
!     La:
!
 &   138.906347E0_realk, 99.90980E0_realk,   3.50000E0_realk,  2.7830455E0_realk,  0.200000E0_realk, &
 &   137.907105E0_realk,  0.09020E0_realk,   5.00000E0_realk,  3.7136460E0_realk,  0.450000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
!
!     Ce
!
 &   139.905433E0_realk, 88.48000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk, &
 &   141.909241E0_realk, 11.08000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk, &
 &   137.905985E0_realk,  0.25000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk, &
 &   135.907140E0_realk,  0.19000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
! 
!     Pr:
!
 &   140.907647E0_realk,100.0000E0_realk,    2.50000E0_realk,  4.275400E0_realk,  -0.058900E0_realk, &
 &    0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &    0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &    0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &    0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &    0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
!
!     Nd:
!
 &   141.907719E0_realk, 27.130000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &   143.910083E0_realk, 23.800000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &   145.913113E0_realk, 17.190000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &   142.909810E0_realk, 12.180000E0_realk,   3.50000E0_realk, -0.065000E0_realk,  -0.630000E0_realk, &
 &   144.912570E0_realk,  8.300000E0_realk,   3.50000E0_realk, -1.065000E0_realk,  -0.330000E0_realk, &
 &   147.916889E0_realk,  5.760000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk/
!
!     Pm:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=61,64) / &
 &   144.912743E0_realk,100.000000E0_realk,   2.50000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
!
!     Sm:
!
 &   151.919728E0_realk, 26.700000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   153.922205E0_realk, 22.700000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   146.914894E0_realk, 15.000000E0_realk,   3.500000E0_realk,-0.814800E0_realk,  -0.259000E0_realk, &
 &   148.917180E0_realk, 13.800000E0_realk,   3.500000E0_realk,-0.671700E0_realk,   0.075000E0_realk, &
 &   147.914819E0_realk, 11.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   149.917273E0_realk,  7.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Eu:
!
 &   152.921225E0_realk, 52.200000E0_realk,   2.500000E0_realk, 1.533000E0_realk,   2.412000E0_realk, &
 &   150.919702E0_realk, 47.800000E0_realk,   2.500000E0_realk, 3.471700E0_realk,   0.903000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
 &     0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk, &
!
!     Gd:
!
 &   157.924019E0_realk, 24.840000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   159.927049E0_realk, 21.860000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   155.922118E0_realk, 20.470000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   156.923956E0_realk, 15.650000E0_realk,   1.500000E0_realk,-0.337260E0_realk,   1.350000E0_realk, &
 &   154.922618E0_realk, 14.800000E0_realk,   1.500000E0_realk,-0.257230E0_realk,   1.270000E0_realk, &
 &   153.920861E0_realk,  2.180000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/ 
!
!     Tb:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=65,70) / &
!
     158.925342E0_realk,100.000000E0_realk,   1.500000E0_realk, 2.014000E0_realk,   1.432000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Dy:
!
 &   163.929171E0_realk, 28.200000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   161.926795E0_realk, 25.500000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   162.928728E0_realk, 24.900000E0_realk,   2.500000E0_realk, 0.672600E0_realk,   2.648000E0_realk, &
 &   160.926930E0_realk, 18.900000E0_realk,   2.500000E0_realk,-0.480300E0_realk,   2.507000E0_realk, &
 &   159.925193E0_realk,  2.340000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   157.924277E0_realk,  0.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Ho:
!
 &   164.930319E0_realk,100.000000E0_realk,   3.500000E0_realk, 4.173000E0_realk,   3.580000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Er:
!
 &   165.930290E0_realk, 33.600000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   167.932368E0_realk, 26.800000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   166.932368E0_realk, 22.950000E0_realk,   3.500000E0_realk,-0.563850E0_realk,   3.565000E0_realk, &
 &   169.935461E0_realk, 14.900000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   163.929198E0_realk,  1.610000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   161.928775E0_realk,  0.140000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
! 
!     Tm:
!
 &   168.934212E0_realk,100.000000E0_realk,   0.500000E0_realk,-0.231600E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Yb:
!
 &   173.938859E0_realk, 31.800000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   171.936378E0_realk, 21.900000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   172.938208E0_realk, 16.120000E0_realk,   2.500000E0_realk,-0.679890E0_realk,   2.800000E0_realk, &
 &   170.936323E0_realk, 14.300000E0_realk,   0.500000E0_realk, 0.493670E0_realk,   0.000000E0_realk, &
 &   175.942564E0_realk, 12.700000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   169.934759E0_realk,  3.050000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/
!
!     Lu:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=71,74) / &
 &   174.940770E0_realk, 97.410000E0_realk,   3.500000E0_realk, 2.232700E0_realk,   3.490000E0_realk, &
 &   175.942679E0_realk,  2.590000E0_realk,   7.000000E0_realk, 3.169200E0_realk,   4.970000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Hf:
!
 &   179.9465457E0_realk,35.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   177.943696E0_realk, 27.297000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   176.943217E0_realk, 18.606000E0_realk,   3.500000E0_realk, 0.793500E0_realk,   3.365000E0_realk, &
 &   178.9458122E0_realk,13.629000E0_realk,   4.500000E0_realk,-0.640900E0_realk,   3.793000E0_realk, &
 &   175.941406E0_realk,  5.206000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   173.940044E0_realk,  0.162000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Ta:
!
 &   180.947992E0_realk, 99.988000E0_realk,   3.500000E0_realk, 2.370500E0_realk,   3.170000E0_realk, &
 &   179.947462E0_realk,  0.012000E0_realk,   8.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     W:
!
 &   183.950928E0_realk, 30.670000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   185.954357E0_realk, 28.600000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   181.948202E0_realk, 26.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   182.950928E0_realk, 14.300000E0_realk,   0.500000E0_realk, 0.11778476,   0.000000E0_realk, &
 &   179.947462E0_realk,  0.162000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/ 
!
!     Re:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=75,80) / &
!
     186.955744E0_realk, 62.600000E0_realk,   2.500000E0_realk, 3.219700E0_realk,   2.070000E0_realk, &
 &   184.952951E0_realk, 37.400000E0_realk,   2.500000E0_realk, 3.187100E0_realk,   2.180000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Os:
!
 &   191.961467E0_realk, 41.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   189.958436E0_realk, 26.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   188.958436E0_realk, 16.100000E0_realk,   1.500000E0_realk, 0.659933E0_realk,   0.856000E0_realk, &
 &   187.955830E0_realk, 13.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   186.955741E0_realk,  1.600000E0_realk,   0.500000E0_realk, 0.06465189E0_realk,   0.000000E0_realk, &
 &   185.953830E0_realk,  1.580000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Ir:
!
 &   192.962917E0_realk, 62.600000E0_realk,   1.500000E0_realk, 0.163700E0_realk,   0.751000E0_realk, &
 &   190.960584E0_realk, 37.400000E0_realk,   1.500000E0_realk, 0.150700E0_realk,   0.816000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Pt:
!
 &   194.964766E0_realk, 33.800000E0_realk,   0.500000E0_realk, 0.609520E0_realk,   0.000000E0_realk, &
 &   193.962655E0_realk, 32.900000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   195.964926E0_realk, 25.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   197.967869E0_realk,  7.200000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   191.961019E0_realk,  0.790000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   189.959917E0_realk,  0.010000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Au:
!   
 &   196.966543E0_realk,100.000000E0_realk,   1.500000E0_realk, 0.148158E0_realk,   0.547000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Hg:
!
 &   201.970617E0_realk, 29.860000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   199.968300E0_realk, 23.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   198.968254E0_realk, 16.870000E0_realk,   0.500000E0_realk, 0.50588549E0_realk,   0.000000E0_realk, &
 &   200.970277E0_realk, 13.180000E0_realk,   1.500000E0_realk,-0.5602257E0_realk ,   0.386000E0_realk, &
 &   197.966743E0_realk,  9.970000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   203.973467E0_realk,  6.870000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/ 
!
!     Tl:
!
      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=81,86) / &
!
     204.974401E0_realk, 70.476000E0_realk,   0.500000E0_realk, 1.63831461E0_realk, 0.000000E0_realk, &
 &   202.972320E0_realk, 29.524000E0_realk,   0.500000E0_realk, 1.62225787E0_realk, 0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
! 
!     Pb:
!   
 &   207.976627E0_realk, 52.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   205.975872E0_realk, 24.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &   206.975872E0_realk, 22.100000E0_realk,   0.500000E0_realk, 0.582583E0_realk,   0.000000E0_realk, &
 &   203.973020E0_realk,  1.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Bi:
!
 &   208.980374E0_realk,100.000000E0_realk,   4.500000E0_realk, 4.110600E0_realk,  -0.516000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Po:
!
 &   208.982404E0_realk,  0.000000E0_realk,   0.500000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     At:
!
 &   209.987126E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
!
!     Rn:
!
 &   222.017571E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk, &
 &     0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/ 
!
!
!
!
!
!
      IF (NIsotop .GT. MAXISO) THEN
         WRITE (LUPRI,'(//,A,2(/,A,I5),A)') ' NISOTOP too large in Isotops.', &
    &         ' Input value:  ',NIsotop, ' Maximum value:',MAXISO, &
    &         ' Program cannot continue.' 
         CALL LSQUIT('MAXISO exceeded in Isotops',lupri)
      END IF
      IF (Atom .GT. MAXCHR) THEN
         WRITE (LUPRI,'(//,A,2(/,A,I5),A)') &
    &         ' Atom too large in Isotops.', &
    &         ' Input value:  ', Atom, &
    &         ' Maximum value:',MAXCHR, &
    &         ' Program cannot continue.'
         CALL LSQUIT('MAXCHR exceeded in Isotops',lupri)
      END IF
!
      IF (Atom .LE. 0) THEN
!        This is a floating orbital, a point charge,
!        or an auxiliary basis set /Mar 2004 hjaaj
         Isotopes = D0
      ELSE IF (TYPE .EQ. 'MASS') THEN
         Isotopes = DATNUC(1,NIsotop,Atom)
      ELSE IF (TYPE .EQ. 'A') THEN
         Isotopes = NINT(DATNUC(1,NIsotop,Atom))
      ELSE IF (TYPE .EQ. 'ABUNDANCE') THEN
         Isotopes = DATNUC(2,NIsotop,Atom)
      ELSE IF (TYPE .EQ. 'SPIN') THEN
         Isotopes = DATNUC(3,NIsotop,Atom)
      ELSE IF (TYPE .EQ. 'MMOM') THEN
         Isotopes = DATNUC(4,NIsotop,Atom)
      ELSE IF (TYPE .EQ. 'GVAL') THEN
         SPIN = DATNUC(3,NIsotop,Atom)
         IF (SPIN .GT. THRESH) THEN
            Isotopes = DATNUC(4,NIsotop,Atom)/SPIN
         ELSE
            Isotopes = D0
         END IF
      ELSE IF (TYPE .EQ. 'LARMOR') THEN
         SPIN = DATNUC(3,NIsotop,Atom)
         IF (SPIN .GT. THRESH) THEN
           Isotopes = ABS(ECHARGE*DATNUC(4,NIsotop,Atom)/(D4*PI*SPIN*DMP))
         ELSE
           Isotopes = D0
         END IF
      ELSE IF (TYPE .EQ. 'QMOM') THEN
         Isotopes = DATNUC(5,NIsotop,Atom)
      ELSE IF (TYPE .EQ. 'NEUTRONS') THEN
         Isotopes = FLOAT(NINT(DATNUC(1,NIsotop,Atom)-Atom))
      ELSE
         WRITE (LUPRI,'(//,3A,/,A)') &
      &       ' Keyword ',TYPE,' unknown in Isotops.', &
      &       ' Program cannot continue.'
         CALL LSQUIT('Illegal keyword in Isotops',lupri)
      END IF
      RETURN
End function Isotopes
!=================!
! Covalent_Radius !       
!=================!
Function CovRad(Charge,lupri)
!     Based on covalent radii and metallic radii in Angstrom.
!     Returns -1 where data is inavailable
Use precision
IMPLICIT NONE
REAL(realk) :: CovRad
INTEGER :: Charge,lupri
REAL(realk), DIMENSION(100) :: RAD
INTEGER :: I
      DATA (RAD(I), I = 1, 100)/ &
             30.0E0_realk,  155.0E0_realk,  160.0E0_realk,  110.0E0_realk, &
      90.0E0_realk,   80.0E0_realk,   70.0E0_realk,   68.0E0_realk,   65.0E0_realk, &
     154.0E0_realk,  190.0E0_realk,  160.0E0_realk,  140.0E0_realk,  110.0E0_realk, &
     110.0E0_realk,  105.0E0_realk,  105.0E0_realk,  190.0E0_realk,  238.0E0_realk, &
     200.0E0_realk,  165.0E0_realk,  145.0E0_realk,  135.0E0_realk,  130.0E0_realk, &
     125.0E0_realk,  125.0E0_realk,  125.0E0_realk,  125.0E0_realk,  125.0E0_realk, &
     140.0E0_realk,  140.0E0_realk,  130.0E0_realk,  120.0E0_realk,  120.0E0_realk, &
     120.0E0_realk,  200.0E0_realk,  255.0E0_realk,  215.0E0_realk,  180.0E0_realk, &
     160.0E0_realk,  145.0E0_realk,  140.0E0_realk,  135.0E0_realk,  130.0E0_realk, &
     130.0E0_realk,  135.0E0_realk,  140.0E0_realk,  155.0E0_realk,  160.0E0_realk, &
     160.0E0_realk,  140.0E0_realk,  140.0E0_realk,  140.0E0_realk,  220.0E0_realk, &
     270.0E0_realk,  220.0E0_realk,  185.0E0_realk,  180.0E0_realk,  180.0E0_realk, &
     180.0E0_realk,  180.0E0_realk,  180.0E0_realk,  200.0E0_realk,  180.0E0_realk, &
     175.0E0_realk,  175.0E0_realk,  175.0E0_realk,  175.0E0_realk,  170.0E0_realk, &
     170.0E0_realk,  170.0E0_realk,  155.0E0_realk,  145.0E0_realk,  140.0E0_realk, &
     135.0E0_realk,  135.0E0_realk,  135.0E0_realk,  135.0E0_realk,  145.0E0_realk, &
     155.0E0_realk,  170.0E0_realk,  175.0E0_realk,  170.0E0_realk,   -100.0E0_realk, &
      -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk, &
      -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk, &
      -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk,   -100.0E0_realk, &
      -100.0E0_realk/
!
      IF (CHARGE .EQ. 0) THEN
!hj      solvent cavity center or floating orbital
        CovRad = -1.0E0_realk
      ELSE IF (CHARGE .LT. 1 .OR. CHARGE .GT. 100) THEN
         WRITE (LUPRI,*) 'ERROR, RADIUS called with NCHARGE =',CHARGE
         CALL LSQUIT('RADIUS called with illegal NCHARGE',-1)
      ELSE
         CovRad = 0.01E0_realk * RAD(CHARGE)
      END IF
End function CovRad
!
End module Fundamental
