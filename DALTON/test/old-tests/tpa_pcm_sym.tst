########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm tpa hf formaldehyde short
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enehf
surf
nuc
nucchg
tes
qrlrve
qrlrve2
omegab
sym
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN RESPONSE
*PCM
.SOLVNT
WATER
.NPCMMT
0
.NESFP
1
.ICESPH
1
*PCMCAV
.CENTER
0.000 0.000 0.000
.RIN
6.00
.AREATS
0.3
.ALPHA
0.529177
**WAVEFUNCTION
.HF
*SCF INPUT
.THRESHOLD
1.0D-8
**RESPONSE
*QUADRATIC
.DIPLEN
.TWO-PHOTON
.ROOTS
 3 3 3 3 
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS


Atomtypes=3 Charge=+2 Angstrom Generators=2 X Y
Charge=6.0 Atoms=1 Basis=STO-3G
C     0.000000     0.000000     0.000000
Charge=8.0 Atoms=1 Basis=STO-3G
O     0.000000     0.000000     1.220000
Charge=1.0 Atoms=1 Basis=STO-3G
H     0.943102     0.000000    -0.544500
END MOLINP

########## Reference Output ########################
START REFOUT


     ************************************************************************
     *************** Dalton - An Electronic Structure Program ***************
     ************************************************************************

    This is output from DALTON 2013.2
   ----------------------------------------------------------------------------
    NOTE:
     
    Dalton is an experimental code for the evaluation of molecular
    properties using (MC)SCF, DFT, CI, and CC wave functions.
    The authors accept no responsibility for the performance of
    the code or for the correctness of the results.
     
    The code (in whole or part) is provided under a licence and
    is not to be reproduced for further distribution without
    the written permission of the authors or their representatives.
     
    See the home page "http://daltonprogram.org" for further information.
     
    If results obtained with this code are published,
    the appropriate citations would be both of:
     
       K. Aidas, C. Angeli, K. L. Bak, V. Bakken, R. Bast,
       L. Boman, O. Christiansen, R. Cimiraglia, S. Coriani,
       P. Dahle, E. K. Dalskov, U. Ekstroem, T. Enevoldsen,
       J. J. Eriksen, P. Ettenhuber, B. Fernandez, L. Ferrighi,
       H. Fliegl, L. Frediani, K. Hald, A. Halkier, C. Haettig,
       H. Heiberg, T. Helgaker, A. C. Hennum, H. Hettema,
       E. Hjertenaes, S. Hoest, I.-M. Hoeyvik, M. F. Iozzi,
       B. Jansik, H. J. Aa. Jensen, D. Jonsson, P. Joergensen,
       J. Kauczor, S. Kirpekar, T. Kjaergaard, W. Klopper,
       S. Knecht, R. Kobayashi, H. Koch, J. Kongsted, A. Krapp,
       K. Kristensen, A. Ligabue, O. B. Lutnaes, J. I. Melo,
       K. V. Mikkelsen, R. H. Myhre, C. Neiss, C. B. Nielsen,
       P. Norman, J. Olsen, J. M. H. Olsen, A. Osted,
       M. J. Packer, F. Pawlowski, T. B. Pedersen, P. F. Provasi,
       S. Reine, Z. Rinkevicius, T. A. Ruden, K. Ruud, V. Rybkin,
       P. Salek, C. C. M. Samson, A. Sanchez de Meras, T. Saue,
       S. P. A. Sauer, B. Schimmelpfennig, K. Sneskov,
       A. H. Steindal, K. O. Sylvester-Hvid, P. R. Taylor,
       A. M. Teale, E. I. Tellgren, D. P. Tew, A. J. Thorvaldsen,
       L. Thoegersen, O. Vahtras, M. A. Watson, D. J. D. Wilson,
       M. Ziolkowski and H. Aagren,
       "The Dalton quantum chemistry program system",
       WIREs Comput. Mol. Sci. 2013. (doi: 10.1002/wcms.1172)
    
    and
    
       Dalton, a Molecular Electronic Structure Program,
       Release DALTON2013.2 (2013), see http://daltonprogram.org
   ----------------------------------------------------------------------------

    Authors in alphabetical order (major contribution(s) in parenthesis):

  Kestutis Aidas,           Vilnius University,           Lithuania   (QM/MM)
  Celestino Angeli,         University of Ferrara,        Italy       (NEVPT2)
  Keld L. Bak,              UNI-C,                        Denmark     (AOSOPPA, non-adiabatic coupling, magnetic properties)
  Vebjoern Bakken,          University of Oslo,           Norway      (DALTON; geometry optimizer, symmetry detection)
  Radovan Bast,             KTH Stockholm,                Sweden      (DALTON installation and execution frameworks)
  Linus Boman,              NTNU,                         Norway      (Cholesky decomposition and subsystems)
  Ove Christiansen,         Aarhus University,            Denmark     (CC module)
  Renzo Cimiraglia,         University of Ferrara,        Italy       (NEVPT2)
  Sonia Coriani,            University of Trieste,        Italy       (CC module, MCD in RESPONS)
  Paal Dahle,               University of Oslo,           Norway      (Parallelization)
  Erik K. Dalskov,          UNI-C,                        Denmark     (SOPPA)
  Thomas Enevoldsen,        Univ. of Southern Denmark,    Denmark     (SOPPA)
  Janus J. Eriksen,         Aarhus University,            Denmark     (PE-MP2/SOPPA, TDA)
  Berta Fernandez,          U. of Santiago de Compostela, Spain       (doublet spin, ESR in RESPONS)
  Lara Ferrighi,            Aarhus University,            Denmark     (PCM Cubic response)
  Heike Fliegl,             University of Oslo,           Norway      (CCSD(R12))
  Luca Frediani,            UiT The Arctic U. of Norway,  Norway      (PCM)
  Bin Gao,                  UiT The Arctic U. of Norway,  Norway      (Gen1Int library)
  Christof Haettig,         Ruhr-University Bochum,       Germany     (CC module)
  Kasper Hald,              Aarhus University,            Denmark     (CC module)
  Asger Halkier,            Aarhus University,            Denmark     (CC module)
  Hanne Heiberg,            University of Oslo,           Norway      (geometry analysis, selected one-electron integrals)
  Trygve Helgaker,          University of Oslo,           Norway      (DALTON; ABACUS, ERI, DFT modules, London, and much more)
  Alf Christian Hennum,     University of Oslo,           Norway      (Parity violation)
  Hinne Hettema,            University of Auckland,       New Zealand (quadratic response in RESPONS; SIRIUS supersymmetry)
  Eirik Hjertenaes,         NTNU,                         Norway      (Cholesky decomposition)
  Maria Francesca Iozzi,    University of Oslo,           Norway      (RPA)
  Brano Jansik              Technical Univ. of Ostrava    Czech Rep.  (DFT cubic response)
  Hans Joergen Aa. Jensen,  Univ. of Southern Denmark,    Denmark     (DALTON; SIRIUS, RESPONS, ABACUS modules, London, and much more)
  Dan Jonsson,              UiT The Arctic U. of Norway,  Norway      (cubic response in RESPONS module)
  Poul Joergensen,          Aarhus University,            Denmark     (RESPONS, ABACUS, and CC modules)
  Joanna Kauczor,           Linkoeping University,        Sweden      (Complex polarization propagator (CPP) module)
  Sheela Kirpekar,          Univ. of Southern Denmark,    Denmark     (Mass-velocity & Darwin integrals)
  Wim Klopper,              KIT Karlsruhe,                Germany     (R12 code in CC, SIRIUS, and ABACUS modules)
  Stefan Knecht,            ETH Zurich,                   Switzerland (Parallel CI and MCSCF)
  Rika Kobayashi,           Australian National Univ.,    Australia   (DIIS in CC, London in MCSCF)
  Henrik Koch,              NTNU,                         Norway      (CC module, Cholesky decomposition)
  Jacob Kongsted,           Univ. of Southern Denmark,    Denmark     (Polarizable embedding, QM/MM)
  Andrea Ligabue,           University of Modena,         Italy       (CTOCD, AOSOPPA)
  Ola B. Lutnaes,           University of Oslo,           Norway      (DFT Hessian)
  Juan I. Melo,             University of Buenos Aires,   Argentina   (LRESC, Relativistic Effects on NMR Shieldings)
  Kurt V. Mikkelsen,        University of Copenhagen,     Denmark     (MC-SCRF and QM/MM)
  Rolf H. Myhre,            NTNU,                         Norway      (Cholesky, subsystems and ECC2)
  Christian Neiss,          Univ. Erlangen-Nuernberg,     Germany     (CCSD(R12))
  Christian B. Nielsen,     University of Copenhagen,     Denmark     (QM/MM)
  Patrick Norman,           Linkoeping University,        Sweden      (Cubic response and complex response in RESPONS)
  Jeppe Olsen,              Aarhus University,            Denmark     (SIRIUS CI/density modules)
  Jogvan Magnus H. Olsen,   Univ. of Southern Denmark,    Denmark     (Polarizable embedding, PE library, QM/MM)
  Anders Osted,             Copenhagen University,        Denmark     (QM/MM)
  Martin J. Packer,         University of Sheffield,      UK          (SOPPA)
  Filip Pawlowski,          Kazimierz Wielki University,  Poland      (CC3)
  Thomas B. Pedersen,       University of Oslo,           Norway      (Cholesky decomposition)
  Patricio F. Provasi,      University of Northeastern,   Argentina   (Analysis of coupling constants in localized orbitals)
  Zilvinas Rinkevicius,     KTH Stockholm,                Sweden      (open-shell DFT, ESR)
  Elias Rudberg,            KTH Stockholm,                Sweden      (DFT grid and basis info)
  Torgeir A. Ruden,         University of Oslo,           Norway      (Numerical derivatives in ABACUS)
  Kenneth Ruud,             UiT The Arctic U. of Norway,  Norway      (DALTON; ABACUS magnetic properties and  much more)
  Pawel Salek,              KTH Stockholm,                Sweden      (DALTON; DFT code)
  Claire C. M. Samson       University of Karlsruhe       Germany     (Boys localization, r12 integrals in ERI)
  Alfredo Sanchez de Meras, University of Valencia,       Spain       (CC module, Cholesky decomposition)
  Trond Saue,               Paul Sabatier University,     France      (direct Fock matrix construction)
  Stephan P. A. Sauer,      University of Copenhagen,     Denmark     (SOPPA(CCSD), SOPPA prop., AOSOPPA, vibrational g-factors)
  Bernd Schimmelpfennig,    Forschungszentrum Karlsruhe,  Germany     (AMFI module)
  Kristian Sneskov,         Aarhus University,            Denmark     (QM/MM, PE-CC)
  Arnfinn H. Steindal,      UiT The Arctic U. of Norway,  Norway      (parallel QM/MM)
  K. O. Sylvester-Hvid,     University of Copenhagen,     Denmark     (MC-SCRF)
  Peter R. Taylor,          VLSCI/Univ. of Melbourne,     Australia   (Symmetry handling ABACUS, integral transformation)
  Andrew M. Teale,          University of Nottingham,     England     (DFT-AC, DFT-D)
  David P. Tew,             University of Bristol,        England     (CCSD(R12))
  Olav Vahtras,             KTH Stockholm,                Sweden      (triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)
  David J. Wilson,          La Trobe University,          Australia   (DFT Hessian and DFT magnetizabilities)
  Hans Agren,               KTH Stockholm,                Sweden      (SIRIUS module, RESPONS, MC-SCRF solvation model)
 --------------------------------------------------------------------------------

     Date and time (Linux)  : Wed Apr  9 11:47:10 2014
     Host name              : rfb-ntf-m0022.bargi.uit.no              

 * Work memory size             :    64000000 =  488.28 megabytes.

 * Directories for basis set searches:
   1) /Users/ast018/dalton/DALTON/test/perl-pid.68818__2014_4_9__11.47
   2) /Users/ast018/dalton/build-release/basis


Compilation information
-----------------------

 Who compiled             | ast018
 Host                     | rfb-ntf-m0022.bargi.uit.no
 System                   | Darwin-13.1.0
 CMake generator          | Unix Makefiles
 Processor                | i386
 64-bit integers          | OFF
 MPI                      | OFF
 Fortran compiler         | /opt/local/bin/gfortran
 Fortran compiler version | GNU Fortran (MacPorts gcc48 4.8.2_0) 4.8.2
 C compiler               | /opt/local/bin/gcc
 C compiler version       | gcc (MacPorts gcc48 4.8.2_0) 4.8.2
 C++ compiler             | /opt/local/bin/g++
 C++ compiler version     | g++ (MacPorts gcc48 4.8.2_0) 4.8.2
 BLAS                     | /usr/lib/libblas.dylib
 LAPACK                   | /usr/lib/liblapack.dylib
 Static linking           | OFF
 Last Git revision        | 4dbbabb174f0d7d8df285d90aa7d832f95fee97e
 Configuration time       | 2014-04-08 10:42:03.275960


   Content of the .dal input file
 ----------------------------------

**DALTON INPUT                                    
.RUN RESPONSE                                     
*PCM                                              
.SOLVNT                                           
WATER                                             
.NPCMMT                                           
0                                                 
.NESFP                                            
1                                                 
.ICESPH                                           
1                                                 
*PCMCAV                                           
.CENTER                                           
0.000 0.000 0.000                                 
.RIN                                              
6.00                                              
.AREATS                                           
0.3                                               
.ALPHA                                            
0.529177                                          
**WAVEFUNCTION                                    
.HF                                               
*SCF INPUT                                        
.THRESHOLD                                        
1.0e-8                                            
**RESPONSE                                        
*QUADRATIC                                        
.DIPLEN                                           
.TWO-PHOTON                                       
.ROOTS                                            
 3 3 3 3                                          
**END OF DALTON INPUT                             


   Content of the .mol file
 ----------------------------

ATOMBASIS                                                                      
                                                                               
                                                                               
Atomtypes=3 Charge=+2 Angstrom Generators=2 X Y                                
Charge=6.0 Atoms=1 Basis=STO-3G                                                
C     0.000000     0.000000     0.000000                                       
Charge=8.0 Atoms=1 Basis=STO-3G                                                
O     0.000000     0.000000     1.220000                                       
Charge=1.0 Atoms=1 Basis=STO-3G                                                
H     0.943102     0.000000    -0.544500                                       


       *******************************************************************
       *********** Output from DALTON general input processing ***********
       *******************************************************************

 --------------------------------------------------------------------------------
   Overall default print level:    0
   Print level for DALTON.STAT:    1

    HERMIT 1- and 2-electron integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed (SIRIUS module)
    Dynamic molecular response properties section will be executed (RESPONSE module)
 --------------------------------------------------------------------------------

 ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT = WATER    **
 Optical and physical constants:
 EPS= 78.390; EPSINF=  1.776; RSOLV=  1.385 A; VMOL=  18.070 ML/MOL;
 TCE=2.57000e-04 1/K; STEN= 71.810 DYN/CM;  DSTEN=  0.6500; CMF=  1.2770


     -----------------------------------
     Input for PCM solvation calculation 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=WATER        EPS   = 78.3900     EPSINF=  1.7760
     RSOLV =  1.3850

     ICESPH =       1     NESFP =       1
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
     POLYG          60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    6.0000


   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************


    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1:                                                                         
 2:                                                                         
    ------------------------------------------------------------------------

  Coordinates are entered in Angstrom and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   6 :
     "/Users/ast018/dalton/build-release/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   8 :
     "/Users/ast018/dalton/build-release/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   1 :
     "/Users/ast018/dalton/build-release/basis/STO-3G"


                         SYMGRP: Point group information
                         -------------------------------

Point group: C2v

   * The point group was generated by:

      Reflection in the yz-plane
      Reflection in the xz-plane

   * Group multiplication table

        |  E   C2z  Oxz  Oyz
   -----+--------------------
     E  |  E   C2z  Oxz  Oyz
    C2z | C2z   E   Oyz  Oxz
    Oxz | Oxz  Oyz   E   C2z
    Oyz | Oyz  Oxz  C2z   E 

   * Character table

        |  E   C2z  Oxz  Oyz
   -----+--------------------
    A1  |   1    1    1    1
    B1  |   1   -1    1   -1
    B2  |   1   -1   -1    1
    A2  |   1    1   -1   -1

   * Direct product table

        | A1   B1   B2   A2 
   -----+--------------------
    A1  | A1   B1   B2   A2 
    B1  | B1   A1   A2   B2 
    B2  | B2   A2   A1   B1 
    A2  | A2   B2   B1   A1 
PCMSPHGEN: SPHERE CENTERS FROM INPUT
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000e+00    0.0000000000e+00    0.0000000000e+00    6.0000000000e+00


                                 Isotopic Masses
                                 ---------------

                           C          12.000000
                           O          15.994915
                           H   _1      1.007825
                           H   _2      1.007825

                       Total mass:    30.010565 amu
                       Natural abundance:  98.633 %

 Center-of-mass coordinates (a.u.):    0.000000    0.000000    1.159649


  Atoms and basis sets
  --------------------

  Number of atom types :    3
  Total number of atoms:    4

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C           1    6.0000    15     5      [6s3p|2s1p]                                        
  O           1    8.0000    15     5      [6s3p|2s1p]                                        
  H           2    1.0000     3     1      [3s|1s]                                            
  ----------------------------------------------------------------------
  total:      4   16.0000    36    12
  ----------------------------------------------------------------------

  Threshold for neglecting AO integrals:  1.00e-12


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   12
  C       :     1  x   0.0000000000    2  y   0.0000000000    3  z   0.0000000000
  O       :     4  x   0.0000000000    5  y   0.0000000000    6  z   2.3054658725
  H   / 1 :     7  x   1.7822044879    8  y   0.0000000000    9  z  -1.0289558751
  H   / 2 :    10  x  -1.7822044879   11  y   0.0000000000   12  z  -1.0289558751


  Symmetry Coordinates
  --------------------

  Number of coordinates in each symmetry:     4    4    3    1

  Symmetry  A1  ( 1)

    1   C     z    3
    2   O     z    6
    3   H     x    [  7  -   10 ]/2
    4   H     z    [  9  +   12 ]/2

  Symmetry  B1  ( 2)

    5   C     x    1
    6   O     x    4
    7   H     x    [  7  +   10 ]/2
    8   H     z    [  9  -   12 ]/2

  Symmetry  B2  ( 3)

    9   C     y    2
   10   O     y    5
   11   H     y    [  8  +   11 ]/2

  Symmetry  A2  ( 4)

   12   H     y    [  8  -   11 ]/2


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H   _1      H   _2
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H   _1:    1.089000    2.000725    0.000000
 H   _2:    1.089000    2.000725    1.886204    0.000000


  Max    interatomic separation is    2.0007 Angstrom (    3.7808 Bohr)
  between atoms    3 and    2, "H   _1" and "O     ".

  Min HX interatomic separation is    1.0890 Angstrom (    2.0579 Bohr)

  Min YX interatomic separation is    1.2200 Angstrom (    2.3055 Bohr)


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.220000
  bond distance:  H   _1     C            1.089000
  bond distance:  H   _2     C            1.089000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H   _1       120.000
  bond angle:     O          C          H   _2       120.000
  bond angle:     H   _1     C          H   _2       120.000




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       1.792803          0.000000    0.000000    1.000000
   IB      13.103106          1.000000    0.000000    0.000000
   IC      14.895908          0.000000    1.000000    0.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         281893.2926          38569.4058          33927.3708 MHz
            9.402948            1.286537            1.131695 cm-1


@  Nuclear repulsion energy :   31.163673729192 Hartree


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:           7    3    2    0


  Symmetry  A1 ( 1)

    1     C        1s         1
    2     C        1s         2
    3     C        2pz        5
    4     O        1s         6
    5     O        1s         7
    6     O        2pz       10
    7     H        1s        11 +   12


  Symmetry  B1 ( 2)

    8     C        2px        3
    9     O        2px        8
   10     H        1s        11 -   12


  Symmetry  B2 ( 3)

   11     C        2py        4
   12     O        2py        9


  No orbitals in symmetry  A2 ( 4)

  Symmetries of electric field:  B1 (2)  B2 (3)  A1 (1)

  Symmetries of magnetic field:  B2 (3)  B1 (2)  A2 (4)


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************



     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************


 Threshold for neglecting two-electron integrals:  1.00e-12
 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014


 MEMORY USED TO GENERATE CAVITY =    432202


 Total number of spheres =    1
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1    0.000000000    0.000000000    0.000000000    3.175062000  126.681817203

 Total number of tesserae =     392
 Surface area =  126.68181720 (A^2)    Cavity volume =  134.07420796 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....

  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   0.07 seconds
 >>>> Total wall time used in HERMIT:   0.09 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          F
     Using EHT  model:          T
     Number of Huckel orbitals each symmetry:    7    3    2    0

 Huckel EHT eigenvalues for symmetry :  1
          -20.704416     -11.377295      -1.495729      -0.964135      -0.593338
           -0.252244      -0.199950

 Huckel EHT eigenvalues for symmetry :  2
           -0.744095      -0.506268      -0.190531

 Huckel EHT eigenvalues for symmetry :  3
           -0.670977      -0.352123

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Wed Apr  9 11:47:11 2014
     Host name              : rfb-ntf-m0022.bargi.uit.no              

 Title lines from ".mol" input file:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

@    Restricted, closed shell Hartree-Fock calculation.

@    Time-dependent Hartree-Fock calculation (random phase approximation).

 Initial molecular orbitals are obtained according to
 ".MOSTART EHT   " input option

@    QM part is embedded in an environment :

@         Solvation model: PCM

     Wave function specification
     ============================
@    Wave function type        >>> HF <<<
@    Number of closed shell electrons          14
@    Number of electrons in active shells       0
@    Total charge of the molecule               2

@    Spin multiplicity and 2 M_S                1         0
@    Total number of symmetries                 4 (point group: C2v)
@    Reference state symmetry                   1 (irrep name : A1 )

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1    2    3    4
                                           |  A1   B1   B2   A2 
                                       --- |  ---  ---  ---  ---
     Total number of orbitals           12 |    7    3    2    0
     Number of basis functions          12 |    7    3    2    0

      ** Automatic occupation of RHF orbitals **

      -- Initial occupation of symmetries is determined from extended Huckel guess.           
      -- Initial occupation of symmetries is :
@    Occupied SCF orbitals               7 |    5    1    1    0

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-08

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.99669
 NUCLEAR APPARENT CHARGE -15.79079
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =    8

 Automatic occupation of symmetries with  14 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    SCF occupation
 -----------------------------------------------------------------------------
     (Precalculated two-electron integrals are transformed to P-supermatrix elements.
      Threshold for discarding integrals :  1.00e-12 )
@  1  -111.294142892     -0.329342268591       2.30e+00  -1.11e+02     5   1   1   0
      Virial theorem: -V/T =      2.030006
@      MULPOP C       1.21; O       0.46; H   _1  0.17; H   _2  0.17; 
 -----------------------------------------------------------------------------
@  2  -111.514696813     -0.333510590827       7.51e-01  -2.21e-01     5   1   1   0
      Virial theorem: -V/T =      2.035325
@      MULPOP C       0.07; O       0.92; H   _1  0.50; H   _2  0.50; 
 -----------------------------------------------------------------------------
@  3  -111.543470257     -0.329772774956       4.39e-01  -2.88e-02     5   1   1   0
      Virial theorem: -V/T =      2.028935
@      MULPOP C       0.56; O       0.39; H   _1  0.52; H   _2  0.52; 
 -----------------------------------------------------------------------------
@  4  -111.560873009     -0.330584549917       1.10e-01  -1.74e-02     5   1   1   0
      Virial theorem: -V/T =      2.030983
@      MULPOP C       0.36; O       0.60; H   _1  0.52; H   _2  0.52; 
 -----------------------------------------------------------------------------
@  5  -111.565569920     -0.330603333656       3.69e-02  -4.70e-03     5   1   1   0
      Virial theorem: -V/T =      2.030395
@      MULPOP C       0.38; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@  6  -111.566150471     -0.330638457315       4.26e-03  -5.81e-04     5   1   1   0
      Virial theorem: -V/T =      2.030171
@      MULPOP C       0.39; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@  7  -111.566155773     -0.330640786916       1.33e-03  -5.30e-06     5   1   1   0
      Virial theorem: -V/T =      2.030214
@      MULPOP C       0.39; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@  8  -111.566156148     -0.330641269667       4.23e-04  -3.76e-07     5   1   1   0
      Virial theorem: -V/T =      2.030202
@      MULPOP C       0.39; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@  9  -111.566156221     -0.330640945191       9.66e-06  -7.23e-08     5   1   1   0
      Virial theorem: -V/T =      2.030205
@      MULPOP C       0.39; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@ 10  -111.566156221     -0.330640961499       2.37e-06  -1.58e-11     5   1   1   0
      Virial theorem: -V/T =      2.030205
@      MULPOP C       0.39; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@ 11  -111.566156221     -0.330640964889       1.55e-07  -1.99e-12     5   1   1   0
      Virial theorem: -V/T =      2.030205
@      MULPOP C       0.39; O       0.59; H   _1  0.51; H   _2  0.51; 
 -----------------------------------------------------------------------------
@ 12  -111.566156221     -0.330640964932       7.16e-09   1.42e-14     5   1   1   0

@ *** DIIS converged in  12 iterations !
@     Converged SCF energy, gradient:   -111.566156220654    7.16e-09
    - total time used in SIRFCK :              0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   14
 Orbital occupations :    5    1    1    0

 Sym       Hartree-Fock orbital energies

1 A1    -21.24419236   -11.75700806    -1.96710566    -1.28470073    -1.11212502
          0.20015224     0.27610554

2 B1     -1.09098267    -0.44093922     0.29510714

3 B2     -1.08694551    -0.26242749

    E(LUMO) :    -0.44093922 au (symmetry 2)
  - E(HOMO) :    -1.08694551 au (symmetry 3)
  ------------------------------------------
    gap     :     0.64600629 au

 >>> Writing SIRIFC interface file <<<

 >>>> CPU and wall time for SCF :       1.022       1.218


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


@    Spin multiplicity:           1
@    Spatial symmetry:            1 ( irrep  A1  in C2v )
@    Total charge of molecule:    2

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

@    Final HF energy:            -111.566156220654                 
@    Nuclear repulsion:            31.163673729192
@    Electronic energy:          -142.399188984914

@    Final gradient norm:           0.000000007163

 
     Date and time (Linux)  : Wed Apr  9 11:47:12 2014
     Host name              : rfb-ntf-m0022.bargi.uit.no              

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species 1  (A1 )
 ------------------------------------------------

    Orbital         1        2        3        4        5        6        7
   1 C   :1s     0.0005   0.9929  -0.1100  -0.1963  -0.0602   0.2179  -0.0305
   2 C   :1s    -0.0068   0.0318   0.2171   0.6406   0.2459  -1.4199   0.1822
   3 C   :2pz   -0.0058   0.0020   0.1441  -0.0947  -0.4750  -0.0428   1.2445
   4 O   :1s     0.9947   0.0002  -0.2211   0.1213  -0.0622  -0.0595   0.0949
   5 O   :1s     0.0242  -0.0061   0.7939  -0.5610   0.3266   0.4288  -0.7420
   6 O   :2pz   -0.0055   0.0019  -0.2460  -0.3804   0.6286  -0.4932   0.7622
   7 H   :1s     0.0003  -0.0061   0.0185   0.1501   0.1698   0.8242   0.4352

 Molecular orbitals for symmetry species 2  (B1 )
 ------------------------------------------------

    Orbital         1        2        3
   1 C   :2px    0.6742  -0.1767   1.1020
   2 O   :2px    0.3443   0.9095  -0.3278
   3 H   :1s     0.2176  -0.3096  -0.8999

 Molecular orbitals for symmetry species 3  (B2 )
 ------------------------------------------------

    Orbital         1        2
   1 C   :2py    0.3159   0.9721
   2 O   :2py    0.8855  -0.5106



 >>>> Total CPU  time used in SIRIUS :      1.03 seconds
 >>>> Total wall time used in SIRIUS :      1.25 seconds

 
     Date and time (Linux)  : Wed Apr  9 11:47:12 2014
     Host name              : rfb-ntf-m0022.bargi.uit.no              


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



                 .------------------------------------------------.
                 | Starting in Dynamic Property Section (RESPONS) |
                 `------------------------------------------------'


 ------------------------------------------------------------------------------
  RESPONSE  -  an MCSCF, MC-srDFT, DFT, and SOPPA response property program
 ------------------------------------------------------------------------------


 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




 Quadratic Response single residue calculation
 ---------------------------------------------

 That is, calculation of second order moments.

 Two-photon transition processes computed (.TWO-PHOTON).
 Both photon frequencies are half the excitation energy.

 Spin of operator A , ISPINA=    0
 Spin of operator B , ISPINB=    0
 Spin of operator C , (Excitation energy) ISPINC=    0

 Print level                                    : IPRSMO =    2
 Threshold for convergence in RSPPP             : THCPP  = 1.000e-03
 Maximum number of iterations in RSPPP          : MAXITP =   60
 Threshold for convergence in RSPLR             : THCLR  = 1.000e-03
 Maximum number of iterations in RSPLR          : MAXITL =   60
 Maximum iterations in optimal orbital algorithm: MAXITO =    5

    1 A OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          ZDIPLEN 

    1 A OPERATORS OF SYMMETRY NO:    2 AND LABELS:

          XDIPLEN 

    1 A OPERATORS OF SYMMETRY NO:    3 AND LABELS:

          YDIPLEN 

    1 B OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          ZDIPLEN 

    1 B OPERATORS OF SYMMETRY NO:    2 AND LABELS:

          XDIPLEN 

    1 B OPERATORS OF SYMMETRY NO:    3 AND LABELS:

          YDIPLEN 

 Integral transformation: Total CPU and WALL times (sec)       0.002       0.005


   SCF energy         :     -111.566156220653525
 -- inactive part     :     -142.399188984914048
 -- nuclear repulsion :       31.163673729192165


                     ***************************************
                     *** RHF response calculation (TDHF) ***
                     ***************************************



 Linear response excitations for quadratic response
 - symmetry of excitation operator    1  ( A1 )
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      13
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      13



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry = 1  ( A1 ); triplet =   F


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   22)
 RSP solution vector no.    1; norm of residual   3.71e-05
 RSP solution vector no.    2; norm of residual   1.31e-05
 RSP solution vector no.    3; norm of residual   9.90e-05

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    2  ( B1 )
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry = 2  ( B1 ); triplet =   F


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   20)
 RSP solution vector no.    1; norm of residual   4.20e-04
 RSP solution vector no.    2; norm of residual   2.08e-04
 RSP solution vector no.    3; norm of residual   1.22e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    3  ( B2 )
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry = 3  ( B2 ); triplet =   F

RSPORT:    2 out of    3 new trial vectors linear dependent

 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   14)
 RSP solution vector no.    1; norm of residual   1.01e-14
 RSP solution vector no.    2; norm of residual   7.30e-15
 RSP solution vector no.    3; norm of residual   3.72e-15

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    4  ( A2 )
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       4
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       3
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       3



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry = 4  ( A2 ); triplet =   F


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   9.12e-16
 RSP solution vector no.    2; norm of residual   6.68e-16
 RSP solution vector no.    3; norm of residual   1.38e-15

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    1  ( A1 )


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      13
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      13


 QRLRVE -- linear response calculation for symmetry  1  ( A1 )
 QRLRVE -- operator label : ZDIPLEN 
 QRLRVE -- operator spin  :   0
 QRLRVE -- frequencies :  0.156891  0.260466  0.472415  0.077257  0.255638
                          0.432130  0.206986  0.331146  0.396124



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A1 ); triplet =   F


 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after    9 linear transformations is     19.68121224
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.01 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after    9 linear transformations is     13.03117779
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.68 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.33115 a.u.
 after    9 linear transformations is    -13.05960917
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -8.05 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.39612 a.u.
 after    9 linear transformations is      9.11183954
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.65 * 10 **  6.0
RSPORT:    5 out of    9 new trial vectors linear dependent

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after   13 linear transformations is     26.40470617
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -8.47 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after   13 linear transformations is     16.01891747
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -1.31 * 10 ** 16.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.33115 a.u.
 after   13 linear transformations is     -7.93150482
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -4.94 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  1,  triplet       F

 GP * SOLUTION vector at frequency     0.39612 a.u.
 after   13 linear transformations is     11.27102370
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -1.40 * 10 ** 16.0

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   26)
 RSP solution vector no.    1; norm of residual   1.05e-15
 RSP solution vector no.    2; norm of residual   5.21e-16
 RSP solution vector no.    3; norm of residual   3.93e-16
 RSP solution vector no.    4; norm of residual   8.23e-16
 RSP solution vector no.    5; norm of residual   5.80e-16
 RSP solution vector no.    6; norm of residual   4.48e-16
 RSP solution vector no.    7; norm of residual   1.14e-15
 RSP solution vector no.    8; norm of residual   6.11e-16
 RSP solution vector no.    9; norm of residual   5.01e-16

 *** RSPCTL MICROITERATIONS CONVERGED

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.156891e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.15689):     12.3942663815    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.260466e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.26047):     17.8351255843    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.472415e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.47241):     26.4047061652    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.772568e-01

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.07726):     11.3829222937    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.255638e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.25564):     17.2063833101    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.432130e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.43213):     16.0189174699    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.206986e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.20699):     13.8309993149    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.331146e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.33115):    -7.93150481524    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A1 )    LABEL   ZDIPLEN     FREQUENCY   0.396124e+00

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.39612):     11.2710237040    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    2  ( B1 )


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12


 QRLRVE -- linear response calculation for symmetry  2  ( B1 )
 QRLRVE -- operator label : XDIPLEN 
 QRLRVE -- operator spin  :   0
 QRLRVE -- frequencies :  0.077257  0.156891  0.255638  0.260466  0.432130
                          0.472415  0.025157  0.177758  0.518615



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 2  ( B1 ); triplet =   F


 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.15689 a.u.
 after    9 linear transformations is      5.40399152
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.74 * 10 **  4.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.25564 a.u.
 after    9 linear transformations is      8.03764290
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -1.07 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.26047 a.u.
 after    9 linear transformations is      8.10167385
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -1.10 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after    9 linear transformations is     15.11529672
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -7.71 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after    9 linear transformations is     25.20404012
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -3.80 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.17776 a.u.
 after    9 linear transformations is      7.17371944
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.80 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.51861 a.u.
 after    9 linear transformations is    -97.17484122
 INERTIA (POS,ZER,NEG) of reduced matrix is   16    0    2
 Determinant of reduced matrix is  6.67 * 10 **  4.0
RSPORT:    6 out of    9 new trial vectors linear dependent

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.15689 a.u.
 after   12 linear transformations is      5.49597203
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -2.97 * 10 ** 13.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.25564 a.u.
 after   12 linear transformations is      8.03780232
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.10 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.26047 a.u.
 after   12 linear transformations is      8.10183400
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.14 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after   12 linear transformations is     15.11744425
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -7.88 * 10 ** 14.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after   12 linear transformations is     25.21493470
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -3.87 * 10 ** 14.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.17776 a.u.
 after   12 linear transformations is      7.17470798
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -2.90 * 10 ** 14.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  2,  triplet       F

 GP * SOLUTION vector at frequency     0.51861 a.u.
 after   12 linear transformations is    -96.78171545
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  6.81 * 10 ** 13.0

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   24)
 RSP solution vector no.    1; norm of residual   6.38e-16
 RSP solution vector no.    2; norm of residual   5.83e-16
 RSP solution vector no.    3; norm of residual   6.42e-16
 RSP solution vector no.    4; norm of residual   4.42e-16
 RSP solution vector no.    5; norm of residual   5.03e-16
 RSP solution vector no.    6; norm of residual   6.14e-16
 RSP solution vector no.    7; norm of residual   6.44e-16
 RSP solution vector no.    8; norm of residual   8.49e-16
 RSP solution vector no.    9; norm of residual   3.10e-16

 *** RSPCTL MICROITERATIONS CONVERGED

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.772568e-01

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.07726):     6.96542750211    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.156891e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.15689):     5.49597203250    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.255638e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.25564):     8.03780231808    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.260466e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.26047):     8.10183400472    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.432130e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.43213):     15.1174442509    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.472415e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.47241):     25.2149346951    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.251572e-01

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.02516):     6.86569374338    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.177758e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.17776):     7.17470798167    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    2  ( B1 )    LABEL   XDIPLEN     FREQUENCY   0.518615e+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.51861):    -96.7817154518    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    3  ( B2 )


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7


 QRLRVE -- linear response calculation for symmetry  3  ( B2 )
 QRLRVE -- operator label : YDIPLEN 
 QRLRVE -- operator spin  :   0
 QRLRVE -- frequencies :  0.260466  0.396124  0.472415  0.025157  0.206986
                          0.177758  0.156891  0.518615  0.331146



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 3  ( B2 ); triplet =   F

RSPORT:    2 out of    9 new trial vectors linear dependent

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  3,  triplet       F

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after    7 linear transformations is      4.42865979
 INERTIA (POS,ZER,NEG) of reduced matrix is   13    0    1
 Determinant of reduced matrix is -4.02 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix. Symmetry  3,  triplet       F

 GP * SOLUTION vector at frequency     0.51861 a.u.
 after    7 linear transformations is      6.12077760
 INERTIA (POS,ZER,NEG) of reduced matrix is   13    0    1
 Determinant of reduced matrix is -4.83 * 10 **  6.0

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   14)
 RSP solution vector no.    1; norm of residual   2.03e-16
 RSP solution vector no.    2; norm of residual   3.97e-16
 RSP solution vector no.    3; norm of residual   3.76e-16
 RSP solution vector no.    4; norm of residual   6.36e-16
 RSP solution vector no.    5; norm of residual   5.81e-16
 RSP solution vector no.    6; norm of residual   4.09e-16
 RSP solution vector no.    7; norm of residual   2.50e-16
 RSP solution vector no.    8; norm of residual   2.06e-16
 RSP solution vector no.    9; norm of residual   2.86e-16

 *** RSPCTL MICROITERATIONS CONVERGED

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.260466e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.26047):     3.40521832472    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.396124e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.39612):     6.57835558540    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.472415e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.47241):     4.42865979371    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.251572e-01

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.02516):     2.80740989356    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.206986e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.20699):     3.14758222929    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.177758e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.17776):     3.04733274578    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.156891e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.15689):     2.98896571257    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.518615e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.51861):     6.12077760262    

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    3  ( B2 )    LABEL   YDIPLEN     FREQUENCY   0.331146e+00

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.33115):     4.02686587677    


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156891    0.313782    4.770825


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260466    0.520931   -1.489717


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    3.196906


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156891    0.313782   20.939698


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260466    0.520931   -0.389093


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    2.952066


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156891    0.313782   -0.232312


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260466    0.520931   -0.131524


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    1.879617


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.077257    0.154514    1.746157


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.255638    0.511276   -9.362958


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.432130    0.864261   -0.769929


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.077257    0.154514    1.746157


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.255638    0.511276   -9.362958


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.432130    0.864261   -0.769929


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.206986    0.413971   -1.245238


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.331146    0.662292   -5.788493


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.396124    0.792249    2.172688


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.206986    0.413971   -1.245238


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.331146    0.662292   -5.788493


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.396124    0.792249    2.172688


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.025157    0.050314    0.047552


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.177758    0.355516    1.307598


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.518615    1.037230    1.149592


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.025157    0.050314    0.047552


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.177758    0.355516    1.307598


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.518615    1.037230    1.149592


       *******************************************************************
       ************ FINAL RESULTS FROM TWO-PHOTON CALCULATION ************
       *******************************************************************

 The two-photon absorption strength for an average molecular  
 orientation is computed according to formulas given by       
 P.R. Monson and W.M. McClain in J. Chem. Phys. 53:29, 1970   
 and W.M. McClain in J. Chem. Phys. 55:2789, 1971.            
 The absorption depends on the light polarization.            
 A monochromatic light source is assumed.                     

 All results are presented in atomic units, except the        
 excitation energy which is given in eV and two-photon cross  
 section which is given in GM. A FWHM of 0.1 eV is assumed.   

   Conversion factors:
      1 a.u. = 1.896788 10^{-50} cm^4 s/photon
      1 GM = 10^{-50} cm^4 s/photon


                  +--------------------------------+
                  | Two-photon transition tensor S |
                  +--------------------------------+
     ---------------------------------------------------------------
     Sym  No  Energy     Sxx     Syy     Szz     Sxy     Sxz     Syz
     ---------------------------------------------------------------
       1   1    8.54    20.9    -0.2     4.8     0.0     0.0     0.0
       1   2   14.18    -0.4    -0.1    -1.5     0.0     0.0     0.0
       1   3   25.71     3.0     1.9     3.2     0.0     0.0     0.0
       2   1    4.20     0.0     0.0     0.0     0.0     1.7     0.0
       2   2   13.91     0.0     0.0     0.0     0.0    -9.4     0.0
       2   3   23.52     0.0     0.0     0.0     0.0    -0.8     0.0
       3   1   11.26     0.0     0.0     0.0     0.0     0.0    -1.2
       3   2   18.02     0.0     0.0     0.0     0.0     0.0    -5.8
       3   3   21.56     0.0     0.0     0.0     0.0     0.0     2.2
       4   1    1.37     0.0     0.0     0.0     0.0     0.0     0.0
       4   2    9.67     0.0     0.0     0.0     1.3     0.0     0.0
       4   3   28.22     0.0     0.0     0.0     1.1     0.0     0.0
     ---------------------------------------------------------------


                    Transition probabilities (a.u.)         
                   -----------------------------------------
                    D  =  2*Df + 4*Dg, Linear   polarization
                    D  = -2*Df + 6*Dg, Circular polarization
                    Df = sum(i,j){ S_ii * S_jj }/30         
                    Dg = sum(i,j){ S_ij * S_ij }/30         

                   Two-photon cross sections                
         ---------------------------------------------------
          sigma  =  8*pi^3*alpha^2*hbar/e^4 * E^2*D   (a.u.)

                             Polarization ratio      
                      -------------------------------
                          R  = (-Df+3*Dg)/(Df+2*Dg)  


                   +-----------------------------------+
                   |   Two-photon absorption summary   |
                   +-----------------------------------+
   ---------------------------------------------------------------------------------
   Sym  No  Energy  Polarization         Df         Dg          D      sigma       R
   ---------------------------------------------------------------------------------
     1   1    8.54   Linear       0.216e+02  0.154e+02  0.105e+03  0.560e+01    0.47
     1   1    8.54   Circular     0.216e+02  0.154e+02  0.490e+02  0.262e+01    0.47
     1   2   14.18   Linear       0.135e+00  0.796e-01  0.588e+00  0.865e-01    0.35
     1   2   14.18   Circular     0.135e+00  0.796e-01  0.208e+00  0.306e-01    0.35
     1   3   25.71   Linear       0.215e+01  0.749e+00  0.729e+01  0.353e+01    0.03
     1   3   25.71   Circular     0.215e+01  0.749e+00  0.196e+00  0.951e-01    0.03
     2   1    4.20   Linear       0.000e+00  0.203e+00  0.813e+00  0.105e-01    1.50
     2   1    4.20   Circular     0.000e+00  0.203e+00  0.122e+01  0.158e-01    1.50
     2   2   13.91   Linear       0.000e+00  0.584e+01  0.234e+02  0.332e+01    1.50
     2   2   13.91   Circular     0.000e+00  0.584e+01  0.351e+02  0.497e+01    1.50
     2   3   23.52   Linear       0.000e+00  0.395e-01  0.158e+00  0.641e-01    1.50
     2   3   23.52   Circular     0.000e+00  0.395e-01  0.237e+00  0.961e-01    1.50
     3   1   11.26   Linear       0.000e+00  0.103e+00  0.413e+00  0.384e-01    1.50
     3   1   11.26   Circular     0.000e+00  0.103e+00  0.620e+00  0.577e-01    1.50
     3   2   18.02   Linear       0.000e+00  0.223e+01  0.894e+01  0.213e+01    1.50
     3   2   18.02   Circular     0.000e+00  0.223e+01  0.134e+02  0.319e+01    1.50
     3   3   21.56   Linear       0.000e+00  0.315e+00  0.126e+01  0.429e+00    1.50
     3   3   21.56   Circular     0.000e+00  0.315e+00  0.189e+01  0.643e+00    1.50
     4   1    1.37   Linear       0.000e+00  0.151e-03  0.603e-03  0.828e-06    1.50
     4   1    1.37   Circular     0.000e+00  0.151e-03  0.904e-03  0.124e-05    1.50
     4   2    9.67   Linear       0.000e+00  0.114e+00  0.456e+00  0.313e-01    1.50
     4   2    9.67   Circular     0.000e+00  0.114e+00  0.684e+00  0.469e-01    1.50
     4   3   28.22   Linear       0.000e+00  0.881e-01  0.352e+00  0.206e+00    1.50
     4   3   28.22   Circular     0.000e+00  0.881e-01  0.529e+00  0.309e+00    1.50
   ---------------------------------------------------------------------------------

 >>>> Total CPU  time used in RESPONSE:   8.80 seconds
 >>>> Total wall time used in RESPONSE:   9.66 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   9.91 seconds
 >>>> Total wall time used in DALTON:  11.09 seconds

 
     Date and time (Linux)  : Wed Apr  9 11:47:22 2014
     Host name              : rfb-ntf-m0022.bargi.uit.no              
END REFOUT
