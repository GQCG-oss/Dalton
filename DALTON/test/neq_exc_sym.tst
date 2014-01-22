########## Test description ########################
START DESCRIPTION
KEYWORDS: dft pcm linear response residue formaldehyde symmetry nonequilibrium short
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
tramom
enedft
tes
nuc
sym
cmass
neqrsp
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN WAVE FUNCTION
.RUN RESPONSE
*PCM
.SOLVNT
WATER
.NPCMMT
0
.NESFP
4
.ICESPH
2
.NEQRSP
*PCMCAV
.INA
1
2
3
4
.RIN 
1.7
1.5
1.2
1.2
.AREATS
0.30
**INTEGRALS
.DIPLEN
**WAVE FUNCTIONS
.DFT
LDA
**RESPONSE
*LINEAR
.SINGLE RESIDUE
.DIPLEN 
.ROOTS
2 2 2 2
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
STO-3G


    3    2 X  Y    1
        6.    1
C     0.000000     0.000000     0.000000
        8.    1
O     0.000000     0.000000     1.220000
        1.    1
H     0.943102     0.000000    -0.544500
END MOLINP

########## Reference Output ########################
START REFOUT


     ************************************************************************
     *************** Dalton - An Electronic Structure Program ***************
     ************************************************************************

    This is output from DALTON (Release Dalton2013 patch 0)
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
       Release DALTON2013.0 (2013), see http://daltonprogram.org
   ----------------------------------------------------------------------------

    Authors in alphabetical order (major contribution(s) in parenthesis):

  Kestutis Aidas,           Vilnius University,           Lithuania   (QM/MM)
  Celestino Angeli,         University of Ferrara,        Italy       (NEVPT2)
  Keld L. Bak,              UNI-C,                        Denmark     (AOSOPPA, non-adiabatic coupling, magnetic properties)
  Vebjoern Bakken,          University of Oslo,           Norway      (DALTON; geometry optimizer, symmetry detection)
  Radovan Bast,             KTH Stockholm                 Sweden      (DALTON installation and execution frameworks)
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
  Filip Pawlowski,          Kazimierz Wielki University   Poland      (CC3)
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

     Date and time (Linux)  : Fri Sep 20 09:38:31 2013
     Host name              : compute-1-3.local                       

 * Work memory size             :    64000000 =  488.28 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/dalton/build-RB-mpi1.4.4-procs4-gnu4.6.2-release/perl-pid.11712__2013_9_20__9.38
   2) /home/arnfinn/dalton/build-RB-mpi1.4.4-procs4-gnu4.6.2-release/basis


Compilation information
-----------------------

 Who compiled             | arnfinn
 Host                     | compute-1-3.local
 System                   | Linux-2.6.18-238.el5
 CMake generator          | Unix Makefiles
 Processor                | x86_64
 64-bit integers          | OFF
 MPI                      | ON
 Fortran compiler         | /global/apps/openmpi/1.4.4-gnu/bin/mpif90
 Fortran compiler version | GNU Fortran (GCC) 4.6.2
 C compiler               | /global/apps/openmpi/1.4.4-gnu/bin/mpicc
 C compiler version       | gcc (GCC) 4.6.2
 C++ compiler             | /global/apps/openmpi/1.4.4-gnu/bin/mpicxx
 C++ compiler version     | unknown
 BLAS                     | -Wl,--start-group;/global/apps/intel/mkl/lib/intel
                          | 64/libmkl_core.so;/global/apps/intel/mkl/lib/intel
                          | 64/libmkl_gf_lp64.so;/global/apps/intel/mkl/lib/in
                          | tel64/libmkl_sequential.so;/usr/lib64/libpthread.s
                          | o;/usr/lib64/libm.so;-fopenmp;-Wl,--end-group
 LAPACK                   | -Wl,--start-group;/global/apps/intel/mkl/lib/intel
                          | 64/libmkl_lapack95_lp64.a;/global/apps/intel/mkl/l
                          | ib/intel64/libmkl_gf_lp64.so;-fopenmp;-Wl,--end-gr
                          | oup
 Static linking           | OFF
 Last Git revision        | fcecae41c6e1063a09a499dd9e7f51b9ee868d86
 Configuration time       | 2013-09-20 09:06:39.448412


   Content of the .dal input file
 ----------------------------------

**DALTON INPUT                                    
.RUN WAVE FUNCTION                                
.RUN RESPONSE                                     
*PCM                                              
.SOLVNT                                           
WATER                                             
.NPCMMT                                           
0                                                 
.NESFP                                            
4                                                 
.ICESPH                                           
2                                                 
.NEQRSP                                           
*PCMCAV                                           
.INA                                              
1                                                 
2                                                 
3                                                 
4                                                 
.RIN                                              
1.7                                               
1.5                                               
1.2                                               
1.2                                               
.AREATS                                           
0.30                                              
**INTEGRALS                                       
.DIPLEN                                           
**WAVE FUNCTIONS                                  
.DFT                                              
LDA                                               
**RESPONSE                                        
*LINEAR                                           
.SINGLE RESIDUE                                   
.DIPLEN                                           
.ROOTS                                            
2 2 2 2                                           
**END OF DALTON INPUT                             


   Content of the .mol file
 ----------------------------

BASIS                                                                          
STO-3G                                                                         
                                                                               
                                                                               
    3    2 X  Y    1                                                           
        6.    1                                                                
C     0.000000     0.000000     0.000000                                       
        8.    1                                                                
O     0.000000     0.000000     1.220000                                       
        1.    1                                                                
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

     ICESPH =       2     NESFP =       4
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =T
     POLYG          60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 


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
  Basis set file used for this atomic type with Z =   6 :
     "/home/arnfinn/dalton/build-RB-mpi1.4.4-procs4-gnu4.6.2-release/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  Basis set file used for this atomic type with Z =   8 :
     "/home/arnfinn/dalton/build-RB-mpi1.4.4-procs4-gnu4.6.2-release/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  Basis set file used for this atomic type with Z =   1 :
     "/home/arnfinn/dalton/build-RB-mpi1.4.4-procs4-gnu4.6.2-release/basis/STO-3G"


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
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000e+00    0.0000000000e+00    0.0000000000e+00    1.7000000000e+00
   2    0.0000000000e+00    0.0000000000e+00    2.3054658725e+00    1.5000000000e+00
   3    1.7822044879e+00    0.0000000000e+00   -1.0289558751e+00    1.2000000000e+00
   4   -1.7822044879e+00    0.0000000000e+00   -1.0289558751e+00    1.2000000000e+00


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

  Basis set used is "STO-3G" from the basis set library.

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


 Default print level:        1

 * Nuclear model: Point charge

 Calculation of one- and two-electron Hamiltonian integrals.

 The following one-electron property integrals are calculated as requested:
          - overlap integrals
          - dipole length integrals

 Center of mass  (bohr):      0.000000000000      0.000000000000      1.159648802225
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************


 Threshold for neglecting two-electron integrals:  1.00e-12
 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014


 MEMORY USED TO GENERATE CAVITY =    432202


 Total number of spheres =    4
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1    0.000000000    0.000000000    0.000000000    2.040000000   25.046551891
   2    0.000000000    0.000000000    1.220000000    1.800000000   22.984715878
   3    0.943102000    0.000000000   -0.544500000    1.440000000    9.281425262
   4   -0.943102000    0.000000000   -0.544500000    1.440000000    9.281425262

 Total number of tesserae =     304
 Surface area =   66.59411829 (A^2)    Cavity volume =   48.28620692 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....

  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   0.10 seconds
 >>>> Total wall time used in HERMIT:   0.35 seconds


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

 
     Date and time (Linux)  : Fri Sep 20 09:38:32 2013
     Host name              : compute-1-3.local                       

 Title lines from ".mol" input file:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

@    Restricted, closed shell Kohn-Sham DFT calculation.

@    Time-dependent Kohn-Sham DFT calculation (TD-DFT).

 Initial molecular orbitals are obtained according to
 ".MOSTART EHT   " input option

@    QM part is embedded in an environment :

@         Solvation model: PCM

     Wave function specification
     ============================
@    For the wave function of type :      >>> KS-DFT <<<
@    Number of closed shell electrons          16
@    Number of electrons in active shells       0
@    Total charge of the molecule               0

@    Spin multiplicity and 2 M_S                1         0
     Total number of symmetries                 4
@    Reference state symmetry                   1
 
     This is a DFT calculation of type: LDA

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1    2    3    4
                                       --- |  ---  ---  ---  ---
     Total number of orbitals           12 |    7    3    2    0
     Number of basis functions          12 |    7    3    2    0

      ** Automatic occupation of RKS orbitals **

      -- Initial occupation of symmetries is determined from extended Huckel guess.           
      -- Initial occupation of symmetries is :
@    Occupied SCF orbitals               8 |    5    2    1    0

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05
 
     This is a DFT calculation of type: LDA

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.97018
 NUCLEAR APPARENT CHARGE -15.78962
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =    5

 Automatic occupation of symmetries with  16 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    SCF occupation
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.964563272998  16.0000008466    8.47e-07
@  1  -111.730265028     -2.014834295276e-02   2.59e+00  -1.12e+02     5   2   1   0
      Virial theorem: -V/T =      1.996185
@      MULPOP C       1.00; O      -0.74; H   _1 -0.13; H   _2 -0.13; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.683589314663  16.0000010604    1.06e-06
@  2  -111.069683994     -4.739886745406e-02   3.47e+00   6.61e-01     5   2   1   0
      Virial theorem: -V/T =      2.003390
@      MULPOP C      -1.90; O       1.43; H   _1  0.24; H   _2  0.24; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.972447269642  16.0000009048    9.05e-07
@  3  -111.918604622     -1.827882238131e-02   1.17e+00  -8.49e-01     5   2   1   0
      Virial theorem: -V/T =      1.999485
@      MULPOP C       0.53; O      -0.74; H   _1  0.10; H   _2  0.10; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.813173332331  16.0000009343    9.34e-07
@  4  -112.031958945     -1.071362942117e-03   1.49e-01  -1.13e-01     5   2   1   0
      Virial theorem: -V/T =      2.007778
@      MULPOP C      -0.07; O      -0.09; H   _1  0.08; H   _2  0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.837346289814  16.0000009364    9.36e-07
@  5  -112.033866467     -2.150639307168e-03   2.04e-02  -1.91e-03     5   2   1   0
      Virial theorem: -V/T =      2.006504
@      MULPOP C      -0.05; O      -0.16; H   _1  0.11; H   _2  0.11; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.834172658732  16.0000009352    9.35e-07
@  6  -112.033907917     -2.019985333266e-03   2.67e-03  -4.15e-05     5   2   1   0
      Virial theorem: -V/T =      2.006688
@      MULPOP C      -0.04; O      -0.16; H   _1  0.10; H   _2  0.10; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.834464847252  16.0000009352    9.35e-07
@  7  -112.033908472     -2.037092760638e-03   4.41e-04  -5.55e-07     5   2   1   0
      Virial theorem: -V/T =      2.006674
@      MULPOP C      -0.04; O      -0.16; H   _1  0.10; H   _2  0.10; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.834502482907  16.0000009351    9.35e-07
@  8  -112.033908486     -2.039563541146e-03   2.16e-05  -1.45e-08     5   2   1   0
      Virial theorem: -V/T =      2.006672
@      MULPOP C      -0.04; O      -0.16; H   _1  0.10; H   _2  0.10; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -13.834504479260  16.0000009351    9.35e-07
@  9  -112.033908486     -2.039693823985e-03   1.55e-07  -3.26e-11     5   2   1   0

@ *** DIIS converged in   9 iterations !
@     Converged SCF energy, gradient:   -112.033908486406    1.55e-07
    - total time used in SIRFCK :              0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    5    2    1    0

 Sym       Kohn-Sham orbital energies

  1    -18.36962387    -9.66197394    -0.90481808    -0.50692421    -0.27783866
         0.36676359     0.53279079

  2     -0.38754503    -0.09325517     0.46802251

  3     -0.27474215     0.02488548

    E(LUMO) :     0.02488548 au (symmetry 3)
  - E(HOMO) :    -0.09325517 au (symmetry 2)
  ------------------------------------------
    gap     :     0.11814065 au

 >>> Writing SIRIFC interface file <<<

 >>>> CPU and wall time for SCF :       2.534       2.599


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


@    Spin multiplicity:           1
@    Spatial symmetry:            1
@    Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

@    Final DFT energy:           -112.033908486406                 
@    Nuclear repulsion:            31.163673729192
@    Electronic energy:          -143.195542521774

@    Final gradient norm:           0.000000155480

 
     Date and time (Linux)  : Fri Sep 20 09:38:34 2013
     Host name              : compute-1-3.local                       

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

    Orbital         1        2        3        4        5        6        7
   1 C   :1s     0.0008   0.9893  -0.1333  -0.1871   0.0449   0.2068   0.1104
   2 C   :1s    -0.0100   0.0485   0.2966   0.5784  -0.1468  -1.2732  -0.7076
   3 C   :2pz   -0.0089   0.0013   0.1775  -0.2835  -0.4156   0.4808  -1.1358
   4 O   :1s     0.9927   0.0007  -0.2127   0.0974  -0.1230  -0.0195  -0.1175
   5 O   :1s     0.0333  -0.0100   0.7200  -0.4013   0.6003   0.1079   0.8568
   6 O   :2pz   -0.0070   0.0043  -0.2264  -0.0660   0.6758  -0.1872  -0.9366
   7 H   :1s     0.0004  -0.0110   0.0348   0.2666   0.1248   0.9083  -0.0869

 Molecular orbitals for symmetry species  2
 ------------------------------------------

    Orbital         1        2        3
   1 C   :2px    0.5918  -0.0483   1.1609
   2 O   :2px    0.2788   0.9222  -0.3534
   3 H   :1s     0.3235  -0.3314  -0.8594

 Molecular orbitals for symmetry species  3
 ------------------------------------------

    Orbital         1        2
   1 C   :2py    0.6087   0.8212
   2 O   :2py    0.6772  -0.7657



 >>>> Total CPU  time used in SIRIUS :      2.54 seconds
 >>>> Total wall time used in SIRIUS :      2.64 seconds

 
     Date and time (Linux)  : Fri Sep 20 09:38:34 2013
     Host name              : compute-1-3.local                       


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




  Linear Response single residue calculation
 -------------------------------------------

 Non-equilibrium PCM solvent model requested    : INERSI =T

 Static dielectric constant                     : EPSTAT = 78.3900
 Optical dielectric constant                    : EPSOL  =  1.7760
 Print level                                    : IPRPP  =   2
 Maximum number of iterations for eigenval.eqs. : MAXITP =  60
 Threshold for convergence of eigenvalue eqs.   : THCPP  = 1.000e-03
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

      2 Excitation energies are calculated for symmetry no.    1

      1 property residues are calculated with labels:

               ZDIPLEN 

      2 Excitation energies are calculated for symmetry no.    2

      1 property residues are calculated with labels:

               XDIPLEN 

      2 Excitation energies are calculated for symmetry no.    3

      1 property residues are calculated with labels:

               YDIPLEN 

      2 Excitation energies are calculated for symmetry no.    4

 Integral transformation: Total CPU and WALL times (sec)       0.002       0.319


   SCF energy         :     -112.033908486406247
 -- inactive part     :     -143.195542521774428
 -- nuclear repulsion :       31.163673729192165


                    *****************************************
                    *** DFT response calculation (TD-DFT) ***
                    *****************************************



 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            2
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      13
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      13
 Electrons in DFTMOMO:   16.00000093514756



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   18)
 RSP solution vector no.    1; norm of residual   4.67e-05
 RSP solution vector no.    2; norm of residual   2.99e-05

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    ZDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  0.60223593     *ENERGY(eV):   11.976358    
@ STATE NO:    2 *TRANSITION MOMENT:  0.64913515     *ENERGY(eV):   16.960758    


  ******************************************************************************
  *** @ Excit. operator sym 1 & ref. state sym 1 => excited state symmetry 1 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.44012308    au
@                      11.976358     eV
@                      96595.849     cm-1
                       1155.5430     kJ / mol

@ Total energy :      -111.59379     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  0.10641827      (Transition moment :  0.60223593     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      5    7  0.119989 -0.035544  0.704583  0.425299  0.024191  0.017044
      9   10  0.220345 -0.026659  0.579400  0.224201  0.061011  0.035350
     11   12 -0.668883 -0.127036  0.949794  0.947563  0.293598  0.278858

@ Overlap diagnostic LAMBDA =    0.8706


 @ Excited state no:    2 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.62329642    au
@                      16.960758     eV
@                      136797.75     cm-1
                       1636.4645     kJ / mol

@ Total energy :      -111.41061     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  0.17509495      (Transition moment :  0.64913515     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      5    6  0.366972 -0.006167  0.410305  0.093488  0.139232  0.057128
      5    7  0.130933 -0.017106  0.704583  0.425299  0.021915  0.015441
      9   10 -0.568911 -0.029786  0.579400  0.224201  0.290657  0.168406
     11   12 -0.162727 -0.057068  0.949794  0.947563  0.011164  0.010603

@ Overlap diagnostic LAMBDA =    0.5452


 Time used in polarization propagator calculation is      2.52 CPU seconds for symmetry 1


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    2

 Number of excitations of this symmetry            2
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       9
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       9
 Electrons in DFTMOMO:   16.00000093514756



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  2; triplet =   F


 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   1.75e-04
 RSP solution vector no.    2; norm of residual   2.66e-04

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    XDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -0.58812794     *ENERGY(eV):   13.279050    
@ STATE NO:    2 *TRANSITION MOMENT:  0.11341229     *ENERGY(eV):   17.902099    


  ******************************************************************************
  *** @ Excit. operator sym 2 & ref. state sym 1 => excited state symmetry 2 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.48799611    au
@                      13.279050     eV
@                      107102.76     cm-1
                       1281.2336     kJ / mol

@ Total energy :      -111.54591     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  0.11253010      (Transition moment : -0.58812794     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      8    6  0.087584 -0.018774  0.725441  0.503705  0.011312  0.008206
      9    6  0.700606  0.020676  0.449403  0.118986  0.462305  0.207762

@ Overlap diagnostic LAMBDA =    0.4571


 @ Excited state no:    2 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.65789007    au
@                      17.902099     eV
@                      144390.18     cm-1
                       1727.2901     kJ / mol

@ Total energy :      -111.37602     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  5.64134088e-03  (Transition moment :  0.11341229     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      9    7  0.702695  0.017428  0.501885  0.301710  0.469591  0.235681

@ Overlap diagnostic LAMBDA =    0.5021


 Time used in polarization propagator calculation is      1.66 CPU seconds for symmetry 2


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    3

 Number of excitations of this symmetry            2
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7
 Electrons in DFTMOMO:   16.00000093514756



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  3; triplet =   F


 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   1.85e-05
 RSP solution vector no.    2; norm of residual   1.83e-06

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    YDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  0.18087511     *ENERGY(eV):   9.3145390    
@ STATE NO:    2 *TRANSITION MOMENT: -0.61568756     *ENERGY(eV):   16.017152    


  ******************************************************************************
  *** @ Excit. operator sym 3 & ref. state sym 1 => excited state symmetry 3 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.34230302    au
@                      9.3145390     eV
@                      75126.828     cm-1
                       898.71644     kJ / mol

@ Total energy :      -111.69161     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  7.46581300e-03  (Transition moment :  0.18087511     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      5   12 -0.704729 -0.043242  0.592639  0.318826  0.437565  0.259318

@ Overlap diagnostic LAMBDA =    0.5943


 @ Excited state no:    2 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.58861951    au
@                      16.017152     eV
@                      129187.05     cm-1
                       1545.4203     kJ / mol

@ Total energy :      -111.44529     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  0.14875246      (Transition moment : -0.61568756     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      4   12 -0.640179 -0.032459  0.669264  0.352450  0.369324  0.247175
     11    6 -0.295411  0.012700  0.411125  0.167306  0.094933  0.039029

@ Overlap diagnostic LAMBDA =    0.6166


 Time used in polarization propagator calculation is      1.63 CPU seconds for symmetry 3


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    4

 Number of excitations of this symmetry            2
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       4
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       3
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       3
 Electrons in DFTMOMO:   16.00000093514756



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  4; triplet =   F

RSPORT:    1 out of    2 new trial vectors linear dependent

 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   1.02e-07
 RSP solution vector no.    2; norm of residual   9.00e-08

 *** RSPCTL MICROITERATIONS CONVERGED


  ******************************************************************************
  *** @ Excit. operator sym 4 & ref. state sym 1 => excited state symmetry 4 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  4
 ---------------------------------------

@ Excitation energy :  0.13701193    au
@                      3.7282844     eV
@                      30070.644     cm-1
                       359.72478     kJ / mol

@ Total energy :      -111.89690     au


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      9   12 -0.709017 -0.052439  0.480781  0.262809  0.431094  0.207262

@ Overlap diagnostic LAMBDA =    0.4808


 @ Excited state no:    2 in symmetry  4
 ---------------------------------------

@ Excitation energy :  0.42490166    au
@                      11.562162     eV
@                      93255.134     cm-1
                       1115.5791     kJ / mol

@ Total energy :      -111.60901     au


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      8   12 -0.706916 -0.010529  0.517868  0.214109  0.484955  0.251142

@ Overlap diagnostic LAMBDA =    0.5178


 Time used in polarization propagator calculation is      0.99 CPU seconds for symmetry 4

 >>>> Total CPU  time used in RESPONSE:   6.81 seconds
 >>>> Total wall time used in RESPONSE:   7.26 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   9.46 seconds
 >>>> Total wall time used in DALTON:  10.49 seconds

 
     Date and time (Linux)  : Fri Sep 20 09:38:42 2013
     Host name              : compute-1-3.local                       
END REFOUT

