########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm qm3 qm3pcm dft residue linear response essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enedft
qm3pcmconv
qm3ele
qm3pol
qm3wdv
qm3pcmene
pcmsol
surf
nucchg
tramom
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON
.RUN RESPONSE
*PCM
.NEQRSP
.SOLVNT
WATER
.NPCMMT
0
.NESFP
8
.ICESPH
2
*PCMCAV
.RIN
1.7
1.5
1.2
1.2
1.5
1.5
1.2
1.2
.INA
1
2
3
4
5
6
7
8
.AREATS
 0.3
.ALPHA
 1.2
 1.2
 1.2
 1.2
 1.2
 1.2
 1.2
 1.2
*QM3
.PRINT
 1
.QM3
.MMPCM
.THRSMP
 1.0D-9
**INTEGRALS
.DIPLEN
.NUCPOT
.NELFLD
.PRINT
 2
**WAVE FUNCTIONS
.DFT
B3LYP
*SCF INPUT
.THRESH
1.0D-12
**RESPONSE
*LINEAR
.DIPLEN
.SINGLE RES
.ROOTS
4
.THCPP
 1.0D-04
.PRINT
5
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS
CH2O PLUS 2 WATERS
------------------------
AtomTypes=5 NoSymmetry Angstrom
        Charge=6.0   Atoms=1    Basis=STO-3G
C           -1.588367    -.770650     .029109 0 1
        Charge=8.0   Atoms=1    Basis=STO-3G
O           -1.657083     .436069    -.009750 0 2
        Charge=1.0   Atoms=2    Basis=STO-3G
H           -.620668   -1.294822      .054251 0 3
H           -2.508043   -1.382001     .040282 0 4
   Charge=-0.669     Atoms=2   Basis=MM
O             .975536    1.507219    -.082590 1 1
O            1.633705   -1.233280    -.078786 2 1
    Charge=0.3345    Atoms=4   Basis=MM
H             .023200    1.300096    -.088063 1 2
H            1.113831    2.040431     .704778 1 3
H            1.600739    -.258639    -.086953 2 2
H            2.303878   -1.461228     .570065 2 3
END MOLINP

########## POTENTIAL.INP ############################
START POTINP
**SYSTP
.NUMMTP
 1
.TYPE
 0
.MODEL
 SPC
.CHARGS
 4
 0.0000
 0.0000
 0.0000
 0.0000
*******
.TYPE
 1-2
.MODEL
 SPC_E01
.ALPISO
 1
 9.501
*******
**TWOIA (i,j=0,1,2,...,N; if i=0 then j.neq.0)
.LJ_A
 2
 0.00000
 0.00000
.LJ_B
 2
 0.000
 0.000
**END OF
END POTINP

########## Reference Output ########################
START REFOUT


   ****************************************************************************
   *************** DALTON2011 - An electronic structure program ***************
   ****************************************************************************

    This is output from DALTON Release 2011 (DEVELOPMENT VERSION)
   ----------------------------------------------------------------------------
    NOTE:
     
    DALTON is an experimental code for the evaluation of molecular
    properties using (MC)SCF, DFT, CI, and CC wave functions.
    The authors accept no responsibility for the performance of
    the code or for the correctness of the results.
     
    The code (in whole or part) is provided under a licence and
    is not to be reproduced for further distribution without
    the written permission of the authors or their representatives.
     
    See the home page "http://daltonprogram.org" for further information.
     
    If results obtained with this code are published,
    an appropriate citation would be:
     
    "Dalton, a molecular electronic structure program,
    Release Dalton2011 (2011), see http://daltonprogram.org"
   ----------------------------------------------------------------------------

    Authors in alphabetical order (major contribution(s) in parenthesis):

  Celestino Angeli,         University of Ferrara,        Italy       (NEVPT2)
  Keld L. Bak,              UNI-C,                        Denmark     (AOSOPPA, non-adiabatic coupling, magnetic properties)
  Vebjoern Bakken,          University of Oslo,           Norway      (DALTON; geometry optimizer, symmetry detection)
  Ove Christiansen,         Aarhus University,            Denmark     (CC module)
  Renzo Cimiraglia,         University of Ferrara,        Italy       (NEVPT2)
  Sonia Coriani,            University of Trieste,        Italy       (CC module, MCD in RESPONS)
  Paal Dahle,               University of Oslo,           Norway      (Parallelization)
  Erik K. Dalskov,          UNI-C,                        Denmark     (SOPPA)
  Thomas Enevoldsen,        SDU - Odense University,      Denmark     (SOPPA)
  Berta Fernandez,          U. of Santiago de Compostela, Spain       (doublet spin, ESR in RESPONS)
  Lara Ferrighi,            Aarhus University,            Denmark     (PCM Cubic response)
  Heike Fliegl,             University of Helsinki,       Finland     (CCSD(R12))
  Luca Frediani,            University of Tromsoe,        Norway      (PCM)
  Bin Gao,                  University of Tromsoe,        Norway      (Gen1Int module)
  Christof Haettig,         Ruhr University Bochum,       Germany     (CC module)
  Kasper Hald,              Aarhus University,            Denmark     (CC module)
  Asger Halkier,            Aarhus University,            Denmark     (CC module)
  Hanne Heiberg,            University of Oslo,           Norway      (geometry analysis, selected one-electron integrals)
  Trygve Helgaker,          University of Oslo,           Norway      (DALTON; ABACUS, ERI, DFT modules, London, and much more)
  Hinne Hettema,            University of Auckland,       New Zealand (quadratic response in RESPONS; SIRIUS supersymmetry)
  Brano Jansik              University of Aarhus          Denmark     (DFT cubic response)
  Hans Joergen Aa. Jensen,  Univ. of Southern Denmark,    Denmark     (DALTON; SIRIUS, RESPONS, ABACUS modules, London, and much more)
  Dan Jonsson,              University of Tromsoe,        Norway      (cubic response in RESPONS module)
  Poul Joergensen,          Aarhus University,            Denmark     (RESPONS, ABACUS, and CC modules)
  Sheela Kirpekar,          SDU - Odense University,      Denmark     (Mass-velocity & Darwin integrals)
  Wim Klopper,              University of Karlsruhe,      Germany     (R12 code in CC, SIRIUS, and ABACUS modules)
  Stefan Knecht,            Univ. of Southern Denmark,    Denmark     (Parallel CI)
  Rika Kobayashi,           ANU Supercomputer Facility,   Australia   (DIIS in CC, London in MCSCF)
  Jacob Kongsted,           Univ. of Southern Denmark,    Denmark     (QM/MM code)
  Henrik Koch,              University of Trondheim,      Norway      (CC module, Cholesky decomposition)
  Andrea Ligabue,           University of Modena,         Italy       (CTOCD, AOSOPPA)
  Ola B. Lutnaes,           University of Oslo,           Norway      (DFT Hessian)
  Juan I. Melo,             University of Buenos Aires,   Argentina   (LRESC, Relativistic Effects on NMR Shieldings)
  Kurt V. Mikkelsen,        University of Copenhagen,     Denmark     (MC-SCRF and QM/MM code)
  Christian Neiss,          Univ. Erlangen-Nuernberg,     Germany     (CCSD(R12))
  Christian B. Nielsen,     University of Copenhagen,     Denmark     (QM/MM code)
  Patrick Norman,           University of Linkoeping,     Sweden      (cubic response and complex response in RESPONS)
  Jeppe Olsen,              Aarhus University,            Denmark     (SIRIUS CI/density modules)
  Jogvan Magnus H. Olsen,   Univ. of Southern Denmark,    Denmark     (Polarizable embedding (PE) module)
  Anders Osted,             Copenhagen University,        Denmark     (QM/MM code)
  Martin J. Packer,         University of Sheffield,      UK          (SOPPA)
  Filip Pawlowski,          Kazimierz Wielki University   Poland      (CC3)
  Thomas B. Pedersen,       University of Oslo,           Norway      (Cholesky decomposition)
  Patricio F. Provasi,      University of Northeastern,   Argentina   (Analysis of coupling constants in localized orbitals)
  Zilvinas Rinkevicius,     KTH Stockholm,                Sweden      (open-shell DFT, ESR)
  Elias Rudberg,            KTH Stockholm,                Sweden      (DFT grid and basis info)
  Torgeir A. Ruden,         University of Oslo,           Norway      (Numerical derivatives in ABACUS)
  Kenneth Ruud,             University of Tromsoe,        Norway      (DALTON; ABACUS magnetic properties and  much more)
  Pawel Salek,              KTH Stockholm,                Sweden      (DALTON; DFT code)
  Claire C.M. Samson        University of Karlsruhe       Germany     (Boys localization, r12 integrals in ERI)
  Alfredo Sanchez de Meras, University of Valencia,       Spain       (CC module, Cholesky decomposition)
  Trond Saue,               CNRS/ULP Toulouse,            France      (direct Fock matrix construction)
  Stephan P. A. Sauer,      University of Copenhagen,     Denmark     (SOPPA(CCSD), SOPPA prop., AOSOPPA, vibrational g-factors)
  Bernd Schimmelpfennig,    Forschungszentrum Karlsruhe,  Germany     (AMFI module)
  Arnfinn H. Steindal,      University of Tromsoe,        Norway      (parallel QM/MM)
  K. O. Sylvester-Hvid,     University of Copenhagen,     Denmark     (MC-SCRF)
  Peter R. Taylor,          VLSCI/Univ. of Melbourne,     Australia   (Symmetry handling ABACUS, integral transformation)
  Andrew M. Teale,          University of Nottingham,     England     (DFT-AC, DFT-D)
  David P. Tew,             University of Bristol,        England     (CCSD(R12))
  Olav Vahtras,             KTH Stockholm,                Sweden      (triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)
  David J. Wilson,          La Trobe University,          Australia   (DFT Hessian and DFT magnetizabilities)
  Hans Agren,               KTH Stockholm,                Sweden      (SIRIUS module, MC-SCRF solvation model)
 --------------------------------------------------------------------------------

     Date and time (Linux)  : Wed Jun  5 20:39:19 2013
     Host name              : c18-12.local                            

 * Work memory size             :    64000000 =  488.28 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/dalton/DALTON/test/perl-pid.25283__2013_6_5__20.39
   2) /home/arnfinn/dalton/DALTON/build-trace/basis


Compilation information
-----------------------

 Who compiled             | arnfinn@stallo-1.local
 System                   | Linux-2.6.32-279.14.1.el6.x86_64
 CMake generator          | Unix Makefiles
 Processor                | x86_64
 64-bit integers          | OFF
 MPI                      | ON
 Fortran compiler         | /global/apps/openmpi/1.6.2/intel/13.0/bin/mpif90
 Fortran compiler version | ifort (IFORT) 13.0.1 20121010
 C compiler               | /global/apps/openmpi/1.6.2/intel/13.0/bin/mpicc
 C compiler version       | icc (ICC) 13.0.1 20121010
 C++ compiler             | /global/apps/openmpi/1.6.2/intel/13.0/bin/mpicxx
 C++ compiler version     | unknown
 BLAS                     | -Wl,--start-group;/global/apps/intel/2013.1/mkl/li
                          | b/intel64/libmkl_core.so;/global/apps/intel/2013.1
                          | /mkl/lib/intel64/libmkl_intel_lp64.so;/global/apps
                          | /intel/2013.1/mkl/lib/intel64/libmkl_sequential.so
                          | ;/usr/lib64/libpthread.so;/usr/lib64/libm.so;-open
                          | mp;-Wl,--end-group
 LAPACK                   | -Wl,--start-group;/global/apps/intel/2013.1/mkl/li
                          | b/intel64/libmkl_lapack95_lp64.a;/global/apps/inte
                          | l/2013.1/mkl/lib/intel64/libmkl_intel_lp64.so;-ope
                          | nmp;-Wl,--end-group
 Static linking           | OFF
 Last Git revision        | c6beed320385b29d2f94c4f311a5082c1d9b67b4
 Configuration time       | 2013-06-04 13:39:05.720816


   Content of the .dal input file
 ----------------------------------

**DALTON                                          
.RUN RESPONSE                                     
*PCM                                              
.NEQRSP                                           
.SOLVNT                                           
WATER                                             
.NPCMMT                                           
0                                                 
.NESFP                                            
8                                                 
.ICESPH                                           
2                                                 
*PCMCAV                                           
.RIN                                              
1.7                                               
1.5                                               
1.2                                               
1.2                                               
1.5                                               
1.5                                               
1.2                                               
1.2                                               
.INA                                              
1                                                 
2                                                 
3                                                 
4                                                 
5                                                 
6                                                 
7                                                 
8                                                 
.AREATS                                           
 0.3                                              
.ALPHA                                            
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
*QM3                                              
.PRINT                                            
 1                                                
.QM3                                              
.MMPCM                                            
.THRSMP                                           
 1.0e-9                                           
**INTEGRALS                                       
.DIPLEN                                           
.NUCPOT                                           
.NELFLD                                           
.PRINT                                            
 2                                                
**WAVE FUNCTIONS                                  
.DFT                                              
B3LYP                                             
*SCF INPUT                                        
.THRESH                                           
1.0e-12                                           
**RESPONSE                                        
*LINEAR                                           
.DIPLEN                                           
.SINGLE RES                                       
.ROOTS                                            
4                                                 
.THCPP                                            
 1.0e-04                                          
.PRINT                                            
5                                                 
**END OF DALTON INPUT                             


   Content of the .mol file
 ----------------------------

ATOMBASIS                                                                      
CH2O PLUS 2 WATERS                                                             
------------------------                                                       
AtomTypes=5 NoSymmetry Angstrom                                                
        Charge=6.0   Atoms=1    Basis=STO-3G                                   
C           -1.588367    -.770650     .029109 0 1                              
        Charge=8.0   Atoms=1    Basis=STO-3G                                   
O           -1.657083     .436069    -.009750 0 2                              
        Charge=1.0   Atoms=2    Basis=STO-3G                                   
H           -.620668   -1.294822      .054251 0 3                              
H           -2.508043   -1.382001     .040282 0 4                              
   Charge=-0.669     Atoms=2   Basis=MM                                        
O             .975536    1.507219    -.082590 1 1                              
O            1.633705   -1.233280    -.078786 2 1                              
    Charge=0.3345    Atoms=4   Basis=MM                                        
H             .023200    1.300096    -.088063 1 2                              
H            1.113831    2.040431     .704778 1 3                              
H            1.600739    -.258639    -.086953 2 2                              
H            2.303878   -1.461228     .570065 2 3                              


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

     ICESPH =       2     NESFP =       8
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =T
     POLYG          60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.7000
     2    0.0000    0.0000    0.0000    1.5000
     3    0.0000    0.0000    0.0000    1.2000
     4    0.0000    0.0000    0.0000    1.2000
     5    0.0000    0.0000    0.0000    1.5000
     6    0.0000    0.0000    0.0000    1.5000
     7    0.0000    0.0000    0.0000    1.2000
     8    0.0000    0.0000    0.0000    1.2000


 Changes of defaults for *QM3   :
 --------------------------------

 +------------------+
 |  WORD: | CHANGE: |
 +------------------+
 |    QM3 |       T |
 |  PRINT |       1 |
 +------------------+
  Settings for determination of induced dipoles:
  Iterative Method is used
  MM/PCM interface active
 +------------------+
  WARNING: QM3 .INTDIR activated since this is a QM/MM/PCM calculation


    *************************************************************************
    *************** Output from MM potential input processing ***************
    *************************************************************************



 ------------------------------------------------------------------------
 | QM-sys type: | Systems:    | Model:  | Electric properties:          |
 ------------------------------------------------------------------------
 |          0   | [   0;   0] | SPC     |     No classical polarization |
 ------------------------------------------------------------------------


 ------------------------------------------------------------------------
 | MM-sys type: | Systems:    | Model:  | Electric properties:          |
 ------------------------------------------------------------------------
 |          1   | [   1;   2] | SPC_E01 | Isotropic polarizability:     |
 |              |             |         | alp( 1)=  9.5010              |
 ------------------------------------------------------------------------




   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************


    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1: CH2O PLUS 2 WATERS                                                      
 2: ------------------------                                                
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
     "/home/arnfinn/dalton/DALTON/build-trace/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   8 :
     "/home/arnfinn/dalton/DALTON/build-trace/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   1 :
     "/home/arnfinn/dalton/DALTON/build-trace/basis/STO-3G"

  Atomic type no.    4
  --------------------
  Nuclear charge:  -0.66900
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  This is an MM atom without basis functions.

  Atomic type no.    5
  --------------------
  Nuclear charge:   0.33450
  Number of symmetry independent centers:    4
  Number of basis sets to read;    2
  This is an MM atom without basis functions.


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1   -3.0015786160e+00   -1.4563174382e+00    5.5008037772e-02    1.7000000000e+00
   2   -3.1314330364e+00    8.2405098160e-01   -1.8424829719e-02    1.5000000000e+00
   3   -1.1728925345e+00   -2.4468589606e+00    1.0251953201e-01    1.2000000000e+00
   4   -4.7395143797e+00   -2.6116033945e+00    7.6121947767e-02    1.2000000000e+00
   5    1.8434958651e+00    2.8482311204e+00   -1.5607248066e-01    1.5000000000e+00
   6    3.0872550190e+00   -2.3305614354e+00   -1.4888396248e-01    1.5000000000e+00
   7    4.3841646100e-02    2.4568253762e+00   -1.6641495175e-01    1.2000000000e+00
   8    2.1048355395e+00    3.8558557669e+00    1.3318373989e+00    1.2000000000e+00


                                 Isotopic Masses
                                 ---------------

                           C          12.000000
                           O          15.994915
                           H           1.007825
                           H           1.007825

                       Total mass:    30.010565 amu
                       Natural abundance:  98.633 %

 Center-of-mass coordinates (a.u.):   -3.067740   -0.312997    0.018175


  Atoms and basis sets
  --------------------

  Number of atom types :    6
  Total number of atoms:   12

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C           1    6.0000    15     5      [6s3p|2s1p]                                        
  O           1    8.0000    15     5      [6s3p|2s1p]                                        
  H           2    1.0000     3     1      [3s|1s]                                            
  O           2   -0.6690     0     0      Point Charge                                       
  H           4    0.3345     0     0      Point Charge                                       
  a           2    0.0000     0     0      Point Charge                                       
  ----------------------------------------------------------------------
  total:     12   16.0000    36    12
  ----------------------------------------------------------------------

  Threshold for neglecting AO integrals:  1.00e-12


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   36
  C       :     1  x  -3.0015786160    2  y  -1.4563174382    3  z   0.0550080378
  O       :     4  x  -3.1314330364    5  y   0.8240509816    6  z  -0.0184248297
  H       :     7  x  -1.1728925345    8  y  -2.4468589606    9  z   0.1025195320
  H       :    10  x  -4.7395143797   11  y  -2.6116033945   12  z   0.0761219478


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H           O           O     
            ------      ------      ------      ------      ------      ------
 C     :    0.000000
 O     :    1.209298    0.000000
 H     :    1.100831    2.018474    0.000000
 H     :    1.104391    2.007988    1.889439    0.000000
 O     :    3.431437    2.843123    3.227697    4.527473    0.000000
 O     :    3.256903    3.690634    2.259133    4.146127    2.818428    0.000000
 H     :    2.626569    1.891038    2.677390    3.690164    0.974615    3.001967
 H     :    3.957345    3.280624    3.815178    5.027183    0.960931    3.406086
 H     :    3.232031    3.331964    2.455252    4.261481    1.873273    0.975233
 H     :    3.989875    4.430026    2.974344    4.841645    3.316946    0.960260
 a  1  :    3.408645    2.807232    3.219816    4.503294    0.065744    2.847367
 a  2  :    3.285589    3.703418    2.294418    4.182574    2.786761    0.065593

            H           H           H           H           a  1        a  2  
            ------      ------      ------      ------      ------      ------
 H     :    0.000000
 H     :    1.538236    0.000000
 H     :    2.217721    2.479847    0.000000
 H     :    3.641364    3.700807    1.540226    0.000000
 a  1  :    0.935667    0.923012    1.906635    3.343454    0.000000
 a  2  :    2.986623    3.363487    0.936414    0.922570    2.815774    0.000000


  Max    interatomic separation is    2.0185 Angstrom (    3.8144 Bohr)
  between atoms    3 and    2, "H     " and "O     ".

  Min HX interatomic separation is    1.1008 Angstrom (    2.0803 Bohr)

  Min YX interatomic separation is    1.2093 Angstrom (    2.2852 Bohr)

  Max QM+MM interatomic separation is    5.0272 Angstrom (    9.5000 Bohr)
  between the QM+MM centers    8 and    4, "H     " and "H     ".




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       1.798871         -0.056827    0.997867   -0.032124
   IB      13.009178          0.998345    0.057081    0.007042
   IC      14.808049         -0.008861    0.031671    0.999459


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         280942.3015          38847.8806          34128.6681 MHz
            9.371226            1.295826            1.138410 cm-1


@  Nuclear repulsion energy :   31.249215315972 Hartree
      QM3 induced dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.020409711860107


    **************************************************************************
    ***************** The MM/MM classical interaction energy *****************
    **************************************************************************


 ------------------------------------------------------------------------
 |  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        |       -0.00626437 |
 |  Epol  = - 1/2*Sum_a[ Pind_a*E^site_a ]          |       -0.00096177 |
 |  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] |        0.00000000 |
 ------------------------------------------------------------------------
 |  E(MM/MM) = Eelec + Epol + Evdw                  |       -0.00722614 |
 ------------------------------------------------------------------------




     ************************************************************************
     *************** The "QM"/MM classical interaction energy ***************
     ************************************************************************


 ------------------------------------------------------------------------
 |  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        |        0.00000000 |
 |  Epol  = - 1/2*Sum_a[ Pind_a*E^(QMclassic)_a ]   |        0.00000000 |
 |  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] |        0.00000000 |
 ------------------------------------------------------------------------
 |  E("QM"/MM) = Eelec + Epol + Evdw                |        0.00000000 |
 ------------------------------------------------------------------------




                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************


 Default print level:        2

 * Nuclear model: Point charge

 Calculation of one- and two-electron Hamiltonian integrals.

 The following one-electron property integrals are calculated as requested:
          - overlap integrals
          - dipole length integrals
          - Electric field at the nuclei
          - Potential energy at the nuclei

 Center of mass  (bohr):     -3.067740315293     -0.312997395263      0.018174673399
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERONE **************************
     ************************************************************************



     ************************************************************************
     ************************** Output from TWOINT **************************
     ************************************************************************


 Threshold for neglecting two-electron integrals:  1.00e-12
 Number of two-electron integrals written:        2907 ( 94.4% )
 Megabytes written:                              0.034


 MEMORY USED TO GENERATE CAVITY =    432042

 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.

 Total number of spheres =   11
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1   -1.588367000   -0.770650000    0.029109000    2.040000000   22.843745166
   2   -1.657083000    0.436069000   -0.009750000    1.800000000   17.320494799
   3   -0.620668000   -1.294822000    0.054251000    1.440000000    4.504430717
   4   -2.508043000   -1.382001000    0.040282000    1.440000000    9.471926497
   5    0.975536000    1.507219000   -0.082590000    1.800000000   16.946491810
   6    1.633705000   -1.233280000   -0.078786000    1.800000000   28.188763663
   7    0.023200000    1.300096000   -0.088063000    1.440000000    4.066167360
   8    1.113831000    2.040431000    0.704778000    1.440000000   11.882478847
   9    1.304620500    0.136969500   -0.080688000    1.471280988    4.271467524
  10   -0.012749934    0.270009434   -0.049134665    1.232658185    0.798536321
  11    0.329445664   -0.588244950   -0.012340260    1.188874163    0.521620751

 Total number of tesserae =     579
 Surface area =  120.81612345 (A^2)    Cavity volume =   96.50955728 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....
 
  ..... DONE GENERATING -Q-  MATRIX .....
 >>>  Time used in Q-MAT      is   0.30 seconds
 >>>> Total CPU  time used in HERMIT:   0.67 seconds
 >>>> Total wall time used in HERMIT:   0.70 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          F
     Using EHT  model:          T
     Number of Huckel orbitals each symmetry:   12

 Huckel EHT eigenvalues for symmetry :  1
          -20.705291     -11.377081      -1.496925      -0.962689      -0.743251
           -0.672270      -0.593413      -0.504635      -0.350830      -0.253214
           -0.198302      -0.193199

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Wed Jun  5 20:39:20 2013
     Host name              : c18-12.local                            

 Title lines from ".mol" input file:
     CH2O PLUS 2 WATERS                                                      
     ------------------------                                                

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

@    Restricted, closed shell Kohn-Sham DFT calculation.

@    Time-dependent Kohn-Sham DFT calculation (TD-DFT).

 Initial molecular orbitals are obtained according to
 ".MOSTART EHT   " input option

@    QM part is embedded in an environment :

@         Solvation model: PCM

@         Model: QM3

     Wave function specification
     ============================
     For the wave function of type :      >>> KS-DFT   <<<
     Number of closed shell electrons          16
     Number of electrons in active shells       0
     Total charge of the molecule               0

     Spin multiplicity and 2 M_S                1         0
     Total number of symmetries                 1
     Reference state symmetry                   1
 
     This is a DFT calculation of type: B3LYP
 Weighted mixed functional:
               HF exchange:    0.20000
                       VWN:    0.19000
                       LYP:    0.81000
                     Becke:    0.72000
                    Slater:    0.80000

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1
                                       --- |  ---
     Occupied SCF orbitals               8 |    8
     Secondary orbitals                  4 |    4
     Total number of orbitals           12 |   12
     Number of basis functions          12 |   12

     Optimization information
     ========================
@    Number of configurations                 1
@    Number of orbital rotations             32
     ------------------------------------------
@    Total number of variables               33

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-12
 
     This is a DFT calculation of type: B3LYP
 Weighted mixed functional:
               HF exchange:    0.20000
                       VWN:    0.19000
                       LYP:    0.81000
                     Becke:    0.72000
                    Slater:    0.80000

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.99085
 NUCLEAR APPARENT CHARGE -15.79248
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.975378044871  16.0000012547    1.25e-06
      QM3 induced Dipole vector converged in  12 iterations.
      Final norm2 of QM3 induced dipole moment vector:   13.241662292657828
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.688237384600798
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.154008077632154
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128059711130834
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126706885498394
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126664581086266
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126675731044768
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126680482362685
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126681854288569
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126682208816078
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126682295680977
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.126682316680574
 MMPCM converged for given density
 No. of iterations:           12
   1  -112.645865418     -4.056810716836e-02    2.41860e+00   -1.13e+02
      Virial theorem: -V/T =      2.004627
@      MULPOP C       0.99; O      -0.72; H      -0.14; H      -0.14; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.758841657192  16.0000013793    1.38e-06
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.040852670815576
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.082429011718478
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.092641306212917
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.094365473198509
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.094641566827926
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.094684876851676
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.094691499416061
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.094692466567628
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.094692595088618
 MMPCM converged for given density
 No. of iterations:            9
   2  -112.409148132     -2.655002900105e-02    2.58352e+00    2.37e-01
      Virial theorem: -V/T =      2.015001
@      MULPOP C      -1.40; O       0.98; H       0.22; H       0.20; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.945237489065  16.0000013047    1.30e-06
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.247060242179387
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.144115822359321
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.130845343274184
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128835970470722
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128524651610388
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128476881582408
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128469744480428
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128468733006708
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.128468605031618
 MMPCM converged for given density
 No. of iterations:            9
   3  -112.860218149     -3.358229838273e-02    8.17454e-01   -4.51e-01
      Virial theorem: -V/T =      2.009390
@      MULPOP C       0.44; O      -0.59; H       0.08; H       0.07; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.866031151391  16.0000013131    1.31e-06
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.092690040469251
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.111670447582247
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.114900961587429
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.115415737318163
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.115496282101365
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.115508662368046
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.115510508670057
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.115510768526119
 MMPCM converged for given density
 No. of iterations:            8
   4  -112.912976030     -2.013543583818e-02    3.84248e-02   -5.28e-02
      Virial theorem: -V/T =      2.014163
@      MULPOP C       0.03; O      -0.18; H       0.08; H       0.07; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.871395063283  16.0000013149    1.31e-06
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.118063592631383
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116807082007326
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116609951901332
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116579117930852
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116574351946044
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116573631705578
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116573527424380
 MMPCM converged for given density
 No. of iterations:            7
   5  -112.913129983     -2.074773907698e-02    7.70486e-03   -1.54e-04
      Virial theorem: -V/T =      2.013787
@      MULPOP C       0.02; O      -0.19; H       0.10; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870742165816  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116277515305010
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116407296440637
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116427686861743
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116430863648249
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116431351250479
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116431424034905
 MMPCM converged for given density
 No. of iterations:            6
   6  -112.913132247     -2.068018558522e-02    2.01874e-03   -2.26e-06
      Virial theorem: -V/T =      2.013841
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870864699476  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116494528436151
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116457755803395
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116451932505411
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116451013826730
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116450869896947
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116450847700309
 MMPCM converged for given density
 No. of iterations:            6
   7  -112.913133038     -2.069724423336e-02    2.19851e-04   -7.91e-07
      Virial theorem: -V/T =      2.013834
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870878973298  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116458352315280
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116454039829781
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453357418978
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453249824789
 MMPCM converged for given density
 No. of iterations:            4
   8  -112.913133103     -2.069924364304e-02    1.62335e-05   -6.49e-08
      Virial theorem: -V/T =      2.013833
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870880120846  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453814467184
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453475322172
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453421711052
 MMPCM converged for given density
 No. of iterations:            3
   9  -112.913133107     -2.069940105393e-02    5.63273e-07   -4.67e-09
      Virial theorem: -V/T =      2.013833
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870880169298  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453440233722
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453423970450
 MMPCM converged for given density
 No. of iterations:            2
  10  -112.913133108     -2.069940787605e-02    7.92254e-08   -1.61e-10
      Virial theorem: -V/T =      2.013833
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870880164192  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453418563726
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419748875
 MMPCM converged for given density
 No. of iterations:            2
  11  -112.913133108     -2.069940711281e-02    3.74246e-09    5.48e-11
      Virial theorem: -V/T =      2.013833
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870880163773  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419721794
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419884980
 MMPCM converged for given density
 No. of iterations:            2
  12  -112.913133108     -2.069940705857e-02    1.94440e-09   -7.11e-14
      Virial theorem: -V/T =      2.013833
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870880163913  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419985015
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419949148
 MMPCM converged for given density
 No. of iterations:            2
  13  -112.913133108     -2.069940707837e-02    3.12236e-12   -9.24e-13
      Virial theorem: -V/T =      2.013833
@      MULPOP C       0.03; O      -0.19; H       0.09; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870880163913  16.0000013143    1.31e-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419943178
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.116453419942232
 MMPCM converged for given density
 No. of iterations:            2
  14  -112.913133108     -2.069940707830e-02    2.16383e-13   -1.42e-14

 *** DIIS converged in  14 iterations !
@    Converged SCF energy, gradient:   -112.913133107583    2.16e-13
   - total time used in SIRFCK :              0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Sym       Kohn-Sham orbital energies

  1    -18.91775721   -10.08361920    -1.01243288    -0.57586966    -0.44370785
        -0.34959224    -0.31992888    -0.16496465     0.07089260     0.39461218
         0.49760331     0.61482162

    E(LUMO) :     0.07089260 au (symmetry 1)
  - E(HOMO) :    -0.16496465 au (symmetry 1)
  ------------------------------------------
    gap     :     0.23585724 au

 >>> Writing SIRIFC interface file <<<

 >>>> CPU and wall time for SCF :      45.070      46.185


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


@    Spin multiplicity:           1
@    Spatial symmetry:            1
@    Total charge of molecule:    0

     QM/MM/PCM calculation converged :

        Electrostatic energy:      -0.007669208901
        Polarization energy:       -0.000833418285
        van der Waals energy:       0.000000000000
        QM/MM/PCM energy:          -0.020699407078

@    Final DFT energy:           -112.913133107583                 
@    Nuclear repulsion:            31.249215315972
@    Electronic energy:          -144.141649016477

@    Final gradient norm:           0.000000000000

 
     Date and time (Linux)  : Wed Jun  5 20:40:06 2013
     Host name              : c18-12.local                            

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

    Orbital         1        2        3        4        5        6        7
   1 C   :1s     0.0007   0.9910  -0.1288  -0.1866  -0.0015   0.0398   0.0001
   2 C   :1s    -0.0089   0.0406   0.2850   0.5871   0.0063  -0.1264  -0.0004
   3 C   :2px    0.0005  -0.0001  -0.0089   0.0120   0.5756   0.0127  -0.0050
   4 C   :2py   -0.0078   0.0011   0.1690  -0.2626   0.0215  -0.4227   0.0177
   5 C   :2pz    0.0003  -0.0000  -0.0056   0.0081   0.0042   0.0119   0.6009
   6 O   :1s     0.9935   0.0005  -0.2150   0.0979  -0.0003  -0.1120  -0.0004
   7 O   :1s     0.0296  -0.0081   0.7390  -0.4141   0.0034   0.5647   0.0021
   8 O   :2px    0.0004  -0.0002   0.0175  -0.0011   0.3297  -0.0297  -0.0064
   9 O   :2py   -0.0064   0.0031  -0.2117  -0.1028   0.0272   0.6802   0.0236
  10 O   :2pz    0.0002  -0.0001   0.0066   0.0034   0.0017  -0.0243   0.6814
  11 H   :1s     0.0003  -0.0088   0.0314   0.2612   0.3170   0.1307   0.0008
  12 H   :1s     0.0003  -0.0088   0.0317   0.2654  -0.3168   0.1466   0.0003

    Orbital         8        9       10
   1 C   :1s    -0.0002  -0.0000   0.2051
   2 C   :1s     0.0000  -0.0002  -1.2684
   3 C   :2px   -0.0892  -0.0074   0.0170
   4 C   :2py    0.0009   0.0262   0.4644
   5 C   :2pz   -0.0002   0.8275  -0.0148
   6 O   :1s     0.0032   0.0001  -0.0225
   7 O   :1s    -0.0144  -0.0002   0.1275
   8 O   :2px    0.9060   0.0065  -0.0011
   9 O   :2py    0.0438  -0.0246  -0.2101
  10 O   :2pz    0.0067  -0.7626   0.0067
  11 H   :1s    -0.3411   0.0004   0.8724
  12 H   :1s     0.3467   0.0000   0.9182



 >>>> Total CPU  time used in SIRIUS :     45.09 seconds
 >>>> Total wall time used in SIRIUS :     46.25 seconds

 
     Date and time (Linux)  : Wed Jun  5 20:40:06 2013
     Host name              : c18-12.local                            

 NOTE:    1 warnings have been issued.
 Check output, result, and error files for "WARNING".


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

 QM/MM response calculation
 Non-equilibrium PCM solvent model requested    : INERSI =T

 Static dielectric constant                     : EPSTAT = 78.3900
 Optical dielectric constant                    : EPSOL  =  1.7760
 Print level                                    : IPRPP  =   5
 Maximum number of iterations for eigenval.eqs. : MAXITP =  60
 Threshold for convergence of eigenvalue eqs.   : THCPP  = 1.000e-04
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

      4 Excitation energies are calculated for symmetry no.    1

      3 property residues are calculated with labels:

               XDIPLEN 
               YDIPLEN 
               ZDIPLEN 

 Integral transformation: Total CPU and WALL times (sec)       0.007       0.011


   SCF energy         :     -112.913133107582794
 -- inactive part     :     -144.282662539261082
 -- nuclear repulsion :       31.249215315972382


                    *****************************************
                    *** DFT response calculation (TD-DFT) ***
                    *****************************************



 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            4
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      32
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      32
 Electrons in DFTMOMO:   16.00000131428623


                         Matrix < abs(phi_p) abs(phi_q) >
                         --------------------------------


               Column   1     Column   2     Column   3     Column   4     Column   5     Column   6     Column   7     Column   8
       1       1.00000001     0.00199563     0.21049475     0.09704217     0.07746212     0.19306334     0.15104992     0.19169152
       2       0.00199563     1.00000055     0.13410930     0.19022437     0.14928021     0.11841411     0.14568113     0.03200435
       3       0.21049475     0.13410930     1.00000005     0.46882440     0.60517410     0.70956035     0.79309361     0.69658954
       4       0.09704217     0.19022437     0.46882440     1.00000004     0.82603197     0.64890286     0.59765571     0.61780758
       5       0.07746212     0.14928021     0.60517410     0.82603197     1.00000004     0.54917202     0.52788194     0.73762333
       6       0.19306334     0.11841411     0.70956035     0.64890286     0.54917202     1.00000000     0.62051646     0.64696297
       7       0.15104992     0.14568113     0.79309361     0.59765571     0.52788194     0.62051646     0.99999999     0.51731898
       8       0.19169152     0.03200435     0.69658954     0.61780758     0.73762333     0.64696297     0.51731898     0.99999998
       9       0.15737825     0.18184850     0.66782058     0.67864879     0.53282948     0.59200082     0.94724205     0.48898593
      10       0.05023827     0.21552930     0.32672900     0.82866151     0.71843100     0.44793287     0.41830817     0.49073935
      11       0.07095010     0.22287843     0.47495003     0.72821385     0.87564171     0.49912794     0.43756232     0.60382383
      12       0.21170599     0.24686412     0.76908519     0.63086162     0.61459841     0.71054325     0.67271817     0.51171038

               Column   9     Column  10     Column  11     Column  12
       1       0.15737825     0.05023827     0.07095010     0.21170599
       2       0.18184850     0.21552930     0.22287843     0.24686412
       3       0.66782058     0.32672900     0.47495003     0.76908519
       4       0.67864879     0.82866151     0.72821385     0.63086162
       5       0.53282948     0.71843100     0.87564171     0.61459841
       6       0.59200082     0.44793287     0.49912794     0.71054325
       7       0.94724205     0.41830817     0.43756232     0.67271817
       8       0.48898593     0.49073935     0.60382383     0.51171038
       9       1.00000000     0.48593365     0.46884117     0.66073088
      10       0.48593365     1.00000006     0.79225371     0.51403404
      11       0.46884117     0.79225371     1.00000008     0.55454844
      12       0.66073088     0.51403404     0.55454844     1.00000001
    ==== End of matrix output ====


                            Matrix < phi_p^2 phi_q^2 >
                            --------------------------


               Column   1     Column   2     Column   3     Column   4     Column   5     Column   6     Column   7     Column   8
       1       1.00000000     0.00000054     0.42517117     0.14436899     0.01409204     0.14000110     0.04222008     0.05259723
       2       0.00000054     1.00000000     0.09981674     0.33713965     0.03364779     0.01866632     0.02442532     0.00056562
       3       0.42517117     0.09981674     1.00000000     0.22194090     0.27665116     0.41516551     0.51994416     0.46482337
       4       0.14436899     0.33713965     0.22194090     1.00000000     0.62694280     0.33309604     0.29378103     0.23872927
       5       0.01409204     0.03364779     0.27665116     0.62694280     1.00000000     0.18797521     0.21578874     0.42243850
       6       0.14000110     0.01866632     0.41516551     0.33309604     0.18797521     1.00000000     0.33620819     0.35120199
       7       0.04222008     0.02442532     0.51994416     0.29378103     0.21578874     0.33620819     1.00000000     0.28941838
       8       0.05259723     0.00056562     0.46482337     0.23872927     0.42243850     0.35120199     0.28941838     1.00000000
       9       0.04438798     0.03811779     0.38916388     0.36959697     0.23400071     0.31428655     0.94181985     0.26134796
      10       0.00940404     0.32148010     0.10177585     0.74605323     0.49563550     0.11372446     0.17058326     0.14329171
      11       0.00967013     0.06664709     0.15722939     0.44556732     0.80145029     0.13671479     0.17038184     0.23716718
      12       0.16396678     0.10666663     0.66820503     0.31075417     0.29499923     0.43963014     0.39705907     0.30036156

               Column   9     Column  10     Column  11     Column  12
       1       0.04438798     0.00940404     0.00967013     0.16396678
       2       0.03811779     0.32148010     0.06664709     0.10666663
       3       0.38916388     0.10177585     0.15722939     0.66820503
       4       0.36959697     0.74605323     0.44556732     0.31075417
       5       0.23400071     0.49563550     0.80145029     0.29499923
       6       0.31428655     0.11372446     0.13671479     0.43963014
       7       0.94181985     0.17058326     0.17038184     0.39705907
       8       0.26134796     0.14329171     0.23716718     0.30036156
       9       1.00000000     0.24007044     0.21007406     0.36759943
      10       0.24007044     1.00000000     0.57709756     0.23828876
      11       0.21007406     0.57709756     1.00000000     0.27086793
      12       0.36759943     0.23828876     0.27086793     1.00000000
    ==== End of matrix output ====

 >>> IN RSPPP:
  THCPP,MAXRM   1.000000000000000e-004         600
  KSYMOP,NGPPP(KSYMOP)            1           3
  LWRK ,LWRK1     63999526    62917438
  KEXCNV,NSIM,LWRK2            4           4     1082605



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F

           0  START configuration VECTORS USING LOWEST DIAGONAL HESSIAN ELEMENTS
           4  START orbital VECTORS

 ** RSPCTL MICROITERATION NUMBER    1
       Electrons: 16.000001( 1.31e-06): LR-DFT*4 evaluation time:       1.5 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.40267e-02  0.00e+00  1.40e-02  7.09e-01    1.51007e-01
         2    1.04260e-01  0.00e+00  1.04e-01  7.09e-01    3.55293e-01
         3    9.02229e-03  0.00e+00  9.02e-03  7.07e-01    4.36525e-01
         4    4.70976e-01  0.00e+00  4.71e-01  7.55e-01    5.31818e-01

 ** RSPCTL MICROITERATION NUMBER    2
       Electrons: 16.000001( 1.31e-06): LR-DFT*4 evaluation time:       1.5 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.03054e-04  0.00e+00  1.03e-04  7.09e-01    1.50894e-01
         2    4.65491e-03  0.00e+00  4.65e-03  7.09e-01    3.49462e-01
         3    1.32931e-03  0.00e+00  1.33e-03  7.07e-01    4.36407e-01
         4    3.73315e-02  0.00e+00  3.73e-02  7.29e-01    4.37184e-01

 ** RSPCTL MICROITERATION NUMBER    3
       Electrons: 16.000001( 1.31e-06): LR-DFT*4 evaluation time:       1.5 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    9.84198e-06  0.00e+00  9.84e-06  7.09e-01    1.50894e-01
         2    1.62865e-04  0.00e+00  1.63e-04  7.09e-01    3.49459e-01
         3    9.91389e-03  0.00e+00  9.91e-03  7.28e-01    4.36267e-01
         4    1.86853e-03  0.00e+00  1.87e-03  7.08e-01    4.36413e-01

 ** RSPCTL MICROITERATION NUMBER    4
       Electrons: 16.000001( 1.31e-06): LR-DFT*3 evaluation time:       1.1 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    8.86598e-06  0.00e+00  8.87e-06  7.09e-01    1.50894e-01
         2    7.23080e-06  0.00e+00  7.23e-06  7.09e-01    3.49459e-01
         3    1.44808e-03  0.00e+00  1.45e-03  7.29e-01    4.36216e-01
         4    1.98095e-04  0.00e+00  1.98e-04  7.08e-01    4.36412e-01

 ** RSPCTL MICROITERATION NUMBER    5
       Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.8 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    9.03441e-06  0.00e+00  9.03e-06  7.09e-01    1.50894e-01
         2    4.71304e-06  0.00e+00  4.71e-06  7.09e-01    3.49459e-01
         3    9.72110e-05  0.00e+00  9.72e-05  7.29e-01    4.36215e-01
         4    1.29840e-05  0.00e+00  1.30e-05  7.08e-01    4.36412e-01

 ** RSPCTL MICROITERATION NUMBER    6
       Electrons: 16.000001( 1.31e-06): LR-DFT*1 evaluation time:       0.4 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    8.82179e-06  0.00e+00  8.82e-06  7.09e-01    1.50894e-01
         2    4.42692e-06  0.00e+00  4.43e-06  7.09e-01    3.49459e-01
         3    3.74365e-05  0.00e+00  3.74e-05  7.29e-01    4.36215e-01
         4    4.45791e-06  0.00e+00  4.46e-06  7.08e-01    4.36412e-01

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   36)
 RSP solution vector no.    1; norm of residual   1.24e-05
 RSP solution vector no.    2; norm of residual   6.24e-06
 RSP solution vector no.    3; norm of residual   5.14e-05
 RSP solution vector no.    4; norm of residual   6.30e-06

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    XDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -3.95725849e-04 *ENERGY(eV):   4.1060393    
@ STATE NO:    2 *TRANSITION MOMENT:  1.65899787e-03 *ENERGY(eV):   9.5092759    
@ STATE NO:    3 *TRANSITION MOMENT:  4.27987055e-02 *ENERGY(eV):   11.870025    
@ STATE NO:    4 *TRANSITION MOMENT:  5.81860300e-03 *ENERGY(eV):   11.875364    

@ Transition operator type:    YDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -9.19588871e-05 *ENERGY(eV):   4.1060393    
@ STATE NO:    2 *TRANSITION MOMENT: -6.43359696e-03 *ENERGY(eV):   9.5092759    
@ STATE NO:    3 *TRANSITION MOMENT: -0.66548695     *ENERGY(eV):   11.870025    
@ STATE NO:    4 *TRANSITION MOMENT: -9.20759902e-02 *ENERGY(eV):   11.875364    

@ Transition operator type:    ZDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  5.03948873e-03 *ENERGY(eV):   4.1060393    
@ STATE NO:    2 *TRANSITION MOMENT: -0.17976795     *ENERGY(eV):   9.5092759    
@ STATE NO:    3 *TRANSITION MOMENT:  2.19058394e-02 *ENERGY(eV):   11.870025    
@ STATE NO:    4 *TRANSITION MOMENT: -4.37888474e-05 *ENERGY(eV):   11.875364    


  ******************************************************************************
  *** @ Excit. operator sym 1 & ref. state sym 1 => excited state symmetry 1 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.15089417    au
@                      4.1060393     eV
@                      33117.443     cm-1
                       396.17259     kJ / mol

@ Total energy :      -112.76224     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  1.57532457e-08  (Transition moment : -3.95725849e-04 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  8.50684702e-10  (Transition moment : -9.19588871e-05 )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  2.55478387e-06  (Transition moment :  5.03948873e-03 )

 Eigenvector for state no.  1

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          29         8(1)   9(1)        0.7077518876        0.0375742150

       31 elements with absolute value less than 7.08e-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      8    9  0.707752  0.037574  0.488986  0.261348  0.449138  0.219622

@ Overlap diagnostic LAMBDA =    0.4890


 @ Excited state no:    2 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.34945947    au
@                      9.5092759     eV
@                      76697.488     cm-1
                       917.50570     kJ / mol

@ Total energy :      -112.56367     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  6.41205461e-07  (Transition moment :  1.65899787e-03 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  9.64302419e-06  (Transition moment : -6.43359696e-03 )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  7.52887516e-03  (Transition moment : -0.17976795     )

 Eigenvector for state no.  2

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          21         6(1)   9(1)        0.7056494895        0.0321398851

       31 elements with absolute value less than 7.06e-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      6    9  0.705649  0.032140  0.592001  0.314287  0.453615  0.268541

@ Overlap diagnostic LAMBDA =    0.5932


 @ Excited state no:    3 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.43621541    au
@                      11.870025     eV
@                      95738.215     cm-1
                       1145.2834     kJ / mol

@ Total energy :      -112.47692     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  5.32685666e-04  (Transition moment :  4.27987055e-02 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  0.12879198      (Transition moment : -0.66548695     )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  1.39549904e-04  (Transition moment :  2.19058394e-02 )

 Eigenvector for state no.  3

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          17         5(1)   9(1)        0.0967565470        0.0010580612
          24         6(1)  12(1)        0.1272163746       -0.0421839814
          25         7(1)   9(1)        0.6677406674        0.1065742950
          31         8(1)  11(1)       -0.1841357707        0.0288508750

       28 elements with absolute value less than 6.68e-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      6   12  0.127216 -0.042184  0.710543  0.439630  0.028696  0.020390
      7    9  0.667741  0.106574  0.947242  0.941820  0.314908  0.298294
      8   11 -0.184136  0.028851  0.603824  0.237167  0.045363  0.027391

@ Overlap diagnostic LAMBDA =    0.8766


 @ Excited state no:    4 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.43641161    au
@                      11.875364     eV
@                      95781.276     cm-1
                       1145.7985     kJ / mol

@ Total energy :      -112.47672     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  9.85014196e-06  (Transition moment :  5.81860300e-03 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  2.46659492e-03  (Transition moment : -9.20759902e-02 )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  5.57868788e-10  (Transition moment : -4.37888474e-05 )

 Eigenvector for state no.  4

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          17         5(1)   9(1)       -0.6999752575       -0.0078165517
          25         7(1)   9(1)        0.0922812537        0.0147534262

       30 elements with absolute value less than 7.00e-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      5    9 -0.699975 -0.007817  0.532829  0.234001  0.479084  0.255270

@ Overlap diagnostic LAMBDA =    0.5384


 Time used in polarization propagator calculation is     11.44 CPU seconds for symmetry 1

 >>>> Total CPU  time used in RESPONSE:  11.45 seconds
 >>>> Total wall time used in RESPONSE:  11.80 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  57.24 seconds
 >>>> Total wall time used in DALTON:  58.91 seconds

 
     Date and time (Linux)  : Wed Jun  5 20:40:18 2013
     Host name              : c18-12.local                            

END REFOUT

