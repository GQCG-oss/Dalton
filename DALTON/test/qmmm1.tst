########## Test description ########################
START DESCRIPTION
KEYWORDS: qmmm dft dipole properties
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
qmmmaniso
qmmmconvergence
qmmmq
qmmmdip
qmmmquad
qmmmelpol
qmmmnucpol
qmmmmulpol
qm3energy
enedft
dipole
dipcompx
dipcompy
dipcompz
excita
excita
OVERRIDE line 9
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON
.RUN PROPERTIES
*QMMM
.QMMM
.MMITER
.PRINT
 1 
.MMPROP
**WAVE FUNCTIONS
.DFT
B3LYP
*SCF INPUT
.THRESHOLD
 1.0D-10
**PROPERTIES
.EXCITA
*EXCITA
.NEXCIT
 2
.DIPSTR
.THRESH
 1.0D-7
**END OF
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS
CH2O PLUS 2 WATERS
------------------------
AtomTypes=3 NoSymmetry Angstrom
        6.0   1    Basis=cc-pVDZ
C           -1.588367    -.770650     .029109 
        8.0   1    Basis=cc-pVDZ
O           -1.657083     .436069    -.009750 
        1.0   2    Basis=cc-pVDZ
H           -.620668   -1.294822      .054251 
H           -2.508043   -1.382001     .040282
END MOLINP

########## POTENTIAL.INP ############################
START POTINP
AU
6 2 2 1
1        1.8434958651        2.8482311204       -0.1560724807       -0.7390827623       -0.1307264223        0.0536972371        0.1276417760       -3.8915470039        0.3042601148        0.1736329543       -4.4593600616        0.4231693516       -4.1915829791        5.4760362505       -0.0866800183       -0.0721680526        5.5822886395       -0.1240156189        5.5133007097
1        0.0438416461        2.4568253762       -0.1664149517        0.3677809440        0.2067300506        0.0345154175       -0.0151419207       -0.0053247942        0.1000149409       -0.0128678147       -0.4962471913       -0.0094448450       -0.5220547903        3.4884652220        0.2984044533       -0.0712037833        1.7393065199       -0.3503076238        1.3605853997
1        2.1048355395        3.8558557669        1.3318373989        0.3713018183       -0.0129776857       -0.1154344939       -0.1761271065       -0.5125105131        0.0301785196        0.0495868218       -0.3480715156        0.2357268893       -0.1533604733        1.1854898700        0.0516528099        0.3136952094        2.3594108855        0.7460384264        2.9489928120
2        3.0872550190       -2.3305614354       -0.1488839625       -0.7392506008        0.1045091928        0.1199578441        0.1050086222       -4.3534010893       -0.2335724574        0.4256499376       -3.8052428657       -0.2036955780       -4.3845193626        5.5552848071        0.0859287441       -0.1208046651        5.4514069390        0.0745648410        5.5662469471
2        3.0249583076       -0.4887568752       -0.1643173557        0.3677014501       -0.0070170988       -0.2093580024       -0.0120115140       -0.5197135585       -0.0029598836       -0.0076502102        0.0139403385        0.0100737937       -0.5191766528        1.5094772028       -0.0367827484       -0.3663623981        3.5546245763        0.1005981894        1.5284304352
2        4.3536984454       -2.7613207262        1.0772667234        0.3715491507       -0.1493905839        0.0323104495       -0.1450960317       -0.2508359582       -0.0754184823        0.2464385390       -0.4973073090       -0.0734587271       -0.2662352700        2.6580864234       -0.2485367227        0.7466547613        1.2237370065       -0.3566527126        2.6080864339
END POTINP

########## Reference Output ########################
START REFOUT


 **************************************************************************************
 ******************** DALTON2011 - An electronic structure program ********************
 **************************************************************************************

    This is output from DALTON Release 2011 (Rev. 0, Mar. 2011)

 --------------------------------------------------------------------------------

    NOTE:
     
    This is an experimental code for the evaluation of molecular
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
    Release DALTON2011 (2011), see http://daltonprogram.org"

 --------------------------------------------------------------------------------

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
  Luca Frediani,            University of Tromsoe,        Norway      (PCM)
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
  Kurt V. Mikkelsen,        University of Copenhagen,     Denmark     (MC-SCRF and QM/MM code)
  Christian B. Nielsen,     University of Copenhagen,     Denmark     (QM/MM code)
  Patrick Norman,           University of Linkoeping,     Sweden      (cubic response and complex response in RESPONS)
  Jeppe Olsen,              Aarhus University,            Denmark     (SIRIUS CI/density modules)
  Anders Osted,             Copenhagen University,        Denmark     (QM/MM code)
  Martin J. Packer,         University of Sheffield,      UK          (SOPPA)
  Thomas B. Pedersen,       University of Oslo,           Norway      (Cholesky decomposition)
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
  Olav Vahtras,             KTH Stockholm,                Sweden      (triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)
  David J. Wilson,          La Trobe University,          Australia   (DFT Hessian and DFT magnetizabilities)
  Hans Agren,               KTH Stockholm,                Sweden      (SIRIUS module, MC-SCRF solvation model)

 --------------------------------------------------------------------------------


     Date and time (Linux)  : Wed Apr  6 09:17:48 2011 
     Host name              : stanley                                 

 * Work memory size             :   100000000 =  762.94 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/jobb/dalton/svn/qmmm_devel/test/2011-04-06T09_17-testjob-pid-973
   2) /home/arnfinn/jobb/dalton/svn/qmmm_devel/basis/


       *******************************************************************
       *********** Output from DALTON general input processing ***********
       *******************************************************************

 --------------------------------------------------------------------------------
   Overall default print level:    0
   Print level for DALTON.ERR :    1

    HERMIT 1- and 2-electron integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed (SIRIUS module)
    Static molecular property section will be executed (ABACUS module)
 --------------------------------------------------------------------------------


 Changes of defaults for *QMMM  :
 --------------------------------

 +------------------+
 |  WORD: | CHANGE: |
 +------------------+
 |   QMMM |       T |
 | MMPROP |       T |
 | MMITER |       T |
 +------------------+


 Induced MM dipoles are solved iteratively
 Max. number of iterations:          100
 Thresshold:   1.00000000000000004E-010



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
  The basis set is "cc-pVDZ" from the basis set library.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/arnfinn/jobb/dalton/svn/qmmm_devel/basis/cc-pVDZ"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/arnfinn/jobb/dalton/svn/qmmm_devel/basis/cc-pVDZ"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Used basis set file for basis set for elements with Z =   1 :
     "/home/arnfinn/jobb/dalton/svn/qmmm_devel/basis/cc-pVDZ"


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 


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

  Number of atom types :    3
  Total number of atoms:    4

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C           1    6.0000    26    14      [9s4p1d|3s2p1d]                                    
  O           1    8.0000    26    14      [9s4p1d|3s2p1d]                                    
  H           2    1.0000     7     5      [4s1p|2s1p]                                        
  ----------------------------------------------------------------------
  total:      4   16.0000    66    38
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00D-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   12
  C       :     1  x  -3.0015786160    2  y  -1.4563174382    3  z   0.0550080378
  O       :     4  x  -3.1314330364    5  y   0.8240509816    6  z  -0.0184248297
  H       :     7  x  -1.1728925345    8  y  -2.4468589606    9  z   0.1025195320
  H       :    10  x  -4.7395143797   11  y  -2.6116033945   12  z   0.0761219478


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.209298    0.000000
 H     :    1.100831    2.018474    0.000000
 H     :    1.104391    2.007988    1.889439    0.000000


  Max interatomic separation is    2.0185 Angstrom (    3.8144 Bohr)
  between atoms    3 and    2, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.209298
  bond distance:  H          C            1.100831
  bond distance:  H          C            1.104391


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H            121.724
  bond angle:     O          C          H            120.357
  bond angle:     H          C          H            117.919




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


@  Nuclear repulsion energy :   31.249215315972


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************


  *********************************** 
  QMMM electrostatic potential: 
  Multipole order                          2
  Anisotropic polarization
  *********************************** 

  ---------------- 
  QMMM information 
  ---------------- 

  MM coordinates in au 
  -------------------- 
   1      1.843496      2.848231     -0.156072
   2      0.043842      2.456825     -0.166415
   3      2.104836      3.855856      1.331837
   4      3.087255     -2.330561     -0.148884
   5      3.024958     -0.488757     -0.164317
   6      4.353698     -2.761321      1.077267

  MM charges 
  ---------- 
   1     -0.739083
   2      0.367781
   3      0.371302
   4     -0.739251
   5      0.367701
   6      0.371549

  MM dipoles (x,y,z) 
  ------------------ 
   1     -0.130726      0.053697      0.127642
   2      0.206730      0.034515     -0.015142
   3     -0.012978     -0.115434     -0.176127
   4      0.104509      0.119958      0.105009
   5     -0.007017     -0.209358     -0.012012
   6     -0.149391      0.032310     -0.145096

  MM quadrupoles (xx,xy,xz,yy,yz,zz) 
  ---------------------------------- 
   1     -3.891547      0.304260      0.173633     -4.459360      0.423169     -4.191583
   2     -0.005325      0.100015     -0.012868     -0.496247     -0.009445     -0.522055
   3     -0.512511      0.030179      0.049587     -0.348072      0.235727     -0.153360
   4     -4.353401     -0.233572      0.425650     -3.805243     -0.203696     -4.384519
   5     -0.519714     -0.002960     -0.007650      0.013940      0.010074     -0.519177
   6     -0.250836     -0.075418      0.246439     -0.497307     -0.073459     -0.266235

  MM polarizabilities in au (xx,xy,xz,yy,yz,zz) 
  --------------------------------------------- 
   1      5.476036     -0.086680     -0.072168      5.582289     -0.124016      5.513301
   2      3.488465      0.298404     -0.071204      1.739307     -0.350308      1.360585
   3      1.185490      0.051653      0.313695      2.359411      0.746038      2.948993
   4      5.555285      0.085929     -0.120805      5.451407      0.074565      5.566247
   5      1.509477     -0.036783     -0.366362      3.554625      0.100598      1.528430
   6      2.658086     -0.248537      0.746655      1.223737     -0.356653      2.608086


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00D-15

 Number of two-electron integrals written:      265159 ( 96.5% )
 Megabytes written:                              3.037

 >>> Time used in TWOINT is   0.20 seconds
 >>>> Total CPU  time used in HERMIT:   0.23 seconds
 >>>> Total wall time used in HERMIT:   0.00 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          T
     Using EHT  model:          F
     Number of Huckel orbitals each symmetry:   12

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
        which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -20.684762     -11.351957      -1.632119      -1.046130      -0.813248
           -0.702067      -0.605910      -0.491826      -0.321033      -0.161573
           -0.131670      -0.108806

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Wed Apr  6 09:17:48 2011 
     Host name              : stanley                                 

 Title lines from ".mol" input file:
     CH2O PLUS 2 WATERS                                                      
     ------------------------                                                

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham DFT calculation.

 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     For the wave function of type :      >>> KS-DFT   <<<
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1
 
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
     Secondary orbitals                 30 |   30
     Total number of orbitals           38 |   38
     Number of basis functions          38 |   38

     Optimization information
     ========================
     Number of configurations                 1
     Number of orbital rotations            240
     ------------------------------------------
     Total number of variables              241

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00D-10
 
     This is a DFT calculation of type: B3LYP
 Weighted mixed functional:
               HF exchange:    0.20000
                       VWN:    0.19000
                       LYP:    0.81000
                     Becke:    0.72000
                    Slater:    0.80000


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------------
         DFT grid generation - Radial Quadrature  : LMG scheme
         DFT grid generation -  partitioning : Original Becke partitioning
         DFT grid generation - Radial integration threshold: 1e-13
         DFT grid generation - Angular polynomials in range [15 35]
         DFT grid generation - Atom:    1*1 points= 21894 compressed from 21894 (117 radial)
         DFT grid generation - Atom:    2*1 points= 21522 compressed from 21522 (117 radial)
         DFT grid generation - Atom:    3*1 points= 20634 compressed from 20634 ( 87 radial)
         DFT grid generation - Atom:    4*1 points= 20634 compressed from 20634 ( 87 radial)
         DFT grid generation - Number of grid points:    84684; grid generation time:      0.1 s
      K-S energy, electrons, error :    -11.716637424232  15.9999999914   -8.63D-09

 Done with induced dipoles in            6  iterations

 Acc. iterations:           6
   1  -114.001203848     -3.382252482250E-02    2.71747D+00   -1.14D+02
      Virial theorem: -V/T =      2.006749
      MULPOP  C       1.08; O      -0.74; H      -0.17; H      -0.17; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.859865843002  15.9999999671   -3.29D-08

 Done with induced dipoles in            7  iterations

 Acc. iterations:          13
   2  -113.462768295      8.621707698802E-03    4.34174D+00    5.38D-01
      Virial theorem: -V/T =      1.998566
      MULPOP  C      -1.67; O       1.17; H       0.28; H       0.21; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -12.009761936348  15.9999999416   -5.84D-08

 Done with induced dipoles in            7  iterations

 Acc. iterations:          20
   3  -114.421416388     -3.502101906112E-02    1.02590D+00   -9.59D-01
      Virial theorem: -V/T =      1.994010
      MULPOP  C       0.36; O      -0.54; H       0.13; H       0.05; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.806144826139  15.9999999442   -5.58D-08

 Done with induced dipoles in            6  iterations

 Acc. iterations:          26
   4  -114.467190768     -2.472441403510E-02    2.46225D-01   -4.58D-02
      Virial theorem: -V/T =      2.009365
      MULPOP  C       0.21; O      -0.26; H       0.07; H      -0.02; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.858469437566  15.9999999431   -5.69D-08

 Done with induced dipoles in            5  iterations

 Acc. iterations:          31
   5  -114.470391098     -2.661344167658E-02    2.53548D-02   -3.20D-03
      Virial theorem: -V/T =      2.005373
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.02; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853504962967  15.9999999431   -5.69D-08

 Done with induced dipoles in            3  iterations

 Acc. iterations:          34
   6  -114.470427534     -2.630591381529E-02    2.52684D-03   -3.64D-05
      Virial theorem: -V/T =      2.005734
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853968549426  15.9999999431   -5.69D-08

 Done with induced dipoles in            2  iterations

 Acc. iterations:          36
   7  -114.470427781     -2.633238260099E-02    1.23203D-04   -2.48D-07
      Virial theorem: -V/T =      2.005700
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853979923873  15.9999999431   -5.69D-08

 Done with induced dipoles in            1  iterations

 Acc. iterations:          37
   8  -114.470427798     -2.633446409214E-02    2.71748D-05   -1.71D-08
      Virial theorem: -V/T =      2.005699
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853977555492  15.9999999431   -5.69D-08

 Done with induced dipoles in            0  iterations

 Acc. iterations:          37
   9  -114.470427815     -2.633403548600E-02    2.51279D-06   -1.71D-08
      Virial theorem: -V/T =      2.005700
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853977610106  15.9999999431   -5.69D-08

 Done with induced dipoles in            0  iterations

 Acc. iterations:          37
  10  -114.470427810     -2.633404392654E-02    7.64732D-07    5.53D-09
      Virial theorem: -V/T =      2.005700
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853977647570  15.9999999431   -5.69D-08

 Done with induced dipoles in            0  iterations

 Acc. iterations:          37
  11  -114.470427809     -2.633404960966E-02    5.79058D-08    1.06D-09
      Virial theorem: -V/T =      2.005700
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853977651043  15.9999999431   -5.69D-08

 Done with induced dipoles in            0  iterations

 Acc. iterations:          37
  12  -114.470427809     -2.633404968259E-02    6.88791D-09    1.08D-10
      Virial theorem: -V/T =      2.005700
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853977652089  15.9999999431   -5.69D-08

 Done with induced dipoles in            0  iterations

 Acc. iterations:          37
  13  -114.470427809     -2.633404972022E-02    4.37177D-10    2.33D-11
      Virial theorem: -V/T =      2.005700
      MULPOP  C       0.16; O      -0.28; H       0.10; H       0.01; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.853977652111  15.9999999431   -5.69D-08

 Done with induced dipoles in            0  iterations

 Acc. iterations:          37
  14  -114.470427809     -2.633404971692E-02    7.61731D-11    4.89D-12

 *** DIIS converged in  14 iterations !
   - total time used in SIRFCK :              0.08 seconds
   - QM/MM times:
     - total time used in QMMMFCK      :       1.29 seconds
     - total time used in QMMM MULPOLES:       0.78 seconds
     - total time used in QMMM_POLARI  :       0.51 seconds
     - MMITER times:
       - total time used in GET_IND_DIPOLES_2:       0.28 seconds
       - total time used in MMPOLARI_ITER2   :       0.28 seconds
       - total time used in F2QMMM           :       0.00 seconds
       - total time used in the iteration    :       0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Sym       Kohn-Sham orbital energies

  1    -19.16731614   -10.28799608    -1.06931232    -0.63545382    -0.49840687
        -0.46174420    -0.40603385    -0.27744815    -0.04832902     0.08479355
         0.16587419     0.20421444     0.45222493

    E(LUMO) :    -0.04832902 au (symmetry 1)
  - E(HOMO) :    -0.27744815 au (symmetry 1)
  ------------------------------------------
    gap     :     0.22911913 au

 >>> Writing SIRIFC interface file <<<


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     QM/MM "QMMM" calculation converged :

     Charge contribution:         -0.014924453275
     Dipole contribution:          0.005065816441
     Quadrupole contribution:     -0.007670309597
     Electronic Pol. energy:      -0.062585214540
     Nuclear pol. energy:          0.059153961354
     Multipole Pol. energy:       -0.005373850100
     Total QM/MM energy:          -0.026334049717

     Final DFT energy:           -114.470427808760                 
     Nuclear repulsion:            31.249215315972
     Electronic energy:          -145.693309075016

     Final gradient norm:           0.000000000076
  -------------------------------------- 
      Output from MM property module     
  ---------------------------------------


  MM total charge:  -1.66533453693773481E-016

  MM total charge dipole moment (x,y,z): 
 -0.11720324286847639       0.74736724825390444       0.99856024493110529     

  MM total permanent dipole moment (x,y,z): 
  1.11274526999999568E-002 -8.43115480999999944E-002 -0.11572617469999998     

  MM total induced dipole moment (x,y,z): 
 -0.23469464574702148       0.40359089589126856       4.44990090830958869E-003

  MM total dipole moment (x,y,z): 
 -0.34077043591549794        1.0666465960451730       0.88728397113941482     

 MM properties skipped since MMITER
  ---------------------------------------


 
     Date and time (Linux)  : Wed Apr  6 09:18:21 2011 
     Host name              : stanley                                 

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

    Orbital         4        5        6        7        8        9       10
   1 C   :1s     0.0081  -0.0022   0.0101   0.0003  -0.0030  -0.0004   0.0229
   2 C   :1s    -0.7005   0.0202   0.1182   0.0018  -0.0036  -0.0027  -0.3898
   3 C   :1s     0.1703   0.0238  -0.1244  -0.0022   0.0386   0.0024  -1.4145
   4 C   :2px    0.0067   0.6488   0.0229   0.0035   0.2689  -0.0075   0.1813
   5 C   :2py    0.2789  -0.0191   0.5661  -0.0176  -0.0135   0.0257   0.3285
   6 C   :2pz   -0.0095   0.0061  -0.0180  -0.5379   0.0010   0.7685  -0.0117
   7 C   :2px   -0.0144  -0.1704   0.0057  -0.0002  -0.1189  -0.0025   0.2835
   8 C   :2py   -0.0432   0.0112  -0.2390   0.0013  -0.0172   0.0066   0.3627
   9 C   :2pz    0.0007  -0.0016   0.0073   0.0505  -0.0004   0.2012  -0.0066
  10 C   :3d2-  -0.0031   0.0097  -0.0000   0.0004  -0.0547   0.0001  -0.0056
  11 C   :3d1-  -0.0007   0.0001  -0.0006  -0.0411  -0.0004  -0.0332   0.0002
  14 C   :3d2+  -0.0125  -0.0002  -0.0184   0.0014   0.0066   0.0010   0.0030
  16 O   :1s     0.3863  -0.0075  -0.2890  -0.0001   0.0357   0.0002  -0.0019
  17 O   :1s     0.0068   0.0111  -0.1285   0.0008   0.0042  -0.0016  -0.0426
  18 O   :2px    0.0104   0.4372   0.0424   0.0074  -0.8224   0.0046  -0.0584
  19 O   :2py    0.1916   0.0648  -0.7370  -0.0236  -0.0097  -0.0199  -0.0982
  20 O   :2pz   -0.0063   0.0010   0.0236  -0.7410  -0.0082  -0.6235   0.0033
  21 O   :2px    0.0059  -0.0193  -0.0202   0.0004  -0.0475   0.0019  -0.0589
  22 O   :2py    0.0099   0.0002   0.0432   0.0006  -0.0146  -0.0045  -0.0417
  23 O   :2pz   -0.0000  -0.0005  -0.0014   0.0232   0.0007  -0.1616   0.0005
  24 O   :3d2-   0.0009  -0.0183  -0.0053  -0.0003   0.0089   0.0001  -0.0029
  25 O   :3d1-   0.0005  -0.0000  -0.0012   0.0242   0.0002  -0.0065   0.0001
  26 O   :3d0    0.0051  -0.0003  -0.0110  -0.0015   0.0022   0.0003  -0.0017
  28 O   :3d2+   0.0068   0.0046  -0.0186  -0.0007  -0.0031   0.0002   0.0040
  29 H   :1s    -0.3761   0.4432  -0.1550  -0.0002   0.4246  -0.0018   0.1855
  30 H   :1s     0.1738  -0.1595   0.0270  -0.0002  -0.0606   0.0037   0.5825
  31 H   :2px    0.0256  -0.0132   0.0083   0.0001  -0.0060  -0.0002   0.0036
  33 H   :2pz    0.0004  -0.0002  -0.0002  -0.0092  -0.0001   0.0246   0.0001
  34 H   :1s    -0.4091  -0.3997  -0.2402   0.0021  -0.4285  -0.0018   0.1501
  35 H   :1s     0.1640   0.1275   0.0481  -0.0007   0.0091   0.0011   1.6435
  36 H   :2px   -0.0229  -0.0068  -0.0116   0.0002  -0.0038  -0.0002  -0.0048
  37 H   :2py   -0.0106  -0.0106   0.0047  -0.0003  -0.0055   0.0009   0.0035
  38 H   :2pz    0.0001   0.0002  -0.0003  -0.0099   0.0000   0.0266  -0.0002



 >>>> Total CPU  time used in SIRIUS :     32.53 seconds
 >>>> Total wall time used in SIRIUS :     33.00 seconds

 
     Date and time (Linux)  : Wed Apr  6 09:18:21 2011 
     Host name              : stanley                                 


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



        *****************************************************************
        ******** Output from **PROPE input processing for ABACUS ********
        *****************************************************************

 QMMM calculation


 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

 Default print level:        0

      Electronic excitation energies 
      Natural orbital connection is used
      for perturbation dependent basis sets.


 Changes of defaults for .EXCITA:
 --------------------------------

 Number of excitation energies:    2    0    0    0    0    0    0    0
 Print level          :    0
 Integral print level :    0
 Threshold            : 1.00E-07
 Maximum iterations   :   60
 Dipole strength

 Center of mass dipole origin  :   -3.067740   -0.312997    0.018175


                 .------------------------------------------------.
                 | Starting in Static Property Section (ABACUS) - |
                 `------------------------------------------------'


 
     Date and time (Linux)  : Wed Apr  6 09:18:21 2011 
     Host name              : stanley                                 
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            3  iterations

 Acc. iterations:          40

 Done with induced dipoles in            6  iterations

 Acc. iterations:          46
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            3  iterations

 Acc. iterations:          49

 Done with induced dipoles in            5  iterations

 Acc. iterations:          54
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            4  iterations

 Acc. iterations:          58

 Done with induced dipoles in            6  iterations

 Acc. iterations:          64
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            4  iterations

 Acc. iterations:          68

 Done with induced dipoles in            5  iterations

 Acc. iterations:          73
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            3  iterations

 Acc. iterations:          76

 Done with induced dipoles in            4  iterations

 Acc. iterations:          80
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            3  iterations

 Acc. iterations:          83

 Done with induced dipoles in            4  iterations

 Acc. iterations:          87
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:          92

 Done with induced dipoles in            5  iterations

 Acc. iterations:          97
       Electrons: 16.000000(-5.69e-08): LR-DFT*2 evaluation time:       2.4 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         102

 Done with induced dipoles in            5  iterations

 Acc. iterations:         107
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         112
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            3  iterations

 Acc. iterations:         115
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         120
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         125
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         130
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         135
       Electrons: 16.000000(-5.69e-08): LR-DFT*1 evaluation time:       1.8 s

 Done with induced dipoles in            5  iterations

 Acc. iterations:         140
 >>> Time used in EXCITA is  33.58 seconds


   ***************************************************************************
   ************************ FINAL RESULTS from ABACUS ************************
   ***************************************************************************


 
     Date and time (Linux)  : Wed Apr  6 09:18:55 2011 
     Host name              : stanley                                 


                             Molecular geometry (au)
                             -----------------------

 C         -3.0015786160           -1.4563174382            0.0550080378
 O         -3.1314330364            0.8240509816           -0.0184248297
 H         -1.1728925345           -2.4468589606            0.1025195320
 H         -4.7395143797           -2.6116033945            0.0761219478





                        Molecular wave function and energy
                        ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0

     Total energy       -114.4704278088 au (Hartrees)
                         -3114.89879178 eV
                           -300542.0619 kJ/mol


                             Relativistic corrections
                             ------------------------

     Darwin correction:                          0.2596110361 au
     Mass-velocity correction:                  -0.3259727291 au

     Total relativistic correction:             -0.0663616930 au (0.0580%)
     Non-relativistic + relativistic energy:  -114.5367895018 au




                                  Dipole moment
                                  -------------

                 au               Debye          C m (/(10**-30)
              1.116817           2.838665           9.468769


                             Dipole moment components
                             ------------------------

                 au               Debye          C m (/(10**-30)

      x      0.12613744         0.32060937         1.06943776
      y     -1.10925318        -2.81944018        -9.40464014
      z      0.03044428         0.07738165         0.25811740


   Units:   1 a.u. =   2.54175 Debye 
            1 a.u. =   8.47835 (10**-30) C m (SI)




                Singlet electronic excitation energies
                --------------------------------------

                =======================================
                 Sym.   Mode   Frequency    Frequency
                ex. st.  No.      (au)          (eV)
                =======================================
                   1        1    0.154260    4.197616
                   1        2    0.321576    8.750517
                ---------------------------------------


                Electric transition dipole moments (in a.u.)
                --------------------------------------------

  Sym.   Mode    Frequency       Velocity/Frequency              Length
 ex. st.  No.      (au)          x       y       z         x       y       z
 ==============================================================================
   1        1     0.154260     0.0011 -0.0012 -0.0239   -0.0000  0.0001 -0.0081
   1        2     0.321576    -0.6462 -0.1780  0.0002   -0.7111 -0.1485 -0.0014
 ------------------------------------------------------------------------------


                               Oscillator strengths
                               --------------------

  Oscillator strengths are dimensionless.

  Sym.   Mode        Frequency     Oscillator-strength 
 ex. st.  No.           (eV)        velocity   length   
 -------------------------------------------------------
   1        1         4.197616        0.0001   0.0000
   1        2         8.750517        0.0963   0.1131


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.209298    0.000000
 H     :    1.100831    2.018474    0.000000
 H     :    1.104391    2.007988    1.889439    0.000000


  Max interatomic separation is    2.0185 Angstrom (    3.8144 Bohr)
  between atoms    3 and    2, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.209298
  bond distance:  H          C            1.100831
  bond distance:  H          C            1.104391


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H            121.724
  bond angle:     O          C          H            120.357
  bond angle:     H          C          H            117.919




 CPU time statistics for ABACUS
 ------------------------------

 EXCITA     00:00:34     100 %

 TOTAL      00:00:34     100 %


   - QM/MM times:
     - total time used in QMMMFIRST:       0.00 seconds
     - total time used in QMMMB2   :       0.00 seconds
 >>>> Total CPU  time used in ABACUS:  33.60 seconds
 >>>> Total wall time used in ABACUS:  34.00 seconds


                   .-------------------------------------------.
                   | End of Static Property Section (ABACUS) - |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  1 minute   6 seconds
 >>>> Total wall time used in DALTON:  1 minute   7 seconds

 
     Date and time (Linux)  : Wed Apr  6 09:18:55 2011 
     Host name              : stanley                                 
END REFOUT

