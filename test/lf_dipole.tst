########## Test description ########################
START DESCRIPTION
PCM test
HF wavefunction with STO-3G basis
Calculation of dipole moment with and without cavity field correction
KEYWORDS: pcm hf dipole localfield essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enehf
dipole
dipolelf
OVERRIDE 7 1.0e-5
nuc
tes
sym
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN WAVE FUNCTION
.RUN PROPERTIES
*PCM
.SOLVNT
WATER
.NPCMMT
0
.NESFP
4
.ICESPH
2
.LOCFLD
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
**WAVE FUNCTIONS
.HF
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
STO-3G
Opt 2.2 paracyclophane chormoph.

    3    0         1
        6.    1
C     0.000000     0.000000     0.000000
        8.    1
O     0.000000     0.000000     1.220000
        1.    2
H     0.943102     0.000000    -0.544500
H    -0.943102     0.000000    -0.544500
END MOLINP

########## Reference Output ########################
START REFOUT


 **************************************************************************************
 ******************** DALTON2011 - An electronic structure program ********************
 **************************************************************************************

    This is output from DALTON Release 2011 (Rev. 0, Dec. 2010)

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
    Release DALTON2011 (2010), see http://daltonprogram.org"

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


     Date and time (Linux)  : Tue Feb 15 12:47:28 2011 
     Host name              : stanley                                 

 * Work memory size             :   100000000 =  762.94 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/jobb/dalton/svn/pure_trunk/test/perl-pid.6105__2011_2_15__12.26
   2) /home/arnfinn/jobb/dalton/svn/pure_trunk/basis/


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

     NON-EQ = F     NEQRSP =F
 POLYG          60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.7000
     2    0.0000    0.0000    0.0000    1.5000
     3    0.0000    0.0000    0.0000    1.2000
     4    0.0000    0.0000    0.0000    1.2000


   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************


    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1: Opt 2.2 paracyclophane chormoph.                                        
 2:                                                                         
    ------------------------------------------------------------------------

  Coordinates are entered in Angstrom and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  Used basis set file for basis set for elements with Z =   6 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  Used basis set file for basis set for elements with Z =   8 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  Used basis set file for basis set for elements with Z =   1 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/STO-3G"


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 
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
                           H           1.007825
                           H           1.007825

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

  Threshold for integrals:  1.00e-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   12
  C       :     1  x   0.0000000000    2  y   0.0000000000    3  z   0.0000000000
  O       :     4  x   0.0000000000    5  y   0.0000000000    6  z   2.3054658725
  H       :     7  x   1.7822044879    8  y   0.0000000000    9  z  -1.0289558751
  H       :    10  x  -1.7822044879   11  y   0.0000000000   12  z  -1.0289558751


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H     :    1.089000    2.000725    0.000000
 H     :    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstrom (    3.7808 Bohr)
  between atoms    3 and    2, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.220000
  bond distance:  H          C            1.089000
  bond distance:  H          C            1.089000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H            120.000
  bond angle:     O          C          H            120.000
  bond angle:     H          C          H            120.000




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


@  Nuclear repulsion energy :   31.163673729192


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************


 Default print level:        1

 * Nuclear model: Point charge

 Calculation of one- and two-electron Hamiltonian integrals.

 Center of mass  (bohr):      0.000000000000      0.000000000000      1.159648802225
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:        1505 ( 48.8% )
 Megabytes written:                              0.021


 MEMORY USED TO GENERATE CAVITY =    432042


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
 >>> Time used in Q-MAT is   0.30 seconds
 >>>> Total CPU  time used in HERMIT:   0.48 seconds
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
          -20.684619     -11.351900      -1.631288      -1.046947      -0.814716
           -0.700208      -0.605487      -0.493669      -0.322892      -0.164075
           -0.130192      -0.105106

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Tue Feb 15 12:47:28 2011 
     Host name              : stanley                                 

 Title lines from ".mol" input file:
     Opt 2.2 paracyclophane chormoph.                                        
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.

 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     For the wave function of type :      >>> HF       <<<
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1

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
     Number of configurations                 1
     Number of orbital rotations             32
     ------------------------------------------
     Total number of variables               33

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.97018
 NUCLEAR APPARENT CHARGE -15.78962
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00e-15
   1  -112.142757101     -2.184941165320e-02    1.57297e+00   -1.12e+02

 Virial theorem: -V/T =      2.002200
      MULPOP C     4.90; O     8.84; H     2.26; 
   2  -112.306432509     -4.221971510745e-04    8.59832e-01   -1.64e-01

 Virial theorem: -V/T =      2.010790
      MULPOP C     6.41; O     7.82; H     1.77; 
   3  -112.355348559     -3.080924397281e-03    4.78592e-02   -4.89e-02

 Virial theorem: -V/T =      2.008440
      MULPOP C     5.92; O     8.21; H     1.86; 
   4  -112.355567191     -2.888297151717e-03    1.11568e-02   -2.19e-04

 Virial theorem: -V/T =      2.008333
      MULPOP C     5.93; O     8.21; H     1.86; 
   5  -112.355580051     -2.869260094137e-03    3.15920e-03   -1.29e-05

 Virial theorem: -V/T =      2.008320
      MULPOP C     5.93; O     8.21; H     1.86; 
   6  -112.355581951     -2.857409074821e-03    6.54862e-04   -1.90e-06

 Virial theorem: -V/T =      2.008327
      MULPOP C     5.93; O     8.21; H     1.86; 
   7  -112.355582031     -2.856183353160e-03    7.06045e-05   -8.02e-08

 Virial theorem: -V/T =      2.008329
      MULPOP C     5.93; O     8.21; H     1.86; 
   8  -112.355582031     -2.856043623204e-03    4.54100e-06   -6.03e-10
 DIIS converged in   8 iterations !
 - total time used in SIRFCK        :       0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Sym       Hartree-Fock orbital energies

  1    -20.29958353   -11.12079213    -1.33320845    -0.80284544    -0.63846097
        -0.54122299    -0.43900852    -0.35348957     0.28385000     0.63726330
         0.77009303     0.90607780

    E(LUMO) :     0.28385000 au (symmetry 1)
  - E(HOMO) :    -0.35348957 au (symmetry 1)
  ------------------------------------------
    gap     :     0.63733957 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -112.355582031491                 
     Nuclear repulsion:            31.163673729192
     Electronic energy:          -143.516399717060

     Final gradient norm:           0.000004541002

 
     Date and time (Linux)  : Tue Feb 15 12:47:31 2011 
     Host name              : stanley                                 

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s    -0.0005  -0.9926   0.1223  -0.1873  -0.0000   0.0268   0.0000
   2 C   :1s     0.0072  -0.0332  -0.2782   0.5808   0.0000  -0.0860  -0.0000
   3 C   :2px    0.0000   0.0000   0.0000  -0.0000   0.5352   0.0000  -0.0000
   4 C   :2py    0.0000  -0.0000   0.0000   0.0000   0.0000   0.0000  -0.6051
   5 C   :2pz    0.0062  -0.0008  -0.1567  -0.2136  -0.0000  -0.4552   0.0000
   6 O   :1s    -0.9943  -0.0002   0.2195   0.1022   0.0000  -0.0895   0.0000
   7 O   :1s    -0.0259   0.0058  -0.7701  -0.4468  -0.0000   0.4791  -0.0000
   8 O   :2px   -0.0000  -0.0000  -0.0000  -0.0000   0.4297   0.0000  -0.0000
   9 O   :2py   -0.0000   0.0000   0.0000  -0.0000  -0.0000  -0.0000  -0.6805
  10 O   :2pz    0.0056  -0.0018   0.1657  -0.1775  -0.0000   0.6779   0.0000
  11 H   :1s    -0.0003   0.0066  -0.0324   0.2612   0.2994   0.1639   0.0000
  12 H   :1s    -0.0003   0.0066  -0.0324   0.2612  -0.2994   0.1639  -0.0000

 Orbital           8        9       10
   1 C   :1s    -0.0000  -0.0000  -0.2009
   2 C   :1s     0.0000   0.0000   1.2706
   3 C   :2px    0.1836   0.0000   0.0000
   4 C   :2py    0.0000  -0.8238   0.0000
   5 C   :2pz    0.0000   0.0000  -0.5035
   6 O   :1s    -0.0000   0.0000   0.0197
   7 O   :1s     0.0000  -0.0000  -0.1040
   8 O   :2px   -0.8807  -0.0000  -0.0000
   9 O   :2py    0.0000   0.7628  -0.0000
  10 O   :2pz    0.0000   0.0000   0.1752
  11 H   :1s     0.3445   0.0000  -0.9049
  12 H   :1s    -0.3445  -0.0000  -0.9049



 >>>> Total CPU  time used in SIRIUS :      2.82 seconds
 >>>> Total wall time used in SIRIUS :      3.00 seconds

 
     Date and time (Linux)  : Tue Feb 15 12:47:31 2011 
     Host name              : stanley                                 


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



        *****************************************************************
        ******** Output from **PROPE input processing for ABACUS ********
        *****************************************************************



 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

 Default print level:        0

      Natural orbital connection is used
      for perturbation dependent basis sets.

 Center of mass dipole origin  :    0.000000    0.000000    1.159649


                 .------------------------------------------------.
                 | Starting in Static Property Section (ABACUS) - |
                 `------------------------------------------------'


 
     Date and time (Linux)  : Tue Feb 15 12:47:31 2011 
     Host name              : stanley                                 
 ----ADDING LOCAL FIELD CONTRIBUTION----


   ***************************************************************************
   ************************ FINAL RESULTS from ABACUS ************************
   ***************************************************************************


 
     Date and time (Linux)  : Tue Feb 15 12:47:32 2011 
     Host name              : stanley                                 


                             Molecular geometry (au)
                             -----------------------

 C          0.0000000000            0.0000000000            0.0000000000
 O          0.0000000000            0.0000000000            2.3054658725
 H          1.7822044879            0.0000000000           -1.0289558751
 H         -1.7822044879            0.0000000000           -1.0289558751





                        Molecular wave function and energy
                        ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0

     Total energy       -112.3555820315 au (Hartrees)
                         -3057.35091079 eV
                           -294989.5351 kJ/mol


                             Relativistic corrections
                             ------------------------

     Darwin correction:                          0.1684676504 au
     Mass-velocity correction:                  -0.2161983088 au

     Total relativistic correction:             -0.0477306584 au (0.0425%)
     Non-relativistic + relativistic energy:  -112.4033126899 au




                                  Dipole moment
                                  -------------

                 au               Debye          C m (/(10**-30)
              0.675590           1.717179           5.727893


                             Dipole moment components
                             ------------------------

                 au               Debye          C m (/(10**-30)

      x      0.00000000         0.00000000         0.00000000
      y     -0.00000000        -0.00000000        -0.00000000
      z     -0.67559026        -1.71717904        -5.72789274


   Units:   1 a.u. =   2.54175 Debye 
            1 a.u. =   8.47835 (10**-30) C m (SI)




                       Local-field corrected dipole moment
                       -----------------------------------

                 au               Debye          C m (/(10**-30)
              0.900539           2.288943           7.635091


                  Local-field corrected dipole moment components
                  ----------------------------------------------

                 au               Debye          C m (/(10**-30)

      x      0.00000000         0.00000000         0.00000000
      y      0.00000000         0.00000000         0.00000000
      z     -0.90053943        -2.28894277        -7.63509125




   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H     :    1.089000    2.000725    0.000000
 H     :    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstrom (    3.7808 Bohr)
  between atoms    3 and    2, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.220000
  bond distance:  H          C            1.089000
  bond distance:  H          C            1.089000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H            120.000
  bond angle:     O          C          H            120.000
  bond angle:     H          C          H            120.000




 CPU time statistics for ABACUS
 ------------------------------




 >>>> Total CPU  time used in ABACUS:   0.44 seconds
 >>>> Total wall time used in ABACUS:   1.00 seconds


                   .-------------------------------------------.
                   | End of Static Property Section (ABACUS) - |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   3.74 seconds
 >>>> Total wall time used in DALTON:   4.00 seconds

 
     Date and time (Linux)  : Tue Feb 15 12:47:32 2011 
     Host name              : stanley                                 
END REFOUT

