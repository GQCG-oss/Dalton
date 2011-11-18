########## Test description ########################
START DESCRIPTION
KEYWORDS: shielding spin-spin coupling pcm dft long
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
gauge_or
dipole
nuc
tes
sym
cmass
shield
Icoupl
Acoupl
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN PROPERTIES
*PCM
.SOLVNT
H2O
.NPCMMT
0
.NESFP
6
.ICESPH
2
*PCMCAV
.INA
1
2
3
4
5
6
.RIN
1.7
1.7
1.3
1.2
1.2
1.2
.AREATS
1.2
**WAVE FUNCTIONS
.DFT
B3LYP
**PROPERTIES
.SHIELD
.SPIN-SPIN
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
STO-3G
Vinyllithium
geometry
    3              
        6.    2   
C1         0.8778403573           -0.0248311095            0.0000000000
C2        -1.0081980992            1.6717602378            0.0000000000
        3.    1  
Li        0.3603461587           -3.6719284437            0.0000000000
        1.    3   
H1         2.7426633015            0.8912639005            0.0000000000
H2        -2.9987984292            1.1226518917            0.0000000000
H3        -0.7303536961            3.7214577549            0.0000000000
END MOLINP

########## Reference Output ########################
START REFOUT


         ****************************************************************
         *********** DALTON - An electronic structure program ***********
         ****************************************************************

    This is output from DALTON (Release 2.0 rev. 0, Mar. 2005)

    Celestino Angeli,         University of Ferrara,        Italy      
    Keld L. Bak,              UNI-C,                        Denmark    
    Vebjoern Bakken,          University of Oslo,           Norway     
    Ove Christiansen,         Aarhus University,            Denmark    
    Renzo Cimiraglia,         University of Ferrara,        Italy      
    Sonia Coriani,            University of Trieste,        Italy      
    Paal Dahle,               University of Oslo,           Norway     
    Erik K. Dalskov,          UNI-C,                        Denmark    
    Thomas Enevoldsen,        SDU - Odense University,      Denmark    
    Berta Fernandez,          U. of Santiago de Compostela, Spain      
    Christof Haettig,         Forschungszentrum Karlsruhe,  Germany    
    Kasper Hald,              Aarhus University,            Denmark    
    Asger Halkier,            Aarhus University,            Denmark    
    Hanne Heiberg,            University of Oslo,           Norway     
    Trygve Helgaker,          University of Oslo,           Norway     
    Hinne Hettema,            University of Auckland,       New Zealand
    Hans Joergen Aa. Jensen,  Univ. of Southern Denmark,    Denmark    
    Dan Jonsson,              KTH Stockholm,                Sweden     
    Poul Joergensen,          Aarhus University,            Denmark    
    Sheela Kirpekar,          SDU - Odense University,      Denmark    
    Wim Klopper,              University of Karlsruhe,      Germany    
    Rika Kobayashi,           ANU Supercomputer Facility,   Australia  
    Jacob Kongsted,           Univ. of Southern Denmark,    Denmark    
    Henrik Koch,              University of Trondheim,      Norway     
    Andrea Ligabue,           University of Modena,         Italy      
    Ola B. Lutnaes,           University of Oslo,           Norway     
    Kurt V. Mikkelsen,        University of Copenhagen,     Denmark    
    Patrick Norman,           University of Linkoeping,     Sweden     
    Jeppe Olsen,              Aarhus University,            Denmark    
    Anders Osted,             Copenhagen University,        Denmark    
    Martin J. Packer,         University of Sheffield,      UK         
    Thomas B. Pedersen,       University of Lund,           Sweden     
    Zilvinas Rinkevicius,     KTH Stockholm,                Sweden     
    Elias Rudberg,            KTH Stockholm,                Sweden     
    Torgeir A. Ruden,         University of Oslo,           Norway     
    Kenneth Ruud,             University of Tromsoe,        Norway     
    Pawel Salek,              KTH Stockholm,                Sweden     
    Alfredo Sanchez de Meras, University of Valencia,       Spain      
    Trond Saue,               University of Strasbourg,     France     
    Stephan P. A. Sauer,      University of Copenhagen,     Denmark    
    Bernd Schimmelpfennig,    Forschungszentrum Karlsruhe,  Germany     
    K. O. Sylvester-Hvid,     University of Copenhagen,     Denmark    
    Peter R. Taylor,          University of Warwick,        UK         
    Olav Vahtras,             KTH Stockholm,                Sweden     
    David J. Wilson,          University of Oslo,           Norway     
    Hans Agren,               KTH Stockholm,                Sweden     

 ---------------------------------------------------------------------

     NOTE:
      
     This is an experimental code for the evaluation of molecular
     properties using (MC)SCF and CC wave functions. The authors
     accept no responsibility for the performance of the code or
     for the correctness of the results.
      
     The code (in whole or part) is provided under a licence and
     is not to be reproduced for further distribution without
     the written permission of the authors or their representatives.
      
     See the home page "http://www.kjemi.uio.no/software/dalton"
     for further information.
      
     If results obtained with this code are published,
     an appropriate citation would be:
      
     "Dalton, a molecular electronic structure program, Release 2.0
     (2005), see http://www.kjemi.uio.no/software/dalton/dalton.html"

     Date and time (Linux)  : Thu Sep 24 00:44:20 2009
     Host name              : stallo-2.local                          

 * Work memory size             :    50000000 =  381.47 megabytes.
 + memory for in-core integrals :    30000000

 * Directories for basis set searches:
   1) /home/ruud/DaltonFix/dalton/test/2009-09-23T17_22-testjob-pid-2787/perl-pid.19196__2009_9_24__0.43
   2) /home/ruud/DaltonFix/dalton/basis/


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

 ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT = H2O      **
 Optical and physical constants:
 EPS= 78.390; EPSINF=  1.776; RSOLV=  1.385 A; VMOL=  18.070 ML/MOL;
 TCE=2.57000e-04 1/K; STEN= 71.810 DYN/CM;  DSTEN=  0.6500; CMF=  1.2770


     -----------------------------------
     Input for PCM solvation calculation 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=H2O          EPS   = 78.3900     EPSINF=  1.7760
     RSOLV =  1.3850

     ICESPH =       2     NESFP =       6
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG          60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.7000
     2    0.0000    0.0000    0.0000    1.7000
     3    0.0000    0.0000    0.0000    1.3000
     4    0.0000    0.0000    0.0000    1.2000
     5    0.0000    0.0000    0.0000    1.2000
     6    0.0000    0.0000    0.0000    1.2000


   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************

    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1: Vinyllithium                                                            
 2: geometry                                                                
    ------------------------------------------------------------------------

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centres:    2
  Used basis set file for basis set for elements with Z =   6 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   3.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   3 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    3
  Used basis set file for basis set for elements with Z =   1 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"


                      SYMADD: Requested addition of symmetry
                      --------------------------------------

 Symmetry test threshold:  5.00e-06

 Symmetry class found: C(s)   

 Symmetry Independent Centres             
 ----------------------------
       6 :      1.95729159     0.00000000     0.00000000  Isotope  1
       6 :     -0.46801060    -0.74395734     0.00000000  Isotope  1
       3 :     -3.32872203     1.57669975     0.00000000  Isotope  1
       1 :      3.57171997    -1.29308872     0.00000000  Isotope  1
       1 :      2.51093358     1.98934417     0.00000000  Isotope  1
       1 :     -0.64226615    -2.81432680     0.00000000  Isotope  1

 The following element was found:   Z        


                         SYMGRP: Point group information
                         -------------------------------

Full point group is: C(s)           
Represented as:      Cs 

   * The point group was generated by:

      Reflection in the xy-plane

   * Group multiplication table

        |  E   Oxy
   -----+----------
     E  |  E 
    Oxy | Oxy   E 

   * Character table

        |  E   Oxy
   -----+----------
    A'  |   1    1
    A"  |   1   -1

   * Direct product table

        | A'   A" 
   -----+----------
    A'  | A' 
    A"  | A"   A' 
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1   -4.6801059946e-01   -7.4395734073e-01    0.0000000000e+00    1.7000000000e+00
   2    1.9572915950e+00    0.0000000000e+00    0.0000000000e+00    1.7000000000e+00
   3   -3.3287220260e+00    1.5766997527e+00    0.0000000000e+00    1.3000000000e+00
   4   -6.4226614841e-01   -2.8143267967e+00    0.0000000000e+00    1.2000000000e+00
   5    2.5109335755e+00    1.9893441650e+00    0.0000000000e+00    1.2000000000e+00
   6    3.5717199730e+00   -1.2930887194e+00    0.0000000000e+00    1.2000000000e+00


                                 Isotopic Masses
                                 ---------------

                           C1         12.000000
                           C2         12.000000
                           Li          7.016005
                           H1          1.007825
                           H2          1.007825
                           H3          1.007825

                       Total mass:    34.039480 amu
                       Natural abundance:  90.435 %

 Center-of-mass coordinates (a.u.):    0.000000    0.000000    0.000000


  Atoms and basis sets
  --------------------

  Number of atom types:     3
  Total number of atoms:    6

  Basis set used is "STO-3G" from the basis set library.

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C2          2    6.0000    15     5      [6s3p|2s1p]                                        
  Li          1    3.0000    15     5      [6s3p|2s1p]                                        
  H3          3    1.0000     3     1      [3s|1s]                                            
  ----------------------------------------------------------------------
  total:      6   18.0000    54    18
  ----------------------------------------------------------------------

  Threshold for integrals:  1.00e-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   18
  C1      :    1  x  -0.4680105995   2  y  -0.7439573407   3  z   0.0000000000
  C2      :    4  x   1.9572915950   5  y   0.0000000000   6  z   0.0000000000
  Li      :    7  x  -3.3287220260   8  y   1.5766997527   9  z   0.0000000000
  H1      :   10  x  -0.6422661484  11  y  -2.8143267967  12  z   0.0000000000
  H2      :   13  x   2.5109335755  14  y   1.9893441650  15  z   0.0000000000
  H3      :   16  x   3.5717199730  17  y  -1.2930887194  18  z   0.0000000000


  Symmetry Coordinates
  --------------------

  Number of coordinates in each symmetry:   12   6

  Symmetry  A'  ( 1)

   1   C1    x    1
   2   C1    y    2
   3   C2    x    4
   4   C2    y    5
   5   Li    x    7
   6   Li    y    8
   7   H1    x   10
   8   H1    y   11
   9   H2    x   13
  10   H2    y   14
  11   H3    x   16
  12   H3    y   17

  Symmetry  A"  ( 2)

  13   C1    z    3
  14   C2    z    6
  15   Li    z    9
  16   H1    z   12
  17   H2    z   15
  18   H3    z   18


   Interatomic separations (in Angstrom):
   --------------------------------------

            C1          C2          Li          H1          H2          H3    
            ------      ------      ------      ------      ------      ------
 C1    :    0.000000
 C2    :    1.342439    0.000000
 Li    :    1.949292    2.919021    0.000000
 H1    :    1.099466    2.027387    2.724012    0.000000
 H2    :    2.139411    1.092723    3.097918    3.040717    0.000000
 H3    :    2.157393    1.094573    3.954756    2.370799    1.825442    0.000000


  Max interatomic separation is    3.9548 Angstrom (    7.4734 Bohr)
  between atoms    6 and    3, "H3    " and "Li    ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C2         C1           1.342439
  bond distance:  Li         C1           1.949292
  bond distance:  H1         C1           1.099466
  bond distance:  H2         C2           1.092723
  bond distance:  H3         C2           1.094573


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     C2         C1         Li           123.897
  bond angle:     C2         C1         H1           111.864
  bond angle:     Li         C1         H1           124.239
  bond angle:     C1         C2         H2           122.605
  bond angle:     C1         C2         H3           124.253
  bond angle:     H2         C2         H3           113.141




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       8.334683          0.967354   -0.253431    0.000000
   IB      43.108403          0.253431    0.967354    0.000000
   IC      51.443086          0.000000    0.000000    1.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

          60635.6583          11723.4454           9824.0414 MHz
            2.022588            0.391052            0.327695 cm-1


  Nuclear repulsion energy :   37.740380557905


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:              15   3


  Symmetry  A' ( 1)

    1     C1       1s         1
    2     C1       1s         2
    3     C1       2px        3
    4     C1       2py        4
    5     C2       1s         6
    6     C2       1s         7
    7     C2       2px        8
    8     C2       2py        9
    9     Li       1s        11
   10     Li       1s        12
   11     Li       2px       13
   12     Li       2py       14
   13     H1       1s        16
   14     H2       1s        17
   15     H3       1s        18


  Symmetry  A" ( 2)

   16     C1       2pz        5
   17     C2       2pz       10
   18     Li       2pz       15

  Symmetries of electric field:  A' (1)  A' (1)  A" (2)

  Symmetries of magnetic field:  A" (2)  A" (2)  A' (1)


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************



     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:        8906 ( 60.6% )
 Megabytes written:                              0.103


 MEMORY USED TO GENERATE CAVITY =    432042


 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING
 Tessera cut in pieces and removed.

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING
 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING

 Total number of spheres =    6
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1   -0.247660542   -0.393685269    0.000000000    2.040000000   19.523217119
   2    1.035754102    0.000000000    0.000000000    2.040000000   20.882040220
   3   -1.761483829    0.834353573    0.000000000    1.560000000   20.499501612
   4   -0.339872607   -1.489277598    0.000000000    1.440000000    9.390637698
   5    1.328728820    1.052715592    0.000000000    1.440000000    9.327821095
   6    1.890072804   -0.684273079    0.000000000    1.440000000    9.350796585

 Total number of tesserae =     126
 Surface area =   88.97401433 (A^2)    Cavity volume =   68.53843516 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....
 
  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   0.09 seconds
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
     Number of Huckel orbitals each symmetry:   15    3

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
        which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -11.354163     -11.345026      -2.488259      -1.491177      -0.982597
           -0.772686      -0.526055      -0.403913      -0.169796      -0.144802
           -0.125633      -0.110256      -0.103046      -0.087557      -0.057634

 Huckel EWMO eigenvalues for symmetry :  2
           -0.570857      -0.260513      -0.112430

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Sep 24 00:44:20 2009
     Host name              : stallo-2.local                          

 Title lines from ".mol" input file:
     Vinyllithium                                                            
     geometry                                                                

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham DFT calculation.


 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     Number of closed shell electrons         18
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                2
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
     Abelian symmetry species          All    1    2
                                       --  --
     Total number of orbitals           18   15    3
     Number of basis functions          18   15    3

      ** Automatic occupation of RKS orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals               9    8    1

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       17.84961
 NUCLEAR APPARENT CHARGE -17.74060
 THEORETICAL -17.77038 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  18 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Angular polynomials in range [15 35]
 Atom:    1*1 points=18676 compressed from 18676 ( 92 radial)
 Atom:    2*1 points=18676 compressed from 18676 ( 92 radial)
 Atom:    3*1 points=18760 compressed from 18760 ( 92 radial)
 Atom:    4*1 points=18150 compressed from 18150 ( 75 radial)
 Atom:    5*1 points=18150 compressed from 18150 ( 75 radial)
 Atom:    6*1 points=18150 compressed from 18150 ( 75 radial)
 Number of grid points:   110562 Grid generation time:       0.1 s
K-S electrons/energy :   17.99998540022957  -11.17140272275581 err:-.15e-04
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -67.91025    67.76686    67.76757   -33.83630    -0.02422
   1  -84.2634350466     -2.421832037074e-02   1.05e+00  -8.43e+01    8  1
 MULPOP C1    4.99; C2    6.07; Li    5.03; H1    0.10; H2    0.09; H3    1.12; 
K-S electrons/energy :   17.99999219999205  -11.24012148049906 err:-.78e-05
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.11382    67.85685    67.85774   -33.83630    -0.03592
   2  -84.2875066919     -3.591973512174e-02   9.28e-01  -2.41e-02    8  1
 MULPOP C1    3.88; C2    6.77; Li    5.68; H1    0.02; H2    0.04; H3    1.05; 
K-S electrons/energy :   17.99999017098086  -11.20589442917609 err:-.98e-05
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.10684    67.86133    67.86238   -33.83630    -0.02787
   3  -84.3429239336     -2.786675964165e-02   4.42e-01  -5.54e-02    8  1
 MULPOP C1    6.16; C2    4.94; Li    5.06; H1    0.04; H2    0.04; H3    1.10; 
K-S electrons/energy :   17.99998813541722  -11.17505508658203 err:-.12e-04
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.03871    67.82931    67.83024   -33.83630    -0.02588
   4  -84.3573512930     -2.588236501144e-02   9.93e-02  -1.44e-02    8  1
 MULPOP C1    3.98; C2    6.89; Li    5.35; H1    0.03; H2    0.04; H3    1.10; 
K-S electrons/energy :   17.99998935703200  -11.19021347711874 err:-.11e-04
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.06954    67.84353    67.84449   -33.83630    -0.02707
   5  -84.3583011584     -2.706519382767e-02   3.25e-03  -9.50e-04    8  1
 MULPOP C1    4.01; C2    6.90; Li    5.31; H1    0.03; H2    0.04; H3    1.09; 
K-S electrons/energy :   17.99998937564979  -11.19013369068265 err:-.11e-04
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.06952    67.84350    67.84447   -33.83630    -0.02708
   6  -84.3583023937     -2.708173739956e-02   2.72e-04  -1.24e-06    8  1
 MULPOP C1    4.01; C2    6.90; Li    5.31; H1    0.03; H2    0.04; H3    1.09; 
K-S electrons/energy :   17.99998937762286  -11.19011671839790 err:-.11e-04
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.06953    67.84350    67.84447   -33.83630    -0.02708
   7  -84.3583024052     -2.708298454817e-02   4.09e-05  -1.15e-08    8  1
 MULPOP C1    4.01; C2    6.90; Li    5.31; H1    0.03; H2    0.04; H3    1.09; 
K-S electrons/energy :   17.99998937793219  -11.19011421008662 err:-.11e-04
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -68.06953    67.84350    67.84447   -33.83630    -0.02708
   8  -84.3583024056     -2.708329210869e-02   3.87e-06  -3.29e-10    8  1
 DIIS converged in   8 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   18
 Orbital occupations :    8    1

 Sym       Kohn-Sham orbital energies

  1     -9.93453682    -9.90160096    -1.83481472    -0.65609189    -0.47317883
        -0.36150688    -0.29201527    -0.13766624     0.06105873     0.11470539
         0.31055858     0.44431161     0.48764224

  2     -0.16418340     0.09282930     0.20547550

    E(LUMO) :     0.06105873 au (symmetry 1)
  - E(HOMO) :    -0.13766624 au (symmetry 1)
  ------------------------------------------
    gap     :     0.19872497 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final DFT energy:            -84.358302405561                 
     Nuclear repulsion:            37.740380557905
     Electronic energy:          -122.071599671357

     Final gradient norm:           0.000003866884

 
     Date and time (Linux)  : Thu Sep 24 00:44:25 2009
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C1  :1s    -0.0023  -0.9906   0.0016   0.1579   0.1495   0.0389  -0.0135
   2 C1  :1s    -0.0093  -0.0445   0.0003  -0.4140  -0.4464  -0.1132   0.0641
   3 C1  :2px   -0.0066  -0.0041   0.0044  -0.1495   0.0139  -0.0527  -0.3601
   4 C1  :2py   -0.0015   0.0010  -0.0039   0.0092   0.1543  -0.2597  -0.2779
   5 C2  :1s     0.9906  -0.0033  -0.0003   0.1962  -0.1087   0.0027   0.0408
   6 C2  :1s     0.0441   0.0100   0.0019  -0.5220   0.3473   0.0006  -0.1297
   7 C2  :2px   -0.0014  -0.0060  -0.0013   0.0776   0.2711   0.3359   0.3068
   8 C2  :2py   -0.0004  -0.0024   0.0013   0.0369   0.1430  -0.3860   0.3549
   9 Li  :1s     0.0002  -0.0004  -0.9869   0.0177   0.0244   0.0289  -0.0325
  10 Li  :1s    -0.0009   0.0062  -0.0483  -0.0164  -0.0320  -0.0454   0.0700
  11 Li  :2px   -0.0012   0.0059  -0.0001  -0.0102  -0.0167  -0.0458   0.0258
  12 Li  :2py    0.0005  -0.0046  -0.0004   0.0096   0.0284  -0.0006  -0.0550
  13 H1  :1s     0.0003   0.0083   0.0026  -0.1191  -0.2777   0.1747   0.3723
  14 H2  :1s    -0.0089  -0.0009   0.0026  -0.1370   0.2684  -0.2397   0.3636
  15 H3  :1s    -0.0091  -0.0006   0.0007  -0.1449   0.2150   0.4273  -0.0677

 Orbital           8        9       10
   1 C1  :1s    -0.0745   0.0435  -0.0052
   2 C1  :1s     0.2746  -0.1953   0.0251
   3 C1  :2px   -0.4348   0.1526  -0.0246
   4 C1  :2py    0.4024  -0.1042  -0.0283
   5 C2  :1s     0.0314  -0.0031   0.0374
   6 C2  :1s    -0.1219   0.0057  -0.1897
   7 C2  :2px    0.1486  -0.0362   0.0513
   8 C2  :2py   -0.0801  -0.0412  -0.0409
   9 Li  :1s    -0.1583  -0.2353  -0.0051
  10 Li  :1s     0.3763   0.8953   0.0139
  11 Li  :2px    0.1955  -0.3623   0.6474
  12 Li  :2py   -0.1815   0.2932   0.7733
  13 H1  :1s    -0.2828   0.0379   0.1220
  14 H2  :1s    -0.2041   0.0255  -0.0052
  15 H3  :1s     0.2573  -0.0417   0.0782

 Molecular orbitals for symmetry species  2
 ------------------------------------------

 Orbital           1        2        3
   1 C1  :2pz    0.5690  -0.2001  -0.8533
   2 C2  :2pz    0.6575   0.5041   0.6107
   3 Li  :2pz    0.1553  -0.8595   0.5206



 >>>> Total CPU  time used in SIRIUS :      4.96 seconds
 >>>> Total wall time used in SIRIUS :      5.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:44:25 2009
     Host name              : stallo-2.local                          


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



        *****************************************************************
        ******** Output from **PROPE input processing for ABACUS ********
        *****************************************************************



 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

 Default print level:        0

      Nuclear magnetic shieldings
      Nuclear spin-spin coupling constants
      Natural orbital connection is used
      for perturbation dependent basis sets.

 Center of mass dipole origin  :    0.000000    0.000000    0.000000

 Center of mass gauge origin   :    0.000000    0.000000    0.000000


                 .------------------------------------------------.
                 | Starting in Static Property Section (ABACUS) - |
                 `------------------------------------------------'


 
     Date and time (Linux)  : Thu Sep 24 00:44:25 2009
     Host name              : stallo-2.local                          
 Symmetry -> DSO by analytical integration.
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.7 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.6 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 17.999989(-1.06e-05): LR-DFT*1 evaluation time:       0.5 s


   ***************************************************************************
   ************************ FINAL RESULTS from ABACUS ************************
   ***************************************************************************


 
     Date and time (Linux)  : Thu Sep 24 00:48:08 2009
     Host name              : stallo-2.local                          


                             Molecular geometry (au)
                             -----------------------

 C1        -0.4680105995           -0.7439573407            0.0000000000
 C2         1.9572915950            0.0000000000            0.0000000000
 Li        -3.3287220260            1.5766997527            0.0000000000
 H1        -0.6422661484           -2.8143267967            0.0000000000
 H2         2.5109335755            1.9893441650            0.0000000000
 H3         3.5717199730           -1.2930887194            0.0000000000





                        Molecular wave function and energy
                        ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0

     Total energy        -84.3583024056 au (Hartrees)
                         -2295.50611292 eV
                           -221482.6912 kJ/mol


                             Relativistic corrections
                             ------------------------

     Darwin correction:                          0.0797216293 au
     Mass-velocity correction:                  -0.1004768550 au

     Total relativistic correction:             -0.0207552257 au (0.0246%)
     Non-relativistic + relativistic energy:   -84.3790576313 au




                                  Dipole moment
                                  -------------

                 au               Debye          C m (/(10**-30)
              1.965473           4.995732          16.663969


                             Dipole moment components
                             ------------------------

                 au               Debye          C m (/(10**-30)

      x     -1.58496208        -4.02857131       -13.43786745
      y      1.16231554         2.95431108         9.85452103


   Units:   1 a.u. =   2.54175 Debye 
            1 a.u. =   8.47835 (10**-30) C m (SI)




  ******************************************************************************
  ************************ ABACUS - CHEMICAL SHIELDINGS ************************
  ******************************************************************************



                 Shielding tensors in symmetry coordinates (ppm)
                 -----------------------------------------------


    Symmetry 1
                       Bz

        C1  z    236.18534171
        C2  z    194.34999156
        Li  z     83.13859889
        H1  z     26.77408030
        H2  z     26.52579434
        H3  z     26.71803797

    Symmetry 2
                       Bx             By

        C1  x     65.20395583     9.71568867
        C1  y    -30.90319212    17.68879144
        C2  x    129.54036528    41.43454235
        C2  y     28.62627911    23.31948802
        Li  x     92.55285111    -8.96971208
        Li  y     -8.31206015    89.08154587
        H1  x     29.12289220     1.13588981
        H1  y     -1.86758072    24.02894561
        H2  x     29.83240256    -0.11371006
        H2  y      2.91972582    25.94159032
        H3  x     30.33292024     2.02328099
        H3  y     -1.30413072    24.07960331




  Chemical shielding for C1    :
  ==============================


  Shielding constant:    106.3594 ppm
  Anisotropy:            194.7390 ppm
  Asymmetry:               0.4007    

  S parameter:           199.8830 ppm
  A parameter:            20.3094 ppm


  Total shielding tensor (ppm):
  -----------------------------

 String on input: "C1      "
 String on input: "C1      "
                    Bx             By             Bz

  C1     x    65.20395583     9.71568867     0.00000000
  C1     y   -30.90319212    17.68879144     0.00000000
  C1     z     0.00000000     0.00000000   236.18534171


  Diamagnetic and paramagnetic contributions (ppm):
  -------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  C1     x   332.9747    18.0723     0.0000     -267.7707    -8.3566     0.0000
  C1     y    10.1708   357.5846     0.0000      -41.0740  -339.8958     0.0000
  C1     z     0.0000     0.0000   376.3721        0.0000     0.0000  -140.1868

  Diamagnetic contribution:    355.643810         Paramagnetic:    -249.284447


  Antisymmetric and traceless symmetric parts (ppm):
  --------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  C1     x     0.0000    20.3094     0.0000      -41.1554   -10.5938     0.0000
  C1     y   -20.3094     0.0000     0.0000      -10.5938   -88.6706     0.0000
  C1     z     0.0000     0.0000     0.0000        0.0000     0.0000   129.8260


  Principal values and axes:
  --------------------------

  C1         15.433871  =  106.36  -  90.93:     0.208190  0.978088  0.000000
  C1         67.458876  =  106.36  -  38.90:     0.978088 -0.208190  0.000000
  C1        236.185342  =  106.36  + 129.83:     0.000000  0.000000  1.000000


  Chemical shielding for C2    :
  ==============================


  Shielding constant:    115.7366 ppm
  Anisotropy:           -154.3941 ppm
  Asymmetry:               0.5275    

  S parameter:           161.3961 ppm
  A parameter:             6.4041 ppm


  Total shielding tensor (ppm):
  -----------------------------

 String on input: "C2      "
 String on input: "C2      "
                    Bx             By             Bz

  C2     x   129.54036528    41.43454235     0.00000000
  C2     y    28.62627911    23.31948802     0.00000000
  C2     z     0.00000000     0.00000000   194.34999156


  Diamagnetic and paramagnetic contributions (ppm):
  -------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  C2     x   326.2754    -1.2779     0.0000     -196.7350    42.7125     0.0000
  C2     y    -2.1933   337.2877     0.0000       30.8195  -313.9682     0.0000
  C2     z     0.0000     0.0000   378.3514        0.0000     0.0000  -184.0014

  Diamagnetic contribution:    347.304820         Paramagnetic:    -231.568205


  Antisymmetric and traceless symmetric parts (ppm):
  --------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  C2     x     0.0000     6.4041     0.0000       13.8038    35.0304     0.0000
  C2     y    -6.4041     0.0000     0.0000       35.0304   -92.4171     0.0000
  C2     z     0.0000     0.0000     0.0000        0.0000     0.0000    78.6134


  Principal values and axes:
  --------------------------

  C2        194.349992  =  115.74  +  78.61:     0.000000  0.000000  1.000000
  C2        140.052630  =  115.74  +  24.32:     0.957803  0.287427  0.000000
  C2         12.807223  =  115.74  - 102.93:    -0.287427  0.957803  0.000000


  Chemical shielding for Li    :
  ==============================


  Shielding constant:     88.2577 ppm
  Anisotropy:             17.0595 ppm
  Asymmetry:               0.0998    

  S parameter:            17.0878 ppm
  A parameter:             0.3288 ppm


  Total shielding tensor (ppm):
  -----------------------------

 String on input: "Li      "
 String on input: "Li      "
                    Bx             By             Bz

  Li     x    92.55285111    -8.96971208     0.00000000
  Li     y    -8.31206015    89.08154587     0.00000000
  Li     z     0.00000000     0.00000000    83.13859889


  Diamagnetic and paramagnetic contributions (ppm):
  -------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  Li     x    99.7164   -14.4586     0.0000       -7.1636     5.4889     0.0000
  Li     y    -7.2919    88.7148     0.0000       -1.0202     0.3668     0.0000
  Li     z     0.0000     0.0000    86.4552        0.0000     0.0000    -3.3166

  Diamagnetic contribution:     91.628801         Paramagnetic:      -3.371136


  Antisymmetric and traceless symmetric parts (ppm):
  --------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  Li     x     0.0000    -0.3288     0.0000        4.2952    -8.6409     0.0000
  Li     y     0.3288     0.0000     0.0000       -8.6409     0.8239     0.0000
  Li     z     0.0000     0.0000     0.0000        0.0000     0.0000    -5.1191


  Principal values and axes:
  --------------------------

  Li         82.003720  =   88.26  -   6.25:     0.633667  0.773606  0.000000
  Li         83.138599  =   88.26  -   5.12:     0.000000  0.000000  1.000000
  Li         99.630677  =   88.26  +  11.37:     0.773606 -0.633667  0.000000


  Chemical shielding for H1    :
  ==============================


  Shielding constant:     26.6420 ppm
  Anisotropy:             -3.9588 ppm
  Asymmetry:               0.8999    

  S parameter:             4.4612 ppm
  A parameter:             1.5017 ppm


  Total shielding tensor (ppm):
  -----------------------------

 String on input: "H1      "
 String on input: "H1      "
                    Bx             By             Bz

  H1     x    29.12289220     1.13588981     0.00000000
  H1     y    -1.86758072    24.02894561     0.00000000
  H1     z     0.00000000     0.00000000    26.77408030


  Diamagnetic and paramagnetic contributions (ppm):
  -------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  H1     x    27.0389     5.4377     0.0000        2.0840    -4.3018     0.0000
  H1     y     3.2228    49.3259     0.0000       -5.0904   -25.2969     0.0000
  H1     z     0.0000     0.0000    22.0205        0.0000     0.0000     4.7536

  Diamagnetic contribution:     32.795105         Paramagnetic:      -6.153132


  Antisymmetric and traceless symmetric parts (ppm):
  --------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  H1     x     0.0000     1.5017     0.0000        2.4809    -0.3658     0.0000
  H1     y    -1.5017     0.0000     0.0000       -0.3658    -2.6130     0.0000
  H1     z     0.0000     0.0000     0.0000        0.0000     0.0000     0.1321


  Principal values and axes:
  --------------------------

  H1         29.149033  =   26.64  +   2.51:     0.997457 -0.071271  0.000000
  H1         26.774080  =   26.64  +   0.13:     0.000000  0.000000  1.000000
  H1         24.002805  =   26.64  -   2.64:     0.071271  0.997457  0.000000


  Chemical shielding for H2    :
  ==============================


  Shielding constant:     27.4333 ppm
  Anisotropy:              4.2784 ppm
  Asymmetry:               0.3637    

  S parameter:             4.3717 ppm
  A parameter:             1.5167 ppm


  Total shielding tensor (ppm):
  -----------------------------

 String on input: "H2      "
 String on input: "H2      "
                    Bx             By             Bz

  H2     x    29.83240256    -0.11371006     0.00000000
  H2     y     2.91972582    25.94159032     0.00000000
  H2     z     0.00000000     0.00000000    26.52579434


  Diamagnetic and paramagnetic contributions (ppm):
  -------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  H2     x    30.9779    10.1936     0.0000       -1.1455   -10.3073     0.0000
  H2     y     5.5618    45.7234     0.0000       -2.6420   -19.7818     0.0000
  H2     z     0.0000     0.0000    22.3354        0.0000     0.0000     4.1903

  Diamagnetic contribution:     33.012238         Paramagnetic:      -5.578975


  Antisymmetric and traceless symmetric parts (ppm):
  --------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  H2     x     0.0000    -1.5167     0.0000        2.3991     1.4030     0.0000
  H2     y     1.5167     0.0000     0.0000        1.4030    -1.4917     0.0000
  H2     z     0.0000     0.0000     0.0000        0.0000     0.0000    -0.9075


  Principal values and axes:
  --------------------------

  H2         25.488448  =   27.43  -   1.94:    -0.307346  0.951598  0.000000
  H2         26.525794  =   27.43  -   0.91:     0.000000  0.000000  1.000000
  H2         30.285545  =   27.43  +   2.85:     0.951598  0.307346  0.000000


  Chemical shielding for H3    :
  ==============================


  Shielding constant:     27.0435 ppm
  Anisotropy:              4.9650 ppm
  Asymmetry:               0.8033    

  S parameter:             5.4730 ppm
  A parameter:             1.6637 ppm


  Total shielding tensor (ppm):
  -----------------------------

 String on input: "H3      "
 String on input: "H3      "
                    Bx             By             Bz

  H3     x    30.33292024     2.02328099     0.00000000
  H3     y    -1.30413072    24.07960331     0.00000000
  H3     z     0.00000000     0.00000000    26.71803797


  Diamagnetic and paramagnetic contributions (ppm):
  -------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  H3     x    41.2680   -12.2075     0.0000      -10.9351    14.2308     0.0000
  H3     y    -9.3121    32.0462     0.0000        8.0080    -7.9666     0.0000
  H3     z     0.0000     0.0000     9.8753        0.0000     0.0000    16.8427

  Diamagnetic contribution:     27.729848         Paramagnetic:      -0.686328


  Antisymmetric and traceless symmetric parts (ppm):
  --------------------------------------------------

                 Bx         By         Bz            Bx         By         Bz

  H3     x     0.0000     1.6637     0.0000        3.2894     0.3596     0.0000
  H3     y    -1.6637     0.0000     0.0000        0.3596    -2.9639     0.0000
  H3     z     0.0000     0.0000     0.0000        0.0000     0.0000    -0.3255


  Principal values and axes:
  --------------------------

  H3         24.058995  =   27.04  -   2.98:    -0.057219  0.998362  0.000000
  H3         26.718038  =   27.04  -   0.33:     0.000000  0.000000  1.000000
  H3         30.353528  =   27.04  +   3.31:     0.998362  0.057219  0.000000


                         +--------------------------------+
                         ! Summary of chemical shieldings !
                         +--------------------------------+

 Definitions from J.Mason, Solid state Nuc.Magn.Res. 2 (1993), 285

 @1atom   shielding       dia      para     skew      span     (aniso     asym)
 @1----------------------------------------------------------------------------
 @1C1      106.3594  355.6438 -249.2844    0.5287  220.7515  194.7390    0.4007
 @1C2      115.7366  347.3048 -231.5682   -0.4018  181.5428  117.9201    1.6186
 @1Li       88.2577   91.6288   -3.3711    0.8712   17.6270   17.0595    0.0998
 @1H1       26.6420   32.7951   -6.1531   -0.0770    5.1462    3.7606    1.1054
 @1H2       27.4333   33.0122   -5.5790    0.5675    4.7971    4.2784    0.3637
 @1H3       27.0435   27.7298   -0.6863    0.1551    6.2945    4.9650    0.8033



                         +--------------------------------+
                         ! Summary of chemical shieldings !
                         +--------------------------------+

 Definitions from Smith, Palke and Grieg, Concepts in Mag.Res. 4 (1992), 107

 @2atom   shielding       dia      para     aniso      asym        S        A
 @2----------------------------------------------------------------------------
 @2C1      106.3594  355.6438 -249.2844  194.7390    0.4007  199.8830   20.3094
 @2C2      115.7366  347.3048 -231.5682 -154.3941    0.5275  161.3961    6.4041
 @2Li       88.2577   91.6288   -3.3711   17.0595    0.0998   17.0878    0.3288
 @2H1       26.6420   32.7951   -6.1531   -3.9588    0.8999    4.4612    1.5017
 @2H2       27.4333   33.0122   -5.5790    4.2784    0.3637    4.3717    1.5167
 @2H3       27.0435   27.7298   -0.6863    4.9650    0.8033    5.4730    1.6637



 *******************************************************************************
 **************** ABACUS - INDIRECT NUCLEAR SPIN-SPIN COUPLINGS ****************
 *******************************************************************************

 Definitions from Smith, Palke, and Grieg, Concepts in Mag. Res. 4 (1992) 107



                       +------------------------------------+
                       ! ABACUS - Final spin-spin couplings !
                       +------------------------------------+



              Indirect spin-spin coupling between C2     and C1    :
              ======================================================


  Mass number atom 1:   13;   Abundance:   1.100 %;   g factor:   1.404824
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     45.0578 Hz
  Anisotropic coupling      :    -46.9707 Hz
  Asymmetry                 :     -0.5363   
  S parameter               :     49.1708 Hz
  A parameter               :      0.3265 Hz
  Isotropic DSO contribution:      0.0907 Hz
  Isotropic PSO contribution:     -8.0955 Hz
  Isotropic SD contribution :      3.5190 Hz
  Isotropic FC contribution :     49.5435 Hz




              Indirect spin-spin coupling between Li     and C1    :
              ======================================================


  Mass number atom 1:    7;   Abundance:  92.500 %;   g factor:   2.170951
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     65.6372 Hz
  Anisotropic coupling      :     10.5347 Hz
  Asymmetry                 :      0.3004   
  S parameter               :     10.6920 Hz
  A parameter               :      0.0758 Hz
  Isotropic DSO contribution:      0.0801 Hz
  Isotropic PSO contribution:     -0.0955 Hz
  Isotropic SD contribution :      0.0206 Hz
  Isotropic FC contribution :     65.6321 Hz



  Mass number atom 1:    6;   Abundance:   7.500 %;   g factor:   0.822047
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     24.8540 Hz
  Anisotropic coupling      :      3.9890 Hz
  Asymmetry                 :      0.3004   
  S parameter               :      4.0486 Hz
  A parameter               :      0.0287 Hz
  Isotropic DSO contribution:      0.0303 Hz
  Isotropic PSO contribution:     -0.0362 Hz
  Isotropic SD contribution :      0.0078 Hz
  Isotropic FC contribution :     24.8521 Hz




              Indirect spin-spin coupling between Li     and C2    :
              ======================================================


  Mass number atom 1:    7;   Abundance:  92.500 %;   g factor:   2.170951
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     -3.1516 Hz
  Anisotropic coupling      :     -3.1661 Hz
  Asymmetry                 :     -0.6198   
  S parameter               :      3.3627 Hz
  A parameter               :      0.2286 Hz
  Isotropic DSO contribution:     -0.0895 Hz
  Isotropic PSO contribution:     -0.0028 Hz
  Isotropic SD contribution :      0.0299 Hz
  Isotropic FC contribution :     -3.0893 Hz



  Mass number atom 1:    6;   Abundance:   7.500 %;   g factor:   0.822047
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     -1.1934 Hz
  Anisotropic coupling      :     -1.1989 Hz
  Asymmetry                 :     -0.6198   
  S parameter               :      1.2733 Hz
  A parameter               :      0.0865 Hz
  Isotropic DSO contribution:     -0.0339 Hz
  Isotropic PSO contribution:     -0.0010 Hz
  Isotropic SD contribution :      0.0113 Hz
  Isotropic FC contribution :     -1.1698 Hz




              Indirect spin-spin coupling between H1     and C1    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     87.4167 Hz
  Anisotropic coupling      :     -5.6935 Hz
  Asymmetry                 :     -0.2823   
  S parameter               :      5.7686 Hz
  A parameter               :      0.6624 Hz
  Isotropic DSO contribution:      0.6077 Hz
  Isotropic PSO contribution:     -1.3886 Hz
  Isotropic SD contribution :      0.2743 Hz
  Isotropic FC contribution :     87.9234 Hz




              Indirect spin-spin coupling between H1     and C2    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :    -12.1053 Hz
  Anisotropic coupling      :     -9.7922 Hz
  Asymmetry                 :     -0.4056   
  S parameter               :     10.0571 Hz
  A parameter               :      0.2893 Hz
  Isotropic DSO contribution:     -0.4907 Hz
  Isotropic PSO contribution:     -1.8449 Hz
  Isotropic SD contribution :     -0.2648 Hz
  Isotropic FC contribution :     -9.5049 Hz




              Indirect spin-spin coupling between H1     and Li    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    7;   Abundance:  92.500 %;   g factor:   2.170951

  Isotropic coupling        :      0.9921 Hz
  Anisotropic coupling      :      6.4064 Hz
  Asymmetry                 :      0.4067   
  S parameter               :      6.5807 Hz
  A parameter               :      0.7526 Hz
  Isotropic DSO contribution:     -0.5486 Hz
  Isotropic PSO contribution:      0.1791 Hz
  Isotropic SD contribution :      0.0093 Hz
  Isotropic FC contribution :      1.3522 Hz



  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    6;   Abundance:   7.500 %;   g factor:   0.822047

  Isotropic coupling        :      0.3757 Hz
  Anisotropic coupling      :      2.4258 Hz
  Asymmetry                 :      0.4067   
  S parameter               :      2.4918 Hz
  A parameter               :      0.2850 Hz
  Isotropic DSO contribution:     -0.2077 Hz
  Isotropic PSO contribution:      0.0678 Hz
  Isotropic SD contribution :      0.0035 Hz
  Isotropic FC contribution :      0.5120 Hz




              Indirect spin-spin coupling between H2     and C1    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :     -0.5925 Hz
  Anisotropic coupling      :     -9.0042 Hz
  Asymmetry                 :     -0.6989   
  S parameter               :      9.7096 Hz
  A parameter               :      0.2528 Hz
  Isotropic DSO contribution:     -0.6451 Hz
  Isotropic PSO contribution:     -1.6689 Hz
  Isotropic SD contribution :     -0.1913 Hz
  Isotropic FC contribution :      1.9128 Hz




              Indirect spin-spin coupling between H2     and C2    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :    131.5705 Hz
  Anisotropic coupling      :     -4.5145 Hz
  Asymmetry                 :     -0.1695   
  S parameter               :      4.5360 Hz
  A parameter               :      0.5129 Hz
  Isotropic DSO contribution:      0.6399 Hz
  Isotropic PSO contribution:     -0.3796 Hz
  Isotropic SD contribution :      0.4905 Hz
  Isotropic FC contribution :    130.8198 Hz




              Indirect spin-spin coupling between H2     and Li    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    7;   Abundance:  92.500 %;   g factor:   2.170951

  Isotropic coupling        :      1.9885 Hz
  Anisotropic coupling      :      2.3939 Hz
  Asymmetry                 :      0.4682   
  S parameter               :      2.4798 Hz
  A parameter               :      0.5518 Hz
  Isotropic DSO contribution:     -0.1158 Hz
  Isotropic PSO contribution:      0.0054 Hz
  Isotropic SD contribution :     -0.0176 Hz
  Isotropic FC contribution :      2.1165 Hz



  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    6;   Abundance:   7.500 %;   g factor:   0.822047

  Isotropic coupling        :      0.7530 Hz
  Anisotropic coupling      :      0.9065 Hz
  Asymmetry                 :      0.4682   
  S parameter               :      0.9390 Hz
  A parameter               :      0.2089 Hz
  Isotropic DSO contribution:     -0.0438 Hz
  Isotropic PSO contribution:      0.0021 Hz
  Isotropic SD contribution :     -0.0067 Hz
  Isotropic FC contribution :      0.8014 Hz




              Indirect spin-spin coupling between H2     and H1    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    1;   Abundance:  99.985 %;   g factor:   5.585694

  Isotropic coupling        :     23.1270 Hz
  Anisotropic coupling      :      8.0201 Hz
  Asymmetry                 :      0.4038   
  S parameter               :      8.2352 Hz
  A parameter               :      0.1671 Hz
  Isotropic DSO contribution:     -3.3795 Hz
  Isotropic PSO contribution:      0.0627 Hz
  Isotropic SD contribution :      0.5229 Hz
  Isotropic FC contribution :     25.9209 Hz




              Indirect spin-spin coupling between H3     and C1    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :    -14.2392 Hz
  Anisotropic coupling      :      9.0661 Hz
  Asymmetry                 :      0.8675   
  S parameter               :     10.1398 Hz
  A parameter               :      0.2467 Hz
  Isotropic DSO contribution:     -0.6848 Hz
  Isotropic PSO contribution:     -2.0213 Hz
  Isotropic SD contribution :     -0.2841 Hz
  Isotropic FC contribution :    -11.2490 Hz




              Indirect spin-spin coupling between H3     and C2    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:   13;   Abundance:   1.100 %;   g factor:   1.404824

  Isotropic coupling        :    118.6129 Hz
  Anisotropic coupling      :     -5.6201 Hz
  Asymmetry                 :     -0.4840   
  S parameter               :      5.8354 Hz
  A parameter               :      0.5521 Hz
  Isotropic DSO contribution:      0.6484 Hz
  Isotropic PSO contribution:     -0.5689 Hz
  Isotropic SD contribution :      0.5063 Hz
  Isotropic FC contribution :    118.0272 Hz




              Indirect spin-spin coupling between H3     and Li    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    7;   Abundance:  92.500 %;   g factor:   2.170951

  Isotropic coupling        :      7.0741 Hz
  Anisotropic coupling      :      1.3013 Hz
  Asymmetry                 :      0.3169   
  S parameter               :      1.3229 Hz
  A parameter               :      0.1029 Hz
  Isotropic DSO contribution:     -0.6395 Hz
  Isotropic PSO contribution:      0.1278 Hz
  Isotropic SD contribution :      0.0364 Hz
  Isotropic FC contribution :      7.5495 Hz



  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    6;   Abundance:   7.500 %;   g factor:   0.822047

  Isotropic coupling        :      2.6787 Hz
  Anisotropic coupling      :      0.4927 Hz
  Asymmetry                 :      0.3169   
  S parameter               :      0.5009 Hz
  A parameter               :      0.0390 Hz
  Isotropic DSO contribution:     -0.2422 Hz
  Isotropic PSO contribution:      0.0484 Hz
  Isotropic SD contribution :      0.0138 Hz
  Isotropic FC contribution :      2.8587 Hz




              Indirect spin-spin coupling between H3     and H1    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    1;   Abundance:  99.985 %;   g factor:   5.585694

  Isotropic coupling        :     14.8891 Hz
  Anisotropic coupling      :     10.1109 Hz
  Asymmetry                 :      0.7435   
  S parameter               :     11.0029 Hz
  A parameter               :      3.9399 Hz
  Isotropic DSO contribution:     -0.6597 Hz
  Isotropic PSO contribution:     -0.4670 Hz
  Isotropic SD contribution :      0.1280 Hz
  Isotropic FC contribution :     15.8878 Hz




              Indirect spin-spin coupling between H3     and H2    :
              ======================================================


  Mass number atom 1:    1;   Abundance:  99.985 %;   g factor:   5.585694
  Mass number atom 2:    1;   Abundance:  99.985 %;   g factor:   5.585694

  Isotropic coupling        :      1.8958 Hz
  Anisotropic coupling      :     22.3313 Hz
  Asymmetry                 :      0.7384   
  S parameter               :     24.2762 Hz
  A parameter               :      6.4820 Hz
  Isotropic DSO contribution:     -3.1966 Hz
  Isotropic PSO contribution:      0.7354 Hz
  Isotropic SD contribution :      0.5007 Hz
  Isotropic FC contribution :      3.8562 Hz




   Interatomic separations (in Angstrom):
   --------------------------------------

            C1          C2          Li          H1          H2          H3    
            ------      ------      ------      ------      ------      ------
 C1    :    0.000000
 C2    :    1.342439    0.000000
 Li    :    1.949292    2.919021    0.000000
 H1    :    1.099466    2.027387    2.724012    0.000000
 H2    :    2.139411    1.092723    3.097918    3.040717    0.000000
 H3    :    2.157393    1.094573    3.954756    2.370799    1.825442    0.000000


  Max interatomic separation is    3.9548 Angstrom (    7.4734 Bohr)
  between atoms    6 and    3, "H3    " and "Li    ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C2         C1           1.342439
  bond distance:  Li         C1           1.949292
  bond distance:  H1         C1           1.099466
  bond distance:  H2         C2           1.092723
  bond distance:  H3         C2           1.094573


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     C2         C1         Li           123.897
  bond angle:     C2         C1         H1           111.864
  bond angle:     Li         C1         H1           124.239
  bond angle:     C1         C2         H2           122.605
  bond angle:     C1         C2         H3           124.253
  bond angle:     H2         C2         H3           113.141




 CPU time statistics for ABACUS
 ------------------------------

 RHSIDE     00:00:01       1 %
 LINRES     00:00:59      27 %
 TRP LR     00:02:37      72 %

 TOTAL      00:03:38     100 %


 >>>> Total CPU  time used in ABACUS:  3 minutes 38 seconds
 >>>> Total wall time used in ABACUS:  3 minutes 43 seconds


                   .-------------------------------------------.
                   | End of Static Property Section (ABACUS) - |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  3 minutes 43 seconds
 >>>> Total wall time used in DALTON:  3 minutes 48 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:48:08 2009
     Host name              : stallo-2.local                          
END REFOUT

