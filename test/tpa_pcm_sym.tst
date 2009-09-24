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
**INTEGRALS
**RESPONSE
*QUADRATIC
.DIPLEN
.TWO-PHOTON
.ROOTS
 3 3 3 3 
*END OF
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

     Date and time (Linux)  : Thu Sep 24 00:48:17 2009
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
  Number of symmetry independent centres:    1
  The basis set is "STO-3G" from the basis set library.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"
  Huckel basis read for this type.

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  The basis set is "STO-3G" from the basis set library.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"
  Huckel basis read for this type.

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    1
  The basis set is "STO-3G" from the basis set library.
  Used basis set file for basis set for elements with Z =   1 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"
  Huckel basis read for this type.


                         SYMGRP: Point group information
                         -------------------------------

Point group: C2v

   * The point group was generated by:

      Reflection in the yz-plane
      Reflection in the xz-plane

   * Group multiplication table

        |  E   C2z  Oxz  Oyz
   -----+--------------------
     E  |  E 
    C2z | C2z   E 
    Oxz | Oxz  Oyz   E 
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
    A1  | A1 
    B1  | B1   A1 
    B2  | B2   A2   A1 
    A2  | A2   B2   B1   A1 
PCMSPHGEN: SPHERE CENTERS FROM INPUT
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000e+00    0.0000000000e+00    0.0000000000e+00    6.0000000000e+00


                                 Isotopic Masses
                                 ---------------

                           C          12.000000
                           O          15.994915
                           H    1      1.007825
                           H    2      1.007825

                       Total mass:    30.010565 amu
                       Natural abundance:  98.633 %

 Center-of-mass coordinates (a.u.):    0.000000    0.000000    1.159649


  Atoms and basis sets
  --------------------

  Number of atom types:     3
  Total number of atoms:    4

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
  C       :    1  x   0.0000000000   2  y   0.0000000000   3  z   0.0000000000
  O       :    4  x   0.0000000000   5  y   0.0000000000   6  z   2.3054658834
  H   / 1 :    7  x   1.7822044964   8  y   0.0000000000   9  z  -1.0289558799
  H   / 2 :   10  x  -1.7822044964  11  y   0.0000000000  12  z  -1.0289558799


  Symmetry Coordinates
  --------------------

  Number of coordinates in each symmetry:    4   4   3   1

  Symmetry  A1  ( 1)

   1   C     z    3
   2   O     z    6
   3   H     x    [  7  -  10 ]/2
   4   H     z    [  9  +  12 ]/2

  Symmetry  B1  ( 2)

   5   C     x    1
   6   O     x    4
   7   H     x    [  7  +  10 ]/2
   8   H     z    [  9  -  12 ]/2

  Symmetry  B2  ( 3)

   9   C     y    2
  10   O     y    5
  11   H     y    [  8  +  11 ]/2

  Symmetry  A2  ( 4)

  12   H     y    [  8  -  11 ]/2


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H    1      H    2
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H    1:    1.089000    2.000725    0.000000
 H    2:    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstrom (    3.7808 Bohr)
  between atoms    3 and    2, "H    1" and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.220000
  bond distance:  H    1     C            1.089000
  bond distance:  H    2     C            1.089000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H    1       120.000
  bond angle:     O          C          H    2       120.000
  bond angle:     H    1     C          H    2       120.000




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       1.792803          0.000000    0.000000    1.000000
   IB      13.103106          1.000000    0.000000    0.000000
   IC      14.895908          0.000000    1.000000    0.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         281893.2924          38569.4057          33927.3707 MHz
            9.402948            1.286537            1.131695 cm-1


  Nuclear repulsion energy :   31.163673581965


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:               7   3   2   0


  Symmetry  A1 ( 1)

    1     C        1s         1
    2     C        1s         2
    3     C        2pz        5
    4     O        1s         6
    5     O        1s         7
    6     O        2pz       10
    7     H        1s        11  +  12


  Symmetry  B1 ( 2)

    8     C        2px        3
    9     O        2px        8
   10     H        1s        11  -  12


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

 Center of mass :      0.000000000000      0.000000000000      1.159648807704
 Operator center:      0.000000000000      0.000000000000      0.000000000000
 Gauge origin   :      0.000000000000      0.000000000000      0.000000000000
 Dipole origin  :      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014


 MEMORY USED TO GENERATE CAVITY =    432042


 Total number of spheres =    1
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1    0.000000000    0.000000000    0.000000000    3.175062000  126.681817203

 Total number of tesserae =     392
 Surface area =  126.68181720 (A^2)    Cavity volume =  134.07420796 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....
 
  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   0.10 seconds
 >>>> Total wall time used in HERMIT:   0.00 seconds


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
          -20.808771     -11.553164      -2.104598      -1.278134      -0.584721
            0.240410       0.392263

 Huckel EHT eigenvalues for symmetry :  2
           -1.053261      -0.439652       0.161628

 Huckel EHT eigenvalues for symmetry :  3
           -0.811063      -0.212037

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Sep 24 00:48:17 2009
     Host name              : stallo-2.local                          

 Title lines from ".mol" input file:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.


     Time-dependent Hartree-Fock calculation (random phase approximation).


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         14
     Number of electrons in active shells      0
     Total charge of the molecule              2

     Spin multiplicity                         1
     Total number of symmetries                4
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species          All    1    2    3    4
                                       --  --  --  --
     Total number of orbitals           12    7    3    2    0
     Number of basis functions          12    7    3    2    0

      ** Automatic occupation of RHF orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals               7    5    1    1    0

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.99669
 NUCLEAR APPARENT CHARGE -15.79079
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  14 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00e-15
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.44278    38.17950    38.17950   -21.79448    -0.33637
   1  -111.152808841     -0.336366680037       3.79e+00  -1.11e+02    5  1  1  0
 MULPOP C    23.42; O     5.43; H     3.15; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.72699    38.32822    38.32822   -21.79448    -0.32975
   2  -111.516982457     -0.329753521569       6.40e-01  -3.64e-01    5  1  1  0
 MULPOP C    18.75; O     4.77; H     3.43; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.54188    38.23295    38.23295   -21.79448    -0.33247
   3  -111.539004080     -0.332469505690       5.21e-01  -2.20e-02    5  1  1  0
 MULPOP C    21.50; O     5.14; H     3.32; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.62132    38.27463    38.27463   -21.79448    -0.33051
   4  -111.559549399     -0.330510284504       1.20e-01  -2.05e-02    5  1  1  0
 MULPOP C    20.28; O     4.99; H     3.38; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.61468    38.27125    38.27125   -21.79448    -0.33058
   5  -111.565209024     -0.330578566265       4.63e-02  -5.66e-03    5  1  1  0
 MULPOP C    20.29; O     5.00; H     3.38; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.60920    38.26844    38.26844   -21.79448    -0.33065
   6  -111.566150765     -0.330649540360       6.72e-03  -9.42e-04    5  1  1  0
 MULPOP C    20.27; O     5.01; H     3.38; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.60984    38.26876    38.26876   -21.79448    -0.33064
   7  -111.566154971     -0.330638386356       2.13e-03  -4.21e-06    5  1  1  0
 MULPOP C    20.24; O     5.00; H     3.38; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.60966    38.26867    38.26867   -21.79448    -0.33064
   8  -111.566156212     -0.330640865400       2.08e-04  -1.24e-06    5  1  1  0
 MULPOP C    20.22; O     5.00; H     3.38; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.60965    38.26867    38.26867   -21.79448    -0.33064
   9  -111.566156221     -0.330640942853       6.04e-05  -8.70e-09    5  1  1  0
 MULPOP C    20.22; O     5.00; H     3.38; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -33.60965    38.26867    38.26867   -21.79448    -0.33064
  10  -111.566156222     -0.330640969195       6.98e-06  -1.42e-09    5  1  1  0
 DIIS converged in  10 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   14
 Orbital occupations :    5    1    1    0

 Sym       Hartree-Fock orbital energies

  1    -21.24419245   -11.75700816    -1.96710589    -1.28470083    -1.11212508
         0.20015242     0.27610540

  2     -1.09098281    -0.44093851     0.29510715

  3     -1.08694639    -0.26242740

    E(LUMO) :    -0.44093851 au (symmetry 2)
  - E(HOMO) :    -1.08694639 au (symmetry 3)
  ------------------------------------------
    gap     :     0.64600789 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    2

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -111.566156222260                 
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -142.399188835030

     Final gradient norm:           0.000006984615

 
     Date and time (Linux)  : Thu Sep 24 00:48:18 2009
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s     0.0005   0.9929  -0.1100  -0.1963  -0.0602   0.2179  -0.0305
   2 C   :1s    -0.0068   0.0318   0.2171   0.6406   0.2459  -1.4199   0.1822
   3 C   :2pz   -0.0058   0.0020   0.1441  -0.0947  -0.4750  -0.0428   1.2445
   4 O   :1s     0.9947   0.0002  -0.2211   0.1213  -0.0622  -0.0595   0.0949
   5 O   :1s     0.0242  -0.0061   0.7939  -0.5610   0.3266   0.4288  -0.7420
   6 O   :2pz   -0.0055   0.0019  -0.2460  -0.3804   0.6286  -0.4932   0.7622
   7 H   :1s     0.0003  -0.0061   0.0185   0.1501   0.1698   0.8242   0.4352

 Molecular orbitals for symmetry species  2
 ------------------------------------------

 Orbital           1        2        3
   1 C   :2px    0.6742  -0.1767  -1.1020
   2 O   :2px    0.3443   0.9095   0.3278
   3 H   :1s     0.2176  -0.3096   0.8999

 Molecular orbitals for symmetry species  3
 ------------------------------------------

 Orbital           1        2
   1 C   :2py    0.3159   0.9721
   2 O   :2py    0.8855  -0.5106



 >>>> Total CPU  time used in SIRIUS :      0.70 seconds
 >>>> Total wall time used in SIRIUS :      1.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:48:18 2009
     Host name              : stallo-2.local                          


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


   SCF energy         :     -111.566156222259778
 -- inactive part     :     -142.399188835030031
 -- nuclear repulsion :       31.163673581965142


                     ***************************************
                     *** RHF response calculation (TDHF) ***
                     ***************************************



 Linear response excitations for quadratic response
 - symmetry of excitation operator    1
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      13
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      13



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   22)
 RSP solution vector no.    1; norm of residual   3.71e-05
 RSP solution vector no.    2; norm of residual   1.31e-05
 RSP solution vector no.    3; norm of residual   9.90e-05

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    2
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  2; triplet =   F


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   20)
 RSP solution vector no.    1; norm of residual   4.20e-04
 RSP solution vector no.    2; norm of residual   2.08e-04
 RSP solution vector no.    3; norm of residual   1.22e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    3
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  3; triplet =   F

RSPORT:    9 out of    3 new trial vectors linear dependent

 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   14)
 RSP solution vector no.    1; norm of residual   8.14e-15
 RSP solution vector no.    2; norm of residual   8.44e-15
 RSP solution vector no.    3; norm of residual   8.33e-15

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    4
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       4
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       3
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       3



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  4; triplet =   F


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   5.54e-16
 RSP solution vector no.    2; norm of residual   4.82e-16
 RSP solution vector no.    3; norm of residual   1.05e-15

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    1


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      13
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      13


 QRLRVE -- linear response calculation for symmetry  1
 QRLRVE -- operator label : ZDIPLEN 
 QRLRVE -- frequencies :  0.156892  0.260465  0.472415  0.077257  0.255638
                          0.432130  0.206986  0.331146  0.396125



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after    9 linear transformations is     19.68372555
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.01 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after    9 linear transformations is     13.03231897
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.68 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.33115 a.u.
 after    9 linear transformations is    -13.05847120
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -8.05 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39612 a.u.
 after    9 linear transformations is      9.11264101
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.65 * 10 **  6.0
RSPORT:   18 out of    9 new trial vectors linear dependent

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after   13 linear transformations is     26.40443728
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -8.47 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after   13 linear transformations is     16.01872997
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -1.31 * 10 ** 16.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.33115 a.u.
 after   13 linear transformations is     -7.93325761
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -4.94 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39612 a.u.
 after   13 linear transformations is     11.27082630
 INERTIA (POS,ZER,NEG) of reduced matrix is   25    0    1
 Determinant of reduced matrix is -1.40 * 10 ** 16.0

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   26)
 RSP solution vector no.    1; norm of residual   6.52e-16
 RSP solution vector no.    2; norm of residual   6.18e-16
 RSP solution vector no.    3; norm of residual   6.37e-16
 RSP solution vector no.    4; norm of residual   1.11e-15
 RSP solution vector no.    5; norm of residual   4.41e-16
 RSP solution vector no.    6; norm of residual   4.90e-16
 RSP solution vector no.    7; norm of residual   9.30e-16
 RSP solution vector no.    8; norm of residual   4.25e-16
 RSP solution vector no.    9; norm of residual   5.91e-16

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.156892e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.15689):     12.3942169019    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.260465e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.26047):     17.8349073326    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.472415e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.47241):     26.4044372848    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.772568e-01
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.07726):     11.3828680519    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.255638e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.25564):     17.2062330810    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.432130e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.43213):     16.0187299662    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.206986e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.20699):     13.8309243723    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.331146e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.33115):    -7.93325760771    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.396125e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.39612):     11.2708263006    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    2


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12


 QRLRVE -- linear response calculation for symmetry  2
 QRLRVE -- operator label : XDIPLEN 
 QRLRVE -- frequencies :  0.077257  0.156892  0.255638  0.260465  0.432130
                          0.472415  0.025157  0.177758  0.518616



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F


 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.15689 a.u.
 after    9 linear transformations is      5.40452767
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.74 * 10 **  4.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.25564 a.u.
 after    9 linear transformations is      8.03761927
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -1.07 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.26047 a.u.
 after    9 linear transformations is      8.10164528
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -1.10 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after    9 linear transformations is     15.11517979
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -7.71 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after    9 linear transformations is     25.20361843
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -3.80 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.17776 a.u.
 after    9 linear transformations is      7.17369554
 INERTIA (POS,ZER,NEG) of reduced matrix is   17    0    1
 Determinant of reduced matrix is -2.80 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51862 a.u.
 after    9 linear transformations is    -97.17079541
 INERTIA (POS,ZER,NEG) of reduced matrix is   16    0    2
 Determinant of reduced matrix is  6.67 * 10 **  4.0
RSPORT:   18 out of    9 new trial vectors linear dependent

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.15689 a.u.
 after   12 linear transformations is      5.49641804
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -2.97 * 10 ** 13.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.25564 a.u.
 after   12 linear transformations is      8.03777866
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.10 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.26047 a.u.
 after   12 linear transformations is      8.10180542
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.14 * 10 ** 15.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.43213 a.u.
 after   12 linear transformations is     15.11732863
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -7.88 * 10 ** 14.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after   12 linear transformations is     25.21451858
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -3.87 * 10 ** 14.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.17776 a.u.
 after   12 linear transformations is      7.17468347
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -2.90 * 10 ** 14.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51862 a.u.
 after   12 linear transformations is    -96.77752466
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  6.81 * 10 ** 13.0

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   24)
 RSP solution vector no.    1; norm of residual   8.94e-16
 RSP solution vector no.    2; norm of residual   5.04e-16
 RSP solution vector no.    3; norm of residual   6.84e-16
 RSP solution vector no.    4; norm of residual   4.86e-16
 RSP solution vector no.    5; norm of residual   4.80e-16
 RSP solution vector no.    6; norm of residual   4.95e-16
 RSP solution vector no.    7; norm of residual   5.19e-16
 RSP solution vector no.    8; norm of residual   7.01e-16
 RSP solution vector no.    9; norm of residual   6.45e-16

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.772568e-01
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.07726):     6.96541538785    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.156892e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.15689):     5.49641803693    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.255638e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.25564):     8.03777866422    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.260465e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.26047):     8.10180541558    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.432130e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.43213):     15.1173286348    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.472415e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.47241):     25.2145185789    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.251574e-01
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.02516):     6.86568169569    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.177758e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.17776):     7.17468347234    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.518616e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.51862):    -96.7775246558    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    3


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7


 QRLRVE -- linear response calculation for symmetry  3
 QRLRVE -- operator label : YDIPLEN 
 QRLRVE -- frequencies :  0.260465  0.396125  0.472415  0.025157  0.206986
                          0.177758  0.156892  0.518616  0.331146



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F

RSPORT:    9 out of    9 new trial vectors linear dependent

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.47241 a.u.
 after    7 linear transformations is      4.42861270
 INERTIA (POS,ZER,NEG) of reduced matrix is   13    0    1
 Determinant of reduced matrix is -4.02 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51862 a.u.
 after    7 linear transformations is      6.12077084
 INERTIA (POS,ZER,NEG) of reduced matrix is   13    0    1
 Determinant of reduced matrix is -4.83 * 10 **  6.0

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   14)
 RSP solution vector no.    1; norm of residual   5.63e-16
 RSP solution vector no.    2; norm of residual   3.35e-16
 RSP solution vector no.    3; norm of residual   3.63e-16
 RSP solution vector no.    4; norm of residual   3.89e-16
 RSP solution vector no.    5; norm of residual   8.42e-16
 RSP solution vector no.    6; norm of residual   7.23e-16
 RSP solution vector no.    7; norm of residual   4.11e-16
 RSP solution vector no.    8; norm of residual   5.01e-16
 RSP solution vector no.    9; norm of residual   4.10e-16

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.260465e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.26047):     3.40522259232    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.396125e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.39612):     6.57839980982    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.472415e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.47241):     4.42861270342    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.251574e-01
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.02516):     2.80741207908    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.206986e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.20699):     3.14758664346    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.177758e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.17776):     3.04733602370    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.156892e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.15689):     2.98897037487    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.518616e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.51862):     6.12077083989    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.331146e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.33115):     4.02687792546    


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156892    0.313784    4.770877


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260465    0.520931   -1.489762


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    3.196850


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156892    0.313784   20.930192


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260465    0.520931   -0.389092


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    2.952140


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156892    0.313784   -0.232313


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260465    0.520931   -0.131525


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    1.879662


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.077257    0.154514    1.746154


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.255638    0.511276   -9.362874


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.432130    0.864261   -0.769917


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.077257    0.154514    1.746153


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.255638    0.511276   -9.362875


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.432130    0.864261   -0.769917


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.206986    0.413972   -1.245257


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.331146    0.662293   -5.788781


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.396125    0.792249    2.172690


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.206986    0.413972   -1.245257


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.331146    0.662293   -5.788778


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.396125    0.792249    2.172690


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.025157    0.050315    0.047545


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.177758    0.355516   -1.307603


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.518616    1.037231    1.149546


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.025157    0.050315    0.047546


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.177758    0.355516   -1.307604


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.518616    1.037231    1.149544


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
   section which is given in GM.                                

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
       4   2    9.67     0.0     0.0     0.0    -1.3     0.0     0.0
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
     1   1    8.54   Linear       0.216e+02  0.154e+02  0.105e+03  0.402e+00    0.47
     1   1    8.54   Circular     0.216e+02  0.154e+02  0.489e+02  0.188e+00    0.47
     1   2   14.18   Linear       0.135e+00  0.796e-01  0.588e+00  0.623e-02    0.35
     1   2   14.18   Circular     0.135e+00  0.796e-01  0.208e+00  0.220e-02    0.35
     1   3   25.71   Linear       0.215e+01  0.749e+00  0.729e+01  0.254e+00    0.03
     1   3   25.71   Circular     0.215e+01  0.749e+00  0.196e+00  0.684e-02    0.03
     2   1    4.20   Linear       0.000e+00  0.203e+00  0.813e+00  0.758e-03    1.50
     2   1    4.20   Circular     0.000e+00  0.203e+00  0.122e+01  0.114e-02    1.50
     2   2   13.91   Linear       0.000e+00  0.584e+01  0.234e+02  0.238e+00    1.50
     2   2   13.91   Circular     0.000e+00  0.584e+01  0.351e+02  0.358e+00    1.50
     2   3   23.52   Linear       0.000e+00  0.395e-01  0.158e+00  0.461e-02    1.50
     2   3   23.52   Circular     0.000e+00  0.395e-01  0.237e+00  0.691e-02    1.50
     3   1   11.26   Linear       0.000e+00  0.103e+00  0.414e+00  0.277e-02    1.50
     3   1   11.26   Circular     0.000e+00  0.103e+00  0.620e+00  0.415e-02    1.50
     3   2   18.02   Linear       0.000e+00  0.223e+01  0.894e+01  0.153e+00    1.50
     3   2   18.02   Circular     0.000e+00  0.223e+01  0.134e+02  0.229e+00    1.50
     3   3   21.56   Linear       0.000e+00  0.315e+00  0.126e+01  0.308e-01    1.50
     3   3   21.56   Circular     0.000e+00  0.315e+00  0.189e+01  0.463e-01    1.50
     4   1    1.37   Linear       0.000e+00  0.151e-03  0.603e-03  0.596e-07    1.50
     4   1    1.37   Circular     0.000e+00  0.151e-03  0.904e-03  0.893e-07    1.50
     4   2    9.67   Linear       0.000e+00  0.114e+00  0.456e+00  0.225e-02    1.50
     4   2    9.67   Circular     0.000e+00  0.114e+00  0.684e+00  0.337e-02    1.50
     4   3   28.22   Linear       0.000e+00  0.881e-01  0.352e+00  0.148e-01    1.50
     4   3   28.22   Circular     0.000e+00  0.881e-01  0.529e+00  0.222e-01    1.50
   ---------------------------------------------------------------------------------

 >>>> Total CPU  time used in RESPONSE:  16.98 seconds
 >>>> Total wall time used in RESPONSE:  36.00 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  17.82 seconds
 >>>> Total wall time used in DALTON:  37.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:48:54 2009
     Host name              : stallo-2.local                          
END REFOUT

