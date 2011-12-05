########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm hf residue symmetry linear response essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enehf
pcmsol
surf
nucchg
tramom
OVERRIDE thr 1.0e-4  1.0e-5
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
.ICESPH
3
*PCMCAV
.AREATS
0.3
**WAVEFUNCTION
.HF
**RESPONSE
*LINEAR
.THCLR
1.0D-4
.DIPLEN
.SINGLE RESIDUE
.ROOTS
3 3 3 3
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS


Atomtypes=3 Charge=+2 Angstrom Generators=2 X Y
Charge=6.0 Atoms=1 Basis=STO-3G  Radius=1.6 Alpha=1.2
C     0.000000     0.000000     0.000000   1.60
Charge=8.0 Atoms=1 Basis=STO-3G  Radius=1.5 Alpha=1.2
O     0.000000     0.000000     1.220000
Charge=1.0 Atoms=1 Basis=STO-3G  Radius=1.2 Alpha=1.2
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

     Date and time (Linux)  : Thu Sep 24 00:57:04 2009
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

     ICESPH =       3     NESFP =       0
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
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
  Number of symmetry independent centres:    1
  The basis set is "STO-3G  Sphere=1" from the basis set library.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"
  Huckel basis read for this type.

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  The basis set is "STO-3G  Sphere=1" from the basis set library.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"
  Huckel basis read for this type.

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    1
  The basis set is "STO-3G  Sphere=1" from the basis set library.
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
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000e+00    0.0000000000e+00    0.0000000000e+00    1.6000000000e+00
   2    0.0000000000e+00    0.0000000000e+00    2.3054658834e+00    1.5000000000e+00
   3    1.7822044964e+00    0.0000000000e+00   -1.0289558799e+00    1.2000000000e+00
   4   -1.7822044964e+00    0.0000000000e+00   -1.0289558799e+00    1.2000000000e+00


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



     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014


 MEMORY USED TO GENERATE CAVITY =    432042

 Tessera cut in pieces and removed.

 Total number of spheres =    4
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1    0.000000000    0.000000000    0.000000000    1.920000000   17.404313423
   2    0.000000000    0.000000000    1.220000000    1.800000000   25.187332839
   3    0.943102000    0.000000000   -0.544500000    1.440000000   11.255487319
   4   -0.943102000    0.000000000   -0.544500000    1.440000000   11.255487319

 Total number of tesserae =     268
 Surface area =   65.10262090 (A^2)    Cavity volume =   45.75393835 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....
 
  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   0.07 seconds
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

 
     Date and time (Linux)  : Thu Sep 24 00:57:04 2009
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

 ESTIMATE OF NUCLEAR CHARGE       15.96782
 NUCLEAR APPARENT CHARGE -15.78910
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  14 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00e-15
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.33014    54.10741    54.10737   -30.92953    -0.48721
   1  -111.303653266     -0.487211104692       3.78e+00  -1.11e+02    5  1  1  0
 MULPOP C    23.42; O     5.43; H     3.15; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.72041    54.29944    54.29936   -30.92953    -0.49033
   2  -111.669295503     -0.490334842105       6.21e-01  -3.66e-01    5  1  1  0
 MULPOP C    18.66; O     4.76; H     3.44; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.61268    54.25841    54.25841   -30.92953    -0.47746
   3  -111.694107287     -0.477462552395       4.48e-01  -2.48e-02    5  1  1  0
 MULPOP C    21.27; O     5.10; H     3.33; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.65574    54.27590    54.27587   -30.92953    -0.48152
   4  -111.710282693     -0.481523646370       1.16e-01  -1.62e-02    5  1  1  0
 MULPOP C    20.14; O     4.97; H     3.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.63986    54.26809    54.26806   -30.92953    -0.48139
   5  -111.715563750     -0.481393150547       4.51e-02  -5.28e-03    5  1  1  0
 MULPOP C    20.18; O     4.98; H     3.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.62756    54.26197    54.26194   -30.92953    -0.48136
   6  -111.716477294     -0.481356510663       7.88e-03  -9.14e-04    5  1  1  0
 MULPOP C    20.13; O     4.98; H     3.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.62787    54.26209    54.26207   -30.92953    -0.48139
   7  -111.716484095     -0.481386440653       2.23e-03  -6.80e-06    5  1  1  0
 MULPOP C    20.09; O     4.98; H     3.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.62733    54.26183    54.26180   -30.92953    -0.48138
   8  -111.716485628     -0.481381694427       1.88e-04  -1.53e-06    5  1  1  0
 MULPOP C    20.08; O     4.98; H     3.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.62733    54.26183    54.26180   -30.92953    -0.48138
   9  -111.716485634     -0.481381191541       5.08e-05  -5.95e-09    5  1  1  0
 MULPOP C    20.08; O     4.98; H     3.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -47.62732    54.26183    54.26180   -30.92953    -0.48138
  10  -111.716485635     -0.481381338699       5.55e-06  -9.80e-10    5  1  1  0
 DIIS converged in  10 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   14
 Orbital occupations :    5    1    1    0

 Sym       Hartree-Fock orbital energies

  1    -21.08782702   -11.61411487    -1.81843564    -1.14185418    -0.97076140
         0.34807001     0.42545034

  2     -0.93988408    -0.29794360     0.44229643

  3     -0.94384861    -0.11091771

    E(LUMO) :    -0.29794360 au (symmetry 2)
  - E(HOMO) :    -0.93988408 au (symmetry 2)
  ------------------------------------------
    gap     :     0.64194049 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    2

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -111.716485634508                 
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -142.398777877774

     Final gradient norm:           0.000005545430

 
     Date and time (Linux)  : Thu Sep 24 00:57:05 2009
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s     0.0005   0.9929  -0.1088   0.1926  -0.0715   0.2178  -0.0333
   2 C   :1s    -0.0068   0.0317   0.2142  -0.6223   0.2853  -1.4190   0.1985
   3 C   :2pz   -0.0058   0.0020   0.1429   0.0704  -0.4767  -0.0304   1.2459
   4 O   :1s     0.9947   0.0002  -0.2225  -0.1222  -0.0533  -0.0591   0.0963
   5 O   :1s     0.0242  -0.0060   0.8000   0.5677   0.2851   0.4238  -0.7504
   6 O   :2pz   -0.0054   0.0019  -0.2375   0.4234   0.6099  -0.4830   0.7639
   7 H   :1s     0.0003  -0.0061   0.0182  -0.1406   0.1779   0.8287   0.4265

 Molecular orbitals for symmetry species  2
 ------------------------------------------

 Orbital           1        2        3
   1 C   :2px    0.6699  -0.1864  -1.1030
   2 O   :2px    0.3562   0.9055   0.3261
   3 H   :1s     0.2154  -0.3113   0.8999

 Molecular orbitals for symmetry species  3
 ------------------------------------------

 Orbital           1        2
   1 C   :2py    0.3049   0.9757
   2 O   :2py    0.8913  -0.5006



 >>>> Total CPU  time used in SIRIUS :      0.54 seconds
 >>>> Total wall time used in SIRIUS :      1.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:57:05 2009
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




  Linear Response single residue calculation
 -------------------------------------------

 Equilibrium PCM solvent model requested        : SOLVNT =T

 Dielectric constant                            : EPSOL  = 78.3900
 Print level                                    : IPRPP  =   2
 Maximum number of iterations for eigenval.eqs. : MAXITP =  60
 Threshold for convergence of eigenvalue eqs.   : THCPP  = 1.000e-03
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

      3 Excitation energies are calculated for symmetry no.    1

      1 property residues are calculated with labels:

               ZDIPLEN 

      3 Excitation energies are calculated for symmetry no.    2

      1 property residues are calculated with labels:

               XDIPLEN 

      3 Excitation energies are calculated for symmetry no.    3

      1 property residues are calculated with labels:

               YDIPLEN 

      3 Excitation energies are calculated for symmetry no.    4


   SCF energy         :     -111.716485634508004
 -- inactive part     :     -142.398777877773796
 -- nuclear repulsion :       31.163673581965142


                     ***************************************
                     *** RHF response calculation (TDHF) ***
                     ***************************************



 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            3
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


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
 (dimension of paired reduced space:   24)
 RSP solution vector no.    1; norm of residual   1.76e-05
 RSP solution vector no.    2; norm of residual   3.23e-06
 RSP solution vector no.    3; norm of residual   2.57e-05

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    ZDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -0.69898419     *ENERGY(eV):   8.3380854    
@ STATE NO:    2 *TRANSITION MOMENT:  -1.0263431     *ENERGY(eV):   14.131640    
@ STATE NO:    3 *TRANSITION MOMENT:  -1.3884031     *ENERGY(eV):   25.555076    


  ******************************************************************************
  *** @ Excit. operator sym 1 & ref. state sym 1 => excited state symmetry 1 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.30641902    au
@                      8.3380854     eV
@                      67251.201     cm-1
                       804.50302     kJ / mol

@ Total energy :      -111.41007     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  9.98065779e-02  (Transition moment : -0.69898419     )


 @ Excited state no:    2 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.51932824    au
@                      14.131640     eV
@                      113979.37     cm-1
                       1363.4961     kJ / mol

@ Total energy :      -111.19716     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  0.36470001      (Transition moment :  -1.0263431     )


 @ Excited state no:    3 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.93913183    au
@                      25.555076     eV
@                      206115.61     cm-1
                       2465.6903     kJ / mol

@ Total energy :      -110.77735     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :   1.2068866      (Transition moment :  -1.3884031     )


 Time used in polarization propagator calculation is      0.23 CPU seconds for symmetry 1


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    2

 Number of excitations of this symmetry            3
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


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
 (dimension of paired reduced space:   22)
 RSP solution vector no.    1; norm of residual   2.14e-06
 RSP solution vector no.    2; norm of residual   1.04e-04
 RSP solution vector no.    3; norm of residual   7.57e-05

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    XDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  3.95340298e-02 *ENERGY(eV):   4.4417624    
@ STATE NO:    2 *TRANSITION MOMENT: -0.92567217     *ENERGY(eV):   13.538906    
@ STATE NO:    3 *TRANSITION MOMENT: -0.18004852     *ENERGY(eV):   23.387834    


  ******************************************************************************
  *** @ Excit. operator sym 2 & ref. state sym 1 => excited state symmetry 2 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.16323178    au
@                      4.4417624     eV
@                      35825.234     cm-1
                       428.56496     kJ / mol

@ Total energy :      -111.55325     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  1.70080927e-04  (Transition moment :  3.95340298e-02 )


 @ Excited state no:    2 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.49754568    au
@                      13.538906     eV
@                      109198.65     cm-1
                       1306.3060     kJ / mol

@ Total energy :      -111.21894     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  0.28422096      (Transition moment : -0.92567217     )


 @ Excited state no:    3 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.85948714    au
@                      23.387834     eV
@                      188635.62     cm-1
                       2256.5832     kJ / mol

@ Total energy :      -110.85700     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  1.85749322e-02  (Transition moment : -0.18004852     )


 Time used in polarization propagator calculation is      0.16 CPU seconds for symmetry 2


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    3

 Number of excitations of this symmetry            3
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


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
 RSP solution vector no.    1; norm of residual   4.44e-07
 RSP solution vector no.    2; norm of residual   4.82e-07
 RSP solution vector no.    3; norm of residual   1.11e-07

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    YDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -0.27554554     *ENERGY(eV):   11.490177    
@ STATE NO:    2 *TRANSITION MOMENT:  0.89229721     *ENERGY(eV):   17.982237    
@ STATE NO:    3 *TRANSITION MOMENT:  0.17704124     *ENERGY(eV):   21.730523    


  ******************************************************************************
  *** @ Excit. operator sym 3 & ref. state sym 1 => excited state symmetry 3 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.42225626    au
@                      11.490177     eV
@                      92674.538     cm-1
                       1108.6337     kJ / mol

@ Total energy :      -111.29423     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  2.13733020e-02  (Transition moment : -0.27554554     )


 @ Excited state no:    2 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.66083509    au
@                      17.982237     eV
@                      145036.54     cm-1
                       1735.0223     kJ / mol

@ Total energy :      -111.05565     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  0.35076876      (Transition moment :  0.89229721     )


 @ Excited state no:    3 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.79858207    au
@                      21.730523     eV
@                      175268.51     cm-1
                       2096.6769     kJ / mol

@ Total energy :      -110.91790     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  1.66869588e-02  (Transition moment :  0.17704124     )


 Time used in polarization propagator calculation is      0.12 CPU seconds for symmetry 3


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    4

 Number of excitations of this symmetry            3
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


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
 RSP solution vector no.    1; norm of residual   2.85e-08
 RSP solution vector no.    2; norm of residual   2.69e-08
 RSP solution vector no.    3; norm of residual   5.15e-09

 *** RSPCTL MICROITERATIONS CONVERGED


  ******************************************************************************
  *** @ Excit. operator sym 4 & ref. state sym 1 => excited state symmetry 4 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  4
 ---------------------------------------

@ Excitation energy :  5.09881416e-02au
@                      1.3874579     eV
@                      11190.604     cm-1
                       133.86935     kJ / mol

@ Total energy :      -111.66550     au


 @ Excited state no:    2 in symmetry  4
 ---------------------------------------

@ Excitation energy :  0.35472185    au
@                      9.6524722     eV
@                      77852.447     cm-1
                       931.32208     kJ / mol

@ Total energy :      -111.36176     au


 @ Excited state no:    3 in symmetry  4
 ---------------------------------------

@ Excitation energy :   1.0425328    au
@                      28.368760     eV
@                      228809.50     cm-1
                       2737.1695     kJ / mol

@ Total energy :      -110.67395     au


 Time used in polarization propagator calculation is      0.03 CPU seconds for symmetry 4

 >>>> Total CPU  time used in RESPONSE:   0.56 seconds
 >>>> Total wall time used in RESPONSE:   1.00 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   1.18 seconds
 >>>> Total wall time used in DALTON:   2.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:57:06 2009
     Host name              : stallo-2.local                          
END REFOUT

