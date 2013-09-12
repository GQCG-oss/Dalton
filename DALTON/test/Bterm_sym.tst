########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm dft water symmetry mcd bterm long
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enedft
tes
nuc
sym
cmass
omegab
Bterm
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN WAVE FUNCTION
.RUN RESPONSE
*PCM
.NPCMMT
0
.EPS
 24.55
.RSOLV
 2.180
*PCMCAV
.AREATS
0.5
**WAVE FUNCTIONS
.DFT
B3LYP
*SCF INPUT
.THRESHOLD
1.0D-8
**RESPONSE
*QUADRATIC
.SINGLE RESIDUE
.MCDBTERM
.ROOTS
 4 4 4 4
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
aug-cc-pVDZ
magnet
nuova base
    2    2 X Y
        8.    1
O      0.0000000000  0.0000000000  -0.1258515023
        1.    1
H      0.0000000000  1.4523500000   0.9986773907
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

     Date and time (Linux)  : Thu Sep 24 00:48:55 2009
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

     -----------------------------------
     Input for PCM solvation calculation 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=INPUT        EPS   = 24.5500     EPSINF=  0.0000
     RSOLV =  2.1800

     ICESPH =       0     NESFP =       0
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
 1: magnet                                                                  
 2: nuova base                                                              
    ------------------------------------------------------------------------

  Atomic type no.    1
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   8 :
     "/home/ruud/DaltonFix/dalton/basis/aug-cc-pVDZ"

  Atomic type no.    2
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   1 :
     "/home/ruud/DaltonFix/dalton/basis/aug-cc-pVDZ"


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
   1    0.0000000000e+00    0.0000000000e+00   -1.2585150230e-01    1.5000000000e+00
   2    0.0000000000e+00    1.4523500000e+00    9.9867739070e-01    1.2000000000e+00
   3    0.0000000000e+00   -1.4523500000e+00    9.9867739070e-01    1.2000000000e+00


                                 Isotopic Masses
                                 ---------------

                           O          15.994915
                           H    1      1.007825
                           H    2      1.007825

                       Total mass:    18.010565 amu
                       Natural abundance:  99.730 %

 Center-of-mass coordinates (a.u.):    0.000000    0.000000    0.000000


  Atoms and basis sets
  --------------------

  Number of atom types:     2
  Total number of atoms:    3

  Basis set used is "aug-cc-pVDZ" from the basis set library.

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  O           1    8.0000    35    23      [10s5p2d|4s3p2d]                                   
  H           2    1.0000    11     9      [5s2p|3s2p]                                        
  ----------------------------------------------------------------------
  total:      3   10.0000    57    41
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00e-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:    9
  O       :    1  x   0.0000000000   2  y   0.0000000000   3  z  -0.1258515023
  H   / 1 :    4  x   0.0000000000   5  y   1.4523500000   6  z   0.9986773907
  H   / 2 :    7  x   0.0000000000   8  y  -1.4523500000   9  z   0.9986773907


  Symmetry Coordinates
  --------------------

  Number of coordinates in each symmetry:    3   2   3   1

  Symmetry  A1  ( 1)

   1   O     z    3
   2   H     y    [  5  -   8 ]/2
   3   H     z    [  6  +   9 ]/2

  Symmetry  B1  ( 2)

   4   O     x    1
   5   H     x    [  4  +   7 ]/2

  Symmetry  B2  ( 3)

   6   O     y    2
   7   H     y    [  5  +   8 ]/2
   8   H     z    [  6  -   9 ]/2

  Symmetry  A2  ( 4)

   9   H     x    [  4  -   7 ]/2


   Interatomic separations (in Angstrom):
   --------------------------------------

            O           H    1      H    2
            ------      ------      ------
 O     :    0.000000
 H    1:    0.972000    0.000000
 H    2:    0.972000    1.537101    0.000000


  Max interatomic separation is    1.5371 Angstrom (    2.9047 Bohr)
  between atoms    3 and    2, "H    2" and "H    1".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  H    1     O            0.972000
  bond distance:  H    2     O            0.972000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     H    1     O          H    2       104.500




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       0.633889          0.000000    1.000000    0.000000
   IB       1.190584          0.000000    0.000000    1.000000
   IC       1.824473          1.000000    0.000000    0.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         797267.3440         424480.0013         277000.0234 MHz
           26.593976           14.159129            9.239726 cm-1


  Nuclear repulsion energy :    9.055004525638


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:              18   7  12   4


  Symmetry  A1 ( 1)

    1     O        1s         1
    2     O        1s         2
    3     O        1s         3
    4     O        1s         4
    5     O        2pz        7
    6     O        2pz       10
    7     O        2pz       13
    8     O        3d0       16
    9     O        3d2+      18
   10     O        3d0       21
   11     O        3d2+      23
   12     H        1s        24  +  25
   13     H        1s        26  +  27
   14     H        1s        28  +  29
   15     H        2py       32  -  33
   16     H        2pz       34  +  35
   17     H        2py       38  -  39
   18     H        2pz       40  +  41


  Symmetry  B1 ( 2)

   19     O        2px        5
   20     O        2px        8
   21     O        2px       11
   22     O        3d1+      17
   23     O        3d1+      22
   24     H        2px       30  +  31
   25     H        2px       36  +  37


  Symmetry  B2 ( 3)

   26     O        2py        6
   27     O        2py        9
   28     O        2py       12
   29     O        3d1-      15
   30     O        3d1-      20
   31     H        1s        24  -  25
   32     H        1s        26  -  27
   33     H        1s        28  -  29
   34     H        2py       32  +  33
   35     H        2pz       34  -  35
   36     H        2py       38  +  39
   37     H        2pz       40  -  41


  Symmetry  A2 ( 4)

   38     O        3d2-      14
   39     O        3d2-      19
   40     H        2px       30  -  31
   41     H        2px       36  -  37

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

 Number of two-electron integrals written:       91307 ( 24.6% )
 Megabytes written:                              1.051


 MEMORY USED TO GENERATE CAVITY =    432042


 Total number of spheres =    3
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1    0.000000000    0.000000000   -0.066597747    1.800000000   24.581233228
   2    0.000000000    0.768550518    0.528477314    1.440000000   11.987539425
   3    0.000000000   -0.768550518    0.528477314    1.440000000   11.987539425

 Total number of tesserae =     152
 Surface area =   48.55631208 (A^2)    Cavity volume =   30.53397368 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....
 
  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   0.06 seconds
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
     Number of Huckel orbitals each symmetry:    4    1    2    0

 Huckel EHT eigenvalues for symmetry :  1
          -20.808625      -2.112474      -0.679872       0.264008

 Huckel EHT eigenvalues for symmetry :  2
           -0.616200

 Huckel EHT eigenvalues for symmetry :  3
           -0.926303       0.099566

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Sep 24 00:48:55 2009
     Host name              : stallo-2.local                          

 Title lines from ".mol" input file:
     magnet                                                                  
     nuova base                                                              

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham DFT calculation.


     Time-dependent Kohn-Sham DFT calculation (TD-DFT).


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         10
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                4
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
     Abelian symmetry species          All    1    2    3    4
                                       --  --  --  --
     Total number of orbitals           41   18    7   12    4
     Number of basis functions          41   18    7   12    4

      ** Automatic occupation of RKS orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals               5    3    1    1    0

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE        9.97104
 NUCLEAR APPARENT CHARGE  -9.59233
 THEORETICAL  -9.59267 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  10 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Angular polynomials in range [15 35]
 Atom:    1*1 points=24126 compressed from 24126 (123 radial)
 Atom:    2*2 points=23238 compressed from 23238 ( 93 radial)
 Number of grid points:    47364 Grid generation time:       0.0 s
K-S electrons/energy :    9.99999963893949   -7.40267055255523 err:-.36e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.63485    26.61972    26.61963   -13.31225    -0.01000
   1  -76.0960993972     -1.000176258944e-02   3.47e+00  -7.61e+01    3  1  1  0
 MULPOP O     7.91; H     2.09; 
K-S electrons/energy :    9.99999967099829   -7.61229073815259 err:-.33e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.72603    26.66724    26.66722   -13.31225    -0.00804
   2  -76.3856340019     -8.037525120855e-03   8.09e-01  -2.90e-01    3  1  1  0
 MULPOP O     7.50; H     2.50; 
K-S electrons/energy :    9.99999967941053   -7.35507633275836 err:-.32e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.64088    26.62576    26.62570   -13.31225    -0.00696
   3  -76.3868224195     -6.961402581739e-03   7.94e-01  -1.19e-03    3  1  1  0
 MULPOP O     7.25; H     2.75; 
K-S electrons/energy :    9.99999967108955   -7.52243666474058 err:-.33e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.70508    26.65498    26.65493   -13.31225    -0.00983
   4  -76.4140555966     -9.834812675809e-03   1.67e-01  -2.72e-02    3  1  1  0
 MULPOP O     7.47; H     2.53; 
K-S electrons/energy :    9.99999967340928   -7.49475472362035 err:-.33e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.69366    26.65021    26.65016   -13.31225    -0.00890
   5  -76.4153542587     -8.895700470863e-03   8.31e-03  -1.30e-03    3  1  1  0
 MULPOP O     7.43; H     2.57; 
K-S electrons/energy :    9.99999967334486   -7.49578246708002 err:-.33e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.69394    26.65031    26.65026   -13.31225    -0.00894
   6  -76.4153573807     -8.935104995775e-03   7.31e-04  -3.12e-06    3  1  1  0
 MULPOP O     7.43; H     2.57; 
K-S electrons/energy :    9.99999967333831   -7.49592125334560 err:-.33e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.69399    26.65033    26.65028   -13.31225    -0.00894
   7  -76.4153574080     -8.936982209827e-03   3.85e-05  -2.73e-08    3  1  1  0
 MULPOP O     7.43; H     2.57; 
K-S electrons/energy :    9.99999967333860   -7.49591458368821 err:-.33e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -26.69398    26.65033    26.65028   -13.31225    -0.00894
   8  -76.4153574081     -8.936834411653e-03   1.75e-06  -1.05e-10    3  1  1  0
 DIIS converged in   8 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   10
 Orbital occupations :    3    1    1    0

 Sym       Kohn-Sham orbital energies

  1    -19.15770709    -1.01435574    -0.39866817    -0.01408856     0.10168961
         0.12975014     0.22734411     0.29590866

  2     -0.32191461     0.11216630     0.32696198     0.89445362     1.18997723
         1.84727001

  3     -0.52976365     0.03359656     0.16210904     0.17626517     0.45234852
         0.53245661

  4      0.28746055     0.95198398     1.75179831     3.31154196

    E(LUMO) :    -0.01408856 au (symmetry 1)
  - E(HOMO) :    -0.32191461 au (symmetry 2)
  ------------------------------------------
    gap     :     0.30782604 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   24.550000

     Final DFT energy:            -76.415357408098                 
     Nuclear repulsion:             9.055004525638
     Electronic energy:           -85.461425099324

     Final gradient norm:           0.000001754066

 
     Date and time (Linux)  : Thu Sep 24 00:48:59 2009
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5
   1 O   :1s     1.0039  -0.0038   0.0030  -0.0051  -0.0211
   2 O   :1s     0.0196   0.8880  -0.3045  -0.2480  -0.1142
   3 O   :1s    -0.0259  -0.1459  -0.1428  -0.1635   0.2592
   4 O   :1s    -0.0101   0.0662  -0.1249  -0.5856   2.5842
   5 O   :2pz    0.0024   0.1543   0.8088  -0.1633  -0.3669
   6 O   :2pz   -0.0041  -0.0863  -0.0713  -0.0312   0.0597
   7 O   :2pz   -0.0016   0.0136   0.0562  -0.1547   1.6590
   8 O   :3d0    0.0000   0.0034   0.0143  -0.0024  -0.0127
  10 O   :3d0    0.0001  -0.0060   0.0105  -0.0031   0.0360
  11 O   :3d2+   0.0007   0.0198   0.0149  -0.0039  -0.0294
  12 H   :1s     0.0014   0.3561   0.3556   0.1015   0.0414
  13 H   :1s     0.0075  -0.2059  -0.1080   0.2451  -1.0633
  14 H   :1s     0.0007  -0.0022   0.0076   0.5880  -0.5215
  15 H   :2py    0.0004  -0.0362  -0.0262   0.0164   0.0312
  16 H   :2pz    0.0004  -0.0209   0.0036   0.0093   0.0220
  17 H   :2py   -0.0027   0.0041  -0.0034  -0.0028   0.2239
  18 H   :2pz   -0.0024   0.0069   0.0138   0.0031   0.0424

 Molecular orbitals for symmetry species  2
 ------------------------------------------

 Orbital           1        2        3
   1 O   :2px    0.9200   0.3478   0.2547
   2 O   :2px   -0.0003   0.0854   0.4494
   3 O   :2px    0.1036  -1.3570   1.1962
   4 O   :3d1+   0.0137   0.0142  -0.0202
   5 O   :3d1+   0.0140  -0.0386   0.0186
   6 H   :2px    0.0261  -0.0079  -0.0104
   7 H   :2px    0.0211   0.2204  -1.1551

 Molecular orbitals for symmetry species  3
 ------------------------------------------

 Orbital           1        2        3
   1 O   :2py    0.7556  -0.2137  -0.4423
   2 O   :2py   -0.1959  -0.0672   0.2408
   3 O   :2py    0.0153  -0.4875   2.4030
   4 O   :3d1-   0.0249  -0.0010  -0.0005
   5 O   :3d1-  -0.0461  -0.0042   0.0758
   6 H   :1s     0.6088   0.0232   0.0455
   7 H   :1s    -0.1713   0.4128  -1.5287
   8 H   :1s     0.0018   2.4695  -2.8936
   9 H   :2py   -0.0214   0.0084   0.0391
  10 H   :2pz   -0.0322   0.0094   0.0277
  11 H   :2py   -0.0232  -0.0485   0.4593
  12 H   :2pz   -0.0159  -0.0661   0.4761

 Molecular orbitals for symmetry species  4
 ------------------------------------------

 Orbital           1        2
   1 O   :3d2-  -0.0339   0.0015
   2 O   :3d2-   0.0636  -1.2230
   3 H   :2px   -0.0123  -0.0394
   4 H   :2px   -1.0827   0.7806



 >>>> Total CPU  time used in SIRIUS :      3.92 seconds
 >>>> Total wall time used in SIRIUS :      4.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:48:59 2009
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

  1 B-frequencies  0.000000e+00

 Spin of operator A , ISPINA=    0
 Spin of operator B , ISPINB=    0
 Spin of operator C , (Excitation energy) ISPINC=    0
 B of Magnetic Circular Dichroism requested    : MCDCAL =T

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

    1 B OPERATORS OF SYMMETRY NO:    2 AND LABELS:

          YANGMOM 

    1 B OPERATORS OF SYMMETRY NO:    3 AND LABELS:

          XANGMOM 

    1 B OPERATORS OF SYMMETRY NO:    4 AND LABELS:

          ZANGMOM 


   SCF energy         :      -76.415357408098018
 -- inactive part     :      -85.461425099324259
 -- nuclear repulsion :        9.055004525637894


                    *****************************************
                    *** DFT response calculation (TD-DFT) ***
                    *****************************************



 Linear response excitations for quadratic response
 - symmetry of excitation operator    1
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      62
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      62



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       1.0 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       1.0 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*3 evaluation time:       0.8 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   30)
 RSP solution vector no.    1; norm of residual   2.56e-04
 RSP solution vector no.    2; norm of residual   3.94e-04
 RSP solution vector no.    3; norm of residual   4.10e-04
 RSP solution vector no.    4; norm of residual   1.87e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    2
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      37
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      37



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   24)
 RSP solution vector no.    1; norm of residual   1.21e-04
 RSP solution vector no.    2; norm of residual   3.61e-04
 RSP solution vector no.    3; norm of residual   6.31e-04
 RSP solution vector no.    4; norm of residual   7.95e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    3
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      52
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      52



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  3; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       1.0 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*3 evaluation time:       0.8 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   30)
 RSP solution vector no.    1; norm of residual   9.34e-05
 RSP solution vector no.    2; norm of residual   5.75e-05
 RSP solution vector no.    3; norm of residual   9.52e-05
 RSP solution vector no.    4; norm of residual   1.98e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response excitations for quadratic response
 - symmetry of excitation operator    4
 - is the operator a triplet operator ?     F


 Perturbation symmetry.     KSYMOP:       4
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      29
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      29



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  4; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       0.9 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   26)
 RSP solution vector no.    1; norm of residual   1.89e-05
 RSP solution vector no.    2; norm of residual   3.28e-05
 RSP solution vector no.    3; norm of residual   4.14e-05
 RSP solution vector no.    4; norm of residual   1.08e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    1


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      62
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      62


 QRLRVE -- linear response calculation for symmetry  1
 QRLRVE -- operator label : ZDIPLEN 
 QRLRVE -- frequencies :  0.260094  0.371727  0.408196  0.458800  0.394814
                          0.468683  0.512613  0.522952  0.316501  0.427812
                          0.443797  0.591275



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.9 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.37173 a.u.
 after   12 linear transformations is     -3.53551236
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -6.11 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.40820 a.u.
 after   12 linear transformations is      5.11862173
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  7.11 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.45880 a.u.
 after   12 linear transformations is     18.43477049
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  3.00 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39481 a.u.
 after   12 linear transformations is      4.50456069
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.05 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46868 a.u.
 after   12 linear transformations is     27.66433151
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  7.98 * 10 ** -3.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51261 a.u.
 after   12 linear transformations is      2.45729548
 INERTIA (POS,ZER,NEG) of reduced matrix is   20    0    4
 Determinant of reduced matrix is  1.73 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.52295 a.u.
 after   12 linear transformations is      6.62655242
 INERTIA (POS,ZER,NEG) of reduced matrix is   20    0    4
 Determinant of reduced matrix is  1.52 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   12 linear transformations is      9.27848109
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  1.32 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   12 linear transformations is     12.77957934
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  8.38 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   12 linear transformations is     -9.17991796
 INERTIA (POS,ZER,NEG) of reduced matrix is   18    0    6
 Determinant of reduced matrix is  1.68 * 10 ** -2.0
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.9 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.37173 a.u.
 after   24 linear transformations is     -3.14691187
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -1.79 * 10 **  9.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.40820 a.u.
 after   24 linear transformations is      6.04562534
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  2.72 * 10 **  8.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.45880 a.u.
 after   24 linear transformations is     19.26486154
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  2.43 * 10 **  7.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39481 a.u.
 after   24 linear transformations is      5.79830790
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -1.00 * 10 **  8.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46868 a.u.
 after   24 linear transformations is   -201.08418507
 INERTIA (POS,ZER,NEG) of reduced matrix is   45    0    3
 Determinant of reduced matrix is -1.40 * 10 **  5.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51261 a.u.
 after   24 linear transformations is      5.73791618
 INERTIA (POS,ZER,NEG) of reduced matrix is   44    0    4
 Determinant of reduced matrix is  1.37 * 10 **  7.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.52295 a.u.
 after   24 linear transformations is      9.47743353
 INERTIA (POS,ZER,NEG) of reduced matrix is   44    0    4
 Determinant of reduced matrix is  8.50 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   24 linear transformations is      9.69218246
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  2.68 * 10 **  8.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   24 linear transformations is     13.11927462
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  1.19 * 10 **  8.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   24 linear transformations is     14.81267388
 INERTIA (POS,ZER,NEG) of reduced matrix is   42    0    6
 Determinant of reduced matrix is  9.89 * 10 **  5.0
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.9 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.37173 a.u.
 after   36 linear transformations is     -3.14661563
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -1.56 * 10 ** 23.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.40820 a.u.
 after   36 linear transformations is      6.04690199
 INERTIA (POS,ZER,NEG) of reduced matrix is   70    0    2
 Determinant of reduced matrix is  2.03 * 10 ** 22.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.45880 a.u.
 after   36 linear transformations is     19.26603741
 INERTIA (POS,ZER,NEG) of reduced matrix is   70    0    2
 Determinant of reduced matrix is  1.43 * 10 ** 21.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39481 a.u.
 after   36 linear transformations is      5.80594707
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -7.89 * 10 ** 21.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46868 a.u.
 after   36 linear transformations is   -170.25487144
 INERTIA (POS,ZER,NEG) of reduced matrix is   69    0    3
 Determinant of reduced matrix is -9.06 * 10 ** 18.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51261 a.u.
 after   36 linear transformations is      5.74186450
 INERTIA (POS,ZER,NEG) of reduced matrix is   68    0    4
 Determinant of reduced matrix is  5.98 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.52295 a.u.
 after   36 linear transformations is      9.48142719
 INERTIA (POS,ZER,NEG) of reduced matrix is   68    0    4
 Determinant of reduced matrix is  3.47 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   36 linear transformations is      9.69271483
 INERTIA (POS,ZER,NEG) of reduced matrix is   70    0    2
 Determinant of reduced matrix is  1.83 * 10 ** 22.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   36 linear transformations is     13.11955629
 INERTIA (POS,ZER,NEG) of reduced matrix is   70    0    2
 Determinant of reduced matrix is  7.56 * 10 ** 21.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   36 linear transformations is     15.07430907
 INERTIA (POS,ZER,NEG) of reduced matrix is   66    0    6
 Determinant of reduced matrix is  2.45 * 10 ** 19.0

 *** THE REQUESTED   12 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   72)
 RSP solution vector no.    1; norm of residual   2.53e-04
 RSP solution vector no.    2; norm of residual   3.54e-05
 RSP solution vector no.    3; norm of residual   1.87e-04
 RSP solution vector no.    4; norm of residual   6.19e-05
 RSP solution vector no.    5; norm of residual   1.54e-04
 RSP solution vector no.    6; norm of residual   2.86e-04
 RSP solution vector no.    7; norm of residual   2.22e-04
 RSP solution vector no.    8; norm of residual   2.47e-04
 RSP solution vector no.    9; norm of residual   1.33e-04
 RSP solution vector no.   10; norm of residual   1.41e-04
 RSP solution vector no.   11; norm of residual   8.76e-05
 RSP solution vector no.   12; norm of residual   3.80e-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.260094e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.26009):     17.0603589390    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.371727e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.37173):    -3.14661562509    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.408196e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.40820):     6.04690198651    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.458800e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.45880):     19.2660374134    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.394814e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.39481):     5.80594706568    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.468683e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.46868):    -170.254871441    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.512613e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.51261):     5.74186449844    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.522952e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.52295):     9.48142719298    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.316501e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.31650):     37.5548954587    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.427812e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.42781):     9.69271482675    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.443797e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.44380):     13.1195562918    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.591275e+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.59128):     15.0743090668    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    2


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      37
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      37


 QRLRVE -- linear response calculation for symmetry  2
 QRLRVE -- operator label : XDIPLEN 
 QRLRVE -- frequencies :  0.397227  0.468547  0.483278  0.394814  0.468683
                          0.512613  0.522952  0.316501  0.427812  0.443797
                          0.591275  0.336845



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.7 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39723 a.u.
 after   12 linear transformations is     12.18015370
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  1.68 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46855 a.u.
 after   12 linear transformations is     22.87359688
 INERTIA (POS,ZER,NEG) of reduced matrix is   20    0    4
 Determinant of reduced matrix is  8.97 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.48328 a.u.
 after   12 linear transformations is     60.08590228
 INERTIA (POS,ZER,NEG) of reduced matrix is   20    0    4
 Determinant of reduced matrix is  9.34 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39481 a.u.
 after   12 linear transformations is     11.62899915
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  1.94 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46868 a.u.
 after   12 linear transformations is     23.02764840
 INERTIA (POS,ZER,NEG) of reduced matrix is   20    0    4
 Determinant of reduced matrix is  9.06 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51261 a.u.
 after   12 linear transformations is     -5.36451809
 INERTIA (POS,ZER,NEG) of reduced matrix is   19    0    5
 Determinant of reduced matrix is -6.38 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.52295 a.u.
 after   12 linear transformations is      2.32025124
 INERTIA (POS,ZER,NEG) of reduced matrix is   19    0    5
 Determinant of reduced matrix is -9.82 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.31650 a.u.
 after   12 linear transformations is      2.91355714
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -6.84 * 10 **  0.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   12 linear transformations is     12.60474913
 INERTIA (POS,ZER,NEG) of reduced matrix is   21    0    3
 Determinant of reduced matrix is -2.02 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   12 linear transformations is     16.72698372
 INERTIA (POS,ZER,NEG) of reduced matrix is   21    0    3
 Determinant of reduced matrix is -1.64 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   12 linear transformations is     -0.92800918
 INERTIA (POS,ZER,NEG) of reduced matrix is   18    0    6
 Determinant of reduced matrix is  9.13 * 10 **  0.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.33685 a.u.
 after   12 linear transformations is      5.50725374
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -3.38 * 10 **  0.0
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.7 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39723 a.u.
 after   24 linear transformations is     12.26267463
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  8.78 * 10 ** 11.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46855 a.u.
 after   24 linear transformations is     24.03918848
 INERTIA (POS,ZER,NEG) of reduced matrix is   44    0    4
 Determinant of reduced matrix is  3.47 * 10 ** 11.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.48328 a.u.
 after   24 linear transformations is     64.79431529
 INERTIA (POS,ZER,NEG) of reduced matrix is   44    0    4
 Determinant of reduced matrix is  3.08 * 10 ** 11.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39481 a.u.
 after   24 linear transformations is     11.71762197
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  1.04 * 10 ** 12.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46868 a.u.
 after   24 linear transformations is     24.19538981
 INERTIA (POS,ZER,NEG) of reduced matrix is   44    0    4
 Determinant of reduced matrix is  3.50 * 10 ** 11.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51261 a.u.
 after   24 linear transformations is     -5.14129615
 INERTIA (POS,ZER,NEG) of reduced matrix is   43    0    5
 Determinant of reduced matrix is -1.83 * 10 ** 12.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.52295 a.u.
 after   24 linear transformations is      2.55926858
 INERTIA (POS,ZER,NEG) of reduced matrix is   43    0    5
 Determinant of reduced matrix is -2.59 * 10 ** 12.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.31650 a.u.
 after   24 linear transformations is      3.12676429
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -5.21 * 10 ** 13.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   24 linear transformations is     13.06893718
 INERTIA (POS,ZER,NEG) of reduced matrix is   45    0    3
 Determinant of reduced matrix is -9.72 * 10 ** 11.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   24 linear transformations is     17.18477509
 INERTIA (POS,ZER,NEG) of reduced matrix is   45    0    3
 Determinant of reduced matrix is -7.00 * 10 ** 11.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   24 linear transformations is     -0.73128703
 INERTIA (POS,ZER,NEG) of reduced matrix is   42    0    6
 Determinant of reduced matrix is  1.39 * 10 ** 13.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.33685 a.u.
 after   24 linear transformations is      5.67951091
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -2.36 * 10 ** 13.0
 Electrons: 10.000000(-3.27e-07): LR-DFT*5 evaluation time:       1.1 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39723 a.u.
 after   29 linear transformations is     12.26268075
 INERTIA (POS,ZER,NEG) of reduced matrix is   56    0    2
 Determinant of reduced matrix is  1.37 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46855 a.u.
 after   29 linear transformations is     24.03925885
 INERTIA (POS,ZER,NEG) of reduced matrix is   54    0    4
 Determinant of reduced matrix is  5.16 * 10 ** 19.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.48328 a.u.
 after   29 linear transformations is     64.79455793
 INERTIA (POS,ZER,NEG) of reduced matrix is   54    0    4
 Determinant of reduced matrix is  4.53 * 10 ** 19.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39481 a.u.
 after   29 linear transformations is     11.71762663
 INERTIA (POS,ZER,NEG) of reduced matrix is   56    0    2
 Determinant of reduced matrix is  1.63 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46868 a.u.
 after   29 linear transformations is     24.19545929
 INERTIA (POS,ZER,NEG) of reduced matrix is   54    0    4
 Determinant of reduced matrix is  5.21 * 10 ** 19.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.51261 a.u.
 after   29 linear transformations is     -5.14123799
 INERTIA (POS,ZER,NEG) of reduced matrix is   53    0    5
 Determinant of reduced matrix is -2.64 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.52295 a.u.
 after   29 linear transformations is      2.55931274
 INERTIA (POS,ZER,NEG) of reduced matrix is   53    0    5
 Determinant of reduced matrix is -3.69 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.31650 a.u.
 after   29 linear transformations is      3.12677674
 INERTIA (POS,ZER,NEG) of reduced matrix is   57    0    1
 Determinant of reduced matrix is -8.48 * 10 ** 21.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   29 linear transformations is     13.06895865
 INERTIA (POS,ZER,NEG) of reduced matrix is   55    0    3
 Determinant of reduced matrix is -1.49 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   29 linear transformations is     17.18480237
 INERTIA (POS,ZER,NEG) of reduced matrix is   55    0    3
 Determinant of reduced matrix is -1.06 * 10 ** 20.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   29 linear transformations is     -0.73119306
 INERTIA (POS,ZER,NEG) of reduced matrix is   52    0    6
 Determinant of reduced matrix is  1.86 * 10 ** 21.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.33685 a.u.
 after   29 linear transformations is      5.67951957
 INERTIA (POS,ZER,NEG) of reduced matrix is   57    0    1
 Determinant of reduced matrix is -3.81 * 10 ** 21.0

 *** THE REQUESTED   12 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   58)
 RSP solution vector no.    1; norm of residual   2.51e-05
 RSP solution vector no.    2; norm of residual   1.03e-05
 RSP solution vector no.    3; norm of residual   1.08e-05
 RSP solution vector no.    4; norm of residual   2.64e-05
 RSP solution vector no.    5; norm of residual   1.03e-05
 RSP solution vector no.    6; norm of residual   2.19e-05
 RSP solution vector no.    7; norm of residual   3.19e-05
 RSP solution vector no.    8; norm of residual   5.74e-05
 RSP solution vector no.    9; norm of residual   1.21e-05
 RSP solution vector no.   10; norm of residual   1.31e-05
 RSP solution vector no.   11; norm of residual   1.15e-05
 RSP solution vector no.   12; norm of residual   4.96e-05

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.397227e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.39723):     12.2626807502    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.39723):    0.434461840477    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.468547e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.46855):     24.0392588451    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.46855):    -3.38537170660    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.483278e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.48328):     64.7945579294    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.48328):    -10.6058501311    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.394814e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.39481):     11.7176266254    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.39481):    0.212451770075    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.468683e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.46868):     24.1954592927    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.46868):    -3.40754799971    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.512613e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.51261):    -5.14123799032    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.51261):     3.26376460752    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.522952e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.52295):     2.55931274472    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.52295):     2.13930485467    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.316501e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.31650):     3.12677673557    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.31650):   -0.914038044129    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.427812e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.42781):     13.0689586543    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.42781):    -1.89613429551    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.443797e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.44380):     17.1848023662    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.44380):    -1.97372918663    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.591275e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.59128):   -0.731193064649    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.59128):     2.57407637381    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.336845e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.33685):     5.67951957286    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.33685):   -0.683355086844    


 QRLRVE -- linear response calculation for symmetry  2
 QRLRVE -- operator label : YANGMOM 
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   5.03e-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   YANGMOM     FREQUENCY   0.000000e+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; YANGMOM  >> (   0.00000):   -9.506284648353e-16

@QRLRVE:  << YANGMOM  ; YANGMOM  >> (   0.00000):    0.563437156969    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    3


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      52
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      52


 QRLRVE -- linear response calculation for symmetry  3
 QRLRVE -- operator label : YDIPLEN 
 QRLRVE -- frequencies :  0.408196  0.458800  0.316501  0.427812  0.443797
                          0.591275  0.336845  0.397227  0.468547  0.483278
                          0.260094  0.371727



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.8 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.40820 a.u.
 after   12 linear transformations is     21.75080902
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.58 * 10 **  0.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.45880 a.u.
 after   12 linear transformations is     76.17918855
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.46 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   12 linear transformations is     30.36267916
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -1.31 * 10 **  0.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   12 linear transformations is     41.34561874
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -5.98 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   12 linear transformations is     -9.65622034
 INERTIA (POS,ZER,NEG) of reduced matrix is   18    0    6
 Determinant of reduced matrix is  2.50 * 10 ** -2.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39723 a.u.
 after   12 linear transformations is     -2.35043901
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -4.29 * 10 ** -1.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46855 a.u.
 after   12 linear transformations is   1001.39927598
 INERTIA (POS,ZER,NEG) of reduced matrix is   23    0    1
 Determinant of reduced matrix is -5.09 * 10 ** -3.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.48328 a.u.
 after   12 linear transformations is    -10.81260143
 INERTIA (POS,ZER,NEG) of reduced matrix is   22    0    2
 Determinant of reduced matrix is  4.85 * 10 ** -2.0
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       2.8 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.40820 a.u.
 after   24 linear transformations is     22.31204708
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -1.45 * 10 ** 10.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.45880 a.u.
 after   24 linear transformations is     79.48783261
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -7.97 * 10 **  8.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   24 linear transformations is     31.08838993
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -1.00 * 10 ** 10.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   24 linear transformations is     42.44842245
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -3.90 * 10 **  9.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   24 linear transformations is     -7.36994794
 INERTIA (POS,ZER,NEG) of reduced matrix is   42    0    6
 Determinant of reduced matrix is  3.80 * 10 **  7.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39723 a.u.
 after   24 linear transformations is      0.71298817
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -4.75 * 10 **  9.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46855 a.u.
 after   24 linear transformations is   4165.71939995
 INERTIA (POS,ZER,NEG) of reduced matrix is   47    0    1
 Determinant of reduced matrix is -6.19 * 10 **  6.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.48328 a.u.
 after   24 linear transformations is     -8.02560404
 INERTIA (POS,ZER,NEG) of reduced matrix is   46    0    2
 Determinant of reduced matrix is  2.17 * 10 **  8.0
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       3.0 s

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.40820 a.u.
 after   36 linear transformations is     22.31250655
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -1.04 * 10 ** 26.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.45880 a.u.
 after   36 linear transformations is     79.48993608
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -4.67 * 10 ** 24.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.42781 a.u.
 after   36 linear transformations is     31.08885372
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -6.67 * 10 ** 25.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.44380 a.u.
 after   36 linear transformations is     42.44908423
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -2.44 * 10 ** 25.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.59128 a.u.
 after   36 linear transformations is     -7.36809836
 INERTIA (POS,ZER,NEG) of reduced matrix is   66    0    6
 Determinant of reduced matrix is  1.10 * 10 ** 23.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.39723 a.u.
 after   36 linear transformations is      0.71576050
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -3.54 * 10 ** 25.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.46855 a.u.
 after   36 linear transformations is   4175.49780614
 INERTIA (POS,ZER,NEG) of reduced matrix is   71    0    1
 Determinant of reduced matrix is -3.47 * 10 ** 22.0

 *** INFO, negative eigenvalues in reduced matrix.

 GP * SOLUTION vector at frequency     0.48328 a.u.
 after   36 linear transformations is     -8.02340689
 INERTIA (POS,ZER,NEG) of reduced matrix is   70    0    2
 Determinant of reduced matrix is  1.14 * 10 ** 24.0

 *** THE REQUESTED   12 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   72)
 RSP solution vector no.    1; norm of residual   2.34e-05
 RSP solution vector no.    2; norm of residual   1.78e-05
 RSP solution vector no.    3; norm of residual   5.21e-05
 RSP solution vector no.    4; norm of residual   2.98e-05
 RSP solution vector no.    5; norm of residual   2.49e-05
 RSP solution vector no.    6; norm of residual   1.78e-05
 RSP solution vector no.    7; norm of residual   4.87e-05
 RSP solution vector no.    8; norm of residual   8.65e-06
 RSP solution vector no.    9; norm of residual   1.54e-05
 RSP solution vector no.   10; norm of residual   2.57e-05
 RSP solution vector no.   11; norm of residual   5.86e-05
 RSP solution vector no.   12; norm of residual   3.69e-05

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.408196e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.40820):     22.3125065472    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.40820):   -0.340685548752    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.458800e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.45880):     79.4899360843    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.45880):    -3.76735417522    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.316501e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.31650):     17.4395631136    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.31650):    0.201009689540    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.427812e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.42781):     31.0888537197    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.42781):   -0.120813070845    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.443797e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.44380):     42.4490842270    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.44380):   -0.573373230039    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.591275e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.59128):    -7.36809836475    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.59128):    -8.96664240621    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.336845e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.33685):     18.8827776445    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.33685):    0.263554944453    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.397227e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.39723):    0.715760495191    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.39723):    -2.91006725932    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.468547e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.46855):     4175.49780614    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.46855):    -473.267733042    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.483278e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.48328):    -8.02340689484    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.48328):     9.57757121008    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.260094e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.26009):     14.9197116682    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.26009):    0.105490842914    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.371727e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.37173):     23.3126962126    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.37173):    0.520316337198    


 QRLRVE -- linear response calculation for symmetry  3
 QRLRVE -- operator label : XANGMOM 
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   1.74e-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   XANGMOM     FREQUENCY   0.000000e+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; XANGMOM  >> (   0.00000):   -1.902991653147e-15

@QRLRVE:  << XANGMOM  ; XANGMOM  >> (   0.00000):     2.07779628658    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    4


 Perturbation symmetry.     KSYMOP:       4
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      29
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      29


 QRLRVE -- linear response calculation for symmetry  4
 QRLRVE -- operator label : ZANGMOM 
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  4; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.5 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       0.4 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   1.44e-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZANGMOM     FREQUENCY   0.000000e+00
 SYMMETRY    4

@QRLRVE:  << ZANGMOM  ; ZANGMOM  >> (   0.00000):     1.14494117422    
 DFT-QR computed in a linearly-scaling fashion.

 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.000000    0.336845   -8.381611

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0
 Excitation energy in au,    moment :                0.336845   -0.737144

 B term contribution:              3.089226
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.000000    0.397227   -7.166815

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0
 Excitation energy in au,    moment :                0.397227   -0.067943

 B term contribution:              0.243469
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.000000    0.468547   -9.871321

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0
 Excitation energy in au,    moment :                0.468547   -0.162248

 B term contribution:              0.800800
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             4    1    0

 omega B, excitation energy, moment :    0.000000    0.483278   -6.976231

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             4    1    0
 Excitation energy in au,    moment :                0.483278   -0.457884

 B term contribution:              1.597151
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.000000    0.336845   -5.468477

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0
 Excitation energy in au,    moment :                0.336845   -0.737144

 B term contribution:              2.015527
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.000000    0.397227   10.999831

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0
 Excitation energy in au,    moment :                0.397227   -0.067943

 B term contribution:             -0.373683
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.000000    0.468547  449.465868

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0
 Excitation energy in au,    moment :                0.468547   -0.162248

 B term contribution:            -36.462432
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             4    1    0

 omega B, excitation energy, moment :    0.000000    0.483278   -6.361964

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             4    1    0
 Excitation energy in au,    moment :                0.483278   -0.457884

 B term contribution:              1.456520
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.000000    0.260094   -9.238084

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0
 Excitation energy in au,    moment :                0.260094   -0.651687

 B term contribution:              3.010169
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.000000    0.371727   -1.045782

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0
 Excitation energy in au,    moment :                0.371727    0.002162

 B term contribution:             -0.001130
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.000000    0.408196    6.718619

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0
 Excitation energy in au,    moment :                0.408196    0.141087

 B term contribution:              0.473954
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             4    2    0

 omega B, excitation energy, moment :    0.000000    0.458800   20.554772

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             4    2    0
 Excitation energy in au,    moment :                0.458800    0.157024

 B term contribution:              1.613796
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.000000    0.260094    2.756298

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0
 Excitation energy in au,    moment :                0.260094   -0.651687

 B term contribution:             -0.898122
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.000000    0.371727    1.412016

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0
 Excitation energy in au,    moment :                0.371727    0.002162

 B term contribution:              0.001526
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.000000    0.408196   -0.155426

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0
 Excitation energy in au,    moment :                0.408196    0.141087

 B term contribution:             -0.010964
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             4    2    0

 omega B, excitation energy, moment :    0.000000    0.458800   -7.028351

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             4    2    0
 Excitation energy in au,    moment :                0.458800    0.157024

 B term contribution:             -0.551810
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.000000    0.394814   -3.141231

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0
 Excitation energy in au,    moment :                0.394814   -0.238922

 B term contribution:              0.375254
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.000000    0.468683  -92.054970

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0
 Excitation energy in au,    moment :                0.468683    0.748664

 B term contribution:            -34.459123
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.000000    0.512613    6.488796

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0
 Excitation energy in au,    moment :                0.512613   -0.203752

 B term contribution:             -0.661051
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             4    3    0

 omega B, excitation energy, moment :    0.000000    0.522952   10.906975

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             4    3    0
 Excitation energy in au,    moment :                0.522952   -0.584103

 B term contribution:             -3.185397
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.000000    0.394814    1.661208

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0
 Excitation energy in au,    moment :                0.394814   -0.238922

 B term contribution:             -0.198449
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.000000    0.468683   -2.145464

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0
 Excitation energy in au,    moment :                0.468683    0.748664

 B term contribution:             -0.803116
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.000000    0.512613   -6.314435

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0
 Excitation energy in au,    moment :                0.512613   -0.203752

 B term contribution:              0.643288
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             4    3    0

 omega B, excitation energy, moment :    0.000000    0.522952   -5.720554

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             4    3    0
 Excitation energy in au,    moment :                0.522952   -0.584103

 B term contribution:              1.670696
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.000000    0.316501    2.981112
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.000000    0.427812    2.398043
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.000000    0.443797   -0.910141
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             4    4    0

 omega B, excitation energy, moment :    0.000000    0.591275  -28.982584
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.000000    0.316501    2.274182
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.000000    0.427812   -4.844422
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.000000    0.443797   10.412323
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             4    4    0

 omega B, excitation energy, moment :    0.000000    0.591275   -1.273233
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.000000    0.316501   -4.255227
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.000000    0.427812    1.014522
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.5 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.000000    0.443797   -0.417551
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       0.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             4    4    0

 omega B, excitation energy, moment :    0.000000    0.591275    7.985711
 No MCD for this SOMOM

 >>>> Total CPU  time used in RESPONSE:  1 minute  35 seconds
 >>>> Total wall time used in RESPONSE:  1 minute  47 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  1 minute  39 seconds
 >>>> Total wall time used in DALTON:  1 minute  52 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:50:47 2009
     Host name              : stallo-2.local                          
END REFOUT

