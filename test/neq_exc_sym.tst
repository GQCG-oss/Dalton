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

     Date and time (Linux)  : Fri Dec 10 15:03:23 2010
     Host name              : stallo-2.local                          

 * Work memory size             :   100000000 =  762.94 megabytes.

 * Directories for basis set searches:
   1) /home/luca/programs/dalton/trunk/test/perl-pid.24302__2010_12_10__15.3
   2) /home/luca/programs/dalton/trunk/basis/


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

     ICESPH =       2     NESFP =       4
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
  Used basis set file for basis set for elements with Z =   6 :
     "/home/luca/programs/dalton/trunk/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   8 :
     "/home/luca/programs/dalton/trunk/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   1 :
     "/home/luca/programs/dalton/trunk/basis/STO-3G"


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
  O       :     4  x   0.0000000000    5  y   0.0000000000    6  z   2.3054658834
  H   / 1 :     7  x   1.7822044964    8  y   0.0000000000    9  z  -1.0289558799
  H   / 2 :    10  x  -1.7822044964   11  y   0.0000000000   12  z  -1.0289558799


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


@  Nuclear repulsion energy :   31.163673581965


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

 The following one-electron property integrals are calculated:
          - overlap integrals
          - dipole length integrals

 Center of mass  (bohr):      0.000000000000      0.000000000000      1.159648807704
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014


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
 >>>> Total CPU  time used in HERMIT:   0.08 seconds
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
          -20.704416     -11.377295      -1.495729      -0.964135      -0.593338
           -0.252244      -0.199950

 Huckel EHT eigenvalues for symmetry :  2
           -0.744095      -0.506268      -0.190531

 Huckel EHT eigenvalues for symmetry :  3
           -0.670977      -0.352123

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Fri Dec 10 15:03:23 2010
     Host name              : stallo-2.local                          

 Title lines from ".mol" input file:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham DFT calculation.


     Time-dependent Kohn-Sham DFT calculation (TD-DFT).


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                4
     Reference state symmetry                  1
 
     This is a DFT calculation of type: LDA

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1    2    3    4
                                       --- |  ---  ---  ---  ---
     Total number of orbitals           12 |    7    3    2    0
     Number of basis functions          12 |    7    3    2    0

      ** Automatic occupation of RKS orbitals **

      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is :
     Occupied SCF orbitals               8 |    5    2    1    0

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

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  16 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Angular polynomials in range [15 35]
 Atom:    1*1 points=18676 compressed from 18676 ( 92 radial)
 Atom:    2*1 points=18280 compressed from 18280 ( 92 radial)
 Atom:    3*2 points=18150 compressed from 18150 ( 75 radial)
 Number of grid points:    55106 Grid generation time:       0.0 s
K-S energy, electrons, error :    -13.964563268538  16.0000008466    8.47e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.68972    60.58642    60.58632   -30.26166    -0.02015
   1  -111.730265027     -2.014834300902e-02   2.59e+00  -1.12e+02    5  2  1  0
      MULPOP C    23.87; O     5.51; H     4.28; 
K-S energy, electrons, error :    -13.683589306844  16.0000010604    1.06e-06
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.27464    60.35152    60.35164   -30.26166    -0.04740
   2  -111.069683961     -4.739886941358e-02   3.47e+00   6.61e-01    5  2  1  0
      MULPOP C    30.44; O     5.92; H     3.78; 
K-S energy, electrons, error :    -13.972447264461  16.0000009048    9.05e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.78057    60.63369    60.63363   -30.26166    -0.01828
   3  -111.918604624     -1.827882230703e-02   1.17e+00  -8.49e-01    5  2  1  0
      MULPOP C    23.43; O     5.39; H     4.35; 
K-S energy, electrons, error :    -13.813173325637  16.0000009343    9.34e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.55603    60.53861    60.53859   -30.26166    -0.00107
   4  -112.031958946     -1.071362887835e-03   1.49e-01  -1.13e-01    5  2  1  0
      MULPOP C    25.38; O     5.57; H     4.19; 
K-S energy, electrons, error :    -13.837346284078  16.0000009364    9.36e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59447    60.55675    60.55673   -30.26166    -0.00215
   5  -112.033866468     -2.150639284714e-03   2.04e-02  -1.91e-03    5  2  1  0
      MULPOP C    25.23; O     5.54; H     4.21; 
K-S energy, electrons, error :    -13.834172652964  16.0000009352    9.35e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.58891    60.55410    60.55408   -30.26166    -0.00202
   6  -112.033907918     -2.019985311012e-03   2.67e-03  -4.15e-05    5  2  1  0
      MULPOP C    25.24; O     5.55; H     4.20; 
K-S energy, electrons, error :    -13.834464062208  16.0000009352    9.35e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.58937    60.55432    60.55430   -30.26166    -0.00204
   7  -112.033908473     -2.037037342962e-03   4.48e-04  -5.54e-07    5  2  1  0
      MULPOP C    25.23; O     5.55; H     4.20; 
K-S energy, electrons, error :    -13.834501696302  16.0000009351    9.35e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.58943    60.55434    60.55432   -30.26166    -0.00204
   8  -112.033908488     -2.039508846991e-03   3.16e-05  -1.50e-08    5  2  1  0
      MULPOP C    25.23; O     5.55; H     4.20; 
K-S energy, electrons, error :    -13.834504384336  16.0000009351    9.35e-07
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.58943    60.55434    60.55432   -30.26166    -0.00204
   9  -112.033908488     -2.039687623228e-03   1.13e-06  -7.13e-11    5  2  1  0
 DIIS converged in   9 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    5    2    1    0

 Sym       Kohn-Sham orbital energies

  1    -18.36962419    -9.66197370    -0.90481816    -0.50692414    -0.27783877
         0.36676366     0.53279077

  2     -0.38754497    -0.09325530     0.46802257

  3     -0.27474219     0.02488546

    E(LUMO) :     0.02488546 au (symmetry 3)
  - E(HOMO) :    -0.09325530 au (symmetry 2)
  ------------------------------------------
    gap     :     0.11814076 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final DFT energy:           -112.033908487735                 
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -143.195542382077

     Final gradient norm:           0.000001128393

 
     Date and time (Linux)  : Fri Dec 10 15:03:24 2010
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s     0.0008   0.9893  -0.1333  -0.1871   0.0449   0.2068   0.1104
   2 C   :1s    -0.0100   0.0485   0.2966   0.5784  -0.1468  -1.2732  -0.7076
   3 C   :2pz   -0.0089   0.0013   0.1775  -0.2835  -0.4156   0.4808  -1.1358
   4 O   :1s     0.9927   0.0007  -0.2127   0.0974  -0.1230  -0.0195  -0.1175
   5 O   :1s     0.0333  -0.0100   0.7200  -0.4013   0.6003   0.1079   0.8568
   6 O   :2pz   -0.0070   0.0043  -0.2264  -0.0660   0.6758  -0.1872  -0.9366
   7 H   :1s     0.0004  -0.0110   0.0348   0.2666   0.1248   0.9083  -0.0869

 Molecular orbitals for symmetry species  2
 ------------------------------------------

 Orbital           1        2        3
   1 C   :2px    0.5918  -0.0483   1.1609
   2 O   :2px    0.2788   0.9222  -0.3534
   3 H   :1s     0.3235  -0.3314  -0.8594

 Molecular orbitals for symmetry species  3
 ------------------------------------------

 Orbital           1        2
   1 C   :2py    0.6087   0.8212
   2 O   :2py    0.6772  -0.7657



 >>>> Total CPU  time used in SIRIUS :      1.15 seconds
 >>>> Total wall time used in SIRIUS :      1.00 seconds

 
     Date and time (Linux)  : Fri Dec 10 15:03:24 2010
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


   SCF energy         :     -112.033908487735161
 -- inactive part     :     -143.195542382077065
 -- nuclear repulsion :       31.163673581965142


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
 Electrons in DFTMOMO:   16.00000093514775



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F

 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*1 evaluation time:       0.1 s

 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   18)
 RSP solution vector no.    1; norm of residual   4.67e-05
 RSP solution vector no.    2; norm of residual   2.99e-05

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    ZDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  0.60223618     *ENERGY(eV):   11.976359    
@ STATE NO:    2 *TRANSITION MOMENT:  0.64913498     *ENERGY(eV):   16.960763    


  ******************************************************************************
  *** @ Excit. operator sym 1 & ref. state sym 1 => excited state symmetry 1 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.44012312    au
@                      11.976359     eV
@                      96595.860     cm-1
                       1155.5431     kJ / mol

@ Total energy :      -111.59379     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  0.10641837      (Transition moment :  0.60223618     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.0635542947


 @ Excited state no:    2 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.62329659    au
@                      16.960763     eV
@                      136797.79     cm-1
                       1636.4650     kJ / mol

@ Total energy :      -111.41061     au

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  0.17509491      (Transition moment :  0.64913498     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.1613756756


 Time used in polarization propagator calculation is      1.15 CPU seconds for symmetry 1


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
 Electrons in DFTMOMO:   16.00000093514775



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s

 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   1.75e-04
 RSP solution vector no.    2; norm of residual   2.66e-04

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    XDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -0.58812808     *ENERGY(eV):   13.279055    
@ STATE NO:    2 *TRANSITION MOMENT:  0.11341246     *ENERGY(eV):   17.902102    


  ******************************************************************************
  *** @ Excit. operator sym 2 & ref. state sym 1 => excited state symmetry 2 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.48799632    au
@                      13.279055     eV
@                      107102.81     cm-1
                       1281.2342     kJ / mol

@ Total energy :      -111.54591     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  0.11253021      (Transition moment : -0.58812808     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.7451204338


 @ Excited state no:    2 in symmetry  2
 ---------------------------------------

@ Excitation energy :  0.65789016    au
@                      17.902102     eV
@                      144390.20     cm-1
                       1727.2904     kJ / mol

@ Total energy :      -111.37602     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  5.64135840e-03  (Transition moment :  0.11341246     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.6758683484


 Time used in polarization propagator calculation is      0.81 CPU seconds for symmetry 2


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
 Electrons in DFTMOMO:   16.00000093514775



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  3; triplet =   F

 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s

 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   1.85e-05
 RSP solution vector no.    2; norm of residual   1.83e-06

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    YDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  0.18087490     *ENERGY(eV):   9.3145409    
@ STATE NO:    2 *TRANSITION MOMENT: -0.61568796     *ENERGY(eV):   16.017150    


  ******************************************************************************
  *** @ Excit. operator sym 3 & ref. state sym 1 => excited state symmetry 3 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.34230310    au
@                      9.3145409     eV
@                      75126.847     cm-1
                       898.71666     kJ / mol

@ Total energy :      -111.69161     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  7.46579734e-03  (Transition moment :  0.18087490     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.0051918674


 @ Excited state no:    2 in symmetry  3
 ---------------------------------------

@ Excitation energy :  0.58861946    au
@                      16.017150     eV
@                      129187.04     cm-1
                       1545.4202     kJ / mol

@ Total energy :      -111.44529     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  0.14875264      (Transition moment : -0.61568796     )


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.0084107380


 Time used in polarization propagator calculation is      0.79 CPU seconds for symmetry 3


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
 Electrons in DFTMOMO:   16.00000093514775



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  4; triplet =   F

 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.2 s
RSPORT:    4 out of    2 new trial vectors linear dependent
 Electrons: 16.000001( 9.35e-07): LR-DFT*1 evaluation time:       0.1 s

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

@ Excitation energy :  0.13701204    au
@                      3.7282872     eV
@                      30070.668     cm-1
                       359.72507     kJ / mol

@ Total energy :      -111.89690     au


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.0057430271


 @ Excited state no:    2 in symmetry  4
 ---------------------------------------

@ Excitation energy :  0.42490158    au
@                      11.562160     eV
@                      93255.118     cm-1
                       1115.5789     kJ / mol

@ Total energy :      -111.60901     au


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)

                         Lambda:     0.0039777580


 Time used in polarization propagator calculation is      0.43 CPU seconds for symmetry 4

 >>>> Total CPU  time used in RESPONSE:   3.19 seconds
 >>>> Total wall time used in RESPONSE:   4.00 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   4.44 seconds
 >>>> Total wall time used in DALTON:   5.00 seconds

 
     Date and time (Linux)  : Fri Dec 10 15:03:28 2010
     Host name              : stallo-2.local                          
END REFOUT

