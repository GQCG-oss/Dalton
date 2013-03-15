########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm polarizability linear response hf formaldehyde essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
surf
enehf
tes
nuc
nucchg
diplen
OVERRIDE thr 1.0e-4
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
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS


Atomtypes=3 Charge=+2 Angstrom Generators=2 X Y
Charge=6.0 Atoms=1 Basis=STO-3G  Radius=1.6 Alpha=1.2
C     0.000000     0.000000     0.000000   
Charge=8.0 Atoms=1 Basis=STO-3G  Radius=1.5 Alpha=1.2
O     0.000000     0.000000     1.220000   
Charge=1.0 Atoms=1 Basis=STO-3G  Radius=1.2 Alpha=1.2
H     0.943102     0.000000    -0.544500   

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


     Date and time (Linux)  : Tue Feb 15 12:47:32 2011 
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
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G  Radius=1.6 Alpha=1.2" from the basis set library.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G  Radius=1.5 Alpha=1.2" from the basis set library.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G  Radius=1.2 Alpha=1.2" from the basis set library.
  Used basis set file for basis set for elements with Z =   1 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/STO-3G"


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
   1    0.0000000000e+00    0.0000000000e+00    0.0000000000e+00    1.9200000000e+00
   2    0.0000000000e+00    0.0000000000e+00    2.3054658725e+00    1.8000000000e+00
   3    1.7822044879e+00    0.0000000000e+00   -1.0289558751e+00    1.4400000000e+00
   4   -1.7822044879e+00    0.0000000000e+00   -1.0289558751e+00    1.4400000000e+00


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

  Threshold for integrals:  1.00e-15


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

         281893.2926          38569.4058          33927.3708 MHz
            9.402948            1.286537            1.131695 cm-1


@  Nuclear repulsion energy :   31.163673729192


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
 >>>> Total CPU  time used in HERMIT:   0.04 seconds
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

 
     Date and time (Linux)  : Tue Feb 15 12:47:32 2011 
     Host name              : stanley                                 

 Title lines from ".mol" input file:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.

     Time-dependent Hartree-Fock calculation (random phase approximation).


 Initial molecular orbitals are obtained according to
 ".MOSTART EHT   " input option.

     Wave function specification
     ============================
     For the wave function of type :      >>> HF       <<<
     Number of closed shell electrons         14
     Number of electrons in active shells      0
     Total charge of the molecule              2

     Spin multiplicity                         1
     Total number of symmetries                4
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1    2    3    4
                                       --- |  ---  ---  ---  ---
     Total number of orbitals           12 |    7    3    2    0
     Number of basis functions          12 |    7    3    2    0

      ** Automatic occupation of RHF orbitals **

      -- Initial occupation of symmetries is determined from extended Huckel guess.           
      -- Initial occupation of symmetries is :
     Occupied SCF orbitals               7 |    5    1    1    0

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
   1  -111.457640446     -0.492839822462       2.26e+00  -1.11e+02    5  1  1  0

 Virial theorem: -V/T =      2.031519
      MULPOP C     5.19; O     5.54; H     3.27; 
   2  -111.677074626     -0.477111282628       6.11e-01  -2.19e-01    5  1  1  0

 Virial theorem: -V/T =      2.035809
      MULPOP C     5.36; O     5.36; H     3.28; 
   3  -111.699766529     -0.486412862884       3.47e-01  -2.27e-02    5  1  1  0

 Virial theorem: -V/T =      2.030291
      MULPOP C     5.36; O     5.23; H     3.41; 
   4  -111.711725712     -0.481213941256       1.04e-01  -1.20e-02    5  1  1  0

 Virial theorem: -V/T =      2.031841
      MULPOP C     5.36; O     5.27; H     3.37; 
   5  -111.715904539     -0.481358836365       3.64e-02  -4.18e-03    5  1  1  0

 Virial theorem: -V/T =      2.031350
      MULPOP C     5.34; O     5.28; H     3.38; 
   6  -111.716482216     -0.481398776177       3.71e-03  -5.78e-04    5  1  1  0

 Virial theorem: -V/T =      2.031174
      MULPOP C     5.34; O     5.28; H     3.38; 
   7  -111.716484876     -0.481376100818       1.76e-03  -2.66e-06    5  1  1  0

 Virial theorem: -V/T =      2.031204
      MULPOP C     5.34; O     5.28; H     3.38; 
   8  -111.716485585     -0.481381576444       3.49e-04  -7.09e-07    5  1  1  0

 Virial theorem: -V/T =      2.031191
      MULPOP C     5.34; O     5.28; H     3.38; 
   9  -111.716485634     -0.481381428263       1.18e-05  -4.89e-08    5  1  1  0

 Virial theorem: -V/T =      2.031193
      MULPOP C     5.34; O     5.28; H     3.38; 
  10  -111.716485634     -0.481381339837       2.21e-06  -9.65e-11    5  1  1  0
 DIIS converged in  10 iterations !
 - total time used in SIRFCK        :       0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   14
 Orbital occupations :    5    1    1    0

 Sym       Hartree-Fock orbital energies

  1    -21.08782657   -11.61411485    -1.81843526    -1.14185401    -0.97076128
         0.34806981     0.42545052

  2     -0.93988395    -0.29794431     0.44229638

  3     -0.94384751    -0.11091784

    E(LUMO) :    -0.29794431 au (symmetry 2)
  - E(HOMO) :    -0.93988395 au (symmetry 2)
  ------------------------------------------
    gap     :     0.64193964 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    2

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -111.716485633675                 
     Nuclear repulsion:            31.163673729192
     Electronic energy:          -142.398778023030

     Final gradient norm:           0.000002214765

 
     Date and time (Linux)  : Tue Feb 15 12:47:33 2011 
     Host name              : stanley                                 

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s     0.0005   0.9929  -0.1088  -0.1926  -0.0715   0.2178  -0.0333
   2 C   :1s    -0.0068   0.0317   0.2142   0.6223   0.2853  -1.4190   0.1985
   3 C   :2pz   -0.0058   0.0020   0.1429  -0.0704  -0.4767  -0.0304   1.2459
   4 O   :1s     0.9947   0.0002  -0.2225   0.1222  -0.0533  -0.0591   0.0963
   5 O   :1s     0.0242  -0.0060   0.8000  -0.5677   0.2851   0.4237  -0.7504
   6 O   :2pz   -0.0054   0.0019  -0.2375  -0.4234   0.6099  -0.4830   0.7639
   7 H   :1s     0.0003  -0.0061   0.0182   0.1406   0.1779   0.8287   0.4265

 Molecular orbitals for symmetry species  2
 ------------------------------------------

 Orbital           1        2        3
   1 C   :2px    0.6699  -0.1864   1.1030
   2 O   :2px    0.3562   0.9055  -0.3261
   3 H   :1s     0.2154  -0.3113  -0.8999

 Molecular orbitals for symmetry species  3
 ------------------------------------------

 Orbital           1        2
   1 C   :2py    0.3049   0.9757
   2 O   :2py    0.8913  -0.5006



 >>>> Total CPU  time used in SIRIUS :      0.67 seconds
 >>>> Total wall time used in SIRIUS :      1.00 seconds

 
     Date and time (Linux)  : Tue Feb 15 12:47:33 2011 
     Host name              : stanley                                 


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




 Linear Response calculation
 ---------------------------

 Equilibrium PCM solvent model requested        : SOLVNT =T

 Dielectric constant                            : EPSOL  = 78.3900
 Print level                                    : IPRLR  =   2
 Maximum number of iterations                   : MAXITL =  60
 Threshold for relative convergence             : THCLR  = 1.000e-04
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

  1 B-frequencies  0.000000e+00

    1 second order properties calculated with symmetry no.    1 and labels:

          ZDIPLEN 

    1 second order properties calculated with symmetry no.    2 and labels:

          XDIPLEN 

    1 second order properties calculated with symmetry no.    3 and labels:

          YDIPLEN 


   SCF energy         :     -111.716485633675163
 -- inactive part     :     -142.398778023029990
 -- nuclear repulsion :       31.163673729192165


                     ***************************************
                     *** RHF response calculation (TDHF) ***
                     ***************************************



 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            0
 Number of response properties of this symmetry    1
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      13
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      13


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : ZDIPLEN 
 RSPLR -- frequencies :   0.000000
 FREQ 2    0.0000000000000000     



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   4.74e-05

 *** RSPCTL MICROITERATIONS CONVERGED


           Final output of second order properties from linear response
           ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< ZDIPLEN  ; ZDIPLEN  >> =  1.211571705126e+01


 Time used in linear response calculation is      0.29 CPU seconds for symmetry 1


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    2

 Number of excitations of this symmetry            0
 Number of response properties of this symmetry    1
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12


 RSPLR -- linear response calculation for symmetry  2
 RSPLR -- operator label : XDIPLEN 
 RSPLR -- frequencies :   0.000000
 FREQ 2    0.0000000000000000     



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   7.54e-06

 *** RSPCTL MICROITERATIONS CONVERGED


           Final output of second order properties from linear response
           ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< XDIPLEN  ; XDIPLEN  >> =  7.327583435349e+00


 Time used in linear response calculation is      0.36 CPU seconds for symmetry 2


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    3

 Number of excitations of this symmetry            0
 Number of response properties of this symmetry    1
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7


 RSPLR -- linear response calculation for symmetry  3
 RSPLR -- operator label : YDIPLEN 
 RSPLR -- frequencies :   0.000000
 FREQ 2    0.0000000000000000     



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   9.77e-05

 *** RSPCTL MICROITERATIONS CONVERGED


           Final output of second order properties from linear response
           ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< YDIPLEN  ; YDIPLEN  >> =  2.908427928051e+00


 Time used in linear response calculation is      0.25 CPU seconds for symmetry 3

 >>>> Total CPU  time used in RESPONSE:   0.91 seconds
 >>>> Total wall time used in RESPONSE:   0.00 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   1.62 seconds
 >>>> Total wall time used in DALTON:   1.00 seconds

 
     Date and time (Linux)  : Tue Feb 15 12:47:33 2011 
     Host name              : stanley                                 
END REFOUT

