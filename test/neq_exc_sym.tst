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


    3    1 X       1
        6.    1
C     0.000000     0.000000     0.000000
        8.    1
O     0.000000     0.000000     1.220000
        1.    1
H     0.943102     0.000000    -0.544500
END MOLINP

########## Reference Output ########################
START REFOUT


    ******************************************************************
    ***********  DALTON - An electronic structure program  ***********
    ******************************************************************

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
    Jakob Kongsted,           Aarhus University,            Denmark    
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

     Date and time (Linux)  : Tue Mar 21 15:41:32 2006
     Host name              : star.chem.uit.no                        

 <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>


 Default print level:        0

    Integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed
    Dynamic molecular property section will be executed

 ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT=WATER    **
 OPTICAL AND PHYSICAL CONSTANTS:
 EPS= 78.390; EPSINF=  1.776; RSOLV=  1.385 A; VMOL=  18.070 ML/MOL;
 TCE= .25700E-03 1/K; STEN= 71.810 DYN/CM;  DSTEN=  0.6500; CMF=  1.2770


     -----------------------------------
     INPUT FOR PCM SOLVATION CALCULATION 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=WATER        EPS   = 78.3900     EPSINF=  1.7760
     RSOLV =  1.3850

     ICESPH =       2     NESFP =       4
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =T
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.7000
     2    0.0000    0.0000    0.0000    1.5000
     3    0.0000    0.0000    0.0000    1.2000
     4    0.0000    0.0000    0.0000    1.2000

Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************


 Default print level:        1

 * Nuclear model: Point charge

 Calculation of one- and two-electron Hamiltonian integrals.

 The following one-electron property integrals are calculated:
          - overlap integrals
          - dipole length integrals


 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

                                                                          
                                                                          

  Coordinates are entered in Angstroms and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A
  Used basis set file for basis set for elements with Z =   6 :
     "/home/lara/programs/main-branch/dalton/basis/STO-3G"
  Used basis set file for basis set for elements with Z =   8 :
     "/home/lara/programs/main-branch/dalton/basis/STO-3G"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/lara/programs/main-branch/dalton/basis/STO-3G"


  Symmetry Operations
  -------------------

  Symmetry operations: 1



                      SYMGRP:Point group information
                      ------------------------------

Point group: Cs 

   * The point group was generated by:

      Reflection in the yz-plane

   * Group multiplication table

        |  E   Oyz
   -----+----------
     E  |  E 
    Oyz | Oyz   E 

   * Character table

        |  E   Oyz
   -----+----------
    A'  |   1    1
    A"  |   1   -1

   * Direct product table

        | A'   A" 
   -----+----------
    A'  | A' 
    A"  | A"   A' 
 ********SPHERES IN SPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.1700000000E+01
   2    0.0000000000E+00    0.0000000000E+00    0.2305465883E+01    0.1500000000E+01
   3    0.1782204496E+01    0.0000000000E+00   -0.1028955880E+01    0.1200000000E+01
   4   -0.1782204496E+01    0.0000000000E+00   -0.1028955880E+01    0.1200000000E+01


                             Isotopic Masses
                             ---------------

                           C          12.000000
                           O          15.994915
                           H    1      1.007825
                           H    2      1.007825

                       Total mass:    30.010565 amu
                       Natural abundance:  98.633 %

 Center-of-mass coordinates (A):    0.000000    0.000000    1.159649


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

  Threshold for integrals:  1.00E-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates:   12

   1   C        x      0.0000000000
   2            y      0.0000000000
   3            z      0.0000000000

   4   O        x      0.0000000000
   5            y      0.0000000000
   6            z      2.3054658834

   7   H    1   x      1.7822044964
   8            y      0.0000000000
   9            z     -1.0289558799

  10   H    2   x     -1.7822044964
  11            y      0.0000000000
  12            z     -1.0289558799


  Symmetry Coordinates
  --------------------

  Number of coordinates in each symmetry:    7   5

  Symmetry  A'  ( 1)

   1   C     y    2
   2   C     z    3
   3   O     y    5
   4   O     z    6
   5   H     x    [  7  -  10 ]/2
   6   H     y    [  8  +  11 ]/2
   7   H     z    [  9  +  12 ]/2

  Symmetry  A"  ( 2)

   8   C     x    1
   9   O     x    4
  10   H     x    [  7  +  10 ]/2
  11   H     y    [  8  -  11 ]/2
  12   H     z    [  9  -  12 ]/2


   Interatomic separations (in Angstroms):
   ---------------------------------------

            C           O           H    1      H    2
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H    1:    1.089000    2.000725    0.000000
 H    2:    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstroms
  between atoms "H    1" and "O     ".


  Bond distances (angstroms):
  ---------------------------

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

   IA    1.792803          0.000000    0.000000    1.000000
   IB   13.103106          1.000000    0.000000    0.000000
   IC   14.895908          0.000000    1.000000    0.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         281893.2924          38569.4057          33927.3707 MHz
            9.402948            1.286537            1.131695 cm-1


  Nuclear repulsion energy :   31.163673581965


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:               9   3


  Symmetry  A' ( 1)

    1     C        1s         1
    2     C        1s         2
    3     C        2py        4
    4     C        2pz        5
    5     O        1s         6
    6     O        1s         7
    7     O        2py        9
    8     O        2pz       10
    9     H        1s        11  +  12


  Symmetry  A" ( 2)

   10     C        2px        3
   11     O        2px        8
   12     H        1s        11  -  12

  Symmetries of electric field:  A" (2)  A' (1)  A' (1)

  Symmetries of magnetic field:  A' (1)  A" (2)  A" (2)


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.01 seconds


 >>> Time used in OVERLA is   0.00 seconds


 >>> Time used in DIPLEN is   0.00 seconds


 >>> Time used in ONEDRV is   0.00 seconds


 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014



 >>> Time used in TWOINT is   0.01 seconds


 MEMORY USED TO GENERATE CAVITY=    432042


 TOTAL NUMBER OF SPHERES=    4
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.000000000    0.000000000    0.000000000    2.040000000   25.046551891
   2    0.000000000    0.000000000    1.220000000    1.800000000   22.984715878
   3    0.943102000    0.000000000   -0.544500000    1.440000000    9.281425262
   4   -0.943102000    0.000000000   -0.544500000    1.440000000    9.281425262

 TOTAL NUMBER OF TESSERAE =     304
 SURFACE AREA=   66.59411829 (A**2)    CAVITY VOLUME=   48.28620692 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.05 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.11 seconds

 >>>> Total CPU  time used in HERMIT:   0.28 seconds
 >>>> Total wall time used in HERMIT:   1.00 seconds

- End of Integral Section


Starting in Wave Function Section -


 *** Output from Huckel module :

     Using EWMO model:          F
     Using EHT  model:          T
     Number of Huckel orbitals each symmetry:    9    3

 Huckel EHT eigenvalues for symmetry :  1
          -20.808771     -11.553164      -2.104598      -1.278134      -0.811063
           -0.584721      -0.212037       0.240410       0.392263

 Huckel EHT eigenvalues for symmetry :  2
           -1.053261      -0.439652       0.161628
 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Tue Mar 21 15:41:33 2006
     Host name              : star.chem.uit.no                        

 Title lines from integral program:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham calculation.


     Time-dependent Kohn-Sham calculation (random phase approximation).


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Number of active orbitals                 0
     Total number of orbitals                 12

     Spin multiplicity                         1
     Total number of symmetries                2
     Reference state symmetry                  1
 
     This is a DFT calculation of type: LDA

     Orbital specifications
     ======================
     Abelian symmetry species           1   2
                                       --  --
     Total number of orbitals           9   3
     Number of basis functions          9   3

      ** Automatic occupation of RKS orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals              6   2

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00E-06

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.97018
 NUCLEAR APPARENT CHARGE -15.78962 THEORETICAL -15.79589 NOT RENORMALIZED
 this is icompcm in icvev 0

  ..... DONE WITH INDUCED NUCLEAR CHARGES .....

 >>> Time used in VNN    is   0.00 seconds



 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  16 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Atom:    1*1 points=18676 compressed from 18676 ( 92 radial)
 Atom:    2*1 points=18206 compressed from 18280 ( 92 radial)
 Atom:    3*2 points=18039 compressed from 18150 ( 75 radial)
 Number of grid points:    54921 Grid generation time:       0.1 s
K-S electrons/energy :   16.00000086596893  -13.89075466270899 err:0.87E-06
   1   -111.738853559538     -0.008423863866   3.69E+00  -1.12E+02    6  2
K-S electrons/energy :   16.00000101615090  -13.74986834238559 err:0.10E-05
   2   -111.702652599968     -0.010498400450   2.09E+00   3.62E-02    6  2
K-S electrons/energy :   16.00000085037999  -14.08065623512538 err:0.85E-06
   3   -111.616381631740     -0.038860487017   2.09E+00   8.63E-02    6  2
K-S electrons/energy :   16.00000093909577  -13.80409388238066 err:0.94E-06
   4   -112.028220098204     -0.000677321270   2.73E-01  -4.12E-01    6  2
K-S electrons/energy :   16.00000093853158  -13.83311992294588 err:0.94E-06
   5   -112.033588055744     -0.001843302221   6.72E-02  -5.37E-03    6  2
K-S electrons/energy :   16.00000093561215  -13.83234091924215 err:0.94E-06
   6   -112.033874517507     -0.001903223067   2.11E-02  -2.86E-04    6  2
K-S electrons/energy :   16.00000093512555  -13.83462390385552 err:0.94E-06
   7   -112.033908389410     -0.002047434197   1.11E-03  -3.39E-05    6  2
K-S electrons/energy :   16.00000093514889  -13.83450041977549 err:0.94E-06
   8   -112.033908487586     -0.002039420477   4.49E-05  -9.82E-08    6  2
K-S electrons/energy :   16.00000093514754  -13.83450477936072 err:0.94E-06
   9   -112.033908487734     -0.002039715185   3.75E-06  -1.48E-10    6  2
K-S electrons/energy :   16.00000093514766  -13.83450448923482 err:0.94E-06
  10   -112.033908487735     -0.002039695267   3.43E-07  -1.44E-12    6  2
 DIIS converged in  10 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    6    2

 Kohn-Sham orbital energies, symmetry 1

       -18.36962379    -9.66197402    -0.90481806    -0.50692423    -0.27783863
        -0.27474214     0.02488548     0.36676356     0.53279078

 Kohn-Sham orbital energies, symmetry 2

        -0.38754504    -0.09325514     0.46802248

    E(LUMO) :     0.02488548 au (symmetry 1)
  - E(HOMO) :    -0.09325514 au (symmetry 2)
  ------------------------------------------
    gap     :     0.11814062 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final DFT energy:           -112.033908487735
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -143.195542374434

     Final gradient norm:           0.000000342568

 
     Date and time (Linux)  : Tue Mar 21 15:41:45 2006
     Host name              : star.chem.uit.no                        

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5        6
   1  C   1s     0.0008   0.9893  -0.1333   0.1871   0.0449   0.0000
   2  C   1s    -0.0100   0.0485   0.2966  -0.5784  -0.1468   0.0000
   3  C   2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.6087
   4  C   2pz   -0.0089   0.0013   0.1775   0.2835  -0.4156   0.0000
   5  O   1s     0.9927   0.0007  -0.2127  -0.0974  -0.1230   0.0000
   6  O   1s     0.0333  -0.0100   0.7200   0.4013   0.6003   0.0000
   7  O   2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.6772
   8  O   2pz   -0.0070   0.0043  -0.2264   0.0660   0.6758   0.0000
   9  H   1s     0.0004  -0.0110   0.0348  -0.2666   0.1248   0.0000

     Molecular orbitals for symmetry species   2

 Orbital          1        2
   1  C   2px    0.5918  -0.0483
   2  O   2px    0.2788   0.9222
   3  H   1s     0.3235  -0.3314



 >>>> Total CPU  time used in SIRIUS :      6.11 seconds
 >>>> Total wall time used in SIRIUS :     12.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 15:41:45 2006
     Host name              : star.chem.uit.no                        

- End of Wave Function Section



  This is output from RESPONSE  -  an MCSCF and SOPPA response property program
 ------------------------------------------------------------------------------



 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




  Linear Response single residue calculation
 -------------------------------------------

 Non-equilibrium PCM solvent model requested    : INERSI =T

 Static dielectric constant                     : EPSTAT = 78.3900
 Optical dielectric constant                    : EPSOL  =  1.7760
 Print level                                    : IPRPP  =   2
 Maximum number of iterations for eigenval.eqs. : MAXITP =  60
 Threshold for convergence of eigenval.eqs.     : THCPP  = 1.000E-03
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

      2 Excitation energies are calculated for symmetry no.    1

      2 property residues are calculated with labels:

               YDIPLEN 
               ZDIPLEN 

      2 Excitation energies are calculated for symmetry no.    2

      1 property residues are calculated with labels:

               XDIPLEN 


   SCF energy         :     -112.033908487735005
 -- inactive part     :     -143.195542374433614
 -- nuclear repulsion :       31.163673581965146


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            2
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      20
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      20



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F

 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.9 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.9 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.9 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.9 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*1 evaluation time:       0.6 s

 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   18)
 RSP solution vector no.    1; norm of residual   1.42E-05
 RSP solution vector no.    2; norm of residual   3.86E-04

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    YDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  0.18087285     *ENERGY(eV):   9.3145380    
@ STATE NO:    2 *TRANSITION MOMENT:  2.53806052E-16 *ENERGY(eV):   11.976357    

@ Transition operator type:    ZDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  4.87550239E-17 *ENERGY(eV):   9.3145380    
@ STATE NO:    2 *TRANSITION MOMENT: -0.60224149     *ENERGY(eV):   11.976357    


@ Reference  state    symmetry  1
@ Excitation operator symmetry  1
@ Excited    state    symmetry  1

@ State no:    1

@ Excitation energy :  0.34230299     au
@                      9.3145380     eV
@                      75126.823     cm-1
                       898.71638     kJ / mol

@ Total energy :      -111.69161     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  7.46562526E-03
 Transition moment   (LENGTH)   :  0.18087285    

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  5.42448091E-34
 Transition moment   (LENGTH)   :  4.87550239E-17

@ State no:    2

@ Excitation energy :  0.44012307     au
@                      11.976357     eV
@                      96595.847     cm-1
                       1155.5429     kJ / mol

@ Total energy :      -111.59379     au

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  1.89010885E-32
 Transition moment   (LENGTH)   :  2.53806052E-16

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  0.10642024    
 Transition moment   (LENGTH)   : -0.60224149    


 Time used in polarization propagator calculation is      5.21 CPU seconds for symmetry 1


 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    2

 Number of excitations of this symmetry            2
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 16.000001( 9.35e-07): LR-DFT*2 evaluation time:       0.9 s
 Electrons: 16.000001( 9.35e-07): LR-DFT*1 evaluation time:       0.6 s

 *** THE REQUESTED    2 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   1.02E-07
 RSP solution vector no.    2; norm of residual   9.00E-08

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    XDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  6.18572249E-17 *ENERGY(eV):   3.7282834    
@ STATE NO:    2 *TRANSITION MOMENT:  2.97840311E-16 *ENERGY(eV):   11.562162    


@ Reference  state    symmetry  1
@ Excitation operator symmetry  2
@ Excited    state    symmetry  2

@ State no:    1

@ Excitation energy :  0.13701190     au
@                      3.7282834     eV
@                      30070.636     cm-1
                       359.72469     kJ / mol

@ Total energy :      -111.89690     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  3.49500576E-34
 Transition moment   (LENGTH)   :  6.18572249E-17

@ State no:    2

@ Excitation energy :  0.42490167     au
@                      11.562162     eV
@                      93255.137     cm-1
                       1115.5792     kJ / mol

@ Total energy :      -111.60901     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  2.51283592E-32
 Transition moment   (LENGTH)   :  2.97840311E-16


 Time used in polarization propagator calculation is      1.85 CPU seconds for symmetry 2

 >>>> Total CPU  time used in RESPONSE:   7.06 seconds
 >>>> Total wall time used in RESPONSE:  15.00 seconds
 >>>> Total CPU  time used in DALTON:  13.45 seconds
 >>>> Total wall time used in DALTON:  28.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 15:42:00 2006
     Host name              : star.chem.uit.no                        
END REFOUT

