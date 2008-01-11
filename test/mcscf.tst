########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm mcscf water energy essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
elactive
dipole
OVERRIDE thr 1.0e-4
nuc
tes
sym
cmass
elactive
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
3
.ICESPH
2
*PCMCAV
.INA
1
2
3
.RIN
1.5
1.2
1.2
.AREATS
0.30
**WAVE FUNCTIONS
.HF
.MP2
.MCSCF
*SCF INPUT
.DOUBLY OCCUPIED
 3 1 1 0
*CONFIGURATION INPUT
.SYMMETRY
 1
.SPIN MUL
 1
.INACTIVE
 1 0 0 0
.ELECTRONS
 8
.CAS SPACE
 4 2 2 0
**PROPERTIES
**END OF INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
6-311G*
magnet
nuova base
    2    2 X Y
        8.    1
O      0.0000000000  0.0000000000  -0.1258515023
        1.    1
H      0.0000000000  1.4523500000  0.9986773907
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

     Date and time (Linux)  : Tue Mar 21 17:51:39 2006
     Host name              : star.chem.uit.no                        

 <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>


 Default print level:        0

    Integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed
    Static molecular property section will be executed

 ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT=WATER    **
 OPTICAL AND PHYSICAL CONSTANTS:
 EPS= 78.390; EPSINF=  1.776; RSOLV=  1.385 A; VMOL=  18.070 ML/MOL;
 TCE= .25700E-03 1/K; STEN= 71.810 DYN/CM;  DSTEN=  0.6500; CMF=  1.2770


     -----------------------------------
     INPUT FOR PCM SOLVATION CALCULATION 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=WATER        EPS   = 78.3900     EPSINF=  1.7760
     RSOLV =  1.3850

     ICESPH =       2     NESFP =       3
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.5000
     2    0.0000    0.0000    0.0000    1.2000
     3    0.0000    0.0000    0.0000    1.2000

Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************



 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

  magnet                                                                  
  nuova base                                                              
  Used basis set file for basis set for elements with Z =   8 :
     "/home/lara/programs/main-branch/dalton/basis/6-311G*"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/lara/programs/main-branch/dalton/basis/6-311G*"


  Symmetry Operations
  -------------------

  Symmetry operations: 2



                      SYMGRP:Point group information
                      ------------------------------

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
 ********SPHERES IN SPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000E+00    0.0000000000E+00   -0.1258515023E+00    0.1500000000E+01
   2    0.0000000000E+00    0.1452350000E+01    0.9986773907E+00    0.1200000000E+01
   3    0.0000000000E+00   -0.1452350000E+01    0.9986773907E+00    0.1200000000E+01


                             Isotopic Masses
                             ---------------

                           O          15.994915
                           H    1      1.007825
                           H    2      1.007825

                       Total mass:    18.010565 amu
                       Natural abundance:  99.730 %

 Center-of-mass coordinates (A):    0.000000    0.000000    0.000000


  Atoms and basis sets
  --------------------

  Number of atom types:     2
  Total number of atoms:    3

  Basis set used is "6-311G*" from the basis set library.

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  O           1    8.0000    31    18      [11s5p1d|4s3p1d]                                   
  H           2    1.0000     5     3      [5s|3s]                                            
  ----------------------------------------------------------------------
  total:      3   10.0000    41    24
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00E-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates:    9

   1   O        x      0.0000000000
   2            y      0.0000000000
   3            z     -0.1258515023

   4   H    1   x      0.0000000000
   5            y      1.4523500000
   6            z      0.9986773907

   7   H    2   x      0.0000000000
   8            y     -1.4523500000
   9            z      0.9986773907


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


   Interatomic separations (in Angstroms):
   ---------------------------------------

            O           H    1      H    2
            ------      ------      ------
 O     :    0.000000
 H    1:    0.972000    0.000000
 H    2:    0.972000    1.537101    0.000000


  Max interatomic separation is    1.5371 Angstroms
  between atoms "H    2" and "H    1".


  Bond distances (angstroms):
  ---------------------------

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

   IA    0.633889          0.000000    1.000000    0.000000
   IB    1.190584          0.000000    0.000000    1.000000
   IC    1.824473          1.000000    0.000000    0.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         797267.3440         424480.0013         277000.0234 MHz
           26.593976           14.159129            9.239726 cm-1


  Nuclear repulsion energy :    9.055004525638


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:              12   4   7   1


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
   10     H        1s        19  +  20
   11     H        1s        21  +  22
   12     H        1s        23  +  24


  Symmetry  B1 ( 2)

   13     O        2px        5
   14     O        2px        8
   15     O        2px       11
   16     O        3d1+      17


  Symmetry  B2 ( 3)

   17     O        2py        6
   18     O        2py        9
   19     O        2py       12
   20     O        3d1-      15
   21     H        1s        19  -  20
   22     H        1s        21  -  22
   23     H        1s        23  -  24


  Symmetry  A2 ( 4)

   24     O        3d2-      14

  Symmetries of electric field:  B1 (2)  B2 (3)  A1 (1)

  Symmetries of magnetic field:  B2 (3)  B1 (2)  A2 (4)


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.01 seconds


 >>> Time used in ONEDRV is   0.01 seconds


 Number of two-electron integrals written:       11044 ( 24.5% )
 Megabytes written:                              0.131



 >>> Time used in TWOINT is   0.05 seconds


 MEMORY USED TO GENERATE CAVITY=    432042


 TOTAL NUMBER OF SPHERES=    3
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.000000000    0.000000000   -0.066597747    1.800000000   24.581233228
   2    0.000000000    0.768550518    0.528477314    1.440000000   11.987539425
   3    0.000000000   -0.768550518    0.528477314    1.440000000   11.987539425

 TOTAL NUMBER OF TESSERAE =     200
 SURFACE AREA=   48.55631208 (A**2)    CAVITY VOLUME=   30.53285203 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.12 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.05 seconds

 >>>> Total CPU  time used in HERMIT:   0.17 seconds
 >>>> Total wall time used in HERMIT:   1.00 seconds

- End of Integral Section


Starting in Wave Function Section -


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

 
     Date and time (Linux)  : Tue Mar 21 17:51:40 2006
     Host name              : star.chem.uit.no                        

 Title lines from integral program:
     magnet                                                                  
     nuova base                                                              

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     MC-SCF optimization.

               Type: complete active space calculation (CAS).

     This is a combination run starting with

               an RHF calculation
               an MP2 calculation 


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons          2
     Number of electrons in active shells      8
     Total charge of the molecule              0

     Number of active orbitals                 8
     Total number of orbitals                 24

     Spin multiplicity                         1
     Total number of symmetries                4
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species           1   2   3   4
                                       --  --  --  --
     Inactive orbitals                  1   0   0   0
     Active orbitals                    4   2   2   0
     Secondary orbitals                 7   2   5   1
     Total number of orbitals          12   4   7   1
     Number of basis functions         12   4   7   1

     Occupied SCF orbitals              3   1   1   0

     Optimization information
     ========================
     Number of configurations            492
     Number of orbital rotations          53
     ---------------------------------------
     Total number of variables           545

     Maximum number of macro iterations      15
     Maximum number of micro iterations     360
     Threshold for gradient            1.00E-05
     Number of initial trial vectors          1
     Number of initial CI iterations          3
     Number of simultaneous trial vectors     1

     This calculation converges to the lowest state
     for the specified symmetry and spin species.

     Maximum number of NEO/NR iterations  24

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE        9.97305
 NUCLEAR APPARENT CHARGE  -9.86625 THEORETICAL  -9.87243 NOT RENORMALIZED
 this is icompcm in icvev 0

  ..... DONE WITH INDUCED NUCLEAR CHARGES .....

 >>> Time used in VNN    is   0.00 seconds



 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy        Error norm    Delta(E)  DIIS dim.
 ---------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00E-15

 >>> Time used in FRMSUP is   0.01 seconds

    1      -75.717172113279    3.53974E+00   -7.57E+01    1
    2      -76.012895061562    4.97183E-01   -2.96E-01    2
    3      -76.027678049804    2.10435E-01   -1.48E-02    3
    4      -76.029739067115    6.99777E-02   -2.06E-03    4
    5      -76.030110158564    1.12229E-02   -3.71E-04    5
    6      -76.030126498174    3.14195E-03   -1.63E-05    6
    7      -76.030128360047    4.28435E-04   -1.86E-06    7
    8      -76.030128390117    5.59025E-05   -3.01E-08    8
    9      -76.030128390504    1.07209E-05   -3.87E-10    9
   10      -76.030128390518    1.16780E-06   -1.48E-11   10
   11      -76.030128390518    8.22769E-08   -1.46E-13   10
 DIIS converged in  11 iterations !


 >>>>> Output from SIRIUS MP2 module <<<<<

 Reference: H.J.Aa.Jensen, P.Jorgensen, H.Agren, and J.Olsen,
            J. Chem. Phys. 88, 3834 (1988); 89, 5354 (1988)



 Check that orbitals are canonical HARTREE-FOCK orbitals

 Number of electrons :   10
 Orbital occupations :    3    1    1    0


 Hartree-Fock electronic energy:         -85.085132916156

 Hartree-Fock total      energy:         -76.030128390518

 Hartree-Fock orbital energies, symmetry 1

       -20.55013880    -1.33874397    -0.56983933     0.14498514     0.59812870
         1.02164275     1.43875677     2.37508250     3.13899157     3.55697683
         5.42351279    51.43285032

 Hartree-Fock orbital energies, symmetry 2

        -0.49955294     1.00283930     3.17787350     5.30999732

 Hartree-Fock orbital energies, symmetry 3

        -0.70288966     0.21350053     0.55399625     1.14905796     2.31299616
         3.88776923     5.55640495

 Hartree-Fock orbital energies, symmetry 4

         3.15309355

    E(LUMO) :     0.14498514
  - E(HOMO) :    -0.49955294
  --------------------------
    gap     :     0.64453809

 ABS SUM OF OFF-DIAGONAL FOCK CORE ELEMENTS ARE :  2.12E-07

 MP2 move   0.098937 electrons to unoccupied HF orbitals


 Hartree-Fock total energy :           -76.0301283905
 + MP2 contribution        :            -0.2230088208
 = MP2 second order energy :           -76.2531372114

 Natural orbital occupation numbers, symmetry 1

         1.99956628     1.98602239     1.97116241     0.02242124     0.01099086
         0.00428229     0.00311040     0.00098523     0.00056418     0.00024913
         0.00015539     0.00002595

 Sum     5.99953575
 RHF     6.00000000
 Diff.  -0.00046425

 Natural orbital occupation numbers, symmetry 2

         1.97467335     0.01818010     0.00322096     0.00049336

 Sum     1.99656777
 RHF     2.00000000
 Diff.  -0.00343223

 Natural orbital occupation numbers, symmetry 3

         1.97004253     0.02413194     0.00474191     0.00083962     0.00061812
         0.00019639     0.00000206

 Sum     2.00057258
 RHF     2.00000000
 Diff.   0.00057258

 Natural orbital occupation numbers, symmetry 4

         0.00332391

 Sum     0.00332391
 RHF     0.00000000
 Diff.   0.00332391

 Time used for MP2 natural orbitals :       0.039 CPU seconds.


        SIRIUS MCSCF optimization (SIROPT)
 ================================================



 <<<<< Output from SIRIUS CI module (CICTL) >>>>>




 (CIST1)  4 lowest diagonal elements:

 Element no. Config.no.    Active energy      Total energy

         1 :         20    -23.8059466955    -76.0288766769
         2 :        131    -22.8436418902    -75.0665718716
         3 :         80    -22.8214819334    -75.0444119148
         4 :        120    -22.5526173979    -74.7755473793


 Convergence threshold for CI optimization :     0.00000500



 *** Reached maximum number of CI iterations:    3
              1 CI roots are not converged.


 CI energies and residuals:
    1      -76.155041847940254       5.71E-02

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  1

   1.988879150   1.975748370   0.024397104   0.011755463

 Symmetry  2

   1.980441312   0.018378999

 Symmetry  3

   1.974758240   0.025641361


 <<< MACRO ITERATION  1 >>>
 --------------------------

 Total MCSCF energy       :      -76.165426189716229       (MACRO    1)
 - Dielec. solvation ener.:       -0.010384341776244

 Norm of total gradient   :        0.146390547946
 -    of CI gradient      :        0.121251874691
 -    of orbital gradient :        0.082025455877
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.08412428
       2       0.02227629
       3       0.17778443
       4      -0.17005695
       5      -0.19116562
       6      -0.29599173
       7      -0.11928063
       8       0.04947570
       9       0.01551645
      10      -0.15967962
      11      -0.06672088
      12      -0.08312596
      13      -0.02381554
      14      -0.05956895
      15       0.01315682
      16      -0.16612778
      17       0.09496828
      18      -0.19340281
      19      -0.14626031
      20      -0.04412329
      21       0.19501324
      22      -0.27505565
      23       0.31814502
      24       0.03589644
      25      -0.22940076
      26      -0.02602562
      27       0.12762517
      28      -0.17965116
      29      -0.05108847
      30      -0.38930377
      31       0.54553916
      32      -0.14187140
      33      -0.02406369
      34      -0.15271335
      35      -0.05778792
      36      -0.21508277
      37      -0.20041971
      38       0.07720945
      39       0.02517799
      40       0.23408260
      41       0.21710193
      42      -0.27866919
      43      -0.03450774
      44      -0.34722747
      45       0.61956457
      46       0.02397196
      47       1.30746383
      48       0.08095189
      49       0.13085786
      50       0.46196655
      51       0.24009690
      52      -0.07474613
      53      -0.17106871
      54       0.04908960
      55       0.01569201
      56      -0.32302722
      57      -0.15926228
      58      -0.10796615
      59      -0.09992379
      60      -0.15495401
      61       0.03347608
      62      -0.25664861
      63      -0.25574390
      64      -0.33599600
      65      -0.19720739
      66      -0.26272855
      67      -0.26153201
      68       0.08701389
      69       0.03104230
      70      -0.09226862
      71      -0.05495812
      72       0.01804865
      73       0.04154139
      74       0.08833367
      75       0.00303065
      76       0.04312205
      77      -0.01640306
      78      -0.01525398
      79      -0.42109032
      80      -0.64100230
      81       0.00195559
      82       0.97842329
      83       0.00939127
      84      -0.03447605
      85      -0.04721786
      86       0.22953555
      87       0.04258555
      88       0.10870850
      89      -0.21795479
      90       0.05098166
      91       0.19978911
      92       0.55476146
      93       0.02164346
      94       0.71133769
      95      -0.01184007
      96      -0.03837298
      97       0.00064479
      98       0.06799477
      99      -0.01699387
     100      -0.16472919
     101       0.01062951
     102      -0.01905400
     103       0.00717214
     104      -0.00578070
     105       0.02721453
     106       0.01838129
     107      -0.04626066
     108       0.00797558
     109       0.04129618
     110       0.06215399
     111       0.00060672
     112       0.00604977
     113       0.00271837
     114      -0.02556275
     115      -0.02258827
     116      -0.00072670
     117       0.00347583
     118       0.00352448
     119      -0.30890251
     120       0.90010402
     121      -0.04764729
     122       0.63382036
     123      -0.00616827
     124       0.02407083
     125      -0.15718624
     126       0.21812114
     127       0.00388420
     128       0.12340282
     129      -0.09165024
     130      -0.10688683
     131       0.07401896
     132      -1.01301415
     133       0.05473003
     134       0.15915947
     135       0.74835138
     136      -0.01320521
     137      -0.27689622
     138       0.23350793
     139      -0.06186037
     140      -0.04814814
     141      -0.09922582
     142      -0.13769657
     143      -0.53988115
     144      -0.27768850
     145       0.43124074
     146      -0.18357513
     147       0.15329523
     148      -0.00354954
     149      -0.00532527
     150       0.05208355
     151      -0.02263574
     152      -0.06785617
     153      -0.00525166
     154       0.19509088
     155       0.13848622
     156      -0.02513524
     157      -0.04942408
     158       0.00935372
     159      -0.02271881
     160       0.03637624
     161       0.15482183
     162       0.57904328
     163      -0.02843481
     164       1.08413307
     165      -0.00691375
     166      -0.02110013
     167      -0.04339308
     168       0.11322868
     169      -0.00084308
     170      -0.06180899
     171       0.00312520
     172       0.03843228
     173      -0.89116187
     174      -0.47410431
     175       0.82375605
     176      -0.36513826
     177       0.07158259
     178       0.00640014
     179      -0.04920101
     180       0.50795459
     181      -0.10151194
     182       0.07486502
     183      -0.14896981
     184      -0.06159962
     185      -0.80417847
     186      -0.23846346
     187       0.49018347
     188      -0.15778101
     189       0.12783601
     190      -0.12566497
     191      -0.26248162
     192       0.13933058
     193      -0.02481226
     194       0.30515478
     195      -0.03397429
     196       0.00597702
     197      -0.04569181
     198       0.01515627
     199       0.01934269
     200       0.01255565
     201      -0.00553272
     202       0.01793429
     203       0.02685339
     204      -0.01042931
     205       0.05073578
     206      -0.04009317
     207      -0.00262867
     208      -0.01233440
     209       0.05404429
     210       0.04407466
     211       0.02381198
     212      -0.00323638
     213       0.00044338
     214      -0.09227562
     215      -0.05881697
     216       0.00398539
     217       0.01757878
     218      -0.00609700
     219       0.29213770
     220      -0.04418355
     221      -0.01080389
     222      -0.02346636
     223       0.07869662
     224       0.07810245
     225      -0.06831661
     226       0.02688643
     227       0.09124156
     228       0.00755907
     229       0.01130295
     230       0.02911524
     231       0.00991204
     232      -0.02718387
     233      -0.00660559
     234      -0.07110536
     235       0.38748754
     236      -0.60553468
     237       0.03877657
     238       0.25703787
     239       0.00518456
     240       0.10298482
     241      -0.17175666
     242       1.35275824
     243      -0.10763706
     244      -0.27634200
     245      -0.04147772
     246      -0.12648344
     247      -0.05408055
     248      -0.08172799
     249       0.03256850
     250       0.12674984
     251       0.10870529
     252       0.16025053
     253      -0.02308404
     254       0.14899583
     255      -0.03478223
     256      -0.13576441
     257       0.01710921
     258      -0.06298882
     259      -0.12834768
     260      -0.13749286
     261      -0.00025396
     262      -0.03067311
     263       0.00207462
     264       0.00332017
     265      -0.00329629
     266       1.17044014
     267      -0.57502215
     268       0.01123327
     269       0.20200439
     270      -0.32902696
     271       0.00818058
     272       0.64889237
     273      -0.27692461
     274      -0.34923580
     275       0.46483407
     276      -0.73910612
     277      -0.11893437
     278      -0.04583218
     279      -0.10299088
     280      -0.40329272
     281      -0.26034398
     282       0.57423769
     283      -0.06114962
     284      -0.14034722
     285       0.12584923
     286       0.18275455
     287       0.22907441
     288      -0.35992445
     289       0.33795938
     290      -0.06114391
     291       0.46683850
     292      -0.98213259
     293      -0.13755218
     294      -0.58304299
     295       0.02901609
     296       0.28437779
     297       0.13678660
     298      -0.22554354
     299      -0.24859268
     300      -0.00710663
     301       0.00107005
     302      -0.07990233
     303      -0.00747247
     304      -0.20359367
     305       0.02266460
     306       1.62112337
     307       0.01003477
     308      -0.23056361
     309      -0.01202377
     310       0.02854563
     311      -0.03096687
     312      -0.44330028
     313      -0.00359339
     314      -1.03030210
     315       0.01689511
     316      -0.03484462
     317      -0.05086417
     318       0.02542798
     319       0.02526530
     320       0.04851944
     321       0.10881477
     322       0.12599319
     323      -0.02307277
     324      -0.03461874
     325      -0.03554001
     326      -0.03074228
     327       0.00609818
     328       0.00365275
     329      -0.29399470
     330       0.55882064
     331      -0.12758485
     332      -0.25554740
     333       0.02020735
     334       0.40050709
     335      -0.23488895
     336      -0.20591344
     337       0.52410622
     338      -0.72214017
     339       0.03706134
     340      -0.04922324
     341      -0.03101777
     342      -0.26340043
     343      -0.33556300
     344       0.31610856
     345       0.01850472
     346      -0.06853137
     347       0.05290624
     348       0.23635097
     349       0.03298883
     350      -0.29562384
     351       0.17200500
     352       0.01257632
     353       0.70646249
     354      -0.11107943
     355      -0.29532693
     356      -0.25561220
     357       0.18974474
     358      -1.22352454
     359       0.03125597
     360      -0.70697906
     361       0.00244145
     362       0.45620357
     363       0.14721147
     364      -0.29257899
     365      -0.35378711
     366       0.00270010
     367      -0.01687639
     368      -0.11866014
     369       0.01631624
     370      -0.29641873
     371       0.06181292
     372       1.75678865
     373       0.00136034
     374      -0.25792029
     375       0.00247778
     376       0.12394217
     377      -0.01956723
     378      -0.47312684
     379       0.00440219
     380      -1.08529006
     381      -0.56329244
     382      -0.24197032
     383      -1.37301473
     384       0.73484439
     385       0.04810163
     386      -0.08579003
     387      -0.02121765
     388       0.01387063
     389       0.04963996
     390       0.11084935
     391       0.06230502
     392       0.22334842
     393      -0.12267023
     394      -0.24540487
     395       0.00692503
     396      -0.01231490
     397       0.00267994
     398      -0.00203561
     399       0.10432954
     400      -0.04162980
     401      -0.00870659
     402       0.03468181
     403      -1.67588338
     404      -0.01995375
     405      -0.02405516
     406       0.01882244
     407      -0.04332528
     408      -0.72344578
     409       0.18595872
     410      -0.09451661
     411      -0.00459695
     412      -0.02575582
     413       1.19378777
     414       0.04025883
     415      -0.00706898
     416       0.02015568
     417      -0.02973140
     418      -0.48161379
     419       0.07619260
     420      -0.03598638
     421      -0.05925608
     422       0.12791780
     423       0.28462663
     424      -0.07479634
     425       0.03617235
     426       0.04286463
     427      -0.04499727
     428      -0.19342812
     429       0.04818802
     430      -0.00051232
     431       0.08021913
     432       0.09273137
     433      -0.32651048
     434       0.04649552
     435      -0.03451547
     436       0.00195308
     437      -0.00118793
     438       0.56221666
     439      -0.00652823
     440       0.00446485
     441      -0.00033934
     442      -0.02501334
     443      -0.12999933
     444      -0.09007960
     445       0.11900787
     446      -0.02425866
     447       0.01728595
     448       0.26967782
     449       0.01157672
     450      -0.01056419
     451      -0.00521567
     452       0.00127367
     453      -0.15183944
     454       0.00891259
     455      -0.00639150
     456      -0.00611036
     457       0.00184604
     458       0.13078716
     459       0.10410771
     460      -0.04108761
     461       0.09577794
     462      -0.52444406
     463       0.25693996
     464       0.01161097
     465      -0.03374821
     466       0.03142551
     467      -0.21246471
     468       0.07800060
     469       0.02454314
     470      -0.12361103
     471       0.08446042
     472      -0.32647616
     473       0.14864473
     474       0.00354933
     475      -0.00702255
     476       0.00184505
     477      -0.57267435
     478       0.26331795
     479       0.15849189
     480       0.31614077
     481       0.32745880
     482       0.10654593
     483      -0.30117793
     484      -0.17733956
     485      -0.01024906
     486       0.06436998
     487      -0.17300759
     488       0.13328986
     489       0.35971531
     490      -0.10545269
     491      -1.40170815
     492       0.17148381
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.02298442
       2      -0.11937261
       3       0.02686569
       4      -0.03290114
       5      -0.00302977
       6       0.04640970
       7      -0.03702844
       8       0.09487895
       9       0.12689267
      10       0.00711502
      11       0.05782531
      12       0.02032826
      13       0.00284071
      14       0.03261565
      15       0.03126685
      16      -0.25379439
      17      -0.04164022
      18      -0.00036983
      19       0.07661993
      20       0.23083321
      21      -0.25721012
      22      -0.23864527
      23      -0.12289638
      24       1.35718608
      25      -0.33464323
      26       0.23788951
      27      -0.56561491
      28       0.07730375
      29      -0.21054129
      30      -0.18236003
      31      -0.35942852
      32       0.02632355
      33       0.15198650
      34       0.03381378
      35       0.03854655
      36       0.12030947
      37      -0.02037858
      38       0.11243404
      39       0.09237424
      40      -0.58907616
      41      -1.01224892
      42       0.04223113
      43      -0.34728804
      44      -0.23719020
      45      -0.39208216
      46      -0.52165453
      47      -1.13592658
      48       0.06333123
      49      -0.10156707
      50      -0.10972510
      51      -1.33451086
      52       0.06061697
      53       0.02686681
      54       0.10494822
      55       0.02777310
      56      -0.13214628
      57       0.15993322
      58       0.09785237
      59       0.02247223
      60       0.11009839
      61       0.08296452
      62      -0.22840060
      63      -0.56481554
      64       0.19781011
      65       0.07253151
      66       0.11114448
      67       0.15658517
      68       0.14394491
      69       0.12343180
      70       0.04751550
      71      -0.11378380
      72       0.03810596
      73      -0.08660066
      74       0.11015004
      75      -0.02363653
      76      -0.14856303
      77      -0.03012222
      78      -0.02502200
      79       0.28357059
      80      -0.83813346
      81       0.07284346
      82      -1.20556874
      83      -0.15367235
      84      -0.12292522
      85       0.03362425
      86      -0.29117219
      87      -0.04086375
      88      -0.16985844
      89      -0.05625600
      90      -0.27108885
      91       0.07711394
      92      -0.91912285
      93      -0.01164662
      94      -1.22710672
      95       0.00629117
      96      -0.02573528
      97       0.03950316
      98      -0.24234692
      99      -0.01687446
     100      -0.00643108
     101       0.15902540
     102       0.10526310
     103       0.00652458
     104       0.01991715
     105       0.01347209
     106       0.01486689
     107       0.01290107
     108      -0.04425679
     109       0.00516669
     110       0.07320021
     111      -0.03894941
     112       0.00265779
     113       0.01133264
     114      -0.05410953
     115       0.00447551
     116       0.00569608
     117       0.01423779
     118      -0.02647587
     119       0.49997936
     120      -0.66040768
     121       0.02610026
     122      -1.32375614
     123      -0.00151687
     124      -0.27084591
     125       0.00411741
     126      -0.39367753
     127       0.01439844
     128      -0.44686724
     129      -0.03914108
     130      -0.05902735
     131       0.92144502
     132       0.19548646
     133      -0.76583432
     134       0.32493301
     135      -0.17686756
     136       0.10652661
     137      -0.27607591
     138      -0.82614369
     139       0.32089647
     140       0.50176428
     141       0.21156163
     142       0.30676517
     143       1.08275201
     144       0.67976950
     145      -0.43689662
     146       0.25243967
     147      -0.69886940
     148       0.16494233
     149       0.35439138
     150      -0.29376529
     151       0.08826163
     152      -0.21876802
     153       0.14820021
     154       0.01278836
     155      -0.05613286
     156       0.01480878
     157      -0.02174387
     158      -0.03449857
     159       0.05193146
     160       0.03128576
     161       0.06092385
     162      -1.15477930
     163       0.03606084
     164      -1.42885977
     165       0.01870504
     166      -0.02017983
     167       0.11733783
     168      -0.28575832
     169       0.00047760
     170       0.01811917
     171       0.08390155
     172      -0.18335914
     173       1.18669944
     174       0.72637070
     175      -0.61734189
     176       0.38438395
     177      -0.49618228
     178       0.26630365
     179       0.58660851
     180      -0.68937663
     181       0.17057325
     182      -0.27401310
     183       0.30544171
     184       0.18153387
     185       1.18655152
     186       0.49165008
     187      -0.94182928
     188       0.43139958
     189      -0.20082474
     190       0.22414659
     191       0.43324536
     192      -0.69210518
     193       0.15321098
     194      -0.30468279
     195       0.20982784
     196      -0.05965081
     197      -0.03384493
     198       0.01707436
     199       0.04854142
     200      -0.01809137
     201       0.01331703
     202      -0.05565938
     203       0.02937917
     204       0.02877073
     205       0.05377915
     206      -0.10734149
     207      -0.05198580
     208       0.00346578
     209      -0.01544591
     210      -0.32781637
     211      -0.04091400
     212       0.03069072
     213       0.00198661
     214       0.06790739
     215       0.01500712
     216       0.00453467
     217       0.01640136
     218      -0.00894319
     219      -1.14373787
     220       0.07608982
     221       0.14930115
     222      -0.06407533
     223       0.04281802
     224       0.06035852
     225      -0.07642742
     226       0.01760943
     227      -0.04072969
     228       0.00678888
     229      -0.04874955
     230      -0.06742455
     231      -0.03823373
     232       0.13921550
     233       0.05258912
     234      -0.05195054
     235       0.12848446
     236      -0.06210832
     237       0.04573632
     238      -0.07742653
     239       0.05288072
     240      -0.01364271
     241      -0.16902794
     242      -1.42035833
     243      -0.04773199
     244       0.09994820
     245      -0.02120785
     246       0.14625942
     247      -0.02939340
     248       0.11826153
     249       0.03292525
     250       0.01038538
     251       0.05428639
     252      -0.03791441
     253      -0.04493639
     254      -0.39446968
     255      -0.04629395
     256      -0.00541134
     257      -0.01646603
     258      -0.00509012
     259      -0.21471808
     260      -0.02658328
     261      -0.01000537
     262      -0.01738838
     263       0.00120053
     264      -0.01745040
     265       0.11829906
     266      -1.30850851
     267       0.24462832
     268       0.56266389
     269      -0.34122287
     270       1.05301447
     271       0.12678209
     272      -0.95187108
     273      -0.91184962
     274      -0.38722179
     275       1.65418278
     276       1.40494541
     277      -0.04092282
     278       0.23216343
     279       0.09660482
     280       0.55616172
     281      -0.24925665
     282      -0.90684050
     283      -0.00918932
     284       0.26317532
     285      -0.11853554
     286      -0.27585621
     287       0.07997688
     288       0.51074737
     289       0.31192186
     290       0.28878015
     291       0.21392288
     292       1.03816060
     293      -0.01744976
     294       0.67259884
     295      -0.01668085
     296      -0.67993114
     297       0.10918700
     298       0.58964036
     299      -0.21115058
     300      -0.19559746
     301       0.05317414
     302       0.16187033
     303       0.02995933
     304       0.51788918
     305       0.18428518
     306      -1.76571507
     307       0.01861829
     308       0.42213629
     309      -0.04637868
     310      -0.23817362
     311      -0.15846973
     312       0.71064425
     313      -0.26264135
     314       0.98551165
     315      -0.07325646
     316       0.06015023
     317       0.03512355
     318      -0.00442926
     319       0.04924220
     320      -0.02605234
     321       0.06651861
     322      -0.10421800
     323      -0.02408159
     324       0.05334418
     325      -0.02405697
     326      -0.02508310
     327       0.01701992
     328      -0.00226028
     329      -0.11874955
     330       0.27938293
     331       0.18471505
     332       0.74003465
     333       0.01470947
     334      -0.77857154
     335       0.24687711
     336      -0.11944163
     337       0.18564085
     338       1.84093671
     339      -0.07253903
     340       0.19162237
     341       0.04155501
     342       0.56606011
     343      -0.10958060
     344      -0.69526075
     345       0.01473255
     346       0.20646650
     347      -0.07690047
     348      -0.47693977
     349      -0.00805573
     350       0.34345371
     351       0.17864517
     352       0.18153051
     353      -0.64376221
     354       1.23748965
     355       0.34431327
     356      -0.10687752
     357       0.16020849
     358       1.29073050
     359       0.04416697
     360       0.85452637
     361      -0.03608179
     362      -0.90870529
     363       0.05743461
     364       0.71031212
     365      -0.17641433
     366      -0.32449813
     367       0.01355040
     368       0.23922275
     369       0.00513776
     370       0.65401526
     371       0.25361003
     372      -1.84459276
     373       0.00482113
     374       0.44557975
     375      -0.00655459
     376      -0.28432104
     377      -0.01052431
     378       0.84796260
     379      -0.13919113
     380       1.00490200
     381       0.65376254
     382      -0.24056566
     383       1.52148428
     384      -0.36218448
     385      -0.16905348
     386       0.10766834
     387       0.03435812
     388      -0.08493696
     389       0.02502558
     390      -0.01206937
     391       0.03511592
     392      -0.32100758
     393      -0.02881006
     394       0.16842090
     395      -0.03491590
     396      -0.09404402
     397      -0.00038741
     398       0.07627109
     399       0.22529614
     400      -0.16711538
     401       0.01574775
     402      -0.04349756
     403       2.35999236
     404      -0.05467549
     405       0.08672100
     406       0.02802238
     407       0.01884032
     408       1.22183332
     409       0.00806501
     410       0.04126027
     411      -0.01216319
     412      -0.03958205
     413      -1.28678809
     414       0.08692431
     415      -0.11126165
     416       0.02373832
     417       0.00625606
     418       0.77517849
     419      -0.14637013
     420       0.15920472
     421      -0.09568423
     422       0.14604042
     423      -0.86681494
     424      -0.03727328
     425       0.00904618
     426       0.04782471
     427       0.01624409
     428       0.39954219
     429       0.02569326
     430      -0.01208975
     431       0.07248348
     432      -0.05200121
     433       0.87791029
     434      -0.07784769
     435       0.04023408
     436       0.03096034
     437      -0.09357003
     438      -1.16099155
     439       0.00454525
     440      -0.01424632
     441      -0.00334963
     442       0.00829427
     443       0.27340683
     444       0.02340048
     445       0.03637246
     446      -0.05706283
     447       0.05015601
     448      -0.36181360
     449       0.02432749
     450      -0.04145101
     451      -0.01073633
     452      -0.00291448
     453       0.46919604
     454      -0.02348811
     455       0.00309660
     456       0.01492168
     457       0.01962302
     458       0.11310904
     459       0.22956747
     460       0.19489713
     461       0.02449212
     462       1.13812403
     463      -0.17452166
     464      -0.06313922
     465       0.07120732
     466      -0.06631957
     467       0.68328173
     468      -0.26779523
     469       0.08560635
     470       0.01687402
     471       0.08559710
     472       1.06704460
     473      -0.52647504
     474      -0.04864015
     475       0.05558140
     476      -0.05462864
     477       0.97450425
     478      -0.30019599
     479       0.17543706
     480       0.06747158
     481       0.09591741
     482       0.05127713
     483      -0.09948673
     484      -0.02485203
     485       0.00197589
     486       0.01579563
     487      -0.18729091
     488       0.03709850
     489      -0.04846289
     490       0.08959938
     491       1.33964905
     492      -0.07503407
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.13072952
       2       0.09259470
       3      -0.02354642
       4       0.34285900
       5       0.37046943
       6       0.46937448
       7       0.27115746
       8      -0.27542078
       9      -0.32034506
      10       0.18151949
      11       0.07111693
      12       0.15996592
      13       0.04635344
      14       0.07224205
      15      -0.14523806
      16      -0.07790380
      17      -0.28834766
      18       0.36053332
      19       0.20120433
      20       0.17742672
      21       0.43460714
      22       0.56311008
      23       0.13757018
      24      -0.53261626
      25       0.87635433
      26       0.48536370
      27       0.81256380
      28       0.33803602
      29       0.45816638
      30       0.79843257
      31      -0.15499510
      32       0.15422795
      33      -0.00687340
      34       0.25039226
      35       0.05654098
      36       0.29310658
      37       0.38398842
      38      -0.33607362
      39      -0.28674468
      40       0.19543628
      41       1.22625554
      42       0.63468261
      43       0.60584223
      44       0.83471884
      45      -0.15223203
      46       0.84731808
      47      -0.02324653
      48      -0.27381353
      49       0.15705609
      50      -0.72181719
      51       3.92989930
      52       0.09279448
      53       0.31843447
      54      -0.35197357
      55      -0.11821869
      56       0.17143072
      57       0.17973518
      58       0.14956061
      59       0.16303687
      60       0.20323178
      61      -0.32312058
      62      -0.18361736
      63       0.32490446
      64       0.49620698
      65       0.26201630
      66       0.37127498
      67       0.31468578
      68      -0.61782181
      69      -0.39285757
      70       0.14036474
      71       0.10563608
      72      -0.01193838
      73       0.01271428
      74      -0.05365795
      75       0.07884958
      76       0.11417716
      77       0.13824763
      78       0.08819650
      79       0.20372534
      80       0.01180795
      81      -0.02863626
      82       0.45174225
      83       0.05218130
      84       0.31029748
      85       0.02725971
      86       0.16951350
      87      -0.04312420
      88       0.43241362
      89       0.45910615
      90       0.82156153
      91      -0.44530889
      92       0.64357875
      93      -0.07554211
      94       1.64478915
      95       0.06627488
      96       0.32145470
      97      -0.10559390
      98       0.52095102
      99       0.09199351
     100       0.62678370
     101      -0.33445937
     102       0.06912781
     103      -0.02789225
     104      -0.00743939
     105      -0.07366549
     106      -0.03118413
     107       0.08334325
     108       0.02336556
     109      -0.10773495
     110      -0.18510937
     111       0.05694594
     112      -0.00125578
     113      -0.02996018
     114       0.18881545
     115       0.03570705
     116      -0.01429139
     117      -0.05977546
     118       0.07124098
     119       0.56831662
     120       0.28580127
     121       0.05057945
     122       0.82769589
     123       0.00957753
     124       0.34738281
     125       0.12299590
     126       0.15550039
     127      -0.04312707
     128       0.72787934
     129       0.23851763
     130       0.86743567
     131       0.05114686
     132      -0.31431940
     133       0.38034265
     134      -0.66853587
     135       0.59676244
     136      -0.12111297
     137      -0.43101235
     138       0.84665063
     139      -0.46558487
     140       0.41721462
     141      -0.28564162
     142      -0.58954845
     143      -0.84727921
     144      -0.30989025
     145       0.13223207
     146      -0.26148145
     147       0.81444244
     148      -0.33383607
     149      -0.59503009
     150       1.25685112
     151      -0.33405652
     152       0.39686473
     153      -0.57013072
     154      -0.83213957
     155      -0.34255823
     156       0.04308759
     157       0.11386050
     158       0.04407741
     159      -0.03234103
     160      -0.12351956
     161      -0.33636025
     162       0.58821004
     163       0.04826185
     164       1.29534233
     165      -0.00673066
     166       0.30189816
     167      -0.06835623
     168       0.52494430
     169      -0.02020220
     170       0.38646304
     171      -0.13642345
     172       0.76718918
     173      -0.41371160
     174      -0.02961774
     175      -0.18970698
     176      -0.01000708
     177       0.47371241
     178      -0.44371010
     179      -0.74632046
     180       0.98184648
     181      -0.37484398
     182       0.27444792
     183      -0.54514430
     184      -0.54180552
     185      -1.30686685
     186      -0.83835779
     187       1.22312668
     188      -0.76005679
     189       0.63941314
     190      -0.49717893
     191      -0.84223492
     192       2.13844583
     193      -0.61090083
     194       0.18952165
     195      -0.82146973
     196      -0.04503385
     197       0.13784272
     198      -0.05201596
     199      -0.16402318
     200       0.02792955
     201      -0.03656515
     202       0.06736926
     203      -0.02576386
     204      -0.03026513
     205      -0.13505952
     206       0.06986098
     207      -0.01134392
     208       0.08228179
     209      -0.13336347
     210       0.10187340
     211       0.02546901
     212      -0.02071130
     213      -0.00273629
     214       0.16973007
     215       0.12062887
     216      -0.01306870
     217      -0.09413599
     218       0.06138890
     219      -0.91437705
     220       0.02954362
     221      -0.06894847
     222       0.08785677
     223      -0.13940787
     224      -0.34554596
     225       0.34468861
     226      -0.12402326
     227      -0.21573314
     228      -0.00652580
     229       0.04872045
     230       0.00835588
     231      -0.01005664
     232      -0.02768112
     233       0.02603663
     234       0.01887372
     235      -0.02739987
     236       1.15096383
     237      -0.07497206
     238      -0.43068630
     239      -0.01077594
     240      -0.21839856
     241       0.31584857
     242       0.41663565
     243       0.24298481
     244       0.41241048
     245       0.06566786
     246       0.06914916
     247       0.13543945
     248       0.15880160
     249      -0.07934016
     250      -0.28106122
     251      -0.22212023
     252      -0.31558665
     253       0.06156130
     254       0.57921902
     255       0.11290683
     256       0.28561699
     257      -0.02120045
     258       0.17030265
     259       0.41256112
     260       0.05990674
     261       0.01220521
     262       0.13053940
     263      -0.02597102
     264      -0.03689741
     265      -0.10959371
     266       0.67710999
     267       0.49061607
     268      -0.61076978
     269      -0.61621890
     270      -0.02403796
     271       0.01406902
     272       0.20836160
     273      -0.02954090
     274      -0.34439260
     275       0.41732576
     276       0.67266490
     277       0.26749834
     278      -0.47921898
     279       0.09698754
     280      -0.43009111
     281       0.46372098
     282       0.59583216
     283       0.16993246
     284      -0.36843135
     285      -0.04400231
     286       0.15491629
     287      -0.48336260
     288      -0.28017835
     289      -0.94643304
     290      -0.23648361
     291      -0.98738339
     292      -0.31973341
     293       0.31643415
     294      -0.25489413
     295      -0.13265794
     296       1.28851992
     297      -0.31321054
     298      -0.56919882
     299       0.45519939
     300       0.19926553
     301      -0.07810519
     302      -0.30871672
     303      -0.06127968
     304      -0.96724657
     305      -0.27526999
     306       0.64085275
     307      -0.05771963
     308      -0.78660200
     309       0.15535459
     310       0.69952223
     311       0.33883051
     312      -0.95234986
     313       0.38189906
     314      -0.02893565
     315       0.08500978
     316      -0.00227527
     317       0.07131160
     318      -0.03723509
     319      -0.08544590
     320      -0.07244783
     321      -0.31008737
     322      -0.09215682
     323       0.05940738
     324      -0.01506974
     325       0.10500708
     326       0.11183171
     327      -0.04821418
     328       0.02778563
     329      -0.09686551
     330      -1.04246625
     331       0.07338845
     332      -0.48778244
     333       0.02235171
     334       0.45532968
     335       0.15162030
     336      -0.48545641
     337      -0.79330730
     338       0.34857705
     339       0.01065477
     340      -0.46250414
     341       0.00452775
     342      -0.58584700
     343       0.67751319
     344       0.65718968
     345      -0.03315125
     346      -0.37920515
     347      -0.05499522
     348       0.21644105
     349      -0.00673797
     350      -0.17822838
     351      -0.52044931
     352      -0.12409169
     353      -0.73818544
     354      -0.16952961
     355      -0.11090845
     356       0.41395723
     357      -0.34232308
     358      -0.46867727
     359      -0.08915184
     360      -0.27716022
     361       0.04280453
     362       1.27402619
     363      -0.30913465
     364      -0.60436891
     365       0.73448931
     366       0.51113252
     367       0.04690356
     368      -0.41387182
     369      -0.04892479
     370      -1.06510773
     371      -0.49020323
     372       0.61049339
     373      -0.02578026
     374      -0.81274313
     375       0.00284922
     376       0.56074002
     377       0.08737111
     378      -1.16874291
     379       0.15356882
     380       0.06006250
     381      -0.35435815
     382       0.55465689
     383      -0.47439646
     384      -0.88304979
     385       0.09596684
     386       0.06418092
     387       0.00648804
     388      -0.08096363
     389      -0.11080008
     390      -0.32152385
     391      -0.17541649
     392      -0.00887593
     393       0.26156048
     394       0.22203017
     395       0.03238622
     396       0.27270126
     397      -0.01136089
     398      -0.15895395
     399      -0.18627377
     400       0.07705652
     401      -0.00512583
     402      -0.08348479
     403      -0.85257342
     404       0.06508914
     405       0.00621258
     406      -0.00269132
     407      -0.01527728
     408      -0.46244899
     409      -0.37937438
     410       0.20096603
     411       0.01040521
     412       0.12764369
     413       0.26202663
     414      -0.18235204
     415       0.17073202
     416      -0.02420603
     417      -0.00028481
     418      -0.54940329
     419       0.01637082
     420      -0.11147211
     421       0.05384418
     422      -0.26562252
     423       1.19583669
     424       0.20694613
     425      -0.11037611
     426      -0.09609430
     427       0.13360548
     428      -0.36284729
     429      -0.12729892
     430       0.00720260
     431      -0.17102149
     432      -0.18004104
     433      -0.96242975
     434       0.08432829
     435      -0.04881506
     436       0.00789240
     437       0.09114240
     438       2.33406242
     439       0.02088968
     440       0.00075253
     441      -0.00462629
     442       0.07302506
     443      -0.68350415
     444       0.13617879
     445      -0.24966771
     446       0.04773817
     447      -0.11589170
     448       0.11814608
     449      -0.16219895
     450       0.15353762
     451       0.03318456
     452       0.00443587
     453      -1.16757716
     454      -0.02287496
     455       0.03423468
     456       0.00287311
     457      -0.01642986
     458      -1.18792473
     459      -0.28126525
     460      -0.03033831
     461       0.05675944
     462      -1.07662619
     463       0.10195487
     464      -0.05129364
     465      -0.12629565
     466       0.00576769
     467      -1.99089724
     468       0.64359837
     469       0.00265824
     470       0.20639653
     471      -0.16628601
     472      -1.57228205
     473       0.59231829
     474       0.01553800
     475      -0.05439319
     476       0.01685069
     477      -1.71213864
     478       0.28066204
     479      -0.32806764
     480      -0.76531755
     481      -0.62499910
     482      -0.36045542
     483       0.72343146
     484       0.30433336
     485       0.09846814
     486      -0.22240361
     487       0.40831575
     488      -0.31405064
     489      -0.66708710
     490       0.10870342
     491       0.36342494
     492      -0.37008103
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.03738403
       2      -0.13648882
       3      -0.12329264
       4      -0.07016123
       5      -0.20349465
       6      -0.12798399
       7      -0.10067985
       8       0.05945522
       9       0.17691071
      10       0.03573827
      11      -0.10299109
      12      -0.12469355
      13      -0.02813496
      14      -0.09981355
      15       0.18951767
      16      -0.12548554
      17      -0.37505107
      18      -0.00607801
      19      -0.15552752
      20      -0.34366449
      21      -0.34374331
      22       0.13778447
      23       0.03549144
      24      -1.12317143
      25      -0.35041099
      26      -0.53998390
      27      -0.16583825
      28      -0.39314638
      29      -0.12157251
      30      -0.17210465
      31      -0.01351217
      32      -0.05421301
      33      -0.31127398
      34      -0.18444580
      35      -0.11946459
      36      -0.04002609
      37       0.00863911
      38       0.03131841
      39       0.18269014
      40      -0.27760430
      41       0.26676067
      42      -0.44894952
      43       0.06168847
      44      -0.23687201
      45      -0.23446284
      46      -0.02882068
      47      -0.21174313
      48      -0.02891875
      49      -0.48324672
      50       0.76968972
      51      -6.40550691
      52      -0.11755510
      53      -0.21833873
      54       0.20230470
      55       0.10387087
      56       0.21040949
      57      -0.10541856
      58      -0.13180589
      59       0.01618108
      60      -0.16723630
      61       0.26771881
      62      -0.41469820
      63       0.25917532
      64      -0.31442753
      65      -0.19366771
      66      -0.08142454
      67      -0.26450861
      68       0.53660319
      69       0.33116792
      70      -0.12675969
      71       0.01916425
      72      -0.07129459
      73       0.03926180
      74      -0.08143981
      75      -0.05363056
      76       0.08594517
      77      -0.10334274
      78      -0.06581153
      79       0.11447770
      80      -0.12337786
      81      -0.00753033
      82      -0.08754419
      83       0.06981669
      84      -0.10683332
      85       0.02324380
      86       0.07442036
      87       0.03648435
      88      -0.53041308
      89      -0.07638196
      90      -1.35736093
      91       0.01094992
      92      -0.19340098
      93       0.07682325
      94      -1.56253592
      95      -0.08501169
      96      -0.55969797
      97       0.05090567
      98      -0.46696714
      99      -0.05896719
     100      -1.19824446
     101       0.09499969
     102      -1.28439890
     103      -0.00940508
     104       0.00025571
     105       0.01108848
     106      -0.03460948
     107      -0.01678318
     108       0.01668024
     109       0.08151638
     110      -0.02297588
     111       0.05451712
     112      -0.03913653
     113       0.01036910
     114      -0.14009869
     115      -0.01919286
     116       0.01526344
     117       0.07891563
     118      -0.09904878
     119       0.36460398
     120      -0.00875380
     121       0.00400376
     122      -0.06515277
     123      -0.02174411
     124      -0.02538593
     125       0.02454707
     126       0.08238717
     127      -0.00364524
     128      -0.57085327
     129       0.04840721
     130      -1.76536619
     131       0.37591083
     132      -0.40675767
     133       0.48685263
     134      -0.70729308
     135       0.05231085
     136       0.00619723
     137       0.42586428
     138      -0.38331994
     139       0.20345815
     140      -0.37099369
     141       0.16498426
     142       0.63159855
     143       0.18823580
     144      -0.09422468
     145      -0.06740033
     146       0.26188768
     147      -0.42471674
     148       0.18670561
     149       0.36650728
     150      -2.13628574
     151       0.65093863
     152      -0.01439510
     153       0.87350853
     154       1.59118856
     155      -0.04526543
     156      -0.02912913
     157      -0.01569769
     158      -0.00714484
     159       0.01812226
     160       0.00745816
     161      -0.16750319
     162       0.08855329
     163      -0.08611096
     164      -1.15003067
     165      -0.04420518
     166      -0.55097568
     167      -0.06263815
     168      -0.41643370
     169       0.04421561
     170      -0.95525157
     171      -0.04825444
     172      -2.06369773
     173      -0.11868344
     174      -0.12950741
     175      -0.00201676
     176       0.16535589
     177      -0.31304567
     178       0.15327921
     179       0.27039904
     180      -1.35505813
     181       0.63983176
     182       0.12855188
     183       0.63647092
     184       0.95385554
     185       1.17108886
     186       0.68063262
     187      -0.93129878
     188       0.66190534
     189      -1.19176807
     190       0.75546311
     191       1.29163404
     192      -3.92175055
     193       1.42754297
     194      -0.26344810
     195       1.74959384
     196       0.83623825
     197       0.05777646
     198      -0.01938699
     199       0.12574193
     200      -0.05709543
     201       0.07004495
     202      -0.00108111
     203       0.00585657
     204      -0.02206617
     205      -0.08750785
     206       0.02888758
     207       0.06094505
     208      -0.00787667
     209       0.04333073
     210       0.15931586
     211       0.06685638
     212      -0.05621738
     213      -0.00522542
     214      -0.22436944
     215       0.02032544
     216      -0.00700534
     217       0.07481344
     218      -0.05267725
     219      -0.08378732
     220      -0.06352831
     221      -0.01338038
     222      -0.00955905
     223       0.04040265
     224       0.01327310
     225       0.09592510
     226       0.05686944
     227       0.25578765
     228      -0.02099178
     229       0.03307742
     230       0.21653327
     231       0.02297344
     232      -0.10637819
     233      -0.05939182
     234       0.04587963
     235      -0.23236501
     236       0.12966802
     237      -0.09523215
     238       0.21879128
     239      -0.10519238
     240      -0.03872115
     241       0.19548847
     242      -0.20190022
     243       0.02286059
     244      -0.31035397
     245       0.01254993
     246      -0.15058407
     247      -0.01704545
     248      -0.56584205
     249      -0.02282011
     250       0.15843379
     251      -0.07348911
     252       0.10136217
     253       0.03296264
     254      -0.25583022
     255       0.02535691
     256      -0.12592892
     257       0.03954591
     258      -0.00076149
     259       0.03367405
     260      -0.05300693
     261       0.01783417
     262      -0.05952245
     263       0.03217209
     264       0.11441402
     265      -0.17687380
     266      -0.36713026
     267       0.73125873
     268      -0.00241239
     269      -0.55490698
     270      -0.44757388
     271      -0.02857478
     272       0.07391575
     273       0.45139583
     274       0.33202876
     275      -0.88245968
     276      -0.95361017
     277       0.01427126
     278       0.18739490
     279      -0.02381914
     280       0.24362233
     281       0.10281059
     282      -0.27002759
     283      -0.09831576
     284       0.45647907
     285       0.09082565
     286      -0.04830123
     287       0.06383463
     288       0.07488736
     289       0.12318688
     290       0.08599849
     291      -0.12244164
     292       0.11663534
     293      -0.01107100
     294      -0.07860650
     295       0.16486881
     296      -1.25042270
     297       0.07026288
     298       0.22641575
     299       0.03418699
     300      -0.11161574
     301      -0.04452277
     302       0.50566433
     303       0.05890291
     304       1.09653212
     305      -0.14793347
     306      -0.13132279
     307      -0.02257964
     308       1.19818894
     309      -0.10304441
     310      -0.78218469
     311      -0.02917199
     312       1.08638152
     313       0.21452052
     314      -0.65227198
     315       0.05251923
     316      -0.10560357
     317      -0.05112993
     318      -0.00149707
     319      -0.07936502
     320       0.05868619
     321       0.05693707
     322       0.21078574
     323       0.01219659
     324      -0.08132257
     325      -0.01010209
     326      -0.00989970
     327       0.01701033
     328      -0.08648151
     329       0.13376266
     330      -0.04135425
     331       0.11276519
     332      -0.40559876
     333      -0.03167900
     334      -0.02062940
     335       0.13409777
     336       0.27580656
     337      -0.73321130
     338      -1.32685866
     339       0.04056682
     340       0.30163389
     341      -0.01987273
     342       0.19992905
     343      -0.01277306
     344      -0.29664705
     345      -0.00477602
     346       0.47826838
     347       0.07578165
     348       0.00933998
     349      -0.02270530
     350       0.20496561
     351       0.07738301
     352      -0.02263891
     353       0.38650132
     354      -1.12228370
     355       0.18558453
     356       0.05038433
     357      -0.17692420
     358       0.08314172
     359      -0.05907867
     360      -0.29287471
     361      -0.00696329
     362      -0.92300945
     363      -0.02687922
     364       0.17784062
     365       0.10952379
     366      -0.22760144
     367      -0.09885423
     368       0.50318577
     369       0.03948136
     370       1.04881642
     371      -0.09487010
     372      -0.33860441
     373      -0.00849263
     374       1.27294229
     375       0.04570942
     376      -0.54797961
     377      -0.10987717
     378       1.14537368
     379       0.18873277
     380      -0.67073276
     381       0.28207462
     382       0.04933002
     383      -0.24600669
     384       1.13685364
     385       0.17698035
     386      -0.14943054
     387      -0.03869580
     388       0.09035614
     389      -0.08155194
     390       0.04559851
     391      -0.00547643
     392       0.31775618
     393       0.03313045
     394      -0.21320147
     395       0.08959145
     396      -0.01146347
     397       0.01040615
     398       0.04519617
     399      -0.15996323
     400       0.06464899
     401       0.01524636
     402       0.04837571
     403       0.08565204
     404       0.07298946
     405      -0.03878243
     406      -0.02222996
     407       0.04490136
     408      -0.23663462
     409       0.00677435
     410      -0.12454835
     411       0.03234528
     412      -0.00856858
     413       0.19604700
     414      -0.03370663
     415       0.07988030
     416      -0.03627085
     417       0.04285920
     418       0.13042860
     419       0.08377981
     420      -0.15452151
     421       0.11422560
     422      -0.06631202
     423      -0.89729703
     424      -0.14007497
     425       0.10762277
     426      -0.02197361
     427      -0.11458022
     428       0.18415304
     429       0.11021489
     430      -0.08932409
     431      -0.04742792
     432       0.13808709
     433       0.54804964
     434       0.00819151
     435       0.03185846
     436      -0.06308769
     437       0.15193476
     438      -3.19346372
     439      -0.04780114
     440       0.03999488
     441       0.01297427
     442      -0.11400077
     443       1.25232709
     444      -0.10982758
     445       0.03768945
     446       0.05679479
     447       0.02169381
     448       0.05159028
     449       0.18855111
     450      -0.10619912
     451      -0.00063193
     452      -0.04569069
     453       1.52519002
     454       0.04395626
     455      -0.01154751
     456      -0.05332815
     457      -0.02087953
     458       2.47352309
     459      -0.19755125
     460      -0.13315294
     461       0.04412893
     462       0.63744379
     463      -0.50452094
     464       0.22145030
     465       0.03477617
     466       0.05694359
     467       3.05268450
     468      -0.80934764
     469      -0.06628918
     470      -0.04452319
     471      -0.05806906
     472       1.00227341
     473      -0.29288493
     474       0.09627993
     475      -0.09192738
     476       0.14407591
     477       2.35650757
     478      -0.41072440
     479      -0.20762595
     480       0.37063882
     481       0.21893530
     482      -0.03589384
     483      -0.17690423
     484      -0.22671015
     485      -0.03349009
     486       0.14047803
     487       0.15401554
     488       0.14953814
     489       0.28413416
     490      -0.06027683
     491      -0.36005370
     492       0.18391633

 Residual norm when dim(red L) =   8
 NEO root     CSF        orbital          total
    1     0.01947279     0.00353513     0.01979108 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  1

   1.987527862   1.975877057   0.024687245   0.013366498

 Symmetry  2

   1.980589323   0.018582169

 Symmetry  3

   1.973942019   0.025427827


 <<< MACRO ITERATION  2 >>>
 --------------------------

 Total MCSCF energy       :      -76.168611452169159       (MACRO    2)
 - Dielec. solvation ener.:       -0.012552130625138

 Norm of total gradient   :        0.057436342990
 -    of CI gradient      :        0.031045877950
 -    of orbital gradient :        0.048322737488
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.03889835
       2      -0.54183428
       3       0.05858407
       4      -0.07169949
       5      -0.06395819
       6      -0.09732073
       7      -0.04549559
       8       0.05287097
       9       0.06252632
      10      -0.05224479
      11      -0.00310124
      12      -0.02738756
      13      -0.00606289
      14       0.00288031
      15      -0.00158822
      16       0.00865133
      17      -0.13537580
      18      -0.07263788
      19      -0.01819827
      20      -1.03662180
      21      -1.77041255
      22      -0.13168837
      23      -0.17968220
      24      -0.86813281
      25      -0.16490394
      26      -0.90332512
      27      -0.26179705
      28      -0.00852316
      29      -0.11121821
      30      -0.17197003
      31       0.01334281
      32      -0.09360311
      33      -0.12605633
      34      -0.01776831
      35       0.01505443
      36      -0.05370711
      37      -0.06679513
      38       0.06910102
      39       0.04768027
      40      -1.52689056
      41      -0.44806958
      42      -0.10392325
      43      -0.17171703
      44      -0.17174178
      45       0.05150206
      46      -0.29028494
      47      -0.10393662
      48       0.11185613
      49       0.03044944
      50      -0.06046864
      51       0.81292747
      52       0.00473284
      53      -0.02940411
      54       0.06116517
      55       0.01729436
      56       0.17021511
      57      -0.03440517
      58      -0.04679064
      59      -0.02905106
      60      -0.01936726
      61       0.03565603
      62      -0.31568826
      63      -0.15185469
      64      -0.06814365
      65      -0.01891284
      66      -0.02775455
      67      -0.02132504
      68       0.05696120
      69       0.03052017
      70      -0.01190850
      71      -0.01211386
      72      -0.04901796
      73      -0.00295720
      74      -0.01228227
      75      -0.01067624
      76      -0.03309637
      77      -0.01642207
      78      -0.00353385
      79       0.32626597
      80       0.34445570
      81       0.03309164
      82      -0.23881547
      83       0.01015059
      84      -0.03769568
      85       0.01444477
      86      -0.01643861
      87       0.02327624
      88       0.10821300
      89      -0.09953773
      90       0.21237326
      91       0.11153748
      92      -0.04892335
      93       0.01199251
      94      -0.16556517
      95       0.00143092
      96       0.05897091
      97       0.02070488
      98      -0.07589251
      99      -0.01998351
     100       0.22464759
     101       0.06355300
     102       0.58036288
     103      -0.00504700
     104      -0.00539965
     105       0.01978059
     106       0.02386635
     107      -0.02634292
     108      -0.00882551
     109       0.01211216
     110       0.07654218
     111      -0.02363007
     112       0.01500980
     113       0.00479731
     114      -0.01069713
     115      -0.00446002
     116      -0.00228749
     117      -0.00547456
     118       0.01364959
     119      -0.01025862
     120       0.02451932
     121       0.01036691
     122      -0.23596853
     123      -0.00168214
     124      -0.04390177
     125       0.00557542
     126      -0.01431033
     127      -0.00267222
     128      -0.02009792
     129      -0.05025798
     130       0.30921006
     131       0.35328389
     132      -0.30616703
     133       0.60893846
     134      -0.41256088
     135       0.11381431
     136       0.06283716
     137       1.20491044
     138      -0.11744684
     139       0.09724108
     140      -0.66962287
     141       0.02867606
     142       0.00508859
     143       0.15751250
     144       0.03279870
     145      -0.08648785
     146       0.08866682
     147      -0.04092795
     148       0.07365400
     149       0.07493368
     150       0.23223792
     151      -0.10557116
     152      -0.11153562
     153      -0.06156137
     154      -0.27843492
     155       0.04855245
     156      -0.01519047
     157      -0.03211605
     158      -0.00101760
     159      -0.00941968
     160       0.03328781
     161       0.06987457
     162      -0.04739987
     163      -0.00525110
     164      -0.29473411
     165       0.01795865
     166       0.05843504
     167       0.01666077
     168      -0.11066026
     169      -0.00549302
     170       0.20816266
     171       0.03589567
     172       0.53823736
     173       0.17807199
     174       0.04498954
     175      -0.12476526
     176      -0.02893651
     177       0.17216338
     178       0.12492246
     179       0.13026044
     180       0.01269699
     181      -0.07928290
     182      -0.09251510
     183       0.03805597
     184      -0.15047951
     185       0.23754969
     186       0.19177290
     187      -0.18122320
     188       0.02832084
     189       0.12487272
     190      -0.03245615
     191      -0.10609777
     192       0.61028943
     193      -0.29487880
     194      -0.03691735
     195      -0.30834872
     196      -0.38115244
     197      -0.03992127
     198       0.01708223
     199       0.00599298
     200       0.01796273
     201      -0.01438841
     202      -0.00550999
     203      -0.05406043
     204       0.01580464
     205       0.03253456
     206       0.04231533
     207       0.05551278
     208      -0.03757561
     209       0.03830756
     210       0.08250404
     211      -0.03444963
     212       0.01427346
     213       0.00434547
     214      -0.00365780
     215      -0.02725764
     216       0.00064759
     217       0.01076392
     218       0.00398263
     219       0.73557917
     220      -0.05085522
     221      -0.03376256
     222      -0.00819610
     223       0.03904285
     224       0.07585535
     225      -0.07853976
     226       0.03392005
     227      -0.00148002
     228       0.02335173
     229      -0.01946286
     230      -0.04016703
     231       0.00620665
     232       0.03107940
     233      -0.04977611
     234       0.11186341
     235      -0.36407291
     236      -1.52942111
     237       0.01278546
     238       0.09927004
     239      -0.01512875
     240       0.06170756
     241      -0.05311839
     242      -0.45905110
     243      -0.05888562
     244      -0.05908028
     245       0.00451713
     246       0.01936395
     247      -0.02908761
     248       0.08042729
     249       0.01595810
     250       0.04759761
     251       0.05607336
     252       0.05377474
     253      -0.01413192
     254      -0.20925456
     255      -0.01850405
     256      -0.07001281
     257       0.02067680
     258      -0.05398971
     259      -0.09954525
     260      -0.00548156
     261      -0.00613123
     262      -0.02628970
     263       0.00559227
     264      -0.00160924
     265       0.05775611
     266      -0.15764983
     267       1.01598707
     268       2.12321366
     269      -0.21990898
     270      -0.44840978
     271       0.03738034
     272      -0.17743448
     273       1.06896096
     274       1.14850175
     275      -1.20894339
     276      -2.08374696
     277      -0.12763135
     278       0.15594026
     279      -0.05289414
     280       0.11312139
     281      -0.07464600
     282      -0.20508627
     283      -0.00967674
     284       0.01289965
     285      -0.05370664
     286       0.07844305
     287       0.07135256
     288       0.09976677
     289       0.21618527
     290       0.07719689
     291       0.21762224
     292       0.34173929
     293      -0.06126416
     294       0.21602139
     295       0.00830433
     296      -0.15811901
     297       0.06106910
     298       0.13685282
     299      -0.09603324
     300      -0.07662533
     301       0.02527821
     302      -0.07236786
     303      -0.00066355
     304      -0.04730400
     305       0.10661230
     306      -0.24559798
     307       0.03870214
     308      -0.12501808
     309      -0.02920131
     310      -0.04756487
     311      -0.07600283
     312      -0.00731107
     313      -0.16721153
     314       0.28820858
     315      -0.04135888
     316       0.01476227
     317      -0.01633750
     318      -0.00812585
     319       0.01955430
     320       0.01337577
     321       0.07013192
     322      -0.04103701
     323      -0.01195278
     324       0.03630175
     325      -0.01410651
     326      -0.04552705
     327       0.00771307
     328       0.01842880
     329       0.32726755
     330       1.44504436
     331      -0.10219483
     332      -0.07280552
     333      -0.03187171
     334      -0.11915100
     335       0.15781968
     336       0.84116026
     337      -0.47639413
     338      -1.87192587
     339       0.01762462
     340       0.10112622
     341      -0.02315687
     342       0.16661731
     343      -0.13184448
     344      -0.17127623
     345       0.00941080
     346       0.01443701
     347       0.02934968
     348       0.00424532
     349       0.00736507
     350       0.04460359
     351       0.10388468
     352       0.06142265
     353       1.39304204
     354      -1.38490144
     355       0.06531883
     356      -0.10199835
     357       0.04405455
     358       0.54693541
     359      -0.00537795
     360       0.27714394
     361       0.00881073
     362      -0.28261730
     363       0.06962694
     364       0.12847796
     365      -0.14414143
     366      -0.14499261
     367       0.00825970
     368      -0.02991266
     369       0.00265187
     370      -0.00147196
     371       0.13101179
     372      -0.25505636
     373       0.02352389
     374      -0.13405296
     375      -0.02280368
     376      -0.08933427
     377       0.00330145
     378       0.07200993
     379      -0.08471879
     380       0.28451286
     381       0.19382379
     382      -0.14040517
     383       0.34612984
     384      -0.15077193
     385      -0.06689662
     386       0.01371398
     387      -0.01293918
     388       0.14604562
     389       0.01594499
     390       0.07173250
     391       0.04006126
     392      -0.08406286
     393      -0.05064299
     394       0.04256352
     395      -0.00415912
     396      -0.06638579
     397       0.00173818
     398       0.04798456
     399       0.02154833
     400      -0.00901966
     401       0.01297623
     402       0.01389331
     403       0.26494085
     404      -0.02810500
     405       0.00769277
     406      -0.00781209
     407       0.03728053
     408       0.12160231
     409       0.05828398
     410      -0.03837895
     411       0.00929623
     412      -0.02167953
     413      -0.30249955
     414       0.12552479
     415      -0.10698521
     416       0.01187856
     417      -0.01969283
     418       0.18994744
     419      -0.05821139
     420       0.07488570
     421      -0.01512690
     422       0.10991952
     423       0.09148101
     424      -0.02502556
     425       0.00814163
     426       0.01665152
     427      -0.01621482
     428       0.12311801
     429      -0.00002940
     430       0.01125030
     431       0.04926036
     432       0.02981709
     433       0.12874588
     434      -0.01642347
     435       0.01028859
     436       0.00586541
     437      -0.04634234
     438       0.17611528
     439       0.00584525
     440      -0.00760757
     441      -0.00015551
     442       0.00622400
     443      -0.19267590
     444      -0.01303007
     445       0.02880139
     446      -0.01729056
     447       0.02254120
     448      -0.09362156
     449       0.00973967
     450      -0.01687959
     451      -0.01213198
     452       0.02594474
     453       0.00001496
     454       0.00135976
     455      -0.00656109
     456       0.00321730
     457      -0.00915950
     458      -0.45653253
     459       0.11379354
     460      -0.04593962
     461      -0.03505982
     462       0.09018258
     463       0.14476679
     464      -0.02018594
     465       0.02956440
     466      -0.01386722
     467      -0.24760097
     468       0.05524229
     469      -0.01130901
     470      -0.07071251
     471       0.06041537
     472       0.18991817
     473       0.00782300
     474      -0.02683537
     475       0.02899338
     476      -0.03505286
     477      -0.05297855
     478       0.01074879
     479       0.07057116
     480       0.16474050
     481       0.11683095
     482       0.10248138
     483      -0.14636880
     484      -0.05732394
     485      -0.04696945
     486       0.07317796
     487      -0.09383926
     488       0.06607104
     489       0.13709599
     490      -0.01486880
     491       0.10068886
     492       0.01984298
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.03958416
       2       0.09560632
       3       0.08496340
       4      -0.07799508
       5      -0.01486452
       6      -0.04332511
       7      -0.08276348
       8       0.17732005
       9       0.18547492
      10       0.10520697
      11       0.03333683
      12      -0.01578675
      13      -0.01358079
      14       0.01823471
      15      -0.02434452
      16       0.33644885
      17       0.47663277
      18      -0.06804269
      19      -0.04579426
      20       0.72962464
      21       0.84744631
      22       0.44499699
      23      -0.02715178
      24       0.82826329
      25       0.52642064
      26       1.36164014
      27       0.36198560
      28       0.04292113
      29      -0.14198747
      30       0.19265690
      31       0.22664065
      32       0.11755458
      33       0.33098079
      34       0.04369374
      35       0.07029383
      36      -0.08509133
      37      -0.06749534
      38       0.20565143
      39       0.12085611
      40       0.78275071
      41       0.66756840
      42       0.00607027
      43      -0.18350841
      44       0.45799630
      45       0.28661243
      46       0.29721916
      47       0.19710783
      48       0.25393579
      49       0.23154558
      50      -0.01675598
      51       2.68497223
      52       0.04164314
      53       0.01045551
      54       0.14595579
      55       0.02059965
      56       0.35588363
      57      -0.08925514
      58      -0.10116191
      59      -0.05935684
      60      -0.09772089
      61       0.07974605
      62       0.17370391
      63       0.70106827
      64      -0.05154956
      65       0.10628731
      66      -0.01255796
      67      -0.02770904
      68       0.12747593
      69       0.06054611
      70       0.00321669
      71       0.06195938
      72      -0.09924798
      73       0.03457019
      74      -0.15176752
      75      -0.01994337
      76      -0.03313808
      77      -0.05306871
      78      -0.00091180
      79      -0.06737270
      80      -0.02748434
      81       0.02876511
      82       0.34093675
      83       0.09028830
      84      -0.01753541
      85      -0.06916202
      86      -0.08027403
      87      -0.02020729
      88      -0.33221031
      89      -0.08593265
      90       0.29126633
      91       0.13073791
      92       0.00901312
      93       0.01769083
      94      -0.16279871
      95      -0.01231833
      96       0.13236817
      97       0.04707921
      98      -0.09627374
      99      -0.04547723
     100       0.27775931
     101       0.22374864
     102       1.19256716
     103      -0.02092637
     104       0.00887221
     105       0.01534191
     106       0.01475265
     107      -0.02491567
     108       0.02604973
     109       0.00750697
     110       0.06788469
     111      -0.07475315
     112       0.02481715
     113       0.01659605
     114      -0.04638481
     115      -0.00633483
     116      -0.00002034
     117      -0.00993071
     118       0.01842096
     119      -0.42482646
     120      -0.21936618
     121       0.00457020
     122       0.35264794
     123      -0.00906581
     124      -0.02805930
     125       0.05400992
     126      -0.05089703
     127       0.05447161
     128      -0.41945404
     129      -0.14440282
     130       0.54756306
     131       0.38856956
     132      -0.22550489
     133       0.75878510
     134      -0.28074467
     135      -0.29427017
     136      -0.53815046
     137      -0.13632605
     138       0.13861167
     139       0.14801645
     140       0.22907860
     141      -0.01307869
     142       0.16445853
     143      -0.21925536
     144       0.07766031
     145      -0.14875702
     146       0.22516513
     147      -0.19004121
     148      -0.05221594
     149      -0.33638323
     150       0.43648569
     151      -0.21998165
     152      -0.07638265
     153      -0.17397761
     154      -0.46402636
     155       0.09382628
     156       0.02861120
     157      -0.11193587
     158      -0.02408572
     159       0.00792759
     160       0.03119895
     161       0.11343333
     162       0.12997362
     163      -0.03895576
     164      -0.27496932
     165       0.03095442
     166       0.11846092
     167       0.00488162
     168      -0.13593791
     169      -0.00059626
     170       0.23098296
     171       0.13200719
     172       1.13483732
     173      -0.29226600
     174       0.05984545
     175      -0.19882509
     176      -0.04310689
     177       0.21688598
     178      -0.12183508
     179      -0.56973060
     180       0.02991267
     181      -0.17873043
     182      -0.08165429
     183       0.01535220
     184      -0.14050892
     185       0.26860775
     186       0.23484165
     187      -0.18528032
     188      -0.02697132
     189       0.12597690
     190      -0.20420045
     191      -0.30752364
     192       1.49136065
     193      -0.62830058
     194       0.05426828
     195      -0.92519472
     196      -0.76974249
     197      -0.10507443
     198       0.04113963
     199       0.04337829
     200       0.02285863
     201      -0.01806040
     202      -0.05042377
     203      -0.10387311
     204       0.05113062
     205       0.02660428
     206       0.11303942
     207       0.09476152
     208      -0.10518435
     209       0.08958346
     210       0.23198436
     211      -0.07870868
     212       0.02381754
     213       0.00496076
     214      -0.05885209
     215      -0.11857443
     216       0.00616688
     217       0.02437880
     218      -0.00288842
     219       0.69634114
     220       0.06496320
     221      -0.15797978
     222       0.01343098
     223      -0.00988562
     224       0.20022720
     225      -0.30234367
     226       0.08270956
     227       0.05186605
     228       0.03237297
     229      -0.05155180
     230      -0.16779764
     231       0.02416680
     232       0.00248081
     233      -0.07676563
     234       0.18764808
     235      -0.26169741
     236       2.32013819
     237       0.02908023
     238       0.05603152
     239      -0.00839422
     240       0.11116500
     241       0.06197930
     242       1.00597167
     243      -0.08138396
     244       0.08334582
     245       0.03601357
     246       0.01671990
     247      -0.01415889
     248       0.28392071
     249       0.02607783
     250       0.04660368
     251       0.05100003
     252       0.00058916
     253      -0.01237878
     254      -0.19571078
     255      -0.03210665
     256      -0.05644494
     257      -0.00843039
     258      -0.13140027
     259      -0.07004591
     260       1.08163581
     261      -0.01756743
     262      -0.06450433
     263       0.01487999
     264      -0.00455538
     265       0.21567561
     266      -0.22885947
     267       0.27422252
     268      -1.91334202
     269       0.13673031
     270      -1.24471393
     271      -0.00674998
     272      -0.08430077
     273       0.88610116
     274      -0.05538185
     275      -2.10285478
     276       1.33958690
     277      -0.11338621
     278       0.36010523
     279      -0.15896701
     280       0.08202089
     281       0.24560300
     282       1.24140433
     283       0.07365212
     284      -0.31966430
     285      -0.24818100
     286      -0.09732687
     287      -0.01425408
     288      -0.03560656
     289       0.07644026
     290      -1.11135835
     291       0.09401954
     292      -0.79107854
     293      -0.06228325
     294       0.26905871
     295       0.01211914
     296      -0.10549213
     297      -0.15498924
     298      -0.33674666
     299       0.12153508
     300       0.99171605
     301       0.09747144
     302      -0.19464962
     303      -0.00148255
     304      -0.18175766
     305       0.29393499
     306      -0.34076287
     307       0.08300806
     308      -0.35231936
     309      -0.05080367
     310      -0.09696509
     311      -0.23530810
     312      -0.15394086
     313      -0.56074988
     314       0.51464739
     315      -0.09606625
     316       0.06311634
     317      -0.00878503
     318      -0.02157579
     319       0.03843966
     320       0.00308707
     321       0.08242336
     322      -0.16413360
     323      -0.00230096
     324       0.10343138
     325      -0.02640689
     326      -0.04791285
     327       0.02282651
     328       0.03329223
     329       0.22452607
     330      -2.41759122
     331      -0.00329950
     332      -0.77074463
     333       0.00652630
     334      -0.03125498
     335      -0.20699343
     336      -0.02694495
     337       0.77472870
     338       1.48910070
     339       0.04729903
     340       0.32428589
     341      -0.04338994
     342       0.21453889
     343      -0.11854896
     344       1.10162617
     345      -0.00721930
     346      -0.31590453
     347      -0.02831762
     348      -0.18617998
     349       0.01634070
     350      -0.08190974
     351       0.12524949
     352      -0.97280564
     353      -1.31412032
     354       0.62583537
     355      -0.86145235
     356       0.38843898
     357      -0.06179795
     358      -0.82698462
     359       0.01662028
     360       0.30629561
     361       0.02097429
     362      -0.26032405
     363       0.05357472
     364      -0.46564376
     365      -0.11965252
     366       1.08576544
     367       0.00927044
     368      -0.10081814
     369       0.02600095
     370      -0.12175459
     371       0.40753427
     372      -0.26438822
     373       0.05577250
     374      -0.35893026
     375      -0.05893293
     376      -0.18556047
     377       0.00549929
     378       0.01124403
     379      -0.27015237
     380       0.52161074
     381      -0.95645922
     382       0.69587400
     383       0.56213512
     384      -0.22004687
     385      -0.06307717
     386       0.01322886
     387       0.01053878
     388       0.19792133
     389       0.05045706
     390       0.21844194
     391       0.07653981
     392      -0.19073505
     393      -0.08022574
     394       0.20432944
     395      -0.06372936
     396      -0.23282587
     397       0.00637915
     398       0.09925142
     399      -0.11143676
     400       0.08941798
     401      -0.05062435
     402      -0.01413384
     403      -0.40444543
     404      -0.02987015
     405      -0.01947785
     406       0.01039916
     407       0.03043652
     408       0.18734388
     409       0.03621639
     410      -0.01136991
     411      -0.00897187
     412      -0.03821965
     413      -0.33965635
     414       0.13341816
     415      -0.15335001
     416       0.01902274
     417       0.01734642
     418       0.30580589
     419      -0.08941914
     420       0.10302947
     421       0.01601146
     422      -0.00355812
     423      -0.23210722
     424       0.12954694
     425      -0.07412900
     426      -0.02068766
     427      -0.03771252
     428      -0.21417947
     429      -0.13526799
     430       0.13015854
     431       0.00639031
     432       0.02233617
     433      -0.69736902
     434      -0.04313186
     435       0.01431535
     436       0.02020179
     437      -0.13749523
     438       0.26515753
     439       0.01735693
     440      -0.03035054
     441       0.00148525
     442      -0.00409910
     443      -0.41341995
     444       0.05807156
     445      -0.00755405
     446      -0.00324134
     447       0.02141565
     448       0.05859921
     449      -0.00398320
     450      -0.02698800
     451      -0.02979398
     452       0.06431295
     453      -0.10954839
     454      -0.00700185
     455      -0.01896259
     456       0.01938982
     457       0.01542331
     458      -0.65672003
     459       0.13700827
     460      -0.21636619
     461      -0.17709969
     462      -0.05630757
     463       0.18508633
     464      -0.07394709
     465       0.05327617
     466      -0.01640086
     467      -0.41917807
     468      -0.06246071
     469      -0.05523159
     470      -0.11803481
     471       0.04599970
     472       0.13788187
     473       0.02217708
     474      -0.04809233
     475       0.09004237
     476      -0.09366849
     477      -0.07606513
     478      -0.13471047
     479       0.03745206
     480      -0.06556017
     481      -0.18088030
     482       0.20533791
     483      -0.01566866
     484       0.17936612
     485      -0.08975286
     486       0.03668866
     487      -0.04747901
     488       0.01666828
     489       0.04715764
     490      -0.01552580
     491      -1.84641791
     492       1.13914347

 Residual norm when dim(red L) =   5
 NEO root     CSF        orbital          total
    1     0.01019588     0.00283885     0.01058371 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  1

   1.988077390   1.976686096   0.023751279   0.012882831

 Symmetry  2

   1.980793425   0.018334407

 Symmetry  3

   1.974985206   0.024489366


 <<< MACRO ITERATION  3 >>>
 --------------------------

 Total MCSCF energy       :      -76.168759959060097       (MACRO    3)
 - Dielec. solvation ener.:       -0.012646485345965

 Norm of total gradient   :        0.010417429969
 -    of CI gradient      :        0.010052631296
 -    of orbital gradient :        0.002732663750
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.03136373
       2      -0.11177987
       3      -0.17477282
       4       0.04286885
       5       0.02279976
       6       0.02769635
       7       0.04347811
       8      -0.13188182
       9      -0.14621616
      10       0.00255432
      11      -0.02686363
      12       0.01358106
      13       0.00600541
      14      -0.01454385
      15      -0.00901956
      16       0.34040092
      17       0.43404079
      18       0.02184437
      19       0.01192032
      20       0.14551416
      21      -0.18059350
      22       0.11861012
      23      -0.04414762
      24      -0.26768076
      25      -0.36349651
      26       0.17315105
      27      -0.30985175
      28      -0.06745583
      29       0.09736514
      30      -0.24731851
      31      -0.26609081
      32       0.04200839
      33       0.03218235
      34      -0.03555251
      35      -0.04597488
      36       0.01315508
      37       0.03305768
      38      -0.14413122
      39      -0.09663770
      40      -0.16292468
      41      -0.14300631
      42      -0.02400013
      43       0.16243063
      44      -0.29378218
      45      -0.40921070
      46      -0.22217412
      47      -0.20087319
      48      -0.27979386
      49      -0.09737862
      50      -0.09728427
      51      -2.21374714
      52      -0.03502101
      53      -0.02137750
      54      -0.15599122
      55      -0.05359286
      56       0.24907856
      57       0.00619011
      58       0.05138837
      59       0.01175917
      60       0.01236795
      61      -0.07786335
      62       0.32081133
      63       0.16062538
      64       0.00775822
      65      -0.02890755
      66      -0.01812318
      67      -0.03534240
      68      -0.13165246
      69      -0.09301512
      70       0.01221863
      71      -0.08775585
      72      -0.09006510
      73      -0.02171651
      74      -0.05274453
      75      -0.00177206
      76       0.05222837
      77       0.02842764
      78       0.00645165
      79      -0.04145321
      80      -0.56762009
      81       0.01638097
      82      -0.62385250
      83      -0.08094208
      84      -0.03312858
      85       0.03532862
      86       0.01569432
      87      -0.01476290
      88      -0.04049917
      89       0.07525047
      90      -0.38010960
      91      -0.08672548
      92      -0.22946459
      93      -0.01254906
      94       0.03294239
      95      -0.00443784
      96      -0.10884853
      97      -0.02273480
      98       0.10721598
      99       0.04143270
     100      -0.32775596
     101      -0.13651820
     102      -1.33004550
     103       0.05234597
     104       0.00025673
     105      -0.01921879
     106      -0.04043661
     107       0.01686073
     108      -0.00801918
     109      -0.00269436
     110      -0.08432013
     111       0.03758748
     112      -0.03041251
     113      -0.00590978
     114       0.01820999
     115      -0.00080161
     116       0.00493298
     117       0.00956263
     118      -0.02088665
     119      -0.59750049
     120      -0.63358539
     121       0.02881578
     122      -0.51496921
     123      -0.03314599
     124      -0.01274007
     125       0.05742843
     126       0.00238000
     127      -0.02323019
     128       0.11039846
     129       0.07013957
     130      -0.53409798
     131       2.02436950
     132      -0.29508202
     133      -0.04527213
     134      -0.26874387
     135      -0.30246608
     136       0.07383090
     137      -0.04397527
     138      -0.39167151
     139      -0.05588103
     140       0.15375574
     141       0.03034165
     142       0.09180714
     143       0.31482195
     144      -0.04476875
     145       0.02930385
     146      -0.08017325
     147      -0.02508544
     148       0.08874814
     149       0.17133456
     150      -0.44746406
     151       0.19213368
     152      -0.02798447
     153       0.13630295
     154       0.44627605
     155       0.17484761
     156      -0.01112975
     157       0.10976060
     158      -0.01537300
     159       0.00445039
     160      -0.04020346
     161      -0.04721297
     162      -0.39290311
     163       0.00248817
     164       0.15925200
     165      -0.02450858
     166      -0.11694254
     167      -0.03343437
     168       0.14786684
     169       0.00796625
     170      -0.31472320
     171      -0.04367642
     172      -1.24083098
     173       0.57595787
     174      -0.02362919
     175       0.06701730
     176      -0.00227488
     177      -0.30161424
     178       0.08641883
     179       0.13481669
     180      -0.18124279
     181       0.16292050
     182      -0.03716646
     183      -0.00684095
     184       0.25110133
     185      -0.16967045
     186      -0.23679812
     187       0.24176448
     188      -0.01636662
     189      -0.22729113
     190       0.15775149
     191       0.38892231
     192      -1.41219172
     193       0.70514243
     194       0.01053569
     195       0.71649144
     196       0.84741943
     197       0.05920416
     198      -0.03357357
     199      -0.00818532
     200      -0.02918715
     201       0.03003695
     202       0.01296969
     203       0.01597196
     204       0.04964835
     205      -0.02514900
     206      -0.07037905
     207       0.00014121
     208      -0.03639664
     209      -0.02962121
     210      -0.11393539
     211       0.01927405
     212      -0.01558182
     213      -0.00157014
     214       0.00669861
     215       0.01438882
     216       0.00099482
     217      -0.01872754
     218      -0.00329741
     219       0.39660924
     220       0.02984129
     221       0.15427723
     222      -0.06947436
     223      -0.02743273
     224      -0.11852216
     225       0.11306643
     226      -0.06344734
     227       0.00211716
     228      -0.03768325
     229       0.04598370
     230       0.06909605
     231      -0.00149932
     232      -0.02425416
     233       0.00670830
     234       0.03644662
     235      -0.13280616
     236       0.28606931
     237      -0.00035241
     238      -0.06203699
     239       0.01625662
     240      -0.04688820
     241       0.03508511
     242      -0.62831131
     243       0.03742324
     244      -0.01056379
     245       0.00449923
     246      -0.04150439
     247       0.01935561
     248      -0.25417221
     249      -0.00937967
     250      -0.03484929
     251      -0.02970470
     252      -0.00097273
     253       0.00765253
     254       0.21411111
     255       0.01070378
     256       0.07044166
     257      -0.01478620
     258       0.08356488
     259       0.08301343
     260      -0.85821355
     261       0.00828120
     262       0.05174435
     263      -0.01200427
     264      -0.00543751
     265      -0.12734930
     266      -0.16895031
     267       0.13873099
     268      -0.07062694
     269       0.54676087
     270      -0.22122596
     271      -0.03934644
     272       0.08325058
     273       0.02846297
     274      -0.07121628
     275      -0.49690841
     276       0.32223192
     277       0.07260256
     278      -0.11362612
     279       0.02830237
     280       0.10067289
     281       0.05234512
     282      -0.84888600
     283      -0.00678595
     284       0.30932898
     285       0.02297321
     286      -0.07185112
     287      -0.01937195
     288       0.03148852
     289      -0.15740796
     290       0.84351443
     291      -0.11180818
     292       0.62430943
     293       0.05300178
     294      -0.12995102
     295      -0.01028764
     296       0.14640641
     297      -0.00752722
     298       0.14906003
     299       0.05102010
     300      -0.82543254
     301      -0.05209180
     302       0.11345558
     303       0.01547888
     304       0.12473522
     305      -0.25182917
     306      -0.08964887
     307      -0.10198145
     308       0.38873907
     309       0.05735940
     310       0.05374440
     311       0.18216591
     312       0.11493332
     313       0.38295350
     314      -0.18015136
     315       0.07361561
     316      -0.03097832
     317       0.01694570
     318       0.03293006
     319      -0.01320446
     320      -0.01667303
     321      -0.05548168
     322       0.11804345
     323       0.01028547
     324      -0.07698772
     325       0.02111649
     326       0.07107868
     327      -0.01537001
     328      -0.04122970
     329       0.11051687
     330      -0.09149882
     331      -0.17396343
     332       0.09041150
     333       0.03113433
     334       0.05032917
     335      -0.29036300
     336       0.10781935
     337       0.61961375
     338      -0.15699170
     339      -0.03527712
     340      -0.14464263
     341       0.04178415
     342       0.02655147
     343       0.05913664
     344      -0.84441064
     345      -0.00926373
     346       0.29064194
     347      -0.00153034
     348       0.01791238
     349       0.01806731
     350       0.08555767
     351      -0.08510566
     352       0.79973650
     353       0.19245588
     354       0.30885960
     355       0.70007772
     356      -0.36139207
     357      -0.03380899
     358       0.47498007
     359       0.01761379
     360      -0.17590651
     361      -0.01474699
     362       0.31055088
     363      -0.02683219
     364       0.15712343
     365       0.07070105
     366      -0.78544646
     367      -0.01882227
     368       0.05082738
     369      -0.00206886
     370       0.08464921
     371      -0.27390920
     372      -0.17479738
     373      -0.03578223
     374       0.43961005
     375       0.03930829
     376       0.10668505
     377       0.00072331
     378      -0.02814289
     379       0.18079834
     380      -0.09560845
     381       0.64958811
     382      -0.35032781
     383      -0.18469883
     384       0.00728608
     385       0.05697738
     386       0.00270312
     387       0.00992149
     388      -0.01514458
     389      -0.00163006
     390      -0.11201429
     391      -0.02781545
     392       0.18197440
     393       0.03808878
     394      -0.10449795
     395      -0.00274683
     396       0.16628375
     397      -0.00282265
     398      -0.09714059
     399       0.07035351
     400      -0.02333867
     401       0.00483159
     402       0.01311765
     403       0.72834128
     404       0.02950605
     405      -0.04292271
     406       0.03038496
     407      -0.02046984
     408      -0.08561532
     409      -0.01940530
     410       0.02470241
     411      -0.02138645
     412       0.02787835
     413       0.19513739
     414      -0.01547550
     415      -0.01192428
     416      -0.00352457
     417       0.05089112
     418      -0.07238730
     419       0.01246475
     420       0.00980377
     421      -0.00922731
     422      -0.11029570
     423      -0.26831681
     424      -0.01599206
     425       0.01607843
     426      -0.01129319
     427       0.01689721
     428       0.23180295
     429       0.03041773
     430      -0.02537272
     431      -0.02912638
     432      -0.00955629
     433       0.35295126
     434       0.01892178
     435      -0.01746613
     436      -0.01243675
     437       0.08101053
     438      -0.48797502
     439      -0.01849034
     440       0.01826245
     441       0.00000712
     442      -0.01269831
     443       0.34682401
     444      -0.00193598
     445      -0.01434396
     446       0.00117558
     447      -0.04205090
     448      -0.11775375
     449      -0.02208331
     450       0.02793982
     451       0.02086093
     452      -0.04035201
     453       0.07974143
     454      -0.00317549
     455       0.00874220
     456      -0.00607646
     457       0.02107726
     458       0.72032417
     459       0.05762262
     460       0.07066566
     461      -0.06346929
     462       0.42550302
     463      -0.41115296
     464       0.03795204
     465      -0.05023142
     466       0.02819143
     467       0.48547208
     468      -0.10453442
     469       0.02369022
     470       0.12462150
     471      -0.06048022
     472       0.21554944
     473      -0.25554142
     474       0.05001017
     475      -0.05558020
     476       0.06216646
     477       0.24612416
     478      -0.06663845
     479      -0.03370422
     480      -0.09750195
     481      -0.05440969
     482      -0.07928471
     483       0.07285284
     484       0.02604880
     485       0.06030645
     486      -0.06682077
     487       0.05964188
     488      -0.04215735
     489      -0.05717245
     490      -0.02535359
     491       1.49973269
     492      -0.94383892
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.01671375
       2      -0.12605449
       3       0.09700288
       4       0.07022313
       5       0.05156545
       6       0.05089686
       7       0.01474456
       8       0.08839009
       9       0.14022022
      10       0.01903481
      11       0.06068442
      12      -0.00145135
      13       0.01008780
      14       0.02247702
      15       0.10236732
      16      -0.07023305
      17       0.12849998
      18      -0.15335716
      19      -0.12072920
      20      -0.11226219
      21      -0.28917208
      22      -0.34361956
      23       0.21443316
      24       0.46074213
      25       0.87110265
      26      -0.11354967
      27       0.75703471
      28       0.06195546
      29      -0.11790942
      30       0.58423233
      31       0.91180098
      32       0.04837689
      33      -0.00656860
      34       0.03351233
      35       0.05536802
      36      -0.10017031
      37      -0.14286836
      38       0.11561064
      39       0.12407680
      40       0.04821451
      41       0.86729793
      42       0.11987236
      43      -0.20358398
      44       0.72620122
      45       1.20986041
      46       0.52234483
      47       0.88289422
      48       0.47574185
      49      -0.15616626
      50       0.79071295
      51       1.80729348
      52       0.04626172
      53       0.03953384
      54       0.23685510
      55       0.15438102
      56      -0.12300401
      57      -0.08433321
      58      -0.11367836
      59      -0.06138933
      60      -0.12462107
      61       0.11972588
      62      -0.06495450
      63      -0.03525303
      64      -0.06994062
      65       0.03032242
      66      -0.04632582
      67      -0.01346574
      68       0.23666330
      69       0.19738682
      70      -0.05671277
      71       0.07796073
      72       0.11045446
      73       0.09765170
      74       0.13149434
      75       0.08432789
      76       0.00727170
      77       0.04678489
      78       0.02819573
      79       0.02533627
      80       1.03572882
      81      -0.05522160
      82       1.97000365
      83       0.05806669
      84       0.11572651
      85      -0.20644060
      86      -0.16169957
      87       0.06191666
      88       0.52051383
      89      -0.05794759
      90       1.22526483
      91       0.07582429
      92       0.89006823
      93      -0.01877239
      94       1.05354729
      95       0.07947037
      96       0.18910259
      97      -0.03520630
      98      -0.11447765
      99      -0.05078094
     100       0.76302043
     101      -0.00636519
     102       2.29215418
     103      -0.13613449
     104       0.01304065
     105       0.02734021
     106       0.07810752
     107      -0.00408338
     108       0.02646339
     109      -0.00450770
     110       0.17114653
     111      -0.03646283
     112       0.07289316
     113      -0.02250853
     114       0.08073369
     115       0.00299343
     116      -0.02147427
     117      -0.02269714
     118       0.06524258
     119      -0.12515588
     120       1.08382186
     121      -0.09939405
     122       1.94063961
     123       0.08479184
     124       0.06658193
     125      -0.32651914
     126      -0.09899645
     127       0.10042558
     128       0.33839417
     129      -0.22009985
     130       1.13296122
     131      -0.94371889
     132      -0.54868649
     133       0.17331271
     134      -0.61069022
     135       1.18505383
     136      -0.15464313
     137      -0.14815892
     138       1.53150386
     139      -0.01280713
     140      -0.09575128
     141      -0.21364142
     142      -1.12583105
     143      -1.64328687
     144       0.19780039
     145      -0.03425197
     146       0.05888405
     147       0.63843544
     148      -0.26232478
     149      -0.39591658
     150       1.19054682
     151      -0.32431005
     152       0.13414550
     153      -0.19221996
     154      -0.80775081
     155       0.00866861
     156       0.03162061
     157      -0.40269209
     158       0.03054058
     159       0.05198823
     160       0.08168998
     161       0.24337969
     162       0.99355027
     163      -0.04596033
     164       1.05346033
     165       0.08619468
     166       0.24743169
     167      -0.01143284
     168      -0.18057558
     169       0.00606507
     170       0.87605505
     171      -0.06032942
     172       2.17517999
     173      -1.77529523
     174       0.26461032
     175      -0.16729191
     176       0.00994693
     177       0.86681507
     178      -0.33135267
     179      -0.36746031
     180       1.31204203
     181      -0.37234696
     182       0.10466460
     183      -0.13114666
     184      -0.97660216
     185      -0.77447378
     186       0.18207144
     187      -0.26777599
     188      -0.06533724
     189       0.93872757
     190      -0.26539731
     191      -0.77332518
     192       1.99558298
     193      -1.08371250
     194       0.04869491
     195      -0.45703446
     196      -1.51097044
     197      -0.12376470
     198       0.03318741
     199      -0.09983287
     200       0.08114143
     201      -0.07638621
     202       0.10756506
     203      -0.11125163
     204      -0.02361780
     205       0.06917644
     206       0.04846237
     207      -0.03830562
     208       0.06638289
     209       0.08810105
     210       0.25794914
     211      -0.00060856
     212       0.04177619
     213      -0.01564575
     214       0.12314614
     215      -0.05927654
     216      -0.01922236
     217      -0.00415989
     218       0.09171182
     219       0.10462583
     220       0.06899222
     221      -0.51156715
     222      -0.01297096
     223       0.00386024
     224       0.21219963
     225      -0.23716867
     226       0.07473853
     227      -0.10513558
     228       0.11206599
     229      -0.10008608
     230       0.13419408
     231       0.02207418
     232      -0.00470014
     233       0.00511782
     234      -0.10137934
     235      -0.04099858
     236       0.12088842
     237      -0.01472694
     238       0.04130786
     239       0.00111489
     240       0.25343502
     241       0.09196556
     242       1.46650630
     243      -0.02075763
     244       0.06590800
     245       0.01148929
     246      -0.04843846
     247       0.00013629
     248       0.29632222
     249      -0.01340862
     250      -0.01253384
     251       0.01981627
     252       0.23752865
     253       0.03114090
     254      -0.41093241
     255      -0.00897708
     256      -0.09208308
     257      -0.00779156
     258      -0.21194947
     259       0.05299213
     260       1.32709221
     261       0.00475476
     262      -0.02147603
     263       0.01650478
     264       0.05926494
     265      -0.04211121
     266       2.03510005
     267       0.28350874
     268      -0.09770936
     269      -0.01954443
     270      -0.09836798
     271      -0.05661675
     272      -0.15089802
     273       0.03973852
     274      -0.30607760
     275      -0.56437681
     276       0.53625502
     277      -0.02795456
     278       0.06137351
     279      -0.20451523
     280      -0.31829390
     281       0.24502247
     282       1.75517289
     283      -0.06809486
     284      -0.64389941
     285      -0.03661865
     286       0.01977497
     287       0.05321808
     288      -0.07057034
     289      -0.04033232
     290      -1.46224931
     291      -0.02806576
     292      -1.09526857
     293      -0.18964690
     294       0.25992045
     295       0.02574827
     296      -0.11523159
     297      -0.02758621
     298      -0.38545879
     299       0.16960155
     300       1.32444871
     301       0.03119527
     302      -0.17752095
     303      -0.05261239
     304      -0.39675998
     305       0.13755931
     306       2.06200933
     307       0.13839501
     308      -0.87844904
     309      -0.05981802
     310       0.13508583
     311      -0.12892053
     312      -0.27529055
     313      -0.11317655
     314      -1.38225527
     315      -0.06898837
     316       0.01234879
     317      -0.03491429
     318      -0.05263803
     319      -0.00260665
     320       0.10247164
     321      -0.00282669
     322      -0.05615642
     323       0.01120769
     324       0.08083865
     325      -0.05878196
     326      -0.11375675
     327       0.00677469
     328       0.08284318
     329       0.08238504
     330      -0.07747774
     331      -0.17147384
     332       0.01015111
     333       0.02597525
     334      -0.05599315
     335      -0.13487100
     336      -0.32626741
     337      -0.01402457
     338       0.42255157
     339       0.12172738
     340       0.21316795
     341      -0.17243055
     342      -0.19432722
     343      -0.13722700
     344       1.82214503
     345       0.00164557
     346      -0.68437126
     347      -0.04218293
     348       0.00582615
     349      -0.05090872
     350      -0.13693835
     351       0.04443865
     352      -1.51823222
     353      -0.49386732
     354       0.06011890
     355      -1.31634446
     356       0.67397144
     357      -0.08314698
     358      -1.19389776
     359       0.03402474
     360       0.27249274
     361       0.00061307
     362      -0.29325293
     363       0.05874254
     364      -0.32139480
     365      -0.10326741
     366       1.44506630
     367       0.04026902
     368      -0.16558619
     369      -0.06115825
     370      -0.41597623
     371       0.12660142
     372       2.11705262
     373       0.03390419
     374      -1.00816690
     375      -0.01359797
     376       0.15083374
     377      -0.04508409
     378      -0.14578402
     379      -0.03387759
     380      -1.61865964
     381      -1.27740111
     382       0.82755980
     383      -1.37276777
     384       0.82802371
     385       0.02387090
     386      -0.07005499
     387      -0.01827102
     388      -0.04730664
     389       0.02926389
     390       0.27240980
     391       0.00190193
     392      -0.06080397
     393      -0.08003570
     394       0.04526699
     395       0.00593689
     396      -0.32229487
     397      -0.00283400
     398       0.14201444
     399       0.01068733
     400       0.01926793
     401      -0.11260140
     402       0.03723589
     403      -2.84305398
     404      -0.04068677
     405       0.00032019
     406       0.05679768
     407       0.03508167
     408       0.32616666
     409       0.11351686
     410      -0.07057728
     411       0.02032825
     412      -0.07269752
     413      -0.38323766
     414       0.00080238
     415       0.01222547
     416       0.05536882
     417      -0.10279691
     418       0.02153772
     419       0.02748559
     420       0.03500061
     421      -0.03715897
     422       0.22962991
     423       1.39298652
     424      -0.00983948
     425      -0.00035015
     426       0.01635030
     427      -0.07335245
     428      -0.74233835
     429       0.00509245
     430       0.03493078
     431       0.00034114
     432       0.09629268
     433      -0.71817929
     434       0.04367582
     435      -0.03214100
     436       0.02188461
     437       0.07789506
     438       2.36750722
     439       0.01879469
     440      -0.01461589
     441      -0.01012087
     442       0.02573088
     443      -0.66015504
     444      -0.04094182
     445       0.10160080
     446       0.06824108
     447       0.02121417
     448       0.15159300
     449       0.01248151
     450      -0.01359437
     451      -0.01409447
     452      -0.00669655
     453      -0.26118141
     454      -0.03989712
     455       0.02983143
     456      -0.03184327
     457      -0.10648160
     458      -1.76244840
     459       0.13377082
     460      -0.24574229
     461       0.08983724
     462      -2.00844035
     463       1.39816413
     464      -0.08423096
     465       0.01779401
     466       0.00824492
     467      -1.58087918
     468       0.90329505
     469       0.04160514
     470      -0.46517335
     471       0.29504125
     472      -1.84719333
     473       1.27742568
     474      -0.08217411
     475      -0.06771765
     476       0.00557201
     477      -1.61917430
     478       1.04285492
     479      -0.06129610
     480       0.00844179
     481      -0.10798202
     482       0.10970974
     483      -0.03958305
     484       0.10935854
     485      -0.07204249
     486       0.03134102
     487       0.01616496
     488       0.04086723
     489       0.25161050
     490      -0.10505742
     491      -2.83162796
     492       1.58790463
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.05365193
       2      -0.05051515
       3       0.20255509
       4      -0.12495081
       5      -0.08774599
       6      -0.07003880
       7      -0.10122687
       8       0.20741194
       9       0.19506070
      10      -0.00265037
      11      -0.00694427
      12      -0.03591985
      13      -0.03064338
      14      -0.00346194
      15      -0.13148530
      16      -0.20327188
      17      -0.17130290
      18      -0.01324322
      19       0.03431506
      20      -0.39882019
      21      -0.78172435
      22      -0.27859130
      23       0.17865266
      24       0.13259391
      25      -0.29366184
      26      -0.49283571
      27      -0.22216447
      28       0.06484209
      29      -0.08813281
      30      -0.27208905
      31      -1.19197226
      32      -0.02397931
      33       0.03059027
      34       0.06325125
      35       0.02508513
      36      -0.03709192
      37      -0.05899112
      38       0.27870141
      39       0.11108146
      40      -0.40928336
      41      -0.62441042
      42      -0.05811346
      43      -0.16676313
      44      -0.31227505
      45      -1.28093709
      46      -0.28089337
      47      -1.31804081
      48       0.06909156
      49       0.49964045
      50      -1.36169991
      51       2.63200937
      52       0.01158325
      53      -0.01123880
      54       0.11027193
      55      -0.08753273
      56      -0.28214824
      57      -0.00919248
      58      -0.05105379
      59      -0.04225251
      60       0.02958543
      61       0.09427873
      62      -0.28370921
      63      -0.11808720
      64       0.03008959
      65      -0.00035322
      66       0.00355053
      67       0.06497640
      68       0.08514440
      69      -0.01821238
      70       0.00367032
      71       0.09073687
      72       0.11236648
      73      -0.00186184
      74       0.07140756
      75      -0.13042878
      76      -0.16541945
      77      -0.14345397
      78      -0.07669797
      79       0.06113156
      80       0.28984592
      81      -0.03016418
      82      -0.28482669
      83       0.09548642
      84      -0.11774518
      85      -0.03676254
      86       0.00652399
      87       0.05075115
      88      -0.82785941
      89      -0.08146941
      90      -1.33980209
      91       0.11073165
      92      -0.63056355
      93       0.05671678
      94      -2.19104632
      95      -0.08865854
      96      -0.04136158
      97       0.14273317
      98      -0.07259628
      99      -0.05122677
     100      -0.58806608
     101       0.38606809
     102      -0.68121473
     103      -0.03324067
     104       0.01144876
     105       0.00570304
     106       0.02588354
     107      -0.01286896
     108      -0.02161106
     109       0.00820738
     110       0.02551757
     111      -0.04125462
     112      -0.03590577
     113       0.05069437
     114      -0.22025914
     115       0.00051998
     116       0.01884559
     117       0.02788036
     118      -0.05980001
     119       0.17482378
     120       0.32300187
     121      -0.07078030
     122      -0.39862600
     123       0.02755129
     124      -0.13115332
     125      -0.12961803
     126       0.03627030
     127       0.07927027
     128      -0.98505684
     129      -0.07141303
     130      -0.78776038
     131      -0.16878166
     132      -0.54429534
     133       0.17415571
     134      -0.46009855
     135       0.39300955
     136       0.09734445
     137       0.12453631
     138      -0.73637258
     139       0.33284288
     140      -0.18180954
     141       0.08992884
     142       0.96558034
     143       0.51016489
     144       0.22375339
     145      -0.08548708
     146       0.37209058
     147      -0.89099701
     148       0.24517309
     149       0.14244626
     150      -1.14730682
     151       0.05314921
     152      -0.10352743
     153      -0.00223516
     154       0.49783507
     155      -0.14921330
     156      -0.00012556
     157       0.00084734
     158       0.00633325
     159       0.01048569
     160      -0.00103310
     161       0.08730050
     162      -0.45311626
     163       0.03608134
     164      -2.52731122
     165      -0.03308176
     166      -0.07225893
     167       0.11949216
     168      -0.07823762
     169      -0.03502536
     170      -0.78927084
     171       0.21404999
     172      -0.90736441
     173       0.15544954
     174       0.24016793
     175      -0.12937687
     176       0.33156374
     177      -0.60949346
     178       0.26509442
     179       0.22438825
     180      -2.01110180
     181       0.24821585
     182      -0.11594914
     183       0.25601580
     184       1.29265391
     185       2.13450379
     186       0.38372742
     187      -0.31052035
     188       0.27511491
     189      -1.17485690
     190       0.21859884
     191       0.21670219
     192      -0.25442007
     193      -0.10539160
     194      -0.21828226
     195      -1.21121722
     196       0.77381495
     197      -0.02374833
     198       0.01762140
     199       0.21329531
     200      -0.05141263
     201       0.03207360
     202      -0.22341945
     203      -0.01609766
     204      -0.06115848
     205       0.04791680
     206       0.05264949
     207      -0.00496236
     208       0.02639148
     209       0.02276788
     210       0.05017378
     211      -0.04036111
     212      -0.02710263
     213       0.03262765
     214      -0.23400818
     215      -0.02872312
     216       0.02024989
     217       0.09960760
     218      -0.12674972
     219      -0.21951126
     220      -0.03114690
     221      -0.12329568
     222       0.04763928
     223       0.07965239
     224       0.08886275
     225      -0.09858318
     226       0.03822388
     227       0.21337044
     228      -0.04883861
     229       0.03593834
     230      -0.42062934
     231       0.00523044
     232       0.03529552
     233       0.01849121
     234      -0.07128883
     235       0.03273827
     236      -0.90794237
     237       0.00359281
     238       0.09634123
     239      -0.01058353
     240       0.06068958
     241      -0.02191197
     242      -0.64075694
     243      -0.05867550
     244      -0.02200137
     245      -0.00712944
     246       0.09322843
     247      -0.03514563
     248       0.25981365
     249       0.01997241
     250       0.10002264
     251       0.05877459
     252      -0.05030411
     253      -0.00925581
     254      -0.11717740
     255      -0.01176463
     256      -0.07286035
     257       0.03413876
     258      -0.08438712
     259      -0.16323867
     260      -0.01617657
     261      -0.02205567
     262      -0.09522153
     263       0.01994916
     264      -0.05854703
     265       0.35339420
     266      -3.66181679
     267       0.26682882
     268       1.09035689
     269      -0.35830601
     270       0.08819472
     271       0.06065532
     272      -0.29706436
     273       0.00234303
     274       0.13866226
     275       0.07577493
     276      -0.84975529
     277      -0.09263874
     278       0.49088358
     279      -0.10717840
     280       0.29483725
     281      -0.07579730
     282      -0.55435423
     283       0.03404863
     284       0.19398758
     285      -0.10023544
     286      -0.17004457
     287       0.01641272
     288       0.13917237
     289       0.28309430
     290       0.23686433
     291       0.14924541
     292       0.35162606
     293      -0.07018463
     294       0.19315464
     295       0.01833284
     296      -0.28137334
     297       0.03241243
     298       0.23575565
     299      -0.09352171
     300      -0.09065672
     301       0.08955491
     302       0.09227367
     303       0.03304118
     304       0.42261827
     305       0.42876561
     306      -3.76599980
     307       0.12525477
     308       0.59521706
     309      -0.06797082
     310      -0.40271732
     311      -0.32772634
     312       0.17068646
     313      -0.86281010
     314       3.51800026
     315      -0.08719170
     316       0.04754969
     317      -0.00415561
     318      -0.01295578
     319       0.03960524
     320      -0.03951204
     321       0.11516644
     322      -0.19765122
     323      -0.03551202
     324       0.06376916
     325      -0.00448936
     326      -0.03064391
     327       0.03363444
     328      -0.02593192
     329       0.04209061
     330       0.81234018
     331      -0.01701208
     332      -0.06710300
     333       0.01579718
     334      -0.24764984
     335       0.13940823
     336      -0.03502652
     337      -0.54890084
     338      -0.31223332
     339       0.07302615
     340       0.47556163
     341      -0.08069029
     342       0.40095989
     343      -0.02876463
     344      -0.43035921
     345       0.02199245
     346       0.22957712
     347      -0.05112376
     348      -0.29661419
     349      -0.01096856
     350       0.04932722
     351       0.13998818
     352       0.25995948
     353       0.24769491
     354      -0.68746926
     355       0.19394066
     356      -0.11692757
     357       0.02938955
     358       0.79749916
     359      -0.03246420
     360       0.27676753
     361       0.04643039
     362      -0.48585712
     363       0.05190089
     364       0.22748811
     365      -0.09425998
     366      -0.34345885
     367      -0.00616863
     368       0.27684988
     369       0.07622022
     370       0.55800292
     371       0.55571025
     372      -3.73482932
     373       0.06440423
     374       0.62974923
     375      -0.08330310
     376      -0.54808013
     377       0.07070170
     378       0.37901230
     379      -0.42105030
     380       3.67938197
     381       0.54204715
     382      -0.33689074
     383       3.33122070
     384      -1.69870406
     385      -0.13847060
     386       0.03809626
     387       0.00247965
     388       0.06345973
     389       0.01391773
     390       0.07449144
     391       0.05706443
     392      -0.32092082
     393      -0.06023466
     394       0.11544951
     395      -0.00403708
     396      -0.10486787
     397       0.00500524
     398       0.02485648
     399      -0.07150994
     400       0.03951778
     401      -0.09328693
     402       0.00698223
     403       1.04114201
     404      -0.05493573
     405       0.05265556
     406       0.00501254
     407       0.01757279
     408       0.36387702
     409       0.01511445
     410      -0.01008722
     411       0.02044930
     412      -0.03120551
     413      -0.28235737
     414       0.03893449
     415      -0.00400848
     416       0.05776491
     417      -0.06660008
     418       0.60013439
     419      -0.03114290
     420       0.02311573
     421      -0.00995290
     422       0.09919520
     423      -1.51339216
     424       0.03186812
     425      -0.03293534
     426       0.01156084
     427       0.00191287
     428       0.51316963
     429      -0.05199356
     430       0.03915640
     431       0.05881913
     432      -0.00294448
     433       0.26236649
     434      -0.11091502
     435       0.07347644
     436       0.00316298
     437      -0.35162503
     438      -3.30333257
     439       0.00714727
     440      -0.01943885
     441       0.00773550
     442      -0.01133635
     443       0.27016542
     444       0.03177920
     445      -0.00531368
     446       0.00044989
     447       0.05715718
     448      -0.14932540
     449       0.03638565
     450      -0.06276265
     451      -0.04192744
     452       0.10776114
     453       0.30783251
     454       0.03482453
     455      -0.05299832
     456       0.05088150
     457       0.04186298
     458       1.63495681
     459       0.01468692
     460      -0.16748927
     461       0.05906512
     462       1.33489706
     463      -0.52757026
     464      -0.03420345
     465       0.13603088
     466      -0.06722912
     467       1.95744561
     468      -1.55967293
     469       0.00854105
     470      -0.18903667
     471       0.10264713
     472       1.56190817
     473      -0.71286058
     474      -0.05259135
     475       0.25741297
     476      -0.19466404
     477       2.68027155
     478      -1.85966841
     479       0.04063143
     480       0.14674806
     481       0.11948676
     482       0.11593711
     483      -0.12648288
     484      -0.06353711
     485      -0.08095953
     486       0.13994225
     487      -0.10051880
     488       0.08355451
     489       0.03336742
     490       0.07091255
     491       0.50418537
     492      -0.17597701
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.00628505
       2       0.20638558
       3       0.12438250
       4      -0.12067817
       5      -0.05856261
       6      -0.10159959
       7      -0.02879618
       8      -0.18593088
       9      -0.32890125
      10      -0.02387777
      11      -0.05887894
      12       0.01701994
      13      -0.00781646
      14      -0.00495699
      15      -0.04789835
      16      -0.26195996
      17      -0.38297240
      18       0.18083265
      19       0.15089061
      20      -0.07144922
      21      -0.01007764
      22       0.26760377
      23      -0.05811043
      24      -0.26337839
      25      -0.28525509
      26       0.03211248
      27      -0.04204219
      28      -0.00321057
      29       0.08063656
      30      -0.19017798
      31       0.96675038
      32      -0.10404133
      33      -0.00423978
      34      -0.02623264
      35      -0.05193625
      36       0.11137211
      37       0.15657783
      38      -0.26006521
      39      -0.24065806
      40       0.10780477
      41       0.01036776
      42      -0.07190218
      43       0.18744621
      44      -0.17003089
      45       0.80623053
      46       0.24480166
      47       1.05157699
      48      -0.62726062
      49      -0.07104996
      50       1.05112720
      51      -2.83863040
      52      -0.04164456
      53      -0.01987815
      54      -0.43475901
      55      -0.13800558
      56      -0.31623680
      57       0.11976899
      58       0.15139154
      59       0.06309771
      60       0.17297830
      61      -0.26816027
      62      -0.39146847
      63      -0.16462966
      64       0.14285074
      65      -0.01319863
      66       0.09996331
      67       0.03502498
      68      -0.45226401
      69      -0.26576769
      70       0.08285084
      71       0.01819020
      72      -0.00838315
      73      -0.09553914
      74      -0.10057461
      75      -0.05131069
      76      -0.02708796
      77       0.01482748
      78       0.00147008
      79      -0.03693533
      80       0.27406845
      81       0.02965517
      82      -0.96264161
      83       0.03227575
      84      -0.14850721
      85       0.18088647
      86       0.24572654
      87      -0.05253362
      88       0.00478292
      89       0.01985863
      90       0.57234387
      91      -0.03441900
      92      -0.27927420
      93       0.02290001
      94       0.65938398
      95      -0.06933418
      96      -0.04030878
      97      -0.02871472
      98       0.10937107
      99       0.07303886
     100       0.36285595
     101      -0.21765468
     102      -0.12775014
     103       0.13487853
     104      -0.00994945
     105      -0.01597033
     106      -0.06448536
     107       0.00343585
     108      -0.02653768
     109       0.01569770
     110      -0.17747623
     111       0.01571550
     112      -0.04134907
     113       0.00099564
     114       0.05728311
     115      -0.02019367
     116       0.01633989
     117      -0.01979768
     118      -0.00333153
     119       0.47006238
     120       0.35271428
     121       0.04717513
     122      -1.19616432
     123      -0.07923822
     124      -0.06967042
     125       0.29831993
     126       0.14129293
     127      -0.04046935
     128       0.12491395
     129       0.26675067
     130       0.35500759
     131      -0.03051568
     132      -0.81995194
     133       0.21815740
     134      -0.70386089
     135       0.48846214
     136      -0.06405357
     137      -0.15179208
     138      -0.93861211
     139      -0.12410414
     140       0.49930416
     141       0.14976343
     142       0.68857424
     143       1.07141113
     144      -0.21168231
     145       0.15640809
     146      -0.04864972
     147      -0.22911572
     148       0.07114300
     149       0.23773881
     150       0.41170407
     151      -0.00099551
     152      -0.09748191
     153       0.03085597
     154      -0.39520364
     155      -0.19012991
     156      -0.02619240
     157       0.47203566
     158      -0.02313689
     159      -0.06773186
     160      -0.07922077
     161      -0.26636669
     162      -0.14285836
     163       0.08533211
     164       0.80028427
     165      -0.11632853
     166      -0.11924893
     167       0.03799036
     168       0.18512685
     169      -0.02076975
     170       0.42998307
     171       0.05091711
     172       0.63657006
     173       0.61516401
     174      -0.32543694
     175       0.38803143
     176      -0.04637747
     177      -0.20623735
     178       0.11338099
     179       0.09576690
     180       0.71427257
     181      -0.13412877
     182       0.02615125
     183      -0.12795069
     184      -0.76420986
     185      -1.00022078
     186      -0.54127858
     187       0.66202752
     188      -0.47294317
     189       0.61842150
     190      -0.47382657
     191       0.35012156
     192       0.70862307
     193       0.94433029
     194       0.39868421
     195       1.23577405
     196      -0.91239163
     197       0.12383526
     198       0.00526429
     199      -0.01789724
     200      -0.03892563
     201       0.03718418
     202       0.04804208
     203       0.13472190
     204      -0.04115994
     205      -0.05408622
     206       0.04550636
     207       0.05386373
     208      -0.06523902
     209      -0.09683154
     210      -0.27094070
     211      -0.02876464
     212      -0.02938622
     213       0.01743301
     214      -0.06280015
     215       0.06978457
     216       0.02772411
     217      -0.04460560
     218      -0.01406536
     219      -0.65942700
     220      -0.09426829
     221       0.54441128
     222       0.11230063
     223       0.04801504
     224      -0.24000222
     225       0.25366779
     226      -0.05672677
     227       0.00945392
     228      -0.07939108
     229       0.05683087
     230       0.08810984
     231      -0.03600487
     232      -0.00068733
     233       0.02088046
     234       0.04210024
     235       0.13018683
     236      -0.20154539
     237       0.05013122
     238      -0.03153902
     239       0.00837751
     240      -0.30908152
     241      -0.15486968
     242       0.31536580
     243      -0.00385572
     244      -0.05336537
     245      -0.01126867
     246       0.08635306
     247      -0.01277070
     248      -0.18372738
     249       0.03958446
     250      -0.02008895
     251       0.01294275
     252      -0.34576929
     253      -0.06285009
     254       0.48192012
     255       0.00433881
     256       0.10827056
     257       0.02277811
     258       0.26908831
     259      -0.12079610
     260      -0.26926137
     261      -0.01229653
     262       0.00483790
     263      -0.03405520
     264      -0.10919863
     265      -0.09611198
     266       2.60362404
     267      -0.07763710
     268       0.26720131
     269      -0.56514685
     270      -0.08626049
     271       0.15433919
     272       0.18193877
     273      -0.28875922
     274      -0.04563239
     275       0.98478692
     276      -0.86996936
     277       0.02292487
     278      -0.22081980
     279       0.12440224
     280       0.10531534
     281      -0.35759410
     282      -0.08970425
     283       0.11673743
     284       0.27851066
     285      -0.04641535
     286       0.15229698
     287      -0.02585719
     288      -0.14949737
     289       0.11887362
     290       0.16364310
     291       0.09806861
     292      -0.55767792
     293       0.25859025
     294      -0.55970478
     295      -0.05974513
     296       0.42736934
     297       0.07222804
     298       0.11757019
     299      -0.26519226
     300      -0.08304979
     301      -0.06613705
     302      -0.34711364
     303       0.10053470
     304      -0.27876279
     305      -0.41514344
     306       2.22682585
     307      -0.33325473
     308      -0.18408153
     309       0.12540059
     310       0.32222573
     311       0.37644905
     312      -0.02694417
     313       0.71281086
     314      -3.09069358
     315       0.05679477
     316       0.00575526
     317       0.04416558
     318       0.05015439
     319       0.02823001
     320      -0.12320741
     321       0.02594785
     322       0.01139540
     323      -0.01847793
     324      -0.06203893
     325       0.06988003
     326       0.12329141
     327      -0.02409814
     328      -0.03873921
     329      -0.06647074
     330       0.09301211
     331       0.27493046
     332      -0.52783112
     333      -0.01151824
     334       0.07225188
     335       0.40587044
     336      -0.13251476
     337      -0.51511256
     338      -0.14594155
     339      -0.11425005
     340      -0.40291573
     341       0.13957894
     342       0.08521629
     343       0.20046331
     344      -0.50455363
     345       0.00830518
     346       0.40490828
     347       0.01002912
     348       0.05952166
     349       0.08623015
     350      -0.04890398
     351      -0.01568045
     352       0.46757174
     353       0.14813982
     354      -0.55414781
     355       0.28956067
     356      -0.07935487
     357       0.14403378
     358      -0.58442346
     359      -0.01921074
     360      -0.57782032
     361      -0.05498245
     362       0.70579906
     363      -0.03669171
     364      -0.08173596
     365       0.12489503
     366      -0.03029203
     367      -0.08209026
     368      -0.53747102
     369       0.11463655
     370      -0.33599606
     371      -0.38901833
     372       2.39959998
     373      -0.09904616
     374      -0.02059673
     375       0.08060674
     376       0.36677472
     377       0.01508408
     378      -0.34671931
     379       0.27326261
     380      -2.96780869
     381      -0.19877828
     382      -0.03489159
     383      -2.45541048
     384       1.23256040
     385      -0.09517079
     386       0.12047181
     387       0.02761119
     388       0.09333482
     389      -0.02621800
     390      -0.30862996
     391       0.01955471
     392      -0.02594491
     393       0.08438144
     394       0.04612475
     395      -0.01145841
     396       0.37938453
     397      -0.00006325
     398      -0.12309826
     399      -0.17131201
     400       0.06032291
     401       0.00621276
     402      -0.09139610
     403       1.47406331
     404       0.04861761
     405       0.03005904
     406      -0.05567188
     407      -0.02014311
     408      -0.51273743
     409      -0.14987582
     410       0.10703250
     411      -0.04030108
     412       0.05027404
     413       0.75825183
     414       0.01623646
     415      -0.01425245
     416       0.00788936
     417       0.09655128
     418      -0.23716138
     419      -0.02640958
     420      -0.08015079
     421       0.03829258
     422      -0.19713608
     423      -0.27073483
     424       0.02566012
     425      -0.02078027
     426      -0.00814664
     427       0.10890741
     428       0.33421921
     429      -0.04675377
     430      -0.02677999
     431       0.04985637
     432      -0.15637614
     433       0.19222962
     434      -0.08475725
     435       0.05574724
     436      -0.04571379
     437      -0.04157106
     438       1.37645866
     439      -0.03302177
     440       0.02446262
     441       0.01640019
     442      -0.05311279
     443      -0.20482739
     444       0.07483551
     445      -0.14326673
     446      -0.13715016
     447       0.01632906
     448       0.23944889
     449      -0.04200548
     450       0.03230350
     451       0.03521219
     452      -0.06143347
     453      -0.22428276
     454       0.02897052
     455      -0.02483478
     456       0.05441551
     457       0.17810990
     458      -1.20155526
     459      -0.20151771
     460       0.10307512
     461      -0.05092954
     462       1.01250191
     463      -0.68131208
     464       0.12278746
     465      -0.04798507
     466      -0.01962265
     467      -1.26883323
     468       0.77032005
     469      -0.03223911
     470       0.43164123
     471      -0.28367089
     472       1.08115953
     473      -0.71290988
     474       0.12535324
     475       0.05014169
     476       0.03501997
     477      -1.64083449
     478       0.78487178
     479       0.15012129
     480      -0.01643799
     481       0.11370843
     482      -0.11785834
     483       0.03801486
     484      -0.11402976
     485       0.08481584
     486      -0.02549155
     487      -0.08116471
     488      -0.05134052
     489      -0.43871369
     490       0.21596131
     491       0.42942912
     492      -0.06629314
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.10300400
       2       0.24222212
       3      -0.13769012
       4       0.21323040
       5       0.20582716
       6       0.19595379
       7       0.11383548
       8       0.12385723
       9       0.54547104
      10      -0.00527286
      11       0.10696319
      12       0.05952543
      13       0.06363593
      14       0.01244986
      15       0.39380315
      16      -0.19262414
      17      -0.30818814
      18      -0.01665795
      19      -0.12034424
      20       0.53354591
      21       0.98454084
      22       0.44032169
      23      -0.24500342
      24       0.32073479
      25       0.16975629
      26       0.53952697
      27      -0.18691011
      28      -0.04667475
      29      -0.09951844
      30       0.11209739
      31      -1.21567038
      32      -0.04988888
      33      -0.10098082
      34      -0.06069911
      35       0.07472995
      36       0.13693055
      37       0.06800251
      38       0.01814314
      39       0.35749786
      40       0.68288417
      41       0.03138189
      42       0.16445929
      43      -0.09344331
      44       0.01113583
      45      -0.86081504
      46      -0.55937022
      47      -1.11790697
      48       1.31803305
      49      -0.72442915
      50      -1.11313326
      51       1.93350689
      52       0.07578043
      53       0.07264127
      54       0.83548702
      55       0.55034863
      56      -0.09160336
      57       0.06961690
      58       0.01062606
      59       0.13060643
      60      -0.12018251
      61       0.36858598
      62      -0.29720290
      63      -0.00478488
      64      -0.06398016
      65       0.08891888
      66       0.07003305
      67      -0.06355506
      68       0.79058816
      69       0.71233355
      70      -0.04968013
      71      -0.17386016
      72      -0.10423003
      73       0.01609325
      74      -0.13225398
      75       0.25452511
      76       0.29078562
      77       0.19123219
      78       0.16786798
      79      -0.19774684
      80       0.49920346
      81       0.14268813
      82       0.15524125
      83      -0.19758298
      84       0.21054212
      85       0.08827137
      86      -0.27291172
      87      -0.11470210
      88       0.71229893
      89       0.06327652
      90       0.38245010
      91      -0.09049461
      92       0.46435607
      93      -0.13035944
      94       1.68338693
      95       0.24612511
      96       0.22113045
      97      -0.16193407
      98      -0.18292483
      99      -0.04907019
     100      -0.29263057
     101      -0.18123428
     102      -0.35286113
     103      -0.02807941
     104      -0.02137649
     105       0.03091243
     106      -0.00362278
     107       0.00446616
     108       0.07757488
     109      -0.03371586
     110       0.11451873
     111       0.03896035
     112       0.14774675
     113      -0.10609201
     114       0.34766032
     115       0.03594085
     116      -0.08114170
     117      -0.00439248
     118       0.11268281
     119       0.54771516
     120       0.68189289
     121       0.14666667
     122       0.48834937
     123      -0.03465293
     124       0.23879386
     125       0.23547007
     126      -0.31066065
     127      -0.14199599
     128       0.90831924
     129      -0.00455855
     130       0.14324244
     131      -0.57569281
     132      -0.53173367
     133      -0.46521811
     134      -0.45982052
     135       0.80267956
     136      -0.44902203
     137      -0.26878921
     138       1.01689571
     139      -0.28281805
     140       0.39307246
     141      -0.03533195
     142      -1.28192983
     143      -0.55038230
     144      -0.01247180
     145      -0.32028006
     146      -0.31871791
     147       0.94466141
     148      -0.36389326
     149      -0.31637131
     150       0.73013010
     151       0.06747514
     152      -0.01755441
     153      -0.11184184
     154       0.30978671
     155       0.05352683
     156       0.05555388
     157      -0.31921544
     158      -0.01898597
     159       0.00881202
     160       0.09057200
     161      -0.10695282
     162       0.06327122
     163      -0.14452671
     164       1.60554846
     165       0.22679994
     166       0.37442416
     167      -0.21994138
     168      -0.35434200
     169       0.06269735
     170      -0.33127692
     171      -0.55860489
     172      -2.25891105
     173       0.32417147
     174       0.10051722
     175      -0.47515904
     176      -0.10611348
     177       0.35619996
     178      -0.26379028
     179      -0.17420526
     180       0.99447687
     181       0.04973701
     182      -0.10001206
     183      -0.10552695
     184       0.18964008
     185      -0.78301287
     186       0.50262440
     187      -1.18177484
     188       0.56883476
     189       0.00243297
     190       1.12025270
     191      -0.79025670
     192      -3.33501986
     193      -1.70572056
     194      -0.73075693
     195      -0.81969093
     196       3.23597849
     197      -0.07910804
     198      -0.05631762
     199      -0.38926488
     200       0.14797632
     201      -0.11703419
     202       0.33442903
     203       0.03957294
     204       0.05704621
     205      -0.07156415
     206      -0.05857415
     207       0.05237582
     208      -0.02951112
     209       0.02556082
     210       0.08232607
     211       0.06189015
     212       0.13378541
     213      -0.11049994
     214       0.64020827
     215       0.07598102
     216      -0.09195706
     217      -0.17819529
     218       0.21866026
     219      -0.50575934
     220       0.08356948
     221      -0.01329915
     222      -0.08978148
     223      -0.07758785
     224       0.04184722
     225       0.03039028
     226       0.03615556
     227      -0.53944355
     228       0.16440603
     229      -0.21436946
     230       0.64760025
     231      -0.01415208
     232       0.00689564
     233      -0.05784007
     234       0.10962850
     235       0.17403669
     236       0.80290192
     237      -0.03197802
     238      -0.11324382
     239      -0.00091183
     240      -0.02212499
     241       0.04848158
     242      -0.56777565
     243       0.10251205
     244       0.08924746
     245       0.03161861
     246      -0.13776396
     247       0.07553518
     248      -0.24757042
     249      -0.05139096
     250      -0.01952605
     251      -0.13966021
     252       0.24805449
     253       0.01386773
     254      -0.56341713
     255      -0.00553467
     256      -0.00001304
     257      -0.06997786
     258      -0.02129805
     259       0.27799104
     260      -0.61949169
     261       0.06319640
     262       0.15303424
     263       0.01818536
     264       0.25410150
     265      -0.27431560
     266      -2.34768143
     267      -0.32984924
     268      -1.05293723
     269      -0.36510019
     270       0.46918067
     271       0.03141583
     272      -0.06889125
     273      -0.32332573
     274      -0.22116616
     275       0.79780704
     276       1.01312230
     277       0.10881924
     278      -0.36478990
     279       0.16249725
     280      -0.16726116
     281       0.19711799
     282      -0.12988452
     283      -0.11472570
     284      -0.44553278
     285       0.17383628
     286       0.04413774
     287      -0.02632356
     288       0.22760475
     289      -0.43153749
     290       0.55323020
     291      -0.23758485
     292       1.19065336
     293       0.05774601
     294       0.71973755
     295       0.04219523
     296      -0.34206068
     297      -0.08739939
     298      -0.15129722
     299       0.17941293
     300      -0.74682861
     301      -0.04592320
     302       0.87636037
     303      -0.22860294
     304      -0.00706338
     305       0.38535229
     306      -1.07356198
     307       0.54377329
     308      -0.16702393
     309      -0.14747013
     310      -0.09548676
     311      -0.35296192
     312      -0.01164623
     313      -0.39666239
     314       3.76157687
     315       0.06228259
     316      -0.06364542
     317      -0.04215562
     318      -0.05780835
     319      -0.12985726
     320       0.14964114
     321      -0.23055400
     322       0.30238602
     323       0.08977374
     324       0.02430785
     325      -0.05753741
     326      -0.10789934
     327      -0.01416622
     328       0.12980240
     329      -0.20325524
     330      -0.82420887
     331       0.35760681
     332       0.08340656
     333      -0.06191885
     334       0.02214102
     335       0.15511789
     336      -0.09910979
     337       0.24864806
     338       0.90199743
     339      -0.07110341
     340      -0.29070005
     341       0.07100365
     342      -0.45244601
     343      -0.12636339
     344       0.09231911
     345      -0.00207876
     346      -0.67327910
     347       0.09735387
     348       0.28395595
     349      -0.00354265
     350       0.27142768
     351      -0.19731484
     352       0.20651526
     353      -0.54455874
     354       0.75646740
     355       0.39686530
     356      -0.18087370
     357      -0.05849106
     358       0.83912256
     359      -0.00010894
     360       0.63806748
     361       0.03096330
     362      -0.51908767
     363      -0.10596805
     364       0.07158626
     365       0.00275422
     366      -0.55142056
     367       0.18999855
     368       0.96236447
     369      -0.32039981
     370      -0.05993642
     371      -0.14444055
     372      -1.55594661
     373       0.09686932
     374      -0.54614807
     375      -0.03855526
     376       0.00267041
     377      -0.19377502
     378       0.42502398
     379      -0.00651998
     380       3.18214275
     381       0.41363643
     382      -0.01817347
     383       2.12496054
     384      -1.21166366
     385       0.27894800
     386      -0.10612984
     387      -0.03736935
     388      -0.14082939
     389      -0.07079349
     390       0.05253047
     391      -0.12293614
     392       0.45510189
     393       0.11041832
     394      -0.10152223
     395       0.04104770
     396      -0.15819948
     397      -0.00548921
     398       0.23454898
     399       0.15612115
     400      -0.15161189
     401       0.18261124
     402      -0.04036325
     403      -0.73232888
     404       0.08508362
     405      -0.04510879
     406      -0.00931346
     407      -0.00578738
     408       0.18697673
     409       0.00097697
     410      -0.06198689
     411       0.00540249
     412       0.10279104
     413      -0.85421283
     414      -0.01165406
     415       0.00222315
     416      -0.08919264
     417       0.08722728
     418      -0.42666622
     419       0.00427389
     420      -0.06500226
     421      -0.00173186
     422      -0.10912613
     423       1.39564444
     424      -0.09540852
     425       0.08941696
     426       0.00670380
     427      -0.03947493
     428      -0.85460171
     429       0.13762208
     430      -0.07673643
     431      -0.11186379
     432       0.04821831
     433      -0.25500941
     434       0.28422710
     435      -0.15257687
     436       0.05599338
     437       0.64364863
     438       1.46078253
     439       0.02746261
     440      -0.01017012
     441      -0.03330517
     442       0.15755058
     443       0.38216470
     444      -0.14044282
     445       0.09765324
     446       0.03871248
     447      -0.15216472
     448      -0.37074771
     449       0.00209421
     450       0.08975058
     451       0.02277202
     452      -0.02503665
     453      -0.05171711
     454      -0.00244448
     455       0.08774759
     456      -0.16163743
     457      -0.37696375
     458       0.97515313
     459       0.03754947
     460       0.35875008
     461       0.09159526
     462      -1.42598314
     463       0.53652955
     464      -0.08172018
     465      -0.08454469
     466       0.14508571
     467       0.17691723
     468       0.47652338
     469      -0.04363261
     470       0.27369095
     471      -0.15687586
     472      -1.97690199
     473       0.89191600
     474      -0.11703481
     475      -0.45898778
     476       0.23515681
     477       0.06235650
     478       0.68589562
     479      -0.13772419
     480      -0.09581497
     481      -0.16278993
     482      -0.13973237
     483       0.16731052
     484       0.07215227
     485       0.05117233
     486      -0.20256675
     487       0.19629550
     488      -0.12229513
     489       0.30075837
     490      -0.28697094
     491       0.88703919
     492      -0.49932527
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.01669443
       2      -0.02228263
       3       0.08108766
       4       0.13958081
       5      -0.10469870
       6       0.03942251
       7       0.05274227
       8      -0.07683339
       9      -0.52537737
      10      -0.02033110
      11      -0.02564461
      12      -0.12762590
      13      -0.04132801
      14      -0.04837063
      15      -0.62197030
      16      -0.34093053
      17      -0.39802859
      18      -0.08392873
      19      -0.14395092
      20       0.36634606
      21       0.90887862
      22      -0.28258309
      23      -0.00216708
      24       0.99417898
      25       0.61929299
      26       0.35096600
      27       0.21862069
      28      -0.16160274
      29       0.18926309
      30       0.90740119
      31       1.05545472
      32       0.05716822
      33      -0.28869810
      34       0.02760732
      35      -0.07009586
      36      -0.07026432
      37      -0.05497123
      38       0.10553620
      39      -0.26775073
      40       0.25790481
      41       0.16616080
      42      -0.15403701
      43       0.12421227
      44       0.68469007
      45       0.66404192
      46       0.23357336
      47       0.14326932
      48      -1.06277376
      49       1.01810394
      50       0.83336156
      51       0.97615569
      52      -0.01419757
      53      -0.06343189
      54      -0.60223689
      55      -0.79152294
      56       0.13786598
      57      -0.09062313
      58      -0.20958343
      59      -0.03788425
      60      -0.17195335
      61      -0.07134076
      62      -0.40546288
      63      -0.32431805
      64      -0.23631835
      65      -0.10574960
      66      -0.14124233
      67      -0.05371911
      68      -0.34400832
      69      -0.83523860
      70      -0.11717082
      71      -0.30927253
      72      -0.21688205
      73       0.13837461
      74       0.10955921
      75       0.10073275
      76       0.09946268
      77      -0.16486392
      78      -0.17348508
      79      -0.14410457
      80      -0.37908879
      81       0.16128857
      82      -0.35328695
      83      -0.27811147
      84       0.19186054
      85      -0.09596302
      86      -0.12266818
      87       0.04520840
      88      -0.35872946
      89       0.00083283
      90      -1.16435423
      91      -0.02769487
      92      -0.54252173
      93      -0.00329251
      94      -1.41176086
      95      -0.07059513
      96      -0.51086317
      97      -0.02140527
      98       0.32897302
      99       0.00806825
     100      -0.28805244
     101       0.55039947
     102       1.00684579
     103      -0.10561863
     104       0.01739661
     105      -0.04618139
     106      -0.01370046
     107       0.00363442
     108      -0.04944162
     109      -0.02752385
     110       0.07487084
     111       0.05567075
     112      -0.14886398
     113       0.14961855
     114      -0.52474825
     115       0.02302681
     116       0.10197118
     117       0.02392172
     118      -0.18624247
     119       0.85228435
     120      -0.47792397
     121       0.18011798
     122      -0.15083719
     123      -0.02446510
     124       0.03617134
     125      -0.16568057
     126      -0.06932542
     127      -0.11657467
     128      -0.56085586
     129      -0.22778560
     130      -1.07571800
     131       0.43157782
     132      -0.06220035
     133      -0.56853739
     134       0.36961119
     135      -0.25834357
     136       0.02502475
     137      -0.14647462
     138      -0.02717334
     139       0.13944472
     140      -0.27288883
     141       0.10344307
     142      -0.19999154
     143       0.07043280
     144       0.34556799
     145       0.16113004
     146      -0.02067349
     147      -0.54392074
     148       0.02854805
     149      -0.08284143
     150      -1.49754812
     151       0.31249225
     152       0.29410568
     153       0.36274592
     154       0.51217777
     155      -0.06612024
     156       0.03964432
     157      -0.38762855
     158      -0.00021304
     159       0.10865105
     160      -0.01656649
     161       0.08773891
     162      -0.72675840
     163      -0.07764738
     164      -1.02556846
     165      -0.00419038
     166      -0.47874490
     167      -0.17472708
     168       0.43343823
     169       0.06609316
     170      -0.01048596
     171       0.40200718
     172       3.17392236
     173       0.72275366
     174       0.32170767
     175       0.01405461
     176      -0.07184660
     177      -0.42156777
     178      -0.03199950
     179       0.14064191
     180      -1.24559380
     181       0.40164922
     182       0.17169341
     183       0.34338289
     184       0.42321923
     185       0.85809265
     186       0.20262402
     187       0.76179554
     188       0.27637507
     189      -0.68500624
     190      -0.89171472
     191       0.59895740
     192       3.39352042
     193       1.25574727
     194       0.33857753
     195      -0.14515909
     196      -5.22600275
     197      -0.00755460
     198      -0.01592147
     199       0.46045472
     200      -0.11034021
     201       0.15095580
     202      -0.70548748
     203      -0.08536014
     204       0.15459023
     205      -0.04420809
     206      -0.24658706
     207      -0.01869351
     208       0.02883129
     209       0.05974954
     210       0.03384280
     211       0.11763665
     212      -0.10565998
     213       0.05978216
     214      -0.52095006
     215      -0.04401388
     216       0.02430199
     217       0.22990794
     218      -0.23429635
     219      -0.38260840
     220       0.13928524
     221      -0.30629015
     222      -0.23562036
     223      -0.04625950
     224       0.16479142
     225      -0.12349037
     226      -0.01630226
     227       0.59071815
     228      -0.05611708
     229       0.26393887
     230      -1.08552751
     231       0.03690643
     232       0.01926219
     233      -0.05498133
     234       0.07222823
     235       0.05000855
     236       0.81012057
     237      -0.11665344
     238       0.11250780
     239      -0.08563967
     240       0.24423923
     241       0.16183053
     242       0.02358094
     243       0.00932707
     244      -0.05874030
     245      -0.04603075
     246      -0.06847332
     247      -0.00891145
     248       0.05353817
     249      -0.04624613
     250      -0.08332259
     251      -0.03109760
     252       0.19230210
     253       0.08342781
     254       0.11566520
     255       0.02738128
     256      -0.16082306
     257       0.01248474
     258      -0.26582972
     259       0.06303829
     260       1.79253851
     261      -0.02251097
     262      -0.02495831
     263       0.04465372
     264       0.07417749
     265       0.30916010
     266       1.58736369
     267      -0.12441059
     268      -1.10042318
     269       0.13416719
     270       0.80718483
     271      -0.08380337
     272      -0.28117686
     273      -0.24537583
     274      -0.31233689
     275       1.59266610
     276       1.74561774
     277       0.03692879
     278       0.04508371
     279       0.04027707
     280       0.09012848
     281       0.20907089
     282       0.37013773
     283      -0.10492331
     284       0.12516571
     285       0.22118900
     286      -0.43755102
     287      -0.00099015
     288      -0.04302074
     289      -0.00688121
     290      -1.39485440
     291      -0.10217386
     292      -0.25893788
     293      -0.29736373
     294      -0.00318520
     295      -0.00576763
     296      -0.22606554
     297      -0.04824806
     298      -0.16244460
     299       0.14123192
     300       1.08605992
     301       0.17552331
     302      -0.21858012
     303      -0.05093761
     304       0.78758496
     305      -0.08500035
     306       0.08743582
     307      -0.23763530
     308       0.84933210
     309      -0.00406654
     310      -0.59549463
     311       0.11310171
     312       0.06132197
     313      -0.12358103
     314      -3.38885709
     315      -0.05795206
     316      -0.00961764
     317      -0.00429297
     318      -0.05871124
     319      -0.02270672
     320       0.08681637
     321       0.05012057
     322      -0.12234504
     323      -0.02896406
     324      -0.00343386
     325      -0.02331630
     326      -0.06741176
     327       0.03227194
     328      -0.19537709
     329      -0.33078317
     330      -0.47349742
     331       0.13990244
     332       0.38538844
     333      -0.09812462
     334      -0.23013868
     335      -0.06803684
     336       0.02520141
     337      -0.09627527
     338       0.91756834
     339       0.07074363
     340       0.27565057
     341      -0.03414084
     342      -0.03543309
     343      -0.10398798
     344       0.90748256
     345      -0.04924205
     346       0.11185734
     347       0.06713578
     348      -0.41086824
     349      -0.15853971
     350      -0.14592988
     351       0.04415258
     352      -1.78719352
     353      -0.52754687
     354       0.57565205
     355      -1.10291311
     356       0.46955579
     357      -0.14247602
     358      -0.21647309
     359      -0.02434425
     360      -0.00094206
     361       0.06139676
     362      -0.30381572
     363      -0.00442254
     364       0.03338385
     365      -0.00883256
     366       0.86045699
     367      -0.05299814
     368      -0.03832296
     369       0.03010974
     370       0.73151176
     371       0.50172539
     372       0.17795195
     373      -0.02591494
     374       0.97920462
     375      -0.06655620
     376      -0.66873649
     377       0.24853144
     378      -0.13308153
     379      -0.15680873
     380      -2.53624138
     381      -0.43903098
     382       0.20774373
     383      -1.22952036
     384       0.89770638
     385       0.10740080
     386      -0.15181968
     387      -0.04027783
     388      -0.16642034
     389      -0.02437641
     390       0.25324013
     391      -0.03563573
     392       0.00992785
     393      -0.06944835
     394      -0.14060970
     395       0.04307585
     396      -0.29544781
     397       0.00896365
     398      -0.17038811
     399       0.46766715
     400      -0.32351451
     401       0.32739513
     402       0.22299978
     403       0.07378223
     404       0.03742036
     405      -0.04844913
     406       0.01318376
     407       0.07006399
     408       0.65259082
     409       0.11902164
     410      -0.13487747
     411       0.09225972
     412      -0.07187293
     413       0.02299488
     414      -0.03438732
     415       0.04476367
     416      -0.09306574
     417       0.00320069
     418       0.15537277
     419       0.01172825
     420       0.11426277
     421      -0.08009713
     422       0.10027957
     423      -0.82183543
     424      -0.00038433
     425       0.01286908
     426      -0.01890281
     427      -0.16526999
     428       0.03345070
     429       0.04314221
     430      -0.00991591
     431      -0.11113608
     432       0.20759799
     433       0.08058757
     434       0.12140304
     435      -0.08148100
     436       0.11658513
     437      -0.14813529
     438      -3.01615643
     439       0.12780262
     440      -0.06002682
     441      -0.00713264
     442      -0.07485188
     443       0.54980662
     444      -0.08507079
     445       0.11599472
     446       0.17604419
     447       0.06802864
     448      -0.07566834
     449       0.05176902
     450      -0.08429324
     451      -0.08218692
     452       0.08290568
     453       0.74272729
     454      -0.19809401
     455       0.07902595
     456      -0.01876222
     457       0.06265379
     458       0.73083425
     459       0.29817282
     460       0.21560081
     461       0.18669222
     462      -0.11912123
     463      -0.18297209
     464      -0.19380042
     465      -0.13044747
     466      -0.13738981
     467       1.67242512
     468      -1.14836021
     469       0.14566896
     470      -0.19571362
     471       0.24378550
     472      -0.06934997
     473      -0.07239287
     474      -0.05051682
     475      -0.00966449
     476      -0.05352446
     477       1.65417382
     478      -1.00840631
     479      -0.18817109
     480       0.09099751
     481       0.03533543
     482       0.04236011
     483      -0.05373193
     484      -0.01924448
     485      -0.08202599
     486       0.05015321
     487       0.13609339
     488       0.10616192
     489       0.34981890
     490      -0.10704029
     491      -1.92700418
     492       0.53669766
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.07901914
       2      -0.05969962
       3       0.21267787
       4      -0.21827376
       5      -0.13939843
       6      -0.08894798
       7      -0.14771797
       8       0.14472847
       9       0.41830617
      10       0.00624366
      11       0.00846146
      12       0.04464991
      13       0.00214059
      14       0.10350633
      15       0.54559353
      16      -0.43812803
      17      -0.42044684
      18      -0.00267284
      19       0.15626306
      20      -0.27853024
      21       0.38412667
      22      -0.48710657
      23       0.14472319
      24      -0.52736834
      25      -0.94236712
      26      -0.37122466
      27       0.00876362
      28       0.27234126
      29      -0.08866661
      30      -1.49868223
      31      -0.91883785
      32       0.08255753
      33      -0.05412798
      34       0.13940423
      35       0.06915403
      36      -0.06532778
      37      -0.07242719
      38       0.07448010
      39       0.16323861
      40      -0.22297907
      41       0.00344639
      42       0.05805920
      43      -0.10433649
      44      -1.01441736
      45      -0.26825090
      46       0.20870787
      47       1.00209929
      48       0.47154258
      49      -0.38797765
      50      -0.86634435
      51      -1.92650358
      52       0.02672119
      53       0.06994614
      54       0.20585364
      55       0.71566322
      56       0.37575607
      57      -0.03234069
      58      -0.02379033
      59      -0.09582376
      60       0.14360623
      61      -0.10936118
      62      -0.41197087
      63      -0.17243466
      64       0.09891150
      65       0.03874926
      66       0.04596548
      67       0.17342049
      68      -0.06313469
      69       0.65627043
      70       0.11201660
      71       0.10249667
      72      -0.08431058
      73      -0.02206055
      74       0.08811403
      75      -0.49404687
      76      -0.42898887
      77      -0.06727757
      78      -0.02241439
      79       0.00141117
      80      -0.22963779
      81      -0.02347212
      82       0.94734446
      83       0.18916398
      84       0.13454540
      85      -0.03868780
      86       0.12189706
      87       0.11818539
      88      -0.04902609
      89      -0.10429308
      90       1.12774823
      91       0.15510860
      92       0.52451078
      93       0.05767768
      94      -0.29617288
      95      -0.11414968
      96       0.46456536
      97       0.48069292
      98      -0.44696579
      99       0.01690299
     100       0.92429399
     101      -0.35619242
     102       0.55796542
     103      -0.07852023
     104       0.01936903
     105       0.03246701
     106       0.05121926
     107      -0.02149446
     108      -0.03990422
     109       0.09039403
     110       0.03939384
     111      -0.13078606
     112       0.04897103
     113      -0.05560283
     114       0.09571953
     115      -0.05651503
     116      -0.04661771
     117       0.03552627
     118       0.14105906
     119       0.92106559
     120      -0.42235887
     121       0.02774193
     122       0.80258984
     123       0.02039571
     124       0.02139512
     125      -0.15690364
     126       0.19484015
     127       0.04799579
     128      -0.14096673
     129       0.00100281
     130       1.73990373
     131       0.49832456
     132      -0.15059243
     133      -0.17852687
     134       0.75214169
     135      -0.33287449
     136       0.26345981
     137       0.11207761
     138       0.38921728
     139       0.03884038
     140      -0.63155937
     141      -0.31012784
     142       0.25202012
     143      -0.74513846
     144      -0.14825277
     145       0.23774277
     146       0.04958190
     147       0.14783661
     148       0.38299390
     149       0.11267174
     150       1.29088036
     151      -0.95212143
     152       0.03482171
     153      -0.65041967
     154      -1.41765111
     155      -0.27119869
     156      -0.07086010
     157       0.05553622
     158       0.10521142
     159      -0.02941121
     160       0.05090712
     161       0.07064295
     162       0.64572861
     163       0.02113215
     164      -0.51934572
     165      -0.10065652
     166       0.27844881
     167       0.46373915
     168      -0.37871620
     169      -0.02825524
     170       0.36200997
     171       0.31850400
     172      -1.11755980
     173      -0.86226294
     174      -0.27560308
     175       0.38477869
     176      -0.27660749
     177       0.40062155
     178       0.25735548
     179      -0.02729781
     180       0.32512284
     181      -0.78678725
     182       0.24916731
     183      -0.42341656
     184      -0.44273415
     185       0.21194929
     186      -0.54948755
     187       0.30663845
     188      -0.69491750
     189       0.77251405
     190      -0.03305236
     191      -0.47077650
     192       0.85335281
     193      -1.23549246
     194       0.57150000
     195      -0.49824222
     196       4.13556484
     197      -0.03254121
     198       0.16685475
     199       0.02995705
     200       0.01028436
     201      -0.12067442
     202       0.69547498
     203      -0.21789068
     204       0.08967788
     205       0.04919356
     206       0.02166245
     207      -0.05160321
     208       0.07583086
     209      -0.01012357
     210       0.17336347
     211      -0.03209736
     212      -0.05670918
     213       0.07603316
     214      -0.24468074
     215      -0.11266921
     216       0.06738849
     217       0.09089107
     218       0.01419696
     219      -0.32939243
     220      -0.00613409
     221      -0.41296670
     222       0.08468680
     223       0.01359235
     224       0.13348463
     225      -0.16821951
     226       0.06925034
     227       0.06650965
     228       0.01536733
     229      -0.08081587
     230       0.62899093
     231       0.04089305
     232       0.01597058
     233      -0.04990027
     234       0.03402003
     235       0.10148585
     236      -0.47538909
     237       0.00781369
     238       0.12203059
     239      -0.04853283
     240      -0.02114544
     241      -0.06853744
     242       0.26664356
     243      -0.09766310
     244      -0.04882275
     245      -0.05584799
     246       0.17285978
     247      -0.09718575
     248       0.64270879
     249       0.01974197
     250       0.24130640
     251       0.11091207
     252      -0.09178008
     253       0.02247379
     254      -0.04272516
     255       0.00655990
     256      -0.02923784
     257       0.06201580
     258       0.00219500
     259      -0.32188148
     260      -1.90717165
     261      -0.05792593
     262      -0.09487896
     263      -0.04133936
     264      -0.17747878
     265       0.02284747
     266      -1.39622583
     267      -0.07831233
     268       0.04358802
     269      -0.15687188
     270       0.26657529
     271      -0.00155089
     272       0.25910253
     273       0.09150634
     274      -0.01085003
     275       0.80555883
     276      -0.87905475
     277      -0.12307629
     278       0.11164497
     279      -0.10586280
     280      -0.29695689
     281      -0.30562931
     282      -0.29847713
     283       0.15040279
     284      -0.01844947
     285      -0.10524222
     286       0.33889398
     287      -0.00233018
     288      -0.28189507
     289       0.48422532
     290       1.41843262
     291       0.23884196
     292      -0.38411158
     293      -0.11448688
     294      -0.50713200
     295       0.01047492
     296       0.44505585
     297       0.08671594
     298       0.14614575
     299      -0.25550483
     300      -0.62589216
     301      -0.12157865
     302      -1.08843635
     303       0.25831247
     304      -1.02785669
     305      -0.11034801
     306       0.01551568
     307      -0.17371729
     308      -1.45905429
     309       0.03192036
     310       1.11377578
     311       0.00064790
     312      -0.43446542
     313       0.20376957
     314       3.19750370
     315      -0.12359043
     316       0.06128422
     317       0.03133924
     318       0.00979179
     319       0.14651079
     320      -0.11116295
     321       0.21566950
     322      -0.36187742
     323      -0.08486810
     324       0.06294955
     325       0.04785051
     326      -0.00264400
     327       0.00957246
     328       0.19489094
     329      -0.44428888
     330       0.82162397
     331       0.02271312
     332       0.13495545
     333       0.02419888
     334       0.25433028
     335      -0.17994254
     336       0.16652333
     337      -0.32181522
     338      -1.55431741
     339       0.10415870
     340       0.15674935
     341      -0.08260593
     342      -0.08794746
     343       0.02615994
     344      -0.70040583
     345       0.02239961
     346       0.23764025
     347      -0.10130827
     348       0.12850599
     349      -0.00810304
     350      -0.29342337
     351       0.23633559
     352       2.09054255
     353       1.13403137
     354      -0.81384808
     355       0.69201837
     356      -0.44493382
     357       0.04248811
     358      -0.30369443
     359      -0.01862144
     360      -0.46983490
     361      -0.00032185
     362       0.43928401
     363       0.09621685
     364      -0.05131754
     365      -0.07194389
     366      -0.49407756
     367      -0.08050691
     368      -1.15762319
     369       0.20978323
     370      -0.75030522
     371      -0.02778489
     372       0.59603016
     373       0.04067381
     374      -1.22854542
     375       0.00666623
     376       0.99868772
     377      -0.03351286
     378      -0.50147700
     379      -0.02014467
     380       2.25855378
     381       0.25456855
     382      -0.38635900
     383       0.56465337
     384      -0.41483893
     385      -0.29517565
     386       0.03318680
     387       0.01476846
     388      -0.03060575
     389       0.08368263
     390       0.05727398
     391       0.12847556
     392      -0.60778611
     393      -0.11747322
     394       0.27538354
     395      -0.04602949
     396      -0.03271600
     397       0.01486602
     398       0.05585064
     399       0.15211838
     400      -0.06524858
     401       0.11511606
     402       0.16906983
     403      -1.65399417
     404      -0.13838901
     405       0.06791220
     406      -0.02681327
     407       0.01597388
     408      -0.41411127
     409      -0.00150819
     410       0.04868396
     411       0.03795236
     412      -0.12043878
     413       0.51861659
     414      -0.03664796
     415       0.08402897
     416      -0.04720455
     417      -0.12893299
     418       0.07797014
     419      -0.06567055
     420       0.10380483
     421      -0.00287153
     422       0.16522316
     423       0.36506439
     424       0.10335338
     425      -0.09481994
     426       0.03986692
     427       0.05629682
     428       0.65427128
     429      -0.15266867
     430       0.10292838
     431       0.10859338
     432      -0.03939253
     433      -0.13603567
     434      -0.32766327
     435       0.16610585
     436      -0.10892631
     437      -0.61630059
     438       1.87294947
     439      -0.16524191
     440       0.09352492
     441       0.03086952
     442      -0.05530289
     443      -1.87729957
     444       0.14280543
     445      -0.08402453
     446      -0.01669454
     447       0.15108835
     448       0.56629402
     449       0.03308759
     450      -0.10016051
     451       0.00736590
     452       0.01369350
     453      -1.17994736
     454       0.22869073
     455      -0.20291052
     456       0.16920224
     457       0.24275770
     458      -2.31414196
     459       0.00386677
     460      -0.13094139
     461       0.05300706
     462      -0.49429596
     463       0.46149766
     464       0.23211801
     465       0.25322915
     466      -0.03863278
     467      -2.16191389
     468       0.62503699
     469       0.04886459
     470      -0.36978431
     471       0.21364242
     472      -0.16709705
     473       0.30208088
     474       0.06470790
     475       0.54582277
     476      -0.29903464
     477      -1.45392160
     478       0.15046369
     479       0.14002872
     480       0.13813640
     481       0.30070450
     482       0.13275885
     483      -0.19465004
     484      -0.16821470
     485      -0.05253597
     486       0.29026797
     487      -0.23158424
     488       0.16110243
     489      -0.08657039
     490       0.21561671
     491       0.72745216
     492      -0.32974317
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.01695333
       2       0.04279095
       3       0.24855727
       4      -0.09111749
       5       0.18722354
       6      -0.07981230
       7      -0.00345779
       8       0.01586170
       9       0.03464585
      10      -0.01864061
      11      -0.06066567
      12       0.10285026
      13      -0.01060122
      14      -0.11445888
      15      -0.18892852
      16      -0.42443509
      17      -0.30753565
      18       0.11763635
      19       0.08429277
      20      -0.32350374
      21       0.00747132
      22      -0.14562895
      23      -0.04468329
      24      -0.01077302
      25       0.98849276
      26      -0.19823893
      27      -0.12477646
      28       0.03016441
      29      -0.16761937
      30       1.58582052
      31       1.51859607
      32      -0.10388913
      33       0.09131635
      34      -0.18887097
      35      -0.04828490
      36       0.04183500
      37       0.10388206
      38       0.08133139
      39       0.17224153
      40      -0.25208165
      41      -0.36590990
      42       0.08424060
      43      -0.10673469
      44       0.87232081
      45       0.66805575
      46      -0.23760528
      47      -1.28436848
      48       0.55106836
      49      -0.52832189
      50       1.78776879
      51       4.88273764
      52      -0.06583914
      53      -0.03195005
      54       0.35542728
      55      -0.27685188
      56       0.41258227
      57       0.03915700
      58       0.25929163
      59       0.02468705
      60       0.06075885
      61       0.20294191
      62      -0.28568253
      63       0.08889563
      64       0.21093312
      65       0.02864366
      66       0.08794242
      67      -0.13175979
      68       0.34796964
      69      -0.08364523
      70       0.06592086
      71       0.39949643
      72       0.10282543
      73       0.01783771
      74       0.00381158
      75       0.50702432
      76       0.20535358
      77       0.35943314
      78       0.27341171
      79      -0.08241061
      80      -0.46428472
      81      -0.19436305
      82       0.45195345
      83       0.32692362
      84      -0.40842593
      85      -0.10984604
      86       0.17146737
      87      -0.02320702
      88      -0.03889000
      89       0.01135304
      90      -0.35969653
      91       0.01339752
      92       0.23049903
      93       0.11458288
      94       0.81869783
      95       0.22193124
      96       0.07940153
      97      -0.78539409
      98       0.34679453
      99      -0.18285221
     100      -0.77110484
     101      -0.13480958
     102      -2.52009332
     103       0.07989839
     104       0.00719315
     105       0.01496011
     106       0.01354015
     107       0.01991980
     108       0.04084166
     109      -0.06053868
     110      -0.10448253
     111       0.03535023
     112       0.14396939
     113      -0.11152848
     114       0.52127461
     115       0.03019994
     116      -0.05810292
     117      -0.17576472
     118       0.03449229
     119       0.49476194
     120      -0.50310206
     121      -0.16890721
     122       0.06117317
     123       0.01582024
     124      -0.21585205
     125       0.02451363
     126       0.07960129
     127       0.15672588
     128      -0.25357955
     129       0.27455998
     130      -1.41924670
     131       0.71312472
     132       0.01693326
     133      -0.04924415
     134       0.72245855
     135      -0.22582871
     136      -0.43043054
     137       0.26452440
     138      -0.40185171
     139       0.36790428
     140       0.05026765
     141       0.27073061
     142       0.77699203
     143       0.07593955
     144       0.05368777
     145      -0.43152421
     146       0.45597753
     147      -0.17355702
     148      -0.59922247
     149       0.05486579
     150      -0.35081519
     151       1.16527603
     152      -0.50497532
     153       0.60340137
     154       1.63072821
     155      -0.05827403
     156      -0.07267703
     157       0.34788315
     158      -0.05483660
     159      -0.03480369
     160      -0.09395257
     161      -0.10627903
     162       0.56641907
     163       0.13538428
     164       0.65943204
     165       0.12608173
     166       0.24952436
     167      -0.25531141
     168       0.07780421
     169      -0.17324617
     170       0.01317210
     171      -0.89079118
     172      -1.76412704
     173      -0.76657311
     174      -0.00504911
     175      -0.40815255
     176       0.26974212
     177       0.21936841
     178      -0.17810961
     179       0.02748593
     180       0.63980949
     181       0.70665916
     182      -0.53967940
     183       0.39185121
     184      -0.13285819
     185      -0.45163167
     186       0.58353541
     187      -1.95125712
     188       0.45340174
     189       0.00197588
     190       0.93434263
     191      -0.03803678
     192      -5.59612965
     193       1.35905658
     194      -1.80078181
     195       1.92094574
     196      -1.11487925
     197       0.08627017
     198      -0.18209452
     199      -0.64630205
     200       0.07979998
     201      -0.03689468
     202      -0.36487967
     203      -0.10289164
     204      -0.07008804
     205       0.00939159
     206       0.28693773
     207       0.02408776
     208      -0.01451146
     209      -0.09288371
     210      -0.05998409
     211      -0.11379532
     212       0.20623964
     213      -0.14903213
     214       0.86304950
     215       0.06915039
     216      -0.06280208
     217      -0.38464006
     218       0.34018864
     219      -0.09878234
     220      -0.20657461
     221       0.18222757
     222       0.23772382
     223       0.15483693
     224      -0.31872617
     225       0.32434493
     226      -0.07238022
     227      -0.81627036
     228       0.06605439
     229      -0.30167253
     230       0.25692225
     231      -0.05660376
     232      -0.05262692
     233      -0.00344305
     234      -0.01383820
     235       0.21030960
     236      -0.72522309
     237       0.15932410
     238      -0.11844669
     239       0.10333273
     240      -0.26257902
     241      -0.18461280
     242      -0.11070000
     243       0.02222540
     244       0.08248926
     245       0.07255233
     246      -0.03354203
     247       0.03471119
     248      -0.44618130
     249       0.06046875
     250      -0.09326267
     251      -0.02090692
     252      -0.28777142
     253      -0.09915768
     254      -0.04795789
     255      -0.03984474
     256       0.16858118
     257      -0.05539661
     258       0.27289870
     259       0.06550646
     260       1.98332409
     261       0.08860905
     262      -0.10242977
     263      -0.01970822
     264      -0.05389210
     265      -0.26988261
     266       2.23060068
     267      -0.14147706
     268       0.74772227
     269      -0.51473919
     270      -0.15667410
     271       0.07424619
     272      -0.23422218
     273       0.39018803
     274       0.19887759
     275      -0.63755872
     276      -1.03403065
     277      -0.03679990
     278       0.45122121
     279      -0.08412784
     280       0.34206052
     281      -0.23574811
     282       0.58176423
     283       0.08533956
     284      -0.25319268
     285      -0.25056579
     286      -0.12062098
     287       0.00925075
     288       0.45186540
     289      -0.19417338
     290      -1.69894308
     291       0.05178752
     292       0.02999415
     293       0.37467110
     294       0.56691515
     295       0.00636218
     296      -0.41410562
     297       0.03075464
     298       0.01092184
     299      -0.13085691
     300       0.72647696
     301      -0.01892339
     302       2.05399172
     303      -0.17526239
     304       0.45694856
     305       0.55225027
     306       0.50576044
     307       0.66343090
     308       1.00768917
     309      -0.02994205
     310      -1.18019419
     311      -0.32716449
     312       1.12735631
     313      -0.59873444
     314      -3.34794717
     315       0.17665861
     316      -0.01427400
     317      -0.02974168
     318       0.09797721
     319      -0.05093027
     320      -0.11942450
     321      -0.17760395
     322       0.40997801
     323       0.07823354
     324      -0.11416084
     325      -0.02415200
     326       0.18124413
     327      -0.03308442
     328      -0.02386069
     329      -0.21925813
     330       0.58034495
     331       0.03204960
     332      -0.04848342
     333       0.14605677
     334      -0.24173554
     335       0.12087148
     336      -0.05940075
     337      -0.14616086
     338      -0.25160552
     339      -0.07299759
     340       0.29218897
     341      -0.00793598
     342       0.57164740
     343       0.04694080
     344       0.66522466
     345       0.06179336
     346      -0.56872592
     347      -0.10910788
     348      -0.18956788
     349       0.23889711
     350       0.43135420
     351      -0.11254475
     352      -2.33916723
     353       0.12973974
     354      -0.24433151
     355      -0.70377178
     356       0.63312200
     357       0.19855314
     358       0.21008297
     359       0.02600627
     360       0.59859643
     361      -0.13835290
     362      -0.44382521
     363       0.03225357
     364       0.00261667
     365      -0.01327646
     366       0.55927292
     367       0.10096521
     368       1.75344436
     369      -0.14843796
     370       0.22415184
     371      -0.33632765
     372      -0.45847637
     373      -0.02174529
     374       0.56324515
     375       0.13653469
     376      -0.92195675
     377      -0.33416512
     378       1.55295973
     379       0.06176583
     380      -2.48049337
     381      -0.47220765
     382       0.45091138
     383       0.26190775
     384      -0.11568126
     385       0.01218450
     386       0.11035670
     387       0.04191968
     388       0.24074382
     389       0.01927480
     390      -0.46196868
     391       0.03324749
     392       0.30450519
     393       0.10354651
     394       0.01741769
     395      -0.03656534
     396       0.50768856
     397      -0.03605315
     398       0.17537199
     399      -0.47895612
     400       0.35992377
     401      -0.32939313
     402      -0.18816207
     403       0.13589656
     404      -0.02658158
     405       0.03550180
     406      -0.05207603
     407      -0.09833022
     408       0.13537925
     409      -0.10613859
     410       0.14695182
     411      -0.12349286
     412       0.14966506
     413      -0.65290468
     414      -0.00304914
     415      -0.01477335
     416       0.05447064
     417      -0.06649864
     418       0.59687096
     419       0.00098205
     420      -0.16443039
     421       0.09975039
     422       0.01320775
     423      -0.13558958
     424      -0.04528751
     425       0.04279255
     426       0.03063416
     427       0.15514353
     428      -0.85285305
     429       0.02531376
     430      -0.00396599
     431       0.14340050
     432      -0.19789346
     433       0.05823490
     434       0.07548001
     435      -0.01851872
     436      -0.04562720
     437       0.60462364
     438       0.33901838
     439       0.02388274
     440      -0.05793080
     441       0.00075487
     442       0.12740999
     443       2.16617286
     444       0.01990466
     445      -0.05569703
     446      -0.19990227
     447      -0.17259940
     448      -0.59969680
     449      -0.13730490
     450       0.20147412
     451       0.09327497
     452      -0.13823921
     453       0.99900481
     454      -0.10913187
     455       0.18237304
     456      -0.15553688
     457      -0.32730641
     458       2.35036097
     459      -0.38429446
     460      -0.27962862
     461      -0.14324915
     462       0.82548576
     463       0.14204080
     464      -0.12978762
     465      -0.06697998
     466       0.15745061
     467       1.06137114
     468       0.07930100
     469      -0.15743233
     470       0.05192725
     471      -0.16525540
     472       1.18447985
     473      -0.24849116
     474      -0.04897100
     475      -0.47362889
     476       0.30375706
     477      -0.04879308
     478       0.64497542
     479       0.19955843
     480      -0.17196254
     481      -0.12088896
     482      -0.12984289
     483       0.14070890
     484       0.05934942
     485       0.13544488
     486      -0.17821530
     487      -0.10623146
     488      -0.12906234
     489      -0.33290152
     490       0.05684636
     491      -0.60172938
     492       1.16259059

 Residual norm when dim(red L) =  15
 NEO root     CSF        orbital          total
    1     0.00003871     0.00001448     0.00004133 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  1

   1.988008599   1.976718935   0.023720661   0.012979398

 Symmetry  2

   1.980681108   0.018399242

 Symmetry  3

   1.975007111   0.024484945


 <<< MACRO ITERATION  4 >>>
 --------------------------

 Total MCSCF energy       :      -76.168770402289368       (MACRO    4)
 - Dielec. solvation ener.:       -0.012753014673001

 Norm of total gradient   :        0.000069249055
 -    of CI gradient      :        0.000035369358
 -    of orbital gradient :        0.000059535201
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.01524520
       2       0.13605217
       3       0.07322544
       4       0.00043594
       5       0.06540312
       6      -0.01039458
       7       0.02096450
       8      -0.04244902
       9      -0.08516562
      10      -0.01724937
      11      -0.00404243
      12       0.00375593
      13       0.00371853
      14      -0.01722698
      15      -0.06029716
      16      -0.76920918
      17      -0.01149464
      18       0.03594267
      19      -0.01885713
      20      -0.36059032
      21       0.31242981
      22      -0.20444478
      23      -0.11576588
      24      -0.58579626
      25       0.12161758
      26      -0.34613626
      27       0.02858019
      28      -0.03416381
      29      -0.01402104
      30       0.04137614
      31      -0.18309681
      32      -0.10205444
      33       0.00991985
      34      -0.02776998
      35       0.00148562
      36       0.00328151
      37       0.02686971
      38      -0.08475411
      39      -0.06127070
      40      -0.58260755
      41      -0.19687937
      42       0.01877804
      43      -0.00203969
      44       0.01495401
      45      -0.13466797
      46      -0.05890511
      47      -0.08373066
      48      -0.05488890
      49      -0.02464153
      50      -0.32394382
      51      -0.47676791
      52      -0.00644130
      53      -0.01480159
      54      -0.01788983
      55      -0.06775286
      56       0.92565463
      57       0.01607830
      58       0.05317042
      59       0.02473085
      60       0.02486368
      61       0.00493175
      62      -0.41924909
      63       0.04261008
      64       0.05234244
      65       0.01664789
      66       0.02728133
      67      -0.01837759
      68       0.00221949
      69      -0.09777360
      70      -0.00993473
      71       0.25577242
      72       0.11661039
      73      -0.05021286
      74       0.01285090
      75       0.07210369
      76       0.07483278
      77      -0.01845039
      78       0.01412409
      79       0.13073781
      80      -0.76837070
      81      -0.18435892
      82       0.45563255
      83       0.16933178
      84      -0.13254027
      85      -0.00327530
      86       0.00109814
      87      -0.00412615
      88       0.08338986
      89       0.03769138
      90       0.06741301
      91      -0.01863427
      92      -0.09145788
      93      -0.05181316
      94       0.21568642
      95       0.01026230
      96      -0.02818169
      97      -0.04540805
      98       0.01146015
      99       0.03452371
     100      -0.09144834
     101       0.09428547
     102       0.53012570
     103       0.03802888
     104      -0.00415771
     105       0.01618578
     106      -0.00097773
     107       0.01702669
     108       0.00090752
     109      -0.02269475
     110      -0.05333902
     111       0.01984152
     112      -0.01348401
     113       0.00666926
     114      -0.00630260
     115       0.01526550
     116       0.00235991
     117       0.01860619
     118      -0.02090228
     119       0.56898641
     120      -1.32825691
     121      -0.19555213
     122       0.25375431
     123       0.05221831
     124      -0.06183116
     125      -0.09669047
     126      -0.02509710
     127       0.12154135
     128       0.13271788
     129       0.06505542
     130       0.03152152
     131       0.42416504
     132      -0.24649299
     133       0.70421903
     134       1.68511118
     135      -0.54666395
     136      -0.02129967
     137       0.07207058
     138       0.20582379
     139      -0.01244456
     140      -0.04679633
     141      -0.09578166
     142      -0.09318006
     143      -0.01416502
     144      -0.04127721
     145       0.01005154
     146       0.17817748
     147      -0.04153540
     148      -0.04643251
     149      -0.04572140
     150       0.10354680
     151      -0.00468488
     152       0.01814701
     153       0.00944776
     154      -0.01376413
     155       0.00043092
     156      -0.02084019
     157       0.08135157
     158      -0.01918007
     159      -0.01750911
     160       0.00109263
     161      -0.05479392
     162       0.07004748
     163       0.00233350
     164       0.11872161
     165       0.01375575
     166       0.01275792
     167      -0.06694318
     168      -0.01155943
     169      -0.00153798
     170      -0.08816230
     171      -0.04693034
     172       0.55375808
     173      -0.42680389
     174      -0.12492686
     175      -0.01918658
     176       0.19110711
     177      -0.00246529
     178      -0.04864897
     179      -0.04572073
     180       0.13094413
     181      -0.01463140
     182       0.01452057
     183      -0.01787303
     184       0.00026838
     185      -0.10421158
     186       0.03755945
     187       0.05687076
     188       0.02656102
     189      -0.05779727
     190      -0.04411173
     191      -0.05189695
     192       0.63778964
     193      -0.33583843
     194       0.13323193
     195      -0.37737152
     196      -0.52492845
     197       0.00556530
     198      -0.03239838
     199      -0.00994391
     200       0.01452243
     201       0.01639788
     202      -0.04410717
     203      -0.08232492
     204      -0.02618656
     205       0.01583855
     206       0.15361219
     207       0.05895914
     208      -0.06059689
     209      -0.02143138
     210      -0.08821326
     211      -0.05278909
     212       0.02669773
     213      -0.02940858
     214       0.06178519
     215       0.05016910
     216      -0.02186896
     217      -0.03400713
     218      -0.00375214
     219       0.04625736
     220      -0.14848361
     221       0.09694794
     222       0.07797246
     223       0.06797028
     224      -0.07699049
     225       0.07983572
     226      -0.01831088
     227      -0.01137686
     228      -0.01039267
     229       0.01751166
     230      -0.06692156
     231      -0.04546843
     232      -0.00261602
     233      -0.04692039
     234      -0.00471390
     235       0.18480277
     236      -0.91610945
     237       0.04363996
     238      -0.08534646
     239       0.01168534
     240      -0.10664910
     241      -0.06896472
     242       0.06502532
     243       0.02070132
     244       0.06975995
     245       0.02999578
     246      -0.03389524
     247       0.02015019
     248      -0.10776367
     249       0.02395531
     250      -0.06927819
     251      -0.02523736
     252      -0.06157550
     253      -0.02890036
     254       0.03686174
     255      -0.01065657
     256       0.03681098
     257       0.00143156
     258       0.08196259
     259       0.06667403
     260      -0.06068421
     261       0.01301957
     262       0.06469552
     263       0.01039346
     264       0.05218445
     265       0.01062436
     266      -0.27171006
     267      -0.41028006
     268       0.86802925
     269      -0.57875498
     270      -0.64271363
     271      -0.06590181
     272       0.21245624
     273       0.57843859
     274       0.34777954
     275      -0.93703883
     276      -1.01510089
     277      -0.02818966
     278       0.09521252
     279       0.01172265
     280       0.05872237
     281      -0.07488826
     282       0.21151758
     283      -0.03287566
     284      -0.06208553
     285      -0.09085170
     286       0.09799138
     287      -0.00551317
     288      -0.07983492
     289      -0.10764683
     290      -0.00781132
     291      -0.00021472
     292      -0.06334961
     293       0.10990436
     294      -0.02728213
     295      -0.01922317
     296      -0.03971668
     297       0.01454278
     298      -0.00286965
     299      -0.00857003
     300      -0.01447950
     301       0.02414818
     302       0.00865492
     303      -0.04995745
     304       0.00571170
     305       0.02548621
     306      -0.07017308
     307      -0.00647607
     308       0.11496537
     309       0.00772721
     310      -0.01699317
     311       0.03257367
     312      -0.17872854
     313       0.04953342
     314       0.29249933
     315       0.02293696
     316      -0.00150335
     317      -0.01384853
     318       0.01294064
     319      -0.03123720
     320       0.01468375
     321      -0.04647779
     322       0.07776985
     323       0.01877665
     324       0.01036389
     325      -0.00803115
     326      -0.01593985
     327      -0.00735094
     328      -0.02460649
     329      -0.21913419
     330       0.69708737
     331      -0.04324061
     332      -0.20971118
     333       0.10348834
     334       0.15804758
     335       0.19437203
     336       0.12211004
     337      -0.16844422
     338      -0.72677970
     339      -0.01733430
     340      -0.03833684
     341      -0.05747858
     342       0.12693036
     343       0.02690333
     344       0.21305361
     345       0.00408271
     346      -0.08237262
     347      -0.07576404
     348       0.11635796
     349       0.05558121
     350      -0.07388753
     351      -0.04868961
     352      -0.04161634
     353       0.27393991
     354      -0.45361978
     355      -0.09757818
     356       0.05678314
     357       0.07828119
     358       0.03535877
     359      -0.00338674
     360      -0.03209307
     361      -0.00712494
     362      -0.04428433
     363      -0.00129000
     364      -0.03728152
     365       0.00002619
     366       0.01550638
     367       0.02049592
     368       0.06501764
     369      -0.04677023
     370      -0.05003771
     371       0.00309121
     372      -0.09668687
     373      -0.00610334
     374       0.04359016
     375       0.00079779
     376      -0.01373144
     377       0.01509430
     378      -0.16746425
     379       0.00653287
     380       0.24543734
     381      -0.11100923
     382       0.06712937
     383      -0.01340608
     384      -0.02499721
     385       0.02435804
     386       0.07013196
     387       0.00789035
     388       0.14712902
     389      -0.00292088
     390      -0.10669202
     391       0.00345894
     392       0.07668401
     393       0.04583599
     394      -0.00161974
     395       0.00208556
     396       0.06271902
     397      -0.00579819
     398      -0.01957944
     399      -0.37486361
     400       0.27674347
     401      -0.22110960
     402      -0.12556462
     403       0.07877736
     404      -0.02447377
     405       0.00625930
     406      -0.01590532
     407      -0.08160079
     408      -0.08010348
     409      -0.04639707
     410       0.04062780
     411      -0.05642083
     412       0.06712977
     413       0.07853846
     414       0.04794842
     415      -0.05889521
     416       0.05147760
     417      -0.01089137
     418       0.23108598
     419       0.00235751
     420      -0.08091798
     421       0.02506456
     422      -0.01281473
     423      -0.11010527
     424      -0.03736950
     425       0.02577083
     426      -0.02049230
     427       0.02857910
     428      -0.05686940
     429       0.05765924
     430      -0.04462897
     431       0.02638253
     432      -0.05127142
     433      -0.08194975
     434       0.03720307
     435      -0.00594729
     436       0.02538743
     437       0.13180514
     438       0.17479649
     439       0.02713481
     440      -0.01948567
     441      -0.00981089
     442       0.00819165
     443       0.05543316
     444      -0.03208747
     445       0.00371201
     446      -0.02753890
     447      -0.07187963
     448       0.06023355
     449      -0.02876047
     450       0.04139707
     451      -0.00042204
     452      -0.01185332
     453       0.04562797
     454      -0.04984335
     455       0.04386070
     456      -0.03202363
     457      -0.05515841
     458       0.08509675
     459      -0.23013143
     460      -0.09836512
     461      -0.16901275
     462       0.13677990
     463       0.21002109
     464      -0.06805742
     465      -0.05263282
     466      -0.01062465
     467       0.06893534
     468       0.08094897
     469      -0.03753518
     470      -0.05990662
     471      -0.02925955
     472       0.17238054
     473       0.07266746
     474      -0.02917034
     475      -0.10494234
     476       0.05137182
     477       0.02500896
     478       0.10146263
     479       0.02495480
     480      -0.07926385
     481      -0.06856393
     482      -0.08633400
     483       0.07520401
     484       0.02962576
     485       0.05466675
     486      -0.09775447
     487       0.00238239
     488      -0.05475136
     489      -0.08198409
     490      -0.00834342
     491      -0.03581665
     492       0.08069283
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1       0.01543204
       2       0.63507506
       3       0.10303135
       4      -0.05317863
       5       0.05076407
       6      -0.05151042
       7       0.01238672
       8       0.02169116
       9      -0.02878659
      10       0.10224834
      11       0.00265639
      12       0.00650101
      13       0.02030503
      14       0.00602474
      15      -0.02012253
      16       0.16938801
      17       0.03152641
      18       0.10801873
      19       0.27683888
      20       0.86682856
      21       1.60148080
      22       0.49689643
      23      -0.28802617
      24       0.46075654
      25       0.16117028
      26       1.19253291
      27      -0.00508053
      28       0.14188068
      29       0.03314762
      30       0.43527891
      31      -0.07602022
      32       0.03379824
      33       0.26902457
      34      -0.01749990
      35       0.04881769
      36       0.14555335
      37      -0.00495985
      38      -0.07169117
      39      -0.02875687
      40       1.35169545
      41      -0.22921598
      42       0.13207785
      43       0.07014141
      44       0.52287968
      45      -0.04249269
      46       0.24854656
      47      -0.03900028
      48       0.02268921
      49       0.07282568
      50      -0.05252585
      51      -0.09110318
      52       0.05042739
      53       0.05430538
      54       0.04299339
      55      -0.02590791
      56       0.36594330
      57      -0.26843478
      58       0.04468647
      59      -0.01621821
      60      -0.24758867
      61       0.01786614
      62       0.05087263
      63       0.19809591
      64       0.04923293
      65       0.18780318
      66       0.13148628
      67      -0.07268466
      68       0.05005753
      69      -0.03899964
      70       0.03930888
      71       0.33823769
      72      -0.05752737
      73       0.07922257
      74       0.04706247
      75       0.00282453
      76      -0.04446756
      77      -0.10604746
      78      -0.01434469
      79      -0.14819023
      80       0.84562732
      81      -0.15761282
      82       0.94743600
      83       0.13091330
      84       0.17472599
      85      -0.04770199
      86      -0.05751551
      87      -0.08867150
      88       0.47824051
      89       0.03020032
      90      -0.35727706
      91      -0.03529238
      92       1.36098199
      93       0.07390850
      94      -0.43649475
      95      -0.12959068
      96      -0.02923798
      97      -0.01491417
      98       0.00472540
      99      -0.00026632
     100      -0.25661695
     101      -0.01399336
     102       0.09006077
     103      -0.01856025
     104      -0.00245288
     105      -0.00634782
     106       0.05310547
     107      -0.06167351
     108       0.05490121
     109      -0.01616070
     110       0.14311865
     111      -0.00193229
     112      -0.01710927
     113       0.01828277
     114      -0.02673450
     115       0.00270039
     116       0.00982094
     117       0.01152143
     118      -0.01571709
     119      -0.00288206
     120       0.55383799
     121      -0.02939497
     122       1.15487854
     123       0.01357724
     124       0.41934074
     125      -0.07245835
     126       0.00348607
     127      -0.17991937
     128       0.81468397
     129       0.29583449
     130      -0.36757161
     131      -0.71541154
     132       0.26914306
     133      -0.46196305
     134      -0.20323358
     135      -0.10264611
     136      -0.23984487
     137      -0.38173758
     138       0.67826158
     139      -0.36453638
     140       0.02700453
     141      -0.09794241
     142      -0.55299674
     143      -1.19992785
     144      -0.11498308
     145      -0.21389204
     146      -0.76759745
     147       1.06070135
     148      -0.07446031
     149      -0.06039994
     150      -0.35792510
     151       0.10435371
     152      -0.02840293
     153       0.12938569
     154       0.24114751
     155      -0.05351536
     156      -0.03767071
     157      -0.02917709
     158       0.02321187
     159      -0.05282013
     160      -0.02589950
     161      -0.11479349
     162       1.30577769
     163       0.06995785
     164      -0.71545968
     165      -0.05228404
     166      -0.04339401
     167       0.04158225
     168      -0.01695863
     169      -0.00303422
     170      -0.27754614
     171      -0.12103534
     172       0.13525412
     173      -0.97775202
     174       0.03318766
     175      -0.31572546
     176      -0.67751426
     177       0.84885291
     178      -0.09751882
     179      -0.21440397
     180      -0.56137426
     181       0.07382530
     182      -0.14898604
     183       0.14631344
     184       0.28406133
     185       0.39128706
     186       0.15330993
     187       0.05778674
     188       0.04086390
     189      -0.28331140
     190       0.05711168
     191      -0.00909344
     192       0.07899497
     193      -0.13897070
     194      -0.06731532
     195      -0.23782336
     196      -0.32042834
     197      -0.00318799
     198       0.00836525
     199       0.01803779
     200       0.00291916
     201       0.02355712
     202      -0.00001781
     203      -0.11761085
     204       0.07629230
     205      -0.09030241
     206       0.25347241
     207       0.03696941
     208       0.03387582
     209      -0.05896657
     210       0.10843323
     211       0.07138952
     212      -0.04864014
     213      -0.01006486
     214      -0.01211459
     215      -0.13584546
     216      -0.00564812
     217      -0.00289367
     218      -0.04880167
     219       0.07989228
     220      -0.03543103
     221       0.01210809
     222       0.18112867
     223      -0.10588352
     224      -0.04227369
     225      -0.01400970
     226       0.06911374
     227       0.01006885
     228      -0.07652297
     229      -0.01731130
     230       0.03878473
     231      -0.02502833
     232      -0.18837822
     233       0.00792640
     234       0.10328425
     235       0.12146746
     236       2.10173391
     237       0.03051314
     238      -0.02270589
     239       0.01142220
     240      -0.22997646
     241      -0.17518028
     242       0.17978935
     243      -0.02204349
     244      -0.01306584
     245       0.07278320
     246       0.32546746
     247       0.04470124
     248       0.05933092
     249       0.03973217
     250      -0.05451381
     251       0.04211265
     252      -0.14119871
     253      -0.02462584
     254       0.02623262
     255      -0.00797277
     256       0.07297141
     257       0.03342270
     258       0.17649040
     259       0.05152688
     260       0.76565245
     261       0.00560971
     262      -0.06148945
     263       0.01944442
     264      -0.01702828
     265       0.00073945
     266      -0.07529235
     267      -0.71016510
     268      -1.87395052
     269       0.24820670
     270       0.24420714
     271      -0.14383252
     272      -0.13670853
     273       0.73812500
     274      -0.30292467
     275       0.39083232
     276       1.64800775
     277      -0.10494594
     278      -0.50536666
     279       0.04141161
     280      -0.90933723
     281       0.00403647
     282       0.18800221
     283       0.09827011
     284       0.03210775
     285      -0.03142811
     286       0.01811614
     287      -0.10868549
     288       0.14191474
     289      -0.07344797
     290      -0.64849455
     291       0.22058006
     292      -0.57478067
     293       0.20521980
     294       0.19206305
     295       0.04828542
     296      -0.18730085
     297      -0.16489866
     298      -0.07957711
     299       0.15589155
     300       0.59283564
     301       0.00398383
     302       0.08448587
     303      -0.01694719
     304       0.09080797
     305       0.01107934
     306      -0.09641168
     307       0.00674047
     308       0.12512051
     309       0.01138973
     310      -0.20340498
     311       0.08378332
     312       0.04080763
     313       0.09588949
     314      -0.02627862
     315       0.00269669
     316       0.00239475
     317       0.00929772
     318      -0.05444499
     319      -0.04265784
     320       0.00598998
     321      -0.00751098
     322      -0.03885459
     323       0.03776653
     324       0.05234969
     325       0.02822820
     326       0.01217180
     327      -0.00111421
     328      -0.01309193
     329      -0.24997298
     330      -1.67958846
     331       0.04196357
     332       0.21883952
     333      -0.00339147
     334      -0.07959135
     335      -0.33647169
     336      -0.33692912
     337       0.26852566
     338       1.54470883
     339      -0.10149938
     340      -0.80758349
     341      -0.02366588
     342      -1.03965895
     343      -0.04120199
     344       0.10662856
     345       0.01773655
     346       0.00846453
     347      -0.02354479
     348       0.14548636
     349       0.06485226
     350       0.14761198
     351      -0.00188502
     352      -0.55795434
     353      -0.94564023
     354       0.94055051
     355      -0.36952402
     356       0.24926931
     357       0.16782316
     358      -0.51180344
     359       0.02683055
     360       0.19671475
     361      -0.01211394
     362      -0.17693117
     363      -0.04075654
     364      -0.17386147
     365       0.06583497
     366       0.68258599
     367       0.01196138
     368       0.11253929
     369       0.01118239
     370       0.08238026
     371      -0.02255527
     372      -0.11972615
     373      -0.02189704
     374       0.08478913
     375       0.00059043
     376      -0.22602042
     377      -0.00880588
     378       0.03472369
     379       0.02375168
     380      -0.04671682
     381      -0.49303517
     382       0.34513262
     383      -0.05559141
     384      -0.02370159
     385       0.28238216
     386      -0.20398735
     387       0.01855276
     388      -0.07051458
     389      -0.04563906
     390       0.03498682
     391       0.02144134
     392      -0.03124916
     393       0.09062633
     394       0.20799441
     395       0.00539860
     396      -0.02340185
     397      -0.00230131
     398       0.04997842
     399      -0.32093819
     400       0.24127216
     401      -0.13938396
     402      -0.17154350
     403      -1.74672656
     404       0.15402962
     405      -0.10822537
     406       0.02008365
     407       0.03535259
     408      -0.04920076
     409      -0.11155338
     410       0.05814085
     411      -0.04935532
     412       0.01945604
     413      -0.34619559
     414      -0.09491596
     415      -0.00201549
     416      -0.03104038
     417      -0.05470249
     418      -1.24611381
     419       0.23243025
     420      -0.13468087
     421       0.06986943
     422       0.16480967
     423       1.52344010
     424       0.10515848
     425      -0.04333035
     426      -0.05738566
     427       0.03374110
     428      -0.04219794
     429      -0.13485365
     430       0.10468700
     431       0.02794235
     432      -0.06687522
     433      -0.18518034
     434      -0.04126883
     435       0.05277067
     436      -0.06554487
     437      -0.05667434
     438      -0.67489724
     439       0.02951260
     440      -0.02431075
     441       0.00054051
     442       0.00331830
     443       0.17867415
     444       0.11559197
     445      -0.14453608
     446      -0.06705711
     447      -0.03743111
     448      -0.09364716
     449       0.05205639
     450      -0.02457861
     451       0.01941678
     452       0.00327382
     453       0.23816657
     454       0.14427925
     455      -0.09555149
     456       0.02626174
     457       0.03016770
     458       0.43360197
     459      -0.16312342
     460      -0.23381047
     461      -0.06774641
     462      -0.92137181
     463       0.58913312
     464       0.18673269
     465       0.10691530
     466       0.08149001
     467       0.33235725
     468      -0.27840080
     469      -0.10507249
     470      -0.16445115
     471       0.02651499
     472      -1.18755429
     473       0.65642493
     474       0.05467180
     475       0.06558940
     476       0.00988701
     477       0.43515722
     478      -0.35863454
     479       0.08108759
     480       0.05601658
     481      -0.20600647
     482       0.19092197
     483      -0.02403007
     484       0.12725617
     485      -0.11081311
     486      -0.02339958
     487      -0.08193943
     488       0.05001099
     489      -0.13131969
     490       0.11008521
     491      -0.72568860
     492       0.58484010
  LINEAR TRANSFORMED CONFIGURATION VECTOR
  **** AFTER PCMLNC **** 

               Column   1
       1      -0.04488490
       2       0.48225292
       3      -0.20234086
       4       0.00246622
       5      -0.19904219
       6      -0.05119963
       7      -0.05630439
       8       0.03590401
       9       0.20881423
      10      -0.02241078
      11      -0.01891638
      12      -0.03365829
      13      -0.02782942
      14       0.00767754
      15       0.16680197
      16       0.30498177
      17      -0.12567436
      18      -0.09013619
      19      -0.25706619
      20       0.07126836
      21       0.54077417
      22       0.01038100
      23       0.02177444
      24      -0.20559980
      25      -0.82093028
      26      -0.33153004
      27      -0.76493511
      28      -0.18962491
      29       0.02031674
      30      -1.04967280
      31       0.32209127
      32       0.04297104
      33      -0.10160001
      34       0.01556944
      35      -0.08350623
      36      -0.10560654
      37      -0.04634531
      38       0.18503204
      39       0.07771101
      40       0.13667396
      41      -0.59134759
      42      -0.26823230
      43      -0.02651679
      44      -1.22276479
      45       0.15306859
      46      -0.97557779
      47       0.17301882
      48       0.02220874
      49      -0.11138452
      50       0.35918489
      51       1.40997779
      52      -0.06900252
      53      -0.09020404
      54      -0.07020089
      55       0.16895231
      56      -0.09726057
      57       0.02329802
      58      -0.12063212
      59      -0.03160692
      60       0.16918797
      61       0.03255908
      62       0.29163079
      63      -0.02961733
      64      -0.16053550
      65      -0.22290327
      66      -0.14582353
      67       0.06661881
      68      -0.10514411
      69       0.20871927
      70      -0.03503706
      71      -0.17821351
      72      -0.07743410
      73      -0.05366762
      74      -0.06362103
      75      -0.19070619
      76      -0.08999272
      77       0.14344083
      78       0.04103821
      79       0.03634115
      80      -0.37154176
      81       0.14126376
      82      -0.83323287
      83      -0.11113893
      84      -0.08629244
      85       0.12900143
      86       0.11256205
      87       0.04172335
      88      -0.61663150
      89      -0.11417292
      90       0.20712162
      91       0.05500068
      92      -1.03973148
      93       0.06264055
      94      -0.43373806
      95       0.08967039
      96      -0.01490303
      97       0.15599329
      98      -0.08373911
      99      -0.04946883
     100       0.20030322
     101      -0.17199559
     102      -1.04382105
     103      -0.00375939
     104      -0.00712465
     105      -0.01745653
     106      -0.05256549
     107       0.02133478
     108      -0.04569671
     109       0.05426238
     110      -0.04389543
     111      -0.01878754
     112       0.03867342
     113      -0.03649986
     114       0.01056556
     115      -0.01957591
     116      -0.01631203
     117      -0.03833211
     118       0.05606819
     119       0.03147556
     120      -0.83851293
     121       0.08407245
     122      -0.83073684
     123       0.06153234
     124      -0.29694717
     125       0.13668062
     126       0.12506796
     127      -0.03800635
     128      -0.90105244
     129      -0.27966560
     130       0.19883063
     131       0.15661596
     132       0.67591341
     133      -0.48464145
     134       0.30807220
     135      -1.07707983
     136       0.62602808
     137       0.13390134
     138      -0.58152593
     139       0.61768711
     140      -0.31847028
     141      -0.04488300
     142       0.63608886
     143       0.97552238
     144      -0.06357790
     145       0.19863495
     146       0.60863108
     147      -0.90486286
     148       0.17368672
     149       0.22668197
     150      -0.13615257
     151       0.02417827
     152       0.03598749
     153      -0.11648596
     154      -0.02483614
     155       0.00573389
     156       0.04255537
     157      -0.09538579
     158      -0.00067798
     159       0.06268886
     160       0.02828257
     161       0.11205255
     162      -1.01459327
     163      -0.07201848
     164      -0.00887104
     165       0.03389004
     166      -0.09537212
     167       0.09612917
     168      -0.02211207
     169       0.06457473
     170       0.17903203
     171       0.19218839
     172      -1.09739749
     173       1.11753412
     174      -0.03989999
     175       0.32865040
     176       0.28333208
     177      -0.71874892
     178       0.24465251
     179       0.32558605
     180       0.11808923
     181       0.00605984
     182       0.15649311
     183      -0.16384094
     184      -0.14013388
     185       0.40709393
     186      -0.17492610
     187      -0.04615843
     188       0.22377411
     189      -0.12144234
     190       0.02035110
     191       0.18497762
     192      -1.04492974
     193       0.92821505
     194      -0.06196320
     195       1.15026985
     196       1.60281417
     197       0.03553465
     198       0.04724288
     199       0.04055349
     200      -0.04860836
     201      -0.04294186
     202       0.09632694
     203       0.01667872
     204       0.07337407
     205       0.02560214
     206      -0.11642463
     207      -0.04229983
     208       0.02310265
     209       0.06927175
     210       0.08303148
     211       0.01886794
     212      -0.02733282
     213       0.05053217
     214      -0.12480472
     215      -0.00639506
     216       0.03788287
     217       0.10765803
     218       0.01795439
     219       0.10300843
     220       0.11915966
     221      -0.14276147
     222      -0.08505551
     223       0.01207344
     224       0.13439340
     225      -0.08915214
     226      -0.05607284
     227       0.05640783
     228       0.03642274
     229       0.03997399
     230       0.06340800
     231       0.07635355
     232       0.02852742
     233       0.04042731
     234       0.00084779
     235       0.04381281
     236      -0.56208784
     237      -0.06806493
     238       0.15372011
     239      -0.05134298
     240       0.20261755
     241       0.07955408
     242      -1.49596336
     243      -0.03831935
     244      -0.18047340
     245      -0.07562471
     246      -0.25091788
     247      -0.04910053
     248      -0.10774513
     249      -0.05124137
     250       0.14840521
     251       0.02435156
     252       0.14722709
     253       0.04476194
     254      -0.05124982
     255       0.00722399
     256      -0.11668776
     257      -0.04573285
     258      -0.24472591
     259      -0.16636825
     260      -1.87628465
     261      -0.03255189
     262      -0.04833057
     263      -0.03663782
     264       0.03330144
     265      -0.01593754
     266       0.32625946
     267      -0.85827461
     268       0.29854442
     269       0.19262247
     270       0.37727473
     271       0.03585019
     272       0.15856075
     273       0.20897436
     274       0.43409221
     275       0.48585170
     276      -0.57349598
     277       0.04392973
     278       0.67139096
     279      -0.07524189
     280       0.72182814
     281      -0.06945330
     282      -1.48860270
     283      -0.01414871
     284       0.19867120
     285       0.06466367
     286       0.04252882
     287       0.12390848
     288      -0.10122926
     289       0.24331595
     290       1.71230284
     291       0.00140199
     292       1.81812596
     293      -0.37007855
     294      -0.11005450
     295       0.05215762
     296       0.23923194
     297      -0.00888763
     298       0.29740033
     299      -0.25393034
     300      -1.81167772
     301      -0.03835339
     302       0.07305723
     303       0.10471399
     304       0.26831318
     305      -0.10687980
     306       0.08659095
     307      -0.04031786
     308      -0.46250187
     309      -0.06028402
     310       0.22894973
     311      -0.14467162
     312       0.16420671
     313      -0.25867086
     314       0.03769343
     315      -0.01345789
     316      -0.02393997
     317      -0.00426208
     318       0.07843056
     319       0.06807724
     320      -0.04001106
     321       0.08128950
     322      -0.00706318
     323      -0.05469300
     324      -0.11429735
     325      -0.01396954
     326       0.02346110
     327       0.01470867
     328       0.05351711
     329      -0.10586316
     330       0.63867074
     331       0.08593007
     332       0.23827157
     333      -0.06284126
     334       0.14781345
     335       0.04740620
     336       0.54361428
     337       0.06283434
     338      -0.81567184
     339       0.08242295
     340       1.07384384
     341       0.09677553
     342       0.87015357
     343      -0.09097983
     344      -1.35538165
     345      -0.04111285
     346       0.19439025
     347       0.07713471
     348       0.03879790
     349      -0.15929494
     350      -0.11007380
     351       0.09789829
     352       1.54746403
     353       1.04022537
     354      -0.18592086
     355       1.31273975
     356      -0.90906728
     357      -0.09886636
     358       1.73726946
     359      -0.02479061
     360      -0.13041180
     361       0.03247863
     362       0.32341579
     363       0.01312362
     364       0.43285999
     365      -0.11288737
     366      -2.00439925
     367      -0.02653569
     368      -0.05466985
     369       0.05718408
     370       0.35401328
     371       0.01197739
     372       0.23096392
     373       0.03844139
     374      -0.23520748
     375      -0.01814868
     376       0.27163576
     377      -0.02791567
     378       0.09784728
     379      -0.05979336
     380       0.19133621
     381       1.55004626
     382      -1.07499834
     383       0.47256865
     384      -0.12187090
     385      -0.22333375
     386       0.01437911
     387      -0.01008033
     388      -0.05700712
     389       0.02339061
     390       0.09563192
     391      -0.02591122
     392      -0.03358746
     393      -0.10627013
     394      -0.26795475
     395      -0.00080816
     396      -0.05677330
     397       0.00994133
     398      -0.07976889
     399       0.11946882
     400      -0.06195141
     401       0.12418503
     402       0.01283373
     403       1.66064615
     404      -0.10175694
     405       0.06166033
     406       0.00663990
     407       0.01360799
     408      -0.13272020
     409       0.22254403
     410      -0.16648490
     411       0.10397025
     412      -0.08367881
     413       0.28518110
     414      -0.11100834
     415       0.09660099
     416      -0.04581243
     417       0.03314710
     418       1.04266448
     419      -0.06434735
     420       0.13399887
     421      -0.00449077
     422      -0.01833407
     423      -1.56373464
     424       0.02435843
     425      -0.00953717
     426       0.04095385
     427      -0.07986942
     428       0.28700293
     429      -0.03929706
     430       0.04597387
     431      -0.05750458
     432       0.12202557
     433       0.54175852
     434      -0.12525406
     435       0.01951587
     436      -0.01192237
     437      -0.25136939
     438       0.00502807
     439      -0.07101295
     440       0.04941002
     441       0.01835016
     442       0.01161796
     443      -0.11449899
     444      -0.15208107
     445       0.18540854
     446       0.11259108
     447       0.13969584
     448       0.02545202
     449       0.05178369
     450      -0.07582399
     451      -0.01458018
     452       0.04309011
     453      -0.35800866
     454      -0.03790401
     455       0.00860224
     456       0.05716678
     457       0.10441466
     458      -0.34977416
     459       0.09894641
     460       0.08720432
     461       0.06726787
     462       0.94243258
     463      -0.67618560
     464      -0.02650462
     465       0.09287407
     466      -0.08657497
     467       0.12414214
     468      -0.15859619
     469       0.01920391
     470       0.05152552
     471      -0.01812845
     472       1.17144198
     473      -0.62532002
     474       0.04822455
     475       0.19523654
     476      -0.11310872
     477       0.04120103
     478      -0.09794628
     479      -0.05286711
     480       0.11300864
     481       0.37964026
     482      -0.05650894
     483      -0.10461745
     484      -0.24264667
     485       0.03850066
     486       0.12441938
     487       0.01452194
     488       0.05621740
     489       0.26078868
     490      -0.09472305
     491       2.53248429
     492      -1.80724478

 Residual norm when dim(red L) =   6
 NEO root     CSF        orbital          total
    1     0.00000636     0.00000622     0.00000890 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  1

   1.988008899   1.976718325   0.023721229   0.012979077

 Symmetry  2

   1.980681197   0.018399249

 Symmetry  3

   1.975007038   0.024484985


 <<< MACRO ITERATION  5 >>>
 --------------------------

 Total MCSCF energy       :      -76.168770402655937       (MACRO    5)
 - Dielec. solvation ener.:       -0.012752890546972

 Norm of total gradient   :        0.000008880220
 -    of CI gradient      :        0.000006350862
 -    of orbital gradient :        0.000006206840

 *** Optimization control: MCSCF converged ***
     Number of macro iterations used            5
     Number of micro iterations used           26
     Total number of CPU seconds used         6.32


               >>> SIRIUS OPTIMIZATION STATISTICS <<<

 
     Date and time (Linux)  : Tue Mar 21 17:51:47 2006
     Host name              : star.chem.uit.no                        


  ITER ITMIC     EMCSCF           GRDNRM        RATIO      STPLNG
 ---------------------------------------------------------------------
    1    6    -76.165426189716   0.1463905479  0.000000   0.1983375306
    2    3    -76.168611452169   0.0574363430  0.931749   0.0250113864
    3   13    -76.168759959060   0.0104174300  1.005999   0.0060317296
    4    4    -76.168770402289   0.0000692491  0.998291   0.0000399989
    5    0    -76.168770402656   0.0000088802  1.006921   0.0000000000


  ITER  INDGCM  GCIMAX      GCINRM     INDGOM  GOBMAX      GOBNRM      GRDNRM
 ------------------------------------------------------------------------------
    1    372    0.027847    0.121252     22    0.040964    0.082025    0.146391
    2    275   -0.008216    0.031046     25   -0.029141    0.048323    0.057436
    3     51   -0.004088    0.010053     25   -0.001874    0.002733    0.010417
    4    134    0.000015    0.000035     24   -0.000040    0.000060    0.000069
    5    132    0.000002    0.000006     20    0.000002    0.000006    0.000009


  ITER ITMIC NCLIN NOLIN   TIMMAC    TIMITR    TIMMIC    TIMLIN    TIMMIC/ITMIC
 ------------------------------------------------------------------------------

    1     6     4     3      2.25      0.03      1.70      1.70      0.28
    2     3     2     2      0.74      0.01      0.61      0.61      0.20
    3    13     8     6      2.25      0.01      2.12      2.11      0.16
    4     4     3     2      0.90      0.01      0.77      0.76      0.19
    5     0     0     0      0.13      0.01      0.00      0.00


 ITER         EMY                 EACTIV              EMCSCF

    1    -61.277934507063    -23.932111866515    -76.165426189716
    2    -61.281056264255    -23.930007582927    -76.168611452169
    3    -61.280857198055    -23.930260801297    -76.168759959060
    4    -61.280818710556    -23.930203202699    -76.168770402289
    5    -61.280818742879    -23.930203294868    -76.168770402656


 ITER         DEPRED              DEACT               RATIO

    1      0.000000000000      0.000000000000      0.000000000000
    2     -0.003418584327     -0.003185262453      0.931748978049
    3     -0.000147621311     -0.000148506891      1.005999000688
    4     -0.000010461105     -0.000010443229      0.998291227801
    5     -0.000000000364     -0.000000000367      1.006921467423


 ITER    BETA           GAMMA             STPLNG              RTRUST

    1      0.20000000  1.00000000      0.198337530597      0.700000000000
    2      0.20000000  1.00000000      0.025011386362      0.700000000000
    3      0.20000000  1.00000000      0.006031729576      0.700000000000
    4      0.20000000  1.00000000      0.000039998864      0.700000000000
    5      0.00000000  0.00000000      0.000000000000      0.700000000000


 Reduced L root no.  1
 ITER         EVAL              EVEC(1)           EVEC(2)           EVEC(3)
 ----------------------------------------------------------------------------
    1   -0.000273057088    0.999214171740   -0.006671963450   -0.037613390345
    2   -0.000011809409    0.999987488846   -0.001052790139   -0.004505908002
    3   -0.000000836887    0.999999272366   -0.000602246686   -0.000860488664
    4   -0.000000000029    0.999999999968   -0.000001839441   -0.000007335519
    5    0.000000000000    0.000000000000    0.000000000000    0.000000000000


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000
     State number:                1

     Final MCSCF energy:          -76.168770402656
     Nuclear repulsion:             9.055004525638
     Electronic energy:           -85.211022037747
     Solvation  energy:            -0.012752890547

     Final gradient norm:           0.000008880220

 
     Date and time (Linux)  : Tue Mar 21 17:51:47 2006
     Host name              : star.chem.uit.no                        


 Occupancies of natural orbitals
 -------------------------------

 Symmetry  1

   2.000000000   1.988008899   1.976718325   0.023721229   0.012979077

 Sum =           6.001427530

 Symmetry  2

   1.980681197   0.018399249

 Sum =           1.999080446

 Symmetry  3

   1.975007038   0.024484985

 Sum =           1.999492024

 Symmetry  4

   No occupied orbitals

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5
   1  O   1s     0.5512  -0.1192  -0.0150   0.0545  -0.1324
   2  O   1s     0.4714  -0.1988  -0.0243   0.0948  -0.2028
   3  O   1s     0.0066   0.5665   0.0612  -0.4599   1.4085
   4  O   1s     0.0002   0.5214  -0.0514  -0.0752  -1.3555
   5  O   2pz    0.0018  -0.0750   0.2501  -0.4536  -0.2385
   6  O   2pz   -0.0003  -0.0859   0.3782  -0.5020  -0.1784
   7  O   2pz    0.0003  -0.1260   0.3763   0.4072   0.4739
   8  O   3d0    0.0000  -0.0045   0.0200   0.0085   0.0184
   9  O   3d2+   0.0001  -0.0101  -0.0135   0.0169   0.0210
  10  H   1s     0.0004   0.0485   0.1199   0.2020  -0.0938
  11  H   1s     0.0001  -0.0104   0.1096   0.2623   0.0504
  12  H   1s     0.0001  -0.0093   0.0212   0.0337   0.0507

     Molecular orbitals for symmetry species   2

 Orbital          1        2
   1  O   2px    0.2945   0.5493
   2  O   2px    0.4234   0.4802
   3  O   2px    0.5050  -0.9662
   4  O   3d1+   0.0189  -0.0253

     Molecular orbitals for symmetry species   3

 Orbital          1        2
   1  O   2py    0.2307   0.4762
   2  O   2py    0.3709   0.5522
   3  O   2py    0.2526  -0.2065
   4  O   3d1-   0.0391   0.0355
   5  H   1s     0.1501  -0.2121
   6  H   1s     0.1716  -0.3316
   7  H   1s     0.0545  -0.0742


 Printout of CI-coefficients larger than 0.05000 for root  1

 (this is the reference state)



  Printout of coefficients in interval  0.3162     to 1.0000    
  =============================================================

 Coefficient of CSF no.        20 is  0.98040186    
 Occupation and spin coupling 
        1    2    5    7
        2    2    2    2


  Printout of coefficients in interval  1.0000E-01 to 0.3162    
  =============================================================
   ( no coefficients )


  Printout of coefficients in interval  5.0000E-02 to 1.0000E-01
  =============================================================

 Coefficient of CSF no.        21 is -6.04641303E-02
 Occupation and spin coupling 
        1    3    5    7
        2    2    2    2

 Coefficient of CSF no.        26 is -5.32506104E-02
 Occupation and spin coupling 
        1    2    6    7
        2    2    2    2

 Coefficient of CSF no.        40 is -6.46889454E-02
 Occupation and spin coupling 
        1    2    5    8
        2    2    2    2

 Coefficient of CSF no.       236 is -6.19816676E-02
 Occupation and spin coupling 
        1    2    5    6    7    8
        2    2    1   -1    1   -1

 Coefficient of CSF no.       268 is  6.85626056E-02
 Occupation and spin coupling 
        1    5    2    3    7    8
        2    2    1   -1    1   -1

 Coefficient of CSF no.       330 is  6.08138104E-02
 Occupation and spin coupling 
        1    7    2    3    5    6
        2    2    1   -1    1   -1


  Norm of printed CI vector ..      0.98410488

   Magnitude of CI coefficients 
  ==============================

   ( Ranges are relative to norm of vector :  1.00E+00 )

  10- 1 to 10- 0         1    0.96118780E+00    0.96118780E+00
  10- 2 to 10- 1        30    0.37244350E-01    0.99843215E+00
  10- 3 to 10- 2        95    0.15237586E-02    0.99995591E+00
  10- 4 to 10- 3       222    0.43795384E-04    0.99999970E+00
  10- 5 to 10- 4       112    0.29517380E-06    0.10000000E+01
  10- 6 to 10- 5        26    0.10092645E-08    0.10000000E+01
  10- 7 to 10- 6         6    0.16297491E-11    0.10000000E+01
  Number of coefficients less than 10^-11 times norm is  0



 >>>> Total CPU  time used in SIRIUS :      6.59 seconds
 >>>> Total wall time used in SIRIUS :      7.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 17:51:47 2006
     Host name              : star.chem.uit.no                        

- End of Wave Function Section



    *****************************************************************
    ******** Output from **PROPE input processing for ABACUS ********
    *****************************************************************



 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

      Natural orbital connection is used
      for perturbation dependent basis sets.


 Default print level:        0


Starting in Static Property Section -


 
     Date and time (Linux)  : Tue Mar 21 17:51:47 2006
     Host name              : star.chem.uit.no                        


 ***************************************************************************
 ************************ FINAL RESULTS FROM ABACUS ************************
 ***************************************************************************


 
     Date and time (Linux)  : Tue Mar 21 17:51:47 2006
     Host name              : star.chem.uit.no                        



                         Molecular geometry (au)
                         -----------------------

 O          0.0000000000            0.0000000000           -0.1258515023
 H    1     0.0000000000            1.4523500000            0.9986773907
 H    2     0.0000000000           -1.4523500000            0.9986773907





                    Molecular wave function and energy
                    ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0


     Total energy        -76.1687704027 au (Hartrees)
                         -2072.65761741 eV
                           -199981.0780 kJ/mol




                         Relativistic corrections
                         ------------------------

     Darwin correction:                          0.1973259369 au
     Mass-velocity correction:                  -0.2489280560 au

     Total relativistic correction:             -0.0516021191 au (0.0677%)
     Non-relativistic + relativistic energy:   -76.2203725218 au




                              Dipole moment
                              -------------


                    1.034668 au           2.629862 Debye




                         Dipole moment components
                         ------------------------

                               au             Debye

                    z      1.03466754      2.62986228

                        1 a.u. =   2.54175 Debye 





   Interatomic separations (in Angstroms):
   ---------------------------------------

            O           H    1      H    2
            ------      ------      ------
 O     :    0.000000
 H    1:    0.972000    0.000000
 H    2:    0.972000    1.537101    0.000000


  Max interatomic separation is    1.5371 Angstroms
  between atoms "H    2" and "H    1".


  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  H    1     O            0.972000
  bond distance:  H    2     O            0.972000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     H    1     O          H    2       104.500




 CPU time statistics for ABACUS
 ------------------------------





 >>>> Total CPU  time used in ABACUS:   0.01 seconds
 >>>> Total wall time used in ABACUS:   0.00 seconds

- End of Static Property Section

 >>>> Total CPU  time used in DALTON:   6.78 seconds
 >>>> Total wall time used in DALTON:   8.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 17:51:47 2006
     Host name              : star.chem.uit.no                        
END REFOUT

