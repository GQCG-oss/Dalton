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
**RESPONSE
*QUADRATIC
.SINGLE RESIDUE
.MCDBTERM
.ROOTS
 4 4 4 4
*END OF INPUT
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

     Date and time (Linux)  : Tue Mar 21 14:59:41 2006
     Host name              : star.chem.uit.no                        

 <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>


 Default print level:        0

    Integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed
    Dynamic molecular property section will be executed

     -----------------------------------
     INPUT FOR PCM SOLVATION CALCULATION 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=INPUT        EPS   = 24.5500     EPSINF=  0.0000
     RSOLV =  2.1800

     ICESPH =       0     NESFP =       0
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 

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
     "/home/lara/programs/main-branch/dalton/basis/aug-cc-pVDZ"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/lara/programs/main-branch/dalton/basis/aug-cc-pVDZ"


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

  Basis set used is "aug-cc-pVDZ" from the basis set library.

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  O           1    8.0000    35    23      [10s5p2d|4s3p2d]                                   
  H           2    1.0000    11     9      [5s2p|3s2p]                                        
  ----------------------------------------------------------------------
  total:      3   10.0000    57    41
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


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.00 seconds


 >>> Time used in ONEDRV is   0.00 seconds


 Number of two-electron integrals written:       91307 ( 24.6% )
 Megabytes written:                              1.051



 >>> Time used in TWOINT is   0.08 seconds


 MEMORY USED TO GENERATE CAVITY=    432042


 TOTAL NUMBER OF SPHERES=    3
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.000000000    0.000000000   -0.066597747    1.800000000   24.581233228
   2    0.000000000    0.768550518    0.528477314    1.440000000   11.987539425
   3    0.000000000   -0.768550518    0.528477314    1.440000000   11.987539425

 TOTAL NUMBER OF TESSERAE =     152
 SURFACE AREA=   48.55631208(A**2)    CAVITY VOLUME=   30.53397368 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.11 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.01 seconds

 >>>> Total CPU  time used in HERMIT:   0.12 seconds
 >>>> Total wall time used in HERMIT:   0.00 seconds

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

 
     Date and time (Linux)  : Tue Mar 21 14:59:41 2006
     Host name              : star.chem.uit.no                        

 Title lines from integral program:
     magnet                                                                  
     nuova base                                                              

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham calculation.


     Time-dependent Kohn-Sham calculation (random phase approximation).


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         10
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Number of active orbitals                 0
     Total number of orbitals                 41

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
     Abelian symmetry species           1   2   3   4
                                       --  --  --  --
     Total number of orbitals          18   7  12   4
     Number of basis functions         18   7  12   4

      ** Automatic occupation of RKS orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals              3   1   1   0

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00E-06

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE        9.97104
 NUCLEAR APPARENT CHARGE  -9.59233 THEORETICAL  -9.59267 NOT RENORMALIZED
 this is icompcm in icvev 0

  ..... DONE WITH INDUCED NUCLEAR CHARGES .....

 >>> Time used in VNN    is   0.00 seconds



 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  10 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Atom:    1*1 points=24126 compressed from 24126 (123 radial)
 Atom:    2*2 points=23192 compressed from 23238 ( 93 radial)
 Number of grid points:    47318 Grid generation time:       0.1 s
K-S electrons/energy :    9.99999963893959   -7.40267067160745 err:-.36E-06
   1    -76.096099516245     -0.010001762589   3.47E+00  -7.61E+01    3  1  1  0
K-S electrons/energy :    9.99999967099836   -7.61229087326312 err:-.33E-06
   2    -76.385634120506     -0.008037525253   8.09E-01  -2.90E-01    3  1  1  0
K-S electrons/energy :    9.99999967941059   -7.35507645041532 err:-.32E-06
   3    -76.386822537107     -0.006961402432   7.94E-01  -1.19E-03    3  1  1  0
K-S electrons/energy :    9.99999967108957   -7.52243679054402 err:-.33E-06
   4    -76.414055717212     -0.009834812663   1.67E-01  -2.72E-02    3  1  1  0
K-S electrons/energy :    9.99999967340932   -7.49475485036522 err:-.33E-06
   5    -76.415354379149     -0.008895700480   8.31E-03  -1.30E-03    3  1  1  0
K-S electrons/energy :    9.99999967334498   -7.49578259379331 err:-.33E-06
   6    -76.415357501193     -0.008935105003   7.31E-04  -3.12E-06    3  1  1  0
K-S electrons/energy :    9.99999967333835   -7.49592138006900 err:-.33E-06
   7    -76.415357528461     -0.008936982218   3.85E-05  -2.73E-08    3  1  1  0
K-S electrons/energy :    9.99999967333863   -7.49591471041138 err:-.33E-06
   8    -76.415357528566     -0.008936834419   1.75E-06  -1.05E-10    3  1  1  0
K-S electrons/energy :    9.99999967333864   -7.49591457765417 err:-.33E-06
   9    -76.415357528569     -0.008936814287   1.59E-07  -2.43E-12    3  1  1  0
 DIIS converged in   9 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   10
 Orbital occupations :    3    1    1    0

 Kohn-Sham orbital energies, symmetry 1

       -19.15770742    -1.01435582    -0.39866831    -0.01408857     0.10168960
         0.12975013     0.22734410     0.29590864

 Kohn-Sham orbital energies, symmetry 2

        -0.32191472     0.11216628     0.32696197     0.89445356     1.18997715
         1.84727002

 Kohn-Sham orbital energies, symmetry 3

        -0.52976372     0.03359655     0.16210902     0.17626516     0.45234851
         0.53245658

 Kohn-Sham orbital energies, symmetry 4

         0.28746054     0.95198394     1.75179832     3.31154184

    E(LUMO) :    -0.01408857 au (symmetry 1)
  - E(HOMO) :    -0.32191472 au (symmetry 2)
  ------------------------------------------
    gap     :     0.30782615 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   24.550000

     Final DFT energy:            -76.415357528569
     Nuclear repulsion:             9.055004525638
     Electronic energy:           -85.461425239920

     Final gradient norm:           0.000000158934

 
     Date and time (Linux)  : Tue Mar 21 15:00:08 2006
     Host name              : star.chem.uit.no                        

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3
   1  O   1s     1.0039  -0.0038   0.0030
   2  O   1s     0.0196   0.8880  -0.3045
   3  O   1s    -0.0259  -0.1459  -0.1428
   4  O   1s    -0.0101   0.0662  -0.1249
   5  O   2pz    0.0024   0.1543   0.8088
   6  O   2pz   -0.0041  -0.0863  -0.0713
   7  O   2pz   -0.0016   0.0136   0.0562
   8  O   3d0    0.0000   0.0034   0.0143
   9  O   3d2+  -0.0003  -0.0075  -0.0025
  10  O   3d0    0.0001  -0.0060   0.0105
  11  O   3d2+   0.0007   0.0198   0.0149
  12  H   1s     0.0014   0.3561   0.3556
  13  H   1s     0.0075  -0.2059  -0.1080
  14  H   1s     0.0007  -0.0022   0.0076
  15  H   2py    0.0004  -0.0362  -0.0262
  16  H   2pz    0.0004  -0.0209   0.0036
  17  H   2py   -0.0027   0.0041  -0.0034
  18  H   2pz   -0.0024   0.0069   0.0138

     Molecular orbitals for symmetry species   2

 Orbital          1
   1  O   2px    0.9200
   2  O   2px   -0.0003
   3  O   2px    0.1036
   4  O   3d1+   0.0137
   5  O   3d1+   0.0140
   6  H   2px    0.0261
   7  H   2px    0.0211

     Molecular orbitals for symmetry species   3

 Orbital          1
   1  O   2py    0.7556
   2  O   2py   -0.1959
   3  O   2py    0.0153
   4  O   3d1-   0.0249
   5  O   3d1-  -0.0461
   6  H   1s     0.6088
   7  H   1s    -0.1713
   8  H   1s     0.0018
   9  H   2py   -0.0214
  10  H   2pz   -0.0322
  11  H   2py   -0.0232
  12  H   2pz   -0.0159



 >>>> Total CPU  time used in SIRIUS :     12.77 seconds
 >>>> Total wall time used in SIRIUS :     27.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 15:00:08 2006
     Host name              : star.chem.uit.no                        

- End of Wave Function Section



  This is output from RESPONSE  -  an MCSCF and SOPPA response property program
 ------------------------------------------------------------------------------



 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




 Quadratic Response single residue calculation
 ---------------------------------------------


 Spin of operator A , ISPINA=    0
 Spin of operator B , ISPINB=    0
 Spin of operator C , (Excitation energy) ISPINC=    0

  1 B-frequencies  0.000000E+00
 B of Magnetic Circular Dichroism requested    : MCDCAL =T

 Print level                                    : IPRSMO =    2
 Threshold for convergence in RSPPP             : THCPP  = 1.000E-03
 Maximum number of iterations in RSPPP          : MAXITP =   60
 Threshold for convergence in RSPLR             : THCLR  = 1.000E-03
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


   SCF energy         :      -76.415357528568862
 -- inactive part     :      -85.461425239920104
 -- nuclear repulsion :        9.055004525637894


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

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.2 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.2 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*3 evaluation time:       2.6 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   30)
 RSP solution vector no.    1; norm of residual   2.56E-04
 RSP solution vector no.    2; norm of residual   3.94E-04
 RSP solution vector no.    3; norm of residual   4.10E-04
 RSP solution vector no.    4; norm of residual   1.87E-04

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

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   24)
 RSP solution vector no.    1; norm of residual   1.21E-04
 RSP solution vector no.    2; norm of residual   3.61E-04
 RSP solution vector no.    3; norm of residual   6.31E-04
 RSP solution vector no.    4; norm of residual   7.95E-04

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

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.2 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.2 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.2 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*3 evaluation time:       2.5 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   30)
 RSP solution vector no.    1; norm of residual   9.34E-05
 RSP solution vector no.    2; norm of residual   5.75E-05
 RSP solution vector no.    3; norm of residual   9.52E-05
 RSP solution vector no.    4; norm of residual   1.98E-04

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

 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*4 evaluation time:       3.1 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.2 s

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   26)
 RSP solution vector no.    1; norm of residual   1.89E-05
 RSP solution vector no.    2; norm of residual   3.28E-05
 RSP solution vector no.    3; norm of residual   4.14E-05
 RSP solution vector no.    4; norm of residual   1.08E-04

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
 QRLRVE -- frequencies :  0.260094  0.371728  0.408196  0.458800  0.394814
                          0.468683  0.512613  0.522952  0.316502  0.427812
                          0.443797  0.591276



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.6 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.37173 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     -3.53553796
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -6.11 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.40820 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      5.11860562
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  7.11 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.45880 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     18.43474531
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  3.00 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.39481 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      4.50455533
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.05 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.46868 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     27.66421209
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  7.98 * 10 ** -3.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.51261 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      2.45732151
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   20    0    4
 DETERMINANT OF REDUCED MATRIX IS  1.73 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   20    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.52295 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      6.62658122
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   20    0    4
 DETERMINANT OF REDUCED MATRIX IS  1.52 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   20    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      9.27846231
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.32 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     12.77956192
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  8.38 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     -9.17992109
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   18    0    6
 DETERMINANT OF REDUCED MATRIX IS  1.68 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   18    0    6
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.6 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.37173 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     -3.14693937
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.79 * 10 **  9.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.40820 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      6.04561007
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  2.72 * 10 **  8.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.45880 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     19.26482564
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  2.43 * 10 **  7.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.39481 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      5.79831121
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.00 * 10 **  8.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.46868 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS   -201.21156101
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   45    0    3
 DETERMINANT OF REDUCED MATRIX IS -1.40 * 10 **  5.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   45    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.51261 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      5.73792867
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   44    0    4
 DETERMINANT OF REDUCED MATRIX IS  1.37 * 10 **  7.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   44    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.52295 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      9.47745215
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   44    0    4
 DETERMINANT OF REDUCED MATRIX IS  8.50 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   44    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      9.69216431
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  2.68 * 10 **  8.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     13.11925669
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.19 * 10 **  8.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     14.81236038
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   42    0    6
 DETERMINANT OF REDUCED MATRIX IS  9.89 * 10 **  5.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   42    0    6
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.6 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.37173 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     -3.14664313
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.56 * 10 ** 23.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.40820 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS      6.04688672
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   70    0    2
 DETERMINANT OF REDUCED MATRIX IS  2.03 * 10 ** 22.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   70    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.45880 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     19.26600152
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   70    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.42 * 10 ** 21.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   70    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.39481 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS      5.80595048
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -7.89 * 10 ** 21.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.46868 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS   -170.34957422
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   69    0    3
 DETERMINANT OF REDUCED MATRIX IS -9.05 * 10 ** 18.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   69    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.51261 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS      5.74187695
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   68    0    4
 DETERMINANT OF REDUCED MATRIX IS  5.98 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   68    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.52295 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS      9.48144579
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   68    0    4
 DETERMINANT OF REDUCED MATRIX IS  3.47 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   68    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS      9.69269668
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   70    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.83 * 10 ** 22.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   70    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     13.11953837
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   70    0    2
 DETERMINANT OF REDUCED MATRIX IS  7.56 * 10 ** 21.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   70    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     15.07398995
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   66    0    6
 DETERMINANT OF REDUCED MATRIX IS  2.45 * 10 ** 19.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   66    0    6

 *** THE REQUESTED   12 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   72)
 RSP solution vector no.    1; norm of residual   2.53E-04
 RSP solution vector no.    2; norm of residual   3.54E-05
 RSP solution vector no.    3; norm of residual   1.87E-04
 RSP solution vector no.    4; norm of residual   6.19E-05
 RSP solution vector no.    5; norm of residual   1.54E-04
 RSP solution vector no.    6; norm of residual   2.86E-04
 RSP solution vector no.    7; norm of residual   2.22E-04
 RSP solution vector no.    8; norm of residual   2.47E-04
 RSP solution vector no.    9; norm of residual   1.33E-04
 RSP solution vector no.   10; norm of residual   1.41E-04
 RSP solution vector no.   11; norm of residual   8.76E-05
 RSP solution vector no.   12; norm of residual   3.80E-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.260094E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.26009):     17.0603627435    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.371728E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.37173):    -3.14664312898    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.408196E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.40820):     6.04688671942    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.458800E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.45880):     19.2660015160    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.394814E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.39481):     5.80595048289    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.468683E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.46868):    -170.349574217    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.512613E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.51261):     5.74187695402    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.522952E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.52295):     9.48144578833    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.316502E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.31650):     37.5548918724    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.427812E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.42781):     9.69269668321    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.443797E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.44380):     13.1195383673    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.591276E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.59128):     15.0739899459    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    2


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      37
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      37


 QRLRVE -- linear response calculation for symmetry  2
 QRLRVE -- operator label : XDIPLEN 
 QRLRVE -- frequencies :  0.397228  0.468548  0.483278  0.394814  0.468683
                          0.512613  0.522952  0.316502  0.427812  0.443797
                          0.591276  0.336845



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.2 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.39723 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     12.18015238
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.68 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.46855 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     22.87364759
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   20    0    4
 DETERMINANT OF REDUCED MATRIX IS  8.97 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   20    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.48328 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     60.08605607
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   20    0    4
 DETERMINANT OF REDUCED MATRIX IS  9.34 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   20    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.39481 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     11.62900059
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.94 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.46868 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     23.02762271
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   20    0    4
 DETERMINANT OF REDUCED MATRIX IS  9.06 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   20    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.51261 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     -5.36452313
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   19    0    5
 DETERMINANT OF REDUCED MATRIX IS -6.38 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   19    0    5

 GP * SOLUTION VECTOR AT FREQUENCY     0.52295 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      2.32026460
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   19    0    5
 DETERMINANT OF REDUCED MATRIX IS -9.82 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   19    0    5

 GP * SOLUTION VECTOR AT FREQUENCY     0.31650 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      2.91355002
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -6.84 * 10 **  0.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     12.60474394
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   21    0    3
 DETERMINANT OF REDUCED MATRIX IS -2.02 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   21    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     16.72698051
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   21    0    3
 DETERMINANT OF REDUCED MATRIX IS -1.64 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   21    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     -0.92803021
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   18    0    6
 DETERMINANT OF REDUCED MATRIX IS  9.13 * 10 **  0.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   18    0    6

 GP * SOLUTION VECTOR AT FREQUENCY     0.33685 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      5.50725112
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.38 * 10 **  0.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.2 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.39723 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     12.26267345
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  8.78 * 10 ** 11.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.46855 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     24.03924798
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   44    0    4
 DETERMINANT OF REDUCED MATRIX IS  3.47 * 10 ** 11.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   44    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.48328 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     64.79450520
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   44    0    4
 DETERMINANT OF REDUCED MATRIX IS  3.08 * 10 ** 11.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   44    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.39481 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     11.71762350
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.04 * 10 ** 12.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.46868 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     24.19537177
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   44    0    4
 DETERMINANT OF REDUCED MATRIX IS  3.50 * 10 ** 11.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   44    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.51261 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     -5.14130113
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   43    0    5
 DETERMINANT OF REDUCED MATRIX IS -1.83 * 10 ** 12.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   43    0    5

 GP * SOLUTION VECTOR AT FREQUENCY     0.52295 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      2.55928230
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   43    0    5
 DETERMINANT OF REDUCED MATRIX IS -2.59 * 10 ** 12.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   43    0    5

 GP * SOLUTION VECTOR AT FREQUENCY     0.31650 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      3.12675733
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -5.21 * 10 ** 13.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     13.06893151
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   45    0    3
 DETERMINANT OF REDUCED MATRIX IS -9.72 * 10 ** 11.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   45    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     17.18477034
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   45    0    3
 DETERMINANT OF REDUCED MATRIX IS -7.00 * 10 ** 11.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   45    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     -0.73130679
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   42    0    6
 DETERMINANT OF REDUCED MATRIX IS  1.39 * 10 ** 13.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   42    0    6

 GP * SOLUTION VECTOR AT FREQUENCY     0.33685 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      5.67950842
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.36 * 10 ** 13.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1
 Electrons: 10.000000(-3.27e-07): LR-DFT*5 evaluation time:       3.7 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.39723 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     12.26267956
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   56    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.37 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   56    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.46855 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     24.03931834
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   54    0    4
 DETERMINANT OF REDUCED MATRIX IS  5.16 * 10 ** 19.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   54    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.48328 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     64.79474784
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   54    0    4
 DETERMINANT OF REDUCED MATRIX IS  4.53 * 10 ** 19.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   54    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.39481 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     11.71762815
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   56    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.63 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   56    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.46868 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     24.19544125
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   54    0    4
 DETERMINANT OF REDUCED MATRIX IS  5.21 * 10 ** 19.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   54    0    4

 GP * SOLUTION VECTOR AT FREQUENCY     0.51261 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     -5.14124298
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   53    0    5
 DETERMINANT OF REDUCED MATRIX IS -2.64 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   53    0    5

 GP * SOLUTION VECTOR AT FREQUENCY     0.52295 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS      2.55932646
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   53    0    5
 DETERMINANT OF REDUCED MATRIX IS -3.69 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   53    0    5

 GP * SOLUTION VECTOR AT FREQUENCY     0.31650 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS      3.12676978
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   57    0    1
 DETERMINANT OF REDUCED MATRIX IS -8.48 * 10 ** 21.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   57    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     13.06895298
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   55    0    3
 DETERMINANT OF REDUCED MATRIX IS -1.49 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   55    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     17.18479762
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   55    0    3
 DETERMINANT OF REDUCED MATRIX IS -1.06 * 10 ** 20.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   55    0    3

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS     -0.73121283
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   52    0    6
 DETERMINANT OF REDUCED MATRIX IS  1.86 * 10 ** 21.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   52    0    6

 GP * SOLUTION VECTOR AT FREQUENCY     0.33685 AU
 AFTER   29 LINEAR TRANSFORMATIONS IS      5.67951709
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   57    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.81 * 10 ** 21.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   57    0    1

 *** THE REQUESTED   12 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   58)
 RSP solution vector no.    1; norm of residual   2.51E-05
 RSP solution vector no.    2; norm of residual   1.03E-05
 RSP solution vector no.    3; norm of residual   1.08E-05
 RSP solution vector no.    4; norm of residual   2.64E-05
 RSP solution vector no.    5; norm of residual   1.03E-05
 RSP solution vector no.    6; norm of residual   2.19E-05
 RSP solution vector no.    7; norm of residual   3.19E-05
 RSP solution vector no.    8; norm of residual   5.74E-05
 RSP solution vector no.    9; norm of residual   1.21E-05
 RSP solution vector no.   10; norm of residual   1.31E-05
 RSP solution vector no.   11; norm of residual   1.15E-05
 RSP solution vector no.   12; norm of residual   4.96E-05

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.397228E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.39723):     12.2626795623    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.39723):    0.434462817481    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.468548E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.46855):     24.0393183440    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.46855):    -3.38537732037    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.483278E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.48328):     64.7947478394    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.48328):    -10.6058751189    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.394814E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.39481):     11.7176281540    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.39481):    0.212453393198    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.468683E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.46868):     24.1954412541    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.46868):    -3.40754258563    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.512613E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.51261):    -5.14124297637    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.51261):     3.26376319031    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.522952E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.52295):     2.55932645874    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.52295):     2.13930213818    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.316502E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.31650):     3.12676978315    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.31650):   -0.914041384220    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.427812E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.42781):     13.0689529833    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.42781):    -1.89613573132    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.443797E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.44380):     17.1847976156    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.44380):    -1.97372939905    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.591276E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.59128):   -0.731212829329    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.59128):     2.57407723861    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.336845E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.33685):     5.67951708771    

@QRLRVE:  << YANGMOM  ; XDIPLEN  >> (   0.33685):   -0.683357176141    


 QRLRVE -- linear response calculation for symmetry  2
 QRLRVE -- operator label : YANGMOM 
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   5.03E-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   YANGMOM     FREQUENCY   0.000000E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; YANGMOM  >> (   0.00000):   -6.450537058168E-17

@QRLRVE:  << YANGMOM  ; YANGMOM  >> (   0.00000):    0.563437652228    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    3


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      52
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      52


 QRLRVE -- linear response calculation for symmetry  3
 QRLRVE -- operator label : YDIPLEN 
 QRLRVE -- frequencies :  0.408196  0.458800  0.316502  0.427812  0.443797
                          0.591276  0.336845  0.397228  0.468548  0.483278
                          0.260094  0.371728



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.6 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.40820 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     21.75077814
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.58 * 10 **  0.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.45880 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     76.17942497
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.46 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     30.36266605
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.31 * 10 **  0.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     41.34564315
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -5.98 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     -9.65628345
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   18    0    6
 DETERMINANT OF REDUCED MATRIX IS  2.50 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   18    0    6

 GP * SOLUTION VECTOR AT FREQUENCY     0.39723 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     -2.35081229
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -4.29 * 10 ** -1.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.46855 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS   1001.51250684
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -5.09 * 10 ** -3.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.48328 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS    -10.81240123
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  4.85 * 10 ** -2.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.6 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.40820 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     22.31201495
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.45 * 10 ** 10.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.45880 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     79.48809378
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -7.96 * 10 **  8.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     31.08837645
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.00 * 10 ** 10.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     42.44844854
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.90 * 10 **  9.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     -7.37000558
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   42    0    6
 DETERMINANT OF REDUCED MATRIX IS  3.80 * 10 **  7.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   42    0    6

 GP * SOLUTION VECTOR AT FREQUENCY     0.39723 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS      0.71269276
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -4.75 * 10 **  9.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.46855 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS   4167.76819611
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   47    0    1
 DETERMINANT OF REDUCED MATRIX IS -6.19 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   47    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.48328 AU
 AFTER   24 LINEAR TRANSFORMATIONS IS     -8.02540874
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   46    0    2
 DETERMINANT OF REDUCED MATRIX IS  2.17 * 10 **  8.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   46    0    2
 Electrons: 10.000000(-3.27e-07): LR-DFT*12 evaluation time:       9.6 s

 GP * SOLUTION VECTOR AT FREQUENCY     0.40820 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     22.31247443
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.04 * 10 ** 26.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.45880 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     79.49019728
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -4.67 * 10 ** 24.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.42781 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     31.08884026
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -6.67 * 10 ** 25.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.44380 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     42.44911034
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.44 * 10 ** 25.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.59128 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     -7.36815599
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   66    0    6
 DETERMINANT OF REDUCED MATRIX IS  1.10 * 10 ** 23.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   66    0    6

 GP * SOLUTION VECTOR AT FREQUENCY     0.39723 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS      0.71546512
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.54 * 10 ** 25.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.46855 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS   4177.55620681
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   71    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.47 * 10 ** 22.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   71    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.48328 AU
 AFTER   36 LINEAR TRANSFORMATIONS IS     -8.02321156
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   70    0    2
 DETERMINANT OF REDUCED MATRIX IS  1.14 * 10 ** 24.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   70    0    2

 *** THE REQUESTED   12 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   72)
 RSP solution vector no.    1; norm of residual   2.34E-05
 RSP solution vector no.    2; norm of residual   1.78E-05
 RSP solution vector no.    3; norm of residual   5.21E-05
 RSP solution vector no.    4; norm of residual   2.98E-05
 RSP solution vector no.    5; norm of residual   2.49E-05
 RSP solution vector no.    6; norm of residual   1.78E-05
 RSP solution vector no.    7; norm of residual   4.87E-05
 RSP solution vector no.    8; norm of residual   8.65E-06
 RSP solution vector no.    9; norm of residual   1.54E-05
 RSP solution vector no.   10; norm of residual   2.57E-05
 RSP solution vector no.   11; norm of residual   5.86E-05
 RSP solution vector no.   12; norm of residual   3.69E-05

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.408196E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.40820):     22.3124744347    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.40820):   -0.340691893908    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.458800E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.45880):     79.4901972804    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.45880):    -3.76741610160    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.316502E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.31650):     17.4395669121    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.31650):    0.201007414075    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.427812E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.42781):     31.0888402567    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.42781):   -0.120822217305    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.443797E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.44380):     42.4491103376    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.44380):   -0.573390204729    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.591276E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.59128):    -7.36815599320    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.59128):    -8.96664354426    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.336845E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.33685):     18.8827847568    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.33685):    0.263552291487    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.397228E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.39723):    0.715465117225    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.39723):    -2.91008849535    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.468548E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.46855):     4177.55620681    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.46855):    -473.506310992    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.483278E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.48328):    -8.02321156490    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.48328):     9.57756626019    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.260094E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.26009):     14.9197134456    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.26009):    0.105489360160    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.371728E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.37173):     23.3127093322    

@QRLRVE:  << XANGMOM  ; YDIPLEN  >> (   0.37173):    0.520312096392    


 QRLRVE -- linear response calculation for symmetry  3
 QRLRVE -- operator label : XANGMOM 
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F

 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   1.74E-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   XANGMOM     FREQUENCY   0.000000E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; XANGMOM  >> (   0.00000):    1.239935532585E-16

@QRLRVE:  << XANGMOM  ; XANGMOM  >> (   0.00000):     2.07779614037    


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

 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s
 Electrons: 10.000000(-3.27e-07): LR-DFT*1 evaluation time:       1.3 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   1.44E-04

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZANGMOM     FREQUENCY   0.000000E+00
 SYMMETRY    4

@QRLRVE:  << ZANGMOM  ; ZANGMOM  >> (   0.00000):     1.14494107236    
 DFT-QR computed in a linearly-scaling fashion.

 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.000000    0.336845   -8.381613

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0
 Excitation energy in au,    moment :                0.336845   -0.737144

 B term contribution:              3.089228
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.000000    0.397228   -7.166795

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0
 Excitation energy in au,    moment :                0.397228   -0.067944

 B term contribution:              0.243469
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.000000    0.468548    9.871271

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0
 Excitation energy in au,    moment :                0.468548    0.162247

 B term contribution:              0.800794
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             4    1    0

 omega B, excitation energy, moment :    0.000000    0.483278    6.976227

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             4    1    0
 Excitation energy in au,    moment :                0.483278    0.457883

 B term contribution:              1.597149
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.000000    0.336845   -5.468479

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0
 Excitation energy in au,    moment :                0.336845   -0.737144

 B term contribution:              2.015529
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.000000    0.397228   10.999923

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0
 Excitation energy in au,    moment :                0.397228   -0.067944

 B term contribution:             -0.373687
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.000000    0.468548 -449.691405

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0
 Excitation energy in au,    moment :                0.468548    0.162247

 B term contribution:            -36.480622
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             4    1    0

 omega B, excitation energy, moment :    0.000000    0.483278    6.361957

 First order moment in a.u. for
 C operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             4    1    0
 Excitation energy in au,    moment :                0.483278    0.457883

 B term contribution:              1.456517
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.000000    0.260094   -9.238086

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0
 Excitation energy in au,    moment :                0.260094   -0.651687

 B term contribution:              3.010172
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.000000    0.371728    1.045784

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0
 Excitation energy in au,    moment :                0.371728   -0.002162

 B term contribution:             -0.001131
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.000000    0.408196   -6.718618

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0
 Excitation energy in au,    moment :                0.408196   -0.141087

 B term contribution:              0.473953
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             4    2    0

 omega B, excitation energy, moment :    0.000000    0.458800  -20.554723

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             4    2    0
 Excitation energy in au,    moment :                0.458800   -0.157023

 B term contribution:              1.613785
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.000000    0.260094    2.756296

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0
 Excitation energy in au,    moment :                0.260094   -0.651687

 B term contribution:             -0.898122
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.000000    0.371728   -1.412016

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0
 Excitation energy in au,    moment :                0.371728   -0.002162

 B term contribution:              0.001526
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.000000    0.408196    0.155424

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0
 Excitation energy in au,    moment :                0.408196   -0.141087

 B term contribution:             -0.010964
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             4    2    0

 omega B, excitation energy, moment :    0.000000    0.458800    7.028364

 First order moment in a.u. for
 C operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             4    2    0
 Excitation energy in au,    moment :                0.458800   -0.157023

 B term contribution:             -0.551808
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.000000    0.394814   -3.141205

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0
 Excitation energy in au,    moment :                0.394814   -0.238923

 B term contribution:              0.375252
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.000000    0.468683   92.103182

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0
 Excitation energy in au,    moment :                0.468683   -0.748664

 B term contribution:            -34.477157
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.000000    0.512613    6.488797

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0
 Excitation energy in au,    moment :                0.512613   -0.203752

 B term contribution:             -0.661053
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             4    3    0

 omega B, excitation energy, moment :    0.000000    0.522952  -10.906978

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             4    3    0
 Excitation energy in au,    moment :                0.522952    0.584102

 B term contribution:             -3.185396
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.000000    0.394814    1.661200

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0
 Excitation energy in au,    moment :                0.394814   -0.238923

 B term contribution:             -0.198449
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.000000    0.468683    2.145468

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0
 Excitation energy in au,    moment :                0.468683   -0.748664

 B term contribution:             -0.803117
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.000000    0.512613   -6.314438

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0
 Excitation energy in au,    moment :                0.512613   -0.203752

 B term contribution:              0.643290
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             4    3    0

 omega B, excitation energy, moment :    0.000000    0.522952    5.720557

 First order moment in a.u. for
 C operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             4    3    0
 Excitation energy in au,    moment :                0.522952    0.584102

 B term contribution:              1.670695
 ------------------------------------------- 
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.000000    0.316502   -2.981123
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.000000    0.427812   -2.398043
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.7 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.000000    0.443797    0.910132
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YANGMOM     2    0
 Excited state no., symmetry, spin:             4    4    0

 omega B, excitation energy, moment :    0.000000    0.591276  -28.982597
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.000000    0.316502   -2.274185
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.000000    0.427812    4.844419
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.000000    0.443797  -10.412326
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.8 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XANGMOM     3    0
 Excited state no., symmetry, spin:             4    4    0

 omega B, excitation energy, moment :    0.000000    0.591276   -1.273230
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.000000    0.316502    4.255232
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.000000    0.427812   -1.014525
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.000000    0.443797    0.417545
 No MCD for this SOMOM
 Electrons: 10.000000(-3.27e-07): QR-DFT/b evaluation time:       1.6 s


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZANGMOM     4    0
 Excited state no., symmetry, spin:             4    4    0

 omega B, excitation energy, moment :    0.000000    0.591276    7.985642
 No MCD for this SOMOM

 >>>> Total CPU  time used in RESPONSE:  4 minutes  8 seconds
 >>>> Total wall time used in RESPONSE:  8 minutes 55 seconds
 >>>> Total CPU  time used in DALTON:  4 minutes 21 seconds
 >>>> Total wall time used in DALTON:  9 minutes 22 seconds

 
     Date and time (Linux)  : Tue Mar 21 15:09:03 2006
     Host name              : star.chem.uit.no                        
END REFOUT

