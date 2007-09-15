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

     Date and time (Linux)  : Tue Mar 21 09:59:34 2006
     Host name              : platina.chem.uit.no                     

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

     ICESPH =       1     NESFP =       1
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    6.0000

Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************


 Default print level:        1

 * Nuclear model: Point charge

 Calculation of one- and two-electron Hamiltonian integrals.


 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

                                                                          
                                                                          

  Coordinates are entered in Angstroms and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

 For atomic type no.  1 the basis set
     STO-3G                                                                          
 from the basis set library will be used.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/luca/programs/dalton/basis/STO-3G"

 Huckel basis read for atomic type no.  1

 For atomic type no.  2 the basis set
     STO-3G                                                                          
 from the basis set library will be used.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/luca/programs/dalton/basis/STO-3G"

 Huckel basis read for atomic type no.  2

 For atomic type no.  3 the basis set
     STO-3G                                                                          
 from the basis set library will be used.
  Used basis set file for basis set for elements with Z =   1 :
     "/home/luca/programs/dalton/basis/STO-3G"

 Huckel basis read for atomic type no.  3


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
SPHGEN: SPHERES CENTERS GIVEN FROM INPUT
 ********SPHERES IN SPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.6000000000E+01


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


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.01 seconds


 >>> Time used in ONEDRV is   0.00 seconds


 Number of two-electron integrals written:        1028 ( 33.4% )
 Megabytes written:                              0.014



 >>> Time used in TWOINT is   0.04 seconds


 MEMORY USED TO GENERATE CAVITY=    432042


 TOTAL NUMBER OF SPHERES=    1
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.000000000    0.000000000    0.000000000    3.175062000  126.681817203

 TOTAL NUMBER OF TESSERAE =     392
 SURFACE AREA=  126.68181720(A**2)    CAVITY VOLUME=  134.07420796 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.11 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.21 seconds

 >>>> Total CPU  time used in HERMIT:   0.32 seconds
 >>>> Total wall time used in HERMIT:   1.00 seconds

- End of Integral Section


Starting in Wave Function Section -


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

 
     Date and time (Linux)  : Tue Mar 21 09:59:35 2006
     Host name              : platina.chem.uit.no                     

 Title lines from integral program:
                                                                             
                                                                             

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

     Number of active orbitals                 0
     Total number of orbitals                 12

     Spin multiplicity                         1
     Total number of symmetries                4
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species           1   2   3   4
                                       --  --  --  --
     Total number of orbitals           7   3   2   0
     Number of basis functions          7   3   2   0

      ** Automatic occupation of RHF orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals              5   1   1   0

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00E-06

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.99669
 NUCLEAR APPARENT CHARGE -15.79079 THEORETICAL -15.79589 NOT RENORMALIZED
 this is icompcm in icvev 0

  ..... DONE WITH INDUCED NUCLEAR CHARGES .....

 >>> Time used in VNN    is   0.00 seconds



 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  14 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00E-15

 >>> Time used in FRMSUP is   0.00 seconds

   1   -111.152808841224     -0.336366680037   3.79E+00  -1.11E+02    5  1  1  0
   2   -111.516982456698     -0.329753521569   6.40E-01  -3.64E-01    5  1  1  0
   3   -111.539004079907     -0.332469505690   5.21E-01  -2.20E-02    5  1  1  0
   4   -111.559549399026     -0.330510284504   1.20E-01  -2.05E-02    5  1  1  0
   5   -111.565209024145     -0.330578566265   4.63E-02  -5.66E-03    5  1  1  0
   6   -111.566150764541     -0.330649540360   6.72E-03  -9.42E-04    5  1  1  0
   7   -111.566154970526     -0.330638386356   2.13E-03  -4.21E-06    5  1  1  0
   8   -111.566156212145     -0.330640865400   2.08E-04  -1.24E-06    5  1  1  0
   9   -111.566156220842     -0.330640942853   6.04E-05  -8.70E-09    5  1  1  0
  10   -111.566156222260     -0.330640969195   6.98E-06  -1.42E-09    5  1  1  0
  11   -111.566156222279     -0.330640961991   4.10E-07  -1.90E-11    5  1  1  0
 DIIS converged in  11 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   14
 Orbital occupations :    5    1    1    0

 Hartree-Fock orbital energies, symmetry 1

       -21.24419226   -11.75700809    -1.96710560    -1.28470071    -1.11212500
         0.20015222     0.27610554

 Hartree-Fock orbital energies, symmetry 2

        -1.09098266    -0.44093924     0.29510711

 Hartree-Fock orbital energies, symmetry 3

        -1.08694543    -0.26242750

    E(LUMO) :    -0.44093924 au (symmetry 2)
  - E(HOMO) :    -1.08694543 au (symmetry 3)
  ------------------------------------------
    gap     :     0.64600619 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    2

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -111.566156222279
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -142.399188842253

     Final gradient norm:           0.000000409815

 
     Date and time (Linux)  : Tue Mar 21 09:59:38 2006
     Host name              : platina.chem.uit.no                     

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5
   1  C   1s     0.0005   0.9929  -0.1100  -0.1963  -0.0602
   2  C   1s    -0.0068   0.0318   0.2171   0.6406   0.2459
   3  C   2pz   -0.0058   0.0020   0.1441  -0.0947  -0.4750
   4  O   1s     0.9947   0.0002  -0.2211   0.1213  -0.0622
   5  O   1s     0.0242  -0.0061   0.7939  -0.5610   0.3266
   6  O   2pz   -0.0055   0.0019  -0.2460  -0.3804   0.6286
   7  H   1s     0.0003  -0.0061   0.0185   0.1501   0.1698

     Molecular orbitals for symmetry species   2

 Orbital          1
   1  C   2px    0.6742
   2  O   2px    0.3443
   3  H   1s     0.2176

     Molecular orbitals for symmetry species   3

 Orbital          1
   1  C   2py    0.3159
   2  O   2py    0.8855



 >>>> Total CPU  time used in SIRIUS :      2.93 seconds
 >>>> Total wall time used in SIRIUS :      3.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 09:59:38 2006
     Host name              : platina.chem.uit.no                     

- End of Wave Function Section



  This is output from RESPONSE  -  an MCSCF and SOPPA response property program
 ------------------------------------------------------------------------------



 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




 Quadratic Response single residue calculation
 ---------------------------------------------

 Two-photon transition process computed         : TWOPHO =T


 Spin of operator A , ISPINA=    0
 Spin of operator B , ISPINB=    0
 Spin of operator C , (Excitation energy) ISPINC=    0

 12 B-frequencies  1.000000E+00  1.000000E+00  1.000000E+00  1.000000E+00  1.000000E+00
  1.000000E+00  1.000000E+00  1.000000E+00  1.000000E+00  1.000000E+00
  1.000000E+00  1.000000E+00

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

    1 B OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          ZDIPLEN 

    1 B OPERATORS OF SYMMETRY NO:    2 AND LABELS:

          XDIPLEN 

    1 B OPERATORS OF SYMMETRY NO:    3 AND LABELS:

          YDIPLEN 


   SCF energy         :     -111.566156222278636
 -- inactive part     :     -142.399188842252983
 -- nuclear repulsion :       31.163673581965142


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

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   22)
 RSP solution vector no.    1; norm of residual   3.71E-05
 RSP solution vector no.    2; norm of residual   1.31E-05
 RSP solution vector no.    3; norm of residual   9.90E-05

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

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   20)
 RSP solution vector no.    1; norm of residual   4.20E-04
 RSP solution vector no.    2; norm of residual   2.08E-04
 RSP solution vector no.    3; norm of residual   1.22E-04

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


 *** THE REQUESTED    3 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   14)
 RSP solution vector no.    1; norm of residual   3.45E-15
 RSP solution vector no.    2; norm of residual   6.31E-16
 RSP solution vector no.    3; norm of residual   4.19E-15

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

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    6)
 RSP solution vector no.    1; norm of residual   4.27E-16
 RSP solution vector no.    2; norm of residual   3.42E-16
 RSP solution vector no.    3; norm of residual   6.21E-16

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
 QRLRVE -- frequencies :  0.156891  0.260466  0.472415  0.077257  0.255638
                          0.432130  0.206985  0.331146  0.396124



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 GP * SOLUTION VECTOR AT FREQUENCY     0.47241 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS     19.67696840
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.01 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.43213 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS     13.02925146
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.68 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.33115 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS    -13.06314498
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -8.05 * 10 **  5.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.39612 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS      9.11046231
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.64 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.47241 AU
 AFTER   13 LINEAR TRANSFORMATIONS IS     26.40472949
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   25    0    1
 DETERMINANT OF REDUCED MATRIX IS -8.47 * 10 ** 15.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   25    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.43213 AU
 AFTER   13 LINEAR TRANSFORMATIONS IS     16.01893246
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   25    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.31 * 10 ** 16.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   25    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.33115 AU
 AFTER   13 LINEAR TRANSFORMATIONS IS     -7.93144321
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   25    0    1
 DETERMINANT OF REDUCED MATRIX IS -4.94 * 10 ** 15.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   25    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.39612 AU
 AFTER   13 LINEAR TRANSFORMATIONS IS     11.27103417
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   25    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.40 * 10 ** 16.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   25    0    1

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   26)
 RSP solution vector no.    1; norm of residual   8.87E-16
 RSP solution vector no.    2; norm of residual   5.88E-16
 RSP solution vector no.    3; norm of residual   4.81E-16
 RSP solution vector no.    4; norm of residual   8.60E-16
 RSP solution vector no.    5; norm of residual   3.49E-16
 RSP solution vector no.    6; norm of residual   3.88E-16
 RSP solution vector no.    7; norm of residual   4.19E-16
 RSP solution vector no.    8; norm of residual   3.41E-16
 RSP solution vector no.    9; norm of residual   5.86E-16

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.156891E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.15689):     12.3942689206    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.260466E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.26047):     17.8351314643    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.472415E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.47241):     26.4047294886    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.772568E-01
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.07726):     11.3829250548    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.255638E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.25564):     17.2063910903    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.432130E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.43213):     16.0189324606    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.206985E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.20699):     13.8310020694    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.331146E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.33115):    -7.93144321182    

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.396124E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.39612):     11.2710341698    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    2


 Perturbation symmetry.     KSYMOP:       2
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      12
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      12


 QRLRVE -- linear response calculation for symmetry  2
 QRLRVE -- operator label : XDIPLEN 
 QRLRVE -- frequencies :  0.077257  0.156891  0.255638  0.260466  0.432130
                          0.472415  0.025157  0.177758  0.518615



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F


 GP * SOLUTION VECTOR AT FREQUENCY     0.15689 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS      5.40398044
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.74 * 10 **  4.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.25564 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS      8.03764421
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.07 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.26047 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS      8.10167489
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.10 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.43213 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS     15.11529641
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -7.71 * 10 **  5.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.47241 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS     25.20402313
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.80 * 10 **  5.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.17776 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS      7.17372068
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   17    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.80 * 10 **  5.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   17    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.51861 AU
 AFTER    9 LINEAR TRANSFORMATIONS IS    -97.17648572
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   16    0    2
 DETERMINANT OF REDUCED MATRIX IS  6.67 * 10 **  4.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   16    0    2

 GP * SOLUTION VECTOR AT FREQUENCY     0.15689 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      5.49595306
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.97 * 10 ** 13.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.25564 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      8.03780373
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.10 * 10 ** 15.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.26047 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      8.10183516
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -1.14 * 10 ** 15.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.43213 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     15.11744677
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -7.88 * 10 ** 14.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.47241 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS     25.21493084
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -3.87 * 10 ** 14.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.17776 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS      7.17470907
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   23    0    1
 DETERMINANT OF REDUCED MATRIX IS -2.90 * 10 ** 14.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   23    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.51861 AU
 AFTER   12 LINEAR TRANSFORMATIONS IS    -96.78293112
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   22    0    2
 DETERMINANT OF REDUCED MATRIX IS  6.81 * 10 ** 13.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   22    0    2

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   24)
 RSP solution vector no.    1; norm of residual   7.99E-16
 RSP solution vector no.    2; norm of residual   3.71E-16
 RSP solution vector no.    3; norm of residual   6.21E-16
 RSP solution vector no.    4; norm of residual   3.66E-16
 RSP solution vector no.    5; norm of residual   3.80E-16
 RSP solution vector no.    6; norm of residual   3.68E-16
 RSP solution vector no.    7; norm of residual   7.32E-16
 RSP solution vector no.    8; norm of residual   6.43E-16
 RSP solution vector no.    9; norm of residual   3.98E-16

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.772568E-01
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.07726):     6.96542824847    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.156891E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.15689):     5.49595305671    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.255638E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.25564):     8.03780373282    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.260466E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.26047):     8.10183516435    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.432130E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.43213):     15.1174467680    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.472415E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.47241):     25.2149308384    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.251572E-01
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.02516):     6.86569447195    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.177758E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.17776):     7.17470907440    

 QRLRVE: SINGLET SOLUTION   LABEL   XDIPLEN     FREQUENCY   0.518615E+00
 SYMMETRY    2

@QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.51861):    -96.7829311218    


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    3


 Perturbation symmetry.     KSYMOP:       3
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:       7
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :       7


 QRLRVE -- linear response calculation for symmetry  3
 QRLRVE -- operator label : YDIPLEN 
 QRLRVE -- frequencies :  0.260466  0.396124  0.472415  0.025157  0.206985
                          0.177758  0.156891  0.518615  0.331146



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F


 GP * SOLUTION VECTOR AT FREQUENCY     0.47241 AU
 AFTER    7 LINEAR TRANSFORMATIONS IS      4.42866313
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   13    0    1
 DETERMINANT OF REDUCED MATRIX IS -4.02 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   13    0    1

 GP * SOLUTION VECTOR AT FREQUENCY     0.51861 AU
 AFTER    7 LINEAR TRANSFORMATIONS IS      6.12077802
 INERTIA (POS,ZER,NEG) OF REDUCED MATRIX IS   13    0    1
 DETERMINANT OF REDUCED MATRIX IS -4.83 * 10 **  6.0


 ***WARNING, negative eigenvalues in reduced matrix, inertia is   13    0    1

 *** THE REQUESTED    9 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   14)
 RSP solution vector no.    1; norm of residual   4.39E-16
 RSP solution vector no.    2; norm of residual   1.93E-16
 RSP solution vector no.    3; norm of residual   1.37E-16
 RSP solution vector no.    4; norm of residual   5.61E-16
 RSP solution vector no.    5; norm of residual   8.20E-16
 RSP solution vector no.    6; norm of residual   4.72E-16
 RSP solution vector no.    7; norm of residual   3.87E-16
 RSP solution vector no.    8; norm of residual   4.49E-16
 RSP solution vector no.    9; norm of residual   3.26E-16

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.260466E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.26047):     3.40521794970    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.396124E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.39612):     6.57835186553    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.472415E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.47241):     4.42866313448    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.251572E-01
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.02516):     2.80740975919    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.206985E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.20699):     3.14758191979    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.177758E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.17776):     3.04733253670    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.156891E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.15689):     2.98896547014    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.518615E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.51861):     6.12077802040    

 QRLRVE: SINGLET SOLUTION   LABEL   YDIPLEN     FREQUENCY   0.331146E+00
 SYMMETRY    3

@QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.33115):     4.02686505261    


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156891    0.313782    4.770822


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260466    0.520931   -1.489719


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    3.196907


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156891    0.313782   20.940095


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260466    0.520931   -0.389093


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    2.952065


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    1    0

 omega B, excitation energy, moment :    0.156891    0.313782   -0.232312


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    1    0

 omega B, excitation energy, moment :    0.260466    0.520931   -0.131524


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    1    0

 omega B, excitation energy, moment :    0.472415    0.944830    1.879614


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.077257    0.154514    1.746157


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.255638    0.511276   -9.362962


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.432130    0.864260   -0.769930


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    2    0

 omega B, excitation energy, moment :    0.077257    0.154514    1.746157


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    2    0

 omega B, excitation energy, moment :    0.255638    0.511276   -9.362962


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    2    0

 omega B, excitation energy, moment :    0.432130    0.864260   -0.769930


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.206985    0.413971   -1.245237


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.331146    0.662292   -5.788484


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      ZDIPLEN     1    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.396124    0.792249    2.172689


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    3    0

 omega B, excitation energy, moment :    0.206985    0.413971   -1.245237


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    3    0

 omega B, excitation energy, moment :    0.331146    0.662292   -5.788484


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      ZDIPLEN     1    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    3    0

 omega B, excitation energy, moment :    0.396124    0.792249    2.172689


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.025157    0.050314    0.047552


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.177758    0.355516   -1.307597


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      YDIPLEN     3    0
 B operator label,  symmetry, spin:      XDIPLEN     2    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.518615    1.037229    1.149602


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             1    4    0

 omega B, excitation energy, moment :    0.025157    0.050314    0.047552


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             2    4    0

 omega B, excitation energy, moment :    0.177758    0.355516   -1.307597


 Second order moment in a.u. for
 A operator label,  symmetry, spin:      XDIPLEN     2    0
 B operator label,  symmetry, spin:      YDIPLEN     3    0
 Excited state no., symmetry, spin:             3    4    0

 omega B, excitation energy, moment :    0.518615    1.037229    1.149602


   *******************************************************************
   ************ FINAL RESULTS FROM TWO-PHOTON CALCULATION ************
   *******************************************************************

   The two-photon absorption strength for an average molecular  
   orientation is computed according to formulas given by       
   P.R. Monson and W.M. MacClain in J. Chem. Phys. 53:29, 1970  
   and W.M. MacClain in J. Chem. Phys. 55:2789, 1971.           
   The absorption depends on the light polarization.            

   A monochromatic light source is assumed.                     


   All results are presented in atomic units, except the        
   excitation energy which is given in eV.                      


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

                             Polarization ratio      
                      -------------------------------
                          R  = (-Df+3*Dg)/(Df+2*Dg)  


                   +-----------------------------------+
                   | Two-photon transition probability |
                   +-----------------------------------+
   ----------------------------------------------------------------------
   Sym  No  Energy  Polarization         Df         Dg          D       R
   ----------------------------------------------------------------------
     1   1    8.54   Linear       0.216E+02  0.154E+02  0.105E+03    0.47
     1   1    8.54   Circular     0.216E+02  0.154E+02  0.490E+02    0.47
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    | 104.78425825   |  0.258489E-18  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |  55.24506807   |  0.136282E-18  |
       +-------------------+----------------+----------------+
       |    circular       |  48.98316632   |  0.120835E-18  |
       +-------------------+----------------+----------------+
     1   2   14.18   Linear       0.135E+00  0.796E-01  0.588E+00    0.35
     1   2   14.18   Circular     0.135E+00  0.796E-01  0.208E+00    0.35
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.58782392   |  0.399666E-20  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   0.26327733   |  0.179004E-20  |
       +-------------------+----------------+----------------+
       |    circular       |   0.20816077   |  0.141530E-20  |
       +-------------------+----------------+----------------+
     1   3   25.71   Linear       0.215E+01  0.749E+00  0.729E+01    0.03
     1   3   25.71   Circular     0.215E+01  0.749E+00  0.196E+00    0.03
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   7.29292536   |  0.163117E-18  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   1.59603528   |  0.356976E-19  |
       +-------------------+----------------+----------------+
       |    circular       |   0.19635744   |  0.439182E-20  |
       +-------------------+----------------+----------------+
     2   1    4.20   Linear       0.000E+00  0.203E+00  0.813E+00    1.50
     2   1    4.20   Circular     0.000E+00  0.203E+00  0.122E+01    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.81308426   |  0.486362E-21  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   1.01635532   |  0.607952E-21  |
       +-------------------+----------------+----------------+
       |    circular       |   1.21962639   |  0.729542E-21  |
       +-------------------+----------------+----------------+
     2   2   13.91   Linear       0.000E+00  0.584E+01  0.234E+02    1.50
     2   2   13.91   Circular     0.000E+00  0.584E+01  0.351E+02    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |  23.37734661   |  0.153107E-18  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |  29.22168326   |  0.191384E-18  |
       +-------------------+----------------+----------------+
       |    circular       |  35.06601991   |  0.229661E-18  |
       +-------------------+----------------+----------------+
     2   3   23.52   Linear       0.000E+00  0.395E-01  0.158E+00    1.50
     2   3   23.52   Circular     0.000E+00  0.395E-01  0.237E+00    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.15807795   |  0.295836E-20  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   0.19759743   |  0.369795E-20  |
       +-------------------+----------------+----------------+
       |    circular       |   0.23711692   |  0.443753E-20  |
       +-------------------+----------------+----------------+
     3   1   11.26   Linear       0.000E+00  0.103E+00  0.413E+00    1.50
     3   1   11.26   Circular     0.000E+00  0.103E+00  0.620E+00    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.41349717   |  0.177542E-20  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   0.51687147   |  0.221928E-20  |
       +-------------------+----------------+----------------+
       |    circular       |   0.62024576   |  0.266314E-20  |
       +-------------------+----------------+----------------+
     3   2   18.02   Linear       0.000E+00  0.223E+01  0.894E+01    1.50
     3   2   18.02   Circular     0.000E+00  0.223E+01  0.134E+02    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   8.93507794   |  0.981946E-19  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |  11.16884743   |  0.122743E-18  |
       +-------------------+----------------+----------------+
       |    circular       |  13.40261692   |  0.147292E-18  |
       +-------------------+----------------+----------------+
     3   3   21.56   Linear       0.000E+00  0.315E+00  0.126E+01    1.50
     3   3   21.56   Circular     0.000E+00  0.315E+00  0.189E+01    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   1.25882084   |  0.197960E-19  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   1.57352605   |  0.247450E-19  |
       +-------------------+----------------+----------------+
       |    circular       |   1.88823126   |  0.296939E-19  |
       +-------------------+----------------+----------------+
     4   1    1.37   Linear       0.000E+00  0.151E-03  0.603E-03    1.50
     4   1    1.37   Circular     0.000E+00  0.151E-03  0.904E-03    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.00060299   |  0.382463E-25  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   0.00075374   |  0.478078E-25  |
       +-------------------+----------------+----------------+
       |    circular       |   0.00090449   |  0.573694E-25  |
       +-------------------+----------------+----------------+
     4   2    9.67   Linear       0.000E+00  0.114E+00  0.456E+00    1.50
     4   2    9.67   Circular     0.000E+00  0.114E+00  0.684E+00    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.45594949   |  0.144386E-20  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   0.56993687   |  0.180482E-20  |
       +-------------------+----------------+----------------+
       |    circular       |   0.68392424   |  0.216578E-20  |
       +-------------------+----------------+----------------+
     4   3   28.22   Linear       0.000E+00  0.881E-01  0.352E+00    1.50
     4   3   28.22   Circular     0.000E+00  0.881E-01  0.529E+00    1.50
 
       Rotationally averaged two-photon transition strengths
                       and rate constants
 
       +-------------------+----------------+----------------+
       |   Polarization    |    DELTA_TP    |   K (0 -> f)   |
       +-------------------+----------------+----------------+
       |  linear (para)    |   0.35242234   |  0.949955E-20  |
       +-------------------+----------------+----------------+
       |  linear (perp)    |   0.44052792   |  0.118744E-19  |
       +-------------------+----------------+----------------+
       |    circular       |   0.52863351   |  0.142493E-19  |
       +-------------------+----------------+----------------+
   ----------------------------------------------------------------------

 >>>> Total CPU  time used in RESPONSE:  20.59 seconds
 >>>> Total wall time used in RESPONSE:  20.00 seconds
 >>>> Total CPU  time used in DALTON:  23.85 seconds
 >>>> Total wall time used in DALTON:  24.00 seconds

 
     Date and time (Linux)  : Tue Mar 21 09:59:58 2006
     Host name              : platina.chem.uit.no                     
END REFOUT

