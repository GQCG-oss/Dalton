########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm polarizability linear response hf formaldehyde essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enehf
tes
nuc
surf
nucchg
diplen
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
*END OF
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS


Atomtypes=3 Charge=+2 Angstrom Generators=2 X Y
Charge=6.0 Atoms=1 Basis=STO-3G  Sphere=1
C     0.000000     0.000000     0.000000   1.60
Charge=8.0 Atoms=1 Basis=STO-3G  Sphere=1
O     0.000000     0.000000     1.220000   1.50
Charge=1.0 Atoms=1 Basis=STO-3G  Sphere=1
H     0.943102     0.000000    -0.544500   1.20

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

     Date and time (Linux)  : Thu Feb  2 12:55:59 2006
     Host name              : platina.chem.uit.no                     

 <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>


 Default print level:        0

    Integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed
    Dynamic molecular property section will be executed

 ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT=WATE    **
 OPTICAL AND PHYSICAL CONSTANTS:
 EPS= 78.390; EPSINF=  1.776; RSOLV=  1.385 A; VMOL=  18.070 ML/MOL;
 TCE= .25700E-03 1/K; STEN= 71.810 DYN/CM;  DSTEN=  0.6500; CMF=  1.2770


     -----------------------------------
     INPUT FOR PCM SOLVATION CALCULATION 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=WATER        EPS   = 78.3900     EPSINF=  1.7760
     RSOLV =  1.3850

     ICESPH =       3     NESFP =       0
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

                                                                          
                                                                          

  Coordinates are entered in Angstroms and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

 For atomic type no.  1 the basis set
     STO-3G  Sphere=1                                                                
 from the basis set library will be used.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/luca/programs/dalton/basis/STO-3G"

 Huckel basis read for atomic type no.  1

 For atomic type no.  2 the basis set
     STO-3G  Sphere=1                                                                
 from the basis set library will be used.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/luca/programs/dalton/basis/STO-3G"

 Huckel basis read for atomic type no.  2

 For atomic type no.  3 the basis set
     STO-3G  Sphere=1                                                                
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
 ********SPHERES IN SPHGEN************
 INDEX        X        Y         Z        R
   1    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.1600000000E+01
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

  label    atoms   charge   prim    cont     basis
  ----------------------------------------------------------------------
  C           1  6.0000      15       5      [6s3p|2s1p]                                        
  O           1  8.0000      15       5      [6s3p|2s1p]                                        
  H           2  1.0000       3       1      [3s|1s]                                            
  ----------------------------------------------------------------------
  total:      4 16.0000      36      12
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



 >>> Time used in TWOINT is   0.03 seconds


 MEMORY USED TO GENERATE CAVITY=    432042

 TESSERA SPEZZATA IN TRONCONI

 TOTAL NUMBER OF SPHERES=    4
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.000000000    0.000000000    0.000000000    1.920000000   17.404313423
   2    0.000000000    0.000000000    1.220000000    1.800000000   25.187332839
   3    0.943102000    0.000000000   -0.544500000    1.440000000   11.255487319
   4   -0.943102000    0.000000000   -0.544500000    1.440000000   11.255487319

 TOTAL NUMBER OF TESSERAE =     268
 SURFACE AREA=   65.10262090(A**2)    CAVITY VOLUME=   45.75393835 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.11 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.07 seconds

 >>>> Total CPU  time used in HERMIT:   0.17 seconds
 >>>> Total wall time used in HERMIT:   0.00 seconds

- End of Integral Section


Starting in Wave Function Section -


 *** Output from Huckel module :

     Using EWMO model:          F
     Using EHT  model:          T
     Number of Huckel orbitals each symmetry:    7    3    2    0

 Huckel EHT eigenvalues for symmetry :  1
          -20.811796     -11.563810      -2.122567      -1.266745      -0.611990
            0.250043       0.404187

 Huckel EHT eigenvalues for symmetry :  2
           -1.021958      -0.446395       0.163030

 Huckel EHT eigenvalues for symmetry :  3
           -0.795010      -0.228090
 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Feb  2 12:55:59 2006
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

 ESTIMATE OF NUCLEAR CHARGE       15.96782
 NUCLEAR APPARENT CHARGE -15.78910 THEORETICAL -15.79589 NOT RENORMALIZED
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

   1   -111.340826376095     -0.485263209759   3.77E+00  -1.11E+02    5  1  1  0
   2   -111.674551584278     -0.489318501668   5.70E-01  -3.34E-01    5  1  1  0
   3   -111.697233198425     -0.477820434081   4.05E-01  -2.27E-02    5  1  1  0
   4   -111.710670810978     -0.481562870194   1.13E-01  -1.34E-02    5  1  1  0
   5   -111.715588794214     -0.481384137683   4.46E-02  -4.92E-03    5  1  1  0
   6   -111.716478518732     -0.481358694633   7.25E-03  -8.90E-04    5  1  1  0
   7   -111.716484433607     -0.481385864171   2.02E-03  -5.91E-06    5  1  1  0
   8   -111.716485628669     -0.481381571469   1.82E-04  -1.20E-06    5  1  1  0
   9   -111.716485633919     -0.481381211729   4.00E-05  -5.25E-09    5  1  1  0
  10   -111.716485634505     -0.481381343792   5.58E-06  -5.86E-10    5  1  1  0
  11   -111.716485634512     -0.481381338316   5.62E-07  -7.01E-12    5  1  1  0
 DIIS converged in  11 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   14
 Orbital occupations :    5    1    1    0

 Hartree-Fock orbital energies, symmetry 1

       -21.08782667   -11.61411484    -1.81843534    -1.14185404    -0.97076131
         0.34806985     0.42545047

 Hartree-Fock orbital energies, symmetry 2

        -0.93988397    -0.29794416     0.44229639

 Hartree-Fock orbital energies, symmetry 3

        -0.94384774    -0.11091781

    E(LUMO) :    -0.29794416 au (symmetry 2)
  - E(HOMO) :    -0.93988397 au (symmetry 2)
  ------------------------------------------
    gap     :     0.64193982 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    2

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -111.716485634512
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -142.398777878161

     Final gradient norm:           0.000000561911

 
     Date and time (Linux)  : Thu Feb  2 12:56:01 2006
     Host name              : platina.chem.uit.no                     

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5
   1  C   1s     0.0005   0.9929  -0.1088   0.1926   0.0715
   2  C   1s    -0.0068   0.0317   0.2142  -0.6223  -0.2853
   3  C   2pz   -0.0058   0.0020   0.1429   0.0704   0.4767
   4  O   1s     0.9947   0.0002  -0.2225  -0.1222   0.0533
   5  O   1s     0.0242  -0.0060   0.8000   0.5677  -0.2851
   6  O   2pz   -0.0054   0.0019  -0.2375   0.4234  -0.6099
   7  H   1s     0.0003  -0.0061   0.0182  -0.1406  -0.1779

     Molecular orbitals for symmetry species   2

 Orbital          1
   1  C   2px    0.6699
   2  O   2px    0.3562
   3  H   1s     0.2154

     Molecular orbitals for symmetry species   3

 Orbital          1
   1  C   2py    0.3049
   2  O   2py    0.8913



 >>>> Total CPU  time used in SIRIUS :      2.28 seconds
 >>>> Total wall time used in SIRIUS :      2.00 seconds

 
     Date and time (Linux)  : Thu Feb  2 12:56:01 2006
     Host name              : platina.chem.uit.no                     

- End of Wave Function Section



  This is output from RESPONSE  -  an MCSCF and SOPPA response property program
 ------------------------------------------------------------------------------



 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




 Linear Response calculation
 ---------------------------

 Equilibrium PCM solvent model requested        : SOLVNT =T

 Dielectric constant                            : EPSOL  = 78.3900
 Print level                                    : IPRLR  =   2
 Maximum number of iterations                   : MAXITL =  60
 Threshold for relative convergence             : THCLR  = 1.000E-04
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

  1 B-frequencies  0.000000E+00

    1 second order properties calculated with symmetry no.    1 and labels:

          ZDIPLEN 

    1 second order properties calculated with symmetry no.    2 and labels:

          XDIPLEN 

    1 second order properties calculated with symmetry no.    3 and labels:

          YDIPLEN 


   SCF energy         :     -111.716485634512381
 -- inactive part     :     -142.398777878161241
 -- nuclear repulsion :       31.163673581965142


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



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   4.74E-05

 *** RSPCTL MICROITERATIONS CONVERGED


       Final output of second order properties from linear response
       ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< ZDIPLEN  ; ZDIPLEN  >> =  1.211570273761E+01


 Time used in linear response calculation is      0.55 CPU seconds for symmetry 1


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



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  2; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   7.54E-06

 *** RSPCTL MICROITERATIONS CONVERGED


       Final output of second order properties from linear response
       ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< XDIPLEN  ; XDIPLEN  >> =  7.327579640227E+00


 Time used in linear response calculation is      0.49 CPU seconds for symmetry 2


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



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  3; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   9.77E-05

 *** RSPCTL MICROITERATIONS CONVERGED


       Final output of second order properties from linear response
       ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< YDIPLEN  ; YDIPLEN  >> =  2.908428615115E+00


 Time used in linear response calculation is      0.33 CPU seconds for symmetry 3

 >>>> Total CPU  time used in RESPONSE:   1.38 seconds
 >>>> Total wall time used in RESPONSE:   1.00 seconds
 >>>> Total CPU  time used in DALTON:   3.84 seconds
 >>>> Total wall time used in DALTON:   3.00 seconds

 
     Date and time (Linux)  : Thu Feb  2 12:56:02 2006
     Host name              : platina.chem.uit.no                     
END REFOUT

