########## Test description ########################
START DESCRIPTION
KEYWORDS: formaldehyde hf nosym pcm linear response short
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enehf
tes
diplen
addlf
symop
diploc
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
.LOCFLD
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
.DIPLOC
**WAVE FUNCTIONS
.HF
**RESPONSE
*LINEAR
.DIPLEN
.DIPLF
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
STO-3G
Opt 2.2 paracyclophane chormoph.

    3    0         1
        6.    1
C     0.000000     0.000000     0.000000
        8.    1
O     0.000000     0.000000     1.220000
        1.    2
H     0.943102     0.000000    -0.544500
H    -0.943102     0.000000    -0.544500
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

     Date and time (Linux)  : Thu Sep 24 00:48:08 2009
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

     ICESPH =       2     NESFP =       4
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
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
 1: Opt 2.2 paracyclophane chormoph.                                        
 2:                                                                         
    ------------------------------------------------------------------------

  Coordinates are entered in Angstrom and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   6 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   8 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    2
  Used basis set file for basis set for elements with Z =   1 :
     "/home/ruud/DaltonFix/dalton/basis/STO-3G"


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 
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
                           H           1.007825
                           H           1.007825

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
  C       :    1  x   0.0000000000   2  y   0.0000000000   3  z   0.0000000000
  O       :    4  x   0.0000000000   5  y   0.0000000000   6  z   2.3054658834
  H       :    7  x   1.7822044964   8  y   0.0000000000   9  z  -1.0289558799
  H       :   10  x  -1.7822044964  11  y   0.0000000000  12  z  -1.0289558799


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H     :    1.089000    2.000725    0.000000
 H     :    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstrom (    3.7808 Bohr)
  between atoms    3 and    2, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  O          C            1.220000
  bond distance:  H          C            1.089000
  bond distance:  H          C            1.089000


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O          C          H            120.000
  bond angle:     O          C          H            120.000
  bond angle:     H          C          H            120.000




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

 Center of mass :      0.000000000000      0.000000000000      1.159648807704
 Operator center:      0.000000000000      0.000000000000      0.000000000000
 Gauge origin   :      0.000000000000      0.000000000000      0.000000000000
 Dipole origin  :      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:        1505 ( 48.8% )
 Megabytes written:                              0.021


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

 >>> Time used in Q-MAT is   0.20 seconds

 >>>> Total CPU  time used in HERMIT:   0.37 seconds
 >>>> Total wall time used in HERMIT:   1.00 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          T
     Using EHT  model:          F
     Number of Huckel orbitals each symmetry:   12

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
        which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -20.684619     -11.351900      -1.631288      -1.046947      -0.814716
           -0.700208      -0.605487      -0.493669      -0.322892      -0.164075
           -0.130192      -0.105106

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Sep 24 00:48:09 2009
     Host name              : stallo-2.local                          

 Title lines from ".mol" input file:
     Opt 2.2 paracyclophane chormoph.                                        
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.


     Time-dependent Hartree-Fock calculation (random phase approximation).


 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species          All    1
                                       --
     Occupied SCF orbitals               8    8
     Secondary orbitals                  4    4
     Total number of orbitals           12   12
     Number of basis functions          12   12

     Optimization information
     ========================
     Number of configurations              1
     Number of orbital rotations          32
     ---------------------------------------
     Total number of variables            33

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05

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

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00e-15
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.70694    60.59333    60.59323   -30.26166    -0.02185
   1  -112.142757100     -2.184941167932e-02    1.57297e+00   -1.12e+02
 MULPOP C     4.90; O    16.64; H     2.26; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.51509    60.51878    60.51878   -30.26166    -0.00042
   2  -112.306432507     -4.221971561265e-04    8.59832e-01   -1.64e-01
 MULPOP C     6.41; O    15.71; H     1.77; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59920    60.55819    60.55817   -30.26166    -0.00308
   3  -112.355348559     -3.080924378981e-03    4.78592e-02   -4.89e-02
 MULPOP C     5.92; O    16.00; H     1.86; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59745    60.55751    60.55748   -30.26166    -0.00289
   4  -112.355567191     -2.888297136007e-03    1.11568e-02   -2.19e-04
 MULPOP C     5.93; O    16.00; H     1.86; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59686    60.55723    60.55720   -30.26166    -0.00287
   5  -112.355580051     -2.869260078207e-03    3.15920e-03   -1.29e-05
 MULPOP C     5.93; O    16.00; H     1.86; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59653    60.55708    60.55705   -30.26166    -0.00286
   6  -112.355581951     -2.857409059320e-03    6.54862e-04   -1.90e-06
 MULPOP C     5.93; O    16.00; H     1.86; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59649    60.55706    60.55703   -30.26166    -0.00286
   7  -112.355582031     -2.856183337382e-03    7.06045e-05   -8.02e-08
 MULPOP C     5.93; O    16.00; H     1.86; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT   -60.59649    60.55706    60.55703   -30.26166    -0.00286
   8  -112.355582032     -2.856043607419e-03    4.54100e-06   -6.03e-10
 DIIS converged in   8 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Sym       Hartree-Fock orbital energies

  1    -20.29958353   -11.12079213    -1.33320845    -0.80284544    -0.63846096
        -0.54122299    -0.43900852    -0.35348957     0.28385000     0.63726330
         0.77009303     0.90607779

    E(LUMO) :     0.28385000 au (symmetry 1)
  - E(HOMO) :    -0.35348957 au (symmetry 1)
  ------------------------------------------
    gap     :     0.63733957 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -112.355582031674                 
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -143.516399570032

     Final gradient norm:           0.000004541003

 
     Date and time (Linux)  : Thu Sep 24 00:48:11 2009
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s    -0.0005  -0.9926   0.1223  -0.1873   0.0000   0.0268   0.0000
   2 C   :1s     0.0072  -0.0332  -0.2782   0.5808   0.0000  -0.0860   0.0000
   3 C   :2px    0.0000   0.0000   0.0000   0.0000   0.5352   0.0000   0.0000
   4 C   :2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000  -0.6051
   5 C   :2pz    0.0062  -0.0008  -0.1567  -0.2136   0.0000  -0.4552   0.0000
   6 O   :1s    -0.9943  -0.0002   0.2195   0.1022   0.0000  -0.0895   0.0000
   7 O   :1s    -0.0259   0.0058  -0.7701  -0.4468   0.0000   0.4791   0.0000
   8 O   :2px    0.0000   0.0000   0.0000   0.0000   0.4297   0.0000   0.0000
   9 O   :2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000  -0.6805
  10 O   :2pz    0.0056  -0.0018   0.1657  -0.1775   0.0000   0.6779   0.0000
  11 H   :1s    -0.0003   0.0066  -0.0324   0.2612   0.2994   0.1639   0.0000
  12 H   :1s    -0.0003   0.0066  -0.0324   0.2612  -0.2994   0.1639   0.0000

 Orbital           8        9       10
   1 C   :1s     0.0000   0.0000   0.2009
   2 C   :1s     0.0000   0.0000  -1.2706
   3 C   :2px    0.1836   0.0000   0.0000
   4 C   :2py    0.0000   0.8238   0.0000
   5 C   :2pz    0.0000   0.0000   0.5035
   6 O   :1s     0.0000   0.0000  -0.0197
   7 O   :1s     0.0000   0.0000   0.1040
   8 O   :2px   -0.8807   0.0000   0.0000
   9 O   :2py    0.0000  -0.7628   0.0000
  10 O   :2pz    0.0000   0.0000  -0.1752
  11 H   :1s     0.3445   0.0000   0.9049
  12 H   :1s    -0.3445   0.0000   0.9049



 >>>> Total CPU  time used in SIRIUS :      1.76 seconds
 >>>> Total wall time used in SIRIUS :      2.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:48:11 2009
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




 Linear Response calculation
 ---------------------------

 Equilibrium PCM solvent model requested        : SOLVNT =T

 Dielectric constant                            : EPSOL  = 78.3900
 Print level                                    : IPRLR  =   2
 Maximum number of iterations                   : MAXITL =  60
 Threshold for relative convergence             : THCLR  = 1.000e-03
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

  1 B-frequencies  0.000000e+00
 ----ADDING LOCAL FIELD CONTRIBUTION----

    6 second order properties calculated with symmetry no.    1 and labels:

          XDIPLEN 
          YDIPLEN 
          ZDIPLEN 
          XDIPLOC 
          YDIPLOC 
          ZDIPLOC 


   SCF energy         :     -112.355582031674473
 -- inactive part     :     -143.516399570032206
 -- nuclear repulsion :       31.163673581965138


                     ***************************************
                     *** RHF response calculation (TDHF) ***
                     ***************************************



 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            0
 Number of response properties of this symmetry    6
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      32
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      32


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : XDIPLEN 
 RSPLR -- frequencies :   0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   1.91e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : YDIPLEN 
 RSPLR -- frequencies :   0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   8.90e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : ZDIPLEN 
 RSPLR -- frequencies :   0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   8.58e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : XDIPLOC 
 RSPLR -- frequencies :   0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   1.90e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : YDIPLOC 
 RSPLR -- frequencies :   0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   8.92e-04

 *** RSPCTL MICROITERATIONS CONVERGED


 RSPLR -- linear response calculation for symmetry  1
 RSPLR -- operator label : ZDIPLOC 
 RSPLR -- frequencies :   0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00e-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   8.54e-04

 *** RSPCTL MICROITERATIONS CONVERGED


           Final output of second order properties from linear response
           ------------------------------------------------------------


@ Spin symmetry of operators: singlet

 Note that minus the linear response function: - << A; B >>(omega) is printed.
 The results are of quadratic accuracy using Sellers formula.

@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES

@ -<< XDIPLEN  ; XDIPLEN  >> =  6.647892730768e+00
@ -<< XDIPLEN  ; YDIPLEN  >> = -4.411291182733e-16
@ -<< XDIPLEN  ; ZDIPLEN  >> = -1.713808339218e-14
@ -<< XDIPLEN  ; XDIPLOC  >> =  9.496466958899e+00
@ -<< XDIPLEN  ; YDIPLOC  >> =  6.144589330454e-16
@ -<< XDIPLEN  ; ZDIPLOC  >> = -2.284810192689e-14
@ -<< YDIPLEN  ; YDIPLEN  >> =  2.548961573109e+00
@ -<< YDIPLEN  ; ZDIPLEN  >> = -3.767602184336e-16
@ -<< YDIPLEN  ; XDIPLOC  >> = -5.163930339203e-16
@ -<< YDIPLEN  ; YDIPLOC  >> =  4.077312570775e+00
@ -<< YDIPLEN  ; ZDIPLOC  >> = -7.042477865380e-16
@ -<< ZDIPLEN  ; ZDIPLEN  >> =  1.153680762600e+01
@ -<< ZDIPLEN  ; XDIPLOC  >> = -2.566606692893e-14
@ -<< ZDIPLEN  ; YDIPLOC  >> =  1.380266650684e-16
@ -<< ZDIPLEN  ; ZDIPLOC  >> =  1.548310723583e+01
@ -<< XDIPLOC  ; XDIPLOC  >> =  1.357527567350e+01
@ -<< XDIPLOC  ; YDIPLOC  >> =  7.753112587007e-16
@ -<< XDIPLOC  ; ZDIPLOC  >> = -3.526240434256e-14
@ -<< YDIPLOC  ; YDIPLOC  >> =  6.522099985054e+00
@ -<< YDIPLOC  ; ZDIPLOC  >> = -2.969001270901e-16
@ -<< ZDIPLOC  ; ZDIPLOC  >> =  2.078619944336e+01


 Time used in linear response calculation is      4.70 CPU seconds for symmetry 1

 >>>> Total CPU  time used in RESPONSE:   5.06 seconds
 >>>> Total wall time used in RESPONSE:   5.00 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   7.20 seconds
 >>>> Total wall time used in DALTON:   8.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:48:16 2009
     Host name              : stallo-2.local                          
END REFOUT

