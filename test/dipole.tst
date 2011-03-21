########## Test description ########################
START DESCRIPTION
HF/6-31G* calculation of dipole moment of ammonia with PCM (water as solvent)
KEYWORDS: pcm hf dipole short
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enedft
nuc
tes
sym
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN PROPERTIES
*PCM
.SOLVNT
WATER
.NPCMMT
0
.ICESPH
1
.NESFP
4
!.LOCFLD
*PCMCAV
.CENTER
 0.318744   0.000000  0.000000
-0.103074   0.476149  0.824714
-0.103074   0.476149 -0.824714
-0.103074  -0.952298  0.000000
.RIN
1.5
1.2
1.2
1.2
.ALPHA
1.2
1.1
1.1
1.1
.AREATS
0.3
**WAVE FUNCTIONS
.DFT
LDA
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
6-31G*
NH3 PCM SCF calculation with symmetry

    2    0
        7.    1
N          0.6023400            0.00000000              0.00000000
        1.    3
H         -0.1947831            0.899791                1.558484
H         -0.1947831            0.899791               -1.558484
H         -0.1947831           -1.799583                0.000000
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

     Date and time (Linux)  : Thu Feb  2 17:50:32 2006
     Host name              : platina.chem.uit.no                     

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

     ICESPH =       1     NESFP =       4
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.3187    0.0000    0.0000    1.5000
     2   -0.1031    0.4761    0.8247    1.2000
     3   -0.1031    0.4761   -0.8247    1.2000
     4   -0.1031   -0.9523    0.0000    1.2000

Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************



 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

  NH3 PCM SCF calculation with symmetry                                   
                                                                          
  Used basis set file for basis set for elements with Z =   7 :
     "/home/luca/programs/dalton/basis/6-31G*"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/luca/programs/dalton/basis/6-31G*"


                      SYMGRP:Point group information
                      ------------------------------

Point group: C1 

   * Character table

        |  E 
   -----+-----
    A   |   1

   * Direct product table

        | A  
   -----+-----
    A   | A  
SPHGEN: SPHERES CENTERS GIVEN FROM INPUT
 ********SPHERES IN SPHGEN************
 INDEX        X        Y         Z        R
   1    0.6023388202E+00    0.0000000000E+00    0.0000000000E+00    0.1500000000E+01
   2   -0.1947816165E+00    0.8997911393E+00    0.1558483478E+01    0.1200000000E+01
   3   -0.1947816165E+00    0.8997911393E+00   -0.1558483478E+01    0.1200000000E+01
   4   -0.1947816165E+00   -0.1799582279E+01    0.0000000000E+00    0.1200000000E+01


                             Isotopic Masses
                             ---------------

                           N          14.003074
                           H           1.007825
                           H           1.007825
                           H           1.007825

                       Total mass:    17.026549 amu
                       Natural abundance:  99.585 %

 Center-of-mass coordinates (A):    0.460792    0.000000    0.000000


  Atoms and basis sets
  --------------------

  Number of atom types:     2
  Total number of atoms:    4

  Basis set used is "6-31G*" from the basis set library.

  label    atoms   charge   prim    cont     basis
  ----------------------------------------------------------------------
  N           1  7.0000      27      14      [10s4p1d|3s2p1d]                                   
  H           1  1.0000       4       2      [4s|2s]                                            
  H           1  1.0000       4       2      [4s|2s]                                            
  H           1  1.0000       4       2      [4s|2s]                                            
  ----------------------------------------------------------------------
  total:      4 10.0000      39      20
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00E-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates:   12

   1   N        x      0.6023400000
   2            y      0.0000000000
   3            z      0.0000000000

   4   H        x     -0.1947831000
   5            y      0.8997910000
   6            z      1.5584840000

   7   H        x     -0.1947831000
   8            y      0.8997910000
   9            z     -1.5584840000

  10   H        x     -0.1947831000
  11            y     -1.7995830000
  12            z      0.0000000000


   Interatomic separations (in Angstroms):
   ---------------------------------------

            N           H           H           H     
            ------      ------      ------      ------
 N     :    0.000000
 H     :    1.041539    0.000000
 H     :    1.041539    1.649428    0.000000
 H     :    1.041539    1.649429    1.649429    0.000000


  Max interatomic separation is    1.6494 Angstroms
  between atoms "H     " and "H     ".


  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  H          N            1.041539
  bond distance:  H          N            1.041539
  bond distance:  H          N            1.041539


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     H          N          H            104.712
  bond angle:     H          N          H            104.712
  bond angle:     H          N          H            104.712




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA    1.813393          0.000000    1.000000    0.000000
   IB    1.813394          0.000000    0.000000    1.000000
   IC    2.741903          1.000000    0.000000    0.000000


 Rotational constants
 --------------------

               A                   B                   C

         278692.4410         278692.3584         184316.8481 MHz
            9.296179            9.296176            6.148148 cm-1


  Nuclear repulsion energy :   11.631995582271


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.01 seconds


 >>> Time used in ONEDRV is   0.01 seconds


 Number of two-electron integrals written:       15646 ( 70.6% )
 Megabytes written:                              0.185



 >>> Time used in TWOINT is   0.07 seconds


 MEMORY USED TO GENERATE CAVITY=    432042


 TOTAL NUMBER OF SPHERES=    4
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.318743975    0.000000000    0.000000000    1.800000000   21.704706410
   2   -0.103073992    0.476148963    0.824713936    1.320000000    9.304256586
   3   -0.103073992    0.476148963   -0.824713936    1.320000000    9.304256586
   4   -0.103073992   -0.952297926    0.000000000    1.320000000    9.304258618

 TOTAL NUMBER OF TESSERAE =     228
 SURFACE AREA=   49.61747820 (A**2)    CAVITY VOLUME=   30.98649979 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.15 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.35 seconds

 >>>> Total CPU  time used in HERMIT:   0.49 seconds
 >>>> Total wall time used in HERMIT:   1.00 seconds

- End of Integral Section


Starting in Wave Function Section -


 *** Output from Huckel module :

     Using EWMO model:          T
     Using EHT  model:          F
     Number of Huckel orbitals each symmetry:    8

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
 which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -15.673821      -1.549901      -0.689979      -0.689979      -0.543845
           -0.196824      -0.155926      -0.155926
 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Feb  2 17:50:33 2006
     Host name              : platina.chem.uit.no                     

 Title lines from integral program:
     NH3 PCM SCF calculation with symmetry                                   
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham calculation.


 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     Number of closed shell electrons         10
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Number of active orbitals                 0
     Total number of orbitals                 20

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1
 
     This is a DFT calculation of type: LDA

     Orbital specifications
     ======================
     Abelian symmetry species           1
                                       --
     Total number of orbitals          20
     Number of basis functions         20

      ** Automatic occupation of RKS orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals              5

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00E-06

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE        9.97273
 NUCLEAR APPARENT CHARGE  -9.86925 THEORETICAL  -9.87243 NOT RENORMALIZED
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
 Atom:    1*1 points=21564 compressed from 21564 (114 radial)
 Atom:    2*1 points=18545 compressed from 18666 ( 81 radial)
 Atom:    3*1 points=18545 compressed from 18666 ( 81 radial)
 Atom:    4*1 points=18551 compressed from 18666 ( 81 radial)
 Number of grid points:    77205 Grid generation time:       0.3 s
K-S electrons/energy :    9.99999901478082   -7.37607695054239 err:-.99E-06
   1    -55.885882874684     -0.012262819154   1.56E+00  -5.59E+01    5
K-S electrons/energy :    9.99999891877401   -7.62625978424713 err:-.11E-05
   2    -56.036423587254     -0.010695592346   7.59E-01  -1.51E-01    5
K-S electrons/energy :    9.99999898335832   -7.42936964809529 err:-.10E-05
   3    -56.056822437272     -0.009132111874   4.64E-01  -2.04E-02    5
K-S electrons/energy :    9.99999895909166   -7.53560894158606 err:-.10E-05
   4    -56.067402157296     -0.011939130256   1.27E-01  -1.06E-02    5
K-S electrons/energy :    9.99999896293748   -7.51257823808632 err:-.10E-05
   5    -56.068281808402     -0.011192329548   1.81E-03  -8.80E-04    5
K-S electrons/energy :    9.99999896291478   -7.51270645170133 err:-.10E-05
   6    -56.068282043162     -0.011182912369   2.02E-04  -2.35E-07    5
K-S electrons/energy :    9.99999896291057   -7.51273147067074 err:-.10E-05
   7    -56.068282045284     -0.011183988028   1.99E-05  -2.12E-09    5
K-S electrons/energy :    9.99999896290967   -7.51273523316053 err:-.10E-05
   8    -56.068282045308     -0.011184033017   3.83E-07  -2.43E-11    5
 DIIS converged in   8 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   10
 Orbital occupations :    5

 Kohn-Sham orbital energies, symmetry 1

       -13.82566873    -0.74770777    -0.38843243    -0.38843082    -0.21030104
         0.06459315     0.14535117     0.14535681     0.64183750     0.64184790

    E(LUMO) :     0.06459315 au (symmetry 1)
  - E(HOMO) :    -0.21030104 au (symmetry 1)
  ------------------------------------------
    gap     :     0.27489419 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final DFT energy:            -56.068282045308
     Nuclear repulsion:            11.631995582271
     Electronic energy:           -67.689093594562	

     Final gradient norm:           0.000000383462

 
     Date and time (Linux)  : Thu Feb  2 17:50:47 2006
     Host name              : platina.chem.uit.no                     

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5
   1  N   1s     0.9942  -0.2019   0.0000   0.0000   0.0736
   2  N   1s     0.0372   0.4055   0.0000   0.0000  -0.1510
   3  N   1s    -0.0225   0.4580   0.0000   0.0000  -0.3051
   4  N   2px   -0.0019  -0.1210   0.0000   0.0000  -0.5407
   5  N   2py    0.0000   0.0000  -0.4859   0.0000   0.0000
   6  N   2pz    0.0000   0.0000   0.0000  -0.4859   0.0000
   7  N   2px    0.0035  -0.0609   0.0000   0.0000  -0.4837
   8  N   2py    0.0000   0.0000  -0.2724   0.0000   0.0000
   9  N   2pz    0.0000   0.0000   0.0000  -0.2724   0.0000
  10  N   3d2-   0.0000   0.0000   0.0331   0.0000   0.0000
  11  N   3d1-   0.0000   0.0000   0.0000  -0.0136   0.0000
  12  N   3d0    0.0000   0.0042  -0.0118   0.0000  -0.0158
  13  N   3d1+   0.0000   0.0000   0.0000   0.0331   0.0000
  14  N   3d2+   0.0000  -0.0073  -0.0068   0.0000   0.0273
  15  H   1s    -0.0003   0.1323  -0.1318  -0.2283   0.0672
  16  H   1s     0.0046   0.0113  -0.0934  -0.1618   0.0507
  17  H   1s    -0.0003   0.1323  -0.1318   0.2283   0.0672
  18  H   1s     0.0046   0.0113  -0.0934   0.1618   0.0507
  19  H   1s    -0.0003   0.1323   0.2637   0.0000   0.0672
  20  H   1s     0.0046   0.0113   0.1868   0.0000   0.0507



 >>>> Total CPU  time used in SIRIUS :     13.67 seconds
 >>>> Total wall time used in SIRIUS :     14.00 seconds

 
     Date and time (Linux)  : Thu Feb  2 17:50:47 2006
     Host name              : platina.chem.uit.no                     

- End of Wave Function Section



    *****************************************************************
    ******** Output from **PROPE input processing for ABACUS ********
    *****************************************************************



 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------



 Default print level:        0


Starting in Static Property Section -


 
     Date and time (Linux)  : Thu Feb  2 17:50:47 2006
     Host name              : platina.chem.uit.no                     


 ***************************************************************************
 ************************ FINAL RESULTS FROM ABACUS ************************
 ***************************************************************************


 
     Date and time (Linux)  : Thu Feb  2 17:50:47 2006
     Host name              : platina.chem.uit.no                     



                         Molecular geometry (au)
                         -----------------------

 N          0.6023400000            0.0000000000            0.0000000000
 H         -0.1947831000            0.8997910000            1.5584840000
 H         -0.1947831000            0.8997910000           -1.5584840000
 H         -0.1947831000           -1.7995830000            0.0000000000





                    Molecular wave function and energy
                    ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0


     Total energy        -56.0682820453 au (Hartrees)
                         -1525.69552143 eV
                           -147207.2534 kJ/mol




                         Relativistic corrections
                         ------------------------

     Darwin correction:                          0.1120853974 au
     Mass-velocity correction:                  -0.1406447082 au

     Total relativistic correction:             -0.0285593108 au (0.0509%)
     Non-relativistic + relativistic energy:   -56.0968413561 au




                              Dipole moment
                              -------------



                    0.913725 au           2.322457 Debye




                         Dipole moment components
                         ------------------------

                               au             Debye

                    x     -0.91372502     -2.32245708
                    y      0.00006776      0.00017224
                    z      0.00000000      0.00000000

                        1 a.u. =   2.54175 Debye 





   Interatomic separations (in Angstroms):
   ---------------------------------------

            N           H           H           H     
            ------      ------      ------      ------
 N     :    0.000000
 H     :    1.041539    0.000000
 H     :    1.041539    1.649428    0.000000
 H     :    1.041539    1.649429    1.649429    0.000000


  Max interatomic separation is    1.6494 Angstroms
  between atoms "H     " and "H     ".


  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  H          N            1.041539
  bond distance:  H          N            1.041539
  bond distance:  H          N            1.041539


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     H          N          H            104.712
  bond angle:     H          N          H            104.712
  bond angle:     H          N          H            104.712




 CPU time statistics for ABACUS
 ------------------------------





 >>>> Total CPU  time used in ABACUS:   0.01 seconds
 >>>> Total wall time used in ABACUS:   0.00 seconds

- End of Static Property Section

 >>>> Total CPU  time used in DALTON:  14.17 seconds
 >>>> Total wall time used in DALTON:  15.00 seconds

 
     Date and time (Linux)  : Thu Feb  2 17:50:47 2006
     Host name              : platina.chem.uit.no                     
END REFOUT

