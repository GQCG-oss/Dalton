########## Test description ########################
START DESCRIPTION
PCM test
HF wavefunction with STO-3G basis
Calculation of dipole moment with and without cavity field correction
KEYWORDS: pcm hf dipole localfield essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enehf
dipole
dipolelf
OVERRIDE 7 1.0e-5
nuc
tes
sym
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
**WAVE FUNCTIONS
.HF
**PROPERTIES
*END OF INPUT
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

     Date and time (Linux)  : Mon Mar 27 14:57:13 2006
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

     ICESPH =       2     NESFP =       4
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
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


 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

  Opt 2.2 paracyclophane chormoph.                                        
                                                                          

  Coordinates are entered in Angstroms and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A
  Used basis set file for basis set for elements with Z =   6 :
     "/home/lara/programs/main-branch/dalton/basis/STO-3G"
  Used basis set file for basis set for elements with Z =   8 :
     "/home/lara/programs/main-branch/dalton/basis/STO-3G"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/lara/programs/main-branch/dalton/basis/STO-3G"


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
                           H           1.007825
                           H           1.007825

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

   7   H        x      1.7822044964
   8            y      0.0000000000
   9            z     -1.0289558799

  10   H        x     -1.7822044964
  11            y      0.0000000000
  12            z     -1.0289558799


   Interatomic separations (in Angstroms):
   ---------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H     :    1.089000    2.000725    0.000000
 H     :    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstroms
  between atoms "H     " and "O     ".


  Bond distances (angstroms):
  ---------------------------

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


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.02 seconds


 >>> Time used in ONEDRV is   0.01 seconds


 Number of two-electron integrals written:        1505 ( 48.8% )
 Megabytes written:                              0.021



 >>> Time used in TWOINT is   0.06 seconds


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

 >>> Time used in PEDRAM is   0.18 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   1.41 seconds

 >>>> Total CPU  time used in HERMIT:   1.77 seconds
 >>>> Total wall time used in HERMIT:   2.00 seconds

- End of Integral Section


Starting in Wave Function Section -


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

 
     Date and time (Linux)  : Mon Mar 27 14:57:15 2006
     Host name              : star.chem.uit.no                        

 Title lines from integral program:
     Opt 2.2 paracyclophane chormoph.                                        
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.


 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Number of active orbitals                 0
     Total number of orbitals                 12

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species           1
                                       --
     Total number of orbitals          12
     Number of basis functions         12

      ** Automatic occupation of RHF orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals              8

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

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00E-15

 >>> Time used in FRMSUP is   0.00 seconds

   1   -112.142757099773     -0.021849411679   1.57E+00  -1.12E+02    8
   2   -112.306432506707     -0.000422197156   8.60E-01  -1.64E-01    8
   3   -112.355348558796     -0.003080924379   4.79E-02  -4.89E-02    8
   4   -112.355567191024     -0.002888297136   1.12E-02  -2.19E-04    8
   5   -112.355580051192     -0.002869260078   3.16E-03  -1.29E-05    8
   6   -112.355581950837     -0.002857409059   6.55E-04  -1.90E-06    8
   7   -112.355582031072     -0.002856183337   7.06E-05  -8.02E-08    8
   8   -112.355582031675     -0.002856043607   4.54E-06  -6.03E-10    8
   9   -112.355582031677     -0.002856039431   6.29E-07  -2.96E-12    8
 DIIS converged in   9 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Hartree-Fock orbital energies, symmetry 1

       -20.29958376   -11.12079225    -1.33320855    -0.80284549    -0.63846100
        -0.54122310    -0.43900865    -0.35348947     0.28384989     0.63726322
         0.77009301     0.90607769

    E(LUMO) :     0.28384989 au (symmetry 1)
  - E(HOMO) :    -0.35348947 au (symmetry 1)
  ------------------------------------------
    gap     :     0.63733936 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -112.355582031677
     Nuclear repulsion:            31.163673581965
     Electronic energy:          -143.516399574212

     Final gradient norm:           0.000000629365

 
     Date and time (Linux)  : Mon Mar 27 14:57:20 2006
     Host name              : star.chem.uit.no                        

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5        6        7
   1  C   1s     0.0005   0.9926  -0.1224   0.1873   0.0000   0.0268   0.0000
   2  C   1s    -0.0072   0.0332   0.2782  -0.5808   0.0000  -0.0860   0.0000
   3  C   2px    0.0000   0.0000   0.0000   0.0000   0.5352   0.0000   0.0000
   4  C   2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000  -0.6051
   5  C   2pz   -0.0062   0.0008   0.1567   0.2136   0.0000  -0.4552   0.0000
   6  O   1s     0.9943   0.0002  -0.2195  -0.1022   0.0000  -0.0895   0.0000
   7  O   1s     0.0259  -0.0058   0.7701   0.4468   0.0000   0.4791   0.0000
   8  O   2px    0.0000   0.0000   0.0000   0.0000   0.4297   0.0000   0.0000
   9  O   2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000  -0.6805
  10  O   2pz   -0.0056   0.0018  -0.1657   0.1775   0.0000   0.6779   0.0000
  11  H   1s     0.0003  -0.0066   0.0324  -0.2612   0.2994   0.1639   0.0000
  12  H   1s     0.0003  -0.0066   0.0324  -0.2612  -0.2994   0.1639   0.0000

 Orbital          8
   1  C   1s     0.0000
   2  C   1s     0.0000
   3  C   2px    0.1836
   4  C   2py    0.0000
   5  C   2pz    0.0000
   6  O   1s     0.0000
   7  O   1s     0.0000
   8  O   2px   -0.8807
   9  O   2py    0.0000
  10  O   2pz    0.0000
  11  H   1s     0.3445
  12  H   1s    -0.3445



 >>>> Total CPU  time used in SIRIUS :      4.85 seconds
 >>>> Total wall time used in SIRIUS :      5.00 seconds

 
     Date and time (Linux)  : Mon Mar 27 14:57:20 2006
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


 
     Date and time (Linux)  : Mon Mar 27 14:57:20 2006
     Host name              : star.chem.uit.no                        
 ----ADDING LOCAL FIELD CONTRIBUTION----


 ***************************************************************************
 ************************ FINAL RESULTS FROM ABACUS ************************
 ***************************************************************************


 
     Date and time (Linux)  : Mon Mar 27 14:57:20 2006
     Host name              : star.chem.uit.no                        



                         Molecular geometry (au)
                         -----------------------

 C          0.0000000000            0.0000000000            0.0000000000
 O          0.0000000000            0.0000000000            2.3054658834
 H          1.7822044964            0.0000000000           -1.0289558799
 H         -1.7822044964            0.0000000000           -1.0289558799





                    Molecular wave function and energy
                    ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0


     Total energy       -112.3555820317 au (Hartrees)
                         -3057.35082404 eV
                           -294989.5383 kJ/mol




                         Relativistic corrections
                         ------------------------

     Darwin correction:                          0.1684676511 au
     Mass-velocity correction:                  -0.2161983092 au

     Total relativistic correction:             -0.0477306581 au (0.0425%)
     Non-relativistic + relativistic energy:  -112.4033126898 au




                              Dipole moment
                              -------------


                    0.675590 au           1.717177 Debye




                         Dipole moment components
                         ------------------------

                               au             Debye

                    x      0.00000000      0.00000000
                    y      0.00000000      0.00000000
                    z     -0.67558959     -1.71717727

                        1 a.u. =   2.54175 Debye 





                   Local-field corrected dipole moment
                   -----------------------------------


                    0.900539 au           2.288940 Debye




              Local-field corrected dipole moment components
              ----------------------------------------------

                               au             Debye

                    x      0.00000000      0.00000000
                    y      0.00000000      0.00000000
                    z     -0.90053851     -2.28894033

                        1 a.u. =   2.54175 Debye 





   Interatomic separations (in Angstroms):
   ---------------------------------------

            C           O           H           H     
            ------      ------      ------      ------
 C     :    0.000000
 O     :    1.220000    0.000000
 H     :    1.089000    2.000725    0.000000
 H     :    1.089000    2.000725    1.886204    0.000000


  Max interatomic separation is    2.0007 Angstroms
  between atoms "H     " and "O     ".


  Bond distances (angstroms):
  ---------------------------

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




 CPU time statistics for ABACUS
 ------------------------------





 >>>> Total CPU  time used in ABACUS:   0.61 seconds
 >>>> Total wall time used in ABACUS:   0.00 seconds

- End of Static Property Section

 >>>> Total CPU  time used in DALTON:   7.23 seconds
 >>>> Total wall time used in DALTON:   7.00 seconds

 
     Date and time (Linux)  : Mon Mar 27 14:57:20 2006
     Host name              : star.chem.uit.no                        
END REFOUT

