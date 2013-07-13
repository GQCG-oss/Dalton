########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm qm3 qm3pcm dft residue linear response essential
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enedft
qm3pcmconv
qm3ele
qm3pol
qm3wdv
qm3pcmene
pcmsol
surf
nucchg
tramom
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON
.RUN RESPONSE
*PCM
.NEQRSP
.SOLVNT
WATER
.NPCMMT
0
.NESFP
8
.ICESPH
2
*PCMCAV
.RIN
1.7
1.5
1.2
1.2
1.5
1.5
1.2
1.2
.INA
1
2
3
4
5
6
7
8
.AREATS
 0.3
.ALPHA
 1.2
 1.2
 1.2
 1.2
 1.2
 1.2
 1.2
 1.2
*QM3
.PRINT
 1
.QM3
.MMPCM
.THRSMP
 1.0D-9
**INTEGRALS
.DIPLEN
.NUCPOT
.NELFLD
.PRINT
 2
**WAVE FUNCTIONS
.DFT
B3LYP
*SCF INPUT
.THRESH
1.0D-12
**RESPONSE
*LINEAR
.DIPLEN
.SINGLE RES
.ROOTS
4
.THCPP
 1.0D-04
.PRINT
5
**END OF DALTON INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS
CH2O PLUS 2 WATERS
------------------------
AtomTypes=5 NoSymmetry Angstrom
        Charge=6.0   Atoms=1    Basis=STO-3G
C           -1.588367    -.770650     .029109 0 1
        Charge=8.0   Atoms=1    Basis=STO-3G
O           -1.657083     .436069    -.009750 0 2
        Charge=1.0   Atoms=2    Basis=STO-3G
H           -.620668   -1.294822      .054251 0 3
H           -2.508043   -1.382001     .040282 0 4
   Charge=-0.669     Atoms=2   Basis=MM
O             .975536    1.507219    -.082590 1 1
O            1.633705   -1.233280    -.078786 2 1
    Charge=0.3345    Atoms=4   Basis=MM
H             .023200    1.300096    -.088063 1 2
H            1.113831    2.040431     .704778 1 3
H            1.600739    -.258639    -.086953 2 2
H            2.303878   -1.461228     .570065 2 3
END MOLINP

########## POTENTIAL.INP ############################
START POTINP
**SYSTP
.NUMMTP
 1
.TYPE
 0
.MODEL
 SPC
.CHARGS
 4
 0.0000
 0.0000
 0.0000
 0.0000
*******
.TYPE
 1-2
.MODEL
 SPC_E01
.ALPISO
 1
 9.501
*******
**TWOIA (i,j=0,1,2,...,N; if i=0 then j.neq.0)
.LJ_A
 2
 0.00000
 0.00000
.LJ_B
 2
 0.000
 0.000
**END OF
END POTINP

########## Reference Output ########################
START REFOUT


   ****************************************************************************
   *************** DALTON2011 - An electronic structure program ***************
   ****************************************************************************

    This is output from DALTON Release 2011 (DEVELOPMENT VERSION)
   ----------------------------------------------------------------------------
    NOTE:
     
    DALTON is an experimental code for the evaluation of molecular
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
    Release Dalton2011 (2011), see http://daltonprogram.org"
   ----------------------------------------------------------------------------

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
  Heike Fliegl,             University of Helsinki,       Finland     (CCSD(R12))
  Luca Frediani,            University of Tromsoe,        Norway      (PCM)
  Bin Gao,                  University of Tromsoe,        Norway      (Gen1Int module)
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
  Juan I. Melo,             University of Buenos Aires,   Argentina   (LRESC, Relativistic Effects on NMR Shieldings)
  Kurt V. Mikkelsen,        University of Copenhagen,     Denmark     (MC-SCRF and QM/MM code)
  Christian Neiss,          Univ. Erlangen-Nuernberg,     Germany     (CCSD(R12))
  Christian B. Nielsen,     University of Copenhagen,     Denmark     (QM/MM code)
  Patrick Norman,           University of Linkoeping,     Sweden      (cubic response and complex response in RESPONS)
  Jeppe Olsen,              Aarhus University,            Denmark     (SIRIUS CI/density modules)
  Anders Osted,             Copenhagen University,        Denmark     (QM/MM code)
  Martin J. Packer,         University of Sheffield,      UK          (SOPPA)
  Filip Pawlowski,          Kazimierz Wielki University   Poland      (CC3)
  Thomas B. Pedersen,       University of Oslo,           Norway      (Cholesky decomposition)
  Patricio F. Provasi,      University of Northeastern,   Argentina   (Analysis of coupling constants in localized orbitals)
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
  David P. Tew,             University of Bristol,        England     (CCSD(R12))
  Olav Vahtras,             KTH Stockholm,                Sweden      (triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)
  David J. Wilson,          La Trobe University,          Australia   (DFT Hessian and DFT magnetizabilities)
  Hans Agren,               KTH Stockholm,                Sweden      (SIRIUS module, MC-SCRF solvation model)
 --------------------------------------------------------------------------------

     Date and time (Linux)  : Fri Mar 22 13:48:42 2013
     Host name              : compute-1-7.local                       

 * Work memory size             :    64000000 =  488.28 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/calculations/test/qm3pcm-merge
   2) /home/arnfinn/dalton/build-merge/basis


   Content of the .dal input file
 ----------------------------------

**DALTON                                          
.RUN RESPONSE                                     
*PCM                                              
.NEQRSP                                           
.SOLVNT                                           
WATER                                             
.NPCMMT                                           
0                                                 
.NESFP                                            
8                                                 
.ICESPH                                           
2                                                 
*PCMCAV                                           
.RIN                                              
1.7                                               
1.5                                               
1.2                                               
1.2                                               
1.5                                               
1.5                                               
1.2                                               
1.2                                               
.INA                                              
1                                                 
2                                                 
3                                                 
4                                                 
5                                                 
6                                                 
7                                                 
8                                                 
.AREATS                                           
 0.3                                              
.ALPHA                                            
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
 1.2                                              
*QM3                                              
.PRINT                                            
 1                                                
.QM3                                              
.MMPCM                                            
.THRSMP                                           
 1.0D-9                                           
**INTEGRALS                                       
.DIPLEN                                           
.NUCPOT                                           
.NELFLD                                           
.PRINT                                            
 2                                                
**WAVE FUNCTIONS                                  
.DFT                                              
B3LYP                                             
*SCF INPUT                                        
.THRESH                                           
1.0D-12                                           
**RESPONSE                                        
*LINEAR                                           
.DIPLEN                                           
.SINGLE RES                                       
.ROOTS                                            
4                                                 
.THCPP                                            
 1.0D-04                                          
.PRINT                                            
5                                                 
**END OF DALTON INPUT                             


   Content of the .mol file
 ----------------------------

ATOMBASIS                                                                      
CH2O PLUS 2 WATERS                                                             
------------------------                                                       
AtomTypes=5 NoSymmetry Angstrom                                                
        Charge=6.0   Atoms=1    Basis=STO-3G                                   
C           -1.588367    -.770650     .029109 0 1                              
        Charge=8.0   Atoms=1    Basis=STO-3G                                   
O           -1.657083     .436069    -.009750 0 2                              
        Charge=1.0   Atoms=2    Basis=STO-3G                                   
H           -.620668   -1.294822      .054251 0 3                              
H           -2.508043   -1.382001     .040282 0 4                              
   Charge=-0.669     Atoms=2   Basis=MM                                        
O             .975536    1.507219    -.082590 1 1                              
O            1.633705   -1.233280    -.078786 2 1                              
    Charge=0.3345    Atoms=4   Basis=MM                                        
H             .023200    1.300096    -.088063 1 2                              
H            1.113831    2.040431     .704778 1 3                              
H            1.600739    -.258639    -.086953 2 2                              
H            2.303878   -1.461228     .570065 2 3                              


       *******************************************************************
       *********** Output from DALTON general input processing ***********
       *******************************************************************

 --------------------------------------------------------------------------------
   Overall default print level:    0
   Print level for DALTON.STAT:    1

    HERMIT 1- and 2-electron integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed (SIRIUS module)
    Dynamic molecular response properties section will be executed (RESPONSE module)
 --------------------------------------------------------------------------------

 ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT = WATER    **
 Optical and physical constants:
 EPS= 78.390; EPSINF=  1.776; RSOLV=  1.385 A; VMOL=  18.070 ML/MOL;
 TCE=2.57000D-04 1/K; STEN= 71.810 DYN/CM;  DSTEN=  0.6500; CMF=  1.2770


     -----------------------------------
     Input for PCM solvation calculation 
     -----------------------------------
     ICOMPCM =       0          SOLVNT=WATER        EPS   = 78.3900     EPSINF=  1.7760
     RSOLV =  1.3850

     ICESPH =       2     NESFP =       8
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
     5    0.0000    0.0000    0.0000    1.5000
     6    0.0000    0.0000    0.0000    1.5000
     7    0.0000    0.0000    0.0000    1.2000
     8    0.0000    0.0000    0.0000    1.2000


 Changes of defaults for *QM3   :
 --------------------------------

 +------------------+
 |  WORD: | CHANGE: |
 +------------------+
 |    QM3 |       T |
 |  PRINT |       1 |
 +------------------+
  Settings for determination of induced dipoles:
  Iterative Method is used
  MM/PCM interface active
 +------------------+


    *************************************************************************
    *************** Output from MM potential input processing ***************
    *************************************************************************



 ------------------------------------------------------------------------
 | QM-sys type: | Systems:    | Model:  | Electric properties:          |
 ------------------------------------------------------------------------
 |          0   | [   0;   0] | SPC     |     No classical polarization |
 ------------------------------------------------------------------------


 ------------------------------------------------------------------------
 | MM-sys type: | Systems:    | Model:  | Electric properties:          |
 ------------------------------------------------------------------------
 |          1   | [   1;   2] | SPC_E01 | Isotropic polarizability:     |
 |              |             |         | alp( 1)=  9.5010              |
 ------------------------------------------------------------------------




   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************


    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1: CH2O PLUS 2 WATERS                                                      
 2: ------------------------                                                
    ------------------------------------------------------------------------

  Coordinates are entered in Angstrom and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   6 :
     "/home/arnfinn/dalton/build-merge/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   8 :
     "/home/arnfinn/dalton/build-merge/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  The basis set is "STO-3G" from the basis set library.
  Basis set file used for this atomic type with Z =   1 :
     "/home/arnfinn/dalton/build-merge/basis/STO-3G"

  Atomic type no.    4
  --------------------
  Nuclear charge:  -0.66900
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  This is an MM atom without basis functions.

  Atomic type no.    5
  --------------------
  Nuclear charge:   0.33450
  Number of symmetry independent centers:    4
  Number of basis sets to read;    2
  This is an MM atom without basis functions.


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1   -3.0015786160D+00   -1.4563174382D+00    5.5008037772D-02    1.7000000000D+00
   2   -3.1314330364D+00    8.2405098160D-01   -1.8424829719D-02    1.5000000000D+00
   3   -1.1728925345D+00   -2.4468589606D+00    1.0251953201D-01    1.2000000000D+00
   4   -4.7395143797D+00   -2.6116033945D+00    7.6121947767D-02    1.2000000000D+00
   5    1.8434958651D+00    2.8482311204D+00   -1.5607248066D-01    1.5000000000D+00
   6    3.0872550190D+00   -2.3305614354D+00   -1.4888396248D-01    1.5000000000D+00
   7    4.3841646100D-02    2.4568253762D+00   -1.6641495175D-01    1.2000000000D+00
   8    2.1048355395D+00    3.8558557669D+00    1.3318373989D+00    1.2000000000D+00


                                 Isotopic Masses
                                 ---------------

                           C          12.000000
                           O          15.994915
                           H           1.007825
                           H           1.007825

                       Total mass:    30.010565 amu
                       Natural abundance:  98.633 %

 Center-of-mass coordinates (a.u.):   -3.067740   -0.312997    0.018175


  Atoms and basis sets
  --------------------

  Number of atom types :    6
  Total number of atoms:   12

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C           1    6.0000    15     5      [6s3p|2s1p]                                        
  O           1    8.0000    15     5      [6s3p|2s1p]                                        
  H           2    1.0000     3     1      [3s|1s]                                            
  O           2   -0.6690     0     0      Point Charge                                       
  H           4    0.3345     0     0      Point Charge                                       
  a           2    0.0000     0     0      Point Charge                                       
  ----------------------------------------------------------------------
  total:     12   16.0000    36    12
  ----------------------------------------------------------------------

  Threshold for neglecting AO integrals:  1.00D-12


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   36
  C       :     1  x  -3.0015786160    2  y  -1.4563174382    3  z   0.0550080378
  O       :     4  x  -3.1314330364    5  y   0.8240509816    6  z  -0.0184248297
  H       :     7  x  -1.1728925345    8  y  -2.4468589606    9  z   0.1025195320
  H       :    10  x  -4.7395143797   11  y  -2.6116033945   12  z   0.0761219478


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           O           H           H           O           O     
            ------      ------      ------      ------      ------      ------
 C     :    0.000000
 O     :    1.209298    0.000000
 H     :    1.100831    2.018474    0.000000
 H     :    1.104391    2.007988    1.889439    0.000000
 O     :    3.431437    2.843123    3.227697    4.527473    0.000000
 O     :    3.256903    3.690634    2.259133    4.146127    2.818428    0.000000
 H     :    2.626569    1.891038    2.677390    3.690164    0.974615    3.001967
 H     :    3.957345    3.280624    3.815178    5.027183    0.960931    3.406086
 H     :    3.232031    3.331964    2.455252    4.261481    1.873273    0.975233
 H     :    3.989875    4.430026    2.974344    4.841645    3.316946    0.960260
 a  1  :    3.408645    2.807232    3.219816    4.503294    0.065744    2.847367
 a  2  :    3.285589    3.703418    2.294418    4.182574    2.786761    0.065593

            H           H           H           H           a  1        a  2  
            ------      ------      ------      ------      ------      ------
 H     :    0.000000
 H     :    1.538236    0.000000
 H     :    2.217721    2.479847    0.000000
 H     :    3.641364    3.700807    1.540226    0.000000
 a  1  :    0.935667    0.923012    1.906635    3.343454    0.000000
 a  2  :    2.986623    3.363487    0.936414    0.922570    2.815774    0.000000


  Max    interatomic separation is    2.0185 Angstrom (    3.8144 Bohr)
  between atoms    3 and    2, "H     " and "O     ".

  Min HX interatomic separation is    1.1008 Angstrom (    2.0803 Bohr)

  Min YX interatomic separation is    1.2093 Angstrom (    2.2852 Bohr)

  Max QM+MM interatomic separation is    5.0272 Angstrom (    9.5000 Bohr)
  between the QM+MM centers    8 and    4, "H     " and "H     ".




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       1.798871         -0.056827    0.997867   -0.032124
   IB      13.009178          0.998345    0.057081    0.007042
   IC      14.808049         -0.008861    0.031671    0.999459


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         280942.3015          38847.8806          34128.6681 MHz
            9.371226            1.295826            1.138410 cm-1


@  Nuclear repulsion energy :   31.249215315972 Hartree
      QM3 induced dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.020409711860107


    **************************************************************************
    ***************** The MM/MM classical interaction energy *****************
    **************************************************************************


 ------------------------------------------------------------------------
 |  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        |       -0.00626437 |
 |  Epol  = - 1/2*Sum_a[ Pind_a*E^site_a ]          |       -0.00096177 |
 |  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] |        0.00000000 |
 ------------------------------------------------------------------------
 |  E(MM/MM) = Eelec + Epol + Evdw                  |       -0.00722614 |
 ------------------------------------------------------------------------




     ************************************************************************
     *************** The "QM"/MM classical interaction energy ***************
     ************************************************************************


 ------------------------------------------------------------------------
 |  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        |        0.00000000 |
 |  Epol  = - 1/2*Sum_a[ Pind_a*E^(QMclassic)_a ]   |        0.00000000 |
 |  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] |        0.00000000 |
 ------------------------------------------------------------------------
 |  E("QM"/MM) = Eelec + Epol + Evdw                |        0.00000000 |
 ------------------------------------------------------------------------




                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************


 Default print level:        2

 * Nuclear model: Point charge

 Calculation of one- and two-electron Hamiltonian integrals.

 The following one-electron property integrals are calculated as requested:
          - overlap integrals
          - dipole length integrals
          - Electric field at the nuclei
          - Potential energy at the nuclei

 Center of mass  (bohr):     -3.067740315293     -0.312997395263      0.018174673399
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERONE **************************
     ************************************************************************



     ************************************************************************
     ************************** Output from TWOINT **************************
     ************************************************************************


 Threshold for neglecting two-electron integrals:  1.00D-12
 Number of two-electron integrals written:        2907 ( 94.4% )
 Megabytes written:                              0.034


 MEMORY USED TO GENERATE CAVITY =    432042

 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.
 Tessera cut in pieces and removed.

 Total number of spheres =   11
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1   -1.588367000   -0.770650000    0.029109000    2.040000000   22.843745166
   2   -1.657083000    0.436069000   -0.009750000    1.800000000   17.320494799
   3   -0.620668000   -1.294822000    0.054251000    1.440000000    4.504430717
   4   -2.508043000   -1.382001000    0.040282000    1.440000000    9.471926497
   5    0.975536000    1.507219000   -0.082590000    1.800000000   16.946491810
   6    1.633705000   -1.233280000   -0.078786000    1.800000000   28.188763663
   7    0.023200000    1.300096000   -0.088063000    1.440000000    4.066167360
   8    1.113831000    2.040431000    0.704778000    1.440000000   11.882478847
   9    1.304620500    0.136969500   -0.080688000    1.471280988    4.271467524
  10   -0.012749934    0.270009434   -0.049134665    1.232658185    0.798536321
  11    0.329445664   -0.588244950   -0.012340260    1.188874163    0.521620751

 Total number of tesserae =     579
 Surface area =  120.81612345 (A^2)    Cavity volume =   96.50955728 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....

  ..... DONE GENERATING -Q-  MATRIX .....
 >>>  Time used in Q-MAT      is   2.18 seconds
 >>>> Total CPU  time used in HERMIT:   4.44 seconds
 >>>> Total wall time used in HERMIT:   4.67 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          F
     Using EHT  model:          T
     Number of Huckel orbitals each symmetry:   12

 Huckel EHT eigenvalues for symmetry :  1
          -20.705291     -11.377081      -1.496925      -0.962689      -0.743251
           -0.672270      -0.593413      -0.504635      -0.350830      -0.253214
           -0.198302      -0.193199

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Fri Mar 22 13:48:47 2013
     Host name              : compute-1-7.local                       

 Title lines from ".mol" input file:
     CH2O PLUS 2 WATERS                                                      
     ------------------------                                                

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

@    Restricted, closed shell Kohn-Sham DFT calculation.

@    Time-dependent Kohn-Sham DFT calculation (TD-DFT).

 Initial molecular orbitals are obtained according to
 ".MOSTART EHT   " input option

@    QM part is embedded in an environment :

@         Solvation model: PCM

@         Model: QM3

     Wave function specification
     ============================
     For the wave function of type :      >>> KS-DFT   <<<
     Number of closed shell electrons          16
     Number of electrons in active shells       0
     Total charge of the molecule               0

     Spin multiplicity and 2 M_S                1         0
     Total number of symmetries                 1
     Reference state symmetry                   1
 
     This is a DFT calculation of type: B3LYP
 Weighted mixed functional:
               HF exchange:    0.20000
                       VWN:    0.19000
                       LYP:    0.81000
                     Becke:    0.72000
                    Slater:    0.80000

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1
                                       --- |  ---
     Occupied SCF orbitals               8 |    8
     Secondary orbitals                  4 |    4
     Total number of orbitals           12 |   12
     Number of basis functions          12 |   12

     Optimization information
     ========================
@    Number of configurations                 1
@    Number of orbital rotations             32
     ------------------------------------------
@    Total number of variables               33

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00D-12
 
     This is a DFT calculation of type: B3LYP
 Weighted mixed functional:
               HF exchange:    0.20000
                       VWN:    0.19000
                       LYP:    0.81000
                     Becke:    0.72000
                    Slater:    0.80000

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.99085
 NUCLEAR APPARENT CHARGE -15.79248
 THEORETICAL -15.79589 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.975378044871  16.0000012547    1.25D-06
      QM3 induced Dipole vector converged in  12 iterations.
      Final norm2 of QM3 induced dipole moment vector:   13.241662292657812
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.646363105050853
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.081849518869684
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.047171361930943
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043718750816638
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043194613844121
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043096838589755
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043077106815990
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043072989862787
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043072115514704
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043071927802125
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043071887143437
 MMPCM converged for given density
 No. of iterations:           12
   1  -112.639529351     -3.423203955690E-02    2.41139D+00   -1.13D+02
      Virial theorem: -V/T =      2.004570
@      MULPOP C       0.99; O      -0.72; H      -0.14; H      -0.14; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.765242249359  16.0000013790    1.38D-06
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.009190575878454
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.020963873629841
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.026331184411254
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.027317183411319
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.027485685307953
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.027514461010933
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.027519442128667
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.027520322241910
 MMPCM converged for given density
 No. of iterations:            8
   2  -112.417759747     -5.284200264228E-03    2.53580D+00    2.22D-01
      Virial theorem: -V/T =      2.014727
@      MULPOP C      -1.38; O       0.94; H       0.24; H       0.20; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.950609959362  16.0000013055    1.31D-06
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.125973962747581
      QM3 induced Dipole vector converged in  10 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.053369399591084
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.044634618617630
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043303188574267
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043089551490429
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043054749332645
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043049014493422
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043048055618050
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.043047892098430
 MMPCM converged for given density
 No. of iterations:            9
   3  -112.848030512     -2.777671776174E-02    8.25875D-01   -4.30D-01
      Virial theorem: -V/T =      2.008969
@      MULPOP C       0.45; O      -0.61; H       0.10; H       0.07; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.870212627640  16.0000013139    1.31D-06
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.021917464306991
      QM3 induced Dipole vector converged in   9 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.033172454108083
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035226928510843
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035566655085391
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035622324889660
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035631501780788
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035633031551232
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035633290675340
 MMPCM converged for given density
 No. of iterations:            8
   4  -112.902375633     -1.008327488218E-02    3.76011D-02   -5.43D-02
      Virial theorem: -V/T =      2.013820
@      MULPOP C       0.03; O      -0.19; H       0.10; H       0.07; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875599499030  16.0000013158    1.32D-06
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.037080243918050
      QM3 induced Dipole vector converged in   8 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036281467570861
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036153106161559
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036132253641701
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036128839641733
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036128275353741
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036128180771217
 MMPCM converged for given density
 No. of iterations:            7
   5  -112.902499831     -1.094028162610E-02    7.11863D-03   -1.24D-04
      Virial theorem: -V/T =      2.013440
@      MULPOP C       0.02; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.874943770423  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.035965772495676
      QM3 induced Dipole vector converged in   7 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036051028038225
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036064839333231
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036067087558026
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036067456150712
 MMPCM converged for given density
 No. of iterations:            5
   6  -112.902504929     -1.084505236531E-02    1.79379D-03   -5.10D-06
      Virial theorem: -V/T =      2.013493
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875051359319  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036100430567094
      QM3 induced Dipole vector converged in   6 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036080585004269
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036077380428492
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036076860546265
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036076775716642
 MMPCM converged for given density
 No. of iterations:            5
   7  -112.902505055     -1.086611026891E-02    2.12941D-04   -1.27D-07
      Virial theorem: -V/T =      2.013488
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875065166005  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036081133021128
      QM3 induced Dipole vector converged in   5 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078524112917
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078102786485
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078034451081
 MMPCM converged for given density
 No. of iterations:            4
   8  -112.902505046     -1.086886214315E-02    1.58447D-05    9.91D-09
      Virial theorem: -V/T =      2.013487
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875066226758  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078340113571
      QM3 induced Dipole vector converged in   4 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078147650972
 MMPCM converged for given density
 No. of iterations:            2
   9  -112.902505045     -1.086906563986E-02    1.42117D-06    8.02D-10
      Virial theorem: -V/T =      2.013487
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875066345980  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078154509796
      QM3 induced Dipole vector converged in   3 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078126957434
 MMPCM converged for given density
 No. of iterations:            2
  10  -112.902505045     -1.086908938105E-02    1.93626D-07    2.12D-10
      Virial theorem: -V/T =      2.013487
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875066332718  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078118350127
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120037223
 MMPCM converged for given density
 No. of iterations:            2
  11  -112.902505045     -1.086908674624E-02    5.68107D-09    4.15D-12
      Virial theorem: -V/T =      2.013487
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875066332111  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120143327
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120299373
 MMPCM converged for given density
 No. of iterations:            2
  12  -112.902505045     -1.086908663503E-02    1.65347D-09   -1.29D-12
      Virial theorem: -V/T =      2.013487
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875066332227  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120362510
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120346061
 MMPCM converged for given density
 No. of iterations:            2
  13  -112.902505045     -1.086908665784E-02    3.20863D-11    0.00D+00
      Virial theorem: -V/T =      2.013487
@      MULPOP C       0.03; O      -0.21; H       0.11; H       0.08; 
 -----------------------------------------------------------------------------
      K-S energy, electrons, error :    -11.875066332229  16.0000013153    1.32D-06
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120343982
      QM3 induced Dipole vector converged in   2 iterations.
      Final norm2 of QM3 induced dipole moment vector:    0.036078120343088
 MMPCM converged for given density
 No. of iterations:            2
  14  -112.902505045     -1.086908665835E-02    3.47038D-13    8.53D-14

 *** DIIS converged in  14 iterations !
@    Converged SCF energy, gradient:   -112.902505044558    3.47D-13
   - total time used in SIRFCK :              0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Sym       Kohn-Sham orbital energies

  1    -18.91321918   -10.08163792    -1.01202420    -0.57188137    -0.44052230
        -0.35083954    -0.31959556    -0.16395483     0.07209758     0.40032444
         0.50229306     0.61656080

    E(LUMO) :     0.07209758 au (symmetry 1)
  - E(HOMO) :    -0.16395483 au (symmetry 1)
  ------------------------------------------
    gap     :     0.23605241 au

 >>> Writing SIRIFC interface file <<<

 >>>> CPU and wall time for SCF :      49.714      50.131


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


@    Spin multiplicity:           1
@    Spatial symmetry:            1
@    Total charge of molecule:    0

     QM/MM/PCM calculation converged :

        Electrostatic energy:      -0.008289821636
        Polarization energy:       -0.000525528099
        van der Waals energy:       0.000000000000
        QM/MM/PCM energy:          -0.010869086658

@    Final DFT energy:           -112.902505044558                 
@    Nuclear repulsion:            31.249215315972
@    Electronic energy:          -144.140851273872

@    Final gradient norm:           0.000000000000

 
     Date and time (Linux)  : Fri Mar 22 13:49:37 2013
     Host name              : compute-1-7.local                       

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

    Orbital         1        2        3        4        5        6        7
   1 C   :1s     0.0007   0.9910  -0.1284  -0.1877  -0.0050   0.0373   0.0001
   2 C   :1s    -0.0089   0.0406   0.2841   0.5911   0.0188  -0.1166   0.0004
   3 C   :2px    0.0005  -0.0001  -0.0095   0.0044   0.5774   0.0157  -0.0069
   4 C   :2py   -0.0078   0.0011   0.1692  -0.2576   0.0180  -0.4248   0.0232
   5 C   :2pz    0.0002  -0.0000  -0.0045   0.0109   0.0066   0.0184   0.5959
   6 O   :1s     0.9935   0.0005  -0.2153   0.0989   0.0015  -0.1101   0.0013
   7 O   :1s     0.0297  -0.0082   0.7406  -0.4198  -0.0048   0.5564  -0.0061
   8 O   :2px    0.0004  -0.0002   0.0176  -0.0038   0.3368  -0.0250  -0.0067
   9 O   :2py   -0.0064   0.0031  -0.2095  -0.1117   0.0231   0.6816   0.0163
  10 O   :2pz    0.0002  -0.0001   0.0080   0.0033   0.0037  -0.0147   0.6859
  11 H   :1s     0.0003  -0.0088   0.0306   0.2514   0.3133   0.1318  -0.0030
  12 H   :1s     0.0003  -0.0088   0.0315   0.2691  -0.3109   0.1497  -0.0014

    Orbital         8        9       10
   1 C   :1s     0.0009   0.0000   0.2048
   2 C   :1s    -0.0052   0.0009  -1.2680
   3 C   :2px   -0.0984  -0.0076   0.0530
   4 C   :2py    0.0046   0.0269   0.4609
   5 C   :2pz   -0.0020   0.8309  -0.0128
   6 O   :1s     0.0034  -0.0002  -0.0229
   7 O   :1s    -0.0158   0.0005   0.1306
   8 O   :2px    0.9045   0.0075  -0.0133
   9 O   :2py    0.0408  -0.0221  -0.2126
  10 O   :2pz    0.0075  -0.7587   0.0060
  11 H   :1s    -0.3391  -0.0021   0.8481
  12 H   :1s     0.3479  -0.0007   0.9404



 >>>> Total CPU  time used in SIRIUS :     49.72 seconds
 >>>> Total wall time used in SIRIUS :     50.20 seconds

 
     Date and time (Linux)  : Fri Mar 22 13:49:37 2013
     Host name              : compute-1-7.local                       

 NOTE:    1 warnings have been issued.
 Check output, result, and error files for "WARNING".


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

 QM/MM response calculation
 Non-equilibrium PCM solvent model requested    : INERSI =T

 Static dielectric constant                     : EPSTAT = 78.3900
 Optical dielectric constant                    : EPSOL  =  1.7760
 Print level                                    : IPRPP  =   5
 Maximum number of iterations for eigenval.eqs. : MAXITP =  60
 Threshold for convergence of eigenvalue eqs.   : THCPP  = 1.000D-04
 Maximum iterations in optimal orbital algorithm: MAXITO =   5

      4 Excitation energies are calculated for symmetry no.    1

      3 property residues are calculated with labels:

               XDIPLEN 
               YDIPLEN 
               ZDIPLEN 

 Integral transformation: Total CPU and WALL times (sec)       0.002       0.118


   SCF energy         :     -112.902505044558382
 -- inactive part     :     -144.282485409391796
 -- nuclear repulsion :       31.249215315972382


                    *****************************************
                    *** DFT response calculation (TD-DFT) ***
                    *****************************************



 >>>>>>>>>> Linear response calculation
 >>>>>>>>>> Symmetry of excitation/property operator(s)    1

 Number of excitations of this symmetry            4
 Number of response properties of this symmetry    0
 Number of C6/C8 properties of this symmetry       0


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      32
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      32
 Electrons in DFTMOMO:   16.00000131525930


                         Matrix < abs(phi_p) abs(phi_q) >
                         --------------------------------


               Column   1     Column   2     Column   3     Column   4     Column   5     Column   6     Column   7     Column   8
       1       1.00000001     0.00199754     0.21078672     0.09829678     0.07893030     0.19221261     0.15190092     0.19131494
       2       0.00199754     1.00000055     0.13378850     0.19067920     0.14964670     0.11823613     0.14468472     0.03398147
       3       0.21078672     0.13378850     1.00000005     0.47119193     0.60888887     0.70949695     0.79405162     0.69606600
       4       0.09829678     0.19067920     0.47119193     1.00000004     0.82451440     0.66055338     0.59950171     0.62382659
       5       0.07893030     0.14964670     0.60888887     0.82451440     1.00000005     0.55552499     0.52976099     0.74566021
       6       0.19221261     0.11823613     0.70949695     0.66055338     0.55552499     1.00000000     0.62080319     0.64825786
       7       0.15190092     0.14468472     0.79405162     0.59950171     0.52976099     0.62080319     0.99999999     0.51960690
       8       0.19131494     0.03398147     0.69606600     0.62382659     0.74566021     0.64825786     0.51960690     0.99999998
       9       0.15654009     0.18265097     0.66635882     0.68291771     0.53601999     0.59117308     0.94591192     0.48968138
      10       0.05090694     0.21527034     0.32755277     0.82573712     0.71542753     0.45531276     0.41720526     0.49622792
      11       0.07061320     0.22270398     0.47259150     0.72697751     0.87373386     0.50023855     0.43484730     0.60665266
      12       0.21145389     0.24693907     0.76789179     0.63440630     0.61812457     0.71156088     0.67161452     0.51400533

               Column   9     Column  10     Column  11     Column  12
       1       0.15654009     0.05090694     0.07061320     0.21145389
       2       0.18265097     0.21527034     0.22270398     0.24693907
       3       0.66635882     0.32755277     0.47259150     0.76789179
       4       0.68291771     0.82573712     0.72697751     0.63440630
       5       0.53601999     0.71542753     0.87373386     0.61812457
       6       0.59117308     0.45531276     0.50023855     0.71156088
       7       0.94591192     0.41720526     0.43484730     0.67161452
       8       0.48968138     0.49622792     0.60665266     0.51400533
       9       1.00000000     0.48779843     0.46839400     0.66073661
      10       0.48779843     1.00000006     0.79089302     0.51715303
      11       0.46839400     0.79089302     1.00000008     0.55472391
      12       0.66073661     0.51715303     0.55472391     1.00000001
    ==== End of matrix output ====


                            Matrix < phi_p^2 phi_q^2 >
                            --------------------------


               Column   1     Column   2     Column   3     Column   4     Column   5     Column   6     Column   7     Column   8
       1       1.00000000     0.00000054     0.42677038     0.14758752     0.01468897     0.13747232     0.04252137     0.05262530
       2       0.00000054     1.00000000     0.09926033     0.34093778     0.03400688     0.01785736     0.02392719     0.00067119
       3       0.42677038     0.09926033     1.00000000     0.22573340     0.28286036     0.41838901     0.52147074     0.46522767
       4       0.14758752     0.34093778     0.22573340     1.00000000     0.62266090     0.34964153     0.29592476     0.24519516
       5       0.01468897     0.03400688     0.28286036     0.62266090     1.00000000     0.19366421     0.21778583     0.43560640
       6       0.13747232     0.01785736     0.41838901     0.34964153     0.19366421     1.00000000     0.33697878     0.35071282
       7       0.04252137     0.02392719     0.52147074     0.29592476     0.21778583     0.33697878     1.00000000     0.29105594
       8       0.05262530     0.00067119     0.46522767     0.24519516     0.43560640     0.35071282     0.29105594     1.00000000
       9       0.04406058     0.03857368     0.38771453     0.37636880     0.23848713     0.31226289     0.93825597     0.25977643
      10       0.00971941     0.32059074     0.10263621     0.74058886     0.48948561     0.11796389     0.16795551     0.14638750
      11       0.00963643     0.06729712     0.15593373     0.44292304     0.79848101     0.13662374     0.16701276     0.23747054
      12       0.16426568     0.10598164     0.66685130     0.31377589     0.29934120     0.44254076     0.39590536     0.30017862

               Column   9     Column  10     Column  11     Column  12
       1       0.04406058     0.00971941     0.00963643     0.16426568
       2       0.03857368     0.32059074     0.06729712     0.10598164
       3       0.38771453     0.10263621     0.15593373     0.66685130
       4       0.37636880     0.74058886     0.44292304     0.31377589
       5       0.23848713     0.48948561     0.79848101     0.29934120
       6       0.31226289     0.11796389     0.13662374     0.44254076
       7       0.93825597     0.16795551     0.16701276     0.39590536
       8       0.25977643     0.14638750     0.23747054     0.30017862
       9       1.00000000     0.24237833     0.21052062     0.36751950
      10       0.24237833     1.00000000     0.57465093     0.24019157
      11       0.21052062     0.57465093     1.00000000     0.26961238
      12       0.36751950     0.24019157     0.26961238     1.00000000
    ==== End of matrix output ====

 >>> IN RSPPP:
  THCPP,MAXRM   1.00000000000000005E-004         600
  KSYMOP,NGPPP(KSYMOP)            1           3
  LWRK ,LWRK1     63999526    62917438
  KEXCNV,NSIM,LWRK2            4           4     1082605



 <<< EXCITATION ENERGIES AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>

 Operator symmetry =  1; triplet =   F

           0  START configuration VECTORS USING LOWEST DIAGONAL HESSIAN ELEMENTS
           4  START orbital VECTORS

 ** RSPCTL MICROITERATION NUMBER    1
       Electrons: 16.000001( 1.32e-06): LR-DFT*4 evaluation time:       1.2 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.38821D-02  0.00D+00  1.39D-02  7.09D-01    1.51503D-01
         2    1.02923D-01  0.00D+00  1.03D-01  7.09D-01    3.57105D-01
         3    9.91618D-03  0.00D+00  9.92D-03  7.07D-01    4.34320D-01
         4    4.69505D-01  0.00D+00  4.70D-01  7.54D-01    5.32177D-01

 ** RSPCTL MICROITERATION NUMBER    2
       Electrons: 16.000001( 1.32e-06): LR-DFT*4 evaluation time:       1.2 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    2.11175D-04  0.00D+00  2.11D-04  7.09D-01    1.51393D-01
         2    4.95820D-03  0.00D+00  4.96D-03  7.09D-01    3.51437D-01
         3    1.88025D-03  0.00D+00  1.88D-03  7.07D-01    4.34183D-01
         4    4.75703D-02  0.00D+00  4.76D-02  7.28D-01    4.39532D-01

 ** RSPCTL MICROITERATION NUMBER    3
       Electrons: 16.000001( 1.32e-06): LR-DFT*4 evaluation time:       1.2 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.47646D-05  0.00D+00  1.48D-05  7.09D-01    1.51393D-01
         2    1.47040D-04  0.00D+00  1.47D-04  7.09D-01    3.51435D-01
         3    8.22088D-04  0.00D+00  8.22D-04  7.07D-01    4.34180D-01
         4    1.69279D-02  0.00D+00  1.69D-02  7.29D-01    4.37672D-01

 ** RSPCTL MICROITERATION NUMBER    4
       Electrons: 16.000001( 1.32e-06): LR-DFT*3 evaluation time:       0.9 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.30591D-05  0.00D+00  1.31D-05  7.09D-01    1.51393D-01
         2    9.49479D-06  0.00D+00  9.49D-06  7.09D-01    3.51435D-01
         3    6.73540D-05  0.00D+00  6.74D-05  7.07D-01    4.34179D-01
         4    1.31555D-03  0.00D+00  1.32D-03  7.29D-01    4.37537D-01

 ** RSPCTL MICROITERATION NUMBER    5
       Electrons: 16.000001( 1.32e-06): LR-DFT*1 evaluation time:       0.4 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.30516D-05  0.00D+00  1.31D-05  7.09D-01    1.51393D-01
         2    8.29402D-06  0.00D+00  8.29D-06  7.09D-01    3.51435D-01
         3    6.20804D-06  0.00D+00  6.21D-06  7.07D-01    4.34179D-01
         4    8.59199D-05  0.00D+00  8.59D-05  7.29D-01    4.37536D-01

 ** RSPCTL MICROITERATION NUMBER    6
       Electrons: 16.000001( 1.32e-06): LR-DFT*1 evaluation time:       0.4 s

      Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue
      ----------------------------------------------------------------
         1    1.30901D-05  0.00D+00  1.31D-05  7.09D-01    1.51393D-01
         2    8.39845D-06  0.00D+00  8.40D-06  7.09D-01    3.51435D-01
         3    5.36640D-06  0.00D+00  5.37D-06  7.07D-01    4.34179D-01
         4    4.19787D-05  0.00D+00  4.20D-05  7.29D-01    4.37536D-01

 *** THE REQUESTED    4 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   34)
 RSP solution vector no.    1; norm of residual   1.85D-05
 RSP solution vector no.    2; norm of residual   1.18D-05
 RSP solution vector no.    3; norm of residual   7.59D-06
 RSP solution vector no.    4; norm of residual   5.76D-05

 *** RSPCTL MICROITERATIONS CONVERGED

@ Transition operator type:    XDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -2.04706688E-03 *ENERGY(eV):   4.1196132    
@ STATE NO:    2 *TRANSITION MOMENT: -1.26730918E-03 *ENERGY(eV):   9.5630241    
@ STATE NO:    3 *TRANSITION MOMENT:  4.29012700E-03 *ENERGY(eV):   11.814620    
@ STATE NO:    4 *TRANSITION MOMENT: -5.24060013E-02 *ENERGY(eV):   11.905969    

@ Transition operator type:    YDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT:  1.19097343E-05 *ENERGY(eV):   4.1196132    
@ STATE NO:    2 *TRANSITION MOMENT:  3.02905036E-03 *ENERGY(eV):   9.5630241    
@ STATE NO:    3 *TRANSITION MOMENT: -3.71519975E-02 *ENERGY(eV):   11.814620    
@ STATE NO:    4 *TRANSITION MOMENT:  0.67670531     *ENERGY(eV):   11.905969    

@ Transition operator type:    ZDIPLEN 
@ STATE NO:    1 *TRANSITION MOMENT: -3.49353378E-03 *ENERGY(eV):   4.1196132    
@ STATE NO:    2 *TRANSITION MOMENT:  0.16972373     *ENERGY(eV):   9.5630241    
@ STATE NO:    3 *TRANSITION MOMENT: -1.12449231E-02 *ENERGY(eV):   11.814620    
@ STATE NO:    4 *TRANSITION MOMENT: -2.16047487E-02 *ENERGY(eV):   11.905969    


  ******************************************************************************
  *** @ Excit. operator sym 1 & ref. state sym 1 => excited state symmetry 1 ***
  ******************************************************************************



 @ Excited state no:    1 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.15139300    au
@                      4.1196132     eV
@                      33226.923     cm-1
                       397.48227     kJ / mol

@ Total energy :      -112.75111     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  4.22939854E-07  (Transition moment : -2.04706688E-03 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  1.43159011E-11  (Transition moment :  1.19097343E-05 )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  1.23181203E-06  (Transition moment : -3.49353378E-03 )

 Eigenvector for state no.  1

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          29         8(1)   9(1)       -0.7077302737       -0.0371243361

       31 elements with absolute value less than 7.08D-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      8    9 -0.707730 -0.037124  0.489681  0.259776  0.449712  0.220216

@ Overlap diagnostic LAMBDA =    0.4897


 @ Excited state no:    2 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.35143468    au
@                      9.5630241     eV
@                      77130.996     cm-1
                       922.69161     kJ / mol

@ Total energy :      -112.55107     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  3.76286399E-07  (Transition moment : -1.26730918E-03 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  2.14964302E-06  (Transition moment :  3.02905036E-03 )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  6.74898545E-03  (Transition moment :  0.16972373     )

 Eigenvector for state no.  2

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          21         6(1)   9(1)       -0.7056491034       -0.0310781332

       31 elements with absolute value less than 7.06D-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      6    9 -0.705649 -0.031078  0.591173  0.312263  0.455046  0.269011

@ Overlap diagnostic LAMBDA =    0.5924


 @ Excited state no:    3 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.43417933    au
@                      11.814620     eV
@                      95291.347     cm-1
                       1139.9376     kJ / mol

@ Total energy :      -112.46833     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  5.32743523E-06  (Transition moment :  4.29012700E-03 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  3.99523397E-04  (Transition moment : -3.71519975E-02 )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  3.66008235E-05  (Transition moment : -1.12449231E-02 )

 Eigenvector for state no.  3

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          17         5(1)   9(1)       -0.7055503354       -0.0080145980

       31 elements with absolute value less than 7.06D-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      5    9 -0.705550 -0.008015  0.536020  0.238487  0.486556  0.260804

@ Overlap diagnostic LAMBDA =    0.5368


 @ Excited state no:    4 in symmetry  1
 ---------------------------------------

@ Excitation energy :  0.43753633    au
@                      11.905969     eV
@                      96028.125     cm-1
                       1148.7515     kJ / mol

@ Total energy :      -112.46497     au

@ Operator type:    XDIPLEN 
@ Oscillator strength (LENGTH)   :  8.01096640E-04  (Transition moment : -5.24060013E-02 )

@ Operator type:    YDIPLEN 
@ Oscillator strength (LENGTH)   :  0.13357403      (Transition moment :  0.67670531     )

@ Operator type:    ZDIPLEN 
@ Oscillator strength (LENGTH)   :  1.36151146E-04  (Transition moment : -2.16047487E-02 )

 Eigenvector for state no.  4

     Orbital operator symmetry = 1
     (only elements abs greater than   10.00 % of max abs value)

      Index(r,s)      r      s        (r s) operator      (s r) operator
      ----------    -----  -----      --------------      --------------
          24         6(1)  12(1)       -0.1288844011        0.0427782882
          25         7(1)   9(1)       -0.6738965145       -0.1076449945
          31         8(1)  11(1)        0.1807040324       -0.0286019793

       29 elements with absolute value less than 6.74D-02 not printed.

 The numbers in parenthesis give the orbital symmetry.

     Configuration operator symmetry = 1
     >> NO ELEMENTS <<


                            PBHT MO Overlap Diagnostic
                            --------------------------

  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.
             J Chem Phys 128, 044118 (2008)


  The dominant contributions:

      I    A    K_IA      K_AI   <|I|*|A|> <I^2*A^2>    Weight   Contrib

      6   12 -0.128884  0.042778  0.711561  0.442541  0.029468  0.020968
      7    9 -0.673897 -0.107645  0.945912  0.938256  0.320641  0.303298
      8   11  0.180704 -0.028602  0.606653  0.237471  0.043809  0.026577

@ Overlap diagnostic LAMBDA =    0.8834


 Time used in polarization propagator calculation is     10.22 CPU seconds for symmetry 1

 >>>> Total CPU  time used in RESPONSE:  10.23 seconds
 >>>> Total wall time used in RESPONSE:  11.55 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  1 minute   4 seconds
 >>>> Total wall time used in DALTON:  1 minute   7 seconds

 
     Date and time (Linux)  : Fri Mar 22 13:49:49 2013
     Host name              : compute-1-7.local                       

END REFOUT

