########## Test description ########################
START DESCRIPTION
KEYWORDS: qm3 dft dipole properties
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enedft
qm3convergence
qm3ele
qm3pol
qm3wdv
qm3energy
dipole
dipcompx
dipcompy
dipcompz
excita
excita
OVERRIDE line 8
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON
.RUN PROPERTIES
*QM3
.PRINT
 1
.QM3
.THRDIP
 1.0D-12
.MAXDIP
 180
.OLDTG
**INTEGRALS
.DIPLEN
.NUCPOT
.NELFLD
.NOSUP
.PRINT
 2
**WAVE FUNCTIONS
.DFT
B3LYP
*SCF INPUT
.THRESH
1.0D-12
**PROPERTIES
.EXCITA
*EXCITA
.NEXCIT
 2
.THRESH
1.0D-9
**END OF
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS
CH2O PLUS 2 WATERS
------------------------
    5    0         1 1.00D-12
        6.0   1    Basis=STO-3G
C           -1.588367    -.770650     .029109 0 1
        8.0   1    Basis=STO-3G
O           -1.657083     .436069    -.009750 0 2
        1.0   2    Basis=STO-3G
H           -.620668   -1.294822      .054251 0 3
H           -2.508043   -1.382001     .040282 0 4
   -0.669     2   Basis=MM
O             .975536    1.507219    -.082590 1 1
O            1.633705   -1.233280    -.078786 2 1
    0.3345    4   Basis=MM
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

     Date and time (Linux)  : Sat Dec 11 20:39:48 2010
     Host name              : stallo-1.local                          

 * Work memory size             :   100000000 =  762.94 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/program/dalton/svn/pure_trunk/test/perl-pid.13663__2010_12_11__20.39
   2) /home/arnfinn/program/dalton/svn/pure_trunk/basis/


       *******************************************************************
       *********** Output from DALTON general input processing ***********
       *******************************************************************

 --------------------------------------------------------------------------------
   Overall default print level:    0
   Print level for DALTON.ERR :    1

    HERMIT 1- and 2-electron integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed (SIRIUS module)
    Static molecular property section will be executed (ABACUS module)
 --------------------------------------------------------------------------------


 Changes of defaults for *QM3   :
 --------------------------------

 +------------------+
 |  WORD: | CHANGE: |
 +------------------+
 |    QM3 |       T |
 | THDISC | 1.0e-12 |
 | MXDIIT |     180 |
 |  OLDTG |       T |
 |  PRINT |       1 |
 +------------------+
  Settings for determination of induced dipoles:
  Iterative Method is used


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
  Number of symmetry independent centres:    1
  The basis set is "STO-3G" from the basis set library.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/arnfinn/program/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  The basis set is "STO-3G" from the basis set library.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/arnfinn/program/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    2
  The basis set is "STO-3G" from the basis set library.
  Used basis set file for basis set for elements with Z =   1 :
     "/home/arnfinn/program/dalton/svn/pure_trunk/basis/STO-3G"

  Atomic type no.    4
  --------------------
  Nuclear charge:  -0.66900
  Number of symmetry independent centres:    2
  This is an MM atom without basis functions.

  Atomic type no.    5
  --------------------
  Nuclear charge:   0.33450
  Number of symmetry independent centres:    4
  This is an MM atom without basis functions.


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 


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

  Number of atom types:     6
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

  Threshold for integrals:  1.00e-12


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   36
  C       :     1  x  -3.0015786302    2  y  -1.4563174451    3  z   0.0550080380
  O       :     4  x  -3.1314330512    5  y   0.8240509855    6  z  -0.0184248298
  H       :     7  x  -1.1728925401    8  y  -2.4468589722    9  z   0.1025195325
  H       :    10  x  -4.7395144021   11  y  -2.6116034068   12  z   0.0761219481


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


  Max interatomic separation is    2.0185 Angstrom (    3.8144 Bohr)
  between atoms    3 and    2, "H     " and "O     ".

  Max QM+MM interatomic separation is    5.0272 Angstrom (    9.5000 Bohr)
  between the QM+MM centers     8 and    4, "H     " and "H     ".




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       1.798871         -0.056827    0.997867   -0.032124
   IB      13.009178          0.998345    0.057081    0.007042
   IC      14.808049         -0.008861    0.031671    0.999459


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         280942.3013          38847.8806          34128.6681 MHz
            9.371226            1.295826            1.138410 cm-1


@  Nuclear repulsion energy :   31.249215168341

 QM3 induced dipole vector converged in  13 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.020409711440023


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

 The following one-electron property integrals are calculated:
          - overlap integrals
          - dipole length integrals
          - Electric field at the nuclei
          - Potential energy at the nuclei

 Center of mass  (bohr):     -3.067740329786     -0.312997396742      0.018174673484
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERONE **************************
     ************************************************************************



     ************************************************************************
     ************************** Output from TWOINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-12

 Number of two-electron integrals written:        2907 ( 94.4% )
 Megabytes written:                              0.034

 >>>> Total CPU  time used in HERMIT:   0.03 seconds
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
     Number of Huckel orbitals each symmetry:   12

 Huckel EHT eigenvalues for symmetry :  1
          -20.705291     -11.377081      -1.496925      -0.962689      -0.743251
           -0.672270      -0.593413      -0.504635      -0.350830      -0.253214
           -0.198302      -0.193199

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Sat Dec 11 20:39:48 2010
     Host name              : stallo-1.local                          

 Title lines from ".mol" input file:
     CH2O PLUS 2 WATERS                                                      
     ------------------------                                                

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham DFT calculation.


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                1
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
     Abelian symmetry species          All |    1
                                       --- |  ---
     Occupied SCF orbitals               8 |    8
     Secondary orbitals                  4 |    4
     Total number of orbitals           12 |   12
     Number of basis functions          12 |   12

     Optimization information
     ========================
     Number of configurations                 1
     Number of orbital rotations             32
     ------------------------------------------
     Total number of variables               33

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-12
 
     This is a DFT calculation of type: B3LYP
 Weighted mixed functional:
               HF exchange:    0.20000
                       VWN:    0.19000
                       LYP:    0.81000
                     Becke:    0.72000
                    Slater:    0.80000


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Angular polynomials in range [15 35]
 Atom:    1*1 points=18676 compressed from 18676 ( 92 radial)
 Atom:    2*1 points=18280 compressed from 18280 ( 92 radial)
 Atom:    3*1 points=18150 compressed from 18150 ( 75 radial)
 Atom:    4*1 points=18150 compressed from 18150 ( 75 radial)
 Number of grid points:    73256 Grid generation time:       0.1 s
K-S energy, electrons, error :    -11.975378041393  16.0000012547    1.25e-06

 QM3 induced Dipole vector converged in  13 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.068478738902117
   1  -112.623001714     -1.770440365548e-02    2.45334e+00   -1.13e+02
      MULPOP C     5.02; O    16.31; H     2.28; 
K-S energy, electrons, error :    -11.746476967977  16.0000013834    1.38e-06

 QM3 induced Dipole vector converged in  13 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.009770909979341
   2  -112.274338902      1.474130243148e-02    2.94671e+00    3.49e-01
      MULPOP C     7.49; O    14.60; H     1.62; 
K-S energy, electrons, error :    -11.950280239732  16.0000013090    1.31e-06

 QM3 induced Dipole vector converged in  13 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.064920294806720
   3  -112.842566660     -1.706188241371e-02    8.61243e-01   -5.68e-01
      MULPOP C     5.61; O    15.92; H     1.81; 
K-S energy, electrons, error :    -11.862966932890  16.0000013128    1.31e-06

 QM3 induced Dipole vector converged in  12 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.035118618121378
   4  -112.900301896     -7.461039424276e-03    6.75709e-02   -5.77e-02
      MULPOP C     6.00; O    15.83; H     1.86; 
K-S energy, electrons, error :    -11.871215579687  16.0000013154    1.32e-06

 QM3 induced Dipole vector converged in  11 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037464088579593
   5  -112.900644179     -8.363423585919e-03    1.32509e-02   -3.42e-04
      MULPOP C     6.00; O    15.81; H     1.82; 
K-S energy, electrons, error :    -11.870525884180  16.0000013146    1.31e-06

 QM3 induced Dipole vector converged in  10 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037301621718752
   6  -112.900659093     -8.300430608597e-03    3.33594e-03   -1.49e-05
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870722893551  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in  10 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037363405920166
   7  -112.900659035     -8.323409781308e-03    4.54921e-04    5.87e-08

!!! Info: DIIS restarted because it was stalled ...
 - energy changed by    5.874e-08
 - gradient changed by -2.881e-03
 - or strange C vector coeff. for current index (=  7) :
    0.000001   -0.000037    0.000492   -0.001449   -0.009559    0.080684
    0.929867

! Backstep to previous Fock matrix because energy increased.
      MULPOP C     6.00; O    15.82; H     1.83; 
   7  Level shift: occupied orbital energies shifted by    0.0000
K-S energy, electrons, error :    -11.871217027125  16.0000013144    1.31e-06

 QM3 induced Dipole vector converged in  10 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037524385612735
   8  -112.900654010     -8.382979259313e-03    6.48038e-03    5.03e-06
      MULPOP C     5.99; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.869827554015  16.0000013147    1.31e-06

 QM3 induced Dipole vector converged in  11 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037072605673195
   9  -112.900650801     -8.215429745224e-03    1.30479e-02    3.21e-06
      MULPOP C     6.00; O    15.81; H     1.83; 
K-S energy, electrons, error :    -11.870752918494  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in  10 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037373699695971
  10  -112.900658917     -8.327223911504e-03    2.26417e-05   -8.12e-06
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754232049  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   7 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374037804598
  11  -112.900658913     -8.327350907024e-03    2.55380e-06    4.39e-09
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754385757  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   6 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107379982
  12  -112.900658912     -8.327376547050e-03    3.89640e-08    8.77e-10
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754387749  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   4 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107680715
  13  -112.900658912     -8.327376694033e-03    1.23229e-08    2.73e-12
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754386325  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   4 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107283510
  14  -112.900658912     -8.327376547890e-03    3.15778e-09   -4.95e-12
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754386665  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   3 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107382257
  15  -112.900658912     -8.327376587732e-03    5.45596e-10    1.28e-12
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754386619  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   3 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107367741
  16  -112.900658912     -8.327376582488e-03    1.89123e-11   -3.13e-13
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754386618  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   2 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107367388
  17  -112.900658912     -8.327376582338e-03    1.17106e-12    8.53e-14
      MULPOP C     6.00; O    15.82; H     1.83; 
K-S energy, electrons, error :    -11.870754386618  16.0000013145    1.31e-06

 QM3 induced Dipole vector converged in   2 iterations.
 Final norm2 of QM3 induced dipole moment vector:    0.037374107367411
  18  -112.900658912     -8.327376582353e-03    3.34857e-13   -2.84e-14
 DIIS converged in  18 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Sym       Kohn-Sham orbital energies

  1    -18.91905091   -10.08389823    -1.01341480    -0.57687160    -0.44434716
        -0.35005810    -0.32119123    -0.16543924     0.06967467     0.39381259
         0.49744572     0.61388889

    E(LUMO) :     0.06967467 au (symmetry 1)
  - E(HOMO) :    -0.16543924 au (symmetry 1)
  ------------------------------------------
    gap     :     0.23511391 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     QM/MM calculation converged :

        Electrostatic energy:      -0.007810375391
        Polarization energy:       -0.000517001192
        van der Waals energy:       0.000000000000
        Total QM/MM energy:        -0.008327376582

     Final DFT energy:           -112.900658911690                 
     Nuclear repulsion:            31.249215168341
     Electronic energy:          -144.141546703449

     Final gradient norm:           0.000000000000

 
     Date and time (Linux)  : Sat Dec 11 20:39:52 2010
     Host name              : stallo-1.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital           1        2        3        4        5        6        7
   1 C   :1s     0.0007   0.9910  -0.1289  -0.1866  -0.0051   0.0400   0.0000
   2 C   :1s    -0.0089   0.0406   0.2852   0.5872   0.0194  -0.1271   0.0005
   3 C   :2px    0.0005  -0.0001  -0.0092   0.0045   0.5760   0.0146  -0.0069
   4 C   :2py   -0.0078   0.0011   0.1690  -0.2633   0.0169  -0.4223   0.0232
   5 C   :2pz    0.0002   0.0000  -0.0046   0.0109   0.0065   0.0182   0.6014
   6 O   :1s     0.9935   0.0005  -0.2149   0.0979   0.0016  -0.1122   0.0012
   7 O   :1s     0.0296  -0.0081   0.7388  -0.4138  -0.0050   0.5655  -0.0059
   8 O   :2px    0.0004  -0.0002   0.0186  -0.0051   0.3295  -0.0220  -0.0067
   9 O   :2py   -0.0064   0.0031  -0.2119  -0.1023   0.0236   0.6807   0.0163
  10 O   :2pz    0.0002  -0.0001   0.0080   0.0030   0.0034  -0.0150   0.6810
  11 H   :1s     0.0003  -0.0088   0.0311   0.2538   0.3180   0.1282  -0.0028
  12 H   :1s     0.0003  -0.0088   0.0320   0.2725  -0.3153   0.1486  -0.0014

 Orbital           8        9       10
   1 C   :1s     0.0008   0.0000   0.2047
   2 C   :1s    -0.0051   0.0009  -1.2661
   3 C   :2px   -0.0893  -0.0076   0.0551
   4 C   :2py    0.0063   0.0267   0.4658
   5 C   :2pz   -0.0019   0.8270  -0.0130
   6 O   :1s     0.0041  -0.0002  -0.0224
   7 O   :1s    -0.0185   0.0004   0.1271
   8 O   :2px    0.9062   0.0074  -0.0146
   9 O   :2py    0.0386  -0.0224  -0.2103
  10 O   :2pz    0.0076  -0.7632   0.0060
  11 H   :1s    -0.3392  -0.0019   0.8467
  12 H   :1s     0.3482  -0.0006   0.9413



 >>>> Total CPU  time used in SIRIUS :      3.95 seconds
 >>>> Total wall time used in SIRIUS :      4.00 seconds

 
     Date and time (Linux)  : Sat Dec 11 20:39:52 2010
     Host name              : stallo-1.local                          

 NOTE:    1 warnings have been issued.
 Check output, result, and error files for "WARNING".


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



        *****************************************************************
        ******** Output from **PROPE input processing for ABACUS ********
        *****************************************************************



 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

 Default print level:        0

      Electronic excitation energies 
      Natural orbital connection is used
      for perturbation dependent basis sets.


 Changes of defaults for .EXCITA:
 --------------------------------

 Number of excitation energies:    2    0    0    0    0    0    0    0
 Print level          :    0
 Integral print level :    0
 Threshold            : 1.00e-09
 Maximum iterations   :   60

 Center of mass dipole origin  :   -3.067740   -0.312997    0.018175


                 .------------------------------------------------.
                 | Starting in Static Property Section (ABACUS) - |
                 `------------------------------------------------'


 
     Date and time (Linux)  : Sat Dec 11 20:39:52 2010
     Host name              : stallo-1.local                          
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*2 evaluation time:       0.5 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*1 evaluation time:       0.3 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*1 evaluation time:       0.3 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*1 evaluation time:       0.3 s
 Electrons: 16.000001( 1.31e-06): LR-DFT*1 evaluation time:       0.3 s
 >>> Time used in EXCITA is   4.94 seconds


   ***************************************************************************
   ************************ FINAL RESULTS from ABACUS ************************
   ***************************************************************************


 
     Date and time (Linux)  : Sat Dec 11 20:39:57 2010
     Host name              : stallo-1.local                          


                             Molecular geometry (au)
                             -----------------------

 C         -3.0015786302           -1.4563174451            0.0550080380
 O         -3.1314330512            0.8240509855           -0.0184248298
 H         -1.1728925401           -2.4468589722            0.1025195325
 H         -4.7395144021           -2.6116034068            0.0761219481





                        Molecular wave function and energy
                        ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0

     Total energy       -112.9006589117 au (Hartrees)
                         -3072.18312003 eV
                           -296420.6374 kJ/mol


                             Relativistic corrections
                             ------------------------

     Darwin correction:                          0.1684489058 au
     Mass-velocity correction:                  -0.2161794680 au

     Total relativistic correction:             -0.0477305622 au (0.0423%)
     Non-relativistic + relativistic energy:  -112.9483894739 au




                                  Dipole moment
                                  -------------

                 au               Debye          C m (/(10**-30)
              0.614666           1.562325           5.211356


                             Dipole moment components
                             ------------------------

                 au               Debye          C m (/(10**-30)

      x      0.07940303         0.20182234         0.67320685
      y     -0.60934621        -1.54880340        -5.16625206
      z      0.01438066         0.03655198         0.12192427


   Units:   1 a.u. =   2.54175 Debye 
            1 a.u. =   8.47835 (10**-30) C m (SI)




                Singlet electronic excitation energies
                --------------------------------------

                 Sym.   Mode   Frequency    Frequency
                ex. st.  No.      (au)          (eV)
                =======================================
                   1        1    0.150234    4.088065
                   1        2    0.349273    9.504191
                ---------------------------------------


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


  Max interatomic separation is    2.0185 Angstrom (    3.8144 Bohr)
  between atoms    3 and    2, "H     " and "O     ".

  Max QM+MM interatomic separation is    5.0272 Angstrom (    9.5000 Bohr)
  between the QM+MM centers     8 and    4, "H     " and "H     ".




 CPU time statistics for ABACUS
 ------------------------------

 EXCITA     00:00:05     100 %

 TOTAL      00:00:05     100 %


 >>>> Total CPU  time used in ABACUS:   4.95 seconds
 >>>> Total wall time used in ABACUS:   5.00 seconds


                   .-------------------------------------------.
                   | End of Static Property Section (ABACUS) - |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   8.93 seconds
 >>>> Total wall time used in DALTON:   9.00 seconds

 
     Date and time (Linux)  : Sat Dec 11 20:39:57 2010
     Host name              : stallo-1.local                          
END REFOUT

