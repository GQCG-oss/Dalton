########## Test description ########################
# 24-aug-07 hjaaj: changed from short to long, uses > 8 min on 2GHz pentium
START DESCRIPTION
KEYWORDS: triplet quadratic response dft pcm long
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enedft
tes
nuc
surf
nucchg
beta
qrlrve2
pcmsol
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON INPUT
.RUN RESPONSE
*PCM
.NEWQR
.SOLVNT
WATER
.ICESPH
2
.NESFP
4
.NPCMMT
0
*PCMCAV
.AREATS
0.4D0
.INA
1
2
3
4
.RIN
1.5
1.7
1.2
1.2
**WAVEFUNCTION
.DFT
B3LYP
*SCF INPUT
.THRESH
1.0D-8
**RESPONSE
!.TRPFLG
.MAXRM
200
*QUADRA
!.MAXIT
!100
!.PRINT
!10
.THCLR
1.0D-4
.APROP
XDIPLEN
.BPROP
ZDIPLEN
.CPROP
XDIPLEN
.ASPIN
1
.BSPIN
0
.CSPIN
1
*END OF INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
3-21G
Calculation of solvation energy

    3  0 0  X  Y       1.0D-15
       8.     1    3    1    1    1
O1  0.000000  0.000000 -2.2800000
       6.     1    3    1    1    1
C2  0.000000  0.000000  0.0000000
        1.    2    2    1    1
H3   1.78680   0.000000  1.140000
H3  -1.78680   0.000000  1.140000

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

     Date and time (Linux)  : Tue Dec 12 20:27:34 2006
     Host name              : platina                                 

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

     NON-EQ = F     NEQRSP =F
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.5000
     2    0.0000    0.0000    0.0000    1.7000
     3    0.0000    0.0000    0.0000    1.2000
     4    0.0000    0.0000    0.0000    1.2000

Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************



 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

  Calculation of solvation energy                                         
                                                                          
  Used basis set file for basis set for elements with Z =   8 :
     "/home/luca/programs/main-branch/dalton/basis/3-21G"
  Used basis set file for basis set for elements with Z =   6 :
     "/home/luca/programs/main-branch/dalton/basis/3-21G"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/luca/programs/main-branch/dalton/basis/3-21G"


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
   1    0.0000000000E+00    0.0000000000E+00   -0.2280000000E+01    0.1500000000E+01
   2    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00    0.1700000000E+01
   3    0.1786800000E+01    0.0000000000E+00    0.1140000000E+01    0.1200000000E+01
   4   -0.1786800000E+01    0.0000000000E+00    0.1140000000E+01    0.1200000000E+01


                             Isotopic Masses
                             ---------------

                           O1         15.994915
                           C2         12.000000
                           H3          1.007825
                           H3          1.007825

                       Total mass:    30.010565 amu
                       Natural abundance:  98.633 %

 Center-of-mass coordinates (A):    0.000000    0.000000   -1.138618


  Atoms and basis sets
  --------------------

  Number of atom types:     3
  Total number of atoms:    4

  Basis set used is "3-21G" from the basis set library.

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  O1          1    8.0000    15     9      [6s3p|3s2p]                                        
  C2          1    6.0000    15     9      [6s3p|3s2p]                                        
  H3          2    1.0000     3     2      [3s|2s]                                            
  ----------------------------------------------------------------------
  total:      4   16.0000    36    22
  ----------------------------------------------------------------------

  Threshold for integrals:  1.00E-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   12
  O1      :    1  x   0.0000000000   2  y   0.0000000000   3  z  -2.2800000000
  C2      :    4  x   0.0000000000   5  y   0.0000000000   6  z   0.0000000000
  H3      :    7  x   1.7868000000   8  y   0.0000000000   9  z   1.1400000000
  H3      :   10  x  -1.7868000000  11  y   0.0000000000  12  z   1.1400000000


   Interatomic separations (in Angstroms):
   ---------------------------------------

            O1          C2          H3          H3    
            ------      ------      ------      ------
 O1    :    0.000000
 C2    :    1.206524    0.000000
 H3    :    2.041901    1.121588    0.000000
 H3    :    2.041901    1.121588    1.891068    0.000000


  Max interatomic separation is    2.0419 Angstroms
  between atoms "H3    " and "O1    ".


  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C2         O1           1.206524
  bond distance:  H3         C2           1.121588
  bond distance:  H3         C2           1.121588


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     O1         C2         H3           122.538
  bond angle:     O1         C2         H3           122.538
  bond angle:     H3         C2         H3           114.923




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       1.802060          0.000000    0.000000    1.000000
   IB      13.122217          1.000000    0.000000    0.000000
   IC      14.924277          0.000000    1.000000    0.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

         280445.1440          38513.2324          33862.8793 MHz
            9.354643            1.284663            1.129544 cm-1


  Nuclear repulsion energy :   31.140735918780


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 Number of two-electron integrals written:       15102 ( 47.0% )
 Megabytes written:                              0.179



 MEMORY USED TO GENERATE CAVITY=    432042


 TOTAL NUMBER OF SPHERES=    4
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1    0.000000000    0.000000000   -1.206524035    1.800000000   22.860798949
   2    0.000000000    0.000000000    0.000000000    2.040000000   24.717903724
   3    0.945533836    0.000000000    0.603262017    1.440000000    9.680866932
   4   -0.945533836    0.000000000    0.603262017    1.440000000    9.680866932

 TOTAL NUMBER OF TESSERAE=     256
 SURFACE AREA=   66.94043654(A**2)    CAVITY VOLUME=   48.54512116 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   0.14 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.53 seconds

 >>>> Total CPU  time used in HERMIT:   0.67 seconds
 >>>> Total wall time used in HERMIT:   1.00 seconds

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
          -20.684793     -11.351633      -1.627151      -1.046455      -0.808323
           -0.702554      -0.608518      -0.491424      -0.320546      -0.162219
           -0.133092      -0.114392
 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Tue Dec 12 20:27:35 2006
     Host name              : platina                                 

 Title lines from integral program:
     Calculation of solvation energy                                         
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Kohn-Sham calculation.


     Time-dependent Kohn-Sham calculation (random phase approximation).


 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     Number of closed shell electrons         16
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Number of active orbitals                 0
     Total number of orbitals                 22

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
     Abelian symmetry species           1
                                       --
     Total number of orbitals          22
     Number of basis functions         22

      ** Automatic occupation of RKS orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals              8

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00E-08

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       15.96948
 NUCLEAR APPARENT CHARGE -15.79043 THEORETICAL -15.79589 NOT RENORMALIZED
 this is icompcm in icvev 0

  ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  16 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------
 Radial Quadrature : LMG scheme
 Space partitioning: Original Becke partitioning
 Radial integration threshold: 1e-13
 Atom:    1*1 points=18821 compressed from 18912 ( 96 radial)
 Atom:    2*1 points=19284 compressed from 19284 ( 96 radial)
 Atom:    3*1 points=18283 compressed from 18406 ( 77 radial)
 Atom:    4*1 points=18283 compressed from 18406 ( 77 radial)
 Number of grid points:    74671 Grid generation time:       0.3 s
K-S electrons/energy :   15.99999884968277  -11.72593556547978 err:-.12E-05
   1   -113.390743852694     -0.019326877785   2.60E+00  -1.13E+02    8
K-S electrons/energy :   15.99999897365772  -11.76206151000974 err:-.10E-05
   2   -112.933540650108     -0.018142730536   4.09E+00   4.57E-01    8
K-S electrons/energy :   15.99999890985154  -11.91058131584606 err:-.11E-05
   3   -113.776960774328     -0.012458708856   8.22E-01  -8.43E-01    8
K-S electrons/energy :   15.99999891568187  -11.76339857845007 err:-.11E-05
   4   -113.803503764256     -0.004704650058   2.05E-01  -2.65E-02    8
K-S electrons/energy :   15.99999891970862  -11.80178250150895 err:-.11E-05
   5   -113.805702888009     -0.005470897700   1.47E-02  -2.20E-03    8
K-S electrons/energy :   15.99999891944522  -11.79894861989234 err:-.11E-05
   6   -113.805715060129     -0.005331001388   2.52E-03  -1.22E-05    8
K-S electrons/energy :   15.99999891941809  -11.79923598040384 err:-.11E-05
   7   -113.805715381689     -0.005353798638   7.98E-05  -3.22E-07    8
K-S electrons/energy :   15.99999891941558  -11.79924729370095 err:-.11E-05
   8   -113.805715382024     -0.005354105026   4.06E-06  -3.35E-10    8
K-S electrons/energy :   15.99999891941536  -11.79924723991375 err:-.11E-05
   9   -113.805715382025     -0.005354120916   8.49E-07  -9.38E-13    8
K-S electrons/energy :   15.99999891941538  -11.79924729484149 err:-.11E-05
  10   -113.805715382025     -0.005354126312   4.43E-08   1.14E-13    8
K-S electrons/energy :   15.99999891941544  -11.79924729917216 err:-.11E-05
  11   -113.805715382025     -0.005354126633   4.02E-10  -2.84E-14    8
 DIIS converged in  11 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   16
 Orbital occupations :    8

 Kohn-Sham orbital energies, symmetry 1

       -19.05344558   -10.21721712    -1.07462042    -0.61727624    -0.48794101
        -0.43424520    -0.39564128    -0.25301270    -0.02631013     0.13403607
         0.21918797     0.28594713     0.68545946

    E(LUMO) :    -0.02631013 au (symmetry 1)
  - E(HOMO) :    -0.25301270 au (symmetry 1)
  ------------------------------------------
    gap     :     0.22670257 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final DFT energy:           -113.805715382025
     Nuclear repulsion:            31.140735918780
     Electronic energy:          -144.941097174172

     Final gradient norm:           0.000000000402

 
     Date and time (Linux)  : Tue Dec 12 20:28:01 2006
     Host name              : platina                                 

     Molecular orbitals for symmetry species   1

 Orbital           1        2        3        4        5        6        7
   1 O1  :1s    -0.9825  -0.0003   0.2112   0.0948   0.0000  -0.1001   0.0000
   2 O1  :1s    -0.1051  -0.0006  -0.2042  -0.0917   0.0000   0.0836   0.0000
   3 O1  :1s     0.0530   0.0013  -0.6363  -0.3615   0.0000   0.4119   0.0000
   4 O1  :2px    0.0000   0.0000   0.0000   0.0000   0.2459   0.0000   0.0000
   5 O1  :2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.4097
   6 O1  :2pz   -0.0045   0.0012  -0.1565   0.0766   0.0000  -0.4247   0.0000
   7 O1  :2px    0.0000   0.0000   0.0000   0.0000   0.2105   0.0000   0.0000
   8 O1  :2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.4154
   9 O1  :2pz    0.0142   0.0056  -0.1389   0.0603   0.0000  -0.4039   0.0000
  10 C2  :1s    -0.0002  -0.9850   0.1229  -0.1737   0.0000   0.0295   0.0000
  11 C2  :1s    -0.0005  -0.1029  -0.1348   0.1982   0.0000  -0.0613   0.0000
  12 C2  :1s    -0.0194   0.0493  -0.1382   0.5180   0.0000   0.0169   0.0000
  13 C2  :2px    0.0000   0.0000   0.0000   0.0000   0.3839   0.0000   0.0000
  14 C2  :2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.3142
  15 C2  :2pz    0.0025   0.0018   0.1688   0.1804   0.0000   0.3285   0.0000
  16 C2  :2px    0.0000   0.0000   0.0000   0.0000   0.2982   0.0000   0.0000
  17 C2  :2py    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.2967
  18 C2  :2pz    0.0153   0.0087  -0.0390   0.1285   0.0000   0.0700   0.0000
  19 H3  :1s     0.0008   0.0022  -0.0229   0.1682   0.1724   0.0749   0.0000
  20 H3  :1s    -0.0006  -0.0141   0.0178   0.0602   0.1053   0.0926   0.0000
  21 H3  :1s     0.0008   0.0022  -0.0229   0.1682  -0.1724   0.0749   0.0000
  22 H3  :1s    -0.0006  -0.0141   0.0178   0.0602  -0.1053   0.0926   0.0000

 Orbital           8        9       10       11       12       13
   1 O1  :1s     0.0000   0.0000   0.0074   0.0000  -0.1268   0.0000
   2 O1  :1s     0.0000   0.0000  -0.0145   0.0000   0.0531   0.0000
   3 O1  :1s     0.0000   0.0000  -0.0109   0.0000   1.5491   0.0000
   4 O1  :2px    0.4699   0.0000   0.0000  -0.1756   0.0000   0.0000
   5 O1  :2py    0.0000   0.3719   0.0000   0.0000   0.0000  -0.0336
   6 O1  :2pz    0.0000   0.0000  -0.0746   0.0000   0.2148   0.0000
   7 O1  :2px    0.5242   0.0000   0.0000  -0.3602   0.0000   0.0000
   8 O1  :2py    0.0000   0.5230   0.0000   0.0000   0.0000  -0.1389
   9 O1  :2pz    0.0000   0.0000  -0.1044   0.0000   0.7275   0.0000
  10 C2  :1s     0.0000   0.0000  -0.1509   0.0000   0.0716   0.0000
  11 C2  :1s     0.0000   0.0000   0.1034   0.0000   0.0185   0.0000
  12 C2  :1s     0.0000   0.0000   1.7360   0.0000  -1.0514   0.0000
  13 C2  :2px   -0.1477   0.0000   0.0000   0.4826   0.0000   0.0000
  14 C2  :2py    0.0000  -0.4358   0.0000   0.0000   0.0000  -1.0498
  15 C2  :2pz    0.0000   0.0000   0.2084   0.0000   0.1403   0.0000
  16 C2  :2px    0.0382   0.0000   0.0000   1.2598   0.0000   0.0000
  17 C2  :2py    0.0000  -0.6094   0.0000   0.0000   0.0000   1.0638
  18 C2  :2pz    0.0000   0.0000   0.6168   0.0000   1.6434   0.0000
  19 H3  :1s    -0.1805   0.0000  -0.1007  -0.0359   0.0198   0.0000
  20 H3  :1s    -0.3102   0.0000  -1.2341  -1.3562  -0.2157   0.0000
  21 H3  :1s     0.1805   0.0000  -0.1007   0.0359   0.0198   0.0000
  22 H3  :1s     0.3102   0.0000  -1.2341   1.3562  -0.2157   0.0000



 >>>> Total CPU  time used in SIRIUS :     24.30 seconds
 >>>> Total wall time used in SIRIUS :     26.00 seconds

 
     Date and time (Linux)  : Tue Dec 12 20:28:01 2006
     Host name              : platina                                 

- End of Wave Function Section



  This is output from RESPONSE  -  an MCSCF and SOPPA response property program
 ------------------------------------------------------------------------------



 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




 CHANGES OF DEFAULTS FOR RSPINP:
 -------------------------------


 USE FOCK TYPE DECOUPLING OF THE TWO-ELECTRON DENSITY MATRIX :
 ADD DV*(FC+FV) INSTEAD OF DV*FC TO E[2] APPROXIMATE ORBITAL DIAGONAL
  AVDIA = T


 Quadratic Response calculation
 ------------------------------

 First hyperpolarizability calculation : HYPCAL= T

 Spin of operator A , ISPINA=    1
 Spin of operator B , ISPINB=    0
 Spin of operator C , ISPINC=    1

  1 B-frequencies  0.000000E+00
  1 C-frequencies  0.000000E+00

 Print level                                    : IPRHYP =   2
 Maximum number of iterations in lin.rsp. solver: MAXITL =  60
 Threshold for convergence of linear resp. eq.s : THCLR  = 1.000E-04
 Maximum iterations in optimal orbital algorithm: MAXITO =   5
 Direct one-index transformation                : DIROIT= T

    1 A OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          XDIPLEN 

    1 B OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          ZDIPLEN 

    1 C OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          XDIPLEN 


   SCF energy         :     -113.805715382024644
 -- inactive part     :     -144.941097174171745
 -- nuclear repulsion :       31.140735918779654


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    1


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:     112
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :     112


 QRLRVE -- linear response calculation for symmetry  1
 QRLRVE -- operator label : ZDIPLEN 
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   F

 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   8.97E-06

 *** RSPCTL MICROITERATIONS CONVERGED

 QRLRVE: SINGLET SOLUTION   LABEL   ZDIPLEN     FREQUENCY   0.000000E+00
 SYMMETRY    1

@QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.00000):     19.7071041996    


 Linear response calculations for quadratic response
 - triplet property operator(s) of symmetry    1


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       T
 Orbital variables.         KZWOPT:     112
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :     112


 QRTRVE -- linear response calculation for symmetry  1
 QRTRVE -- operator label : XDIPLEN 
 QRTRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry =  1; triplet =   T

 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s
 Electrons: 15.999999(-1.08e-06): LR-DFT*1 evaluation time:       1.7 s

 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00E-04
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   5.53E-05

 *** RSPCTL MICROITERATIONS CONVERGED
 
 ======================================================================
 >>>>>>>>    L I N E A R   R E S P O N S E   F U N C T I O N   <<<<<<<<
 ======================================================================

 The -<<A;B>>(omega_b) functions from vectors generated
 in a *QUADRA calculation of <<A;B,C>>(omega_b,omega_c)

 Note: the accuracy of off-diagonal elements will be linear
 in the convergence threshold THCLR =  1.00E-04

 All zero because spin symmetry of A and B is different.
 ISPINA and ISPINB :    1   0


  Results from quadratic response calculation
 --------------------------------------------

 DFT-QR computed in a linearly-scaling fashion.

 Electrons: 15.999999(-1.08e-06): QR-DFT/b evaluation time:       2.2 s

@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,X) =   -125.80876395

 >>>> Total CPU  time used in RESPONSE:  25.16 seconds
 >>>> Total wall time used in RESPONSE:  27.00 seconds
 >>>> Total CPU  time used in DALTON:  50.13 seconds
 >>>> Total wall time used in DALTON:  54.00 seconds

 
     Date and time (Linux)  : Tue Dec 12 20:28:28 2006
     Host name              : platina                                 
END REFOUT

