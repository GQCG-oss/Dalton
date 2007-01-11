########## Test description ########################
START DESCRIPTION
KEYWORDS: pcm hf ecd 1-fluoroethanol medium
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
dipole
nuc
tes
sym
cmass
gauge_or
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
5
.ICESPH
2
*PCMCAV
.RIN
1.7
2.0
1.4
1.6   
1.2
.INA
1
2
3
4
5
.AREATS
0.70
**WAVE FUNCTION
.HF
**PROPERTIES
.ECD
*EXCITA
.NEXCIT
 3
**END OF INPUT
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
cc-pVDZ


    4    0          
        6.    2
C         -0.1077810630            0.0016311448           -0.0277533647
C         -0.0509885280           -0.0063400396            2.8152173831
        9.    1
F          2.3195193068            0.0622555881           -0.9275026109
        8.    1
O         -1.2729658590           -2.1533608644           -0.8603501924
        1.    5
H         -1.0158622780            1.6763033629           -0.7793105516
H         -1.2942535134           -2.1745591590           -2.6533440159
H         -1.9586539793           -0.0071208309            3.5572943234
H          0.9181052495           -1.6875960201            3.4664685991
H          0.9386804771            1.6487921657            3.5023462928

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

     Date and time (Linux)  : Fri Mar 24 13:25:24 2006
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

     ICESPH =       2     NESFP =       5
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG 60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.7000
     2    0.0000    0.0000    0.0000    2.0000
     3    0.0000    0.0000    0.0000    1.4000
     4    0.0000    0.0000    0.0000    1.6000
     5    0.0000    0.0000    0.0000    1.2000

Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************



 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

                                                                          
                                                                          
  Used basis set file for basis set for elements with Z =   6 :
     "/home/lara/programs/main-branch/dalton/basis/cc-pVDZ"
  Used basis set file for basis set for elements with Z =   9 :
     "/home/lara/programs/main-branch/dalton/basis/cc-pVDZ"
  Used basis set file for basis set for elements with Z =   8 :
     "/home/lara/programs/main-branch/dalton/basis/cc-pVDZ"
  Used basis set file for basis set for elements with Z =   1 :
     "/home/lara/programs/main-branch/dalton/basis/cc-pVDZ"


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
   1   -0.1077810630E+00    0.1631144800E-02   -0.2775336470E-01    0.1700000000E+01
   2   -0.5098852800E-01   -0.6340039600E-02    0.2815217383E+01    0.2000000000E+01
   3    0.2319519307E+01    0.6225558810E-01   -0.9275026109E+00    0.1400000000E+01
   4   -0.1272965859E+01   -0.2153360864E+01   -0.8603501924E+00    0.1600000000E+01
   5   -0.1015862278E+01    0.1676303363E+01   -0.7793105516E+00    0.1200000000E+01


                             Isotopic Masses
                             ---------------

                           C          12.000000
                           C          12.000000
                           F          18.998403
                           O          15.994915
                           H           1.007825
                           H           1.007825
                           H           1.007825
                           H           1.007825
                           H           1.007825

                       Total mass:    64.032443 amu
                       Natural abundance:  97.504 %

 Center-of-mass coordinates (A):    0.302504   -0.528873    0.143931


  Atoms and basis sets
  --------------------

  Number of atom types:     4
  Total number of atoms:    9

  Basis set used is "cc-pVDZ" from the basis set library.

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C           2    6.0000    26    14      [9s4p1d|3s2p1d]                                    
  F           1    9.0000    26    14      [9s4p1d|3s2p1d]                                    
  O           1    8.0000    26    14      [9s4p1d|3s2p1d]                                    
  H           5    1.0000     7     5      [4s1p|2s1p]                                        
  ----------------------------------------------------------------------
  total:      9   34.0000   139    81
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00E-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates:   27

   1   C        x     -0.1077810630
   2            y      0.0016311448
   3            z     -0.0277533647

   4   C        x     -0.0509885280
   5            y     -0.0063400396
   6            z      2.8152173831

   7   F        x      2.3195193068
   8            y      0.0622555881
   9            z     -0.9275026109

  10   O        x     -1.2729658590
  11            y     -2.1533608644
  12            z     -0.8603501924

  13   H        x     -1.0158622780
  14            y      1.6763033629
  15            z     -0.7793105516

  16   H        x     -1.2942535134
  17            y     -2.1745591590
  18            z     -2.6533440159

  19   H        x     -1.9586539793
  20            y     -0.0071208309
  21            z      3.5572943234

  22   H        x      0.9181052495
  23            y     -1.6875960201
  24            z      3.4664685991

  25   H        x      0.9386804771
  26            y      1.6487921657
  27            z      3.5023462928


   Interatomic separations (in Angstroms):
   ---------------------------------------

            C           C           F           O           H           H     
            ------      ------      ------      ------      ------      ------
 C     :    0.000000
 C     :    1.504741    0.000000
 F     :    1.370253    2.344677    0.000000
 O     :    1.369216    2.343528    2.233818    0.000000
 H     :    1.083713    2.161409    1.962376    2.031586    0.000000
 H     :    1.910707    3.181761    2.427374    0.948945    2.271067    0.000000
 H     :    2.135043    1.083181    3.280087    2.624212    2.511720    3.498628
 H     :    2.124336    1.083185    2.610347    2.578298    3.043698    3.453212
 H     :    2.134482    1.083335    2.594986    3.278374    2.490710    4.012564

            H           H           H     
            ------      ------      ------
 H     :    0.000000
 H     :    1.763676    0.000000
 H     :    1.766185    1.765676    0.000000


  Max interatomic separation is    4.0126 Angstroms
  between atoms "H     " and "H     ".


  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C          C            1.504741
  bond distance:  F          C            1.370253
  bond distance:  O          C            1.369216
  bond distance:  H          C            1.083713
  bond distance:  H          O            0.948945
  bond distance:  H          C            1.083181
  bond distance:  H          C            1.083185
  bond distance:  H          C            1.083335


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     C          C          F            109.193
  bond angle:     C          C          O            109.170
  bond angle:     C          C          H            112.214
  bond angle:     F          C          O            109.258
  bond angle:     F          C          H            105.604
  bond angle:     O          C          H            111.303
  bond angle:     C          C          H            110.111
  bond angle:     C          C          H            109.260
  bond angle:     C          C          H            110.057
  bond angle:     H          C          H            109.000
  bond angle:     H          C          H            109.218
  bond angle:     H          C          H            109.171
  bond angle:     C          O          H            109.690




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA   52.248587         -0.216669    0.209580    0.953483
   IB   56.452522          0.841008    0.536036    0.073287
   IC   96.963569         -0.495741    0.817766   -0.292401


 Rotational constants
 --------------------

               A                   B                   C

           9672.5871           8952.2839           5212.0504 MHz
            0.322643            0.298616            0.173855 cm-1


  Nuclear repulsion energy :  135.226836327282


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in HUCKEL is   0.13 seconds


 >>> Time used in ONEDRV is   0.16 seconds


 Number of two-electron integrals written:     5495998 ( 99.6% )
 Megabytes written:                             62.932



 >>> Time used in TWOINT is   7.28 seconds


 MEMORY USED TO GENERATE CAVITY=    432042

 TESSERA SPEZZATA IN TRONCONI

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
IT IS VALUABLE ALMOST ONLY FOR TESTING

 TOTAL NUMBER OF SPHERES=    5
 SPHERE             CENTER  (X,Y,Z) (A)               RADIUS (A)      AREA(A*A)
   1   -0.057035282    0.000863165   -0.014686448    2.040000000    3.740433668
   2   -0.026981967   -0.003355004    1.489748876    2.400000000   50.281695884
   3    1.227436751    0.032944238   -0.490813242    1.680000000   17.505819464
   4   -0.673624520   -1.139509491   -0.455277713    1.920000000   25.650918549
   5   -0.537571164    0.887061534   -0.412393382    1.440000000    8.646852065

 TOTAL NUMBER OF TESSERAE=     245
 SURFACE AREA=  105.82571963(A**2)    CAVITY VOLUME=   92.06411535 (A**3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATING CAVITY .....

 >>> Time used in PEDRAM is   7.64 seconds


  ..... DONE GENERATING -Q-  MATRIX .....

 >>> Time used in Q-MAT  is   0.20 seconds

 >>>> Total CPU  time used in HERMIT:   7.85 seconds
 >>>> Total wall time used in HERMIT:  11.00 seconds

- End of Integral Section


Starting in Wave Function Section -


 *** Output from Huckel module :

     Using EWMO model:          T
     Using EHT  model:          F
     Number of Huckel orbitals each symmetry:   25

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
        which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -26.384283     -20.685288     -11.354910     -11.348708      -1.986117
           -1.595436      -1.446837      -1.023120      -0.866682      -0.805631
           -0.772840      -0.692346      -0.623528      -0.614481      -0.579041
           -0.540562      -0.455981      -0.219901      -0.219072      -0.162851
           -0.128588      -0.122676      -0.118179      -0.113677      -0.106764
 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Fri Mar 24 13:25:36 2006
     Host name              : star.chem.uit.no                        

 Title lines from integral program:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.


 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     Number of closed shell electrons         34
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Number of active orbitals                 0
     Total number of orbitals                 81

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species           1
                                       --
     Total number of orbitals          81
     Number of basis functions         81

      ** Automatic occupation of RHF orbitals **
      -- Initial occupation of symmetries is determined from Huckel guess.                    
      -- Initial occupation of symmetries is : --

     Occupied SCF orbitals             17

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00E-06

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       33.94200
 NUCLEAR APPARENT CHARGE -33.55690 THEORETICAL -33.56627 NOT RENORMALIZED
 this is icompcm in icvev 0

  ..... DONE WITH INDUCED NUCLEAR CHARGES .....

 >>> Time used in VNN    is   0.00 seconds



 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Automatic occupation of symmetries with  34 electrons.

 Iter     Total energy      Solvation energy  Error norm  Delta(E)    HF occupation
 ----------------------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00E-15

 >>> Time used in FRMSUP is   1.15 seconds

   1   -252.188970817199     -0.016197485737   3.48E+00  -2.52E+02   17
   2   -252.835208674230     -0.009564151352   1.69E+00  -6.46E-01   17
   3   -252.972507643799     -0.011121715162   3.11E-01  -1.37E-01   17
   4   -252.978092679027     -0.011340397237   7.78E-02  -5.59E-03   17
   5   -252.978518757816     -0.011681471661   1.92E-02  -4.26E-04   17
   6   -252.978559691654     -0.011745396176   5.26E-03  -4.09E-05   17
   7   -252.978563412086     -0.011779504211   1.76E-03  -3.72E-06   17
   8   -252.978563741345     -0.011788102748   6.26E-04  -3.29E-07   17
   9   -252.978563785550     -0.011788995502   1.01E-04  -4.42E-08   17
  10   -252.978563787079     -0.011789167400   1.83E-05  -1.53E-09   17
  11   -252.978563787119     -0.011789175530   5.68E-06  -4.03E-11   17
  12   -252.978563787120     -0.011789167303   1.17E-06  -8.12E-13   17
  13   -252.978563787120     -0.011789165586   2.68E-07   2.23E-14   17
 DIIS converged in  13 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   34
 Orbital occupations :   17

 Hartree-Fock orbital energies, symmetry 1

       -26.28287786   -20.56613836   -11.35729927   -11.22167329    -1.61021825
        -1.38734416    -1.02405080    -0.85795006    -0.75328999    -0.72505410
        -0.68785060    -0.62232755    -0.56693054    -0.55035065    -0.53830711
        -0.51671866    -0.46872537     0.19812997     0.22853080     0.25326087
         0.26886597     0.27881928

    E(LUMO) :     0.19812997 au (symmetry 1)
  - E(HOMO) :    -0.46872537 au (symmetry 1)
  ------------------------------------------
    gap     :     0.66685534 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -252.978563787120
     Nuclear repulsion:           135.226836327282
     Electronic energy:          -388.193610948815

     Final gradient norm:           0.000000267823

 
     Date and time (Linux)  : Fri Mar 24 13:27:33 2006
     Host name              : star.chem.uit.no                        

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5        6        7
   1  C   1s    -0.0002  -0.0001  -1.0013  -0.0026   0.0067  -0.0049   0.0032
   2  C   1s    -0.0005  -0.0004  -0.0073  -0.0001  -0.2418   0.2613  -0.4762
   3  C   1s     0.0039   0.0041   0.0090   0.0023   0.1433  -0.1183   0.1367
   4  C   2px    0.0005  -0.0002  -0.0011   0.0000  -0.1212  -0.1104  -0.0002
   5  C   2py    0.0000  -0.0003   0.0011   0.0000   0.0158  -0.1518  -0.0422
   6  C   2pz   -0.0003  -0.0003   0.0011   0.0002   0.0511  -0.0485  -0.1993
   7  C   2px    0.0012  -0.0004   0.0003   0.0001   0.0425   0.0372   0.0079
   8  C   2py   -0.0005  -0.0014  -0.0002   0.0001  -0.0110   0.0597   0.0065
   9  C   2pz   -0.0001  -0.0001   0.0010   0.0008  -0.0140   0.0212   0.0248
  10  C   3d2-   0.0000  -0.0001  -0.0001   0.0000  -0.0033   0.0189  -0.0004
  11  C   3d1-   0.0000   0.0000   0.0000   0.0000  -0.0033   0.0199   0.0060
  12  C   3d0   -0.0001   0.0001   0.0002  -0.0002   0.0087  -0.0039  -0.0179
  13  C   3d1+  -0.0001   0.0000   0.0000   0.0000   0.0157   0.0125  -0.0004
  14  C   3d2+   0.0001   0.0001   0.0000   0.0000  -0.0211  -0.0174  -0.0018
  15  C   1s     0.0000   0.0000  -0.0023   1.0019   0.0014   0.0002   0.0209
  16  C   1s     0.0001   0.0000   0.0003   0.0106  -0.0310   0.0562  -0.6117
  17  C   1s    -0.0010  -0.0010  -0.0028  -0.0127  -0.0043  -0.0255   0.1846
  18  C   2px    0.0000   0.0000   0.0000   0.0000  -0.0058  -0.0084   0.0025
  19  C   2py    0.0000   0.0000   0.0000   0.0000   0.0008  -0.0103  -0.0047
  20  C   2pz    0.0002   0.0002   0.0003   0.0002   0.0255  -0.0303   0.1058
  21  C   2px   -0.0002   0.0000   0.0000   0.0000  -0.0005   0.0061   0.0009
  22  C   2py    0.0002   0.0002  -0.0001   0.0001  -0.0005   0.0091   0.0010
  23  C   2pz    0.0004   0.0003   0.0015   0.0008   0.0030   0.0136  -0.0456
  24  C   3d2-   0.0000   0.0001   0.0000   0.0000   0.0004  -0.0015   0.0007
  25  C   3d1-   0.0000   0.0000   0.0001   0.0000  -0.0004   0.0027   0.0021
  26  C   3d0   -0.0001  -0.0001   0.0000  -0.0002  -0.0061   0.0046  -0.0138
  27  C   3d1+   0.0001   0.0000  -0.0001   0.0000   0.0018   0.0016  -0.0006
  28  C   3d2+   0.0000  -0.0001   0.0000   0.0000   0.0015   0.0012   0.0005
  29  F   1s     1.0006   0.0001  -0.0001   0.0000   0.0074  -0.0003   0.0004
  30  F   1s     0.0025   0.0003  -0.0006  -0.0001  -0.8923  -0.2037   0.1633
  31  F   1s    -0.0042  -0.0009   0.0005   0.0000  -0.0034   0.0032   0.0242
  32  F   2px   -0.0013   0.0001   0.0004   0.0000   0.0895  -0.0260   0.1240
  33  F   2py    0.0000   0.0000   0.0000   0.0000   0.0040  -0.0177  -0.0139
  34  F   2pz    0.0005   0.0000  -0.0001  -0.0001  -0.0327  -0.0015  -0.1057
  35  F   2px    0.0020   0.0002  -0.0009   0.0000  -0.0118  -0.0036  -0.0018
  36  F   2py    0.0001   0.0002   0.0000   0.0000  -0.0002   0.0010   0.0004
  37  F   2pz   -0.0007   0.0000   0.0002   0.0002   0.0047   0.0017   0.0093
  38  F   3d2-   0.0000   0.0000   0.0000   0.0000  -0.0005   0.0012   0.0005
  39  F   3d1-   0.0000   0.0000   0.0000   0.0000   0.0001  -0.0004   0.0000
  40  F   3d0    0.0001   0.0000  -0.0001   0.0000   0.0025  -0.0010   0.0005
  41  F   3d1+   0.0001   0.0000  -0.0001   0.0000   0.0045  -0.0003   0.0048
  42  F   3d2+  -0.0001   0.0000   0.0002   0.0000  -0.0061   0.0014  -0.0037
  43  O   1s     0.0001   1.0008   0.0000   0.0000   0.0030  -0.0106  -0.0010
  44  O   1s     0.0003   0.0039  -0.0006  -0.0002  -0.1042   0.8226   0.1892
  45  O   1s    -0.0010  -0.0052   0.0006   0.0000  -0.0025  -0.0720  -0.0097
  46  O   2px    0.0000   0.0008  -0.0001   0.0000  -0.0251   0.0489  -0.0365
  47  O   2py    0.0000   0.0016  -0.0003   0.0000  -0.0266   0.0971  -0.0831
  48  O   2pz    0.0000  -0.0013  -0.0001  -0.0001   0.0009  -0.0452  -0.1399
  49  O   2px   -0.0005  -0.0013   0.0004   0.0000   0.0011  -0.0098  -0.0002
  50  O   2py   -0.0002  -0.0023   0.0007   0.0000  -0.0011  -0.0198   0.0028
  51  O   2pz   -0.0001   0.0004  -0.0002   0.0004  -0.0065   0.0420   0.0401
  52  O   3d2-   0.0000  -0.0001   0.0002   0.0000  -0.0038   0.0059  -0.0037
  53  O   3d1-   0.0000   0.0000   0.0000   0.0000  -0.0008   0.0030  -0.0043
  54  O   3d0    0.0000   0.0002   0.0000   0.0000   0.0016   0.0002   0.0046
  55  O   3d1+   0.0000   0.0000   0.0000   0.0000  -0.0010   0.0016  -0.0024
  56  O   3d2+   0.0000   0.0001  -0.0001   0.0000   0.0008  -0.0042   0.0029
  57  H   1s    -0.0004  -0.0005  -0.0007  -0.0003  -0.0509   0.0550  -0.1577
  58  H   1s     0.0002   0.0002  -0.0010   0.0002   0.0219  -0.0264   0.0758
  59  H   2px   -0.0002   0.0001   0.0001   0.0000  -0.0042   0.0021  -0.0055
  60  H   2py    0.0000   0.0001   0.0000   0.0000   0.0062  -0.0073   0.0094
  61  H   2pz    0.0000   0.0000  -0.0001   0.0001  -0.0021   0.0011  -0.0055
  62  H   1s     0.0000   0.0000  -0.0005   0.0004  -0.0289   0.2843   0.1789
  63  H   1s    -0.0002  -0.0001   0.0001   0.0000   0.0108  -0.1398  -0.0864
  64  H   2px    0.0000   0.0000  -0.0001   0.0000  -0.0010   0.0017   0.0002
  65  H   2py   -0.0001   0.0001  -0.0002   0.0000  -0.0008   0.0032  -0.0009
  66  H   2pz    0.0000  -0.0007  -0.0002   0.0001  -0.0049   0.0381   0.0190
  67  H   1s     0.0000   0.0001   0.0001   0.0001  -0.0050   0.0150  -0.2230
  68  H   1s     0.0001  -0.0002  -0.0003   0.0012  -0.0005  -0.0060   0.1065
  69  H   2px    0.0000   0.0000  -0.0001  -0.0005  -0.0007   0.0010  -0.0177
  70  H   2py    0.0000   0.0000   0.0000   0.0000   0.0002  -0.0014  -0.0002
  71  H   2pz    0.0000   0.0001   0.0001   0.0001   0.0012  -0.0012   0.0079
  72  H   1s     0.0001   0.0001   0.0002   0.0002  -0.0049   0.0152  -0.2194
  73  H   1s    -0.0001  -0.0001  -0.0004   0.0012   0.0009  -0.0056   0.1025
  74  H   2px    0.0000   0.0000   0.0000   0.0002   0.0001  -0.0019   0.0083
  75  H   2py    0.0000   0.0000   0.0000  -0.0004  -0.0009   0.0003  -0.0152
  76  H   2pz    0.0001   0.0001   0.0001   0.0001   0.0004  -0.0013   0.0074
  77  H   1s     0.0000  -0.0001   0.0001   0.0001  -0.0047   0.0069  -0.2253
  78  H   1s    -0.0001   0.0002  -0.0003   0.0012   0.0019  -0.0077   0.1052
  79  H   2px    0.0000   0.0000   0.0000   0.0002  -0.0001  -0.0004   0.0093
  80  H   2py    0.0000   0.0000   0.0000   0.0004   0.0010  -0.0004   0.0158
  81  H   2pz    0.0001   0.0000   0.0001   0.0001   0.0004  -0.0015   0.0077

 Orbital          8        9       10       11       12       13       14
   1  C   1s    -0.0074  -0.0025  -0.0049   0.0003   0.0011   0.0030  -0.0010
   2  C   1s    -0.5463  -0.0168   0.0436  -0.0545  -0.1220   0.0034   0.0012
   3  C   1s     0.2448   0.0075   0.0383  -0.0055   0.0007  -0.0541   0.0174
   4  C   2px   -0.0102   0.0670  -0.4464  -0.0286   0.1770  -0.1001  -0.0258
   5  C   2py   -0.0616   0.4198   0.1042  -0.1135  -0.0886  -0.2134  -0.0160
   6  C   2pz    0.1996  -0.1023  -0.0111  -0.2705  -0.0349  -0.2267   0.0272
   7  C   2px    0.0303  -0.0302   0.1529   0.0089  -0.0466   0.0422   0.0246
   8  C   2py   -0.0034  -0.1008  -0.0203   0.0215   0.0489   0.0954   0.0012
   9  C   2pz   -0.0516   0.0142  -0.0175   0.0370   0.0068   0.0322  -0.0098
  10  C   3d2-  -0.0036  -0.0085   0.0174  -0.0099   0.0233   0.0221   0.0035
  11  C   3d1-   0.0005  -0.0062  -0.0072   0.0130  -0.0174   0.0113   0.0144
  12  C   3d0    0.0155   0.0151   0.0044   0.0115  -0.0096  -0.0078   0.0006
  13  C   3d1+   0.0086  -0.0023   0.0140  -0.0187  -0.0259   0.0303   0.0118
  14  C   3d2+  -0.0049   0.0172  -0.0283  -0.0085   0.0018  -0.0107  -0.0180
  15  C   1s    -0.0129  -0.0044  -0.0009  -0.0051  -0.0002  -0.0021   0.0005
  16  C   1s     0.4006   0.0755  -0.0270   0.1270   0.0163   0.0049   0.0006
  17  C   1s    -0.1293   0.0029   0.0057   0.0121  -0.0009   0.0286  -0.0056
  18  C   2px    0.0068   0.0166  -0.1493  -0.0548   0.0245   0.2970   0.4333
  19  C   2py   -0.0159   0.1334   0.0358  -0.0333  -0.1304  -0.2111   0.4280
  20  C   2pz    0.1478   0.0790  -0.0149   0.2318   0.0890   0.1850  -0.0326
  21  C   2px   -0.0024  -0.0043   0.0220  -0.0024  -0.0222  -0.0320  -0.0851
  22  C   2py    0.0062  -0.0322  -0.0102   0.0201  -0.0074   0.0234  -0.0668
  23  C   2pz   -0.0241  -0.0205  -0.0040  -0.0454  -0.0199  -0.0396   0.0064
  24  C   3d2-   0.0010   0.0033  -0.0002  -0.0012  -0.0029  -0.0078   0.0157
  25  C   3d1-   0.0021  -0.0051  -0.0043   0.0031  -0.0004   0.0022   0.0084
  26  C   3d0   -0.0132  -0.0051   0.0009  -0.0134  -0.0044  -0.0099   0.0010
  27  C   3d1+   0.0009   0.0003   0.0074  -0.0029  -0.0054   0.0077   0.0095
  28  C   3d2+   0.0012  -0.0017   0.0052   0.0011  -0.0020  -0.0083  -0.0152
  29  F   1s     0.0020  -0.0003   0.0010   0.0005   0.0003   0.0002   0.0010
  30  F   1s     0.1870  -0.0447   0.1308  -0.0075  -0.0140  -0.0051  -0.0075
  31  F   1s     0.0277  -0.0161   0.0391  -0.0078  -0.0129   0.0016  -0.0126
  32  F   2px    0.3914  -0.1727   0.4384  -0.2760  -0.3499   0.1805   0.0570
  33  F   2py   -0.0343   0.3691   0.2723  -0.3247   0.6102   0.3445  -0.0336
  34  F   2pz   -0.0410  -0.0496  -0.3877  -0.6211  -0.2423   0.1997  -0.2245
  35  F   2px   -0.0010  -0.0013   0.0059  -0.0035  -0.0108   0.0072   0.0010
  36  F   2py    0.0018  -0.0011   0.0036  -0.0068   0.0238   0.0182  -0.0028
  37  F   2pz   -0.0036   0.0007  -0.0034  -0.0108  -0.0107   0.0168  -0.0158
  38  F   3d2-   0.0009  -0.0078  -0.0048   0.0050  -0.0061  -0.0024   0.0006
  39  F   3d1-  -0.0002   0.0032   0.0016  -0.0015   0.0024   0.0010   0.0002
  40  F   3d0    0.0054  -0.0021   0.0015  -0.0064  -0.0037   0.0014  -0.0003
  41  F   3d1+   0.0044  -0.0001   0.0103   0.0079   0.0004   0.0001   0.0029
  42  F   3d2+  -0.0096   0.0034  -0.0080   0.0042   0.0047  -0.0014  -0.0011
  43  O   1s     0.0018  -0.0011   0.0000   0.0014  -0.0018  -0.0015  -0.0002
  44  O   1s     0.1762  -0.0335   0.0065  -0.1202   0.1607  -0.1170   0.0739
  45  O   1s    -0.0196  -0.0293   0.0051  -0.0531   0.0970  -0.0365   0.0353
  46  O   2px   -0.1026  -0.1838  -0.1559   0.1792   0.0239   0.1985   0.1166
  47  O   2py   -0.2169  -0.1328   0.2113   0.1517  -0.2346   0.4703  -0.0680
  48  O   2pz   -0.1160  -0.4773   0.2014  -0.1191   0.3444  -0.1599   0.1651
  49  O   2px   -0.0021   0.0136  -0.0082   0.0009  -0.0020   0.0021  -0.0001
  50  O   2py    0.0074   0.0170  -0.0099  -0.0031  -0.0065   0.0161  -0.0120
  51  O   2pz    0.0310   0.0698  -0.0255   0.0093  -0.0118  -0.0013  -0.0048
  52  O   3d2-  -0.0083  -0.0072  -0.0027   0.0060   0.0001   0.0099   0.0007
  53  O   3d1-  -0.0008  -0.0097   0.0027  -0.0039   0.0059  -0.0065   0.0034
  54  O   3d0    0.0093   0.0120  -0.0060  -0.0034  -0.0042  -0.0045  -0.0038
  55  O   3d1+  -0.0006  -0.0044   0.0013  -0.0034   0.0028  -0.0019   0.0004
  56  O   3d2+   0.0062   0.0007  -0.0100  -0.0018   0.0055  -0.0076   0.0021
  57  H   1s    -0.3097   0.2590   0.2537  -0.0123  -0.2081  -0.0586  -0.0031
  58  H   1s     0.1372  -0.1035  -0.0880   0.0033   0.0677   0.0161   0.0014
  59  H   2px   -0.0098   0.0084  -0.0005  -0.0008  -0.0043  -0.0036  -0.0014
  60  H   2py    0.0164  -0.0084  -0.0112  -0.0006   0.0067  -0.0001  -0.0008
  61  H   2pz   -0.0057   0.0062   0.0073  -0.0030  -0.0052  -0.0032  -0.0001
  62  H   1s     0.1833   0.4264  -0.1843   0.0535  -0.2489   0.0761  -0.1355
  63  H   1s    -0.0713  -0.1738   0.0735  -0.0252   0.1064  -0.0624   0.0558
  64  H   2px   -0.0016  -0.0025  -0.0051   0.0084  -0.0002   0.0038   0.0040
  65  H   2py   -0.0046   0.0014   0.0031   0.0059  -0.0081   0.0133  -0.0023
  66  H   2pz    0.0175   0.0355  -0.0144   0.0053  -0.0181   0.0044  -0.0089
  67  H   1s     0.2199   0.0527   0.0990   0.1843   0.0167  -0.2100  -0.4407
  68  H   1s    -0.0840  -0.0213  -0.0358  -0.0740  -0.0053   0.0620   0.1159
  69  H   2px    0.0162   0.0043   0.0040   0.0091   0.0008  -0.0048  -0.0138
  70  H   2py   -0.0001   0.0027   0.0003  -0.0009  -0.0018  -0.0042   0.0088
  71  H   2pz   -0.0042   0.0002  -0.0030  -0.0015  -0.0003   0.0086   0.0087
  72  H   1s     0.2337  -0.0254  -0.1021   0.1391   0.1455   0.3797  -0.1600
  73  H   1s    -0.0912   0.0056   0.0388  -0.0478  -0.0602  -0.1269   0.0514
  74  H   2px   -0.0089   0.0014   0.0001  -0.0054  -0.0037  -0.0045   0.0132
  75  H   2py    0.0147   0.0012  -0.0057   0.0071   0.0045   0.0134   0.0015
  76  H   2pz   -0.0031   0.0027   0.0033   0.0018  -0.0007  -0.0030   0.0020
  77  H   1s     0.2146   0.1566  -0.0480   0.0876  -0.0428   0.0469   0.5729
  78  H   1s    -0.0830  -0.0578   0.0258  -0.0443   0.0331  -0.0002  -0.1572
  79  H   2px   -0.0081  -0.0044  -0.0019  -0.0037   0.0018   0.0049  -0.0064
  80  H   2py   -0.0131  -0.0062   0.0023  -0.0049  -0.0012  -0.0059  -0.0168
  81  H   2pz   -0.0035  -0.0029   0.0016   0.0033   0.0024   0.0003  -0.0115

 Orbital         15       16       17
   1  C   1s    -0.0007  -0.0035  -0.0076
   2  C   1s     0.1085  -0.1203  -0.0316
   3  C   1s     0.0254   0.0218   0.0642
   4  C   2px   -0.0297   0.0026  -0.1220
   5  C   2py    0.1918   0.1222   0.2826
   6  C   2pz    0.2537  -0.2716  -0.3399
   7  C   2px    0.0146  -0.0609   0.0519
   8  C   2py   -0.0687  -0.0588  -0.0976
   9  C   2pz   -0.0613   0.0772   0.0974
  10  C   3d2-  -0.0216   0.0046  -0.0253
  11  C   3d1-  -0.0250  -0.0140  -0.0142
  12  C   3d0    0.0248  -0.0248  -0.0238
  13  C   3d1+   0.0074   0.0012   0.0040
  14  C   3d2+  -0.0060   0.0300  -0.0282
  15  C   1s     0.0019  -0.0016  -0.0005
  16  C   1s     0.0290  -0.0497  -0.0864
  17  C   1s    -0.0039  -0.0246  -0.0651
  18  C   2px    0.2955   0.1561   0.0102
  19  C   2py   -0.3721  -0.0452  -0.0897
  20  C   2pz   -0.2950   0.2983   0.3182
  21  C   2px   -0.0605   0.0061  -0.0181
  22  C   2py    0.0675   0.0008   0.0280
  23  C   2pz    0.0602  -0.0520  -0.0408
  24  C   3d2-  -0.0138  -0.0014  -0.0046
  25  C   3d1-  -0.0129  -0.0051  -0.0108
  26  C   3d0    0.0133  -0.0151  -0.0160
  27  C   3d1+   0.0087   0.0016   0.0027
  28  C   3d2+  -0.0119  -0.0063  -0.0006
  29  F   1s     0.0009  -0.0023   0.0009
  30  F   1s    -0.0090   0.0121  -0.0182
  31  F   1s    -0.0154   0.0280  -0.0202
  32  F   2px    0.0532  -0.1795   0.1979
  33  F   2py   -0.1541  -0.2669  -0.1128
  34  F   2pz   -0.3463   0.1580   0.2076
  35  F   2px    0.0011  -0.0140   0.0181
  36  F   2py   -0.0095  -0.0149  -0.0104
  37  F   2pz   -0.0273   0.0161   0.0187
  38  F   3d2-   0.0000   0.0016  -0.0015
  39  F   3d1-  -0.0005  -0.0013   0.0001
  40  F   3d0    0.0002  -0.0018  -0.0004
  41  F   3d1+   0.0032  -0.0012   0.0012
  42  F   3d2+  -0.0015   0.0036  -0.0019
  43  O   1s     0.0018   0.0000  -0.0010
  44  O   1s    -0.0271   0.1048   0.0651
  45  O   1s    -0.0360   0.0525   0.0451
  46  O   2px   -0.1295  -0.6626   0.4150
  47  O   2py   -0.2147   0.1155  -0.4208
  48  O   2pz   -0.0951   0.2460   0.1725
  49  O   2px   -0.0014  -0.0435   0.0325
  50  O   2py   -0.0022  -0.0010  -0.0428
  51  O   2pz    0.0099  -0.0080  -0.0065
  52  O   3d2-  -0.0064  -0.0097   0.0013
  53  O   3d1-   0.0013  -0.0005   0.0017
  54  O   3d0    0.0069  -0.0052  -0.0053
  55  O   3d1+   0.0007   0.0057  -0.0052
  56  O   3d2+   0.0021  -0.0075   0.0061
  57  H   1s     0.1595   0.1247   0.5006
  58  H   1s    -0.0348  -0.0218  -0.0980
  59  H   2px    0.0015   0.0057   0.0065
  60  H   2py   -0.0009  -0.0030  -0.0137
  61  H   2pz    0.0078  -0.0011   0.0024
  62  H   1s     0.1110  -0.2033  -0.1622
  63  H   1s    -0.0182   0.0773   0.0524
  64  H   2px   -0.0024  -0.0199   0.0127
  65  H   2py   -0.0056   0.0043  -0.0119
  66  H   2pz    0.0078  -0.0120  -0.0095
  67  H   1s    -0.3979  -0.0514   0.0672
  68  H   1s     0.0944   0.0480  -0.0207
  69  H   2px   -0.0142   0.0010   0.0045
  70  H   2py   -0.0082  -0.0017  -0.0021
  71  H   2pz    0.0029   0.0047   0.0066
  72  H   1s     0.3892   0.1961   0.1818
  73  H   1s    -0.1021  -0.0463   0.0171
  74  H   2px   -0.0032  -0.0017  -0.0045
  75  H   2py    0.0093   0.0080   0.0053
  76  H   2pz   -0.0125   0.0022  -0.0004
  77  H   1s    -0.2661   0.1253  -0.0035
  78  H   1s     0.0546  -0.0317   0.0096
  79  H   2px    0.0136  -0.0010   0.0001
  80  H   2py    0.0039  -0.0071  -0.0033
  81  H   2pz    0.0013   0.0035   0.0063



 >>>> Total CPU  time used in SIRIUS :    113.58 seconds
 >>>> Total wall time used in SIRIUS :    118.00 seconds

 
     Date and time (Linux)  : Fri Mar 24 13:27:33 2006
     Host name              : star.chem.uit.no                        

- End of Wave Function Section



    *****************************************************************
    ******** Output from **PROPE input processing for ABACUS ********
    *****************************************************************



 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

      Electronic excitation energies 
      Electronic circular dichroism
      Natural orbital connection is used
      for perturbation dependent basis sets.


 Default print level:        0



 Changes of defaults for EXCITA:
 -------------------------------

 Number of excitation energies:    3    0    0    0    0    0    0    0
 Print level in EXCITA        :    0
 Integral print level in EXCITA        :    0
 Threshold in EXCITA          : 1.00E-04
 Maximum iterations in EXCITA :   60

 Center of mass gauge origin:    0.302504   -0.528873    0.143931

Starting in Static Property Section -


 
     Date and time (Linux)  : Fri Mar 24 13:27:33 2006
     Host name              : star.chem.uit.no                        

 >>> Time used in EXCITA is  2 minutes 37 seconds



 ***************************************************************************
 ************************ FINAL RESULTS FROM ABACUS ************************
 ***************************************************************************


 
     Date and time (Linux)  : Fri Mar 24 13:31:02 2006
     Host name              : star.chem.uit.no                        



                         Molecular geometry (au)
                         -----------------------

 C         -0.1077810630            0.0016311448           -0.0277533647
 C         -0.0509885280           -0.0063400396            2.8152173831
 F          2.3195193068            0.0622555881           -0.9275026109
 O         -1.2729658590           -2.1533608644           -0.8603501924
 H         -1.0158622780            1.6763033629           -0.7793105516
 H         -1.2942535134           -2.1745591590           -2.6533440159
 H         -1.9586539793           -0.0071208309            3.5572943234
 H          0.9181052495           -1.6875960201            3.4664685991
 H          0.9386804771            1.6487921657            3.5023462928





                    Molecular wave function and energy
                    ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0


     Total energy       -252.9785637871 au (Hartrees)
                         -6883.89670076 eV
                           -664195.1238 kJ/mol




                         Relativistic corrections
                         ------------------------

     Darwin correction:                          0.6428392565 au
     Mass-velocity correction:                  -0.8098212371 au

     Total relativistic correction:             -0.1669819806 au (0.0660%)
     Non-relativistic + relativistic energy:  -253.1455457677 au




                              Dipole moment
                              -------------



                    0.761783 au           1.936260 Debye




                         Dipole moment components
                         ------------------------

                               au             Debye

                    x     -0.60130378     -1.52836158
                    y      0.42770972      1.08712954
                    z     -0.18924058     -0.48100152

                        1 a.u. =   2.54175 Debye 






                Singlet electronic excitation energies
                --------------------------------------

               Sym.   Mode   Frequency    Frequency
              ex. st.  No.      (au)          (eV)
              ---------------------------------------
                 1        1    0.390777   10.633579
                 1        2    0.432720   11.774898
                 1        3    0.460931   12.542565



                Electric transition dipole moments (in a.u.)
                --------------------------------------------

  Sym.   Mode    Frequency       Velocity/Frequency              Length       
 ex. st.  No.      (au)          x       y       z         x       y       z  
 ------------------------------------------------------------------------------
   1        1     0.390777      0.000   0.000   0.000    -0.073  -0.171   0.103
   1        2     0.432720      0.000   0.000   0.000    -0.405   0.723  -0.263
   1        3     0.460931      0.000   0.000   0.000     0.001  -0.903  -0.115



                Magnetic transition dipole moments (in a.u.)
                --------------------------------------------

  Sym.   Mode    Frequency        Conventional                 London       
 ex. st.  No.      (eV)          x       y       z         x       y       z  
 ------------------------------------------------------------------------------
   1        1    10.633579      0.253   0.152  -0.220     0.227   0.117  -0.173
   1        2    11.774898     -0.181   0.074  -0.150    -0.201   0.024  -0.124
   1        3    12.542565      1.056  -0.154   0.082     1.015  -0.136   0.090



                   Oscillator and Scalar Rotational Strengths
                   ------------------------------------------

  Units: 10**(-40) (esu**2)*(cm**2) (rotational strength)
         dimensionless              (oscillator strength)

  Sym.   Mode     Frequency       Oscillator-strength      Rotational-strength   
 ex. st.  No.        (eV)          velocity  length     velocity  length   London
 -------------------------------------------------------------------------------
   1        1     10.633579        0.0000   0.0118       0.000  -31.644  -25.686
   1        2     11.774898        0.0000   0.2183       0.000   78.550   61.905
   1        3     12.542565        0.0000   0.2545       0.000   61.475   53.517


   Interatomic separations (in Angstroms):
   ---------------------------------------

            C           C           F           O           H           H     
            ------      ------      ------      ------      ------      ------
 C     :    0.000000
 C     :    1.504741    0.000000
 F     :    1.370253    2.344677    0.000000
 O     :    1.369216    2.343528    2.233818    0.000000
 H     :    1.083713    2.161409    1.962376    2.031586    0.000000
 H     :    1.910707    3.181761    2.427374    0.948945    2.271067    0.000000
 H     :    2.135043    1.083181    3.280087    2.624212    2.511720    3.498628
 H     :    2.124336    1.083185    2.610347    2.578298    3.043698    3.453212
 H     :    2.134482    1.083335    2.594986    3.278374    2.490710    4.012564

            H           H           H     
            ------      ------      ------
 H     :    0.000000
 H     :    1.763676    0.000000
 H     :    1.766185    1.765676    0.000000


  Max interatomic separation is    4.0126 Angstroms
  between atoms "H     " and "H     ".


  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C          C            1.504741
  bond distance:  F          C            1.370253
  bond distance:  O          C            1.369216
  bond distance:  H          C            1.083713
  bond distance:  H          O            0.948945
  bond distance:  H          C            1.083181
  bond distance:  H          C            1.083185
  bond distance:  H          C            1.083335


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     C          C          F            109.193
  bond angle:     C          C          O            109.170
  bond angle:     C          C          H            112.214
  bond angle:     F          C          O            109.258
  bond angle:     F          C          H            105.604
  bond angle:     O          C          H            111.303
  bond angle:     C          C          H            110.111
  bond angle:     C          C          H            109.260
  bond angle:     C          C          H            110.057
  bond angle:     H          C          H            109.000
  bond angle:     H          C          H            109.218
  bond angle:     H          C          H            109.171
  bond angle:     C          O          H            109.690




 CPU time statistics for ABACUS
 ------------------------------

 RHSIDE     00:00:44      22 %
 EXCITA     00:02:37      78 %

 TOTAL      00:03:22     100 %



 >>>> Total CPU  time used in ABACUS:  3 minutes 22 seconds
 >>>> Total wall time used in ABACUS:  3 minutes 29 seconds

- End of Static Property Section

 >>>> Total CPU  time used in DALTON:  5 minutes 23 seconds
 >>>> Total wall time used in DALTON:  5 minutes 38 seconds

 
     Date and time (Linux)  : Fri Mar 24 13:31:02 2006
     Host name              : star.chem.uit.no                        
END REFOUT

