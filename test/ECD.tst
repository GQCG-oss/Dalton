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

     Date and time (Linux)  : Thu Sep 24 00:50:52 2009
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
    Static molecular property section will be executed (ABACUS module)
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

     ICESPH =       2     NESFP =       5
     OMEGA = 40.0000     RET   =  0.2000     FRO   =  0.7000

     IPRPCM=       0

     NON-EQ = F     NEQRSP =F
 POLYG          60

     INPUT FOR CAVITY DEFINITION 
     ---------------------------
     ATOM         COORDINATES           RADIUS 
     1    0.0000    0.0000    0.0000    1.7000
     2    0.0000    0.0000    0.0000    2.0000
     3    0.0000    0.0000    0.0000    1.4000
     4    0.0000    0.0000    0.0000    1.6000
     5    0.0000    0.0000    0.0000    1.2000


   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************

    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1:                                                                         
 2:                                                                         
    ------------------------------------------------------------------------

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centres:    2
  Used basis set file for basis set for elements with Z =   6 :
     "/home/ruud/DaltonFix/dalton/basis/cc-pVDZ"

  Atomic type no.    2
  --------------------
  Nuclear charge:   9.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   9 :
     "/home/ruud/DaltonFix/dalton/basis/cc-pVDZ"

  Atomic type no.    3
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centres:    1
  Used basis set file for basis set for elements with Z =   8 :
     "/home/ruud/DaltonFix/dalton/basis/cc-pVDZ"

  Atomic type no.    4
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centres:    5
  Used basis set file for basis set for elements with Z =   1 :
     "/home/ruud/DaltonFix/dalton/basis/cc-pVDZ"


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 
 ********SPHERES IN PCMSPHGEN************
 INDEX        X        Y         Z        R
   1   -1.0778106300e-01    1.6311448000e-03   -2.7753364700e-02    1.7000000000e+00
   2   -5.0988528000e-02   -6.3400396000e-03    2.8152173831e+00    2.0000000000e+00
   3    2.3195193068e+00    6.2255588100e-02   -9.2750261090e-01    1.4000000000e+00
   4   -1.2729658590e+00   -2.1533608644e+00   -8.6035019240e-01    1.6000000000e+00
   5   -1.0158622780e+00    1.6763033629e+00   -7.7931055160e-01    1.2000000000e+00


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

 Center-of-mass coordinates (a.u.):    0.302504   -0.528873    0.143931


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

  Threshold for integrals:  1.00e-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   27
  C       :    1  x  -0.1077810630   2  y   0.0016311448   3  z  -0.0277533647
  C       :    4  x  -0.0509885280   5  y  -0.0063400396   6  z   2.8152173831
  F       :    7  x   2.3195193068   8  y   0.0622555881   9  z  -0.9275026109
  O       :   10  x  -1.2729658590  11  y  -2.1533608644  12  z  -0.8603501924
  H       :   13  x  -1.0158622780  14  y   1.6763033629  15  z  -0.7793105516
  H       :   16  x  -1.2942535134  17  y  -2.1745591590  18  z  -2.6533440159
  H       :   19  x  -1.9586539793  20  y  -0.0071208309  21  z   3.5572943234
  H       :   22  x   0.9181052495  23  y  -1.6875960201  24  z   3.4664685991
  H       :   25  x   0.9386804771  26  y   1.6487921657  27  z   3.5023462928


   Interatomic separations (in Angstrom):
   --------------------------------------

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


  Max interatomic separation is    4.0126 Angstrom (    7.5826 Bohr)
  between atoms    9 and    6, "H     " and "H     ".


  Bond distances (Angstrom):
  --------------------------

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

   IA      52.248587         -0.216669    0.209580    0.953483
   IB      56.452522          0.841008    0.536036    0.073287
   IC      96.963569         -0.495741    0.817766   -0.292401


 Rotational constants
 --------------------

               A                   B                   C

           9672.5871           8952.2839           5212.0504 MHz
            0.322643            0.298616            0.173855 cm-1


  Nuclear repulsion energy :  135.226836327282


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************



     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00e-15

 Number of two-electron integrals written:     5495979 ( 99.6% )
 Megabytes written:                             62.932


 >>> Time used in TWOINT is   3.22 seconds


 MEMORY USED TO GENERATE CAVITY =    432042

 Tessera cut in pieces and removed.

 ** WARNING ** A VERY POOR TESSELATION HAS BEEN CHOSEN
 IT IS VALUABLE ALMOST ONLY FOR TESTING

 Total number of spheres =    5
 Sphere             Center  (X,Y,Z) (A)               Radius (A)      Area (A^2)
   1   -0.057035282    0.000863165   -0.014686448    2.040000000    3.740433668
   2   -0.026981967   -0.003355004    1.489748876    2.400000000   50.281695884
   3    1.227436751    0.032944238   -0.490813242    1.680000000   17.505819464
   4   -0.673624520   -1.139509491   -0.455277713    1.920000000   25.650918549
   5   -0.537571164    0.887061534   -0.412393382    1.440000000    8.646852065

 Total number of tesserae =     245
 Surface area =  105.82571963 (A^2)    Cavity volume =   92.06411535 (A^3)

          THE SOLUTE IS ENCLOSED IN ONE CAVITY

 ..... DONE GENERATION CAVITY .....

 >>> Time used in PEDRAM is   3.32 seconds

 
  ..... DONE GENERATING -Q-  MATRIX .....
 >>>> Total CPU  time used in HERMIT:   3.40 seconds
 >>>> Total wall time used in HERMIT:   3.00 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


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

 
     Date and time (Linux)  : Thu Sep 24 00:50:55 2009
     Host name              : stallo-2.local                          

 Title lines from ".mol" input file:
                                                                             
                                                                             

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

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species          All    1
                                       --
     Occupied SCF orbitals              17   17
     Secondary orbitals                 64   64
     Total number of orbitals           81   81
     Number of basis functions          81   81

     Optimization information
     ========================
     Number of configurations              1
     Number of orbital rotations        1088
     ---------------------------------------
     Total number of variables          1089

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00e-05

          -------------------------------------
          ---- POLARISABLE CONTINUUM MODEL ----
          ----      UNIVERSITY OF PISA     ----
          -------------------------------------

 ESTIMATE OF NUCLEAR CHARGE       33.94200
 NUCLEAR APPARENT CHARGE -33.55690
 THEORETICAL -33.56627 NOT RENORMALIZED

 ..... DONE WITH INDUCED NUCLEAR CHARGES .....


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------

 Precalculated two-electron integrals are transformed to P-supermatrix elements.
 Threshold for discarding integrals :  1.00e-15

 >>> Time used in FORMSUP is   1.29 seconds

PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.22818   225.31524   225.31425  -112.71685    -0.01620
   1  -252.188970817     -1.619748573660e-02    3.48289e+00   -2.52e+02
 MULPOP C    10.46; F     9.56; O     8.58; H     5.39; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.71634   225.56549   225.56542  -112.71685    -0.00956
   2  -252.835208674     -9.564151351569e-03    1.69062e+00   -6.46e-01
 MULPOP C    12.42; F     9.12; O     8.13; H     4.33; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.53315   225.47256   225.47204  -112.71685    -0.01112
   3  -252.972507644     -1.112171516193e-02    3.10615e-01   -1.37e-01
 MULPOP C    11.60; F     9.38; O     8.36; H     4.66; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52499   225.46824   225.46777  -112.71685    -0.01134
   4  -252.978092679     -1.134039723665e-02    7.78036e-02   -5.59e-03
 MULPOP C    11.63; F     9.37; O     8.36; H     4.65; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.53224   225.47152   225.47105  -112.71685    -0.01168
   5  -252.978518758     -1.168147166113e-02    1.91808e-02   -4.26e-04
 MULPOP C    11.63; F     9.37; O     8.37; H     4.63; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52968   225.47018   225.46971  -112.71685    -0.01175
   6  -252.978559692     -1.174539617573e-02    5.26369e-03   -4.09e-05
 MULPOP C    11.62; F     9.37; O     8.37; H     4.64; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52961   225.47011   225.46964  -112.71685    -0.01178
   7  -252.978563412     -1.177950421103e-02    1.76493e-03   -3.72e-06
 MULPOP C    11.62; F     9.37; O     8.37; H     4.64; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52949   225.47004   225.46957  -112.71685    -0.01179
   8  -252.978563741     -1.178810274756e-02    6.25957e-04   -3.29e-07
 MULPOP C    11.62; F     9.37; O     8.37; H     4.64; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52944   225.47001   225.46954  -112.71685    -0.01179
   9  -252.978563786     -1.178899550159e-02    1.01464e-04   -4.42e-08
 MULPOP C    11.62; F     9.37; O     8.37; H     4.64; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52944   225.47001   225.46955  -112.71685    -0.01179
  10  -252.978563787     -1.178916740022e-02    1.83276e-05   -1.53e-09
 MULPOP C    11.62; F     9.37; O     8.37; H     4.64; 
PCMFCK: PCMEE, PCMEN, PCMNE, PCMNN, ESOLT  -225.52944   225.47001   225.46955  -112.71685    -0.01179
  11  -252.978563787     -1.178917552987e-02    5.67530e-06   -4.31e-11
 DIIS converged in  11 iterations !


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   34
 Orbital occupations :   17

 Sym       Hartree-Fock orbital energies

  1    -26.28287800   -20.56613859   -11.35729965   -11.22167230    -1.61021837
        -1.38734423    -1.02405064    -0.85795003    -0.75329004    -0.72505419
        -0.68785063    -0.62232756    -0.56693052    -0.55035043    -0.53830698
        -0.51671865    -0.46872535     0.19813000     0.22853084     0.25326088
         0.26886604     0.27881933

    E(LUMO) :     0.19813000 au (symmetry 1)
  - E(HOMO) :    -0.46872535 au (symmetry 1)
  ------------------------------------------
    gap     :     0.66685535 au

 >>> Writing SIRIFC interface file <<<


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     SOLVATION MODEL: polarizable continuum model (PCM),
          dielectric constant =   78.390000

     Final HF energy:            -252.978563787121                 
     Nuclear repulsion:           135.226836327282
     Electronic energy:          -388.193610938872

     Final gradient norm:           0.000005675296

 
     Date and time (Linux)  : Thu Sep 24 00:51:45 2009
     Host name              : stallo-2.local                          

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

 Orbital          13       14       15       16       17       18       19
   1 C   :1s     0.0030   0.0010  -0.0007  -0.0035  -0.0076   0.0377   0.0297
   2 C   :1s     0.0034  -0.0012   0.1085  -0.1203  -0.0316  -0.0147   0.0218
   3 C   :1s    -0.0541  -0.0174   0.0254   0.0218   0.0642  -1.0067  -0.7103
   4 C   :2px   -0.1001   0.0258  -0.0297   0.0026  -0.1220   0.0822   0.0252
   5 C   :2py   -0.2134   0.0160   0.1918   0.1222   0.2826  -0.2258  -0.1674
   6 C   :2pz   -0.2267  -0.0272   0.2537  -0.2716  -0.3399   0.0867  -0.0456
   7 C   :2px    0.0422  -0.0246   0.0146  -0.0609   0.0519   0.2436   0.1574
   8 C   :2py    0.0954  -0.0012  -0.0687  -0.0588  -0.0976  -0.2340  -0.0766
   9 C   :2pz    0.0322   0.0098  -0.0613   0.0772   0.0974   0.1804   0.0756
  10 C   :3d2-   0.0221  -0.0035  -0.0216   0.0046  -0.0253   0.0006  -0.0050
  11 C   :3d1-   0.0113  -0.0144  -0.0250  -0.0140  -0.0142  -0.0015  -0.0122
  12 C   :3d0   -0.0078  -0.0006   0.0248  -0.0248  -0.0238   0.0030   0.0080
  13 C   :3d1+   0.0303  -0.0118   0.0074   0.0012   0.0040  -0.0058  -0.0006
  14 C   :3d2+  -0.0107   0.0180  -0.0060   0.0300  -0.0282   0.0085   0.0063
  15 C   :1s    -0.0021  -0.0005   0.0019  -0.0016  -0.0005   0.0544  -0.0675
  16 C   :1s     0.0049  -0.0006   0.0290  -0.0497  -0.0864  -0.0921   0.0868
  17 C   :1s     0.0286   0.0056  -0.0039  -0.0246  -0.0651  -1.3091   1.6971
  18 C   :2px    0.2970  -0.4333   0.2955   0.1561   0.0102   0.0228   0.0279
  19 C   :2py   -0.2111  -0.4280  -0.3721  -0.0452  -0.0897  -0.0549  -0.0811
  20 C   :2pz    0.1850   0.0326  -0.2950   0.2983   0.3182  -0.1519   0.1555
  21 C   :2px   -0.0320   0.0851  -0.0605   0.0061  -0.0181  -0.0096   0.0548
  22 C   :2py    0.0234   0.0668   0.0675   0.0008   0.0280  -0.0967  -0.2833
  23 C   :2pz   -0.0396  -0.0064   0.0602  -0.0520  -0.0408  -0.3260   0.2350
  24 C   :3d2-  -0.0078  -0.0157  -0.0138  -0.0014  -0.0046  -0.0004   0.0010
  25 C   :3d1-   0.0022  -0.0084  -0.0129  -0.0051  -0.0108   0.0081   0.0115
  26 C   :3d0   -0.0099  -0.0010   0.0133  -0.0151  -0.0160   0.0066  -0.0111
  28 C   :3d2+  -0.0083   0.0152  -0.0119  -0.0063  -0.0006  -0.0024  -0.0022
  30 F   :1s    -0.0051   0.0075  -0.0090   0.0121  -0.0182  -0.0500  -0.0088
  31 F   :1s     0.0016   0.0126  -0.0154   0.0280  -0.0202  -0.0170   0.0194
  32 F   :2px    0.1805  -0.0570   0.0532  -0.1795   0.1979  -0.0344   0.0225
  33 F   :2py    0.3445   0.0336  -0.1541  -0.2669  -0.1128   0.0390   0.0251
  34 F   :2pz    0.1997   0.2245  -0.3463   0.1580   0.2076  -0.0100   0.0799
  35 F   :2px    0.0072  -0.0010   0.0011  -0.0140   0.0181  -0.0228  -0.0161
  36 F   :2py    0.0182   0.0028  -0.0095  -0.0149  -0.0104   0.0293   0.0158
  37 F   :2pz    0.0168   0.0158  -0.0273   0.0161   0.0187  -0.0230   0.0198
  43 O   :1s    -0.0015   0.0002   0.0018   0.0000  -0.0010   0.0089   0.0121
  44 O   :1s    -0.1170  -0.0739  -0.0271   0.1048   0.0651  -0.1310  -0.0886
  45 O   :1s    -0.0365  -0.0353  -0.0360   0.0525   0.0451  -0.2699  -0.2853
  46 O   :2px    0.1985  -0.1166  -0.1295  -0.6626   0.4150  -0.0020  -0.0112
  47 O   :2py    0.4703   0.0680  -0.2147   0.1155  -0.4208   0.0440  -0.0365
  48 O   :2pz   -0.1599  -0.1651  -0.0951   0.2460   0.1725   0.1077   0.1909
  49 O   :2px    0.0021   0.0001  -0.0014  -0.0435   0.0325  -0.0006   0.0168
  50 O   :2py    0.0161   0.0120  -0.0022  -0.0010  -0.0428   0.1052   0.0824
  51 O   :2pz   -0.0013   0.0048   0.0099  -0.0080  -0.0065   0.1219   0.1990
  57 H   :1s    -0.0586   0.0031   0.1595   0.1247   0.5006  -0.0185   0.0620
  58 H   :1s     0.0161  -0.0014  -0.0348  -0.0218  -0.0980   1.0702   0.6750
  60 H   :2py   -0.0001   0.0008  -0.0009  -0.0030  -0.0137  -0.0165  -0.0159
  61 H   :2pz   -0.0032   0.0001   0.0078  -0.0011   0.0024   0.0020  -0.0101
  62 H   :1s     0.0761   0.1355   0.1110  -0.2033  -0.1622   0.0904   0.0505
  63 H   :1s    -0.0624  -0.0558  -0.0182   0.0773   0.0524   0.6406   0.7949
  64 H   :2px    0.0038  -0.0040  -0.0024  -0.0199   0.0127  -0.0005  -0.0034
  65 H   :2py    0.0133   0.0023  -0.0056   0.0043  -0.0119   0.0010  -0.0041
  66 H   :2pz    0.0044   0.0089   0.0078  -0.0120  -0.0095   0.0142   0.0194
  67 H   :1s    -0.2100   0.4407  -0.3979  -0.0514   0.0672   0.0039  -0.0132
  68 H   :1s     0.0620  -0.1159   0.0944   0.0480  -0.0207   0.7529  -0.6928
  69 H   :2px   -0.0048   0.0138  -0.0142   0.0010   0.0045   0.0107  -0.0091
  72 H   :1s     0.3797   0.1600   0.3892   0.1961   0.1818  -0.0166  -0.0824
  73 H   :1s    -0.1269  -0.0514  -0.1021  -0.0463   0.0171   0.3929  -1.2463
  74 H   :2px   -0.0045  -0.0132  -0.0032  -0.0017  -0.0045  -0.0025   0.0076
  75 H   :2py    0.0134  -0.0015   0.0093   0.0080   0.0053   0.0033  -0.0159
  76 H   :2pz   -0.0030  -0.0020  -0.0125   0.0022  -0.0004  -0.0043   0.0022
  77 H   :1s     0.0469  -0.5729  -0.2661   0.1253  -0.0035   0.0517   0.0633
  78 H   :1s    -0.0002   0.1572   0.0546  -0.0317   0.0096   0.8172  -0.5080
  79 H   :2px    0.0049   0.0064   0.0136  -0.0010   0.0001  -0.0057   0.0031
  80 H   :2py   -0.0059   0.0168   0.0039  -0.0071  -0.0033  -0.0076   0.0093
  81 H   :2pz    0.0003   0.0115   0.0013   0.0035   0.0063  -0.0032   0.0048



 >>>> Total CPU  time used in SIRIUS :     45.47 seconds
 >>>> Total wall time used in SIRIUS :     50.00 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:51:45 2009
     Host name              : stallo-2.local                          


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
      Electronic circular dichroism
      Natural orbital connection is used
      for perturbation dependent basis sets.


 Changes of defaults for EXCITA:
 -------------------------------

 Number of excitation energies:    3    0    0    0    0    0    0    0
 Print level in EXCITA        :    0
 Integral print level in EXCITA        :    0
 Threshold in EXCITA          : 1.00e-04
 Maximum iterations in EXCITA :   60

 Center of mass dipole origin  :    0.302504   -0.528873    0.143931

 Center of mass gauge origin   :    0.302504   -0.528873    0.143931


                 .------------------------------------------------.
                 | Starting in Static Property Section (ABACUS) - |
                 `------------------------------------------------'


 
     Date and time (Linux)  : Thu Sep 24 00:51:45 2009
     Host name              : stallo-2.local                          

 >>> Time used in EXCITA is  1 minute  35 seconds



   ***************************************************************************
   ************************ FINAL RESULTS from ABACUS ************************
   ***************************************************************************


 
     Date and time (Linux)  : Thu Sep 24 00:53:48 2009
     Host name              : stallo-2.local                          


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

     Darwin correction:                          0.6428392553 au
     Mass-velocity correction:                  -0.8098212351 au

     Total relativistic correction:             -0.1669819798 au (0.0660%)
     Non-relativistic + relativistic energy:  -253.1455457669 au




                                  Dipole moment
                                  -------------

                 au               Debye          C m (/(10**-30)
              0.761784           1.936263           6.458677


                             Dipole moment components
                             ------------------------

                 au               Debye          C m (/(10**-30)

      x     -0.60130351        -1.52836090        -5.09806321
      y      0.42771116         1.08713321         3.62628604
      z     -0.18924239        -0.48100612        -1.60446370


   Units:   1 a.u. =   2.54175 Debye 
            1 a.u. =   8.47835 (10**-30) C m (SI)




                Singlet electronic excitation energies
                --------------------------------------

                 Sym.   Mode   Frequency    Frequency
                ex. st.  No.      (au)          (eV)
                =======================================
                   1        1    0.390777   10.633580
                   1        2    0.432720   11.774899
                   1        3    0.460931   12.542562
                ---------------------------------------


                Electric transition dipole moments (in a.u.)
                --------------------------------------------

  Sym.   Mode    Frequency       Velocity/Frequency              Length
 ex. st.  No.      (au)          x       y       z         x       y       z
 ==============================================================================
   1        1     0.390777     0.0000  0.0000  0.0000   -0.0729 -0.1708  0.1032
   1        2     0.432720     0.0000  0.0000  0.0000   -0.4053  0.7234 -0.2630
   1        3     0.460931     0.0000  0.0000  0.0000    0.0010 -0.9028 -0.1152
 ------------------------------------------------------------------------------


                Magnetic transition dipole moments (in a.u.)
                --------------------------------------------

  ( mu_B*<0|l_i|n>, where mu_B = 0.5 is the Bohr magneton)

 Symm.  Mode   Frequency         Conventional                  London
 ex.st.  No.     (eV)        x        y       z         x        y        z
 ==============================================================================
   1      1    10.6336    0.25309  0.15191 -0.22017   0.22685  0.11745 -0.17328
   1      2    11.7749   -0.18116  0.07438 -0.14971  -0.20123  0.02385 -0.12354
   1      3    12.5426    1.05595 -0.15374  0.08187   1.01545 -0.13613  0.09004
 ------------------------------------------------------------------------------


                   Oscillator and Scalar Rotational Strengths
                   ------------------------------------------

  Units: 10**(-40) (esu**2)*(cm**2) (rotational strength)
         dimensionless              (oscillator strength)

 Symm.  Mode   Frequency    Oscillator-strength      Rotational-strength   
 ex.st.  No.     (eV)        velocity  length     velocity  length   London
 ==========================================================================
   1      1    10.6336        0.0000   0.0118       0.000  -31.644  -25.686
   1      2    11.7749        0.0000   0.2183       0.000   78.550   61.905
   1      3    12.5426        0.0000   0.2545       0.000   61.478   53.520
 --------------------------------------------------------------------------


   Interatomic separations (in Angstrom):
   --------------------------------------

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


  Max interatomic separation is    4.0126 Angstrom (    7.5826 Bohr)
  between atoms    9 and    6, "H     " and "H     ".


  Bond distances (Angstrom):
  --------------------------

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

 RHSIDE     00:00:18      16 %
 EXCITA     00:01:35      84 %

 TOTAL      00:01:53     100 %


 >>>> Total CPU  time used in ABACUS:  1 minute  53 seconds
 >>>> Total wall time used in ABACUS:  2 minutes  3 seconds


                   .-------------------------------------------.
                   | End of Static Property Section (ABACUS) - |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  2 minutes 42 seconds
 >>>> Total wall time used in DALTON:  2 minutes 56 seconds

 
     Date and time (Linux)  : Thu Sep 24 00:53:48 2009
     Host name              : stallo-2.local                          
END REFOUT

