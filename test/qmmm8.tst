########## Test description ########################
START DESCRIPTION
KEYWORDS: qmmm hf dipole properties
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
qmmmaniso
qmmmconvergence
qmmmq
qmmmdip
qmmmquad
qmmmoct
qmmmelpol
qmmmnucpol
qmmmmulpol
qm3energy
enehf
dipole
dipcompx
dipcompy
dipcompz
excita
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON
.RUN PROPERTIES
*QMMM
.QMMM
.PRINT
 1 
.MMPROP
**WAVE FUNCTIONS
.HF 
*SCF INPUT
.THRESHOLD
 1.0D-10
**PROPERTIES
.EXCITA
*EXCITA
.NEXCIT
 1
.DIPSTR
.THRESH
 1.0D-7
**END OF
END DALINP

########## MOLECULE.INP ############################
START MOLINP
ATOMBASIS
Acrolein plus 2 water molecules
------------------------
AtomTypes=3 NoSymmetry Angstrom
        6.0   3    Basis=cc-pVDZ
C             -0.145335   -0.546770    0.000607
C              1.274009   -0.912471   -0.000167
C              1.630116   -2.207690   -0.000132
        8.0   1    Basis=cc-pVDZ
O             -0.560104    0.608977    0.000534
        1.0   4    Basis=cc-pVDZ
H             -0.871904   -1.386459    0.001253
H              2.004448   -0.101417   -0.000710
H              0.879028   -3.000685    0.000484
H              2.675323   -2.516779   -0.000673
END MOLINP

########## POTENTIAL.INP ############################
START POTINP
AU
10 3 2 1
1       -6.2900497831       -0.1950764279       -0.0007861261       -0.7423294438        0.0303539300        0.3280749969        0.0000894127       -3.9516312016       -0.0561791973        0.0008348984       -4.5778807726        0.0000430036       -5.0206878337       -0.0150265431       -0.0987395310       -0.0000662324        0.0267925411       -0.0002404801        0.0219905017        0.0316173819       -0.0000537894        0.2063402594        0.0001754735        1.5931619660        0.0800886126       -0.0012113893        2.5255643422        0.0013486011        3.3670040387
1       -4.7314868231        0.7808745191        0.0006406172        0.3699305781       -0.1005787187       -0.0557286568       -0.0000906399       -0.5774205281       -0.0533068952       -0.0000246837       -0.6043997755       -0.0000510833       -0.5595888419       -0.1619181713        0.0850921795        0.0000929568       -0.1146181832        0.0000703947       -0.1289664836        0.0249611871       -0.0000856861       -0.0002388295       -0.0003063695        0.7929865198        0.1542303414        0.0001503753        0.6015087761        0.0002890322        0.5925152779
1       -7.6330082202        1.0331680746       -0.0016062672        0.3723988657        0.0911510444       -0.0729658327        0.0000575297       -0.5583105371        0.0466697048        0.0000093367       -0.6224134008        0.0000240992       -0.5582982578        0.2225477804        0.0892957698       -0.0000119240        0.0847835272        0.0000788789        0.1297854415       -0.0229222402        0.0000709112       -0.0215875620        0.0002947454        0.7209493065       -0.1783490996        0.0000586528        0.6420261847        0.0000623389        0.5751449171
1       -5.5107683031        0.2928990456       -0.0000727545        0.0000000000       -0.1154554696       -0.1098625035       -0.0001131928        0.6931173508        0.3993825482        0.0004427892        0.4811735442        0.0003651060        0.2335453026        1.1456775180        0.1549332625        0.0007929257        0.4485264365        0.0002257656        0.0913068047        0.1384100831        0.0003985145       -0.0411780894        0.0001923295        3.4974039172        2.1355046659        0.0022286584        1.8455550714        0.0014941130        1.4124398311
1       -6.9615290017        0.4190458234       -0.0011961966        0.0000000000        0.0925537150       -0.1281941897        0.0000470915        0.5497923459       -0.4072157699        0.0001632249        0.6322989679       -0.0002450231        0.2419436597       -1.0679926401        0.1883003226       -0.0006460315       -0.4727878386        0.0000929861       -0.1029722061        0.3945907516       -0.0002891323       -0.0214685564       -0.0002611650        2.6916145417       -2.2463816416        0.0003377505        2.5542051973       -0.0013254304        1.4293331897
2        3.2924641584        4.4245310497       -0.0014078460       -0.7424407021       -0.2840371815        0.1671869101        0.0013661730       -4.4188113058        0.2801417448        0.0002126277       -4.1129942280        0.0039064597       -5.0206529493        0.0330216712        0.0760732530       -0.0003899377        0.0087745562        0.0009840165       -0.1813726316       -0.1218544851       -0.0016261873        0.1033223765        0.0025668497        2.2823444229       -0.4207398269       -0.0006457588        1.8324399300       -0.0069062619        3.3666855109
2        1.5905409054        3.7261752013       -0.0015495754        0.3699635356        0.1039146227        0.0490621032        0.0000401115       -0.6455428297       -0.0041180935        0.0001453595       -0.5362599824        0.0001224056       -0.5596286658        0.0810557396        0.0035896134        0.0001532545       -0.0090165544       -0.0002204117        0.0757377031        0.2842212967        0.0003840935        0.1046774555        0.0011004803        0.8132771942        0.1387943320        0.0004103386        0.5814638105       -0.0000004496        0.5927954021
2        3.0850893929        6.2328855678        0.0078499223        0.3724771666        0.0058194611       -0.1165815150       -0.0005822607       -0.5562677027        0.0450154899        0.0002165578       -0.6244100291       -0.0004104801       -0.5583164804       -0.2048195849       -0.1193833899       -0.0004155498       -0.0283899604        0.0001497946       -0.0584957773       -0.0907038156        0.0007682659       -0.1181777112       -0.0013985302        0.4993246061       -0.0196042351        0.0000571000        0.8612729797        0.0020132268        0.5750737753
2        2.4415025319        4.0753531255       -0.0014787107        0.0000000000        0.1562106237        0.0289341838       -0.0001592549        0.9309922644        0.2280573432       -0.0002276453        0.2429875998       -0.0003906527        0.2330344430       -0.9479164549       -0.3465723180        0.0000503448       -0.2208226146       -0.0005070175       -0.0203355542       -0.7822646206       -0.0024267913       -0.1000074397       -0.0013521451        4.2946372081        1.2694080541        0.0004092117        1.0566337570       -0.0044943748        1.4133256363
2        3.1887767756        5.3287083088        0.0032210382        0.0000000000        0.0501665747       -0.1493857957       -0.0008283012        0.2176177849       -0.1663997493       -0.0007688598        0.9642929866        0.0038606681        0.2416297737        0.3960654724        0.3304373259        0.0011603890        0.0383098210       -0.0006969579        0.0790638970        1.2632196844        0.0054054256        0.0723757018        0.0006068808        0.6174280837       -0.4402546796       -0.0002937434        4.6220075983        0.0174865911        1.4307100148
END POTINP

########## Reference Output ########################
START REFOUT


   ****************************************************************************
   *************** DALTON2011 - An electronic structure program ***************
   ****************************************************************************

    This is output from DALTON Release 2011 (Rev. 0, Mar. 2011)
   ----------------------------------------------------------------------------
    NOTE:
     
    This is an experimental code for the evaluation of molecular
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
    Release DALTON2011 (2011), see http://daltonprogram.org"
 --------------------------------------------------------------------------------
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
  Luca Frediani,            University of Tromsoe,        Norway      (PCM)
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
  Kurt V. Mikkelsen,        University of Copenhagen,     Denmark     (MC-SCRF and QM/MM code)
  Christian B. Nielsen,     University of Copenhagen,     Denmark     (QM/MM code)
  Patrick Norman,           University of Linkoeping,     Sweden      (cubic response and complex response in RESPONS)
  Jeppe Olsen,              Aarhus University,            Denmark     (SIRIUS CI/density modules)
  Anders Osted,             Copenhagen University,        Denmark     (QM/MM code)
  Martin J. Packer,         University of Sheffield,      UK          (SOPPA)
  Thomas B. Pedersen,       University of Oslo,           Norway      (Cholesky decomposition)
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
  Olav Vahtras,             KTH Stockholm,                Sweden      (triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)
  David J. Wilson,          La Trobe University,          Australia   (DFT Hessian and DFT magnetizabilities)
  Hans Agren,               KTH Stockholm,                Sweden      (SIRIUS module, MC-SCRF solvation model)
 --------------------------------------------------------------------------------

     Date and time (Linux)  : Thu Apr  7 15:36:38 2011 
     Host name              : stanley                                 

 * Work memory size             :   100000000 =  762.94 megabytes.

 * Directories for basis set searches:
   1) /home/arnfinn/jobb/dalton/svn/pure_trunk/test/2011-04-07T15_29-testjob-pid-10167
   2) /home/arnfinn/jobb/dalton/svn/pure_trunk/basis/


       *******************************************************************
       *********** Output from DALTON general input processing ***********
       *******************************************************************

 --------------------------------------------------------------------------------
   Overall default print level:    0
   Print level for DALTON.STAT:    1

    HERMIT 1- and 2-electron integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed (SIRIUS module)
    Static molecular property section will be executed (ABACUS module)
 --------------------------------------------------------------------------------


 Changes of defaults for *QMMM  :
 --------------------------------

 +------------------+
 |  WORD: | CHANGE: |
 +------------------+
 |   QMMM |       T |
 | MMPROP |       T |
 +------------------+



   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************


    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1: Acrolein plus 2 water molecules                                         
 2: ------------------------                                                
    ------------------------------------------------------------------------

  Coordinates are entered in Angstrom and converted to atomic units.
          - Conversion factor : 1 bohr = 0.52917721 A

  Atomic type no.    1
  --------------------
  Nuclear charge:   6.00000
  Number of symmetry independent centers:    3
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Used basis set file for basis set for elements with Z =   6 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/cc-pVDZ"

  Atomic type no.    2
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Used basis set file for basis set for elements with Z =   8 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/cc-pVDZ"

  Atomic type no.    3
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    4
  Number of basis sets to read;    2
  The basis set is "cc-pVDZ" from the basis set library.
  Used basis set file for basis set for elements with Z =   1 :
     "/home/arnfinn/jobb/dalton/svn/pure_trunk/basis/cc-pVDZ"


                         SYMGRP: Point group information
                         -------------------------------

Point group: C1 


                                 Isotopic Masses
                                 ---------------

                           C          12.000000
                           C          12.000000
                           C          12.000000
                           O          15.994915
                           H           1.007825
                           H           1.007825
                           H           1.007825
                           H           1.007825

                       Total mass:    56.026215 amu
                       Natural abundance:  96.446 %

 Center-of-mass coordinates (a.u.):    0.973773   -1.393790    0.000425


  Atoms and basis sets
  --------------------

  Number of atom types :    3
  Total number of atoms:    8

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  C           3    6.0000    26    14      [9s4p1d|3s2p1d]                                    
  O           1    8.0000    26    14      [9s4p1d|3s2p1d]                                    
  H           4    1.0000     7     5      [4s1p|2s1p]                                        
  ----------------------------------------------------------------------
  total:      8   30.0000   132    76
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00D-15


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:   24
  C       :     1  x  -0.2746433464    2  y  -1.0332455534    3  z   0.0011470638
  C       :     4  x   2.4075280908    5  y  -1.7243202870    6  z  -0.0003155843
  C       :     7  x   3.0804727920    8  y  -4.1719294689    9  z  -0.0002494438
  O       :    10  x  -1.0584431615   11  y   1.1507997464   12  z   0.0010091138
  H       :    13  x  -1.6476597673   14  y  -2.6200277935   15  z   0.0023678268
  H       :    16  x   3.7878577518   17  y  -0.1916503544   18  z  -0.0013417055
  H       :    19  x   1.6611221762   20  y  -5.6704728374   21  z   0.0009146274
  H       :    22  x   5.0556277659   23  y  -4.7560230271   24  z  -0.0012717857


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           C           C           O           H           H     
            ------      ------      ------      ------      ------      ------
 C     :    0.000000
 C     :    1.465700    0.000000
 C     :    2.431231    1.343281    0.000000
 O     :    1.227919    2.383018    3.568007    0.000000
 H     :    1.110397    2.197637    2.633349    2.019650    0.000000
 H     :    2.195429    1.091490    2.139278    2.661125    3.150355    0.000000
 H     :    2.659139    2.125241    1.092234    3.885970    2.381489    3.110036
 H     :    3.440501    2.130137    1.089951    4.498704    3.722962    2.506800

            H           H     
            ------      ------
 H     :    0.000000
 H     :    1.860334    0.000000


  Max interatomic separation is    4.4987 Angstrom (    8.5013 Bohr)
  between atoms    8 and    4, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C          C            1.465700
  bond distance:  C          C            1.343281
  bond distance:  O          C            1.227919
  bond distance:  H          C            1.110397
  bond distance:  H          C            1.091490
  bond distance:  H          C            1.092234
  bond distance:  H          C            1.089951


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     C          C          O            124.190
  bond angle:     C          C          H            116.421
  bond angle:     O          C          H            119.389
  bond angle:     C          C          C            119.821
  bond angle:     C          C          H            117.558
  bond angle:     C          C          H            122.621
  bond angle:     C          C          H            121.182
  bond angle:     C          C          H            121.847
  bond angle:     H          C          H            116.971




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA      10.696865         -0.668732    0.743503    0.000230
   IB     108.794831          0.743503    0.668732   -0.000587
   IC     119.491696          0.000590    0.000221    1.000000


 Rotational constants
 --------------------

 The molecule is planar.

               A                   B                   C

          47245.5253           4645.2483           4229.4069 MHz
            1.575941            0.154949            0.141078 cm-1


@  Nuclear repulsion energy :  102.664846860515


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************


  *********************************** 
  QMMM electrostatic potential: 
  Multipole order                          3
  Anisotropic polarization
  *********************************** 

  ---------------- 
  QMMM information 
  ---------------- 

  MM coordinates in au 
  -------------------- 
   1     -6.290050     -0.195076     -0.000786
   2     -4.731487      0.780875      0.000641
   3     -7.633008      1.033168     -0.001606
   4     -5.510768      0.292899     -0.000073
   5     -6.961529      0.419046     -0.001196
   6      3.292464      4.424531     -0.001408
   7      1.590541      3.726175     -0.001550
   8      3.085089      6.232886      0.007850
   9      2.441503      4.075353     -0.001479
  10      3.188777      5.328708      0.003221

  MM charges 
  ---------- 
   1     -0.742329
   2      0.369931
   3      0.372399
   4      0.000000
   5      0.000000
   6     -0.742441
   7      0.369964
   8      0.372477
   9      0.000000
  10      0.000000

  MM dipoles (x,y,z) 
  ------------------ 
   1      0.030354      0.328075      0.000089
   2     -0.100579     -0.055729     -0.000091
   3      0.091151     -0.072966      0.000058
   4     -0.115455     -0.109863     -0.000113
   5      0.092554     -0.128194      0.000047
   6     -0.284037      0.167187      0.001366
   7      0.103915      0.049062      0.000040
   8      0.005819     -0.116582     -0.000582
   9      0.156211      0.028934     -0.000159
  10      0.050167     -0.149386     -0.000828

  MM quadrupoles (xx,xy,xz,yy,yz,zz) 
  ---------------------------------- 
   1     -3.951631     -0.056179      0.000835     -4.577881      0.000043     -5.020688
   2     -0.577421     -0.053307     -0.000025     -0.604400     -0.000051     -0.559589
   3     -0.558311      0.046670      0.000009     -0.622413      0.000024     -0.558298
   4      0.693117      0.399383      0.000443      0.481174      0.000365      0.233545
   5      0.549792     -0.407216      0.000163      0.632299     -0.000245      0.241944
   6     -4.418811      0.280142      0.000213     -4.112994      0.003906     -5.020653
   7     -0.645543     -0.004118      0.000145     -0.536260      0.000122     -0.559629
   8     -0.556268      0.045015      0.000217     -0.624410     -0.000410     -0.558316
   9      0.930992      0.228057     -0.000228      0.242988     -0.000391      0.233034
  10      0.217618     -0.166400     -0.000769      0.964293      0.003861      0.241630

                     MM octupoles (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yyz,yzz,zzz) 
  ------------------------------------------------------------------------------
   1     -0.015027     -0.098740     -0.000066      0.026793     -0.000240      0.021991      0.031617     -0.000054      0.206340      0.000175
   2     -0.161918      0.085092      0.000093     -0.114618      0.000070     -0.128966      0.024961     -0.000086     -0.000239     -0.000306
   3      0.222548      0.089296     -0.000012      0.084784      0.000079      0.129785     -0.022922      0.000071     -0.021588      0.000295
   4      1.145678      0.154933      0.000793      0.448526      0.000226      0.091307      0.138410      0.000399     -0.041178      0.000192
   5     -1.067993      0.188300     -0.000646     -0.472788      0.000093     -0.102972      0.394591     -0.000289     -0.021469     -0.000261
   6      0.033022      0.076073     -0.000390      0.008775      0.000984     -0.181373     -0.121854     -0.001626      0.103322      0.002567
   7      0.081056      0.003590      0.000153     -0.009017     -0.000220      0.075738      0.284221      0.000384      0.104677      0.001100
   8     -0.204820     -0.119383     -0.000416     -0.028390      0.000150     -0.058496     -0.090704      0.000768     -0.118178     -0.001399
   9     -0.947916     -0.346572      0.000050     -0.220823     -0.000507     -0.020336     -0.782265     -0.002427     -0.100007     -0.001352
  10      0.396065      0.330437      0.001160      0.038310     -0.000697      0.079064      1.263220      0.005405      0.072376      0.000607

  MM polarizabilities in au (xx,xy,xz,yy,yz,zz) 
  --------------------------------------------- 
   1      1.593162      0.080089     -0.001211      2.525564      0.001349      3.367004
   2      0.792987      0.154230      0.000150      0.601509      0.000289      0.592515
   3      0.720949     -0.178349      0.000059      0.642026      0.000062      0.575145
   4      3.497404      2.135505      0.002229      1.845555      0.001494      1.412440
   5      2.691615     -2.246382      0.000338      2.554205     -0.001325      1.429333
   6      2.282344     -0.420740     -0.000646      1.832440     -0.006906      3.366686
   7      0.813277      0.138794      0.000410      0.581464     -0.000000      0.592795
   8      0.499325     -0.019604      0.000057      0.861273      0.002013      0.575074
   9      4.294637      1.269408      0.000409      1.056634     -0.004494      1.413326
  10      0.617428     -0.440255     -0.000294      4.622008      0.017487      1.430710


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 Threshold for neglecting two-electron integrals:  1.00D-15

 Number of two-electron integrals written:     4216362 ( 98.5% )
 Megabytes written:                             48.284

 >>>  Time used in TWOINT     is   3.05 seconds
 >>>> Total CPU  time used in HERMIT:   3.17 seconds
 >>>> Total wall time used in HERMIT:   3.20 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          T
     Using EHT  model:          F
     Number of Huckel orbitals each symmetry:   24

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
        which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -20.684561     -11.357400     -11.351264     -11.346036      -1.724840
           -1.477979      -1.046134      -1.031570      -0.744448      -0.729231
           -0.672044      -0.551742      -0.540543      -0.516178      -0.438809
           -0.319925      -0.220785      -0.162868      -0.131432      -0.129120
           -0.123435      -0.114530      -0.093199      -0.085427

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Thu Apr  7 15:36:41 2011 
     Host name              : stanley                                 

 Title lines from ".mol" input file:
     Acrolein plus 2 water molecules                                         
     ------------------------                                                

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     Restricted, closed shell Hartree-Fock calculation.

 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option.

     Wave function specification
     ============================
     For the wave function of type :      >>> HF       <<<
     Number of closed shell electrons         30
     Number of electrons in active shells      0
     Total charge of the molecule              0

     Spin multiplicity                         1
     Total number of symmetries                1
     Reference state symmetry                  1

     Orbital specifications
     ======================
     Abelian symmetry species          All |    1
                                       --- |  ---
     Occupied SCF orbitals              15 |   15
     Secondary orbitals                 61 |   61
     Total number of orbitals           76 |   76
     Number of basis functions          76 |   76

     Optimization information
     ========================
     Number of configurations                 1
     Number of orbital rotations            915
     ------------------------------------------
     Total number of variables              916

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00D-10


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy       Solvation energy    Error norm    Delta(E)
 -----------------------------------------------------------------------------
     (Precalculated two-electron integrals are transformed to P-supermatrix elements.
      Threshold for discarding integrals :  1.00D-15 )
 >>>  Time used in FORMSUP    is   0.39 seconds
   1  -189.984636021     -4.864808855971E-02    3.51976D+00   -1.90D+02
      Virial theorem: -V/T =      2.006263
      MULPOP  C       0.87; C       0.36; C       0.49; O      -0.93; H      -0.23; H      -0.23; H      -0.17; H      -0.18; 
 -----------------------------------------------------------------------------
   2  -190.480702010     -5.297194642815E-03    2.43947D+00   -4.96D-01
      Virial theorem: -V/T =      2.002406
      MULPOP  C      -0.48; C      -0.30; C      -0.16; O       0.40; H       0.11; H       0.15; H       0.14; H       0.14; 
 -----------------------------------------------------------------------------
   3  -190.789823866     -3.367834484884E-02    5.36328D-01   -3.09D-01
      Virial theorem: -V/T =      1.998888
      MULPOP  C       0.39; C      -0.15; C       0.03; O      -0.52; H       0.04; H       0.07; H       0.06; H       0.07; 
 -----------------------------------------------------------------------------
   4  -190.803881180     -3.130005668284E-02    1.01501D-01   -1.41D-02
      Virial theorem: -V/T =      2.002766
      MULPOP  C       0.32; C      -0.18; C       0.04; O      -0.44; H       0.05; H       0.08; H       0.06; H       0.07; 
 -----------------------------------------------------------------------------
   5  -190.804903590     -3.211373428889E-02    2.84205D-02   -1.02D-03
      Virial theorem: -V/T =      2.001816
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.45; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
   6  -190.805010950     -3.237128381527E-02    1.03919D-02   -1.07D-04
      Virial theorem: -V/T =      2.001853
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
   7  -190.805029214     -3.249412278894E-02    4.16588D-03   -1.83D-05
      Virial theorem: -V/T =      2.001860
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
   8  -190.805031951     -3.254332442942E-02    1.26082D-03   -2.74D-06
      Virial theorem: -V/T =      2.001864
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
   9  -190.805032215     -3.255360179244E-02    3.28197D-04   -2.64D-07
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  10  -190.805032239     -3.255662281233E-02    1.58708D-04   -2.32D-08
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  11  -190.805032242     -3.255727287602E-02    4.01227D-05   -3.23D-09
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  12  -190.805032242     -3.255754794212E-02    1.28175D-05   -4.26D-10
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  13  -190.805032242     -3.255760376667E-02    3.24394D-06   -2.77D-11
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  14  -190.805032242     -3.255761735465E-02    7.31315D-07   -2.64D-12
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  15  -190.805032242     -3.255761892699E-02    1.27688D-07   -1.71D-13
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  16  -190.805032242     -3.255761815637E-02    4.47081D-08    2.84D-13
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  17  -190.805032242     -3.255761802714E-02    1.44637D-08    6.54D-13
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  18  -190.805032242     -3.255761789452E-02    2.56115D-09   -7.39D-13
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  19  -190.805032242     -3.255761786084E-02    1.28186D-09    9.38D-13
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  20  -190.805032242     -3.255761785172E-02    3.81600D-10   -2.84D-14
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  21  -190.805032242     -3.255761785062E-02    2.00707D-10   -6.54D-13
      Virial theorem: -V/T =      2.001866
      MULPOP  C       0.33; C      -0.20; C       0.02; O      -0.46; H       0.06; H       0.10; H       0.07; H       0.08; 
 -----------------------------------------------------------------------------
  22  -190.805032242     -3.255761785066E-02    8.88223D-11   -5.68D-14

 *** DIIS converged in  22 iterations !
   - total time used in SIRFCK :              0.00 seconds


 *** SCF orbital energy analysis ***
    (incl. solvent contribution)

 Only the five lowest virtual orbital energies printed in each symmetry.

 Number of electrons :   30
 Orbital occupations :   15

 Sym       Hartree-Fock orbital energies

  1    -20.58391570   -11.36231562   -11.28783936   -11.26469396    -1.40917986
        -1.09493158    -0.89758420    -0.79510472    -0.68711392    -0.67727552
        -0.61179802    -0.56310536    -0.55845963    -0.45709971    -0.40513161
         0.07084221     0.18659463     0.20486436     0.21738597     0.25255840

    E(LUMO) :     0.07084221 au (symmetry 1)
  - E(HOMO) :    -0.40513161 au (symmetry 1)
  ------------------------------------------
    gap     :     0.47597381 au

 >>> Writing SIRIFC interface file <<<

 >>> Transforming 2-el. integrals acc. to .FINAL TRANSFORMATION =  4 <<<

 >>>> CPU and wall time for SCF :      34.400      34.627


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


     Spin multiplicity:           1
     Spatial symmetry:            1
     Total charge of molecule:    0

     QM/MM "QMMM" calculation converged :

     Charge contribution:         -0.025284994968
     Dipole contribution:          0.006484510340
     Quadrupole contribution:     -0.008011935055
     Octuple contribution:        -0.003435218159
     Electronic Pol. energy:      -0.076373706621
     Nuclear pol. energy:          0.073837929366
     Multipole Pol. energy:        0.000225797247
     Total QM/MM energy:          -0.032557617851

     Final HF energy:            -190.805032242244                 
     Nuclear repulsion:           102.664846860515
     Electronic energy:          -293.437321484909

     Final gradient norm:           0.000000000089
  -------------------------------------- 
      Output from MM property module     
  ---------------------------------------


  MM total charge:   1.00000119296339562E-010
  The MM region is charged 

  MM total charge dipole moment (x,y,z): 
 -0.63044801060850464        1.2336355263258616       3.61824895377539863E-003

  MM total permanent dipole moment (x,y,z): 
  3.00986018000000202E-002 -5.94602993999999579E-002 -1.73331100000000018E-004

  MM total induced dipole moment (x,y,z): 
 -3.09596839541735371E-002  4.56630235296047390E-002  6.02883911214753047E-005

  MM total dipole moment (x,y,z): 
 -0.63130909276267821        1.2198382504554663       3.50520624489687391E-003


 Molecular polarizability of the MM region

              Column   1    Column   2    Column   3
       1     18.13017497    0.68391670    0.00161344
       2      0.68391670   17.11580373    0.00996866
       3      0.00161344    0.00996866   14.64811272

 Isotropic polarizability 
   16.631363807498062     

 Isotropic polarizability pr. pol. site
   1.6631363807498061     

 Isotropic OPTROT (beta)
  1.70977176812031277E-005

  ---------------------------------------


 
     Date and time (Linux)  : Thu Apr  7 15:37:16 2011 
     Host name              : stanley                                 

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species  1
 ------------------------------------------

    Orbital        11       12       13       14       15       16       17
   1 C   :1s     0.0052  -0.0014   0.0000  -0.0066   0.0000   0.0000  -0.0211
   2 C   :1s    -0.0379  -0.1170   0.0001  -0.0423   0.0000  -0.0000   0.0806
   3 C   :1s    -0.0874  -0.0030   0.0000   0.0605  -0.0000  -0.0000   0.6048
   4 C   :2px   -0.0776  -0.1773   0.0004  -0.3266   0.0000  -0.0003  -0.1466
   5 C   :2py    0.3789   0.0567   0.0001  -0.0986   0.0000  -0.0001  -0.1439
   6 C   :2pz   -0.0001   0.0004   0.4786   0.0002  -0.0180  -0.4598   0.0001
   7 C   :2px    0.0234   0.0423  -0.0000   0.1087  -0.0000  -0.0001  -0.2820
   8 C   :2py   -0.1618  -0.0145  -0.0000  -0.0260  -0.0000  -0.0000  -0.3103
   9 C   :2pz    0.0000  -0.0000  -0.0347  -0.0001  -0.0175  -0.1443   0.0002
  10 C   :3d2-  -0.0267   0.0001   0.0000   0.0513  -0.0000   0.0000   0.0106
  11 C   :3d1-   0.0000   0.0000   0.0430  -0.0000  -0.0303   0.0242  -0.0000
  13 C   :3d1+   0.0000  -0.0000  -0.0077   0.0000   0.0274  -0.0239  -0.0000
  14 C   :3d2+  -0.0211   0.0128  -0.0000  -0.0363   0.0000  -0.0000   0.0059
  16 C   :1s     0.0518   0.0660  -0.0000  -0.1073  -0.0000   0.0000  -0.0324
  17 C   :1s    -0.0326   0.0474  -0.0000  -0.1834   0.0000  -0.0001   0.1781
  18 C   :2px    0.0396   0.2365  -0.0000   0.3580   0.0003  -0.0002   0.0314
  19 C   :2py   -0.4262   0.3324  -0.0002  -0.0685   0.0001  -0.0001   0.0516
  20 C   :2pz    0.0000  -0.0001   0.2080  -0.0002   0.5514  -0.2867  -0.0000
  21 C   :2px   -0.0199  -0.0335   0.0000  -0.0390   0.0000  -0.0001   0.0663
  22 C   :2py    0.1310  -0.0698   0.0000   0.0942  -0.0000  -0.0000   0.2466
  23 C   :2pz   -0.0000   0.0000  -0.0241  -0.0000   0.0217  -0.1912  -0.0001
  24 C   :3d2-  -0.0084   0.0303  -0.0000   0.0200  -0.0000  -0.0000   0.0049
  25 C   :3d1-   0.0000  -0.0000  -0.0041  -0.0000  -0.0309  -0.0251  -0.0000
  27 C   :3d1+   0.0000  -0.0000  -0.0134   0.0000   0.0008   0.0307   0.0000
  28 C   :3d2+  -0.0240  -0.0144  -0.0000  -0.0099   0.0000   0.0000  -0.0034
  29 C   :1s     0.0021  -0.0018   0.0000  -0.0012  -0.0000  -0.0000  -0.0559
  30 C   :1s     0.0677  -0.0604   0.0000   0.0175  -0.0000  -0.0000   0.1770
  31 C   :1s    -0.0108   0.0113  -0.0000   0.0926  -0.0000   0.0000   1.6762
  32 C   :2px   -0.3811  -0.3775   0.0003  -0.0945   0.0003   0.0003  -0.0035
  33 C   :2py    0.3514  -0.3035   0.0002   0.0792   0.0001   0.0001  -0.2404
  34 C   :2pz    0.0001   0.0004   0.1227   0.0001   0.5235   0.4943   0.0001
  35 C   :2px    0.0861   0.0640  -0.0000  -0.0150   0.0000   0.0002  -0.0588
  36 C   :2py   -0.0738   0.0607  -0.0000   0.0438  -0.0000   0.0001  -0.4725
  37 C   :2pz   -0.0000  -0.0001  -0.0092   0.0000   0.0243   0.2983   0.0001
  38 C   :3d2-  -0.0024   0.0319  -0.0000   0.0150   0.0000  -0.0000   0.0067
  39 C   :3d1-  -0.0000  -0.0000   0.0100  -0.0000   0.0304  -0.0124  -0.0000
  40 C   :3d0    0.0016  -0.0007   0.0000  -0.0030   0.0000  -0.0000  -0.0214
  41 C   :3d1+   0.0000  -0.0000  -0.0028   0.0000  -0.0110   0.0079  -0.0000
  42 C   :3d2+  -0.0167  -0.0101   0.0000  -0.0135  -0.0000   0.0000   0.0032
  44 O   :1s    -0.1016  -0.0105   0.0000  -0.0124  -0.0000   0.0000  -0.0323
  45 O   :1s    -0.0369  -0.0086   0.0000   0.0282   0.0000   0.0000   0.0386
  46 O   :2px    0.1529  -0.2662   0.0006   0.6893  -0.0002   0.0002   0.0230
  47 O   :2py   -0.3958  -0.1407   0.0003   0.2512  -0.0001   0.0001   0.0712
  48 O   :2pz    0.0000   0.0006   0.6698  -0.0005  -0.3298   0.3996  -0.0000
  49 O   :2px    0.0011  -0.0150   0.0000   0.0810  -0.0000   0.0001   0.0655
  50 O   :2py    0.0032  -0.0073   0.0000   0.0264  -0.0000   0.0000   0.0505
  51 O   :2pz    0.0000   0.0000   0.0188  -0.0001  -0.0363   0.1500  -0.0000
  53 O   :3d1-   0.0000  -0.0000  -0.0253   0.0000   0.0077   0.0034  -0.0000
  57 H   :1s    -0.2536   0.0090  -0.0000   0.4207  -0.0000   0.0000   0.0224
  58 H   :1s     0.0566   0.0122  -0.0000  -0.0862   0.0000  -0.0000  -0.8779
  59 H   :2px   -0.0128  -0.0003   0.0000   0.0049   0.0000  -0.0000  -0.0080
  60 H   :2py   -0.0018  -0.0021   0.0000   0.0101  -0.0000  -0.0000  -0.0003
  61 H   :2pz    0.0000   0.0000   0.0085  -0.0000   0.0009  -0.0149   0.0000
  62 H   :1s    -0.2425   0.4634  -0.0003   0.1523  -0.0000  -0.0000   0.0562
  63 H   :1s     0.0678  -0.1486   0.0001  -0.0133  -0.0000   0.0000  -0.1465
  64 H   :2px    0.0117  -0.0101   0.0000  -0.0016   0.0000  -0.0000  -0.0019
  65 H   :2py    0.0025  -0.0148   0.0000  -0.0057   0.0000  -0.0000   0.0004
  66 H   :2pz   -0.0000   0.0000   0.0052   0.0000   0.0118  -0.0080   0.0000
  67 H   :1s     0.0358   0.4845  -0.0003   0.0553  -0.0000  -0.0000  -0.0322
  68 H   :1s    -0.0144  -0.1218   0.0001   0.0147  -0.0000   0.0000  -1.3766
  70 H   :2py    0.0079   0.0131  -0.0000   0.0022   0.0000   0.0000  -0.0147
  71 H   :2pz    0.0000  -0.0000   0.0032   0.0000   0.0117   0.0149   0.0000
  72 H   :1s    -0.4195  -0.3323   0.0002  -0.1305   0.0000   0.0000  -0.0601
  73 H   :1s     0.1210   0.0900  -0.0000   0.0287  -0.0000   0.0000  -0.9614
  74 H   :2px    0.0164   0.0080  -0.0000   0.0024   0.0000   0.0000   0.0135
  75 H   :2py   -0.0004  -0.0113   0.0000  -0.0025   0.0000   0.0000  -0.0048
  76 H   :2pz   -0.0000  -0.0000   0.0024  -0.0000   0.0118   0.0167  -0.0000



 >>>> Total CPU  time used in SIRIUS :     34.46 seconds
 >>>> Total wall time used in SIRIUS :     34.69 seconds

 
     Date and time (Linux)  : Thu Apr  7 15:37:16 2011 
     Host name              : stanley                                 


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



        *****************************************************************
        ******** Output from **PROPE input processing for ABACUS ********
        *****************************************************************

 QMMM calculation


 The following molecular properties will be calculated in this run:
 ------------------------------------------------------------------

 Default print level:        0

      Electronic excitation energies 
      Natural orbital connection is used
      for perturbation dependent basis sets.


 Changes of defaults for .EXCITA:
 --------------------------------

 Number of excitation energies:    1    0    0    0    0    0    0    0
 Print level          :    0
 Integral print level :    0
 Threshold            : 1.00E-07
 Maximum iterations   :   60
 Dipole strength

 Center of mass dipole origin  :    0.973773   -1.393790    0.000425


                 .------------------------------------------------.
                 | Starting in Static Property Section (ABACUS) - |
                 `------------------------------------------------'


 
     Date and time (Linux)  : Thu Apr  7 15:37:16 2011 
     Host name              : stanley                                 

 TRACTL_1: Integral transformation abandoned,
 the required MO integrals are already available.
 >>>  Time used in EXCITA     is  12.85 seconds


   ***************************************************************************
   ************************ FINAL RESULTS from ABACUS ************************
   ***************************************************************************


 
     Date and time (Linux)  : Thu Apr  7 15:37:29 2011 
     Host name              : stanley                                 


                             Molecular geometry (au)
                             -----------------------

 C         -0.2746433464           -1.0332455534            0.0011470638
 C          2.4075280908           -1.7243202870           -0.0003155843
 C          3.0804727920           -4.1719294689           -0.0002494438
 O         -1.0584431615            1.1507997464            0.0010091138
 H         -1.6476597673           -2.6200277935            0.0023678268
 H          3.7878577518           -0.1916503544           -0.0013417055
 H          1.6611221762           -5.6704728374            0.0009146274
 H          5.0556277659           -4.7560230271           -0.0012717857





                        Molecular wave function and energy
                        ----------------------------------

     Spin multiplicity  1     State number       1     Total charge       0

     Total energy       -190.8050322422 au (Hartrees)
                         -5192.06904154 eV
                           -500958.5349 kJ/mol


                             Relativistic corrections
                             ------------------------

     Darwin correction:                          0.3812031093 au
     Mass-velocity correction:                  -0.4768191229 au

     Total relativistic correction:             -0.0956160136 au (0.0501%)
     Non-relativistic + relativistic energy:  -190.9006482558 au




                                  Dipole moment
                                  -------------

                 au               Debye          C m (/(10**-30)
              1.865132           4.740692          15.813246


                             Dipole moment components
                             ------------------------

                 au               Debye          C m (/(10**-30)

      x      0.91613569         2.32858450         7.76732181
      y     -1.62462676        -4.12938907       -13.77415929
      z     -0.00005385        -0.00013686        -0.00045652


   Units:   1 a.u. =   2.54175 Debye 
            1 a.u. =   8.47835 (10**-30) C m (SI)




                Singlet electronic excitation energies
                --------------------------------------

                 Sym.   Mode   Frequency    Frequency
                ex. st.  No.      (au)          (eV)
                =======================================
                   1        1    0.184456    5.019298
                ---------------------------------------


                Electric transition dipole moments (in a.u.)
                --------------------------------------------

  Sym.   Mode    Frequency       Velocity/Frequency              Length
 ex. st.  No.      (au)          x       y       z         x       y       z
 ==============================================================================
   1        1     0.184456     0.0001 -0.0000  0.0535    0.0001 -0.0000  0.0567
 ------------------------------------------------------------------------------


                               Oscillator strengths
                               --------------------

  Oscillator strengths are dimensionless.

  Sym.   Mode        Frequency     Oscillator-strength 
 ex. st.  No.           (eV)        velocity   length   
 -------------------------------------------------------
   1        1         5.019298        0.0004   0.0004


   Interatomic separations (in Angstrom):
   --------------------------------------

            C           C           C           O           H           H     
            ------      ------      ------      ------      ------      ------
 C     :    0.000000
 C     :    1.465700    0.000000
 C     :    2.431231    1.343281    0.000000
 O     :    1.227919    2.383018    3.568007    0.000000
 H     :    1.110397    2.197637    2.633349    2.019650    0.000000
 H     :    2.195429    1.091490    2.139278    2.661125    3.150355    0.000000
 H     :    2.659139    2.125241    1.092234    3.885970    2.381489    3.110036
 H     :    3.440501    2.130137    1.089951    4.498704    3.722962    2.506800

            H           H     
            ------      ------
 H     :    0.000000
 H     :    1.860334    0.000000


  Max interatomic separation is    4.4987 Angstrom (    8.5013 Bohr)
  between atoms    8 and    4, "H     " and "O     ".


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  C          C            1.465700
  bond distance:  C          C            1.343281
  bond distance:  O          C            1.227919
  bond distance:  H          C            1.110397
  bond distance:  H          C            1.091490
  bond distance:  H          C            1.092234
  bond distance:  H          C            1.089951


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     C          C          O            124.190
  bond angle:     C          C          H            116.421
  bond angle:     O          C          H            119.389
  bond angle:     C          C          C            119.821
  bond angle:     C          C          H            117.558
  bond angle:     C          C          H            122.621
  bond angle:     C          C          H            121.182
  bond angle:     C          C          H            121.847
  bond angle:     H          C          H            116.971




 CPU time statistics for ABACUS
 ------------------------------

 EXCITA     00:00:13      99 %

 TOTAL      00:00:13     100 %


 >>>> Total CPU  time used in ABACUS:  12.93 seconds
 >>>> Total wall time used in ABACUS:  12.96 seconds


                   .-------------------------------------------.
                   | End of Static Property Section (ABACUS) - |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:  50.58 seconds
 >>>> Total wall time used in DALTON:  50.87 seconds

 
     Date and time (Linux)  : Thu Apr  7 15:37:29 2011 
     Host name              : stanley                                 
END REFOUT

