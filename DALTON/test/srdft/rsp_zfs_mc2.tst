########## Test description ########################
START DESCRIPTION
CAS/DZP calculation of zero-field splitting (Benzene 3Sigmau)
KEYWORDS: cas short zfs
END DESCRIPTION

########## Check list ##############################
START CHECKLIST
enemc
zfs
END CHECKLIST

########## DALTON.INP ##############################
START DALINP
**DALTON
.RUN RESP
**WAVE FUN
.HF
.MP2
.MCSCF
*SCF INPUT
.DOUBLY OCCUPIED
 6 5 4 3 1 1 1 0
*CONFIGU
.SYMMETRY
 2
.ELECTRONS
 6
.SPIN MULT
 3
.INACTIVE
 6 5 4 3 0 0 0 0
.CAS SPACE
 0 0 0 0 2 2 1 1
*OPTIMI
.DETERM
*ORBITAL
.NOSUPSYM
**RESPONSE
*ESR
.ZFS
*END OF
END DALINP

########## MOLECULE.INP ############################
START MOLINP
BASIS
DZP(Dunning)
Test1 i test-suite: Geometrioptimering med symmetri, med beregning av
egenskaper i foerste og siste punkt. Generelt kontraherte basissett
    2     
        6.    6
C1     2.62337000000000    -0.00000000000000     0.00000000000000
C2     1.31168500000000     2.27190506352598     0.00000000000000
C3    -1.31168500000000     2.27190506352598     0.00000000000000
C4    -2.62337000000000     0.00000000000000     0.00000000000000
C5    -1.31168500000000    -2.27190506352598     0.00000000000000
C6     1.31168500000000    -2.27190506352598     0.00000000000000
        1.    6
H1     4.65135000000000    -0.00000000000000     0.00000000000000
H2     2.32567500000000     4.02818726189275     0.00000000000000
H3    -2.32567500000000     4.02818726189275     0.00000000000000
H4    -4.65135000000000     0.00000000000000     0.00000000000000
H5    -2.32567500000000    -4.02818726189275     0.00000000000000
H6     2.32567500000000    -4.02818726189275     0.00000000000000
END MOLINP

########## Reference Output ########################
START REFOUT


    ******************************************************************
    ***********  DALTON - An electronic structure program  ***********
    ******************************************************************

               This is output from DALTON (Release 1.1, September 2000)

                          Principal authors:

            Trygve Helgaker,     University of Oslo,        Norway 
            Hans Joergen Jensen, SDU - Odense University,   Denmark
            Poul Joergensen,     Aarhus University,         Denmark
            Jeppe Olsen,         Aarhus University,         Denmark
            Kenneth Ruud,        San Diego Superc. Center,  USA    
            Hans Agren,         KTH Stockholm,             Sweden 

                          Contributors:

            Keld L. Bak,         UNI-C,                     Denmark
            Vebjoern Bakken,     University of Oslo,        Norway 
            Ove Christiansen,    University of Lund,        Sweden 
            Paal Dahle,          University of Oslo,        Norway 
            Erik K. Dalskov,     UNI-C,                     Denmark
            Thomas Enevoldsen,   SDU - Odense University,   Denmark
            Berta Fernandez,     U.of Santiago de Compostela,Spain 
            Hanne Heiberg,       University of Oslo,        Norway 
            Hinne Hettema,       University of Auckland,    NZ     
            Dan Jonsson,         San Diego Superc. Center,  USA    
            Sheela Kirpekar,     SDU - Odense University,   Denmark
            Rika Kobayashi,      Cambridge University,      England
            Henrik Koch,         SDU - Odense University,   Denmark
            Kurt V.Mikkelsen,    University of Copenhagen,  Denmark
            Patrick Norman,      University of Linkoeping,  Sweden 
            Martin J. Packer,    University of Sheffield,   UK     
            Torgeir A. Ruden,    University of Oslo,        Norway 
            Trond Saue,          University of Toulouse,    France 
            Stephan P. A. Sauer, University of Copenhagen,  Denmark
            Peter R.Taylor,      San Diego Superc. Center,  USA    
            Olav Vahtras,        PDC, Stockholm,            Sweden 
------------------------------------------------------------------------



     NOTE:
      
     This is an experimental code for the evaluation of molecular
     properties using (MC)SCF wave functions. The authors accept no
     responsibility for the performance of the code or for the
     correctness of the results.
      
     The code (in whole or part) is provided under a license and
     is not to be reproduced for further distribution without
     the written permission of the authors or their representatives.
     See licence agreement for further information.
      
     If results obtained with this code are published, an
     appropriate citation would be:
      
     T. Helgaker, H.J.Aa. Jensen, P. Joergensen, J. Olsen, K. Ruud,H. Agren,
     K.L. Bak, V. Bakken, O. Christiansen,P. Dahle, E.K. Dalskov, T. Enevoldsen,
     B. Fernandez,H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
     H. Koch, K.V. Mikkelsen, P. Norman, M.J. Packer, T.A. Ruden, T. Saue,
     S.P.A. Sauer, P.R. Taylor, and O. Vahtras:
     "Dalton release 1.1 (2000), an electronic structure program"

 
     Date and time (Linux)  : Thu Jun  7 16:20:54 2001
     Host name              : ghiaurov.pdc.kth.se                     

 <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>


 Default print level:        0

    Integral sections will be executed
    Wave function sections will be executed
    Dynamic molecular property section will be executed
    Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************



 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

  Test1 i test-suite: Geometrioptimering med symmetri, med beregning av   
  egenskaper i foerste og siste punkt. Generelt kontraherte basissett     


                  SYMADD: Requested addition of symmetry
                  --------------------------------------

 Symmetry threshold:  0.50E-05

 Symmetry class found: D(6h)  

 The following elements were found:   X  Y  Z  


  Symmetry Operations
  -------------------

  Symmetry operations: 3



                      SYMGRP:Point group information
                      ------------------------------

Full group is:  D(6h)          
Represented as: D2h

   * The point group was generated by:

      Reflection in the yz-plane
      Reflection in the xz-plane
      Reflection in the xy-plane

   * Group multiplication table

        |  E   C2z  C2y  C2x   i   Oxy  Oxz  Oyz
   -----+----------------------------------------
     E  |  E 
    C2z | C2z   E 
    C2y | C2y  C2x   E 
    C2x | C2x  C2y  C2z   E 
     i  |  i   Oxy  Oxz  Oyz   E 
    Oxy | Oxy   i   Oyz  Oxz  C2z   E 
    Oxz | Oxz  Oyz   i   Oxy  C2y  C2x   E 
    Oyz | Oyz  Oxz  Oxy   i   C2x  C2y  C2z   E 

   * Character table

        |  E   C2z  C2y  C2x   i   Oxy  Oxz  Oyz
   -----+----------------------------------------
    Ag  |   1    1    1    1    1    1    1    1
    B3u |   1   -1   -1    1   -1    1    1   -1
    B2u |   1   -1    1   -1   -1    1   -1    1
    B1g |   1    1   -1   -1    1    1   -1   -1
    B1u |   1    1   -1   -1   -1   -1    1    1
    B2g |   1   -1    1   -1    1   -1    1   -1
    B3g |   1   -1   -1    1    1   -1   -1    1
    Au  |   1    1    1    1   -1   -1   -1   -1

   * Direct product table

        | Ag   B3u  B2u  B1g  B1u  B2g  B3g  Au 
   -----+----------------------------------------
    Ag  | Ag 
    B3u | B3u  Ag 
    B2u | B2u  B1g  Ag 
    B1g | B1g  B2u  B3u  Ag 
    B1u | B1u  B2g  B3g  Au   Ag 
    B2g | B2g  B1u  Au   B3g  B3u  Ag 
    B3g | B3g  Au   B1u  B2g  B2u  B1g  Ag 
    Au  | Au   B3g  B2g  B1u  B1g  B2u  B3u  Ag 


  Atoms and basis sets
  --------------------

  Number of atom types:     2
  Total number of atoms:   12

  Basis set used is "DZP(Dunning)" from the basis set library.

  label    atoms   charge   prim    cont     basis   
  ----------------------------------------------------------------------
  C1          4       6      29      15      [9s5p1d|4s2p1d]                              
  C6          2       6      29      15      [9s5p1d|4s2p1d]                              
  H1          4       1       7       5      [4s1p|2s1p]                                  
  H6          2       1       7       5      [4s1p|2s1p]                                  
  ----------------------------------------------------------------------
  total:     12      42     216     120
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for integrals:  1.00D-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates: 36


   1   C1   1   x      1.3116850000
   2            y      2.2719050635
   3            z      0.0000000000

   4   C1   2   x     -1.3116850000
   5            y      2.2719050635
   6            z      0.0000000000

   7   C1   3   x      1.3116850000
   8            y     -2.2719050635
   9            z      0.0000000000

  10   C1   4   x     -1.3116850000
  11            y     -2.2719050635
  12            z      0.0000000000

  13   C6   1   x      2.6233700000
  14            y      0.0000000000
  15            z      0.0000000000

  16   C6   2   x     -2.6233700000
  17            y      0.0000000000
  18            z      0.0000000000

  19   H1   1   x      2.3256750000
  20            y      4.0281872619
  21            z      0.0000000000

  22   H1   2   x     -2.3256750000
  23            y      4.0281872619
  24            z      0.0000000000

  25   H1   3   x      2.3256750000
  26            y     -4.0281872619
  27            z      0.0000000000

  28   H1   4   x     -2.3256750000
  29            y     -4.0281872619
  30            z      0.0000000000

  31   H6   1   x      4.6513500000
  32            y      0.0000000000
  33            z      0.0000000000

  34   H6   2   x     -4.6513500000
  35            y      0.0000000000
  36            z      0.0000000000



  Symmetry Coordinates
  --------------------

  Number of coordinates in each symmetry:   6  6  6  6  4  4  2  2


  Symmetry 1

   1   C1    x    [  1  -   4  +   7  -  10 ]/4
   2   C1    y    [  2  +   5  -   8  -  11 ]/4
   3   C6    x    [ 13  -  16 ]/2
   4   H1    x    [ 19  -  22  +  25  -  28 ]/4
   5   H1    y    [ 20  +  23  -  26  -  29 ]/4
   6   H6    x    [ 31  -  34 ]/2


  Symmetry 2

   7   C1    x    [  1  +   4  +   7  +  10 ]/4
   8   C1    y    [  2  -   5  -   8  +  11 ]/4
   9   C6    x    [ 13  +  16 ]/2
  10   H1    x    [ 19  +  22  +  25  +  28 ]/4
  11   H1    y    [ 20  -  23  -  26  +  29 ]/4
  12   H6    x    [ 31  +  34 ]/2


  Symmetry 3

  13   C1    x    [  1  -   4  -   7  +  10 ]/4
  14   C1    y    [  2  +   5  +   8  +  11 ]/4
  15   C6    y    [ 14  +  17 ]/2
  16   H1    x    [ 19  -  22  -  25  +  28 ]/4
  17   H1    y    [ 20  +  23  +  26  +  29 ]/4
  18   H6    y    [ 32  +  35 ]/2


  Symmetry 4

  19   C1    x    [  1  +   4  -   7  -  10 ]/4
  20   C1    y    [  2  -   5  +   8  -  11 ]/4
  21   C6    y    [ 14  -  17 ]/2
  22   H1    x    [ 19  +  22  -  25  -  28 ]/4
  23   H1    y    [ 20  -  23  +  26  -  29 ]/4
  24   H6    y    [ 32  -  35 ]/2


  Symmetry 5

  25   C1    z    [  3  +   6  +   9  +  12 ]/4
  26   C6    z    [ 15  +  18 ]/2
  27   H1    z    [ 21  +  24  +  27  +  30 ]/4
  28   H6    z    [ 33  +  36 ]/2


  Symmetry 6

  29   C1    z    [  3  -   6  +   9  -  12 ]/4
  30   C6    z    [ 15  -  18 ]/2
  31   H1    z    [ 21  -  24  +  27  -  30 ]/4
  32   H6    z    [ 33  -  36 ]/2


  Symmetry 7

  33   C1    z    [  3  +   6  -   9  -  12 ]/4
  34   H1    z    [ 21  +  24  -  27  -  30 ]/4


  Symmetry 8

  35   C1    z    [  3  -   6  -   9  +  12 ]/4
  36   H1    z    [ 21  -  24  -  27  +  30 ]/4


   Interatomic separations (in Angstroms):
   ---------------------------------------

            C1 1        C1 2        C1 3        C1 4        C6 1        C6 2

   C1 1    0.000000
   C1 2    1.388228    0.000000
   C1 3    2.404481    2.776455    0.000000
   C1 4    2.776455    2.404481    1.388228    0.000000
   C6 1    1.388228    2.404481    1.388228    2.404481    0.000000
   C6 2    2.404481    1.388228    2.404481    1.388228    2.776455    0.000000
   H1 1    1.073161    2.137438    3.376770    3.849616    2.137438    3.376770
   H1 2    2.137438    1.073161    3.849616    3.376770    3.376770    2.137438
   H1 3    3.376770    3.849616    1.073161    2.137438    2.137438    3.376770
   H1 4    3.849616    3.376770    2.137438    1.073161    3.376770    2.137438
   H6 1    2.137438    3.376770    2.137438    3.376770    1.073161    3.849616
   H6 2    3.376770    2.137438    3.376770    2.137438    3.849616    1.073161

            H1 1        H1 2        H1 3        H1 4        H6 1        H6 2

   H1 1    0.000000
   H1 2    2.461388    0.000000
   H1 3    4.263250    4.922777    0.000000
   H1 4    4.922777    4.263250    2.461388    0.000000
   H6 1    2.461388    4.263250    2.461388    4.263250    0.000000
   H6 2    4.263250    2.461388    4.263250    2.461388    4.922777    0.000000




  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2                           distance
                  ------     ------                           --------
  bond distance:    C1 2       C1 1                           1.388228
  bond distance:    C1 4       C1 3                           1.388228
  bond distance:    C6 1       C1 1                           1.388228
  bond distance:    C6 1       C1 3                           1.388228
  bond distance:    C6 2       C1 2                           1.388228
  bond distance:    C6 2       C1 4                           1.388228
  bond distance:    H1 1       C1 1                           1.073161
  bond distance:    H1 2       C1 2                           1.073161
  bond distance:    H1 3       C1 3                           1.073161
  bond distance:    H1 4       C1 4                           1.073161
  bond distance:    H6 1       C6 1                           1.073161
  bond distance:    H6 2       C6 2                           1.073161


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3                   angle
                  ------     ------     ------                   -----
  bond angle:       C1 2       C1 1       C6 1                 120.000
  bond angle:       C1 2       C1 1       H1 1                 120.000
  bond angle:       C6 1       C1 1       H1 1                 120.000
  bond angle:       C1 1       C1 2       C6 2                 120.000
  bond angle:       C1 1       C1 2       H1 2                 120.000
  bond angle:       C6 2       C1 2       H1 2                 120.000
  bond angle:       C1 4       C1 3       C6 1                 120.000
  bond angle:       C1 4       C1 3       H1 3                 120.000
  bond angle:       C6 1       C1 3       H1 3                 120.000
  bond angle:       C1 3       C1 4       C6 2                 120.000
  bond angle:       C1 3       C1 4       H1 4                 120.000
  bond angle:       C6 2       C1 4       H1 4                 120.000
  bond angle:       C1 1       C6 1       C1 3                 120.000
  bond angle:       C1 1       C6 1       H6 1                 120.000
  bond angle:       C1 3       C6 1       H6 1                 120.000
  bond angle:       C1 2       C6 2       C1 4                 120.000
  bond angle:       C1 2       C6 2       H6 2                 120.000
  bond angle:       C1 4       C6 2       H6 2                 120.000


  Nuclear repulsion energy :  204.624364126512


  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:         26  26  19  19   9   9   6   6


  Symmetry  Ag ( 1)

    1     C1       1s         1  +   2  +   3  +   4
    2     C1       1s         5  +   6  +   7  +   8
    3     C1       1s         9  +  10  +  11  +  12
    4     C1       1s        13  +  14  +  15  +  16
    5     C1       2px       17  -  18  +  19  -  20
    6     C1       2py       21  +  22  -  23  -  24
    7     C1       2px       29  -  30  +  31  -  32
    8     C1       2py       33  +  34  -  35  -  36
    9     C1       3d2-      41  -  42  -  43  +  44
   10     C1       3d0       49  +  50  +  51  +  52
   11     C1       3d2+      57  +  58  +  59  +  60
   12     C6       1s        61  +  62
   13     C6       1s        63  +  64
   14     C6       1s        65  +  66
   15     C6       1s        67  +  68
   16     C6       2px       69  -  70
   17     C6       2px       75  -  76
   18     C6       3d0       85  +  86
   19     C6       3d2+      89  +  90
   20     H1       1s        91  +  92  +  93  +  94
   21     H1       1s        95  +  96  +  97  +  98
   22     H1       2px       99  - 100  + 101  - 102
   23     H1       2py      103  + 104  - 105  - 106
   24     H6       1s       111  + 112
   25     H6       1s       113  + 114
   26     H6       2px      115  - 116


  Symmetry  B3u( 2)

   27     C1       1s         1  -   2  +   3  -   4
   28     C1       1s         5  -   6  +   7  -   8
   29     C1       1s         9  -  10  +  11  -  12
   30     C1       1s        13  -  14  +  15  -  16
   31     C1       2px       17  +  18  +  19  +  20
   32     C1       2py       21  -  22  -  23  +  24
   33     C1       2px       29  +  30  +  31  +  32
   34     C1       2py       33  -  34  -  35  +  36
   35     C1       3d2-      41  +  42  -  43  -  44
   36     C1       3d0       49  -  50  +  51  -  52
   37     C1       3d2+      57  -  58  +  59  -  60
   38     C6       1s        61  -  62
   39     C6       1s        63  -  64
   40     C6       1s        65  -  66
   41     C6       1s        67  -  68
   42     C6       2px       69  +  70
   43     C6       2px       75  +  76
   44     C6       3d0       85  -  86
   45     C6       3d2+      89  -  90
   46     H1       1s        91  -  92  +  93  -  94
   47     H1       1s        95  -  96  +  97  -  98
   48     H1       2px       99  + 100  + 101  + 102
   49     H1       2py      103  - 104  - 105  + 106
   50     H6       1s       111  - 112
   51     H6       1s       113  - 114
   52     H6       2px      115  + 116


  Symmetry  B2u( 3)

   53     C1       1s         1  +   2  -   3  -   4
   54     C1       1s         5  +   6  -   7  -   8
   55     C1       1s         9  +  10  -  11  -  12
   56     C1       1s        13  +  14  -  15  -  16
   57     C1       2px       17  -  18  -  19  +  20
   58     C1       2py       21  +  22  +  23  +  24
   59     C1       2px       29  -  30  -  31  +  32
   60     C1       2py       33  +  34  +  35  +  36
   61     C1       3d2-      41  -  42  +  43  -  44
   62     C1       3d0       49  +  50  -  51  -  52
   63     C1       3d2+      57  +  58  -  59  -  60
   64     C6       2py       71  +  72
   65     C6       2py       77  +  78
   66     C6       3d2-      81  -  82
   67     H1       1s        91  +  92  -  93  -  94
   68     H1       1s        95  +  96  -  97  -  98
   69     H1       2px       99  - 100  - 101  + 102
   70     H1       2py      103  + 104  + 105  + 106
   71     H6       2py      117  + 118


  Symmetry  B1g( 4)

   72     C1       1s         1  -   2  -   3  +   4
   73     C1       1s         5  -   6  -   7  +   8
   74     C1       1s         9  -  10  -  11  +  12
   75     C1       1s        13  -  14  -  15  +  16
   76     C1       2px       17  +  18  -  19  -  20
   77     C1       2py       21  -  22  +  23  -  24
   78     C1       2px       29  +  30  -  31  -  32
   79     C1       2py       33  -  34  +  35  -  36
   80     C1       3d2-      41  +  42  +  43  +  44
   81     C1       3d0       49  -  50  -  51  +  52
   82     C1       3d2+      57  -  58  -  59  +  60
   83     C6       2py       71  -  72
   84     C6       2py       77  -  78
   85     C6       3d2-      81  +  82
   86     H1       1s        91  -  92  -  93  +  94
   87     H1       1s        95  -  96  -  97  +  98
   88     H1       2px       99  + 100  - 101  - 102
   89     H1       2py      103  - 104  + 105  - 106
   90     H6       2py      117  - 118


  Symmetry  B1u( 5)

   91     C1       2pz       25  +  26  +  27  +  28
   92     C1       2pz       37  +  38  +  39  +  40
   93     C1       3d1-      45  +  46  -  47  -  48
   94     C1       3d1+      53  -  54  +  55  -  56
   95     C6       2pz       73  +  74
   96     C6       2pz       79  +  80
   97     C6       3d1+      87  -  88
   98     H1       2pz      107  + 108  + 109  + 110
   99     H6       2pz      119  + 120


  Symmetry  B2g( 6)

  100     C1       2pz       25  -  26  +  27  -  28
  101     C1       2pz       37  -  38  +  39  -  40
  102     C1       3d1-      45  -  46  -  47  +  48
  103     C1       3d1+      53  +  54  +  55  +  56
  104     C6       2pz       73  -  74
  105     C6       2pz       79  -  80
  106     C6       3d1+      87  +  88
  107     H1       2pz      107  - 108  + 109  - 110
  108     H6       2pz      119  - 120


  Symmetry  B3g( 7)

  109     C1       2pz       25  +  26  -  27  -  28
  110     C1       2pz       37  +  38  -  39  -  40
  111     C1       3d1-      45  +  46  +  47  +  48
  112     C1       3d1+      53  -  54  -  55  +  56
  113     C6       3d1-      83  +  84
  114     H1       2pz      107  + 108  - 109  - 110


  Symmetry  Au ( 8)

  115     C1       2pz       25  -  26  -  27  +  28
  116     C1       2pz       37  -  38  -  39  +  40
  117     C1       3d1-      45  -  46  +  47  -  48
  118     C1       3d1+      53  +  54  -  55  -  56
  119     C6       3d1-      83  -  84
  120     H1       2pz      107  - 108  - 109  + 110

  Symmetries of electric field:  B3u(2)  B2u(3)  B1u(5)

  Symmetries of magnetic field:  B3g(7)  B2g(6)  B1g(4)


 ************************************************************************
 ************************** Output from HERINT **************************
 ************************************************************************


 >>> Time used in ONEDRV is   0.17 seconds


 >>> Time used in OVERLA is   0.17 seconds


 Number of two-electron integrals written:     3519734 ( 13.4% )
 Megabytes written:                             40.308



 >>> Time used in TWOINT is  1 minute   1 second 



 ***********************************************************************
 ************************* Output from FORMSUP *************************
 ***********************************************************************

  NOSYM :  F
 OLDSUP :  F

 Threshold for discarding integrals :  1.00D-15

 >>> Time used in FRMSUP is   2.49 seconds

 >>>> Total CPU  time used in HERMIT:  1 minute  12 seconds
 >>>> Total wall time used in HERMIT:  1 minute  36 seconds

- End of Integral Section
     Starting Wave Function Section -

 ************************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program  *Apr 96*
 ************************************************************************

 
     Date and time (Linux)  : Thu Jun  7 16:22:31 2001
     Host name              : ghiaurov.pdc.kth.se                     

 Title lines from integral program:
     Test1 i test-suite: Geometrioptimering med symmetri, med beregning av   
     egenskaper i foerste og siste punkt. Generelt kontraherte basissett     

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

     MC-SCF optimization.


     Multi-configurational response calculation.

               Type: complete active space calculation (CAS).

     This is a combination run starting with

               a RHF calculation
               an MP2 calculation 


 Initial molecular orbitals are obtained according to
 ".MOSTART HUCKEL" input option.

     Wave function specification
     ============================
     Number of closed shell electrons         36
     Number of electrons in active shells      6
     Number of active orbitals                 6
     Total number of orbitals                120

     Spin multiplicity                         3
     Total number of symmetries                8
     Reference state symmetry                  2

     Orbital specifications
     ======================
     Abelian symmetry species           1   2   3   4   5   6   7   8
                                       --  --  --  --  --  --  --  --
     Inactive orbitals                  6   5   4   3   0   0   0   0
     Active orbitals                    0   0   0   0   2   2   1   1
     Secondary orbitals                20  21  15  16   7   7   5   5
     Total number of orbitals          26  26  19  19   9   9   6   6
     Number of basis functions         26  26  19  19   9   9   6   6

     Occupied HF orbitals               6   5   4   3   1   1   1   0

     Optimization information
     ========================
     Number of determinants               61
     Number of orbital rotations         371
     ---------------------------------------
     Total number of variables           432

     Maximum number of macro iterations      15
     Maximum number of micro iterations     360
     Threshold for gradient            1.00D-05
     Number of initial trial vectors          1
     Number of initial CI iterations          3
     Number of simultaneous trial vectors     1

     This calculation converges to the lowest state
     for the specified symmetry and spin species.

     Maximum number of NEO/NR iterations  24


 >>>>> DIIS optimization of Hartree-Fock <<<<<

 C1-DIIS algorithm; max error vectors =   10

 Iter      Total energy        Error norm    Delta(E)
 ----------------------------------------------------
    1     -229.642523158134    6.45828D+00   -2.30D+02
    2     -230.667740242564    1.02041D+00   -1.03D+00
    3     -230.729614726100    4.54243D-01   -6.19D-02
    4     -230.742700319767    4.45098D-02   -1.31D-02
    5     -230.742971712235    6.69432D-03   -2.71D-04
    6     -230.742978564004    1.04234D-03   -6.85D-06
    7     -230.742978728368    1.23290D-04   -1.64D-07
    8     -230.742978730670    1.81687D-05   -2.30D-09
    9     -230.742978730727    2.08412D-06   -5.75D-11
   10     -230.742978730728    2.58135D-07   -3.98D-13
 DIIS converged in  10 iterations !


 >>>>> Output from SIRIUS MP2 module <<<<<

 Reference: H.J.Aa.Jensen, P. Jorgensen, H. Agren, and J. Olsen,
            J. Chem. Phys. 88, 3834 (1988); 89, 5354 (1988)



 Check that orbitals are canonical HARTREE-FOCK orbitals

 Number of electrons :   42
 Orbital occupations :    6    5    4    3    1    1    1    0


 Hartree-Fock electronic energy:        -435.367342857240

 Hartree-Fock total      energy:        -230.742978730728

 Hartree-Fock orbital energies, symmetry 1

       -11.23958388   -11.23776004    -1.15319745    -0.82567709    -0.71295443
        -0.49488732     0.27821347     0.30214828     0.48841679     0.53656354
         0.58301245     0.67112306     0.84462407     1.13500616     1.20634814
         1.47490678     1.94188893     2.03030330     2.13537087     2.61744348
         2.66587427     3.05640619     3.45356451     3.75452724    23.84171550
        24.33647541

 Hartree-Fock orbital energies, symmetry 2

       -11.23900570   -11.23715237    -1.01705185    -0.64528901    -0.58942202
         0.32386753     0.33194337     0.40626243     0.57530615     0.65626603
         0.85493793     0.99447916     1.17963703     1.19521837     1.31900884
         1.76138423     2.09975251     2.29357382     2.39851424     2.75281121
         3.29917753     3.33043696     3.47841033     4.30976009    24.17047257
        24.46517687

 Hartree-Fock orbital energies, symmetry 3

       -11.23900570    -1.01705185    -0.61885520    -0.58942202     0.32386753
         0.40626243     0.41760498     0.65626603     0.85493793     1.17963703
         1.31900884     1.81439508     2.09975251     2.29357382     2.39851424
         2.57848887     3.33043696     3.47841033    24.17047257

 Hartree-Fock orbital energies, symmetry 4

       -11.23776004    -0.82567709    -0.49488732     0.27821347     0.48841679
         0.53656354     0.65186125     0.84462407     1.05560072     1.13500616
         1.47490678     2.13537087     2.42945564     2.61744348     2.66587427
         3.05640619     3.65616446     3.75452724    24.33647541

 Hartree-Fock orbital energies, symmetry 5

        -0.50132573     0.12889581     0.47806230     0.62541958     1.44314981
         1.61396000     2.11673651     2.48881752     2.78109635

 Hartree-Fock orbital energies, symmetry 6

        -0.33557206     0.33100525     0.52857507     0.71835379     1.57452506
         2.11550063     2.51520423     2.55071606     3.00618409

 Hartree-Fock orbital energies, symmetry 7

        -0.33557206     0.52857507     1.43274787     1.57452506     2.51520423
         2.55071606

 Hartree-Fock orbital energies, symmetry 8

         0.12889581     0.62541958     1.61396000     2.11673651     2.78109635
         2.91981348

    E(LUMO) :     0.12889581
  - E(HOMO) :    -0.33557206
  --------------------------
    gap     :     0.46446787

 ABS SUM OF OFF-DIAGONAL FOCK CORE ELEMENTS ARE :  1.22D-06

 MP2 move   0.575443 electrons to unoccupied HF orbitals


 Hartree-Fock total energy :          -230.7429787307
 + MP2 contribution        :            -0.8537540638
 = MP2 second order energy :          -231.5967327945

 Natural orbital occupation numbers, symmetry 1

         1.99947446     1.99946452     1.98550570     1.97477195     1.97018289
         1.95746631     0.02294036     0.01759627     0.01574271     0.01162312
         0.00839856     0.00743648     0.00285439     0.00201397     0.00190234
         0.00162118     0.00112228     0.00094197     0.00057840     0.00055689
         0.00032459     0.00030692     0.00023031     0.00014833     0.00013495
         0.00006073

 Sum    11.98340059
 RHF    12.00000000
 Diff.  -0.01659941

 Natural orbital occupation numbers, symmetry 2

         1.99947136     1.99946120     1.98144905     1.97008862     1.96439667
         0.02396292     0.01801011     0.01607124     0.01300085     0.00825307
         0.00496167     0.00260422     0.00218115     0.00153081     0.00110194
         0.00090486     0.00058720     0.00043249     0.00036744     0.00031616
         0.00031102     0.00024100     0.00013102     0.00010492     0.00008946
         0.00004558

 Sum    10.01007602
 RHF    10.00000000
 Diff.   0.01007602

 Natural orbital occupation numbers, symmetry 3

         1.99947136     1.98144905     1.96439667     1.95916253     0.02396292
         0.01607124     0.01101581     0.00825307     0.00496167     0.00218115
         0.00153081     0.00150182     0.00110194     0.00056769     0.00043249
         0.00036744     0.00031102     0.00013102     0.00010492

 Sum     7.97697461
 RHF     8.00000000
 Diff.  -0.02302539

 Natural orbital occupation numbers, symmetry 4

         1.99946452     1.97477195     1.95746631     0.02294036     0.01574271
         0.01521031     0.00743648     0.00285439     0.00196556     0.00190234
         0.00112228     0.00094197     0.00055689     0.00032459     0.00031545
         0.00023031     0.00014833     0.00006105     0.00006073

 Sum     6.00351654
 RHF     6.00000000
 Diff.   0.00351654

 Natural orbital occupation numbers, symmetry 5

         1.95252249     0.06578104     0.00960531     0.00554273     0.00478633
         0.00310540     0.00088294     0.00062909     0.00046898

 Sum     2.04332431
 RHF     2.00000000
 Diff.   0.04332431

 Natural orbital occupation numbers, symmetry 6

         1.91900630     0.02870783     0.00937900     0.00344204     0.00260850
         0.00124400     0.00059893     0.00044666     0.00025967

 Sum     1.96569291
 RHF     2.00000000
 Diff.  -0.03430709

 Natural orbital occupation numbers, symmetry 7

         1.91900630     0.00937900     0.00758628     0.00344204     0.00124400
         0.00059893

 Sum     1.94125654
 RHF     2.00000000
 Diff.  -0.05874346

 Natural orbital occupation numbers, symmetry 8

         0.06578104     0.00478633     0.00310540     0.00088294     0.00073380
         0.00046898

 Sum     0.07575849
 RHF     0.00000000
 Diff.   0.07575849

 Time used for MP2 natural orbitals :      -2.179 CPU seconds.


        SIRIUS MCSCF optimization (SIROPT)
 ================================================



 <<<<< Output from SIRIUS CI module (CICTL) >>>>>




 (CIST1)  4 lowest diagonal elements:

 Element no. Config.no.    Active energy      Total energy

         1 :          2     -6.1495461916   -230.5231153754
         2 :         46     -6.1495461916   -230.5231153754
         3 :          7     -5.7394614275   -230.1130306112
         4 :         30     -5.7157869408   -230.0893561245


 (CIST1) number of start vectors is increased from  1 to  2 because of degeneracy


 Convergence threshold for CI optimization :     0.00000500



 *** Reached maximum number of CI iterations:    3
              1 CI roots are not converged.


 CI energies and residuals:
    1     -230.629340727407200       3.20D-02

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  5

   1.916272966   0.519708245

 Symmetry  6

   1.479994995   0.084226916

 Symmetry  7

   1.470430470

 Symmetry  8

   0.529366408


 <<< MACRO ITERATION  1 >>>
 --------------------------

 Total MCSCF energy       :     -230.629340727407400       (MACRO    1)

 Norm of total gradient   :        0.366166118259
 -    of CI gradient      :        0.063918160446
 -    of orbital gradient :        0.360544165015

 Residual norm when dim(red L) =   5
 NEO root     CSF        orbital          total
    1     0.04483281     0.01024069     0.04598753 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  5

   1.909162345   0.533808192

 Symmetry  6

   1.467105450   0.088992399

 Symmetry  7

   1.457435462

 Symmetry  8

   0.543496152


 <<< MACRO ITERATION  2 >>>
 --------------------------

 Total MCSCF energy       :     -230.669478339613200       (MACRO    2)

 Norm of total gradient   :        0.042080485378
 -    of CI gradient      :        0.028911120086
 -    of orbital gradient :        0.030576369716

 Residual norm when dim(red L) =   5
 NEO root     CSF        orbital          total
    1     0.00286314     0.00160581     0.00328271 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  5

   1.911256445   0.531791256

 Symmetry  6

   1.468565427   0.088028512

 Symmetry  7

   1.460294311

 Symmetry  8

   0.540064048


 <<< MACRO ITERATION  3 >>>
 --------------------------

 Total MCSCF energy       :     -230.669912157510900       (MACRO    3)

 Norm of total gradient   :        0.003312033057
 -    of CI gradient      :        0.002913515780
 -    of orbital gradient :        0.001575115478

 Residual norm when dim(red L) =  16
 NEO root     CSF        orbital          total
    1     0.00002672     0.00000443     0.00002708 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  5

   1.910746318   0.535566584

 Symmetry  6

   1.464746758   0.088628059

 Symmetry  7

   1.464708647

 Symmetry  8

   0.535603635


 <<< MACRO ITERATION  4 >>>
 --------------------------

 Total MCSCF energy       :     -230.669917814369500       (MACRO    4)

 Norm of total gradient   :        0.000028912898
 -    of CI gradient      :        0.000026729899
 -    of orbital gradient :        0.000011021261

 Residual norm when dim(red L) =   4
 NEO root     CSF        orbital          total
    1     0.00000504     0.00000608     0.00000789 converged

 (NEONEX) NEO vector is converged.

   <<< OUTPUT FROM SIRCNO >>>    Keyword = FD+NO 


 Occupations of CAS natural orbitals:

 Symmetry  5

   1.910745849   0.535575640

 Symmetry  6

   1.464736007   0.088630838

 Symmetry  7

   1.464714921

 Symmetry  8

   0.535596744


 <<< MACRO ITERATION  5 >>>
 --------------------------

 Total MCSCF energy       :     -230.669917814688100       (MACRO    5)

 Norm of total gradient   :        0.000007893530
 -    of CI gradient      :        0.000005037628
 -    of orbital gradient :        0.000006077016

 *** Optimization control: MCSCF converged ***
     Number of macro iterations used            5
     Number of micro iterations used           23
     Total number of CPU seconds used         0.00


               >>> SIRIUS OPTIMIZATION STATISTICS <<<

 
     Date and time (Linux)  : Thu Jun  7 16:32:02 2001
     Host name              : ghiaurov.pdc.kth.se                     


  ITER ITMIC     EMCSCF           GRDNRM        RATIO      STPLNG
 ---------------------------------------------------------------------
    1    4   -230.629340727407   0.3661661183  0.000000   0.3603662769
    2    3   -230.669478339613   0.0420804854  0.923964   0.0282942486
    3   14   -230.669912157511   0.0033120331  1.002426   0.0061726921
    4    2   -230.669917814370   0.0000289129  1.000379   0.0000286766
    5    0   -230.669917814688   0.0000078935  0.997507   0.0000000000


  ITER  INDGCM  GCIMAX      GCINRM     INDGOM  GOBMAX      GOBNRM      GRDNRM
 ------------------------------------------------------------------------------
    1     45    0.026695    0.063918    368    0.127809    0.360544    0.366166
    2     45   -0.012718    0.028911    367    0.012216    0.030576    0.042080
    3     30    0.001590    0.002914    348    0.000604    0.001575    0.003312
    4     29    0.000009    0.000027    367    0.000005    0.000011    0.000029
    5     35   -0.000002    0.000005    350   -0.000001    0.000006    0.000008


  ITER ITMIC NCLIN NOLIN   TIMMAC    TIMITR    TIMMIC    TIMLIN    TIMMIC/ITMIC
 ------------------------------------------------------------------------------

    1     4     1     3      0.00      0.00      0.00      0.00      0.00
    2     3     2     2      0.00      0.00      0.00      0.00      0.00
    3    14     9     6      0.00      0.00      0.00      0.00      0.00
    4     2     2     1      0.00      0.00      0.00      0.00      0.00
    5     0     0     0      0.00      0.00      0.00      0.00


 ITER         EMY                 EACTIV              EMCSCF

    1   -428.997933310287     -6.255771543633   -230.629340727407
    2   -428.969840396548     -6.324002069577   -230.669478339613
    3   -428.977448756307     -6.316827527716   -230.669912157511
    4   -428.976923426650     -6.317358514232   -230.669917814370
    5   -428.976925555991     -6.317356385209   -230.669917814688


 ITER         DEPRED              DEACT               RATIO

    1      0.000000000000      0.000000000000      0.000000000000
    2     -0.043440686041     -0.040137612206      0.923963589527
    3     -0.000432768159     -0.000433817898      1.002425637570
    4     -0.000005654716     -0.000005656859      1.000378933562
    5     -0.000000000319     -0.000000000319      0.997506691850


 ITER    BETA           GAMMA             STPLNG              RTRUST

    1      0.20000000  1.00000000      0.360366276859      0.700000000000
    2      0.20000000  1.00000000      0.028294248642      0.700000000000
    3      0.20000000  1.00000000      0.006172692080      0.700000000000
    4      0.20000000  1.00000000      0.000028676606      0.700000000000
    5      0.00000000  0.00000000      0.000000000000      0.700000000000


 Reduced L root no.  1
 ITER         EVAL              EVEC(1)           EVEC(2)           EVEC(3)
 ----------------------------------------------------------------------------
    1   -0.003457295773    0.997412798098   -0.070753603176   -0.004646049609
    2   -0.000034620344    0.999983989094   -0.002122472383   -0.005133998632
    3   -0.000000452377    0.999999237958   -0.001000999709   -0.000226407965
    4   -0.000000000026    0.999999999984   -0.000004662606   -0.000003028779
    5    0.000000000000    0.000000000000    0.000000000000    0.000000000000


                    >>> FINAL RESULTS FROM SIRIUS <<<

     Spin multiplicity:           3
     Spatial symmetry:            2
     State number:                1

     Final MCSCF energy:         -230.669917814688
     Nuclear repulsion:           204.624364126512
     Electronic energy:          -435.294281941200

     Final gradient norm:           0.000007893530

 
     Date and time (Linux)  : Thu Jun  7 16:32:02 2001
     Host name              : ghiaurov.pdc.kth.se                     


 Occupancies of natural orbitals
 -------------------------------

 Symmetry  1

   2.000000000   2.000000000   2.000000000   2.000000000   2.000000000
   2.000000000

 Sum =          12.000000000

 Symmetry  2

   2.000000000   2.000000000   2.000000000   2.000000000   2.000000000

 Sum =          10.000000000

 Symmetry  3

   2.000000000   2.000000000   2.000000000   2.000000000

 Sum =           8.000000000

 Symmetry  4

   2.000000000   2.000000000   2.000000000

 Sum =           6.000000000

 Symmetry  5

   1.910745849   0.535575640

 Sum =           2.446321490

 Symmetry  6

   1.464736007   0.088630838

 Sum =           1.553366845

 Symmetry  7

   1.464714921

 Sum =           1.464714921

 Symmetry  8

   0.535596744

 Sum =           0.535596744

     Molecular orbitals for symmetry species   1

 Orbital          1        2        3        4        5        6
   1  C1  1s     0.2452  -0.1735  -0.0569  -0.0309  -0.0054   0.0014
   2  C1  1s     0.1785  -0.1263  -0.0756  -0.0414  -0.0070   0.0019
   3  C1  1s     0.0022  -0.0013   0.2000   0.1097   0.0210  -0.0045
   4  C1  1s    -0.0001   0.0004   0.1042   0.0858   0.0248  -0.0027
   5  C1  2px   -0.0001   0.0003  -0.0362  -0.1586   0.1113   0.2899
   6  C1  2py   -0.0002  -0.0002  -0.0626   0.1350   0.1927  -0.0162
   7  C1  2px    0.0000  -0.0003   0.0038   0.0146   0.0118   0.0497
   8  C1  2py    0.0000  -0.0003   0.0066   0.0175   0.0204   0.0059
   9  C1  3d2-  -0.0001   0.0000  -0.0034   0.0073  -0.0004   0.0023
  10  C1  3d0   -0.0001   0.0001  -0.0103  -0.0024   0.0006  -0.0006
  11  C1  3d2+   0.0000  -0.0001   0.0019   0.0117   0.0002  -0.0202
  12  C6  1s     0.2453   0.3469  -0.0569   0.0618  -0.0054  -0.0028
  13  C6  1s     0.1785   0.2526  -0.0756   0.0828  -0.0070  -0.0038
  14  C6  1s     0.0022   0.0027   0.2000  -0.2194   0.0210   0.0089
  15  C6  1s    -0.0001  -0.0008   0.1042  -0.1715   0.0248   0.0053
  16  C6  2px   -0.0003   0.0000  -0.0723  -0.0752   0.2225  -0.2618
  17  C6  2px    0.0000   0.0008   0.0076  -0.0449   0.0235  -0.0600
  18  C6  3d0   -0.0001  -0.0001  -0.0103   0.0049   0.0006   0.0012
  19  C6  3d2+  -0.0001   0.0000  -0.0039  -0.0009  -0.0005  -0.0241
  20  H1  1s     0.0002  -0.0002   0.0386   0.0581   0.1073   0.0788
  21  H1  1s     0.0000   0.0002  -0.0024   0.0211   0.0776   0.0787
  22  H1  2px    0.0000   0.0000  -0.0025  -0.0052  -0.0042   0.0015
  23  H1  2py   -0.0001   0.0000  -0.0043  -0.0041  -0.0072  -0.0048
  24  H6  1s     0.0002   0.0003   0.0386  -0.1161   0.1073  -0.1576
  25  H6  1s     0.0000  -0.0004  -0.0024  -0.0421   0.0776  -0.1573
  26  H6  2px   -0.0001  -0.0001  -0.0050   0.0124  -0.0084   0.0068

     Molecular orbitals for symmetry species   2

 Orbital          1        2        3        4        5
   1  C1  1s     0.1734   0.2454  -0.0395  -0.0348   0.0062
   2  C1  1s     0.1262   0.1787  -0.0525  -0.0469   0.0083
   3  C1  1s     0.0015   0.0018   0.1358   0.1250  -0.0211
   4  C1  1s    -0.0007  -0.0010   0.1128   0.1427  -0.0401
   5  C1  2px    0.0002   0.0001   0.0928   0.0774   0.0098
   6  C1  2py   -0.0003   0.0001  -0.0740   0.1341   0.1912
   7  C1  2px    0.0006   0.0003  -0.0332   0.0106   0.0339
   8  C1  2py   -0.0005   0.0005   0.0286   0.0183   0.0034
   9  C1  3d2-  -0.0001   0.0000  -0.0088   0.0123   0.0125
  10  C1  3d0   -0.0001  -0.0001  -0.0056  -0.0022  -0.0002
  11  C1  3d2+  -0.0001   0.0000  -0.0095  -0.0071   0.0098
  12  C6  1s     0.3469  -0.2453  -0.0790   0.0348   0.0123
  13  C6  1s     0.2525  -0.1786  -0.1050   0.0469   0.0166
  14  C6  1s     0.0030  -0.0018   0.2716  -0.1250  -0.0422
  15  C6  1s    -0.0014   0.0010   0.2256  -0.1427  -0.0802
  16  C6  2px   -0.0002  -0.0001  -0.0354  -0.1548   0.3409
  17  C6  2px   -0.0002  -0.0005   0.0163  -0.0211   0.0398
  18  C6  3d0   -0.0001   0.0001  -0.0112   0.0022  -0.0004
  19  C6  3d2+  -0.0001   0.0000  -0.0058  -0.0142   0.0119
  20  H1  1s     0.0002   0.0002   0.0398   0.1356   0.0792
  21  H1  1s     0.0001  -0.0002   0.0087   0.1000   0.0736
  22  H1  2px    0.0000   0.0000  -0.0017  -0.0050  -0.0029
  23  H1  2py    0.0000   0.0000  -0.0048  -0.0086  -0.0035
  24  H6  1s     0.0003  -0.0002   0.0796  -0.1356   0.1584
  25  H6  1s     0.0002   0.0002   0.0174  -0.1000   0.1473
  26  H6  2px   -0.0001   0.0000  -0.0101   0.0099  -0.0090

     Molecular orbitals for symmetry species   3

 Orbital          1        2        3        4
   1  C1  1s     0.3004  -0.0684   0.0000   0.0107
   2  C1  1s     0.2186  -0.0909   0.0000   0.0144
   3  C1  1s     0.0026   0.2352   0.0000  -0.0365
   4  C1  1s    -0.0012   0.1954   0.0000  -0.0694
   5  C1  2px   -0.0003  -0.0740   0.2690   0.1912
   6  C1  2py   -0.0001   0.0073  -0.1553   0.2305
   7  C1  2px   -0.0005   0.0286   0.0305   0.0034
   8  C1  2py    0.0000  -0.0002  -0.0176   0.0379
   9  C1  3d2-   0.0000  -0.0007  -0.0065   0.0047
  10  C1  3d0   -0.0001  -0.0097   0.0000  -0.0004
  11  C1  3d2+   0.0001   0.0088  -0.0113  -0.0125
  12  C6  2py    0.0004   0.1356   0.3106  -0.1006
  13  C6  2py    0.0009  -0.0497   0.0353   0.0320
  14  C6  3d2-  -0.0002  -0.0146  -0.0131   0.0170
  15  H1  1s     0.0003   0.0689   0.0000   0.1372
  16  H1  1s     0.0002   0.0150   0.0000   0.1276
  17  H1  2px    0.0000  -0.0048   0.0045  -0.0035
  18  H1  2py   -0.0001  -0.0073  -0.0026  -0.0070
  19  H6  2py    0.0000   0.0011   0.0052  -0.0009

     Molecular orbitals for symmetry species   4

 Orbital          1        2        3
   1  C1  1s     0.3005  -0.0536   0.0025
   2  C1  1s     0.2188  -0.0717   0.0033
   3  C1  1s     0.0023   0.1900  -0.0077
   4  C1  1s    -0.0007   0.1485  -0.0046
   5  C1  2px    0.0002   0.1350  -0.0162
   6  C1  2py   -0.0001  -0.0028   0.2712
   7  C1  2px    0.0003   0.0175   0.0059
   8  C1  2py    0.0006   0.0348   0.0566
   9  C1  3d2-   0.0000  -0.0033   0.0228
  10  C1  3d0   -0.0001  -0.0042  -0.0010
  11  C1  3d2+   0.0000  -0.0073  -0.0023
  12  C6  2py    0.0004   0.2366  -0.2993
  13  C6  2py   -0.0001  -0.0045  -0.0463
  14  C6  3d2-  -0.0001  -0.0159   0.0189
  15  H1  1s     0.0003   0.1006   0.1365
  16  H1  1s    -0.0003   0.0365   0.1363
  17  H1  2px    0.0000  -0.0041  -0.0048
  18  H1  2py    0.0000  -0.0100  -0.0040
  19  H6  2py    0.0000   0.0028  -0.0043

     Molecular orbitals for symmetry species   5

 Orbital          1        2
   1  C1  2pz    0.2592  -0.2442
   2  C1  2pz    0.0764  -0.2045
   3  C1  3d1-  -0.0106  -0.0154
   4  C1  3d1+  -0.0061   0.0190
   5  C6  2pz    0.2592   0.4884
   6  C6  2pz    0.0764   0.4090
   7  C6  3d1+  -0.0123   0.0076
   8  H1  2pz    0.0052  -0.0062
   9  H6  2pz    0.0052   0.0124

     Molecular orbitals for symmetry species   6

 Orbital          1        2
   1  C1  2pz    0.2034   0.4467
   2  C1  2pz    0.0856   0.3244
   3  C1  3d1-  -0.0123   0.0070
   4  C1  3d1+   0.0129   0.0040
   5  C6  2pz    0.4068  -0.4467
   6  C6  2pz    0.1712  -0.3244
   7  C6  3d1+  -0.0085  -0.0080
   8  H1  2pz    0.0045   0.0079
   9  H6  2pz    0.0090  -0.0079

     Molecular orbitals for symmetry species   7

 Orbital          1
   1  C1  2pz    0.3523
   2  C1  2pz    0.1483
   3  C1  3d1-  -0.0014
   4  C1  3d1+  -0.0123
   5  C6  3d1-   0.0200
   6  H1  2pz    0.0078

     Molecular orbitals for symmetry species   8

 Orbital          1
   1  C1  2pz    0.4230
   2  C1  2pz    0.3542
   3  C1  3d1-  -0.0013
   4  C1  3d1+   0.0154
   5  C6  3d1-   0.0279
   6  H1  2pz    0.0107


 Printout of CI-coefficients larger than 0.05000 for root  1

 (this is the reference state)



  Printout of coefficients in interval  0.3162E+00 to 0.1000E+01
  ==============================================================

  Coefficient of determinant        2 is 0.66020646E+00
  alpha-string:   1   3   5   6
   beta-string:   1   3

  Coefficient of determinant       46 is 0.66019053E+00
  alpha-string:   1   2   3   5
   beta-string:   1   5


  Printout of coefficients in interval  0.1000E+00 to 0.3162E+00
  ==============================================================

  Coefficient of determinant       26 is -.12299576E+00
  alpha-string:   1   2   5   6
   beta-string:   1   2

  Coefficient of determinant       30 is -.19292282E+00
  alpha-string:   1   3   4   5
   beta-string:   3   5

  Coefficient of determinant       40 is -.12299543E+00
  alpha-string:   1   2   3   6
   beta-string:   1   6


  Printout of coefficients in interval  0.5000E-01 to 0.1000E+00
  ==============================================================

  Coefficient of determinant        1 is -.55832699E-01
  alpha-string:   1   2   3   4
   beta-string:   1   3

  Coefficient of determinant        4 is 0.55835947E-01
  alpha-string:   1   4   5   6
   beta-string:   1   3

  Coefficient of determinant        7 is 0.54074305E-01
  alpha-string:   1   3   5   6
   beta-string:   2   3

  Coefficient of determinant        9 is 0.64465145E-01
  alpha-string:   1   4   5   6
   beta-string:   2   3

  Coefficient of determinant       13 is 0.74889046E-01
  alpha-string:   2   3   5   6
   beta-string:   1   4

  Coefficient of determinant       21 is 0.64463767E-01
  alpha-string:   1   2   3   4
   beta-string:   5   6

  Coefficient of determinant       22 is -.54071125E-01
  alpha-string:   1   3   5   6
   beta-string:   5   6

  Coefficient of determinant       47 is -.55834224E-01
  alpha-string:   1   2   4   5
   beta-string:   1   5

  Coefficient of determinant       48 is -.55834198E-01
  alpha-string:   1   3   4   6
   beta-string:   1   5

  Coefficient of determinant       50 is -.54071795E-01
  alpha-string:   1   2   3   5
   beta-string:   2   5

  Coefficient of determinant       51 is 0.64982850E-01
  alpha-string:   1   2   4   5
   beta-string:   2   5

  Coefficient of determinant       54 is -.54072296E-01
  alpha-string:   1   2   3   5
   beta-string:   3   6

  Coefficient of determinant       56 is 0.64982956E-01
  alpha-string:   1   3   4   6
   beta-string:   3   6


  Norm of printed CI vector ..  0.98572956E+00

   Magnitude of CI coefficients 
  ==============================

   ( Ranges are relative to norm of vector :  1.00E+00 )

  10- 1 to 10- 0         5    0.93919915E+00    0.93919915E+00
  10- 2 to 10- 1        28    0.60371848E-01    0.99957100E+00
  10- 3 to 10- 2        16    0.42900066E-03    0.10000000E+01
  10- 6 to 10- 5         4    0.24567789E-10    0.10000000E+01
  10- 7 to 10- 6         7    0.69533804E-12    0.10000000E+01
  10- 8 to 10- 7         1    0.56007590E-14    0.10000000E+01
  Number of coefficients less than 10^-11 times norm is             0



 >>>> Total CPU  time used in SIRIUS :    154.58 seconds
 >>>> Total wall time used in SIRIUS :    573.00 seconds

 
     Date and time (Linux)  : Thu Jun  7 16:32:03 2001
     Host name              : ghiaurov.pdc.kth.se                     

- End of Wave Function Section


  This is output from RESPONS  -  an MCSCF and SOPPA response property program 
 ------------------------------------------------------------------------------



 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>



 ********* ESRINP ********

 *ESR   input ignored because no operators requested.


 MCSCF energy         :     -230.669917814688100
 -- inactive part     :     -428.976925555991400
 --   active part     :       -6.317356385208711
 -- nuclear repulsion :      204.624364126512100


 Output from ZFSDRV
 ------------------



 Trace-less zero field splitting tensor (cm-1)
 ---------------------------------------------


               Column   1     Column   2     Column   3
       1      -0.05304940     0.00000000     0.00000000
       2       0.00000000    -0.05304472     0.00000000
       3       0.00000000     0.00000000     0.10609412


 ZFS energy eigenvalues (cm-1)
 -----------------------------


               Column   1
       1      -0.10609412
       2       0.05304472
       3       0.05304940


 Zero field splitting paramters
 ------------------------------

@ZFS parameter D =   0.159141 cm-1 =    4770.93 MHz
@ZFS parameter E =   0.000002 cm-1 =       0.07 MHz


 ZFS tensor eigenvalues and cosines
 ----------------------------------

            cm-1             MHz       direction cosines
       -0.053049    -1590.381025  1.0000  0.0000  0.0000
       -0.053045    -1590.240750  0.0000  1.0000  0.0000
        0.106094     3180.621775  0.0000  0.0000  1.0000

 >>>> Total CPU  time used in RESPONSE:  5 minutes 50 seconds
 >>>> Total wall time used in RESPONSE:  6 minutes 27 seconds
 >>>> Total CPU  time used in DALTON:  9 minutes 36 seconds
 >>>> Total wall time used in DALTON: 17 minutes 37 seconds

 
     Date and time (Linux)  : Thu Jun  7 16:38:31 2001
     Host name              : ghiaurov.pdc.kth.se                     
END REFOUT
