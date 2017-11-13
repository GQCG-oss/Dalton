!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
C
C  /* Deck herctl */
#ifndef PRG_DIRAC
      SUBROUTINE HERCTL(WORK,LWORK)
C
C     The code for calculation of one- and two-electron integrals
C     was written by T. Helgaker in 1984 at the University of Oslo.
C
C     General contractions were implemented by T. Helgaker at the
C     University of Aarhus in Feb-Mar 1988.
C
C     Symmetry was implemented by P. R. Taylor and T. Helgaker at
C     NASA Ames in Apr 1988.
C
C     The supermatrix code is written by O. Kvalheim at the University
C     of Bergen.
C
C     Spin-orbit integrals by O. Vahtras and T. Helgaker at the
C     University of Oslo, Nov 1989.
C
C     Cartesian and spherical moments integrals by T. Helgaker,
C     University of Oslo, Sep and Oct 1990
C
C     Integrals for indirect spin-spin coupling tensors by T. Helgaker
C     and O. Vahtras at the University of Oslo, Feb 1991.
C
C     Half-derivative overlap integrals for NACMES by T. Helgaker,
C     at University of Aarhus, Jun 1991
C
C     Overlap matrix differentiated with respect to external magetic
C     field, K. Ruud & T. Helgaker, Oct 1991
C
C     Electronic angular momentum around fixed center or nuclei,
C     K. Ruud & T. Helgaker, Oct 1991
C
C     One-electron contribution to the magnetic moment of molecules,
C     K. Ruud & T. Helgaker, Nov. 1991
C
C     Kinetic energy integrals, K. Ruud, Nov. 1991
C
C     Cosine and sine integrals, T. Helgaker, Jun. 1993
C
C     Mass-velocity and Darwin integrals
C     S. Kirpekar & H.J.Aa. Jensen, Jul. 1993
C
C     Magnetic field derivative integrals of dipole length
C     K.Ruud, Aug.-93
C

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qm3.h"
#include "maxorb.h"
C
#include "cbiher.h"
#include "cbihr1.h"
#include "cbihr2.h"
#include "cbihrs.h"
#include "cbieri.h"
#include "gnrinf.h"
#include "exeinf.h"
#include "inftap.h"
#include "huckel.h"
#include "r12int.h"
#include "nuclei.h"
#include "pcmlog.h"
C
Cholesky
#include "ccdeco.h"
Cholesky
C
C inftra : USE_INTSORT
C
#include "inftra.h"
C

      LOGICAL SET, FEXIST
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('HERCTL')
      CALL GETTIM(TIMHER,WALHER)
C
C     *************************
C     ***** Input Section *****
C     *************************
C
      TIMSTR = SECOND()
      CALL HERINP(WORK,LWORK)
      IF (TSTINP) THEN
         WRITE (LUPRI,'(/15X,A)')
     &        '*** End of input test for HERMIT ***'
         CALL QUIT('*** End of input test for HERMIT ***')
      END IF
      TIMINP = SECOND() - TIMSTR
      CALL FLSHFO(LUPRI)
C
C     ************************************
C     ***** Calculation of Integrals *****
C     ************************************
C

      CALL GPINQ('AOPROPER','EXIST',FEXIST)
      IF (FEXIST) THEN
         IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ',' ',
     &                                  IDUMMY,.FALSE.)
         CALL GPCLOSE(LUPROP,'DELETE')
      END IF
C
      IF (QM3) THEN
         QM3LO1 = .TRUE.
         QM3LO2 = .TRUE.
         LONUPO = .FALSE.
         LOELFD = .FALSE.
         CALL GPOPEN(LUQM3E,'ELFDMM','UNKNOWN',' ',' ',
     &               IDUMMY,.FALSE.)
         CALL GPOPEN(LUQM3P,'POTMM','UNKNOWN',' ',' ',
     &               IDUMMY,.FALSE.)
         FORQM3 = .TRUE.
C        make sure QM3 integrals are evaluated for relevant properties;
C        is reset to .false. before return (is used in other routines
C        later to signal when integrals are for QM3 and when not).
      ELSE
         FORQM3 = .FALSE.
      END IF
C
      RUNQM3 = .FALSE.
C
      IF (DORLM) CALL GPOPEN(LUSOL,
     &     'AOSOLINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)

C     Delete old MO integrals integrals
      IF (USE_INTSORT) THEN
         IF (LUORDA .LE. 0) CALL DAOPEN(LUORDA,'AOORDINT')
                            CALL DARMOV(LUORDA)
         IF (LUMINT .LE. 0) CALL DAOPEN(LUMINT,'MOTWOINT')
                            CALL DARMOV(LUMINT)
      ELSE
         CALL GPINQ('AOSRTINT.DA','EXIST',FEXIST)
         IF (FEXIST) THEN
            CALL SYSTEM('rm -f AOSRTINT.DA*')
         END IF
         CALL GPINQ('MOTWOINT','EXIST',FEXIST)
         IF (FEXIST) THEN
            CALL GPOPEN(LUINTM,'MOTWOINT',
     &           'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
            CALL GPCLOSE(LUINTM,'DELETE')
         END IF
      END IF

C     Delete any existing old 2-el. AO integrals
C     (they may have been generated with the MAKE_AOTWOINT routine,
C      even if RUNTWO is .false. now ! This will for example be
C      the case if DIRCAL and DOMP2))

      IF (LUINTA  .LE. 0) CALL GPOPEN(LUINTA,'AOTWOINT',
     &        'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
Casm
C     Let's give a chance to the user who konows what he does
C
C     CALL GPCLOSE(LUINTA,'DELETE')
      IF (RMAOTWO) THEN
         CALL GPCLOSE(LUINTA,'DELETE')
      ELSE
         CALL GPCLOSE(LUINTA,'KEEP')
      END IF
C
      IF (DOSRIN) THEN
         IF (LUSRINT .LE. 0) CALL GPOPEN(LUSRINT,'AO2SRINT',
     &        'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL GPCLOSE(LUSRINT,'DELETE')
      END IF
      IF (LUSUPM .LE. 0) CALL GPOPEN(LUSUPM,'AOSUPINT',
     &        'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPCLOSE(LUSUPM,'DELETE')

C     IF (RUNTWO .or. RUNERI) :
C        Open LUINTA and LUSRINT for new 2-el. integrals
C        to be calculated in the CALL HERINT call below

      IF ((RUNTWO .OR. RUNERI) .AND. .NOT. DCCR12) THEN
         CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         IF (DOSRIN) CALL
     &      GPOPEN(LUSRINT,'AOSR2INT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      END IF
C
      CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ',' ',IDUMMY,.FALSE.)
      CALL MOLLAB('INTEGRAL',LUONEL,LUPRI)
      CALL HERINT(WORK,LWORK)

      IF (PCM) THEN
         CALL IEFMAT(WORK,LWORK)
      END IF
C
      IF ( (QM3) .AND. (.NOT.INTDIR) ) THEN
        DO 752 KI = 1, ISYTP
        IF (MDLWRD(KI)(1:3) .EQ. 'SPC') THEN
          IF ( .NOT. (LONUPO) ) THEN
            CALL QUIT('Error, NUCPOT integrals not computed')
          END IF
        END IF
  752   CONTINUE
C
        DO 753 KI = 1, ISYTP
        IF (MDLWRD(KI)(1:5) .EQ. 'SPC_E') THEN
          IF ( .NOT. (LOELFD) ) THEN
            CALL QUIT('Error, NELFLD integrals not computed')
          END IF
        END IF
  753   CONTINUE
      END IF

#if defined(BUILD_GEN1INT)
C     test suite of Gen1Int interface, added by Bin Gao, Oct. 7, 2011
      if (TEST_GEN1INT)
     &  call gen1int_host_test(LWORK, WORK, LUPRI, IPRDEF)
#endif

C
C     **********************************
C     ***** P Supermatrix Elements *****
C     **********************************
C
C     Integrals on AOSUPINT
C     HJAaJ Oct 2003: CALL FORMSUP has been moved to SIRIUS and
C     is only done when relevant ...
C
C     **********************************
C     ***** Integral sorting       *****
C     **********************************
C
C     Integrals on AOORDINT
C
      IF (USE_INTSORT .AND. CHOINT) THEN
         WRITE(LUPRI,'(/A)')
     &      ' * INFO: .SORT INT ignored when .CHOLESKY'
         USE_INTSORT = .FALSE.
      END IF
      IF (USE_INTSORT) THEN
         CALL GPINQ('AOTWOINT','EXIST',FEXIST)
         IF (FEXIST) THEN
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL SORTAO(WORK,LWORK)
            CALL TIMER('SORTAO',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         ELSE
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/A)')
     &      ' * WARNING: .SORT INT requested but no integrals to sort!'
         END IF
      END IF
C
C     ******************************
C     ***** End of Calculation *****
C     ******************************
C
      CALL GPCLOSE(LUONEL,'KEEP')
      IF (LUPROP  .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      IF (LUQM3E  .GT. 0) CALL GPCLOSE(LUQM3E,'KEEP')
      IF (LUQM3P  .GT. 0) CALL GPCLOSE(LUQM3P,'KEEP')
      IF (LUSOL   .GT. 0) CALL GPCLOSE(LUSOL ,'KEEP')
C
      IF (LUINTA  .GT. 0) CALL GPCLOSE(LUINTA ,'KEEP')
      IF (LUSRINT .GT. 0) CALL GPCLOSE(LUSRINT,'KEEP')
C
      IF (RUNTWO .OR. RUNERI) FTRCTL = .TRUE.
C     ..transformation control, old AO/MO files cannot be used any more.
C
      CALL GETTIM(TEND,WEND)
      TIMHER = TEND - TIMHER
      WALHER = WEND - WALHER
      CALL TIMTXT(' Total CPU  time used in HERMIT:',TIMHER,LUPRI)
      CALL TIMTXT(' Total wall time used in HERMIT:',WALHER,LUPRI)
      CALL QEXIT('HERCTL')
      RETURN
      END
#endif
C  /* Deck herinp */
      SUBROUTINE HERINP(WORK,LWORK)
C
C     --- General Input for HERMIT ---
C
C     Trygve Helgaker, extensively extended by others 
C
C     950602-vebjorn:
C     Added flag HRINPC to ensure that HERMIT input processing is done only once.
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
C     NTABLE has been increased 118 (WK/UniKA/07-25-2005).
      PARAMETER (NDIR = 9, NTABLE = 146)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D4 = 4.0D0, D10 = 1.0D1)
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7
      DIMENSION WORK(LWORK)
C
#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "hrunit.h"
#include "cbieri.h"
#else
C
C inftra : USE_INTSORT
C
#include "dummy.h"
#include "gnrinf.h"
#include "cbihrs.h"
#include "inftap.h"
#include "infinp.h"
#include "inftra.h"
#endif
C
#include "cbiher.h"
#include "cbihr1.h"
#include "nuclei.h"
#include "ccom.h"
#include "orgcom.h"

#include "efield.h"
#include "cbisol.h"
#include "qm3.h"
#include "qmmm.h"
#include "huckel.h"
#include "drw2el.h"
#include "elweak.h"
#include "r12int.h"
#include "codata.h"
#include "denfit.h"
#include "infpar.h"

C     Keywords for R12 have been added (WK/UniKA/04-11-2002).
cLig <> added RANGMO and RPSO
C Bin Gao, December 17, 2009; adds integrals S2MBRA, S2MKET, and S2MMIX
      DATA TABDIR /'*END OF', '*READIN', '*ONEINT', '*TWOINT',
     &             '*SUPINT', '*ER2INT', '*SORINT', '*DENFIT',
     &             '*QM3INP'/
C     adds test suite of Gen1Int interface as term 239
C     Bin Gao, Oct. 3, 2011
!                    101        102        103        104
      DATA TABLE  /'.PRINT ', '.INPTES', '.NOSUP ', '.SPIN-O',
!                    105        106        107        108
     &             '.DIPLEN', '.NO HAM', '.SOTEST', '.DIPVEL',
!                    109        110        111        112
     &             '.QUADRU', '.PHASEO', '.SECMOM', '.SUPONL',
!                    113        114        115        116
     &             '.CARMOM', '.SPHMOM', '.FC    ', '.PSO   ',
!                    117        118        119        120
     &             '.SD    ', '.DSO   ', '.POINTS', '.SELECT',
!                    121        122        123        124
     &             '.QUASUM', '.SD+FC ', '.PROPRI', '.HDO   ',
!                    125        126        127        128
     &             '.S1MAG ', '.S2MAG ', '.ANGMOM', '.ANGLON',
!                    129        130        131        132
     &             '.LONMOM', '.MAGMOM', '.S1MAGT', '.MGMOMT',
!                    133        134        135        136
     &             '.KINENE', '.S2MAGT', '.DSUSNL', '.DSUSLL',
!                    137        138        139        140
     &             '.DSUSLH', '.DIASUS', '.DSUTST', '.NSTNOL',
!                    141        142        143        144
     &             '.NSTLON', '.NST   ', '.NSNLTS', '.NSLTST',
!                    145        146        147        148
     &             '.NELFLD', '.NSTTST', '.EFGCAR', '.EFGSPH',
!                    149        150        151        152
     &             '.S1MAGL', '.S1MAGR', '.HDOBR ', '.S1MLT ',
!                    153        154        155        156
     &             '.HDOBRT', '.S1MRT ', '.NUCPOT', '.NPOTST',
!                    157        158        159        160
     &             '.MGMO2T', '.MGMTHR', '.HBDO  ', '.SUSCGO',
!                    161        162        163        164
     &             '.NSTCGO', '.EXPIKR', '.MASSVE', '.DARWIN',
!                    165        166        167        168
     &             '.CM-1  ', '.CM-2  ', '.SQHDOL', '.SQHDOR',
!                    169        170        171        172
     &             '.NOTWO ', '.GAUGEO', '.DIPORG', '.NO2SO ',
!                    173        174        175        176
     &             '.S1ELE ', '.S1ELB ', '.ONEELD', '.THETA ',
!                    177        178        179        180
     &             '.NUCMOD', '.SORT I', '.DIPGRA', '.QUAGRA',
!                    181        182        183        184
     &             '.OCTGRA', '.ROTSTR', '.THIRDM', '.SOFIEL',
!                    185        186        187        188
     &             '.SOMAGM', '.DEROVL', '.DERHAM', '.ELGDIA', 
!                    189        190        191        192
     &             '.ELGDIL', '.MNF-SO', '.DPTOVL', '.DPTPOT', 
!                    193        194        195        196
     &             '.XDDXR3', '.AD2DAR', '.PVIOLA', '.WEINBG', 
!                    197        198        199        200
     &             '.FINDPT', '.FNPROP', '.PVP   ', '.1ELPOT',
!                    201        202        203        204
     &             '.QDBINT', '.QDBTST', '.R12   ', '.R12GTG',
!                    205        206        207        208
     &             '.AUXBAS', '.NOTV12', '.R12INT', '.U12INT',
!                    209        210        211        212
     &             '.U21INT', '.DCCR12', '.RANGMO', '.RPSO  ',
!                    213        214        215        216
     &             '.DPTPXP', '.NOPICH', '.OZ-KE ', '.PSO-KE', 
!                    217        218        219        220
     &             '.DNS-KE', '.SD-KE ', '.FC-KE ', '.DSO-KE', 
!                    221        222        223        224
     &             '.PSO-OZ', '.SQHD2O', '.R12STG', '.R12RSG',
!                    225        226        227        228
     &             '.R12RGG', '.R12ERF', '.R12RRF', '.DIPLOC',
!                    229        230        231        232
     &             '.2-EL D', '.EFBDER', '.EFBDR2', '.MAGQDP',
!                    233        234        235        236
     &             '.MQDPTS', '.ORB-OR', '.DERAM ', '.DIPANH',
!                    237        238        239        240
     &             '.NEWTRA', '.KINADI', '.GENINT', '.S2MBRA',
!                    241        242        243        244
     &             '.S2MKET', '.S2MMIX', '.LRINTS', '.NORMAO',
!                    245        246
     &             '.TWOBYT', '.LFDIPL'/
C
C
      CALL QENTER('HERINP')
C
      RMAOTWO = .TRUE.
C
      IPRDEF = IPRUSR + 1
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
C     Check if input has been processed earlier.
C     HRINPC variable makes sure HERMIT input is only read once
C
      IF (HRINPC) GOTO 1000
      HRINPC = .TRUE.
C
C     *** Initialize ***
C
      ! CALL HERINI - Apr. 2011: HERINI now called in EXEHER /hjaaj
      !               so Hermit variables are also initialized if only
      !               ABACUS

      !Initialize defaults for density fitting
      CALL DENSFIT_INI
C
C     If WESTA - calculate integrals for WESTA program
C
      IF (WESTA) THEN
         SQHDOR = .TRUE.
         SQHD2O = .TRUE.
         DIPLEN = .TRUE.
         DIPVEL = .TRUE.
         ANGMOM = .TRUE.
         ONEPRP = .TRUE.
      END IF
C
      CALL TITLER('Output from **INTEGRALS input processing (HERMIT)',
     &   '*',118)
C
C     **** Find Hermit input *****
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',ERR=910,END=920) WORD1
      CALL UPCASE(WORD1)
      IF ((WORD1 .EQ. '**HERMI') .OR. (WORD1 .EQ. '*HERMIT') .OR.
     &    (WORD1 .EQ. '**INTEG')) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//A)')
     &   'WARNING **INTEGRALS input: error reading Dalton input file'
  920 CONTINUE
         WORD  = '**END O'
         WORD1 = '**END O'
         WRITE (LUPRI,'(/A)')
     &   ' - Using defaults, no **INTEGRALS input found'
         GOTO 300
  930 CONTINUE
C
      CALL TITLER('Output from HERMIT input processing','*',118)
C
C     ***** Process input for COMMON  /CBIHER/  *****
C
  100 READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124,125,126,127,128,129,130,
     &                131,132,133,134,135,136,137,138,139,140,
     &                141,142,143,144,145,146,147,148,149,150,
     &                151,152,153,154,155,156,157,158,159,160,
     &                161,162,163,164,165,166,167,168,169,170,
     &                171,172,173,174,175,176,177,178,179,180,
     &                181,182,183,184,185,186,187,188,189,190,
     &                191,192,193,194,195,196,197,198,199,200,
     &                201,202,203,204,205,206,207,208,209,210,
     &                211,212,213,214,215,216,217,218,219,220,
     &                221,222,223,224,225,226,227,228,229,230,
     &                231,232,233,234,235,236,237,238,239,240,
     &                241,242,243,244,245,246), I
            END IF
         END DO
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)')
     &         'ERROR: Keyword ',WORD,' not recognized under ',WORD1
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword under '//WORD1)
  101    CONTINUE
            READ (LUCMD,*) IPRDEF
            GO TO 100
  102    CONTINUE
            TSTINP = .TRUE.
            GO TO 100
  103    CONTINUE
            SUPMAT = .FALSE.
            GO TO 100
  104    CONTINUE
            SPNORB = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  105    CONTINUE
            DIPLEN = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  106    CONTINUE
            HAMILT = .FALSE.
            SUPMAT = .FALSE.
            GO TO 100
  107    CONTINUE
            SOTEST = .TRUE.
            GO TO 100
  108    CONTINUE
            DIPVEL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  109    CONTINUE
            QUADRU = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  110    CONTINUE
            READ (LUCMD, *,IOSTAT=IOS) (ORIGIN(I),I=1,3)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error in reading (ORIGIN(I),I=1,3)!'
              WRITE(LUPRI,*) '(ORIGIN(I),I=1,3):', (ORIGIN(I),I=1,3)
              CALL QUIT('Error in reading (ORIGIN(I),I=1,3)!')
            ENDIF
            GO TO 100
  111    CONTINUE
            SECMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  112    CONTINUE
            SUPMAT = .TRUE.
            HAMILT = .FALSE.
            NOTWO  = .TRUE.
            ONEPRP = .FALSE.
            GO TO 100
  113    CONTINUE
            CARMOM = .TRUE.
            READ (LUCMD,*) IORCAR_input
            IORCAR = MAX(IORCAR,IORCAR_input)
            ONEPRP = .TRUE.
            GO TO 100
  114    CONTINUE
            SPHMOM = .TRUE.
            READ (LUCMD,*) IORSPH
            ONEPRP = .TRUE.
            GO TO 100
  115    CONTINUE
            FERMI = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  116    CONTINUE
            PSO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  117    CONTINUE
            SPIDIP = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  118    CONTINUE
            DSO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  119    CONTINUE
            READ (LUCMD,*) NPQUAD
            GO TO 100
  120    CONTINUE
            READ (LUCMD, *) NPATOM
            IF (NPATOM .GT. MXCENT) THEN
               WRITE (LUPRI,'(//2A/A,I8/A,I6)')
     &             'ERROR: Too many atoms selected under ',WORD,
     &             'ERROR: Number of atoms selected:',NPATOM,
     &             'ERROR: Number of atoms allowed: ',MXCENT
               CALL QUIT('Error in '//WORD1)
            END IF
            READ (LUCMD, *) (IPATOM(I),I=1,NPATOM)
            ALLATM = .FALSE.
            GO TO 100
  121    CONTINUE
            TRIANG = .FALSE.
            GO TO 100
  122    CONTINUE
            SDFC = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  123    CONTINUE
            PROPRI = .TRUE.
            GO TO 100
  124    CONTINUE
            HDO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  125    CONTINUE
            S1MAG = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  126    CONTINUE
            S2MAG = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  127    CONTINUE
            ANGMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  128    CONTINUE
            ANGLON = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  129    CONTINUE
            LONMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  130    CONTINUE
            MAGMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  131    CONTINUE
            S1MAG  = .TRUE.
            S1MAGT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  132    CONTINUE
            MGMOMT = .TRUE.
            LONMOM = .TRUE.
            HAMILT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  133    CONTINUE
            KINENE = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  134    CONTINUE
            S2MAG  = .TRUE.
            S2MAGT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  135    CONTINUE
            DSUSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  136    CONTINUE
            DSUSLL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  137    CONTINUE
            DSUSLH = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  138    CONTINUE
            DIASUS = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  139    CONTINUE
            DSUTST = .TRUE.
            DSUSLL = .TRUE.
            ANGLON = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  140    CONTINUE
            NUCSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  141    CONTINUE
            NUCSLO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  142    CONTINUE
            NUCSHI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  143    CONTINUE
            NSNLTS = .TRUE.
            NUCSLO = .TRUE.
            PSO    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  144    CONTINUE
            NSLTST = .TRUE.
            NELFLD = .TRUE.
            NUCSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  145    CONTINUE
            NELFLD = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  146    CONTINUE
            NSTTST = .TRUE.
            NUCSLO = .TRUE.
            NUCSNL = .TRUE.
            NUCSHI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  147    CONTINUE
            EFGCAR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  148    CONTINUE
            EFGSPH = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  149    CONTINUE
            S1MAGL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  150    CONTINUE
            S1MAGR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  151    CONTINUE
            HDOBR  = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  152    CONTINUE
            S1MLT  = .TRUE.
            S1MAGL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  153    CONTINUE
            HDOBR  = .TRUE.
            HDOBRT = .TRUE.
            DIPVEL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  154    CONTINUE
            S1MRT  = .TRUE.
            S1MAGR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  155    CONTINUE
            NUCPOT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  156    CONTINUE
            NUCPOT = .TRUE.
            HAMILT = .TRUE.
            NPOTST = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  157    CONTINUE
            MGMO2T = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  158    CONTINUE
            READ (LUCMD,*) PRTHRS
            GO TO 100
  159    CONTINUE
            HBDO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  160    CONTINUE
            SUSCGO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  161    CONTINUE
            NSTCGO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  162    CONTINUE
            EXPIKR = .TRUE.
            ONEPRP = .TRUE.
            READ (LUCMD,*) (EXPKR(I),I=1,3)
            GOTO 100
  163    CONTINUE
            MASSVL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  164    CONTINUE
            DARWIN = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  165    CONTINUE
            READ (LUCMD,'(A7)') FIELD1
            IF (.NOT. ((FIELD1 .EQ. 'X-FIELD')
     &          .OR.   (FIELD1 .EQ. 'Y-FIELD')
     &          .OR.   (FIELD1 .EQ. 'Z-FIELD')
     &          .OR.   (FIELD1 .EQ. 'XYZ-ALL'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD1,
     &               '" illegal'
               CALL QUIT('Illegal field directions for CM-1 integrals')
            END IF
            CM1    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  166    CONTINUE
            READ (LUCMD,'(A7)') FIELD2
C Modified by Bin Gao, December 14, 2009
C adding XYZ-ALL
            IF (.NOT. ((FIELD2 .EQ. 'X-FIELD')
     &          .OR. (FIELD2 .EQ. 'Y-FIELD')
     &          .OR. (FIELD2 .EQ. 'Z-FIELD')
     &          .OR. (FIELD2 .EQ. 'XYZ-ALL'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD2,
     &               '" illegal'
               CALL QUIT('Illegal field directions for CM-2 integrals')
            END IF
            CM2    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  167    CONTINUE ! .SQHDOL - Bra differentiated half-derivative overlap matrix
            SQHDOL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  168    CONTINUE ! .SQHDOR - Ket differentiated half-derivative overlap matrix
            SQHDOR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  169    CONTINUE
            NOTWO  = .TRUE.
            SUPMAT = .FALSE.
            GO TO 100
  170    CONTINUE
            READ (LUCMD, *, IOSTAT=IOS) (GAGORG(I),I=1,3)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error in reading (GAGORG(I),I=1,3)!'
              WRITE(LUPRI,*) '(GAGORG(I),I=1,3):', (GAGORG(I),I=1,3)
              CALL QUIT('Error in reading (GAGORG(I),I=1,3)!')
            ENDIF
            GO TO 100
  171    CONTINUE
            READ (LUCMD, *, IOSTAT=IOS) (DIPORG(I),I=1,3)
            IF (IOS.NE.0) THEN
              WRITE(LUPRI,*) 'Error in reading (DIPORG(I),I=1,3)!'
              WRITE(LUPRI,*) '(DIPORG(I),I=1,3):', (DIPORG(I),I=1,3)
              CALL QUIT('Error in reading (DIPORG(I),I=1,3)!')
            ENDIF
            GO TO 100
  172    CONTINUE
c ach
            no2so=.true.
            GOTO 100
  173    CONTINUE
            S1ELE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  174    CONTINUE
            S1ELB  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  175    CONTINUE
            ONEELD = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  176    CONTINUE
            THETA  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  177    CONTINUE
C&&& Nuclear model
            READ (LUCMD, *,IOSTAT=IOS) INUC
            IF (IOS.NE.0) THEN
              CALL QUIT('Error in reading nuclear model INUC  !')
            ENDIF
            IF    (INUC.EQ.1) THEN
              GAUNUC = .FALSE.
            ELSEIF(INUC.EQ.2) THEN
              GAUNUC = .TRUE.
            ELSE
              WRITE(LUPRI,'(//3A,I5)')
     &         '*** ERROR ',WORD1,' *** Unknown nuclear model:',INUC
              CALL QUIT('*** ERROR '//WORD1//': Unknown nuclear model')
            ENDIF
            GO TO 100
!      .SORT Integrals
  178    CONTINUE
            SUPMAT = .FALSE.
            USE_INTSORT = .TRUE.
            NEWTRA = .TRUE.
            GO TO 100
  179    CONTINUE
            DPLGRA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  180    CONTINUE
            QUAGRA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  181    CONTINUE
            OCTGRA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  182    CONTINUE
            ROTSTR = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100                   
  183    CONTINUE
            THRMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  184    CONTINUE
            SOFLD  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  185    CONTINUE
            SOMM   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  186    CONTINUE
            DEROVL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  187    CONTINUE
            DERHAM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100                            
  188    CONTINUE
            ELGDIA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  189    CONTINUE
            ELGDIL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  190    CONTINUE
            MNF_SO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  191    CONTINUE
            DPTOVL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  192    CONTINUE
            DPTPOT = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  193    CONTINUE
            XDDXR3 = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100        
  194    CONTINUE
            AD2DAR = .TRUE.
            READ (LUCMD, *) DARFAC
            GO TO 100                            
  195    CONTINUE
            PVIOLA = .TRUE.
            SPNORB = .TRUE.
            ONEPRP = .TRUE. 
            GO TO 100                            
  196    CONTINUE
            READ(LUCMD,*) BGWEIN
            BGWEIL = .TRUE.
            BGWEIN = D1 - D4*BGWEIN
            GOTO 100   
  197    CONTINUE
            ONEPRP = .TRUE. 
            DPTPOT = .TRUE.
            FINDPT = .TRUE.
            READ (LUCMD, *) DPTFAC
            GOTO 100 
  198    CONTINUE
            READ (LUCMD, *) X10FNP
            PVFINN = D10 ** X10FNP
            GOTO 100
  199    CONTINUE
            PVPINT = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100                            
  200    CONTINUE
            POTENE = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  201    CONTINUE
C Modified by Bin Gao, December 14, 2009: adds XYZ-ALL
            READ (LUCMD,'(A7)') FIELD3
            IF (.NOT. ((FIELD3 .EQ. 'XX-FGRD')
     &          .OR.   (FIELD3 .EQ. 'XY-FGRD')
     &          .OR.   (FIELD3 .EQ. 'XZ-FGRD')
     &          .OR.   (FIELD3 .EQ. 'YY-FGRD')
     &          .OR.   (FIELD3 .EQ. 'YZ-FGRD')
     &          .OR.   (FIELD3 .EQ. 'ZZ-FGRD')
     &          .OR.   (FIELD3 .EQ. 'XYZ-ALL'))) THEN
               WRITE (LUPRI,'(/4A/)') 'ERROR: Field direction "',FIELD3,
     &               '" illegal for input ',WORD
               CALL QUIT('Illegal field direction for QDBINT integrals')
            END IF
            QDBINT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  202    CONTINUE
            QDBTST = .TRUE.
            THETA  = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  203    CONTINUE
C           R12 calculation requested with linear r12 terms (WK/UniKA/04-11-2002).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            GOTO 100
  204    CONTINUE
C           R12 calculation requested with Gaussian-type geminals (WK/UniKA/03-24-2005).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
C           SLATER = .TRUE.
            READ(LUCMD,*) NTOGAM
            IF (NTOGAM .LE. 0 .OR. NTOGAM .GT. MAXGAM)
     &      CALL QUIT('Invalid number of Gaussian-damped terms.')
            DO NTGAM = 1, NTOGAM
              READ(LUCMD,*) GAMAB(NTGAM), GAMAX(NTGAM)
            END DO
            GOTO 100
  205    CONTINUE
C           Calculation with multiple basis sets requested (WK/UniKA/04-11-2002).
            CALL QUIT('ERROR, .AUXBAS has been renamed to .R12AUX'//
     &        ' and moved to *MOLBAS')
            GOTO 100
  206    CONTINUE
C           No calculation of 1/r12 integrals (WK/UniKA/04-11-2002).
            NOTV12 = .TRUE.
            V12INT = .FALSE.
            GOTO 100
  207    CONTINUE
C           Calculation of integrals of the operator r12 (WK/UniKA/04-11-2002).
            R12INT = .TRUE.
            GOTO 100
  208    CONTINUE
C           Calculation of integrals of the operator U12 (WK/UniKA/04-11-2002).
            U12INT = .TRUE.
            GOTO 100
  209    CONTINUE
C           Calculation of integrals of the operator U21 (WK/UniKA/04-11-2002).
            U21INT = .TRUE.
            GOTO 100
  210    CONTINUE
C           Integrals are computed for DIRCCR12 program (WK/UniKA/25-11-2002).
            DCCR12 = .TRUE.
            SUPMAT = .FALSE.
            GOTO 100
cLig added what to do for RANGMO and RPSO options
  211    CONTINUE
            RANGMO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  212    CONTINUE
            RPSO   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  213    CONTINUE
            PXPINT = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  214    CONTINUE
            NOPICH = .TRUE.
            GO TO 100
  215    CONTINUE
            OZKE   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  216    CONTINUE
            PSOKE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  217    CONTINUE
            DNSKE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100

  218    CONTINUE
            SDKE   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  219    CONTINUE
            FCKE   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  220    CONTINUE
            DSOKE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  221    CONTINUE
            PSOOZ = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  222    CONTINUE ! .SQH2DO - Ket second differentiated half-derivative overlap matrix
            SQHD2O = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  223    CONTINUE
C           R12 calculation requested with fitted Slater-type geminal (WK/UniKA/03-24-2005).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
            SLATER = .TRUE.
            READ(LUCMD,*) ALPSTG
            READ(LUCMD,*) NTOGAM
            IF (NTOGAM .LE. 0 .OR. NTOGAM .GT. MAXGAM)
     &      CALL QUIT('Invalid number of Gaussians.')
            CALL STGFIT(ALPSTG,NTOGAM,GAMAB,GAMAX)
            GOTO 100
  224    CONTINUE
C           R12 calculation requested with fitted "r12-times-Slater" geminal (WK/UniKA/03-24-2005).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
            READ(LUCMD,*) ALPSTG
            READ(LUCMD,*) NTOGAM
            IF (NTOGAM .LE. 0 .OR. NTOGAM .GT. MAXGAM)
     &      CALL QUIT('Invalid number of Gaussians.')
            CALL RSTGFT(ALPSTG,NTOGAM,GAMAB,GAMAX)
            GOTO 100
  225    CONTINUE
C           R12 calculation requested with Gaussian-damped r12 integrals (WK/UniKA/03-24-2005).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
            READ(LUCMD,*) NTOGAM
            IF (NTOGAM .LE. 0 .OR. NTOGAM .GT. MAXGAM)
     &      CALL QUIT('Invalid number of Gaussians.')
            DO NTGAM = 1, NTOGAM
              READ(LUCMD,*) GAMAB(NTGAM), GAMAX(NTGAM)
            END DO
            GOTO 100
  226    CONTINUE
C           R12 calculation requested with fitted ERFC-type geminal (WK/UniKA/03-29-2005).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
            SLATER = .TRUE.
            READ(LUCMD,*) ALPSTG
            READ(LUCMD,*) NTOGAM
            IF (NTOGAM .LE. 0 .OR. NTOGAM .GT. MAXGAM)
     &      CALL QUIT('Invalid number of Gaussians.')
            CALL ERFFIT(ALPSTG,NTOGAM,GAMAB,GAMAX)
            GOTO 100
  227    CONTINUE
C           R12 calculation requested with fitted "r12-times-ERFC" geminal (WK/UniKA/03-29-2005).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
            READ(LUCMD,*) ALPSTG
            READ(LUCMD,*) NTOGAM
            IF (NTOGAM .LE. 0 .OR. NTOGAM .GT. MAXGAM)
     &      CALL QUIT('Invalid number of Gaussians.')
            CALL RRFFIT(ALPSTG,NTOGAM,GAMAB,GAMAX)
            GOTO 100
  228    CONTINUE
            DIPLEN=.TRUE.
            ONEPRP=.TRUE.
            GO TO 100
  229    CONTINUE
            DAR2EL = .TRUE.
            GO TO 100
  230    CONTINUE
            EFBDER = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  231    CONTINUE
            EFB2DR = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  232    CONTINUE
            MAGQDP = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  233    CONTINUE
            MAGQDP = .TRUE.
            ANGMOM = .TRUE.
            ONEPRP = .TRUE.
            MQDPTS = .TRUE.
            GOTO 100
  234    CONTINUE
            ORBORB = .TRUE.
            GOTO 100
  235    CONTINUE
            DERAM  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  236    CONTINUE
            DIPANH = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
!        .NEWTRA ("new transformation" from 1988 (!)
!        this option is without calling SORTAO
  237    CONTINUE
            NEWTRA = .TRUE.
            GO TO 100
  238    CONTINUE
            KINADI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
  239    CONTINUE
#if defined(BUILD_GEN1INT)
            TEST_GEN1INT = .true.
#endif
            goto 100
C copied from linsca abacus/herdrv.F, Bin Gao, December 17, 2009
  240    CONTINUE
            S2MBRA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  241    CONTINUE
            S2MKET = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  242    CONTINUE
            S2MMIX = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
  243    CONTINUE
c jim-dbg : setting all integrals for LRESC module
            ONEPRP = .TRUE.
            DOLRINTS=.TRUE.
            FERMI = .TRUE.  ! 115.Fermi 
            KINENE = .TRUE. ! 133.Kinenerg
            SOFLD  = .TRUE. ! 184.somf
            DPTOVL = .TRUE. ! 191.
            SPIDIP = .TRUE. ! 117.Sd
            ANGMOM = .TRUE. ! 127.
            DARWIN = .TRUE. ! 164
            MASSVL = .TRUE. ! 163
            NSTCGO = .TRUE. ! 161.Dia
            PSO = .TRUE.    ! 116.Pso
            DIPVEL = .TRUE. ! 108.DipVel
            OZKE   = .TRUE. ! 215.Lkin
            PSOKE  = .TRUE. ! 216.PsoKin
            DNSKE  = .TRUE. ! 217.DiaKin
            PSOOZ  = .TRUE. ! 221. PSO-OZ for AngPso          
Cxx  el problema es que el origen de gauge debe estar en el atomo q elegimos, LRATOM. 
c       1 - Leer INP file y buscar LRATOM 
c       2 - hacer que GAGORG = CORD(LRATOM)
c Por ahora GAGORG = (0.0 0.0 0.0)
c  ACLARARLO EN OUT FILE !!
cx            write(lupri,*)' At INTEGRAL SECT to do GAUGEO on LRATOM' 
cx            write(lupri,*)' Calling lrscinp for SELECT' 
cx            CALL LRSCINP('.SELECT')
cx            write(lupri,*)' Despues. LRATOM ', LRATOM 
cx            GAGORG(1) = 0.000000
cx            GAGORG(2) = 0.000000
cx            GAGORG(3) = 0.000000
c =================================================
c -------------------------------------------------
css           x    READ (LUCMD,*) (LRGAUG(IS), IS = 1, 3)
css           x    DO ICENT = 1, NUCIND
c             x     NAME =  NAMEX(3*ICENT)(1:4)
css           x       WRITE (LUPRI,'(2X,A,3X," : ",3(A1,2X,A,F15.10))')
css     &     x       NAME, '1' , 'x' , CORD(1,ICENT),
css     &     x             '2' , 'y' , CORD(2,ICENT),
css     &     x             '3' , 'z' , CORD(3,ICENT)
css           x    ENDDO
cs            READ (LUCMD, *, IOSTAT=IOS) (GAGORG(I),I=1,3)
cs            IF (IOS.NE.0) THEN
cs              WRITE(LUPRI,*) 'Error in reading (GAGORG(I),I=1,3)!'
cs              WRITE(LUPRI,*) '(GAGORG(I),I=1,3):', (GAGORG(I),I=1,3)
cs              CALL QUIT('Error in reading (GAGORG(I),I=1,3)!')
cs            ENDIF
            GO TO 100
  244    CONTINUE
            RMAOTWO = .FALSE.
            GO TO 100 
  245    CONTINUE   ! .TWOBYTEPACKING
            TWOBYTEPACKING = .TRUE.
            GO TO 100 
  246    CONTINUE
            LFDIPLN = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 300
      ELSE
         WRITE (LUPRI,'(/3A/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in '//WORD1)
      END IF
  300 CONTINUE

      IF (.NOT. BGWEIL) BGWEIN = WEINBG 
cLig <>
      NMRISS = FERMI .OR. PSO .OR. SPIDIP .OR. DSO .OR. SDFC 
     &               .OR. RPSO
      IF (TSTINP) WRITE (LUPRI,'(/A/)') ' @@ Input test run only !!!'
      WRITE (LUPRI,'(/A,I5)') ' Default print level:    ',IPRDEF
      IF(GAUNUC) THEN
        WRITE(LUPRI,'(/A)') ' * Nuclear model: '//
     &    'Gaussian charge distribution.'
      ELSE
        WRITE(LUPRI,'(/A)') ' * Nuclear model: Point charge'
      ENDIF
      IF (SPNORB .AND. SOTEST) THEN
         WRITE (LUPRI,'(/A)')
     *    ' FATAL ERROR in input: .SPIN-ORBIT and .SOTEST cannot '/
     *    /'both be specified.'
          CALL QUIT('Error in **INTEGRALS input.')
      END IF
      IF (HAMILT) THEN
          IF (NOTWO) THEN
             WRITE (LUPRI,'(/A)')
     &          ' Calculation of one-electron Hamiltonian integrals.'
          ELSE
             WRITE (LUPRI,'(/A)')
     &    ' Calculation of one- and two-electron Hamiltonian integrals.'
          END IF
      ELSE
         IF (NOTWO) WRITE (LUPRI,'(/A)')
     &          ' Two-electron integrals not calculated.'
      END IF
      IF (R12INT) WRITE (LUPRI,'(A)')
     &    ' Calculation of two-electron integrals over r12.'
      IF (U12INT) WRITE (LUPRI,'(A)')
     &    ' Calculation of two-electron integrals over [T1,r12].'
      IF (U21INT) WRITE (LUPRI,'(A)')
     &    ' Calculation of two-electron integrals over [T2,r12].'
      IF (DAR2EL) WRITE (LUPRI,'(/A)')
     &     ' Calculate two-electron Darwin integrals'
      IF (AD2DAR) WRITE (LUPRI,'(/A/A,F18.14/)')
     &    ' Two-electron Darwin integrals are added to standard',
     &    ' repulsion integrals with perturbation parameter:',DARFAC
      IF (SPNORB .AND. .NOT.no2so) WRITE (LUPRI,'(/A)')
     &     ' Calculate two-electron spin-orbit integrals'
      IF (ORBORB) WRITE (LUPRI,'(/A)')
     &     ' Calculate two-electron orbit-orbit integrals'
      IF (TWOBYTEPACKING .AND. NBASIS .LE. 255) WRITE (LUPRI,'(/A)')
     &     ' (Two byte packing of indices on files requested.)'
      IF (FINDPT) THEN
         WRITE (LUPRI,'(/A/A,F18.14/A,F18.14,A/)')
     &    ' A direct relativistic perturbation is added to',
     &    ' Hamiltonian and metric with perturbation parameter:',DPTFAC,
     &    ' (perturbation parameter * c^{-2} = ',DPTFAC*ALPHAC**2,')'
         IF (NODTOT .GT. 0) CALL PARQUIT('FINDPT')
         IF (DIRCAL) CALL QUIT(
     &    'Direct calculation (.DIRECT) is not compatible with .FINDPT')
      END IF
      IF (ONEPRP) THEN
         WRITE (LUPRI,'(/A)')
     &      ' The following one-electron property integrals are'
     &    //' calculated as requested:'
                     WRITE (LUPRI,'(10X,A)') '- overlap integrals'
         IF (DIPLEN) WRITE (LUPRI,'(10X,A)') '- dipole length integrals'
         IF (DIPVEL) WRITE (LUPRI,'(10X,A)')
     &                         '- dipole velocity integrals'
C
         IF (QUADRU) WRITE (LUPRI,'(10X,A)')
     &                         '- quadrupole moment integrals'
         IF (THETA)  WRITE (LUPRI,'(10X,A)')
     &                         '- traceless quadrupole moment integrals'
         IF (SECMOM) WRITE (LUPRI,'(10X,A)')'- second moment integrals'
C
         IF (THRMOM) WRITE (LUPRI,'(10X,A)')'- third moment integrals'
C
         IF (CARMOM) THEN
            IF (IORCAR .GT. 0) THEN
               WRITE (LUPRI,'(10X,A,I2,A)')
     &          '- Cartesian multipole moment integrals of orders',
     &          ABS(IORCAR),' and lower'
            ELSE
               WRITE (LUPRI,'(10X,A,I2)')
     &          '- Cartesian multipole moment integrals of order',
     &           ABS(IORCAR)
            END IF
         END IF
         IF (SPHMOM) THEN
            IF (IORSPH .GT. 0) THEN
               WRITE (LUPRI,'(10X,A,I2,A)')
     &          '- Spherical multipole moment integrals of orders',
     &          ABS(IORSPH),' and lower'
            ELSE
               WRITE (LUPRI,'(10X,A,I2)')
     &          '- Spherical multipole moment integrals of order',
     &           ABS(IORSPH)
            END IF
         END IF
C
         IF (KINENE) WRITE (LUPRI,'(10X,A)')
     &      '- electronic kinetic energy'
         IF (KINADI) WRITE (LUPRI,'(10X,A)')
     &      '- adiabatic kinetic energy from nuclear derivatives'
         IF (MASSVL) WRITE (LUPRI,'(10X,A)')
     &       '- mass velocity integrals'
         IF (DARWIN) WRITE (LUPRI,'(10X,A)')
     &       '- 1-electron Darwin integrals'
         IF (SPNORB) WRITE (LUPRI,'(10X,A)')
     &                         '- spatial spin-orbit integrals'
c ach
         IF (no2so) WRITE (LUPRI,'(10X,A)')
     &                        'but no two-electron spin-orbit integrals'
         IF (NMRISS) THEN
            IF (FERMI) THEN
               WRITE (LUPRI,'(10X,A)')'- Fermi contact integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (Dirac delta function integrals)'
            END IF
            IF (PSO) THEN
               WRITE (LUPRI,'(10X,A)')
     &              '- paramagnetic spin-orbit integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (nuclear moment - electron orbit coupling)'
            END IF
            IF (SPIDIP) THEN
               WRITE (LUPRI,'(10X,A)')'- spin-dipole integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (electron spin - nuclear moment coupling)'
            END IF
            IF (DSO) THEN
               WRITE (LUPRI,'(10X,A)')
     &              '- diamagnetic spin-orbit integrals'
               WRITE (LUPRI,'(10X,A)')
     &              '  (indirect nuclear dipole - dipole coupling)'
            END IF
            IF (SDFC) THEN
               WRITE (LUPRI,'(10X,A)')
     &             '- spin-dipole + Fermi contact integrals'
               WRITE (LUPRI,'(10X,A)')
     &             '  (electron spin - nuclear magnetic field coupling)'
            END IF
cLig
            IF (RPSO) THEN
               WRITE (LUPRI,'(10X,A)')
     &              '- CTOCD-DZ diamagnetic shielding integrals'
            END IF
cLig
         END IF
         IF (HDO) WRITE (LUPRI,'(10X,A)')
     &       '- antisymmetrized half-derivative overlap integrals'
         IF (S1MAG) WRITE (LUPRI,'(10X,A)')
     &       '- first magnetic derivatives of overlap integrals'
         IF (S1MAGT) WRITE (LUPRI,'(10X,A)')
     &       '- test of first magnetic derivative of overlap integrals'
         IF (S2MAG) WRITE (LUPRI,'(10X,A)')
     &       '- second magnetic derivatives of overlap integrals'
         IF (S2MAGT) WRITE (LUPRI,'(10X,A)')
     &      '- test of second magnetic derivatives of overlap integrals'
         IF (ANGMOM) WRITE (LUPRI,'(10X,A)')
     &      '- electronic angular momentum around the origin'
C copied from linsca abacus/herdrv.F, Bin Gao, December 17, 2009
         IF (S2MBRA) WRITE (LUPRI,'(10X,A)')
     &       '- second magnetic derivatives of overlap integrals, '//
     &       'bra part'
         IF (S2MKET) WRITE (LUPRI,'(10X,A)')
     &       '- second magnetic derivatives of overlap integrals, '//
     &       'ket part'
         IF (S2MMIX) WRITE (LUPRI,'(10X,A)')
     &       '- second magnetic derivatives of overlap integrals, '//
     &       'mixed bra and ket part'
         IF (ANGLON) WRITE (LUPRI,'(10X,A)')
     &      '- electronic angular momentum around the nuclei'
         IF (LONMOM) WRITE (LUPRI,'(10X,A)')
     &      '- London orbital contribution to angular momentum'
         IF (MAGMOM) WRITE (LUPRI,'(10X,A)')
     &      '- one-electron contribution to magnetic moment'
         IF (MGMOMT) WRITE (LUPRI,'(10X,A)')
     &      '- test of London contribution to angular momentum'
         IF (DSUSNL) WRITE (LUPRI,'(10X,A)')
     &   '- Magnetic susceptibility without London orbital contribution'
         IF (DSUSLL) WRITE (LUPRI,'(10X,A)')
     &      '- Angular London orbital contribution to magnetic susc.'
         IF (DSUSLH) WRITE (LUPRI,'(10X,A)')
     &      '- London orbital contribution to magnetic susceptibility'
         IF (DIASUS) WRITE (LUPRI,'(10X,A)')
     &      '- Magnetic susceptibility integrals'
         IF (DSUTST) WRITE (LUPRI,'(10X,A)')
     &     '- Test of London orbital contr. to magnetic susc. integrals'
         IF (NUCSNL) WRITE (LUPRI,'(10X,A)')
     &     '- Nuclear shieldings without London orbital contribution'
         IF (NUCSLO) WRITE (LUPRI,'(10X,A)')
     &     '- London orbital contribution to nuclear shieldings'
         IF (NUCSHI) WRITE (LUPRI,'(10X,A)')
     &     '- Nuclear shielding tensor integrals'
         IF (NSNLTS) WRITE (LUPRI,'(10X,A)')
     &     '- Test of London orbital contribution to nuclear shieldings'
         IF (ELGDIA) WRITE (LUPRI,'(10X,A)')
     &     '- Diamagnetic one-electron spin-orbit (no-London)'
         IF (ELGDIL) WRITE (LUPRI,'(10X,A)')
     &     '- Diamagnetic one-electron spin-orbit (London)'
         IF (MNF_SO) WRITE (LUPRI,'(10X,A)')
     &     '- Mean field spin-orbit integrals'
         IF (NELFLD) WRITE (LUPRI,'(10X,A)')
     &     '- Electric field at the nuclei'
         IF (NSNLTS) WRITE(LUPRI,'(10X,A)')
     &     '- Test of non-London orbital contr. to nuclear shieldings'
         IF (NSTTST) WRITE (LUPRI,'(10X,A)')
     &     '- Test of nuclear shielding tensor integrals'
         IF (EFGCAR) WRITE (LUPRI,'(10X,A)')
     &            '- Cartesien electric field gradient integrals'
         IF (EFGSPH) WRITE (LUPRI,'(10X,A)')
     &            '- Spherical electric field gradient integrals'
         IF (S1MAGL) WRITE (LUPRI,'(10X,A)')
     &        '- Bra-differentiated overlap matrix with respect to B'
         IF (S1MAGR) WRITE (LUPRI,'(10X,A)')
     &        '- Ket-differentiated overlap matrix with respect to B'
         IF (HBDO) WRITE (LUPRI,'(10X,A)')
     &        '-Half B-differentiated overlap matrix'
         IF (HDOBR) WRITE (LUPRI,'(10X,A)')
     &        '- Ket-differentiated hdo-integrals with respect to B'
         IF (S1MLT) WRITE (LUPRI,'(10X,A)')
     &        '- Test of bra-diff. overlap matrix with respect to B'
         IF (S1MRT) WRITE (LUPRI,'(10X,A)')
     &        '- Test of ket-diff. overlap matrix with respect to B'
         IF (HDOBRT) WRITE (LUPRI,'(10X,A)')
     &        '- Test of ket-diff. hdo-integrals with respect to B'
         IF (SQHDOL) WRITE (LUPRI,'(10X,A)')
     &      '- Bra differentiated half-derivative overlap matrix'
         IF (SQHDOR) WRITE (LUPRI,'(10X,A)')
     &      '- Ket differentiated half-derivative overlap matrix'
         IF (SQHD2O) WRITE (LUPRI,'(10X,A)')
     &      '- Ket second differentiated half-derivative overlap matrix'
         IF (NUCPOT) WRITE (LUPRI,'(10X,A)')
     &            '- Potential energy at the nuclei'
         IF (NPOTST) WRITE (LUPRI,'(10X,A)')
     &            '- Test of nuclear potential energy'
         IF (MGMO2T) WRITE (LUPRI,'(10X,A)')
     &            '- Test of two-electron part of magnetic moment'
         IF (SUSCGO) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic magnetizability using common gauge origin'
         IF (NSTCGO) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic shielding tensor using common gauge origin'
         IF (EXPIKR) WRITE (LUPRI,'(10X,A)')
     &      '- Cosine and sine integals'
         IF (DPLGRA) WRITE (LUPRI,'(10X,A)')
     &      '- Dipole gradient integrals'
         IF (DIPANH) WRITE (LUPRI,'(10X,A)')
     &        '- Anharmonic corrections to dipole gradients'
         IF (QUAGRA) WRITE (LUPRI,'(10X,A)')
     &      '- Quadrupole gradient integrals'
         IF (OCTGRA) WRITE (LUPRI,'(10X,A)')
     &      '- Octupole gradient integrals'
         IF (ROTSTR) WRITE (LUPRI,'(10X,A)')
     &      '- Rotational strength integrals'
         IF (SOFLD) WRITE (LUPRI,'(10X,A)')
     &      '- Magnetic-field correction to one-electron SO integrals'
         IF (SOMM) WRITE (LUPRI,'(10X,A)')
     &      '- Magnetic-moment correction to one-electron SO integrals'
         IF (POTENE) WRITE (LUPRI,'(10X,A)')
     &      '- Potential energy Hamiltonian integrals'   
         IF (PVPINT) WRITE (LUPRI,'(10X,A)')
     &      '- pVp integrals of the Douglas-Kroll Hamiltonian'
         IF (PXPINT) WRITE (LUPRI,'(10X,A)')
     &      '- small component dipole length integrals for DPT'
         IF (DKTRAN) WRITE (LUPRI,'(10X,A)')
     &      '- The second order Douglas-Kroll-Hess transformation'//
     &      ' will be applied (DKH2)'
         IF (CM1) THEN
            WRITE (LUPRI,'(10X,A)')
     &         '- First order magnetic derivative of electric field'
            WRITE (LUPRI,'(12X,A,A1,A)')
     &         'Electric field applied in ',FIELD1(1:1),'-direction'
         END IF
         IF (CM2) THEN
            WRITE (LUPRI,'(10X,A)')
     &         '- Second order magnetic derivative of electric field'
            WRITE (LUPRI,'(12X,A,A1,A)')
     &         'Electric field applied in ',FIELD2(1:1),'-direction'
         END IF
         IF (QDBINT) THEN
            WRITE (LUPRI,'(10X,A)')
     &         '- First order magnetic derivative of electric '//
     &         'field gradient'
Cajt qdb            WRITE (LUPRI,'(12X,A,A1,A)')
Cajt qdb     &         'Electric field gradient applied in ',FIELD3(1:2),
            WRITE (LUPRI,'(12X,A,A7,A)')
     &         'Electric field gradient applied in ',FIELD3(1:7),
Cajt qdb
     &         '-direction'
         END IF
         IF (EFBDER) WRITE (LUPRI,'(10X,A)')
     &        '- 1.order London orbital correction to electric field at'
     &        //' a position in space'
         IF (EFB2DR) WRITE (LUPRI,'(10X,A)')
     &        '- 2.order London orbital correction to electric field at'
     &        //' a position in space'
         IF (QDBTST) WRITE (LUPRI,'(10X,A)') 
     &           'Magnetic-field derivative integrals of electric '//
     &           'field gradients will be tested'
         IF (DEROVL) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of overlap integrals'
         IF (DERAM) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of angular momentum integrals'
         IF (DERHAM) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of one-electron Hamiltonian '//
     &        'integrals'
         IF (MAGQDP) WRITE (LUPRI,'(10X,A)')
     &        '- Magnetic quadrupole integrals calculated'
         IF (MQDPTS) WRITE (LUPRI,'(10X,A)')
     &        '- Test of magnetic quadrupole integrals calculated'
         IF (PSOKE) WRITE (LUPRI,'(10X,A)')
     &      '- kinetic energy correction to paramagnetic spin-orbit'
         IF (DNSKE) WRITE (LUPRI,'(10X,A)')
     &      '- kinetic energy correction to diamagnetic shielding'
         IF (OZKE) WRITE (LUPRI,'(10X,A)')
     &      '- kinetic energy correction to orbital Zeeman'
         IF (SDKE) WRITE (LUPRI,'(10X,A)')
     &      '- kinetic energy correction to spin-dipole operator'
         IF (FCKE) WRITE (LUPRI,'(10X,A)')
     &      '- kinetic energy correction to Fermi Contact operator'
         IF (DSOKE) WRITE (LUPRI,'(10X,A)')
     &      '- kinetic energy correction to diamagnetic spin-orbit '//
     &        'operator'
         IF (PSOOZ) WRITE (LUPRI,'(10X,A)')
     &      '- orbital Zeeman correction to paramagnetic spin-orbit'//
     &        ' integral'
         IF (DPTPOT) WRITE (LUPRI,'(10X,A,A)')
     &        '- small component potential energy for DPT'
         IF (DPTOVL) WRITE (LUPRI,'(10X,A)')
     &        '- small component overlap integrals for DPT'
         IF (XDDXR3) WRITE (LUPRI,'(10X,A)')
     &        '-d/dB d/dK   < 1/Rk^3 > type of integrals. '
         IF (PROPRI) WRITE (LUPRI,'(/A)')
     &      ' All one-electron property integrals are printed.'
         IF (S1ELE) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     &        ' Type A'
         IF (S1ELB) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     &        ' Type B'
         IF (ONEELD) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of one-electron'//
     &       ' Hamiltonian integrals'
         IF (PVIOLA) WRITE (LUPRI,'(10X,A)')
     &       '- parity violation electroweak interaction'
cLig added message to display for RPSO and RAGNMO optiion
         IF (RANGMO) WRITE (LUPRI,'(10X,A)')
     &      '- Product between (r - r_go) and l_go for the'//
     &      '  computation of CTOCD-DZ magnetizability in an'//
     &      '  analytical way'
cLig
         IF (LFDIPLN) WRITE (LUPRI,'(10X,A)')
     &     '- Effective dipole operator - only possible with PE library'
      END IF

      WRITE (LUPRI,'(4(/A,3F20.12))')
     &    ' Center of mass  (bohr):', (CMXYZ(I),I=1,3),
     &    ' Operator center (bohr):', (ORIGIN(I),I=1,3),
     &    ' Gauge origin    (bohr):', (GAGORG(I),I=1,3),
     &    ' Dipole origin   (bohr):', (DIPORG(I),I=1,3)

      IF (EXPIKR) WRITE (LUPRI,'(/A,3F20.15)')
     &    ' Wave numbers for exp(ikr):', (EXPKR(I),I=1,3)
      IF (SOTEST) WRITE (LUPRI,'(/A)')
     *    ' Test of spatial spin-orbit integrals.'
      IF (.NOT.HAMILT) WRITE (LUPRI,'(/A)')
     *    ' Ordinary (field-free non-relativistic) Hamiltonian '/
     *     /'integrals not calculated.'
      IF (MGMO2T) WRITE (LUPRI,'(/A,1P,D10.2)')
     &    ' Threshold for testing two-electron part of magnetic moment:'
     &    ,PRTHRS
      IF (NMRISS) THEN
         IF (ALLATM) THEN
            WRITE (LUPRI,'(/2A)')
     &         ' Integrals for all indirect spin-spin',
     &         ' coupling and/or shielding tensors are calculated.'
         ELSE
            WRITE (LUPRI,'(/2A/)')
     &         ' Indirect spin-spin integrals involving the following',
     &         ' nuclei are calculated:'
            WRITE (LUPRI,'(10X,20I3)') (IPATOM(I),I = 1, NPATOM)
         END IF
         IF (DSO) THEN
            WRITE (LUPRI,'(/2A,I3)')
     &        ' Number of integration points for diamagnetic',
     &        ' spin-orbit integrals: ',NPQUAD
            IF (.NOT.TRIANG) WRITE (LUPRI,'(A)')
     &        ' Integrals for symmetry related coupling tensors'
     &           //' JAB and JBA calculated.'
         END IF
      END IF
C
C     **** Process input for various program sections  *****
C
      CALL ER2INI
  301 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         READ (LUCMD, '(A7)') WORD
         GO TO 301
      ELSE IF (PROMPT .EQ. '*') THEN
         IF (WORD(1:2) .EQ. '**') GO TO 399
         CALL UPCASE(WORD)
         DO I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (399,302,303,304,305,306,307,308,309), I
            END IF
         END DO
         WRITE (LUPRI,'(/3A/)') ' Input module ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input modules',LUPRI)
         CALL QUIT('Illegal input directive in '//WORD1)
      ELSE
         WRITE (LUPRI,'(/3A/)') ' Prompter "',PROMPT,
     &      '" illegal or out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input modules',LUPRI)
         CALL QUIT('Error in prompt in '//WORD1//', input section.')
      END IF

!     DATA TABDIR /'*END OF', '*READIN', '*ONEINT', '*TWOINT',  ! TABDIR(1:4)
!    &             '*SUPINT', '*ER2INT', '*SORINT', '*DENFIT',  ! TABDIR(5:8)
!    &             '*QM3INP'/                                   ! TABDIR(9)

C 302   CALL REAINP(WORD,RELCAL,TSTINP)
  302   CALL QUIT('Input ERROR, *READIN has been moved to **DALTON'//
     &     ' section as *MOLBAS')
      GO TO 301
  303   CALL HR1INP(WORD)
      GO TO 301
  304   CALL HR2INP(WORD)
      GO TO 301
  305   CALL HRSINP(WORD)
      GO TO 301
  306   CALL ER2INP(WORD)
      GO TO 301
  307   CALL SRTINP(WORD)
      GO TO 301
        !Density-fitting input keywords
  308   CALL DENSFIT_INP(WORD)
      GO TO 301
  309   CALL QUIT('ERROR, *QM3 has been moved to **DALTON section')
      GO TO 301
C
  399 CONTINUE
C
      CALL HR1INP(WORD)
      CALL HR2INP(WORD)
      CALL HRSINP(WORD)
      CALL ER2INP(WORD)
      CALL SRTINP(WORD)
 1000 CONTINUE
C
      IF (QMMM) CALL READMM(WORK,LWORK)

      CALL SETDCH
      CALL GPCLOSE(LUCMD,'KEEP')
      IF (THRSUP .LT. 0.0D0) THRSUP = THRS
C
      IF (DORLM .AND. .NOT. CAVUSR)
     &   WRITE(LUPRI,'(/A,3F12.6)') ' Cavity center (center of mass):',
     &                               (CAVORG(I), I = 1, 3)
C
      IF (TESTIN) THEN
         CALL QUIT('*** End of input test for Molecule specification')
      END IF
C
      CALL QEXIT('HERINP')
      RETURN
C
      END
C  /* Deck hr1inp */
      SUBROUTINE HR1INP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "orgcom.h"
#include "cbiher.h"
#include "cbihr1.h"
C
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.SOLVEN', '.NOT AL', '.CAVORG',
     &            '.FNMC  ', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      CALL QENTER('HR1INP')
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         GOTO 199
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHR1/
C
      RUNONE = .TRUE.
      IPRONE = IPRDEF
      DORLM  = .FALSE.
      ALLRLM = .TRUE.
      CAVUSR = .FALSE.
      FNMC   = .FALSE.
C
      NEWDEF = WORD .EQ. '*ONEINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ONEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ONEINP.')
    1          CONTINUE
                  RUNONE = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRONE
                  IF (IPRONE .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  DORLM = .TRUE.
                  READ (LUCMD,*) LMAX
               GO TO 100
    4          CONTINUE
                  ALLRLM = .FALSE.
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) (CAVORG(I),I = 1, 3)
                  CAVUSR = .TRUE.
               GO TO 100
    6          CONTINUE
                  FNMC = .TRUE.
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
    9          CONTINUE
               GO TO 100
   10          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ONEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ONEINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) GOTO 199
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for ONEINP:',1)
         IF (.NOT.RUNONE) THEN
            WRITE (LUPRI,'(A)') ' No one-electron integrals calculated.'
         ELSE
            IF (IPRONE .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level in ONEINT:',IPRONE
         END IF
         IF (DORLM) THEN
            WRITE (LUPRI,'(A/A,I2)')
     &         ' One-electron RLM integrals calculated.',
     &         ' Maximum L quantum number: ', LMAX
            IF (ALLRLM) THEN
               WRITE (LUPRI,'(A)') ' All symmetries saved on file.'
            ELSE
               WRITE (LUPRI,'(A)')
     &            ' Only totally symmetric integrals saved on file.'
            END IF
            IF (CAVUSR) WRITE(LUPRI,'(A,3F15.10)')
     &         ' User supplied cavity center',(CAVORG(I),I=1,3)
         END IF
         IF (FNMC) WRITE (LUPRI,'(/,2A)') ' Finite nuclear mass',
     &      ' correction by modifying the kinetic energy integrals.'
      END IF
 199  CALL QEXIT('HR1INP')
      RETURN
      END
C  /* Deck hr2inp */
      SUBROUTINE HR2INP(WORD)
#ifdef ENABLE_EFS
      use efs_interface
#endif

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D1 = 1.0D0,D0 = 0.0D0)
      PARAMETER (NTABLE = 16)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "cbiher.h"
#include "cbihr2.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.PANAS ', '.RETURN', '.SOFOCK',
     &            '.TIME  ', '.ICEDIF', '.IFTHRS', '.THRFAC', '.ERF   ',
     &            '.ERFEXP', '.DOSRIN', '.KSPT2M', '.ERFGAU', '.COMLAM',
     &            '.EFS   '/
      DATA SET/.FALSE./
C
      CALL QENTER('HR2INP')
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         GOTO 199
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHR2/
C
      RUNTWO = .NOT.NOTWO
      IPRTWO = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      RTNTWO = .FALSE.
      TKTIME = .FALSE.
      SOFOCK = .FALSE.
      ICDIFF = 1
      IEDIFF = 1
      IFTHRS = 20
      USRSCR = .FALSE.
      THRFAC(1) = D1
      THRFAC(2) = D1
C
C     Put ERFEXP in /GNRINF/ to awoid TKTIME conflict with twosta.h
C
C      ERFEXP = .FALSE.
C
      NEWDEF = WORD .EQ. '*TWOINT'
      ICHANG = 0
      IPRSUM = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in HR2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in HR2INP.')
    1          CONTINUE
                  RUNTWO = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, *) IPRTWO,
     &                     IPRNTA, IPRNTB, IPRNTC, IPRNTD
                  IPRSUM = IPRNTA + IPRNTB + IPRNTC + IPRNTD
                  IF (IPRTWO .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*,ERR=35) PANAS
                  GOTO 36
 35               PANAS = 0.25D0
                  BACKSPACE LUCMD
 36               CONTINUE 
C     
C     We cannot use new integral code for Panas correction
C
                  SEGBAS = .FALSE.
               GO TO 100
    4          CONTINUE
                  RTNTWO = .TRUE.
               GO TO 100
    5          CONTINUE
C&&&& SOFOCK - construction of Fock matrices in SO-basis
                  SOFOCK = .TRUE.
               GO TO 100
    6          CONTINUE
                  TKTIME = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& ICEDIF Separate screening of Coulomb and exchange contributions
C&&&& in direct SCF
                  READ (LUCMD,*) ICDIFF,IEDIFF
               GO TO 100
    8          CONTINUE
C&&&& IFTHRS Screening threshold in direct construction of Fock matrices
                  READ (LUCMD,*) IFTHRS
                  USRSCR = .TRUE.
               GO TO 100
    9          CONTINUE
C&&& THRFAC: Factors to multiply LL-integral threshold for SL- and SS - integrals
C&&& This option only used in DIRAC
                 READ(LUCMD,*) THRFAC(1),THRFAC(2)
               GO TO 100
   10          CONTINUE
C&&& ERF  : Value of \chi in separation of two-electron operator
C&&&        Vee=Vsr+Vlr ; Vlr=erf(\chi*r_12)/(r_12)
               READ (LUCMD,*,ERR=37) CHIVAL
               GO TO 100
 37            CALL QUIT('Error reading CHIVAL for *TWOINT/.ERF')
   11          CONTINUE
C&&& ERFEXP: Value of \chi in separation of two-electron operator
C&&&         NEW VERSION March 2010 mnp+hjaj
               ERFEXP(0) = .TRUE.
               ERFEXP(2) = .TRUE.
               READ (LUCMD,*,ERR=39) CHIVAL
               GO TO 100
 39            CALL QUIT('Error reading CHIVAL for *TWOINT/.ERFEXP')
   12          CONTINUE
C&&& DOSRIN : Calculate short-range 2-elec. integrals needed for
C&&&          computing the short-range Hartree term (U_sr) in MC-srDFT
               DOSRIN = .TRUE.
               GO TO 100
C&&& KSPT2M
   13          CONTINUE
               READ (LUCMD,*,ERR=38) CHIVAL 
               CHI1ST = .TRUE.
               RUNTWO = .TRUE.
               GO TO 100
 38            CALL QUIT('Error reading CHIVAL for *TWOINT/.MU ')
   14          CONTINUE
C&&& ERFGAU: Value of \chi in separation of two-electron operator
C&&&         Vee=Vsr+Vlr ; Vlr=erf(\chi*r_12)+N*exp(\chi)/(r_12)
               ERFEXP(0) = .TRUE.
               ERFEXP(1) = .TRUE.
               READ (LUCMD,*,ERR=141) CHIVAL
               GO TO 100
  141          CALL QUIT('Error reading CHIVAL for *TWOINT/.ERFGAU')
   15          CONTINUE

C K. Sharkas and J. Toulouse beg
C COMLAM : COMplement LAMbda: value of lambda for linear separation of the electron-electron interaction: Vee= lambda/r_12 + (1-lambda)/r_12
C K. Sharkas, A. Savin, H. J. Aa. Jensen, J. Toulouse, J. Chem. Phys. 137, 044104 (2012)
C For now, range separation must also be activated (with very large mu) when using the linear separation.
C Example of input:
C**INTEGRALS
C*TWOINT
C.DOSRIN
C.COMLAM
C 0.25 (value of lambda)
C.ERF
C 10000000000 (large value of mu)
C
C See in srdft.F for the associated lambda-dependent density functionals.

               COMLAM = .TRUE.
               READ (LUCMD,*,ERR=151) VLAMBDA
               GO TO 100
  151          CALL QUIT('Error reading Lambda for *TWOINT/.COMLAM')
C K. Sharkas and J. Toulouse end
   16          CONTINUE
C&&& EFS
#ifdef ENABLE_EFS
               use_efs = .true.
#else
               call quit('.EFS not available in this version')
#endif
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized in HR2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in HR2INP.')
            END IF
      END IF
  300 CONTINUE
      ICEDIF = ICDIFF + 2*IEDIFF
C     IF (ICHANG .EQ. 0.AND..NOT.(DIRCAL.OR.CHI1ST)) GOTO 199
      IF (NEWDEF) THEN
         CALL HEADER('Set-up from HR2INP:',1)
         IF (.NOT.(RUNTWO.OR.DIRCAL)) THEN
            WRITE (LUPRI,'(A)') ' No two-electron integrals calculated.'
         ELSE
            WRITE (LUPRI,'(A,I5)') ' Print level in TWOINT:',IPRTWO
            IF (IPRSUM .GT. 0) THEN
                 WRITE (LUPRI,'(A,4I3)')
     &                ' Extra output for the following shells:',
     &                 IPRNTA, IPRNTB, IPRNTC, IPRNTD
                IF (RTNTWO) WRITE (LUPRI,'(A)')
     &               ' Program will exit TWOINT after these shells.'
            END IF
            IF (TKTIME) WRITE (LUPRI,'(/,2A)') ' Detailed timing for',
     &         ' integral calculation will be provided.'
            IF (PANAS .NE. 0.0D0) WRITE (LUPRI,'(/,A,F10.5)')
     &           ' Coulomb integrals screened with a factor of',PANAS
            IF (CHIVAL .NE. -1.0D0) THEN
               IF (CHI1ST) THEN
                  DIRCAL = .TRUE.
                  WRITE (LUPRI,'(2A,E15.5)') 
     &              ' Two-electron integrals calculated',
     &              ' for KS-PT2 with mu = ',CHIVAL
                  WRITE (LUPRI,'(A)') 
     &              ' DFT optimization is run .DIRECT' 
               ELSE
                  IF (ERFEXP(1)) THEN
                     WRITE (LUPRI,410)
                  ELSE IF (ERFEXP(2)) THEN
                     WRITE (LUPRI,420)
                  ELSE
                     WRITE (LUPRI,400)
                  ENDIF
                  WRITE (LUPRI,'(A,F10.5)')
     &            '              with the coupling parameter : ',CHIVAL
                  IF (COMLAM) THEN
                     WRITE (LUPRI,'(2A,E15.5)')
     &                 ' Two-electron integrals calculated',
     &                 ' for linear coupling with Lambda = ',VLAMBDA
                  ENDIF
            ENDIF
            END IF
            IF (DIRCAL) THEN
              IF(SOFOCK) THEN
                WRITE(LUPRI,'(A)')
     &            ' * Direct calculation of Fock matrices in SO-basis.' 
              ELSE 
                WRITE(LUPRI,'(A)')
     &            ' * Direct calculation of Fock matrices in AO-basis.' 
              ENDIF

              IF (USRSCR) THEN
                FCKTHR = -IFTHRS
                FCKTHR = 10.0D0**FCKTHR
                WRITE(LUPRI,'(2A,1P,D9.2)') 
     &            ' * Screening threshold in direct Fock ',
     &            'matrix construction: ',FCKTHR
                IF(IFTHRS.GE.16) WRITE(LUPRI,'(4X,A)') 
     &             '---> WARNING : Screening turned off !'
              ELSE
                WRITE(LUPRI,'(A)')
     &    ' * Program controlled screening thresholds used for this.'
              END IF
              IF(ICDIFF.EQ.1) WRITE(LUPRI,'(A)')
     &    ' * Separate density screening of Coulomb integral batches'
              IF(IEDIFF.EQ.1) WRITE(LUPRI,'(A)')
     &    ' * Separate density screening of exchange integral batches'
            ENDIF
            IF (RELCAL .AND. (THRFAC(1).NE.D1.OR.THRFAC(2).NE.D1)) THEN
              WRITE(LUPRI,'(A,2(/3X,A,1P,D10.3))')
     +           ' * Threshold factors for omitting integrals:',
     +           'SL-integrals: ',THRFAC(1),
     +           'SS-integrals: ',THRFAC(2)
            END IF
         END IF
 400  FORMAT(/' DFT-hybrid : Using a Erf type two-elec. operator')
 410  FORMAT(/' DFT-hybrid : Using a Erf+Gau type two-elec. operator')
 420  FORMAT(/' DFT-hybrid : Using a Erf+Exp type two-elec. operator')
      END IF
 199  CALL QEXIT('HR2INP')
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck hrsinp */
      SUBROUTINE HRSINP(WORD)
C
C     Input for FRMSUP routine (FoRM SUPermatrix 2-el integrals)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "gnrinf.h"
#include "cbiher.h"
#include "cbihrs.h"
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NOSYMM', '.ONLY J', '.THRESH',
     &            'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      CALL QENTER('HRSINP')
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         GOTO 199
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHRS/
C
      RUNSUP = SUPMAT
      IPRSUP = IPRDEF
      NOSSUP = .FALSE.
      ONLY_J = .FALSE.
      THRSUP = -1.0D0
C
      NEWDEF = WORD .EQ. '*SUPINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in SUPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in SUPINP.')
    1          CONTINUE
                  RUNSUP = .FALSE.
                  SUPMAT = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRSUP
                  IF (IPRSUP .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  NOSSUP = .TRUE.
               GO TO 100
    4          CONTINUE
                  ONLY_J = .TRUE.
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) THRSUP
                  IF (THR_REDFAC .GT. 0.0D0) THEN
                    ICHANG = ICHANG + 1
                    WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &           ' thresholds multiplied with general factor',THR_REDFAC
                   THRSUP = THRSUP*THR_REDFAC
                  END IF
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
    9          CONTINUE
               GO TO 100
   10          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in SUPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in SUPINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) GOTO 199
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for *SUPINT:',1)
         IF (IPRSUP .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level                        :',IPRSUP
         IF (ONLY_J) THEN
            WRITE(LUPRI,'(A)') ' Only Coulomb (no exchange)'
         ELSE
            WRITE(LUPRI,'(A)') ' P-supermatrix integrals,'//
     &         ' (if requested for DFT: with scaled exchange)'
         END IF
         IF (THRSUP .NE. -1.0D0) WRITE (LUPRI,'(A,D12.2)')
     &         ' Threshold for supermatrix integrals:',THRSUP
         IF (NOSSUP) THEN
            WRITE (LUPRI,'(A)') ' No symmetry used in SUPINT.'
         ELSE
            WRITE (LUPRI,'(A)') ' Symmetry used in SUPINT.'
         END IF
      END IF
 199  CALL QEXIT('HRSINP')
      RETURN
      END
C  /* Deck herint */
      SUBROUTINE HERINT(WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0)
C
#include "gnrinf.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbihr2.h"
#include "cbihrs.h"
#include "cbieri.h"
#include "orgcom.h"
#include "nuclei.h"
#include "huckel.h"
#include "inftap.h"
#include "efield.h"
#include "r12int.h"
#include "drw2el.h"
#include "qm3.h"
#include "pcmlog.h"
C
      DIMENSION   RLMORI(3)
      CHARACTER*8 LABELT(3), LABELS(6)
C
      DIMENSION WORK(LWORK)
#include "memint.h"
C
C     Control routine for calculation of undifferentiated one- and
C     two-electron Hamiltonian integrals and transformation of
C     two-electron integrals to P supermatrix elements.
C
      IF (SKIP) RETURN
      CALL QENTER('HERINT')
C
      IF (IPRDEF .EQ. 1) CALL TITLER('Output from HERINT','*',126)
C
      I2TYP = 0
C
C     **********************************
C     ***** One-Electron Integrals *****
C     **********************************
C
C     Both Hamiltionian and property one-electron integrals.
C
      IF (RUNONE) THEN
         IF (IPRDEF .GT. 1) CALL TITLER('Output from HERONE','*',126)

!        Open standard property integral file AOPROPER if LUPROP .le. 0
!        (if calculating integrals for DKH, then LUPROP is already opened
!         as AODKHINT)
         IF (LUPROP .LE. 0)
     &   CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
C           ***** These integrals are used to modify the Hamiltonian ****
C           ***** thus calculated first of all if requested          ****
C           ***** DPT potential energy  integrals                    ****
C
         IF (ONEPRP .AND. DPTPOT) THEN
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL PR1INT('DPTPO1 ',WORK,LWORK,IDUMMY,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL FLSHFO(LUPRI)
            CALL PR1INT('DPTPO2 ',WORK,LWORK,IDUMMY,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('DPTPOT',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         END IF
C
C        ***********************************************
C        ***** Overlap, kinetic energy and nuclear *****
C        ***** attraction integrals                *****
C        ***********************************************
C
C        Contract any primitive Douglas-Kroll integrals if necessary
C     
         IF (DKTRAN .AND. .NOT. PVPINT) THEN
            CALL DKH_CONT(PROPRI,IPRONE,WORK,LWORK)
         END IF
C
C        *****************************************
C        ***** R(l,m) integrals for solvent  *****
C        *****************************************
C
         IF (DORLM) THEN
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL RLMSET_CBISOL(LMAX)
            CALL RLMNUC(WORK,LWORK,IPRONE)
            CALL TIMER('RLMNUC',TIMSTR,TIMEND)
            CALL DCOPY(3,DIPORG,1,RLMORI,1)
            CALL DCOPY(3,CAVORG,1,DIPORG,1)
            DO 100 IORDER = 0, LMAX
               CALL PR1INT('SOLVENT',WORK,LWORK,IORDER,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
  100       CONTINUE
            CALL DCOPY(3,RLMORI,1,DIPORG,1)
            CALL TIMER('RLMINT',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         ELSE
! some /cbisol/ variables are used for allocation, also when DORLM false /hjaaj
            CALL RLMSET_CBISOL(0)
         END IF
C
C        *****************************************
C        Integrals written on LUONEL -
C        Some previous information has been written in READIN
C        NOTE: must be called after CALL RLMSET_CBISOL
C              because /CBISOL/ is used for allocation in ONEDRV
C        *****************************************
C
         IF (.NOT. ONLYOV) THEN
            IF (HAMILT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL ONEDRV(WORK,LWORK,IPRONE,.FALSE.,0,.FALSE.,.TRUE.,
     &                     .TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,PCM)
               CALL TIMER('ONEDRV',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
         END IF
C
C        *******************************************
C        ***** One-electron property integrals *****
C        *******************************************
C
         IF (DOHUCKEL) THEN
C
C           ***** Overlap integrals for Huckel *****
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL PR1INT('HUCKEL  ',WORK,LWORK,IDUMMY,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('HUCKEL',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         END IF
      END IF ! RUNONE
      IF (RUNONE .OR. ONEPRP) THEN
C
C           ***** Overlap integrals *****
C
! HJAaJ Jan 2011: Always calculate overlap integrals
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL PR1INT('OVERLAP ',WORK,LWORK,IDUMMY,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('OVERLAP',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
C
            IF (ONLYOV) THEN
              CALL QUIT('Overlap integrals calculated')
            END IF
C
C           ***** Modified overlap integrals *****
C                 for population analysis /hjaaj Mar 2006
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL PR1INT('POPOVLP',WORK,LWORK,IDUMMY,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('POPOVLP',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
C
C           R**2 and R**4 integrals /hjaaj Nov 2017
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            IORDER = 2
            CALL PR1INT('R2     ',WORK,LWORK,IORDER,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('R2     ',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
            CALL TIMER('START ',TIMSTR,TIMEND)
            IORDER = 4
            CALL PR1INT('R4     ',WORK,LWORK,IORDER,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('R4     ',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
C
C           ***** Dipole length integrals *****
C
! HJAaJ Jan 2011: Always calculate dipole and quadrupole integrals
!                 and kinetic energy integrals
!                 (For wave function analysis in sirpop.F and RESPONS)
!           IF (DIPLEN) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DIPLEN ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DIPLEN',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
!           END IF
C
C           ***** Quadrupole integrals *****
C
!           IF (QUADRU) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('QUADRUP',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('QUADRUP',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
!           END IF
C
C           ***** Kinetic energy integrals *****
C
!           IF (KINENE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('KINENER',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('KINENE',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
!           END IF

            IF (KINADI) THEN ! reduced mass diagonal correction
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('KINADIA',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('KINADI',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF

C        ------------------------------------------------------

         IF (ONEPRP) THEN

C        ------------------------------------------------------

C
C           ***** Dipole velocity integrals *****
C
            IF (DIPVEL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DIPVEL ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DIPVEL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Effective dipole integrals *****
C           - only when polarizable environment turned on
C             using PE library - otherwise they are just
C             standard dipole integrals

            IF (LFDIPLN) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('LFDIPLN',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('LFDIPL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Quadrupole integrals *****
C                 (Other than 'QUADRUP'
C
            IF (THETA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('THETA  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('THETA ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Second moments integrals *****
C
            IF (SECMOM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SECMOM ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SECMOM',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Third moments integrals *****
C
            IF (THRMOM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('THRMOM ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('THRMOM',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Potential energy integrals *****
C
            IF (POTENE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('POTENER',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('POTENE',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Douglas-Kroll-Hess pVp integrals *****
C
            IF (PVPINT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PVPINT ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PVPINT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Mass-velocity integrals *****
C
            IF (MASSVL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('MASSVEL',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('MASSVL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Darwin integrals *****
C
            IF (DARWIN) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DARWIN ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DARWIN',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Spatial one-electron spin-orbit integrals *****
C
            IF (SPNORB) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SPNORB ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SPNORB',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Cartesian moments integrals *****
C
            IF (CARMOM) THEN
               ISTR = 0
               IF (IORCAR .LT. 0) ISTR = ABS(IORCAR)
               DO 200 IORDER = ISTR, ABS(IORCAR)
                  CALL TIMER('START ',TIMSTR,TIMEND)
                  CALL PR1INT('CARMOM ',WORK,LWORK,IORDER,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  CALL TIMER('CARMOM',TIMSTR,TIMEND)
                  CALL FLSHFO(LUPRI)
  200          CONTINUE
            END IF
C
C           ***** Spherical moments integrals *****
C
            IF (SPHMOM) THEN
               ISTR = 0
               IF (IORSPH .LT. 0) ISTR = ABS(IORSPH)
               DO 300 IORDER = ISTR, ABS(IORSPH)
                  CALL TIMER('START ',TIMSTR,TIMEND)
                  CALL PR1INT('SPHMOM ',WORK,LWORK,IORDER,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  CALL TIMER('SPHMOM',TIMSTR,TIMEND)
                  CALL FLSHFO(LUPRI)
  300          CONTINUE
            END IF
C
C           ***** Fermi contact integrals *****
C
            IF (FERMI) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('FERMI C',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('FERMI ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Paramagnetic spin-orbit integrals *****
C
            IF (PSO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PSO    ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PSO   ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** spin-dipole integrals *****
C
            IF (SPIDIP) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SPIN-DI',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SPIN-D',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** spin-dipole + Fermi contact integrals *****
C
            IF (SDFC) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SDFC   ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SDFC  ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** diamagnetic spin-orbit integrals *****
C
            IF (DSO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DSO    ',WORK,LWORK,IDUMMY,
     &                     NPQUAD,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DSO   ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Half-derivative overlap integrals *****
C
            IF (HDO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('HDO    ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('HDO   ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
            IF (SQHD2O) THEN ! second derivative
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SQHD2OR',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SQHD2OR',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Bra-half derivative overlap integrals ****
C
            IF (SQHDOL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SQHDOL ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SQHDOL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Ket-half derivative overlap integrals *****
C
            IF (SQHDOR) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SQHDOR ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SQHDOR',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** First magnetic derivative of overlap integrals *****
C
            IF (S1MAG) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S1MAG  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S1MAG ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           **** Bra-differentiated overlap integrals (with respect to B) ****
C
            IF (S1MAGL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S1MAGL ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S1MAGL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           **** Ket-differentiated overlap integrals (with respect to B) ****
C
            IF (S1MAGR) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S1MAGR ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S1MAGR',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Half B-differentiated overlap matrix *****
C
            IF (HBDO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('HBDO   ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('HBDO  ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Ket-differentiated hdo-integrals with respect to B *****
C
            IF (HDOBR) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('HDOBR  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('HDOBR ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           *** Test of ket-differentiated hdo-integrals with respect to B ***
C
            IF (HDOBRT) THEN
               LABELT(1) = 'XDIPVEL '
               LABELT(2) = 'YDIPVEL '
               LABELT(3) = 'ZDIPVEL '
               CALL HDBTST(WORK,LWORK,IPRINT,LABELT,NATOM,ORIGIN)
            END IF
C
C           ***** Test of first magnetic derivative of overlap integrals *****
C
            IF (S1MAGT) THEN
               LABELT(1) = 'dS/dBX  '
               LABELT(2) = 'dS/dBY  '
               LABELT(3) = 'dS/dBZ  '
               CALL MG1TST(WORK,LWORK,IPRONE,'OVERLAP ',LABELT,ORIGIN,
     &                     'SM1 ')
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Test of bra-diff. overlap matrix with respect to B *****
C
            IF (S1MLT) THEN
               LABELT(1) = 'd<S|/dBX'
               LABELT(2) = 'd<S|/dBY'
               LABELT(3) = 'd<S|/dBZ'
               CALL MG1TST(WORK,LWORK,IPRONE,'OVERLAP ',LABELT,ORIGIN,
     &                     'S1ML')
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Test of ket-diff. overlap matrix with respect to B *****
C     Need to be rewritten for both operator center and gauge origin
C     K.Ruud, jul.-94
C
            IF (S1MRT) THEN
               WRITE (LUPRI,'(A)') 'This option is currently disabled'//
     &              ' (conflicting origin definitions)'
               CALL QUIT('Unimplemented module S1MRT')
               LABELT(1) = 'd|S>/dBX'
               LABELT(2) = 'd|S>/dBY'
               LABELT(3) = 'd|S>/dBZ'
               CALL MG1TST(WORK,LWORK,IPRONE,'OVERLAP ',LABELT,ORIGIN,
     &                     'S1MR')
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Second magnetic derivatives of overlap integrals *****
C
            IF (S2MAG) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S2MAG  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S2MAG ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Test of second magnetic derivatives of overlap integrals ****
C
            IF (S2MAGT) THEN
               LABELS(1) = 'dS/dB2XX'
               LABELS(2) = 'dS/dB2XY'
               LABELS(3) = 'dS/dB2XZ'
               LABELS(4) = 'dS/dB2YY'
               LABELS(5) = 'dS/dB2YZ'
               LABELS(6) = 'dS/dB2ZZ'
               CALL MGTST2(WORK,LWORK,IPRONE,'OVERLAP ',LABELS)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Second magnetic derivatives of overlap integrals *****
Cdj bra-diff. part 
C
            IF (S2MBRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S2MBRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S2MBRA',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Second magnetic derivatives of overlap integrals *****
Cdj ket-diff. part 
C
            IF (S2MKET) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S2MKET ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S2MKET',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Second magnetic derivatives of overlap integrals *****
Cdj mixed bra-diff. ket-diff. part 
C
            IF (S2MMIX) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('S2MMIX ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('S2MMIX',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Electronic angular momentum around the origin *****
C
            IF (ANGMOM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('ANGMOM ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ANGMOM',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Electronic angular momentum around the nuclei *****
C
            IF (ANGLON) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('ANGLON ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ANGLON',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** London orbital contribution to angular momentum *****
C
            IF (LONMOM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('LONMOM ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('LONMOM',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** One-electron contribution to magnetic moment *****
C
            IF (MAGMOM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('MAGMOM ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('MAGMOM',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Test of London contribution to angular momentum *****
C
            IF (MGMOMT) THEN
               LABELT(1) = 'XLONMOM '
               LABELT(2) = 'YLONMOM '
               LABELT(3) = 'ZLONMOM '
               CALL MG1TST(WORK,LWORK,IPRONE,'ONEHAMIL',LABELT,ORIGIN,
     &                     'LN  ')
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Diamagnetic magnetizability using common gauge origin *****
C
            IF (SUSCGO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DSUSCGO',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SUSCGO',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     **** Diamagnetic shielding tensor using common gauge origin ****
C
            IF (NSTCGO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('NSTCGO ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('NSTCGO',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     **** Diamagnetic shielding tensor using common gauge origin ****
C    
            IF (ELGDIA) THEN
                  CALL TIMER('START ',TIMSTR,TIMEND) 
               CALL PR1INT('ELGDIAN',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ELGDIA',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     **** Diamagnetic shielding tensor using common gauge origin ****
C
            IF (ELGDIL) THEN
                  CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('ELGDIAL',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ELGDIL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ** Diamagnetic contribution to magnetizability, no london contribution **
C
            IF (DSUSNL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DSUSNOL',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DSUSNL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ** Magnetizability, angular part of London orbital contribution **
C
            IF (DSUSLL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DSUSLAN',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DSUSLL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Magnetizability, London orbital contribution *****
C
            IF (DSUSLH) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DSUSLH ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DSUSLH',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Diamagnetic part of magnetizability *****
C
            IF (DIASUS) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DIASUS ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DIASUS',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     **** Test of London orbital contribution to magnetizability ****
C
            IF (DSUTST) THEN
               LABELS(1) = 'XXDSUSLL'
               LABELS(2) = 'XYDSUSLL'
               LABELS(3) = 'XZDSUSLL'
               LABELS(4) = 'YYDSUSLL'
               LABELS(5) = 'YZDSUSLL'
               LABELS(6) = 'ZZDSUSLL'
               LABELT(1) = 'XANGLON '
               LABELT(2) = 'YANGLON '
               LABELT(3) = 'ZANGLON '
               CALL SUSTST(WORK,LWORK,IPRONE,LABELT,LABELS)
            END IF
C
C           ***** Kinetic energy correction to orbital Zeeman *****
C
            IF (OZKE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('OZKE   ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('OZKE  ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Kinetic energy correction to paramagnetic SO integrals *****
C
            IF (PSOKE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PSOKE  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PSOKE ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Kinetic energy correction to diamagnetic shielding *****
C
            IF (DNSKE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DNSKE  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DNSKE ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Kinetic energy correction to spin-dipole integrals *****
C
            IF (SDKE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SDKE   ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SDKE  ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Kinetic energy correction to Fermi contact integrals *****
C
            IF (FCKE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('FCKE   ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('FCKE  ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Kinetic energy correction to diamagnetic spin-orbit *****
C
            IF (DSOKE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DSOKE  ',WORK,LWORK,IDUMMY,
     &                     NPQUAD,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DSOKE ',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ****  Magnetic-field corrected spin-orbit  ****
C
            IF (PSOOZ) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PSOOZ  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PSO-OZ ',TIMSTR,TIMEND)
            END IF
C
C     ***** Potential energy at the individual nuclei *****
C
            IF (NUCPOT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('NUCPOT ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('NUCPOT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Test of potential energy ****
C
            IF (NPOTST) THEN
               CALL NPTST(WORK,LWORK,NATOM)
            END IF
C
C     ***** Electric field at the individual nuclei *****
C
            IF (NELFLD) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('NEFIELD',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('NELFLD',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Electric field gradient at the individual nuclei, cartesian *****
C
            IF (EFGCAR) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('ELFGRDC',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ELFGRD',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Electric field gradient at the individual nuclei, spherical *****
C
            IF (EFGSPH) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('ELFGRDS',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ELFGRD',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Nuclear shielding without London orbital contribution *****
C
            IF (NUCSNL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('NUCSNLO',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('NSTNLO',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** London orbital contribution to nuclear shielding integrals *****
C
            IF (NUCSLO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('NUCSLO ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('NUCSLO',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Nuclear shielding tensor integrals *****
C
            IF (NUCSHI) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('NUCSHI ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('NUCSHI',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C     ***** Test of London orbital contribution to nuclear shieldings *****
C
            IF (NSNLTS) THEN
               CALL NS1TST(WORK,LWORK,IPRINT,DOATOM,NATOM)
            END IF
C
C     ***** Test of non-London orbital contribution to nuclear shieldings *****
C
            IF (NSLTST) THEN
               CALL NSTST2(WORK,LWORK,IPRINT,DOATOM,NATOM)
            END IF
C
C     ***** Test of nuclear shielding tensor integrals *****
C
            IF (NSTTST) THEN
               CALL NSTST3(WORK,LWORK,IPRINT,DOATOM,NATOM)
            END IF
C
C           ***** Cosine and sine integrals *****
C
            IF (EXPIKR) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL DCOPY(3,ORIGIN,1,RLMORI,1)
               CALL DCOPY(3,EXPKR,1,ORIGIN,1)
               CALL PR1INT('EXPIKR ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL DCOPY(3,RLMORI,1,ORIGIN,1)
               CALL TIMER('EXPIKR',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ****  First order magnetic derivative of electric field  ****
C
            IF (CM1) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               IF (FIELD1.EQ.'XYZ-ALL') THEN
                 FIELD1 = 'X-FIELD'
                 CALL PR1INT('CM1    ',WORK,LWORK,IDUMMY,
     &                       IDUMMY,TRIANG,PROPRI,IPRONE)
                 FIELD1 = 'Y-FIELD'
                 CALL PR1INT('CM1    ',WORK,LWORK,IDUMMY,
     &                       IDUMMY,TRIANG,PROPRI,IPRONE)
                 FIELD1 = 'Z-FIELD'
                 CALL PR1INT('CM1    ',WORK,LWORK,IDUMMY,
     &                       IDUMMY,TRIANG,PROPRI,IPRONE)
               ELSE    
                 CALL PR1INT('CM1    ',WORK,LWORK,IDUMMY,
     &                       IDUMMY,TRIANG,PROPRI,IPRONE)
               END IF
               CALL TIMER('CM1   ',TIMSTR,TIMEND)
            END IF
C
C           **** 1.order B derivative of electric field at a point ****
C
            IF (EFBDER) THEN
               CALL TIMER('START',TIMSTR,TIMEND)
               CALL PR1INT('EFIELB1',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('EFBDER',TIMSTR,TIMEND)
            END IF
C
C           **** 2.order B derivative of electric field at a point ****
C
            IF (EFB2DR) THEN
               CALL TIMER('START',TIMSTR,TIMEND)
               CALL PR1INT('EFIELB2',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('EFBDR2',TIMSTR,TIMEND)
            END IF
C
C           ****  Second order magnetic derivative of electric field  ****
C
            IF (CM2) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
C Modified by Bin Gao, December 14, 2009
C copied from linsca
               IF (FIELD2 .EQ. 'XYZ-ALL') THEN
                  FIELD2 = 'X-FIELD'
                  CALL PR1INT('CM2    ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD2 = 'Y-FIELD'
                  CALL PR1INT('CM2    ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD2 = 'Z-FIELD'
                  CALL PR1INT('CM2    ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
               ELSE
                  CALL PR1INT('CM2    ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
               END IF
               CALL TIMER('CM2 ',TIMSTR,TIMEND)
            END IF
C
C           ****  Second order magnetic derivative of electric field  ****
C
            IF (QDBINT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
Cajt qdb               CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
Cajt qdb     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               IF (FIELD3 .EQ. 'XYZ-ALL') THEN
                  FIELD3 = 'XX-FGRD'
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD3 = 'XY-FGRD'
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD3 = 'YY-FGRD'
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD3 = 'XZ-FGRD'
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD3 = 'YZ-FGRD'
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
                  FIELD3 = 'ZZ-FGRD'
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
               ELSE
                  CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                        IDUMMY,TRIANG,PROPRI,IPRONE)
               END IF
Cajt qdb
               CALL TIMER('QDBINT',TIMSTR,TIMEND)
            END IF
C
C           ***** Test of London contribution to angular momentum *****
C
            IF (QDBTST) THEN
               LABELT(1) = FIELD3(1:2)//'-QDB X'
               LABELT(2) = FIELD3(1:2)//'-QDB Y'
               LABELT(3) = FIELD3(1:2)//'-QDB Z'
               LABELS(1) = FIELD3(1:2)//'THETA '
               CALL MG1TST(WORK,LWORK,IPRONE,LABELS,LABELT,ORIGIN,
     &                     '    ')
               CALL FLSHFO(LUPRI)
            END IF
C
C           ****  First geometrical derivative of dipole integrals  ****
C
            IF (DPLGRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DPLGRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DPLGRA',TIMSTR,TIMEND)
            END IF
C
C           ****  Second geometrical derivative of dipole integrals  ****
C
            IF (DIPANH) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DIPANH ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DIPANH',TIMSTR,TIMEND)
            END IF
C
C           ****  First geometrical derivative of quadrupole integrals  ****
C
            IF (QUAGRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('QUAGRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('QUAGRA',TIMSTR,TIMEND)
            END IF
C
C           ****  First geometrical derivative of octupole integrals  ****
C
            IF (OCTGRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('OCTGRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('OCTGRA ',TIMSTR,TIMEND)
            END IF
C
C           **** Magnetic quadrupole integrals ****
C
            IF (MAGQDP) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('MAGQDP ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('MAGQDP',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Test of magnetic quadrupole integrals *****
C
            IF (MQDPTS) THEN
               LABELT(1) = 'XXMAGQDP'
               LABELT(2) = 'YXMAGQDP'
               LABELT(3) = 'ZXMAGQDP'
               LABELS(1) = 'XANGMOM '
               CALL MG1TST(WORK,LWORK,IPRONE,LABELS,LABELT,ORIGIN,
     &                     'MQDP')
               LABELT(1) = 'XYMAGQDP'
               LABELT(2) = 'YYMAGQDP'
               LABELT(3) = 'ZYMAGQDP'
               LABELS(1) = 'YANGMOM '
               CALL MG1TST(WORK,LWORK,IPRONE,LABELS,LABELT,ORIGIN,
     &                     'MQDP')
               LABELT(1) = 'XZMAGQDP'
               LABELT(2) = 'YZMAGQDP'
               LABELT(3) = 'ZZMAGQDP'
               LABELS(1) = 'ZANGMOM '
               CALL MG1TST(WORK,LWORK,IPRONE,LABELS,LABELT,ORIGIN,
     &                     'MQDP')
               CALL FLSHFO(LUPRI)
            END IF
C
C           ****  Rotational strength integrals  ****
C
            IF (ROTSTR) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('ROTSTR ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('ROTSTR ',TIMSTR,TIMEND)
            END IF
C
C           ****  Magnetic-field corrected spin-orbit  ****
C
            IF (SOFLD) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SOFIELD',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SOFIELD',TIMSTR,TIMEND)
            END IF
C
C           ****  ??? integrals  ****
C
            IF (SOMM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('SOMAGMO',WORK,LWORK,IDUMMY,
     &                     NPQUAD,TRIANG,PROPRI,IPRONE)
               CALL TIMER('SOMAGMO',TIMSTR,TIMEND)
            END IF
C
C           ****   First electric derivative of overlap integrals    ****
C           ****   Type A   ****
            IF (S1ELE) THEN
               CALL PR1INT('S1ELE  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
            END IF
C     
C           ****Mean-field spin-orbit integrals                     ****
C           
            IF (MNF_SO) THEN
               CALL AMFIINTEGRALS(PROPRI,IPRONE,WORK,LWORK)
            END IF
C
C           ****   First electric derivative of overlap integrals    ****
C           ****   Type B   ****
            IF (S1ELB) THEN
               CALL PR1INT('S1ELB  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
            END IF
C
C           ****First electric deriv. of one-electron Ham. integrals****
C
            IF (ONEELD) THEN
               CALL PR1INT('ONEELD ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
            END IF
C
C           **** First geometrical derivatives of integrals ***
C
            IF (DEROVL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DEROVLP',WORK,LWORK,IDUMMY,IDUMMY,
     &                     TRIANG,PROPRI,IPRONE)
               CALL TIMER('DEROVL',TIMSTR,TIMEND)
            END IF
C
C           **** First geometrical derivatives of integrals ***
C
            IF (DERAM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DERANGM',WORK,LWORK,IDUMMY,IDUMMY,
     &                     TRIANG,PROPRI,IPRONE)
               CALL TIMER('DERAM ',TIMSTR,TIMEND)
            END IF
C
C           **** First geometr. deriv. of 1el hamiltonian integrals ***
C
            IF (DERHAM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DERHAMI',WORK,LWORK,IDUMMY,IDUMMY,
     &                     TRIANG,PROPRI,IPRONE)
               CALL TIMER('DERHAM',TIMSTR,TIMEND)
            END IF
C
C           ***** DPT overlap integrals ********************************
C
            IF (DPTOVL) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DPTOVL ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DPTOVL',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** Parity violating electroweak interaction *****
C
            IF (PVIOLA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PVIOLA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PVIOLA',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
C
C           ***** DPT pso-lookalikes SQUARED *****
C
            IF (XDDXR3) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('XDDXR3 ',WORK,LWORK,IDUMMY,
     &              IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('XDDXR3',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END  IF
cLig call to PR1INT form RANGMO and RPSO options
C
C           ***** (r-r')l' integrals, for calculation of
C           magetizability  *****
C
            IF (RANGMO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('RANGMO ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('RANGMO ',TIMSTR,TIMEND) 
               CALL FLSHFO(LUPRI)  
            END IF
C
C           ***** (r-r')l'/|r-R_I|**3 integrals, for shieldings *****
C
            IF (RPSO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('RPSO   ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE) 
               CALL TIMER('RPSO  ',TIMSTR,TIMEND)  
               CALL FLSHFO(LUPRI) 
            END IF
C
C           ***** pXp integrals, for DPT correction to dipole *****
C
            IF (PXPINT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PXPINT ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PXPINT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('PRPINT ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('PRPINT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
cLig
         END IF ! (ONEPRP)
      END IF  ! (RUNONE .OR. ONEPRP)
      RUNONE = .TRUE. ! enable one-electron integrals for next iteration if e.g. geometry optimization
C
C     Determine number of different integrals in R12 method (WK/UniKA/19-11-2002).
      NOPP12 = 0
      IF (R12EIN) THEN
         NOPP12 = NOPP12 + 1
      ELSE
         IF (V12INT) NOPP12 = NOPP12 + 1
         IF (R12INT) NOPP12 = NOPP12 + 1
         IF (U12INT) NOPP12 = NOPP12 + 1
         IF (U21INT) NOPP12 = NOPP12 + 1
      END IF
      IF (NOPP12 .EQ. 0)
     &CALL QUIT('Wrong number of integral types in HERINT')
C
C     Determine pointers to various integrals (WK/UniKA/26-11-2002).
      IADV12 = 0
      IADR12 = 0
      IADU12 = 0
      IADU21 = 0
      IF (R12EIN) THEN
         IADV12 = 1
      ELSE
         IF (V12INT) IADV12 = 1
         IF (R12INT) IADR12 = 1
         IF (U12INT) IADU12 = 1
         IF (U21INT) IADU21 = 1
      END IF
      IADR12 = IADV12 + IADR12
      IADU12 = IADR12 + IADU12
      IADU21 = IADU12 + IADU21
C
C     **********************************
C     ***** Two-Electron Integrals *****
C     **********************************
C
C     Cauchy-Scwartz integrals for screening 
C
      IF(DIRCAL) THEN
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL GABGEN(0,0,0,0,0,.FALSE.,IPRDEF,WORK,LWORK)
C       CALL GABGEN(IJOB,ITYPE,IGTYP,MAXDIF,JATOM,NOCONT,WORK,LWORK)
        CALL TIMER('GABGEN',TIMSTR,TIMEND)
      ENDIF
C
C     Integrals on LUINTA
C
      IF (RUNTWO) THEN
         I2TYP = 0
C
C        ********************************************
C        ***** Two-electron repulsion integrals *****
C        ********************************************
C
         IF (HAMILT) THEN
            IF (RUNERI) THEN
               IF (IPRDEF .GT. 1) CALL TITLER('Output from ERI 2INT',
     &            '*',126)
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL MEMGET('REAL',KCCFBT,MXCONT*MXPRIM,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KINDXB,8*MXSHEL*MXCONT,WORK,KFREE,
     &                     LFREE)
               CALL ER2INT(WORK(KCCFBT),WORK(KINDXB),WORK(KFREE),LFREE)
               CALL MEMREL('HERINT.ER2INT',WORK,KWORK,KWORK,KFREE,LFREE)
               CALL TIMER('ERI 2INT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            ELSE
               IF (IPRDEF .GT. 1) CALL TITLER('Output from TWOINT',
     &            '*',126)
               NPAO = MXSHEL*MXAOVC
               CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)
               CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                     WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),0,
     &                     .FALSE.,IPRDEF)
               CALL MEMREL('HERINT.PAOVEC',WORK,KWORK,KJORBS,KFREE,
     &                     LFREE)
               CALL TIMER('START ',TIMSTR,TIMEND)
               IRNTYP = 0
               CALL TWOINT(WORK(KFREE),LFREE,DUMMY,DUMMY,DUMMY,1,IDUMMY,
     &                     IDUMMY,DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,
     &                     .TRUE.,.TRUE.,.FALSE.,TKTIME,
     &                     IPRTWO,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                     RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &                     WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &                     IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                     .FALSE.,.false.)
               CALL TIMER('TWOINT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
               CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
               IF (CHI1ST) CHIVAL = -1.0D0
            END IF
         END IF
C
C        ******************************************************
C        ***** Short-range two-electron integrals (MCDFT) *****
C        ******************************************************
C
         IF(DOSRIN) THEN
            IF (PANAS .EQ. D0 .AND. CHIVAL .EQ. -1.d0) 
     &      CALL QUIT('HERINT ERROR : Short-range integrals '//
     &                'specified with the unmodified 2e-operator!')
            SRINTS = .TRUE.
            NPAO = MXSHEL*MXAOVC
            CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)
            CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                  WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),0,
     &                  .FALSE.,IPRDEF)
            CALL MEMREL('HERINT.PAOVEC',WORK,KWORK,KJORBS,KFREE,
     &                  LFREE)
            CALL TIMER('START ',TIMSTR,TIMEND)
            IRNTYP = 0
            CALL TWOINT(WORK(KFREE),LFREE,DUMMY,DUMMY,DUMMY,1,IDUMMY,
     &                  IDUMMY,DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,
     &                  .TRUE.,.TRUE.,.FALSE.,TKTIME,
     &                  IPRTWO,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                  RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &                  WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                  .FALSE.,.false.)
            CALL TIMER('TWOINT short-range',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            SRINTS = .FALSE.
         ENDIF
C
C        *****************************************************
C        ***** Spatial two-electron spin-orbit integrals *****
C        *****************************************************
C
         IF (SPNORB.and..not.no2so) THEN
            NPAO = MXSHEL*MXAOVC
            CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                  WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),0,
     &                  .FALSE.,IPRDEF)
            CALL MEMREL('HERINT.PAOVEC',WORK,KWORK,KJORBS,KFREE,LFREE)
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL GETTIM(SO1CPU,SO1WAL)
            IRNTYP = - 2
            CALL TWOINT(WORK(KFREE),LFREE,DUMMY,DUMMY,DUMMY,1,IDUMMY,
     &                  IDUMMY,
     &                  DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,.TRUE.,.TRUE.,
     &                  .FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &                  IPRNTD,RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &                  WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                  .FALSE.,.false.)
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            CALL GETTIM(SO2CPU,SO2WAL)
            SOCPU=SO2CPU-SO1CPU
            SOWAL=SO2WAL-SO1WAL
            CALL TIMER('TWOINT 2-el spin-orbit',TIMSTR,TIMEND)
            WRITE(LUPRI,'(2(/A),2(/A,F7.0,A))')
     &      '   Two-electron spin-orbit integrals',
     &      '   =================================',
     &      '   Spin-orbit 2-electron CPU  time ',SOCPU,' seconds',
     &      '   Spin-orbit 2-electron wall time ',SOWAL,' seconds'
            CALL FLSHFO(LUPRI)
         END IF
C
C        *****************************************
C        ***** Two-electron Darwin integrals *****
C        *****************************************
C
         IF (DAR2EL) THEN
            NPAO = MXSHEL*MXAOVC
            CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                  WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),0,
     &                  .FALSE.,IPRDEF)
            CALL MEMREL('HERINT.PAOVEC',WORK,KWORK,KJORBS,KFREE,LFREE)
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL GETTIM(DA1CPU,DA1WAL)
            IRNTYP = 0
            DO2DAR = .TRUE.
            IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
            CALL GPOPEN(LUINTA,'AO2ELDAR','UNKNOWN',' ',' ',
     &                  IDUMMY,.FALSE.)
            CALL TWOINT(WORK(KFREE),LFREE,DUMMY,DUMMY,DUMMY,1,IDUMMY,
     &                  IDUMMY,
     &                  DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,.TRUE.,.TRUE.,
     &                  .FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &                  IPRNTD,RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &                  WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                  .FALSE.,.false.)
            DO2DAR = .FALSE.
            IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            CALL GETTIM(DA2CPU,DA2WAL)
            DACPU=DA2CPU-DA1CPU
            DAWAL=DA2WAL-DA1WAL
            CALL TIMER('TWOINT 2-el Darwin',TIMSTR,TIMEND)
            WRITE(LUPRI,'(2(/A),2(/A,F5.0,A))')
     &      '   Two-electron Darwin integrals',
     &      '   =================================',
     &      '   Spin-orbit Darwin CPU  time ',DACPU,' seconds',
     &      '   Spin-orbit Darwin wall time ',DAWAL,' seconds'
            CALL FLSHFO(LUPRI)
         END IF
C
C        **********************************************
C        ***** Two-electron orbit-orbit integrals *****
C        **********************************************
C
         IF (ORBORB) THEN
            NPAO = MXSHEL*MXAOVC
            CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)     
            CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                  WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),0,
     &                  .FALSE.,IPRDEF)
            CALL MEMREL('HERINT.PAOVEC',WORK,KWORK,KJORBS,KFREE,LFREE)
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL GETTIM(DA1CPU,DA1WAL)
            IRNTYP = 0
            BPH2OO = .TRUE.
            IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
            CALL GPOPEN(LUINTA,'2-ORBORB','UNKNOWN',' ',' ',
     &                  IDUMMY,.FALSE.)
            CALL TWOINT(WORK(KFREE),LFREE,DUMMY,DUMMY,DUMMY,1,IDUMMY,
     &                  IDUMMY,
     &                  DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,.TRUE.,.TRUE.,
     &                  .FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &                  IPRNTD,RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &                  WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                  .FALSE.,.false.)
            BPH2OO = .FALSE.
            IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            CALL GETTIM(DA2CPU,DA2WAL)
            DACPU=DA2CPU-DA1CPU
            DAWAL=DA2WAL-DA1WAL
            CALL TIMER('TWOINT',TIMSTR,TIMEND)
            WRITE(LUPRI,'(2(/A),2(/A,F5.0,A))')
     &      '   Two-electron Darwin integrals',
     &      '   =================================',
     &      '   Orbit-orbit CPU  time ',DACPU,' seconds',
     &      '   Orbit-orbit wall time ',DAWAL,' seconds'
            CALL FLSHFO(LUPRI)
         END IF
C
C        *************************************
C        ***** Test spin-orbit integrals *****
C        *************************************
C
         IF (SOTEST) THEN
            CALL SOCHK2(WORK,LWORK,IPRDEF)
            CALL FLSHFO(LUPRI)
         END IF
C
C        ********************************************************
C        ***** Test of two-electron part of magnetic moment *****
C        ********************************************************
C
         IF (MGMO2T) THEN
            CALL MM2TST(WORK,LWORK,PRTHRS,IPRINT)
         END IF
      END IF
C
      CALL QEXIT('HERINT')
      RETURN
      END
C
      SUBROUTINE MAKE_AOTWOINT(WORK,LWORK)
C
C 11-jan-2007 Hans Joergen Aa. Jensen
C
C Calculate two-electron integrals
C
#include "implicit.h"
      DIMENSION WORK(LWORK)
#include "priunit.h"
C
C Used from common blocks:
C  GNRINF: DIRCAL, DOSRIN
C  CBIHER: ONEPRP
C  CBIHR1: RUNONE
C  CBIHR2: RUNTWO
C  INFTAP: LUINTA, LUSRINT
#include "mxcent.h"
#include "gnrinf.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbihr2.h"
#include "inftap.h"
C
      LOGICAL AOEXIST, DIRCAL_SAV, RUNONE_SAV, ONEPRP_SAV, RUNTWO_SAV !, DOSRIN_SAV
C
      CALL QENTER('MAKE_AOTWOINT')
      CALL GPINQ('AOTWOINT','EXIST',AOEXIST)
      IF (.NOT. AOEXIST) THEN
         RUNONE_SAV = RUNONE
         ONEPRP_SAV = ONEPRP
         RUNTWO_SAV = RUNTWO
!        DOSRIN_SAV = DOSRIN
         DIRCAL_SAV = DIRCAL
         RUNONE = .FALSE.
         ONEPRP = .FALSE.
         RUNTWO = .TRUE.
!        DOSRIN = DOSRIN .AND. .NOT.DIRCAL
         DIRCAL = .FALSE.
C       ... hjaaj: AOSR2INT is only needed if not DIRCAL (not direct)
C       /Jun 09 -- not true any more - hjaaj April 2010 (parallel
C       activates DIRCAL !!!)
C        Open LUINTA/LUSRINT for new integrals
         WRITE (LUPRI,'(/A)') '---> (Re)generating AOTWOINT'
         CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         IF (DOSRIN) THEN
            WRITE (LUPRI,'(A)') '---> and (re)generating AOSR2INT'
            CALL GPOPEN(LUSRINT,'AOSR2INT','UNKNOWN',' ',' ',
     &                  IDUMMY,.FALSE.)
         END IF
         IF (LUPROP .LT. 0)
     &   CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ',' ',IDUMMY,.FALSE.)
         CALL HERINT(WORK,LWORK)
         RUNONE = RUNONE_SAV
         ONEPRP = ONEPRP_SAV
         RUNTWO = RUNTWO_SAV
!        DOSRIN = DOSRIN_SAV
         DIRCAL = DIRCAL_SAV
      END IF

      CALL QEXIT('MAKE_AOTWOINT')
      RETURN
      END
#endif /* ifndef PRG_DIRAC */
C  /* Deck setdch */
      SUBROUTINE SETDCH
C
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "dorps.h"
C
      DO 100 IREP = 0, MAXREP
         DOREPS(IREP) = .TRUE.
  100 CONTINUE
      RETURN
      END
C  /* Deck srtinp */
      SUBROUTINE SRTINP(WORD)
#include "implicit.h"
#include "priunit.h"
C
C inftra : USE_INTSORT
C
#include "cbisor.h"
#include "inftra.h"
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'xxxxxxx', 'xxxxxx ', '.THRQ  ', '.PRINT ', '.IO PRI',
     &            '.INTSYM', '.DELAO ', '.KEEP  ', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
C
C     Initialize /INTSRT/
C
      THRQ2   = 1.0D-15
      ISPRINT = 0
      ISPRFIO = 0
      ISNTSYM = 1
      DELAO  = .FALSE.
      DO I=1,8
         ISKEEP(I) = 0
      END DO
C
      NEWDEF = WORD .EQ. '*SORINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in *SORINT.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in *SORINT.')
    1          CONTINUE
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
                  READ (LUCMD, *) THRQ2
                  IF (THRQ2.NE.1D-15) ICHANG=ICHANG+1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD, *) ISPRINT
                  IF (ISPRINT.NE.0) ICHANG=ICHANG+1
               GO TO 100
    5          CONTINUE
                  READ (LUCMD, *) ISPRFIO
                  IF (ISPRFIO.NE.0) ICHANG=ICHANG+1
               GO TO 100
    6          CONTINUE
                  READ (LUCMD, *) ISNTSYM
                  IF (ISNTSYM .NE. 1) ICHANG = ICHANG + 1
    7          CONTINUE
                  DELAO = .TRUE.
                  ICHANG=ICHANG+1
               GO TO 100
    8          CONTINUE
                  READ (LUCMD,*) (ISKEEP(I), I=1,8)
               GO TO 100
    9          CONTINUE
               GO TO 100
   10          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in *SORINT.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *SORINT.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) RETURN
      IF (USE_INTSORT .AND. NEWDEF) THEN
         CALL HEADER('Changes of defaults for SORTAO:',1)
         IF (THRQ2 .NE. 1.0D-15) WRITE (LUPRI,'(A,D12.2)')
     &        ' Threshold for sorted integrals:', THRQ2
         IF (ISPRINT .NE. 0) WRITE(LUPRI,'(A,I5)')
     &        ' Print level in SORTAO',ISPRINT
         IF (ISPRFIO .NE. 0) WRITE(LUPRI,'(A,I5)')
     &        ' Print level in fast I/O routines',ISPRFIO
         IF (ISNTSYM.NE.1) WRITE(LUPRI,'(A,I5)')
     &        ' Integral symmetry ',ISNTSYM
         IF (DELAO) WRITE(LUPRI,'(A)')
     &        ' Delete basic file of two-electron integrals'
         IF (ISUM(8,ISKEEP,1) .NE.0) WRITE (LUPRI,'(A,8I2)')
     &        ' Integrals to be kept in symmetries with "0":',
     &        (ISKEEP(I),I=1,8)
      END IF
      RETURN
      END
C/* Deck herini */
      SUBROUTINE HERINI
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "gnrinf.h"
#include "inftra.h"
#include "cbirea.h"
#include "cbisol.h"
#include "cbiher.h"      
#include "drw2el.h"
#include "orgcom.h"
#include "elweak.h"
#include "r12int.h"
#include "infinp.h"
#include "qm3.h"
C
Cholesky
#include "ccdeco.h"
Cholesky
      LOGICAL FIRST_CALL
      SAVE    FIRST_CALL
      DATA    FIRST_CALL /.TRUE./
C
      IF (.NOT. FIRST_CALL) RETURN
      FIRST_CALL = .FALSE.
C
C     Initialize /CBIHER/
C
      USE_INTSORT = .FALSE. ! in inftra.h
      IF(RELCAL) THEN
        HAMILT = .FALSE.
        ONEPRP = .FALSE.
        DIPVEL = .FALSE.
        DIRAC_HER = .TRUE.
        SUPMAT = .FALSE.
      ELSE
        HAMILT = .TRUE.
        ONEPRP = .FALSE.
        DIPVEL = .FALSE.
        DIRAC_HER = .FALSE.
        SUPMAT = .NOT. DIRCAL .AND. .NOT. USE_INTSORT
      ENDIF
      SKIP   = .FALSE.
      TSTINP = .FALSE.
      ONEPRP = .FALSE.
      HAMILT = .TRUE.
      DIPLEN = .FALSE.
      QUADRU = .FALSE.
      SPNORB = .FALSE.
      no2so  = .false.
      SOTEST = .FALSE.
      NOTWO  = DIRCAL
      TWOBYTEPACKING = .FALSE.
      SECMOM = .FALSE.
      CARMOM = .FALSE.
      SPHMOM = .FALSE.
      FERMI  = .FALSE.
      PSO    = .FALSE.
      SPIDIP = .FALSE.
      DSO    = .FALSE.
      SDFC   = .FALSE.
      PROPRI = .FALSE.
      HDO    = .FALSE.
      S1MAG  = .FALSE.
      S2MAG  = .FALSE.
      ANGMOM = .FALSE.
      ANGLON = .FALSE.
      LONMOM = .FALSE.
      MAGMOM = .FALSE.
      S1MAGT = .FALSE.
      MGMOMT = .FALSE.
      KINENE = .FALSE.
      S2MAGT = .FALSE.
      DSUSNL = .FALSE.
      DSUSLL = .FALSE.
      DSUSLH = .FALSE.
      DIASUS = .FALSE.
      DSUTST = .FALSE.
      NUCSNL = .FALSE.
      NUCSLO = .FALSE.
      NUCSHI = .FALSE.
      NSTTST = .FALSE.
      NSLTST = .FALSE.
      NELFLD = .FALSE.
      NSNLTS = .FALSE.
      EFGCAR = .FALSE.
      EFGSPH = .FALSE.
      S1MAGL = .FALSE.
      S1MAGR = .FALSE.
      HDOBR  = .FALSE.
      S1MLT  = .FALSE.
      S1MRT  = .FALSE.
      HDOBRT = .FALSE.
      NUCPOT = .FALSE.
      NPOTST = .FALSE.
      MGMO2T = .FALSE.
      HBDO   = .FALSE.
      SUSCGO = .FALSE.
      NSTCGO = .FALSE.
      MASSVL = .FALSE.
      DARWIN = .FALSE.
      CM1    = .FALSE.
      CM2    = .FALSE.
      QDBINT = .FALSE.
      SQHDOL = .FALSE.
      SQHDOR = .FALSE.
      S1ELE  = .FALSE.
      S1ELB  = .FALSE.
      ONEELD = .FALSE.
      THETA  = .FALSE.
      DPLGRA = .FALSE.
      QUAGRA = .FALSE.
      OCTGRA = .FALSE.
      ROTSTR = .FALSE.
      THRMOM = .FALSE.
      SOFLD  = .FALSE.
      SOMM   = .FALSE.
      DEROVL = .FALSE.
      DERAM  = .FALSE.
      DERHAM = .FALSE.   
      ELGDIA = .FALSE.
      ELGDIL = .FALSE.
      MNF_SO = .FALSE.
      PVPINT = .FALSE.
      POTENE = .FALSE.
      PXPINT = .FALSE.
      NOPICH = .FALSE.
      PRTHRS = 1.0D-10
      NPQUAD = 40
      ALLATM = .TRUE.
      TRIANG = .TRUE.
      EXPIKR = .FALSE.
      DO2DAR = .FALSE.
      DAR2EL = .FALSE.
      AD2DAR = .FALSE.
      ORBORB = .FALSE.
      FINDPT = .FALSE.
      NFIELD = 0
      DPTINT = .FALSE.
      NO2DPT = .FALSE.
      DPTFAC = 0.D0
      DARFAC = 0.D0
      S4CENT = .FALSE.
      DPTOVL = .FALSE.
      DPTPOT = .FALSE.
      XDDXR3 = .FALSE.
      LFDIPLN = .FALSE.
C     Initialize R12 logical variables R12CAL (WK/UniKA/04-11-2002).
      R12DIA = .TRUE.
      R12SVD = .FALSE.
      NOAUXB = .FALSE.
      R12PRP = .FALSE.
      NOTONE = .FALSE.
      NOTTWO = .FALSE.
      NOTTRE = .FALSE.
      IANCC2 = 0
      IAPCC2 = 0
      R12LEV = 0.D0
      R12RST = .FALSE.
      R12CAL = .FALSE.
      R12XXL = .FALSE.
      R12CBS = .FALSE.
      R12OLD = .FALSE.
      R12SQR = .FALSE.
      R12TRA = .FALSE.
      R12INT = .FALSE.
      R12EIN = .FALSE.
      R12ECO = .FALSE.
      R12EOR = .FALSE.
      SLATER = .FALSE.
      R12NOA = .FALSE.
      R12NOB = .FALSE.
      R12NOP = .FALSE.
      R12HYB = .TRUE.
      NORXR  = .FALSE.
      MAGQDP = .FALSE.
      MQDPTS = .FALSE.
      U12INT = .FALSE.
      U21INT = U12INT
      V12INT = .TRUE.
      NOTV12 = .NOT. V12INT
      NOTR12 = .NOT. R12INT
      NOTU12 = .NOT. U12INT
      ANTICO = .FALSE.
      DIRSCF = .FALSE.
      U12DIR = .FALSE.
      R12DIR = .FALSE.
      V12DIR = .FALSE.
      DCCR12 = .FALSE.
      STEPIV = .FALSE.
      CUSP12 = .FALSE.
      DUMR12 = .FALSE.
      ONEAUX = .FALSE.
      INTGAC = 0
      INTGAD = 0
      COMBSS = .FALSE.
      LAUXBS = .FALSE.
      LOOPDP = .FALSE.
      VCLTHR =  0D0
      SVDTHR =  1D-12
      DO I = 1,8
        NRXR12(I) = 0
        LOCFRO(I) = 0
      END DO
cLig setting RANGMO and RPSO false
      RANGMO = .FALSE.
      RPSO   = .FALSE.
cLig
      DOLRINTS=.FALSE. ! jim-dgb DOLRINT init to FALSE
C
      CALL DZERO(ORIGIN,3)
chj may 09: GAGORG and DIPORG now initialized in READIN
chj   CALL DZERO(GAGORG,3)
chj   CALL DZERO(DIPORG,3)
      CALL DZERO(EXPKR ,3)
C
Cholesky
      IF (CHOINT) NOTWO = .TRUE.
Cholesky
C
CTSAUE
C    Initialize nuclear model in NUCLEI
C      1 - Point nucleus
C      2 - Gaussian nucleus
C
      IF(RELCAL) THEN
        GAUNUC = .TRUE.
      ELSE
        GAUNUC = .FALSE.
      ENDIF
C     
C     Initialize /CBISOL/ (10-Dec-92 th+hjaaj)
C     -- not used in Hermit; SOLVNT must be false
C        in order to skip solvent modules in ONEDRV
C        and cavity center in VIBMAS
C
      SOLVNT = .FALSE.
C
      PVFINN = 1.D0
      BGWEIL = .FALSE.
      PVIOLA = .FALSE.
C
C     For .R12AUX the cartesian moments CMxxxxxx up to order two
C     are needed.
      IF (LMULBS) CALL SET_CARMOM_2
C
      RETURN
      END
C/* Deck amfiin */
      SUBROUTINE AMFIINTEGRALS(PROPRI,IPRINT,WORK,LWORK)
C
C     Interface to B.Schimmelpfennig's atomic mean-field integral code
C     Written by K.Ruud in northern Norway in January 2001.
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
C
      LOGICAL BREIT, PROPRI
      DIMENSION WORK(LWORK)
#include "gnrinf.h"
#include "inftap.h"
#include "nuclei.h"
#include "cbirea.h"
C
C     Open Mean-field input file generate in READIN
C
      LUTMP2   = -1
      LUMNFPRP = -1
      CALL GPOPEN(LUTMP2,'AMFI_MNF.INP','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      CALL GPOPEN(LUMNFPRP,'AOPROPER.MNF','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
C
C     For the time being, no Breit-integrals (but soon....)
C
      BREIT = .NOT.DKTRAN
CBS  SO FAR ONLY BREIT-PAULI  !!!   Not no-pair
C
C     We loop over symmetry-independent centers
C
      REWIND (LUMNFPRP)
      DO IATOM = 1, NUCIND
         CALL AMFI(LUTMP2,LUMNFPRP,BREIT,GAUNUC,GNUEXP(IATOM),
     &             WORK,LWORK)
      END DO
C
C     We should now have all the atomic mean-field SO integrals, 
C     now generate them in symmetry-orbital basis
C
      REWIND (LUMNFPRP)
      CALL SYMTRAFO(LUMNFPRP,PROPRI,IPRINT,WORK,LWORK)
C
      CALL GPCLOSE(LUTMP2,'KEEP')
      CALL GPCLOSE(LUMNFPRP,'DELETE')
      RETURN
      END
C/* Deck rstgft */
      subroutine rstgft(alpha,i,ggan,ggas)
C
C R12 calculation requested with fitted "r12-times-Slater" geminal (WK/UniKA/03-24-2005).
C
#include "implicit.h"
      dimension ggax(6,6),ggac(6,6),ggas(6),ggan(6)
      character*5 fn
      data ggax /
     &0.5230D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.3383D+00,0.2652D+01,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.2670D+00,0.1503D+01,0.7928D+01,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.2272D+00,0.1077D+01,0.4270D+01,0.1885D+02,0.0000D+00,0.0000D+00,
     &0.2011D+00,0.8520D+00,0.2941D+01,0.9710D+01,0.3980D+02,0.0000D+00,
     &0.1824D+00,0.7118D+00,0.2252D+01,0.6474D+01,0.1966D+02,0.7792D+02 
     & /
       data ggac /
     &0.6498D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4806D+00,0.3307D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.3880D+00,0.3184D+00,0.1768D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.3259D+00,0.3108D+00,0.1755D+00,0.1102D+00,0.0000D+00,0.0000D+00,
     &0.2804D+00,0.3026D+00,0.1785D+00,0.1100D+00,0.7450D-01,0.0000D+00,
     &0.2454D+00,0.2938D+00,0.1815D+00,0.1128D+00,0.7502D-01,0.5280D-01 
     & /
c-----data ggax /---old fit without weight function exp(-2x)------
c    1 0.1760d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,
c    2 0.1069d0, 0.4323d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,
c    3 0.0801d0, 0.2359d0, 0.9192d0, 0.0000d0, 0.0000d0, 0.0000d0,
c    4 0.0654d0, 0.1644d0, 0.4662d0, 1.7980d0, 0.0000d0, 0.0000d0,
c    5 0.0561d0, 0.1273d0, 0.3080d0, 0.8643d0, 3.3200d0, 0.0000d0,
c    6 0.0495d0, 0.1047d0, 0.2289d0, 0.5475d0, 1.5300d0, 5.8660d0
c    & /
c      data ggac /
c    1 0.2811d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,
c    2 0.1036d0, 0.3999d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,
c    3 0.0454d0, 0.2352d0, 0.3692d0, 0.0000d0, 0.0000d0, 0.0000d0,
c    4 0.0220d0, 0.1459d0, 0.2781d0, 0.3006d0, 0.0000d0, 0.0000d0,
c    5 0.0114d0, 0.0929d0, 0.2126d0, 0.2601d0, 0.2354d0, 0.0000d0,
c    6 0.0062d0, 0.0603d0, 0.1620d0, 0.2259d0, 0.2212d0, 0.1828d0
c----& /
      if (i .ge. 1 .and. i .le. 6) then
       do k=1,i
        ggas(k)=ggax(k,i)*alpha**2
       enddo
       write(fn,'(a4,i1)') 'sin.',i
       LUCFN=-1
       CALL GPOPEN(LUCFN,FN,'NEW',' ','FORMATTED',IDUMMY,.FALSE.)
       write(LUCFN,'(a,i1,a)') 'set output ''rsg',i,'.ps'''
       write(LUCFN,'(a)') 'set term postscript'
       write(LUCFN,'(a,f25.20,a)') 'rsg(x)=x*exp(-x*',alpha,')'
       do k=1,i
        ggan(k)=ggac(k,i)
       enddo
       do k=1,i
        if (k.eq.1) then
         if (k.lt.i) then
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*x*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*x*exp(-x**2*',ggas(k),')'
         endif
        else
         if (k.lt.i) then
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*x*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*x*exp(-x**2*',ggas(k),')'
         endif
        endif
       enddo
       write(LUCFN,'(a,i10,a)') 'set xrange [0:',nint(5d0/alpha),']'
       write(LUCFN,'(a,F10.1,a)') 
     &    'set yrange [0:',nint(5d0/alpha)/10.0d0,']'
       write(LUCFN,'(a)') 'plot fit(x), rsg(x)'
       call gpclose(LUCFN,'KEEP')
      else
       call quit('This STG fit is not possible. Sorry.')
      end if
      end
C/* Deck stgfit */
      subroutine stgfit(alpha,i,ggan,ggas)
C
C R12 calculation requested with fitted Slater-type geminal (WK/UniKA/03-24-2005).
C
#include "implicit.h"
      dimension ggax(6,6),ggac(6,6),ggas(6),ggan(6)
      character*5 fn
      data ggax /
     &0.6853D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4254D+00,0.4520D+01,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.3303D+00,0.2321D+01,0.1628D+02,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.2783D+00,0.1591D+01,0.7637D+01,0.4574D+02,0.0000D+00,0.0000D+00,
     &0.2447D+00,0.1225D+01,0.4924D+01,0.1988D+02,0.1127D+03,0.0000D+00,
     &0.2209D+00,0.1004D+01,0.3622D+01,0.1216D+02,0.4587D+02,0.2544D+03
     &/
       data ggac /
     &0.7354D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.5640D+00,0.3102D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4683D+00,0.3087D+00,0.1529D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4025D+00,0.3090D+00,0.1570D+00,0.8898D-01,0.0000D+00,0.0000D+00,
     &0.3532D+00,0.3072D+00,0.1629D+00,0.9321D-01,0.5619D-01,0.0000D+00,
     &0.3144D+00,0.3037D+00,0.1681D+00,0.9811D-01,0.6024D-01,0.3726D-01
     &/
      if (i .ge. 1 .and. i .le. 6) then
       do k=1,i
        ggas(k)=ggax(k,i)*alpha**2
       enddo
       write(fn,'(a4,i1)') 'gin.',i
       LUCFN=-1
       CALL GPOPEN(LUCFN,FN,'NEW',' ','FORMATTED',IDUMMY,.FALSE.)
       write(LUCFN,'(a,i1,a)') 'set output ''stg',i,'.ps'''
       write(LUCFN,'(a)') 'set term postscript'
       write(LUCFN,'(a,f25.20,a)') 'stg(x)=exp(-x*',alpha,')'
       do k=1,i
        ggan(k)=ggac(k,i)
       enddo
       do k=1,i
        if (k.eq.1) then
         if (k.lt.i) then
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*exp(-x**2*',ggas(k),')'
         endif
        else
         if (k.lt.i) then
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*exp(-x**2*',ggas(k),')'
         endif
        endif
       enddo
       write(LUCFN,'(a,i10,a)') 'set xrange [0:',nint(3d0/alpha),']'
       write(LUCFN,'(a,i10,a)') 'set yrange [0:1]'
       write(LUCFN,'(a)') 'plot fit(x), stg(x)'
       call gpclose(LUCFN,'KEEP')
      else
       call quit('This STG fit is not possible. Sorry.')
      end if
      end
C/* Deck erffit */
      subroutine erffit(alpha,i,ggan,ggas)
C
C R12 calculation requested with fitted ERFC-type geminal (WK/UniKA/03-29-2005).
C
#include "implicit.h"
      dimension ggax(6,6),ggac(6,6),ggas(6),ggan(6)
      character*5 fn
      data ggax /
     &0.1701D+01,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.1370D+01,0.8312D+01,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.1256D+01,0.4702D+01,0.2877D+02,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.1198D+01,0.3495D+01,0.1406D+02,0.8056D+02,0.0000D+00,0.0000D+00,
     &0.1161D+01,0.2888D+01,0.9417D+01,0.3575D+02,0.1991D+03,0.0000D+00,
     &0.1136D+01,0.2521D+01,0.7177D+01,0.2232D+02,0.8204D+02,0.4515D+03 
     &/
       data ggac /
     &0.7658D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.6232D+00,0.2676D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.5466D+00,0.2623D+00,0.1308D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4948D+00,0.2604D+00,0.1327D+00,0.7585D-01,0.0000D+00,0.0000D+00,
     &0.4563D+00,0.2577D+00,0.1363D+00,0.7886D-01,0.4774D-01,0.0000D+00,
     &0.4260D+00,0.2543D+00,0.1393D+00,0.8244D-01,0.5094D-01,0.3157D-01
     &/
      if (i .ge. 1 .and. i .le. 6) then
       do k=1,i
        ggas(k)=ggax(k,i)*alpha**2
       enddo
       write(fn,'(a4,i1)') 'cin.',i
       LUCFN=-1
       CALL GPOPEN(LUCFN,FN,'NEW',' ','FORMATTED',IDUMMY,.FALSE.)
       write(LUCFN,'(a,i1,a)') 'set output ''etg',i,'.ps'''
       write(LUCFN,'(a)') 'set term postscript'
       write(LUCFN,'(a,f25.20,a)') 'etg(x)=1.0-erf(x*',alpha,')'
       do k=1,i
        ggan(k)=ggac(k,i)
       enddo
       do k=1,i
        if (k.eq.1) then
         if (k.lt.i) then
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*exp(-x**2*',ggas(k),')'
         endif
        else
         if (k.lt.i) then
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*exp(-x**2*',ggas(k),')'
         endif
        endif
       enddo
       write(LUCFN,'(a,i10,a)') 'set xrange [0:',nint(3d0/alpha),']'
       write(LUCFN,'(a,i10,a)') 'set yrange [0:1]'
       write(LUCFN,'(a)') 'plot fit(x), etg(x)'
       call gpclose(LUCFN,'KEEP')
      else
       call quit('This STG fit is not possible. Sorry.')
      end if
      end
C/* Deck rrffit */
      subroutine rrffit(alpha,i,ggan,ggas)
C
C  R12 calculation requested with fitted "r12-times-ERFC" geminal (WK/UniKA/03-29-2005).
C
#include "implicit.h"
      dimension ggax(6,6),ggac(6,6),ggas(6),ggan(6)
      character*5 fn
      data ggax /
     &0.1492D+01,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.1266D+01,0.5257D+01,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.1185D+01,0.3351D+01,0.1461D+02,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.1143D+01,0.2643D+01,0.8308D+01,0.3411D+02,0.0000D+00,0.0000D+00,
     &0.1117D+01,0.2269D+01,0.6008D+01,0.1808D+02,0.7171D+02,0.0000D+00,
     &0.1099D+01,0.2037D+01,0.4815D+01,0.1239D+02,0.3603D+02,0.1405D+03 
     &/
       data ggac /
     &0.6939D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.5564D+00,0.2815D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4835D+00,0.2677D+00,0.1493D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.4350D+00,0.2602D+00,0.1460D+00,0.9297D-01,0.0000D+00,0.0000D+00,
     &0.3993D+00,0.2535D+00,0.1468D+00,0.9191D-01,0.6278D-01,0.0000D+00,
     &0.3716D+00,0.2472D+00,0.1478D+00,0.9346D-01,0.6282D-01,0.4442D-01
     &/
      if (i .ge. 1 .and. i .le. 6) then
       do k=1,i
        ggas(k)=ggax(k,i)*alpha**2
       enddo
       write(fn,'(a4,i1)') 'din.',i
       LUCFN=-1
       CALL GPOPEN(LUCFN,FN,'NEW',' ','FORMATTED',IDUMMY,.FALSE.)
       write(LUCFN,'(a,i1,a)') 'set output ''reg',i,'.ps'''
       write(LUCFN,'(a)') 'set term postscript'
       write(LUCFN,'(a,f25.20,a)') 'reg(x)=x-x*erf(x*',alpha,')'
       do k=1,i
        ggan(k)=ggac(k,i)
       enddo
       do k=1,i
        if (k.eq.1) then
         if (k.lt.i) then
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*x*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(a,f25.20,a,f25.20,a)') 
     &    'fit(x)=',ggan(k),'*x*exp(-x**2*',ggas(k),')'
         endif
        else
         if (k.lt.i) then
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*x*exp(-x**2*',ggas(k),') + \\'
         else
           write(LUCFN,'(f25.20,a,f25.20,a)') 
     &    ggan(k),'*x*exp(-x**2*',ggas(k),')'
         endif
        endif
       enddo
       write(LUCFN,'(a,i10,a)') 'set xrange [0:',nint(30/alpha),']'
       write(LUCFN,'(a,F10.1,a)') 
     &    'set yrange [0:',nint(5d0/alpha)/10.0d0,']'
       write(LUCFN,'(a)') 'plot fit(x), reg(x)'
       call gpclose(LUCFN,'KEEP')
      else
       call quit('This STG fit is not possible. Sorry.')
      end if
      end
C -- end of herdrv.F --
