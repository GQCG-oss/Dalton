C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: herdrv.F,v 1.1.1.1 2001-02-08 13:33:19 hjj Exp $
C
#include <single.h>
#ifdef UNDEF
/* Comdeck notes */
950602-vebjorn:
      Added flag HRINPC to ensure that HERMIT input processing is done
      only once.
950505-hjaaj: let .NOTWO and .NOHAM imply .NOSUP
 new option .SUPONLy for just SUPMAT .and. .NOT.HAMILT
    (for generating AOSUPINT file from existing AOTWOINT)
    Kenneth: will you put .SUPONLy into HERMIT manual ?
    Done by Kenneth before jan.-96
#endif
C  /* Deck herctl */
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
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <mxorb.h>
C
#include <cbiher.h>
#include <cbihr1.h>
#include <cbihr2.h>
#include <cbihrs.h>
#include <cbieri.h>
#include <gnrinf.h>
#include <exeinf.h>
#include <inftap.h>
#include <huckel.h>
C
      LOGICAL SET
C
      DIMENSION WORK(LWORK)
C
      CALL GETTIM(TIMHER,WALHER)
C
C     *************************
C     ***** Input Section *****
C     *************************
C
      TIMSTR = SECOND()
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
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
      IF (ONEPRP .OR. ADDSTO) CALL GPOPEN(LUPROP,
     &     'AOPROPER','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF (DORLM) CALL GPOPEN(LUSOL,
     &     'AOSOLINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF (RUNTWO .OR. ERITWO .OR. RUNSUP) CALL GPOPEN (LUINTA,
     &     'AOTWOINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF (RUNSUP) CALL GPOPEN(LUSUPM,'AOSUPINT','NEW',' ',' ',IDUMMY,
     &                        .FALSE.)
      CALL HERINT(WORK,LWORK)
C
C     ******************************
C     ***** End of Calculation *****
C     ******************************
C
      CALL GPCLOSE(LUONEL,'KEEP')
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
C
      IF (RUNTWO .OR. RUNERI) FTRCTL = .TRUE.
C
C     ..transformation control, old AO/MO files cannot be used any more.
C
      IF (LUSUPM .GT. 0) CALL GPCLOSE(LUSUPM,'KEEP')
      IF (LUSOL  .GT. 0) CALL GPCLOSE(LUSOL,'KEEP')
      CALL GETTIM(TEND,WEND)
      TIMHER = TEND - TIMHER
      WALHER = WEND - WALHER
      CALL TIMTXT('>>>> Total CPU  time used in HERMIT:',TIMHER,LUPRI)
      CALL TIMTXT('>>>> Total wall time used in HERMIT:',WALHER,LUPRI)
      RETURN
      END
C  /* Deck herinp */
      SUBROUTINE HERINP(WORK,LWORK)
C
C <<< General Input for HERMIT >>>
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <efield.h>
      PARAMETER (NDIR = 7, NTABLE = 90)
      PARAMETER (D0 = 0.0D0)
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     &          WORD1*7
      DIMENSION WORK(LWORK)
C
#include <gnrinf.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbihrs.h>
#include <cbisol.h>
#include <nuclei.h>
#include <orgcom.h>
#include <huckel.h>
#include <ccom.h>
#include <inftra.h>
#include <inftap.h>
#include <cbisor.h>
#include <ecpinf.h>
C
      DATA TABDIR /'*END OF', '*READIN', '*ONEINT', '*TWOINT',
     &             '*SUPINT', '*ER2INT', '*SORINT'/
      DATA TABLE  /'.PRINT ', '.INPTES', '.NOSUP ', '.SPIN-O',
     &             '.DIPLEN', '.NO HAM', '.SOTEST', '.DIPVEL',
     &             '.QUADRU', '.PHASEO', '.SECMOM', '.SUPONL',
     &             '.CARMOM', '.SPHMOM', '.FC    ', '.PSO   ',
     &             '.SD    ', '.DSO   ', '.POINTS', '.SELECT',
     &             '.QUASUM', '.SD+FC ', '.PROPRI', '.HDO   ',
     &             '.S1MAG ', '.S2MAG ', '.ANGMOM', '.ANGLON',
     &             '.LONMOM', '.MAGMOM', '.S1MAGT', '.MGMOMT',
     &             '.KINENE', '.S2MAGT', '.DSUSNL', '.DSUSLL',
     &             '.DSUSLH', '.DIASUS', '.DSUTST', '.NSTNOL',
     &             '.NSTLON', '.NST   ', '.NSNLTS', '.NSLTST',
     &             '.NELFLD', '.NSTTST', '.EFGCAR', '.EFGSPH',
     &             '.S1MAGL', '.S1MAGR', '.HDOBR ', '.S1MLT ',
     &             '.HDOBRT', '.S1MRT ', '.NUCPOT', '.NPOTST',
     &             '.MGMO2T', '.MGMTHR', '.HBDO  ', '.SUSCGO',
     &             '.NSTCGO', '.EXPIKR', '.MASSVE', '.DARWIN',
     &             '.CM-1  ', '.CM-2  ', '.SQHDOL', '.SQHDOR',
     &             '.NOTWO ', '.GAUGEO', '.DIPORG', '.GFACDI',
     &             '.S1ELE ', '.S1ELB ', '.ONEELD', '.THETA ',
     &             '.NUCMOD', '.SORT I', '.DIPGRA', '.QUAGRA',
     &             '.OCTGRA', '.ROTSTR', '.THIRDM', '.SOFIEL',
     &             '.SOMAGM', '.DEROVL', '.DERHAM', '.ELGDIA',
     &             '.ELGDIL', '.MNF-SO'/
C
C     Check if input has been processed earlier.
C
      IPRDEF = IPRUSR + 1
      IF (HRINPC) THEN
         CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         GOTO 1000
      END IF
      HRINPC = .TRUE.
C
C     Initialize effective core potential parameters
C
      ECP = .FALSE.
      NTYECP = 0
C
C     Initialize /CBIHER/
C
CTSAUE
      RNINTS = NEWTRA
      IF(RELCAL) THEN
        HAMILT = .FALSE.
        ONEPRP = .FALSE.
        DIPVEL = .FALSE.
        DIRAC  = .TRUE.
        SUPMAT = .FALSE.
      ELSE
        HAMILT = .TRUE.
        ONEPRP = .FALSE.
        DIPVEL = .FALSE.
        DIRAC  = .FALSE.
        SUPMAT = .NOT. DIRCAL .AND. .NOT. RNINTS
      ENDIF
      SKIP   = .FALSE.
      TSTINP = .FALSE.
      ONEPRP = .FALSE.
      HAMILT = .TRUE.
      DIPLEN = .FALSE.
      QUADRU = .FALSE.
      SPNORB = .FALSE.
      SOTEST = .FALSE.
      NOTWO  = DIRCAL
      SECMOM = .FALSE.
      CARMOM = .FALSE.
      SPHMOM = .FALSE.
      OCTUPO = .FALSE.
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
      SQHDOL = .FALSE.
      SQHDOR = .FALSE.
      S1ELE  = .FALSE.
      S1ELB  = .FALSE.
      ONEELD = .FALSE.
      GFACDI = .FALSE.
      THETA  = .FALSE.
      DPLGRA = .FALSE.
      QUAGRA = .FALSE.
      OCTGRA = .FALSE.
      ROTSTR = .FALSE.
      THRMOM = .FALSE.
      SOFLD  = .FALSE.
      SOMM   = .FALSE.
      DEROVL = .FALSE.
      DERHAM = .FALSE.   
      ELGDIA = .FALSE.
      ELGDIL = .FALSE.
      MNF_SO = .FALSE.
      PRTHRS = 1.0D-10
      NPQUAD = 40
      ALLATM = .TRUE.
      TRIANG = .TRUE.
      EXPIKR = .FALSE.
      CALL DZERO(ORIGIN,3)
      CALL DZERO(GAGORG,3)
      CALL DZERO(DIPORG,3)
      CALL DZERO(EXPKR ,3)
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
      CALL TITLER('Output from HERMIT input processing','*',118)
C
C     **** Find Hermit input *****
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',ERR=910,END=920) WORD
      IF ((WORD .EQ. '**HERMI') .OR. (WORD .EQ. '*HERMIT') .OR.
     &    (WORD .EQ. '**INTEG')) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         NWARN = NWARN + 1
         WRITE (LUPRI,*) 'WARNING HERINP: error reading Dalton input'
  920 CONTINUE
         CALL REAINI(IPRDEF,RELCAL,TSTINP)
         WORD = '**END O'
         GOTO 1
  930 CONTINUE
      WORD1 = WORD
C
C     ***** Process input for COMMON  /CBIHER/  *****
C
  100 READ (LUCMD, '(A7)') WORD
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124,125,126,127,128,129,130,
     &                131,132,133,134,135,136,137,138,139,140,
     &                141,142,143,144,145,146,147,148,149,150,
     &                151,152,153,154,155,156,157,158,159,160,
     &                161,162,163,164,165,166,167,168,169,170,
     &                171,172,173,174,175,176,177,178,179,180,
     &                181,182,183,184,185,186,187,188,189,190), I
            END IF
  99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in HERINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in HERINP.')
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
            READ (LUCMD, *) (ORIGIN(I),I=1,3)
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
            READ (LUCMD,*) IORCAR
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
               WRITE (LUPRI,'(/A,/A,I3,/A,I3)')
     &             ' Too many atoms selected.',
     &             ' Number of atoms selected:    ',NPATOM,
     &             ' Number of atoms allowed:     ',MXCENT
               CALL QUIT('Error in HERINP')
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
 126     CONTINUE
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
 129     CONTINUE
            LONMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 130     CONTINUE
            MAGMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 131     CONTINUE
            S1MAG  = .TRUE.
            S1MAGT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 132     CONTINUE
            MGMOMT = .TRUE.
            LONMOM = .TRUE.
            HAMILT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 133     CONTINUE
            KINENE = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 134     CONTINUE
            S2MAG  = .TRUE.
            S2MAGT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 135     CONTINUE
            DSUSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 136     CONTINUE
            DSUSLL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 137     CONTINUE
            DSUSLH = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 138     CONTINUE
            DIASUS = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 139     CONTINUE
            DSUTST = .TRUE.
            DSUSLL = .TRUE.
            ANGLON = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 140     CONTINUE
            NUCSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 141     CONTINUE
            NUCSLO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 142     CONTINUE
            NUCSHI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 143     CONTINUE
            NSNLTS = .TRUE.
            NUCSLO = .TRUE.
            PSO    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 144     CONTINUE
            NSLTST = .TRUE.
            NELFLD = .TRUE.
            NUCSNL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 145     CONTINUE
            NELFLD = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 146     CONTINUE
            NSTTST = .TRUE.
            NUCSLO = .TRUE.
            NUCSNL = .TRUE.
            NUCSHI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 147     CONTINUE
            EFGCAR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 148     CONTINUE
            EFGSPH = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 149     CONTINUE
            S1MAGL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 150     CONTINUE
            S1MAGR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 151     CONTINUE
            HDOBR  = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 152     CONTINUE
            S1MLT  = .TRUE.
            S1MAGL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 153     CONTINUE
            HDOBR  = .TRUE.
            HDOBRT = .TRUE.
            DIPVEL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 154     CONTINUE
            S1MRT  = .TRUE.
            S1MAGR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 155     CONTINUE
            NUCPOT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 156     CONTINUE
            NUCPOT = .TRUE.
            HAMILT = .TRUE.
            NPOTST = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 157     CONTINUE
            MGMO2T = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 158     CONTINUE
            READ (LUCMD,*) PRTHRS
            GO TO 100
 159     CONTINUE
            HBDO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 160     CONTINUE
            SUSCGO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 161     CONTINUE
            NSTCGO = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 162     CONTINUE
            EXPIKR = .TRUE.
            ONEPRP = .TRUE.
            READ (LUCMD,*) (EXPKR(I),I=1,3)
            GOTO 100
 163     CONTINUE
            MASSVL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 164     CONTINUE
            DARWIN = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 165     CONTINUE
            READ (LUCMD,'(A7)') FIELD1
            IF (.NOT. ((FIELD1 .EQ. 'X-FIELD')
     &          .OR. (FIELD1 .EQ. 'Y-FIELD')
     &          .OR. (FIELD1 .EQ. 'Z-FIELD'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD1,
     &               '" illegal'
               CALL QUIT('Illegal field directions for CM-1 integrals')
            END IF
            CM1    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 166     CONTINUE
            READ (LUCMD,'(A7)') FIELD2
            IF (.NOT. ((FIELD2 .EQ. 'X-FIELD')
     &          .OR. (FIELD2 .EQ. 'Y-FIELD')
     &          .OR. (FIELD2 .EQ. 'Z-FIELD'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD2,
     &               '" illegal'
               CALL QUIT('Illegal field directions for CM-2 integrals')
            END IF
            CM2    = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 167     CONTINUE
            SQHDOL = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 168     CONTINUE
            SQHDOR = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 169     CONTINUE
            NOTWO = .TRUE.
            SUPMAT = .FALSE.
            GO TO 100
 170     CONTINUE
            READ (LUCMD, *) (GAGORG(I),I=1,3)
            GO TO 100
 171     CONTINUE
            READ (LUCMD, *) (DIPORG(I),I=1,3)
            GO TO 100
 172     CONTINUE
            GFACDI = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 173     CONTINUE
            S1ELE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 174     CONTINUE
            S1ELB  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 175     CONTINUE
            ONEELD = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 176     CONTINUE
            THETA  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 177     CONTINUE
C&&& Nuclear model
            READ (LUCMD, *) INUC
            IF    (INUC.EQ.1) THEN
              GAUNUC = .FALSE.
            ELSEIF(INUC.EQ.2) THEN
              GAUNUC = .TRUE.
            ELSE
              WRITE(LUPRI,'(A,I5)')
     &         '*** ERROR *** Unknown nuclear model:',INUC
              CALL QUIT('*** ERROR *** HERINP: Unknown nuclear model')
            ENDIF
            GO TO 100
 178     CONTINUE
            SUPMAT = .FALSE.
            RNINTS = .TRUE.
            NEWTRA = .TRUE.
            GO TO 100
 179     CONTINUE
            DPLGRA = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 180     CONTINUE
            QUAGRA = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 181     CONTINUE
            OCTGRA = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 182     CONTINUE
            ROTSTR = .TRUE.
            PROPRI = .TRUE.
            GO TO 100                   
 183     CONTINUE
            THRMOM = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 184     CONTINUE
            SOFLD  = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 185     CONTINUE
            SOMM   = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 186     CONTINUE
            DEROVL = .TRUE.
            PROPRI = .TRUE.
            GO TO 100
 187     CONTINUE
            DERHAM = .TRUE.
            PROPRI = .TRUE.
            GO TO 100                            
 188     CONTINUE
            ELGDIA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 189     CONTINUE
            ELGDIL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 190     CONTINUE
            MNF_SO = .TRUE.
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 199
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in HERINP.')
      END IF
 199  CONTINUE
      NMRISS = FERMI .OR. PSO .OR. SPIDIP .OR. DSO .OR. SDFC
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
          CALL QUIT('Error in HERMIT input.')
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
      IF (ONEPRP) THEN
         WRITE (LUPRI,'(/A)')
     &  ' The following one-electron property integrals are calculated:'
                     WRITE (LUPRI,'(10X,A)') '- overlap integrals'
         IF (DIPLEN) WRITE (LUPRI,'(10X,A)') '- dipole length integrals'
         IF (DIPVEL) WRITE (LUPRI,'(10X,A)')
     &                         '- dipole velocity integrals'
         IF (QUADRU) WRITE (LUPRI,'(10X,A)')
     &                         '- quadrupole moment integrals'
         IF (THETA)  WRITE (LUPRI,'(10X,A)')
     &                         '- traceless quadrupole moment integrals'
         IF (SECMOM) WRITE (LUPRI,'(10X,A)')'- second moments integrals'
         IF (THRMOM) WRITE (LUPRI,'(10X,A)')'- third moment integrals'   
         IF (SPNORB) WRITE (LUPRI,'(10X,A)')
     &                                  '- spatial spin-orbit integrals'
         IF (OCTUPO) WRITE (LUPRI,'(10X,A)')
     &                                    '- octupole moment integrals'
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
         IF (KINENE) WRITE (LUPRI,'(10X,A)')
     &      '- electronic kinetic energy'
         IF (MASSVL) WRITE (LUPRI,'(10X,A)')
     &       '- mass velocity integrals'
         IF (DARWIN) WRITE (LUPRI,'(10X,A)')
     &       '- 1-electron Darwin integrals'
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
         END IF
         IF (HDO) WRITE (LUPRI,'(10X,A)')
     &                            '- half-derivative overlap integrals'
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
     &     '- Electric field at the nucleus'
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
     &        '- Test og ket-diff. hdo-integrals with respect to B'
         IF (SQHDOL) WRITE (LUPRI,'(10X,A)')
     &        '- Bra differentiated half-derivative overlap matrix'
         IF (SQHDOR) WRITE (LUPRI,'(10X,A)')
     &        '- Ket differentiated half-derivative overlap matrix'
         IF (NUCPOT) WRITE (LUPRI,'(10X,A)')
     &            '- Potenial energy at the nuclei'
         IF (NPOTST) WRITE (LUPRI,'(10X,A)')
     &            '- Test of nuclear potential energy'
         IF (MGMO2T) WRITE (LUPRI,'(10X,A)')
     &            '- Test of two-electron part of magnetic moment'
         IF (SUSCGO) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic magnetizability using common gauge origin'
         IF (NSTCGO) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic shielding tensor using common gauge origin'
         IF (GFACDI) WRITE (LUPRI,'(10X,A)')
     &      '- Diamagnetic London contribution to rotational g-factors'
         IF (EXPIKR) WRITE (LUPRI,'(10X,A)')
     &      '- Cosine and sine integals'
         IF (DPLGRA) WRITE (LUPRI,'(10X,A)')
     &      '- Dipole gradient integrals'
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
         IF (DEROVL) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of overlap integrals'
         IF (DERHAM) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of one-electron Hamiltonian '//
     &        'integrals'                                                   
         IF (PROPRI) WRITE (LUPRI,'(/A)')
     &      ' All one-electron property integrals are printed.'
         IF (S1ELE) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     $        'Type A'
         IF (S1ELB) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     $        'Type B'
         IF (ONEELD) THEN
            WRITE (LUPRI,'(10X,A)')
     $       '- first electric derivatives of one-electron'
            WRITE (LUPRI,'(12X,A)') 'Hamiltonian integrals'
         ENDIF
      END IF
      IF (DNRM2(3,ORIGIN,1) .GT. D0) WRITE (LUPRI,'(/A,3F20.15)')
     &    ' Operator center:', (ORIGIN(I),I=1,3)
      IF (DNRM2(3,GAGORG,1) .GT. D0) WRITE (LUPRI,'(/A,3F20.15)')
     &    ' Gauge origin   :', (GAGORG(I),I=1,3)
      IF (DNRM2(3,DIPORG,1) .GT. D0) WRITE (LUPRI,'(/A,3F20.15)')
     &    ' Dipole origin  :', (DIPORG(I),I=1,3)
      IF (EXPIKR) WRITE (LUPRI,'(/A,3F20.15)')
     &    ' Wave numbers for exp(ikr):', (EXPKR(I),I=1,3)
      IF (SOTEST) WRITE (LUPRI,'(/A)')
     *    ' Test of spatial spin-orbit integrals.'
      IF (.NOT.HAMILT) WRITE (LUPRI,'(/A)')
     *    ' Ordinary (field-free non-relativistic) Hamiltonian '/
     *     /'integrals not calculated.'
      IF (SUPMAT) WRITE (LUPRI,'(/A/A)') ' In this run '//
     *     ' precalculated two-electron integrals are transformed',
     *     ' to P-supermatrix elements.'
      IF (MGMO2T) THEN
         WRITE (LUPRI,'(/A,D12.6)')
     &    ' Threshold for testing two-electron part of magnetic moment:'
     &    ,PRTHRS
      END IF
      IF (NMRISS) THEN
         IF (ALLATM) THEN
            WRITE (LUPRI,'(/2A)')
     &         ' Integrals for all indirect spin-spin',
     &         ' coupling tensors are calculated.'
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
  200 CONTINUE
      CALL REAINI(IPRDEF,RELCAL,TSTINP)
 201  PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         IF (WORD(1:2) .EQ. '**') GO TO 1
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7), I
            END IF
  210    CONTINUE
         WRITE (LUPRI,'(/3A/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in HERINP.')
      ELSE
         WRITE (LUPRI,'(/3A/)') ' Prompter "',PROMPT,
     *      '" illegal or out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in HERINP, error in prompt.')
      END IF
    2   CALL REAINP(WORD,IPRDEF,RELCAL,TSTINP)
      GO TO 201
    3   CALL HR1INP(WORD)
      GO TO 201
    4   CALL HR2INP(WORD)
      GO TO 201
    5   CALL HRSINP(WORD)
      GO TO 201
    6   CALL ER2INP(WORD)
      GO TO 201
    7   CALL SRTINP(WORD)
      GO TO 201
C
    1 CONTINUE
C
      CALL HR1INP(WORD)
      CALL HR2INP(WORD)
      CALL HRSINP(WORD)
      CALL ER2INP(WORD)
      CALL SRTINP(WORD)
 1000 CONTINUE
C
C     Read geometry, orbital spec., etc. (MOLECULE format).
C     Parameter .TRUE. in call indicates that LUONEL must be written.
C
      CALL READIN(WORK,LWORK,.TRUE.)
      CALL SETDCH
      CALL GPCLOSE(LUCMD,'KEEP')
      IF (TSTINP) CALL GPCLOSE(LUONEL,' ')
      IF (THRSUP .LT. 0.0D0) THRSUP = THRS
      IF (DORLM .AND. .NOT. CAVUSR) THEN
         DO 234 I = 1, MXCENT
            ISOTOP(I) = 1
 234     CONTINUE
         KGEOM = 1
         KMASS = KGEOM + 3*NATOMS
         KNAT  = KMASS + NATOMS
         KLAST = KNAT + (NATOMS + 1)/IRAT
         IF (KLAST .GT. LWORK) CALL STOPIT('HERDRV','CMMASS',KLAST,
     &                                     LWORK)
         CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),ISOTOP,IPRDEF)
         CALL DCOPY(3,CMXYZ,1,CAVORG,1)
         WRITE(LUPRI,'(/A,3F12.6)') ' Cavity center (center of mass):',
     &                               (CAVORG(I), I = 1, 3)
      END IF
C
      RETURN
C
      END
C  /* Deck hr1inp */
      SUBROUTINE HR1INP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include <orgcom.h>
#include <cbiher.h>
#include <cbihr1.h>
C
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.SOLVEN', '.NOT AL', '.CAVORG',
     &            'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
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
      CALL DZERO(CAVORG,3)
C
      NEWDEF = WORD .EQ. '*ONEINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
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
      IF (ICHANG .EQ. 0) RETURN
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
      END IF
      RETURN
      END
C  /* Deck hr2inp */
      SUBROUTINE HR2INP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0,D0 = 0.0D0)
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include <cbiher.h>
#include <cbihr2.h>
#include <gnrinf.h>
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.PANAS ', '.RETURN', '.SOFOCK',
     &            '.TIME  ', '.ICEDIF', '.IFTHRS', '.THRFAC', '.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
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
      NEWDEF = WORD .EQ. '*TWOINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
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
     *            '" not recognized in HR2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in HR2INP.')
    1          CONTINUE
                  RUNTWO = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(5I5)') IPRTWO,
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
C&&&& Screening threshold in direct construction of Fock matrices
                  READ (LUCMD,*) IFTHRS
                  USRSCR = .TRUE.
               GO TO 100
    9          CONTINUE
C&&& THRFAC: Factors to multiply LL-integral threshold for SL- and SS - integrals
C&&& This option only used in DIRAC
                 READ(LUCMD,*) THRFAC(1),THRFAC(2)
               GO TO 100
   10          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in HR2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in HR2INP.')
            END IF
      END IF
  300 CONTINUE
      ICEDIF = ICDIFF + 2*IEDIFF
      IF (ICHANG .EQ. 0.AND.(.NOT.DIRCAL)) RETURN
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for HR2INP:',1)
         IF (.NOT.(RUNTWO.OR.DIRCAL)) THEN
            WRITE (LUPRI,'(A)') ' No two-electron integrals calculated.'
         ELSE
            IF (IPRTWO .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level in TWOINT:',IPRTWO
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
            IF (DIRCAL) THEN
              IF(SOFOCK) THEN
                WRITE(LUPRI,'(/1X,A)')
     &            '* Direct calculation of Fock matrices in SO-basis.' 
              ELSE 
                WRITE(LUPRI,'(/1X,A)')
     &            '* Direct calculation of Fock matrices in AO-basis.' 
              ENDIF
              IF (USRSCR) THEN
                FCKTHR = -IFTHRS
                FCKTHR = 10.0D0**FCKTHR
                WRITE(LUPRI,'(1X,A,A,1P,E8.2)') 
     &            '* Screening threshold in direct Fock ',
     &            'matrix construction: ',FCKTHR
                IF(IFTHRS.GE.16) WRITE(LUPRI,'(4X,A)') 
     &             '---> WARNING : Screening turned off !'
              ELSE
                WRITE(LUPRI,'(1X,A)')
     &        '* Program controlled screening thresholds used for this.'
              END IF
              IF(ICDIFF.EQ.1) WRITE(LUPRI,'(1X,A)')
     &    '* Separate density screening of Coulomb integral batches'
              IF(IEDIFF.EQ.1) WRITE(LUPRI,'(1X,A)')
     &    '* Separate density screening of exchange integral batches'
            ENDIF
            IF(RELCAL) THEN
              WRITE(LUPRI,'(1X,A,2(/3X,A,E9.3))')
     +           '* Threshold factors for integrals:',
     +           'SL-integrals: ',THRFAC(1),
     +           'SS-integrals: ',THRFAC(2)
            ENDIF
         END IF
      END IF
      RETURN
      END
C  /* Deck hrsinp */
      SUBROUTINE HRSINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include <cbiher.h>
#include <cbihrs.h>
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NOSYMM', '.OLDSUP', '.THRESH',
     &            'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
C
C     Initialize /CBIHRS/
C
      RUNSUP = SUPMAT
      IPRSUP = IPRDEF
      NOSSUP = .FALSE.
      OLDSUP = .FALSE.
      THRSUP = -1.0D0
C
      NEWDEF = WORD .EQ. '*SUPINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
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
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRSUP
                  IF (IPRSUP .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  NOSSUP = .TRUE.
               GO TO 100
    4          CONTINUE
                  OLDSUP = .TRUE.
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) THRSUP
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
      IF (ICHANG .EQ. 0) RETURN
      IF (RUNSUP. AND. NEWDEF) THEN
         CALL HEADER('Changes of defaults for SUPINT:',1)
         IF (IPRSUP .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level in SUPINT:',IPRSUP
         IF (THRSUP .NE. -1.0D0) WRITE (LUPRI,'(A,D12.2)')
     &         ' Threshold for supermatrix integrals:', THRSUP
         IF (OLDSUP) THEN
            WRITE (LUPRI,'(A)')
     &         ' Old format for supmatrix integral file.'
            IF (NOSSUP) THEN
               NOSSUP = .FALSE.
               WRITE (LUPRI,'(/A/)' )
     &           ' Note: NOSSUP has been set .FALSE. to conform with'
     &           //' .OLDSUP in input.'
            END IF
         END IF
         IF (NOSSUP) THEN
            WRITE (LUPRI,'(A)') ' No symmetry used in SUPINT.'
         ELSE
            WRITE (LUPRI,'(A)') ' Symmetry used in SUPINT.'
         END IF
      END IF
      RETURN
      END
C  /* Deck herint */
      SUBROUTINE HERINT(WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <dummy.h>
C
#include <gnrinf.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbihr2.h>
#include <cbihrs.h>
#include <cbieri.h>
#include <orgcom.h>
#include <nuclei.h>
#include <huckel.h>
#include <ccom.h>
#include <cbisor.h>
C
C
      DIMENSION   RLMORI(3)
      CHARACTER*8 LABELT(3), LABELS(6)
C
      DIMENSION WORK(LWORK)
#include <memint.h>
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
      IF (RUNONE) THEN
         IF (IPRDEF .GT. 1) CALL TITLER('Output from HERONE','*',126)
C
C        ***********************************************
C        ***** Overlap, kinetic energy and nuclear *****
C        ***** attraction integrals                *****
C        ***********************************************
C
C        Integrals written on LUONEL -
C        Some previous information has been written in READIN
C
         IF (HAMILT) THEN
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL ONEDRV(WORK,LWORK,IPRONE,.FALSE.,0,.FALSE.,.TRUE.,
     &                  .TRUE.,.FALSE.)
            CALL TIMER('ONEDRV',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         END IF
C
C        *****************************************
C        ***** R(l,m) integrals for solvent  *****
C        *****************************************
C
         IF (DORLM) THEN
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL RLMNUC(WORK,LWORK,LMAX,IPRONE)
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
         END IF
C
C        *******************************************
C        ***** One-electron property integrals *****
C        *******************************************
C
         IF (ONEPRP .OR. ADDSTO) THEN
C
C           ***** Overlap integrals *****
C
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL PR1INT('OVERLAP ',WORK,LWORK,IDUMMY,
     &                  IDUMMY,TRIANG,PROPRI,IPRONE)
            CALL TIMER('OVERLAP',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
C
C           ***** Dipole length integrals *****
C
            IF (DIPLEN) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DIPLEN ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DIPLEN',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
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
C           ***** Quadrupole integrals *****
C
            IF (QUADRU) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('QUADRUP',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('QUADRUP',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            END IF
            IF (THETA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('THETA  ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('THETA',TIMSTR,TIMEND)
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
C           ***** Kinetic energy integrals *****
C
            IF (KINENE) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('KINENER',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('KINENE',TIMSTR,TIMEND)
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
C           ***** Octupole integrals *****
C
            IF (OCTUPO) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('OCTUPOL',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('OCTUPO',TIMSTR,TIMEND)
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
C           ***** Diamagnetic contribution to molecular g-factors *****
C
            IF (GFACDI) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('GFACDIA',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('GFACDI',TIMSTR,TIMEND)
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
               CALL PR1INT('CM1    ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('CM1   ',TIMSTR,TIMEND)
            END IF
C
C           ****  Second order magnetic derivative of electric field  ****
C
            IF (CM2) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('CM2    ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('CM2 ',TIMSTR,TIMEND)
            END IF
C
C           ****  First geometricl derivative of dipole integrals  ****
C
            IF (DPLGRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DPLGRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DPLGRA',TIMSTR,TIMEND)
            END IF                                                  
C
C           ****  First geometricl derivative of quadrupole integrals  ****
C
            IF (QUAGRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('QUAGRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('QUAGRA',TIMSTR,TIMEND)
            END IF                                               
C
C           ****  First geometricl derivative of octupole integrals  ****
C
            IF (OCTGRA) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('OCTGRA ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('OCTGRA ',TIMSTR,TIMEND)
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
               CALL AMFIIN(PROPRI,IPRONE,WORK,LWORK)
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
C           **** First geometr. deriv. of 1el hamiltonian integrals ***
C
            IF (DERHAM) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DERHAMI',WORK,LWORK,IDUMMY,IDUMMY,
     &                     TRIANG,PROPRI,IPRONE)
               CALL TIMER('DERHAM',TIMSTR,TIMEND)
            END IF                                        
         END IF
      END IF
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
         IF (IPRDEF .GT. 1) CALL TITLER('Output from TWOINT','*',126)
         I2TYP = 0
C
C        ********************************************
C        ***** Two-electron repulsion integrals *****
C        ********************************************
C
         IF (HAMILT) THEN
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
            IRNTYP = 0
            CALL TWOINT(WORK(KFREE),LFREE,DUMMY,DUMMY,DUMMY,1,IDUMMY,
     &                  IDUMMY,
     &                  DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,.TRUE.,.TRUE.,
     &                  .FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &                  IPRNTD,RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &                  WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                  .FALSE.)
            CALL TIMER('TWOINT',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
         END IF
C
C        *****************************************************
C        ***** Spatial two-electron spin-orbit integrals *****
C        *****************************************************
C
         IF (SPNORB) THEN
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
     &                  .FALSE.)
            CALL GETTIM(SO2CPU,SO2WAL)
            SOCPU=SO2CPU-SO1CPU
            SOWAL=SO2WAL-SO1WAL
            CALL TIMER('TWOINT',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            WRITE(LUPRI,'(2(/A),2(/A,F5.0,A))')
     &      '   Two-electron spin-orbit integrals',
     &      '   =================================',
     &      '   Spin-orbit 2-electron CPU  time ',SOCPU,' seconds',
     &      '   Spin-orbit 2-electron wall time ',SOWAL,' seconds'
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
C     ********************************
C     ***** New integral program *****
C     ********************************
C
      IF (RUNERI) THEN
         IF (IPRDEF .GT. 1) CALL TITLER('Output from ER2INT','*',126)
C
C        ********************************************
C        ***** Two-electron repulsion integrals *****
C        ********************************************
C
         IF (HAMILT) THEN
            CALL TIMER('START ',TIMSTR,TIMEND)
            IRNTYP = 0
            IF (DISTST) THEN
               CALL ER2DIS(WORK,LWORK)
            ELSE
               CALL MEMGET('REAL',KCCFBT,MXCONT*MXPRIM  ,
     &                     WORK,KFREE,LFREE)     
               CALL MEMGET('INTE',KINDXB,8*MXSHEL*MXCONT,
     &                     WORK,KFREE,LFREE)     
               CALL ER2INT(WORK(KCCFBT),WORK(KINDXB),WORK(KFREE),LFREE)
               CALL MEMREL('HERINT.ER2INT',WORK,KWORK,KWORK,
     &                     KFREE,LFREE)
            END IF
            CALL TIMER('ERIDRV',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         END IF
      END IF
C
C     **********************************
C     ***** P Supermatrix Elements *****
C     **********************************
C
C     Integrals on LUSUPM
C
      IF (RUNSUP) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL FRMSUP(WORK,LWORK,NOSSUP,OLDSUP,THRSUP,IPRSUP)
         CALL TIMER('FRMSUP',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
      END IF
C
C     **********************************
C     ***** Integral sorting       *****
C     **********************************
C
C     Integrals on LU?
C
      IF (RNINTS) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL SORTAO(WORK,LWORK)
         CALL TIMER('SORTAO',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
      END IF
C
      CALL QEXIT('HERINT')
      RETURN
      END
C  /* Deck setdch */
      SUBROUTINE SETDCH
C
#include <implicit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
C
#include <symmet.h>
#include <dorps.h>
C
      DO 100 IREP = 0, MAXREP
         DOREPS(IREP) = .TRUE.
  100 CONTINUE
      RETURN
      END
C  /* Deck srtinp */
      SUBROUTINE SRTINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <cbisor.h>
      PARAMETER (NTABLE = 10)
C
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.SORTST', '.COMND ', '.THRQ  ', '.PRINT ', '.IO PRI',
     &            '.INTSYM', '.DELAO ', '.KEEP  ', 'XXXXXXX', 'XXXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF' .AND. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
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
      SORTST = 0.D0 ! not used
      COMND  = ' '
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
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(A40)') COMND
                  IF (COMND.NE.' ') ICHANG=ICHANG+1
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
     *            '" not recognized in SUPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in SUPINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) RETURN
      IF (RNINTS. AND. NEWDEF) THEN
         CALL HEADER('Changes of defaults for SORTAO:',1)
         IF (COMND  .NE. ' ') WRITE (LUPRI,'(A)')
     &        ' COMND ',COMND
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
C/* Deck amfiin */
      SUBROUTINE AMFIIN(PROPRI,IPRINT,WORK,LWORK)
C
C     Interface to B.Schimmelpfennig's atomic mean-field integral code
C     Written by K.Ruud in northern Norway in January 2001.
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
C
      LOGICAL BREIT, PROPRI
      DIMENSION WORK(LWORK)
#include <inftap.h>
#include <nuclei.h>
#include <cbirea.h>
C
C     Open Mean-field input file generate in READIN
C
      LUTMP2 = 0
      MNFPRP = 0
      CALL GPOPEN(LUTMP2,'MNF.INP','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      CALL GPOPEN(MNFPRP,'AOPROPER.MNF','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
C
C     For the time being, no Breit-integrals (but soon....)
C
      BREIT = .TRUE.
CBS  SO FAR ONLY BREIT-PAULI  !!!   Not no-pair
C
C     We loop over symmetry-independent center
C
      REWIND (MNFPRP)
      DO IATOM = 1, NUCIND
CBS      print *, 'lwork'
         CALL AMFI(LUTMP2,MNFPRP,OLDNRM,BREIT,GAUNUC,GNUEXP(IATOM),
     &             WORK,LWORK)
      END DO
C
C     We should now have all the atomic mean-field SO integrals, 
C     now generate them in symmetry-orbital basis
C
      REWIND (MNFPRP)
CBS   print *, 'lwork',lwork
         CALL SYMTRAFO(MNFPRP,LUPROP,PROPRI,IPRINT,WORK,LWORK)
      RETURN
      END


