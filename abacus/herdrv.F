C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
#include <single.h>
C
#ifdef UNDEF
c===========================================================================
c/* Comdeck notes */
c951129-hjaaj: force NEWTRA true if .INTSORT (later changed to .SORT Integrals)
c950602-vebjorn:
c  Added flag HRINPC to ensure that HERMIT input processing is done only once.
c950505-hjaaj: let .NOTWO and .NOHAM imply .NOSUP
c new option .SUPONLy for just SUPMAT .and. .NOT.HAMILT
c    (for generating AOSUPINT file from existing AOTWOINT)
c===========================================================================
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
#include <maxorb.h>
C
#include <cbiher.h>
#include <cbihr1.h>
#include <cbihr2.h>
#include <cbihrs.h>
#include <cbisor.h>
#include <cbieri.h>
#include <gnrinf.h>
#include <exeinf.h>
#include <inftap.h>
#include <huckel.h>
#include <r12int.h>
#include <nuclei.h>
C
      LOGICAL SET, FEXIST
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
      IF (ONEPRP .OR. DOHUCKEL .OR. DKTRAN) CALL GPOPEN(LUPROP,
     &     'AOPROPER','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF (DORLM) CALL GPOPEN(LUSOL,
     &     'AOSOLINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF ((RUNTWO .OR. RUNERI) .AND. .NOT. DCCR12) THEN
C        Delete old MO integrals and AOSUP integrals
         IF (RNINTS) THEN
            IF (LUORDA .LE. 0) CALL DAOPEN(LUORDA,'AOORDINT')
                               CALL DARMOV(LUORDA)
            IF (LUMINT .LE. 0) CALL DAOPEN(LUMINT,'MOTWOINT')
                               CALL DARMOV(LUMINT)
         ELSE
            CALL GPINQ('AOSRTINT.DA','EXIST',FEXIST)
            IF (FEXIST) THEN
               CALL SYSTEM('rm -f AOSRTINT.DA*')
               LUINTM = -1
            END IF
            CALL GPINQ('MOTWOINT','EXIST',FEXIST)
            IF (FEXIST) THEN
               CALL GPOPEN(LUINTM,'MOTWOINT',
     &              'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
               CALL GPCLOSE(LUINTM,'DELETE')
            END IF
         END IF
         CALL GPOPEN(LUSUPM,'AOSUPINT',
     &           'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL GPCLOSE(LUSUPM,'DELETE')
C        Open LUINTA for new integrals
         CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      END IF
C
      CALL HERINT(WORK,LWORK)
C
C     **********************************
C     ***** P Supermatrix Elements *****
C     **********************************
C
C     Integrals on AOSUPINT
C     HJAaJ Oct 2003: CALL FRMSUP has been moved to SIRIUS and
C     is only done when relevant ...
C
C     **********************************
C     ***** Integral sorting       *****
C     **********************************
C
C     Integrals on AOORDINT
C
      IF (RNINTS) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL SORTAO(WORK,LWORK)
         CALL TIMER('SORTAO',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
      END IF
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
C     ..transformation control, old AO/MO files cannot be used any more.
C
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
C     tuh, extensively extended by others 
C
C <<< General Input for HERMIT >>>
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <efield.h>
C     NTABLE has been increased from 102 to 105 (WK/UniKA/04-11-2002).
cLig <> changed NTABLE from 111 to 112
C     NTABLE has been increased 114 (WK/UniKA/09-03-2004).
      PARAMETER (NDIR = 8, NTABLE = 123)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D4 = 4.0D0, D10 = 1.0D1)
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
#include <drw2el.h>
#include <ccom.h>
#include <inftra.h>
#include <inftap.h>
#include <infinp.h>
#include <cbisor.h>
#include <ecpinf.h>
#include <elweak.h>
#include <r12int.h>
#include <codata.h>
#include <denfit.h>
C     Keywords '.R12INT', '.R12EIN' and '.MULBAS' have been added (WK/UniKA/04-11-2002).
cLig <> added RANGMO and RPSO
      DATA TABDIR /'*END OF', '*READIN', '*ONEINT', '*TWOINT',
     &             '*SUPINT', '*ER2INT', '*SORINT', '*DENFIT'/
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
     &             '.NOTWO ', '.GAUGEO', '.DIPORG', '.NO2SO ',
     &             '.S1ELE ', '.S1ELB ', '.ONEELD', '.THETA ',
     &             '.NUCMOD', '.SORT I', '.DIPGRA', '.QUAGRA',
     &             '.OCTGRA', '.ROTSTR', '.THIRDM', '.SOFIEL',
     &             '.SOMAGM', '.DEROVL', '.DERHAM', '.ELGDIA', 
     &             '.ELGDIL', '.MNF-SO', '.DPTOVL', '.DPTPOT', 
     &             '.XDDXR3', '.AD2DAR', '.PVIOLA', '.WEINBG', 
     &             '.FINDPT', '.FNPROP', '.PVP   ', '.1ELPOT',
     &             '.QDBINT', '.QDBTST', '.R12   ', '.R12EXP',
     &             '.AUXBAS', '.NOTV12', '.R12INT', '.U12INT',
     &             '.U21INT', '.DCCR12', '.RANGMO', '.RPSO  ',
     &             '.DPTPXP', '.NOPICH', '.OZ-KE ', '.PSO-KE', 
     &             '.DNS-KE', '.SD-KE ', '.FC-KE ', '.DSO-KE', 
     &             '.PSO-OZ', '.SQHD2O', 'xxxxxxx'/
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
C     *** Initialize ***
C
      CALL HERINI
      !Initialize defaults for density fitting
      CALL DFIINI
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
         CALL ER2INI
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
C              Labels 203-209 have been added (WK/UniKA/04-11-2002).
cLig <> added go to 211 and 212
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
     &                221,222,223), I
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
c ach
            no2so=.true.
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
            ONEPRP = .TRUE.
            GO TO 100
 180     CONTINUE
            QUAGRA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 181     CONTINUE
            OCTGRA = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 182     CONTINUE
            ROTSTR = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100                   
 183     CONTINUE
            THRMOM = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 184     CONTINUE
            SOFLD  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 185     CONTINUE
            SOMM   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 186     CONTINUE
            DEROVL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 187     CONTINUE
            DERHAM = .TRUE.
            ONEPRP = .TRUE.
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
            ONEPRP = .TRUE.
            GO TO 100
 191     CONTINUE
            DPTOVL = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 192     CONTINUE
            DPTPOT = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 193     CONTINUE
            XDDXR3 = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100        
 194     CONTINUE
            AD2DAR = .TRUE.
            READ (LUCMD, *) DARFAC
            GO TO 100                            
 195     CONTINUE
            PVIOLA = .TRUE.
            SPNORB = .TRUE.
            ONEPRP = .TRUE. 
            GO TO 100                            
 196     CONTINUE
            READ(LUCMD,*) BGWEIN
            BGWEIL = .TRUE.
            BGWEIN = D1 - D4*BGWEIN
            GOTO 100   
 197     CONTINUE
            IF (DIRCAL) CALL QUIT('No direct calculation with FINDPT')
            ONEPRP = .TRUE. 
            DPTPOT = .TRUE.
            FINDPT = .TRUE.
            READ (LUCMD, *) DPTFAC
            GOTO 100 
 198     CONTINUE
            READ (LUCMD, *) X10FNP
            PVFINN = D10 ** X10FNP
            GOTO 100
 199     CONTINUE
            PVPINT = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100                            
 200     CONTINUE
            POTENE = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 201     CONTINUE
            READ (LUCMD,'(A7)') FIELD3
            IF (.NOT. ((FIELD3 .EQ. 'XX-FGRD')
     &          .OR.   (FIELD3 .EQ. 'XY-FGRD')
     &          .OR.   (FIELD3 .EQ. 'XZ-FGRD')
     &          .OR.   (FIELD3 .EQ. 'YY-FGRD')
     &          .OR.   (FIELD3 .EQ. 'YZ-FGRD')
     &          .OR.   (FIELD3 .EQ. 'ZZ-FGRD'))) THEN
               WRITE (LUPRI,'(/,3A,/)') ' Field direction "',FIELD1,
     &               '" illegal'
               CALL QUIT('Illegal field direction for QDBINT integrals')
            END IF
            QDBINT = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 202     CONTINUE
            QDBTST = .TRUE.
            THETA  = .TRUE.
            ONEPRP = .TRUE.
            GOTO 100
 203     CONTINUE
C           R12 calculation requested with linear r12 terms (WK/UniKA/04-11-2002).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            GOTO 100
 204     CONTINUE
C           R12 calculation requested with Gaussian-damped r12 integrals (WK/UniKA/04-11-2002).
            R12CAL = .TRUE.
            CARMOM = .TRUE.
            IORCAR = MAX(IORCAR,2)
            ONEPRP = .TRUE.
            R12EOR = .TRUE.
            READ(LUCMD,*) GAMAC
            GAMAD = GAMAC
C           CALL QUIT('The Gaussian-damped R12 method is '//
C    &                'not available in this version. Sorry.')
            GOTO 100
 205     CONTINUE
C           Calculation with multiple basis sets requested (WK/UniKA/04-11-2002).
            LMULBS = .TRUE.
            GOTO 100
 206     CONTINUE
C           No calculation of 1/r12 integrals (WK/UniKA/04-11-2002).
            NOTV12 = .TRUE.
            V12INT = .FALSE.
            GOTO 100
 207     CONTINUE
C           Calculation of integrals of the operator r12 (WK/UniKA/04-11-2002).
            R12INT = .TRUE.
            GOTO 100
 208     CONTINUE
C           Calculation of integrals of the operator U12 (WK/UniKA/04-11-2002).
            U12INT = .TRUE.
            GOTO 100
 209     CONTINUE
C           Calculation of integrals of the operator U21 (WK/UniKA/04-11-2002).
            U21INT = .TRUE.
            GOTO 100
 210     CONTINUE
C           Integrals are computed for DIRCCR12 program (WK/UniKA/25-11-2002).
            DCCR12 = .TRUE.
            SUPMAT = .FALSE.
            GOTO 100
cLig >> added what to do for RANGMO and RPSO options
 211     CONTINUE
            RANGMO = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 212     CONTINUE
            RPSO   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 213     CONTINUE
            PXPINT = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 214     CONTINUE
            NOPICH = .TRUE.
            GO TO 100
cLig <<
 215     CONTINUE
            OZKE   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 216     CONTINUE
            PSOKE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 217     CONTINUE
            DNSKE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100

 218     CONTINUE
            SDKE   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 219     CONTINUE
            FCKE   = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 220     CONTINUE
            DSOKE  = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 221     CONTINUE
            PSOOZ = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
C        .SQH2DO
 222     CONTINUE
            SQHD2O = .TRUE.
            ONEPRP = .TRUE.
            GO TO 100
 223     CONTINUE
            GO TO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 299
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in HERINP.')
      END IF
 299  CONTINUE
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
      IF (R12INT) WRITE (LUPRI,'(A)')
     &    ' Calculation of two-electron integrals over r12.'
      IF (U12INT) WRITE (LUPRI,'(A)')
     &    ' Calculation of two-electron integrals over [T1,r12].'
      IF (U21INT) WRITE (LUPRI,'(A)')
     &    ' Calculation of two-electron integrals over [T2,r12].'
      IF (AD2DAR) WRITE (LUPRI,'(/A/A,F18.14/)')
     &    ' Two-electron Darwin integrals are added to standard',
     &    ' repulsion integrals with perturbation parameter:',DARFAC
      IF (FINDPT) WRITE (LUPRI,'(/A/A,F18.14/A,F18.14,A/)')
     &    ' A direct relativistic perturbation is added to',
     &    ' Hamiltonian and metric with perturbation parameter:',DPTFAC,
     &    ' (perturbation parameter * c^{-2} = ',DPTFAC*ALPHAC**2,')'
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
     &                         '- spatial spin-orbit integrals'
c ach
         IF (no2so) WRITE (LUPRI,'(10X,A)')
     &                        'but no two-electron spin-orbit integrals'
         IF (OCTUPO) WRITE (LUPRI,'(10X,A)')
     &                         '- octupole moment integrals'
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
cLig >>
            IF (RPSO) THEN
               WRITE (LUPRI,'(10X,A)')
     &              '- CTOCD-DZ diamagnetic shielding integrals'
            END IF
cLig <<
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
     &        '- Test of ket-diff. hdo-integrals with respect to B'
         IF (SQHDOL) WRITE (LUPRI,'(10X,A)')
     &      '- Bra differentiated half-derivative overlap matrix'
         IF (SQHDOR) WRITE (LUPRI,'(10X,A)')
     &      '- Ket differentiated half-derivative overlap matrix'
         IF (SQHD2O) WRITE (LUPRI,'(10X,A)')
     &      '- Ket second differentiated half-derivative overlap matrix'
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
         IF (POTENE) WRITE (LUPRI,'(10X,A)')
     &      '- Potential energy Hamiltonian integrals'   
         IF (PVPINT) WRITE (LUPRI,'(10X,A)')
     &      '- pVp integrals of the Douglas-Kroll Hamiltonian'
         IF (PXPINT) WRITE (LUPRI,'(10X,A)')
     &      '- small component dipole lengths integrals for DPT'
         IF (DKTRAN) WRITE (LUPRI,'(10X,A)')
     &      '- The second order Douglas-Kroll transformation'//
     &      ' will be applied'
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
            WRITE (LUPRI,'(12X,A,A1,A)')
     &         'Electric field gradient applied in ',FIELD3(1:2),
     &         '-direction'
         END IF
         IF (QDBTST) WRITE (LUPRI,'(10X,A)') 
     &           'Magnetic-field derivative integrals of electric '//
     &           'field gradients will be tested'
         IF (DEROVL) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of overlap integrals'
         IF (DERHAM) WRITE (LUPRI,'(10X,A)')
     &        '- Geometrical derivatives of one-electron Hamiltonian '//
     &        'integrals'
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
     &        'Type A'
         IF (S1ELB) WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of overlap integrals,'//
     &        'Type B'
         IF (ONEELD) THEN
            WRITE (LUPRI,'(10X,A)')
     &       '- first electric derivatives of one-electron'
            WRITE (LUPRI,'(12X,A)') 'Hamiltonian integrals'
         ENDIF
         IF (PVIOLA) WRITE (LUPRI,'(10X,A)')
     &       '- parity violation electroweak interaction'
cLig >> added message to display for RPSO and RAGNMO optiion
         IF (RANGMO) WRITE (LUPRI,'(10X,A)')
     &      '- Product between (r - r_go) and l_go for the'//
     &      '  computation of CTOCD-DZ magnetizability in an'//
     &      '  analytical way'
cLig <<
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
      IF (MGMO2T) THEN
         WRITE (LUPRI,'(/A,D12.6)')
     &    ' Threshold for testing two-electron part of magnetic moment:'
     &    ,PRTHRS
      END IF
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
  300 CONTINUE
      CALL REAINI(IPRDEF,RELCAL,TSTINP)
      CALL ER2INI
 301  PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 300
      ELSE IF (PROMPT .EQ. '*') THEN
         IF (WORD(1:2) .EQ. '**') GO TO 1
         DO 310 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO (1,2,3,4,5,6,7,8), I
            END IF
  310    CONTINUE
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
      GO TO 301
    3   CALL HR1INP(WORD)
      GO TO 301
    4   CALL HR2INP(WORD)
      GO TO 301
    5   CALL HRSINP(WORD)
      GO TO 301
    6   CALL ER2INP(WORD)
      GO TO 301
    7   CALL SRTINP(WORD)
      GO TO 301
        !Density-fitting input keywords
    8   CALL DFIINP(WORD)
      GO TO 301
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
C
      IF (DORLM .AND. .NOT. CAVUSR)
     &   WRITE(LUPRI,'(/A,3F12.6)') ' Cavity center (center of mass):',
     &                               (CAVORG(I), I = 1, 3)
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
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NOSYMM', 'xxxxxxx', '.THRESH',
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
#include <maxorb.h>
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
#include <efield.h>
#include <r12int.h>
#include <drw2el.h>
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
C     Contract any primitive Douglas-Kroll integrals if necessary
C     
         IF (DKTRAN .AND. .NOT. PVPINT) THEN
            CALL DKCONT(PROPRI,IPRONE,WORK,LWORK)
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
         IF (ONEPRP) THEN
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
C           ***** Douglas-Kroll pVp integrals *****
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
            IF (SQHD2O) THEN
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
C           ****  Second order magnetic derivative of electric field  ****
C
            IF (CM2) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('CM2    ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('CM2 ',TIMSTR,TIMEND)
            END IF
C
C           ****  Second order magnetic derivative of electric field  ****
C
            IF (QDBINT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('QDBINT ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
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
C           ***** DPT potential energy  integrals ***********************
C
            IF (DPTPOT) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL PR1INT('DPTPO1 ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL FLSHFO(LUPRI)
               CALL PR1INT('DPTPO2 ',WORK,LWORK,IDUMMY,
     &                     IDUMMY,TRIANG,PROPRI,IPRONE)
               CALL TIMER('DPTPOT',TIMSTR,TIMEND)
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
cLig >> call to PR1INT form RANGMO and RPSO options
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
cLig <<
         END IF
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
     &                  .TRUE.,.FALSE.,.FALSE.,.FALSE.)
            CALL TIMER('ONEDRV',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
         END IF
      END IF
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
         IF (IPRDEF .GT. 1) CALL TITLER('Output from TWOINT','*',126)
         I2TYP = 0
C
C        ********************************************
C        ***** Two-electron repulsion integrals *****
C        ********************************************
C
         IF (HAMILT) THEN
            IF (RUNERI) THEN
               CALL TIMER('START ',TIMSTR,TIMEND)
               CALL MEMGET('REAL',KCCFBT,MXCONT*MXPRIM,WORK,KFREE,LFREE)
               CALL MEMGET('INTE',KINDXB,8*MXSHEL*MXCONT,WORK,KFREE,
     &                     LFREE)
               CALL ER2INT(WORK(KCCFBT),WORK(KINDXB),WORK(KFREE),LFREE)
               CALL MEMREL('HERINT.ER2INT',WORK,KWORK,KWORK,KFREE,LFREE)
               CALL TIMER('ERIDRV',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
            ELSE
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
     &                     .FALSE.)
               CALL TIMER('TWOINT',TIMSTR,TIMEND)
               CALL FLSHFO(LUPRI)
               CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            END IF
         END IF
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
     &                  .FALSE.)
            CALL MEMREL('HERINT.TWOINT',WORK,KWORK,KWORK,KFREE,LFREE)
            CALL GETTIM(SO2CPU,SO2WAL)
            SOCPU=SO2CPU-SO1CPU
            SOWAL=SO2WAL-SO1WAL
            CALL TIMER('TWOINT',TIMSTR,TIMEND)
            WRITE(LUPRI,'(2(/A),2(/A,F5.0,A))')
     &      '   Two-electron spin-orbit integrals',
     &      '   =================================',
     &      '   Spin-orbit 2-electron CPU  time ',SOCPU,' seconds',
     &      '   Spin-orbit 2-electron wall time ',SOWAL,' seconds'
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
C  /* Deck setdch */
      SUBROUTINE SETDCH
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
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
      DATA TABLE /'xxxxxxx', 'xxxxxx ', '.THRQ  ', '.PRINT ', '.IO PRI',
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
      IF (RNINTS. AND. NEWDEF) THEN
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
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <nuclei.h>
#include <gnrinf.h>
#include <inftra.h>
#include <cbisol.h>
#include <cbisor.h>
#include <cbiher.h>      
#include <ecpinf.h>
#include <drw2el.h>
#include <orgcom.h>
#include <elweak.h>
#include <r12int.h>
#include <infinp.h>
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
C#ifdef VAR_INT64
C        SUPMAT = .FALSE.
C#else
        SUPMAT = .NOT. DIRCAL .AND. .NOT. RNINTS
C#endif
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
      AD2DAR = .FALSE.
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
C     Initialize R12 logical variables R12CAL-LMULBS (WK/UniKA/04-11-2002).
      R12DIA = .TRUE.
      R12SVD = .FALSE.
      NOAUXB = .FALSE.
      R12PRP = .FALSE.
      NOTONE = .FALSE.
      NOTTWO = .FALSE.
      R12LEV = 0.D0
      R12RST = .FALSE.
      R12CAL = .FALSE.
      R12XXL = .FALSE.
      R12OLD = .FALSE.
      R12SQR = .FALSE.
      R12TRA = .FALSE.
      R12INT = .FALSE.
      R12EIN = .FALSE.
      R12ECO = .FALSE.
      R12EOR = .FALSE.
      R12NOA = .FALSE.
      R12NOB = .FALSE.
      R12NOP = .FALSE.
      R12HYB = .TRUE.
      NORXR  = .FALSE.
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
      LMULBS = .FALSE.
      LOOPDP = .FALSE.
      VCLTHR =  0D0
      SVDTHR =  1D-12
cLig >> setting RANGMO and RPSO false
      RANGMO = .FALSE.
      RPSO   = .FALSE.
cLig <<
C
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
      PVFINN = 1.D0
      BGWEIL = .FALSE.
      PVIOLA = .FALSE.
C
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
#include <gnrinf.h>
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
      BREIT = .NOT.DKTRAN
CBS  SO FAR ONLY BREIT-PAULI  !!!   Not no-pair
C
C     We loop over symmetry-independent center
C
      REWIND (MNFPRP)
      DO IATOM = 1, NUCIND
         CALL AMFI(LUTMP2,MNFPRP,OLDNRM,BREIT,GAUNUC,GNUEXP(IATOM),
     &             WORK,LWORK)
      END DO
C
C     We should now have all the atomic mean-field SO integrals, 
C     now generate them in symmetry-orbital basis
C
      REWIND (MNFPRP)
      CALL SYMTRAFO(MNFPRP,LUPROP,PROPRI,IPRINT,WORK,LWORK)
C
      CALL GPCLOSE(LUTMP2,'KEEP')
      CALL GPCLOSE(MNFPRP,'DELETE')
      RETURN
      END
C  /* Deck dkint */
      SUBROUTINE DKINTF(WORK,LWORK)
C
C     Driver for calculating Douglas-Kroll integrals in the primitive
C     basis and transforming them to the appropriate form before 
C     writing them to file
C
C     K.Ruud, September-00
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL HUCKELBKP
C
#include <inftap.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbihr2.h>
#include <cbieri.h>
#include <cbihrs.h>
#include <gnrinf.h>
#include <huckel.h>
C
      CALL GETTIM(TIMHER,WALHER)
C
      CALL HERINI
      HAMILT = .FALSE.
      PVPINT = .TRUE.
      POTENE = .TRUE.
      KINENE = .TRUE.
      PROPRI = .FALSE.
      HUCKELBKP = DOHUCKEL
      DOHUCKEL = .FALSE.
      ONEPRP = .TRUE.
      RUNONE = .TRUE.
      RUNSUP = .FALSE.
      RUNTWO = .FALSE.
C
C     Process basis set, although now as a decontracted set without
C     normalization of the primitive orbitals (UNCONT .true.)
C     (READIN will be called again in HERINP and generate the
C     proper contracted basis and the info necessary to
C     transform Douglas-Kroll integrals to contracted basis)
C      
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
      CALL DKREIN(WORK,LWORK)
C
      TIMSTR = SECOND()
      TIMINP = SECOND() - TIMSTR
      CALL FLSHFO(LUPRI)
      CALL GPOPEN(LUPROP,'DKINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL HERINT(WORK,LWORK)
      CALL GPCLOSE(LUONEL,'DELETE')
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
C
C     We now calculate all the Douglas-Kroll quantities in the primitive
C     basis
C
      CALL DKDRIV(WORK,LWORK,IPRDEF)
C
      DOHUCKEL = HUCKELBKP
      RDINPC = .FALSE.
      RUNSUP = SUPMAT
      PVPINT = .FALSE.
      RUNTWO = .NOT.NOTWO
      HAMILT = .TRUE.
      CALL GETTIM(TEND,WEND)
      TIMHER = TEND - TIMHER
      WALHER = WEND - WALHER
      CALL TIMTXT('>>>> Total CPU  time used in DKINT:',TIMHER,LUPRI)
      CALL TIMTXT('>>>> Total wall time used in DKINT:',WALHER,LUPRI)
      RETURN
      END
C/* Deck dkrein */
      SUBROUTINE DKREIN(WORK,LWORK)
C
C     Read specification of primitive basis (i.e. decontracted)
C     for calculation of Douglas-Kroll integrals.
C
#include <implicit.h>
#include <mxcent.h>
      DIMENSION WORK(LWORK)
      LOGICAL UNCONTBKP
#include <maxaqn.h>
#include <ccom.h>
#include <mxsymm.h>
#include <cbirea.h>
C
C     Define cbirea.h for call of READIN below
C
      THRS = 1.0D-15
C
      LUMLCL = 0
      MAXPRI = MAXPRD
      BIGVC  = .FALSE.
      DIRAC  = .FALSE.
      GENCON = .TRUE.
      BASIS  = .FALSE.
      ATOMBA = .FALSE.
      OLDNRM = .FALSE.
      TOLLRN = 0.0D0
      ZCMVAL = 1.0D0
      LCMMAX = -1
C
C     We also need two variables from cbiher.h
C
      CALL RDCBIH(IPREAD, INPTST)
C
C     Now call READIN to get the specification of the primitive basis.
C     Note that UNCONT = .true., to force decontracted specification
C     for Douglas-Kroll integrals.
C
      UNCONTBKP = UNCONT
      UNCONT = .TRUE.
      CALL READIN(WORK,LWORK,.TRUE.)
      UNCONT = UNCONTBKP
      CALL SETDCH
      RETURN
      END
C/* Deck dkdriv */
      SUBROUTINE DKDRIV(WORK,LWORK,IPRINT)
C
C     Driver routine for calculation of the various Douglas-Kroll matrices
C     in the primitive basis
C
C     K.Ruud, September 2000
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      DIMENSION WORK(LWORK)
#include <shells.h>
#include <symmet.h>
#include <inftap.h>
#include <ibtfun.h>
C
      CALL GPOPEN(LUPROP,'DKINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
      NBAST = 0
      DO 400 ISHELL = 1, KMAX
         DO 410 KB = 0,MAXREP
            IF (IBTAND(KB,ISTBAO(ISHELL)) .EQ. 0)
     *            NBAST = NBAST + KHKT(ISHELL)
  410    CONTINUE
  400 CONTINUE
      NNBASX = NBAST*(NBAST + 1)/2
      N2BASX = NBAST*NBAST
C
C     Memory allocation
C
      KOVERL = 1
      KXMAT  = KOVERL + NNBASX
      KEIGVC = KXMAT  + N2BASX
      KSQMT1 = KEIGVC + N2BASX
      KSQMT2 = KSQMT1 + N2BASX
      KTRMT1 = KSQMT2 + N2BASX
      KTRMT2 = KTRMT1 + NNBASX
      KE1    = KTRMT2 + NNBASX
      KA1    = KE1    + NBAST
      KEK    = KA1    + NBAST
      KD1    = KEK    + NBAST
      KWEWMT = KD1    + NBAST
      KVMAT  = KWEWMT + N2BASX
      KPVPMT = KVMAT  + NNBASX
      KLAST  = KPVPMT + NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DKDRIV','DKDRI1',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
C     Call the driver routine
C     
      CALL DKDRI1(WORK(KOVERL),WORK(KXMAT),WORK(KEIGVC),WORK(KSQMT1),
     &            WORK(KSQMT2),WORK(KTRMT1),WORK(KTRMT2),WORK(KWEWMT),
     *            WORK(KE1),WORK(KA1),WORK(KEK),WORK(KD1),WORK(KVMAT),
     &            WORK(KPVPMT),WORK(KLAST),LWRK,NBAST,NNBASX,N2BASX,
     &            IPRINT)
C
      RETURN
      END
C/* Deck dkdri1 */
      SUBROUTINE DKDRI1(OVERLP,XMAT,EIGVEC,SQMAT1,SQMAT2,TRMAT1,TRMAT2,
     &                  WEWMAT,E1DIAG,A1DIAG,EKDIAG,D1DIAG,VMAT,PVPMAT,
     &                  WORK,LWORK,NBAST,NNBASX,N2BASX,IPRINT)
C
C     K.Ruud, Sept.-2000
C
#include <implicit.h>
#include <codata.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, CVEL2 = CVEL*CVEL, DM1 = -D1)
      CHARACTER*8 RTNLBL(2)
      LOGICAL FNDLAB
      DIMENSION OVERLP(NNBASX), XMAT(NBAST,NBAST), EIGVEC(NBAST,NBAST),
     &          SQMAT1(NBAST,NBAST), SQMAT2(NBAST,NBAST),
     &          TRMAT1(NNBASX), TRMAT2(NNBASX), E1DIAG(NBAST),
     &          A1DIAG(NBAST), EKDIAG(NBAST), D1DIAG(NBAST),
     &          WEWMAT(N2BASX), VMAT(NNBASX), PVPMAT(NNBASX)
      DIMENSION WORK(LWORK)
#include <dkinf.h>      
#include <priunit.h>
#include <inftap.h>
#include <aosotr.h>
#include <symmet.h>
C
C     Allocate temporary storage
C
      KSCR1 = 1
      KSCR2 = KSCR1 + N2BASX
      KLAST = KSCR2 + N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DKDRI1','START ',KLAST,LWORK)
C
C     Read in overlap matrix
C
      IF (FNDLAB('OVERLAP ',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,OVERLP)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Overlap matrix read from LUPROP in DKDRI1',-1)
            CALL OUTPAK(OVERLP,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: Overlap matrix'//
     *      ' not found on LUPROP.'
         CALL QUIT('HUCDR1 error: property not found on LUPROP')
      END IF
C
C     Diagonalize the overlap matrix
C
      CALL DCOPY(NNBASX,OVERLP,1,TRMAT1,1)
      CALL DUNIT(EIGVEC,NBAST)
      CALL JACO(TRMAT1,EIGVEC,NBAST,NBAST,NBAST,WORK(KSCR1),WORK(KSCR2))
C
      CALL DZERO(XMAT,N2BASX)
      DO I = 1, NBAST
         II = I*(I + 1)/2
         E1DIAG(I) = D1/SQRT(TRMAT1(II))
      END DO
      DO I = 1, NBAST
         CONST =E1DIAG(I)
         DO J = 1, NBAST
            XMAT(J,I) = EIGVEC(J,I)*CONST
         END DO
      END DO
C
C     Read in the kinetic energy and transform to orthonormal basis
C     Note that we reuse OVERLP to save memory
C     
      REWIND (LUPROP)
      IF (FNDLAB('KINENERG',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,OVERLP)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Kinetic energy matrix read from LUPROP '//
     &                  'in DKDRI1',-1)
            CALL OUTPAK(OVERLP,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: Kinetic energy'//
     *      ' matrix not found on LUPROP.'
         CALL QUIT('HUCDR1 error: property not found on LUPROP')
      END IF
C     
      CALL UTHU(OVERLP,TRMAT1,XMAT,WORK(KSCR1),NBAST,NBAST)
C
C     Diagonalize the kinetic energy matrix in orthonormal basis
C
      CALL DUNIT(EIGVEC,NBAST)
      CALL JACO(TRMAT1,EIGVEC,NBAST,NBAST,NBAST,WORK(KSCR1),WORK(KSCR2))
C
C     Extract kinetic energy and construct Epi, Api and Dpi
C     
      DO I = 1, NBAST
         II = I*(I + 1)/2
         DIAG = D2*TRMAT1(II)
         EKDIAG(I) = DIAG
         E1DIAG(I) = CVEL*SQRT(DIAG + CVEL2)
         A1DIAG(I) = SQRT((CVEL2 + E1DIAG(I))/(D2*E1DIAG(I)))
         D1DIAG(I) = CVEL*A1DIAG(I)/(CVEL2 + E1DIAG(I))
         E1DIAG(I) = E1DIAG(I) - CVEL2
      END DO
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Diagonal E1 matrix',-1)
         CALL OUTPUT(E1DIAG,1,NBAST,1,1,NBAST,1,1,LUPRI)
         CALL HEADER('Diagonal A1 matrix',-1)
         CALL OUTPUT(A1DIAG,1,NBAST,1,1,NBAST,1,1,LUPRI)
         CALL HEADER('Diagonal Ek matrix',-1)
         CALL OUTPUT(EKDIAG,1,NBAST,1,1,NBAST,1,1,LUPRI)
         CALL HEADER('Diagonal D1 matrix',-1)
         CALL OUTPUT(D1DIAG,1,NBAST,1,1,NBAST,1,1,LUPRI)
      END IF
C
C
C     Read in nuclear attraction and pVp integrals from file LUPROP (DKINT)
C     and transform to orthonormal basis. Then we transform them to the 
C     diagonal basis of the kinetic operator.
C
      REWIND (LUPROP)
      IF (FNDLAB('POTENERG',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,VMAT)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Potential energy matrix read from LUPROP '//
     &                  'in DKDRI1',-1)
            CALL OUTPAK(VMAT,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: potential energy'//
     *      ' matrix not found on LUPROP.'
         CALL QUIT('DKDRI1 error: property not found on LUPROP')
      END IF
      CALL UTHU(VMAT,OVERLP,XMAT,WORK(KSCR1),NBAST,NBAST)
      CALL UTHU(OVERLP,VMAT,EIGVEC,WORK(KSCR1),NBAST,NBAST)
C     
      REWIND (LUPROP)
      IF (FNDLAB('pVpINTEG',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,PVPMAT)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('pVp integrals read from LUPROP '//
     &                  'in DKDRI1',-1)
            CALL OUTPAK(PVPMAT,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: pVp integral'//
     *      ' matrix not found on LUPROP.'
         CALL QUIT('DKDRI1 error: property not found on LUPROP')
      END IF
C
      CALL UTHU(PVPMAT,WORK(KSCR2),XMAT,WORK(KSCR1),NBAST,NBAST)
      CALL UTHU(WORK(KSCR2),PVPMAT,EIGVEC,WORK(KSCR1),NBAST,NBAST)
C
C     Create second-order Douglas-Kroll contributions
C
      CALL DZERO(WEWMAT,N2BASX)
      CALL WWEWEW(VMAT,PVPMAT,E1DIAG,A1DIAG,EKDIAG,TRMAT1,TRMAT2,
     &            SQMAT1,SQMAT2,WEWMAT,WORK(KSCR1),NBAST,NNBASX,N2BASX,
     &            IPRINT)
C
C     Backtransform the second-order Douglas-Kroll matrix to 
C     orthonormal basis and make it square
C
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           WEWMAT,NBAST,0.D0,
     &           SQMAT1,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           EIGVEC,NBAST,0.D0,
     &           WEWMAT,NBAST)
      CALL DSCAL(N2BASX,-1.0D0,WEWMAT,1)
C
C     -1/2WWE - 1/2 EWW - WEW
C
C
C     We now backtransform the various matrices to primitive basis
C     and construct the full effective Douglas-Kroll Hamiltonian
C
      CALL DZERO(SQMAT1,N2BASX)
      DO I = 1, NBAST
         SQMAT1(I,I) = E1DIAG(I) - CVEL2
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           SQMAT1,NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT2,NBAST,
     &           EIGVEC,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     AHA
C
      CALL DZERO(SQMAT1,N2BASX)
      DO I = 1, NBAST
         SQMAT1(I,I) = A1DIAG(I)
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           SQMAT1,NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT2,NBAST,
     &           EIGVEC,NBAST,0.D0,
     &           SQMAT1,NBAST)
      CALL DSPTSI(NBAST,OVERLP,WORK(KSCR1))
      CALL DGEMM('T','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           WORK(KSCR1),NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,-1.D0,
     &           SQMAT2,NBAST,
     &           SQMAT1,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     DWD
C
      CALL DZERO(SQMAT1,N2BASX)
      DO I = 1, NBAST
         SQMAT1(I,I) = D1DIAG(I)
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           SQMAT1,NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT2,NBAST,
     &           EIGVEC,NBAST,0.D0,
     &           SQMAT1,NBAST)
      CALL DSPTSI(NBAST,WORK(KSCR2),WORK(KSCR1))
      CALL DGEMM('T','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           WORK(KSCR1),NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,-1.D0,
     &           SQMAT2,NBAST,
     &           SQMAT1,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     Pack and backtransform the final effect Hamiltonian to primitive basis
C
      CALL DSITSP(NBAST,WEWMAT,TRMAT1)
      CALL DGEINV(NBAST,XMAT,SQMAT1,WORK(KSCR1),WORK(KSCR2),INFO)
      CALL UTHU(TRMAT1,TRMAT2,SQMAT1,WORK(KSCR1),NBAST,NBAST)
C
C     Write the final effective Hamiltonian to DKINT and keep until we
C     know what kind of contraction coefficients we have
C
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('Final Douglas-Kroll effective one-electron'//
     &               ' potential (SO basis)',-1)
         CALL OUTPAK(TRMAT2,NBAST,1,LUPRI)
      END IF
C
C     We transform to AO basis in order to more easily contract it
C
      IF (MAXREP .GT. 0) THEN
         IATOM = 0
         CALL TRSOAO(TRMAT2,WORK,LWORK,NBAST,IATOM,IPRINT)
         CALL DCOPY(NNBASX,WORK,1,TRMAT2,1)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Final Douglas-Kroll effective one-electron'//
     &                 ' potential (AO basis)',-1)
            CALL OUTPAK(TRMAT2,NBAST,1,LUPRI)
         END IF
      END IF
      NPRIM  = NBAST
      NNPRMX = NNBASX
      REWIND(LUPROP)
      CALL WRTPRO(TRMAT2,NNBASX,'DKINTPRI',RTNLBL)
      CALL GPCLOSE(LUPROP,'KEEP')
C
C     End of DKDRIV
C     
      RETURN
      END
C/* Deck wwewew */
      SUBROUTINE WWEWEW(VMAT,PVPMAT,E1DIAG,A1DIAG,EKDIAG,VTMP,GTMP,
     &                  SQMAT1,SQMAT2,WEWMAT,RDIAG,NBAST,NNBASX,N2BASX,
     &                  IPRINT)
C
C     Subroutine for calculating second-order Douglas-Kroll contributions
C
#include <implicit.h>
#include <codata.h>
      PARAMETER (CVEL2 = CVEL*CVEL, DP5 = 0.50D0)
#include <priunit.h>
      DIMENSION VMAT(NNBASX), PVPMAT(NNBASX), E1DIAG(NBAST),
     &          A1DIAG(NBAST), EKDIAG(NBAST), WEWMAT(NBAST,NBAST),
     &          SQMAT1(NBAST,NBAST), SQMAT2(NBAST,NBAST),
     &          VTMP(NNBASX), GTMP(NNBASX), RDIAG(NBAST)
C
      DO I = 1, NBAST
         E1DIAG(I) = E1DIAG(I) + CVEL2
         RDIAG(I) = CVEL/(CVEL2 + E1DIAG(I))
      END DO
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
            VTMP(IJ) = VMAT(IJ)/(E1DIAG(I) + E1DIAG(J))
            GTMP(IJ) = PVPMAT(IJ)/(E1DIAG(I) + E1DIAG(J))
            SQMAT1(I,J) = A1DIAG(I)*RDIAG(I)*GTMP(IJ)*A1DIAG(J)
     &                   *A1DIAG(J)
            SQMAT1(J,I) = A1DIAG(J)*RDIAG(J)*GTMP(IJ)*A1DIAG(I)
     &                   *A1DIAG(I)
            SQMAT2(I,J) = RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                  - GTMP(IJ)*A1DIAG(J)*RDIAG(J)/EKDIAG(I)
            SQMAT2(J,I) = RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                  - GTMP(IJ)*A1DIAG(I)*RDIAG(I)/EKDIAG(J)
            IJ = IJ + 1
         END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,0.D0,
     &           WEWMAT,NBAST)
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('1.transformation of W-matrix')
         CALL OUTPUT(WEWMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
           SQMAT1(I,J) = A1DIAG(I)*VTMP(IJ)*A1DIAG(J)*A1DIAG(J)*RDIAG(J)
           SQMAT1(J,I) = A1DIAG(J)*VTMP(IJ)*A1DIAG(I)*A1DIAG(I)*RDIAG(I)
           SQMAT2(I,J) = -EKDIAG(I)*RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                 + GTMP(IJ)*A1DIAG(J)*RDIAG(J)
           SQMAT2(J,I) = -EKDIAG(J)*RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                 + GTMP(IJ)*A1DIAG(I)*RDIAG(I)
           IJ = IJ + 1
        END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,1.D0,
     &           WEWMAT,NBAST)
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('2.transformation of W-matrix')
         CALL OUTPUT(WEWMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     Make 1/2 EW*W + 1/2 W*WE
C
      DO I = 1, NBAST
         DO J = 1, NBAST
            WEWMAT(I,J) = DP5*WEWMAT(I,J)*(E1DIAG(I)+E1DIAG(J))
         END DO
      END DO
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('3.transformation of W-matrix')
         CALL OUTPUT(WEWMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     Start constructing WEW
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
           SQMAT1(I,J) = A1DIAG(I)*RDIAG(I)*GTMP(IJ)*A1DIAG(J)
     &                  *A1DIAG(J)*E1DIAG(J)
           SQMAT1(J,I) = A1DIAG(J)*RDIAG(J)*GTMP(IJ)*A1DIAG(I)
     &                  *A1DIAG(I)*E1DIAG(I)
           SQMAT2(I,J) = RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                 - GTMP(IJ)*A1DIAG(J)*RDIAG(J)/EKDIAG(I)
           SQMAT2(J,I) = RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                 - GTMP(IJ)*A1DIAG(I)*RDIAG(I)/EKDIAG(J)
           IJ = IJ + 1
        END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
            SQMAT1(I,J) = A1DIAG(I)*VTMP(IJ)*A1DIAG(J)*A1DIAG(J)
     &                   *RDIAG(J)*E1DIAG(J)
            SQMAT1(J,I) = A1DIAG(J)*VTMP(IJ)*A1DIAG(I)*A1DIAG(I)
     &                   *RDIAG(I)*E1DIAG(I)
            SQMAT2(I,J) = -EKDIAG(I)*RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                  + GTMP(IJ)*A1DIAG(J)*RDIAG(J)
            SQMAT2(J,I) = -EKDIAG(J)*RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                  + GTMP(IJ)*A1DIAG(I)*RDIAG(I)
            IJ = IJ + 1
         END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     The entire matrix is now constructed in the diagonal T basis
C
      RETURN
      END
C/* Deck dkcont */
      SUBROUTINE DKCONT(PROPRI,IPRINT,WORK,LWORK)
C
C     This subroutine reads in the Douglas-Kroll potential in the primitive
C     basis and contracts it with the contraction coefficients of the 
C     basis set used in the calculation
C
C     K.Ruud, September 2000
C     
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <dkinf.h>      
      LOGICAL FNDLAB, ANTI, PROPRI
      CHARACTER*8 LABINT
      DIMENSION WORK(LWORK)
#include <shells.h>
#include <symmet.h>
C
C     
      NBAST = ISUM(MAXREP+1,NAOS,1)
      NNBASX = NBAST*(NBAST + 1)/2
C
      KWMAT  = 1
      KCONTC = KWMAT + NPRIM*NPRIM 
      KLAST  = KCONTC + NPRIM*NBAST
      KTMP   = KLAST + NPRIM*NPRIM
      LWRK   = LWORK - KLAST + 1
      IF (KTMP .GT. LWORK) CALL STOPIT('DKCONT','PRI2CN',KLAST,LWORK)
      LUDKIN = -9000
      CALL GPOPEN(LUDKIN,'DKINT','OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      IF (FNDLAB('DKINTPRI',LUDKIN)) THEN
         IF (NNPRMX .GT. LWORK) CALL STOPIT('DKCONT','WMAT  ',
     &                                      NNPRMX,LWORK)
         CALL READT(LUDKIN,NNPRMX,WORK(KWMAT))
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Primitive Douglas-Kroll integrals read from '//
     &                  'LUDKIN in DKCONT',-1)
            CALL OUTPAK(WORK(KWMAT),NPRIM,1,LUPRI)
         END IF
         CALL DCOPY(NNPRMX,WORK(KWMAT),1,WORK(KLAST),1)
         CALL DSPTSI(NPRIM,WORK(KLAST),WORK(KWMAT))
         NOPTYP = 1
         LABINT = 'DKPOTINT'
         ANTI   = .FALSE.
         INTREP = 0
         CALL PRI2CN(WORK(KWMAT),WORK(KCONTC),WORK(KLAST),
     &               LWRK,ANTI,NOPTYP,NBAST,NPRIM,LABINT,
     &               INTREP,PROPRI,IPRINT)
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKCONT error: Primitive Douglas'//
     &      '-Kroll integrals not found on LUDKIN.'
         CALL QUIT('DKCONT error: Integrals not found on LUDKIN')
      END IF
      CALL GPCLOSE(LUDKIN,'DELETE')
      RETURN
      END
C/* Deck pri2cn */
      SUBROUTINE PRI2CN(AOINT,CONTC,WORK,LWORK,ANTI,
     &                  NOPTYP,NBAST,NPRIM,LABINT,INTREP,PROPRI,IPRINT)
C
C     Subroutine to transform primitive AO's to contracted SO basis
C
#include <implicit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <pi.h>
      PARAMETER (D2 = 2.0D0, DP5 = 0.50D0, DP75 = 0.75D0, D4 = 4.0D0)
#include <priunit.h>
      LOGICAL ANTI, PROPRI
      CHARACTER*8 RTNLBL(2), LABINT(NOPTYP)
      DIMENSION AOINT(NPRIM,NPRIM), CONTC(NPRIM,NBAST), WORK(LWORK)
      DIMENSION INTREP(NOPTYP)
#include <pgroup.h>
#include <nuclei.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <ibtfun.h>
C
      ICNT   = 0
      IADR   = 0
      IORBA  = 0
      IF (NOPTYP.GT. 1) CALL QUIT('More than one integral type in'//
     &                            'PRI2CN')
      NELMNT = NBAST*(NBAST + 1)/2
C
      LUCNMT = -9000
      CALL GPOPEN(LUCNMT,'CNTMAT','OLD',' ','FORMATTED',IDUMMY,.FALSE.)
      READ (LUCNMT,*) ((CONTC(I,J),I=1,NPRIM),J=1,NBAST)
      CALL GPCLOSE(LUCNMT,'DELETE')
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Contraction matrix read from '//
     &               'LUCNMT in PRI2CN',-1)
         CALL OUTPUT(CONTC,1,NPRIM,1,NBAST,NPRIM,NBAST,1,LUPRI)
      END IF
      CALL DGEMM('T','N',NBAST,NPRIM,NPRIM,1.0D0,CONTC,NPRIM,
     &           AOINT,NPRIM,0.0D0,WORK,NBAST)
      CALL DGEMM('N','N',NBAST,NPRIM,NPRIM,1.0D0,WORK,NBAST,
     &           CONTC,NPRIM,0.0D0,AOINT,NBAST)
      CALL DSITSP(NBAST,AOINT,WORK)
      CALL DCOPY(NELMNT,WORK,1,AOINT,1)
C
C     First we transform to SO basis, then for consistency with WRTUND,
C     we need to symmetry pack the integrals
C
      IF (IPRINT .GE. 6) THEN
         CALL AROUND('Contracted Douglas-Kroll '//
     &               'integrals in AO basis')
         CALL OUTPAK(AOINT,NBAST,1,LUPRI)
      END IF
      IF (MAXREP .GT. 0) THEN
         IATOM = -1
         CALL TRAOSO(AOINT,WORK,LWORK,NBAST,IATOM,IPRINT)
         CALL DCOPY(NELMNT,WORK,1,AOINT,1)
         CALL PKSYM1(AOINT,AOINT,NAOS,MAXREP+1,2)
      END IF
C
C     Write integrals to file AOPROPER and print them if people are
C     interested
C     
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      IF (ANTI) THEN
         RTNLBL(2)='ANTISYMM'
      ELSE
         RTNLBL(2)='SYMMETRI'
      END IF
C     
      DO 600 I = 1, NOPTYP
         IF (PROPRI .OR. IPRINT .GE. 4) THEN
            CALL AROUND('Integrals of operator: '//LABINT(I))
            WRITE (LUPRI,'(A,2X,A3,A1,I2,A1)')
     &           ' Symmetry of operator:',
     &           REP(INTREP(I)),'(',(INTREP(I) + 1),')'
            CALL OUTPAK(AOINT,NBAST,1,LUPRI)
         END IF
         CALL WRTPRO(AOINT,NELMNT,LABINT(I),RTNLBL)
 600  CONTINUE
C     
      RETURN
      END
C/* Deck rdcbih */
      SUBROUTINE RDCBIH(IPDF, INTST)
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiher.h>
      LOGICAL INTST
      IPDF  = IPRDEF
      INTST = TSTINP
      RETURN
      END
