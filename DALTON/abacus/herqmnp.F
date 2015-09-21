!...
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C  /* Deck qmnpinp */
      SUBROUTINE QMNPINP(WORD)
C
C Purpose:
C   reads "*QMNPMM" input group and sets up various flags for
C   QM/NP/MM type embedding
C
C Input:
C   WORD - DALTON input group keyword
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "gnrinf.h"
#include "qmnpmm.h"

#ifdef VAR_MPI
#include "iprtyp.h"
#endif

C
      PARAMETER ( NTABLE = 13)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      DATA TABLE /'.QMNP  ', '.QMNPMM', '.QMMM  ', '.NPPOIN',
     &            '.MMGAUS', '.PRINT ', '.MQITER', '.CMXPOL',
     &            '.NONPCA', '.MMPOLA', '.MMCAPA', '.XXXXXX',
     &            '.DAMPED'/
C
      CALL QENTER('QMNPINP')
C
      CALL SET_QMNPMM()
C
      NEWDEF = (WORD .EQ. '*QMNPMM')
      ICHANG = 0
C
      IF (NEWDEF) THEN
        WORD1 = WORD
C
  100   CONTINUE
C
        READ (LUCMD, '(A7)') WORD
        CALL UPCASE(WORD)
        PROMPT = WORD(1:1)
C
        IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
          GO TO 100
        ELSE IF (PROMPT .EQ. '.') THEN
          ICHANG = ICHANG + 1
          DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
              GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13), I
            END IF
  200     CONTINUE
          IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
          END IF
C
          WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,
     *                           '" not recognized for *QMNPMM'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword for *QMNPMM')
C
C         ".QMNP "
    1     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .TRUE.
            DOMMSUB = .FALSE.
            DONPCAP = .TRUE.
            DONPPOL = .TRUE.
          GO TO 100
C
C         ".QMNPMM"
    2     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .TRUE.
            DOMMSUB = .TRUE.
            DONPCAP = .TRUE.
            DONPPOL = .TRUE.
          GO TO 100
C
C         ".QMMM  "
    3     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .FALSE.
            DOMMSUB = .TRUE.
          GO TO 100
C
C         ".NPPOIN"
    4     CONTINUE
            NPMQGAU = .FALSE.
          GO TO 100
C
C         ".MMGAUS"
    5     CONTINUE
            MMMQGAU = .TRUE.
          GO TO 100
C
C         ".PRINT "
    6     CONTINUE
            READ(LUCMD,*) IPRTLVL
          GO TO 100
C
C         ".MQITER"
    7     CONTINUE
            MQITER = .TRUE.
          GO TO 100
C
C         ".CMXPOL"
    8     CONTINUE
            CMXPOL = .TRUE.
          GO TO 100
C
C         ".NONPCA"
    9     CONTINUE
            DONPCAP = .FALSE.
          GO TO 100
C
C         ".MMPOLA"
   10     CONTINUE
            DOMMPOL = .TRUE.
          GO TO 100
C
C         ".MMCAPA"
   11     CONTINUE
            CMXPOL = .TRUE.
          GO TO 100
C
C         ".XXXXXX"
   12     CONTINUE
          GO TO 100
C         ".DAMPED"
 13       CONTINUE
            NOVDAMP = .FALSE.
          GO TO 100

C
        ELSE IF (PROMPT .EQ. '*') THEN
          GO TO 300
        ELSE
          WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,
     *                             '" not recognized for *QMNPMM'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal prompt for *QMNPMM')
        END IF
      END IF
C
  300 CONTINUE
C
      CALL QEXIT('QMNPINP')
      RETURN
      END
C  /* Deck set_qmnpmm */
      SUBROUTINE SET_QMNPMM()
C
C Purpose:
C   sets initial values of all parameters related to
C   QM/NP/MM type embedding
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      DONPSUB = .FALSE.
      DOMMSUB = .FALSE.
      DONPCAP = .FALSE.
      DOMMCAP = .FALSE.
      DONPPOL = .FALSE.
      DOMMPOL = .FALSE.
      NPMQGAU = .TRUE.
      MMMQGAU = .FALSE.
      MQITER  = .FALSE.
      CMXPOL  = .FALSE.
      NOVDAMP = .TRUE.
      TNPBLK  = 0
      TMMBLK  = 0
      IPRTLVL = 1
      TNPATM  = 0
      TMMATM  = 0
      TPOLATM = 0
      TNPFF   = 0
      TMMFF   = 0
C
      CALL DZERO(NPCORD,3*MXNPATM)
      CALL DZERO(MMCORD,3*MXMMATM)
      CALL DZERO(NPCHRG,MAXBLK)
      CALL DZERO(MMCHRG,MAXBLK)
      CALL IZERO(NPFTYP,MXNPATM)
      CALL IZERO(MMFTYP,MXMMATM)
      CALL IZERO(NPATOM,MAXBLK)
      CALL IZERO(MMATOM,MAXBLK)
      CALL DZERO(NPFPOL,MXNPFF)
      CALL DZERO(NPFCAP,MXNPFF)
      CALL DZERO(NPFOMG1,MXNPFF)
      CALL DZERO(NPFGAM1,MXNPFF)
      CALL DZERO(NPFOMG2,MXNPFF)
      CALL DZERO(NPFGAM2,MXNPFF)
      CALL DZERO(NPFFAC,MXNPFF)
      CALL DZERO(MMFM0,MXMMFF)
      CALL DZERO(MMFPOL,MXMMFF)
      CALL DZERO(MMMOL,MXMMATM)
      CALL IZERO(MMSKIP,MXMMATM)
C
      ENSOLQNP = 0.0D0
      EESOLQNP = 0.0D0
      ENSOLMNP = 0.0D0
      EESOLMNP = 0.0D0
C
      ENSOLQMM = 0.0D0
      EESOLQMM = 0.0D0
      ENSOLMMM = 0.0D0
      EESOLMMM = 0.0D0
C
      RETURN
      END
C  /* Deck read_pot_qmnpmm */
      SUBROUTINE READ_POT_QMNPMM()
C
C Purpose:
C   reads NP and MM subsystems geometry, total charges, and
C   force field data from POTENTIAL.INP file
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      PARAMETER (XFACT = 1.0D0/XTANG)
C
      CHARACTER UNITS*2, NPWORD*2, FFWORD*4, NPLBL(MXNPATM)*2
      CHARACTER MMWORD*2, MMLBL(MXMMATM)*2
C
      CALL QENTER('READ_POT_QMNPMM')
C     Open POTENTIAL.INP file
      LUQMNP=-1
      CALL GPOPEN(LUQMNP,'POTENTIAL.INP','OLD',' ',
     &            'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUQMNP)
C     Read geometry input units, number of NP and MM
C     separate subsystems
      READ(LUQMNP,*) UNITS, TNPBLK, TMMBLK
      CALL UPCASE(UNITS)
C     Check input consistency
      IF ((UNITS.NE.'AA').AND.(UNITS.NE.'AU')) THEN
         CALL QUIT('Unknown units in POTENTIAL.INP')
      ENDIF
      IF ((DONPSUB).AND.(TNPBLK.LT.1)) THEN
         CALL QUIT('NP system is missing in POTENTIAL.INP')
      END IF
      IF ((DOMMSUB).AND.(TMMBLK.LT.1)) THEN
         CALL QUIT('MM system is missing in POTENTIAL.INP')
      END IF
      IF (TNPBLK.GT.MAXBLK) THEN
        WRITE(LUPRI,'(/2X,A)')
     &        'Maximum number of NP subsystems exceeded!'
        WRITE(LUPRI,'(2X,A,I3,2X,A,I3)') 'Input TNPBLK=',
     &        TNPBLK, 'Maximum allowed:', MAXBLK
        CALL QUIT('Increase MAXBLK in qmnpmm.h')
      END IF
      IF (TMMBLK.GT.MAXBLK) THEN
        WRITE(LUPRI,'(/2X,A)')
     &        'Maximum number of MM subsystems exceeded!'
        WRITE(LUPRI,'(2X,A,I3,2X,A,I3)') 'Input TNPBLK=',
     &        TNPBLK, 'Maximum allowed:', MAXBLK
        CALL QUIT('Increase MAXBLK in qmnpmm.h')
      END IF
C     Read NP subsystems data
      IF (DONPSUB) THEN
         ISTART = 0
         DO I=1,TNPBLK
C           Read NP subsystem header
            READ(LUQMNP,*) NPWORD, NPCHRG(I), NPATOM(I)
C           Check input consistency
            CALL UPCASE(NPWORD)
            IF (NPWORD.NE.'NP') THEN
               WRITE(LUPRI,'(/2X,A,I2,1X,A)')
     &        'Incorrect NP subsystem I=', I,
     &        'header in POTENTIAL.INP file!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
            IF (NPATOM(I).LE.0) THEN
               WRITE(LUPRI,'(/2X,A,I5,A,I2,A)')
     &        'Incorrect number of atoms ', NPATOM(I),
     &        'in NP subsystem', I, '!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
C           Read NP subsystem data
            DO J=1,NPATOM(I)
               JOFF = ISTART + J
               READ(LUQMNP,*) NPLBL(JOFF), NPCORD(1,JOFF),
     &              NPCORD(2,JOFF), NPCORD(3,JOFF),
     &              NPFTYP(JOFF)
            END DO
            ISTART = ISTART + NPATOM(I)
            TNPATM = TNPATM + NPATOM(I)
            IF (TNPATM.GT.MXNPATM) THEN
               WRITE(LUPRI,'(/2X,A)')
     &         'Maximum number of NP atoms exceeded!'
               WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')
     &         'Input TNPATM=', TNPATM,
     &         'Maximum allowed:', MXNPATM
               CALL QUIT('Increase MXNPATM in qmnpmm.h')
            END IF
         END DO
C        Convert to atomic units if neeeded
         IF (UNITS.EQ.'AA') THEN
            CALL DSCAL(3*TNPATM,XFACT,NPCORD,1)
         END IF
C        Read force field data
         READ(LUQMNP,*) FFWORD, TNPFF
         CALL UPCASE(FFWORD)
C        Check NP force field data consistency
         IF (FFWORD.NE.'NPFF') THEN
            WRITE(LUPRI,'(/2X,A,A)') 'Corrupted NP force field ',
     &           'header in POTENTIAL.INP file!'
            CALL QUIT('Corrupted POTENTIAL.INP')
         END IF
         IF (TNPFF.GT.MXNPFF) THEN
            WRITE(LUPRI,'(/2X,A)')
     &            'Maximum number of NP force field types exceeded!'
            WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')
     &            'Input TNPFF=', TNPFF,
     &            'Maximum allowed:', MXNFF
            CALL QUIT('Increase MXNPFF in qmnpmm.h')
         END IF
C        Read force field data
         DO I=1,TNPFF
            READ(LUQMNP,*) NPFPOL(I), NPFCAP(I), NPFOMG1(I),
     &                     NPFGAM1(I), NPFOMG2(I), NPFGAM2(I),
     &                     NPFFAC(I)
         END DO
C        Check force field definitions
         DO I=1,TNPATM
            IF (NPFTYP(I).GT.TNPFF) THEN
               WRITE(LUPRI,'(/2X,A,I4,1X,A)')
     &         'Unknown force field requested for atom', I,
     &         'in NP region!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
         END DO
         IF (IPRTLVL.GE.5) THEN
            CALL PRINT_NPREGION(NPLBL)
         END IF
      END IF
C     MM subsystems input
      IF (DOMMSUB) THEN
         ISTART = 0
         DO I=1,TMMBLK
C           Read MM subsystem header
            READ(LUQMNP,*) MMWORD,  MMATOM(I)
C           Check input consistency
            CALL UPCASE(MMWORD)
            IF (MMWORD.NE.'MM') THEN
               WRITE(LUPRI,'(/2X,A,I2,1X,A)')
     &        'Incorrect MM subsystem I=', I,
     &        'header in POTENTIAL.INP file!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
            IF (MMATOM(I).LE.0) THEN
               WRITE(LUPRI,'(/2X,A,I5,A,I2,A)')
     &        'Incorrect number of atoms ', MMATOM(I),
     &        'in MM subsystem', I, '!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
C           Read MM subsystem data
            DO J=1,MMATOM(I)
               JOFF = ISTART + J
               READ(LUQMNP,*) MMLBL(JOFF), MMMOL(JOFF),
     &                        MMCORD(1,JOFF), MMCORD(2,JOFF),
     &                        MMCORD(3,JOFF), MMFTYP(JOFF)
            END DO
            ISTART = ISTART + MMATOM(I)
            TMMATM = TMMATM + MMATOM(I)
            IF (TMMATM.GT.MXMMATM) THEN
               WRITE(LUPRI,'(/2X,A)')
     &         'Maximum number of MM atoms exceeded!'
               WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')
     &         'Input TMMATM=', TMMATM,
     &         'Maximum allowed:', MXMMATM
               CALL QUIT('Increase MXMMATM in qmnpmm.h')
            END IF
         END DO
C        Convert to atomic units if neeeded
         IF (UNITS.EQ.'AA') THEN
            CALL DSCAL(3*TMMATM,XFACT,MMCORD,1)
         END IF
C        Read force field data
         READ(LUQMNP,*) FFWORD, TMMFF
         CALL UPCASE(FFWORD)
C        Check NP force field data consistency
         IF (FFWORD.NE.'MMFF') THEN
            WRITE(LUPRI,'(/2X,A,A)') 'Corrupted MM force field ',
     &           'header in POTENTIAL.INP file!'
            CALL QUIT('Corrupted POTENTIAL.INP')
         END IF
C        Read force field data
         DO I=1,TMMFF
            IF ((.NOT.DOMMCAP).AND.(.NOT.DOMMPOL)) THEN
               READ(LUQMNP,*) MMFM0(I)
            END IF
            IF ((.NOT.DOMMCAP).AND.DOMMPOL) THEN
               READ(LUQMNP,*) MMFM0(I),MMFPOL(I)
            END IF
         END DO
C        Check force field definitions
         DO I=1,TMMATM
            IF (MMFTYP(I).GT.TMMFF) THEN
               WRITE(LUPRI,'(/2X,A,I4,1X,A)')
     &         'Unknown force field requested for atom', I,
     &         'in MM region!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
C           Set polarization centers
            IF (MMFPOL(MMFTYP(I)) .GT. 1.0D-6) THEN
               TPOLATM = TPOLATM+1
               MMSKIP(I) = 1
C              works only for one non-metallic MM region
            END IF
         END DO
         IF (IPRTLVL.GE.5) THEN
            CALL PRINT_MMREGION(MMLBL)
         END IF
      END IF
C     Close POTENTIAL.INP file
      CALL GPCLOSE(LUQMNP,'KEEP')
C
      CALL QEXIT('READ_POT_QMNPMM')
C
      RETURN
      END
C  /* Deck print_npregion */
      SUBROUTINE PRINT_NPREGION(ATMLBL)
C
C Purpose:
C     prints detailed information about nanoparticles
C     in QM/NP/MM embedding
C
C Input:
C  ATMLBL - List of atom labels for NP region
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      CHARACTER ATMLBL(MXNPATM)*2
C
      IF (DOMMSUB) THEN
         CALL TITLER('QM/NP/MM Embedding: Nanoparticle(s) data','*',103)
      ELSE
        CALL TITLER('QM/NP Embedding: Nanoparticle(s) data','*',103)
      END IF
C
      IF (TNPBLK.GT.1) THEN
         WRITE(LUPRI,'(/2X,A,I3,1X,A)') 'Nanoparticle region contains',
     &         TNPBLK, 'separate nanoparticles.'
      ELSE
        WRITE(LUPRI,'(/2X,A,A)') 'Nanoparticle region contains single',
     &        ' nanoparticle.'
      END IF
      ISTART = 0
      DO I=1,TNPBLK
         WRITE(LUPRI, '(/,2X,A,I3,1X,A,F6.3,A)') 'Nanoparticle',I,
     &         'with charge', NPCHRG(I), '. Coordinates in au.'
         WRITE(LUPRI, '(2X,A)')
     &         '--------------------------------------------------'
         WRITE(LUPRI, '(2X,A)')
     &         'Atom    Coord. X    Coord. Y    Coord. Z   FF Type'
         WRITE(LUPRI, '(2X,A)')
     &         '=================================================='
         DO J=1,NPATOM(I)
            JOFF = ISTART+J
            WRITE(LUPRI,'(3X,A,1X,F11.5,1X,F11.5,1X,F11.5,3X,I4)')
     &            ATMLBL(JOFF), NPCORD(1,JOFF), NPCORD(2,JOFF),
     &            NPCORD(3,JOFF), NPFTYP(JOFF)
         END DO
         WRITE(LUPRI, '(2X,A)')
     &         '=================================================='
         ISTART = ISTART + NPATOM(I)
      END DO
      WRITE(LUPRI,'(/2X,A,A)') 'Static force field(s) data for ',
     &      'nanoparticle region.'
      WRITE(LUPRI, '(2X,A)')
     &      '-------------------------------------'
      WRITE(LUPRI, '(2X,A)')
     &      'FF Type Polarizability    Capacitance'
      WRITE(LUPRI, '(2X,A)')
     &      '====================================='
      DO I=1,TNPFF
         WRITE(LUPRI,'(3X,I2,5X,F11.5,5X,F11.5)') I, NPFPOL(I),
     &         NPFCAP(I)
      END DO
      WRITE(LUPRI, '(2X,A)')
     &      '====================================='
      WRITE(LUPRI,'(/2X,A,A)') 'Dynamic force field(s) data for ',
     &      'nanoparticle region.'
      WRITE(LUPRI, '(2X,A,A)')
     &      '---------------------------------------------------',
     &      '------'
      WRITE(LUPRI, '(2X,A,A)')
     &      'FF Type Omega(1)  Gamma(1)  Omega(2)  Gamma(2)',
     &      '  Factor'
      WRITE(LUPRI, '(2X,A,A)')
     &      '===================================================',
     &      '======'
      DO I=1,TNPFF
         WRITE(LUPRI,'(3X,I2,3X,5(1X,F9.5))') I, NPFOMG1(I), NPFGAM1(I),
     &         NPFOMG2(I), NPFGAM2(I), NPFFAC(I)
      END DO
      WRITE(LUPRI, '(2X,A,A)')
     &      '===================================================',
     &      '======'
C
      RETURN
      END
C  /* Deck print_mmregion */
      SUBROUTINE PRINT_MMREGION(ATMLBL)
C
C Purpose:
C     prints detailed information about MM region
C     in QM/NP/MM embedding
C
C Input:
C  ATMLBL - List of atom labels for MM region
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      CHARACTER ATMLBL(MXNPATM)*2
C
      IF (DOMMSUB) THEN
         CALL TITLER('QM/NP/MM Embedding: MM region(s) data','*',103)
      END IF
C
      IF (TMMBLK.GT.1) THEN
         WRITE(LUPRI,'(/2X,A,I3,1X,A)') 'MM region contains',
     &         TMMBLK, 'separate non-metallic subregions.'
      END IF
      ISTART = 0
      DO I=1,TMMBLK
         WRITE(LUPRI, '(/,2X,A,I3,1X,A)') 'MM subregion',I,
     &         '. Coordinates in au.'
         WRITE(LUPRI, '(2X,A)')
     &         '--------------------------------------------------'
         WRITE(LUPRI, '(2X,A)')
     &         'Atom    Coord. X    Coord. Y    Coord. Z   FF Type'
         WRITE(LUPRI, '(2X,A)')
     &         '=================================================='
         DO J=1,MMATOM(I)
            JOFF = ISTART+J
            WRITE(LUPRI,'(3X,A,1X,F11.5,1X,F11.5,1X,F11.5,3X,I4)')
     &            ATMLBL(JOFF), MMCORD(1,JOFF), MMCORD(2,JOFF),
     &            MMCORD(3,JOFF), MMFTYP(JOFF)
         END DO
         WRITE(LUPRI, '(2X,A)')
     &         '=================================================='
C        works only for one MM region
         WRITE(LUPRI, '(2X,A,I4)') 'Number of polarizable centers   : ',
     &         TPOLATM
         WRITE(LUPRI, '(2X,A,I4)') 'Number of nonpolarizable centers: ',
     &         MMATOM(I)-TPOLATM
         ISTART = ISTART + NPATOM(I)
      END DO
      WRITE(LUPRI,'(/2X,A)') 'Force field(s) data for MM region'
C
      IF ((.NOT.DOMMCAP).AND.(.NOT.DOMMPOL)) THEN
         WRITE(LUPRI, '(2X,A)') '--------------------------------'
         WRITE(LUPRI, '(2X,A)') 'FF Type MM Charge Polarizability'
         WRITE(LUPRI, '(2X,A)') '================================'
         DO I=1,TMMFF
            WRITE(LUPRI,'(3X,I2,3X,1X,F9.5,4X,F9.5,1X,I2)') I, MMFM0(I),
     &            MMFPOL(I), MMSKIP(I)
         END DO
         WRITE(LUPRI, '(2X,A)') '================================'
      END IF
      RETURN
      END
C  /* Deck getdim_relmat */
      SUBROUTINE GETDIM_RELMAT(IMATDIM, SQFLAG)
C
C Purpose:
C     determines Relay matrix dimensions.
C
C Input:
C   SQFLAG - Request total size of Relay matrix
C Output:
C  IMATDIM - Size of Relay matrix or it's dimension
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      LOGICAL SQFLAG
C
      IMATDIM = 0
C     Add nanoparticle contribution
      IF (DONPSUB) THEN
         IMATDIM = IMATDIM + 3*TNPATM
         IF (DONPCAP) THEN
            IMATDIM = IMATDIM + TNPATM
         END IF
      END IF
C     Add Lagrangian term
      IF (DONPCAP.OR.DOMMCAP) THEN
         IMATDIM = IMATDIM + 1
      END IF
C     Get requested dimmension
      IF ((.NOT.MQITER).AND.SQFLAG) THEN
         IMATDIM = IMATDIM*IMATDIM
      END IF
C
      RETURN
      END
C  /* Deck getdim_relmat */
      SUBROUTINE GETDIM_MMMAT(IMATDIM, SQFLAG)
C
C Purpose:
C     determines Relay matrix dimensions for MM region.
C
C Input:
C   SQFLAG - Request total size of Relay matrix
C Output:
C  IMATDIM - Size of Relay matrix or it's dimension
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      LOGICAL SQFLAG
C
      IMATDIM = 0
C     Add molecular contribution
      IF (DOMMSUB) THEN
         IMATDIM = IMATDIM + 3*TPOLATM
      END IF
C     Get requested dimmension
      IF ((.NOT.MQITER).AND.SQFLAG) THEN
         IMATDIM = IMATDIM*IMATDIM
      END IF
C
      RETURN
      END
C  /* Deck getdim_relmat */
      SUBROUTINE COMP_RELMAT(FMAT,WORK,LWORK)
C
C Purpose:
C     Computes Relay matrix and inverts it or estimates initial
C     MQ vector values for iterative MQ vector determination
C     algorithm.
C
C Input:
C   WORK  - Temporary memory array
C   LWORK - Size of temporary memory array
C Output:
C if .not. MQITER
C  FMAT   - Inverted Relay matrix
C else
C  FMAT   - Initial guess of real part of MQ vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(*), WORK(LWORK)
C
      KFREE = 1
      LFREE = LWORK
C     Initialize arrays
      CALL GETDIM_RELMAT(IDIM,.TRUE.)
      CALL DZERO(FMAT,IDIM)
C     Reset matrix dimension parameter
      CALL GETDIM_RELMAT(IDIM,.FALSE.)
      IF (.NOT.MQITER) THEN
C        Compute polarizabilty dependent terms
         IF (DONPPOL.OR.DOMMPOL) THEN
            CALL GET_AMAT(FMAT,IDIM)
         END IF
C        Compute capacitancy dependent term
         IF (DONPCAP.OR.DOMMCAP) THEN
           CALL GET_CMAT(FMAT,IDIM)
           CALL GET_MMAT(FMAT,IDIM)
           CALL GET_QLAG(FMAT,IDIM)
         END IF
      ELSE
C       Conjugated gradient method via paradiso solver
      END IF
C     Print Relay matrix
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed Relay matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
      END IF
C     Invert Relay matrix
      IF (.NOT.MQITER) THEN
         CALL MEMGET('INTE',KIPIV,IDIM,WORK,KFREE,LFREE)
         CALL DGETRF(IDIM,IDIM,FMAT,IDIM,WORK(KIPIV),IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('LU factorization failed in COMP_RELMAT')
         END IF
         CALL DGETRI(IDIM,FMAT,IDIM,WORK(KIPIV),WORK(KFREE),LFREE,
     &               IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('Inversion failed in COMP_RELMAT')
         END IF
         CALL MEMREL('COMP_RELMAT',WORK,1,1,KFREE,LFREE)
         IF (IPRTLVL.GE.15) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Inverted Relay matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
         END IF
      END IF
C
      RETURN
      END
C  /* Deck getdim_mmrelmat */
      SUBROUTINE COMP_MMRELMAT(FMAT,WORK,LWORK)
C
C Purpose:
C     Computes Relay matrix and inverts it or estimates initial
C     MQ vector values for iterative MQ vector determination
C     algorithm. (MM region)
C
C Input:
C   WORK  - Temporary memory array
C   LWORK - Size of temporary memory array
C Output:
C if .not. MQITER
C  FMAT   - Inverted Relay matrix
C else
C  FMAT   - Initial guess of real part of MQ vector
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(*), WORK(LWORK)
C
      KFREE = 1
      LFREE = LWORK
C     Initialize arrays
      CALL GETDIM_MMMAT(IDIM,.TRUE.)
      CALL DZERO(FMAT,IDIM)
C     Reset matrix dimension parameter
      CALL GETDIM_MMMAT(IDIM,.FALSE.)
      IF (.NOT.MQITER) THEN
C        Compute polarizabilty dependent terms
         CALL GET_MM_AMAT(FMAT,IDIM)
      ELSE
C       Conjugated gradient method via paradiso solver
      END IF
C     Print Relay matrix
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed Relay (MM) matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
      END IF
C     Invert Relay matrix
      IF (.NOT.MQITER) THEN
         CALL MEMGET('INTE',KIPIV,IDIM,WORK,KFREE,LFREE)
         CALL DGETRF(IDIM,IDIM,FMAT,IDIM,WORK(KIPIV),IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('LU factorization failed in COMP_MMRELMAT')
         END IF
         CALL DGETRI(IDIM,FMAT,IDIM,WORK(KIPIV),WORK(KFREE),LFREE,
     &               IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('Inversion failed in COMP_MMRELMAT')
         END IF
         CALL MEMREL('COMP_MMRELMAT',WORK,1,1,KFREE,LFREE)
         IF (IPRTLVL.GE.15) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Inverted Relay(MM) matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
         END IF
      END IF
C
      RETURN
      END
C  /* Deck get_amat */
      SUBROUTINE GET_AMAT(FMAT,IDIM)
C
C Purpose:
C     Computes A matrix component of Relay matrix for NP and
C     MM regions.
C
C Input:
C   IDIM - Dimension of Relay matrix
C Output:
C   FMAT - Relay matrix with A matrix contribution added up
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(IDIM,IDIM)
C
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0)
C
      DIMENSION RIJ(3)
C     Set diagonal components
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            RPOL = D1/NPFPOL(NPFTYP(I))
            DO J=1,3
               JOFF = (I-1)*3+J
               FMAT(JOFF,JOFF) = RPOL
            END DO
         END DO
      END IF
C     Set off-diagonal components
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            DO J=I+1,TNPATM
C              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RAD2 = RAD*RAD
               RAD51 = D1/(RAD2*RAD2*RAD)
               IF (NPMQGAU) CALL GET_GG_AFACT(FACTA,FACTB,I,J,RAD)
C              Distribute interaction tensor
               DO K=1,3
                  KOFF = (I-1)*3+K
                  DO L=1,3
                     LOFF = (J-1)*3+L
                     RVAL = D3*RIJ(K)*RIJ(L)
                     IF (K.EQ.L) RVAL = RVAL-RAD2
                     RVAL = RVAL*RAD51
                     IF (NPMQGAU) THEN
                        RVAL = RVAL*FACTA-FACTB*RIJ(K)*RIJ(L)
                     END IF
                     FMAT(KOFF,LOFF) = -RVAL
                     FMAT(LOFF,KOFF) = -RVAL
                  END DO
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck get_amat */
      SUBROUTINE GET_MM_AMAT(FMAT,IDIM)
C
C Purpose:
C     Computes A matrix component of Relay matrix for MM region.
C
C Input:
C   IDIM - Dimension of Relay matrix
C Output:
C   FMAT - Relay matrix with A matrix contribution added up
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(IDIM,IDIM)
C
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0)
C
      DIMENSION RIJ(3)
C     Set diagonal components
       IOFF = 1
       DO I=1,TMMATM
         IF (MMSKIP(I) .EQ. 0) CYCLE
         RPOL = D1/MMFPOL(MMFTYP(I))
         DO J=1,3
            JOFF = (IOFF-1)*3+J
            FMAT(JOFF,JOFF) = RPOL
         END DO
         IOFF = IOFF+1
      END DO
C     Set off-diagonal components
      IOFF = 1
      DO I=1,TMMATM
         IF (MMSKIP(I) .EQ. 0) CYCLE
         JOFF = IOFF+1
         DO J=I+1,TMMATM
            IF (MMSKIP(J) .EQ. 0) CYCLE
C           Compute distance dependent parameters
            RIJ(1) = MMCORD(1,I)-MMCORD(1,J)
            RIJ(2) = MMCORD(2,I)-MMCORD(2,J)
            RIJ(3) = MMCORD(3,I)-MMCORD(3,J)
            RAD = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
            RAD2 = RAD*RAD
            RAD51 = D1/(RAD2*RAD2*RAD)
C           Distribute interaction tensor
            DO K=1,3
               KOFF = (IOFF-1)*3+K
               DO L=1,3
                  LOFF = (JOFF-1)*3+L
                  RVAL = D3*RIJ(K)*RIJ(L)
                  IF (K.EQ.L) RVAL = RVAL-RAD2
                  RVAL = RVAL*RAD51
                  FMAT(KOFF,LOFF) = -RVAL
                  FMAT(LOFF,KOFF) = -RVAL
               END DO
            END DO
            JOFF = JOFF+1
         END DO
         IOFF = IOFF+1
      END DO
C
      RETURN
      END
C  /* Deck get_gg_afact */
      SUBROUTINE GET_GG_AFACT(FACTA,FACTB,IATM,JATM,RAD)
C
C Purpose:
C     Determines damping factors for T(2) operator for
c     Gaussian/Gaussian dipole model.
C
C Input:
C   IATM - I-th atom in NP region
C   JATM - J-th atom in NP region
C   RAD  - Distance between I and J atoms
C Output:
C   FACTA  - Scalling factor
C   FACTB  - Additional factor
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "pi.h"
#include "qmnpmm.h"
C
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0)
      PARAMETER (D13 = 1.0D0/3.0D0)
C     Get polarizabilities
      RIPOL = NPFPOL(NPFTYP(IATM))/D3
      RJPOL = NPFPOL(NPFTYP(JATM))/D3
C     Get damping radius
      RDIM = SQRT(D2)/SQRTPI
      RIM  = (RDIM*RIPOL)**D13
      RJM  = (RDIM*RJPOL)**D13
      RIJM = SQRT(RIM*RIM+RJM*RJM)
      RVAL = RAD/RIJM
C     Compute factors
      FACTA = DERF(RVAL)-D2*RVAL*DEXP(-RVAL*RVAL)/SQRTPI
      FACTB = D4*DEXP(-RVAL*RVAL)
      FACTB = FACTB/(RAD*RAD*RIJM*RIJM*RIJM*SQRTPI)
C
      RETURN
      END
C  /* Deck get_cmat */
      SUBROUTINE GET_CMAT(FMAT,IDIM)
C
C Purpose:
C     Computes C matrix component of Relay matrix for NP and
C     MM regions.
C
C Input:
C   IDIM   - Dimension of Relay matrix
C Output:
C   FMAT   - Relay matrix with C matrix contribution added up.
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(IDIM,IDIM)
C
      PARAMETER (D1 = 1.0D0)
C
      DIMENSION RIJ(3)
C     Set diagonal components
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
C        Fix me: MM shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            FMAT(IOFF,IOFF) = -D1/NPFCAP(NPFTYP(I))
         END DO
      END IF
C     Set off-diagonal components
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
C        Fix me: MM shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            DO J=I+1,TNPATM
               JOFF = ISTART+J
C              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD  = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RVAL = D1/RAD
               IF (NPMQGAU) THEN
                  CALL GET_GG_CFACT(FACT,I,J,RAD)
                  RVAL = FACT*RVAL
               END IF
               FMAT(IOFF,JOFF) = -RVAL
               FMAT(JOFF,IOFF) = -RVAL
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck get_gg_cfact */
      SUBROUTINE GET_GG_CFACT(FACT,IATM,JATM,RAD)
C
C Purpose:
C     Determines damping factors for T(0) operator for
c     Gaussian/Gaussian dipole model.
C
C Input:
C   IATM - I-th atom in NP region
C   JATM - J-th atom in NP region
C   RAD  - Distance between I and J atoms
C Output:
C   FACT  - Scalling factor
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "pi.h"
#include "qmnpmm.h"
C
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0)
      ! PARAMETER (D21 = SQRT(D2)) - sqrt not allowed in PARAMETR for PGI compilers

      D21 = SQRT(D2)
C     Get capacitancies
      RICAP = D21*NPFCAP(NPFTYP(IATM))/SQRTPI
      RJCAP = D21*NPFCAP(NPFTYP(JATM))/SQRTPI
C     Get damping radius & scalling factor
      RIJM = SQRT(RICAP*RICAP+RJCAP*RJCAP)
      FACT = DERF(RAD/RIJM)
C
      RETURN
      END
C  /* Deck get_mmat */
      SUBROUTINE GET_MMAT(FMAT,IDIM)
C
C Purpose:
C     Computes M matrix component of Relay matrix for NP and
C     MM regions.
C
C Input:
C   IDIM - Dimension of Relay matrix
C Output:
C   FMAT - Relay matrix with M matrix contribution added up
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(IDIM,IDIM)
C
      PARAMETER (D1 = 1.0D0)
C
      DIMENSION RIJ(3)
C     Set off-diagonal components
      IF (DONPPOL.AND.DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
         DO I=1,TNPATM
            DO J=1,TNPATM
               IF (I.EQ.J) CYCLE
C              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD  = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RAD3 = D1/(RAD*RAD*RAD)
C              Get damping factor
               IF (NPMQGAU) CALL GET_GG_MFACT(FACT,I,J,RAD)
C              Distribute contributions
               DO M=1,3
                  IOFF = (I-1)*3+M
                  JOFF = ISTART+J
                  RVAL = RIJ(M)*RAD3
                  IF (NPMQGAU) RVAL = RVAL*FACT
                  FMAT(IOFF,JOFF) = -RVAL
                  FMAT(JOFF,IOFF) = -RVAL
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck get_gg_mfact */
      SUBROUTINE GET_GG_MFACT(FACT,IATM,JATM,RAD)
C
C Purpose:
C     Determines damping factors for T(1) operator for
c     Gaussian/Gaussian dipole model.
C
C Input:
C   IATM - I-th atom in NP region
C   JATM - J-th atom in NP region
C   RAD  - Distance between I and J atoms
C Output:
C   FACT  - Scalling factor
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "pi.h"
#include "qmnpmm.h"
C
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0)
      ! PARAMETER (D21 = SQRT(D2)) - sqrt not allowed in PARAMETR for PGI compilers
      PARAMETER (D13 = 1.0D0/3.0D0)

      D21 = SQRT(D2)
C     Get polarizabilities
      RIPOL = NPFPOL(NPFTYP(IATM))/D3
C     Get damping radius
      RDIM = D21/SQRTPI
      RIM  = (RDIM*RIPOL)**D13
C     Get j-th capacitancy
      RJCAP = RDIM*NPFCAP(NPFTYP(JATM))
C     Get damping radius & scalling factor
      RIJM = SQRT(RIM*RIM+RJCAP*RJCAP)
      RADX = RAD/RIJM
      FACT = DERF(RADX)-D2*RADX*DEXP(-RADX*RADX)/SQRTPI
C
      RETURN
      END
C  /* Deck get_qlag */
      SUBROUTINE GET_QLAG(FMAT,IDIM)
C
C Purpose:
C     Sets charge constrain in Relay matrix for NP and
C     MM regions.
C
C Input:
C   IDIM   - Dimension of Relay matrix
C Output:
C   FMAT   - Relay matrix with C matrix contribution added up.
C
C Last updated: 22/03/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "qmnpmm.h"
C
      DIMENSION FMAT(IDIM,IDIM)
C
      PARAMETER (D1 = 1.0D0)
C     Set diagonal components
      IF (DONPCAP.OR.DOMMCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
C        Fix me: MM shift
         IF (DONPCAP) THEN
            DO I=1,TNPATM
               IOFF = ISTART+I
               FMAT(IOFF,IDIM) = D1
               FMAT(IDIM,IOFF) = D1
            END DO
         END IF
      END IF
C
      RETURN
      END
C  /* Deck write_relmat */
      SUBROUTINE WRITE_RELMAT(FMAT)
C
C Purpose:
C     Stores Relay matrix in binary file.
C
C Input:
C  FMAT   - Inverted Relay matrix
C
C Last updated: 16/08/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "inftap.h"
C
      DIMENSION FMAT(*)

C     determine dimensions
      CALL GETDIM_RELMAT(IDIM,.TRUE.)
C     write inverted Relay matrix
      LUQMNP = -1
      CALL GPOPEN(LUQMNP,'QMMMNP','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUQMNP)
      CALL WRTIEF(FMAT, IDIM, 'QQMNPMAT', LUQMNP)
      CALL GPCLOSE(LUQMNP,'KEEP')
C
      RETURN
      END
C  /* Deck read_relmat */
      SUBROUTINE READ_RELMAT(FMAT)
C
C Purpose:
C     Reads Relay matrix in binary file.
C
C Input:
C  FMAT   - Inverted Relay matrix
C
C Last updated: 16/08/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "inftap.h"
C
      DIMENSION FMAT(*)
C
      LOGICAL FNDLAB
C
C     determine dimensions
      CALL GETDIM_RELMAT(IDIM,.TRUE.)
C     read inverted Relay matrix
      LUQMNP = -1
      CALL GPOPEN(LUQMNP,'QMMMNP','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUQMNP)
      IF (FNDLAB('QQMNPMAT',LUQMNP)) THEN
        CALL READT(LUQMNP,IDIM,FMAT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMNP file.')
      ENDIF
      CALL GPCLOSE(LUQMNP,'KEEP')
C
      RETURN
      END
C  /* Deck comp_dampvmat */
      SUBROUTINE COMP_DAMPVMAT(FCAO, MQVEC, WORK, LWORK)
C
C Purpose:
C     Computes damped potential matrix in AO basis.
C
C Input:
C  FMAT   - Inverted Relay matrix
C
C Last updated: 16/08/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "dummy.h"
#include "maxorb.h"
#include "aovec.h"
#include "shells.h"
#include "primit.h"
C
      DOUBLE PRECISION MQVEC 
C
      DIMENSION FCAO(*), MQVEC(*), WORK(LWORK)
C
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0)
      PARAMETER (D13 = 1.0D0/3.0D0)
C
C     determine dimensions
      CALL GETDIM_RELMAT(IDIM,.FALSE.)
C
      KFREE = 1
      LFREE = LWORK
C     allocate and set gaussian broadening paramenters
      CALL MEMGET('REAL',KRDVEC,TNPATM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRQVEC,TNPATM,WORK,KFREE,LFREE)
      CALL SET_DAMPARAM(WORK(KRDVEC),WORK(KRQVEC))
      IF (IPRTLVL.GE.15) THEN 
         WRITE(LUPRI,'(/,2X,A)') '*** Computed MQ vector ***'
         CALL OUTPUT(MQVEC,1,IDIM,1,1,IDIM,1,1,LUPRI)
      END IF     
C
#ifdef VAR_MPI
           call mpixbcast (QMCMM_WORK  , 1, 'INTEGER', 0)
           iprint = 0
           call mpixbcast (iprint, 1, 'INTEGER', 0)
#endif

#ifdef ENABLE_VPOTDAMP
           call vpotdamped(kmax,nhkt,nuco,nrco,jstrt,cent,priccf,priexp,
     &                fcao,
     &                npcord,
     &                mqvec(3*tnpatm+1), work(krqvec),
     &                mqvec            , work(krdvec),
     &                tnpatm)
#else
            call quit('VPOTDAMP not compiled in this version')
#endif

      CALL MEMREL('COMP_DAMPVMAT',WORK,1,1,KFREE,LFREE)
C
      RETURN
      END
C  /* Deck set_damparam */
      SUBROUTINE SET_DAMPARAM(RDVEC, RQVEC)
C
C Purpose:
C     Sets damping parameters vectors for dipoles and charges.
C
C Input:
C  FMAT   - Inverted Relay matrix
C
C Last updated: 16/08/2013 by Z. Rinkevicius.
C
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "dummy.h"
#include "pi.h"
C
      DIMENSION RDVEC(*), RQVEC(*)
C
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (D13 = 1.0D0/3.0D0)
C
      RDIM = SQRT(D2)/SQRTPI
      DO I=1,TNPATM
         RIPOL = NPFPOL(NPFTYP(I))/D3
         RDVEC(I) = (RDIM*RIPOL)**D13
      END DO
C
      DO I=1,TNPATM
         RQVEC(I) = RDIM*NPFCAP(NPFTYP(I))
      END DO
C
      RETURN
      END

