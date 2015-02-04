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

module qmcmm

   public read_pot_qmnpmm
   public qmnpinp
   public getdim_relmat
   public getdim_mmmat
   public comp_relmat
   public read_relmat
   public write_relmat
   public comp_mmrelmat
   public comp_dampvmat

   private

contains

      SUBROUTINE QMNPINP(WORD)
!
! Purpose:
!   reads "*QMNPMM" input group and sets up various flags for
!   QM/NP/MM type embedding
!
! Input:
!   WORD - DALTON input group keyword
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "gnrinf.h"
#include "qmnpmm.h"

#ifdef VAR_MPI
#include "iprtyp.h"
#endif

!
      PARAMETER ( NTABLE = 13)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
!
      DATA TABLE /'.QMNP  ', '.QMNPMM', '.QMMM  ', '.NPPOIN',           &
             '.MMGAUS', '.PRINT ', '.MQITER', '.CMXPOL',                &
             '.NONPCA', '.MMPOLA', '.MMCAPA', '.XXXXXX',                &
             '.DAMPED'/
!
      CALL QENTER('QMNPINP')
!
      CALL SET_QMNPMM()
!
      NEWDEF = (WORD .EQ. '*QMNPMM')
      ICHANG = 0
!
      IF (NEWDEF) THEN
        WORD1 = WORD
!
  100   CONTINUE
!
        READ (LUCMD, '(A7)') WORD
        CALL UPCASE(WORD)
        PROMPT = WORD(1:1)
!
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
!
          WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,                     &
                            '" not recognized for *QMNPMM'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword for *QMNPMM')
!
!         ".QMNP "
    1     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .TRUE.
            DOMMSUB = .FALSE.
            DONPCAP = .TRUE.
            DONPPOL = .TRUE.
          GO TO 100
!
!         ".QMNPMM"
    2     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .TRUE.
            DOMMSUB = .TRUE.
            DONPCAP = .TRUE.
            DONPPOL = .TRUE.
          GO TO 100
!
!         ".QMMM  "
    3     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .FALSE.
            DOMMSUB = .TRUE.
          GO TO 100
!
!         ".NPPOIN"
    4     CONTINUE
            NPMQGAU = .FALSE.
          GO TO 100
!
!         ".MMGAUS"
    5     CONTINUE
            MMMQGAU = .TRUE.
          GO TO 100
!
!         ".PRINT "
    6     CONTINUE
            READ(LUCMD,*) IPRTLVL
          GO TO 100
!
!         ".MQITER"
    7     CONTINUE
            MQITER = .TRUE.
          GO TO 100
!
!         ".CMXPOL"
    8     CONTINUE
            CMXPOL = .TRUE.
          GO TO 100
!
!         ".NONPCA"
    9     CONTINUE
            DONPCAP = .FALSE.
          GO TO 100
!
!         ".MMPOLA"
   10     CONTINUE
            DOMMPOL = .TRUE.
          GO TO 100
!
!         ".MMCAPA"
   11     CONTINUE
            CMXPOL = .TRUE.
          GO TO 100
!
!         ".XXXXXX"
   12     CONTINUE
          GO TO 100
!         ".DAMPED"
 13       CONTINUE
            NOVDAMP = .FALSE.
          GO TO 100

!
        ELSE IF (PROMPT .EQ. '*') THEN
          GO TO 300
        ELSE
          WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,                      &
                              '" not recognized for *QMNPMM'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal prompt for *QMNPMM')
        END IF
      END IF
!
  300 CONTINUE
!
      CALL QEXIT('QMNPINP')
      END
      SUBROUTINE SET_QMNPMM()
!
! Purpose:
!   sets initial values of all parameters related to
!   QM/NP/MM type embedding
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
!
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
!
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
!
      ENSOLQNP = 0.0D0
      EESOLQNP = 0.0D0
      ENSOLMNP = 0.0D0
      EESOLMNP = 0.0D0
!
      ENSOLQMM = 0.0D0
      EESOLQMM = 0.0D0
      ENSOLMMM = 0.0D0
      EESOLMMM = 0.0D0
!
      END
      SUBROUTINE READ_POT_QMNPMM()
!
! Purpose:
!   reads NP and MM subsystems geometry, total charges, and
!   force field data from POTENTIAL.INP file
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      PARAMETER (XFACT = 1.0D0/XTANG)
!
      CHARACTER UNITS*2, NPWORD*2, FFWORD*4, NPLBL(MXNPATM)*2
      CHARACTER MMWORD*2, MMLBL(MXMMATM)*2
!
      CALL QENTER('READ_POT_QMNPMM')
!     Open POTENTIAL.INP file
      LUQMNP=-1
      CALL GPOPEN(LUQMNP,'POTENTIAL.INP','OLD',' ',                     &
             'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUQMNP)
!     Read geometry input units, number of NP and MM
!     separate subsystems
      READ(LUQMNP,*) UNITS, TNPBLK, TMMBLK
      CALL UPCASE(UNITS)
!     Check input consistency
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
        WRITE(LUPRI,'(/2X,A)')                                          &
         'Maximum number of NP subsystems exceeded!'
        WRITE(LUPRI,'(2X,A,I3,2X,A,I3)') 'Input TNPBLK=',               &
         TNPBLK, 'Maximum allowed:', MAXBLK
        CALL QUIT('Increase MAXBLK in qmnpmm.h')
      END IF
      IF (TMMBLK.GT.MAXBLK) THEN
        WRITE(LUPRI,'(/2X,A)')                                          &
         'Maximum number of MM subsystems exceeded!'
        WRITE(LUPRI,'(2X,A,I3,2X,A,I3)') 'Input TNPBLK=',               &
         TNPBLK, 'Maximum allowed:', MAXBLK
        CALL QUIT('Increase MAXBLK in qmnpmm.h')
      END IF
!     Read NP subsystems data
      IF (DONPSUB) THEN
         ISTART = 0
         DO I=1,TNPBLK
!           Read NP subsystem header
            READ(LUQMNP,*) NPWORD, NPCHRG(I), NPATOM(I)
!           Check input consistency
            CALL UPCASE(NPWORD)
            IF (NPWORD.NE.'NP') THEN
               WRITE(LUPRI,'(/2X,A,I2,1X,A)')                           &
         'Incorrect NP subsystem I=', I,                                &
         'header in POTENTIAL.INP file!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
            IF (NPATOM(I).LE.0) THEN
               WRITE(LUPRI,'(/2X,A,I5,A,I2,A)')                         &
         'Incorrect number of atoms ', NPATOM(I),                       &
         'in NP subsystem', I, '!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
!           Read NP subsystem data
            DO J=1,NPATOM(I)
               JOFF = ISTART + J
               READ(LUQMNP,*) NPLBL(JOFF), NPCORD(1,JOFF),              &
               NPCORD(2,JOFF), NPCORD(3,JOFF),                          &
               NPFTYP(JOFF)
            END DO
            ISTART = ISTART + NPATOM(I)
            TNPATM = TNPATM + NPATOM(I)
            IF (TNPATM.GT.MXNPATM) THEN
               WRITE(LUPRI,'(/2X,A)')                                   &
          'Maximum number of NP atoms exceeded!'
               WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')                         &
          'Input TNPATM=', TNPATM,                                      &
          'Maximum allowed:', MXNPATM
               CALL QUIT('Increase MXNPATM in qmnpmm.h')
            END IF
         END DO
!        Convert to atomic units if neeeded
         IF (UNITS.EQ.'AA') THEN
            CALL DSCAL(3*TNPATM,XFACT,NPCORD,1)
         END IF
!        Read force field data
         READ(LUQMNP,*) FFWORD, TNPFF
         CALL UPCASE(FFWORD)
!        Check NP force field data consistency
         IF (FFWORD.NE.'NPFF') THEN
            WRITE(LUPRI,'(/2X,A,A)') 'Corrupted NP force field ',       &
            'header in POTENTIAL.INP file!'
            CALL QUIT('Corrupted POTENTIAL.INP')
         END IF
         IF (TNPFF.GT.MXNPFF) THEN
            WRITE(LUPRI,'(/2X,A)')                                      &
             'Maximum number of NP force field types exceeded!'
            WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')                            &
             'Input TNPFF=', TNPFF,                                     &
             'Maximum allowed:', MXNFF
            CALL QUIT('Increase MXNPFF in qmnpmm.h')
         END IF
!        Read force field data
         DO I=1,TNPFF
            READ(LUQMNP,*) NPFPOL(I), NPFCAP(I), NPFOMG1(I),            &
                      NPFGAM1(I), NPFOMG2(I), NPFGAM2(I),               &
                      NPFFAC(I)
         END DO
!        Check force field definitions
         DO I=1,TNPATM
            IF (NPFTYP(I).GT.TNPFF) THEN
               WRITE(LUPRI,'(/2X,A,I4,1X,A)')                           &
          'Unknown force field requested for atom', I,                  &
          'in NP region!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
         END DO
         IF (IPRTLVL.GE.5) THEN
            CALL PRINT_NPREGION(NPLBL)
         END IF
      END IF
!     MM subsystems input
      IF (DOMMSUB) THEN
         ISTART = 0
         DO I=1,TMMBLK
!           Read MM subsystem header
            READ(LUQMNP,*) MMWORD,  MMATOM(I)
!           Check input consistency
            CALL UPCASE(MMWORD)
            IF (MMWORD.NE.'MM') THEN
               WRITE(LUPRI,'(/2X,A,I2,1X,A)')                           &
         'Incorrect MM subsystem I=', I,                                &
         'header in POTENTIAL.INP file!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
            IF (MMATOM(I).LE.0) THEN
               WRITE(LUPRI,'(/2X,A,I5,A,I2,A)')                         &
         'Incorrect number of atoms ', MMATOM(I),                       &
         'in MM subsystem', I, '!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
!           Read MM subsystem data
            DO J=1,MMATOM(I)
               JOFF = ISTART + J
               READ(LUQMNP,*) MMLBL(JOFF), MMMOL(JOFF),                 &
                         MMCORD(1,JOFF), MMCORD(2,JOFF),                &
                         MMCORD(3,JOFF), MMFTYP(JOFF)
            END DO
            ISTART = ISTART + MMATOM(I)
            TMMATM = TMMATM + MMATOM(I)
            IF (TMMATM.GT.MXMMATM) THEN
               WRITE(LUPRI,'(/2X,A)')                                   &
          'Maximum number of MM atoms exceeded!'
               WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')                         &
          'Input TMMATM=', TMMATM,                                      &
          'Maximum allowed:', MXMMATM
               CALL QUIT('Increase MXMMATM in qmnpmm.h')
            END IF
         END DO
!        Convert to atomic units if neeeded
         IF (UNITS.EQ.'AA') THEN
            CALL DSCAL(3*TMMATM,XFACT,MMCORD,1)
         END IF
!        Read force field data
         READ(LUQMNP,*) FFWORD, TMMFF
         CALL UPCASE(FFWORD)
!        Check NP force field data consistency
         IF (FFWORD.NE.'MMFF') THEN
            WRITE(LUPRI,'(/2X,A,A)') 'Corrupted MM force field ',       &
            'header in POTENTIAL.INP file!'
            CALL QUIT('Corrupted POTENTIAL.INP')
         END IF
!        Read force field data
         DO I=1,TMMFF
            IF ((.NOT.DOMMCAP).AND.(.NOT.DOMMPOL)) THEN
               READ(LUQMNP,*) MMFM0(I)
            END IF
            IF ((.NOT.DOMMCAP).AND.DOMMPOL) THEN
               READ(LUQMNP,*) MMFM0(I),MMFPOL(I)
            END IF
         END DO
!        Check force field definitions
         DO I=1,TMMATM
            IF (MMFTYP(I).GT.TMMFF) THEN
               WRITE(LUPRI,'(/2X,A,I4,1X,A)')                           &
          'Unknown force field requested for atom', I,                  &
          'in MM region!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
!           Set polarization centers
            IF (MMFPOL(MMFTYP(I)) .GT. 1.0D-6) THEN
               TPOLATM = TPOLATM+1
               MMSKIP(I) = 1
!              works only for one non-metallic MM region
            END IF
         END DO
         IF (IPRTLVL.GE.5) THEN
            CALL PRINT_MMREGION(MMLBL)
         END IF
      END IF
!     Close POTENTIAL.INP file
      CALL GPCLOSE(LUQMNP,'KEEP')
!
      CALL QEXIT('READ_POT_QMNPMM')
!
      END
      SUBROUTINE PRINT_NPREGION(ATMLBL)
!
! Purpose:
!     prints detailed information about nanoparticles
!     in QM/NP/MM embedding
!
! Input:
!  ATMLBL - List of atom labels for NP region
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      CHARACTER ATMLBL(MXNPATM)*2
!
      IF (DOMMSUB) THEN
         CALL TITLER('QM/NP/MM Embedding: Nanoparticle(s) data','*',103)
      ELSE
        CALL TITLER('QM/NP Embedding: Nanoparticle(s) data','*',103)
      END IF
!
      IF (TNPBLK.GT.1) THEN
         WRITE(LUPRI,'(/2X,A,I3,1X,A)') 'Nanoparticle region contains', &
          TNPBLK, 'separate nanoparticles.'
      ELSE
        WRITE(LUPRI,'(/2X,A,A)') 'Nanoparticle region contains single', &
         ' nanoparticle.'
      END IF
      ISTART = 0
      DO I=1,TNPBLK
         WRITE(LUPRI, '(/,2X,A,I3,1X,A,F6.3,A)') 'Nanoparticle',I,      &
          'with charge', NPCHRG(I), '. Coordinates in au.'
         WRITE(LUPRI, '(2X,A)')                                         &
          '--------------------------------------------------'
         WRITE(LUPRI, '(2X,A)')                                         &
          'Atom    Coord. X    Coord. Y    Coord. Z   FF Type'
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
         DO J=1,NPATOM(I)
            JOFF = ISTART+J
            WRITE(LUPRI,'(3X,A,1X,F11.5,1X,F11.5,1X,F11.5,3X,I4)')      &
             ATMLBL(JOFF), NPCORD(1,JOFF), NPCORD(2,JOFF),              &
             NPCORD(3,JOFF), NPFTYP(JOFF)
         END DO
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
         ISTART = ISTART + NPATOM(I)
      END DO
      WRITE(LUPRI,'(/2X,A,A)') 'Static force field(s) data for ',       &
       'nanoparticle region.'
      WRITE(LUPRI, '(2X,A)')                                            &
       '-------------------------------------'
      WRITE(LUPRI, '(2X,A)')                                            &
       'FF Type Polarizability    Capacitance'
      WRITE(LUPRI, '(2X,A)')                                            &
       '====================================='
      DO I=1,TNPFF
         WRITE(LUPRI,'(3X,I2,5X,F11.5,5X,F11.5)') I, NPFPOL(I),         &
          NPFCAP(I)
      END DO
      WRITE(LUPRI, '(2X,A)')                                            &
       '====================================='
      WRITE(LUPRI,'(/2X,A,A)') 'Dynamic force field(s) data for ',      &
       'nanoparticle region.'
      WRITE(LUPRI, '(2X,A,A)')                                          &
       '---------------------------------------------------',           &
       '------'
      WRITE(LUPRI, '(2X,A,A)')                                          &
       'FF Type Omega(1)  Gamma(1)  Omega(2)  Gamma(2)',                &
       '  Factor'
      WRITE(LUPRI, '(2X,A,A)')                                          &
       '===================================================',           &
       '======'
      DO I=1,TNPFF
         WRITE(LUPRI,'(3X,I2,3X,5(1X,F9.5))') I, NPFOMG1(I), NPFGAM1(I),&
          NPFOMG2(I), NPFGAM2(I), NPFFAC(I)
      END DO
      WRITE(LUPRI, '(2X,A,A)')                                          &
       '===================================================',           &
       '======'
!
      END
      SUBROUTINE PRINT_MMREGION(ATMLBL)
!
! Purpose:
!     prints detailed information about MM region
!     in QM/NP/MM embedding
!
! Input:
!  ATMLBL - List of atom labels for MM region
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      CHARACTER ATMLBL(MXNPATM)*2
!
      IF (DOMMSUB) THEN
         CALL TITLER('QM/NP/MM Embedding: MM region(s) data','*',103)
      END IF
!
      IF (TMMBLK.GT.1) THEN
         WRITE(LUPRI,'(/2X,A,I3,1X,A)') 'MM region contains',           &
          TMMBLK, 'separate non-metallic subregions.'
      END IF
      ISTART = 0
      DO I=1,TMMBLK
         WRITE(LUPRI, '(/,2X,A,I3,1X,A)') 'MM subregion',I,             &
          '. Coordinates in au.'
         WRITE(LUPRI, '(2X,A)')                                         &
          '--------------------------------------------------'
         WRITE(LUPRI, '(2X,A)')                                         &
          'Atom    Coord. X    Coord. Y    Coord. Z   FF Type'
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
         DO J=1,MMATOM(I)
            JOFF = ISTART+J
            WRITE(LUPRI,'(3X,A,1X,F11.5,1X,F11.5,1X,F11.5,3X,I4)')      &
             ATMLBL(JOFF), MMCORD(1,JOFF), MMCORD(2,JOFF),              &
             MMCORD(3,JOFF), MMFTYP(JOFF)
         END DO
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
!        works only for one MM region
         WRITE(LUPRI, '(2X,A,I4)') 'Number of polarizable centers   : ',&
          TPOLATM
         WRITE(LUPRI, '(2X,A,I4)') 'Number of nonpolarizable centers: ',&
          MMATOM(I)-TPOLATM
         ISTART = ISTART + NPATOM(I)
      END DO
      WRITE(LUPRI,'(/2X,A)') 'Force field(s) data for MM region'
!
      IF ((.NOT.DOMMCAP).AND.(.NOT.DOMMPOL)) THEN
         WRITE(LUPRI, '(2X,A)') '--------------------------------'
         WRITE(LUPRI, '(2X,A)') 'FF Type MM Charge Polarizability'
         WRITE(LUPRI, '(2X,A)') '================================'
         DO I=1,TMMFF
            WRITE(LUPRI,'(3X,I2,3X,1X,F9.5,4X,F9.5,1X,I2)') I, MMFM0(I),&
             MMFPOL(I), MMSKIP(I)
         END DO
         WRITE(LUPRI, '(2X,A)') '================================'
      END IF
      END
      SUBROUTINE GETDIM_RELMAT(IMATDIM, SQFLAG)
!
! Purpose:
!     determines Relay matrix dimensions.
!
! Input:
!   SQFLAG - Request total size of Relay matrix
! Output:
!  IMATDIM - Size of Relay matrix or it's dimension
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      LOGICAL SQFLAG
!
      IMATDIM = 0
!     Add nanoparticle contribution
      IF (DONPSUB) THEN
         IMATDIM = IMATDIM + 3*TNPATM
         IF (DONPCAP) THEN
            IMATDIM = IMATDIM + TNPATM
         END IF
      END IF
!     Add Lagrangian term
      IF (DONPCAP.OR.DOMMCAP) THEN
         IMATDIM = IMATDIM + 1
      END IF
!     Get requested dimmension
      IF ((.NOT.MQITER).AND.SQFLAG) THEN
         IMATDIM = IMATDIM*IMATDIM
      END IF
!
      END
      SUBROUTINE GETDIM_MMMAT(IMATDIM, SQFLAG)
!
! Purpose:
!     determines Relay matrix dimensions for MM region.
!
! Input:
!   SQFLAG - Request total size of Relay matrix
! Output:
!  IMATDIM - Size of Relay matrix or it's dimension
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      LOGICAL SQFLAG
!
      IMATDIM = 0
!     Add molecular contribution
      IF (DOMMSUB) THEN
         IMATDIM = IMATDIM + 3*TPOLATM
      END IF
!     Get requested dimmension
      IF ((.NOT.MQITER).AND.SQFLAG) THEN
         IMATDIM = IMATDIM*IMATDIM
      END IF
!
      END
      SUBROUTINE COMP_RELMAT(FMAT,WORK,LWORK)
!
! Purpose:
!     Computes Relay matrix and inverts it or estimates initial
!     MQ vector values for iterative MQ vector determination
!     algorithm.
!
! Input:
!   WORK  - Temporary memory array
!   LWORK - Size of temporary memory array
! Output:
! if .not. MQITER
!  FMAT   - Inverted Relay matrix
! else
!  FMAT   - Initial guess of real part of MQ vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(*), WORK(LWORK)
      integer, allocatable :: ipiv(:)
!
      KFREE = 1
      LFREE = LWORK
!     Initialize arrays
      CALL GETDIM_RELMAT(IDIM,.TRUE.)
      CALL DZERO(FMAT,IDIM)
!     Reset matrix dimension parameter
      CALL GETDIM_RELMAT(IDIM,.FALSE.)
      IF (.NOT.MQITER) THEN
!        Compute polarizabilty dependent terms
         IF (DONPPOL.OR.DOMMPOL) THEN
            CALL GET_AMAT(FMAT,IDIM)
         END IF
!        Compute capacitancy dependent term
         IF (DONPCAP.OR.DOMMCAP) THEN
           CALL GET_CMAT(FMAT,IDIM)
           CALL GET_MMAT(FMAT,IDIM)
           CALL GET_QLAG(FMAT,IDIM)
         END IF
      ELSE
!       Conjugated gradient method via paradiso solver
      END IF
!     Print Relay matrix
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed Relay matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
      END IF
!     Invert Relay matrix
      IF (.NOT.MQITER) THEN
         allocate(ipiv(idim))
         CALL DGETRF(IDIM,IDIM,FMAT,IDIM,IPIV,IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('LU factorization failed in COMP_RELMAT')
         END IF
         CALL DGETRI(IDIM,FMAT,IDIM,IPIV,WORK(KFREE),LFREE,             &
                IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('Inversion failed in COMP_RELMAT')
         END IF
         deallocate(ipiv)
         IF (IPRTLVL.GE.15) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Inverted Relay matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
         END IF
      END IF
!
      END
      SUBROUTINE COMP_MMRELMAT(FMAT,WORK,LWORK)
!
! Purpose:
!     Computes Relay matrix and inverts it or estimates initial
!     MQ vector values for iterative MQ vector determination
!     algorithm. (MM region)
!
! Input:
!   WORK  - Temporary memory array
!   LWORK - Size of temporary memory array
! Output:
! if .not. MQITER
!  FMAT   - Inverted Relay matrix
! else
!  FMAT   - Initial guess of real part of MQ vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(*), WORK(LWORK)
      integer, allocatable :: ipiv(:)
!
      KFREE = 1
      LFREE = LWORK
!     Initialize arrays
      CALL GETDIM_MMMAT(IDIM,.TRUE.)
      CALL DZERO(FMAT,IDIM)
!     Reset matrix dimension parameter
      CALL GETDIM_MMMAT(IDIM,.FALSE.)
      IF (.NOT.MQITER) THEN
!        Compute polarizabilty dependent terms
         CALL GET_MM_AMAT(FMAT,IDIM)
      ELSE
!       Conjugated gradient method via paradiso solver
      END IF
!     Print Relay matrix
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed Relay (MM) matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
      END IF
!     Invert Relay matrix
      IF (.NOT.MQITER) THEN
         allocate(ipiv(idim))
         CALL DGETRF(IDIM,IDIM,FMAT,IDIM,IPIV,IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('LU factorization failed in COMP_MMRELMAT')
         END IF
         CALL DGETRI(IDIM,FMAT,IDIM,IPIV,WORK(KFREE),LFREE,             &
                IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('Inversion failed in COMP_MMRELMAT')
         END IF
         deallocate(ipiv)
         IF (IPRTLVL.GE.15) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Inverted Relay(MM) matrix ***'
            CALL OUTPUT(FMAT,1,IDIM,1,IDIM,IDIM,IDIM,1,LUPRI)
         END IF
      END IF
!
      END
      SUBROUTINE GET_AMAT(FMAT,IDIM)
!
! Purpose:
!     Computes A matrix component of Relay matrix for NP and
!     MM regions.
!
! Input:
!   IDIM - Dimension of Relay matrix
! Output:
!   FMAT - Relay matrix with A matrix contribution added up
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(IDIM,IDIM)
!
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0)
!
      DIMENSION RIJ(3)
!     Set diagonal components
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            RPOL = D1/NPFPOL(NPFTYP(I))
            DO J=1,3
               JOFF = (I-1)*3+J
               FMAT(JOFF,JOFF) = RPOL
            END DO
         END DO
      END IF
!     Set off-diagonal components
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            DO J=I+1,TNPATM
!              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RAD2 = RAD*RAD
               RAD51 = D1/(RAD2*RAD2*RAD)
               IF (NPMQGAU) CALL GET_GG_AFACT(FACTA,FACTB,I,J,RAD)
!              Distribute interaction tensor
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
      END
      SUBROUTINE GET_MM_AMAT(FMAT,IDIM)
!
! Purpose:
!     Computes A matrix component of Relay matrix for MM region.
!
! Input:
!   IDIM - Dimension of Relay matrix
! Output:
!   FMAT - Relay matrix with A matrix contribution added up
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(IDIM,IDIM)
!
      PARAMETER (D1 = 1.0D0, D3 = 3.0D0)
!
      DIMENSION RIJ(3)
!     Set diagonal components
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
!     Set off-diagonal components
      IOFF = 1
      DO I=1,TMMATM
         IF (MMSKIP(I) .EQ. 0) CYCLE
         JOFF = IOFF+1
         DO J=I+1,TMMATM
            IF (MMSKIP(J) .EQ. 0) CYCLE
!           Compute distance dependent parameters
            RIJ(1) = MMCORD(1,I)-MMCORD(1,J)
            RIJ(2) = MMCORD(2,I)-MMCORD(2,J)
            RIJ(3) = MMCORD(3,I)-MMCORD(3,J)
            RAD = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
            RAD2 = RAD*RAD
            RAD51 = D1/(RAD2*RAD2*RAD)
!           Distribute interaction tensor
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
!
      END
      SUBROUTINE GET_GG_AFACT(FACTA,FACTB,IATM,JATM,RAD)
!
! Purpose:
!     Determines damping factors for T(2) operator for
!     Gaussian/Gaussian dipole model.
!
! Input:
!   IATM - I-th atom in NP region
!   JATM - J-th atom in NP region
!   RAD  - Distance between I and J atoms
! Output:
!   FACTA  - Scalling factor
!   FACTB  - Additional factor
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "pi.h"
#include "qmnpmm.h"
!
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0)
      PARAMETER (D13 = 1.0D0/3.0D0)
!     Get polarizabilities
      RIPOL = NPFPOL(NPFTYP(IATM))/D3
      RJPOL = NPFPOL(NPFTYP(JATM))/D3
!     Get damping radius
      RDIM = SQRT(D2)/SQRTPI
      RIM  = (RDIM*RIPOL)**D13
      RJM  = (RDIM*RJPOL)**D13
      RIJM = SQRT(RIM*RIM+RJM*RJM)
      RVAL = RAD/RIJM
!     Compute factors
      FACTA = DERF(RVAL)-D2*RVAL*DEXP(-RVAL*RVAL)/SQRTPI
      FACTB = D4*DEXP(-RVAL*RVAL)
      FACTB = FACTB/(RAD*RAD*RIJM*RIJM*RIJM*SQRTPI)
!
      END
      SUBROUTINE GET_CMAT(FMAT,IDIM)
!
! Purpose:
!     Computes C matrix component of Relay matrix for NP and
!     MM regions.
!
! Input:
!   IDIM   - Dimension of Relay matrix
! Output:
!   FMAT   - Relay matrix with C matrix contribution added up.
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(IDIM,IDIM)
!
      PARAMETER (D1 = 1.0D0)
!
      DIMENSION RIJ(3)
!     Set diagonal components
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
!        Fix me: MM shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            FMAT(IOFF,IOFF) = -D1/NPFCAP(NPFTYP(I))
         END DO
      END IF
!     Set off-diagonal components
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
!        Fix me: MM shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            DO J=I+1,TNPATM
               JOFF = ISTART+J
!              Compute distance dependent parameters
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
      END
      SUBROUTINE GET_GG_CFACT(FACT,IATM,JATM,RAD)
!
! Purpose:
!     Determines damping factors for T(0) operator for
!     Gaussian/Gaussian dipole model.
!
! Input:
!   IATM - I-th atom in NP region
!   JATM - J-th atom in NP region
!   RAD  - Distance between I and J atoms
! Output:
!   FACT  - Scalling factor
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "pi.h"
#include "qmnpmm.h"
!
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0, D4 = 4.0D0)
      PARAMETER (D21 = 1.41421356237309504880D0)
!     Get capacitancies
      RICAP = D21*NPFCAP(NPFTYP(IATM))/SQRTPI
      RJCAP = D21*NPFCAP(NPFTYP(JATM))/SQRTPI
!     Get damping radius & scalling factor
      RIJM = SQRT(RICAP*RICAP+RJCAP*RJCAP)
      FACT = DERF(RAD/RIJM)
!
      END
      SUBROUTINE GET_MMAT(FMAT,IDIM)
!
! Purpose:
!     Computes M matrix component of Relay matrix for NP and
!     MM regions.
!
! Input:
!   IDIM - Dimension of Relay matrix
! Output:
!   FMAT - Relay matrix with M matrix contribution added up
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(IDIM,IDIM)
!
      PARAMETER (D1 = 1.0D0)
!
      DIMENSION RIJ(3)
!     Set off-diagonal components
      IF (DONPPOL.AND.DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
         DO I=1,TNPATM
            DO J=1,TNPATM
               IF (I.EQ.J) CYCLE
!              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD  = SQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RAD3 = D1/(RAD*RAD*RAD)
!              Get damping factor
               IF (NPMQGAU) CALL GET_GG_MFACT(FACT,I,J,RAD)
!              Distribute contributions
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
      END
      SUBROUTINE GET_GG_MFACT(FACT,IATM,JATM,RAD)
!
! Purpose:
!     Determines damping factors for T(1) operator for
!     Gaussian/Gaussian dipole model.
!
! Input:
!   IATM - I-th atom in NP region
!   JATM - J-th atom in NP region
!   RAD  - Distance between I and J atoms
! Output:
!   FACT  - Scalling factor
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "pi.h"
#include "qmnpmm.h"
!
      PARAMETER (D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (D21 = 1.41421356237309504880D0)
      PARAMETER (D13 = 1.0D0/3.0D0)

!     Get polarizabilities
      RIPOL = NPFPOL(NPFTYP(IATM))/D3
!     Get damping radius
      RDIM = D21/SQRTPI
      RIM  = (RDIM*RIPOL)**D13
!     Get j-th capacitancy
      RJCAP = RDIM*NPFCAP(NPFTYP(JATM))
!     Get damping radius & scalling factor
      RIJM = SQRT(RIM*RIM+RJCAP*RJCAP)
      RADX = RAD/RIJM
      FACT = DERF(RADX)-D2*RADX*DEXP(-RADX*RADX)/SQRTPI
!
      END
      SUBROUTINE GET_QLAG(FMAT,IDIM)
!
! Purpose:
!     Sets charge constrain in Relay matrix for NP and
!     MM regions.
!
! Input:
!   IDIM   - Dimension of Relay matrix
! Output:
!   FMAT   - Relay matrix with C matrix contribution added up.
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      DIMENSION FMAT(IDIM,IDIM)
!
      PARAMETER (D1 = 1.0D0)
!     Set diagonal components
      IF (DONPCAP.OR.DOMMCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
!        Fix me: MM shift
         IF (DONPCAP) THEN
            DO I=1,TNPATM
               IOFF = ISTART+I
               FMAT(IOFF,IDIM) = D1
               FMAT(IDIM,IOFF) = D1
            END DO
         END IF
      END IF
!
      END
      SUBROUTINE WRITE_RELMAT(FMAT)
!
! Purpose:
!     Stores Relay matrix in binary file.
!
! Input:
!  FMAT   - Inverted Relay matrix
!
! Last updated: 16/08/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "inftap.h"
!
      DIMENSION FMAT(*)

!     determine dimensions
      CALL GETDIM_RELMAT(IDIM,.TRUE.)
!     write inverted Relay matrix
      LUQMNP = -1
      CALL GPOPEN(LUQMNP,'QMMMNP','UNKNOWN','SEQUENTIAL','UNFORMATTED', &
             IDUMMY,.FALSE.)
      REWIND(LUQMNP)
      CALL WRTIEF(FMAT, IDIM, 'QQMNPMAT', LUQMNP)
      CALL GPCLOSE(LUQMNP,'KEEP')
!
      END
      SUBROUTINE READ_RELMAT(FMAT)
!
! Purpose:
!     Reads Relay matrix in binary file.
!
! Input:
!  FMAT   - Inverted Relay matrix
!
! Last updated: 16/08/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "inftap.h"
!
      DIMENSION FMAT(*)
!
      LOGICAL FNDLAB
!
!     determine dimensions
      CALL GETDIM_RELMAT(IDIM,.TRUE.)
!     read inverted Relay matrix
      LUQMNP = -1
      CALL GPOPEN(LUQMNP,'QMMMNP','UNKNOWN','SEQUENTIAL','UNFORMATTED', &
             IDUMMY,.FALSE.)
      REWIND(LUQMNP)
      IF (FNDLAB('QQMNPMAT',LUQMNP)) THEN
        CALL READT(LUQMNP,IDIM,FMAT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMNP file.')
      ENDIF
      CALL GPCLOSE(LUQMNP,'KEEP')
!
      END




   subroutine comp_dampvmat(fmat, mqvec)
      ! computes damped potential matrix in ao basis

      implicit none

      real(8), intent(in)    :: fmat(*) ! inverted relay matrix
      real(8), intent(inout) :: mqvec(*)

#include "priunit.h"
#include "qmnpmm.h"
#include "maxorb.h"
#include "aovec.h"
#include "shells.h"
#include "primit.h"

      real(8), allocatable :: rdvec(:)
      real(8), allocatable :: rqvec(:)
      integer              :: idimension
      integer              :: iprint

      call getdim_relmat(idimension, .false.)

      ! allocate and set gaussian broadening paramenters

      allocate(rdvec(tnpatm))
      allocate(rqvec(tnpatm))

      call set_damparam(rdvec, rqvec)

      if (iprtlvl > 14) then
         write(lupri, '(/,2x,a)') '*** Computed MQ vector ***'
         call output(mqvec, 1, idimension, 1, 1, idimension, 1, 1, lupri)
      end if

#ifdef VAR_MPI
      call mpixbcast(QMCMM_WORK, 1, 'INTEGER', 0)
      iprint = 0
      call mpixbcast(iprint, 1, 'INTEGER', 0)
#endif

#ifdef ENABLE_VPOTDAMP
      call vpotdamped(kmax, nhkt, nuco, nrco, jstrt, cent, priccf, priexp, &
                      fmat,                                                &
                      npcord,                                              &
                      mqvec(3*tnpatm+1), rqvec,                            &
                      mqvec            , rdvec,                            &
                      tnpatm)
#else
      call quit('VPOTDAMP not compiled in this version')
#endif

      deallocate(rdvec)
      deallocate(rqvec)

   end subroutine


   pure subroutine set_damparam(rdvec, rqvec)
      ! sets damping parameters vectors for dipoles and charges
      ! see Eqs 14 and 15 in J. Phys. Chem. C 112, 40 (2008)

      implicit none

      real(8), intent(inout) :: rdvec(*)
      real(8), intent(inout) :: rqvec(*)

! tnpatm, npfpol, npfcap
#include "qmnpmm.h"

! sqrtpi
#include "pi.h"

      real(8) :: rdim
      real(8) :: ripol
      integer :: i

      rdim = sqrt(2.0d0)/sqrtpi
      do i = 1, tnpatm
         ripol = npfpol(npftyp(i))/3.0d0
         rdvec(i) = (rdim*ripol)**(1.0d0/3.0d0)
         rqvec(i) = rdim*npfcap(npftyp(i))
      end do

   end subroutine

end module
