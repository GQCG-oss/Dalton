!
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
!
!
! FILE: sirius/sirstex.F
!
C#######################################################################
       SUBROUTINE STXINP(WORD,INPERR,ALLOPT)
C
C Copyright Nov 1990 Hans Joergen Aa. Jensen
C
#include "implicit.h"
      PARAMETER (NTABLE = 10)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12

      LOGICAL ALLOPT
C
#include "stex.h"
#include "inftap.h"
#include "inforb.h"
#include "priunit.h"
C
      DATA TABLE /'.XAS   ','.XES   ','.SHAKE ','.PRINT ','.OPEN S',    ! 5
     &            '.COEFFI','.AO    ','.AUGER ','.AUGERT','.XXXXXX'/    !10
C
C Initialize
C
      XAS    = .FALSE.
      XES    = .FALSE.
      SHAKE  = .FALSE.
      STEXMO = .TRUE.
      AUGER  = .FALSE.
      RPATST = .FALSE.
      IPRSTX = 0
      CALL IZERO(NHLSTX,8)
      CALL IZERO(IHLSTX,NOPMAX*8)
      CALL DZERO(CJSTX,NOPMAX*8)
      CALL DZERO(CKSTX,NOPMAX*8)
C
C     ***** PROCESS INPUT FOR STXINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &   LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,                           ! 5
     &                   106,107,108,109,100), II                       !10
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI
     &            )
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword  for "*STEX  "')
C
C *** Option 1 >XAS<   X-ray absorption spectroscopy
  101 CONTINUE
      XAS    = .TRUE.
      GO TO 100
C *** Option 2 >XES<   X-ray emission spectrosopy 
  102 CONTINUE
      XES    = .TRUE.
      GO TO 100
C *** Option 3 >SHAKE< Electron shake-up/off spectroscopy
  103 CONTINUE
      SHAKE  = .TRUE.
      GO TO 100
C *** Option 4 >PRINT< Print level in stex routines
  104 CONTINUE
      READ(LUCMD,*)  IPRSTX
      GO TO 100
C *** Option 5 >OPEN S< Open shells in stex calculation
  105 CONTINUE
      READ(LUCMD,*) (NHLSTX(ISYM), ISYM=1,NSYM)
      READ(LUCMD,*) ((IHLSTX(I,ISYM),I=1,NHLSTX(ISYM)),ISYM=1,NSYM)
      GO TO 100
C *** Option 6 >COEFF< Coefficients for modified Fock matrices
  106 CONTINUE
      READ(LUCMD,*) ((CJSTX(I,ISYM),I=1,NHLSTX(ISYM)),ISYM=1,NSYM)
      READ(LUCMD,*) ((CKSTX(I,ISYM),I=1,NHLSTX(ISYM)),ISYM=1,NSYM)
      GO TO 100
C *** Option 7 >AO   < Save stex matrices in AO basis
  107 CONTINUE
      STEXMO = .FALSE.
      GO TO 100
C *** Option 8 >AUGER< Auger orbital pairs
  108 CONTINUE
      AUGER = .TRUE.
      READ(LUCMD,*) NAUGER,(IAUGER(I),JAUGER(I),I=1,NAUGER)
      GO TO 100
C *** Option 9 >AUGERTEST   < Build RPA matrix
  109 CONTINUE
      AUGER  = .TRUE.
      RPATST = .TRUE.
      GO TO 100
C
C
      END
C#######################################################################
      SUBROUTINE STXCTL(WORK,LWORK)
#include "implicit.h"
      DIMENSION WORK(LWORK)
C
C  Driver routine for static exchange calculations
C
C Global
C    NCMOT
C
#include "stex.h"
#include "inforb.h"
#include "infpri.h"
#include "priunit.h"
C
C Local
C
C
      CALL QENTER('STXCTL')
      CALL TIMER('START ',STXSTA,STXEND)
      CALL TITLER('Static exchange calculation','*',105)
      IF (XAS) CALL HEADER('X-ray absorption spectrum',lupri)
      IF (XES) CALL HEADER('X-ray emission spectrum',lupri)
      IF (SHAKE) CALL HEADER('Shake-up/shake-off spectrum',lupri)
      IF (AUGER) THEN
         CALL AUGCTL(WORK,LWORK)
         GO TO 9000
      END IF
      IF (IPRSTX.GE.0) THEN
         WRITE(LUPRI,'(7X,A,I5)') 
     &   'Print level in STEX routines ', IPRSTX
         DO ISYM=1,NSYM
            IF (NHLSTX(ISYM).GT.0) THEN
               WRITE(LUPRI,'(7X,A)')
     &            'Static exchange Hamiltonian calculated for'
               WRITE(LUPRI,'(7X,A,10(I4,A,I1,A))')
     &            'orbital (symmetry) ',
     &            (IHLSTX(I,ISYM), ' (',ISYM,')', I=1,NHLSTX(ISYM))
               WRITE(LUPRI,'(7X,A,10F4.1)')
     &            'J coefficients     ', 
     &            (CJSTX(I,ISYM),I=1,NHLSTX(ISYM))
               WRITE(LUPRI,'(7X,A,10F4.1)')
     &            'K coefficients     ', 
     &            (CKSTX(I,ISYM),I=1,NHLSTX(ISYM))
            END IF
         END DO
      END IF
      CALL FLSHFO(LUPRI)
C
C  Allocate memory
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET2('REAL','CMO' ,KCMO ,NCMOT   ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DEN' ,KDEN ,2*N2BASX,WORK,KFREE,LFREE) !maybe a multiple of
      CALL MEMGET2('REAL','FSTX',KFSTX,2*N2BASX,WORK,KFREE,LFREE) !n2basx
C
C  Augment basis set
C
      CALL STXAUG(WORK(KCMO),WORK,KFREE,LFREE)
C
C  Build modified densities
C
      CALL STXDEN(WORK(KCMO),WORK(KDEN),WORK,KFREE,LFREE)
C
C  Build STEX Fock matrix
C
      CALL STXFCK(WORK(KDEN),WORK(KFSTX),WORK,KFREE,LFREE)
C
C  Process Fock matrix, calculate primitive spectrum
C
C     CALL STXEIG()
C
C  Stietjes imageing for continuous spectrum
C
C     CALL STXSTM()
C
      CALL TIMER('STXCTL',STXSTA,STXEND)
      CALL MEMREL('STXCTL',WORK,1,1,KFREE,LFREE)
 9000 CALL QEXIT('STXCTL')
      END
C
C#######################################################################
      SUBROUTINE STXAUG(CMO,WORK,KFREE,LFREE)
#include "implicit.h"
C
C   Input: KFREE, LFREE  start and length of free workspace
C
C
C   Output: CMO, molecular orbitals defined on a smaller basis projected 
C           on the current (augmented), from STXAUGORB
C           Note: the CMO from the last call to STXAUGORB is returned
C
      INTEGER KFREE, LFREE
      DIMENSION CMO(*), WORK(*)
C
#include "inforb.h"
#include "inftap.h"      
#include "stex.h"
#include "priunit.h"
C
      CHARACTER*8 FILE 
C
      CALL QENTER('STXAUG')
C
C  Prepare MO interface file
C
C  1. # basis functions, # occupied
C  2. # initial state energy, nuclear energy
C  3. # initial state occupied orbitals
C  4. # initial state occupied orbital energies
C  5. # initial state occupied orbital indices
C  6. # initial state energy, nuclear energy
C  7. # initial state occupied orbitals
C  8. # initial state occupied orbital energies
C  9. # initial state occupied orbital indices
C 10. # MO overlap matrix
C
      LUSTEX = -1
      CALL GPOPEN(LUSTEX,'UNIT2','NEW','SEQUENTIAL','UNFORMATTED',0,
     &   .FALSE.)
      REWIND LUSTEX
      IF (STEXMO) THEN
         WRITE(LUSTEX) NORBT, NOCCT
      ELSE
         WRITE(LUSTEX) NBAST, NOCCT
      END IF
C
C Augment ground state orbitals
C
      FILE = 'SIRIFC1'
cvc
c ground state -> icall=0
      icall=0
      CALL STXAUGORB(FILE,EMC,ENUC,CMO,WORK,KFREE,LFREE,icall)
      IF (IPRSTX.GT.1) THEN
         CALL HEADER('Augmented initial state',LUPRI)
         DO ISYM=1,NSYM
            ICMOI=ICMO(ISYM)+1
            NBASI=NBAS(ISYM)
            NOCCI=NOCC(ISYM)
            CALL OUTPUT(CMO(ICMOI),1,NBASI,1,NOCCI,NBASI,NOCCI,1,LUPRI)
         END DO
      END IF
cvc
C
C Append occupied gs orbitals to file LUSTEX
C
      CALL STXAPPOCC(CMO,EMC,ENUC,LUSTEX,WORK,KFREE,LFREE)
C
C Append gs orbital energies to file LUSTEX
C
      CALL STXAPPENE(FILE,LUSTEX,WORK,KFREE,LFREE)
C
C Augment final state orbitals
C
      FILE = 'SIRIFC2'
cvc
c final ionic state -> icall=1
      icall=1
      CALL STXAUGORB(FILE,EMC,ENUC,CMO,WORK,KFREE,LFREE,icall)
      IF (IPRSTX.GT.1) THEN
         CALL HEADER('Augmented final state',LUPRI)
         DO ISYM=1,NSYM
            ICMOI=ICMO(ISYM)+1
            NBASI=NBAS(ISYM)
            NOCCI=NOCC(ISYM)
            CALL OUTPUT(CMO(ICMOI),1,NBASI,1,NOCCI,NBASI,NOCCI,1,LUPRI)
         END DO
      END IF
C
C Append final state occupied orbitals to file LUSTEX
C
      CALL STXAPPOCC(CMO,EMC,ENUC,LUSTEX,WORK,KFREE,LFREE)
C
C Append occupied final orbital energies to file LUSTEX
C
      CALL STXAPPENE(FILE,LUSTEX,WORK,KFREE,LFREE)
C
C Append molecular overlap to LUSTEX
C
      CALL STXMOVLAP(LUSTEX,WORK,KFREE,LFREE)
      CALL GPCLOSE(LUSTEX,'KEEP')
      CALL QEXIT('STXAUG')
      END
C#######################################################################
      SUBROUTINE STXAUGORB(FILE,EMC,ENUC,CMO,WORK,KFREE,LFREE,icall)
#include "implicit.h"
C
C   Input: KFREE, LFREE  start and length of free workspace
C
C
C   Output: CMO, molecular orbitals defined on a smaller basis projected 
C           on the current (augmented), from STXAUGORB
C          
C
      INTEGER KFREE, LFREE
      DIMENSION CMO(*), WORK(*)
      CHARACTER FILE*(*)
C
C  Global variables
C  LUSIFC, LBSIFC
C
#include "stex.h"
#include "inftap.h"
#include "inforb.h"
#include "infpri.h"
#include "maxorb.h"
#include "infinp.h"
#include "priunit.h"
C
C  Local variables
C
      DIMENSION MULD2H_(8,8) ,
     &   NFRO_(8), NRHF_(8), NISH_(8), NASH_(8), NORB_(8), NBAS_(8)
      LOGICAL OLFLAG
      PARAMETER (D1=1.0D0)
      LOGICAL DEBUG
      CHARACTER*8 LABEL
C
C  External functions
C
      LOGICAL FNDLAB
      EXTERNAL FNDLAB
C
C--------------------------------------------------------
C
C  Find old orbitals and dimensions from interface file
C
C
      CALL QENTER('STXAUGORB')
      DEBUG=IPRSTX.GT.10
      LUSIFC=-1
      CALL GPOPEN(LUSIFC,FILE,'OLD','SEQUENTIAL','UNFORMATTED',
     &     0,.FALSE.)
      REWIND LUSIFC
      IF (FNDLAB(LBSIFC,LUSIFC) ) THEN
         READ (LUSIFC) POTNUC_, EMY_, EACTIV_, EMCSCF_
         READ (LUSIFC) 
     *      NISHT_,NASHT_,NOCCT_,NORBT_,NBAST_,NCONF_,NWOPT_,NWOPH_,
     *      NCDETS_,NCMOT_,NNASHX_,NNASHY_,NNORBT_,N2ORBT_,
     *      NSYM_,MULD2H_, NRHF_,NFRO_,
     *      NISH_,NASH_,NORB_,NBAS_
         CALL MEMGET2('REAL','CMO',KCMO,NCMOT_,WORK,KFREE,LFREE)
         CALL READT(LUSIFC,NCMOT_,WORK(KCMO))
      ELSE
         CALL QUIT('ERROR (STXAUGORB) MO''s not found on SIRIFC')
      END IF
      CALL GPCLOSE (LUSIFC,'KEEP')
C
C  Return mc and nuclear energies
C
      EMC=EMCSCF_
      ENUC=POTNUC_
C
C   Consider trial mo:s from H1DIAG for spanning virtual space
C
      REWIND LUIT1
      LABEL='NEWORB  '
      if (debug) 
     &   write(lupri,*)'augorb(trial)fndlab label luit1',label,luit1
      IF (FNDLAB(LABEL,LUIT1)) THEN
         CALL READT(LUIT1,NCMOT,CMO)
         if (debug) call prorb(cmo,.false.,lupri)
      ELSE
         CALL QUIT('STXAUGORB: LABEL "'//LABEL//'" NOT FOUND')
      END IF
C
C   Need AO overlap matrix
C
      CALL MEMGET2('REAL','AOTR',KAOTR,NNBAST,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','AOSQ',KAOSQ,N2BAS(1),WORK,KFREE,LFREE)
      CALL RDONEL('OVERLAP ',.TRUE.,WORK(KAOTR),NNBAST)
C
C  Copy old orbitals to CMO with new dimensions
C
      IAO = 0
      IMO = 0
      IAO_ = 0
      IMO_ = 0
      IOFF = 0
      IOFF_ = 0
      DO ISYM=1, NSYM
         NAO = NBAS(ISYM)
         NAO_ = NBAS_(ISYM)
         NMO = NORB(ISYM)
         NMO_ = NORB_(ISYM)
         NLEN = NAO*NMO
         NLEN_ = NAO_*NMO_
C
C   Unitary transformation to bring trial to actual
C
         CALL DSPTGE(NAO,WORK(KAOTR+IIBAS(ISYM)),WORK(KAOSQ))
         CALL UAUG(
     &      WORK(KCMO+IOFF_),NAO_,NMO_,
     &      CMO(1+IOFF),NAO,NMO,
     &      WORK(KAOSQ),
     &      WORK(KFREE),LFREE
     &      )
         CALL DZERO(CMO(1+IOFF),NAO*NMO_)
         CALL MCOPY(NAO_,NMO_,WORK(KCMO+IOFF_),NAO_,CMO(1+IOFF),NAO)
         IF (DEBUG) THEN
            call STXoutful(cmo(1+ioff),nao,nmo,'add orbitals')
         END IF
         IOFF = IOFF + NAO*NMO
         IOFF_ = IOFF_+ NAO_*NMO_
         IF (DEBUG) THEN
            print *,'isym',isym
            print *,'nao ',nao,nao_
            print *,'nmo ',nmo,nmo_
            print *,'nlen',nlen,nlen_
            print *,'ioff',ioff,ioff_
         END IF
      END DO
C
C Print orbitals
C
      IF (IPRSTX.GT.10) THEN
CVC         CALL HEADER('Augmented orbitals from '//FILE,LUPRI)
         CALL HEADER('Augmented orbitals from ... ',LUPRI)
         CALL PRORB(CMO,.FALSE.,LUPRI)
         CALL FLSHFO(LUPRI)
      END IF
C
C Orthonormalize new orbitals using Gram-Schmidt
C
      OLFLAG = FLAG(32)
      FLAG(32) = .FALSE.
      CALL MEMGET2('REAL','S',KS,N2BASX,WORK,KFREE,LFREE)
      CALL ORTHO(CMO,WORK(KS),WORK(KFREE),LFREE)
      FLAG(32) = OLFLAG
      IF (IPRSTX.GT.10) THEN
CVC         CALL HEADER('Renormalized orbitals from '//FILE,LUPRI)
         CALL HEADER('Renormalized orbitals from ... ',LUPRI)
         CALL PRORB(CMO,.FALSE.,LUPRI)
         CALL FLSHFO(LUPRI)
      END IF
c
cvc
      if(icall.eq.0) then

c  scrive su aug.mo gli orbitali molecolari dello stato fondamentale
C  su base aumentata
c
        luaumo=-1
        CALL GPOPEN(luaumo,'AUMO','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &  0,.FALSE.)
!D      write(lupri,*)' write on ', luaumo,
!D   &            ' the gs MOs on the augmented ', 
!D   &            'basis set'
        CALL Waugmo(CMO,luaumo)
        CALL FLSHFO(luaumo)
        CALL GPCLOSE (luaumo,'KEEP')
!D      write(lupri,*)' UNIT ',luaumo,' is closed'
      endif
cvc
C
C Save new orbitals on file
C
      CALL APPORB(FILE,CMO,.FALSE.)
      CALL MEMREL('STXAUGORB',WORK,KCMO,KCMO,KFREE,LFREE)
      CALL QEXIT('STXAUGORB')
      END
C#######################################################################
      SUBROUTINE Waugmo(CMO,IOUT)
C
C Written 8-Sep-1998 by VC
C
C Purpose:
C  Write augmented molecular orbital coefficients etc. on unit IOUT.
C
C Input:
C  CMO, MO orbital coefficients (symmetry blocked)
C  IOUT, output file unit
C
#include "implicit.h"
      DIMENSION CMO(*)
C
C Used from common blocks:
C   INFINP : CENT,TYPE,SUPSYM,?
C   INFORB : NSYM,...
C   INFIND : ISSMO(),
C
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
#include "priunit.h"
C
      rewind iout
      write(iout) nsym
!D    write(lupri,*) ' nsym ',nsym
      write(iout) (norb(i),i=1,nsym)
!D    write(lupri,*) ' norb ',(norb(i),i=1,nsym)
      write(iout) (nbas(i),i=1,nsym)
!D    write(lupri,*) ' nbas ',(nbas(i),i=1,nsym)
      nbast=0
      ncmot=0
      do i=1,nsym
        nbast=nbast+nbas(i)
        ncmot=ncmot+norb(i)*nbas(i)
      enddo
!D    write(lupri,*) ' nbast,ncmot ', nbast,ncmot
      write(iout) (cent(i),i=1,nbast)
!D    write(lupri,*) ' vector cent '
      write(iout) (type(i),i=1,nbast)
!D    write(lupri,*) ' vector type '
      write(iout) (cmo(i),i=1,ncmot)
!D    write(lupri,*) ' vector cmo '
c     write(lupri,1) (cmo(i),i=1,(nbas(1)*5+1),nbas(1))
c   1 format(5f12.5)
C
C *** End of subroutine Waugmo
C
      RETURN
      END
C#######################################################################
      SUBROUTINE STXAPPOCC(CMO,EMC,ENUC,LUSTEX,WORK,KFREE,LFREE)
#include "implicit.h"
      REAL*8 CMO(*), WORK(*)
      INTEGER KFREE,LFREE
#include "inforb.h"
#include "stex.h"
      CALL QENTER('STXAPPOCC')
C
      WRITE(LUSTEX) EMC, ENUC
C
C Extract occupied orbitals
C
      CALL MEMGET2('REAL','MOOCC',KMOOCC,NBAST*NOCCT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KMOOCC),NBAST*NOCCT)
      IF (IPRSTX.GT.10) THEN
         call header('STXAPPOCC input:occupied orbitals  ',LUPRI)
         call prorb(cmo,.true.,LUPRI)
      END IF
      JOFF = 0
      DO ISYM=1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NOCC(ISYM)
         NORBI = NORB(ISYM)
         ICMOI = ICMO(ISYM)
         CALL MCOPY(NBASI,NOCCI,CMO(ICMOI+1),NBASI,
     &                          WORK(KMOOCC+JOFF),NBAST)
         JOFF = JOFF + NBAST*NOCC(ISYM) + NBAS(ISYM)
      END DO
      CALL MEMCHK('STXAPPOCC:after MCOPY loop',WORK,KMOOCC)
      IF (IPRSTX.GT.10) THEN
         call STXoutful(
     &      work(kmoocc),nbast,nocct,'STXAPPOCC:Extracted occupied'
     &      )
      END IF
C
C Transform orbitals to SIRIFC1 basis
C
      IF (STEXMO) THEN
         CALL LTRA('SIRIFC1 ',WORK(KMOOCC),WORK(KFREE),LFREE)
         NDIMT=NORBT
      ELSE
         NDIMT=NBAST
      END IF
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(WORK(KMOOCC),NDIMT,NOCCT,
     &               'STXAPPOCC:Orbitals written to LUSTEX')
      END IF
      CALL WRITT(LUSTEX,NDIMT*NOCCT,WORK(KMOOCC))
      CALL MEMREL('STXAPPOCC',WORK,KMOOCC,KMOOCC,KFREE,LFREE)
      CALL QEXIT('STXAPPOCC')
      END
C#######################################################################
      SUBROUTINE STXAPPENE(FILE,LUSTEX,WORK,KFREE,LFREE)
#include "implicit.h"
      REAL*8 WORK(*)
      INTEGER LUSTEX,KFREE,LFREE
      CHARACTER FILE*(*)
#include "inforb.h"
#include "inftap.h"
#include "stex.h"
#include "maxorb.h"
#include "maxash.h"
#include "infind.h"
      DIMENSION MULD2H_(8,8) ,
     &   NFRO_(8), NRHF_(8), NISH_(8), NASH_(8), NORB_(8), NBAS_(8)
      LOGICAL FNDLAB, DEBUG
      EXTERNAL FNDLAB
      PARAMETER (DH = 0.5D0)
      CALL QENTER('STXAPPENE')
      DEBUG=IPRSTX.GT.10
      IF (DEBUG) PRINT *, '--- STXAPPENE called'
C
C Read orbital energies from SIRIFC for open-shell case
C we have to read FC, FV and diagonalize FC + 2*FV
C
      CALL GPOPEN(LUSIFC,FILE,'OLD','SEQUENTIAL','UNFORMATTED',
     &     0,.FALSE.)
      REWIND LUSIFC
      IF (FNDLAB(LBSIFC,LUSIFC) ) THEN
         READ (LUSIFC) ! 1. E...
         READ (LUSIFC) ! 2. NASHT...
     *      NISHT_,NASHT_,NOCCT_,NORBT_,NBAST_,NCONF_,NWOPT_,NWOPH_,
     *      NCDETS_,NCMOT_,NNASHX_,NNASHY_,NNORBT_,N2ORBT_,
     *      NSYM_,MULD2H_, NRHF_,NFRO_,
     *      NISH_,NASH_,NORB_,NBAS_,
     *      NELMN1_, NELMX1_, NELMN3_, NELMX3_, MCTYPE_,
     *      NAS1_, NAS2_, NAS3_
         READ (LUSIFC) ! 3. CMO
         READ (LUSIFC) ! 4. REF
         READ (LUSIFC) ! 5. DV
         CALL MEMGET2('REAL','F',KF,N2ORBT_,WORK,KFREE,LFREE)
         CALL READT(LUSIFC,N2ORBT_,WORK(KF))
      END IF
      CALL GPCLOSE(LUSIFC,'KEEP')
      IF (IPRSTX.GT.10) THEN
         CALL HEADER('STXAPPENE:F',LUPRI)
         CALL OUTPTB(WORK(KF),NORB_,NSYM_,1,LUPRI)
      END IF
      CALL MEMGET2('REAL','EVEC',KEVEC,NORB_(1)*NBAS_(1),
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVAL',KEVAL,NORB_(1)*NBAS_(1),
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SCR ',KSCR,NORB_(1),WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','ISCR',ISCR,NORB_(1),WORK,KFREE,LFREE)
C
C Extract Fock matrix diagonal
C
      IOFF  = 0
      IIOFF = 0
      I2OFF = 0
      DO ISYM=1,NSYM_
         if (debug) 
     &      print *,'isym ioff iioff i2off',isym,ioff,iioff,i2off
         NORBI=NORB_(ISYM)
         NBASI=NBAS_(ISYM)
         NISHI=NISH_(ISYM)
         NOCCI=NISH_(ISYM)+NASH_(ISYM)
         if (debug) print *,'norbi nbasi nocci ',norbi,nbasi,nocci
         CALL DCOPY(NOCCI,WORK(KF+I2OFF),NORBI+1,WORK(KEVAL+IOFF),1)
!
! inactive diagonals are a factor 2 two large
!
         CALL DSCAL(NISHI,DH,WORK(KEVAL+IOFF),1)
         IOFF  =  IOFF + NOCCI
         IIOFF = IIOFF + IROW(NORBI+1)
         I2OFF = I2OFF + NORBI*NORBI
      END DO
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(WORK(KEVAL),NOCCT_,1,
CVC     &       'Eigenvalues to LUSTEX from file '//FILE)
     &       'Eigenvalues to LUSTEX from file ...')
      END IF
!
! Append energies to LUSTEX
!
      CALL WRITT(LUSTEX,NOCCT_,WORK(KEVAL))
!
! Append orbital indices to LUSTEX
!
      IF (IPRSTX.GT.10) THEN
         WRITE(*,*) 'Orbital indices',NOCCT_,(I,I=1,NOCCT_)
      END IF
      WRITE(LUSTEX) NOCCT_, (I, I=1, NOCCT_) !?
      CALL QEXIT('STXAPPENE')
      END
C#######################################################################
      SUBROUTINE STXMOVLAP(LUSTEX,WORK,KFREE,LFREE)
#include "implicit.h"
      REAL*8 WORK(*)
      INTEGER LUSTEX,KFREE,LFREE
#include "inforb.h"
#include "stex.h"
      CALL QENTER('STXMOVLAP')
      IF (STEXMO) THEN
         NDIMT=NORBT
      ELSE
         NDIMT=NBAST
      END IF
C
C Allocate gound and final state orbitals
C
      CALL MEMGET2('REAL','GS',KGS,NDIMT*NOCCT,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','CH',KCH,NDIMT*NOCCT,WORK,KFREE,LFREE)
C
C Read ground state orbitals
C
      REWIND LUSTEX
      READ(LUSTEX)  !# basis functions, #occupied

      READ(LUSTEX)  !mc and nuclear energies
      CALL READT(LUSTEX,NDIMT*NOCCT,WORK(KGS))
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(WORK(KGS),NDIMT,NOCCT,'GS orbitals')
      END IF
      READ(LUSTEX)  ! orbital energies
      READ(LUSTEX)  ! orbital indices
C
C Read final state orbitals
C
      READ(LUSTEX)  !mc and nuclear energies
      CALL READT(LUSTEX,NDIMT*NOCCT,WORK(KCH))
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(WORK(KCH),NDIMT,NOCCT,'CH orbitals')
      END IF
      READ(LUSTEX)  ! orbital energies
      READ(LUSTEX)  ! orbital indices
      CALL MEMGET2('REAL','MOVLP',KMOVLP,N2OCCX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KMOVLP),N2OCCX)
      IF (STEXMO) THEN
C
C The core hole orbitals is already the MO overlap
C
         DO ISYM=1,NSYM
            CALL MCOPY(
     &         NOCC(ISYM),NOCC(ISYM),
     &         WORK(KCH),NORBT,
     &         WORK(KMOVLP+IOCC(ISYM)*(NOCCT+1)),NOCCT
     &         )
         END DO
      ELSE
C
C Read AO overlap matrix
C
         CALL MEMGET2('REAL','AOTR',KAOTR,NNBAST,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','AOSQ',KAOSQ,N2BAST,WORK,KFREE,LFREE)
         CALL RDONEL('OVERLAP ',.TRUE.,WORK(KAOTR),NNBAST)
         IF (IPRSTX.GT.0) THEN
            CALL HEADER('AO overlap from RDONEL',LUPRI)
            CALL OUTPKB(WORK(KAOTR),NBAS,NSYM,1,LUPRI)
         END IF
         ITROFF = 0
         ISQOFF = 0
         IMOOFF = 0
         IOVOFF = 0
         DO ISYM=1,NSYM
            IBASI = IBAS(ISYM)
            NNBASI = NNBAS(ISYM)
            N2BASI = N2BAS(ISYM)
            NBASI = NBAS(ISYM)
            NOCCI = NOCC(ISYM)
            CALL DSPTGE(NBASI,WORK(KAOTR+ITROFF),WORK(KAOSQ+ISQOFF))
   !
   ! Form Smo = c(gs)T Sao c(ch)
   !
            CALL STXOUTFUL(WORK(KAOSQ+ISQOFF),NBASI,NBASI,'AO overlap')
            CALL DGEMM('N','N',NBASI,NOCCI,NBASI,1.D0,
     &                 WORK(KAOSQ+ISQOFF),NBASI,
     &                 WORK(KCH+IMOOFF),NBAST,0.D0,
     &                 WORK(KFREE),NBASI)
            CALL DGEMM('T','N',NOCCI,NOCCI,NBASI,1.D0,
     &                 WORK(KGS+IMOOFF),NBAST,
     &                 WORK(KFREE),NBASI,0.D0,
     &                 WORK(KMOVLP+IOVOFF),NOCCT)
            ITROFF = ITROFF + NNBASI
            ISQOFF = ISQOFF + N2BASI
            IMOOFF = IMOOFF + NBAST*NOCCI + NBASI
            IOVOFF = IOVOFF + NOCCT*NOCCI + NOCCI
         END DO
      END IF
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(WORK(KMOVLP),NOCCT,NOCCT,'MO overlap')
      END IF
      CALL WRITT(LUSTEX,N2OCCX,WORK(KMOVLP))
      CALL MEMREL('STXMOVLAP',WORK,KGS,KGS,KFREE,LFREE)
      CALL QEXIT('STXMOVLAP')
      END
C#######################################################################
      SUBROUTINE STXDEN(CMO,DEN,WORK,KFREE,LFREE)
#include "implicit.h"
C
C Build modified densities for stex Hamiltonian
C
C Input:  MO coefficients, CMO
C Output: Coulomb and exchange density matrices, DC, DX
C
      REAL*8 CMO, DEN, WORK
      INTEGER KFREE, LFREE
C
C Global
C
#include "stex.h"
#include "inforb.h"
#include "infpri.h"
#include "priunit.h"
C
C Dimensions
C
      DIMENSION CMO(*), DEN(N2BASX,*), WORK(*)
C
C Local
C
      LOGICAL GETDC, GETDV
      PARAMETER (D2 = 2.0D0, DH = 0.5D0)
      LOGICAL DEBUG
C
      CALL QENTER('STXDEN')
      DEBUG=IPRSTX.GT.10
      IF (DEBUG) PRINT *,'---STXDEN called'
C
C Get normal density matrix
C
      GETDC = .TRUE.
      GETDV = .FALSE.
      KDC=1
      KDX=1+N2BASX
      IF (IPRSTX.GT.10) THEN
         CALL HEADER('STXDEN|Scaled MO coefficients',LUPRI)
         SQRT2=SQRT(D2)
         SQRTH=1/SQRT2
         CALL DSCAL(NCMOT,SQRT2,CMO,1)
         CALL PRORB(CMO,.FALSE.,LUPRI)
         CALL DSCAL(NCMOT,SQRTH,CMO,1)
      END IF
      CALL FCKDEN(GETDC,GETDV,DEN(1,1),DUM,CMO,DUM,WORK(KFREE),LFREE)
      CALL DCOPY(N2BASX,DEN(1,1),1,DEN(1,2),1)
      IF (IPRSTX.GT.1) THEN
         CALL HEADER('Density matrix from FCKDEN',LUPRI)
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C Construct densities corrections
C
C
C Loop over holes
C
      CALL MEMGET2('REAL','DCORR',KDCORR,N2BASX,WORK,KFREE,LFREE)
      ICOFF=0
      IDOFF=0
      CALL DZERO(WORK(KDCORR),N2BASX)
      DO ISYM=1,NSYM
         NBASI=NBAS(ISYM)
         NORBI=NORB(ISYM)
         DO I=1,NHLSTX(ISYM)
            IHL = IHLSTX(I,ISYM)
            ICST = ICOFF + NBASI*(IHL-1) + 1
            CALL DGEMM('N','N',NBASI,NBASI,1,1.D0,
     &                 CMO(ICST),NBASI,
     &                 CMO(ICST),1,0.D0,
     &                 WORK(KDCORR+IDOFF),NBAST)
            IF (IPRSTX.GT.10) THEN
               WRITE(LUPRI,'(/7X,A,I5,A,I5)')
     &            'Scaled Density correction, symmetry',ISYM,
     &            '  orbital',IHL
               KST = KDCORR+IDOFF
               KEN = KDCORR+IDOFF+NBASI*NBAST+NBASI
               CALL DSCAL(KEN-KST+1,D2,WORK(KST),1)
               CALL OUTPUT(WORK(KST),
     &                     1,NBASI,1,NBASI,NBAST,NBAST,1,LUPRI)
               CALL DSCAL(KEN-KST+1,DH,WORK(KST),1)
            END IF
C
C Form modified densities
C
            
            CJFAC=-CJSTX(I,ISYM)
            CKFAC= CKSTX(I,ISYM)*2
            CALL DAXPY(N2BASX,CJFAC,WORK(KDCORR),1,DEN(1,1),1)
            CALL DAXPY(N2BASX,CKFAC,WORK(KDCORR),1,DEN(1,2),1)
            CALL DZERO(WORK(KDCORR),N2BASX)
         END DO
         ICOFF=ICOFF+NBASI*NORBI
         IDOFF=IDOFF+NBASI*NBAST+NBASI
      END DO
C
C Final densities
C
      IF (IPRSTX.GT.1) THEN
         CALL HEADER('Final density matrices, coulomb ',LUPRI)
         CALL OUTPUT(DEN(1,1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Final density matrices, exchange ',LUPRI)
         CALL OUTPUT(DEN(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL MEMREL('STXDEN',WORK,KDCORR,KDCORR,KFREE,LFREE)
      CALL QEXIT('STXDEN')
      END
C#######################################################################
      SUBROUTINE STXFCK(DEN,FSTX,WORK,KFREE,LFREE)
#include "implicit.h"
C
C Purpose: given Coulomb and exchange density matrices calculate 
C STEX Fock matrix
C
C Input: density matrices, DC,DX
C Output: Stex Fock matrix, FSTX
C
      REAL*8 DEN(*), FSTX(*), WORK(*)
      INTEGER KFREE, LFREE
C
C
C
C Global
C
#include "stex.h"
#include "infpri.h"
#include "inforb.h"
#include "gnrinf.h"
C
C Local
C
      PARAMETER (D1=1.0D0, NDMAX=10)
      LOGICAL STXDIR
      INTEGER IFCTYP(NDMAX), ISYMDM(NDMAX)
      LOGICAL DEBUG
C
      DEBUG=IPRSTX.GT.10 .or. .true.
      IF (DEBUG) PRINT *,'---STXFCK called'
      CALL QENTER('STXFCK')
C
C Build two-electron STEX Fock matrix
C
c     print *,'stxfck: allocated memory ',kfree
      NDMAT = 2
      IFCTYP(1) = 11
      IFCTYP(2) = 12
      ISYMDM(1) = 1
      ISYMDM(2) = 1
      STXDIR = DIRCAL
      CALL DZERO(FSTX,NDMAT*N2BASX)
      CALL SIRFCK(FSTX,DEN,NDMAT,ISYMDM,IFCTYP,STXDIR,WORK(KFREE),LFREE) 
      IC = 1
      IX = 1 + N2BASX
      IF (IPRSTX.GT.10) THEN
       CALL STXOUTFUL(FSTX(IC),NBAST,NBAST,'Two-electron STEX coulomb')
       CALL STXOUTFUL(FSTX(IX),NBAST,NBAST,'Two-electron STEX exchange')
      END IF
      CALL DAXPY(N2BASX,D1,FSTX(IX),1,FSTX(IC),1)
      IF (DEBUG) THEN
       CALL STXOUTFUL(FSTX,NBAST,NBAST,'Two-electron STEX fock matrix')
      END IF
C
C ***** one-electron part in WRK(KH1AO)
C
      CALL MEMGET2('REAL','H1AO',KH1AO,N2BASX,WORK,KFREE,LFREE)
      CALL SIRH1(WORK(KH1AO),WORK(KFREE),LFREE)
      CALL MEMGET2('REAL','TMP',KTMP,NNBASX,WORK,KFREE,LFREE)
      CALL PKSYM1(WORK(KTMP),WORK(KH1AO),NBAS,NSYM,-1)
      CALL DSPTGE(NBAST,WORK(KTMP),WORK(KH1AO))
      IF (DEBUG) CALL STXOUTFUL(
     &   WORK(KH1AO),NBAST,NBAST,'One-electron fock matrix')
C
C  Add one+two 
C
      CALL DAXPY(N2BASX,D1,WORK(KH1AO),1,FSTX,1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(FSTX,NBAST,NBAST,'Total STEX fock matrix')
      END IF
C
C  Prepare interface file STEXFOCK
C
      LUSTXF = -1
      if (debug) write(2,*)'stxfck:open UNIT1 unit ',lustxf
      CALL GPOPEN(LUSTXF,'UNIT1','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &   0,.FALSE.)
      CALL STXAPPFCK(LUSTXF,FSTX,WORK,KFREE,LFREE)
      CALL STXAPPPRP(LUSTXF,WORK,KFREE,LFREE)
      CALL GPCLOSE(LUSTXF,'KEEP')
C
      CALL MEMREL('STXFCK',WORK,KH1AO,KH1AO,KFREE,LFREE)
      CALL QEXIT('STXFCK')
      END
C#######################################################################
      SUBROUTINE STXAPPFCK(LU,F,WORK,KFREE,LFREE)
#include "implicit.h"
      INTEGER LU, KFREE, LFREE
      REAL*8 F(*), WORK(*)
#include "stex.h"
#include "inforb.h"
#include "inftap.h"
#include "dummy.h"
#include "codata.h"
#include "priunit.h"
      CHARACTER*8 FILE
      INTEGER NSTEXF,ICALC,IHOLE,IHOL1,ISPIN,ICOUP,ICH
      LOGICAL FNDLAB,DEBUG
      EXTERNAL FNDLAB
      CALL QENTER('STXAPPFCK')
      DEBUG=IPRSTX.GT.10

      if (debug) print *,'ENTER STXAPPFCK---'
      if (debug) write(LUPRI,*)'appfck unit ',lu
      IF (STEXMO) THEN
         NDIMT=NORBT
         NNDIMX=NNORBX
      ELSE
         NDIMT=NBAST
         NNDIMX=NNBASX
      END IF
!
! Get energies, mo:s from SIRIFC1
!
      FILE='SIRIFC1 '
      LUSIFC=-1
      CALL GPOPEN(LUSIFC,FILE,'OLD','SEQUENTIAL','UNFORMATTED',
     &     0,.FALSE.)
      REWIND LUSIFC
      IF (FNDLAB(LBSIFC,LUSIFC) ) THEN
         READ (LUSIFC) POTNUC1, EMY1, EACTIV1, EMCSCF1
         READ (LUSIFC)
     &      NISHT1,NASHT1,NOCCT1,NORBT1,NBAST1
      ELSE
         CALL QUIT('ERROR (STXAPPFCK) MO''s not found on SIRIFC')
      END IF
      CALL GPCLOSE (LUSIFC,'KEEP')
!
! Get energies, mo:s from SIRIFC2
!
      FILE='SIRIFC2'
      LUSIFC=-1
      CALL GPOPEN(LUSIFC,FILE,'OLD','SEQUENTIAL','UNFORMATTED',
     &     0,.FALSE.)
      REWIND LUSIFC
      IF (FNDLAB(LBSIFC,LUSIFC) ) THEN
         READ (LUSIFC) POTNUC2, EMY2, EACTIV2, EMCSCF2
         READ (LUSIFC)
     &      NISHT2,NASHT2,NOCCT2,NORBT2,NBAST2
      ELSE
         CALL QUIT('ERROR (STXAPPFCK) MO''s not found on SIRIFC')
      END IF
      CALL GPCLOSE (LUSIFC,'KEEP')
      IF (IPRSTX.GT.10) THEN
         WRITE(*,*)' STXAPPFCK: Energies read from LUSIFC1-2'
         WRITE(*,*)' GS ', EMCSCF1
         WRITE(*,*)' CH ', EMCSCF2
      END IF
!
! Find gs orbital with maximum overlap with core hole orbital
! from mo overlap matrix
!
      LUSTEX = -1
      CALL GPOPEN (LUSTEX,'UNIT2','OLD','SEQUENTIAL','UNFORMATTED',
     &   0,.FALSE.)
      REWIND LUSTEX
      READ(LUSTEX) !nbast,nocct
      READ(LUSTEX) !gs energies
      READ(LUSTEX) !gs orbitals
      READ(LUSTEX) !gs orbital energies
      READ(LUSTEX) !gs orbital indices 
      READ(LUSTEX) !ch energies 
      READ(LUSTEX) !ch orbitals
      READ(LUSTEX) !ch orbital energies
      READ(LUSTEX) !ch orbital indices 
      CALL MEMGET2('REAL','MOV',KMOV,N2OCCX,WORK,KFREE,LFREE)
      CALL READT(LUSTEX,N2OCCX,WORK(KMOV))
      CALL GPCLOSE(LUSTEX, 'KEEP')
C
!
! Get open shell index of excited state
!
      ICH = 0
      DO ISYM=1,NSYM
         IF (NHLSTX(ISYM).EQ.1) THEN
            ICH = IOCC(ISYM) + IHLSTX(1,ISYM)
         END IF
      END DO
      if (debug) print *, 'open shell ich = ',ich
cvc
c  nasht2=2 for core-excited state orbitals
cvc
!D    write(LUPRI,*) ' nasht2 =',nasht2
c     IF (NASHT2.EQ.1.or.NASHT2.EQ.2) THEN
      IF (NASHT2.ne.0) THEN
!
! Get corresponding orbital index in ground state
!
         CALL GETCMX(ICH,WORK(KMOV),NOCCT,INGSCH,OVMAX)
         WRITE(LUPRI,'(7X,A/7X,A,I3,A,F6.4,A)')
     &      'Core hole orbital has maxixmum overlap with ',
     &      'ground state orbital ',INGSCH,
     &      ' (',OVMAX,')'
      ELSE IF (NASHT2.EQ.0) THEN
         WRITE(LUPRI,'(7X,A,I3,A,F6.4,A)')
     &      'Frozen orbital calculaton : orbital hole',
     &       ICH
      ELSE
         CALL QUIT('STXAPPFCK: NASHT2 =/= ZERO OR ONE')
      END IF
!
! Write to file
!
      CALL MEMGET2('REAL','TRI',KTRI,NNBASX,WORK,KFREE,LFREE)
      CALL DGETSP(NBAST,F,WORK(KTRI))
!
! Transform to MO
!
      IF (STEXMO) THEN
         IFSYM = 1
         CALL PRPTRA(
     &      'SIRIFC1 ',WORK(KTRI),IFSYM,WORK(KFREE),LFREE,
     &      .FALSE.
     &      )
         if (iprstx.gt.0) then
            call header('Static exchange Fock matrix in mo basis',lupri)
            call outpak(work(ktri),norbt,1,lupri)
         end if
      END IF
      IF (XAS) THEN
!
! 1. Number of Fock matrices
!
         NSTEXF = 1
         WRITE(LU) NSTEXF
!
! 2. IP, calculation type, spin, coupling, ihol2
!
         !
         !  Ionization potential: relaxed - scf energy difference
         !                        frozen  - open shell orbital energy
         IF (NASHT2.GT.0) THEN
            EIP = EMCSCF2-EMCSCF1
         ELSE
            LUMO=-1
            CALL GPOPEN(LUMO,'UNIT2','OLD','SEQUENTIAL','UNFORMATTED',
     &         0,.FALSE.)
            READ(LUMO)
            READ(LUMO)
            READ(LUMO)
            READ(LUMO) (EIP, I=1,ICH)
            CALL GPCLOSE(LUMO,'KEEP')
            if (debug) print *,' i ich eip ',i,ich,eip
            EIP=-EIP
         END IF
          
         ! ?
          
         ISPIN=0               
          
         ! Calculation type :  Closed shell ground state
          
         ICALC = 1        
          
         ! Singlet coupling
          
         ICOUP=0
         IF (IPRSTX.GT.10) THEN
            WRITE(*,*)'Ionization potential ',EIP
            WRITE(*,*)'Calculation type     ',ICALC
            WRITE(*,*)'Spin                 ',0
            WRITE(*,*)'Coupling             ',0
         END IF
         WRITE(LU) EIP, ICALC, ISPIN, ICOUP, ICH,
     &      NOCCT, DUMMY, DUMMY, DUMMY, DUMMY
!
! 3. ihole, ihol1 !???
!
         IHOLE = 0
         IHOL1 = 0
         IF (IPRSTX.GT.10) THEN
            WRITE(*,*)'IHOLE                ',IHOLE
            WRITE(*,*)'IHOL1                ',IHOL1
         END IF
         WRITE(LU) IHOLE,IHOL1 
!
! 4. fock matrix
!
         CALL WRITT(LU,NNDIMX,WORK(KTRI))
      END IF
      CALL QEXIT('STXAPPFCK')
      END
C#######################################################################
      SUBROUTINE STXEIG
      PRINT *,'---STXEIG called'
      END
C#######################################################################
      SUBROUTINE STXSTM
      PRINT *,'---STXSTM called'
      END
C#######################################################################
      BLOCK DATA STXBD
#include "implicit.h"
#include "stex.h"
      DATA DOSTEX/.FALSE./
      END

#ifdef INCLUDE_NOT_USED_CODE
C  /* Deck fckds3 */
C#######################################################################
      SUBROUTINE STXFCKDSX(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C     970303-tsaue : index permutation
C     941011-hjaaj: renamed from FCKDI1 to FCKDS3
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.50D00)
      INTEGER P, Q, R, S
#include "inforb.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST,*),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include "stex.h"
      LOGICAL DEBUG

C
      DEBUG=IPRSTX.GT.10
      IF (DEBUG) THEN
        call STXoutful(dmat(1,1,1),nbast,nbast,'in stxfckdsx density 1')
        call STXoutful(dmat(1,1,2),nbast,nbast,'in stxfckdsx density 2')
      END IF
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = DP5*GINT
         IF (R.EQ.S)              GINT = DP5*GINT
         IF (P.EQ.R .AND. S.EQ.Q) GINT = DP5*GINT
C           coulomb:
         FADD = GINT*(DMAT(R,S,1)+DMAT(S,R,1))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q,1)+DMAT(Q,P,1))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
C           exchange:
         GINT = DP5*GINT
         FMAT(P,R) = FMAT(P,R) - GINT*DMAT(Q,S,2)
         FMAT(R,P) = FMAT(R,P) - GINT*DMAT(S,Q,2)
         FMAT(P,S) = FMAT(P,S) - GINT*DMAT(Q,R,2)
         FMAT(S,P) = FMAT(S,P) - GINT*DMAT(R,Q,2)
         FMAT(Q,R) = FMAT(Q,R) - GINT*DMAT(P,S,2)
         FMAT(R,Q) = FMAT(R,Q) - GINT*DMAT(S,P,2)
         FMAT(Q,S) = FMAT(Q,S) - GINT*DMAT(P,R,2)
         FMAT(S,Q) = FMAT(S,Q) - GINT*DMAT(R,P,2)
         !print *,'fockmatrix after integral',p,q,r,s,gint
         !call STXoutful(fmat,nbast,nbast)
  100 CONTINUE
      RETURN
      END
C#######################################################################
C  /* Deck fckdsj */
      SUBROUTINE FCKDSJ(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen October 1994: general routine for
C symmetric singlet Fock matrix from symmetric singlet density matrix
C (based on FCKDS3)
C
C This subroutine adds derivative two-electron integrals to
C Fock matrices. The Fock matrices are assumed
C to be square matrices in full dimension without symmetry reduction
C in size. Remember to zero out the fock matrices before starting
C to accumulate.
C
C tsaue Feb 1997: indices permuted for exchange
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.50D00, DP25 = 0.25D0)
      INTEGER P, Q, R, S
#include "inforb.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)

C
      DO 100 INT = 1, LENGTH
         DINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DP5*DINT
         IF (R.EQ.S)              DINT = DP5*DINT
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DP5*DINT
         EINT = DP25*DINT
C        ... because DMAT is multiplied by extra 2 in SIRFCK
C            in order to get full Coulomb contribution:
C            DINT*(DMAT(R,S) + DMAT(S,R)) = 2*DINT*DMAT(R,S)
C            Compare general code in FCKDS3.
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(S,R)
C        FMAT(P,R) = FMAT(P,R) - EINT*DMAT(Q,S)
C        FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(P,S)
C        FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(P,R)
C        FMAT(P,S) = FMAT(P,S) - EINT*DMAT(Q,R)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(Q,P)
C        ... DMAT(i,j) = DMAT(j,i) has been used above
  100 CONTINUE
      RETURN
      END
#endif   /* INCLUDE_NOT_USED_CODE */
C#######################################################################
      SUBROUTINE STXAPPPRP(LUSTXF,WORK,KFREE,LFREE)
#include "implicit.h"
#include "priunit.h"
      INTEGER LUSTXF, KFREE, LFREE
      REAL*8 WORK(*)
#include "inftap.h"
#include "inforb.h"
#include "dummy.h"
#include "stex.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"
#include "dipole.h"
      LOGICAL FNDLB2,DEBUG
      EXTERNAL FNDLB2
      CHARACTER*8 RTNLBL(2)
      PARAMETER (NPROP = 9)
      CHARACTER*8 PROP(NPROP)
      INTEGER SYMPRP(NPROP)
      LOGICAL ANTISY(NPROP)
      DATA PROP/'XDIPLEN ', 'YDIPLEN ', 'ZDIPLEN',
     &          'XDIPVEL ', 'YDIPVEL ', 'ZDIPVEL',
     &          'XANGMOM ', 'YANGMOM ', 'ZANGMOM'/
      DATA ANTISY/
     &   .FALSE., .FALSE., .FALSE.,
     &   .TRUE., .TRUE., .TRUE.,
     &   .TRUE., .TRUE., .TRUE./

      CALL QENTER('STXAPPPRP')
      CALL FLSHFO(LUPRI)
      DEBUG=IPRSTX.GT.10
      DO I=1,3
         SYMPRP(I) = ISYMAX(I,1)+1
         SYMPRP(I+3) = SYMPRP(I)
         SYMPRP(I+6) = ISYMAX(I,2)+1
      END DO
      if (debug) write(LUPRI,*)'appprp:lustxf',lustxf
      IF (STEXMO) THEN
         NDIMT=NORBT
         NNDIMT=NNORBT
         NNDIMX=NNORBX
      ELSE
         NDIMT=NBAST
         NNDIMT=NNBAST
         NNDIMX=NNBASX
      END IF
!
! Calculate integrals if not found
!
      LUPROP=-1
      CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &   0,.FALSE.)
      REWIND LUPROP
      IF (.NOT.FNDLB2('XDIPLEN ',RTNLBL,LUPROP)) THEN
         if (debug) write(LUPRI,*)'calling pr1int DIPLEN'
         CALL PR1INT('DIPLEN ',WORK(KFREE),LFREE,IDUMMY,
     &               IDUMMY,.TRUE.,DEBUG,0)
      END IF
      REWIND LUPROP
      IF (.NOT.FNDLB2('XDIPVEL ',RTNLBL,LUPROP)) THEN
         if (debug) write(LUPRI,*)'calling pr1int DIPVEL'
         CALL PR1INT('DIPVEL ',WORK(KFREE),LFREE,IDUMMY,
     &               IDUMMY,.TRUE.,DEBUG,0)
      END IF
      IF (.NOT.FNDLB2('XANGMOM ',RTNLBL,LUPROP)) THEN
         if (debug) write(LUPRI,*)'calling pr1int DIPVEL'
         CALL PR1INT('ANGMOM ',WORK(KFREE),LFREE,IDUMMY,
     &               IDUMMY,.TRUE.,DEBUG,0)
      END IF
!
! Copy overlap and h1 integrals to LUSTXF
!
      CALL MEMGET2('REAL','AONNT',KAONNT,NNBAST,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','AONNX',KAONNX,NNBASX,WORK,KFREE,LFREE)

      CALL RDONEL('OVERLAP ',.TRUE.,WORK(KAONNT),NNBAST)
      CALL DZERO(WORK(KAONNX),NNBASX)
      if (debug) call header('STXAPPPRP:nt2nx for overlap',LUPRI)
      CALL NT2NX(WORK(KAONNT),WORK(KAONNX))
      IF (STEXMO) THEN
         CALL PRPTRA(
     &      'SIRIFC1 ',WORK(KAONNX),1,WORK(KFREE),LFREE,.FALSE.
     &      )
         CALL WRITT(LUSTXF,NNORBX,WORK(KAONNX))
      ELSE
         CALL WRITT(LUSTXF,NNBASX,WORK(KAONNX))
      END IF
!
! Copy property integrals to LUSTXF
!
      REWIND LUPROP
      CALL MEMGET2('REAL','PRP',KPRP,NNBASX,WORK,KFREE,LFREE)
      if (debug) write(LUPRI,*) 'kprp,kfree,lfree',kprp,kfree,lfree
      DO IPROP=1,NPROP
         IF (FNDLB2(PROP(IPROP),RTNLBL,LUPROP)) THEN
            CALL READT(LUPROP,NNBASX,WORK(KPRP))
            IF (STEXMO) THEN
               if (debug) 
     &            call header('Transforming '//prop(iprop),LUPRI)
               CALL PRPTRA('SIRIFC1 ',WORK(KPRP),SYMPRP(IPROP),
     &                     WORK(KFREE),LFREE,ANTISY(IPROP))
               CALL WRITT(LUSTXF,NNORBX,WORK(KPRP))
            ELSE
               CALL WRITT(LUSTXF,NNBASX,WORK(KPRP))
            END IF
         ELSE
            CALL QUIT('STXAPPPRP: ERROR READING INTEGRALS FROM LUPROP:'
     &              //PROP(IPROP))
         END IF
      END DO
      CALL GPCLOSE(LUPROP,'KEEP')
!
! Skip force
!
      CALL DZERO(WORK(KPRP),NNDIMX)
      CALL WRITT(LUSTXF,NNDIMX,WORK(KPRP))
      CALL WRITT(LUSTXF,NNDIMX,WORK(KPRP))
      CALL WRITT(LUSTXF,NNDIMX,WORK(KPRP))
!
! Skip torque
!
      CALL WRITT(LUSTXF,NNDIMX,WORK(KPRP))
      CALL WRITT(LUSTXF,NNDIMX,WORK(KPRP))
      CALL WRITT(LUSTXF,NNDIMX,WORK(KPRP))
!
! Append nuclear dipole moments  (need to initialize common dipole.h)
!
      CALL MEMGET2('REAL','CSTRA',KCSTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','SCTRA',KSCTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL DIPNUC(WORK(KCSTRA),WORK(KSCTRA),0,.FALSE.)
      CALL WRITT(LUSTXF,3,DIPMN)
      IF (IPRSTX.GT.0) THEN
         CALL STXOUTFUL(DIPMN,1,3,'Nuclear dipole moment')
      END IF
!
! Append nuclear quadrupole moments
!
      WRITE(LUSTXF) DUMMY,DUMMY,DUMMY
!
      CALL MEMREL('STXAPPPRP',WORK,KPRP,KPRP,KFREE,LFREE)
      CALL QEXIT('STXAPPPRP')
      END
C#######################################################################
      SUBROUTINE STXOUTFUL(A,M,N,TEXT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(M,N)
      CHARACTER TEXT*(*)
      CALL HEADER(TEXT,LUPRI)
      CALL OUTPUT(A,1,M,1,N,M,N,1,LUPRI)
      CALL FLSHFO(LUPRI)
      END
C#######################################################################
      SUBROUTINE APPORB(LABEL3,CMO,REWIT1)
C
C  19-Oct-1989; 3-Aug-1990 (REWIT1; hjaaj)
C
C  If (REWIT1) then
C     Write molecular orbitals with label 'LABEL3  ' at
C     beginning of LUIT1
C  else
C     Overwrite 'LABEL3  ' on LUIT1
C  end if
C
#include "implicit.h"
      CHARACTER*8 LABEL3
      DIMENSION CMO(NCMOT)
      LOGICAL   REWIT1
C
#include "dummy.h"
C
C Used from common blocks:
C  INFORB : NCMOT
C  INFTAP : LUIT1
C
#include "inforb.h"
#include "inftap.h"
C
      LOGICAL FNDLAB
      CHARACTER*8 MOLBL(3),LABTBL(2)
      DATA MOLBL /'********','        ','        '/
      DATA LABTBL/'APPORB  ','EODATA  '/
C
      CALL QENTER('APPORB')
C
C     Write orbitals to LUIT1
C
      CALL GETDAT(MOLBL(2),MOLBL(3))
      NCMOT4 = MAX(4,NCMOT)
      IF (REWIT1) THEN
         CALL NEWIT1
      ELSE
         REWIND LUIT1
         IF (FNDLAB(LABTBL(1),LUIT1)) THEN
            BACKSPACE LUIT1
         ELSE
            REWIND LUIT1
            IF (FNDLAB(LABTBL(2),LUIT1)) BACKSPACE LUIT1
         END IF
      END IF
      WRITE (LUIT1) MOLBL(1),MOLBL(2),LABTBL(1),LABEL3
C     ...label 'APPORB  '
      CALL WRITT(LUIT1,NCMOT4,CMO)
      WRITE (LUIT1) MOLBL,LABTBL(2)
C     ...label 'EODATA  '
      REWIND LUIT1
C
      CALL QEXIT('APPORB')
      RETURN
      END
C#######################################################################
      SUBROUTINE LTRA(LABEL,CMOCC,WORK,LWORK)
#include "implicit.h"
      REAL*8  CMOCC(*), WORK(*)
      INTEGER  LWORK
      CHARACTER*8 LABEL
#include "inftap.h"
#include "inforb.h"
#include "stex.h"
      LOGICAL FNDLAB, DEBUG
      EXTERNAL FNDLAB
      CALL QENTER('LTRA')
      DEBUG=IPRSTX.GT.10
      IF (DEBUG) THEN
         print *,'ENTER LTRA---'
         call STXoutful(cmocc,nbast,nocct,'LTRA input CMOCC')
      END IF
C
C Transform input mo:s by mo:s defined by LABEL
C
      KFREE=1
      LFREE=LWORK
C     OPEN(LUIT1,FORM='UNFORMATTED',STATUS='OLD')
      REWIND LUIT1
      IF (FNDLAB(LABEL,LUIT1)) THEN
         CALL MEMGET2('REAL','CMO',KCMO,NCMOT,WORK,KFREE,LFREE)
         CALL READT(LUIT1,NCMOT,WORK(KCMO))
         IF (IPRSTX.GT.10) THEN
            call prorb(work(kcmo),.false.,LUPRI)
         END IF 
      ELSE
         CALL QUIT('LTRA: LABEL "'//LABEL//'" NOT FOUND')
      END IF
C     CLOSE(LUIT1)
C
C
C Get AO overlap
C
      CALL MEMGET2('REAL','AOTR',KAOTR,NNBAST,WORK,KFREE,LFREE)
      CALL RDONEL('OVERLAP ',.TRUE.,WORK(KAOTR),NNBAST)
      IF (DEBUG) THEN
         call header('LTRA:overlap',LUPRI)
         call outpkb(WORK(KAOTR),NBAS,NSYM,1,LUPRI)
      END IF
C
C Form SAO CM
C
      CALL MEMGET2('REAL','AOSQ',KAOSQ,N2BAST,WORK,KFREE,LFREE)
      DO ISYM=1,NSYM
         NBASI=NBAS(ISYM)
         NORBI=NORB(ISYM)
         ICMOI=ICMO(ISYM)
         ITROFF=IIBAS(ISYM)
         ISQOFF=I2BAS(ISYM)
         CALL DSPTGE(NBASI,WORK(KAOTR+ITROFF),WORK(KAOSQ+ISQOFF))
         CALL DGEMM('N','N',NBASI,NORBI,NBASI,1.D0,
     &              WORK(KAOSQ+ISQOFF),NBASI,
     &              WORK(KCMO+ICMOI),NBASI,0.D0,
     &              WORK(KFREE),NBASI)
         CALL DCOPY(NBASI*NORBI,WORK(KFREE),1,WORK(KCMO+ICMOI),1)
      END DO
      IF (DEBUG) THEN
         call header('LTRA:overlap square packed',LUPRI)
         call outptb(WORK(KAOSQ),NBAS,NSYM,1,LUPRI)
         call header('after   S*cmo',LUPRI)
         call prorb(work(kcmo),.false.,LUPRI)
      END IF
      
C
C Left transform input matrix
C
      CALL MEMGET2('REAL','MO',KMO,NORBT*NOCCT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KMO),NORBT*NOCCT)
      if (debug) write(LUPRI,*)'kmo kfree lfree',kmo,kfree,lfree
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(CMOCC,NBAST,NOCCT,'LTRA input AO BASIS')
      END IF
      DO ISYM=1,NSYM
         if (debug) write(LUPRI,*) 'isym',isym
         IA = IOCC(ISYM)*NBAST+IBAS(ISYM)
         IM = IOCC(ISYM)*NORBT+IORB(ISYM)
         if (debug) write(LUPRI,*) ' ia im ',ia,im
         CALL DGEMM('T','N',NORB(ISYM),NOCC(ISYM),NBAS(ISYM),1.D0,
     &              WORK(KCMO+ICMO(ISYM)),NBAS(ISYM),
     &              CMOCC(IA+1),NBAST,0.D0,
     &              WORK(IM+KMO),NORBT)
      END DO
      IF (IPRSTX.GT.10) THEN
         CALL STXOUTFUL(WORK(KMO),NORBT,NOCCT,'LTRA output MO BASIS')
      END IF
      CALL DCOPY(NORBT*NOCCT,WORK(KMO),1,CMOCC,1)

      CALL MEMREL('LTRA',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('LTRA')
      END
C#######################################################################
      SUBROUTINE PRPTRA(LABEL,PRPTR,PRPSYM,WORK,LWORK,ANTI)
#include "implicit.h"
      REAL*8  PRPTR(*), WORK(*)
      INTEGER  LWORK, PRPSYM
      CHARACTER*8 LABEL
#include "inftap.h"
#include "inforb.h"
#include "stex.h"
      LOGICAL FNDLAB
      EXTERNAL FNDLAB
      LOGICAL ANTI, DEBUG
      CALL QENTER('PRPTRA')
C
C Transform input mo:s by mo:s defined by LABEL
C
      DEBUG=IPRSTX.GT.10
      KFREE=1
      LFREE=LWORK
C     OPEN(LUIT1,FORM='UNFORMATTED',STATUS='OLD')
      REWIND LUIT1
      if (debug) write(LUPRI,*)'fndlab label luit1 ',label,luit1
      IF (FNDLAB(LABEL,LUIT1)) THEN
         CALL MEMGET2('REAL','CMO',KCMO,NCMOT,WORK,KFREE,LFREE)
         CALL READT(LUIT1,NCMOT,WORK(KCMO))
         if (debug) call prorb(work(kcmo),.false.,LUPRI)
      ELSE
         CALL QUIT('PRPTRA: LABEL "'//LABEL//'" NOT FOUND')
      END IF
C     CLOSE(LUIT1)
C
C Form square and transform
C
      CALL MEMGET2('REAL','PRPSQ',KPRPSQ,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','PRPMO',KPRPMO,N2ORBX,WORK,KFREE,LFREE)
      IF (ANTI) THEN
         CALL DAPTGE(NBAST,PRPTR,WORK(KPRPSQ))
      ELSE
         CALL DSPTGE(NBAST,PRPTR,WORK(KPRPSQ))
      END IF
      CALL DZERO(WORK(KPRPMO),N2ORBX)
      IF (DEBUG) THEN
         print *,'operator symmetry ',prpsym
         call STXoutful(work(kprpsq),nbast,nbast,'prptra before')
      END IF
      DO ISYM=1,NSYM
         JSYM=MULD2H(PRPSYM,ISYM)
         NBASI=NBAS(ISYM)
         NORBI=NORB(ISYM)
         ICMOI=ICMO(ISYM)
         NBASJ=NBAS(JSYM)
         NORBJ=NORB(JSYM)
         ICMOJ=ICMO(JSYM)
         IF (DEBUG) THEN
            write(LUPRI,*)
     &         'isym jsym nbasi nbasj norbi norbj icmoi icmoj ',
     &          isym,jsym,nbasi,nbasj,norbi,norbj,icmoi,icmoj
         END IF
         CALL UTHV(
     &      WORK(KCMO+ICMOI), WORK(KPRPSQ), WORK(KCMO+ICMOJ),
     &      ISYM, JSYM, NBASI, NBASJ, WORK(KPRPMO), WORK(KFREE)
     &      )
      END DO
      IF (DEBUG) THEN
         call STXoutful(work(kprpmo),norbt,norbt,'prptra after ')
      END IF
      IF (ANTI) THEN
         CALL DGETAP(NORBT,WORK(KPRPMO),PRPTR)
      ELSE
         CALL DGETSP(NORBT,WORK(KPRPMO),PRPTR)
      END IF
      CALL MEMREL('PRPTRA',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('PRPTRA')
      END
C#######################################################################
      SUBROUTINE NT2NX(ANT,ANX)
      REAL*8 ANT(*), ANX(*)
#include "inforb.h"
#include "maxorb.h"
#include "maxash.h"
#include "infind.h"
      CALL QENTER('NT2NX')
      CALL DCOPY(NNBAS(1),ANT,1,ANX,1)
      DO ISYM=2,NSYM
         DO IB=1,NBAS(ISYM)
            CALL DCOPY(
     &         IB,
     &         ANT(IIBAS(ISYM)+IROW(IB)+1),1,
     &         ANX(IROW(IBAS(ISYM)+IB)+IBAS(ISYM)+1),1
     &         )
         END DO
      END DO
      CALL QEXIT('NT2NX')
      END
C#######################################################################
      SUBROUTINE GETCMX(NCOL,A,NDIM,NROW,VAL)
#include "implicit.h"
      INTEGER NCOL,NROW,NDIM
      REAL*8 A(NDIM,*), VAL
      NROW=IDAMAX(NDIM,A(1,NCOL),1)
      VAL=A(NROW,NCOL)
      END
C#######################################################################
      SUBROUTINE UAUG(CMO1,NAO1,NMO1,CMO,NAO,NMO,S,WORK,LWORK)
#include "implicit.h"
      DIMENSION CMO1(NAO1,NMO1), CMO(NAO,NMO), WORK(LWORK), S(NAO,NAO)
      LOGICAL DEBUG
      data debug/.false./
      IF (NAO1.EQ.NAO .OR. NAO1.EQ.0) RETURN
C
C CMO1 is original MO
C
C CMO  is trial MO is augmented basis (A B) where A is the dimension of CMO1
C
C   1. Project CMO1 out of CMO
C   2. Select linearly independent components to span the virtual space
C
C
C 1. Form MO overlap <1|A>, <1|B>
C
C
      DEBUG=.TRUE.
      CALL QENTER('UAUG')
      !IF (DEBUG) THEN
         CALL STXOUTFUL(CMO1,NAO1,NMO1,'UAUG:Initial MO')
         CALL STXOUTFUL(CMO,NAO,NMO,'UAUG:Trial MO (augmented)')
      !ENDIF 
      NAO2=NAO-NAO1
      NMO2=NMO-NMO1
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET2('REAL','1A',K1A,NMO1*NMO1,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','1B',K1B,NMO1*NMO2,WORK,KFREE,LFREE)
      LTMP=NAO*MAX(NMO1,NMO2)
      CALL MEMGET2('REAL','TMP',KTMP,LTMP,WORK,KFREE,LFREE)
C
C  <1|A>
C
C  SAO  => S(1:NAO1,1:NAO)
C  CMOA => CMO(1:NAO,1:NMO1)
C  S1A  =  MATMUL( TRANSPOSE(CMO1) , SAO, CMOA )
C
      CALL DGEMM('N','N',NAO1,NMO1,NAO,1.D0,
     &           S,NAO,
     &           CMO,NAO,0.D0,
     &           WORK(KTMP),NAO1)
      CALL DGEMM('T','N',NMO1,NMO1,NAO1,1.D0,
     &           CMO1,NAO1,
     &           WORK(KTMP),NAO1,0.D0,
     &           WORK(K1A),NMO1)
      IF (DEBUG) THEN
         CALL HEADER('SX',2)
         CALL OUTPUT(S,1,NAO1,1,NAO,NAO,NAO,1,2)
         CALL STXOUTFUL(WORK(KTMP),NAO1,NMO1,'SCA')
         CALL STXOUTFUL(WORK(K1A),NMO1,NMO1,'UAUG: <1|A> ')
         CALL MEMCHK('UAUG 1A',WORK,K1A)
      END IF
C
C  <1|B>
C
C  SAO => S(1:NAO1,1:NAO)
C  CMOB => CMO(1:NAO,NMO1+1:NMO)
C  S1B = MATMUL( TRANSPOSE(CMO1) , SAO, CMOB )
C
      CALL DGEMM('N','N',NAO1,NMO2,NAO,1.D0,
     &           S,NAO,
     &           CMO(1,NMO1+1),NAO,0.D0,
     &           WORK(KTMP),NAO1)
      CALL DGEMM('T','N',NMO1,NMO2,NAO1,1.D0,
     &           CMO1,NAO1,
     &           WORK(KTMP),NAO1,0.D0,
     &           WORK(K1B),NMO1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(K1B),NMO1,NMO2,'UAUG: <1|B> ')
         CALL MEMCHK('UAUG 1B',WORK,K1A)
      END IF
C
C Form MO overlap ov projected basis (1-P1) ( |A> |B> )
C
C           S = ( 1-<A|1><1|A>      - <A|1><1|B> )
C               (  -<B|1><1|A>    1 - <B|1><1|B> )
C
      CALL MEMGET2('REAL','S',KS,NMO*NMO,WORK,KFREE,LFREE)
      CALL DUNIT(WORK(KS),NMO)
      ISAA = KS
      ISBA = KS + NMO1
      ISAB = KS + NMO1*NMO
      ISBB = KS + NMO1*NMO + NMO1
      CALL DGEMM('T','N',NMO1,NMO1,NMO1,-1.D0,
     &           WORK(K1A),NMO1,
     &           WORK(K1A),NMO1,1.D0,
     &           WORK(ISAA),NMO)
      CALL DGEMM('T','N',NMO1,NMO2,NMO1,-1.D0,
     &           WORK(K1A),NMO1,
     &           WORK(K1B),NMO1,1.D0,
     &           WORK(ISAB),NMO)
      CALL DGEMM('T','N',NMO2,NMO1,NMO1,-1.D0,
     &           WORK(K1B),NMO1,
     &           WORK(K1A),NMO1,1.D0,
     &           WORK(ISBA),NMO)
      CALL DGEMM('T','N',NMO2,NMO2,NMO1,-1.D0,
     &           WORK(K1B),NMO1,
     &           WORK(K1B),NMO1,1.D0,
     &           WORK(ISBB),NMO)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KS),NMO,NMO,'projected overlap')
         CALL MEMCHK('UAUG KS',WORK,K1A)
      END IF
C
C Diagonalize MO overlap
C
      CALL MEMGET2('REAL','EIG',KEIG,NMO,WORK,KFREE,LFREE)
      CALL DSYEV(
     &   'V','U',NMO,WORK(KS),NMO,WORK(KEIG),WORK(KFREE),LFREE,INFO
     &   ) 
      IF (INFO.NE.0) THEN
         CALL QUIT('UAUG:DSYEV diagonalization failed')
      END IF
C
C WORK(KS) now contains the eigenvectors
C
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KEIG),1,NMO1,
     &      'Eigenvalues of rejected vectors')
         CALL STXOUTFUL(WORK(KEIG+NMO1),1,NMO2,
     &      'Eigenvalues of selected vectors')
         CALL STXOUTFUL(WORK(KS),NMO,NMO,'Eigenvectors')
      END IF
C
C  Projected trial vectors
C
C  CMOA => CMO(:,1:NMO1)
C  CMOB => CMO(:,NMO+1:NMO2)
C
C  CMOA = CMOA - MATMUL(TRANSPOSE(CMO1), S1A)
C  CMOB = CMOB - MATMUL(TRANSPOSE(CMO1), S1B)
C
      CALL DGEMM('T','N',NAO1,NMO1,NMO1,-1.D0,
     &           CMO1,NMO1,
     &           WORK(K1A),NMO1,1.D0,
     &           CMO(1,1),NAO)
      CALL DGEMM('T','N',NAO1,NMO2,NMO1,-1.D0,
     &           CMO1,NMO1,
     &           WORK(K1B),NMO1,1.D0,
     &           CMO(1,NMO1+1),NAO)
      IF (DEBUG) THEN
         CALL STXOUTFUL(CMO,NAO,NMO,'Projected mo')
      END IF
C
C   Transform the projected trial vectors with the selected
C   overlap eigenvectors
C
C  CMOB = MATMUL(CMO,S(1,NMO1+1:NMO))
C
      CALL DGEMM('N','N',NAO,NMO2,NMO,1.D0,
     &           CMO,NAO,
     &           WORK(ISAB),NMO,0.D0,
     &           WORK(KTMP),NAO)
      CALL DCOPY(NAO*NMO2,WORK(KTMP),1,CMO(1,NMO1+1),1)
      CALL DZERO(CMO,NAO*NMO1)
      CALL MCOPY(NAO1,NMO1,CMO1,NAO1,CMO,NAO)
      IF (DEBUG) THEN
         CALL STXOUTFUL(CMO,NAO,NMO,'UAUG:Final MO')
      END IF
      CALL MEMREL('UAUG',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('UAUG')
      END
C#######################################################################
      SUBROUTINE UAUGOLD(CMO1,NAO1,NMO1,CMO,NAO,NMO,S,WORK,LWORK)
#include "implicit.h"
      DIMENSION CMO1(NAO1,NMO1), CMO(NAO,NMO), WORK(LWORK), S(NAO,NAO)
      LOGICAL DEBUG
      data debug/.true./
      IF (NAO1.EQ.NAO .OR. NAO1.EQ.0) RETURN
C
C CMO1 is original MO
C
C CMO  is trial MO is augmented basis (A B) where A is the dimension of CMO1
C
C Set up the unitary transformation defined by
C
C        U CMO(A) = CMO1
C
C conditions on U
C
C   U(A,A) = <A|1>
C   U(B,A) = <B|1>
C   U(A,B)+ = - U(B,B)**(-1) <B|1><1|A>
C   U(B,B) = SQRT( 1 - <B|1><1|B>)
C
C      
C Purpose: define virtual space in the augmented basis
C          so that linear dependencies is avoided
C
      CALL QENTER('UAUG')
      KFREE = 1
      LFREE = LWORK
      NAO2 = NAO - NAO1
      NMO2 = NMO - NMO1
      IF (DEBUG) THEN
         CALL STXOUTFUL(CMO1,NAO1,NMO1,'UAUG input: CMO1')
         CALL STXOUTFUL(CMO,NAO,NMO,'UAUG input: CMO')
         CALL STXOUTFUL(S,NAO,NAO,'UAUG input: S')
      END IF
C
C 1. Form MO overlap <A|1>, <1|B>
C
C
      CALL MEMGET2('REAL','A1',KA1,NMO1*NMO1,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','1B',K1B,NMO1*NMO2,WORK,KFREE,LFREE)
      NTMP=NAO*MAX(NMO1,NMO2)
      CALL MEMGET2('REAL','TMP',KTMP,NTMP,WORK,KFREE,LFREE)
C
C  <A|1>
C
C  CMOA => CMO(1:NAO,1:NMO1)
C  SAO  => S(1:NAO,1:NAO1)
C  SA1  =  MATMUL( TRANSPOSE(CMOA) , SAO, CMO1 )
C
      CALL DGEMM('N','N',NAO,NMO1,NAO1,1.D0,
     &           S,NAO,
     &           CMO1,NAO1,0.D0,
     &           WORK(KTMP),NAO)
      CALL DGEMM('T','N',NMO1,NMO1,NAO,1.D0,
     &           CMO,NAO,
     &           WORK(KTMP),NAO,0.D0,
     &           WORK(KA1),NMO1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KA1),NMO1,NMO1,'UAUG: <A|1> ')
         CALL MEMCHK('UAUG A1',WORK,KA1)
      END IF
C
C  <1|B>
C
C  CMOB => CMO(1:NAO,NMO1+1:NMO)
C  SAO => S(1:NAO1,1:NAO)
C  S1B = MATMUL( TRANSPOSE(CMO1) , SAO, CMOB )
C
      CALL DGEMM('N','N',NAO1,NMO2,NAO,1.D0,
     &           S,NAO,
     &           CMO(1,NMO1+1),NAO,0.D0,
     &           WORK(KTMP),NAO1)
      CALL DGEMM('T','N',NMO1,NMO2,NAO1,1.D0,
     &           CMO1,NAO1,
     &           WORK(KTMP),NAO1,0.D0,
     &           WORK(K1B),NMO1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(K1B),NMO1,NMO2,'UAUG: <1|B> ')
         CALL MEMCHK('UAUG 1B',WORK,KA1)
      END IF
C
C 2.  Form U(B,B)
C 
C
C    SBB = MATMUL (TRANPOSE(S1B), S1B))
C
      CALL MEMGET2('REAL','SBB',KSBB,NMO2*NMO2,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','UBB',KUBB,NMO2*NMO2,WORK,KFREE,LFREE)
      CALL DGEMM('T','N',NMO2,NMO2,NMO1,1.D0,
     &           WORK(K1B),NMO1,
     &           WORK(K1B),NMO1,0.D0,
     &           WORK(KSBB),NMO2)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KSBB),NMO2,NMO2,'UAUG: <B|1><1|B>')
         CALL MEMCHK('UAUG SBB',WORK,KA1)
      END IF
C
C     UBB = SQRT(1 - SBB) diagnoalize
C
      D1=1.0D0
      CALL DSCAL(NMO2*NMO2,-D1,WORK(KSBB),1)
      DO I=0,NMO2*NMO2,NMO2+1
         WORK(KSBB+I) = WORK(KSBB+I) + D1
      END DO
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KSBB),NMO2,NMO2,'UAUG: (1-SBB)')
      END IF
      IF (DEBUG) CALL DCOPY(NMO2*NMO2,WORK(KSBB),1,WORK(KTMP),1)
      CALL MATFUN('SQRT',NMO2,WORK(KSBB),WORK(KUBB),WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL DGEMM('N','N',NMO2,NMO2,NMO2,1.D0,
     &              WORK(KUBB),NMO2,
     &              WORK(KUBB),NMO2,0.D0,
     &              WORK(KSBB),NMO2)
         CALL STXOUTFUL(WORK(KSBB),NMO2,NMO2,' UBB*UBB ')
      END IF
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KUBB),NMO2,NMO2,'UAUG: UBB')
         CALL MEMCHK('UAUG UBB',WORK,KA1)
      END IF
C
C  Transform trial MO:s with U(B,B)
C
C  CMO(:,NMO1+1:NMO) = MATMUL( CMO(:,NMO+1:NMO), UBB )
C
      CALL DGEMM('N','N',NAO,NMO2,NMO2,1.D0,
     &           CMO(1,NMO1+1),NAO,
     &           WORK(KUBB),NMO2,0.D0,
     &           WORK(KTMP),NAO)
      CALL DCOPY(NAO*NMO2, WORK(KTMP),1,CMO(1,NMO1+1),1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(
     &      CMO(1,NMO1+1),NAO,NMO2,
     &      ' First transformation of trial virtuals: UBB'
     &      )
         CALL MEMCHK('UAUG: CMO*UBB',WORK,KA1)
      END IF
C
C 3.  Form U(A,B)
C
C   The transpose
C   <A|1><1|B> forms right-hand side vectors for the solution of
C   U(A,B)
C
      CALL MEMGET2('REAL','A11B',KA11B,NMO1*NMO2,WORK,KFREE,LFREE)
      CALL DGEMM('N','N',NMO1,NMO2,NMO1,1.D0,
     &           WORK(KA1),NMO1,
     &           WORK(K1B),NMO1,0.D0,
     &           WORK(KA11B),NMO1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KA11B),NMO1,NMO2,'UAUG: <A|1><1|B>')
         CALL MEMCHK('UAUG UAB RHS',WORK,KA1)
      END IF
      CALL STX_TRANSP(NMO1,NMO2,WORK(KA11B),WORK(KTMP))
      KB11A = KA11B
      CALL DCOPY(NMO1*NMO2,WORK(KTMP),1,WORK(KB11A),1)
      CALL DSCAL(NMO1*NMO2,-D1,WORK(KB11A),1)
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KB11A),NMO2,NMO1,'UAUG: -<B|1><1|A>')
         CALL MEMCHK('UAUG TRSP RHS',WORK,KA1)
      END IF
C              UBB UAB+ = - <B|1><1|A>
C     
      CALL MEMGET2('INTE','IPIV',KIPIV,NMO2,WORK,KFREE,LFREE)
      LUBB = KTMP
      CALL DCOPY(NMO2*NMO2,WORK(KUBB),1,WORK(LUBB),1)
      CALL DGETRF(!NOTE LU-DECOMPOSES UBB
     &   NMO2, NMO2, WORK(LUBB), NMO2, 
     &   WORK(KIPIV), INFO
     &   )
      IF (INFO.NE.0) THEN
         CALL QUIT('UAUG: SGETRF FACTORIZATION FAILED')
      END IF
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(LUBB),NMO2,NMO2,'UAUG: LU(UBB)')
         CALL MEMCHK('UAG LU(UBB)',WORK,KA1)
      END IF
      CALL DGETRS(  
     &   'N', NMO2, NMO1, WORK(LUBB), NMO2, 
     &   WORK(KIPIV), WORK(KB11A), NMO2,INFO
     &   )
C
C On exit the rhs is replaced by solution vectors
C
      KUABT = KB11A
      IF (INFO.NE.0) THEN
         CALL QUIT('UAUG: SGETRS SOLVER FAILED')
      END IF
      IF (DEBUG) THEN
         CALL DGEMM('N','N',NMO2,NMO1,NMO2,1.D0,
     &              WORK(KUBB),NMO2,
     &              WORK(KUABT),NMO2,0.D0,
     &              WORK(KTMP),NMO2)
         CALL STXOUTFUL(WORK(KTMP),NMO2,NMO1,'UAUG:UBB*UAB(T)')
      END IF
C
C Transpose back the solution
C
      KUAB = KTMP
      CALL STX_TRANSP(NMO2,NMO1,WORK(KUABT),WORK(KUAB))
      IF (DEBUG) THEN
         CALL STXOUTFUL(WORK(KUAB),NMO1,NMO2,'UAUG: UAB')
         CALL MEMCHK('UAG UAB',WORK,KA1)
      END IF
C
C   Transform trial MO:s with U(A,B)
C
C     CMO(:,NMO1+1:NMO) += MATMUL( CMO(:,1:NMO1), UAB )
C
      CALL DGEMM('N','N',NAO,NMO2,NMO1,1.D0,
     &           CMO(1,1),NAO,
     &           WORK(KUAB),NMO1,1.D0,
     &           CMO(1,NMO1+1),NAO)
      IF (DEBUG) THEN
         CALL STXOUTFUL(
     &      CMO(1,NMO1+1),NAO,NMO2,
     &      ' Final transformation of trial virtuals: UBB'
     &      )
      END IF
      IF (DEBUG) THEN
      !
      ! Test if matrix is unitary
      !
         CALL MEMGET2('REAL','U',KU,NMO*NMO,WORK,KFREE,LFREE)
         CALL MEMGET2('REAL','1',K1,NMO*NMO,WORK,KFREE,LFREE)
         IAA = KU
         IBA = KU + NMO1
         IAB = KU + NMO1*NMO
         IBB = KU + NMO1*NMO + NMO1
         CALL MCOPY(NMO1,NMO1,WORK(KA1),NMO1,WORK(IAA),NMO)
         CALL MCOPY(NMO1,NMO2,WORK(KUAB),NMO1,WORK(IAB),NMO)
         KB1=KTMP
         CALL STX_TRANSP(NMO1,NMO2,WORK(K1B),WORK(KB1))
         CALL MCOPY(NMO2,NMO1,WORK(KB1),NMO2,WORK(IBA),NMO)
         CALL MCOPY(NMO2,NMO2,WORK(KUBB),NMO2,WORK(IBB),NMO)
         CALL DGEMM('T','N',NMO,NMO,NMO,1.D0,
     &              WORK(KU),NMO,
     &              WORK(KU),NMO,0.D0,
     &              WORK(K1),NMO)
         CALL STXOUTFUL(WORK(KU),NMO,NMO,'UAUG: U')
         CALL STXOUTFUL(WORK(K1),NMO,NMO,'UAUG: U+U = 1 ?')
         CALL DGEMM('N','T',NMO,NMO,NMO,1.D0,
     &              WORK(KU),NMO,
     &              WORK(KU),NMO,0.D0,
     &              WORK(K1),NMO)
         CALL STXOUTFUL(WORK(K1),NMO,NMO,'UAUG: UU+ = 1 ?')
      END IF
      CALL MEMREL('UAUG:END',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('UAUG')
      END
C#######################################################################
      SUBROUTINE MATFUN(FUNC,N,A,FA,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
C
C  Matrix functions
C  On input: FUNC (implemented 'SQRT')
C            A  matrix
C            N  dimension
C  On output:
C            A  eigenvectors
C           FA  results
C  I/O  : 

      DIMENSION A(N,N), FA(N,N), WORK(LWORK)
      CHARACTER*(*) FUNC
      CALL QENTER('MATFUN')
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET2('REAL','EIG',KEIG,N,WORK,KFREE,LFREE)
      CALL DSYEV(
     &   'V','U',N,A,N,WORK(KEIG),WORK(KFREE),LFREE,INFO
     &   )
      IF (INFO.NE.0) THEN
         CALL QUIT('MATFUN: DSYEV DIAGONALIZATION FAILED')
      END IF
      call STXoutful(work(keig),n,1,'matfun:eigenvalues')
      IF (FUNC(1:4).EQ.'SQRT') THEN
         CALL DCOPY(N*N,A,1,FA,1)
         DO I=1,N
            EIG = WORK(KEIG+I-1)
C            print *,'eig ',eig
            IF (EIG.LT.0D0) THEN
               WRITE(0,*)
     &            '***WARNING***   MATFUN(SQRT):NEGATIVE EIGENVALUES ',
     &            EIG
               WRITE(0,*) 'ZERO WILL BE ASSUMED'
               WRITE(LUPRI,*)
     &            '***WARNING***   MATFUN(SQRT):NEGATIVE EIGENVALUES ',
     &            EIG
               WRITE(LUPRI,*) 'ZERO WILL BE ASSUMED'
               EIG = 0D0
            ELSE
               EIG = SQRT(EIG)
            END IF
            CALL DSCAL(N,EIG,FA(1,I),1)
         END DO
         CALL MEMGET2('REAL','TMP',KTMP,N*N,WORK,KFREE,LFREE)
         CALL DGEMM('N','T',N,N,N,1.D0,
     &              A,N,
     &              FA,N,0.D0,
     &              WORK(KTMP),N)
         CALL DCOPY(N*N,WORK(KTMP),1,FA,1)
      ELSE
CVC         CALL QUIT('MATFUN: '//FUNC//' NOT IMPLEMENTED')
         CALL QUIT('MATFUN: ...  NOT IMPLEMENTED')
      END IF
      CALL MEMREL('MATFUN',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('MATFUN')
      END
C#######################################################################
      SUBROUTINE STX_TRANSP(N,M,A,AT)
#include "implicit.h"
      DIMENSION A(N,M), AT(M,N)
C
C     AT = TRANSPOSE(A)
C
      DO I = 1,N
         DO J = 1,M
            AT(J,I) = A(I,J)
         END DO
      END DO
      END
C#######################################################################
      SUBROUTINE AUGCTL(WORK,LWORK)
#include "implicit.h"
      DIMENSION WORK(LWORK)
cvc
      dimension over(2485),sss(100),vgs(4900)
cvc
#include "inforb.h"
#include "stex.h"
#include "inftap.h"
#include "maxorb.h"
#include "maxash.h"
#include "infind.h"
#include "priunit.h"
      CHARACTER*8 STRCOU, STREXC, LABEL
      INTEGER IJLIST(2,NOPMAX)
      LOGICAL FNDLAB
      LOGICAL DEBUG/.TRUE./
C
      CALL QENTER('AUGCTL')
      CALL HEADER('Auger calculation',LUPRI)
      WRITE(LUPRI,*)'   Number of orbital pairs ', NAUGER
      IF (NAUGER.GT.NOPMAX) THEN
         WRITE(LUPRI,*)' Specified number of pairs ',NAUGER
         WRITE(LUPRI,*)' Maximum   number of pairs ',NOPMAX
         CALL QUIT('AUGER:NAUGER > NOPMAX - INCREASE NOPMAX in stex.h')
      END IF
      WRITE(LUPRI,*)'   Orbital pairs           '
C
C Set up index array for SCA_VW
C
      DO I=1,NAUGER
         WRITE(LUPRI,*) 
     &          '                           ', IAUGER(I),JAUGER(I)
         IJLIST(1,I) = IAUGER(I)
         IJLIST(2,I) = JAUGER(I)
      END DO
C
C Allocate 
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET2('REAL','VWMAT',KVWMAT,N2BASX*2*NAUGER,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','CMO',KCMO,NCMOT,WORK,KFREE,LFREE)
C
C Augment
C
      CALL STXAUG(WORK(KCMO),WORK,KFREE,LFREE)
C
c legge gli orbitali gs (aumentati)
c
        luaumo=-1
        CALL GPOPEN(luaumo,'AUMO','OLD','SEQUENTIAL','UNFORMATTED',
     &  0,.FALSE.)
!D    write(lupri,*)' read from ', luaumo,' the augmented gs MOs',nnbast 
c
      rewind luaumo
c
c     read(luaumo) nsym
      read(luaumo) 
c     write(lupri,*) ' nsym ',nsym
c     read(luaumo) (norb(i),i=1,nsym)
      read(luaumo) 
c     write(lupri,*) ' norb ',(norb(i),i=1,nsym)
c     read(luaumo) (nbas(i),i=1,nsym)
      read(luaumo) 
c     write(lupri,*) ' nbas ',(nbas(i),i=1,nsym)
c     nbast=0
c     norbt=0
c     ncmot=0
c     do i=1,nsym
c       icmo(i)=ncmot
c       nbast=nbast+nbas(i)
c       norbt=norbt+norb(i)
c       ncmot=ncmot+norb(i)*nbas(i)
c     enddo
c     write(LUPRI,*) ' nbast,norbt,ncmot ', nbast,norbt,ncmot
c     read(luaumo) (cent(i),i=1,nbast)
      read(luaumo) 
c     read(luaumo) (type(i),i=1,nbast)
      read(luaumo) 
      read(luaumo) (vgs(i),i=1,nnbast)
C
      CALL GPCLOSE (luaumo,'KEEP')
!D    write(lupri,*)' UNIT ',luaumo,' is closed'
c
      CALL HEADER('AUGCTL: Orbitali gs aumentati ',LUPRI)
      CALL PRORB(vgs,.true.,LUPRI)
C
c legge gli orbitali per calcolo Auger (aumentati)
c
      REWIND LUIT1
!D    write(LUPRI,*) 'LUIT1 =',LUIT1
      LABEL='SIRIFC2 '
      IF (RPATST) THEN
         write (LUPRI,*) 'RPA TEST: CMO unit matrix => AO integrals!'
         DO I = 1,NSYM
            CALL DUNIT(WORK(KCMO+ICMO(I)),NBAS(I))
         END DO
      ELSE IF (FNDLAB(LABEL,LUIT1)) THEN
            CALL HEADER('AUGCTL: Trovata la label',LUPRI)
            write (LUPRI,*) LABEL
         CALL READT(LUIT1,NCMOT,WORK(KCMO))
c
c test aggiunto da HJ giu99 per calcolo usando orbitals 2a1
c di hondo nesbet
c
c       call hj(WORK(KCMO),NORBT,NBAST)
c
cvc     IF (DEBUG) THEN
           CALL HEADER('AUGCTL: Orbitali per  matrici Auger',LUPRI)
cvc 
cvc   TRUE=  stampa solo orbitali occupati
cvc   FALSE= stampa tutti gli orbitali
        CALL PRORB(WORK(KCMO),.TRUE.,LUPRI)
cvc        CALL PRORB(WORK(KCMO),.FALSE.,LUPRI)
cvc     END IF
!D    write(LUPRI,*)' nauger, nbast =',nauger,nbast
      nom=0
      DO I=1,NAUGER
         WRITE(LUPRI,*) ' orbitale DALTON ', IAUGER(I)
         if(IAUGER(I).gt.nom) nom=IAUGER(I)
         j1=kcmo+(iauger(i)-1)*nbast+1
         j2=j1+nbast-1
         WRITE(LUPRI,99) (WORK(j),j=j1,j2)
         WRITE(LUPRI,*) ' orbitale DALTON', JAUGER(I)
         if(JAUGER(I).gt.nom) nom=JAUGER(I)
         j1=kcmo+(jauger(i)-1)*nbast+1
         j2=j1+nbast-1
         WRITE(LUPRI,99) (WORK(j),j=j1,j2)
      END DO
   99 format(5f12.5)
      read(5,*) iflagorb
      WRITE(LUPRI,*) ' iflagorb',iflagorb
      call quit('ERROR in AUGCTL')
 8888 continue
cvc
      ELSE
         CALL QUIT('AUGCTL : LABEL "'//LABEL//'" NOT FOUND')
      END IF
cvc
C
C Get AO overlap
C
      nnn=nbast*(nbast+1)/2
      nom=nom+1
!D    write(LUPRI,*) ' nbast,nnn,nom =',nbast,nnn,nom
      CALL RDONEL('OVERLAP ',.TRUE.,over,nnn)
!D    call header('matrice over AO ',LUPRI)
!D    call outpkb(over,NBAS,NSYM,1,LUPRI)
      do jjj=1,nom
        do iii=1,jjj
          sij=0.d0
          ml=0
          do m=1,nbast
            do l=1,m
              ml=ml+1
              il=(iii-1)*nbast+l
              im=(iii-1)*nbast+m
              mj=(jjj-1)*nbast+m
              lj=(jjj-1)*nbast+l
              sij=sij+vgs(im)*over(ml)*vgs(lj)
              if(m.ne.l) then
                sij=sij+vgs(il)*over(ml)*vgs(mj)
              endif
            enddo
          enddo
          ind=(jjj-1)*jjj/2+iii
          sss(ind)=sij
        enddo
      enddo
!D    write(LUPRI,*) ' overlap orbitali gs'
!D    call outpkb(sss,nom,NSYM,1,LUPRI)
      do jjj=1,nom
        do iii=1,jjj
          sij=0.d0
          ml=0
          do m=1,nbast
            do l=1,m
              ml=ml+1
              il=(iii-1)*nbast+l
              im=(iii-1)*nbast+m
              mj=(jjj-1)*nbast+m
              lj=(jjj-1)*nbast+l
              sij=sij+work(kcmo-1+im)*over(ml)*work(kcmo-1+lj)
              if(m.ne.l) then
                sij=sij+work(kcmo-1+il)*over(ml)*work(kcmo-1+mj)
              endif
            enddo
          enddo
          ind=(jjj-1)*jjj/2+iii
          sss(ind)=sij
        enddo
      enddo
!D    write(LUPRI,*) ' overlap orbitali auger'
!D    call outpkb(sss,nom,NSYM,1,LUPRI)
      do jjj=1,nom
        do iii=1,nom
          sij=0.d0
          ml=0
          do m=1,nbast
            do l=1,m
              ml=ml+1
              il=(iii-1)*nbast+l
              im=(iii-1)*nbast+m
              mj=(jjj-1)*nbast+m
              lj=(jjj-1)*nbast+l
              sij=sij+vgs(im)*over(ml)*work(kcmo-1+lj)
              if(m.ne.l) then
                sij=sij+vgs(il)*over(ml)*work(kcmo-1+mj)
              endif
            enddo
          enddo
          ind=(jjj-1)*nom+iii
          sss(ind)=sij
        enddo
      enddo
      write(LUPRI,*) ' overlap orbitali gs / orbitali auger'
      do jjj=1,nom
        write(LUPRI,999) jjj,(sss(iii),iii=jjj,nom*nom,nom)
      enddo
  999 format(i3,9f8.3)
      CALL QUIT('Hardwired STOP in AUGCTL')
cvc
C
C Build matrices
C
      CALL HEADER('AUGCTL: costruisce le matrici',LUPRI)
      CALL SCA_VW(
     &   NAUGER,IJLIST,WORK(KVWMAT),WORK(KCMO),WORK,KFREE,LFREE
     &   )
C
C Write to file AUGER
C
      STRCOU(1:8)='(IJ|**) '
      STREXC(1:8)='(I*|*J) '
      CALL MEMGET2('REAL','VWMO',KVWMO,N2ORBX,WORK,KFREE,LFREE)
CVC
      LUAUG=-1
!D    write(LUPRI,*) ' apre unita',luaug
CVC
      CALL GPOPEN(LUAUG,'AUGER','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &   0,.FALSE.)
CVC
      REWIND LUAUG
CVC
!D    write(LUPRI,*) 'scrive su unita',luaug,nauger,' coppie di matrici'
CVC
      WRITE(LUAUG)NAUGER,NORBT
CVC
      ICOFF = 0
      IXOFF = N2BASX
      DO N=1,NAUGER
         I=IJLIST(1,N)
         J=IJLIST(2,N)
!D    write(LUPRI,*) ' scrive su unita',luaug,' indici matrici ',i,j
         WRITE(LUAUG)I,J
         IF (IPRSTX.GT.10) THEN
            WRITE(STRCOU(2:2),'(I1)') I
            WRITE(STRCOU(3:3),'(I1)') J
            WRITE(STREXC(2:2),'(I1)') I
            WRITE(STREXC(6:6),'(I1)') J
           CALL STXOUTFUL(WORK(KVWMAT+ICOFF),NBAST,NBAST,STRCOU//'(AO)')
           CALL STXOUTFUL(WORK(KVWMAT+IXOFF),NBAST,NBAST,STREXC//'(AO)')
         END IF
         IF (STEXMO) THEN
            ISYMI=ISMO(I)
            ISYMJ=ISMO(J)
            IVWSYM=MULD2H(ISYMI,ISYMJ)
            CALL STXGENTRA(
     &         'SIRIFC1 ',WORK(KVWMAT+ICOFF),WORK(KVWMO),IVWSYM,
     &         WORK(KFREE),LFREE
     &         )
!D    write(LUPRI,*) ' scrive matrice exch1 (ij|**) base MO su unita',
!D   $           luaug,'   n2orbx =',n2orbx
            CALL WRITT(LUAUG,N2ORBX,WORK(KVWMO))
            IF (IPRSTX.GT.10) THEN
               CALL STXOUTFUL(WORK(KVWMO),NORBT,NORBT,STRCOU//'(MO)')
cvc
            else
cvc stampa ridotta alle prime 10 colonne
               CALL STXOUTFUL(WORK(KVWMO),NORBT,10,STRCOU//'(MO)')
cvc
            END IF
            CALL STXGENTRA(
     &         'SIRIFC1 ',WORK(KVWMAT+IXOFF),WORK(KVWMO),IVWSYM,
     &         WORK(KFREE),LFREE
     &         )
!D    write(LUPRI,*) ' scrive matrice exch2 (i*|*j)base MO su unita',
!D   $           luaug,'   n2orbx =',n2orbx
            CALL WRITT(LUAUG,N2ORBX,WORK(KVWMO))
            IF (IPRSTX.GT.10) THEN
               CALL STXOUTFUL(WORK(KVWMO),NORBT,NORBT,STREXC//'(MO)')
            else
cvc stampa ridotta alle prime 10 colonne
               CALL STXOUTFUL(WORK(KVWMO),NORBT,10,STREXC//'(MO)')
            END IF
         ELSE
!D    write(LUPRI,*) ' scrive matrice exch1 (ij|**) base AO su unita',
!D   &   luaug
            CALL WRITT(LUAUG,N2BASX,WORK(KVWMAT+ICOFF))
!D    write(LUPRI,*) ' scrive matrice exch2 (i*|*j)base AO su unita',
!D   &   luaug
            CALL WRITT(LUAUG,N2BASX,WORK(KVWMAT+IXOFF))
         END IF
         ICOFF = ICOFF + 2*N2BASX
         IXOFF = IXOFF + 2*N2BASX
      END DO
CVC
!D    write(lupri,*) ' chiude unita',luaug
      CALL GPCLOSE(LUAUG,'KEEP')
CVC
      CALL MEMREL('AUGCTL',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('AUGCTL')
      END
C#######################################################################
      SUBROUTINE AUGAB(CMO,VWMAT,A,B,VWMO,WORK,LWORK)
C
C Set up A and B matrices, check with E(2) explicitly (*LINEAR.ABCHK) 
C
C Input: (IP|**), (I*|P*) for all nonredundant rotations I->P
C
C Transform to matrices to MO and select (for non-diagonal elements)
C
C A(IP,JQ) = (PI|JQ) - 0.5*(PQ|JI)
C B(IP,JQ) = (PI|QJ) - 0.5*(PJ|QI)
C
C Note: the factor -0.5 is included in the exchange matrices
C
#include "implicit.h"
#include "inforb.h"
#include "maxorb.h"
#include "infvar.h"
#include "stex.h"
      DIMENSION CMO(NCMOT), A(NWOPT,NWOPT), B(NWOPT,NWOPT)
      DIMENSION VWMAT(NBAST,NBAST,2,NWOPT), VWMO(NORBT,NORBT,2)
      INTEGER I,J,P,Q,PRPSYM
      LOGICAL DEBUG
      PARAMETER (DH=0.5D0)
      CALL QENTER('AUGAB')
      DEBUG=IPRSTX.GT.10
      CALL DZERO(A,NWOPT*NWOPT)
      CALL DZERO(B,NWOPT*NWOPT)
      CALL DZERO(VWMO,2*N2ORBX)
      PRPSYM=1
C
C Loop over all nonredundant rotations
C
      DO IG=1,NWOPT
         I=JWOP(1,IG)
         P=JWOP(2,IG)
         IF (DEBUG) THEN
            PRINT *,'=================='
            PRINT *,' IG (I,P) ',IG,I,P
            PRINT *,'=================='
            CALL STXOUTFUL(VWMAT(1,1,1,IG),NBAST,NBAST,'VWMAT 1')
            CALL STXOUTFUL(VWMAT(1,1,2,IG),NBAST,NBAST,'VWMAT 2')
         END IF
C
C Transform (IP|**) and (I*|P*) to mo basis
C
         DO ISYM=1,NSYM
            JSYM=MULD2H(PRPSYM,ISYM)
            NBASI=NBAS(ISYM)
            NORBI=NORB(ISYM)
            ICMOI=ICMO(ISYM)
            print *,'isym nbasi norbi icmoi ',isym,nbasi,norbi,icmoi
            NBASJ=NBAS(JSYM)
            NORBJ=NORB(JSYM)
            ICMOJ=ICMO(JSYM)
            print *,'jsym nbasj norbj icmoj ',jsym,nbasj,norbj,icmoj
            CALL UTHV(
     &        CMO(1+ICMOI), VWMAT(1,1,1,IG), CMO(1+ICMOJ),
     &        ISYM, JSYM, NBASI, NBASJ, VWMO(1,1,1), WORK
     &        )
            CALL UTHV(
     &        CMO(1+ICMOI), VWMAT(1,1,2,IG), CMO(1+ICMOJ),
     &        ISYM, JSYM, NBASI, NBASJ, VWMO(1,1,2), WORK
     &        )
         END DO
         IF (DEBUG) THEN
            CALL STXOUTFUL(VWMO(1,1,1),NORBT,NORBT,'VWMO 1')
            CALL STXOUTFUL(VWMO(1,1,2),NORBT,NORBT,'VWMO 2')
         END IF
C
C  Select the matrix elements that contribute to RPA matrices A and B 
C
C  this scaling fixes the RPA comparison
C
         call dscal(n2orbx, 4.0d0,vwmo(1,1,1),1)
         call dscal(n2orbx, 4.0d0,vwmo(1,1,2),1)
         DO JG=1,IG
            J = JWOP(1,JG)
            Q = JWOP(2,JG)
            IF (DEBUG) PRINT *,' JG (J,Q) ',JG,J,Q
            A(IG,JG) = VWMO(J,Q,1) + VWMO(J,Q,2)
            A(JG,IG) = A(IG,JG)
            B(IG,JG) = - VWMO(Q,J,1) - VWMO(Q,J,2)
            B(JG,IG) = B(IG,JG)
         END DO
         IF (DEBUG) THEN
            CALL STXOUTFUL(A,NWOPT,NWOPT,'Partial A matrix')
            CALL STXOUTFUL(B,NWOPT,NWOPT,'Partial B matrix')
         END IF
      END DO
      CALL QEXIT('AUGAB')
      END
C#######################################################################
      SUBROUTINE STXGENTRA(LABEL,PRPAO,PRPMO,PRPSYM,WORK,LWORK)
#include "implicit.h"
      REAL*8  PRPAO(*), PRPMO(*), WORK(*)
      INTEGER  LWORK, PRPSYM
      CHARACTER*8 LABEL
#include "inftap.h"
#include "inforb.h"
#include "stex.h"
      LOGICAL FNDLAB
      EXTERNAL FNDLAB
      LOGICAL DEBUG
      CALL QENTER('STXGENTRA')
C
C Transform  PRPAO by mo:s defined by LABEL: output PRPMO
C
      DEBUG=IPRSTX.GT.10
      KFREE=1
      LFREE=LWORK
C     OPEN(LUIT1,FORM='UNFORMATTED',STATUS='OLD')
      REWIND LUIT1
      if (debug)
     &   write(LUPRI,*)'STXGENTRA fndlab: label luit1 ',label,luit1
      IF (FNDLAB(LABEL,LUIT1)) THEN
         CALL MEMGET2('REAL','CMO',KCMO,NCMOT,WORK,KFREE,LFREE)
         CALL READT(LUIT1,NCMOT,WORK(KCMO))
         if (debug) call prorb(work(kcmo),.false.,LUPRI)
      ELSE
         CALL QUIT('STXGENTRA: LABEL "'//LABEL//'" NOT FOUND')
      END IF
C     CLOSE(LUIT1)
C
C Form square and transform
C
      CALL DZERO(PRPMO,N2ORBX)
      DO ISYM=1,NSYM
         JSYM=MULD2H(PRPSYM,ISYM)
         NBASI=NBAS(ISYM)
         NORBI=NORB(ISYM)
         ICMOI=ICMO(ISYM)
         NBASJ=NBAS(JSYM)
         NORBJ=NORB(JSYM)
         ICMOJ=ICMO(JSYM)
         IF (DEBUG) THEN
            write(LUPRI,*)
     &         'isym jsym nbasi nbasj norbi norbj icmoi icmoj ',
     &          isym,jsym,nbasi,nbasj,norbi,norbj,icmoi,icmoj
         END IF
         CALL UTHV(
     &      WORK(KCMO+ICMOI), PRPAO, WORK(KCMO+ICMOJ),
     &      ISYM, JSYM, NBASI, NBASJ, PRPMO, WORK(KFREE)
     &      )
      END DO
      IF (DEBUG) THEN
         call STXoutful(prpmo,norbt,norbt,'STXgentra after ')
      END IF
      CALL MEMREL('STXGENTRA',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('STXGENTRA')
      END
! -- end of sirius/sirstex.F --
