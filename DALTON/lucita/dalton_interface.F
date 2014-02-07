!
!...   Copyright (c) 2013 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2013 (2013), see http://daltonprogram.org"
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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig, 2001        *
* Dalton adaptation by Jeppe Olsen, Hans Joergen Aa. Jensen and       *
* Stefan Knecht May 08 - Dec 10.                                      *
*                                                                     *
***********************************************************************

      SUBROUTINE DMPINT(LUINT)
*
* Dump integrals in WORK(KINT1),WORK(KINT2) on file LUINT
*
      use lucita_energy_types
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "glbbas.inc"
#include "wrkspc.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      CALL REWINO(LUINT)
*.1 : One-electron integrals
      WRITE(LUINT,'(E22.15)')
     &     (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
      WRITE(LUINT,'(E22.15)')
     &     (WORK(KINT2-1+INT2),INT2=1,NINT2)
*.3 : Core energy
      WRITE(LUINT,'(E22.15)') ECORE_ORIG
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_CMOAO(CMO)
*
* Obtain AO-MO transformation matrix
*
* Jeppe Olsen, November 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "clunit.inc"
#include "orbinp.inc"
*. Output
      DIMENSION CMO(*)

!     write(lupri,*) 'GET_CMOAO: ENVIRO = ',ENVIRO !hjaaj DEBUG
      IF(ENVIRO(1:6).EQ.'DALTON') THEN
        CALL GET_CMOAO_DALTON(CMO,NMOS_ENV(1),NAOS_ENV(1),NSMOB)
      ELSE IF(ENVIRO(1:6).EQ.'MOLCAS') THEN
*.
        CALL QUIT('MOLCAS environment not available in this version')
      ELSE IF(ENVIRO(1:5).EQ.'LUCIA' ) THEN
*. Read in from LUCIA 1e file : unit 91
        LU91 = 91
        CALL GET_CMOAO_LUCIA(CMO,NMOS_ENV,NAOS_ENV,LU91)
      ELSE IF(ENVIRO(1:4).EQ.'NONE') THEN
        WRITE(lupri,*) ' GET_CMOAO, Warning : Called with ENVIRO = NONE'
        WRITE(lupri,*) ' No coefficients read in '
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_CMOAO_DALTON(CMO,NBAS,NMO,NSM)
*
* Obtain MO-AO expansion matrix from SIRIUS/DALTON file SIRGEOM
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
*. Input
      INTEGER NBAS(*), NMO(*)
*. Output
      DIMENSION CMO(*)
*

!     write(lupri,*) 'GET_CMOAO_DALTON: NSM  = ',NSM         !hjaaj DEBUG
!     write(lupri,*) 'GET_CMOAO_DALTON: NBAS = ',NBAS(1:NSM) !hjaaj DEBUG
!     write(lupri,*) 'GET_CMOAO_DALTON: NMO  = ',NMO (1:NSM) !hjaaj DEBUG
      NTEST  = 0
      ITAP30 = 16
      OPEN(ITAP30,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND ITAP30
      CALL MOLLAB('TRCCINT ',ITAP30,6)
*. Skip record containing dimensions of orbitals
      READ (ITAP30) NSYM
*. And skip record containing eigenvalues etc
      READ (ITAP30)
C     READ (ITAP30) NSYMHF,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYMHF),
C    *              (NORB(I),I=1,NSYMHF),(NBAS(I),I=1,NSYMHF),
C    *              POTNUC,EMCSCF
C
C
C     READ (ITAP30) (WRK(KEIGVL+I-1),I=1,NORBT),
C    *              (IWRK(KEIGSY+I-1),I=1,NORBT)
*. And then the MO-AO expansion matrix
      NCOEF = 0
      DO ISM = 1, NSM
        NCOEF = NCOEF + NMO(ISM)*NBAS(ISM)
      END DO
      READ (ITAP30) (CMO(I),I=1,NCOEF)
      CLOSE(ITAP30,STATUS='KEEP')
*
      NTEST = 000 !
!     NTEST = 100 ! hjaaj DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) '  MO - AO expansion matrix '
        WRITE(lupri,*) '============================='
        WRITE(lupri,*)
        CALL APRBLM2(CMO,NBAS,NMO,NSM,0)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_CMOAO_LUCIA(CMO,NMOS,NAOS,LUH)
*
* Obtain CMOAO expansion matrix from LUCIA formatted file LUH
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
*. Input
      INTEGER NMOS(*),NAOS(*)
*. Output
      DIMENSION CMO(*)
*
* Structure of file
* 1 : Number of syms
* 2 : NMO's per sym
* 3 : NAO's per SYM
* 4 : Number of elements in CMOAO
* Note : CMOAO and property integrals written in form
*     given by ONEEL_MAT_DISC
*
* Jeppe Olsen, Feb. 98
*
      WRITE(lupri,*)  ' GET_CMOAO_LUCIA, LUH = ', LUH
      CALL REWINO(LUH)
*. skip Number of orbital symmetries
      READ(LUH,*)
*. skip Number of MO's per symmetry
      READ(LUH,*)
*. skip Number of AO's per symmetry
      READ(LUH,*)
*. skip read Length of CMO-AO expansion
      READ(LUH,*)
*. read CMO-AO expansion matrix
      CALL ONEEL_MAT_DISC(CMO,1,NSMOB,NAOS,NMOS,LUH,1)
C          ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) ' MO-AO transformation read in '
        CALL PRHONE(CMO,NMOS,1,NSMOB,0)
C            PRHONE(C,NFUNC,M,NSM,IPACK)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_CMOMO(CMOMO)
*
* Obtain MO-MO transformation matrix CMOMO for transforming to
* final set of orbitals
*
* Output matrix CMOMO is returned in symmetry packed form
*
*. Density matrix is assumed in place
*
* Type of final orbitals is provided by the keyword
* keywords ITRACI_CR, ITRACI_CN
*
* ITRACI_CR : COMP => Rotate all orbitals
*             REST => Rotalte only inside orbital subspaces
*
* ITRACI_CN : NATU => Transform to natural orbitals
* ITRACI_CR : CANO => Transform to canonical orbitals
*
* Jeppe Olsen, February 1998 ( from FINMO)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
      INTEGER*8 KMAT1, KMAT2, KMAT3, KMAT4
      ! for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Output
      DIMENSION CMOMO(*)
*
      NTEST = 100
      IF(NTEST.GE.1) THEN
        WRITE(lupri,*)
        WRITE(lupri,*) ' ===================='
        WRITE(lupri,*) ' GET_CMOMO in action'
        WRITE(lupri,*) ' ===================='
        WRITE(lupri,*)
      END IF

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETCMO')
      CALL MEMMAN(KMAT1,NTOOB**2,'ADDL  ',2,'MAT1  ')
      CALL MEMMAN(KMAT2,NTOOB**2,'ADDL  ',2,'MAT2  ')
      CALL MEMMAN(KMAT3,NTOOB**2,'ADDL  ',2,'MAT3  ')
      CALL MEMMAN(KMAT4,2*NTOOB**2,'ADDL  ',2,'MAT4  ')
*
*. Matrix defining final orbitals
*
      IF(ITRACI_CN(1:4).EQ.'CANO' ) THEN
*. Construct FI+FA in WORK(KMAT1)
        CALL COPVEC(WORK(KINT1O),WORK(KMAT1),NINT1)
        CALL FIFAM(WORK(KMAT1))
      ELSE IF(ITRACI_CN(1:4).EQ.'NATU' ) THEN
*. Symmetry order density matrix
        CALL TYPE_TO_SYM_REO_MAT(WORK(KRHO1),WORK(KMAT2))
*. Pack to triangular form
        CALL TRIPAK_BLKM(WORK(KMAT2),WORK(KMAT1),1,NTOOBS,NSMOB)
*. multiply by minus one to get natural orbitals
*. with largest occupations first
        ONEM = -1.0D0
        LDIM = 0
        DO ISM = 1, NSMOB
          LDIM = LDIM + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        END DO
        CALL SCALVE(WORK(KMAT1),ONEM,LDIM)
        IF(NTEST.GE.100) THEN
          WRITE(lupri,*) ' Packed density matrix ( times - 1 )'
          CALL APRBLM2(WORK(KMAT1),NACOBS,NACOBS,NSMOB,1)
        END IF
      END IF
*
* Diagonalize
*
      IF(ITRACI_CR(1:4).EQ.'REST') THEN
*. Diagonalize symmetry-type blocks
        CALL DIAG_BLKS(WORK(KMAT1),CMOMO,NGSOB,NTOOBS,MXPOBS,
     &                 NSMOB,NGAS,WORK(KMAT3),WORK(KMAT4))
*. Reorder to assure max diag dominance
        IREO = 1
        IF(IREO.NE.0) THEN
          WRITE(lupri,*) ' CMOMO reordered to assure max. diag. dom.'
          DO ISM = 1, NSMOB
            IF(ISM.EQ.1) THEN
              IOFF = 1
            ELSE
              IOFF = IOFF + NTOOBS(ISM-1)**2
            END IF
            L  = NTOOBS(ISM)
            CALL GET_DIAG_DOM(CMOMO(IOFF),WORK(KMAT1),L,WORK(KMAT2))
            CALL COPVEC(WORK(KMAT1),CMOMO(IOFF),L*L)
          END DO
        END IF
      ELSE IF (ITRACI_CR(1:4).EQ.'COMP') THEN
*. Diagonalize symmetry blocks
        CALL DIAG_BLKS(WORK(KMAT1),CMOMO,NACOBS,NTOOBS,MXPOBS,
     &                 NSMOB,1,WORK(KMAT3),WORK(KMAT4))
*. Reorder to assure max diag dominance
        IREO = 1
        IF(IREO.NE.0) THEN
          WRITE(lupri,*) ' CMOMO reordered to assure max. diag. dom.'
          DO ISM = 1, NSMOB
            IF(ISM.EQ.1) THEN
              IOFF = 1
            ELSE
              IOFF = IOFF + NTOOBS(ISM-1)**2
            END IF
            L  = NTOOBS(ISM)
            CALL GET_DIAG_DOM(CMOMO(IOFF),WORK(KMAT1),L,WORK(KMAT2))
            CALL COPVEC(WORK(KMAT1),CMOMO(IOFF),L*L)
          END DO
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
         WRITE(lupri,*) ' Output set of MO''s '
         CALL APRBLM2(CMOMO,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETCMO')
      RETURN
      END
***********************************************************************
*
* Obtain property integrals with LABEL LABEL from LU91,
* LUCIA format
*
* Jeppe Olsen, Feb.98

      SUBROUTINE GET_H1AO(LABEL,H1AO,IHSM,NBAS)
*
* Obtain 1 electron integrals with label LABEL
*
* Jeppe Olsen, Feb.98
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
      INTEGER*8 KLSCR
      ! for addressing of WORK
*
#include "mxpdim.inc"
#include "crun.inc"
#include "orbinp.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
*
      CHARACTER*8 LABEL
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_H1A')
*
      write(lupri,*) 'GET_H1AO: ENVIRO = ',ENVIRO !hjaaj DEBUG
      IF(ENVIRO(1:6).EQ.'DALTON') THEN
        LSCR = NTOOB**2
        CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'GTH1SC')
        CALL GET_H1AO_DALTON(LABEL,H1AO,IHSM,WORK(KLSCR),NBAS,NSMOB)
C            GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS,NSM)
      ELSE IF (ENVIRO(1:5).EQ.'LUCIA') THEN
        LU91 = 91
        CALL GET_H1AO_LUCIA(LABEL,H1AO,LU91)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT_H1A')
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS,NSM)
*
*. Obtain one-electron integrals in ao basis from dalton
*
* Label of integrals LABEL from FILE AORPROPER
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
*. Input
      CHARACTER*8 LABEL
      DIMENSION NBAS(*)
#include "multd2h.inc"
*. output
      DIMENSION H1AO(*)
*. Scratch
      DIMENSION SCR(*)
*
      LOGICAL FNDLAB
*
      write(lupri,*) 'GET_H1AO_DALTON: LABEL = ',LABEL !hjaaj DEBUG
      NTEST =   02
      IF(NTEST.GE.2) THEN
        WRITE(lupri,*) ' Fetching one-electron integrals with Label ',
     &  LABEL
        WRITE(lupri,*) ' IHSM NSM', IHSM,NSM
      END IF
*
*. Number of elements in : Complete lower half array
*                          Symmetry restricted complete matrix
*                          Symmetry restricted lower half matrix
*-- I am not completely sure about the input format of the integrals
      NBAST = 0
      DO ISM = 1, NSM
       NBAST = NBAST + NBAS(ISM)
      END DO
      NINT01 = NBAST*(NBAST+1)/2
C     write(lupri,*) ' IHSM = ', IHSM
*
      NINT10 = 0
      DO IRSM = 1, NSM
       ICSM = MULTD2H(IHSM,IRSM)
       NINT10 = NINT10 + NBAS(IRSM)*NBAS(ICSM)
      END DO
*
      NINT11 = 0
      DO IRSM = 1, NSM
       ICSM = MULTD2H(IHSM,IRSM)
       IF(IRSM.GT.ICSM) THEN
        NINT11 = NINT11 + NBAS(IRSM)*NBAS(ICSM)
       ELSE IF(IRSM.EQ.ICSM) THEN
        NINT11 = NINT11 + NBAS(IRSM)*(NBAS(IRSM)+1)/2
       END IF
      END DO
*
*. Read in integrals, assumed in complete lower half format
*
         LUPRP = 15
         OPEN (LUPRP,STATUS='OLD',FORM='UNFORMATTED',FILE='AOPROPER')
         REWIND (LUPRP)
         IF (FNDLAB(LABEL,LUPRP)) THEN
C           write(lupri,*) ' Label obtained'
            READ(LUPRP) (SCR(I),I=1,NINT01)
C           write(lupri,*) 'integrals readin'
C           call prsym(scr,NBAST)
C           CALL READT(LUPRP,NBAST*(NBAST+1)/2,WRK(KSCR2))
         ELSE
            WRITE(lupri,*)'Property label "',LABEL,'" not found on file'
            Call Abend2( 'Wrong input or integrals not generated' )
         ENDIF
        CLOSE(LUPRP,STATUS='KEEP')
*
C        WRITE(lupri,*) ' Number of symmetry apdapted integrals',NINT10
*
*. Transfer integrals to symmetry adapted form, complete form
*
         IBINT = 1
*. Loop over symmetry blocks
         DO IRSM = 1, NSM
           ICSM = MULTD2H(IHSM,IRSM)
           NR = NBAS(IRSM)
           NC = NBAS(ICSM)
*. Offsets
           IBR = 1
           DO ISM = 1, IRSM - 1
             IBR = IBR + NBAS(ISM)
           END DO
           IBC = 1
           DO ISM = 1, ICSM - 1
             IBC = IBC + NBAS(ISM)
           END DO
*. Complete block, stored in usual column wise fashion
           DO ICORB = 1, NC
             DO IRORB = 1, NR
               ICABS = IBC + ICORB -1
               IRABS = IBR + IRORB -1
               ICRMX = MAX(ICABS,IRABS)
               ICRMN = MIN(ICABS,IRABS)
               H1AO(IBINT-1 + (ICORB-1)*NR+IRORB) =
     &         SCR(ICRMX*(ICRMX-1)/2+ICRMN)
             END DO
           END DO
           IBINT = IBINT + NR*NC
         END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) ' One-electron integrals obtained from AOPROPER'
        CALL PRSYM(SCR,NBAST)
*
        WRITE(lupri,*) ' One-electron integrals in packed form'
        CALL PRHONE(H1AO,NBAS,IHSM,NSM,0)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_H1AO_LUCIA(LABEL,H1,LUH)
*
*
* Obtain property integrals with LABEL LABEL from LU91,
* LUCIA format
*
* Jeppe Olsen, Feb.98
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
#include "wrkspc.inc"
C     CHARACTER*1 XYZ(3)
C     DATA XYZ/'X','Y','Z'/
      CHARACTER*8 LABEL, LABEL2, LABELX
*. Output
      DIMENSION H1(*)
*
* Structure of file
* 1 : Number of syms
* 2 : NMO's per sym
* 3 : NAO's per SYM
* 4 : Number of elements in CMOAO
* 4 : CMOAO-expansion matrix (in symmetry packed form)
* 5 : Number of property AO lists
*     Loop over number of properties
*     Label, offset and length of each proprty list
*
*     Property integrals for prop1,prop2 ...
*
* Note : CMOAO and property integrals written in form
*     given by ONEEL_MAT_DISC
*
* Jeppe Olsen, Feb. 98
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETH1A')
*
*. DIPOLE => DIPLEN
      IF(LABEL(1:6).EQ.'DIPOLE') THEN
        LABELX = 'DIPLEN  '
      ELSE
        LABELX = LABEL
      END IF
*
      CALL REWINO(LUH)
*. Skip Number of orbital symmetries
      READ (LUH,*) NSMOB
*. Skip Number of MO's per symmetry
      READ (LUH,*) (NMOS_ENV(ISM),ISM=1,NSMOB)
*. Skip Number of AO's per symmetry
      READ (LUH,*) (NAOS_ENV(ISM),ISM=1,NSMOB)
*. Length of CMO-AO expansion
      READ(LUH,*) LENGTH
*. And skip
      DO IJ = 1, LENGTH
        READ(LUH,'(E22.15)')
      END DO
*. Total number of properties ( 3 for each rank1, 6 for each rank 2)
      READ(LUH,*) NPROP_COMP
      IFOUND = 0
      WRITE(lupri,*) ' NPROP_COMP = ', NPROP_COMP
      DO IPROP_COM = 1, NPROP_COMP
        READ(LUH,'(A,I6,I6)') LABEL2,IOFF,LENGTH
        IF(LABEL2.EQ.LABELX) THEN
          IOFFA = IOFF
          LENGTHA = LENGTH
          IFOUND = 1
        END IF
      END DO
      IF(IFOUND.EQ.0) THEN
        WRITE(lupri,*) ' Label not found on file 91'
        WRITE(lupri,'(2A)' ) ' Label = ', LABELX
        Call Abend2( ' Label not found on file 91' )
      END IF
*. Skip to start of integrals
      WRITE(lupri,*) ' IOFFA, LENGTHA ', IOFFA,LENGTHA
      DO IJ = 1, IOFFA - 1
        READ(LUH,*)
      END DO
*. and read
      CALL SYM_FOR_OP(LABEL,IXYZSYM,IOPSM)
      CALL ONEEL_MAT_DISC(H1,IOPSM,NSMOB,
     &                    NAOS_ENV,NAOS_ENV,LUH,1)
C          ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) ' Property integrals read in '
        CALL PRHONE(H1,NAOS_ENV,IOPSM,NSMOB,0)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETH1A')
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_ORB_DIM_ENV(ECORE_ENV)
*
* Obtain number of orbitals and basis functions from the
* programming environment.
* results stored in NAOS_ENV, NMOS_ENV
*
* Obtain environments CORE energy, ECORE_ENV
*
* Jeppe Olsen, December 97
*
#ifdef VAR_MPI
      use dalton_mpi
#endif
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "parluci.h"
! potnuc on infinp.h could also deliver ecore_env
!
      NTEST = 000
!     NTEST = 9999 ! debug

      call izero(naos_env,MXPOBS)
      call izero(nmos_env,MXPOBS)
      if(luci_myproc .eq. luci_master)then
        CALL GETOBS_DALTON(ECORE_ENV,NAOS_ENV,NMOS_ENV,NSMOB)
      end if
#ifdef VAR_MPI
      if(luci_nmproc .gt. 1)then
        call dalton_mpi_bcast(nsmob,     luci_master, mpi_comm_world)
        call dalton_mpi_bcast(ecore_env, luci_master, mpi_comm_world)
        call dalton_mpi_bcast(naos_env,  luci_master, mpi_comm_world)
        call dalton_mpi_bcast(nmos_env,  luci_master, mpi_comm_world)
      end if
#endif
!
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) ' From GET_ORB_FROM_ENV : '
        WRITE(lupri,*) ' ======================='
        WRITE(lupri,*) ' NSMOB   ',NSMOB
        WRITE(lupri,*) ' NAOS_ENV'
        CALL IWRTMA(NAOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(lupri,*) ' NMOS_ENV'
        CALL IWRTMA(NMOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(lupri,*) ' ECORE_ENV=', ECORE_ENV
      END IF
!
      END
***********************************************************************
*
      SUBROUTINE GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW)
*
*. Obtain Property integrals in MO basis for operator with
*  label LABEL.
*
* If ILOW = 1, only the elements below the diagonal are
* obtained.
*
* Jeppe Olsen, June 1997
*              September 97 : ILOW added
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
*. Input
      DIMENSION NMO(*),NBAS(*)
#include "multd2h.inc"
*. Output
      DIMENSION H(*)
*. Scratch
      DIMENSION SCR(*)
*. Scratch should atleaest be of length  **
*
      NTEST = 000
*. Integrals in AO basis, neglect symmetry
      NBAST = 0
      NMOT = 0
      DO ISM = 1, NSM
        NBAST = NBAST + NBAS(ISM)
        NMOT  = NMOT  + NMO(ISM)
      END DO
C?    WRITE(lupri,*) ' Total number of basis functions ',NBAST
      LINTMX = NBAST*NBAST
*
      KLH1AO = 1
      KLFREE = KLH1AO + LINTMX
*
      KLC = KLFREE
      KLFREE = KLC + LINTMX
*
      KLSCR = KLFREE
*. Currently only DALTON route is working
      IDALTON = 1
      IF(IDALTON.EQ.1) THEN
C?      WRITE(lupri,*) ' Dalton route in action'
*. Obtain AO property integrals
C            GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS)
C            GET_H1AO(LABEL,H1AO,IHSM,NBAS)
        CALL GET_H1AO(LABEL,SCR(KLH1AO),IHSM,NBAS)
C       CALL GET_H1AO_DALTON(LABEL,SCR(KLH1AO),IHSM,
C    &       SCR(KLSCR),NBAS,NSM)
*. Obtain MO-AO transformation matrix
        CALL GET_CMOAO(SCR(KLC))
*. Transform from AO to MO basis
C            TRAH1(NBAS,NORB,NSYM,HAO,C,HMO,IHSM,SCR)
        CALL TRAH1(NBAS,NMO,NSM,SCR(KLH1AO),SCR(KLC),H,IHSM,
     &             SCR(KLSCR))
      END IF
*
      IF(NTEST .GE. 100 ) THEN
        WRITE(lupri,*) 'electron integrals in MO basis, full format '
        CALL PRHONE(H,NMO,IHSM,NSM,0)
      END IF
      IF(ILOW.EQ.1) THEN
*. Complete to lower half form
        IOFF_IN = 1
        IOFF_OUT = 1
        DO ISM = 1, NSM
          JSM = MULTD2H(ISM,IHSM)
          IF(ISM.EQ.JSM) THEN
*. Copy lower half
            LDIM = NMO(ISM)
            NELMNT_IN = LDIM * LDIM
            NELMNT_OUT = LDIM * (LDIM + 1)/2
            CALL COPVEC(H(IOFF_IN),SCR(KLSCR),NELMNT_IN)
            SIGN = 1.0D0
            CALL TRIPAK_LUCI(SCR(KLSCR),H(IOFF_OUT),1,LDIM,LDIM,SIGN)
            IOFF_IN = IOFF_IN + NELMNT_IN
            IOFF_OUT = IOFF_OUT + NELMNT_OUT
          ELSE IF(ISM.LT.JSM) THEN
*. Just skip block in input matrix
            LIDIM = NMO(ISM)
            LJDIM = NMO(JSM)
            IOFF_IN = IOFF_IN + LIDIM*LJDIM
          ELSE IF(ISM.GT.JSM) THEN
*. Copy block to block
            LIDIM = NMO(ISM)
            LJDIM = NMO(JSM)
            NELMNT = LIDIM*LJDIM
C           CALL TRPMAT(H(IOFF_IN),LIDIM,LJDIM,H(IOFF_OUT))
            CALL COPVEC(H(IOFF_IN),H(IOFF_OUT),NELMNT)
            IOFF_IN = IOFF_IN + NELMNT
            IOFF_OUT = IOFF_OUT + NELMNT
          END IF
        END DO
      END IF
*. The one-electron integrals reside in a NMOT X NMOT matrix.
*. Zero trivial integrals
      IF(ILOW.EQ.1) THEN
        NELMNT = IOFF_OUT-1
      ELSE
        LENGTH = 0
        DO ISM = 1, NSM
          JSM = MULTD2H(ISM,IHSM)
          NELMNT = NELMNT + NMO(ISM)*NMO(JSM)
        END DO
        IFREE = NELMNT + 1
      END IF
C?    WRITE(lupri,*) ' GET_PROP : NELMNT= ', NELMNT
      ZERO = 0.0D0
      NZERO = NMOT*NMOT - NELMNT
      IFREE = NELMNT + 1
      CALL SETVEC(H(IFREE),ZERO,NZERO)

      IF(NTEST .GE. 50 ) THEN
        WRITE(lupri,*) 'electron integrals in MO basis '
        CALL PRHONE(H,NMO,IHSM,NSM,ILOW)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GETH1(H,ISM,ITP,JSM,JTP)
*
* One-electron integrals over orbitals belonging to
* given OS class
*
*
* The orbital symmetries  are used to obtain the total
* symmetry of the one-electron integrals.
* It is therefore assumed that ISM, JSM represents 
*   a correct symmetry block
* of the integrals
*
* Jeppe Olsen, Version of fall 97
*              Summer of 98 : CC options added
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
#include "mxpdim.inc"
#include "wrkspc.inc"
*.Global pointers
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cc_exc.inc"
*.Output
      DIMENSION H(*)
*
      NI = NOBPTS(ITP,ISM)
      NJ = NOBPTS(JTP,JSM)
*
      IF(ICC_EXC.EQ.0) THEN
*
* Normal one-electron integrals
*
        IJ = 0
        DO J = 1, NJ
          DO I = 1, NI
            IJ = IJ+1
            H(IJ) = GETH1E(I,ITP,ISM,J,JTP,JSM)
          END DO
        END DO
      ELSE
*
* Single excitation coefficients dressed up as integrals
* taken from KCC
C           GET_SX_BLK(HBLK,H,IGAS,ISM,JGAS,JSM)
*. Note : WORK(KCC1) not perfect choice
!      CALL GET_SX_BLK(H,WORK(KCC1),ITP,ISM,JTP,JSM)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(lupri,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
        CALL WRTMT_LU(H,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END
***********************************************************************
      FUNCTION GETH1E(IORB,ITP,ISM,JORB,JTP,JSM)
*
* One-electron integral for active
* orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
*
* The orbital symmetries are used to obtain the
* total symmetry of the operator
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "multd2h.inc"
#include "intform.inc"
*
      IJSM = MULTD2H(ISM,JSM)
      IF(IH1FORM.EQ.1) THEN
*. Normal integrals, lower triangular packed
        IF(IJSM.EQ.1) THEN
          GETH1E =
     &    GTH1ES(IREOTS,WORK(KPINT1),WORK(KINT1),IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        ELSE
          GETH1E =
     &    GTH1ES(IREOTS,WORK(KPGINT1(IJSM)),WORK(KINT1),IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        END IF
      ELSE
*. Integrals are in full blocked form
        GETH1E =
     &  GTH1ES(IREOTS,WORK(KPGINT1A(IJSM)),WORK(KINT1),IBSO,MXPNGAS,
     &         IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,0)
      END IF
*
      RETURN
      END
***********************************************************************
      FUNCTION GETH1I(IORB,JORB)
*
* Obtain one -electron integral H(IORB,JOB)
*
* Interface from EXPHAM to LUCIA
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "priunit.h"
*
#include "mxpdim.inc"
#include "orbinp.inc"
*
      ISM = ISMFTO(IORB)
      ITP = ITPFSO(IREOTS(IORB))
      IREL = IORB - IOBPTS(ITP,ISM) + 1
*
      JSM = ISMFTO(JORB)
      JTP = ITPFSO(IREOTS(JORB))
      JREL = JORB - IOBPTS(JTP,JSM) + 1
*
      GETH1I = GETH1E(IREL,ITP,ISM,JREL,JTP,JSM)
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(lupri,*) ' GETH1I : IORB JORB ', IORB, JORB
        WRITE(lupri,*) ' ISM ITP IREL ', ISM,ITP,IREL
        WRITE(lupri,*) ' JSM JTP JREL ', JSM,JTP,JREL
        WRITE(lupri,*) ' GETH1I = ', GETH1I
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GETINCN2(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,I2INDX,
     &                  ICOUL)
*
* Obtain integrals
*
*     ICOUL = 0 :
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 :
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 2 :  XINT(IL,JK) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0
*
*
* Version for integrals stored in INTLST
*
* If type equals zero, all integrals of given type are fetched
* ( added aug8, 98)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
*. Pair of orbital indices ( symmetry ordered ) => address in symmetry packed
*. matrix
      Dimension I2INDX(*)
*.Output
      DIMENSION XINT(*)
*. Local scratch
      DIMENSION IJARR(MXPORB)
*
      IF(ITP.GE.1) THEN
        iOrb=NOBPTS(ITP,ISM)
      ELSE
        IORB = NTOOBS(ISM)
      END IF
*
      IF(JTP.GE.1) THEN
        jOrb=NOBPTS(JTP,JSM)
      ELSE
        JORB = NTOOBS(JSM)
      END IF
*
      IF(KTP.GE.1) THEN
        kOrb=NOBPTS(KTP,KSM)
      ELSE
        KORB = NTOOBS(KSM)
      END IF
*
      IF(LTP.GE.1) THEN
        lOrb=NOBPTS(LTP,LSM)
      ELSE
        LORB = NTOOBS(LSM)
      END IF
*
*. Offsets relative to start of all orbitals, symmetry ordered
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO

*
*     Collect Coulomb terms
*
      ijblk = max(ism,jsm)*(max(ism,jsm)-1)/2 + min(ism,jsm)
      klblk = max(ksm,lsm)*(max(ksm,lsm)-1)/2 + min(ksm,lsm)
*
      IF(IJBLK.GT.KLBLK) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.EQ.KLBLK) THEN
       IJRELKL = 0
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.LT.KLBLK) THEN
       IJRELKL = -1
       IBLOFF = IJKLOF(MAX(KSM,LSM),MIN(KSM,LSM),MAX(ISM,JSM))
      END IF
*
      itOrb=NTOOBS(iSm)
      jtOrb=NTOOBS(jSm)
      ktOrb=NTOOBS(kSm)
      ltOrb=NTOOBS(lSm)
*
      If(ISM.EQ.JSM) THEN
       IJPAIRS = ITORB*(ITORB+1)/2
      ELSE
       IJPAIRS = ITORB*JTORB
      END IF
*
      IF(KSM.EQ.LSM) THEN
        KLPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KLPAIRS = KTORB*LTORB
      END IF
*
      iInt=0
      Do lJeppe=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
        Do jJeppe=jMin,jOff+jOrb-1
*
*
*. Set up array IJ*(IJ-1)/2
          IF(IJRELKL.EQ.0) THEN
            DO II = IOFF,IOFF+IORB-1
              IJ = I2INDX((JJEPPE-1)*NTOOB+II)
              IJARR(II) = IJ*(IJ-1)/2
            END DO
          END IF
*
          Do kJeppe=kOff,kOff+kOrb-1
            iMin = iOff
            kl = I2INDX(KJEPPE+(LJEPPE-1)*NTOOB)
            If(IKSM.ne.0) iMin = kJeppe
            IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &               + (KJEPPE-KOFF)*Jorb*Iorb
     &               + (JJEPPE-JOFF)*Iorb
            ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)
                IINT = (KJEPPE-KOFF)*JORB*LORB*IORB
     &               + (JJEPPE-JOFF)     *LORB*IORB
     &               + (LJEPPE-LOFF)          *IORB
            END IF
*
            IF(IJRELKL.EQ.1) THEN
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) > (ksm,lsm)
              IJKL0 = IBLOFF-1+(kl-1)*ijPairs
              IJ0 = (JJEPPE-1)*NTOOB
              Do iJeppe=iMin,iOff+iOrb-1
                  ijkl = ijkl0 + I2INDX(IJEPPE+IJ0)
                  iInt=iInt+1
                  Xint(iInt) = Intlst(ijkl)
              End Do
            END IF
*
*. block (ISM JSM !ISM JSM)
            IF(IJRELKL.EQ.0) THEN
              IJ0 = (JJEPPE-1)*NTOOB
              KLOFF = KL*(KL-1)/2
              IJKL0 = (KL-1)*IJPAIRS-KLOFF
              Do iJeppe=iMin,iOff+iOrb-1
                ij = I2INDX(IJEPPE+IJ0   )
                If ( ij.ge.kl ) Then
C                 ijkl=ij+(kl-1)*ijPairs-klOff
                  IJKL = IJKL0 + IJ
                Else
                  IJOFF = IJARR(IJEPPE)
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
                iInt=iInt+1
                Xint(iInt) = Intlst(iblOff-1+ijkl)
              End Do
            END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) < (ksm,lsm)
            IF(IJRELKL.EQ.-1) THEN
              ijkl0 = IBLOFF-1+KL - KLPAIRS
              IJ0 = (JJEPPE-1)*NTOOB
              Do iJeppe=iMin,iOff+iOrb-1
                IJKL = IJKL0 + I2INDX(IJEPPE + IJ0)*KLPAIRS
                iInt=iInt+1
                Xint(iInt) = Intlst(ijkl)
              End Do
            END IF
*
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
      IF(ISM.EQ.LSM) THEN
       ILPAIRS = ITORB*(ITORB+1)/2
      ELSE
       ILPAIRS = ITORB*LTORB
      END IF
*
      IF(KSM.EQ.JSM) THEN
        KJPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KJPAIRS = KTORB*JTORB
      END IF
*
        ilblk = max(ism,lsm)*(max(ism,lsm)-1)/2 + min(ism,lsm)
        kjblk = max(ksm,jsm)*(max(ksm,jsm)-1)/2 + min(ksm,jsm)
        IF(ILBLK.GT.KJBLK) THEN
          ILRELKJ = 1
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.EQ.KJBLK) THEN
          ILRELKJ = 0
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.LT.KJBLK) THEN
          ILRELKJ = -1
          IBLOFF = IJKLOF(MAX(KSM,JSM),MIN(KSM,JSM),MAX(ISM,LSM))
        END IF
*
        iInt=0
        Do lJeppe=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lJeppe
*
          IF(ILRELKJ.EQ.0) THEN
           DO II = IOFF,IOFF+IORB-1
             IL = I2INDX(II+(LJEPPE-1)*NTOOB)
             IJARR(II) = IL*(IL-1)/2
           END DO
          END IF
*
          Do jJeppe=jMin,jOff+jOrb-1
            Do kJeppe=kOff,kOff+kOrb-1
              KJ = I2INDX(KJEPPE+(JJEPPE-1)*NTOOB)
              KJOFF = KJ*(KJ-1)/2
              iMin = iOff
*
              IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                  IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &                  + (KJEPPE-KOFF)*Jorb*Iorb
     &                  + (JJEPPE-JOFF)*Iorb
              ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)
                IINT = (KJEPPE-KOFF)*JORB*LORB*IORB
     &               + (JJEPPE-JOFF)     *LORB*IORB
     &               + (LJEPPE-LOFF)          *IORB
              END IF
*
              If(IKSM.ne.0) iMin = kJeppe
*
              IF(ILRELKJ.EQ.1) THEN
                ILKJ0 = IBLOFF-1+( kj-1)*ilpairs
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE + IL0)
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.0) THEN
                IL0 = (LJEPPE-1)*NTOOB
                ILKJ0 = (kj-1)*ilPairs-kjOff
                Do iJeppe=iMin,iOff+iOrb-1
                  IL = I2INDX(IJEPPE + IL0 )
                  If ( il.ge.kj ) Then
C                     ilkj=il+(kj-1)*ilPairs-kjOff
                      ILKJ = IL + ILKJ0
                    Else
                      ILOFF = IJARR(IJEPPE)
                      ilkj=kj+(il-1)*kjPairs-ilOff
                    End If
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(iBLoff-1+ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.-1) THEN
                ILKJ0 = IBLOFF-1+KJ-KJPAIRS
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE+ IL0)*KJPAIRS
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
              END IF
*
            End Do
          End Do
        End Do
      End If
*
      RETURN
      END
***********************************************************************
      SUBROUTINE LGETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,ICOUL)

*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cc_exc.inc"
#include "crun.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*
      CALL QENTER('LGETINT')
      NTEST = 00
*
      IF(NTEST.GE.5)
     &WRITE(lupri,*) ' LGETINT : ICC_EXC and ICOUL = ', ICC_EXC, ICOUL

* =======================
* Usual/Normal  integrals
* =======================
*
*. Integrals in core in internal LUCIA format
      CALL GETINCN2(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &              IXCHNG,IKSM,JLSM,WORK(KINT2),
     &              WORK(KPINT2),NSMOB,WORK(KINH1),ICOUL)

      IF(NTEST.NE.0) THEN
        IF(ITP.EQ.0) THEN
          NI = NTOOBS(ISM)
        ELSE
          NI = NOBPTS(ITP,ISM)
        END IF
*
        IF(KTP.EQ.0) THEN
          NK = NTOOBS(KSM)
        ELSE
          NK = NOBPTS(KTP,KSM)
        END IF
*
        IF(IKSM.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
*
        IF(JTP.EQ.0) THEN
          NJ = NTOOBS(JSM)
        ELSE
          NJ = NOBPTS(JTP,JSM)
        END IF
*
        IF(LTP.EQ.0) THEN
          NL = NTOOBS(LSM)
        ELSE
          NL = NOBPTS(LTP,LSM)
        END IF
*
        IF(JLSM.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(lupri,*) ' 2 electron integral block for TS blocks '
        WRITE(lupri,*) ' Ixchng :', IXCHNG
        WRITE(lupri,*) ' After GETINC '
        WRITE(lupri,'(1H ,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMT_LU(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('LGETINT')
      RETURN
      END
***********************************************************************

      SUBROUTINE GETOBS_DALTON(ECORE_ENV,NAOS_ENV,NMOS_ENV,NSYM_ENV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "priunit.h"
*. Scratch
      DIMENSION NBAS(8), NOCC(8), NORB(8)
*. Output
      DIMENSION NAOS_ENV(*), NMOS_ENV(*)

C
C     Read AO and MO information on file SIRIFC written from SIRIUS.
C
      ITAP30 = 16
      OPEN(ITAP30,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND ITAP30
      CALL MOLLAB('TRCCINT ',ITAP30,lupri)
      READ (ITAP30) NSYM,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYM),
     *              (NORB(I),I=1,NSYM),(NBAS(I),I=1,NSYM),
     *              POTNUC,EMCSCF

!     transfer to output variables
      NSYM_ENV  = NSYM
      ECORE_ENV = POTNUC
      CALL ICOPY(NSYM,NORB,1,NMOS_ENV,1)
      CALL ICOPY(NSYM,NBAS,1,NAOS_ENV,1)

#ifdef LUCI_DEBUG

      write(lupri,*) 'GETOBS_DALTON' !hjaaj DEBUG
      WRITE(lupri,*) ' Number of basis functions per sym '
      CALL IWRTMA(NBAS,1,8,1,8)
*
      WRITE(lupri,*) ' Norb as delivered from environment '
      CALL IWRTMA(NORB,1,8,1,8)
*
      WRITE(lupri,*) ' NOCC as delivered from DALTON (discarded)'
      CALL IWRTMA(NOCC,1,8,1,8)
*.
      WRITE(lupri,*) ' from DALTON: NORBT, NBAST, NCMOT :',
     &   NORBT,NBAST,NCMOT

#endif

      END
***********************************************************************
       SUBROUTINE GETOBS_LUCIA(NAOS_ENV,NMOS_ENV)
*
* Obtain info on orbital dimensions from LU91 - LUCIA format
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      INTEGER NMOS_ENV(*),NAOS_ENV(*)
*
      LUH = 91
      CALL REWINO(LUH)
*.
      READ(LUH,*) NSMOB
*.
      READ(LUH,*) (NMOS_ENV(ISM),ISM=1, NSMOB)
*
      READ(LUH,*) (NAOS_ENV(ISM),ISM=1, NSMOB)
*
      RETURN
      END
***********************************************************************
      FUNCTION GTIJKL(I,J,K,L)
*
* Obtain  integral (I J ! K L )
* where I,J,K and l refers to active orbitals in
* Type ordering
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*.GLobal pointers
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
      IF(INTIMP .EQ. 2 ) THEN
*. LUCAS ordering
        I12S = 0
        I34S = 0
        I1234S = 1
        GTIJKL = GIJKLL(IREOTS(1),WORK(KPINT2),WORK(KLSM2),
     &                  WORK(KINT2),ISMFTO,IBSO,NACOB,NSMOB,
     &                  NOCOBS,I,J,K,L)
      ELSE IF (INTIMP.EQ.1.OR.INTIMP.EQ.5.or.INTIMP.eq.6) THEN
!       sirius or dirac integral format
        GTIJKL = GMIJKL(I,J,K,L,WORK(KINT2),WORK(KPINT2))
      END IF

      END

***********************************************************************
      FUNCTION GMIJKL(IORB,JORB,KORB,LORB,INTLST,IJKLOF)
*
* Obtain integral (IORB JORB ! KORB LORB) MOLCAS version
* Integrals assumed in core
*
* Version for integrals stored in INTLST
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
*. The orbital list corresponds to type ordered indices, reform to
*. symmetry ordering
*
      IABS = IREOTS(IORB)
      ISM  = ISMFTO(IORB)
      IOFF = IBSO(ISM)
*
      JABS = IREOTS(JORB)
      JSM  = ISMFTO(JORB)
      JOFF = IBSO(JSM)
*
      KABS = IREOTS(KORB)
      KSM  = ISMFTO(KORB)
      KOFF = IBSO(KSM)
*
      LABS = IREOTS(LORB)
      LSM  = ISMFTO(LORB)
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(lupri,*) ' GMIJKL at your service '
        WRITE(lupri,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(lupri,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(lupri,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(lupri,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
        call flshfo(lupri)
      END IF
*
      If ( jSm.gt.iSm .or. ( iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JABS - JOFF + 1
        J = IABS - IOFF + 1
      Else
        iSym=iSm
        jSym=jSm
        I = IABS - IOFF + 1
        J = JABS - JOFF + 1
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( lSm.gt.kSm  .or. ( kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LABS -LOFF + 1
        L = KABS - KOFF + 1
      Else
        kSym=kSm
        lSym=lSm
        K = KABS - KOFF + 1
        L = LABS -LOFF + 1
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(lupri,*) ' i j k l ',i,j,k,l
        write(lupri,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
        call flshfo(lupri)
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(lupri,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )then
        write(lupri,*) ' ijkl ', ijkl
        write(lupri,*) ' iblOff       ', iblOff
        write(lupri,*) ' iblOff-1+ijkl', iblOff-1+ijkl
        call flshfo(lupri)
      end if
*
      GMIJKL = Intlst(iblOff-1+ijkl)
      If( ntest .ge. 100 )
     & write(lupri,*) ' GMIJKL ', GMIJKL
*
      END
***********************************************************************
      SUBROUTINE GTJK(RJ,RK,NTOOB)
*
* Interface routine for obtaining Coulomb (RJ) and
* Exchange integrals (RK)
*
* Ordering of intgrals is in the internal order
      IMPLICIT REAL*8(A-H,O-Z)
*
*.CRUN
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,NP1,NP2,NQ,INCORE,MXCIV,ICISTR,
C    &            NOCSF,IDIAG
#include "mxpdim.inc"
#include "parluci.h"
#include "crun.inc"
*.Output
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)

      IF(INTIMP.EQ.1.OR.INTIMP.EQ.5.or.INTIMP.eq.6) THEN
*. Interface to SIRIUS
        CALL GTJKS(RJ,RK,NTOOB)
      ELSE
*. Interface to LUCAS integrals
        CALL GTJKL(RJ,RK,NTOOB)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(luwrt,*) ' RJ and RK from GTJK '
        CALL WRTMT_LU(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
        CALL WRTMT_LU(RK,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GTJKL(RJ,RK,NTOOB)
*
* Obtain Coulomb  integrals (II!JJ)
*        exchange integrals (IJ!JI)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*
      DO 100 IORB = 1, NTOOB
        DO 50 JORB = 1, NTOOB
          RJ(IORB,JORB) = GTIJKL(IORB,IORB,JORB,JORB)
          RK(IORB,JORB) = GTIJKL(IORB,JORB,JORB,IORB)
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(lupri,*) ' RJ and RK from GTJK '
        CALL WRTMT_LU(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
        CALL WRTMT_LU(RK,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
***********************************************************************

* Working on EXPHAM
* some known problems :
*     1 : if CSF are used diagonal is not delivered to H0mat
      SUBROUTINE GTJKS(J,K,NORB)
*
* Obtain Coulomb and Exchange integrals
* from complete integral list stored in core
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "priunit.h"
      REAL*8           J(NORB,NORB),K(NORB,NORB)
*
      DO 200 IORB = 1, NORB
        DO 100 JORB = 1, NORB
!         write(lupri,*) 'iorb,jorb ==> j',iorb,jorb
          J(IORB,JORB) = GTIJKL(IORB,IORB,JORB,JORB)
!         write(lupri,*) 'iorb,jorb ==> k',iorb,jorb
          K(IORB,JORB) = GTIJKL(IORB,JORB,JORB,IORB)
  100   CONTINUE
  200 CONTINUE
*
      END
***********************************************************************
      subroutine hello_dalton_lucita
************************************************************************
*                                                                      *
*     Print the program banner, date and time of execution             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*     Modified, Timo Fleig, Dec 2001, for DIRAC                        *
*                           Aug 2004                                   *
*                           Aug 2006                                   *
*               HJAaJ, May 2008, for DALTON                            *
*                                                                      *
************************************************************************
#include "priunit.h"
      Character*8   Fmt
      Character*70  Line,StLine
      integer    :: lpaper = 120 
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        StLine(i:i)='*'
      End Do
      left=(lPaper-lLine)/2
      Write(Fmt,'(A,I3.3,A)') '(',left,'X,A)'
*----------------------------------------------------------------------*
*     Print the program header                                         *
*----------------------------------------------------------------------*
      write(lupri,'(/1x,a)') StLine
      write(lupri,'( 1x,a)') StLine,' ',
     &   'D A L T O N - L U C I T A',
     &   'An interface section for LUCITA under DALTON', ' ',
     &   'Authors: J. Olsen, Univ. Aarhus',
     &   '         H. J. Aa. Jensen, Univ. Southern Denmark ',
     &   '         S. Knecht, Univ. Southern Denmark',' ',
!    &   'Based on LUCITA-DIRAC interface',
!    &   '    Author: Timo Fleig, Univ. Toulouse', ' ',
     &   'Using LUCIA version 1999', ' ',
     &   '    Author: J. Olsen, Lund/Aarhus',' ',
     &   'Parallelization of LUCITA, Duesseldorf/Odense:    ',
     &   '  S. Knecht, Univ. Southern Denmark', ' ',
     &   'Citations:',
     &   '  J. Olsen, P. Joergensen, J. Simons,             ',
     &   '          Chem. Phys. Lett. 169 (1990) 463        ',
     &   '  S. Knecht, H. J. Aa. Jensen and T. Fleig,       ',
     &   '          J. Chem. Phys., 128 (2008) 014108       ',' ',
     &   StLine, Stline, Stline

      end
***********************************************************************
      Function I2EAD(IORB,JORB,KORB,LORB)
*
* Find adress of integral in LUCIA order
*
      IMPLICIT REAL*8           (A-H,O-Z)
*
#include "mxpdim.inc"

#include "glbbas.inc"
*
#include "wrkspc.inc"
*
      I2EAD = I2EADS(IORB,JORB,KORB,LORB,WORK(KPINT2))
*
      RETURN
      END
***********************************************************************
      FUNCTION I2EADS(IORB,JORB,KORB,LORB,IJKLOF)
*
* Obtain address of integral (IORB JORB ! KORB LORB) in MOLCAS order
* IORB JORB KORB LORB corresponds to SYMMETRY ordered indices !!
* Integrals assumed in core
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "priunit.h"
*
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
      IABS = IORB
      ISM = ISMFTO(IREOST(IORB))
      IOFF = IBSO(ISM)
*
      JABS = JORB
      JSM = ISMFTO(IREOST(JORB))
      JOFF = IBSO(JSM)
*
      KABS = KORB
      KSM = ISMFTO(IREOST(KORB))
      KOFF = IBSO(KSM)
*
      LABS = LORB
      LSM = ISMFTO(IREOST(LORB))
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(lupri,*) ' I2EADS at your service '
        WRITE(lupri,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(lupri,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(lupri,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(lupri,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
      END IF
*
      If ( jSm.gt.iSm .or. ( iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JABS - JOFF + 1
        J = IABS - IOFF + 1
      Else
        iSym=iSm
        jSym=jSm
        I = IABS - IOFF + 1
        J = JABS - JOFF + 1
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( lSm.gt.kSm  .or. ( kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LABS -LOFF + 1
        L = KABS - KOFF + 1
      Else
        kSym=kSm
        lSym=lSm
        K = KABS - KOFF + 1
        L = LABS -LOFF + 1
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(lupri,*) ' i j k l ',i,j,k,l
        write(lupri,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(lupri,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )
     & write(lupri,*) ' ijkl ', ijkl
*
      I2EADS = iblOff-1+ijkl
      If( ntest .ge. 100 ) then
        write(lupri,*) 'i j k l ', i,j,k,l
        write(lupri,*) ' ibloff ijkl ',ibloff,ijkl
        write(lupri,*) ' I2EADS  = ', I2EADS
      END IF
*
      END

***********************************************************************
      subroutine intim(citask_id,int1_mcscf,int2_mcscf,
     &                 update_ijkl_from_env,
     &                 orbital_trial_vector,
     &                 ci_trial_vector,
     &                 mcscf_provides_integrals)
*
* Interface to external integrals
*
* If NOINT .ne. 0, only pointers are constructed
* Jeppe Olsen, Winter of 1991
*
* Version : Fall 97
*
*
#ifdef VAR_MPI
      use sync_coworkers
#endif
      use lucita_cfg
      use lucita_energy_types
      use lucita_mcscf_ci_cfg
#include "implicit.h"
      ! for addressing of WORK
      INTEGER*8 current_free_mem, kintdal_interface
      character (len=12), intent(in) :: citask_id
      real(8), intent(inout)         :: int1_mcscf(*)
      real(8), intent(inout)         :: int2_mcscf(*)
      logical, intent(inout)         :: update_ijkl_from_env
      logical, intent(in)            :: orbital_trial_vector
      logical, intent(in)            :: ci_trial_vector
      logical, intent(in)            :: mcscf_provides_integrals
#include "priunit.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "lucinp.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "parluci.h"
#include "oper.inc"
*./CINTFO/
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      character (len=25)   :: file_name
      logical              :: update_ijkl, integral_available
      integer              :: ngsh_lucita_tmp(max_number_of_gas_spaces,
     &                                        max_number_of_ptg_irreps)
      integral_available = .false.

!     pointers for symmetry blocks of integrals
      CALL INTPNT(WORK(KPINT1),WORK(KLSM1),
     &            WORK(KPINT2),WORK(KLSM2))

!     pointer for orbital indices for symmetry blocked matrices
      CALL ORBINH1(WORK(KINH1),NTOOBS,NTOOB,NSMOB)

!
!     load one-electron integrals and two-electron integrals
      if(citask_id.eq.'return CIdia' .or. 
     &   citask_id.eq.'perform Davi' .or. 
     &   citask_id.eq.'ijkl resort ' .or. 
     &   citask_id.eq.'return sigma' )then

        IF(INCORE.EQ.1)THEN
          write(file_name,'(a17,i3,a1,i4)')
     &    'IJKL_REOD_LUCITA-',1,'.',luci_myproc
          do i = 18, 20
            if(file_name(i:i).eq. ' ') file_name(i:i) = '0'
          end do
          do i = 22, 25
            if(file_name(i:i).eq. ' ') file_name(i:i) = '0'
          end do
          LUINT_LUCITA = -1
!
!         check for existing file with reordered integrals -
!         if available read from it...
          if(citask_id.eq.'perform Davi' .or. 
     &       citask_id.eq.'ijkl resort ')then
            if(luci_myproc .eq. luci_master)then
              inquire(file=file_name,exist=integral_available)
            end if
          end if
!
          CALL GPOPEN(LUINT_LUCITA,file_name,'UNKNOWN',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
          is_2e_int_active = i12

          update_ijkl =
     &    update_ijkl_from_env.or.ci_trial_vector.or.
     &    orbital_trial_vector

!         default for reordering task
          if(citask_id.eq.'ijkl resort ') update_ijkl = .true.

          if(update_ijkl)then

            if(luci_myproc .eq. luci_master)then
              if(.not.integral_available)then
                CALL MEMMAN(current_free_mem,0,'SFREEM',2,'SHOWMM')
                LSCR = current_free_mem - 1000
                idum = 0
!               set local marker + allocate space for lucita-dalton integral interface
                call memman(idum,idum,'MARK  ',idum,'dalint')
                call memman(kintdal_interface,lscr,'ADDL  ',2,'dalint')

                CALL LUCITA_GETINT_DALTON(einact_mc2lu,WORK(KINT1),
     &                                    WORK(KINT2),
     &                                    WORK(kintdal_interface),LSCR,
     &                                    int1_mcscf,int2_mcscf,
     &                                    mcscf_provides_integrals,1)
                ecore = ecore_orig + einact_mc2lu
!               release memory
                CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'dalint')
              else
!               read 1e- and 2e-integrals from file
!               -----------------------------------
                CALL MOLLAB('LUCIIJKL',LUINT_LUCITA,LUPRI)
                READ (LUINT_LUCITA) nint1, nint2, ecore, ecore_orig ,
     &                              einact_mc2lu
                read(LUINT_LUCITA) 
     &          ((ngsh_lucita_tmp(i,j),i=1,max_number_of_gas_spaces),
     &                                 j=1,max_number_of_ptg_irreps)
                CALL READT(LUINT_LUCITA, nint1, work(kint1))
                if(is_2e_int_active .gt. 1)
     &          CALL READT(LUINT_LUCITA, nint2, work(kint2))
!               check consistency in active spaces for integral run 
!               and actual run
                do i = 1, max_number_of_gas_spaces
                  do j = 1, max_number_of_ptg_irreps
                    if(ngsh_lucita_tmp(i,j) /= ngsh_lucita(i,j))then
                    write(lupri,*) '*** error in intim: active spaces'//
     &              ' in integral run and present CI run do not match.'
                    write(lupri,*) 'integral run ',ngsh_lucita_tmp(:,:)
                    write(lupri,*) 'actual CI    ',ngsh_lucita(:,:)
                    call quit('inconsistency in active spaces for the'//
     &              ' integral CI-run and the present CI-run.')
                    end if
                  end do
                end do
              end if
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
              write(lupri,*) 
     &        '*** master   ',luci_myproc,'reports 1-ints:',
     &                        nint1
              call wrtmatmn(work(kint1),1,nint1,1,nint1,lupri)
              write(lupri,*) 
     &        '*** master   ',luci_myproc,'reports 2-ints:',
     &                        nint2
              call wrtmatmn(work(kint2),1,nint2,1,nint2,lupri)
              write(lupri,*) 
     &        '*** master   ',luci_myproc,'reports ecore  ',
     &                        ecore, ecore_orig, einact_mc2lu
#undef LUCI_DEBUG
#endif
            end if ! luci_master only
#ifdef VAR_MPI
            if(luci_nmproc .gt. 1)then
!             set sync_ctrl option
              sync_ctrl_array(3) = update_ijkl
!             synchronize the co-workers with the 1-/2-electron integrals
              call sync_coworkers_cfg(work(kint1),work(kint2))

#ifdef LUCI_DEBUG
              if(luci_myproc .eq. luci_master+1)then
                print *, '*** co-worker',luci_myproc,'reports 1-ints:',
     &                        nint1
                call wrtmatmn(work(kint1),1,nint1,1,nint1,lupri)
                print *, '*** co-worker',luci_myproc,'reports 2-ints:',
     &                        nint2
                call wrtmatmn(work(kint2),1,nint2,1,nint2,lupri)
              end if ! debug print for co-worker-id luci_master+1
!#undef LUCI_DEBUG
#endif
#ifdef STEFAN_PLEASE_FIX
<<<<<<< HEAD
=======
        ELSE
          CALL QUIT('LUCITA out-of-core not implemented yet, sorry!')
        END IF
*
* ===============
* . DIRAC format
* ===============
*
      else if (INTIMP.eq.6.and.NOINT.eq.0) then
        call quit('*** error in intim: dirac integral interface
     & missing ***')
!       call rdint_dirac(KINT1,KINT2,DOUGRP,IPRNT)
      END IF
*
*. Change one-electron integrals to inactive fock matrix
      IF(NOINT.EQ.0) THEN

        call dcopy(nint1,work(kint1),1,work(kint1o),1)

        IF(IUSE_PH.EQ.1) THEN
           CALL FI(WORK(KINT1),ECORE_HEX,1)
        ELSE
           ECORE_HEX = 0.0D0
        END IF
      END IF

      ECORE_ORIG = ECORE
      ECORE      = ECORE + ECORE_HEX
#ifdef LUCI_DEBUG
      if(luci_myproc.eq.luci_master)then
        write(lupri,'(/2x,e15.8)') 'Updated core energy: ',ECORE
      end if
#endif
*
      END
***********************************************************************
*
*  Create usual input file for LUCITA by writing keywords and
*  values to unit ..
*
      subroutine luciwrt(IPLOCAL,NSHPGS,MXPNGAS,MXPIRR,MXPICI,
     &                   NGAS,NSYM,NCISPC,NSEQ,TYPELUC,MAXIT,
     &                   MXSEQ,IRESDIM,LCSBLK,NACTEL,MSCOMB,SCOMB,
     &                   IRREP,NROOT,MULTIP,MS2,NELGS,
     &                   CALSIZ,IVFLEV,
     &                   IUSE_PH,IDEMOL,IRSTRT,MXCIVL,
     &                   IPRSTR,IPRCIX,IPRORB,IPRDIA,
     &                   IPRXT,IPRRSP,IPRDEN,IPROCC,IPRNCIV)
*
      implicit real*8 (A-H,O-Z)
*
#include "priunit.h"
#include "units.inc"
#include "parluci.h"
*
      logical      MSCOMB
      character*2  TYPELUC(MXSEQ)
      character*3  CALSIZ
      character*72 HESTR
      CHARACTER    LUCIAIN_MPIN*7, NLUCIAIN_MPI*11
      INTEGER      LUCIAIN_MPIL
      dimension    NSHPGS(MXPNGAS,MXPIRR),NELGS(MXPNGAS,2,MXPICI)
*
*-----------------------------------------------------------------*
*  Print input for LUCITA input file           
*-----------------------------------------------------------------*
*
!     IPLOCAL = 5
      if (IPLOCAL.ge.5) then
         write(lupri,*) '==========================================='
         write(lupri,*) '  Input written to file LUCIAIN     '
         write(lupri,*) '==========================================='
*
         write(lupri,'(A13)') '* &LUCIA &END'
*
         if (CALSIZ.eq.'NOR') then
           write(lupri,'(A6)') 'MEGACI'
         else
           write(lupri,'(A6)') 'TERACI'
         end if
*
         write(lupri,'(A6)') 'DALTON'
*
         write(lupri,'(A3)') 'D2H'
*
         write(lupri,'(A6)') 'NIRREP'
         write(lupri,'(1X,I3)') NSYM
*
         write(lupri,'(A6)') 'IREFSM'
         write(lupri,'(1X,I3)') IRREP
*
         write(lupri,'(A5)') 'GASSH'
         write(lupri,'(1X,I3)') NGAS
*
         HESTR = ' '
         write (HESTR,'(A,I1,A)') '(1X,',NSYM-1,'(I3,A1),I3)'
         if (NSYM.eq.1) write (HESTR,'(A)') '(I4)'
         do IGAS = 1,NGAS,1
            write(lupri,HESTR)
     &      (NSHPGS(IGAS,ISYM),',',ISYM=1,NSYM-1),NSHPGS(IGAS,NSYM)
         end do
*
         write(lupri,'(A6)') 'NACTEL'
         write(lupri,'(1X,I3)') NACTEL
*
         write(lupri,'(A6)') 'GASSPC'
         write(lupri,'(1X,I3)') NCISPC
         do ISPC=1,NCISPC,1
            write(lupri,'(1X,2I3)')
     &           (NELGS(IGAS,1,ISPC),NELGS(IGAS,2,ISPC),IGAS=1,NGAS)
         end do
*
         write(lupri,'(A6)') 'SEQUEN'
         do ICISPC=1,NCISPC,1
            write(lupri,'(1X,I3)') NSEQ
            do ISEQ=1,NSEQ,1
               if (TYPELUC(ISEQ).eq.'CI') then
                  write(lupri,'(1X,A3,I3)') 'CI,',MAXIT
               else if (TYPELUC(ISEQ).eq.'PT') then
                  write(lupri,'(1X,A6,I3)') 'PERTU,',MAXIT
               else if (TYPELUC(ISEQ).eq.'VF') then
                  write(lupri,'(1X,A8,I3)') 'VECFREE,',IVFLEV
               end if
            end do
         end do
*
         write(lupri,'(A6)') 'RESTRT'
         write(lupri,'(1X,I3)') IRSTRT
*
         write(lupri,'(A5)') 'MULTS'
         write(lupri,'(1X,I3)') MULTIP
*
         write(lupri,'(A3)') 'MS2'
         write(lupri,'(1X,I3)') MS2
*
         if (MSCOMB) then
            write(lupri,'(A6)') 'MSCOMB'
            write(lupri,'(F4.1)') SCOMB
         end if
*
         write(lupri,'(A5)') 'ROOTS'
         write(lupri,'(1X,I3)') NROOT

         write(lupri,'(A5)') 'MXCIV'
         write(lupri,'(1X,I3)') MXCIVL
*
         write(lupri,'(A5)') 'DENSI'
         write(lupri,'(1X,I3)') IDEMOL
*
         write(lupri,'(A5)') 'NOCSF'
*
         write(lupri,'(A6)') 'IPRSTR'
         write(lupri,'(1X,I3)') IPRSTR
         write(lupri,'(A6)') 'IPRCIX'
         write(lupri,'(1X,I3)') IPRCIX
         write(lupri,'(A6)') 'IPRORB'
         write(lupri,'(1X,I3)') IPRORB
         write(lupri,'(A6)') 'IPRDIA'
         write(lupri,'(1X,I3)') IPRDIA
         write(lupri,'(A5)') 'IPRXT'
         write(lupri,'(1X,I3)') IPRXT
         write(lupri,'(A6)') 'IPRRSP'
         write(lupri,'(1X,I3)') IPRRSP
         write(lupri,'(A6)') 'IPRDEN'
         write(lupri,'(1X,I3)') IPRDEN
         write(lupri,'(A6)') 'IPROCC'
         write(lupri,'(1X,I3)') IPROCC
         write(lupri,'(A6)') 'PRNCIV'
         write(lupri,'(1X,I3,A)') IPRNCIV
*
         write(lupri,'(A6)') 'RESDIM'
         write(lupri,'(1X,I8)') IRESDIM
*
         write(lupri,'(A6)') 'LCSBLK'
         write(lupri,'(I20)') LCSBLK
*
         if (IUSE_PH.eq.1) then
            write(lupri,'(A6)') 'USE_PH'
            write(lupri,'(A6)') 'ADVICE'
         end if
*
         if (CALSIZ.ne.'NOR') then
            write(lupri,'(A6)') 'CLSSEL'
         end if
*
         write(lupri,'(A12)') 'END OF INPUT'
*
      end if
*
*--------------------------------------------------
* Write to file LUCIAIN
*--------------------------------------------------
      LUCIAIN_MPIN = "LUCIAIN"
      IF (LUCI_MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUCIAIN_MPI,'(A7,A1,I1)') LUCIAIN_MPIN,'.',LUCI_MYPROC
         LUCIAIN_MPIL=9
      ELSEIF (LUCI_MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUCIAIN_MPI,'(A7,A1,I2)') LUCIAIN_MPIN,'.',LUCI_MYPROC
         LUCIAIN_MPIL=10
      ELSEIF (LUCI_MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUCIAIN_MPI,'(A7,A1,I3)') LUCIAIN_MPIN,'.',LUCI_MYPROC
         LUCIAIN_MPIL=11
      ELSE
         CALL QUIT("LUCI_NMPROC.GT.1000! EXTEND LUCINOD.F MODULE")
      ENDIF

      open(unit=LUINP,file=NLUCIAIN_MPI(1:LUCIAIN_MPIL),
     &     status='UNKNOWN',form='FORMATTED')
         Rewind(LUINP)
         write(LUINP,'(A13)') '* &LUCIA &END'
*
         if (CALSIZ.eq.'NOR') then
           write(LUINP,'(A6)') 'MEGACI'
         else
           write(LUINP,'(A6)') 'TERACI'
         end if
*
         write(LUINP,'(A6)') 'DALTON'
*
         write(LUINP,'(A3)') 'D2H'
*
         write(LUINP,'(A6)') 'NIRREP'
         write(LUINP,'(1X,I3)') NSYM
*
         write(LUINP,'(A6)') 'IREFSM'
         write(LUINP,'(1X,I3)') IRREP
*
         write(LUINP,'(A5)') 'GASSH'
         write(LUINP,'(1X,I3)') NGAS
*
         HESTR = ' '
         write (HESTR,'(A,I1,A)') '(1X,',NSYM-1,'(I3,A1),I3)'
         if (NSYM.eq.1) write (HESTR,'(A)') '(1X,I3)'
         do IGAS = 1,NGAS,1
            write(LUINP,HESTR)
     &      (NSHPGS(IGAS,ISYM),',',ISYM=1,NSYM-1),NSHPGS(IGAS,NSYM)
         end do
*
         write(LUINP,'(A6)') 'NACTEL'
         write(LUINP,'(I4)') NACTEL
*
         write(LUINP,'(A6)') 'GASSPC'
         write(LUINP,'(I4)') NCISPC
         do ISPC=1,NCISPC,1
            write(LUINP,'(1X,2I3)')
     &           (NELGS(IGAS,1,ISPC),NELGS(IGAS,2,ISPC),IGAS=1,NGAS)
         end do
*
         write(LUINP,'(A6)') 'SEQUEN'
         do ICISPC=1,NCISPC,1
            write(LUINP,'(1X,I3)') NSEQ
            do ISEQ=1,NSEQ,1
               if (TYPELUC(ISEQ).eq.'CI') then
                  write(LUINP,'(1X,A3,I3)') 'CI,',MAXIT
               else if (TYPELUC(ISEQ).eq.'PT') then
                  write(LUINP,'(1X,A6,I3)') 'PERTU,',MAXIT
               else if (TYPELUC(ISEQ).eq.'VF') then
                  write(LUINP,'(1X,A8,I3)') 'VECFREE,',IVFLEV
               end if
            end do
         end do
*
         if (IRSTRT.eq.1) then
           write(LUINP,'(A6)') 'RESTRT'
         end if
*
         write(LUINP,'(A5)') 'MULTS'
         write(LUINP,'(1X,I3)') MULTIP
*
         write(LUINP,'(A3)') 'MS2'
         write(LUINP,'(1X,I3)') MS2
*
         if (MSCOMB) then
            write(LUINP,'(A6)') 'MSCOMB'
            write(LUINP,'(F4.1)') SCOMB
         end if
*
         write(LUINP,'(A5)') 'ROOTS'
         write(LUINP,'(1X,I3)') NROOT

         write(LUINP,'(A5)') 'MXCIV'
         write(LUINP,'(1X,I3)') MXCIVL
*
         write(LUINP,'(A5)') 'DENSI'
         write(LUINP,'(1X,I3)') IDEMOL
*
         write(LUINP,'(A5)') 'NOCSF'
*
         write(LUINP,'(A6)') 'IPRSTR'
         write(LUINP,'(1X,I3)') IPRSTR
         write(LUINP,'(A6)') 'IPRCIX'
         write(LUINP,'(1X,I3)') IPRCIX
         write(LUINP,'(A6)') 'IPRORB'
         write(LUINP,'(1X,I3)') IPRORB
         write(LUINP,'(A6)') 'IPRDIA'
         write(LUINP,'(1X,I3)') IPRDIA
         write(LUINP,'(A5)') 'IPRXT'
         write(LUINP,'(1X,I3)') IPRXT
         write(LUINP,'(A6)') 'IPRRSP'
         write(LUINP,'(1X,I3)') IPRRSP
         write(LUINP,'(A6)') 'IPRDEN'
         write(LUINP,'(1X,I3)') IPRDEN
         write(LUINP,'(A6)') 'IPROCC'
         write(LUINP,'(1X,I3)') IPROCC
         write(LUINP,'(A6)') 'PRNCIV'
         write(LUINP,'(1X,I3)') IPRNCIV
*
         write(LUINP,'(A6)') 'RESDIM'
         write(LUINP,'(1X,I8)') IRESDIM
*
         write(LUINP,'(A6)') 'LCSBLK'
         write(LUINP,'(I20)') LCSBLK
*
         if (IUSE_PH.eq.1) then
            write(LUINP,'(A6)') 'USE_PH'
            write(LUINP,'(A6)') 'ADVICE'
         end if
*
         if (CALSIZ.ne.'NOR') then
            write(LUINP,'(A6)') 'CLSSEL'
         end if
*
         write(LUINP,'(A12)') 'END OF INPUT'
*
      close (unit=LUINP)
*
      end
***********************************************************************
*  Conversion routine from 'dummy' DALTON input for LUCITA to
*  true LUCITA input. This routine creates the required
*  default settings and processes the input already given by
*  DALTON.
*************************************************************
*
      subroutine dalton2lucita(IRREP,NROOT,MULTIP,MXHOL1,MXELR3,INACOB,
     &                         IDEMOL,IRSTRT,MXCIVL,IMAXCIITER,
     &                         TYPE,WAFFCT,CALSIZ,PRILUC,NINOBS,
     &                         NGAS,NCISPC,IPRT)
*
      use lucita_cfg
      implicit real*8 (A-H,O-Z)
#include "priunit.h"
#include "mxpdim.inc"
#include "lucinp.inc"
#include "SysDef.inc"
*
#include "inforb.h"
#include "parluci.h"
*
      logical MSCOMB
*
      character*2 TYPELUC(MXSEQ)
      character*3 CALSIZ
      character*6 TYPE, WAFFCT, PRILUC
*
      dimension NSHPGS(MXPNGAS,MXPIRR),NELGS(MXPNGAS,2,MXPICI)
      dimension IHELP(4),NINOBS(MXPIRR)
*
*
!     NSYM = NSMOB ! transfer number of symmetries from lucinp.inc
!     print *, '*** note from stefan: is it not dangerous 
!    &to re-set NSYM here in dalton2lucita??? NSYM is global in sirius!'
!     print *,'NSMOB is',NSMOB
!
!     set IRESDIM independent from ci type
      IRESDIM = 100 000
*
*------------------------------------------------------------------*
*  Create default settings for LUCITA
*------------------------------------------------------------------*
*
*------------------------------------------------------------------*
*  Import required information from Dalton
*------------------------------------------------------------------*
      call set_nracte_wf(type,WAFFCT,NACRAS)
!     print *, ' stefan: define NACRAS for proper FCI/SDCI!!!***'
*
*.................................................................
*  1) Inactive space as GAS 1 (if stated in input)
*.................................................................
*
*
*     inactive shells and total number of inactive orbitals
      nsum    =  0
      iretish = -1
      if(inacob .eq. 1)then
        call icopy(nsym,nish_lucita,1,ninobs,1)
        nsum = 0
        do isym = 1, nsym
          nsum = nsum + ninobs(isym)
        end do
      end if
      if(nsum.gt.0) iretish = 0
*
*.................................................................
*  2) TYPE dependencies (keywords depending on type of CI chosen)
*.................................................................
*         ................
      if (TYPE.eq.'FCI   ') then
*         ................
*
*  NGSSH, NGAS
        if (INACOB.eq.0) then
          NGAS = 1
          do ISYM=1,NSYM,1
            NSHPGS(1,ISYM) = NORB(ISYM)
          end do
        else if (INACOB.eq.1) then
          NGAS = 2
          do ISYM=1,NSYM,1
            NSHPGS(1,ISYM) = NINOBS(ISYM)
            NSHPGS(2,ISYM) = NORB(ISYM) - NINOBS(ISYM)
          end do
        end if
*
*  NACTEL
        NTOCCH = 0
        NTOCCR = 0
        if (NACTEL.eq.-1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NTOCCH = NTOCCH + 2*NOCC(ISYM)
            end do
            NACTEL = NTOCCH
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NTOCCR = NTOCCR + 2*NISH(ISYM)
            end do
            NOCCIM = NTOCCR + NACRAS
            NACTEL = NOCCIM
          end if
          write(lupri,'(A42,I3/)')
     &      ' Number of active electrons ........... ',NACTEL
        end if
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
         write(lupri,*) 'Number of active electrons exceeds '
         write(lupri,*) 'the orbital space.'
         write(lupri,*) 'Consider Paulis famous principle and restart!'
         Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        NCISPC = 1
        ICISPC = 1
        if (INACOB.eq.0) then
          NELGS(1,1,ICISPC) = NACTEL
          NELGS(1,2,ICISPC) = NACTEL
        else if (INACOB.eq.1) then
          NELGS(1,1,ICISPC) = 2*NSUM
          NELGS(1,2,ICISPC) = 2*NSUM
          NELGS(2,1,ICISPC) = NACTEL
          NELGS(2,2,ICISPC) = NACTEL
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(lupri,*) 'Illegal spin multiplicity given.'
           write(lupri,*) 'Number of active electrons: ',NACTEL
           write(lupri,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(lupri,*) 'Illegal spin multiplicity given.'
           write(lupri,*) 'Number of active electrons: ',NACTEL
           write(lupri,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(lupri,*) 'Illegal spin multiplicity given.'
              write(lupri,*) 'Compare with number of active'
              write(lupri,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(lupri,*) 'Illegal spin multiplicity given.'
              write(lupri,*) 'Compare with number of active'
              write(lupri,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        MAXIT = IMAXCIITER
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
*              ..............
      else if (TYPE.eq.'SDCI  ') then
*              ..............
*
*  NGSSH, NGAS
        if (INACOB.eq.0) then
          if (WAFFCT.eq.'HF_SCF') then
            NGAS = 2
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NOCC(ISYM)
              NSHPGS(2,ISYM) = NSSH(ISYM)
              if (NSHPGS(2,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),2,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=3,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do
              end if
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              NGAS = 2
              do ISYM=1,NSYM,1
                NSHPGS(1,ISYM) = NASH(ISYM)
                NSHPGS(2,ISYM) = NSSH(ISYM)
                if (NSHPGS(2,ISYM).gt.MXTSOB) then
                  call gassplit(NSHPGS(1,ISYM),2,MXPNGAS,
     &                          MXTSOB,NGAS,IPRT)
                else
                  do IGAS=3,NGAS,1
                    NSHPGS(IGAS,ISYM) = 0
                  end do
                end if
              end do
            else if (IRETISH.eq.0) then
              NGAS = 3
              do ISYM=1,NSYM,1
                NSHPGS(1,ISYM) = NISH(ISYM)
                NSHPGS(2,ISYM) = NASH(ISYM)
                NSHPGS(3,ISYM) = NSSH(ISYM)
                if (NSHPGS(3,ISYM).gt.MXTSOB) then
                  call gassplit(NSHPGS(1,ISYM),3,MXPNGAS,
     &                          MXTSOB,NGAS,IPRT)
                else
                  do IGAS=4,NGAS,1
                    NSHPGS(IGAS,ISYM) = 0
                  end do
                end if
              end do
            else
              write(lupri,*) 'IRETISH has unallowed value.'
              write(lupri,*) 'IRETISH = ',IRETISH
              Call Abend2('Quitting.')
>>>>>>> origin/master
#endif  /* STEFAN_PLEASE_FIX */
            end if
#endif   /* VAR_MPI */
            
            if(.not.orbital_trial_vector.and..not.ci_trial_vector)then
!             put 1e- and 2e-integrals to file
!             --------------------------------
              WRITE (LUINT_LUCITA) ('********',I=1,3), 'LUCIIJKL'
              WRITE(LUINT_LUCITA) nint1, nint2, ecore, ecore_orig ,
     &                            einact_mc2lu
              write(LUINT_LUCITA) 
     &        ((ngsh_lucita(i,j),i=1,max_number_of_gas_spaces),
     &                           j=1,max_number_of_ptg_irreps)
              CALL WRITT(LUINT_LUCITA, nint1, work(kint1))
              if(is_2e_int_active .gt. 1)
     &        CALL WRITT(LUINT_LUCITA, nint2, work(kint2))
              WRITE (LUINT_LUCITA) ('********',I=1,3), 'EODATA  '
              update_ijkl_from_env = .false.
            end if
!           1e-integrals from MCSCF environment, 2e-integrals from file
            if(ci_trial_vector)then
              if(is_2e_int_active .gt. 1)then
                CALL MOLLAB('LUCIIJKL',LUINT_LUCITA,LUPRI)
                READ (LUINT_LUCITA) nint1, nint2, ecore, ecore_orig ,
     &                              einact_mc2lu
                read(LUINT_LUCITA) 
     &          ((ngsh_lucita_tmp(i,j),i=1,max_number_of_gas_spaces),
     &                                 j=1,max_number_of_ptg_irreps)
                CALL READT(LUINT_LUCITA, nint1, work(kint2))
                CALL READT(LUINT_LUCITA, nint2, work(kint2))
              end if
#ifdef STEFAN_PLEASE_FIX
<<<<<<< HEAD
=======
            end do
          end if
        end if
*
*  NACTEL
        NTOCCH = 0
        NTOCCR = 0
        if (NACTEL.eq.-1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NTOCCH = NTOCCH + 2*NOCC(ISYM)
            end do
            NACTEL = NTOCCH
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NTOCCR = NTOCCR + 2*NISH(ISYM)
            end do
            NOCCIM = NTOCCR + NACRAS
            NACTEL = NOCCIM
          end if
          write(lupri,'(A42,I3/)')
     &      ' Number of active electrons ........... ',NACTEL
        end if
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
         write(lupri,*) 'Number of active electrons exceeds '
         write(lupri,*) 'the orbital space.'
         write(lupri,*) 'Consider Paulis famous principle and restart!'
         Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        NCISPC = 1                  ! might well be made more flexible
        NOCCEL = 0
        if (INACOB.eq.0) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NOCC(ISYM)
            end do
            NIMN = max(0,NOCCEL-2)
            NIMX = NOCCEL
            NEMN = NOCCEL
            NEMX = NOCCEL
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NIMN
              NELGS(1,2,ICISPC) = NIMX
              do IGS=2,NGAS-1,1
                NELGS(IGS,1,ICISPC) = NEMN - 2
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              NOCCEL = NACRAS
            else
              do ISYM=1,NSYM,1
                NOCCEL = NOCCEL + 2*NISH(ISYM)
              end do
>>>>>>> origin/master
#endif  /* STEFAN_PLEASE_FIX */
            end if
          else
!           read 1e- and 2e-integrals from file
!           -----------------------------------
            CALL MOLLAB('LUCIIJKL',LUINT_LUCITA,LUPRI)
            READ (LUINT_LUCITA) nint1, nint2, ecore, ecore_orig ,
     &                          einact_mc2lu
            read(LUINT_LUCITA) 
     &      ((ngsh_lucita_tmp(i,j),i=1,max_number_of_gas_spaces),
     &                             j=1,max_number_of_ptg_irreps)
            CALL READT(LUINT_LUCITA, nint1, work(kint1))
            if(is_2e_int_active .gt. 1)
     &      CALL READT(LUINT_LUCITA, nint2, work(kint2))
          end if
          CALL GPCLOSE(LUINT_LUCITA, 'KEEP')
        ELSE
          CALL QUIT('LUCITA out-of-core not implemented yet, sorry!')
        END IF

        call dcopy(nint1,work(kint1),1,work(kint1o),1)

        IF(IUSE_PH.EQ.1) CALL FI(WORK(KINT1),ECORE_HEX,1)

        ECORE_ORIG = ECORE
        ECORE      = ECORE + ECORE_HEX

#ifdef STEFAN_PLEASE_FIX
<<<<<<< HEAD
#ifdef LUCI_DEBUG
        if(luci_myproc.eq.luci_master)then
          write(lupri,'(/2x,a,e15.8)') 'Updated core energy: ',ECORE
        end if
#endif
=======
!       print *, ' shells per gas space in dalton2lucita:ngas,nsym',
!    &  ngas, nsym
!       do j = 1, ngas
!         print *, ' nshpgs per gas: ',
!    &    (nshpgs(j,i),i=1,nsym)
!       end do

*  NACTEL
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
         write(lupri,*) 'Number of active electrons exceeds '
         write(lupri,*) 'the orbital space.'
         write(lupri,*) 'Consider Paulis famous principle and restart!'
         Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        do ICISPC=1,NCISPC,1
          do IGAS=1,NGAS,1
            do i = 1, 2
              nelgs(igas,i,icispc) = ngso_lucita(igas,i)
            end do
          end do
        end do
        if (NACTEL.ne.NELGS(NGAS,2,NCISPC)) then
          write(lupri,*) 'I am sure that fascinating discussions about'
          write(lupri,*) 'the energy of such a wave function exist,'
          write(lupri,*) 'but I am just a dumb program, so I will stop.'
          write(lupri,*)
          write(lupri,*) 'Number of active electrons does not match'
          write(lupri,*) 'total number of electrons in GAS spaces.'
          Call Abend2( 'Quitting.' )
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(lupri,*) 'Illegal spin multiplicity given.'
           write(lupri,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(lupri,*) 'Illegal spin multiplicity given.'
           write(lupri,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(lupri,*) 'Illegal spin multiplicity given.'
              write(lupri,*) 'Compare with number of active'
              write(lupri,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(lupri,*) 'Illegal spin multiplicity given.'
              write(lupri,*) 'Compare with number of active'
              write(lupri,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        MAXIT = IMAXCIITER
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
*              ................
      else if (TYPE.eq.'RASCI ') then
*              ................
*
*  NGSSH, NGAS
*
!       synchronize module lucita_orbital_spaces (initialized after input ci/mcscf reading) 
!       and internal lucita common blocks

        do isym = 1, nsym
          nshpgs(1,isym) = nas1_lucita(isym)
          nshpgs(2,isym) = nas2_lucita(isym)
          nshpgs(3,isym) = nas3_lucita(isym)
        end do

        ngas         = 3
        ngas_tmp_ras = 3
        do isym = 1, nsym
          if(nshpgs(ngas_tmp_ras,isym) .gt. mxtsob)then
            call gassplit(nshpgs(1,isym),ngas_tmp_ras,mxpngas,
     &                    mxtsob,ngas,iprt)
          end if
        end do
*
*  NACTEL
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
         write(lupri,*) 'Number of active electrons exceeds '
         write(lupri,*) 'the orbital space.'
         write(lupri,*) 'Consider Paulis famous principle and restart!'
         Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        ICISPC = 1
        NCISPC = 1
        MXELR1 = 0
        MXELR2 = 0
        do ISYM=1,NSYM,1
          MXELR1 = MXELR1 + 2*NSHPGS(1,ISYM)
          MXELR2 = MXELR2 + 2*NSHPGS(2,ISYM)
        end do
*
        NFMN = 0
        NFMX = 0
        NIMX = MXELR1 + NFMX
        NIMN = NIMX - MXHOL1
        NEMN = NACTEL + NFMX
        NEMX = NACTEL + NFMX
        NAMX = NEMX
        NAMN = NAMX - MXELR3
*
        do ICISPC=1,NCISPC,1
          NELGS(1,1,ICISPC) = NIMN
          NELGS(1,2,ICISPC) = NIMX
          NELGS(2,1,ICISPC) = NAMN
          NELGS(2,2,ICISPC) = NAMX
          NELGS(3,1,ICISPC) = NEMN
          NELGS(3,2,ICISPC) = NEMX
        end do
        NACTEL = NACTEL+NFMX
        if (NACTEL.ne.NELGS(NGAS,2,NCISPC)) then
          write(lupri,*) 'I am sure that fascinating discussions about'
          write(lupri,*) 'the energy of such a wave function exist,'
          write(lupri,*) 'but I am just a dumb program, so I will stop.'
          write(lupri,*)
          write(lupri,*) 'Number of active electrons does not match'
          write(lupri,*) 'total number of electrons in GAS spaces.'
          Call Abend2( 'Quitting.' )
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(lupri,*) 'Illegal spin multiplicity given.'
           write(lupri,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(lupri,*) 'Illegal spin multiplicity given.'
           write(lupri,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(lupri,*) 'Illegal spin multiplicity given.'
              write(lupri,*) 'Compare with number of active electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(lupri,*) 'Illegal spin multiplicity given.'
              write(lupri,*) 'Compare with number of active electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
>>>>>>> origin/master
#endif  /* STEFAN_PLEASE_FIX */
      else

      end if ! citask_id switch

      END
***********************************************************************

      SUBROUTINE PUTINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM)
*
* Put integrals in permanent integral list
*
* Jeppe Olsen, Jan. 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*. Specific input
      DIMENSION XINT(*)
*
      CALL QENTER('PUTIN')
*. Offset and number of integrals
*
      IF(ITP.EQ.0) THEN
        NI = NTOOBS(ISM)
      ELSE
        NI = NOBPTS(ITP,ISM)
      END IF
*
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      IF(JTP.EQ.0) THEN
        NJ = NTOOBS(JSM)
      ELSE
        NJ = NOBPTS(JTP,JSM)
      END IF
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      IF(KTP.EQ.0) THEN
        NK = NTOOBS(KSM)
      ELSE
        NK = NOBPTS(KTP,KSM)
      END IF
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      IF(LTP.EQ.0) THEN
        NL = NTOOBS(LSM)
      ELSE
        NL = NOBPTS(LTP,LSM)
      END IF
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO
*
      INT_IN = 0
      DO LOB = LOFF,LOFF+NL-1
       DO KOB = KOFF,KOFF+NK-1
        DO JOB = JOFF,JOFF+NJ-1
         DO IOB = IOFF,IOFF+NI-1
C?         WRITE(6,*) ' IOB, JOB, KOB, LOB', IOB,JOB,KOB,LOB
           INT_OUT = I2EAD(IOB,JOB,KOB,LOB)
           INT_IN = INT_IN + 1
C?         WRITE(6,*) ' INT_OUT, INT_IN ', INT_OUT,INT_IN
C?         WRITE(6,*) ' KINT2-1+INT_OUT = ',KINT2-1+INT_OUT
           WORK(KINT2-1+INT_OUT) = XINT(INT_IN)
         END DO
        END DO
       END DO
      END DO
*
      CALL QEXIT('PUTIN')
      END
      SUBROUTINE MAKE_LUCITA_INTEGRALS(CMO,WORK,LWORK)
C
C     Dec 2010 : get EMY and get FCAC and H2AC matrices in Dalton
C
      use lucita_cfg
#include "implicit.h"
      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX
#include "maxorb.h"
#include "infinp.h"

      DIMENSION CMO(*), WORK(LWORK)
 
C Used from include files:
C  priunit: LUPRI
C  inforb: NCMOT, nsym
C  inftap: LUINTM, FNINTM
C  mxpdim: MXPOBS
C  orbinp: ntoobs
#include "priunit.h"
#include "inforb.h"
#include "inftap.h"


      LOGICAL CLOSE_MOTWOINT, DISKH2
C
C     -------------------------------------------------------
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
C
      IF (LUINTM .LE. 0) THEN
          CALL GPOPEN(LUINTM,FNINTM,'OLD',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
          CLOSE_MOTWOINT = .TRUE.
      ELSE
          CLOSE_MOTWOINT = .FALSE.
      END IF

!     allocate scratch memory and read-in 1/2-electron integrals from dalton-sirius files
      call memget('REAL',kfcac,nnashx   ,WORK,kfree,lfree)
      call memget('REAL',kh2ac,nnashx**2,WORK,kfree,lfree)

      DISKH2 = .FALSE.
      CALL CIINTS(DISKH2,CMO,EMY,WORK(kfcac),
     &            WORK(kh2ac),WORK(kfree),lfree)

      IF (CLOSE_MOTWOINT) THEN
         CALL GPCLOSE(LUINTM,'KEEP')
      END IF


#ifdef LUCI_DEBUG
      WRITE(LUPRI,'(//A)') 'MAKE_LUCITA_INTEGRALS: FCAC matrix'
      CALL OUTPAK(work(kfcac),NASHT,-1,LUPRI)
      WRITE(LUPRI,'(//A)') 'MAKE_LUCITA_INTEGRALS: H2AC matrix'
      CALL OUTPUT(work(kh2ac),1,NNASHX,1,NNASHX,NNASHX,
     &            NNASHX,-1,LUPRI)
#endif

!     Write one- and two-electron integrals to file FCIDUMP
 
      if (lucita_cfg_fci_dump)
     &   call write_fcidump_file(work(kfcac),NASHT,work(kh2ac),EMY)

      LUINT_LUCITA = -1
      CALL GPOPEN(LUINT_LUCITA,'INTEGRALS_LUCITA','NEW',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'LUCIINTS'
      WRITE (LUINT_LUCITA) EMY, NNASHX, NSYM, DUM, DUM
      CALL WRITT(LUINT_LUCITA, NNASHX, WORK(KFCAC))
      CALL WRITT(LUINT_LUCITA, NNASHX*NNASHX, WORK(KH2AC))
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'EODATA  '
      CALL GPCLOSE(LUINT_LUCITA, 'KEEP')
!     release scratch space
      call memrel('lucita-dalton-1/2int-write.done',WORK,
     &            kfrsav,kfrsav,kfree,lfree)

      END ! SUBROUTINE MAKE_LUCITA_INTEGRALS
***********************************************************************

      SUBROUTINE LUCITA_GETINT_DALTON(EMY,X1INT,X2INT,
     &                                tmp_dal_work,LWORK,
     &                                int1_mcscf,
     &                                int2_mcscf,
     &                                mcscf_provides_integrals,iintden)
C
C     Dec 2010 : get EMY and get FCAC and H2AC matrices in Dalton
C
#include "implicit.h"

      DIMENSION X1INT(*), X2INT(*), tmp_dal_work(lwork)
      real(8), intent(inout) :: int1_mcscf(*)
      real(8), intent(inout) :: int2_mcscf(*)
      integer, intent(in)    :: iintden
      logical, intent(in)    :: mcscf_provides_integrals
 
C Used from include files:
C  priunit: LUPRI
C  inforb: NCMOT, nsym
C  inftap: LUINTM, FNINTM
C  mxpdim: MXPOBS
C  orbinp: ntoobs
C  cgas: ngas
#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "oper.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

      LOGICAL CLOSE_MOTWOINT, DISKH2
!     scratch space
      integer, allocatable :: mult_ptg(:,:)
      integer              :: is2e_active
C
C     -------------------------------------------------------
C
!     allocate scratch matrix
      allocate(mult_ptg(nsym,nsym))

!     set multiplication table
      mult_ptg = 0
      call set_ptg_multiplication_table(mult_ptg,nsym)
      
      is2e_active = i12

      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK

!     write(lupri,*) 'mcscf_provides_integrals = ',
!    &                mcscf_provides_integrals

      if(.not.mcscf_provides_integrals)then
!       allocate scratch memory and read-in 1/2-electron integrals from dalton-sirius files
        call memget('REAL',kfcac,nnashx   ,tmp_dal_work,kfree,lfree)
        call memget('REAL',kh2ac,nnashx**2,tmp_dal_work,kfree,lfree)

        LUINT_LUCITA = -1
        CALL GPOPEN(LUINT_LUCITA,'INTEGRALS_LUCITA','OLD',' ',
     &                  'UNFORMATTED',IDUMMY,.FALSE.)
        CALL MOLLAB('LUCIINTS',LUINT_LUCITA,LUPRI)
        READ (LUINT_LUCITA) EMY, NNASHX_x, NSYM_x
        IF (NNASHX .NE. NNASHX_x) THEN
           call quit('NNASHX in common .ne. NNASHX on INTEGRALS_LUCITA')
        END IF
        IF (NSYM .NE. NSYM_x) THEN
           call quit('NSYM in common .ne. NSYM on INTEGRALS_LUCITA')
        END IF
        CALL READT(LUINT_LUCITA, NNASHX, tmp_dal_WORK(KFCAC))
        CALL READT(LUINT_LUCITA, NNASHX*NNASHX, tmp_dal_WORK(KH2AC))
        CALL GPCLOSE(LUINT_LUCITA, 'KEEP')

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
        WRITE(LUPRI,'(//A)') 'LUCITA_GETINT_DALTON: FCAC matrix'
        WRITE(LUPRI,'(A,2i8)') 
     &              'dimension: NASHT,NASHT * NASHT',NASHT,NASHT**2
        CALL OUTPAK(tmp_dal_work(kfcac),NASHT,-1,LUPRI)
        WRITE(LUPRI,'(//A)') 'LUCITA_GETINT_DALTON: H2AC matrix'
        WRITE(LUPRI,'(A,2i8)') 
     &              'dimension: NNASHX, NNASHX**2',NNASHX,NNASHX**2
        CALL OUTPUT(tmp_dal_work(kh2ac),1,NNASHX,1,NNASHX,NNASHX,
     &              NNASHX,-1,LUPRI)
#endif
!#undef LUCI_DEBUG

!--------------------------------------------------------------
!     transfer 1e- and 2e-integrals from Dalton to LUCITA order
!--------------------------------------------------------------
        call sirint2luint(tmp_dal_work(kfcac),tmp_dal_work(kh2ac),
     &                    dummy,x1int,x2int,mult_ptg,adsxa,nasht,nnashx,
     &                    ntoobs,nsym,mxpobs,nint1,nint2,ireost,ireots,
     &                    ngas,is2e_active,iintden,-1)
!       release scratch space
        call memrel('lucita-dalton-1/2int-interface.done',tmp_dal_work,
     &              kfrsav,kfrsav,kfree,lfree)

      else ! mcscf_provides_integrals
!       write(lupri,*) 'distribute ints from mc environment...'
        call sirint2luint(int1_mcscf,int2_mcscf,
     &                    dummy,x1int,x2int,mult_ptg,adsxa,nasht,nnashx,
     &                    ntoobs,nsym,mxpobs,nint1,nint2,ireost,ireots,
     &                    ngas,is2e_active,iintden,-1)
      end if

      deallocate(mult_ptg)

      END
***********************************************************************

      subroutine lucita_srdft_h1_adapt(X1INT,x1int_scr,int1_mcscf)
C
#include "implicit.h"

      real(8), intent(inout) :: x1int(*)
      real(8), intent(inout) :: x1int_scr(*)
      real(8), intent(inout) :: int1_mcscf(*)
 
C Used from include files:
C  priunit: LUPRI
C  inforb: NCMOT, nsym
C  inftap: LUINTM, FNINTM
C  mxpdim: MXPOBS
C  orbinp: ntoobs
C  cgas: ngas
#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

!     scratch space
      integer, allocatable :: mult_ptg(:,:)
      integer              :: is2e_active
      integer              :: iintden = 1
C
C     -------------------------------------------------------
C
!     allocate scratch matrix
      allocate(mult_ptg(nsym,nsym))

!     set multiplication table
      mult_ptg = 0
      call set_ptg_multiplication_table(mult_ptg,nsym)
      
      is2e_active = 1 ! no 2-electron integrals

!-----------------------------------------------------------------
!     transfer 1e- srDFT contributions from Dalton to LUCITA order
!-----------------------------------------------------------------
      call sirint2luint(x1int,dummy,              
     &                  dummy,x1int_scr,dummy,mult_ptg,
     &                  adsxa,nasht,nnashx,
     &                  ntoobs,nsym,mxpobs,nint1,nxxxx,ireost,ireots,
     &                  ngas,is2e_active,iintden,-1)

!     scale 1-electron integrals with the sr-dft contributions
      call daxpy(nint1,1.0d0,x1int_scr,1,int1_mcscf,1)

      deallocate(mult_ptg)

      END
***********************************************************************

      subroutine lucita_putdens_generic(rho1lu,rho2lu,rho1gen,rho2gen,
     &                                  pv_full_scratch,
     &                                  isrho2_active,iintden,rhotype,
     &                                  state_id)
!
!     Sep 2011 : driver for sorting the 1-/2-particle density matrices from 
!                LUCITA to a generic format (here: SIRIUS)
!
#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
#endif
#include "implicit.h"
#include "parluci.h"

      real(8), intent(in)    :: rho1lu(*)
      real(8), intent(in)    :: rho2lu(*)
      real(8), intent(out)   :: rho1gen(*)
      real(8), intent(out)   :: rho2gen(*)
      real(8), intent(inout) :: pv_full_scratch(*)
      integer, intent(in)    :: iintden
      integer, intent(in)    :: isrho2_active
      integer, intent(in)    :: rhotype
      integer, intent(in)    :: state_id
 
! Used from include files:
!  priunit: LUPRI
!  inforb: NCMOT, nsym
!  inftap: LUINTM, FNINTM
!  mxpdim: MXPOBS
!  orbinp: ntoobs

#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

!     scratch space
      integer, allocatable :: mult_ptg(:,:)
      integer              :: noint_tmp
      character (len=25)   :: file_name
!
!     -------------------------------------------------------
!
!     allocate scratch matrix
      allocate(mult_ptg(nsym,nsym))

!---------------------------------------------------------------------
!     transfer 1e- and 2e-density matrices from LUCITA to SIRIUS order
!---------------------------------------------------------------------

!     set multiplication table
      mult_ptg = 0
      call set_ptg_multiplication_table(mult_ptg,nsym)

      call dzero(rho1gen,nnashx**1)
      if(isrho2_active .gt. 1) call dzero(rho2gen,nnashx**2)

      call sirint2luint(rho1gen,rho2gen,pv_full_scratch,
     &                  rho1lu,rho2lu,mult_ptg,adsxa,nasht,nnashx,
     &                  ntoobs,nsym,mxpobs,nint1,nint2,ireost,
     &                  ireots,ngas,isrho2_active,iintden,
     &                  rhotype)

      deallocate(mult_ptg)

!     put 1e- and 2e-density matrices to file
!     ---------------------------------------

      write(file_name,'(a17,i3,a1,i4)') 
     &      'DENSITIES_LUCITA-',state_id,'.',luci_myproc
      do i = 18, 20
        if(file_name(i:i).eq. ' ') file_name(i:i) = '0'
      end do
      do i = 22, 25
        if(file_name(i:i).eq. ' ') file_name(i:i) = '0'
      end do

      LUINT_LUCITA = -1
      CALL GPOPEN(LUINT_LUCITA,file_name,'UNKNOWN',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'LUCIDENS'
        CALL WRITT(LUINT_LUCITA, NNASHX, rho1gen)
      if(isrho2_active .gt. 1) 
     &CALL WRITT(LUINT_LUCITA, NNASHX*NNASHX, rho2gen)
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'EODATA  '
      CALL GPCLOSE(LUINT_LUCITA, 'KEEP')

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(LUPRI,'(//A)') 
     & 'lucita_putdens_generic: 1-particle density matrix'
      WRITE(LUPRI,'(A,2i8)') 
     &            'dimension: NASHT,NASHT * NASHT',NASHT,NASHT**2
      CALL OUTPAK(rho1gen,NASHT,-1,LUPRI)
      WRITE(LUPRI,'(//A)') 
     & 'lucita_putdens_generic: 2-particle density matrix'
      WRITE(LUPRI,'(A,2i8)') 
     &            'dimension: NNASHX, NNASHX**2',NNASHX,NNASHX**2
      CALL OUTPUT(rho2gen,1,NNASHX,1,NNASHX,NNASHX,
     &            NNASHX,-1,LUPRI)
#endif
!#undef LUCI_DEBUG

      END
***********************************************************************

      subroutine lucita_spinputdens_1p(rho1slu,rho2lu,rho1gen,rho2gen,
     &                                 pv_full_scratch,
     &                                 isrho2_active,iintden,rhotype,
     &                                 state_id,iab)
!
!     Aug 2012 : driver for sorting the 1-particle spin-density matrices from 
!                LUCITA to a generic format (here: SIRIUS)
!
#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
      use mpi
#else
#include "mpif.h"
#endif
#endif
#include "implicit.h"
#include "parluci.h"

      real(8), intent(in)    :: rho1slu(*)
      real(8), intent(in)    :: rho2lu(*)
      real(8), intent(out)   :: rho1gen(*)
      real(8), intent(out)   :: rho2gen(*)
      real(8), intent(inout) :: pv_full_scratch(*)
      integer, intent(in)    :: iintden
      integer, intent(in)    :: isrho2_active
      integer, intent(in)    :: rhotype
      integer, intent(in)    :: state_id
      integer, intent(in)    :: iab
 
! Used from include files:
!  priunit: LUPRI
!  inforb: NCMOT, nsym
!  inftap: LUINTM, FNINTM
!  mxpdim: MXPOBS
!  orbinp: ntoobs

#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

!     scratch space
      integer, allocatable :: mult_ptg(:,:)
      integer              :: noint_tmp
      character (len=25)   :: file_name
!
!     -------------------------------------------------------
!
!     allocate scratch matrix
      allocate(mult_ptg(nsym,nsym))

!---------------------------------------------------------------------
!     transfer 1e-spin density matrix from LUCITA to Sirius order
!---------------------------------------------------------------------

!     set multiplication table
      mult_ptg = 0
      call set_ptg_multiplication_table(mult_ptg,nsym)

      call dzero(rho1gen,nnashx**1)
      if(isrho2_active .gt. 1) call dzero(rho2gen,nnashx**2)

      call sirint2luint(rho1gen,rho2gen,pv_full_scratch,
     &                  rho1slu,rho2lu,mult_ptg,adsxa,nasht,nnashx,
     &                  ntoobs,nsym,mxpobs,nint1,nint2,ireost,
     &                  ireots,ngas,isrho2_active,iintden,
     &                  rhotype)

      deallocate(mult_ptg)

!     put 1e- and 2e-density matrices to file
!     ---------------------------------------

      select case(iab)
        case(1) ! alpha spin
          write(file_name,'(a17,i3,a1,i4)') 
     &       'SPINDENSITY-A-1p-',state_id,'.',luci_myproc
        case(2) ! beta  spin
          write(file_name,'(a17,i3,a1,i4)') 
     &        'SPINDENSITY-B-1p-',state_id,'.',luci_myproc
        case default ! ??? ask einstein or dirac
          stop ' wrong program/hamiltonian - spin is a good quantum num'
      end select

      do i = 18, 20
        if(file_name(i:i).eq. ' ') file_name(i:i) = '0'
      end do
      do i = 22, 25
        if(file_name(i:i).eq. ' ') file_name(i:i) = '0'
      end do

      LUINT_LUCITA = -1
      CALL GPOPEN(LUINT_LUCITA,file_name,'UNKNOWN',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'LUCIDENS'
        CALL WRITT(LUINT_LUCITA, NNASHX, rho1gen)
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'EODATA  '
      CALL GPCLOSE(LUINT_LUCITA, 'KEEP')

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(LUPRI,'(//A,i2)') 
     & 'lucita_putdens_generic: 1-particle spin density matrix :',iab
      WRITE(LUPRI,'(A,2i8)') 
     &            'dimension: NASHT, triangular packed',NASHT
      CALL OUTPAK(rho1gen,NASHT,-1,LUPRI)
#undef LUCI_DEBUG
#endif

      END
***********************************************************************

      subroutine lucita_putdens_ensemble(rho1_ensemble)
!
#include "implicit.h"

      real(8), intent(in)    :: rho1_ensemble(*)
 
! Used from include files:
!  priunit: LUPRI
!  inforb: NCMOT, nsym
!  inftap: LUINTM, FNINTM
!  mxpdim: MXPOBS
!  orbinp: ntoobs

#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

!     scratch space
      character (len=25)   :: file_name
!
!     -------------------------------------------------------
!
!     put 1e-density matrix to file
!     -----------------------------

      LUINT_LUCITA = 99
      write(file_name,'(a16)') 
     &      'ENSEMBLE-DENSITY'

         open(file=file_name(1:16),unit=luint_lucita,status='replace',
     &        form='unformatted',action='write',position='rewind')

      WRITE (LUINT_LUCITA) ('********',I=1,3), 'LUCIDENS'
        CALL WRITT(LUINT_LUCITA, NNASHX, rho1_ensemble)
      WRITE (LUINT_LUCITA) ('********',I=1,3), 'EODATA  '

      close(luint_lucita,status='keep')

      END
***********************************************************************

      subroutine sirint2luint(h1_in,h2_in,pv_full_scratch,h1_out,h2_out,
     &                        mult_ptg,return_sym_of_jorb,nr_actorb_tot,
     &                        nnashx,nr_orb_per_sym,nr_sym,
     &                        max_orb_per_sym,nr_1int,nr_2int,
     &                        ireost,ireots,nr_gas,is2e_active,iintden,
     &                        rhotype)
!
!     purpose: iintden == 1: reorder integrals from dalton-sirius format to internal
!                            lucita format.
!                            1-electron integrals: on output symmetry-reduced list
!                            2-electron integrals: on output symmetry-reduced list
!              iintden == 2: reorder particle-density matrices from to internal lucita format to 
!                            dalton-sirius format.
!                            rhotype controls the format:
!                            ==> 1: rho2(ij,kl) density matrix
!                            ==> 2: rho2(i,j,k,l) used for response
!                            ==> 3: rho2(ij,kl) symmetrized transition density matrix
!
!     -----------------------------------------------------------------
      implicit none
#include "priunit.h"
      integer :: nr_actorb_tot
      integer :: nnashx
      integer :: nr_sym
      integer :: max_orb_per_sym
      integer :: nr_1int
      integer :: nr_2int
      integer :: return_sym_of_jorb(max_orb_per_sym,max_orb_per_sym*2)
      integer :: nr_orb_per_sym(max_orb_per_sym)
      integer :: mult_ptg(nr_sym,nr_sym)
      integer :: iintden
      integer :: nr_gas
      integer :: is2e_active
      integer :: ireost(*)
      integer :: ireots(*)
      integer :: rhotype
!     real(8) :: h1_in(nnashx)
!     real(8) :: h2_in(nnashx,nnashx)
      real(8) :: h1_in(*)
      real(8) :: h2_in(*)
      real(8) :: pv_full_scratch(nr_actorb_tot,nr_actorb_tot,
     &                           nr_actorb_tot,nr_actorb_tot)
      real(8) :: h1_out(*)
      real(8) :: h2_out(*)
!     ------------------------------------------------------------------

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      if(iintden .eq. 2)then
        write(lupri,*) '  1-el density matrix'
        call wrtmatmn(h1_out,1,nr_actorb_tot**2,1,
     &                nr_actorb_tot**2,lupri)
        if(is2e_active .gt. 1)then
         write(lupri,*) '  2-el density matrix'
         call wrtmatmn(h2_out,1,nr_actorb_tot**2*(nr_actorb_tot**2+1)/2,
     &                 1,nr_actorb_tot**2*(nr_actorb_tot**2+1)/2,lupri)
        end if
      end if
#endif

!     1-electron integrals/density matrix
!     -----------------------------------
      call sir1int2lu1int(h1_in,h1_out,return_sym_of_jorb,nnashx,
     &                    nr_orb_per_sym,nr_sym,max_orb_per_sym,
     &                    nr_1int,nr_actorb_tot,ireost,iintden,rhotype)
      
!     2-electron integrals/density matrix
!     -----------------------------------
      if(is2e_active .gt. 1)then
        if(iintden .eq. 1)then
          call sir2int2lu2int(h2_in,h2_out,mult_ptg,nr_actorb_tot,
     &                        nnashx,nr_orb_per_sym,nr_sym,
     &                        max_orb_per_sym,nr_2int)
        else
          call lu2dens2sir2dens(h2_in,h2_out,pv_full_scratch,
     &                          nr_actorb_tot,ireost,ireots,
     &                          nr_sym,nr_gas,rhotype)
        end if
      end if

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      if(iintden .eq. 1)then
        write(lupri,*) '  symmetry reduced 1-el'
        call wrtmatmn(h1_out,1,nr_1int,1,nr_1int,lupri)
        if(is2e_active .gt. 1)then
          write(lupri,*) '  symmetry reduced 2-el'
          call wrtmatmn(h2_out,1,nr_2int,1,nr_2int,lupri)
        end if
      else
        write(lupri,*) '  1-el density matrix'
        call outpak(h1_in,nr_actorb_tot,-1,lupri)
        if(is2e_active .gt. 1)then
          write(lupri,*) '  2-el density matrix'
          call output(h2_in,1,nnashx,1,nnashx,nnashx,
     &                nnashx,-1,lupri)
        end if
      end if
#undef LUCI_DEBUG
#endif

      end

***********************************************************************
      subroutine sir1int2lu1int(h1_in,h1_out,return_sym_of_jorb,nnashx,
     &                          nr_orb_per_sym,nr_sym,max_orb_per_sym,
     &                          nr_1int,nr_actorb_tot,ireost,iintden,
     &                          rhotype)
!
!     purpose: reorder integrals from dalton-sirius format to internal
!              lucita format.
!
!              1-electron integrals: on output symmetry-reduced list
!
!              note from stefan: i am not sure whether my code actually 
!              works for non-totally symmetric 1-electron operators as
!              well, i.e. if we do property calculations. 
!              this should be checked with due care.
!
!     -----------------------------------------------------------------
      implicit none
#include "priunit.h"
      integer :: nnashx
      integer :: nr_sym
      integer :: max_orb_per_sym
      integer :: nr_1int
      integer :: return_sym_of_jorb(max_orb_per_sym,max_orb_per_sym*2)
      integer :: nr_orb_per_sym(max_orb_per_sym)
      integer :: nr_actorb_tot
      integer :: ireost(*)
      integer :: iintden
      integer :: rhotype
      real(8) :: h1_in(*)
      real(8) :: h1_out(*)
!     -----------------------------------------------------------------
      integer :: isym, jsym, ijsym
      integer :: offset_ij, offset_internal, isorb, jsorb
      integer :: nr_act_isym, nr_act_jsym, orb_tmp
      integer :: i_index, j_index
      integer :: loop_counter_i
      character (len=90) :: error_message
!     -----------------------------------------------------------------

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      if(rhotype > 0)then
        write(lupri,*) ' distributing 1p-dens mat elements from h1_out:'
        call wrtmatmn(h1_out,nr_actorb_tot,nr_actorb_tot,
     &                       nr_actorb_tot,nr_actorb_tot,lupri)
      else
        write(lupri,*) ' distributing 1e-ints from h1_in:'
        do isym = 1, nnashx
          write(lupri,*) ' h1_in(',isym,') = ',h1_in(isym)
        end do
      end if
!#undef LUCI_DEBUG
#endif

!     1-electron integrals/density matrix elements
!     --------------------------------------------
      orb_tmp         = 0
      offset_ij       = 0
      offset_internal = 1
      ijsym           = 1

!     insert check for symmetry of 1-electron operator - quit if not
!     totally symmetric for the time being
      if(ijsym .gt. 1)then
        error_message = '*** error in integral resorting for dalton'//
     &                  ' 1-e ints. operator is not totally symmetric.'
        call quit(error_message)
      end if

      do isym = 1, nr_sym


        nr_act_isym = nr_orb_per_sym(isym)
        jsym        = return_sym_of_jorb(isym,ijsym)

#ifdef LUCI_DEBUG
        write(lupri,*) ' isym, jsym',isym, jsym
        write(lupri,*) ' nr_act_isym', nr_act_isym
#endif

        if(jsym >= isym)then

          nr_act_jsym = nr_orb_per_sym(jsym)

#ifdef LUCI_DEBUG
          write(lupri,*) ' nr_act_jsym', nr_act_jsym
#endif

          do jsorb = 1, nr_act_jsym
            loop_counter_i = jsorb
            if(rhotype > 0) loop_counter_i = nr_act_isym
            do isorb = 1, loop_counter_i
!             calculate offset on dalton 1e-array (lower triangular matrix)
              i_index   =  (orb_tmp+isorb)
              j_index   = ((orb_tmp+jsorb)*(orb_tmp+jsorb-1)/2)
!             no triangular packing if rhotype > 1              
              if(rhotype > 0) j_index   = ((orb_tmp+jsorb-1))
!    &                                  * ( orb_tmp+jsorb  ))
     &                                  * nr_actorb_tot
              offset_ij = i_index + j_index
  
              if(iintden .eq. 1)then ! 1e integral
                h1_out(offset_internal) = h1_in(offset_ij)
              else ! 1p density matrix element
!
!               LUCITA: type-symmetry ordering <==> SIRIUS (or other programs) symmetry-type
!                       array ireost takes care ot that...
!
                offset_internal  = 
     &          (ireost(isorb+orb_tmp)-1) * nr_actorb_tot + 
     &           ireost(jsorb+orb_tmp) 
                if(rhotype > 1) offset_internal = 
     &          (ireost(jsorb+orb_tmp)-1) * nr_actorb_tot + 
     &           ireost(isorb+orb_tmp) 
     
                h1_in(offset_ij) = h1_out(offset_internal)
              end if
  
#ifdef LUCI_DEBUG
       write(lupri,*)'offset_ij, offset_int,jsorb,isorb,h1,switch',
     & offset_ij, offset_internal,jsorb,isorb,
     & h1_out(offset_internal), h1_in(offset_ij), iintden
#endif

!             count the number of non-vanishing 1-electron integrals
              offset_internal = offset_internal + 1
            end do
          end do
        end if ! jsym >= isym

!       total # of active orbitals for each isym
        orb_tmp = orb_tmp + nr_act_isym
      end do

!     final # of 1-el ints
      nr_1int = offset_internal - 1

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      if(rhotype <= 0)then
        write(lupri,*) '  final 1e- list: # of ints =', nr_1int
        call wrtmatmn(h1_out,1,nr_1int,1,nr_1int,lupri)
      end if
#undef LUCI_DEBUG
#endif

!     symmetrize 1p-dens matrix
      if(rhotype >  0)then

!       write(lupri,*) '  temp 1p-tdm: '
!       call wrtmatmn(h1_in,nr_actorb_tot,nr_actorb_tot,
!    &                nr_actorb_tot,nr_actorb_tot,lupri)

        call trpad(h1_in,1.0d0,nr_actorb_tot)

!       write(lupri,*) '  temp 1p-tdm: part 2'
!       call wrtmatmn(h1_in,nr_actorb_tot,nr_actorb_tot,
!    &                nr_actorb_tot,nr_actorb_tot,lupri)

        if(rhotype == 1) call dscal(nr_actorb_tot**2,0.5d0,h1_in,1)

!       write(lupri,*) '  temp 1p-tdm: part 3'
!       call wrtmatmn(h1_in,nr_actorb_tot,nr_actorb_tot,
!    &                nr_actorb_tot,nr_actorb_tot,lupri)

!       pack symmetrized 1p-dens matrix (h1_in) as triangular matrix h1_out
        call dsitsp(nr_actorb_tot,h1_in,h1_out)
!       save triangular matrix in h1_in
        call dzero(h1_in,nr_actorb_tot**2)
        call dcopy((nr_actorb_tot*(nr_actorb_tot+1))/2,h1_out,1,h1_in,1)
      end if

      end
***********************************************************************

      subroutine sir2int2lu2int(h2_in,h2_out,mult_ptg,nr_actorb_tot,
     &                          nnashx,nr_orb_per_sym,nr_sym,
     &                          max_orb_per_sym,nr_2int)
!
!     purpose: reorder integrals from dalton-sirius format to internal
!              lucita format.
!
!              2-electron integrals: on output symmetry-reduced list
!
!     -----------------------------------------------------------------
      implicit none
#include "priunit.h"
      integer :: nr_actorb_tot
      integer :: nnashx
      integer :: nr_sym
      integer :: max_orb_per_sym
      integer :: nr_2int
      integer :: nr_orb_per_sym(max_orb_per_sym)
      integer :: mult_ptg(nr_sym,nr_sym)
      real(8) :: h2_in(nnashx,nnashx)
      real(8) :: h2_out(*)
!     -----------------------------------------------------------------
      integer :: isym, jsym, ksym, lsym, llsym, ijsym, klsym, ijklsym
      integer :: offset_internal, offset_ijkl
      integer :: nr_act_isym, nr_act_jsym, nr_act_ksym, nr_act_lsym
      integer :: ioff, joff, koff, loff
      integer :: ij_res, kl_res, ijkl_res
      integer :: i_first, j_first, j_last, l_last
      integer :: i, j, k, l, ij, kl
!     -----------------------------------------------------------------

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      print *,' h2_in:'
      do isym = 1, nnashx
        do jsym = 1, nnashx
          print *,' h2_in(',isym,jsym,') = ',
     &              h2_in(isym,jsym)
        end do
      end do
      print *,' mult_ptg:'
      do isym = 1, nr_sym
        do jsym = 1, nr_sym
          print *,' mult_ptg(',isym,jsym,') = ',
     &              mult_ptg(isym,jsym)
        end do
      end do
#endif

!     2-electron integrals/2-particle density matrix elements
!     -------------------------------------------------------
      offset_internal = 1
      ioff            = 1
!     4-index loop
      do isym = 1, nr_sym
        nr_act_isym = nr_orb_per_sym(isym)
        joff = 1
        do jsym = 1, isym
          nr_act_jsym = nr_orb_per_sym(jsym)
          koff = 1
          do ksym = 1, isym
            nr_act_ksym = nr_orb_per_sym(ksym)
            llsym = ksym

            if(ksym .eq. isym) llsym = jsym

            loff  = 1
            do lsym = 1, llsym
              nr_act_lsym = nr_orb_per_sym(lsym)

              ijsym   = mult_ptg(isym,jsym)
!             print *, ' isym, jsym',isym,jsym
              klsym   = mult_ptg(ksym,lsym)
!             print *, ' ksym, lsym',ksym,lsym
              ijklsym = mult_ptg(ijsym,klsym)

!             restrict loop to non-vanishing symmetry combinations of (ij|kl)

              if(ijklsym .eq. 1)then

#ifdef LUCI_DEBUG
                print *, ' check current symmetries'
                print('(/a,4i6)'), 'isym,jsym,ksym,lsym',
     &                              isym,jsym,ksym,lsym
#endif

!               index restrictions
                if(isym .eq. jsym)then
                  ij_res = 1
                else
                  ij_res = 0
                end if
                if(ksym .eq. lsym)then
                  kl_res = 1
                else
                  kl_res = 0
                end if
 
!               particle symmetry
                if(isym.eq.ksym.and.jsym.eq.lsym)then
                  ijkl_res = 1
                else
                  ijkl_res = 0
                end if
!
! K
!  L
!   I
!    J
!               loop over non-redundant output indices including symmetry offsets
                do k = koff, koff + nr_act_ksym - 1

                  if(kl_res .eq. 1)then
                    l_last = k
                  else
                    l_last = loff + nr_act_lsym - 1
                  end if

                  do l = loff, l_last

                    if(ijkl_res .eq. 1)then
                      i_first = k
                    else
                      i_first = ioff
                    end if

                    do i = i_first, ioff + nr_act_isym - 1

                      if(ij_res .eq. 1)then
                        j_last = i
                      else
                        j_last = joff + nr_act_jsym - 1
                      end if

                      if(ijkl_res .eq. 1 .and. i .eq. k)then
                        j_first = l
                      else
                        j_first = joff
                      end if

                      do j = j_first, j_last

!                       indices on SIRIUS array
                        ij                      = (i*(i-1))/2 + j
                        kl                      = (k*(k-1))/2 + l

!                       pick element
                        h2_out(offset_internal) = h2_in(ij,kl)

#ifdef LUCI_DEBUG
                        print ('(/2x,a)'),
     &        '-------------------------------------------------------'
                        print ('(2x,a,2i4,a,2i4,a)'),'next integral:'//
     &                    ' (ij|kl) = (',i,j,'|',k,l,')'
                        print *, 'ij, kl =',ij, kl
                        print *, 'offset_internal =',offset_internal
                        print *, 'h2_out(offset_internal) =',
     &                            h2_out(offset_internal)
                        print('(2x,a/)'),
     &        '-------------------------------------------------------'
#endif
                        offset_internal         = offset_internal + 1
                      end do ! j
                    end do ! i
                  end do ! l
                end do ! k
              end if ! non-vanishing non-redundant (wrt symmetry) integrals
              loff = loff + nr_act_lsym
            end do ! lsym
            koff = koff + nr_act_ksym
          end do ! ksym
          joff = joff + nr_act_jsym
        end do ! jsym
        ioff = ioff + nr_act_isym
      end do ! isym

!     final # of 2-el ints
      nr_2int = offset_internal - 1

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(lupri,*) '  final 2e- list: # of ints =', nr_2int
      call wrtmatmn(h2_out,1,nr_2int,1,nr_2int,lupri)
#undef LUCI_DEBUG
#endif

      end
***********************************************************************

      subroutine lu2dens2sir2dens(pv_sirius,pv_lu,pv_full_scratch,nasht,
     &                            ireost,ireots,nr_sym,nr_gas,rhotype)
!
!     purpose: reorder 2e-RDM from dalton-sirius format to internal
!              lucita format.
!
!              2e-RDM: on output symmetry-reduced list
!
!     -----------------------------------------------------------------
      implicit none
#include "priunit.h"
      integer, intent(in)    :: nasht
      integer, intent(in)    :: nr_sym
      integer, intent(in)    :: nr_gas
      integer, intent(in)    :: rhotype
      integer, intent(in)    :: ireost(*)
      integer, intent(in)    :: ireots(*)
      real(8), intent(out)   :: pv_sirius((nasht*(nasht+1))/2,
     &                                  (nasht*(nasht+1))/2)
      real(8), intent(in)    :: pv_lu(*)
      real(8), intent(inout) :: pv_full_scratch(nasht,nasht,nasht,nasht)
!     -----------------------------------------------------------------
      integer                :: nnashx
      integer                :: n2ashx
      integer                :: kl_ST_max, kl_ST_min, kl_ST_index
      integer                :: l, k, k_ST, l_ST
      integer                :: i, j, ijkl,klij, ij,kl
      real(8), allocatable   :: temp_mat(:)
      logical                :: do_reorder
      real(8)                :: value
!     -----------------------------------------------------------------

      nnashx = (nasht*(nasht+1))/2
      n2ashx = nasht**2

!     expand the lower-triangle 2-particle density matrix (lu format) to the full one
      call dsptsi(n2ashx,pv_lu,pv_full_scratch)

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(lupri,*) ' final rho2 (full)'
      call wrtmatmn(pv_full_scratch,n2ashx,n2ashx,n2ashx,n2ashx,lupri)
#endif

      if(rhotype > 1)then
      call symmetrize_transition_density_matrix(pv_full_scratch,nasht,1)
      end if

#ifdef LUCI_DEBUG
      write(lupri,*) ' final rho2 (full) - symm attempt'
      call wrtmatmn(pv_full_scratch,n2ashx,n2ashx,n2ashx,n2ashx,lupri)
      do i = 1, nasht
        do j = 1, nasht
          do k = 1, nasht
            do l = 1, nasht
              write(lupri,'(a,i2,a1,i2,a1,i2,a1,i2,a,f16.10)') 
     &                       ' pv(',i,',',j,',',k,',',l,') = ',
     &                         pv_full_scratch(i,j,k,l)
            end do
          end do
        end do
      end do
#endif

      do_reorder = (nr_gas .gt. 1 .and. nr_sym .gt. 1)
      if (do_reorder) allocate(temp_mat(n2ashx))

      do l = 1, nasht
        l_ST = ireots(l)
        do k = 1, l
           k_ST = ireots(k)
           kl_ST_max = max(k_ST,l_ST)
           kl_ST_min = min(k_ST,l_ST)
           kl_ST_index = (kl_ST_max*(kl_ST_max-1))/2 + kl_ST_min
           if(do_reorder)then
!            reorder from type-symmetry ordering to a symmetry-type one
             call reormt(pv_full_scratch(1,1,k,l),temp_mat,
     &                   nasht,nasht,ireost,ireost)
!            distribute elements on lower triangle
             call dgetsp(nasht,temp_mat,pv_sirius(1,kl_ST_index))
           else
!            distribute elements on lower triangle
             call dgetsp(nasht,pv_full_scratch(1,1,k,l),
     &                   pv_sirius(1,kl_ST_index))
           end if
        end do
      end do


      if(do_reorder) deallocate(temp_mat)

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(lupri,*) ' final PV in SIRIUS format'
      call output(pv_sirius,1,nnashx,1,nnashx,nnashx,nnashx,-1,lupri)
#undef LUCI_DEBUG
#endif

      end
***********************************************************************

      subroutine set_ptg_multiplication_table(mult_ptg,nr_ptg_irreps)
!
!     purpose: set point-group multiplication table for 
!              dalton -> lucita integral transfer.
!
!     -----------------------------------------------------------------
      implicit none
#include "multd2h.inc"
      integer :: nr_ptg_irreps, i, j
      integer :: mult_ptg(nr_ptg_irreps,nr_ptg_irreps)
!     -----------------------------------------------------------------

      do j = 1, nr_ptg_irreps
         do i = 1, nr_ptg_irreps
            mult_ptg(i,j) = multd2h(i,j)
         end do
      end do

      end
***********************************************************************

      subroutine write_fcidump_file(AMATRIX,nrow,BMATRIX,corenergy)
!     -----------------------------------------------------------------
!     Written by Katharina Boguslawski and Pawel Tecmer  
!     ETH Zurich, March 2013
!
!     purpose: write one- and two-electron integrals for a given active
!              space to disk
!     
!     filename: FCIDUMP
!     -----------------------------------------------------------------
  
      use lucita_cfg

#include "implicit.h"
#include "maxorb.h"
#include "infinp.h"
#include "priunit.h"
#include "inforb.h"
#include "inftap.h"
#include "clunit.inc"
!#include "orbinp.inc"
      
      character(len=30) :: form1
      character(len=30) :: form2
      character(len=30) :: form3
      DIMENSION AMATRIX(*)
      DIMENSION BMATRIX(*)
      INTEGER NDUMMY, start
      INTEGER NTRIANGLE, norbtot, noccend, noccendi, noccendj
      integer, allocatable :: mult_ptg(:,:)
      integer nactiveorb(NSYM)
      integer :: offseti, offsetj, offsetk, offsetl
      real*8 :: corenergy
!     -----------------------------------------------------------------
      integer :: ISYM, KSYM, syml, JSYM, ijsym, klsym, ijklsym
      integer :: i, j, k, l
!     -----------------------------------------------------------------

      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX
             
!     allocate memory for multiplication table to write only
!     symmetry-allowed integrals
      allocate(mult_ptg(NSYM,NSYM))
      mult_ptg = 0
      call set_ptg_multiplication_table(mult_ptg,NSYM)

      LUFCID = -1

!     define printing format
      form1="(F17.12, 4X, I6, I6, I6, I6)"
      form2="(A11,I3,A7,I2,A5,I2,A1)"
      form3="(F17.8, 4X, I6, I6, I6, I6)"


!     FCIDUMP info
      write(lupri,*) 
      write(lupri,*)" FCIDUMP file details "
      write(lupri,*)"======================" 
      write(lupri,*) 
         

!     calculate number of active electrons
      norbtot = 0
      do KSYM=1,NSYM,1
         norbtot = norbtot + NOCC(KSYM) - NISH(KSYM)
         nactiveorb(KSYM) = NASH(KSYM)
      end do

!     Print header of FCIDUMP file
!     using MOLPRO format
      if (LUFCID .le. 0) then
      CALL GPOPEN(LUFCID,'FCIDUMP','new',' ',
     &                'FORMATTED',IDUMMY,.false.)
      end if
      write(LUFCID,form2) ' &FCI NORB=', norbtot , 
     & ',NELEC=', lucita_cfg_nr_active_e, ',MS2=',
     &  lucita_cfg_is_spin_multiplett-1, ','
      write(LUFCID,"(A)",advance='no') '  ORBSYM=' 
      do ISYM=1,NSYM,1
         if (nactiveorb(ISYM).ne.0) then
            do i=1,nactiveorb(ISYM),1
            write(LUFCID,"(I1,A1)",advance='no') ISYM,','
         end do
         end if
      end do


      write(lupri,form2) ' &FCI NORB=', norbtot , 
     & ',NELEC=', lucita_cfg_nr_active_e, ',MS2=',
     &  lucita_cfg_is_spin_multiplett-1, ','


      write(LUFCID,*) 
      write(LUFCID,"(A7,I1)") '  ISYM=',lucita_cfg_ptg_symmetry
      write(LUFCID,"(A5)") ' &END'

!     NOW:
!     Print two-electron integrals:
!     integrals are sorted as KL|IJ
!     order is different as in MOLPRO, where integrals are sorted as
!     IJ|KL
!     Integrals are stored in a set of triangular matrices, define index
!     which accesses matrix elements
      IELE = 0
!     Determine number of elements in one triangular matrix:
      NTRIANGLE = norbtot*(norbtot+1)/2

! KL IJ
!     define orbital offset for K
      offsetk = 0
!     sum over all irreducible representations
      do KSYM=1,NSYM,1
         if (nactiveorb(KSYM).ne.0) then
            do k=1,nactiveorb(KSYM),1 ! iterate through all elements
            offsetl = 0 !define orbital offset for L
!     restrict summation to prevent double counting
            do syml=1,KSYM,1
            if (nactiveorb(syml).ne.0) then
!     set upper summation bound for orbital index (prevent double
!     counting):
!     if not the same irrep l goes from 1 to number of orbitals
               if (KSYM.eq.syml) then
                  noccend = k
               else
                  noccend = nactiveorb(syml)
               end if
               do l=1,noccend,1
!     orbital offset for I
                  offseti = 0
!     restrict summation to prevent double counting for both irrep ISYM and
!     orbital indices i
                  do ISYM=1,KSYM,1
                     if (nactiveorb(ISYM).ne.0) then
                        ijksyml = 0
                        if (ISYM.eq.KSYM) then
                           noccendi = k
                        else
                           noccendi = nactiveorb(ISYM)
                        end if
                        do i=1,noccendi,1
!     orbital offset J                        
                        offsetj = 0
!     double counting problem: irrep of J must be smaller or equal to
!     irrep of I
      do JSYM=1,ISYM,1
!     Collect only integrals which are allowed by symmetry, all others
!     are neglected.
!     Two cases have to be distinguished: IJ|KL  and IK|JL
          if (nactiveorb(JSYM).ne.0) then
              if(ISYM.eq.JSYM.and.KSYM.eq.syml) then ! first case
                 ijsym   = mult_ptg(ISYM,JSYM)
                 ksyml   = mult_ptg(KSYM,syml)
                 ijksyml = mult_ptg(ijsym,ksyml)
              else ! second case
                 ijsym   = mult_ptg(ISYM,KSYM)
                 ksyml   = mult_ptg(JSYM,syml)
                 ijksyml = mult_ptg(ijsym,ksyml)
            
              end if

              if(ijksyml .eq. 1) then

!     Set matric index: first determine first index of block KL
              IELE = (offsetk+k)*(offsetk+k-1)/2*NTRIANGLE
     &              +(offsetl+l-1)*NTRIANGLE+1
!     Prevent double printing of symmetry redundant indices:
!     set upper summation index
!     if IJKL in same irrep, restrict j to at most i
              if (JSYM.eq.ISYM.and.KSYM.eq.syml) then !.and.ISYM.eq.KSYM) then
                  noccendj = i
!redundant    else if (JSYM.eq.ISYM.and.KSYM.eq.syml.and.ISYM.ne.KSYM) then
!                 noccendj = i
!     if LJ in irrep1 and IK in irrep2
              else if (syml.eq.JSYM.and.KSYM.eq.ISYM
     &                 .and.syml.ne.ISYM) then
!             second restriction to prevent double printing, J<=L in KL IJ
                  if (k.eq.i) then
                      noccendj = l
                  else ! otherwise all J are needed
                      noccendj = nactiveorb(JSYM)
                  end if
              else
                  noccendj = nactiveorb(JSYM)
              end if
              IELE = IELE+(offseti+i)*(offseti+i-1)/2+offsetj
              do j=1,noccendj,1
!                check for redundant integrals
                 if (JSYM.eq.ISYM.and.KSYM.eq.syml
     &               .and.ISYM.eq.KSYM) then
                     if (k.eq.i.and.l.eq.i.and.j.lt.i) then

                     else if (k.eq.i.and.j.lt.l) then
                      
                     else
                         write(LUFCID,form1) BMATRIX(IELE),
     &                        i+offseti, j+offsetj, k+offsetk, l+offsetl
                     end if
                 else
                     write(LUFCID,form1) BMATRIX(IELE),
     &                     i+offseti, j+offsetj, k+offsetk, l+offsetl
                 end if
                 IELE = IELE + 1
               end do ! do j
               offsetj = offsetj + nactiveorb(JSYM)

              else ! if not symmetry-allowed go to next irrep and update
!             orbital offset J
                   offsetj = offsetj + nactiveorb(JSYM)
              end if ! if ijksyml
          end if ! if nocc(jsym)
       end do ! do jsym

                        end do ! do i
                        offseti = offseti + nactiveorb(ISYM) !update orbital
!                       offset I
                     end if ! if nocc(isym)
                  end do ! do isym
               end do ! do l
               offsetl = offsetl + nactiveorb(syml)
               end if ! if nocc(syml)
            end do ! do syml
            end do ! do k
         end if ! if nocc(ksym)
            offsetk = offsetk + nactiveorb(KSYM)
      end do ! do ksym


!     Print one-electron integrals:
!     Reset matrix index
      IELE = 1
!     number of dummy indices to be ignored in one-electron
!     integrals because of symmetry ISYM < actual ISYM
      NDUMMY = 0
      
      do ISYM=1,NSYM,1
         if (nactiveorb(ISYM).ne.0) then
            IELE = IELE + NDUMMY
            do i=1,nactiveorb(ISYM),1 ! only loop through same irrep
               do j=1,i,1
                  write(LUFCID,form1) AMATRIX(IELE), i+NDUMMY, j+NDUMMY,
     &                                0, 0
                  IELE = IELE + 1
               end do
!              add orbital offset if irrep changes
               if (i.lt.nactiveorb(ISYM)) then
                   IELE = IELE + NDUMMY
               end if
            end do
!           update orbital offset
            NDUMMY = NDUMMY + nactiveorb(ISYM)
         end if
      end do

      write(LUFCID,form3) POTNUC + corenergy, 0,0 ,0,0
 
      call GPCLOSE(LUFCID,'KEEP')
 
      deallocate(mult_ptg)

      return     
 
      end 
