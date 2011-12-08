!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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
#include <single.h>
C
      SUBROUTINE C3FOCK(IBEQC,FDTI,VECB,VECC,
     *                 ZYMB,ZYMC,KZYVB,KZYVC,
     *                 ISYMA,ISYMB,ISYMC,
     *                 CMO,FC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C To calculate Fock matrix (FDTI) with doubly one-index transformed
C integrals to be used in direct cubic RPA.
C March 1997: tsaue Just one call to SIRFCK, big speedup !
C
#include <implicit.h>
C
      DIMENSION FDTI(N2ORBX)
      DIMENSION VECB(*),VECC(*)
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT)
      DIMENSION CMO(*),FC(*),WRK(*),ISYMDM(3)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
#include <infvar.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <wrkrsp.h>
#include <infrsp.h>
C
      CALL QENTER('C3FOCK')
      IF (8*N2BASX .GT. LWRK) THEN
         WRITE (LUPRI,'(A)') ' Due to lack of available memory, I '//
     &     'use a slower C3FOCK routine'
         CALL C3FCKM(IBEQC,FDTI,VECB,VECC,
     *               ZYMB,ZYMC,KZYVB,KZYVC,
     *               ISYMA,ISYMB,ISYMC,
     *               CMO,FC,MJWOP,WRK,LWRK)
      ELSE
         CALL C3FCKO(IBEQC,FDTI,VECB,VECC,
     *               ZYMB,ZYMC,KZYVB,KZYVC,
     *               ISYMA,ISYMB,ISYMC,
     *               CMO,FC,MJWOP,WRK(1),LWRK)
      END IF
      CALL QEXIT('C3FOCK')
      RETURN
      END
C
      SUBROUTINE C3FCKO(IBEQC,FDTI,VECB,VECC,
     *                 ZYMB,ZYMC,KZYVB,KZYVC,
     *                 ISYMA,ISYMB,ISYMC,
     *                 CMO,FC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C To calculate Fock matrix (FDTI) with doubly one-index transformed
C integrals to be used in direct cubic RPA.
C March 1997: tsaue Just one call to SIRFCK, big speedup !
C
#include <implicit.h>
#include <dummy.h>
#include <thrzer.h>
C
      DIMENSION FDTI(N2ORBX)
      DIMENSION VECB(*),VECC(*)
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT)
      DIMENSION CMO(*),FC(*),WRK(*),ISYMDM(3),IFCTYP(3)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
#include <infvar.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0 )
C
C Allocate work space for Fock matrices
C Corresponding Fock matrices in equations
C
C   F-1 <--> F12 + F21
C   F-2 <--> F1
C   F-3 <--> F2
C
      CALL QENTER('C3FCKO')
      KF1AO = 1
      KF2AO = KF1AO + N2BASX
      KF3AO = KF2AO + N2BASX
      KD1AO = KF3AO + N2BASX
      KD2AO = KD1AO + N2BASX
      KD3AO = KD2AO + N2BASX
      KFREE = KD3AO + N2BASX
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('C3FCKO',LFREE,LWRK)
C
C     Unpack the response vectors
C
      CALL GTZYMT(1,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL GTZYMT(1,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
C
C
C Construct density matrices in AO basis
C
C
      CALL DZERO(WRK(KF1AO),6*N2BASX)
C
C Alternative construction need to assign KTMP, KBCZYM, ISYMBC
C Further scale D1AO with factor 2
C
      ISYMDM(1) = MULD2H(ISYMB,ISYMC)
      ISYMDM(2) = ISYMB
      ISYMDM(3) = ISYMC
C
      KBCZYM = KF1AO
      KTMP   = KF2AO
C
      CALL ZYTRA1(ISYMB,ZYMB,CMO,WRK(KTMP),WRK(KD2AO))
      CALL ZYTRA1(ISYMC,ZYMC,CMO,WRK(KTMP),WRK(KD3AO))
C
      CALL DZERO(WRK(KTMP),N2BASX)
      CALL ZYMUL2(ISYMB,ISYMC,ZYMB,ZYMC,WRK(KBCZYM))
      CALL ZYMUL2(ISYMC,ISYMB,ZYMC,ZYMB,WRK(KBCZYM))
      CALL ZYTRA2(ISYMDM(1),WRK(KBCZYM),CMO,WRK(KTMP),WRK(KD1AO))
C
C
C Transpose and scale by 2 to agree with convention used by SIRFCK
C Further scale by 2 due to permutations. D2AO and D3AO scaled in one call.
C
      CALL DSCAL(3*N2BASX,D4,WRK(KD1AO),1)
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'D12 in C3FOCK'
         CALL OUTPUT(WRK(KD1AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D1 in C3FOCK'
         CALL OUTPUT(WRK(KD2AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D2 in C3FOCK'
         CALL OUTPUT(WRK(KD3AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C Construct Fock matrices in AO basis
C
      NFMATX = 1
C
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
      JD1AO = KD1AO
      DO I = 1,3
         IX = 10 * MATSYM(NBAST,NBAST,WRK(JD1AO),THRZER)
C        INTEGER FUNCTION MATSYM(N,NDIM,AMAT,THRZER)
C
         IF (IX .EQ. 30) THEN
C           zero density matrix, do nothing !
            IFCTYP(I) = 0
CCCC     ELSE IF ("triplet Fock matrix") THEN
C           only exchange ! hjaaj 990505: triplet not implemented!
CCCC        IFCTYP(I) = IX + 2
         ELSE IF (IX .EQ. 20) THEN
C           Only exchange if antisymmetric density matrix
            IFCTYP(I) = IX + 2
         ELSE
C           Coulomb+exchange
            IFCTYP(I) = IX + 3
         END IF
         JD1AO = JD1AO + N2BASX
      END DO
C
      CALL DZERO(WRK(KF1AO),3*N2BASX)
      IF (IBEQC.EQ.1) THEN
        NFMATX = 2
        CALL SIRFCK(WRK(KF1AO),WRK(KD1AO),NFMATX,ISYMDM,IFCTYP,
     *            DIRFCK,WRK(KFREE),LFREE)
        CALL DCOPY(N2BASX,WRK(KF2AO),1,WRK(KF3AO),1)
      ELSE
        NFMATX = 3
        CALL SIRFCK(WRK(KF1AO),WRK(KD1AO),NFMATX,ISYMDM,IFCTYP,
     *            DIRFCK,WRK(KFREE),LFREE)
      ENDIF
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F12 in AO basis in C3FOCK'
         CALL OUTPUT(WRK(KF1AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F1 in AO basis in C3FOCK'
         CALL OUTPUT(WRK(KF2AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F2 in AO basis in C3FOCK'
         CALL OUTPUT(WRK(KF3AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C Transform Fock matrices to MO basis
C Reuse memory used for density matrices
C
      KF2MO = KD1AO
      KF3MO = KF2MO + N2ORBX
C
      CALL FAOMO(ISYMDM(1),WRK(KF1AO),FDTI,CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(2),WRK(KF2AO),WRK(KF2MO),CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(3),WRK(KF3AO),WRK(KF3MO),CMO,WRK(KFREE),LFREE)
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F12 in MO basis in C3FOCK'
         CALL OUTPUT(FDTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F1 in MO basis in C3FOCK'
         CALL OUTPUT(WRK(KF2MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F2 in MO basis in C3FOCK'
         CALL OUTPUT(WRK(KF3MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Evaluate FDTI see formula
C
C 2*F2 + [k2,Fc] --> F-3
C 2*F1 + [k1,Fc] --> F-2
C
      KSYMOP = ISYMC
      CALL FCKOIN(1,FC,DUMMY,ZYMC,WRK(KF3MO),DUMMY)
      KSYMOP = ISYMB
      CALL FCKOIN(1,FC,DUMMY,ZYMB,WRK(KF2MO),DUMMY)
C
C F12  + [k2,F-3] --> FDTI
C FDTI + [k1,F-2] --> FDTI
C
      CALL OITH1(ISYMC,ZYMC,WRK(KF2MO),FDTI,ISYMDM(2))
      CALL OITH1(ISYMB,ZYMB,WRK(KF3MO),FDTI,ISYMDM(3))
C
C Scale with 1/2 ( from definition of E3 )
C
      CALL DSCAL(NORBT*NORBT,1/D2,FDTI,1)
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(//A)')'Final result in C3FCKO'
         CALL OUTPUT(FDTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL QEXIT('C3FCKO')
      RETURN
      END
C
      SUBROUTINE C3FCKM(IBEQC,FDTI,VECB,VECC,
     *                 ZYMB,ZYMC,KZYVB,KZYVC,
     *                 ISYMA,ISYMB,ISYMC,
     *                 CMO,FC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C To calculate Fock matrix (FDTI) with doubly one-index transformed
C integrals to be used in direct cubic RPA.
C
#include <implicit.h>
#include <dummy.h>
#include <thrzer.h>
C
      DIMENSION FDTI(N2ORBX)
      DIMENSION VECB(*),VECC(*)
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT)
      DIMENSION CMO(*),FC(*),WRK(*)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
#include <infvar.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0 )
C
      NFMATX = 1
C
C Allocate work space for Fock matrices
C Corresponding Fock matrices in equations
C
C   F-1 <--> F12 + F21
C   F-2 <--> F1
C   F-3 <--> F2
C
      CALL QENTER('C3FCKM')
      KFAO  = 1
      KDAO  = KFAO  + N2BASX
      KFRE1 = KDAO  + N2BASX
      KFREE = KFRE1 + N2BASX
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('C3FCKM',LFREE,LWRK)
C
      KFREE = KFRE1
      LFREE = LWRK - KFREE
C
C     Unpack the response vectors
C
      CALL GTZYMT(1,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL GTZYMT(1,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
C
C
C Construct density matrices in AO basis
C
C
      CALL DZERO(WRK(KFAO),3*N2BASX)
      CALL DZERO(FDTI,N2ORBX)
      CALL CDENS2(ISYMB,ISYMC,CMO,ZYMB,ZYMC,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),FDTI)
      CALL CDENS2(ISYMC,ISYMB,CMO,ZYMC,ZYMB,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),FDTI)
      CALL DSCAL(N2BASX,D2,WRK(KDAO),1)
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'D12 in C3FOCK'
         CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      ISYMDM = MULD2H(ISYMB,ISYMC)
      CALL DZERO(WRK(KFAO),N2BASX)
Chj
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
      IX = 10 * MATSYM(NBAST,NBAST,WRK(KDAO),THRZER)
C
      IF (IX .EQ. 30) THEN
C        zero density matrix, do nothing !
         IFCTYP = 0
      ELSE IF (IX .EQ. 20) THEN
C        Only exchange if antisymmetric density matrix
         IFCTYP = IX + 2
      ELSE
C        Coulomb+exchange
         IFCTYP = IX + 3
      END IF
      IF (IFCTYP .NE. 0) THEN
         CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *               DIRFCK,WRK(KFREE),LFREE)
      END IF
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F12 in AO basis in C3FOCK'
         CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL FAOMO(ISYMDM,WRK(KFAO),FDTI,CMO,WRK(KFREE),LFREE)
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F12 in MO basis in C3FOCK'
         CALL OUTPUT(FDTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Second contribution
C
      CALL DZERO(WRK(KFAO),2*N2BASX)
      CALL CDENS1(ISYMB,CMO,ZYMB,WRK(KDAO),WRK(KFREE),LFREE)
      CALL DSCAL(N2BASX,D4,WRK(KDAO),1)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'D1 in C3FOCK'
         CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      ISYMDM = ISYMB
Chj
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
      IX = 10 * MATSYM(NBAST,NBAST,WRK(KDAO),THRZER)
C
      IF (IX .EQ. 30) THEN
C        zero density matrix, do nothing !
         IFCTYP = 0
      ELSE IF (IX .EQ. 20) THEN
C        Only exchange if antisymmetric density matrix
         IFCTYP = IX + 2
      ELSE
C        Coulomb+exchange
         IFCTYP = IX + 3
      END IF
      IF (IFCTYP .NE. 0) THEN
         CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *               DIRFCK,WRK(KFREE),LFREE)
      END IF
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'F1 in AO basis in C3FOCK'
         CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'F1 in MO basis in C3FOCK'
         CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      KSYMOP = ISYMB
      CALL FCKOIN(1,FC,DUMMY,ZYMB,WRK(KDAO),DUMMY)
      IF (IBEQC .EQ. 1) THEN
C
C     If B and C vector identical, all that remains to do is to scale
C     this term by a factor of two.
C
         CALL DSCAL(N2ORBX,D2,WRK(KDAO),1)
      ELSE
         CALL OITH1(ISYMC,ZYMC,WRK(KDAO),FDTI,ISYMDM)
C
C     Third contribution
C
         CALL DZERO(WRK(KFAO),2*N2BASX)
         CALL CDENS1(ISYMC,CMO,ZYMC,WRK(KDAO),WRK(KFREE),LFREE)
         CALL DSCAL(N2BASX,D4,WRK(KDAO),1)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'D2 in C3FOCK'
            CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C
         ISYMDM = ISYMC
Chj
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
         IX = 10 * MATSYM(NBAST,NBAST,WRK(KDAO),THRZER)
C
         IF (IX .EQ. 30) THEN
C           zero density matrix, do nothing !
            IFCTYP = 0
         ELSE IF (IX .EQ. 20) THEN
C           Only exchange if antisymmetric density matrix
            IFCTYP = IX + 2
         ELSE
C           Coulomb+exchange
            IFCTYP = IX + 3
         END IF
         IF (IFCTYP .NE. 0) THEN
            CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *                  DIRFCK,WRK(KFREE),LFREE)
         END IF
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F2 in AO basis in C3FOCK'
            CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C
         CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F2 in MO basis in C3FOCK'
            CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
C
         KSYMOP = ISYMC
         CALL FCKOIN(1,FC,DUMMY,ZYMC,WRK(KDAO),DUMMY)
         CALL OITH1(ISYMB,ZYMB,WRK(KDAO),FDTI,ISYMDM)
      END IF
C
C     All contributions to FDTI now finished, we 
C     scale with 1/2 ( from definition of E3 )
C
      CALL DSCAL(NORBT*NORBT,1/D2,FDTI,1)
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(//A)')'Final result in C3FCKM'
         CALL OUTPUT(FDTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL QEXIT('C3FCKM')
      RETURN
      END
C
      SUBROUTINE C4FOCK(IBCDEQ,FTTI,VECB,VECC,VECD,
     *                 ZYMB,ZYMC,ZYMD,KZYVB,KZYVC,KZYVD,
     *                 ISYMA,ISYMB,ISYMC,ISYMD,
     *                 CMO,FC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C To calculate Fock matrix (FTTI) with triply one-index transformed
C integrals to be used in direct cubic RPA.
C March 1997: tsaue Just one call to SIRFCK, big speedup !
C July 1997:  kruud. Interface routine; further path depends on available
C                    memory
C
#include <implicit.h>
C
      DIMENSION FTTI(N2ORBX)
      DIMENSION VECB(*),VECC(*),VECD(*)
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT),ZYMD(NORBT,NORBT)
      DIMENSION CMO(*),FC(*),WRK(*)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
#include <infvar.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <wrkrsp.h>
#include <infrsp.h>
C
      IF (14*N2BASX .GT. LWRK) THEN
         WRITE (LUPRI,'(A)') ' Due to potential lack of available '//
     &     'memory, I use a slower C4FOCK routine'
         CALL C4FCKM(IBCDEQ,FTTI,VECB,VECC,VECD,
     *               ZYMB,ZYMC,ZYMD,KZYVB,KZYVC,KZYVD,
     *               ISYMA,ISYMB,ISYMC,ISYMD,
     *               CMO,FC,MJWOP,WRK,LWRK)
      ELSE
         CALL C4FCKO(IBCDEQ,FTTI,VECB,VECC,VECD,
     *               ZYMB,ZYMC,ZYMD,KZYVB,KZYVC,KZYVD,
     *               ISYMA,ISYMB,ISYMC,ISYMD,
     *               CMO,FC,MJWOP,WRK,LWRK)
      END IF
C
      RETURN
      END
C
      SUBROUTINE C4FCKO(IBCDEQ,FTTI,VECB,VECC,VECD,
     *                  ZYMB,ZYMC,ZYMD,KZYVB,KZYVC,KZYVD,
     *                  ISYMA,ISYMB,ISYMC,ISYMD,
     *                  CMO,FC,MJWOP,WRK,LWRK)
#include <implicit.h>
#include <dummy.h>
C
      DIMENSION FTTI(N2ORBX)
      DIMENSION VECB(*),VECC(*),VECD(*)
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT),ZYMD(NORBT,NORBT)
      DIMENSION CMO(*),FC(*),WRK(*),ISYMDM(7),IFCTYP(7)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
#include <infvar.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, D6 = 6.0D0 , D8 = 8.0D0 ,
     &           D12 = 12.0D0)
C Allocate work space for Fock matrices
C Corresponding Fock matrices in equations
C
C   F-1 <--> F123 + permutations
C   F-2 <--> F12 + F21
C   F-3 <--> F13 + F31
C   F-4 <--> F23 + F32
C   F-5 <--> F1
C   F-6 <--> F2
C   F-7 <--> F3
C
      KF1AO = 1
      KF2AO = KF1AO + N2BASX
      KF3AO = KF2AO + N2BASX
      KF4AO = KF3AO + N2BASX
      KF5AO = KF4AO + N2BASX
      KF6AO = KF5AO + N2BASX
      KF7AO = KF6AO + N2BASX
      KD1AO = KF7AO + N2BASX
      KD2AO = KD1AO + N2BASX
      KD3AO = KD2AO + N2BASX
      KD4AO = KD3AO + N2BASX
      KD5AO = KD4AO + N2BASX
      KD6AO = KD5AO + N2BASX
      KD7AO = KD6AO + N2BASX
      KFREE = KD7AO + N2BASX
      LFREE = LWRK  - KFREE
C
      IF (IPRRSP.GT.7) WRITE(LUPRI,'(//A,2(/A,I10))')
     *' Memory allocation in C4FCKO',
     *' Available workspace: ', LWRK,
     *' Allocated          : ', KFREE
C
C     Unpack the response vectors
C
      CALL GTZYMT(1,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL GTZYMT(1,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
      CALL GTZYMT(1,VECD,KZYVD,ISYMD,ZYMD,MJWOP)
C
      CALL DZERO(WRK(KF1AO),14*N2BASX)
C
C Alternativ construction need to assign KTMP, K..ZYM, K...ZYM
C
      KBCZYM  = KF1AO
      KBCDZYM = KF2AO
      KTMP    = KF3AO
      KBDZYM  = KBCZYM
      KCDZYM  = KBCZYM
C
      ISYMDM(1) = ISYMA
      ISYMDM(2) = MULD2H(ISYMB,ISYMC)
      ISYMDM(3) = MULD2H(ISYMB,ISYMD)
      ISYMDM(4) = MULD2H(ISYMC,ISYMD)
      ISYMDM(5) = ISYMB
      ISYMDM(6) = ISYMC
      ISYMDM(7) = ISYMD
      ISYMBCD   = MULD2H(ISYMDM(2),ISYMD)
C
      CALL ZYTRA1(ISYMB,ZYMB,CMO,WRK(KTMP),WRK(KD5AO))
      CALL ZYTRA1(ISYMC,ZYMC,CMO,WRK(KTMP),WRK(KD6AO))
      CALL ZYTRA1(ISYMD,ZYMD,CMO,WRK(KTMP),WRK(KD7AO))
C
      CALL DZERO(WRK(KBCZYM),3*N2ORBX)
C
      CALL ZYMUL2(ISYMB,ISYMD,ZYMB,ZYMD,WRK(KBDZYM))
      CALL ZYMUL2(ISYMD,ISYMB,ZYMD,ZYMB,WRK(KBDZYM))
      CALL ZYMUL3(ISYMDM(3),ISYMC,WRK(KBDZYM),ZYMC,WRK(KBCDZYM))
      CALL ZYTRA2(ISYMDM(3),WRK(KBDZYM),CMO,WRK(KTMP),WRK(KD3AO))
C
      CALL DZERO(WRK(KBDZYM),N2ORBX)
      CALL DZERO(WRK(KTMP),N2BASX)
      CALL ZYMUL2(ISYMB,ISYMC,ZYMB,ZYMC,WRK(KBCZYM))
      CALL ZYMUL2(ISYMC,ISYMB,ZYMC,ZYMB,WRK(KBCZYM))
      CALL ZYMUL3(ISYMDM(2),ISYMD,WRK(KBCZYM),ZYMD,WRK(KBCDZYM))
      CALL ZYTRA2(ISYMDM(2),WRK(KBCZYM),CMO,WRK(KTMP),WRK(KD2AO))
C
      CALL DZERO(WRK(KCDZYM),N2ORBX)
      CALL DZERO(WRK(KCDZYM),N2ORBX)
      CALL ZYMUL2(ISYMC,ISYMD,ZYMC,ZYMD,WRK(KCDZYM))
      CALL ZYMUL2(ISYMD,ISYMC,ZYMD,ZYMC,WRK(KCDZYM))
      CALL ZYMUL3(ISYMDM(4),ISYMB,WRK(KCDZYM),ZYMB,WRK(KBCDZYM))
      CALL ZYTRA2(ISYMDM(4),WRK(KCDZYM),CMO,WRK(KTMP),WRK(KD4AO))
C
      CALL DZERO(WRK(KTMP),N2BASX)
      CALL ZYTRA1(ISYMBCD,WRK(KBCDZYM),CMO,WRK(KTMP),WRK(KD1AO))
C
C Transpose and scale by 2 to agree with convention used by SIRFCK
C Further scale by 3 due to permutations. D2AO -> D7AO scaled in one call
C
      CALL DSCAL(N2BASX,D8,WRK(KD1AO),1)
      CALL DSCAL(3*N2BASX,D12,WRK(KD2AO),1)
      CALL DSCAL(3*N2BASX,D6,WRK(KD5AO),1)
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'D123 in C4FOCK'
         CALL OUTPUT(WRK(KD1AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D12 in C4FOCK'
         CALL OUTPUT(WRK(KD2AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D13 in C4FOCK'
         CALL OUTPUT(WRK(KD3AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D23 in C4FOCK'
         CALL OUTPUT(WRK(KD4AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D1 in C4FOCK'
         CALL OUTPUT(WRK(KD5AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D2 in C4FOCK'
         CALL OUTPUT(WRK(KD6AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'D3 in C4FOCK'
         CALL OUTPUT(WRK(KD7AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C Construct Fock matrices in AO basis
C
      DO I = 1,7
        IFCTYP(I) = 3
      ENDDO
      NFMATX = 7
      IF (IBCDEQ.EQ.24) THEN
         NFMATX    = 5
         IFCTYP(3) = 0
         IFCTYP(4) = 0
         IFCTYP(6) = 0
         IFCTYP(7) = 0
      ELSE IF (IBCDEQ.EQ.2) THEN
         IFCTYP(4) = 0
         IFCTYP(6) = 0
      ELSE IF (IBCDEQ.EQ.3) THEN
         NFMATX    = 6
         IFCTYP(4) = 0
         IFCTYP(7) = 0
      ELSE IF (IBCDEQ.EQ.4) THEN
         NFMATX    = 6
         IFCTYP(3) = 0
         IFCTYP(7) = 0
      ENDIF
C
      CALL DZERO(WRK(KF1AO),7*N2BASX)
      CALL SIRFCK(WRK(KF1AO),WRK(KD1AO),NFMATX,ISYMDM,IFCTYP,
     *            DIRFCK,WRK(KFREE),LFREE)
C
      IF (IBCDEQ.EQ.24) THEN
         CALL DCOPY(N2BASX,WRK(KF2AO),1,WRK(KF3AO),1)
         CALL DCOPY(N2BASX,WRK(KF2AO),1,WRK(KF4AO),1)
         CALL DCOPY(N2BASX,WRK(KF5AO),1,WRK(KF6AO),1)
         CALL DCOPY(N2BASX,WRK(KF5AO),1,WRK(KF7AO),1)
      ELSE IF (IBCDEQ.EQ.2) THEN
         CALL DCOPY(N2BASX,WRK(KF3AO),1,WRK(KF4AO),1)
         CALL DCOPY(N2BASX,WRK(KF5AO),1,WRK(KF6AO),1)
      ELSE IF (IBCDEQ.EQ.3) THEN
         CALL DCOPY(N2BASX,WRK(KF2AO),1,WRK(KF4AO),1)
         CALL DCOPY(N2BASX,WRK(KF5AO),1,WRK(KF7AO),1)
      ELSE IF (IBCDEQ.EQ.4) THEN
         CALL DCOPY(N2BASX,WRK(KF2AO),1,WRK(KF3AO),1)
         CALL DCOPY(N2BASX,WRK(KF6AO),1,WRK(KF7AO),1)
      END IF
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F123 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF1AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F12 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF2AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F13 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF3AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F23 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF4AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F1 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF5AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F2 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF6AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F3 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KF7AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C Transform Fock matrices to MO basis
C Reuse memory used for density matrices
C
      KF2MO = KD1AO
      KF3MO = KF2MO + N2ORBX
      KF4MO = KF3MO + N2ORBX
      KF5MO = KF4MO + N2ORBX
      KF6MO = KF5MO + N2ORBX
      KF7MO = KF6MO + N2ORBX
C
      CALL FAOMO(ISYMDM(1),WRK(KF1AO),FTTI,CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(2),WRK(KF2AO),WRK(KF2MO),CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(3),WRK(KF3AO),WRK(KF3MO),CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(4),WRK(KF4AO),WRK(KF4MO),CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(5),WRK(KF5AO),WRK(KF5MO),CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(6),WRK(KF6AO),WRK(KF6MO),CMO,WRK(KFREE),LFREE)
      CALL FAOMO(ISYMDM(7),WRK(KF7AO),WRK(KF7MO),CMO,WRK(KFREE),LFREE)
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F123 in MO basis in C4FOCK'
         CALL OUTPUT(FTTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F12 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KF2MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F13 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KF3MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F23 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KF4MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F1 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KF5MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F2 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KF6MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')'F3 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KF7MO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Evaluate FTTI see formula
C
C 3*F3 + [k3,Fc] --> F-7
C 3*F2 + [k2,Fc] --> F-6
C 3*F1 + [k1,Fc] --> F-5
C
      KSYMOP = ISYMD
      CALL FCKOIN(1,FC,DUMMY,ZYMD,WRK(KF7MO),DUMMY)
      KSYMOP = ISYMC
      CALL FCKOIN(1,FC,DUMMY,ZYMC,WRK(KF6MO),DUMMY)
      KSYMOP = ISYMB
      CALL FCKOIN(1,FC,DUMMY,ZYMB,WRK(KF5MO),DUMMY)
C
C F23 + [k2,F-7] --> F-4
C F-4 + [k3,F-6] --> F-4
C
      CALL OITH1(ISYMC,ZYMC,WRK(KF7MO),WRK(KF4MO),ISYMDM(7))
      CALL OITH1(ISYMD,ZYMD,WRK(KF6MO),WRK(KF4MO),ISYMDM(6))
C
C F13 + [k1,F-7] --> F-3
C F-3 + [k3,F-5] --> F-3
C
      CALL OITH1(ISYMB,ZYMB,WRK(KF7MO),WRK(KF3MO),ISYMDM(7))
      CALL OITH1(ISYMD,ZYMD,WRK(KF5MO),WRK(KF3MO),ISYMDM(5))
C
C F12 + [k2,F-5] --> F-2
C F-2 + [k1,F-6] --> F-2
C
      CALL OITH1(ISYMC,ZYMC,WRK(KF5MO),WRK(KF2MO),ISYMDM(5))
      CALL OITH1(ISYMB,ZYMB,WRK(KF6MO),WRK(KF2MO),ISYMDM(6))
C
C F123 + [k1,F-4] --> FTTI
C FTTI + [k2,F-3] --> FTTI
C FTTI + [k3,F-2] --> FTTI
C
      CALL OITH1(ISYMB,ZYMB,WRK(KF4MO),FTTI,ISYMDM(4))
      CALL OITH1(ISYMC,ZYMC,WRK(KF3MO),FTTI,ISYMDM(3))
      CALL OITH1(ISYMD,ZYMD,WRK(KF2MO),FTTI,ISYMDM(2))
C
C Scale with 1/6 ( from definition of E4, without minus sign )
C
      CALL DSCAL(NORBT*NORBT,1/D6,FTTI,1)
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(//A)')'Final result in C4FOCK'
         CALL OUTPUT(FTTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
C
      SUBROUTINE C4FCKM(IBCDEQ,FTTI,VECB,VECC,VECD,
     *                 ZYMB,ZYMC,ZYMD,KZYVB,KZYVC,KZYVD,
     *                 ISYMA,ISYMB,ISYMC,ISYMD,
     *                 CMO,FC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C To calculate Fock matrix (FTTI) with triply one-index transformed
C integrals to be used in direct cubic RPA.
C July  1997: kruud: Memory efficient version of C4FOCK
C
#include <implicit.h>
#include <dummy.h>
C
      DIMENSION FTTI(N2ORBX)
      DIMENSION VECB(*),VECC(*),VECD(*)
      DIMENSION ZYMB(NORBT,NORBT),ZYMC(NORBT,NORBT),ZYMD(NORBT,NORBT)
      DIMENSION CMO(*),FC(*),WRK(*)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
#include <infvar.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, D6 = 6.0D0, D3 = 3.0D0)
C
C Allocate work space for Fock matrices
C Corresponding Fock matrices in equations
C
C   F-1 <--> F123 + permutations
C   F-2 <--> F12 + F21
C   F-3 <--> F13 + F31
C   F-4 <--> F23 + F32
C   F-5 <--> F1
C   F-6 <--> F2
C   F-7 <--> F3
C
      IFCTYP = 3
      NFMATX = 1
C
      KFAO  = 1
      KDAO  = KFAO  + N2BASX
      KFRE1 = KDAO  + N2BASX
      KFRE2 = KFRE1 + N2BASX
      KFREE = KFRE2 + N2BASX
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('C4FCKM',LFREE,LWRK)
C
      IF (IPRRSP.GT.7) WRITE(LUPRI,'(//A,2(/A,I10))')
     *' Memory allocation in C4FCKM',
     *' Available workspace: ', LWRK,
     *' Allocated          : ', KFREE
C
      KFREE = KFRE1
      LFREE = LWRK - KFREE
C
C     Unpack the response vectors
C
      CALL GTZYMT(1,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL GTZYMT(1,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
      CALL GTZYMT(1,VECD,KZYVD,ISYMD,ZYMD,MJWOP)
C
C     F123 contribution
C
      CALL DZERO(WRK(KFAO),4*N2BASX)
      CALL CDENS3(ISYMB,ISYMC,ISYMD,CMO,ZYMB,ZYMC,ZYMD,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL CDENS3(ISYMB,ISYMD,ISYMC,CMO,ZYMB,ZYMD,ZYMC,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL CDENS3(ISYMC,ISYMB,ISYMD,CMO,ZYMC,ZYMB,ZYMD,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL CDENS3(ISYMC,ISYMD,ISYMB,CMO,ZYMC,ZYMD,ZYMB,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL CDENS3(ISYMD,ISYMB,ISYMC,CMO,ZYMD,ZYMB,ZYMC,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL CDENS3(ISYMD,ISYMC,ISYMB,CMO,ZYMD,ZYMC,ZYMB,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
C
      CALL DSCAL(N2BASX,D2,WRK(KDAO),1)
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'D123 in C4FOCK'
         CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      ISYMDM = ISYMA
      CALL DZERO(WRK(KFAO),N2BASX)
      CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *            DIRFCK,WRK(KFREE),LFREE)
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F123 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL FAOMO(ISYMDM,WRK(KFAO),FTTI,CMO,WRK(KFREE),LFREE)
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(//A)')'F123 in MO basis in C4FOCK'
         CALL OUTPUT(FTTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     F12 matrix
C
      CALL DZERO(WRK(KFAO),4*N2BASX)
      CALL CDENS2(ISYMB,ISYMC,CMO,ZYMB,ZYMC,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL CDENS2(ISYMC,ISYMB,CMO,ZYMC,ZYMB,
     *            WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
      CALL DSCAL(N2BASX,D6,WRK(KDAO),1)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'D12 in C4FOCK'
         CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      ISYMDM = MULD2H(ISYMB,ISYMC)
      CALL DZERO(WRK(KFAO),N2BASX)
      CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *            DIRFCK,WRK(KFREE),LFREE)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'F12 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'F12 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      IF (IBCDEQ .EQ. 24) THEN
         CALL DSCAL(N2BASX,D3,WRK(KDAO),1)
      END IF
      CALL OITH1(ISYMD,ZYMD,WRK(KDAO),FTTI,ISYMDM)
C
C     F13 matrix
C
      IF (IBCDEQ .NE. 24) THEN
         CALL DZERO(WRK(KFAO),4*N2BASX)
         CALL CDENS2(ISYMB,ISYMD,CMO,ZYMB,ZYMD,
     *               WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
         CALL CDENS2(ISYMD,ISYMB,CMO,ZYMD,ZYMB,
     *               WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
         CALL DSCAL(N2BASX,D6,WRK(KDAO),1)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'D13 in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         ISYMDM = MULD2H(ISYMB,ISYMD)
         CALL DZERO(WRK(KFAO),N2BASX)
         CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *        DIRFCK,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F13 in AO basis in C4FOCK'
            CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C
         CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F13 in MO basis in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
         CALL OITH1(ISYMC,ZYMC,WRK(KDAO),FTTI,ISYMDM)
C     
C     F23 contribution
C
         CALL DZERO(WRK(KFAO),4*N2BASX)
         CALL CDENS2(ISYMC,ISYMD,CMO,ZYMC,ZYMD,
     *        WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
         CALL CDENS2(ISYMD,ISYMC,CMO,ZYMD,ZYMC,
     *        WRK(KDAO),WRK(KFAO),WRK(KFRE1),WRK(KFRE2))
         CALL DSCAL(N2BASX,D6,WRK(KDAO),1)
         IF (IPRRSP .GT. 20) THEN
            WRITE(LUPRI,'(//A)')'D23 in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         ISYMDM = MULD2H(ISYMC,ISYMD)
         CALL DZERO(WRK(KFAO),N2BASX)
         CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *        DIRFCK,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F23 in AO basis in C4FOCK'
            CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F23 in MO basis in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
         CALL OITH1(ISYMB,ZYMB,WRK(KDAO),FTTI,ISYMDM)
      END IF
C
C     F1 contribution
C
      CALL DZERO(WRK(KFAO),2*N2BASX)
      CALL CDENS1(ISYMB,CMO,ZYMB,WRK(KDAO),
     *            WRK(KFREE),LFREE)
      CALL DSCAL(N2BASX,D6,WRK(KDAO),1)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'D1 in C4FOCK'
         CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      ISYMDM = ISYMB
      CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *            DIRFCK,WRK(KFREE),LFREE)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'F1 in AO basis in C4FOCK'
         CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
      IF (IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)')'F1 in MO basis in C4FOCK'
         CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      KSYMOP = ISYMB
      CALL FCKOIN(1,FC,DUMMY,ZYMB,WRK(KDAO),DUMMY)
      CALL DZERO(WRK(KFAO),N2BASX)
      CALL OITH1(ISYMD,ZYMD,WRK(KDAO),WRK(KFAO),ISYMDM)
      IF (IBCDEQ .EQ. 24) THEN
         CALL DSCAL(N2ORBX,D3,WRK(KFAO),1)
      END IF
      CALL OITH1(ISYMC,ZYMC,WRK(KFAO),FTTI,MULD2H(ISYMB,ISYMD))
C
      CALL DZERO(WRK(KFAO),N2BASX)
      CALL OITH1(ISYMC,ZYMC,WRK(KDAO),WRK(KFAO),ISYMDM)
      IF (IBCDEQ .EQ. 24) THEN
         CALL DSCAL(N2ORBX,D3,WRK(KFAO),1)
      END IF
      CALL OITH1(ISYMD,ZYMD,WRK(KFAO),FTTI,MULD2H(ISYMB,ISYMC))
C
C     F2 contribution
C
      IF (IBCDEQ .NE. 24) THEN
         CALL DZERO(WRK(KFAO),2*N2BASX)
         CALL CDENS1(ISYMC,CMO,ZYMC,WRK(KDAO),
     *        WRK(KFREE),LFREE)
         CALL DSCAL(N2BASX,D6,WRK(KDAO),1)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'D2 in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         ISYMDM = ISYMC
         CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *        DIRFCK,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F2 in AO basis in C4FOCK'
            CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F2 in MO basis in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
C     
         KSYMOP = ISYMC
         CALL FCKOIN(1,FC,DUMMY,ZYMC,WRK(KDAO),DUMMY)
         CALL DZERO(WRK(KFAO),N2BASX)
         CALL OITH1(ISYMD,ZYMD,WRK(KDAO),WRK(KFAO),ISYMDM)
         CALL OITH1(ISYMB,ZYMB,WRK(KFAO),FTTI,MULD2H(ISYMC,ISYMD))
C     
         CALL DZERO(WRK(KFAO),N2BASX)
         CALL OITH1(ISYMB,ZYMB,WRK(KDAO),WRK(KFAO),ISYMDM)
         CALL OITH1(ISYMD,ZYMD,WRK(KFAO),FTTI,MULD2H(ISYMB,ISYMC))
C     
C     F3 contribution
C     
         CALL CDENS1(ISYMD,CMO,ZYMD,WRK(KDAO),WRK(KFREE),LFREE)
         CALL DSCAL(N2BASX,D6,WRK(KDAO),1)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'D3 in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         ISYMDM = ISYMD
         CALL SIRFCK(WRK(KFAO),WRK(KDAO),NFMATX,ISYMDM,IFCTYP,
     *        DIRFCK,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F3 in AO basis in C4FOCK'
            CALL OUTPUT(WRK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
C     
         CALL FAOMO(ISYMDM,WRK(KFAO),WRK(KDAO),CMO,WRK(KFREE),LFREE)
         IF (IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(//A)')'F3 in MO basis in C4FOCK'
            CALL OUTPUT(WRK(KDAO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
C     
         KSYMOP = ISYMD
         CALL FCKOIN(1,FC,DUMMY,ZYMD,WRK(KDAO),DUMMY)
         CALL DZERO(WRK(KFAO),N2BASX)
         CALL OITH1(ISYMC,ZYMC,WRK(KDAO),WRK(KFAO),ISYMDM)
         CALL OITH1(ISYMB,ZYMB,WRK(KFAO),FTTI,MULD2H(ISYMC,ISYMD))
C     
         CALL DZERO(WRK(KFAO),N2BASX)
         CALL OITH1(ISYMB,ZYMB,WRK(KDAO),WRK(KFAO),ISYMDM)
         CALL OITH1(ISYMC,ZYMC,WRK(KFAO),FTTI,MULD2H(ISYMB,ISYMD))
      END IF
C     
C     All contributions now gathered in FTTI, and we finish by 
C     scaling with 1/6 ( from definition of E4, without minus sign )
C     
      CALL DSCAL(NORBT*NORBT,1/D6,FTTI,1)
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(//A)')'Final result in C4FOCK'
         CALL OUTPUT(FTTI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE FAOMO(ISYMF,FAO,FMO,CMO,WRK,LWRK)
C
#include <implicit.h>
C
      DIMENSION FAO(N2BASX),FMO(N2ORBX),CMO(*),WRK(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infinp.h>
#include <wrkrsp.h>
C
      CALL DZERO(FMO,N2ORBX)
      DO 100 ISYM=1,NSYM
         JSYM   = MULD2H(ISYM,ISYMF)
         NORBI  = NORB(ISYM)
         NORBJ  = NORB(JSYM)
         IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 100
         CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     *              FAO,NBAS,NBAST,FMO,NORB,NORBT,WRK,LWRK)
 100  CONTINUE
      RETURN
      END
