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
*=======================================================================
      SUBROUTINE CC_R12NO(T1AM,T2AM,WORK,LWORK)
*-----------------------------------------------------------------------
* Purpose: Calculate semi-natural virtual orbitals for R12
*
* Note: T1AM, T2AM will be overwritten! T2AM is expected to contain
*       the integrals (ia|jb) on entry!
*
* Christian Neiss, autumn 2005
*-----------------------------------------------------------------------
C
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <inftap.h>
#include <r12int.h>
C
      INTEGER KFOCKD,KYMAT,KEIVAL,KEIVEC,KT2SQ,KSCR,KFCKR12
      INTEGER KCMO,KFCKHF,LUFCK,KOFF1,KOFF2,KOFF3,KCMOTR
      INTEGER LWORK,LWRK1,LWRK2,KEND1,KEND2
      INTEGER MATZ,KFV1,KFV2,IERR,ISYM,IDUMMY
      INTEGER NFCKR12,IFCKR12(8),ICOUNT,IVIR1(8),KFOCKX,KCMOX
      INTEGER NSYMX,NORBTSX,NBASTX,NLAMDSX,NRHFSX(8),NORBSX(8)
      INTEGER NR12xVIR,IR12xVIR(8),KFCKMIX,LUNIT
      LOGICAL LOCDBG
#if defined (SYS_CRAY)
      REAL WORK(LWORK),T1AM(NT1AMX),T2AM(NT2AMX)
#else
      DOUBLE PRECISION WORK(LWORK),T1AM(NT1AMX),T2AM(NT2AMX)
#endif
C
      PARAMETER (LOCDBG = .FALSE.)       
C
      CALL QENTER('CC_R12NO')
      IF (LOCDBG) WRITE(LUPRI,*)'Entered CC_R12NO'
C
      NFCKR12 = 0
      ICOUNT  = 0
      DO ISYM = 1, NSYM
        IFCKR12(ISYM) = NFCKR12
        IVIR1(ISYM)   = ICOUNT
        NFCKR12 = NFCKR12 + NRXR12(ISYM)*NRXR12(ISYM)
        ICOUNT  = ICOUNT + NVIR(ISYM)
      END DO
C
      KFOCKD = 1
      KCMO   = KFOCKD + NORBTS
      KFCKHF = KCMO  + NLAMDS
      KYMAT  = KFCKHF + N2BAST
      KEIVAL = KYMAT + NMATAB(1)
      KEIVEC = KEIVAL + NVIRT 
      KFCKR12= KEIVEC + NMATAB(1)
      KEND1  = KFCKR12 + NFCKR12
      LWRK1  = LWORK - KEND1
      IF (LWRK1.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
C
      !read canonical orbital energies and CMO matrix:
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)  
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      READ (LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP') 
      IF (FROIMP .OR. FROEXP)
     *  CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'vir-vir block of CMO matrix:'
        DO ISYM = 1, NSYM
          KOFF1 = KCMO + IGLMVI(ISYM,ISYM)
          CALL OUTPUT(WORK(KOFF1),1,NBAS(ISYM),1,NVIR(ISYM),
     &                NBAS(ISYM),NVIR(ISYM),1,LUPRI)
        END DO
      END IF
C
      !read SCF AO-Fock matrix:
      LUFCK = -1
      CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *           IDUMMY,.FALSE.)
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(1))
      CALL GPCLOSE(LUFCK,'KEEP' )
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'AO-Fock matrix read in:'
        DO ISYM = 1, NSYM
          WRITE(LUPRI,*)' Symmetry block: ',ISYM
          CALL OUTPUT(WORK(KFCKHF+IAODIS(ISYM,ISYM)),1,NBAS(ISYM),
     &                1,NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),1,LUPRI)
        END DO
      END IF
C
      !make MP2-Guess:
      CALL CCSD_GUESS(T1AM,T2AM,WORK(KFOCKD),IPRINT)
C
      !make 2C-E combination of T2AM:
      KT2SQ = KEND1
      KEND2 = KT2SQ + NT2SQ(1)
      LWRK2 = LWORK - KEND2
      IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
      CALL CC_T2SQ(T2AM,WORK(KT2SQ),1)
      CALL CCRHS_T2TR(WORK(KT2SQ),WORK(KEND2),LWRK2,1)
C     CALL DSCAL(NT2SQ(1),2.0D0,WORK(KT2SQ),1)
C
      !calculate virtual-virtual block of density matrix:
      CALL CC_YI(WORK(KYMAT),WORK(KT2SQ),1,T2AM,1,
     &           WORK(KEND2),LWRK2)
C     CALL CC_PRSQ(IDUMMY,WORK(KT2SQ),1,0,1)
C     CALL CC_PRP(IDUMMY,T2AM,1,0,1)
C     WRITE(LUPRI,*)'YMAT in CC_R12NO:'
C     DO ISYM = 1, NSYM
C       WRITE(LUPRI,*)'symmetry block ',ISYM
C       CALL OUTPUT(WORK(KYMAT+IMATAB(ISYM,ISYM)),1,NVIR(ISYM),
C    &              1,NVIR(ISYM),NVIR(ISYM),NVIR(ISYM),
C    &              1,LUPRI)
C     END DO
C
      !allocate work space for transformed CMO matrix:
      KCMOTR = KEND1
      KEND1  = KCMOTR + NLAMDS
      LWRK1  = LWORK - KEND1
      IF (LWRK1.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
      CALL DCOPY(NLAMDS,WORK(KCMO),1,WORK(KCMOTR),1)
C
      !symmetrize and diagonalize:
C     CALL DSCAL(NMATAB(1),0.5D0,WORK(KYMAT),1)
      DO ISYM = 1, NSYM
        KSCR  = KEND1
        KEND2 = KSCR + NVIR(ISYM)*NVIR(ISYM)
        LWRK2 = LWORK - KEND2
        IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
C        
        CALL TRSREC(NVIR(ISYM),NVIR(ISYM),
     &              WORK(KYMAT+IMATAB(ISYM,ISYM)),WORK(KSCR)) 
        CALL DAXPY(NVIR(ISYM)*NVIR(ISYM),1.0D0,WORK(KSCR),1,
     &             WORK(KYMAT+IMATAB(ISYM,ISYM)),1)
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'vir-vir block of density matrix before diag.:'
          CALL OUTPUT(WORK(KYMAT+IMATAB(ISYM,ISYM)),1,NVIR(ISYM),
     &                1,NVIR(ISYM),NVIR(ISYM),NVIR(ISYM),1,LUPRI)
        END IF
C
        MATZ = 1
        KFV1  = KEND1
        KFV2  = KFV1 + NVIR(ISYM)
        KEND2 = KFV2 + NVIR(ISYM)
        LWRK2 = LWORK - KEND2
        IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
        CALL RS(NVIR(ISYM),NVIR(ISYM),WORK(KYMAT+IMATAB(ISYM,ISYM)),
     &          WORK(KEIVAL+IVIR1(ISYM)),
     &          MATZ,WORK(KEIVEC+IMATAB(ISYM,ISYM)),
     &          WORK(KFV1),WORK(KFV2),IERR)
        IF ( IERR.NE.0 ) THEN
          WRITE(LUPRI,'(/A,I5)')
     *    ' EIGENVALUE PROBLEM NOT CONVERGED IN RS, IERR =',IERR
          CALL QUIT(' CCRED: EIGENVALUE EQUATION NOT CONVERGED ')
        END IF
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'Eigenvalues before reorder: symmetry=',ISYM
          CALL OUTPUT(WORK(KEIVAL+IVIR1(ISYM)),1,NVIR(ISYM),
     &                1,1,NVIR(ISYM),1,1,LUPRI)
          WRITE(LUPRI,*)'Eigenvectors before reorder:'
          CALL OUTPUT(WORK(KEIVEC+IMATAB(ISYM,ISYM)),1,NVIR(ISYM),
     &                1,NVIR(ISYM),NVIR(ISYM),NVIR(ISYM),1,LUPRI)
        END IF
        CALL RGORD(NVIR(ISYM),NVIR(ISYM),WORK(KEIVAL+IVIR1(ISYM)),
     &             WORK(KFV1),WORK(KEIVEC+IMATAB(ISYM,ISYM)),.TRUE.)
        WRITE(LUPRI,*)'R12-NATVIR: Nat. occ. numbers after ordering: '//
     &                'symmetry=',ISYM
        IF (NVIR(ISYM).EQ.0) THEN
          WRITE(LUPRI,*) 'This block is empty'
        ELSE
          CALL OUTPUT(WORK(KEIVAL+IVIR1(ISYM)),1,MIN(10,NVIR(ISYM)),
     &                1,1,NVIR(ISYM),1,1,LUPRI)
        END IF
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'Eigenvectors after reorder:'
          CALL OUTPUT(WORK(KEIVEC+IMATAB(ISYM,ISYM)),1,NVIR(ISYM),
     &                1,NVIR(ISYM),NVIR(ISYM),NVIR(ISYM),1,LUPRI)
          CALL FLSHFO(LUPRI)
        END IF
C
        !calculate CMO-matrix for natural virtual orbitals:
        KOFF1 = KCMO + IGLMVI(ISYM,ISYM)
        KOFF2 = KEIVEC+IMATAB(ISYM,ISYM)
        KOFF3 = KCMOTR + IGLMVI(ISYM,ISYM)
        CALL DGEMM('N','N',NBAS(ISYM),NRXR12(ISYM),NVIR(ISYM),1.0D0,
     &             WORK(KOFF1),MAX(1,NBAS(ISYM)),
     &             WORK(KOFF2),MAX(1,NVIR(ISYM)),
     &             0.0D0,WORK(KOFF3),MAX(1,NBAS(ISYM)))
C
      END DO
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'vir-vir block of CMO matrix with natural '//
     &                'virtual R12-orbitals:'
        DO ISYM = 1, NSYM
          KOFF1 = KCMOTR + IGLMVI(ISYM,ISYM)
          CALL OUTPUT(WORK(KOFF1),1,NBAS(ISYM),1,NVIR(ISYM),
     &                NBAS(ISYM),NVIR(ISYM),1,LUPRI)
        END DO
      END IF
C
      !transform Fock matrix into subspace of additional r12-orbitals 
      !(which are chosen from the natural vir. orbials just generated):
      DO ISYM = 1, NSYM
        KOFF1 = KCMOTR + IGLMVI(ISYM,ISYM)
        KOFF2 = KFCKHF + IAODIS(ISYM,ISYM)
        KOFF3 = KEND1
        CALL DGEMM('N','N',NBAS(ISYM),NRXR12(ISYM),NBAS(ISYM),1.0D0,
     &             WORK(KOFF2),MAX(1,NBAS(ISYM)),
     &             WORK(KOFF1),MAX(1,NBAS(ISYM)),0.0D0,
     &             WORK(KOFF3),MAX(1,NBAS(ISYM)))
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'CMO part used:'
          WRITE(LUPRI,*)'symmetry block: ',ISYM
          CALL OUTPUT(WORK(KOFF1),1,NBAS(ISYM),
     &                1,NRXR12(ISYM),NBAS(ISYM),NVIR(ISYM),
     &                1,LUPRI)
          WRITE(LUPRI,*)'AO-Fock part used:'
          WRITE(LUPRI,*)'symmetry block: ',ISYM
          CALL OUTPUT(WORK(KOFF2),1,NBAS(ISYM),
     &                1,NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),
     &                1,LUPRI)
          WRITE(LUPRI,*)'half-transformed R12-Fock-part:'
          WRITE(LUPRI,*)'symmetry block: ',ISYM
          CALL OUTPUT(WORK(KOFF3),1,NBAS(ISYM),
     &                1,NRXR12(ISYM),NBAS(ISYM),NRXR12(ISYM),
     &                1,LUPRI)
        END IF
        CALL DGEMM('T','N',NRXR12(ISYM),NRXR12(ISYM),NBAS(ISYM),1.0D0,
     &             WORK(KOFF3),MAX(1,NBAS(ISYM)),
     &             WORK(KOFF1),MAX(1,NBAS(ISYM)),0.0D0,
     &             WORK(KFCKR12+IFCKR12(ISYM)),MAX(1,NRXR12(ISYM)))
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'R12-Fock-part before diagonalization:'
          WRITE(LUPRI,*)'symmetry block: ',ISYM
          CALL OUTPUT(WORK(KFCKR12+IFCKR12(ISYM)),1,NRXR12(ISYM),
     &                1,NRXR12(ISYM),NRXR12(ISYM),NRXR12(ISYM),
     &                1,LUPRI)
        END IF
        !diagonalize:
        MATZ = 1
        KFV1  = KEND1
        KFV2  = KFV1 + NRXR12(ISYM)
        KEND2 = KFV2 + NRXR12(ISYM)
        LWRK2 = LWORK - KEND2
        IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
        CALL RS(NRXR12(ISYM),NRXR12(ISYM),WORK(KFCKR12+IFCKR12(ISYM)),
     &          WORK(KEIVAL+IVIR1(ISYM)),
     &          MATZ,WORK(KEIVEC+IMATAB(ISYM,ISYM)),
     &          WORK(KFV1),WORK(KFV2),IERR)
        IF ( IERR.NE.0 ) THEN
          WRITE(LUPRI,'(/A,I5)')
     *    ' EIGENVALUE PROBLEM NOT CONVERGED IN RS, IERR =',IERR
          CALL QUIT(' CCRED: EIGENVALUE EQUATION NOT CONVERGED ')
        END IF
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'Fock-Eigenvalues after RS: symmetry=',ISYM
          CALL OUTPUT(WORK(KEIVAL+IVIR1(ISYM)),1,NRXR12(ISYM),
     &                1,1,NVIR(ISYM),1,1,LUPRI)
          WRITE(LUPRI,*)'Eigenvectors after RS:'
          CALL OUTPUT(WORK(KEIVEC+IMATAB(ISYM,ISYM)),1,NRXR12(ISYM),
     &                1,NRXR12(ISYM),NRXR12(ISYM),NRXR12(ISYM),1,LUPRI)
        END IF
        CALL RGORD(NRXR12(ISYM),NRXR12(ISYM),WORK(KEIVAL+IVIR1(ISYM)),
     &             WORK(KFV1),WORK(KEIVEC+IMATAB(ISYM,ISYM)),.FALSE.)
        WRITE(LUPRI,*)
        WRITE(LUPRI,*)'R12-NATVIR: Fock-Eigenvalues after ordering: '//
     &                'symmetry=',ISYM
        IF (NRXR12(ISYM).EQ.0) THEN
          WRITE(LUPRI,*) 'This block is empty'
        ELSE IF (NVIR(ISYM).LT.NRXR12(ISYM)) THEN
          WRITE(LUPRI,*) 'You specified more virtual R12 orbitals '//
     &                   'than there are virtual orbitals in this '//
     &                   'symmery:'
          WRITE(LUPRI,*) 'NVIR(',ISYM,') = ',NVIR(ISYM)
          WRITE(LUPRI,*) 'NRXR12(',ISYM,') = ',NRXR12(ISYM)
          CALL QUIT('Too many virtual R12 orbitals')
        ELSE
          CALL OUTPUT(WORK(KEIVAL+IVIR1(ISYM)),1,NRXR12(ISYM),
     &                1,1,NVIR(ISYM),1,1,LUPRI)
        END IF
        IF (LOCDBG) THEN
          WRITE(LUPRI,*)'Eigenvectors after reorder:'
          CALL OUTPUT(WORK(KEIVEC+IMATAB(ISYM,ISYM)),1,NRXR12(ISYM),
     &                1,NRXR12(ISYM),NRXR12(ISYM),NRXR12(ISYM),1,LUPRI)
          CALL FLSHFO(LUPRI)
        END IF
 
        !calculate semi-natural virtual R12-part of CMO-matrix :
        KOFF1 = KCMOTR + IGLMVI(ISYM,ISYM)
        KOFF2 = KEIVEC+IMATAB(ISYM,ISYM)
        KOFF3 = KEND1
        CALL DGEMM('N','N',NBAS(ISYM),NRXR12(ISYM),NRXR12(ISYM),1.0D0,
     &             WORK(KOFF1),MAX(1,NBAS(ISYM)),
     &             WORK(KOFF2),MAX(1,NRXR12(ISYM)),
     &             0.0D0,WORK(KOFF3),MAX(1,NBAS(ISYM)))
        CALL DCOPY(NBAS(ISYM)*NRXR12(ISYM),WORK(KOFF3),1,WORK(KOFF1),1)
      END DO
C
      WRITE(LUPRI,*)
C
      IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,*)'R12-NATVIR: MO coefficient matrix of '//
     &                'semi-natural virtual R12-orbitals:'
        DO ISYM = 1, NSYM
          IF (NRXR12(ISYM).GT.0) THEN
            WRITE(LUPRI,*) 'symmetry=',ISYM
            KOFF1 = KCMOTR + IGLMVI(ISYM,ISYM)
            CALL OUTPUT(WORK(KOFF1),1,NBAS(ISYM),1,NRXR12(ISYM),
     &                  NBAS(ISYM),NRXR12(ISYM),1,LUPRI)
          END IF
        END DO
      END IF
C
      !generate mixed r12-vir Fock-matrix and write it to file:
      NR12xVIR = 0
      DO ISYM = 1, NSYM
        IR12xVIR(ISYM) = NR12xVIR
        NR12xVIR = NR12xVIR + NVIR(ISYM)*NRXR12(ISYM)
      END DO
C
      KFCKMIX = KEND1
      KEND2   = KFCKMIX + NR12xVIR
      LWRK2   = LWORK - KEND2
      IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
      DO ISYM = 1, NSYM
        KSCR  = KEND2
        KEND2 = KSCR + NBAS(ISYM)*NRXR12(ISYM)
        LWRK2   = LWORK - KEND2
        IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
        KOFF1 = KCMOTR + IGLMVI(ISYM,ISYM)
        KOFF2 = KFCKHF + IAODIS(ISYM,ISYM)
        KOFF3 = KSCR
        CALL DGEMM('N','N',NBAS(ISYM),NRXR12(ISYM),NBAS(ISYM),1.0D0,
     &             WORK(KOFF2),MAX(1,NBAS(ISYM)),
     &             WORK(KOFF1),MAX(1,NBAS(ISYM)),0.0D0,
     &             WORK(KOFF3),MAX(1,NBAS(ISYM)))
        KOFF1 = KCMO + IGLMVI(ISYM,ISYM)
        CALL DGEMM('T','N',NVIR(ISYM),NRXR12(ISYM),NBAS(ISYM),1.0D0,
     &             WORK(KOFF1),MAX(1,NBAS(ISYM)),
     &             WORK(KOFF3),MAX(1,NBAS(ISYM)),0.0D0,
     &             WORK(KFCKMIX+IR12xVIR(ISYM)),MAX(1,NVIR(ISYM)))
      END DO
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'R12 x VIR FOCK MATRIX:'
        DO ISYM = 1, NSYM
          WRITE(LUPRI,*) 'Symmetry block: ',ISYM
          CALL OUTPUT(WORK(KFCKMIX+IR12xVIR(ISYM)),1,NVIR(ISYM),
     &                1,NRXR12(ISYM),NVIR(ISYM),NRXR12(ISYM),
     &                1,LUPRI)
        END DO
      END IF
      LUNIT = -1
      CALL GPOPEN(LUNIT,'R12FVIR','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE (LUNIT) (WORK(KFCKMIX+I-1), I=1, NR12xVIR)
      CALL GPCLOSE(LUNIT,'KEEP')
C
      !merge results to SIRIUS interface file:
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('FULLBAS ',LUSIFC,LUPRI)
      READ (LUSIFC) NSYMX,NORBTSX,NBASTX,NLAMDSX,(NRHFSX(I),I=1,NSYM),
     &              (NORBSX(I),I=1,NSYM)
      KFOCKX = KEND1
      KCMOX  = KFOCKX + NORBTSX
      KEND2  = KCMOX + NLAMDSX
      LWRK2  = LWORK - KEND2
      IF (LWRK2.LT.0) CALL QUIT ('Insufficient memory in CC_R12NO')
      READ (LUSIFC) (WORK(KFOCKX+I-1), I=1,NORBTSX)
      READ (LUSIFC) (WORK(KCMOX+I-1),I=1,NLAMDSX)
C     CALL GPCLOSE(LUSIFC,'KEEP')
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'Orbital energies read in:'
        WRITE(LUPRI,*) (WORK(KFOCKX+I-1),I=1,NORBTSX)
      END IF
C
      KOFF2 = KFOCKX 
      DO ISYM = 1, NSYM
        KOFF1 = KEIVAL + IVIR1(ISYM)
        KOFF2 = KOFF2 + NRHFSA(ISYM)
        CALL DCOPY(NRXR12(ISYM),WORK(KOFF1),1,WORK(KOFF2),1)
        KOFF2 = KOFF2 + NRXR12(ISYM) + NRHFSA(ISYM)
C       CALL DCOPY(NRXR12(ISYM),WORK(KOFF1),1,WORK(KOFF2),1)
        KOFF2 = KOFF2 + NVIR(ISYM) + NORB2(ISYM)
      END DO
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'Orbital energies write out:'
        WRITE(LUPRI,*) (WORK(KFOCKX+I-1),I=1,NORBTSX)
      END IF
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'MO-coefficient matrix read in:'
        KOFF1 = KCMOX
        DO ISYM = 1, NSYM
          WRITE(LUPRI,*)' Symmetry block: ',ISYM
          CALL OUTPUT(WORK(KOFF1),1,MBAS1(ISYM)+MBAS2(ISYM),
     &                1,NORBSX(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
     &                NORBSX(ISYM),1,LUPRI)
          KOFF1 = KOFF1 + NORBSX(ISYM)*(MBAS1(ISYM)+MBAS2(ISYM))
        END DO
      END IF   
      KOFF2 = KCMOX
C     KOFF3 = KCMOX
      DO ISYM = 1, NSYM
        KOFF1 = KCMOTR + IGLMVI(ISYM,ISYM)
        KOFF2 = KOFF2 + NRHFSA(ISYM)*(MBAS1(ISYM)+MBAS2(ISYM))
C       KOFF3 = KOFF3 + NRHFSB(ISYM)*(MBAS1(ISYM)+MBAS2(ISYM)) + 
C    &          NRHFSA(ISYM)*(MBAS1(ISYM)+MBAS2(ISYM))
        DO I = 1, NRXR12(ISYM)
          CALL DCOPY(MBAS1(ISYM),WORK(KOFF1),1,WORK(KOFF2),1)
          KOFF2 = KOFF2 + MBAS1(ISYM)+MBAS2(ISYM)
C         CALL DCOPY(MBAS1(ISYM),WORK(KOFF1),1,WORK(KOFF3),1)
          KOFF1 = KOFF1 + MBAS1(ISYM)
C         KOFF3 = KOFF3 + MBAS1(ISYM)+MBAS2(ISYM)
        END DO
        KOFF2 = KOFF2 + (NORB1(ISYM)+NORB2(ISYM))*
     &                  (MBAS1(ISYM)+MBAS2(ISYM))
C       KOFF3 = KOFF3 + NORB2(ISYM)*(MBAS1(ISYM)+MBAS2(ISYM))
      END DO
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'MO-coefficient matrix write out:'
        KOFF1 = KCMOX
        DO ISYM = 1, NSYM
          WRITE(LUPRI,*)' Symmetry block: ',ISYM
          CALL OUTPUT(WORK(KOFF1),1,MBAS1(ISYM)+MBAS2(ISYM),
     &                1,NORBSX(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
     &                NORBSX(ISYM),1,LUPRI)
          KOFF1 = KOFF1 + NORBSX(ISYM)*(MBAS1(ISYM)+MBAS2(ISYM))
        END DO
      END IF
C     CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
C    &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('FULLBAS ',LUSIFC,LUPRI)
      READ (LUSIFC) 
      WRITE (LUSIFC) (WORK(KFOCKX+I-1), I=1,NORBTSX)
      WRITE (LUSIFC) (WORK(KCMOX+I-1),I=1,NLAMDSX)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C     Compute new matrix <ij|r12**2|kl>
      CALL R12AUX1(WORK(KCMOX),WORK(KEND2),LWRK2)
C
      IF (LOCDBG) WRITE(LUPRI,*)'Leaving CC_R12NO'
      CALL QEXIT('CC_R12NO')
      RETURN
      END
*=======================================================================

