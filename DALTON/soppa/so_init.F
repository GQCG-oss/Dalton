C
C  /* Deck so_init */
      SUBROUTINE SO_INIT(WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Stephan P. A. Sauer, April 2006: triplet 2p-2h pointers added
C
C     PURPOSE: Set the pointer-arrays IIJDEN, IABDEN, and IAIDEN
C              and the dimensions NIJDEN, NABDEN, and NAIDEN and
C              NIJDET, NABDET, and NAIDET as well as NSOO, NTOO, NSVV,
C              NTVV, NT2AMT1, NT2AMT2, NT2AMTT, IT2AMT1, IT2AMT2 and IT2AMT3.
C
#include "implicit.h"
C
      PARAMETER (IBIG = -100 000 000)
C
#include "priunit.h"
C
#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"
C Need DOCCSD parameter from gnring.h
#include "gnrinf.h"
C Need irat
#include "iratdef.h"

      DIMENSION WORK(*)
#ifdef VAR_MPI
#include "maxorb.h"
#include "aovec.h"
      LOGICAL so_get_herdir, herdir, so_get_direct,
     &        get_mxcall
      get_mxcall = .false.
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_INIT')
C
C-------------------------------------------
C     Initializes memory statistics routine.
C-------------------------------------------
C
      CALL SO_MEMMAX ('START',LWORK)
C
C---------------------------------------------------------------
C     If only RPA and if SOPPA is calculating the MP2 amplitude without
C     first doing CCSD calculations, then initialize some pointers and sort AO-integrals
C     into distributions.
C---------------------------------------------------------------
C
      IF (.NOT. DOCCSD) THEN
         CALL RP_INIT(WORK,LWORK)
C
#ifdef VAR_MPI
C     In pure aorpa calculations, eridi has not been initialized
         herdir = so_get_herdir()
         get_mxcall = so_get_direct()
#endif
      END IF
C
#ifdef VAR_MPI
C
C In AORPA we cannot be sure that the number of integral distributions
C is known in advance. This is needed in order to allocate memory for
C load-distribution later on
      if (get_mxcall) then
         if (herdir) then
            call herdi1(work(1), lwork, iprint)
         else
            kccfbt  = 1
            kindxbt = kccfbt + mxprim*mxcont
            kend3   = kindxbt + 8*mxshel*mxcont/irat
            lwork3 = lwork - kend3
            call eridi1(kodcl1, kodcl2, kodbc1, kodbc2,
     &                  krdbc1, krdbc2, kodpp1, kodpp2,
     &                  krdpp1, krdpp2, kfree,  lfree,
     &                  kend3, work(kccfbt), work(kindxbt),
     &                  work(kend3), lwork3, iprint)
         endif
         get_mxcall = .false.
      endif
#endif

      CALL DZERO(SOTIME,LSOTIM)
      CALL DZERO(SOWTIM,LSOWTI)
      CALL DZERO(SOORWC,LORWCI)
C
C-----------------------------------------------------------------
C     Calculate number of elements and offsets for Density matrix.
C-----------------------------------------------------------------
C
      DO 100 ISYMIJ = 1,NSYM
C
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
C
         DO 200 ISYMI = 1,NSYM
C
            ISYMJ = MULD2H(ISYMI,ISYMIJ)
C
C
            IIJDEN(ISYMI,ISYMJ) = ICOUN1
            IABDEN(ISYMI,ISYMJ) = ICOUN2
CRF I flipped the indicies on IAIDEN!
            IAIDEN(ISYMJ,ISYMI) = ICOUN3
C
            ICOUN1 = ICOUN1 + NRHF(ISYMI) * NRHF(ISYMJ)
            ICOUN2 = ICOUN2 + NVIR(ISYMI) * NVIR(ISYMJ)
C            ICOUN3 = ICOUN3 + NVIR(ISYMI) * NRHF(ISYMJ)
            ICOUN3 = ICOUN3 + NRHF(ISYMI) * NVIR(ISYMJ)
Cend-Pi
C
  200    CONTINUE
C
         NIJDEN(ISYMIJ) = ICOUN1
         NABDEN(ISYMIJ) = ICOUN2
         NAIDEN(ISYMIJ) = ICOUN3
C
  100 CONTINUE
C
      NIJDET = 0
      NABDET = 0
      NAIDET = 0
C
      DO 300 ISYMIJ = 1,NSYM
C
         NIJDET = NIJDET + NIJDEN(ISYMIJ)
         NABDET = NABDET + NABDEN(ISYMIJ)
         NAIDET = NAIDET + NAIDEN(ISYMIJ)
C
  300 CONTINUE
C
C========================================================================
C     Calculate number of elements and offsets for triplet 2p-2h vectors.
C========================================================================
C
C-----------------------
C     Initialize arrays.
C-----------------------
C
      DO I = 1, 8
         NSOO(I)    = 0
         NTOO(I)    = 0
         NSVV(I)    = 0
         NTVV(I)    = 0
         NT2AMT1(I) = 0
         NT2AMT2(I) = 0
         NT2AMT3(I) = 0
         NT2AMTT(I) = 0
         DO J = 1, 8
            IT2AMT1(I,J) = IBIG
            IT2AMT2(I,J) = IBIG
            IT2AMT3(I,J) = IBIG
            ISOO(I,J)    = IBIG
            ITOO(I,J)    = IBIG
            ISVV(I,J)    = IBIG
            ITVV(I,J)    = IBIG
         END DO
      END DO
C
C--------------------------------------------------------------------------
C     Find no. of vir/vir-pairs with a>=b (NSVV(ISYM)) and a>b (NTVV(ISYM))
C     and no. of occ/occ-pairs with i>=j (NSOO(ISYM)) and i>j (NTOO(ISYM))
C--------------------------------------------------------------------------
C
      DO ISYMIJ = 1, NSYM
C
         ICOUNSO = 0
         ICOUNTO = 0
         ICOUNSV = 0
         ICOUNTV = 0
C
         DO ISYMI = 1, NSYM
C
            ISYMJ = MULD2H(ISYMI,ISYMIJ)
C
            IF (ISYMI .EQ. ISYMJ) THEN
C
               ISOO(ISYMI,ISYMJ) = ICOUNSO
               ITOO(ISYMI,ISYMJ) = ICOUNTO
               ISVV(ISYMI,ISYMJ) = ICOUNSV
               ITVV(ISYMI,ISYMJ) = ICOUNTV
C
               ICOUNSO = ICOUNSO + NRHF(ISYMI) * (NRHF(ISYMI) + 1) / 2
               ICOUNTO = ICOUNTO + NRHF(ISYMI) * (NRHF(ISYMI) - 1) / 2
               ICOUNSV = ICOUNSV + NVIR(ISYMI) * (NVIR(ISYMI) + 1) / 2
               ICOUNTV = ICOUNTV + NVIR(ISYMI) * (NVIR(ISYMI) - 1) / 2
C
            ELSE IF (ISYMI .GT. ISYMJ) THEN
C
               ISOO(ISYMI,ISYMJ) = ICOUNSO
               ITOO(ISYMI,ISYMJ) = ICOUNTO
               ISVV(ISYMI,ISYMJ) = ICOUNSV
               ITVV(ISYMI,ISYMJ) = ICOUNTV
C
               ISOO(ISYMJ,ISYMI) = ICOUNSO
               ITOO(ISYMJ,ISYMI) = ICOUNTO
               ISVV(ISYMJ,ISYMI) = ICOUNSV
               ITVV(ISYMJ,ISYMI) = ICOUNTV
C               
               ICOUNSO = ICOUNSO + NRHF(ISYMI) * NRHF(ISYMJ)
               ICOUNTO = ICOUNTO + NRHF(ISYMI) * NRHF(ISYMJ)
               ICOUNSV = ICOUNSV + NVIR(ISYMI) * NVIR(ISYMJ)
               ICOUNTV = ICOUNTV + NVIR(ISYMI) * NVIR(ISYMJ)
C
            END IF
C
         END DO
C
         NSOO(ISYMIJ) = ICOUNSO
         NTOO(ISYMIJ) = ICOUNTO
         NSVV(ISYMIJ) = ICOUNSV
         NTVV(ISYMIJ) = ICOUNTV
C
      END DO
C
C-----------------------------------------
C     Find no. of 2p-2h triplet operators:
C     T1: for a>b and i>j
C     T2: for a>=b and i>j
C     T3: for a>b and i>=j
C-----------------------------------------
C
      DO ISYM = 1, NSYM
C
         DO ISYMAB = 1, NSYM
C
            ISYMIJ = MULD2H(ISYM,ISYMAB)
C
            NT2AMT1(ISYM) = NT2AMT1(ISYM) + NTOO(ISYMIJ)*NTVV(ISYMAB)
            NT2AMT2(ISYM) = NT2AMT2(ISYM) + NTOO(ISYMIJ)*NSVV(ISYMAB)
            NT2AMT3(ISYM) = NT2AMT3(ISYM) + NSOO(ISYMIJ)*NTVV(ISYMAB)
C
         END DO
C
         NT2AMTT(ISYM) = NT2AMT1(ISYM) + NT2AMT2(ISYM) + NT2AMT3(ISYM)
C
      END DO
C
C------------------------------------------------------------------------
C     Find offsets for of 2p-2h triplet operators: IT2AMT?(ISYMTR,ISYMAB)
C     for a given symmetry ISYMTR and
C     for a given symmetry ISYMAB of the pair of virtual orbitals ab
C     T1: for a>b and i>j
C     T2: for a>=b and i>j
C     T3: for a>b and i>=j
C------------------------------------------------------------------------
C
      DO ISYMTR = 1, NSYM
C
         ICOUNT1 = 0
         ICOUNT2 = 0
         ICOUNT3 = 0
C
         DO ISYMIJ = 1, NSYM
C
            ISYMAB = MULD2H(ISYMTR,ISYMIJ)
C
            IT2AMT1(ISYMIJ,ISYMAB) = ICOUNT1
            IT2AMT2(ISYMIJ,ISYMAB) = NT2AMT1(ISYMTR) + ICOUNT2
            IT2AMT3(ISYMIJ,ISYMAB) = NT2AMT1(ISYMTR) + NT2AMT2(ISYMTR)
     &                               + ICOUNT3
C
            ICOUNT1 = ICOUNT1 + NTOO(ISYMIJ) * NTVV(ISYMAB)
            ICOUNT2 = ICOUNT2 + NTOO(ISYMIJ) * NSVV(ISYMAB)
            ICOUNT3 = ICOUNT3 + NSOO(ISYMIJ) * NTVV(ISYMAB)
C
         END DO
C
      END DO
CPFP**
      IF (IPRSOP .GE. 5) THEN
         WRITE(LUPRI,'(//,1X,A,/)')
     &        'Counter and offset arrays for triplet 2p-2h vectors:'
         WRITE(LUPRI,'(A,8I4)') 'NSOO ',(NSOO(I),I=1,8)
         WRITE(LUPRI,'(A,8I4)') 'NTOO ',(NTOO(I),I=1,8)
         WRITE(LUPRI,'(A,8I4)') 'NSVV ',(NSVV(I),I=1,8)
         WRITE(LUPRI,'(A,8I4)') 'NTVV ',(NTVV(I),I=1,8)
C
         WRITE(LUPRI,'(8(A,8I4,/))') ('ISOO ',(ISOO(I,J),J=1,8),I=1,8)
         WRITE(LUPRI,'(8(A,8I4,/))') ('ITOO ',(ITOO(I,J),J=1,8),I=1,8)
         WRITE(LUPRI,'(8(A,8I4,/))') ('ISVV ',(ISVV(I,J),J=1,8),I=1,8)
         WRITE(LUPRI,'(8(A,8I4,/))') ('ITVV ',(ITVV(I,J),J=1,8),I=1,8)
C
         WRITE(LUPRI,'(A,8I4)') 'NT2AMT1 ',(NT2AMT1(I),I=1,8)
         WRITE(LUPRI,'(A,8I4)') 'NT2AMT2 ',(NT2AMT2(I),I=1,8)
         WRITE(LUPRI,'(A,8I4)') 'NT2AMT3 ',(NT2AMT3(I),I=1,8)
         WRITE(LUPRI,'(A,8I4)') 'NT2AMTT ',(NT2AMTT(I),I=1,8)
C
         WRITE(LUPRI,'(8(A,8I4,/))')
     &        ('IT2AMT1 ',(IT2AMT1(I,J),J=1,8),I=1,8)
         WRITE(LUPRI,'(8(A,8I4,/))')
     &        ('IT2AMT2 ',(IT2AMT2(I,J),J=1,8),I=1,8)
         WRITE(LUPRI,'(8(A,8I4,/))')
     &        ('IT2AMT3 ',(IT2AMT3(I,J),J=1,8),I=1,8)
         CALL FLUSH(LUPRI)
      END IF
Cend-PFP**
C
C---------------------------------------------------------------------
C     Set no. of 2p-2h operators N2P2HOP for singlet and triplet case:
C---------------------------------------------------------------------
C
      DO ISYMTR = 1, NSYM
C
         IF (TRIPLET) THEN
C
            N2P2HOP(ISYMTR) = NT2AMTT(ISYMTR)
C
         ELSE
C
            N2P2HOP(ISYMTR) = NT2AM(ISYMTR)
C
         END IF
C
      END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_INIT')
C
      RETURN
      END
