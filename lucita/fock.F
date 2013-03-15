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
C 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE FI(FIMAT,ECC,IDOH2)
*
*. Construct inactive fockmatrix + core-core interaction energy.
*. I.e. add contributions from all orbitals
*  that belong to hole orbital spaces ( as defined by IPHGAS).
*
* Note that this is a more general definition of the
* Inactive Fockmatrix than usually used.
*
*. On input FIMAT should be the inactive Fock matrix, in symmetry packed form
*
* Jeppe Olsen
*
* Revision : Dec 97 : General hole spaces
      ImplICIT REAL*8(A-H,O-Z)
*
      DIMENSION FIMAT(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "lucinp.inc"
#include "cgas.inc"
*
*
      CALL FIH(FIMAT,ECC,IBSO,NSMOB,ITPFSO,IPHGAS,NTOOBS,NTOOB,IREOST,
     &         IDOH2)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE FIFAM(FIFA)
*
*. Construct inactive + active fock matrix
*
*. On input FIFAM Should be the inactive Fock matrix, in symmetry packed form
*
* Jeppe Olsen
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION FIFA(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "glbbas.inc"
#include "lucinp.inc"
*
*
      CALL FIFAMS(FIFA,WORK(KRHO1),IBSO,NSMOB,
     &            NTOOBS,NACOB,NTOOB,IREOST)
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE FIFAMS(FIFA,RHO1,IOBSM,NSMOB,LOBSM,NACOB,
     &                  NORBT,ISTOB)
*
* Update inactive fock matrix with active contributions
*
*     FIFA(I,J) = FIFA(I,J) + sum(k,l) ((ij!ab)-0.5*(ib!ja))
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FIFA(*),RHO1(NACOB,NACOB)
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
*
*.  Assume spatial symmetric fock matrix
      IJSM = 1
      IJ = 0
      DO ISM = 1, NSMOB
        CALL SYMCOM(2,6,ISM,JSM,IJSM)
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        IF(JSM.NE.0) THEN
          DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
            DO J = IOBSM(JSM),I
C?            write(6,*) ' I J ', I,J
              IP = ISTOB(I)
              JP = ISTOB(J)
C?            write(6,*) ' IP JP ', IP,JP
               IJ= IJ + 1
               DO IA = 1, NACOB
                 DO IB = 1, NACOB
                   FIFA(IJ) = FIFA(IJ)
     &           + RHO1(IA,IB)
     &           *(GTIJKL(IP,JP,IA,IB)-0.5*GTIJKL(IP,IB,IA,JP))
                 END DO
               END DO
            END DO
          END DO
        END IF
      END DO
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' FI + FA in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*)
       ISYM = 1
       CALL APRBLM2(FIFA,LOBSM,LOBSM,NSMOB,ISYM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE FIH(FI,ECC,IOBSM,NSMOB,ITPFSO,IPHGAS,LOBSM,NORBT,ISTOB,
     &               IDOH2)
*
* construct inactive fock matrix
*
*     FI(I,J) = FI(I,J) + sum(h) (2(ij!hh)-(ih!jh))
*
* where h is summed over all hole orbitals (as declaed by IPHGAS)
* Note that this is a more general definition of the Inactive
* Fock matrix than usually used.
* (Normal realization : see FIS )
*
* Jeppe Olsen ( I admit )
*
* Dec 97
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FI(*)
      INTEGER IOBSM(*),LOBSM(*),ISTOB(*)
      INTEGER ITPFSO(*), IPHGAS(*)

*
      NTEST = 0
C?    WRITE(6,*) ' FIH: IPHGAS and ITPFSO '
C?    CALL IWRTMA(IPHGAS,1,2,1,2)
C?    CALL IWRTMA(ITPFSO,1,3,1,3)
*
* Core-Core energy
*
      ECC = 0.0D0
      IJSM = 1
*. One-electron part
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE
          IIOFF = IIOFF + LOBSM(ISM-1)*(LOBSM(ISM-1)+1)/2
        END IF
        II = IIOFF-1
        DO I = IOBSM(ISM),IOBSM(ISM)+LOBSM(ISM)-1
          II = II + (I-IOBSM(ISM)+1)
          IF(IPHGAS(ITPFSO(I)).EQ.2) ECC = ECC + 2*FI(II)
        END DO
      END DO
C?    WRITE(6,*) ' one-electron part to ECC ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IOBSM(ISM), IOBSM(ISM) + LOBSM(ISM)-1
          DO J = IOBSM(JSM), IOBSM(JSM) + LOBSM(JSM)-1
              IP = ISTOB(I)
              JP = ISTOB(J)
              IF(IPHGAS(ITPFSO(I)).EQ.2.AND.IPHGAS(ITPFSO(J)).EQ.2)
     &        ECC = ECC +2*GTIJKL(IP,IP,JP,JP)-GTIJKL(IP,JP,JP,IP)
          END DO
          END DO
        END DO
        END DO
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
      IF(IDOH2.NE.0) THEN
        IJSM = 1
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IF(JSM.NE.0) THEN
            DO I = IOBSM(ISM),IOBSM(ISM) + LOBSM(ISM)-1
              DO J = IOBSM(JSM),I
                IP = ISTOB(I)
                JP = ISTOB(J)
                IJ= IJ + 1
                DO KSYM = 1, NSMOB
                  DO K = IOBSM(KSYM),IOBSM(KSYM)-1+LOBSM(KSYM)
                    KP = ISTOB(K)
                    IF(IPHGAS(ITPFSO(K)).EQ.2) FI(IJ) = FI(IJ)
     &            + 2.0D0*GTIJKL(IP,JP,KP,KP)-GTIJKL(IP,KP,KP,JP)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
*
      IF(NTEST.NE.0) THEN
*
       WRITE(6,*) ' FI in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*)
       ISYM = 1
       CALL APRBLM2(FI,LOBSM,LOBSM,NSMOB,ISYM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE FOCK_MAT(F,I12)
*
* Construct Fock matrix
*
* F(i,j) = SUM(K) H(i,K) * RHO1(j,K)
*          + SUM(M,K,L) I  (i M K L ) * RHO2( j M K L )
*
* Helsingfors, december 11 (1996)
* (after the EFG Winter School)
*
* Unless I12 = 2, only one-electron part is calculated
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLINT, KLDEN, KLFBLK, IIOFF, IDOFF
!               for addressing of WORK
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cgas.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Output
      DIMENSION F(*)
*
      NTEST = 1
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'FOO   ')
*
      ONE = 1.0D0
      ZERO = 0.0D0
*. Largest set of orbitals with given symmetry and type
CTF
* Using local MXTSOB_L (MXTSOB is now a common parameter!)
      MXTSOB_L = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IGAS,ISM))
      END DO
      END DO
C?    WRITE(6,*) 'MXTSOB_L = ', MXTSOB_L
*. Allocate scratch space for 2-electron integrals and
*. two-electron densities
      MX4IBLK = MXTSOB_L ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT ')
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN ')
*. And a block of F
      MX2IBLK = MXTSOB_L ** 2
      CALL MEMMAN(KLFBLK,MX2IBLK,'ADDL  ',2,'KLFBL ')
*.
*
      ONE = 1.0D0
      DO IJSM = 1, NSMOB
        ISM = IJSM
        JSM = IJSM
        NIJS = NOCOBS(IJSM)
*
        IF(IJSM.EQ.1) THEN
         IFOFF = 1
        ELSE
         IFOFF = IFOFF+NOCOBS(IJSM-1)**2
        END IF
*
        DO JGAS = 1, NGAS
          IF(JGAS.EQ.1) THEN
            IJ = 1
          ELSE
            IJ = IJ + NOBPTS(JGAS-1,JSM)
          END IF
          NJ = NOBPTS(JGAS,IJSM)
          DO IGAS = 1, NGAS
            IF(NTEST.GE.1000) THEN
              WRITE(6,*)
     &        ' Fock matrix for ISM IGAS JGAS',ISM,IGAS,JGAS
            END IF
            NI = NOBPTS(IGAS,ISM)
            IF(IGAS.EQ.1) THEN
              II = 1
            ELSE
              II = II + NOBPTS(IGAS-1,ISM)
            END IF
*
*  =======================
*. block F(ijsm,igas,jgas)
*  =======================
*
            CALL SETVEC(WORK(KLFBLK),ZERO,NI*NJ)
* 1 : One-electron part
            DO KGAS = 1, NGAS
              KSM = IJSM
              NK = NOBPTS(KGAS,KSM)
*. blocks of one-electron integrals and one-electron density
              CALL GETD1(WORK(KLDEN),JSM,JGAS,KSM,KGAS)
              CALL GETH1(WORK(KLINT),ISM,IGAS,KSM,KGAS)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*)
     &          ' 1-e ints for ISM IGAS KGAS ',ISM,IGAS,KGAS
                CALL WRTMT_LU(WORK(KLINT),NI,NK,NI,NK)
                WRITE(6,*)
     &          ' 1-e densi for ISM JGAS KGAS ',ISM,JGAS,KGAS
                CALL WRTMT_LU(WORK(KLDEN),NJ,NK,NJ,NK)
              END IF
*. And then a matrix multiply( they are pretty much in fashion
*. these days )
              CALL MATML7(WORK(KLFBLK),WORK(KLINT),WORK(KLDEN),
     &                    NI,NJ,NI,NK,NJ,NK,ONE,ONE,2)
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated block '
                 CALL WRTMT_LU(WORK(KLFBLK),NI,NJ,NI,NJ)
               END IF

            END DO
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One-electron contributions'
              WRITE(6,*) ' =========================='
              CALL WRTMT_LU(WORK(KLFBLK),NI,NJ,NI,NJ)
            END IF
            IF(I12.EQ.2) THEN
*. 2 : Two-electron part
            DO KSM = 1, NSMOB
            DO LSM = 1, NSMOB
*. Obtain MSM
              CALL  SYMCOM(3,1,KSM,LSM,KLSM)
              CALL  SYMCOM(3,1,KLSM,ISM,IKLSM)
              IMKLSM = 1
              CALL  SYMCOM(2,1,IKLSM,MSM,IMKLSM)
*
              DO MGAS = 1, NGAS
              DO KGAS = 1, NGAS
              DO LGAS = 1, NGAS
                NM = NOBPTS(MGAS,MSM)
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)

*. Blocks of density matrix and integrals : (K L ! I M),D2(K L, J M)
                IXCHNG = 0
                ICOUL  = 1
C                CALL LGETINT(WORK(KLINT),
C    &               ISM,IGAS,MSM,MGAS,KSM,KGAS,LSM,LGAS,
C    &               IXCHNG,0,0,ICOUL)
                CALL LGETINT(WORK(KLINT),
     &               KGAS,KSM,LGAS,LSM,IGAS,ISM,MGAS,MSM,
     &               IXCHNG,0,0,ICOUL)

                CALL GETD2 (WORK(KLDEN),
     &               KSM,KGAS,LSM,LGAS,JSM,JGAS,MSM,MGAS,ICOUL)
C               GETINT(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
C    &                     LTYP,LSM,IXCHNG,0,0,ICOUL)
                NKL = NK*NL
                DO M = 1, NM
                  IIOFF = KLINT + (M-1)*NKL*NI
                  IDOFF = KLDEN + (M-1)*NKL*NJ
                  CALL MATML7(WORK(KLFBLK),WORK(IIOFF),WORK(IDOFF),
     &                        NI,NJ,NKL,NI,NKL,NJ,ONE,ONE,1)
                END DO
              END DO
              END DO
              END DO
            END DO
            END DO
            END IF
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' One- + two-electron contributions'
              WRITE(6,*) ' ================================='
              CALL WRTMT_LU(WORK(KLFBLK),NI,NJ,NI,NJ)
            END IF
*. Block has been constructed , transfer to -complete-
*. symmetry blocked Fock matrix
            DO J = 1, NJ
              DO I = 1, NI
C?              WRITE(6,*) 'IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1',
C?   &                      IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1
                F(IFOFF-1+(J+IJ-1-1)*NIJS + I+II-1 ) =
     &          WORK(KLFBLK-1+(J-1)*NI+I)
              END DO
            END DO
*
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from FOO '
        WRITE(6,*) ' ================'
        CALL APRBLM2(F,NOCOBS,NOCOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'FOO   ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PTFOCK(LU0,LUN,N,ISM,ISPC)
*
* Perturbation expansion of general fock matrix
* through arbitrary order.
*
* It is assumed that this calculation is preceded  by
* a call to the perturbation routine to obtain the
* wave function corrections to the neutral state.
*
* Input
*       LUN : File containing wave function corrections
*       LU0 : File containing reference wave funcrtion
*         N : Max order of expansion
*      ISM : Symmetry of reference state
*      ISPC : Space of referencestate
*
* Jeppe, June 98
*
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 KLVEC1, KLVEC2, KLDEN1, KLDEN1P, KLDEN2, KLFN, KLFSCR,
     &          KLFNS, KLDNS, KLSCR4, KLSCR5, KLSCR6
!               for addressing of WORK
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "lucinp.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      NTEST = 500
*
      WRITE(6,*)
      WRITE(6,*) ' ======================== '
      WRITE(6,*) ' PTFOCK is now in CONTROL '
      WRITE(6,*) ' ======================== '
      WRITE(6,*)

* a bit on files :
* LUSC36 is LUN.
* Two additional scratch files to be used are  LUSC1 and LUSC2
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PTFOCK')
*
*
*     ========================
* 1 : Local memory allocation
*     ========================
*
*. Allocate space for two vector chunks
      CALL MEMMAN(KLVEC1,LBLOCK,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,LBLOCK,'ADDL  ',2,'VEC2  ')
* one-body Density matrices through order n
      NMAT = N+1
      LENGTH = NMAT * NTOOB ** 2
      CALL MEMMAN(KLDEN1,LENGTH,'ADDL  ',2,'DENN1 ')
*. in symmetry packed form
      CALL MEMMAN(KLDEN1P,LENGTH,'ADDL  ',2,'DENN1P')
*. and Two-body density matrices
       LENGTH = NMAT * NTOOB ** 2*(NTOOB**2+1)/2
      CALL MEMMAN(KLDEN2,LENGTH,'ADDL  ',2,'DENN2 ')
*. Fock matrices through order n
      LENGTH = NMAT * 2*NINT1
      CALL MEMMAN(KLFN,LENGTH,'ADDL  ',2,'F(N)  ')
*. A scratch matrix ( not a nice thing to say about a matrix )
      LENGTH =  2*NINT1
      CALL MEMMAN(KLFSCR,LENGTH,'ADDL  ',2,'FSCR  ')
*. Space for Fock matrices, density matrix belonging to a
*  given symmetry
      LENGTH = NMAT * 2*NINT1
      CALL MEMMAN(KLFNS,LENGTH,'ADDL  ',2,'F(N)S ')
      LENGTH = NMAT * 2*NINT1
      CALL MEMMAN(KLDNS,LENGTH,'ADDL  ',2,'D(N)S ')
*.
      CALL MEMMAN(KLSCR4,2*NTOOB**2,'ADDL  ',2,'SCR4  ')
      CALL MEMMAN(KLSCR5,2*NTOOB**2,'ADDL  ',2,'SCR5  ')
      CALL MEMMAN(KLSCR6,2*NTOOB**2,'ADDL  ',2,'SCR6  ')
*
* ===============================================
* 2 : Construct density matrices through order N
* ===============================================
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
*. No print in density matrices
      IPRDEN_SAVE = IPRDEN
      IPRDEN = 0
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLSCR4),ZERO,LRHO1)
      I12_SAVE = I12
      I12 = 2
      DO K = 0, N
        CALL PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
     &       WORK(KLDEN1+(K-0)*LRHO1),
     &       WORK(KLDEN2+(K-0)*LRHO2),LUSC1,LUSC2)
      END DO
      IPRDEN = IPRDEN_SAVE
      I12 = I12_SAVE
*
      WRITE(6,*) ' Memtest 1 : '
      CALL LMEMCHK('PTFOCK')
*
* =====================================================
* 3 : Construct for all orders Fock matrix with unpartioned Hamiltonian
* =====================================================
*
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB **2*(NTOOB**2+1)/2
      LFOCK = 2*NINT1
*
      ONE = 1.0D0
      ONEM = -1.0D0
      DO K = 0, N
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLFN+(K-0)*LFOCK),ZERO,LFOCK)
*. Full Hamiltonian with K order density
        CALL COPVEC(WORK(KLDEN1+(K-0)*LRHO1),WORK(KRHO1),LRHO1)
        CALL COPVEC(WORK(KLDEN2+(K-0)*LRHO2),WORK(KRHO2),LRHO2)
        CALL FOCK_MAT(WORK(KLFSCR),2)
        CALL COPVEC(WORK(KLFSCR),WORK(KLFN+(K-0)*LFOCK),LFOCK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 'Correction to Fock matrix of order =',K
          CALL APRBLM2(WORK(KLFN+(K-0)*LFOCK),
     &                 NOCOBS,NOCOBS,NSMOB,0)
        END IF
*
      END DO
*
      IF(NTEST.GE.100) THEN
* Accumulate corrections to Fock matrix
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLSCR4),ZERO,LFOCK)
        ONE = 1.0D0
        DO K = 0, N
          CALL VECSUM(WORK(KLSCR4),WORK(KLSCR4),
     &         WORK(KLFN+(K-0)*LFOCK),ONE,ONE,LFOCK)
        END DO
*
        WRITE(6,*) '  sum(k) Fock(k) '
        WRITE(6,*) ' =============== '
        CALL APRBLM2(WORK(KLSCR4),NOCOBS,NOCOBS,NSMOB,0)
      END IF
*
*.Finito
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'PTFOCK')
*
      RETURN
      END
