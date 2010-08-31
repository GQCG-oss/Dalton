C
C  /* Deck rp_charge */
      SUBROUTINE RP_CHARGE(ISYMTR,NEXCI,SECMAT,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1998
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate RPA charge radius (second moment of charge)
C              of ground and excited states and molecular orbitals.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION WORK(LWORK)
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <soppinf.h>
#include <maxorb.h>
CSPAS:15/03-2006: merge with Dalton-2.0  
C#include <cbilrs.h>
#include <cbiexc.h>
CKeinSPASmehr
#include <secmom.h>
C
      DIMENSION SECMAT(3,MXNEXI,NSYM)
      CHARACTER*8 LABEL,RTNLBL(2)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_CHARGE')
C
C----------------------------------------------
C     Set symmetry of xx, yy, and zz operators.
C----------------------------------------------
C
      KSYMOP = 1
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LPRP1  = N2BST(KSYMOP)
      LTR1E  = NT1AM(ISYMTR)
      LTR1D  = NT1AM(ISYMTR)
C
      KPRP1   = 1
      KTR1E   = KPRP1  + LPRP1
      KTR1D   = KTR1E  + LTR1E
      KAIE    = KTR1D  + LTR1D
      KAID    = KAIE   + LTR1E
      KEND1   = KAID   + LTR1D
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_CHARGE',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_CHARGE',' ',KEND1,LWORK)
C
C-------------------------------------
C     Loop over cartesian coordinates.
C-------------------------------------
C
      DO 400 ICAR = 1, 3
C
         IF (ICAR .EQ. 1) LABEL = 'XXSECMOM'
         IF (ICAR .EQ. 2) LABEL = 'YYSECMOM'
         IF (ICAR .EQ. 3) LABEL = 'ZZSECMOM'
C
C-------------------------------------------
C        Get MO one-electron propety matrix.
C-------------------------------------------
C
         CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL,KSYMOP,RTNLBL,
     &                  WORK(KEND1),LWORK1)
C
C-------------------------------------------------------------
C        Extract second moments of charge for each orbital and 
C        the SCF ground state.
C-------------------------------------------------------------
C
         GROUND = ZERO
C
         DO 50 ISYMI = 1,NSYM
C
            KOFF1 = KPRP1 + IFCRHF(ISYMI,ISYMI) 
C
            DO 40  I = 1,NRHF(ISYMI)
C
               KOFF2 = KOFF1 + I + (I-1) * NORB(ISYMI) - 1
C
               SECOMO (ICAR,I,ISYMI) = WORK(KOFF2)
C
               GROUND = GROUND + WORK(KOFF2) * TWO
C
   40       CONTINUE
C
            ISYMA = ISYMI
            KOFF3 = KPRP1 + IFCVIR(ISYMA,ISYMA) + NRHF(ISYMA)
C
            DO 41  A = 1,NVIR(ISYMA)
C
               KOFF4 = KOFF3 + A + (A-1) * NORB(ISYMA) - 1
C
               SECVMO (ICAR,A,ISYMA) = WORK(KOFF4)
C
   41       CONTINUE
C
   50    CONTINUE
C
         SECGR (ICAR) = GROUND
C
C-----------------------------------------
C        Open files with solution vectors.
C-----------------------------------------
C
         CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
         CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
C
C--------------------------------------------------------
C        Loop over excitations of the specified symmetry.
C--------------------------------------------------------
C
         DO 300 IEXCI = 1,NEXCI
C
C---------------------------------
C           Read solution vectors.
C---------------------------------
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,IEXCI)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,IEXCI)
C
C------------------------------------------------------------------
C           Calculate excited state minus ground state one-electron
C           property.
C------------------------------------------------------------------
C
            DO 100 ISYMI = 1, NSYM
C
               ISYMJ = ISYMI
               ISYMA = MULD2H(ISYMJ,ISYMTR)
C
               NVIRA = MAX(NVIR(ISYMA),1)
               NORBJ = MAX(NORB(ISYMJ),1)
C
               KOFF1 = KTR1E + IT1AM(ISYMA,ISYMJ)
               KOFF2 = KTR1D + IT1AM(ISYMA,ISYMJ)
               KOFF3 = KPRP1 + IFCRHF(ISYMJ,ISYMI) 
               KOFF4 = KAIE  + IT1AM(ISYMA,ISYMI)
               KOFF5 = KAID  + IT1AM(ISYMA,ISYMI)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                    ONE,WORK(KOFF1),NVIRA,WORK(KOFF3),NORBJ,ZERO,
     &                    WORK(KOFF4),NVIRA)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                    ONE,WORK(KOFF2),NVIRA,WORK(KOFF3),NORBJ,ZERO,
     &                    WORK(KOFF5),NVIRA)
C
  100       CONTINUE
C
            DOTP1 = DDOT(LTR1E,WORK(KAIE),1,WORK(KTR1E),1)
            DOTP2 = DDOT(LTR1D,WORK(KAID),1,WORK(KTR1D),1)
C
            DO 200 ISYMA = 1, NSYM
C
               ISYMB = ISYMA
               ISYMI = MULD2H(ISYMB,ISYMTR)
C
               NORBA = MAX(NORB(ISYMA),1)
               NVIRB = MAX(NVIR(ISYMB),1)
               NVIRA = MAX(NVIR(ISYMA),1)
C
               KOFF1 = KPRP1 + IFCVIR(ISYMA,ISYMB) + NRHF(ISYMA)
               KOFF2 = KTR1E + IT1AM(ISYMB,ISYMI)
               KOFF3 = KTR1D + IT1AM(ISYMB,ISYMI)
               KOFF4 = KAIE  + IT1AM(ISYMA,ISYMI)
               KOFF5 = KAID  + IT1AM(ISYMA,ISYMI)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,WORK(KOFF1),NORBA,WORK(KOFF2),NVIRB,ZERO,
     &                    WORK(KOFF4),NVIRA)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,WORK(KOFF1),NORBA,WORK(KOFF3),NVIRB,ZERO,
     &                    WORK(KOFF5),NVIRA)
C
  200       CONTINUE
C
            DOTP3 = DDOT(LTR1E,WORK(KAIE),1,WORK(KTR1E),1)
            DOTP4 = DDOT(LTR1D,WORK(KAID),1,WORK(KTR1D),1)
C
            DOTP = DOTP3 + DOTP4 - DOTP1 - DOTP2
C
            SECMAT(ICAR,IEXCI,ISYMTR) = DOTP
C
  300    CONTINUE
C
C-----------------------------------------
C        Close files with solution vectors.
C-----------------------------------------
C
         CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
         CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
  400 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_CHARGE')
C
      RETURN
      END
