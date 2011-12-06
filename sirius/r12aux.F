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
C  /* Deck r12aux */  
      SUBROUTINE R12AUX(WORK,LWORK)
#include <implicit.h>
#include <priunit.h>   
C
      PARAMETER (D0 = 0D0, D1 = 1D0)
      DIMENSION WORK(LWORK)
#include <iratdef.h>
#include <dummy.h>
#include <thrldp.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <gnrinf.h>
#include <cbirea.h>
#include <r12int.h>
#include "incore.h"
      LOGICAL DELEMO, FOUND, LDUMMY, TEMP_DIRFCK
      INTEGER NBASF(8)
#include <memint.h>
C
      CALL QENTER('R12AUX')
      TEMP_DIRFCK = DIRFCK
      DIRFCK = .TRUE.
      AOSAVE = .FALSE.
      IF (.NOT. LMULBS) THEN
         DO ISYM = 1, NSYM
            MBAS1(ISYM) = NBAS(ISYM)
            MBAS2(ISYM) = 0
         END DO
      END IF
C
C     ***** ALLOCATE MEM AND READ MOLECULAR ORBITALS *****
C
      NOCCT_R12 = 0
      DO ISYM = 1, NSYM
         NOCCT_R12 = NOCCT_R12 + NOCC(ISYM) + NRXR12(ISYM)
      END DO
      NCMOTI  = NCMOT
      IF (BOYORB.OR.PIPORB) THEN
         NCMOT = 0
         NORBT = 0
         DO ISYM = 1, NSYM
cwmk
cwmk        Bug fixed (WK/UniKA/17-08-2005).
cwmk        NORB1(ISYM) = NORB1(ISYM) + NOCC(ISYM)
cwmk
            NORBT = NORBT + NORB1(ISYM)
            NCMOT  = NCMOT + (NORB1(ISYM)+NORB2(ISYM))*NBAS(ISYM)
            IF (NRXR12(ISYM) .NE. 0) CALL QUIT
     &      ('Local MP2-R12 not implemented with .R12ORBa. Sorry.')
         END DO
         NNBOYM = MAX(NORB1(1) + NORB2(1), NBAS(1))
         N2BOYX = NNBOYM * NNBOYM
         NNBOYT = NNBOYM * (NNBOYM + 1) / 2
      ELSE
         NNBOYM = NBAST + NOCCT_R12
         N2BOYX = NNBOYM * NNBOYM
         NNBOYT = NNBOYM * (NNBOYM + 1) / 2
      END IF
c      NMOCCT = 2 * NOCCT * NOCCT
      NMORBT = 2 * NORBT * NORBT
      NMOCCT = 2 * NOCCT_R12 * NOCCT_R12
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KCMO , NCMOT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFVAO,     0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDCAO,N2BOYX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVAO,     0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDV  ,     0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAUX ,N2BOYX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,MAX(NNBOYT,N2BOYX),WORK,KFREE,LFREE)
      IF (R12PRP) THEN
         CALL MEMGET('REAL',KSMO,MAX(7*NMORBT,NNBOYT),WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KSMO,MAX(7*NMOCCT,NNBOYT),WORK,KFREE,LFREE)
      ENDIF
      CALL MEMGET('REAL',KFCAO,N2BOYX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAUX,NNBOYT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFQQ,NOCCT_R12**4,WORK,KFREE,LFREE)
      IF (R12PRP) THEN
         CALL MEMGET('REAL',KFQQA,NVIRT*NOCCT**3,WORK,KFREE,LFREE)
c         CALL MEMGET('REAL',KFQQA,NORBT*NOCCT**3,WORK,KFREE,LFREE)
      ENDIF
      CALL MEMGET('WORK',KWRK,LWRK,WORK,KFREE,LFREE)
      NCMOT  = NCMOTI
      IF (BOYORB.OR.PIPORB) THEN
         OPEN(99,FILE='LOCMO',FORM='FORMATTED')
         KL = KCMO
         DO IJ = 1, NBAS(1)*NORB1(1)
            READ(99,'(D30.20)') WORK(KL)
            KL = KL + 1 
         END DO
         CLOSE(99)
      ELSE
         JRDMO = 9
         CALL READMO(WORK(KCMO),JRDMO,WORK(KFREE),LFREE)
      END IF
C
C     ***** INITIALIZE AUXILIARY BASIS SET *****               
C
      ISTC = KAUX
      DO ISYM = 1, NSYM
         NORBI = MBAS1(ISYM)
         NBASI = NBAS(ISYM)
         NAUXI = MBAS2(ISYM)
         NORB2(ISYM) = NAUXI
         IF (NAUXI .NE. 0) THEN
           IJ = ISTC
           DO J = 1, NAUXI
             DO I = 1, NBASI
               IF (I - NORBI .EQ. J) THEN
                 WORK(IJ) = D1
               ELSE
                 WORK(IJ) = D0
               END IF 
               IJ = IJ + 1
             END DO
           END DO
         END IF
         ISTC = ISTC + NAUXI*NBASI  
      END DO
C
C     ***** ORTHOGONALIZE AUXILIARY BASIS SET *****
C
      LAUXBS = .TRUE.
      IF (R12ECO .OR. R12CBS) THEN
          CALL DELMO(WORK(KAUX),WORK(KWRK),LWRK,THROVL,CMAXMO,DELEMO)
          CALL AUXORT(WORK(KCMO),WORK(KAUX),WORK(KWRK),LWRK)
      END IF
      LAUXBS = .TRUE.
      CALL DELMO(WORK(KAUX),WORK(KWRK),LWRK,THROVL,CMAXMO,DELEMO)
      LAUXBS = .FALSE.
C
C     **** CONSTRUCT FOCK MATRIX IN AUXILIARY BASIS ****
C
      CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDCAO),WORK(KDVAO),
     *               WORK(KCMO),DUMMY,WORK(KWRK),LWRK)

      R12TRA = .TRUE.
      NORXR = NORXR .OR. R12HYB
      IF (R12ECO) THEN
         MBSMAX = 6 
      ELSE
         IF (R12NOB .OR. NORXR) THEN
            MBSMAX = 5
         ELSE
            MBSMAX = 6
         END IF
      END IF
      CALL FCKMAO(.TRUE.,EMCMY,WORK(KFCAO),WORK(KFVAO),
     *     WORK(KDCAO),WORK(KDVAO),WORK(KDV),WORK(KCMO),WORK(KWRK),LWRK)

      WRITE(LUPRI,'(/A,I1,A)') ' Computation of exchange matrix'//
     *                         ' done [',MBSMAX,']'
      MBSMAX = 4
      R12TRA = .FALSE.
      CALL DCOPY(N2BASX,WORK(KFCAO),1,WORK(KDCAO),1)
      CALL DGETSP(NBAST,WORK(KDCAO),WORK(KFCAO))
      IF (NSYM .GT. 1) CALL PKSYM1(WORK(KFCAO),WORK(KFCAO),NBAS,NSYM,2)
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXAOX','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,112) (WORK(KFCAO + J), J = 0, NNBAST - 1) 
      CALL GPCLOSE(LUMULB,'KEEP')
      CALL UTHUB(WORK(KFCAO),WORK(KFAUX),WORK(KAUX),WORK(KDCAO),NSYM,
     *           NBAS,NORB2)   

C
C     **** FORM CANONICAL AUXILIARY ORBITALS ****
C
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXBAS','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,111) NSYM
      ISTF = KFAUX
      ISTC = KAUX
      WRITE(LUPRI,'(/A/A)') ' AUXILIARY BASIS LOW EIGENVALUES',
     *                      ' ISYM         EIGENVALUE'
      DO ISYM = 1, NSYM
         NORBI = NORB1(ISYM)
         NAUXI = NORB2(ISYM)
         NBASI = NBAS(ISYM)
         WRITE (LUMULB,111) ISYM,NORBI,NAUXI,NBASI
         IF (NAUXI .NE. 0) THEN         
            IF (R12NOB .OR. NORXR) THEN
               DO I=1,NAUXI
                  WORK(KWRK + I - 1) = 0D0
               END DO
            ELSE
               CALL JACO(WORK(ISTF),WORK(ISTC),NAUXI,NAUXI,NBASI,
     *                   WORK(KWRK),WORK(KWRK+NAUXI))
               DO I=1,NAUXI
                  WORK(KWRK + I - 1) = WORK(ISTF + I*(I+1)/2 - 1)
               END DO
               CALL ORDER(WORK(ISTC),WORK(KWRK),NAUXI,NBASI)
            END IF
            WRITE(LUPRI,402) ISYM, WORK(KWRK)
  402       FORMAT(I4,F20.12)
            IJ = ISTC - 1
            DO I=1,NAUXI
               WRITE(LUMULB,111) I
               WRITE(LUMULB,112) WORK(KWRK + I - 1)
               WRITE(LUMULB,112) (WORK(IJ + J), J = 1, NBASI)  
               IJ = IJ + NBASI
            END DO
         END IF
         ISTF   = ISTF + NAUXI*(NAUXI+1)/2
         ISTC   = ISTC + NAUXI*NBASI  
      END DO
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM FOCK MATRIX IN ORTHOGONAL BASIS ****
C
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      CALL DCOPY(N2BASX,WORK(KFCAO),1,WORK(KSMAT),1)
      ISMO = KFCAO
      ISMO1 = KCMO
      ISMO2 = KAUX
      NNBASF = 0
      DO ISYM = 1, NSYM
         NORBI = NORB1(ISYM)
         NAUXI = NORB2(ISYM)
         NBASI = NBAS(ISYM) 
         CALL DCOPY(NORBI*NBASI,WORK(ISMO1),1,WORK(ISMO),1)
         ISMO = ISMO + NORBI*NBASI
         CALL DCOPY(NAUXI*NBASI,WORK(ISMO2),1,WORK(ISMO),1)
         ISMO = ISMO + NAUXI*NBASI
         ISMO1 = ISMO1 + NORBI*NBASI   
         ISMO2 = ISMO2 + NAUXI*NBASI
         NBASF(ISYM) = NORBI + NAUXI
         NNBASF = NNBASF + NBASF(ISYM)*(NBASF(ISYM)+1)/2
      END DO
      CALL UTHUB(WORK(KSMAT),WORK(KSMO),WORK(KFCAO),WORK(KDCAO),NSYM,
     *           NBAS,NBASF)    
      WRITE(LUMULB,112) (WORK(KSMO + J), J = 0, NNBASF - 1) 
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM OVERLAP IN ORTHOGONAL BASIS ****
C
      CALL GPOPEN(LUMULB,'AUXOVL','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      FOUND = .TRUE.
      CALL RDONEL('OVERLAP ',FOUND,WORK(KSMAT),NNBAST)
      CALL UTHUB(WORK(KSMAT),WORK(KSMO),WORK(KFCAO),WORK(KDCAO),NSYM,
     *           NBAS,NBASF)    
      WRITE(LUMULB,112) (WORK(KSMO + J), J = 0, NNBASF - 1) 
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM <ij|r12**2|kl> IN ORTHOGONAL BASIS ****
C
      ISMO = KFCAO
      ISMO1 = KCMO
      CALL DZERO(WORK(ISMO),NOCCT_R12*NBAST)
      DO ISYM = 1, NSYM
c         NOCCI = NOCC(ISYM)
         NOCCI = NOCC(ISYM) + NRXR12(ISYM)
         NORBI = NORB1(ISYM)
         NBASI = NBAS(ISYM) 
         NVIRI = NVIR(ISYM) 
         DO INOC = 1, NOCCI 
           CALL DCOPY(NBASI,WORK(ISMO1),1,WORK(ISMO),1)
           ISMO = ISMO + NBAST
           ISMO1 = ISMO1 + NBASI
         END DO
         ISMO1 = ISMO1 + (NORBI - NOCCI) * NBASI
         ISMO = ISMO + NBASI
      END DO 
      ISMO  = KFCAO + NOCCT_R12*NBAST
      ISMO2 = KCMO  + (NOCC(1)+NRXR12(1))*NBAS(1)
      CALL DZERO(WORK(ISMO),NVIRT*NBAST)
      DO ISYM = 1, NSYM
         NOCCI = NOCC(ISYM) + NRXR12(ISYM)
         NORBI = NORB1(ISYM)
         NBASI = NBAS(ISYM) 
         NVIRI = NVIR(ISYM)
         DO INOC = 1, NVIRI 
           CALL DCOPY(NBASI,WORK(ISMO2),1,WORK(ISMO),1)
           ISMO = ISMO + NBAST
           ISMO2 = ISMO2 + NBASI
         END DO
         ISMO = ISMO + NBASI
         IF (ISYM .LT. NSYM)
     &   ISMO2 = ISMO2 + (NOCC(ISYM+1) + NRXR12(ISYM+1)) * NBAS(ISYM+1)
      END DO 
      IF (R12PRP) THEN 
         CALL R12QV(.FALSE.,.FALSE.,WORK(KFQQA),WORK(KSMO+6*NMORBT),
     *              WORK(KSMO+5*NMORBT),
     *              WORK(KSMO+4*NMORBT),WORK(KSMO+3*NMORBT),
     *              WORK(KSMO+2*NMORBT),WORK(KSMO+NMORBT),WORK(KSMO), 
     *              WORK(KFCAO),WORK(KFCAO),WORK(KWRK),LWRK,
     *              NBAST,NORBT,NOCCT_R12,NVIRT,.FALSE.)
         LUNIT = -1
         CALL GPOPEN(LUNIT,'AUXQA12','UNKNOWN','SEQUENTIAL',
     &                      'FORMATTED',IDUMMY,LDUMMY)
         WRITE(LUNIT,112) (WORK(KFQQA + J), J = 0, NVIRT*NOCCT**3 - 1)
         CALL GPCLOSE(LUNIT,'KEEP')
      ENDIF

      CALL R12QQ(WORK(KFQQ),WORK(KSMO+6*NMOCCT),WORK(KSMO+5*NMOCCT),
     *           WORK(KSMO+4*NMOCCT),WORK(KSMO+3*NMOCCT),
     *           WORK(KSMO+2*NMOCCT),WORK(KSMO+NMOCCT),WORK(KSMO), 
     *           WORK(KFCAO),WORK(KFCAO),WORK(KWRK),LWRK,
     *           NBAST,NOCCT_R12)
      CALL GPOPEN(LUMULB,'AUXQ12','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,112) (WORK(KFQQ + J), J = 0, NOCCT_R12**4 - 1)
      CALL GPCLOSE(LUMULB,'KEEP')
      CALL MEMREL('R12AUX',WORK,1,1,KFREE,LFREE)
C
      DIRFCK = TEMP_DIRFCK
      CALL QEXIT('R12AUX')
      RETURN
  111 FORMAT(10I5)
  112 FORMAT(4E30.20)
      END
C  /* Deck auxort */
      SUBROUTINE AUXORT(CMO,AUX,SCRA,LSCRA)
#include <implicit.h>
#include <priunit.h>
      DIMENSION CMO(*), AUX(*), SCRA(*)
      PARAMETER (D1 = 1.0D0)
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <r12int.h>
      KOVLP = 1
      KSMOS = KOVLP + NNBAST
      KSCR1 = KSMOS
      DO ISYM = 1,NSYM
        KSCR1 = KSCR1 + NBAS(ISYM) * (NORB1(ISYM) + NORB2(ISYM))
      END DO
      KSCR2 = KSCR1 + NBAST 
      IF (KSCR2 .GT. LSCRA) CALL ERRWRK('AUXORT',KSCR2+NBASMA,LSCRA)
      CALL RDONEL('OVERLAP ',.TRUE.,SCRA(KOVLP),NNBAST)
      ISMO = KSMOS
      ISMO1 = 1
      ISMO2 = 1
      DO ISYM = 1,NSYM
        ISSYM = KOVLP + IIBAS(ISYM)
        NBASI = NBAS(ISYM)
        NORBI = NORB1(ISYM)
        NAUXI = NORB2(ISYM)
        NORBN = NORBI + NAUXI
        IF (NORBN .NE. 0) THEN        
          CALL DCOPY(NORBI*NBASI,CMO(ISMO1),1,SCRA(ISMO),1)
          CALL DCOPY(NAUXI*NBASI,AUX(ISMO2),1,SCRA(ISMO+NORBI*NBASI),1)
          CALL AUXNRM(SCRA(ISSYM),SCRA(ISMO),NBASI,
     *                NORBI,NAUXI,SCRA(KSCR1))
          CALL DCOPY(NAUXI*NBASI,SCRA(ISMO+NORBI*NBASI),1,AUX(ISMO2),1)
C         CALL AROUND('Auxiliary basis with CMOs projected out')
C         CALL OUTPUT(AUX(ISMO2),1,NAUXI,1,NBASI,NAUXI,NBASI,1,LUPRI)
          ISMO1 = ISMO1 + NORBI*NBASI
          ISMO2 = ISMO2 + NAUXI*NBASI  
          ISMO = ISMO + (NORBI + NAUXI) * NBASI
        END IF
      END DO
      RETURN
      END
C  /* Deck auxnrm */
      SUBROUTINE AUXNRM(S,VC,N,M1,M2,W)                
C
C     PROJECT OUT M1 VECTORS
C
#include <implicit.h>
      DIMENSION S(*), VC(*), W(*)
      PARAMETER ( D0 = 0D0, D1 = 1D0)
      IF (M1. LE. 0 .OR. M2 .LE. 0) RETURN
      IVCI = M1 * N + 1
      DO I = 1, M2
        CALL MPAPV(N,S,VC(IVCI),W)
        IVCJ = 1
        DO J = 1, M1
          T = DDOT(N,VC(IVCJ),1,W,1)
          DO K = 0, N - 1
             VC(IVCI + K) = VC(IVCI + K) - T * VC(IVCJ + K)  
          ENDDO
          IVCJ = IVCJ + N
        ENDDO
        IVCI = IVCI + N
      ENDDO   
      RETURN
      END
C
C  /* Deck r12aux1 */  
      SUBROUTINE R12AUX1(CMO,WORK,LWORK)
C
C     The matrix <ij|r12**2|kl> must be recomputed if canonical virtual
C     orbitals have been replaced by natural orbitals (WK/UniKA/22-11-2005).
C
#include <implicit.h>
#include <priunit.h>   
      LOGICAL FOUND
      DIMENSION WORK(LWORK), CMO(*), NBASF(8)
#include <iratdef.h>
#include <dummy.h>
#include <thrldp.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <gnrinf.h>
#include <r12int.h>
#include <memint.h>
C
      CALL QENTER('R12AUX1')
C
C     ***** ALLOCATE MEMORY *****
C
      KFREE = 1
      LFREE = LWORK
      NMOCCT = 2 * NOCCT * NOCCT
      KBNOB = (NBAST + NOCCT) ** 2
      CALL MEMGET('REAL',KFMO,NOCCT*NBAST,WORK,KFREE,LFREE)
      LENKSM = MAX(7*NMOCCT,KBNOB)
      CALL MEMGET('REAL',KSMO,LENKSM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFQQ,NOCCT**4,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAO,KBNOB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDCAO,KBNOB,WORK,KFREE,LFREE)
      CALL MEMGET('WORK',KWRK,LWRK,WORK,KFREE,LFREE)
C
C     **** COPY MOLECULAR ORBITALS ****
C
      IFMO = KFMO
      ICMA = 1     
      CALL DZERO(WORK(IFMO),NOCCT*NBAST)
      DO ISYM = 1, NSYM
         NOCCI = NOCC(ISYM)
         NORBI = NORB1(ISYM) + NORB2(ISYM)
         NBASI = NBAS(ISYM) 
         DO I = 1, NOCCI 
            CALL DCOPY(NBASI,CMO(ICMA),1,WORK(IFMO),1)
            IFMO = IFMO + NBAST
            ICMA = ICMA + NBASI
         END DO
         ICMA = ICMA + NORBI * NBASI
         IFMO = IFMO + NBASI
      END DO 
C
C     **** FORM <ij|r12**2|kl> IN ORTHOGONAL BASIS ****
C
      KSMO1 = KSMO  + NMOCCT
      KSMO2 = KSMO1 + NMOCCT
      KSMO3 = KSMO2 + NMOCCT
      KSMO4 = KSMO3 + NMOCCT
      KSMO5 = KSMO4 + NMOCCT
      KSMO6 = KSMO5 + NMOCCT
      CALL R12QQ(WORK(KFQQ),WORK(KSMO6),WORK(KSMO5),
     *           WORK(KSMO4),WORK(KSMO3),
     *           WORK(KSMO2),WORK(KSMO1),WORK(KSMO), 
     *           WORK(KFMO),WORK(KFMO),WORK(KWRK),LWRK,
     *           NBAST,NOCCT)
      LUMULB = -1
      CALL GPOPEN(LUMULB,'AUXQ12','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,'(4E30.20)') (WORK(KFQQ + J), J = 0, NOCCT**4 - 1)
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** CONSTRUCT EXCHANGE MATRIX IN AUXILIARY BASIS ****
C
      LUMULB = -1
      CALL GPOPEN(LUMULB,'AUXAOX','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      READ(LUMULB,112) (WORK(KFCAO + J), J = 0, NNBAST - 1)
      CALL GPCLOSE(LUMULB,'KEEP')
      NNBASF = 0
      DO ISYM = 1, NSYM
         NORBI = NORB1(ISYM)
         NAUXI = NORB2(ISYM)
         NBASF(ISYM) = NORBI + NAUXI + NOCC(ISYM)
         NNBASF = NNBASF + NBASF(ISYM)*(NBASF(ISYM)+1)/2
      END DO
      CALL UTHUB(WORK(KFCAO),WORK(KSMO),CMO,WORK(KDCAO),NSYM,
     *           NBAS,NBASF)
      LUMULB = -1
      CALL GPOPEN(LUMULB,'AUXFCKN','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,112) (WORK(KSMO + J), J = 0, NNBASF - 1)
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM OVERLAP IN ORTHOGONAL BASIS ****
C
      LUMULB = -1
      FOUND = .TRUE.
      CALL GPOPEN(LUMULB,'AUXOVLN','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      CALL RDONEL('OVERLAP ',FOUND,WORK(KFCAO),NNBAST)
      CALL UTHUB(WORK(KFCAO),WORK(KSMO),CMO,WORK(KDCAO),NSYM,
     *           NBAS,NBASF)    
      WRITE(LUMULB,112) (WORK(KSMO + J), J = 0, NNBASF - 1) 
      CALL GPCLOSE(LUMULB,'KEEP')
C
  112 FORMAT(4E30.20)
      CALL MEMREL('R12AUX1',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('R12AUX1')
      RETURN
      END
!  --- end of r12aux.F ---
