*=====================================================================*
      SUBROUTINE LAMMATS(XLAMDP,XLAMDH,T1AM,ISYMT,NOT1AM,LRES,
     &                   NGLMDS,IGLMRHS,IGLMVIS,ICMO,WORK,LWORK)
*=====================================================================*
C
C     PURPOSE:
C             Calculate transformation matrices which include
C             all (i.e. frozen + active) orbitals
C             When response calculation (lres) then 
C             calculate Lambda bar 
C
C             NOT1AM - assume T1 amplitudes are zero
C                      --> just order CMO into Lambda matrices
C
C     C. Haettig, spring 2004
C     added virtual blocks, Christof Haettig, spring 2005
*----------------------------------------------------------------------*
       implicit none
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <r12int.h>
C
      LOGICAL LOCDBG, LRES, NOT1AM
      PARAMETER ( LOCDBG = .FALSE. ) 
C
      INTEGER NGLMDS(8), IGLMRHS(8,8), IGLMVIS(8,8),
     &        LWORK, KCMO, LWRK1, ISYM, ISYMI, ISYMA, ISYMJ, KOFF6,
     &        ISYMP, ISYMB, NBASP, NVIRB, ISYMT, NVIRA, KOFF4, KOFF5,
     &        KOFF1, KOFF2, KOFF3, KEND, ICOUNT, ICOUNT2, ISYM1, ISYM2
      INTEGER NCMO(8), ICMO(8,8)
#if defined (SYS_CRAY)
      REAL XLAMDH(*),XLAMDP(*),WORK(*),T1AM(*),ONE,CMO(NLAMDS)
#else
      DOUBLE PRECISION  XLAMDH(*),XLAMDP(*),WORK(*),T1AM(*),ONE,
     &                  CMO(NLAMDS)
#endif
      PARAMETER (ONE = 1.0d0)
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB(LABEL,LUSIFC,LUPRI)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (CMO(I), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
      
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'CMOS out of LAMMATS:'
        DO ISYM = 1, NSYM
          ISYMI = MULD2H(ISYM,ISYMT)
          CALL OUTPUT(CMO,1,NBAS(ISYM),1,NRHFS(ISYMI),
     &                NBAS(ISYM),NRHFS(ISYMI),1,LUPRI)
        END DO
      END IF
C
C---------------------------------------
C     Reorder the MO-coefficient matrix.
C---------------------------------------
C
      IF (LRES) THEN
        CALL DZERO(XLAMDH,NGLMDS(ISYMT)) 
        CALL DZERO(XLAMDP,NGLMDS(ISYMT)) 
      ELSE
        DO ISYM = 1,NSYM
         KOFF1  = ICMO(ISYM,ISYM) + 1
         KOFF2  = IGLMRHS(ISYM,ISYM) + 1
         CALL DCOPY(NBAS(ISYM)*NRHFS(ISYM),CMO(KOFF1),1,XLAMDH(KOFF2),1)

         KOFF1  = ICMO(ISYM,ISYM) + NBAS(ISYM)*NRHFS(ISYM) + 1
         KOFF2  = IGLMVIS(ISYM,ISYM) + 1
         CALL DCOPY(NBAS(ISYM)*NVIRS(ISYM),CMO(KOFF1),1,XLAMDH(KOFF2),1)
        END DO
        CALL DCOPY(NGLMDS(1),XLAMDH,1,XLAMDP,1)
      END IF
C
C-------------------------------------------
C     Calculate the transformation matrices.
C-------------------------------------------
C
      IF (.NOT. NOT1AM) THEN
       DO ISYMP = 1,NSYM
 
         ISYMI = MULD2H(ISYMP,ISYMT)
         ISYMB = MULD2H(ISYMI,ISYMT)
 
         NBASP = MAX(NBAS(ISYMP),1)
         NVIRB = MAX(NVIR(ISYMB),1)
 
         KOFF1 = ICMO(ISYMP,ISYMB) + NBAS(ISYMP)*NRHFS(ISYMB) + 1
         KOFF2 = IT1AM(ISYMB,ISYMI) + 1
         KOFF3 = IGLMRHS(ISYMP,ISYMI) + NBAS(ISYMP)*NRHFFR(ISYMI) + 1
 
         CALL DGEMM('N','N',NBAS(ISYMP),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,CMO(KOFF1),NBASP,T1AM(KOFF2),NVIRB,
     *              ONE,XLAMDH(KOFF3),NBASP)
 
         ISYMA = MULD2H(ISYMP,ISYMT)
         ISYMJ = MULD2H(ISYMA,ISYMT)
 
         NBASP = MAX(NBAS(ISYMP),1)
         NVIRA = MAX(NVIR(ISYMA),1)
 
         KOFF4 = ICMO(ISYMP,ISYMJ) + NBAS(ISYMP)*NRHFFR(ISYMJ) + 1
         KOFF5 = IT1AM(ISYMA,ISYMJ) + 1
         KOFF6 = IGLMVIS(ISYMP,ISYMA) + 1
 
         CALL DGEMM('N','T',NBAS(ISYMP),NVIR(ISYMA),NRHF(ISYMJ),
     *              -ONE,CMO(KOFF4),NBASP,T1AM(KOFF5),NVIRA,
     *              ONE,XLAMDP(KOFF6),NBASP)
       END DO
      END IF
C
C-------------------------------------------
C     Print the matrices:
C-------------------------------------------
C
      IF (LOCDBG) THEN
C
         IF (LRES) THEN
           CALL AROUND('Lambda Particle bar matrix in LAMMATS')
         ELSE
           CALL AROUND('Lambda Particle matrix in LAMMATS')
         END IF

         DO ISYM = 1,NSYM
            ISYMI = MULD2H(ISYM,ISYMT)
            WRITE(LUPRI,1) ISYM,ISYMI
            WRITE(LUPRI,2)
            WRITE(LUPRI,3)
            IF (NRHF(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,4)
            ELSE
              KOFF1 = 1 + IGLMRHS(ISYM,ISYMI)
              CALL OUTPUT(XLAMDP(KOFF1),1,NBAS(ISYM),1,NRHFS(ISYMI),
     *                    NBAS(ISYM),NRHFS(ISYMI),1,LUPRI)
            END IF
         END DO
C
         IF (LRES) THEN
           CALL AROUND('Lambda Hole bar matrix in LAMMATS')
         ELSE
           CALL AROUND('Lambda Hole matrix in LAMMATS')
         END IF

         DO ISYM = 1,NSYM
            ISYMI = MULD2H(ISYM,ISYMT)
            WRITE(LUPRI,1) ISYM, ISYMI
            WRITE(LUPRI,7)
            WRITE(LUPRI,8)
            IF (NRHF(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,4)
            ELSE
              KOFF1 = 1 + IGLMRHS(ISYM,ISYMI)
              CALL OUTPUT(XLAMDH(KOFF1),1,NBAS(ISYM),1,NRHFS(ISYMI),
     *                    NBAS(ISYM),NRHFS(ISYMI),1,LUPRI)
            END IF
          END DO
C
      END IF
C
      RETURN
C
    1 FORMAT(/,/,7X,'Symmetry number :',2I5)
    2 FORMAT(/,/,7X,'Lambda particle occupied part')
    3 FORMAT(7X,'-----------------------------')
    4 FORMAT(/,/,7X,'This symmetry is empty')
    5 FORMAT(/,/,7X,'Lambda particle virtual part')
    6 FORMAT(7X,'----------------------------')
    7 FORMAT(/,/,7X,'Lambda hole occupied part')
    8 FORMAT(7X,'-------------------------')
    9 FORMAT(/,/,7X,'Lambda hole virtual part')
   10 FORMAT(7X,'------------------------')
C
      END
*======================================================================*
