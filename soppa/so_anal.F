C  /* Deck so_anal */
      SUBROUTINE SO_ANAL(TR1E,TR1D,LTR1E,TR2E,TR2D,LTR2E,THR1,THR2,
     &                   ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1997, based on CC_PRAM by Ove Christiansen.
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Analysis of eigenvectors in terms of most important
C              orbital excitations which are printed.
C
#include <implicit.h>
#include <priunit.h>
#include <soppinf.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      PARAMETER (ZERO = 0.0D+00, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION TR1E(LTR1E), TR1D(LTR1E)
      DIMENSION TR2E(LTR2E), TR2D(LTR2E)
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ANAL')
C
C---------------------------
C     Write header of table.
C---------------------------
C
      WRITE(LUPRI,'(A)')
     *     ' |=====================================================|'
      WRITE(LUPRI,'(1X,A)')
     *     '|          | Symmetry|  Orbital index  '
     *     //'|   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|          | index   |   a   b   i   j '
     *     //'|              |'
      WRITE(LUPRI,'(A)')
     *     ' |=====================================================|'
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
C---------------------------------------
C     Loop through one excitation part.
C---------------------------------------
C
  1   CONTINUE
      N1 = 0
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            MI = IORB(ISYMI) + I
C
            DO 120 A=1,NVIR(ISYMA)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
               IF (ABS(TR1E(NAI)) .GT. THR1 ) THEN
C
                  WRITE(LUPRI,9990) ISYMA,ISYMI,A,I,TR1E(NAI)
C
                  N1 = N1 + 1
cKeld                  SUMOFP = SUMOFP + TR1E(NAI)*TR1E(NAI)
C
               ENDIF
C
               IF (ABS(TR1D(NAI)) .GT. THR1 ) THEN
C
                  WRITE(LUPRI,9991) ISYMA,ISYMI,A,I,TR1D(NAI)
C
                  N1 = N1 + 1
cKeld                  SUMOFP = SUMOFP + TR1E(NAI)*TR1E(NAI)
C
               ENDIF
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      IF ((N1 .LT. 1) .AND. (THR1 .GT. 1.0D-4) ) THEN
         THR1 = THR1/5.0D0
         GOTO 1
      ENDIF
C
      WRITE(LUPRI,'(A)')
     *     ' |-----------------------------------------------------|'
C
C-------------------------------------------
C     Loop through Double excitation vector.
C-------------------------------------------
C
 2    CONTINUE
C
      N2 = 0
C
      DO 200 ISYMAI = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 240 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 250 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 260 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ)).OR.(ISYMAI.LT.ISYMBJ))
     *                          GOTO 260
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           IF (ABS(TR2E(NAIBJ)) .GT. THR2 ) THEN
C
                              WRITE(LUPRI,9992) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,TR2E(NAIBJ)
                              N2 = N2 + 1
C
cKeld                            SUMOFP = SUMOFP + CAM(NAIBJ)*CAM(NAIBJ)
C
                           ENDIF
C
                           IF (ABS(TR2D(NAIBJ)) .GT. THR2 ) THEN
C
                              WRITE(LUPRI,9993) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,TR2D(NAIBJ)
                              N2 = N2 + 1
C
cKeld                           SUMOFP = SUMOFP + CAM(NAIBJ)*CAM(NAIBJ)
C
                           ENDIF
C
  260                   CONTINUE
C
  250                CONTINUE
C
  240             CONTINUE
C
  230          CONTINUE
C
  220       CONTINUE
C
  210    CONTINUE
C
  200 CONTINUE
C
      IF ((N2 .LT. 1) .AND. (THR2 .GT. 1.0D-4) ) THEN
         THR2 = THR2/5D00
         GOTO 2
      ENDIF
C
      WRITE(LUPRI,'(A)')
     *     ' |=====================================================|'
C
cKeld      WRITE(LUPRI,'(//10X,A,8X,F10.4)')
cKeld     *     'Norm of Printed Amplitude Vector : ',SQRT(SUMOFP)
      WRITE(LUPRI,'(/1X,A53,1X,F8.6)')
     *  'Printed all single excitation amplitudes greater than',THR1
      WRITE(LUPRI,'(1X,A53,1X,F8.6)')
     *  'Printed all double excitation amplitudes greater than',THR2
C
 9990 FORMAT(1X,'|    exci. | ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',
     *       1x, F10.6,'  |')
 9991 FORMAT(1X,'| de-exci. | ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',
     *       1x, F10.6,'  |')
 9992 FORMAT(1X,'|    exci. | ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',1x, F10.6,'  |')
 9993 FORMAT(1X,'| de-exci. | ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',1x, F10.6,'  |')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ANAL')
C
      RETURN
      END
