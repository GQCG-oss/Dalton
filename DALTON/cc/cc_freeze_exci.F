C  /* Deck cc_freeze_exci */
      SUBROUTINE CC_FREEZE_exci(CAM,ISYMTR,
     &           MAXCORE, MAXION,
     &           NRHFCORE,IRHFCORE,NVIRION,IVIRION,
     &           LBOTH)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     09-7-2015 Sonia Coriani
C
C     Purpose: Project out specific excitations 
C              from a trial vector (by zeroing 
C              specific elements)
C     Ex1: zero all ai and aibj elements where i and j
C     are valence orbitals (CORE-VALENCE SEPARATION)
C     Ex2: zero all a an b elements that do not correspond
C     to a specific virtual orbitals
C
C Based on cc_pram()
! CAM is the vector analyzed, of symmetry ISYMTR
! LBOTH checks if both CVS and IONISATION are requested 
! Control is passed via argument list, not via common block
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      Implicit none

#include "ccsdsym.h"
      Double precision CAM(*)
      Integer MAXCORE, NRHFCORE(8),IRHFCORE(MAXCORE,8)
      Integer MAXION,NVIRION(8),IVIRION(MAXION,8)
      integer ISYMTR,ISYMAI,ISYMI,ISYMA,ISYMJ,ISYMB,ISYMBJ
      Double precision TWO, THR1, THR2, zero
      PARAMETER (TWO = 2.0D0,zero=0.0d0)
      Logical LOCDBG, ikeep, LBOTH
      Parameter (Locdbg = .false.)
      Integer AA, II, MA, MI, JJ, BB, NBJ, NAI, MJ, MB
      Integer KAIBJ, NAIBJ, INDEX
C
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "priunit.h"
Cholesky
#include "maxorb.h"
#include "ccdeco.h"
C
      LOGICAL CCSEFF
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CCSEFF = CCS .OR. (CHOINT.AND.CC2)
C
      THR1 = 1.0D-2
      THR2 = 1.0D-2
C
C------------------------------------------
C     Loop through single excitation part.
C------------------------------------------
C
      if (locdbg) then
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
      end if
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
      ikeep=.false.
      DO 100 ISYMA = 1,NSYM
         ISYMI = MULD2H(ISYMAI,ISYMA)
         DO 110 I = 1,NRHF(ISYMI)
            MI = IORB(ISYMI) + I
            DO 120 A=1,NVIR(ISYMA)
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
               ikeep = .false.
               IF (LBOTH) THEN
                 do ii = 1, NRHFCORE(ISYMI)
                  IF (I==IRHFCORE(II,ISYMI)) THEN
                     do aa = 1, NVIRION(ISYMA)
                        IF (A==IVIRION(AA,ISYMA)) THEN
                           ikeep = .true.
                           go to 333
                        END IF
                     end do
                  END IF
                 end do
               ELSE
                 do ii = 1, NRHFCORE(ISYMI)
                   IF (I==IRHFCORE(II,ISYMI)) THEN
                      ikeep = .true.
                      go to 333
                   end if
                 end do
                 do aa = 1, NVIRION(ISYMA)
                   IF (A==IVIRION(AA,ISYMA)) THEN
                    !write(lupri,*)'IONISATION ONLY'
                    ikeep = .true.
                    go to 333
                   END IF
                 end do
               end if
  333          continue
               if (.not.ikeep) CAM(NAI) = zero
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through double excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCSEFF .OR. CCP2 )) THEN
C
      if (locdbg) then
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
      end if
C
      ikeep = .false.
      DO 200 ISYMAI = 1,NSYM
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
         IF (ISYMAI.lt.ISYMBJ) GO TO 200
         DO 210 ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO 220 ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO 230 J = 1,NRHF(ISYMJ)
                  MJ = IORB(ISYMJ) + J
                  DO 240 B = 1,NVIR(ISYMB)
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
                     DO 250 I = 1,NRHF(ISYMI)
                        MI = IORB(ISYMI) + I
                        DO 260 A = 1,NVIR(ISYMA)
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
                           IF ((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ))
     *                          GOTO 260
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSE
                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
         
                           ENDIF
                           KAIBJ = NAIBJ + NT1AM(ISYMTR)  !same vector contains singles and doubles
                                                          !skip the singles
                           ikeep = .false.
                           if (LBOTH) then
                              do ii = 1, nrhfcore(isymi)
                                 if (i==IRHFCORE(II,ISYMI)) then
                                   do aa = 1, NVIRION(ISYMA)
                                     IF (A.EQ.IVIRION(AA,ISYMA)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                   do bb = 1, NVIRION(ISYMB)
                                     IF (B==IVIRION(BB,ISYMB)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                 end if
                              end do
                              do jj = 1, nrhfcore(isymj)
                                 if (j==IRHFCORE(jj,ISYMJ)) then
                                   do aa = 1, NVIRION(ISYMA)
                                     IF (A.EQ.IVIRION(AA,ISYMA)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                   do bb = 1, NVIRION(ISYMB)
                                     IF (B==IVIRION(BB,ISYMB)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                 end if
                              end do
                           else
                             do ii = 1, nrhfcore(isymi)
                              if (i==IRHFCORE(II,ISYMI)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do aa = 1, nvirion(isyma)
                              if (a==IVIRION(aa,ISYMA)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do jj = 1, nrhfcore(isymj)
                              if (j==IRHFCORE(JJ,ISYMJ)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do bb = 1, nvirion(isymb)
                              if (b==IVIRION(bb,ISYMB)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                           end if
!========================================
  444                      continue
                           if (.not.ikeep) CAM(KAIBJ) = zero
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F15.9,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | ',I12,' | ',1x,F15.9,'  |')

      RETURN
      END
!---

C  /* Deck cc_pram0 */
      SUBROUTINE CC_PRAM0(CAM,PT1,ISYMTR,LGRS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! THIS ROUTINE IS BUGGED FOR SYMMETRY. WRONG STORAGE ASSUMED! SONIA
!
C     30-5-1995 Ove Christiansen
C
C     Purpose: Writes out vector:
C              %T1 and %T2 and ||T1||/||T2||
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
C
      PARAMETER (TWO = 2.0D0, THPRT = 1.0D-9, HUNDRED = 100.0D0)
C
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "priunit.h"
Cholesky
#include "maxorb.h"
#include "ccdeco.h"
C
      LOGICAL CCSEFF, LOCDBG
      parameter (locdbg=.false.)
Cholesky
C
C
      LOGICAL LGRS
      DIMENSION CAM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
Cholesky
      CCSEFF = CCS .OR. (CHOINT.AND.CC2)
Cholesky
C
C------------------------
C     Add up the vectors.
C------------------------
C
      C1NOSQ = 0.0D0
      C2NOSQ = 0.0D0
      KC1 = 1
      KC2 = KC1 + NT1AM(ISYMTR)
      C1NOSQ = DDOT(NT1AM(ISYMTR),CAM(KC1),1,CAM(KC1),1)
Chol  IF (.NOT. CCS) C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      IF (.NOT. CCSEFF)
     &   C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      CNOSQ  = C1NOSQ + C2NOSQ
C
      IF (.NOT. (CCSEFF.OR.CCP2) .AND. CNOSQ.NE.0.0D0) THEN
C
         WRITE(LUPRI,'(//10X,A)')
     *     'CC_PRAM:Overall Contribution of the Different Components'
         WRITE(LUPRI,'(10X,A//)')
     *     '--------------------------------------------------------'
         WRITE(LUPRI,'(/10X,A,10X,F10.4,A)')
     *     'Single Excitation Contribution : ',
     *     (C1NOSQ/CNOSQ)*HUNDRED,' %'
         WRITE(LUPRI,'(/10X,A,10X,F10.4,A)')
     *     'Double Excitation Contribution : ',
     *     (C2NOSQ/CNOSQ)*HUNDRED,' %'
         WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *     '||T1||/||T2||                  : ',
     *      SQRT(C1NOSQ/C2NOSQ)
         IF (LGRS) WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *     'Tau1 diagnostic                : ',
     *      SQRT(C1NOSQ/(TWO*DFLOAT(NRHFT)))
         PT1 = (C1NOSQ/CNOSQ)*HUNDRED
      ELSE
         PT1 = HUNDRED
      ENDIF
      WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *  'Norm of Total Amplitude Vector : ',SQRT(CNOSQ)
C
      CALL FLSHFO(LUPRI)
C
C----------------------------------------------
C     Initialize threshold etc from Printlevel.
C----------------------------------------------
C
      NL = MAX(1,2*IPRINT)
C
      CNOSQ = MAX(CNOSQ,THPRT)
C
      THR1 = SQRT(C1NOSQ/CNOSQ)/NL
      THR2 = SQRT(C2NOSQ/CNOSQ)/NL
      THR1 = MAX(THR1,1.0D-02)
      THR2 = MAX(THR2,1.0D-02)
      SUMOFP = 0.0D00
C
      IF (DEBUG.or.locdbg) THR1 = 1.0D-2
      IF (DEBUG.or.locdbg) THR2 = 1.0D-2
C
C-------------------------------------
C     Loop until a few is Printed out.
C-------------------------------------
C
C
C---------------------------------------
C     Loop through One excitation part.
C---------------------------------------
C
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
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
               IF (ABS(CAM(NAI)) .GE. THR1 ) THEN
C
                  WRITE(LUPRI,9990) ISYMA,ISYMI,A,I,NAI,CAM(NAI)
C
                  N1 = N1 + 1
                  SUMOFP = SUMOFP + CAM(NAI)*CAM(NAI)
C
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF ((N1.LT.1).AND.(SQRT(C1NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR1 = THR1/5.0D0
         GOTO 1
      ENDIF
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through Double excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCSEFF .OR. CCP2 )) THEN
C
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
C

 2    CONTINUE
      N2 = 0
C
      DO 200 ISYMAI = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)

         if (isymai.lt.isymbj) go to 200
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

                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                           ENDIF
C
                           KAIBJ = NAIBJ + NT1AM(ISYMTR)

                           IF (ABS(CAM(KAIBJ)) .GT. THR2 ) THEN
C
                              WRITE(LUPRI,9991) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
     *                                      CAM(KAIBJ)
                              N2 = N2 + 1
C
                              SUMOFP = SUMOFP + CAM(KAIBJ)*CAM(KAIBJ)
C
                           ENDIF
C
  260                   CONTINUE  !loop A
  250                CONTINUE     !loop I
  240             CONTINUE        !loop B
  230          CONTINUE           !loop J
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF ((N2.LT.1).AND.(SQRT(C2NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR2 = THR2/5D00
         GOTO 2
      ENDIF
C
      ENDIF
C
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      WRITE(LUPRI,'(//10X,A,8X,F10.4)')
     *     'Norm of Printed Amplitude Vector : ',SQRT(SUMOFP)
      WRITE(LUPRI,'(/10X,A43,F9.6)')
     *     'Printed all single excitations greater than',THR1
      IF (.NOT. (CCSEFF.OR.CCP2)) THEN
         WRITE(LUPRI,'(/10X,A43,F9.6)')
     *     'Printed all double excitations greater than',THR2
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F10.6,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | ',I12,' | ',1x,F10.6,'  |')
C
      RETURN
      END


C  /* Deck cc_freeze_start */
      SUBROUTINE CC_FREEZE_start(CAM,ISYMTR,
     &           MAXCORE, MAXION,
     &           NRHFCORE,IRHFCORE,NVIRION,IVIRION,
     &           LBOTH)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     24/8/2016, Eirik & Sonia
C
C     Purpose: allow automated selection of start vectors
C              for core excitations and ionizations
C              Achieved by setting all elements of the 
C              diagonal of the Fock matrix non refering to
C              the selected core/diffuse orbital to a huge
C              number so that it will be discarded by the 
C              FNDM3 routine later on (the routine picks up
C              the indices of the lowest energy eigenvalue).
C    In other words, we push the valence eigenvalues above the core 
C              specific elements)
C
C Based on cc_pram()
! CAM is the vector analyzed, of symmetry ISYMTR
! LBOTH checks if both CVS and IONISATION are requested 
! Control is passed via argument list, not via common block
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      Implicit none

#include "ccsdsym.h"
      Double precision CAM(*)
      Integer MAXCORE, NRHFCORE(8),IRHFCORE(MAXCORE,8)
      Integer MAXION,NVIRION(8),IVIRION(MAXION,8)
      integer ISYMTR,ISYMAI,ISYMI,ISYMA,ISYMJ,ISYMB,ISYMBJ
      Double precision TWO, THR1, THR2, crazy
      PARAMETER (TWO = 2.0D0,crazy=1.0d+14)
      Logical LOCDBG, ikeep, LBOTH
      Parameter (Locdbg = .false.)
      Integer AA, II, MA, MI, JJ, BB, NBJ, NAI, MJ, MB
      Integer KAIBJ, NAIBJ, INDEX
C
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "priunit.h"
Cholesky
#include "maxorb.h"
#include "ccdeco.h"
C
      LOGICAL CCSEFF
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CCSEFF = CCS .OR. (CHOINT.AND.CC2)
C
      THR1 = 1.0D-2
      THR2 = 1.0D-2
C
C------------------------------------------
C     Loop through single excitation part.
C------------------------------------------
C
      if (locdbg) then
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
      end if
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
      ikeep=.false.
      DO 100 ISYMA = 1,NSYM
         ISYMI = MULD2H(ISYMAI,ISYMA)
         DO 110 I = 1,NRHF(ISYMI)
            MI = IORB(ISYMI) + I
            DO 120 A=1,NVIR(ISYMA)
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
               ikeep = .false.
               IF (LBOTH) THEN
                 do ii = 1, NRHFCORE(ISYMI)
                  IF (I==IRHFCORE(II,ISYMI)) THEN
                     do aa = 1, NVIRION(ISYMA)
                        IF (A==IVIRION(AA,ISYMA)) THEN
                           ikeep = .true.
                           go to 333
                        END IF
                     end do
                  END IF
                 end do
               ELSE
                 do ii = 1, NRHFCORE(ISYMI)
                   IF (I==IRHFCORE(II,ISYMI)) THEN
                      ikeep = .true.
                      go to 333
                   end if
                 end do
                 do aa = 1, NVIRION(ISYMA)
                   IF (A==IVIRION(AA,ISYMA)) THEN
                    ikeep = .true.
                    go to 333
                   END IF
                 end do
               end if
  333          continue
               if (.not.ikeep) CAM(NAI) = crazy
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through double excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCSEFF .OR. CCP2 )) THEN
C
      if (locdbg) then
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
      end if
C
      ikeep = .false.
      DO 200 ISYMAI = 1,NSYM
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
         IF (ISYMAI.lt.ISYMBJ) GO TO 200
         DO 210 ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO 220 ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO 230 J = 1,NRHF(ISYMJ)
                  MJ = IORB(ISYMJ) + J
                  DO 240 B = 1,NVIR(ISYMB)
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
                     DO 250 I = 1,NRHF(ISYMI)
                        MI = IORB(ISYMI) + I
                        DO 260 A = 1,NVIR(ISYMA)
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
                           IF ((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ))
     *                          GOTO 260
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSE
                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
         
                           ENDIF
                           KAIBJ = NAIBJ + NT1AM(ISYMTR)  !same vector contains singles and doubles
                                                          !skip the singles
                           ikeep = .false.
                           if (LBOTH) then
                              do ii = 1, nrhfcore(isymi)
                                 if (i==IRHFCORE(II,ISYMI)) then
                                   do aa = 1, NVIRION(ISYMA)
                                     IF (A.EQ.IVIRION(AA,ISYMA)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                   do bb = 1, NVIRION(ISYMB)
                                     IF (B==IVIRION(BB,ISYMB)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                 end if
                              end do
                              do jj = 1, nrhfcore(isymj)
                                 if (j==IRHFCORE(jj,ISYMJ)) then
                                   do aa = 1, NVIRION(ISYMA)
                                     IF (A.EQ.IVIRION(AA,ISYMA)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                   do bb = 1, NVIRION(ISYMB)
                                     IF (B==IVIRION(BB,ISYMB)) THEN
                                        ikeep = .true.
                                        go to 444
                                     END IF
                                   end do
                                 end if
                              end do
                           else
                             do ii = 1, nrhfcore(isymi)
                              if (i==IRHFCORE(II,ISYMI)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do aa = 1, nvirion(isyma)
                              if (a==IVIRION(aa,ISYMA)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do jj = 1, nrhfcore(isymj)
                              if (j==IRHFCORE(JJ,ISYMJ)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do bb = 1, nvirion(isymb)
                              if (b==IVIRION(bb,ISYMB)) then
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                           end if
!========================================
  444                      continue
                           if (.not.ikeep) CAM(KAIBJ) = crazy
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F15.9,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | ',I12,' | ',1x,F15.9,'  |')

      RETURN
      END
!---
