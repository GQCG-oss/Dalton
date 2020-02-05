C  /* Deck cc_freeze_tripletexci */
      SUBROUTINE CC_freeze_tripletexci(CAM1,CAMP,CAMM,ISYMTR,
     &           MAXCORE, MAXION,
     &           NRHFCORE,IRHFCORE,NVIRION,IVIRION,
     &           LBOTH)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Purpose: Project out specific excitations 
C              from a trial vector (by zeroing 
C              specific elements)
C     Ex1: zero all ai and aibj elements where i and j
C     are valence orbitals (CORE-VALENCE SEPARATION)
C     Ex2: zero all a an b elements that do not correspond
C     to a specific virtual orbitals
C     Sonia and Eirik 2016
! Based on cc_pram3()
! CAM is the vector analyzed, of symmetry ISYMTR
! LBOTH checks if both CVS and IONISATION are requested 
! Control is passed via argument list, not via common block
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      Implicit none
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      Double precision CAM1(*), CAMP(*), CAMM(*)
      Integer MAXCORE, NRHFCORE(8),IRHFCORE(MAXCORE,8)
      Integer MAXION,NVIRION(8),IVIRION(MAXION,8)
      integer ISYMTR,ISYMAI,ISYMI,ISYMA,ISYMJ,ISYMB,ISYMBJ
      Double precision TWO, THR1, THR2, zero, thprt
      PARAMETER (TWO = 2.0D0,zero=0.0d0, THPRT = 1.0D-9)
      Logical LOCDBG, ikeep, LBOTH
      Parameter (Locdbg = .false.)
      Integer AA, II, MA, MI, JJ, BB, NBJ, NAI, MJ, MB
      Integer KAIBJ, NAIBJ, INDEX
      double precision c1nosq, C2MNOSQ, C2PNOSQ, cnosq
      double precision pt1, ptm, ptp, sumofp, ddot
      integer nl, n1,n2

C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_FREEZE_TRIPLETEXCI')
C
      CALL FLSHFO(LUPRI)
C

C---------------------------------------
C     Loop through single excitation part.
C---------------------------------------
C
!      WRITE(LUPRI,'(//A)')
!     *     ' +=============================================='
!     *    //'===============================+'
!      WRITE(LUPRI,'(1X,A)')
!     *     '| symmetry|  orbital index  |   Excitation Numbers'
!     *     //'             |   Amplitude  |'
!      WRITE(LUPRI,'(1X,A)')
!     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
!     *     //'     NAIBJ    |              |'
!      WRITE(LUPRI,'(A)')
!     *     ' +=============================================='
!     *    //'===============================+'
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
      N1 = 0
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
               if (.not.ikeep) CAM1(NAI) = zero
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C
      CALL FLSHFO(LUPRI)
C--------------------------------------------
C     Loop through Doublee excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCS .OR. CCP2 )) THEN
C
!      WRITE(LUPRI,'(A)')
!     *     ' +----------------------------------------------'
!     *    //'-------------------------------+'
C

      N2 = 0
C
      ikeep = .false.
      DO 200 ISYMAI = 1,NSYM
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
         IF (ISYMAI.LT.ISYMBJ) CYCLE
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
     *                             + INDEX(NAI,NBJ)
                           ELSE
                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                           ENDIF
!***
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
                           if (.not.ikeep) CAMP(NAIBJ) = zero
                           if (.not.ikeep) CAMM(NAIBJ) = zero
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
C
      ENDIF
C
      CALL QEXIT('CC_FREEZE_TRIPLETEXCI')
C
      RETURN
      END

C  /* Deck cc_freeze_tripletstart */
      SUBROUTINE CC_freeze_tripletstart(CAM1,CAMP,CAMM,ISYMTR,
     &           MAXCORE, MAXION,
     &           NRHFCORE,IRHFCORE,NVIRION,IVIRION,
     &           LBOTH)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Filter dia for start vectors
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
C     Sonia and Eirik 2016
C
! Based on cc_pram3()
! CAM is the vector analyzed, of symmetry ISYMTR
! LBOTH checks if both CVS and IONISATION are requested 
! Control is passed via argument list, not via common block
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      Implicit none
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      Double precision CAM1(*), CAMP(*), CAMM(*)
      Integer MAXCORE, NRHFCORE(8),IRHFCORE(MAXCORE,8)
      Integer MAXION,NVIRION(8),IVIRION(MAXION,8)
      integer ISYMTR,ISYMAI,ISYMI,ISYMA,ISYMJ,ISYMB,ISYMBJ
      Double precision TWO, THR1, THR2, crazy, thprt
      PARAMETER (TWO = 2.0D0,crazy=1.0d9, THPRT = 1.0D-9)
      Logical LOCDBG, ikeep, LBOTH
      Parameter (Locdbg = .false.)
      Integer AA, II, MA, MI, JJ, BB, NBJ, NAI, MJ, MB
      Integer KAIBJ, NAIBJ, INDEX
      double precision c1nosq, C2MNOSQ, C2PNOSQ, cnosq
      double precision pt1, ptm, ptp, sumofp, ddot
      integer nl, n1,n2

C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_FREEZE_TRIPLETSTART')
C
C      CALL FLSHFO(LUPRI)
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
      N1 = 0
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
               if (.not.ikeep) CAM1(NAI) = crazy
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through Doublee excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCS .OR. CCP2 )) THEN
C
      N2 = 0
C
      ikeep = .false.
      DO 200 ISYMAI = 1,NSYM
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
         IF (ISYMAI.LT.ISYMBJ) CYCLE
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
     *                             + INDEX(NAI,NBJ)
                           ELSE
                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                           ENDIF
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
                           if (.not.ikeep) CAMP(NAIBJ) = crazy
                           if (.not.ikeep) CAMM(NAIBJ) = crazy
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
C
      ENDIF
C
C
      CALL QEXIT('CC_FREEZE_TRIPLETSTART')
C
      RETURN
      END
