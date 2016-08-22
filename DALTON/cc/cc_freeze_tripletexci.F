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
C------------------------
C     Add up the vectors.
C------------------------
C
!      C1NOSQ = DDOT(NT1AM(ISYMTR),CAM1,1,CAM1,1)
!      IF (.NOT. CCS) THEN
!        C2PNOSQ = DDOT(NT2AM(ISYMTR),CAMP,1,CAMP,1)
!        C2MNOSQ = DDOT(NT2AM(ISYMTR),CAMM,1,CAMM,1)
!      ELSE
!        C2PNOSQ = 0.0D0
!        C2MNOSQ = 0.0D0
!      END IF
C
!      CNOSQ  = C1NOSQ + C2PNOSQ + C2MNOSQ
C
!      IF (CNOSQ .EQ. 0.0D0) THEN
!         PT1 = 0.0D0
!         PTP = 0.0D0
!         PTM = 0.0D0
!      ELSE
!         PT1 = (C1NOSQ /CNOSQ)*100.0D0
!         PTP = (C2PNOSQ/CNOSQ)*100.0D0
!         PTM = (C2MNOSQ/CNOSQ)*100.0D0
!      END IF
C
!      IF (.NOT. CCS .AND. CNOSQ .NE. 0.0D0) THEN
!        WRITE(LUPRI,'(//5X,A)')
!     *     'CC_PRAM:Overall Contribution of the Different Components'
!        WRITE(LUPRI,'(5X,A//)')
!     *     '--------------------------------------------------------'
!        WRITE(LUPRI,'(/5X,A,5X,F10.4,A)')
!     *     'Single Excitation Contribution      : ', PT1,' %'
!        WRITE(LUPRI,'(/5X,A,5X,F10.4,A,F10.4,A)')
!     *     'Double Excitation Contribution (+/-): ',
!     *     PTP,' %    /',PTM,' % '
!        WRITE(LUPRI,'(/5X,A,5X,F10.4)')
!     *     '||T1||/||T2||                       : ',
!     *      SQRT(C1NOSQ/(C2PNOSQ+C2MNOSQ))
!!        IF (LGRS) WRITE(LUPRI,'(/5X,A,5X,F10.4)')
!!     *     'Tau1 diagnostic                     : ',
!!     *      SQRT(C1NOSQ/(TWO*DFLOAT(NRHFT)))
!      END IF
!C
!      WRITE(LUPRI,'(/5X,A,5X,F10.4)')
!     *  'Norm of Total Amplitude Vector : ',SQRT(CNOSQ)
C
      CALL FLSHFO(LUPRI)
C
C----------------------------------------------
C     Initialize threshold etc from Printlevel.
C----------------------------------------------
C
!      NL = MAX(1,2*IPRINT)
!      CNOSQ = MAX(CNOSQ,THPRT)
!      THR1 = SQRT(C1NOSQ/CNOSQ)/NL
!      THR2 = SQRT((C2PNOSQ+C2MNOSQ)/CNOSQ)/NL
!      THR1 = MAX(THR1,1.0D-03)
!      THR2 = MAX(THR2,1.0D-03)
!      SUMOFP = 0.0D00
!C
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
!  1   CONTINUE
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
!                           write(lupri,*)
!     &    'BOTH REQ I:', I, ' II:', II, ' A:', ' AA:', AA
!                           call flshfo(lupri)
                           ikeep = .true.
                           go to 333
                        END IF
                     end do
                  END IF
                 end do
               ELSE
                 do ii = 1, NRHFCORE(ISYMI)
                   IF (I==IRHFCORE(II,ISYMI)) THEN
                       !write(lupri,*)'CORE ONLY'
                       !call flshfo(lupri)
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
               !write(lupri,*)'ikeep=', ikeep
               if (.not.ikeep) CAM1(NAI) = zero
               !WRITE(LUPRI,9990) ISYMA,ISYMI,A,I,NAI,CAM1(NAI)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
!      IF ((N1.LT.1).AND.(SQRT(C1NOSQ/CNOSQ).GT.1.0D-3)) THEN
!         THR1 = THR1/5.0D0
!         GOTO 1
!      ENDIF
C
      CALL FLSHFO(LUPRI)
!      if (.false.) then
C
C--------------------------------------------
C     Loop through Doublee excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCS .OR. CCP2 )) THEN
C
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
C

! 2    CONTINUE
      N2 = 0
C
      ikeep = .false.
      DO 200 ISYMAI = 1,NSYM
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
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
!Sonia: can't remember why we changed this
!                           IF (((ISYMAI.EQ.ISYMBJ).AND.
!     *                         (NAI .LT. NBJ)).OR.(ISYMAI.LT.ISYMBJ))
!     *                          GOTO 260
                           IF ((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ))
     *                          GOTO 260
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
!                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
!     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                               NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                           ENDIF

!                           KAIBJ = NAIBJ + NT1AM(ISYMTR)  !same vector contains singles and doubles
!                                                          !skip the singles
C
!                           IF (ABS(CAMP(NAIBJ)) .GT. THR2 ) THEN
!                              WRITE(LUPRI,9991) ISYMA,ISYMB,ISYMI,ISYMJ,
!     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
!     *                                      CAMP(NAIBJ)
!                              N2 = N2 + 1
!                              SUMOFP = SUMOFP + CAMP(NAIBJ)*CAMP(NAIBJ)
!                           ENDIF
C
!                           IF (ABS(CAMM(NAIBJ)) .GT. THR2 ) THEN
!                              WRITE(LUPRI,9992) ISYMA,ISYMB,ISYMI,ISYMJ,
!     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
!     *                                      CAMM(NAIBJ)
!                              N2 = N2 + 1
!                              SUMOFP = SUMOFP + CAMM(NAIBJ)*CAMM(NAIBJ)
!                           ENDIF
C
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
                             !write(lupri,*)'Get in here'
                             do ii = 1, nrhfcore(isymi)
                              if (i==IRHFCORE(II,ISYMI)) then
                              !write(lupri,*)'loop ii',ii, 'I keep'
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do aa = 1, nvirion(isyma)
                              if (a==IVIRION(aa,ISYMA)) then
                              !write(lupri,*)'loop aa',aa, 'I keep'
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do jj = 1, nrhfcore(isymj)
                              if (j==IRHFCORE(JJ,ISYMJ)) then
                              !write(lupri,*)'loop jj',jj, 'I keep'
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                             do bb = 1, nvirion(isymb)
                              if (b==IVIRION(bb,ISYMB)) then
                              !write(lupri,*)'loop bb',bb, 'I keep'
                                 ikeep = .true.
                                 go to 444
                              end if
                             end do
                           end if
!========================================
  444                      continue
                           if (.not.ikeep) CAMP(NAIBJ) = zero
                           if (.not.ikeep) CAMM(NAIBJ) = zero
!                           if (abs(CAM(KAIBJ)).gt.thr2) then
!                              WRITE(LUPRI,9991) ISYMA,ISYMB,ISYMI,ISYMJ,
!     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
!     *                                      CAM(KAIBJ)
                           !end if
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
!      IF ((N2.LT.1).AND.(SQRT((C2PNOSQ+C2MNOSQ)/CNOSQ).GT.1.0D-3)) THEN
!         THR2 = THR2/5D00
!         GOTO 2
!      ENDIF
C
      ENDIF
C
      !end if !false sonia
!      WRITE(LUPRI,'(A)')
!     *     ' +=============================================='
!     *    //'===============================+'
C
!      WRITE(LUPRI,'(//10X,A,8X,F10.4)')
!     *     'Norm of Printed Amplitude Vector : ',SQRT(SUMOFP)
!      WRITE(LUPRI,'(//10X,A43,1X,F8.6)')
!     *     'Printed all single excitations greater than',THR1
!      IF (.NOT. (CCS.OR.CCP2)) THEN
!         WRITE(LUPRI,'(//10X,A43,1X,F8.6)')
!     *     'Printed all double excitations greater than',THR2
!      ENDIF
C
! 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
!     *       ' | ',12x,' | ',1x, F10.6,'  |')
! 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
!     *       I3,1X,I3,1X,I3,1X,I3,' | ',
!     *       I8,1x,I8,' | (+)',I9,' | ',1x,F10.6,'  |')
! 9992 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
!     *       I3,1X,I3,1X,I3,1X,I3,' | ',
!     *       I8,1x,I8,' | (-)',I9,' | ',1x,F10.6,'  |')
C
      CALL QEXIT('CC_FREEZE_TRIPLETEXCI')
C
      RETURN
      END
