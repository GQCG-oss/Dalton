!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
***********************************************************************

      SUBROUTINE MICDV6(VEC1,VEC2,C2,LU1,LU2,RNRM,EIG,nfinal_vec,
     &                  FINEIG,root_converged,root_residual,
     &                  MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,IROOTHOMING,LUWRTOUT,c_state_of_interest,
     &                  doensemble,weight)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version
*
* From MICDV4, Jeppe Olsen, April 1997
*              Modified : Oct 97 : root homing added
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scratch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "parluci.h"
       DIMENSION VEC1(*),VEC2(*),C2(*)
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
       LOGICAL, intent(in) :: doensemble
       real*8 , intent(in) :: weight(nroot)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
*      IROOTHOMING : Do roothoming, i.e. select the
*      eigenvectors in iteration n+1 as the approximations
*      with largest overlap with the previous space
*
      real(8), intent(out) :: FINEIG(nroot)
      real(8), intent(out) :: root_residual(nroot)
      integer, intent(out) :: root_converged(nroot)
      integer, intent(out) :: nfinal_vec
      integer, intent(in)  :: c_state_of_interest

      LOGICAL CONVER,LUCEND
      REAL*8 INPRDD, INPROD
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*     New characters used for timing in send of VEC from LU3
      CHARACTER CPULUSTEP*12, WALLLUSTEP*12
*     New characters used for timing in cpback of sigma vec and calc of
*     aproj(ij)
      CHARACTER CPUFSTEP*12, WALLFSTEP*12
*     New characters used for timing in step3
      CHARACTER CPUTSTEP*12, WALLTSTEP*12
*     XJEP is used for ROOTHOMING
      real(8), allocatable ::  xjep(:)
      integer, allocatable :: ixjep(:)
      character (len=13)   :: convergence

!     in the Davidson we perform only sigma vector calculations
      ISIGDEN = 1

      lucend = .FALSE.

      MYNEW_ID   = 0
      MYNEW_COMM = 0
C
      WRITE(LUWRTOUT,'(/10X,A)')
     &' **************************************************'
      WRITE(LUWRTOUT,'(10X,A)')
     &'     entering MICDV6 (sequential solver routine)   '
      WRITE(LUWRTOUT,'(10X,A )')
     &' **************************************************'
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
         MAXVEC = MIN(MAXVEC,2)
       END IF
*
!      IPRT = 9999 ! IPRT = 0000  DEBUG
       IPRT = 0000
       lbl  = 0
c      IPRT = IPRTXX
#ifdef LUCI_DEBUG
         WRITE(LUWRTOUT,*) ' Initial set of eigenvectors '
         CALL REWINE(LU1,-1)
         DO IROOT = 1, NROOT
           WRITE(LUWRTOUT,1066) IROOT, LBLK, LBL
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
#endif
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(LUWRTOUT,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(LUWRTOUT,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(LUWRTOUT,*)' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(LUWRTOUT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(LUWRTOUT,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(LUWRTOUT,*) ' Enforced stop on MICDV6 '
         Call Abend1( 20 )
       END IF
*
      if(iroothoming.eq.1)then
        write(luwrtout,*) ' Root homing performed '
        allocate(xjep( maxvec*3))
        allocate(ixjep(maxvec*3))
      end if
*
      KAPROJ = 1
      KFREE  = KAPROJ+ MAXVEC*(MAXVEC+1)/2
      CONVER = .FALSE.
!     set convergence thresholds
!     norm   
      THRQC  = E_CONV
!     THRQC  = 1.0d-04
!     energy change
!     THRQE  = E_CONV * 1.0d-10
      THRQE  = E_CONV
*
* ===================
*.Initial iteration
* ===================
      ITER = 1
      CALL REWINE(LU1,-1)
      CALL REWINE(LU2,-1)
*
      DO IVEC = 1,NINVEC
        CALL REWINE(LU5,-1)
        CALL REWINE(LU6,-1)
        CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
*       Check timings for sigma vector generation
        CALL GETTIM(CPUSIG1,WALLSIG1)

        CALL SIGDEN_CI(VEC1,VEC2,C2,LU5,LU6,cdummy,sdummy,ISIGDEN,-1,
     &                 xxs2dum)

        CALL GETTIM(CPUSIG2,WALLSIG2)
        WALLTID = SECTID(WALLSIG2-WALLSIG1)
        if(ivec .eq. 1)then
          write(luwrtout,9777) WALLTID
        else
          write(luwrtout,9778) WALLTID
        end if
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
        CALL REWINE(LU6,-1)
        CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
        CALL REWINE(LU2,-1)
          DO JVEC = 1, IVEC
             CALL REWINE(LU5,-1)
             IJ = IVEC*(IVEC-1)/2 + JVEC
             APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
          END DO ! loop over IVEC
       END DO ! loop over NINVEC
*
       IF( IPRT .GE.5 ) THEN
         WRITE(LUWRTOUT,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRTOUT,*) 
     &   ' Eigenvalues of initial iteration (with shift)'
         WRITE(LUWRTOUT,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
        WRITE(luwrt,*) ' Eigenvalues of initial iteration (no shift)'
        WRITE(luwrt,'(5F18.13)')
     &  ( EIG(1,IROOT),IROOT=1,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
C     TRAVCD(VEC1,VEC2,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)
*
       IF(IPRT.GE.600) THEN
         WRITE(LUWRTOUT,*) ' Initial set of eigenvectors '
         CALL REWINE(LU1,-1)
         DO IROOT = 1, NROOT
           WRITE(LUWRTOUT,1066) IROOT, LBLK, LBL
 1066 FORMAT(/'lets see IROOT =', 3I4/)
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(LUWRTOUT,*) ' Initial set of sigma vectors '
         CALL REWINE(LU2,-1)
         DO IROOT = 1, NROOT
           WRITE(LUWRTOUT,1067) IROOT
 1067 FORMAT(/'lets see IROOT again=', I4/)
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       CALL dzero(APROJ,NINVEC*(NINVEC+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO
*
       NVEC = NROOT
CJAN20 IF (MAXIT .EQ. 1 ) GOTO  901
 12    CONTINUE 
       IF (MAXIT .EQ. 1 ) GOTO  1001
*
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        CALL GETTIM(CPUITR1,WALLITR1)
        write(LUWRTOUT,*)
        write(LUWRTOUT,'(A21,3X,I3)') ' Info from iteration ',ITER
        write(LUWRTOUT,*) '_______________________'
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       CALL REWINE(LU1,-1)
       CALL REWINE(LU2,-1)
*
       DO 100 IROOT = 1, NROOT
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINE(LU7,-1)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
*
         CALL REWINE(LU7,-1)
         CALL REWINE(LU5,-1)
         FACTOR = - EIGAPR
         CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU2,LU7,LU5,0,LBLK)
*
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(LUWRTOUT,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU5,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
         WRITE(LUWRTOUT,'(A19,I8,1p,1E15.5,0p,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
*
         CALL FLSHFO(LUWRTOUT)
*
!        write(LUWRTOUT,*) 'RNORM, THRQC, abs, THRQE',
!    &                      RNORM, THRQC, 
!    &                      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)),
!    &                      THRQE

!hjajj/s IF(RNORM.LT. THRQC .OR.
!    &      (ITER.GT.2.AND.
!    &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.THRQE)) THEN
         IF(RNORM.LT. THRQC)THEN
            root_converged(IROOT) = 0
         ELSE
            root_converged(IROOT) = iroot
            CONVER = .FALSE.
         END IF
         if(c_state_of_interest /= 0)then
           if(root_converged(abs(c_state_of_interest)) == 0)then ! desired root converged - we are done!
             conver = .true. 
           end if
         end if
         !> do not converge any undesired root
         if(doensemble)then
           if(abs(weight(iroot)) == 0.0d0)then
             root_converged(IROOT) = 0
           end if
         end if
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU6
CSK         print*,'NP1,NP2,NQ ',NP1,NP2,NQ
         IF( root_converged(IROOT) .gt. 0) THEN
           IF(IPRT.GE.3) THEN
             WRITE(LUWRTOUT,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINE(LUDIA,-1)
           CALL REWINE(LU5,-1)
           CALL REWINE(LU6,-1)
           CALL H0M1TD(LU6,LUDIA,LU5,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
C     H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
           IF ( IPRT  .GE. 600) THEN
             WRITE(LUWRTOUT,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU6,1,LBLK)
           END IF
CSK        print*,'IOLSTM ',IOLSTM
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU5
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
             CALL REWINE(LUDIA,-1)
*
             CALL H0M1TD(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
             GAMMA = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
CSK             WRITE(LUWRTOUT,*) ' Gamma = ', Gamma
             IF(IPRT.GE.1000) WRITE(LUWRTOUT,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU7,LU5,0,LBLK))
CSK              write(LUWRTOUT,*) ' VNORM = ', VNORM
              IF(IPRT.GE.1000) write(LUWRTOUT,*) ' VNORM = ', VNORM
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
             IF(IPRT.GE.5) WRITE(LUWRTOUT,*) ' Olsen Correction active '
CSK                WRITE(LUWRTOUT,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                FACTOR = -(DELTA/GAMMA)
                IF(IPRT.GE.5) WRITE(LUWRTOUT,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(LUWRTOUT,*)   DELTA,GAMMA,FACTOR
CSK                WRITE(LUWRTOUT,*) ' DELTA,GAMMA,FACTOR'
CSK                WRITE(LUWRTOUT,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU6,LU5,LU7,1,LBLK)
                CALL COPVCD(LU7,LU6,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(LUWRTOUT,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU6,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
*.. Vectors on LU1
*
           CALL REWINE( LU1 ,LBLK)
           DO IVEC = 1, NROOT
               CALL REWINE(LU6,-1)
               WORK(IVEC) = -INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
           END DO
           IF(IPRT.GE.1000) THEN
             Write(LUWRTOUT,*) ' Overlap with vectors on LU1'
             call WRTMT_LU(work,1,nroot,1,nroot)
           END IF
           ONE = 1.0D0
           CALL REWINE(LU1,-1)
           CALL MVCSMD(LU1,WORK,LU7,LU5,VEC1,VEC2,NROOT,1,LBLK)
           CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU6,LU5,1,LBLK)
           IF(IPRT.GE.1000) THEN
             write(LUWRTOUT,*) ' orthogonalized to vectors on LU1'
             CALL WRTVCD(VEC1,LU5,1,LBLK)
           END IF
*.. Vectors on LU3
           IF(NVEC+IADD-1-NROOT.GT.0) THEN
             CALL REWINE( LU3 ,LBLK)
             DO IVEC = 1, NVEC+IADD-1-NROOT
                 CALL REWINE(LU6,-1)
                 WORK(IVEC) = -INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
             END DO
             ONE = 1.0D0
             CALL REWINE(LU3,-1)
             CALL MVCSMD(LU3,WORK,LU7,LU6,VEC1,VEC2,NVEC+IADD-1-NROOT,
     &            1,LBLK)
             CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU5,LU6,1,LBLK)
           ELSE
             CALL REWINE(LU3,-1)
             CALL COPVCD(LU5,LU6,VEC1,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU6,LBLK)
           CALL SCLVCD(LU6,LU3,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU6,LU7,FACTOR,VEC1,1,LBLK)
             WRITE(LUWRTOUT,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU7,1,LBLK)
           END IF
*
         END IF
  100 CONTINUE
*
      IF( CONVER ) THEN
         GOTO  1001
      END IF
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
*
* ====================================================
*  2 : Optimal combination of new and old directions
* ====================================================
*
*  2.1: Multiply new directions with matrix
*
      IF(IPRT.GE.1000) THEN
       IF (luci_MYPROC.EQ.luci_MASTER) THEN
        WRITE(LUWRTOUT,*) ' Vectors on LU3'
        WRITE(LUWRTOUT,*) 'NVEC-NROOT+IADD =',NVEC-NROOT+IADD
        CALL REWINE(LU3,-1)
        DO IVEC = 1, NVEC-NROOT+IADD
          CALL WRTVCD(VEC1,LU3,0,LBLK)
        END DO
       END IF
       CALL REWINE(LU3,-1)
      END IF
*
      CALL SKPVCD(LU3,NVEC-NROOT,VEC1,1,LBLK)
      CALL SKPVCD(LU4,NVEC-NROOT,VEC1,1,LBLK)
*
*     
 150  DO IVEC = 1, IADD
        CALL REWINE(LU5,LBLK)
        CALL REWINE(LU6,LBLK)
        CALL COPVCD(LU3,LU5,VEC1,0,LBLK)

!       write(luwrtout,*) 'new trial vec ...'
!       CALL WRTVCD(VEC1,LU5,1,LBLK)
*Check timings for sigma vector generation
        CALL GETTIM(CPUSIG1,WALLSIG1)
        CALL SIGDEN_CI(VEC1,VEC2,C2,LU5,LU6,cdummy,sdummy,ISIGDEN,-1,
     &                 xxs2dum)
        CALL GETTIM(CPUSIG2,WALLSIG2)
        WALLTID = SECTID(WALLSIG2-WALLSIG1)
        if(ivec .eq. 1)then
          WRITE(LUWRTOUT,9400) WALLTID
        else
          WRITE(LUWRTOUT,9401) WALLTID
        end if
        CALL REWINE(LU6,LBLK)
        CALL COPVCD(LU6,LU4,VEC1,0,LBLK)
*2.2 Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO JVEC = 1, NROOT
          CALL REWINE(LU6,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
!         write(LUWRT,*) 'ij, aproj(ij):',ij,aproj(ij)
        END DO

        CALL REWINE(LU3,LBLK)
C       DO JVEC = NROOT+1, NVEC+IADD
        DO JVEC = NROOT+1, NVEC+IVEC
         CALL REWINE(LU6,LBLK)
         IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
         APROJ(IJ) = INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
!        write(LUWRT,*) 'ij, aproj(ij):',ij,aproj(ij)
        END DO
      END DO
*     /\ End do over new trial vectors
*
*. 2.3 Diagonalize projected matrix
      NVEC = NVEC + IADD
!     write(LUWRT,*) 'Dimension of projected matrix:',NVEC
!     call wrtmatmn(aproj,1,NVEC*(NVEC+1)/2,1,NVEC*(NVEC+1)/2,luwrtout)
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NVEC,0,1)
*. Test : transform projected matrix
C TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
C     CALL TRAN_SYM_BLOC_MAT(APROJ,AVEC,1,NVEC,XJEP(1000),XJEP(1))
C     WRITE(LUWRTOUT,*) ' Explicitly transformed matrix'
C     CALL PRSYM(XJEP(1000),NVEC)

      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(LUWRTOUT,*) 
     &  ' Var and Pert solution, energy and coefficients'
        WRITE(LUWRTOUT,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(LUWRTOUT,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
*
      IF(IROOTHOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first
*
*. Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
       END DO
       WRITE(LUWRTOUT,*)
     & ' Norm of projections to previous vector space '
       CALL WRTMT_LU(XJEP,1,NVEC,1,NVEC)
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL SCALVE(XJEP,ONEM,NVEC)
       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
       IF(NSORT.LT.NVEC) THEN
         WRITE(LUWRTOUT,*) ' Warning : Some elements lost in sorting '
         WRITE(LUWRTOUT,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
       IF(IPRT.GE.3) THEN
         WRITE(LUWRTOUT,*) ' New roots choosen as vectors '
         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
       END IF
*. Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
       END DO
       CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
       END DO
       DO INEW = 1, NVEC
         WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
       END DO
*
       IF(IPRT.GE.3) THEN
         WRITE(LUWRTOUT,*) ' Reordered WORK and AVEC arrays '
         CALL PRSYM(WORK,NVEC)
         CALL WRTMT_LU(AVEC,NVEC,NVEC,NVEC,NVEC)
       END IF
*
      END IF
*     ^ End of root homing procedure
      DO IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
      END DO
*
!     IPRT = 6

      IF(IPRT .GE. 3 ) THEN
        WRITE(LUWRTOUT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(LUWRTOUT,'(5F18.13)')
     &  ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
        WRITE(LUWRTOUT,'(A)') ' Norm of Residuals (Previous it) '
        WRITE(LUWRTOUT,'(5F18.13)')
     &  ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
      END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRTOUT,*) ' Projected matrix'
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRTOUT,*) ' eigen pairs '
        WRITE(LUWRTOUT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMT_LU(AVEC,NVEC,NROOT,NVEC,NROOT)
      END IF

!     IPRT = 0
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset
*
*
* case 1 : Only NROOT vectors can be stored
*          save current eigenvector approximations
* Case 2 : Atleast 2*NROOT eigenvectors can be saved
*          Current eigenvactor approximations+
*          vectors allowing generation of previous approxs.
*
*
*
*  check timing of part3 in loop iteration
      CALL GETTIM(CPUSTEP3,WALLSTEP3)

#ifdef debug_blabla
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU5,LBLK)
        DO 320 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  320   CONTINUE
*. Transfer C vectors to LU1
        CALL REWINE( LU5,LBLK)
        CALL REWINE( LU1,LBLK)
        DO 411 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  411   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 329 IROOT = 1, NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  329   CONTINUE
*
* Transfer HC's to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 400 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  400   CONTINUE
        NVEC = NROOT
*
        CALL DZERO(AVEC,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
*
      END IF
*
#endif

!     IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
!     activating the IF... part above does NOT work as it is...
!     TODO: transfer the davidson subroutine from kr-ci (dirac) which 
!           has a working increasing subspace dimension routine.
!           it's basically the code you find above in #ifdef ...
*
*. 3.1 : Current wave function approximations, collect on LU7
!#ifdef blabla_debug
!       IPRT = 1000
        IF(IPRT.GE.1000) THEN
        write(LUWRTOUT,*) ' I am going to reset '
        write(LUWRTOUT,*) ' nroot, nvec ', nroot,nvec
        END IF
        IF(IPRT.GE.1000) THEN
          WRITE(LUWRTOUT,*) ' Initial vectors on LU1'
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
          WRITE(LUWRTOUT,*) ' Initial vectors on LU3'
          CALL REWINE(LU3,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF
*
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU1 to LU5
          CALL REWINE(LU1,-1)
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of c reset, part 1'
*. and add LU3 to LU5
          CALL REWINE(LU3,-1)
          ONE = 1.0D0
          CALL MVCSMD2(LU3,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of c reset, part 2'
C              MVCSMD2(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
C    &         VEC1,VEC2,NVEC,IREW,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
          CALL REWINE(LU5,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
CSK          WRITE(LUWRTOUT,*) ' SCALE = ', SCALE
          WORK(IROOT) = SCALE
*. scale LU5 => LU7
          CALL REWINE(LU5,-1)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of c reset, part 3'
        END DO
*. Transfer C vectors from LU7 to LU1
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU1,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU1,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(LUWRTOUT,*) ' end of c reset, part 4'
        IF(IPRT.GE.1000) THEN
          WRITE(LUWRTOUT,*) ' Reset C-vectors on LU1 '
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
        END IF
*. 3.2 : corresponding sigma vectors
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU2 to LU5
          CALL REWINE(LU2,-1)
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of s reset, part 1'
*. and add LU4 to LU5
          CALL REWINE(LU4,-1)
          CALL MVCSMD2(LU4,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of s reset, part 2'
          SCALE  = WORK(IROOT)
*. scale LU5 => LU7
          CALL REWINE(LU5,-1)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of s reset, part 3'
        END DO
*. Transfer sigma  vectors from LU7 to LU2
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU2,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU2,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(LUWRTOUT,*) ' end of s reset, part 4'
        NNVEC = NROOT

        IF(3*NROOT.LE.MAXVEC) THEN
*
*. Orthogonalize the
*. last set of correction vectors to the current
*. eigenvectors on LU1, and save on LU2
*. Overlap with root approximations
*. Start of last set of trial vectors
          ISTART = NVEC-NROOT-IADD+1
          CALL SKPVCD(LU3,ISTART-1,VEC1,1,LBLK)
*
          CALL REWINE(LU5,-1)
          DO JVEC = 1, IADD
*. Orthogonalize to vectors on LU1
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU3,LU7,VEC1,0,LBLK)
            IF(IPRT.GE.1000) THEN
            write(LUWRTOUT,*) ' Initial vector on LU7'
            CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
            CALL REWINE(LU1,-1)
            DO IROOT = 1, NROOT
              CALL REWINE(LU7,-1)
              WORK(IROOT+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU1,0,LBLK)
            END DO
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' c, part1 finito'
*. And to trial vectors on LU5
            CALL REWINE(LU5,-1)
            DO KVEC = 1, JVEC-1
              CALL REWINE(LU7,-1)
              WORK(NROOT+KVEC+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU5,0,LBLK)
            END DO
            WORK(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' c, part2 finito'
*
            ONE = 1.0D0
            CALL MVCSMD2(LU1,WORK(1+(JVEC-1)*2*NROOT),ONE ,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            ONE = 1.0D0
            CALL MVCSMD2(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000) THEN
              write(LUWRTOUT,*) ' c, part4 finito'
              write(LUWRTOUT,*) ' Vector after sec ort'
              CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
*
            FACTOR = INPRDD(VEC1,VEC1,LU7,LU7,1,LBLK)
            SCALE = 1.0D0/SQRT(FACTOR)
            CALL SCALVE(WORK((JVEC-1)*2*NROOT+1),SCALE,
     &           NROOT+JVEC)
            CALL REWINE(LU7,-1)
            CALL SCLVCD(LU7,LU5,SCALE,VEC1,0,LBLK)
C                SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
           IF(IPRT.GE.1000)
     &       write(LUWRTOUT,*) ' c, part6 finito'
         END DO
*        /\ End of loop over orthogonalized directions
*. Transfer modified directions to LU3
         CALL REWINE(LU3,-1)
         CALL REWINE(LU5,-1)
         DO JVEC =1, IADD
           CALL COPVCD(LU5,LU3,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000) THEN
           write(LUWRTOUT,*) ' c, part7 finito'
           WRITE(LUWRTOUT,*) ' Additional trial vectors on LU3'
           CALL REWINE(LU3,-1)
           DO JVEC = 1, IADD
            CALL WRTVCD(VEC1,LU3,0,LBLK)
           END DO
         END IF
* Sigma vectors corresponding to orthogonalized directions
         CALL SKPVCD(LU4,ISTART-1,VEC1,1,LBLK)
         CALL REWINE(LU5,-1)
         DO JVEC = 1, IADD
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU4,LU7,VEC1,0,LBLK)
*
            FACT = WORK(NROOT+JVEC+(JVEC-1)*2*NROOT)
            CALL MVCSMD2(LU2,WORK(1+(JVEC-1)*2*NROOT),FACT,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' s, part 1 finito '
            ONE = 1.0D0
            CALL MVCSMD2(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' s, part 2 finito '
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU7,LU5,VEC1,0,LBLK)
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' s, part 3 finito '
         END DO
*        /\ End of loop over orthogonalized directions
*. Copy back to LU4
         CALL REWINE(LU4,-1)
         CALL REWINE(LU5,-1)
         DO JVEC = 1, IADD
            CALL COPVCD(LU5,LU4,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000)
     &   write(LUWRTOUT,*) ' s, part 4 finito '
C         WRITE(LUWRTOUT,*) 'I was in 3.1 setting NNVEC=NROOT+IADD', 
C     &luci_myproc,(NROOT + IADD)
         NNVEC = NROOT + IADD
       END IF
       IF(IPRT.GE.1000) THEN
         WRITE(LUWRTOUT,*) ' Additional sigma vectors on LU4'
         CALL REWINE(LU4,-1)
         DO JVEC = 1, IADD
           CALL WRTVCD(VEC1,LU4,0,LBLK)
         END DO
       END IF
*
*      /\ End if more than NROOT vectors to be reset
       NVEC = NNVEC
*       write(LUWRTOUT,*)'NVEC(1)',NVEC NVEC = IADD (iadd can change!!!)+ NROOT

!      END IF ! reset

*. New subspace
*. Calculate subspace Hamiltonian from actual vectors
*. on disc
       IF(IPRT.GE.1000) write(LUWRTOUT,*) ' Subspace hamiltonian'
       CALL REWINE(LU1,-1)
       CALL REWINE(LU3,-1)
CSK       ididaproj = 0
CSK       write(LUWRTOUT,*)'NVEC = ', NVEC
       DO IVEC = 1, NVEC
*
         CALL REWINE(LU5,-1)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINE(LU2,-1)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
         END DO
         CALL REWINE(LU4,-1)
         DO JVEC = NROOT+1,IVEC
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
         END DO
CSK        ididaproj = ididaproj + 1
       END DO
       if (IPRT.ge.2) then
         write(LUWRTOUT,*) ' Reset hamiltonian'
         call prsym(aproj,nvec)
       end if
!      IPRT = 0

!#endif
CSK       write(LUWRTOUT,*)'actual dimension of aproj',MAXVEC*(MAXVEC+1)/2
CSK       write(LUWRTOUT,*)'calculated elements of aproj',ididaproj
*  check timing of part3 in loop iteration
       CALL GETTIM(CPUSTEP4,WALLSTEP4)
       WALLTSTEP = SECTID(WALLSTEP4-WALLSTEP3)
       WRITE(LUWRTOUT,9600) WALLTSTEP
*
 998   CONTINUE
*  Timing of this iteration
       CALL GETTIM(CPUITR2,WALLITR2)
       WALLTID = SECTID(WALLITR2-WALLITR1)
       WRITE(LUWRTOUT,9300) WALLTID
*. End of resetting business
 999  CONTINUE
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
* ( End of loop over iterations )
C
C
      IF( .NOT. CONVER ) THEN
*       CONVERGENCE WAS NOT OBTAINED
        IF(IPRT .GE. 2 )
     &  WRITE(LUWRTOUT,1170) MAXIT
 1170   FORMAT('   Convergence was not obtained in ',I3,' iterations')
      ELSE
*       CONVERGENCE WAS OBTAINED
        ITER = ITER - 1
        IF (IPRT .GE. 2 )
     &  WRITE(LUWRTOUT,1180) ITER
 1180   FORMAT(1X,' Convergence was obtained in ',I3,' iterations')
      END IF
*
      do iroot = 1, nroot
        write(LUWRTOUT,'(/1x,a)')    '-------------------'
        write(LUWRTOUT,'( 1x,a,i6)') 'root number =',IROOT
        write(LUWRTOUT,'(1x,a/)')    '-------------------'
        do i = 1, iter
          write(LUWRTOUT,'( 1x,a9,i4,f20.10,4x,1p,e10.3)')
     & 'iteration',i,eig(i,iroot)+eigshf,rnrm(i,iroot)
        end do
      end do

      write(LUWRTOUT,'(/1x,a)') '**********************************'//
     &           '******************************'
      write(LUWRTOUT,*) '  iter   root         energy          '//
     &           'residual     convergence'
      write(LUWRTOUT,'(1x,a )') '**********************************'//
     &           '******************************'
      DO 1601 IROOT = 1, NROOT
         FINEIG(IROOT)        = EIG(ITER,IROOT)+EIGSHF
         root_residual(IROOT) = rnrm(iter,iroot)
         if(root_converged(iroot) .eq. 0) then
           convergence = ' converged'
         else
           convergence = 'NOT converged'
         end if
         WRITE(LUWRTOUT,'(i7,i7,f20.10,4x,1p,E10.3,4x,a)')
     &         ITER,IROOT,FINEIG(IROOT),root_residual(IROOT),
     &         convergence
 1601 CONTINUE
      nfinal_vec = nvec
*
      if(iroothoming.eq.1)then
        deallocate(ixjep)
        deallocate( xjep)
      end if

      CALL FLSHFO(LUWRTOUT)
*
      RETURN
 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION            : ',A)
 9400 FORMAT(/' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9401 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9600 FORMAT(/' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A)
 9777 FORMAT(/' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9778 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
      END
***********************************************************************

      subroutine blubb_micdv4(VEC1,VEC2,C2,LU1,LU2,RNRM,EIG,nfinal_vec,
     &                  FINEIG,root_converged,root_residual,
     &                  MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,lupri)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Version using H0 + Lambda V as Sigma routine
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "parluci.h"
       DIMENSION VEC1(*),VEC2(*), C2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2                               
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       real(8), intent(out) :: FINEIG(nroot)
       LOGICAL CONVER,RTCNV(NROOT)
       real(8), intent(out) :: root_residual(nroot)
       integer, intent(out) :: root_converged(nroot)
       integer, intent(out) :: nfinal_vec
       REAL*8 INPRDD
       CHARACTER SECTID*12, CPUTID*12, WALLTID*12,FILELAB*8
*
       call qenter('MICDV4')
*
       WRITE(LUWRT,*) '             ::::::::::::::::::::::::   '
       WRITE(LUWRT,*) '               Entering MICDV4 (rel)    '
       WRITE(LUWRT,*) '             ::::::::::::::::::::::::   '
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(LUWRT,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(LUWRT,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(LUWRT,*) ' Sorry MICDV4 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(LUWRT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(LUWRT,*) ' Raise MXCIV_CI to be at least 2 * Nroot '
         WRITE(LUWRT,*) ' Enforced stop on MICDV4 '
         STOP 20
       END IF
*
!      IPRT = 5000 ! debug
       IPRT = 0010

       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINE(LU1,-1)
       CALL REWINE(LU2,-1)
       DO 10 IVEC = 1,NINVEC
         CALL REWINE(LU3,-1)
         CALL REWINE(LU4,-1)
         CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
        CALL SIGDEN_CI(VEC1,VEC2,C2,LU3,LU4,cdummy,sdummy,1,-1,
     &                 xxs2dum)
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINE(LU4,-1)
         CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINE(LU2,-1)
         DO 8 JVEC = 1, IVEC
           CALL REWINE(LU3,-1)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU3,0,LBLK)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.3 ) THEN
           WRITE(LUWRT,*) ' INITIAL PROJECTED MATRIX  '
           CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL DCOPY(NINVEC*(NINVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
       CALL EIGEN_luci(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*

       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRT,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(LUWRT,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(LUWRT,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        CALL GETTIM(CPUITR1,WALLITR1)
        write(LUWRT,'(/A,I6/A)') ' Info from iteration',ITER,
     &                             '___________________________'
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         EIGAPR = EIG(ITER-1,IROOT)
*
         CALL REWINE(LU1,-1)
         CALL REWINE(LU2,-1)
         EIGAPR = EIG(ITER-1,IROOT)
         DO 60 IVEC = 1, NVEC
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           IF(IVEC.EQ.1) THEN
             CALL REWINE( LU3 ,-1)
*                 SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
             CALL SCLVCD(LU2,LU3,FACTOR,VEC1,0,LBLK)
           ELSE
             CALL REWINE(LU3,-1)
             CALL REWINE(LU4,-1)
C                 VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU2,LU3,0,LBLK)
           END IF
C
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINE(LU3,-1)
           CALL REWINE(LU4,-1)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(LUWRT,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
         END IF
*  Strange place to put convergence but ....
C                      INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
*
         WRITE(LUWRT,'(A19,7X,I3,3X,1E18.13,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
         CALL FLSHFO(LUWRT)
*
         IF(RNORM.LT. e_conv)THEN
            root_converged(IROOT) = 0
         ELSE
            root_converged(IROOT) = iroot
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new
*        direction
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU3
         IF( root_converged(IROOT) /= 0) THEN
           IF(IPRT.GE.3) THEN
             WRITE(LUWRT,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINE(LUDIA,-1)
           CALL REWINE(LU3,-1)
           CALL REWINE(LU4,-1)
*. Assuming diagonal preconditioner
           IPRECOND = 1
           CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
C               H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND)
           IF ( IPRT  .GE. 600) THEN
             WRITE(LUWRT,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* Current eigen-vector on LU5
             CALL REWINE(LU1,-1)
             DO 66 IVEC = 1, NVEC
               FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
               IF(IVEC.EQ.1) THEN
                 IF(NVEC.EQ.1) THEN
                   CALL REWINE( LU5 ,-1)
                   CALL SCLVCD(LU1,LU5,FACTOR,VEC1,0,LBLK)
                 ELSE
                   CALL REWINE( LU4 ,-1)
                   CALL SCLVCD(LU1,LU4,FACTOR,VEC1,0,LBLK)
                 END IF
               ELSE
                 CALL REWINE(LU5,-1)
                 CALL REWINE(LU4,-1)
                 CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                 CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
               END IF
   66        CONTINUE
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(LUWRT,*) '  (current  X ) '
               CALL WRTVCD(VEC1,LU5,1,LBLK)
             END IF
* (H0 - E )-1  * X on LU4
             CALL REWINE(LU5,-1)
             CALL REWINE(LU4,-1)
             CALL REWINE(LUDIA,-1)
*
             CALL H0M1TD(LU4,LUDIA,LU5,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPRDD(VEC1,VEC2,LU5,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU5,LU4,1,LBLK))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(LUWRT,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU5,LU3,1,LBLK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(LUWRT,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(LUWRT,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
                CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(LUWRT,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU3,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
           CALL REWINE( LU1 ,LBLK)
           DO 80 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             WORK(IVEC) = INPRDD(VEC1,VEC2,LU1,LU3,0,LBLK)
C?       WRITE(6,*) ' MICDV4 : Overlap ', WORK(IVEC)
   80      CONTINUE
*
           CALL REWINE(LU1,LBLK)
           DO 82 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             CALL REWINE(LU4,LBLK)
             CALL VECSMD(VEC1,VEC2,-WORK(IVEC),1.0D0,LU1,LU3,
     &                   LU4,0,LBLK)
             CALL COPVCD(LU4,LU3,VEC1,1,LBLK)
   82      CONTINUE
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(LUWRT,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU3,LU4,FACTOR,VEC1,1,LBLK)
             WRITE(LUWRT,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,NVEC,VEC1,1,LBLK)
      DO 150 IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU1,LU3,VEC1,0,LBLK)

        write(luwrt,*) 'new c trial vec ...'
        CALL WRTVCD(VEC1,LU3,1,LBLK)

        CALL SIGDEN_CI(VEC1,VEC2,C2,LU3,LU4,cdummy,sdummy,1,-1,
     &                 xxs2dum)

        CALL REWINE(LU4,LBLK)
        CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO 140 JVEC = 1, NVEC+IVEC
          CALL REWINE(LU4,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU4,0,LBLK)
          write(LUWRT,*) 'ij, aproj(ij):',ij,aproj(ij)
  140   CONTINUE
  150 CONTINUE
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      write(LUWRT,*) 'Dimension of projected matrix:',NVEC
      call wrtmatmn(aproj,1,NVEC*(NVEC+1)/2,1,NVEC*(NVEC+1)/2,luwrt)
      CALL DCOPY(NVEC*(NVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
      CALL EIGEN_luci(WORK(KAPROJ),AVEC,NVEC,0,1)
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(LUWRT,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(LUWRT,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(LUWRT,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRT,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset      
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU5,LBLK)
        DO 320 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  320   CONTINUE
*. Transfer C vectors to LU1
        CALL REWINE( LU5,LBLK)
        CALL REWINE( LU1,LBLK)
        DO 411 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  411   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 329 IROOT = 1, NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  329   CONTINUE
*
* Transfer HCs to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 400 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  400   CONTINUE
        NVEC = NROOT
*
        CALL DZERO(AVEC,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
*
      END IF
*
*  Timing of this iteration
      CALL GETTIM(CPUITR2,WALLITR2)
      CPUTID = SECTID(CPUITR2-CPUITR1)
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRT,9300) CPUTID,WALLTID
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
 
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(LUWRT,1170) MAXIT
 1170    FORMAT(/'  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(LUWRT,1180) ITER
 1180    FORMAT(/,'  Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 0 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(LUWRT,*)
          WRITE(LUWRT,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(LUWRT,*)
     &    '============================================'
          WRITE(LUWRT,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          IF (RTCNV(IROOT)) THEN
             WRITE(6,1190) IROOT,FINEIG(IROOT)+EIGSHF
          ELSE
             WRITE(6,1191) IROOT,FINEIG(IROOT)+EIGSHF
          END IF
 1190 FORMAT(' The final eigenvalue',I5,F22.10,' (converged)')
 1191 FORMAT(' The final eigenvalue',I5,F22.10,' (NOT converged)')
          IF(IPRT.GE.400) THEN
            WRITE(LUWRT,1200)
 1200       FORMAT(' The final approximation to eigenvector ')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(LUWRT,1300)
 1300     FORMAT(/' Summary of iterations',
     +           /' ---------------------')
          WRITE(LUWRT,1310)
 1310     FORMAT
     &    (/,'  Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(LUWRT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        CALL REWINE(LU1,LBLK)
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
        DO IROOT = 1, NROOT
          FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           root_residual(IROOT) = rnrm(iter,iroot)
          IF (RTCNV(IROOT)) THEN
             WRITE(6,1190) IROOT,FINEIG(IROOT)
          ELSE
             WRITE(6,1191) IROOT,FINEIG(IROOT)
          END IF
        END DO
        write(LUWRT,*)
        write(LUWRT,*) '*++++++++++++++++++++++++++++++++++++++++++++*'
        write(LUWRT,*)
      END IF
      nfinal_vec = nvec
*
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(LUWRT,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      call qexit('MICDV4')
C
 1030 FORMAT(/,3X,7F15.8,/,(2X,7F15.8))
 1120 FORMAT(/,I6,7F15.8,/,(5X,7F15.8))
 9300 FORMAT(' >>>  CPU (WALL) TIME IN ITERATION: ',A,'(',A,')')
      END

      SUBROUTINE MICDV4_ENLMD(VEC1,VEC2,C2,LU1,LU2,RNRM,EIG,nfinal_vec,
     &                  FINEIG,root_converged,root_residual,
     &                  MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,lupri)
*
* Davidson algorithm , requires two blocks in core
* Multi root version
*
* Jeppe Olsen Winter of 1991
*
* Updated to allow general preconditioner, October 1993
*
* Version using H0 + Lambda V as Sigma routine
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*), C2(*)
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2                               
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       real(8), intent(out) :: FINEIG(nroot)
       real(8), intent(out) :: root_residual(nroot)
       integer, intent(out) :: root_converged(nroot)
       integer, intent(out) :: nfinal_vec

       LOGICAL CONVER
       REAL*8 INPRDD, INPROD
*
      iprt = 6
!     in the Davidson we perform only sigma vector calculations
      ISIGDEN = 1

       IPICO = 0
       IF(IPICO.NE.0) THEN
C?       WRITE(lupri,*) ' Perturbative solver '
         MAXVEC = MIN(MAXVEC,2)
       ELSE IF(IPICO.EQ.0) THEN
C?       WRITE(lupri,*) ' Variational  solver '
       END IF
*
      WRITE(lupri,'(/10X,A)')
     &' **************************************************'
      WRITE(lupri,'(10X,A)')
     &'     entering MICDV4 (sequential solver routine)   '
      WRITE(lupri,'(10X,A )')
     &' **************************************************'
 
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(lupri,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(lupri,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(lupri,*) ' Sorry MICDV4 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(lupri,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(lupri,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(lupri,*) ' Enforced stop on MICDV4 '
         STOP 20
       END IF
*
       KAPROJ = 1
       KFREE  = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST   = 1.0D-8
       CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
       ITER = 1
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       DO 10 IVEC = 1,NINVEC
         CALL REWINO(LU3)
         CALL REWINO(LU4)
         CALL COPVCD(LU1,LU3,VEC1,0,LBLK)

         CALL SIGDEN_CI(VEC1,VEC2,C2,LU3,LU4,cdummy,sdummy,ISIGDEN,-1,
     &                  xxs2dum)

*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
         CALL REWINO(LU4)
         CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Projected matrix
         CALL REWINO(LU2)
         DO 8 JVEC = 1, IVEC
           CALL REWINO(LU3)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU3,0,LBLK)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.3 ) THEN
         WRITE(lupri,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(lupri,'(A,I4)') ' Eigenvalues of initial iteration '
         WRITE(lupri,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
         WRITE(lupri,*) ' Initial set of eigen values (no shift) '
         CALL WRTMAT(EIG(1,1),1,NROOT,MAXIT,NROOT)
       END IF
       NVEC = NINVEC

       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)

*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       CALL dzero(APROJ,NROOT*(NROOT+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO

       IF (MAXIT .EQ. 1 ) GOTO  901
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        IF(IPRT  .GE. 10 ) THEN
         WRITE(lupri,*) ' Info from iteration .... ', ITER
        END IF
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT

         EIGAPR = EIG(ITER-1,IROOT)
*
         CALL REWINO(LU1)
         CALL REWINO(LU2)
         DO 60 IVEC = 1, NVEC
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           IF(IVEC.EQ.1) THEN
             CALL REWINO( LU3 )
*                 SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
             CALL SCLVCD(LU2,LU3,FACTOR,VEC1,0,LBLK)
           ELSE
             CALL REWINO(LU3)
             CALL REWINO(LU4)
C                 VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU2,LU3,0,LBLK)
           END IF
C
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
           CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU1,LU4,0,LBLK)
   60    CONTINUE
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(lupri,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU4,1,LBLK)
         END IF
*  Strange place to put convergence but ....
C                      INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU4,LU4,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
         WRITE(lupri,'(A19,I8,1p,1E15.5,0p,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
*
         IF(RNORM.LT. e_conv)THEN
            root_converged(IROOT) = 0
         ELSE
            root_converged(IROOT) = iroot
            CONVER = .FALSE.
         END IF

         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU3
         IF( root_converged(IROOT) .gt. 0 ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(lupri,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINO(LUDIA)
           CALL REWINO(LU3)
           CALL REWINO(LU4)
*. Assuming diagonal preconditioner
           IPRECOND = 1
           CALL H0M1TD(LU3,LUDIA,LU4,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
C               H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG,IPRECOND)
           IF ( IPRT  .GE. 600) THEN
             WRITE(lupri,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* Current eigen-vector on LU5
             CALL REWINO(LU1)
             DO 66 IVEC = 1, NVEC
               FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
               IF(IVEC.EQ.1) THEN
                 IF(NVEC.EQ.1) THEN
                   CALL REWINO( LU5 )
                   CALL SCLVCD(LU1,LU5,FACTOR,VEC1,0,LBLK)
                 ELSE
                   CALL REWINO( LU4 )
                   CALL SCLVCD(LU1,LU4,FACTOR,VEC1,0,LBLK)
                 END IF
               ELSE
                 CALL REWINO(LU5)
                 CALL REWINO(LU4)
                 CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU4,LU1,LU5,0,LBLK)
                 CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
               END IF
   66        CONTINUE
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(lupri,*) '  (current  X ) '
               CALL WRTVCD(VEC1,LU5,1,LBLK)
             END IF
* (H0 - E )-1  * X on LU4
             CALL REWINO(LU5)
             CALL REWINO(LU4)
             CALL REWINO(LUDIA)
*
             CALL H0M1TD(LU4,LUDIA,LU5,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPRDD(VEC1,VEC2,LU5,LU4,1,LBLK)
* is X an eigen vector for (H0 - 1 ) - 1
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU5,LU4,1,LBLK))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(lupri,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU5,LU3,1,LBLK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(lupri,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(lupri,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU3,LU4,LU5,1,LBLK)
                CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(lupri,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU3,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
           CALL REWINE( LU1 ,LBLK)
           DO 80 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             WORK(IVEC) = INPRDD(VEC1,VEC2,LU1,LU3,0,LBLK)
C?       WRITE(lupri,*) ' MICDV4 : Overlap ', WORK(IVEC)
   80      CONTINUE
*
           CALL REWINE(LU1,LBLK)
           DO 82 IVEC = 1,NVEC+IADD-1
             CALL REWINE(LU3,LBLK)
             CALL REWINE(LU4,LBLK)
             CALL VECSMD(VEC1,VEC2,-WORK(IVEC),1.0D0,LU1,LU3,
     &                   LU4,0,LBLK)
             CALL COPVCD(LU4,LU3,VEC1,1,LBLK)
   82      CONTINUE
           IF ( IPRT  .GE. 600 ) THEN
             WRITE(lupri,*) '   Orthogonalized (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU3,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU3,LBLK)
           CALL SCLVCD(LU3,LU1,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU3,LU4,FACTOR,VEC1,1,LBLK)
             WRITE(lupri,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : Optimal combination of new and old directions
*
*  2.1: Multiply new directions with matrix
      CALL SKPVCD(LU1,NVEC,VEC1,1,LBLK)
      CALL SKPVCD(LU2,NVEC,VEC1,1,LBLK)
      DO 150 IVEC = 1, IADD
        CALL REWINE(LU3,LBLK)
        CALL COPVCD(LU1,LU3,VEC1,0,LBLK)
        CALL SIGDEN_CI(VEC1,VEC2,C2,LU3,LU4,cdummy,sdummy,ISIGDEN,-1,
     &                 xxs2dum)
        CALL REWINE(LU4,LBLK)
        CALL COPVCD(LU4,LU2,VEC1,0,LBLK)
*. Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO 140 JVEC = 1, NVEC+IVEC
          CALL REWINE(LU4,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU4,0,LBLK)
  140   CONTINUE
  150 CONTINUE
*. Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN_luci(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(lupri,*) ' Var and Pert solution, energy and coefficients'
        WRITE(lupri,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(lupri,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(lupri,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(lupri,'(5F18.13)')
     &   ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
         WRITE(lupri,'(A)') ' Norm of Residuals (Previous it) '
         WRITE(lupri,'(5F18.13)')
     &   ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(lupri,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(lupri,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMAT(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset      
*
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU5,LBLK)
        DO 320 IROOT = 1, NROOT
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU3,LU3,1,LBLK)
          CALL REWINE(LU3,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
          WORK(IROOT) = SCALE
          CALL SCLVCD(LU3,LU5,SCALE,VEC1,0,LBLK)
  320   CONTINUE
*. Transfer C vectors to LU1
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 411 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU1,VEC1,0,LBLK)
  411   CONTINUE
*. corresponding sigma vectors
        CALL REWINE (LU5,LBLK)
        CALL REWINE (LU2,LBLK)
        DO 329 IROOT = 1, NROOT
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU3,LU4,VEC1,VEC2,NVEC,1,LBLK)
*
          CALL REWINE(LU3,LBLK)
          CALL SCLVCD(LU3,LU5,WORK(IROOT),VEC1,0,LBLK)
  329   CONTINUE
*
* Transfer HC's to LU2
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU5,LBLK)
        DO 400 IVEC = 1,NROOT
          CALL COPVCD(LU5,LU2,VEC1,0,LBLK)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
*
      END IF
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
 
* ( End of loop over iterations )
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(lupri,1170) MAXIT
 1170   FORMAT('   Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(lupri,1180) ITER
 1180   FORMAT(1X,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          WRITE(lupri,*)
          WRITE(lupri,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(lupri,*)
     &    '============================================'
          WRITE(lupri,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          root_residual(IROOT) = rnrm(iter,iroot)
          WRITE(lupri,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(lupri,1200)
 1200       FORMAT(/'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(lupri,1300)
 1300     FORMAT(/' Summary of iterations',
     +           /' ---------------------')
          WRITE(lupri,1310)
 1310     FORMAT
     &    (/' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(lupri,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           root_residual(IROOT) = rnrm(iter,iroot)
 1601   CONTINUE
      END IF
      nfinal_vec = nvec
*
      IF(IPRT .le. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(lupri,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
C
      RETURN
 1030 FORMAT(/2X,7F15.8,/,(2X,7F15.8))
 1120 FORMAT(/2X,I3,7F15.8,/,(5X,7F15.8))
      END

      SUBROUTINE CLASS_TRUNC(NCLS,ICLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                       E_CONV,ICLS_A,N_TRN_CLS,E_TRUNC,W_TRUNC,
     &                       IPRNT)
*
* Decide which classes of parameters that can be eliminated
*
* Jeppe Olsen, Jan 97
*              March '97 updated
*
*. Note in current version all energy contributions are
*. positive.
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS_L(NCLS)
      DIMENSION CLS_CT(NCLS),CLS_ET(NCLS),CLS_C(NCLS),CLS_E(NCLS)
*. Output
      INTEGER ICLS_A(NCLS)
*. Giving  the truncated classes
*  Additional output is
*              N_TRN_CLS : Number of truncated classes
*              E_TRUNC   : Estimated error of eliminating these classes
*                          (compared to expansion with only largest
*                           coefficients included).
*                          E_TRUNC is thus the error that should be
*                          added to the estmated error arising from
*                          eliminating small terms
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' Welcome to TRUNC_CLASS '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*)
     &  ' Required threshold for convergence of energy',E_CONV
      END IF
*
*. First find total energy correction and
*. largest class contribution
*
      E_TOT = 0.0D0
      E_CLS_MAX= 0.0D0
      DO JCLS = 1, NCLS
       E_TOT = E_TOT + CLS_ET(JCLS)
       E_CLS_MAX = MAX( E_CLS_MAX,CLS_ET(JCLS))
      END DO
C?    WRITE(6,*) ' E_CLS_MAX = ', E_CLS_MAX
*
*. The truncation is done in two steps, first an overall
*  threshold for deleting classes is constructed,
*  and all classes with energy contributions less than this
*  threshold are excluded.
*. Then one takes a pass through remaining classes and
*  eliminates untill the quota is filled
*
*. Threshold for including determinants, start with E_CLS_MAX
*. and decrease until quota is filled
*
      ILOOP = 0
      STEP = 1.2
 1000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
          E_THRES =  E_CLS_MAX
        ELSE
          E_THRES = E_THRES/STEP
        END IF
        E_ELI = 0.0D0
        N_ELI = 0
        DO JCLS = 1, NCLS
          IF(ABS(CLS_ET(JCLS)).LE.E_THRES) THEN
            E_ELI = E_ELI + CLS_ET(JCLS)
            N_ELI = N_ELI + 1
          END IF
        END DO
        IF(ILOOP.EQ.10000) THEN
          WRITE(6,*) ' Loop count exceeded 10000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          Call Abend2( 'CLASS_TRUNC: Iloop.eq.10000' )
        END If
      IF(E_ELI.GT. E_CONV) GOTO 1000
      IF(NTEST.GE.10) WRITE(6,*)
     & ' Overall Threshold for eliminating classes',E_THRES
*
      N_PAS_CLS = 0
      L_PAS_CLS = 0
      N_TRN_CLS = 0
      L_TRN_CLS = 0
      N_ACT_CLS = 0
      L_ACT_CLS = 0
*
      E_TRUNC = 0.0D0
      E_TRUNCT = 0.0D0
      W_TRUNCT= 0.0D0
      W_TRUNC = 0.0D0
*
*. Eliminate classes with energy contribution less than E_THRES
*
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).NE.0.0D0.AND
     &     .ABS(CLS_ET(JCLS)).LE. E_THRES) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        END IF
      END DO
      IF(NTEST.GE.10) THEN
      WRITE(6,*)
     &  ' Number of classes with contributions less than E_THRES',
     &   N_TRN_CLS
      WRITE(6,*) ' Energy contributions from these classes ',
     &   E_TRUNCT
      END IF
* Eliminate remaining classes until thres hold is obtained
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).EQ.0.0D0) THEN
*. Passive class, no contribution before truncation
          N_PAS_CLS = N_PAS_CLS + 1
          L_PAS_CLS = L_PAS_CLS + ICLS_L(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1
     &  .AND.ABS(E_TRUNCT+CLS_ET(JCLS)).LT.E_CONV) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1) THEN
*. Class is active
          N_ACT_CLS = N_ACT_CLS + 1
          L_ACT_CLS = L_ACT_CLS + ICLS_L(JCLS)
           ICLS_A(JCLS) = 1
        END IF
      END DO
*. Correct for missing minus in first order correction (not here, not now )
      E_TRUNC  =  E_TRUNC
      E_TRUNCT =  E_TRUNCT
*
      IF(NTEST.GE.10) THEN
      WRITE(6,*)
      WRITE(6,*) ' Estimated complete energy contribution of ',
     &            '  eliminated classes ', -E_TRUNCT
      WRITE(6,*) ' Estimated truncated energy contribution of ',
     &            '  eliminated classes ', -E_TRUNC
      WRITE(6,*) ' Truncation error in eliminated classes was ',
     &              -E_TRUNC+E_TRUNCT
      WRITE(6,*) ' Estimated truncated weight contribution of ',
     &            '  eliminated classes ', W_TRUNC
      WRITE(6,*) ' Estimated energy contribution without trunc',
     &             -E_TOT
      WRITE(6,*) ' Energy contribution of of active classes ',
     &             -E_TOT+E_TRUNCT
      END IF

C?    WRITE(6,*)
C?    WRITE(6,*) '  Class      Number     Dimension  '
C?    WRITE(6,*) ' ================================= '
C?    WRITE(6,'(A,5X,I5,5X,I10)') 'Passive  ', N_PAS_CLS,L_PAS_CLS
C?    WRITE(6,'(A,5X,I5,5X,I10)') 'Truncated', N_TRN_CLS,L_TRN_CLS
C?    WRITE(6,'(A,5X,I5,5X,I10)') 'Active   ', N_ACT_CLS,L_ACT_CLS
*
*. Eliminate classes with energy contributions less than E_TOT* FAC2
*  (temporary elimination, these classes may be invoked later in
*  the iterative sequence). Only the active classes are examined.
*
      FAC2 = 0.1
      ILOOP = 0
      STEP = 1.2
      E_TEMP_TRUNC =  ABS(E_TOT)*FAC2
 2000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
          E_THRES2 =   E_TEMP_TRUNC
        ELSE
          E_THRES2 = E_THRES2/STEP
        END IF
        E_ELI2 = 0.0D0
        N_ELI2 = 0
        DO JCLS = 1, NCLS
          IF(ICLS_A(JCLS).EQ.1.AND.
     &       ABS(CLS_ET(JCLS)).LE.E_THRES2) THEN
            E_ELI2 = E_ELI2 + ABS(CLS_ET(JCLS))
            N_ELI2 = N_ELI2 + 1
          END IF
        END DO
        IF(ILOOP.EQ.10000) THEN
          WRITE(6,*) ' Loop count exceeded 10000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          Call Abend2( 'CLASS_TRUNC: Iloop.eq.10000' )
        END If
      IF(E_ELI2.GT. E_TEMP_TRUNC) GOTO 2000
      IF(NTEST.GE.10) THEN
      WRITE(6,*)
     & ' Temporary elimination of classes with total contribution'
      WRITE(6,*) ' less than ',E_TEMP_TRUNC
      WRITE(6,*) ' gives threshold for temporary elimination',E_THRES2
      END IF
*
*. Eliminate classes with energy contribution less than E_THRES2
*
      N_TRN_CLS2 = 0
      L_TRN_CLS2 = 0
      E_TRUNC2 = 0
      W_TRUNC2 = 0
      DO JCLS = 1, NCLS
        IF(ICLS_A(JCLS).EQ.1.AND.CLS_CT(JCLS).NE.0.0D0.AND
     &     .ABS(CLS_ET(JCLS)).LE. E_THRES2) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          N_TRN_CLS2= N_TRN_CLS2+ 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          L_TRN_CLS2= L_TRN_CLS2+ ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNC2 = E_TRUNC2 + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNC2 = W_TRUNC2+ CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        END IF
      END DO
      IF(NTEST.GE.1) THEN
      WRITE(6,*)
      WRITE(6,'(A,F25.12)')
     &  ' Energy contributions from eliminated classes ', -E_TRUNC
      WRITE(6,'(A,F25.12)')
     &  ' Norm of eliminated classes                   ',SQRT(W_TRUNC)
      END IF
*
      N_ACT_CLS = N_ACT_CLS - N_TRN_CLS2
      L_ACT_CLS = L_ACT_CLS - L_TRN_CLS2
*
      E_TRUNC = - E_TRUNC
      E_TRUNCT = -  E_TRUNCT
*
      if (NTEST.ge.5) then
      WRITE(6,*)
      WRITE(6,*) '             Number     Dimension  '
      WRITE(6,*) ' ================================= '
      WRITE(6,'(A,5X,I5,5X,I10)') 'Passive  ', N_PAS_CLS,L_PAS_CLS
      WRITE(6,'(A,5X,I5,5X,I10)') 'Truncated', N_TRN_CLS,L_TRN_CLS
      WRITE(6,'(A,5X,I5,5X,I10)') '(Temp)   ',
     &N_TRN_CLS2,L_TRN_CLS2
      WRITE(6,'(A,5X,I5,5X,I10)') 'Active   ', N_ACT_CLS,L_ACT_CLS
*
      WRITE(6,*)
      WRITE(6,*) ' Information about classes '
      WRITE(6,*) ' =========================='
      WRITE(6,*)
      WRITE(6,*)
     & ' Class    Dimension      E          E(Trunc)      C',
     & '         C(Trunc)   Active'
      WRITE(6,*)
     & ' =============================================================',
     & '================'
       DO JCLS = 1, NCLS
         IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
     &   THEN
           WRITE(6,'(2X,I4,1X,I10,3X,4(E12.6,1X),1X,I2)')
     &     JCLS,ICLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
     &     CLS_CT(JCLS),CLS_C(JCLS),ICLS_A(JCLS)
         END IF
       END DO
      end if
*
      RETURN
      END
***********************************************************************

      SUBROUTINE CLS_TO_BASE(CLS_E,EBASC,CLS_C,CBASC,NCLS,NSPC,
     &                       IBASSPC,IPRNT)
*
* Class info => Base space info
* for energy abd wf correction
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBASSPC(NCLS)
      DIMENSION CLS_C(NCLS),CLS_E(NCLS)
*. Output
      DIMENSION EBASC(*),CBASC(*)
*
      ZERO = 0.0D0
      CALL SETVEC(EBASC,NSPC,ZERO)
      CALL SETVEC(CBASC,NSPC,ZERO)
*
      DO ICLS = 1, NCLS
        ISPC = IBASSPC(ICLS)
        EBASC(ISPC) = EBASC(ISPC) + CLS_E(ICLS)
        CBASC(ISPC) = CBASC(ISPC) + CLS_C(ICLS)
      END DO
*
      NTESTL = 000
      NTEST = max(IPRNT,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =============================================='
        WRITE(6,*) ' Contribution  to energy and wf per base space '
        WRITE(6,*) ' =============================================='
        WRITE(6,*)
        WRITE(6,'(A)') '  Class         Energy          wf '
        WRITE(6,'(A)') ' ==========================================='
        DO ISPC = 1, NSPC
          WRITE(6,'(2X,I3,3X,E12.6,2X,E12.6)')
     &          ISPC,EBASC(ISPC),CBASC(ISPC)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE DIAG_BLKS(A,X,LGAS,NOBPSM,MXPOBS,NSMOB,NGAS,
     &                     SCR1,SCR2)
*
* A matrix A and an orbital partitioning LGAS is given.
* Diagonalize diagonal blocks of A
*
*. A is assumed to contain only active orbitals and is assumed to be in
* packed form and symmetry ordered
*
*
* Input
* =====
* A : Input matrix
* LGAS : Orbital partitioning
* NOBPSM : Number of orbitals per symmetry
* MXPOBS : Max number of orbital symmetries
* NSMOB  : Number of orbital symmetries
* NGAS   : Number of orbital partitionings
*
* Output
* ======
* X : Eigenvector expansion, sorted according to eigenvalues in
*     each subspace
* A : Correspinding eigenvalues
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION A(*),X(*)
      DIMENSION LGAS(MXPOBS,*)
      DIMENSION NOBPSM(*)
*. Scratch : Number of orbitals **2 ( atmost )
      DIMENSION SCR1(*), SCR2(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DIAG_BLKS'
        WRITE(6,*) ' ========='
        WRITE(6,*) ' NSMOB NGAS ', NSMOB,NGAS
      END IF
*
      DO ISMOB = 1, NSMOB
*. offset for symmetryblocks in matrices
        IF(ISMOB.EQ.1) THEN
          IOFFMTP = 1
          IOFFMTC = 1
        ELSE
          IOFFMTP = IOFFMTP + NOBPSM(ISMOB-1)*(NOBPSM(ISMOB-1)+1)/2
          IOFFMTC = IOFFMTC + NOBPSM(ISMOB-1) ** 2
        END IF
*. Zero symmetry block of eigenvector matrix to avoid interactions
*. between different blocks
        LOBPS = NOBPSM(ISMOB)
        ZERO = 0.0D0
        CALL SETVEC(X(IOFFMTC),ZERO,LOBPS**2)
*. Loop over subbloks, extract,  diagonalize, and expand
        DO ITPOB = 1, NGAS
          IF(ITPOB.EQ.1) THEN
            IOFFOB=1
          ELSE
            IOFFOB = IOFFOB + LGAS(ISMOB,ITPOB-1)
          END IF
          LOB = LGAS(ISMOB,ITPOB)
*. Extract
          IJ2 = 0
          DO IOB = IOFFOB,IOFFOB+LOB-1
            JOBMX = IOB
            DO JOB = IOFFOB,JOBMX
              IJ1 = IOFFMTP -1 + IOB*(IOB-1)/2+JOB
              IJ2 = IJ2 + 1
              SCR1(IJ2) = A(IJ1)
            END DO
          END DO
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Extracted block of matrix for  ISMOB,ITPOB = ',
     &      ISMOB,ITPOB
            CALL PRSYM(SCR1,LOB)
          END IF
*, Diagonalize
C         CALL EIGEN_LUCI(WORK(KMAT1-1+IOFFP),WORK(KMAT2-1+IOFFC),LORB,0,1)
          CALL EIGEN_LUCI(SCR1, SCR2, LOB,0,1)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Corresponding eigenvalues and eigenvectors'
            WRITE(6,*) (SCR1(I*(I+1)/2),I=1, LOB)
            WRITE(6,*)
            CALL WRTMT_LU(SCR2,LOB,LOB,LOB,LOB)
          END IF
*. Expand eigenvector to full symmetry block
          IJ2 = 0
          DO JOB = IOFFOB,IOFFOB+LOB-1
            DO IOB = IOFFOB,IOFFOB+LOB-1
              IJ1 = IOFFMTC -1 + (JOB-1)*LOBPS + IOB
              IJ2 = IJ2 + 1
              X(IJ1) = SCR2(IJ2)
            END DO
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Complete eigenvector matrix'
        CALL APRBLM2(X,NOBPSM,NOBPSM,NSMOB,0)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE EKTPERT(F,S,NDIM,NORD,EN,C,
     &                   VEC1,VEC2,VEC3,AMAT1,AMAT2,AMAT3)
*
* Perturbation expansion of generalized eigenvalue problem
* Special version for EKT problem where there are singularities in
* the zero order matrices.
*
* Ordering the matrices so the occupied orbitals come first, and
* then the virtual orbitals, the  zero order matrices have the form
*
*     (x *                     )
*     ( x*               0     )
*     (**x                     )
*     (                        )
*     (                        )
*     (                        )
*     (                        )
*     (                        )
*     (  0              0      )
*     (                        )
*     (                        )
*     (                        )
*
*
*
* The matrices in the  eigenvalue problem FC = E SC
* are separated into orders
*
* F = Sum(k=0,NORD) F(K)
* S = Sum(K=0,NORD) S(K)
*
*. Obtain corrections to energy and wave functions
*
*. The normalization condition used is C(K)T S(0) C(0) = 0
*
* The energy corrections become
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
* and the wave function corrections in the occupied orbital space are
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(I=0,N-K)E(N-K-L)S(L)C(K))
*
* Whereas they read in the virtual space
*
* C(N-2) = (F(2)-E(0)S(2))-1 (-Sum(K=0,N-3)F(N-K)C(K)
*                             +Sum(K=0,N-3)Sum(I=0,N-K)E(N-K-L)S(L)C(K)
*                             + -(F(2)-E0*S(2)) C(N-2)OCC )
* Where only the
*
* The zero order matrices F(0),S(0) are assumed diagonal
*
* Jeppe and Dage, Oct 22 1997

*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION F(NDIM**2,*),S(NDIM**2,*)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
      DIMENSION AMAT1(NDIM,NDIM),AMAT2(NDIM,NDIM)
      DIMENSION AMAT3(NDIM,NDIM)
*
      NTEST = 1

*. Zero order wavefunction in virtual space
* obtained by solving - in the virtual space, the
* equations
*
* (F(2)-E(0)S(2))(virt,virt) C0(virt) = -(F(2)-E(0)S(2))C0(occ))(virt)
*
*
*. Obtain the number of occupied orbitals by examining S(0)
      NOCC = 0
      DO I = 1, NDIM
       IF(S((I-1)*NDIM+I,1) .NE. 0) NOCC = NOCC + 1
      END DO
      NVIRT = NDIM - NOCC
      WRITE(6,*) ' Number of occupied orbitals ', NOCC
      WRITE(6,*) ' Number of virtual  orbitals ', NVIRT
*
*. Zero order energy
*
C          MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(F,C,VEC1,NDIM,NDIM,0)
      C0FC0 = INPROD(VEC1,C,NDIM)
*
      CALL MATVCB(S,C,VEC1,NDIM,NDIM,0)
      C0SC0 = INPROD(VEC1,C,NDIM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'C(0)T F(0) C(0) = ', C0FC0
        WRITE(6,*) 'C(0)T S(0) C(0) = ', C0SC0
      END IF
      E0 = C0FC0/C0SC0
      EN(0) = E0
*. Save diagonal of F(0) - E(0)C(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = F((I-1)*NDIM+I,1)-E0*S((I-1)*NDIM+I,1)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Zero order diagonal '
        CALL WRTMT_LU(VEC3,1,NDIM,1,NDIM)
      END IF
*
*. The virtual-virtual  part of the operator F(2)-E(0)S(2)
*. will be used to solve linear equations. Set up inverse
*
      ONE = 1.0D0
      CALL VECSUM(AMAT1,F(1,3),S(1,3),ONE,-E0,NDIM**2)
      DO IOCC = 1, NOCC
        DO J = 1, NDIM
          AMAT1(IOCC,J) = 0.0D0
          AMAT1(J,IOCC) = 0.0D0
        END DO
        AMAT1(IOCC,IOCC) = 0.001
      END DO
C?    write(6,*) ' Input matrix to INVMAT'
C?    CALL WRTMT_LU(AMAT1,NDIM,NDIM,NDIM,NDIM)
      CALL COPVEC(AMAT1,AMAT3,NDIM**2)
      CALL INVMAT(AMAT1,AMAT2,NDIM,NDIM)
*     ^ Returns inverse matrix in AMAT1
      WRITE(6,*) ' Inverse matrix obtained by INVMAT'
      CALL WRTMT_LU(AMAT1,NDIM,NDIM,NDIM,NDIM)
*. Inverse matrix by diagonalization
      CALL COPVEC(AMAT3,AMAT1,NDIM**2)
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
      CALL INVERT_BY_DIAG(AMAT1,AMAT2,AMAT3,VEC1,NDIM)
*. Zero occ-occ part
      DO IOCC = 1, NOCC
        DO J = 1, NDIM
          AMAT1(IOCC,J) = 0.0D0
          AMAT1(J,IOCC) = 0.0D0
        END DO
      END DO
*. -(F(2)-E(0)S(2)) Will be used in the future, reconstruct in AMAT2
      CALL VECSUM(AMAT2,F(1,3),S(1,3),ONE,-E0,NDIM**2)
      ONEM = -1.0D0
      CALL SCALVE(AMAT2,ONEM,NDIM**2)
*. Obtain virtual part of C0
* (F(2)-E(0)S(2))(virt,virt) C0(virt) = -(F(2)-E(0)S(2))C0(occ))(virt)
      CALL MATVCB(AMAT2,C(1,1),VEC1,NDIM,NDIM,0)
      CALL MATVCB(AMAT1,VEC1,VEC2,NDIM,NDIM,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Virtual part of C0 '
        CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
      END IF
      CALL COPVEC(VEC2(1+NOCC),C(1+NOCC,1),NVIRT)
C?    WRITE(6,*) ' complete zero order correction vector '
C?    CALL WRTMT_LU(C(1,1),1,NDIM,1,NDIM)
*. And zero it once again ( COnstructed below )
      ZERO = 0.0D0
      CALL SETVEC(C(1+NOCC,1),ZERO,NVIRT)
*
*. And then start the iterations
      DO IORD = 1, NORD
*
* =================================================
* The (IORD-2) wf corrections in the virtual space
* =================================================
*
        IF(IORD.GE.2) THEN
          ZERO = 0.0D0
          CALL SETVEC(VEC2,ZERO,NDIM)
*. Note : Only the occupied part of C(N-2) is included in RHS of
*         expression. The virtual part was carefully zeroed !
          DO K = 0, IORD -2
            CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            ONE = 1.0D0
            ONEM = -1.0D0
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
          END DO
*
          DO K = 0, IORD -2
            DO L = 0, IORD -K
              CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
              CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
            END DO
          END DO
*. Multiply with (E(2)-E0S(2))-1
          CALL MATVCB(AMAT1,VEC2,VEC1,NDIM,NDIM,0)
*. And save
          CALL COPVEC(VEC1(1+NOCC),C(1+NOCC,IORD-2+1),NVIRT)
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' occ+virtual part for order = ',IORD-2
            CALL WRTMT_LU(C(1,IORD-2+1),1,NDIM,1,NDIM)
          END IF
*
        END IF
*       ^ End of construction of virtual part of C(VIRT,IORD-2)
*
*  =================
*. Energy correction
*  =================
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      - Sum(I=0,N-1)Sum(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
        EN(IORD) = 0.0D0
        DO I = 1, IORD
          CALL MATVCB(F(1,I+1),C(1,IORD-I+1),VEC1,NDIM,NDIM,0)
          EN(IORD) = EN(IORD)+INPROD(C,VEC1,NDIM)
        END DO
C?      write(6,*) ' First term to En ', EN(IORD)
        DO I = 0,IORD-1
          DO J = 1, IORD-I
            CALL MATVCB(S(1,J+1),C(1,I+1),VEC1,NDIM,NDIM,0)
            EN(IORD) = EN(IORD)-EN(IORD-I-J)*INPROD(C,VEC1,NDIM)
C?          write(6,*) ' second term to EN: I J EN ',I,J,EN(IORD)
          END DO
        END DO
        EN(IORD) = EN(IORD)/C0SC0
        WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ===========================================
*. Wave function corrections, occupied of IORD
*  ===========================================
*
* The occupied part obtained from
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(L=0,N-K)E(N-K-L)S(L)C(K))
        ZERO = 0.0D0
        CALL SETVEC(VEC2,ZERO,NDIM)
        DO K = 0, IORD -1
          CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
        END DO
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
*
        DO K = 0, IORD -1
          DO L = 0, IORD -K
            CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
C?          write(6,*) ' second term to rhs, K,L = ', K,L
C?          CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
          END DO
        END DO
*. Check overlap with S(0) times zero order state ( should be zero )
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
        CALL MATVCB(S(1,1),C(1,1),VEC1,NDIM,NDIM,0)
        C0SSC0 = INPROD(VEC1,VEC1,NDIM)
        OVLAP = INPROD(VEC1,VEC2,NDIM)
        WRITE(6,*) ' OVLAP = ', OVLAP
        FACTOR = -OVLAP/C0SSC0
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vector before ortho'
          CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
        END IF
        CALL VECSUM(VEC2,VEC2,VEC1  ,ONE,FACTOR,NDIM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vector before DIAVC2'
          CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
        END IF
*. Multiply with (F(0)-E(0)S(0))-1
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        JEPZAP = 1
        IF(IORD.EQ.1.AND.JEPZAP.EQ.1) THEN
          WRITE(6,*) ' First order correction zapped '
          ZERO = 0.0D0
          CALL SETVEC(VEC1,ZERO,NOCC)
        END IF
*
        CALL COPVEC(VEC1(1),C(1,IORD+1),NOCC)
*. The virtual part is still not known
        ZERO = 0.0D0
        CALL SETVEC(C(1+NOCC,IORD+1),ZERO,NVIRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Eigenfunction correction ', IORD
          CALL WRTMT_LU(C(1,IORD+1),1,NDIM,1,NDIM)
        END IF
      END DO
*
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order        Correction '
      WRITE(6,*) ' ==============================='
      DO IORD = 1, NORD
        WRITE(6,'(3X,I3,E25.13)')IORD,EN(IORD)
      END DO
*
      ETOT = E0
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,'(A,E25.13)') ' Zero-order energy ', E0
      WRITE(6,'(A,E25.13)') ' Sum(K=0,NORD) E(K) ', ETOT
*
      RETURN
      END
***********************************************************************

      SUBROUTINE H0DIAG(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,NROOT,
     &           EIGVAL,EIGVEC,SCR,NTESTG,ECORE)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           approximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Obtain the lowest NROOT eigenvectors
*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* NROOT : Number of roots to be obtained
*
* ======
* Output
* ======
* EIGVAL(*) : Converged eigen values
* EIGVEC(IROOT,*) : Complete eigenvector IROOT
*
* Note : The NROOT eigenpairs to be obtained are assumed
*        to be ' slightly ' perturbed eigensolutions
*        of PHP.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL HPQTVM
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*)
* Output
      DIMENSION EIGVAL(*),EIGVEC(NP1DM+NP2DM+NQDM,*)
* Scratch
      DIMENSION SCR(*)
*.SCR Should atleast be dimensioned ??????
      LOGICAL CONVER
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
* There are two routes :2 Iterative diagonalization of complete matrix
*                       1 complete diagonalizations of partitioned
*                         matrices.
*.
* ========
* Route 1 :
* ========
* The Q-space can be partitioned into the P -space
* to give the effective eigenvalue equation
*
* (PHP - PHQ  (QHQ-E)**-1 QHP ) VP = E VP
*
* This leads to a simple iterative scheme

*
* Note : Only NROOT = 1 tested
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
      CALL QENTER('H0DIA')
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0DIAG speaking '
        WRITE(6,*) ' =============== '
      END IF
*
C?    write(6,*) ' QHQ as delivered '
C?    call WRTMT_LU(QHQ,1,NQDM,1,NQDM)
C?    write(6,*) ' PHP as delivered '
C?    call PRSYM(PHP,NP1DM)
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
*. A bit of memory
*
      I12 = 2
      IF(I12.EQ.1. OR. I12. EQ.3 ) THEN
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. A PQ matrix
      KLPQ = KLFREE
      KLFREE = KLFREE + NPDM * NQDM
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
*
*. Initial eigenvalues
      CALL COPVEC(PHP,SCR(KLPP1),NPDM*(NPDM+1)/2)
      CALL EIGEN_LUCI(SCR(KLPP1),SCR(KLPP2),NPDM,0,1)
*. Extract eigenvalues
      CALL XTRCDI(SCR(KLPP1),EIGVAL,NROOT,1)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Initial set of eigenvalues '
        CALL WRTMT_LU(EIGVAL,1,NROOT,1,NROOT)
      END IF
*. Largest allowed number of iterations
      MAXIT = 5
      DO 1000 IROOT = 1, NROOT
        CONVER = .FALSE.
        EINI = EIGVAL(IROOT)
        DO 900 ITER = 1, MAXIT
          IF(NTEST.GE.2) WRITE(6,*) ' Info from iteration ', ITER
*. Current eigenvalue and eigenvector
          E = EIGVAL(IROOT)
          CALL COPVEC(SCR(KLPP2+(IROOT-1)*NPDM),SCR(KLV1),NPDM)
* ==============================
* HPP - PHQ (QHQ - E) **-1 * QHP
* ==============================
*. QHP in KLPQ
C         TRPMAT(XIN,NROW,NCOL,XOUT)
          CALL TRPMAT(PHQ,NP1DM,NQDM,SCR(KLPQ))
*.Multiply with (QHQ-E)**-1
          DO 30 IP1 = 1, NP1DM
C           DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
            IOFF = KLPQ + (IP1-1)*NQDM
            CALL DIAVC3_LUCI(SCR(IOFF),SCR(IOFF),QHQ,-E,NQDM,XDUMMY)
   30     CONTINUE
*. Multiply with PHQ
          CALL MATML4(SCR(KLPP1),PHQ,SCR(KLPQ),NP1DM,NP1DM,
     &                NP1DM,NQDM,NQDM,NP1DM,0)
*.
C?      write(6,*) ' PHQ (QHQ-E)-1 QHP matrix '
C?      CALL WRTMT_LU(SCR(KLPP1),NP1DM,NP1DM,NP1DM,NP1DM)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          SIGNTP = 1.0
          CALL TRIPAK_LUCI(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,SIGNTP)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL PRSYM(SCR(KLPP1),NPDM)
          END IF
*.Diagonalize
          CALL EIGEN_LUCI(SCR(KLPP1),SCR(KLPP2),NPDM,0,1)
*. Extract eigenvalues
           EIGVAL(IROOT) = SCR(KLPP1-1+IROOT*(IROOT+1)/2)
           IF(NTEST.GE.2)
     &     WRITE(6,*) ' Eigenvalue ', EIGVAL(IROOT)
           IF(NTEST.GE.10) THEN
             WRITE(6,*) ' P-space eigenvector '
             CALL WRTMT_LU(SCR(KLPP2+(IROOT-1)*NPDM),
     &            1,NPDM,1,NPDM)
           END IF
*. Check for convergence
           EVALDF = ABS(E-EIGVAL(IROOT))
           EVECOV = SQRT(INPROD(SCR(KLV1),SCR(KLPP2+(IROOT-1)*NPDM),
     &                          NPDM) )
           IF(EVALDF.LT.1.0D-7.AND.EVECOV.GT.0.999D0) CONVER = .TRUE.
           IF(CONVER) GOTO 901
  900    CONTINUE
  901    CONTINUE
*. P-part of eigenvector
         CALL COPVEC(SCR(KLPP2+(IROOT-1)*NPDM),EIGVEC(1,IROOT),
     &               NPDM)
*. Obtain Q part of eigenvector
*.    -(QHQ-E)**-1 HQP XP

         CALL MATML4(SCR(KLV1),PHQ,EIGVEC(1,IROOT),
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
*
C             DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3_LUCI(EIGVEC(NPDM+1,IROOT),SCR(KLV1),QHQ,
     &               -EIGVAL(IROOT),NQDM,XDUMMY)
         CALL SCALVE(EIGVEC(NPDM+1,IROOT),-1.0D0,NQDM)
*. Normalize
         XNORM = INPROD(EIGVEC(1,IROOT),EIGVEC(1,IROOT),NPQDM)
         SCALE = 1.0D0/SQRT(XNORM)
         CALL SCALVE(EIGVEC(1,IROOT),SCALE,NPQDM)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Initial and final eigenvalue ',
     &     EINI,EIGVAL(IROOT)
           WRITE(6,*)
     &   ' Part of eigenvector in Q space',
     &     SQRT(ABS(XNORM-1.0D0)/XNORM)
C?         WRITE(6,*) ' Eigenvector in PQ space '
C?         CALL WRTMT_LU(EIGVEC(1,IROOT),1,NPQDM,1,NPQDM)
         END IF
*
 1000 CONTINUE
       END IF
       IF( I12. EQ.2 .OR. I12 .EQ.3 ) THEN
*
*. Iterative scheme
*
*
*. Initial eigenvalues and eigenvectors
*
       KLPP1 = 1
       KLFREE = KLPP1 + NP1DM*(NP1DM+1)/2
       KLPP2 = KLFREE
       KLFREE = KLPP2 + NP1DM*NP1DM
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LU4 = 37
        LU5 = 38
        LUDIA = 39
        CALL COPVEC(PHP,SCR(KLPP1),NP1DM*(NP1DM+1)/2)
        CALL EIGEN_LUCI(SCR(KLPP1),SCR(KLPP2),NP1DM,0,1)
*. Extract eigenvalues and eigenvectors on LU1
        CALL XTRCDI(SCR(KLPP1),EIGVAL,NROOT,1)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Initial set of eigenvalues '
          CALL WRTMT_LU(EIGVAL,1,NROOT,1,NROOT)
        END IF
*. Eigenvectors on LU1
        CALL SETVEC(EIGVEC(1,1),0.0D0,NPQDM)
        CALL REWINE(LU1,-1)
        DO 510 IROOT = 1, NROOT
          CALL COPVEC(SCR(KLPP2+(IROOT-1)*NP1DM),EIGVEC(1,1),NP1DM)
          CALL TODSC_LUCI(EIGVEC(1,1),NPQDM,-1,LU1)
  510   CONTINUE
*
*. Iterations
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LUDIA)
*. Davidson CI diagonalization
        MINST = 1
        NBLK = 1
        INICI = -1
        MAXCIT = 15
        IPRTCI = MAX(NTEST-2,0)
        MXVCCI = 3 * NROOT
        SHIFT = 0.0D0
        IDIAG = 1
        ICISTR = 1
        IDIAG = 1
        THRES= 1.0D-8
* FIXME !!!
C  Only 33 arguments in this call, 70 are required.
C
c       CALL CIEIG5(HPQTVM,INICI,EIGVAL,SCR(KLV1),SCR(KLV2),
c    &            MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,NPQDM ,
c    &            NBLK,NROOT,MXVCCI,MAXCIT,LU1,IPRTCI,
c    &            DUMMY,0,DUMMY,IDUMMY,
c    &            0,0,0,DUMMY,ECORE,ICISTR,NPQDM,IDIAG,DUMMY,THRES)
C       CALL CIEIG5(HPQTVM,INICI,ENOT,SCR(KLV1),SCR(KLV2),
C    &            MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,NPQDM ,
C    &            NBLK,NROOT,MXVCCI,MAXCIT,LU1,IPRTCI,
C    &            DUMMY,0,DUMMY,IDUMMY,
C    &            0,0,0,DUMMY,0.0D0,1,NPQDM,1,DUMMY)
        CALL REWINE(LU1,-1)
        DO 1286 JROOT = 1, NROOT
         CALL FRMDSC_LUCI(EIGVEC(1,JROOT),NPQDM,-1,LU1,IMZERO,IAMPACK)
 1286   CONTINUE
      END IF
*
      CALL QEXIT('H0DIA')
*
      RETURN
      END
***********************************************************************

      SUBROUTINE H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,
     &                    IOCA,IOCB,NGAS,MXPNGAS,INTH0SPC,NELFTP)
*
* Set up INTH0SPC : Division of CI space, so only determinants
*                   belonging to the same space  have nonvanishing
*                   matrix elements of H0
*
* =====
* Input
* =====
*
* IH0SPC : ne. 0 : Interacting subspaces have been defined
*          .eq.0 : Interacting subspaces not defined, let
*                  evrything interact
* NPTSPC : Number of subspaces defined
* IOCPTSPC : Allowed occumulated occupation of each subspace
* NOCTPA :  Number of alpha occupation types
* NOCTPB : Number of beta occupation types
* IOCA : Occupation  of alpha string
* IOCB : Occupation  of beta string
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      DIMENSION IOCPTSPC(2,MXPNGAS,*)
      DIMENSION IOCA(MXPNGAS,*),IOCB(MXPNGAS,*)
      DIMENSION NELFTP(*)
*. Output
      DIMENSION INTH0SPC(NOCTPA,NOCTPB)
*
      IF(IH0SPC.EQ.0) THEN
*. All interactions allowed
        IONE = 1
        CALL ISETVC(INTH0SPC,IONE,NOCTPA*NOCTPB)
      ELSE
*. Explicit construction of matrix giving partitioning of
*  subspaces
        IZERO = 0
        CALL ISETVC(INTH0SPC,IZERO,NOCTPA*NOCTPB)
*
        DO ISPC = 1, NPTSPC
          DO IATP = 1, NOCTPA
            DO IBTP = 1, NOCTPB
              IAMOKAY = 1
              IEL = 0
C?            WRITE(6,*) ' ISPC IATP IBTP ', ISPC,IATP,IBTP
              DO IGAS = 1, NGAS
               IEL = IEL
     &             + NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
C?             WRITE(6,*) ' IGAS IEL ', IGAS,IEL
C?             WRITE(6,*)
C?   &          ' Limits :',IOCPTSPC(1,IGAS,ISPC),IOCPTSPC(2,IGAS,ISPC)
               IF(IEL.LT.IOCPTSPC(1,IGAS,ISPC).OR.
     &            IEL.GT.IOCPTSPC(2,IGAS,ISPC)    ) IAMOKAY = 0
              END DO
C?            WRITE(6,*) ' IAMOKAY = ', IAMOKAY
*. Allowed
              IF(IAMOKAY.EQ.1.AND.INTH0SPC(IATP,IBTP).EQ.0) THEN
                INTH0SPC(IATP,IBTP) = ISPC
              END IF
            END DO
          END DO
        END DO
      END IF
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Output from  H0INTSPC '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' Output from H0INTSPC '
        WRITE(6,*)
        CALL IWRTMA(INTH0SPC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS

*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVER
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL QENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3_LUCI(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          SIGNTP = 1.0
          CALL TRIPAK_LUCI(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,SIGNTP)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK_LUCI(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM,SIGNTP)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMT_LU(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
C               INVMAT(A,B,MATDIM,NDIM)
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================

         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3_LUCI(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LUDIA = 37
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LUDIA)
*. Initial Guess
        CALL REWINE(LU1,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LU1)
*. Right hand side
        CALL REWINE(LU2,-1)
        CALL TODSC_LUCI(RHS,NPQDM,-1,LU2)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM,LU1,LU2,LU3,LUDIA,SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(SCR(KLV1),NPQDM,-1,LU1,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      CALL QEXIT('H0LNS')
      RETURN
      END
***********************************************************************

      SUBROUTINE H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
* Disc version
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  LUOUT       LUDIA        LUIN
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal on file LUDIA
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements on LUDIA corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , September 1993
*
*
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0
* LUIN  : File Containing input vector
* LBLK : Defines format of files
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
COLD  REAL * 8  INPROD
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
!     work is h0scr in calling routine
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10)
     & write(6,*) ' H0M1TV , NPQDM = ', NPQDM
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
*. Obtain subspace components of input vector
C            GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,NTESTG)
        IZERO = 0
        CALL GATVCD(LUIN,LBLK,NPQDM,IPNTR,WORK(KLV1),VEC1,
     &              IZERO)
*. Solve linear equations in subspace
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
      END IF
*
*. Calculate inverse diagonal and scatter results from subspace,
*. Write to file LUOUT
C     DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
C    &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
      CALL DMTVDS(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &            IPNTR,WORK(KLV2),NPQDM,LBLK,XH0PSX)
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TD '
        WRITE(6,*) ' ========================= '
*. Note : works only if result vector is first file on LUOUT
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
        WRITE(6,*) ' Overlap between input and output vector',
     &               XH0PSX
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal Diag
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements in DIAG corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , May 1990

*
*
* =====
* Input
* =====
* DIAG : Diagonal of matrix
* VECIN : Input vector
* NVAR : Dimension of full matrix
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
* ======
* Output
* ======
* VECUT : Output vector (you guessed ?? ), can occupy same space
*         as VECIN or DIAG
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      REAL * 8  INPROD
*
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
!               work is h0scr in calling routine
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
C?    write(6,*) ' H0M1TV , NPQDM = ', NPQDM
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
        CALL GATVEC(WORK(KLV1),VECIN,IPNTR,NPQDM)
* X(T)(DIAG+SHIFT)-1 X in subspace, for later subtraction
        CALL GATVEC(WORK(KLV2),DIAG,IPNTR,NPQDM)
        CALL DIAVC3_LUCI(WORK(KLV2),WORK(KLV1),
     &       WORK(KLV2),SHIFT,NPQDM,X1)
       ELSE
         X1 = 0.0D0
       END IF
*
      CALL DIAVC3_LUCI(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
*
      IF(NPQDM .NE. 0 ) THEN
C                H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
C    &           X,RHS,S,SCR)
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
C?     write(6,*) ' KLPHP KLPHQ KLQHQ ',KLPHP,KLPHQ,KLQHQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
         X3 = INPROD(WORK(KLV1),WORK(KLV2),NPQDM)
         CALL SCAVEC(VECUT,WORK(KLV2),IPNTR,NPQDM)
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(6,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3


*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TV '
        WRITE(6,*) ' ========================= '
        CALL WRTMT_LU(VECUT,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE HINTV(LURHS,LUX,SHIFT,SHIFT_DIA,VECIN,VECOUT,
     &                LBLK,LUPROJ,LUPROJ2)
*
* Solve  (H+Shift)X = RHS
*
* Where H is matrix rep of operator defined by /OPER/ in space defined by
* /CANDS/
*
* If ICISTR.EQ.1 VECIN contains RHS, else RHS is assumed  on LURHS
* Output : solution is on LUX
*
* Jeppe Olsen, Winter of 1996
*
* Jan. 98 : SHIFT_DIA added
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLCLBTx, KLCLEBTx, KLCI1BTx, KLCIBTx, KLCIOIO, KLCBLTPx
!               for addressing of WORK
      REAL*8 INPROD , INPRDD
      LOGICAL CONVER
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
#include "mxpdim.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
#include "wrkspc.inc"
#include "clunit.inc"
* SCRATCH files used : LUSC3,LUSC34,LUSC35,LUSC37
#include "oper.inc"
#include "crun.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "stinf.inc"
#include "csm.inc"
*. Max number of iterations is picked from MAXCIT in crun
      EXTERNAL H0TVM
      DIMENSION ERROR(100)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HINTV ')
*
      NTEST = 0
*
*
* 1 : Construct diagonal  on LUSC3
* ================================
*
*. Use type of H0 as type of zero order operator  in all spaces
       IF(IAPR.NE.0) THEN
*. Not exact Hamiltonian
         IF(IH0INSPC(1).EQ.1.OR.IH0INSPC(1).EQ.3) THEN
*. Mp operator
           I12 = 1
           CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
           CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
         ELSE IF (IH0INSPC(1).EQ.2.OR.IH0INSPC(1).EQ.4) THEN
*. EN diagonal
           I12 = 2
         END IF
       END IF
       ECOREX = SHIFT_DIA
*
* Partitioning and blockstructure of CI vector
*
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NTTS = MXNTTS
      CALL MEMMAN(KLCLBTx ,  NTTS,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLCLEBTx,  NTTS,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLCI1BTx,  NTTS,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLCIBTx ,8*NTTS,'ADDL  ',1,'CIBT  ')
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ISSPC,WORK(KLCIOIO))
*
      CALL MEMMAN(KLCBLTPx,NSMST,'ADDL  ',2,'CBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTPx),0)
*
*. Batches  of C vector
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KLCBLTPx),WORK(KNSTSO(IATP)),
     &              WORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
     &              ISMOST(1,ISSM),NBATCH,WORK(KLCLBTx),
     &              WORK(KLCLEBTx),WORK(KLCI1BTx),WORK(KLCIBTx),0,
     &              ITTSS_ORD)
*. Number of BLOCKS
        NBLOCK = IFRMR(WORK(KLCI1BTx),1,NBATCH)
     &         + IFRMR(WORK(KLCLBTx),1,NBATCH) - 1
C?      WRITE(6,*) ' HINTV : Number of blocks ', NBLOCK
        CALL GASDIAT(VECIN,LUSC3,ECOREX,ICISTR,I12,
     &               WORK(KLCBLTPx),NBLOCK,WORK(KLCIBTx),IDUMMY)
C      CALL GASDIAT(VECIN,ISSM,ISSPC,LUSC3,ECOREX,ICISTR,I12)
*. Clean up time
       IF(IH0INSPC(1).EQ.1.OR.IH0INSPC(1).EQ.3) THEN
*. MP operator
         CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
         CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
       END IF
*
* 2 : Solve linear set of equations
* ==================================
*
      ZERO = 0.0D0
      IF(LBLK.GT.0 ) THEN
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        Call Abend2(' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1')
*. Two vectors can be in core
*. Initial Guess on LUX
        NDIM = LBLK
        CALL SETVEC(VECOUT,ZERO,NDIM)
        CALL REWINE(LUX,-1)
        CALL TODSC_LUCI(VECOUT,NDIM,-1,LUX)
*. Right hand side on LUSC34
        CALL REWINE(LUSC34,-1)
        CALL TODSC_LUCI(VECIN,NDIM,-1,LUSC34)
*. Max number of its and convergence thresholds are picked up from
* corresponding eigenvalue info
        CONVER = .FALSE.
        TEST = SQRT(THRES_E) * SQRT(INPROD(VECIN,VECIN,NDIM))
        ILNPRT = MAX(NTEST-10,0)
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXITLE = ',MXITLE
        CALL MINGCG(MV8,LUX,LUSC34,LUSC35,LUSC3,VECIN,VECOUT,
     &              MXIT_LOC,CONVER,TEST,SHIFT,ERROR,NDIM,
     &              LUPROJ,LUPROJ2,ILNPRT)
        CALL REWINE(LUX,-1)
        CALL FRMDSC_LUCI(VECOUT,NDIM,-1,LUX,IMZERO,IAMPACK)
*
         IF(NTEST.GE.5) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(VECOUT,1,NDIM,1,NDIM)
         END IF
*
      ELSE IF(LBLK.LE.0)   THEN
*
*. Use path allowing us to work with segments of vectors
*
*. Initial guess on LUX
        CALL SETVCD(LUSC3,LUX,VECOUT,ZERO,1,LBLK)
*. (Right hand side vector is assumed in place)
*. convergence threshold is picked up from
* corresponding eigenvalue info
        TEST =
     &  SQRT(THRES_E) * SQRT(INPRDD(VECIN,VECOUT,LURHS,LURHS,1,-1))
        ILNPRT = NTEST
        SHIFT2 = 0.0D0
        CONVER = .FALSE.
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXIT_LOC ',MXIT_LOC
        CALL MICGCG(H0TVM,LUX,LURHS,LUSC34,LUSC35,LUSC37,LUSC3,
     &              VECIN,VECOUT,MXIT_LOC,
     &              CONVER,TEST,SHIFT2,ERROR,NDIM,LUPROJ,LUPROJ2,
     &              ILNPRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Solution to linear set of Equations '
          CALL WRTVCD(VECIN,LUX,1,LBLK)
        END IF
*
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HINTV ')
      RETURN
      END
***********************************************************************

      SUBROUTINE HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                   IABSPC,JABSPC,IABOP,JABOP,IIF,JDOH2,
     &                   IDOH2,IMZERO,IDIAG)
*
* Decide upon the treatment of matrix element
*
* <IASM IATP IBSM IBTP | H(apr) | JASM, JATP JBSM, JBTP>
*
* and do preparations (IIF = 1 )
* or counteract preparations (IIF = 2)
*
* Jeppe Olsen, The last day of January 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      IF((IABSPC.NE.JABSPC).OR.
     &   (IABSPC.EQ.JABSPC.AND.(IABOP.EQ.1.OR.IABOP.EQ.2).AND.
     &    (JASM.NE.IASM.OR.IATP.NE.JATP.OR.IBTP.NE.JBTP))) THEN
*. Zero
        IMZERO = 1
      ELSE
*. Not zero !
        IMZERO = 0
      END IF
*. Diagonal approximation?
        IF(IABOP.EQ.1.OR.IABOP.EQ.2) THEN
          IDIAG = 1
        ELSE
          IDIAG = 0
        END IF
*. Moller Plesset or normal operator ?
        IF(IABOP.EQ.1.OR.IABOP.EQ.3.OR.IABOP.EQ.5)THEN
          IMP = 1
        ELSE
          IMP = 0
        END IF
C     END IF
*. Two - or one- electron operator
      IF(IABOP.EQ.1.OR.IABOP.EQ.3) THEN
        IDOH2 = 0
      ELSE
C       IDOH2 = JDOH2
        IDOH2 = 1
      END IF
*. Put MP integrals in place ( Or put good old one-electron integrals
*. back where they belong
      IF(IMP.EQ.1) THEN
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
         WRITE(6,*) ' HMATAPR speaking '
         WRITE(6,*) ' ================='
         WRITE(6,*) ' | IASM IATP IBSM IBTP > :',
     &              '|',IASM,IATP,IBSM,IBTP,'>'
         WRITE(6,*) ' | JASM JATP JBSM JBTP > :',
     &              '|',JASM,JATP,JBSM,JBTP,'>'
         WRITE(6,*) ' IABSPC,JABSPC :', IABSPC,JABSPC
         WRITE(6,*) ' IABOP ', IABOP
         WRITE(6,*) ' Results : IMP IDIAG IDOH2 IMZERO : ',
     &                IMP,IDIAG,IDOH2, IMZERO
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE HTV(VECIN,VECOUT,LUIN,LUOUT)
*
* Full operator times vector
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "oper.inc"
*
* complete operator in action
      IPERTOP = 0
      I12     = 2
      ISIGDEN = 1
*
      call quit('*** HTV is disabled in this lucita version. ***')
      CALL SIGDEN_CI(VECIN,VECOUT,C2,LUIN,LUOUT,cdummy,sdummy,
     &               ISIGDEN,-1,xxs2dum)
*
      RETURN
      END
***********************************************************************

      SUBROUTINE INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
*
* Invert symmetric  - hopefully nonsingular - matrix A
* by diagonalization
*
* Jeppe Olsen, Oct 97 to check INVMAT
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output matrix
      DIMENSION A(*)
*. Scratch matrices and vector
      DIMENSION B(*),SCR(*),VEC(*)
*
      NTEST = 01
*. Reform a to symmetric packed form
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
      SIGNTP = 1.0
      CALL TRIPAK_LUCI(A,SCR,1,NDIM,NDIM,SIGNTP)
*. Diagonalize
      CALL EIGEN_LUCI(SCR,B,NDIM,0,1)
      CALL COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMT_LU(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements - without safety at the moment
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode inactivated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL WRTMT_LU(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE MATPERT(H0,V,NDIM,NORD,EN,C,VEC1,VEC2,VEC3,ECORE)
*
* Perturbation expansion of simple  eigenvalue problem
*
* Explicit matrix version
*
*
*. Obtain corrections to energy and wawe functions
*
*. The normalization condition used is C(K)T  C(0) = 0
*
* The energy corrections are
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
*
*
* Jeppe Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION H0(NDIM**2),V(NDIM**2)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
*
*. Zero order energy
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(H0,C,VEC1,NDIM,NDIM,0)
      E0   = INPROD(VEC1,C,NDIM)
*
      WRITE(6,*) 'E0  = ', E0
      EN(0) = E0
*. Save diagonal of H0 - E(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = H0((I-1)*NDIM+I)-E0
      END DO
C?    WRITE(6,*) ' Zero order diagonal '
C?    CALL WRTMT_LU(VEC3,1,NDIM,1,NDIM)
*. And then start the iterations
      DO IORD = 1, NORD
*
*  =================
*. Energy correction
*  =================
*
* E(n) =  C(0)T V C(N-1)
        CALL MATVCB(V,C(1,IORD+1-1),VEC1,NDIM,NDIM,0)
        EN(IORD) = INPROD(C,VEC1,NDIM)
C?      WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ==========================
*. Wave function corrections
*  ==========================
*
* C(N) = (H(0)-E(0))-1 (-VC(N-1)
*                           +Sum(K=1,N)E(K)C(N-K))
        CALL MATVCB(V,C(1,IORD+1-1),VEC2,NDIM,NDIM,0)
        ONEM = -1.0D0
        CALL SCALVE(VEC2,ONEM,NDIM)
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
*
        ONE = 1.0D0
        DO K = 1, IORD
          CALL VECSUM(VEC2,VEC2,C(1,IORD+1-K),ONE,EN(K),NDIM)
        END DO
*. Check overlap with zero order state ( should be zero )
        OVLAP = INPROD(C(1,1),VEC2,NDIM)
        FACTOR = -OVLAP
C?      WRITE(6,*) ' OVLAP = ',OVLAP
        CALL VECSUM(VEC2,VEC2,C(1,1),ONE,FACTOR,NDIM)
*. Multiply with (H0(0)-E(0))-1
C            DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        ZERO = 0.0D0
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        CALL COPVEC(VEC1,C(1,IORD+1),NDIM)
*
C?      WRITE(6,*) ' Eigenfunction correction ', IORD
C?      CALL WRTMT_LU(C(1,IORD+1),1,NDIM,1,NDIM)
      END DO
*
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order             Correction '
      WRITE(6,*) ' ===================================='
      DO IORD = 1, NORD
        WRITE(6,'(3X,I3,E20.8)')IORD,EN(IORD)
      END DO
*
      ETOT = E0 + ECORE
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,*) ' Zero-order energy ', E0 + ECORE
      WRITE(6,*) ' Sum(K=0,NORD) E(K) ', ETOT
*
      END
***********************************************************************

      SUBROUTINE MICGCG(MV8,LU1,LU2,LU3,LU4,LU5,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,LUPROJ2,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
* The vector used to orthogonalize is saved on LUPROJ2
*
* Version using blocks of vectors
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPRDD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CALL QENTER('MICGC')
      NTEST = 005
      NTEST = MAX(NTEST,IPRT)
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' Welcome to MICGCG '
        WRITE(6,*) ' =================='
        WRITE(6,*)
*
C?    WRITE(6,*) ' NTEST ,LU1,LU2,LU3 = ', NTEST,LU1,LU2,LU3
      END IF
      CONVER = .FALSE.
      ITER = 1
*
      LBLK = -1
*
      ONE = 1.0D0
      ONEM = -1.0D0
      ZERO = 0.0D0
*. Overlap between LUPROJ and LUPROJ2
      IF(LUPROJ.GT.0) THEN
        X12 = INPRDD(VEC1,VEC2,LUPROJ,LUPROJ2,1,LBLK)
      ELSE
        X12 = 0.0D0
      END IF
C?    WRITE(6,*) ' MICGCG : X12 = ', X12
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X on LU2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector on LU1 '
        CALL WRTVCD(VEC1,LU1,1,LBLK)
        WRITE(6,*) ' Vector on LU2 '
        CALL WRTVCD(VEC1,LU2,1,LBLK)
      END IF
      CALL MV8(VEC1,VEC2,LU1,LU3)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector on LU3 '
        CALL WRTVCD(VEC1,LU3,1,LBLK)
      END IF
*
C          VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LU3,LU4,1,LBLK)
      CALL COPVCD(LU4,LU2,VEC1,1,LBLK)
*
*
      RNORM = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
      ERROR(1) = SQRT(RNORM)
      IF(ERROR(1).LE.TEST) THEN
*. Convergence in one shot- you are lucky or have
* supplied a vaninshing RHS
        NITER = 0
        CONVER = .TRUE.
        GOTO 1001
      END IF
*
*. Preconditioner H times initial residual, H * R on LU4
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Diagonal and input to diagonal '
       CALL WRTVCD(VEC1,LUDIA,1,LBLK)
       CALL WRTVCD(VEC1,LU2  ,1,LBLK)
       WRITE(6,*) ' SHIFT = ', W
      END IF
      CALL DMTVCD(VEC1,VEC2,LUDIA,LU2,LU4,W,1,1,LBLK)
      IF(LUPROJ.NE.0) THEN
        OVLAP = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
        FACTOR = -OVLAP/X12
        CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU4,LUPROJ2,LU3,1,LBLK)
        CALL COPVCD(LU3,LU4,VEC1,1,LBLK)
        OVLAP2 = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
        WRITE(6,*) ' Updated overlap of trial vector ', OVLAP2
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPRDD(VEC1,VEC2,LU2,LU4,1,LBLK)
*. P = RHO * H*R on LU3
      RHO = 1.0D0
      CALL SCLVCD(LU4,LU3,RHO,VEC1,1,LBLK)
*.S = AP on LU4
      CALL MV8(VEC1,VEC2,LU3,LU4)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*
* Vectors on files :
*     X on LU1
*     R on LU2
*     P on LU3
*  S=AP on LU4
*     H on LUDIA

        NITER = NITER + 1
       IF ( NTEST .GE. 10 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPRDD(VEC1,VEC2,LU3,LU4,1,LBLK)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S on LU2
        CALL VECSMD(VEC1,VEC2,ONE,-A,LU2,LU4,LU5,1,LBLK)
        CALL COPVCD(LU5,LU2,VEC1,1,LBLK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Residual on LU2 '
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*.    new residual has been obtained , check for convergence
        RNORM = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
        RNORME = MAX(RNORM,0.0D0)
        ERROR(ITER+1) = SQRT(RNORME)
*.    X = X + A * P
C?      WRITE(6,*) ' MICGCG : A = ', A
        CALL VECSMD(VEC1,VEC2,ONE,A,LU1,LU3,LU5,1,LBLK)
        CALL COPVCD(LU5,LU1,VEC1,1,LBLK)
*
        IF( SQRT(RNORME) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.H * R on LU4
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Diagonal and input to diagonal '
             CALL WRTVCD(VEC1,LUDIA,1,LBLK)
             CALL WRTVCD(VEC1,LU2  ,1,LBLK)
             WRITE(6,*) ' SHIFT = ', W
           END IF
*
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU2,LU4,W,1,1,LBLK)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Preconditioner times residual '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
           END IF
           IF(LUPROJ.NE.0) THEN
             OVLAP = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
             FACTOR = -OVLAP/X12
             CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU4,LUPROJ2,LU5,1,LBLK)
             CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
             OVLAP2 = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
C?           WRITE(6,*) ' Updated overlap of trial vector ', OVLAP2
*. Overlap between X and LUPROJ
             OVLAP3 = INPRDD(VEC1,VEC2,LUPROJ,LU1,1,LBLK)
             WRITE(6,*) ' Overlap between LU1 and LUPROJ ', OVLAP3
           END IF
*. GAMMA = <R!H!R>
           GAMMA = INPRDD(VEC1,VEC2,LU2,LU4,1,LBLK)
           B = GAMMA/C
*. P = RHO*(H*R + B*P) on LU3
           RHO = 1.0D0
           CALL VECSMD(VEC1,VEC2,ONE,B,LU4,LU3,LU5,1,LBLK)
           CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*.S = AP on LU4
           CALL MV8(VEC1,VEC2,LU3,LU4)
*.End of prepations for next iteration
        END IF
*
 1000 CONTINUE
 1001 CONTINUE
      IF(NTEST .GT. 0 ) THEN
*
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(/'  convergence was obtained in...',I3,' iterations',
     +        /'  norm of residual..............',1P,E15.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT+1)
 1020  FORMAT(/' convergence was not obtained in',I3,'iterations',
     +        /' norm of residual...............',1P,E15.8)
      END IF
*
      END IF
*
      IF(NTEST.GE. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(/' solution to set of linear equations')
       CALL WRTVCD(VEC1,LU1,1,LBLK)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call WRTMT_LU(vec2,1,nvar,1,nvar)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(/10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
*
      CALL QEXIT('MICGC')
      RETURN
      END
***********************************************************************

      SUBROUTINE MINDV4(MV7,
     &                  VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &                  IOLSEN,IPICO)
*
* Davidson algorithm , requires two vectors in core
* Multi root version
*
* Allows updating of preconditioning matrix so this is
* the current eigenvector approximation
* is an eigenvector for the preconditioner
*
* Jeppe Olsen Sept 89
*             Jan  92 : MV7 entry
*
* Input :
* =======
*        MV7 : Name of routine performing matrix*vector calculation
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    complete vector
*        LU2,LU3   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
* On input LU1 is supposed to hold initial guess to eigenvectors
*
* IOLSEN : Use inverse iteration modified Davidson
* IPICO  : Use perturbation estimate of new vector instead of
*          variational method
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
*
       EXTERNAL MV7
*
       WRITE(6,*) ':::::::::::::::::::'
       WRITE(6,*) '  Entering MINDV4  '
       WRITE(6,*) ':::::::::::::::::::'
       IOLSTM = IOLSEN
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Variational'
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Perturbational'
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Normal Davidson, Variational '
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Normal Davidson, Perturbational'
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' SORRY MINDV2 WOUNDED , MAXVEC .LT. 2*NROOT '
         Call Abend2( ' ENFORCED STOP IN MINDV2' )
       END IF
*
       IF(IPICO.NE.0) THEN
         MAXVEC = 2*NROOT
       END IF
*
       ISIGDEN = 1
       KAPROJ  = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-7
       CONVER = .FALSE.
       DO 1234 MACRO = 1,1
*
*.   INITAL ITERATION
       ITER = 1
       CALL REWINE( LU1 ,-1)
       CALL REWINE( LU2 ,-1)
       DO 10 IVEC = 1,NINVEC
         CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL quit('*** MINDV4 is disabled in this lucita version. ***')
         CALL SIGDEN_CI(VEC1,VEC2,C2,0,0,cdummy,sdummy,ISIGDEN,-1,
     &                  xxs2dum)
         CALL TODSC_LUCI(VEC2,NVAR,-1  ,LU2)
*        PROJECTED MATRIX
         CALL REWINE( LU1,-1)
         DO 8 JVEC = 1, IVEC
           CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*  DIAGONALIZE INITIAL PROJECTED MATRIX
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF( IPRT  .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Initial set of eigenvalues '
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF),IROOT=1,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
** LOOP OVER ITERATIONS
*
 1000 CONTINUE
      IF(IPRT  .GE. 10 ) THEN
       WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
      END IF


        ITER = ITER + 1
*
** 1          NEW DIRECTION TO BE INCLUDED
*
*   1.1 : R = H*X - EIGAPR*X
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         CALL SETVEC(VEC1,0.0D0,NVAR)
*
         CALL REWINE( LU2,-1)
         DO 60 IVEC = 1, NVEC
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   60    CONTINUE
*
         EIGAPR = EIG(ITER-1,IROOT)
         CALL REWINE( LU1,-1)
         DO 50 IVEC = 1, NVEC
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   50    CONTINUE
           IF ( IPRT  .GE.600 ) THEN
             WRITE(6,*) '  ( HX - EX ) '
             CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
           END IF
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
         RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
         RNRM(ITER-1,IROOT) = RNORM
         IF(RNORM.LT. TEST ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
*.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
           CALL REWINE( LUDIA,-1)
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL H0M1TV(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
              CALL REWINE(LU3,-1)
              CALL TODSC_LUCI(VEC1,NVAR,-1,LU3)
* Current eigen vector
              CALL REWINE( LU1,-1)
              CALL SETVEC(VEC1,0.0D0,NVAR)
              DO 59 IVEC = 1, NVEC
                CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
                FACTOR = AVEC((IROOT-1)*NVEC+ IVEC)
                CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   59         CONTINUE
              IF ( IPRT  .GE. 600 ) THEN
                WRITE(6,*) ' And X  '
                CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
              END IF
              CALL TODSC_LUCI(VEC1,NVAR,-1,LU3)
* (H0 - E )-1  * X
              CALL REWINE( LUDIA,-1)
              CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
              CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
     &                   H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
              CALL TODSC_LUCI(VEC2,NVAR,-1,LU3)
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPROD(VEC2,VEC1,NVAR)
* is X an eigen vector for (H0 - 1 ) - 1
              CALL VECSUM(VEC2,VEC1,VEC2,GAMMA,-1.0D0,NVAR)
              VNORM = SQRT(MAX(0.0D0,INPROD(VEC2,VEC2,NVAR)))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                CALL REWINE(LU3,-1)
                CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU3,IMZERO,IAMPACK)
                DELTA = INPROD(VEC1,VEC2,NVAR)
                CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSUM(VEC1,VEC1,VEC2,FACTOR,1.0D0,NVAR)
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) '  Modified new trial vector '
                  CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
                END IF
              ELSE
                IF(IPRT.GT.0) WRITE(6,*)
     &          ' Inverse correction switched of'
                CALL REWINE(LU3,-1)
                CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
              END IF
            END IF
*. 1.3 ORTHOGONALIZE TO ALL PREVIOUS VECTORS
           XNRMI =    INPROD(VEC1,VEC1,NVAR)
           CALL REWINE( LU1 ,-1)

           DO 80 IVEC = 1,NVEC+IADD-1
             CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
   80      CONTINUE
*. 1.4 Normalize vector and check for linear dependency
           SCALE = INPROD(VEC1,VEC1,NVAR)
           IF(ABS(SCALE)/XNRMI .LT. 1.0D-10) THEN
*. Linear dependency
             IADD = IADD - 1
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(6,*) '  Trial vector linear dependent so OUT !!! '
             END IF
           ELSE
             C1NRM = SQRT(SCALE)
             FACTOR = 1.0D0/SQRT(SCALE)
             CALL SCALVE(VEC1,FACTOR,NVAR)
*
             CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU1)
             IF ( IPRT  .GE.600 ) THEN
               WRITE(6,*) 'ORTHONORMALIZED (D-E)-1 *( HX - EX ) '
               CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
             END IF
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
*
*  2.1: MULTIPLY NEW DIRECTION WITH MATRIX
       CALL REWINE( LU1,-1)
       CALL REWINE( LU2,-1)
       DO 110 IVEC = 1, NVEC
         CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
         CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
  110  CONTINUE
*
      DO 150 IVEC = 1, IADD
        CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        call quit('*** MINDV4 is disabled in this lucita version. ***')
        CALL SIGDEN_CI(VEC1,VEC2,C2,0,0,cdummy,sdummy,ISIGDEN,-1,
     &                  xxs2dum)
        CALL TODSC_LUCI(VEC2,NVAR,-1  ,LU2)
*   AUGMENT PROJECTED MATRIX
        CALL REWINE( LU1,-1)
        DO 140 JVEC = 1, NVEC+IVEC
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
  140   CONTINUE
  150 CONTINUE
*  DIAGONALIZE PROJECTED MATRIX
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' PROJECTED MATRIX AND EIGEN PAIRS '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMT_LU(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  PERHAPS RESET OR ASSEMBLE CONVERGED EIGENVECTORS
*
  901 CONTINUE
*
      IPULAY = 1
      IF(IPULAY.EQ.1 .AND. MAXVEC.EQ.3 .AND.NVEC.GE.2.
     &   .AND. .NOT.CONVER) THEN
* Save trial vectors : 1 -- current trial vector
*                      2 -- previous trial vector orthogonalized
        CALL REWINE( LU3,-1)
        CALL REWINE( LU1,-1)
*. Current trial vector
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2200 IVEC = 1, NVEC
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
         CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2200   CONTINUE
        SCALE = INPROD(VEC1,VEC1,NVAR)
        SCALE  = 1.0D0/SQRT(SCALE)
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU3)
* Previous trial vector orthonormalized
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        SCALE2 = INPROD(VEC2,VEC2,NVAR)
        SCALE2 = 1.0D0/SQRT(SCALE2)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC_LUCI(VEC2,NVAR,-1,LU3)
*
        CALL REWINE( LU1,-1)
        CALL REWINE( LU3,-1)
        DO 2411 IVEC = 1,2
          CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC1,NVAR,-1,  LU1)
 2411   CONTINUE
*. Corresponding sigma vectors
        CALL REWINE ( LU3,-1)
        CALL REWINE( LU2,-1)
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2250 IVEC = 1, NVEC
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
          CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2250   CONTINUE
*
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC_LUCI(VEC1,NVAR,-1,  LU3)
* Sigma vector corresponding to second vector on LU1
        CALL REWINE(LU2,-1)
        CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC_LUCI(VEC2,NVAR,-1,LU3)
*
        CALL REWINE( LU2,-1)
        CALL REWINE( LU3,-1)
        DO 2400 IVEC = 1,2
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC2,NVAR,-1  ,LU2)
 2400   CONTINUE
        NVEC = 2
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 2410 IROOT = 1,NVEC
          AVEC((IROOT-1)*NVEC+IROOT) = 1.0D0
 2410   CONTINUE
*.Projected hamiltonian
       CALL REWINE( LU1 ,-1)
       DO 2010 IVEC = 1,NVEC
         CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL REWINE( LU2,-1)
         DO 2008 JVEC = 1, IVEC
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
 2008    CONTINUE
 2010  CONTINUE
      END IF
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU3,-1)
        DO 320 IROOT = 1, NROOT
          CALL REWINE( LU1,-1)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 200 IVEC = 1, NVEC
            CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  200     CONTINUE
*
          SCALE = INPROD(VEC1,VEC1,NVAR)
          SCALE  = 1.0D0/SQRT(SCALE)
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU3)
  320   CONTINUE
        CALL REWINE( LU1,-1)
        CALL REWINE( LU3,-1)
        DO 411 IVEC = 1,NROOT
          CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC1,NVAR,-1,  LU1)
  411   CONTINUE
* CORRESPONDING SIGMA VECTOR
        CALL REWINE ( LU3,-1)
        DO 329 IROOT = 1, NROOT
          CALL REWINE( LU2,-1)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 250 IVEC = 1, NVEC
            CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  250     CONTINUE
*
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC_LUCI(VEC1,NVAR,-1,  LU3)
  329   CONTINUE
* PLACE C IN LU1 AND HC IN LU2
        CALL REWINE( LU2,-1)
        CALL REWINE( LU3,-1)
        DO 400 IVEC = 1,NROOT
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC2,NVAR,-1  ,LU2)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
C
      END IF
C
C     IF( ITER .LT. MAXIT .AND. .NOT. CONVER) GOTO 1000
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*. Place first eigenvector in vec1
      CALL REWINE(LU1,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
* ( End of loop over iterations )
*
*
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT(/'  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(/' Convergence was obtained in ',I3,' iterations')
        END IF
*. Final eigenvalues
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
*
      IF ( IPRT .GT. 1 ) THEN
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(/'The final approximation to eigenvector')
            CALL REWINE( LU1,-1)
            CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
          END IF
          WRITE(6,1300)
 1300     FORMAT(/' Summary of iterations',
     +           /' ---------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (/' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,F14.7,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),
     &                 RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
 1234 CONTINUE
C
      RETURN
 1030 FORMAT(/2X,7F15.8,/,(2X,7F15.8))
 1120 FORMAT(/2X,I3,7F15.8,/,(5X,7F15.8))
      END
***********************************************************************

      SUBROUTINE MINGCG(MV8,LU1,LU2,LU3,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPROD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CONVER = .FALSE.
      ITER = 1
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X
      CALL REWINE(LU1,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINE(LU2,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
      CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-1.0D0,NVAR)
*
      RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
      ERROR(1) = RNORM
      CALL REWINE(LU2,-1)
      CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU2)
*. Preconditioner H times initial vector , H * R
*.H * R
      CALL REWINE(LUDIA,-1)
      CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
      CALL DIAVC2(VEC2,VEC1,VEC2,W,NVAR)
      IF(LUPROJ.NE.0) THEN
        CALL REWINE(LUPROJ,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL REWINE(LU2,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPROD(VEC1,VEC2,NVAR)
*. P = RHO * H*R
      RHO = 1.0D0
      CALL SCALVE(VEC2,RHO,NVAR)
      CALL REWINE(LU3,-1)
      CALL TODSC_LUCI(VEC2,NVAR,-1  ,LU3)
      CALL COPVEC(VEC2,VEC1,NVAR)
*.S = AP
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINE (LU3,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*.    P is assumed in VEC1 and S = A*P in VEC2

        NITER = NITER + 1
       IF ( NTEST .GE. 10 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPROD(VEC1,VEC2,NVAR)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S
        CALL REWINE(LU2,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-A,NVAR)
        CALL REWINE(LU2,-1)
        CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU2)
*.    new residual has been obtained , check for convergence
        RNORM = INPROD(VEC1,VEC1,NVAR)
        ERROR(ITER+1) = SQRT(RNORM)
*.    X = X + A * P
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        CALL REWINE(LU3,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC2,VEC1,1.0D0,A,NVAR)
        CALL REWINE(LU1,-1)
        CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU1)
*
        IF( SQRT(RNORM) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.       H * R
           CALL REWINE(LU2,-1)
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           CALL REWINE(LUDIA,-1)
           CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL DIAVC2(VEC1,VEC2,VEC1 ,W,NVAR)
           IF(LUPROJ.NE.0) THEN
             CALL REWINE(LUPROJ,-1)
             CALL FRMDSC_LUCI(VEC2,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
             CALL REWINE(LU2,-1)
             CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
           END IF
           GAMMA = INPROD(VEC1,VEC2,NVAR)
           B = GAMMA/C
*.       P = RHO*(H*R + B*P)
           CALL REWINE(LU3,-1)
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,B,NVAR)
*.       Define next RHO
           RHO = 1.0D0
           CALL SCALVE(VEC1,RHO,NVAR)
           CALL REWINE(LU3,-1)
           CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU3)
*.       S = MATRIX * P
           CALL MV8(VEC1,VEC2,0,0)
           CALL REWINE(LU3,-1)
           CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
*.End of prepations for next iteration
        END IF
*
*
 1000 CONTINUE
 1001 CONTINUE
      IF(NTEST .GT. 0 ) THEN
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(/'  convergence was obtained in...',I3,' iterations',
     +        /'  norm of residual..............',F13.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT +1 )
 1020  FORMAT(/' convergence was not obtained in',I3,'iterations',
     +        /' norm of residual...............',F13.8)
      END IF
      END IF
C
      IF(NTEST.GT. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(/' solution to set of linear equations')
       CALL REWINE(LU1,-1)
       CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
       CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call WRTMT_LU(vec2,1,nvar,1,nvar)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(/10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
C
      RETURN
      END
***********************************************************************

      SUBROUTINE MIXHONE(H1,H2,NREPTP,IREPTP,NOBTP,NSMOB)
*
* Replace selected type blocks of H1 with the corresponding blocks
* in H2
*
*. H1 and H2 are assumed to be in symmetry order !
*. -and total symmetric
*
*     Jeppe Olsen, March 14 1996 ( Still snowing in Lund )
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "orbinp.inc"
*. Specific input
      DIMENSION IREPTP(*)
      DIMENSION H2(*)
*. Input and output
      DIMENSION H1(*)
*
      DO ISMOB = 1, NSMOB
        IF (ISMOB.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISMOB-1)*(NTOOBS(ISMOB-1)+1)/2
        END IF
*. Loop over types for given symmetry
        DO ITP = 1, NOBTP
          IF(ITP.EQ.1) THEN
           IOBOFF = 1
          ELSE
            IOBOFF = IOBOFF + NOBPTS(ITP-1,ISMOB)
          END IF
          DO JTP = 1, ITP
            IF(JTP.EQ.1) THEN
             JOBOFF = 1
            ELSE
              JOBOFF = JOBOFF + NOBPTS(JTP-1,ISMOB)
            END IF
*. Number of elements in this type-type block
            LIOB = NOBPTS(ITP,ISMOB)
            LJOB = NOBPTS(JTP,ISMOB)
*
*. Should this block of H1 be replaced by corresponding block of H2
            IF(ITP.EQ.JTP) THEN
              IMOVE = 0
              DO KTP = 1, NREPTP
                IF(IREPTP(KTP).EQ.ITP) IMOVE = 1
              END DO
*
              IF(IMOVE.EQ.1) THEN
C?              WRITE(6,*) ' Block transfer ISMOB ITP JTP ',
C?   &          ISMOB,ITP,JTP
                DO IOB = IOBOFF,IOBOFF+LIOB-1
                  DO JOB = JOBOFF, IOB
                    H1(IOFF-1+IOB*(IOB-1)/2+JOB)
     &            = H2(IOFF-1+IOB*(IOB-1)/2+JOB)
                  END DO
                END DO
              END IF
*
            END IF
          END DO
        END DO
      END DO
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' MIXHONE in action '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' NSMOB NOBTP ', NSMOB,NOBTP
        WRITE(6,*) ' Types to be changed '
        CALL IWRTMA(IREPTP,1,NREPTP,1,NREPTP)
        WRITE(6,*) ' output H1 and H2 '
C       APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        CALL APRBLM2(H1,NTOOBS,NTOOBS,NSMOB,1)
        CALL APRBLM2(H2,NTOOBS,NTOOBS,NSMOB,1)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
* Transfer one-electron matrix H between memory and disc file in
* LUCIA format
*
* IFT = 1 => From disc ( read)
* IFT = 2 => To   disc (write)
*
* Note : File LUH is supposed to be at start of correct integral block
*
* Jeppe Olsen, Feb. 98
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER NRPSM(NSM),NCPSM(*)
      DIMENSION H(*)
*
*. Order of integrals are
*
*     Loop over Symmetry of row index => Symmetry of column  index
*      Loop over columns in symmetry block
*        Loop over rows in symmetry block
*        End of loop over rows in symmetry block
*      End of Loop over columns in symmetry block
*     End of loop over symmetry of row index
*
* Each symmetry block is thus given in complete form
* Note all integrals are in a single record
*
*. Length of list
C              NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,IPACK)
      LENGTH = NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,0)
*. and read/write
      WRITE(6,*) ' ONEEL, LUH = ', LUH
      IF(IFT.EQ.1) THEN
        DO IJ = 1, LENGTH
          READ(LUH,'(E22.15)') H(IJ)
        END DO
      ELSE IF (IFT.EQ.2) THEN
        DO IJ = 1, LENGTH
          WRITE(LUH,'(E22.15)') H(IJ)
        END DO
      END IF
*
      RETURN
      END
