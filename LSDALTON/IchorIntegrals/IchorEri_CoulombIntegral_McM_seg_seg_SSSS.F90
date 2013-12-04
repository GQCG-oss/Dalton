module IchorCoulombIntegral_seg_seg_SSSS_mod
use IchorprecisionModule
use IchorCommonModule
private 
public :: IchorCoulombIntegral_seg_seg_SSSS
CONTAINS
  subroutine IchorCoulombIntegral_seg_seg_SSSS(nPrimP,nPrimQ,nPasses,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,TABFJW,&
       & reducedExponents,integralPrefactor,PQorder,CDAB)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses
!    integer,intent(in) :: nTABFJW1,nTABFJW2
    real(realk),intent(in) :: pcent(3,nPrimP,nPasses),qcent(3,nPrimQ)
    real(realk),target,intent(in) :: QpreExpFac(nPrimQ)
    real(realk),target,intent(in) :: PpreExpFac(nPrimP,nPasses)
!    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    real(realk),intent(in) :: TABFJW(0:3,0:1200)
    real(realk),intent(inout) :: CDAB(nPasses)
    real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
    logical :: PQorder
    !
    real(realk) :: RJ000,squaredDistance
    real(realk) :: px,py,pz,qx,qy,qz,pqx,pqy,pqz,tmpP
    integer :: nPrimPQ,offset,IPNT,ipq,iPrimQ,iPrimP,IPASSP
    REAL(REALK), PARAMETER :: D2 = 2.0E0_realk,D4 = 4.0E0_realk,D1 = 1.0E0_realk
    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
    REAL(REALK), PARAMETER :: SQRPIH = 1.77245385090551602730E00_realk/D2
    Real(realk), parameter :: PI=3.14159265358979323846E0_realk
    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI, D12 = 12.0E0_realk
    REAL(REALK), PARAMETER :: D100 = 100.0E0_realk,TENTH = 1.0E-2_realk
    REAL(REALK), PARAMETER :: COEF2 = 0.5E0_realk, COEF3=-D1/6.0E0_realk
    REAL(REALK), PARAMETER :: D36 = 36E0_realk,HALF=0.5E0_realk
    real(realk) :: WVAL,R,RWVAL,REXPW,GVAL,WDIFF,W2,W3,tmp
!ifdef VAR_DEBUG
    IF(.NOT.PQorder)THEN
       call Ichorquit('IchorCoulombIntegral_seg_seg_SSSS expect to get PQ ordering',-1)
    ENDIF
    !nested OMP parallel or GPU/MIC
    !reducedExponents, integralPrefactor as input
    !I would like to avoid the IF statement inside the loop (not good from a gpu point of view)
    !Add nPassP as well so that nPasses is big!! far greater than nThreads/nGPUcores
    !precalc squaredDistance and sort -> alot of work.
    nPrimPQ = nPrimP*nPrimQ

!All cores have 2*nPrimP*nPrimQ : RJ000, squaredDistance
  
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IPASSP,iPrimQ,&
!$OMP qx,qy,qz,offset,iPrimP,pqx,pqy,pqz,squaredDistance,&
!$OMP ipq,WVAL,IPNT,WDIFF,W2,W3,R,RJ000,REXPW,RWVAL,GVAL,&
!$OMP tmp,tmpQ) Shared(Qcent,Pcent,reducedExponents,TABFJW,&
!$OMP integralPrefactor,QpreExpFac,PpreExpFac,CDAB,nPasses,&
!$OMP nPrimQ,nPrimP,nPrimPQ)
    DO IPASSP = 1,nPasses
       !build squaredDistance(nPrimP,nPrimQ,nPasses)
       tmp = 0.0E0_realk
       DO iPrimP=1, nPrimP
          px = -Pcent(1,iPrimP,iPassP)
          py = -Pcent(2,iPrimP,iPassP)
          pz = -Pcent(3,iPrimP,iPassP)
          TMPP = PpreExpFac(iPrimP,iPassP)
          DO iPrimQ=1, nPrimQ
             pqx = px + Qcent(1,iPrimQ)
             pqy = py + Qcent(1,iPrimQ)
             pqz = pz + Qcent(1,iPrimQ)
             !(nPrimP,nPrimQ)
             squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz
             !(nPrimP,nPrimQ)*(nPrimP,nPrimQ)
             WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
             IF (WVAL .LT. D12) THEN
                IPNT = NINT(D100*WVAL)
                WDIFF = WVAL - TENTH*IPNT
                W2    = WDIFF*WDIFF
                W3    = W2*WDIFF
                W2    = W2*COEF2
                W3    = W3*COEF3
                RJ000 = TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3
             !  12 < WVAL <= (2J+36) 
             ELSE IF (WVAL.LE.D36) THEN
                REXPW = HALF*EXP(-WVAL)
                RWVAL = D1/WVAL
                GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
                RJ000 = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
                !  (2J+36) < WVAL 
             ELSE
                RWVAL = PID4/WVAL
                RJ000 = SQRT(RWVAL)
             END IF
             !scale RJ000
             tmp = tmp + QpreExpFac(iPrimQ)*integralPrefactor(IPrimQ,iPrimP)*RJ000*TMPP
          ENDDO
       ENDDO
       CDAB(iPassP) = tmp
    ENDDO
!$OMP END PARALLEL DO
  end subroutine IchorCoulombIntegral_seg_seg_SSSS

END MODULE IchorCoulombIntegral_seg_seg_SSSS_mod
