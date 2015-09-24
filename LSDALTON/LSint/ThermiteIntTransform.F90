!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE ThermiteIntTransform_module
use precision
use memory_handling
use ls_util
!!
use infpar_module
private
logical,save :: DoThermiteIntTransform
integer,save :: nT1,nTMO1,nT2,nTMO2,nTA,nTMaxN,TITThreadID,nTITthreads
real(realk),save,pointer :: TCMO1(:,:)
real(realk),save,pointer :: TCMO2(:,:)
real(realk),save,pointer :: TmpAlphaCD(:)
integer,save,pointer :: TOrbToFull(:,:),iLocalTIT(:),iLocalTIT2(:)
integer,save,pointer :: TITGindexToLocal(:)
!$OMP THREADPRIVATE(TITThreadID)
public :: InitThermiteIntTransform, FreeThermiteIntTransform,&
     & ThermiteIntTransform_AOtoMOinternal,DoThermiteIntTransform,&
     & iLocalTIT,TOrbToFull,DF3centerTrans3,TITThreadID,&
     & ThermiteIntTransform_alloc_TmpArray,iLocalTIT2,&
     & Initthermiteintthreadid,nTITthreads,&
     & ThermiteIntTransform_AOtoMOinternalFinal,&
     & TITGindexToLocal,InitThermiteIntTransform1,&
     & nullThermiteIntTransform

CONTAINS
subroutine nullThermiteIntTransform() 
implicit none
nullify(TITGindexToLocal)
nullify(TCMO1)
nullify(TCMO2)
nullify(TOrbToFull)
nullify(TmpAlphaCD)
nullify(iLocalTIT)  
nullify(iLocalTIT2)  
end subroutine nullThermiteIntTransform

subroutine InitThermiteIntTransform1(GindexToLocal,nAux) 
implicit none
integer :: nAux
integer :: GindexToLocal(nAux)
!
integer :: I
call mem_alloc(TITGindexToLocal,nAux)
do I=1,nAux
   TITGindexToLocal(I) = GindexToLocal(I)
enddo
end subroutine InitThermiteIntTransform1

subroutine InitThermiteIntTransform(n1,nMO1,n2,nMO2,C1,C2,nA,MaxN,nthreads)
  implicit none
  integer,intent(in) :: n1,nMO1,n2,nMO2,nA,MaxN,nthreads
  real(realk),target :: C1(n1,nMO1),C2(n2,nMO2)
  
  nTA = nA
  nTMaxN = MaxN
  nT1 = n1
  nTMO1 = nMO1
  nT2 = n2
  nTMO2 = nMO2
  TCMO1 => C1
  TCMO2 => C2
!  call mem_alloc(TCMO1,nT1,nTMO1)
!  call mem_alloc(TCMO2,nT2,nTMO2)
!  call dcopy(nT1*nTMO1,C1,1,TCMO1,1)
!  call dcopy(nT2*nTMO2,C2,1,TCMO2,1)
  call mem_alloc(TOrbToFull,MaxN,nthreads)
  nullify(TmpAlphaCD)
  DoThermiteIntTransform = .TRUE.
  call mem_alloc(iLocalTIT,nthreads)  
  iLocalTIT = 0
  call mem_alloc(iLocalTIT2,nthreads)  
  iLocalTIT2 = 0
  nTITthreads = nthreads 
end subroutine InitThermiteIntTransform

subroutine InitThermiteIntThreadID(threadID)
  implicit none
  integer,intent(in) :: ThreadID
  TITThreadID = ThreadID+1
end subroutine InitThermiteIntThreadID

subroutine ThermiteIntTransform_alloc_TmpArray(use_bg_buf)
  implicit none
  logical :: use_bg_buf
  integer(kind=long) :: n8
  n8 = nTMaxN*nT1*nTMO2
  IF(use_bg_buf)THEN
     call mem_pseudo_alloc(TmpAlphaCD,n8)
  ELSE
     call mem_alloc(TmpAlphaCD,n8)
  ENDIF
end subroutine ThermiteIntTransform_alloc_TmpArray

!TmpArray is batch of  AO integrals. FullAlphaCD is used to allocate MO integrals
subroutine ThermiteIntTransform_AOtoMOInternal(FullAlphaCD,nA,n1,n2,TmpArray,m1,m2,m3)
implicit none
integer,intent(in) :: nA,n1,n2,m1,m2,m3
real(realk),intent(inout) :: TmpArray(m1,m2,m3,nTITthreads)
real(realk),intent(inout) :: FullAlphaCD(nA,n1,n2)
!local variables
integer :: M,N,K
#ifdef VAR_LSDEBUG
IF(n1.NE.nTMO1)call lsquit('dimmismatch1 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(n2.NE.nTMO2)call lsquit('dimmismatch2 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(m1.NE.nTMaxN)call lsquit('dimmismatch3 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(m2.NE.nT1)call lsquit('dimmismatch4 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(m3.NE.nT2)call lsquit('dimmismatch5 in ThermiteIntTransform_AOtoMOInternal',-1)
#endif
!have to be omp critical because only 1 TmpAlphaCD array - and this would be overwritten
!$OMP CRITICAL (AOtoMOInternal)
M = nTMaxN*nT1   !rows of Output Matrix
N = nTMO2        !columns of Output Matrix
K = m3           !summation dimension
#ifdef VAR_OMP
IF(nTITthreads.GT.1)THEN
   IF(TITThreadID.EQ.1)THEN
      call dgemm_TS('N','N',M,N,K,1.0E0_realk,TmpArray,M,TCMO2,&
           & K,0.0E0_realk,TmpAlphaCD,M)
   ELSE
      call dgemm_TS('N','N',M,N,K,1.0E0_realk,TmpArray(1,1,1,TITThreadID),M,TCMO2,&
           & K,0.0E0_realk,TmpAlphaCD,M)
   ENDIF
ELSE
   call dgemm('N','N',M,N,K,1.0E0_realk,TmpArray,M,TCMO2,&
        & K,0.0E0_realk,TmpAlphaCD,M)
ENDIF
#else
call dgemm('N','N',M,N,K,1.0E0_realk,TmpArray(1,1,1,TITThreadID),M,TCMO2,&
     & K,0.0E0_realk,TmpAlphaCD,M)
#endif
!TmpAlphaCD(nTMaxN*nT1,nTMO2)
call DF3centerTrans2(nTMO1,nTMO2,nT1,nTMaxN,nTA,TCMO1,TmpAlphaCD,FullAlphaCD)
!$OMP END CRITICAL (AOtoMOInternal)
iLocalTIT(TITThreadID) = 0
iLocalTIT2(TITThreadID) = 0
! NO OpenMP parallelization because this is called from a OpenMP parallized part
DO M=1,m3
   DO N=1,m2
      DO K=1,m1
         TmpArray(K,N,M,TITThreadID) = 0.0E0_realk
      ENDDO
   ENDDO
ENDDO
end subroutine ThermiteIntTransform_AOtoMOInternal

subroutine DF3centerTrans2(nvirt,nocc,nbast,nAuxA,nAFull,Cvirt,AlphaCD2,FullAlphaCD)
  implicit none
  integer,intent(in) :: nvirt,nocc,nbast,nAuxA,nAFull
  real(realk),intent(in) :: Cvirt(nbast,nvirt),AlphaCD2(nAuxA,nBast,nocc)
  real(realk),intent(inout) :: FullAlphaCD(nAFull,nvirt,nocc)
  !
  integer ::ADIAG,IDIAG,ALPHAAUX,ALPHA,A
  real(realk) :: TMP
! NO OpenMP parallelization because this is called from a OpenMP parallized part
  do IDIAG = 1,nocc
     do ADIAG = 1,nvirt
        do ALPHAAUX = 1,iLocalTIT2(TITThreadID)
           A = TOrbToFull(ALPHAAUX,TITThreadID)
           TMP = 0.0E0_realk
           do ALPHA = 1,nbast
              TMP = TMP + Cvirt(ALPHA,ADIAG)*AlphaCD2(ALPHAAUX,ALPHA,IDIAG)
           enddo
           FullAlphaCD(A,ADIAG,IDIAG) = TMP
        enddo
     enddo
  enddo
end subroutine DF3centerTrans2

!TmpArray is batch of  AO integrals. FullAlphaCD is used to allocate MO integrals
subroutine ThermiteIntTransform_AOtoMOInternalFinal(FullAlphaCD,nA,n1,n2,TmpArray,m1,m2,m3)
implicit none
integer,intent(in) :: nA,n1,n2,m1,m2,m3
real(realk),intent(inout) :: TmpArray(m1,m2,m3,nTITthreads)
real(realk),intent(inout) :: FullAlphaCD(nA,n1,n2)
!local variables
integer :: M,N,K
#ifdef VAR_LSDEBUG
IF(nA.NE.nTA)call lsquit('dimmismatch0 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(n1.NE.nTMO1)call lsquit('dimmismatch1 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(n2.NE.nTMO2)call lsquit('dimmismatch2 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(m1.NE.nTMaxN)call lsquit('dimmismatch3 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(m2.NE.nT1)call lsquit('dimmismatch4 in ThermiteIntTransform_AOtoMOInternal',-1)
IF(m3.NE.nT2)call lsquit('dimmismatch5 in ThermiteIntTransform_AOtoMOInternal',-1)
#endif
M = nTMaxN*nT1   !rows of Output Matrix
N = nTMO2        !columns of Output Matrix
K = m3           !summation dimension
call dgemm('N','N',M,N,K,1.0E0_realk,TmpArray(:,:,:,TITThreadID),M,TCMO2,&
     & K,0.0E0_realk,TmpAlphaCD,M)
call DF3centerTrans2OMP(nTMO1,nTMO2,nT1,nTMaxN,nTA,TCMO1,TmpAlphaCD,FullAlphaCD)
end subroutine ThermiteIntTransform_AOtoMOInternalFinal

subroutine DF3centerTrans2OMP(nvirt,nocc,nbast,nAuxA,nAFull,Cvirt,AlphaCD2,FullAlphaCD)
  implicit none
  integer,intent(in) :: nvirt,nocc,nbast,nAuxA,nAFull
  real(realk),intent(in) :: Cvirt(nbast,nvirt),AlphaCD2(nAuxA,nBast,nocc)
  real(realk),intent(inout) :: FullAlphaCD(nAFull,nvirt,nocc)
  !
  integer ::ADIAG,IDIAG,ALPHAAUX,ALPHA,A,TID,nA
  real(realk) :: TMP
  TID = TITThreadID 
  nA = iLocalTIT2(TID)
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
  !$OMP PRIVATE(ADIAG,IDIAG,ALPHAAUX,ALPHA,TMP,A) &
  !$OMP SHARED(nvirt,nocc,nbast,Cvirt,AlphaCD2,&
  !$OMP        FullAlphaCD,iLocalTIT2,TOrbToFull,TID,nA) 
  do IDIAG = 1,nocc
     do ADIAG = 1,nvirt
        do ALPHAAUX = 1,iLocalTIT2(TID)
           A = TOrbToFull(ALPHAAUX,TID)
           TMP = 0.0E0_realk
           do ALPHA = 1,nbast
              TMP = TMP + Cvirt(ALPHA,ADIAG)*AlphaCD2(ALPHAAUX,ALPHA,IDIAG)
           enddo
           FullAlphaCD(A,ADIAG,IDIAG) = TMP
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine DF3centerTrans2OMP

subroutine FreeThermiteIntTransform(use_bg_buf)
  implicit none
  logical :: use_bg_buf
  IF(associated(TmpAlphaCD))THEN
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(TmpAlphaCD)
     ELSE
        call mem_dealloc(TmpAlphaCD)
     ENDIF
!     call mem_dealloc(TCMO1)
!     call mem_dealloc(TCMO2)
     call mem_dealloc(TOrbToFull)
     call mem_dealloc(iLocalTIT)  
     call mem_dealloc(iLocalTIT2)  
  ENDIF
  IF(associated(TITGindexToLocal))THEN
     call mem_dealloc(TITGindexToLocal)
  ENDIF
  call nullThermiteIntTransform() 
end subroutine FreeThermiteIntTransform

END MODULE ThermiteIntTransform_module
