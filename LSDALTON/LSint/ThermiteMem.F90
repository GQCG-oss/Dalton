!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE ThermiteMem_module
use precision
use memory_handling
integer(kind=long) :: allocIntmaxTUVdim
integer(kind=long) :: allocIntmaxWORK
integer(kind=long) :: allocIntmaxOrb
integer(kind=long) :: allocIntmaxPrim
!1=LHS,2=RHS,3=PASS,4=FTUV,5=PassF
integer(kind=long) :: nODint(5) 
integer(kind=long) :: nODrealk(5)
integer(kind=long) :: nODint_tp(5) 
integer(kind=long) :: nODrealk_tp(5)
integer(kind=long) :: INTWORKLENGTH
! PASSQ and PASSF is always private, FTUV always shared
! LHS always private 
! but RHS can be both private and shared 
real(realk),pointer :: INTWORK(:)
integer,pointer :: ODLHSintbuffer(:)
integer,pointer :: ODRHSintbuffer(:)
integer,pointer :: ODRHSintbuffer_tp(:)
integer,pointer :: ODPASSintbuffer(:)
integer,pointer :: ODFTUVintbuffer(:)
integer,pointer :: ODPASFintbuffer(:)
real(realk),pointer :: ODLHSrealkbuffer(:)
real(realk),pointer :: ODRHSrealkbuffer(:)
real(realk),pointer :: ODRHSrealkbuffer_tp(:)
real(realk),pointer :: ODPASSrealkbuffer(:)
real(realk),pointer :: ODFTUVrealkbuffer(:)
real(realk),pointer :: ODPASFrealkbuffer(:)
integer(kind=long) :: LODLHSintbuffer !Last nonused element (start being 1)
integer(kind=long) :: LODRHSintbuffer
integer(kind=long) :: LODLHSrealkbuffer
integer(kind=long) :: LODRHSrealkbuffer
integer(kind=long) :: LODRHSintbuffer_tp
integer(kind=long) :: LODRHSrealkbuffer_tp
integer(kind=long) :: LODPASSintbuffer
integer(kind=long) :: LODPASSrealkbuffer
integer(kind=long) :: LODFTUVintbuffer
integer(kind=long) :: LODFTUVrealkbuffer
integer(kind=long) :: LODPASFintbuffer
integer(kind=long) :: LODPASFrealkbuffer
integer(kind=long) :: LINTWORK
integer(kind=long) :: LINTWORKALLOC
integer(kind=long) :: LMAXWORKALLOC
!$OMP THREADPRIVATE(ODLHSintbuffer,ODLHSrealkbuffer,LODLHSintbuffer,&
!$OMP LODLHSrealkbuffer,LODRHSintbuffer_tp,LODRHSrealkbuffer_tp,&
!$OMP allocIntmaxTUVdim,allocIntmaxOrb,allocIntmaxPrim,allocIntmaxWORK,&
!$OMP ODPASSintbuffer,ODPASSrealkbuffer,ODPASFintbuffer,ODPASFrealkbuffer,&
!$OMP ODRHSintbuffer_tp,ODRHSrealkbuffer_tp,nODint_tp,nODrealk_tp,&
!$OMP LODPASSintbuffer,LODPASSrealkbuffer,LODPASFintbuffer,LODPASFrealkbuffer,&
!$OMP INTWORK,INTWORKLENGTH,LINTWORK,LINTWORKALLOC,LMAXWORKALLOC)

INTERFACE mem_ODpointer_alloc
  MODULE PROCEDURE real_setODpointer_1dim, int_setODpointer_1dim
END INTERFACE
CONTAINS
#ifdef VAR_LSDEBUGINT
subroutine set_allocIntmaxTUVdim(maxPrim,maxTUVdim,maxOrb)
implicit none
integer :: maxTUVdim,maxOrb,maxPrim
!
allocIntmaxPrim = maxPrim
allocIntmaxTUVdim = maxTUVdim
allocIntmaxOrb = maxOrb
end subroutine set_allocIntmaxTUVdim
#endif

subroutine init_workmem(maxWORK)
implicit none
integer(kind=long),intent(in) :: maxWORK
integer(kind=long) :: nsize
INTWORKLENGTH = maxWork+1
LINTWORK = 1
LMAXWORKALLOC = 0
LINTWORKALLOC = 0
call mem_alloc(INTWORK,maxwork)
nsize = SIZE(INTWORK)*mem_realsize
call mem_allocated_mem_intwork(nsize)
end subroutine init_workmem

subroutine free_workmem()
implicit none
integer(kind=long) :: nsize
nsize = SIZE(INTWORK)*mem_realsize
call mem_deallocated_mem_intwork(nsize)
INTWORKLENGTH = 0
LINTWORKALLOC = 0
LINTWORK = 1
CALL MEM_DEALLOC(INTWORK)
end subroutine free_workmem

SUBROUTINE INIT_BUFCOUNTERS(ielec)
implicit none
integer,intent(in) :: ielec
if(mem_InsideOMPsection)THEN
   nODint_tp(ielec) = 0
   nODrealk_tp(ielec) = 0
else
   nODint(ielec) = 0
   nODrealk(ielec) = 0
   IF(MOD(ielec,2).NE.0)then!1,3,5 = LHS,PASSQ,PASSF always private
      nODint_tp(ielec) = 0
      nODrealk_tp(ielec) = 0
   ENDIF
endif
END SUBROUTINE INIT_BUFCOUNTERS

SUBROUTINE SET_BUFCOUNTERS(ielec,nint,nrealk)
implicit none
integer,intent(in) :: ielec
integer(kind=long),intent(in) :: nint,nrealk 
if(mem_InsideOMPsection)THEN
   nODint_tp(ielec) = nint
   nODrealk_tp(ielec) = nrealk
else
   nODint(ielec) = nint
   nODrealk(ielec) = nrealk
   IF(MOD(ielec,2).NE.0)then!1,3,5 = LHS,PASSQ,PASSF always private
      nODint_tp(ielec) = nint
      nODrealk_tp(ielec) = nrealk
   ENDIF
endif
END SUBROUTINE SET_BUFCOUNTERS

SUBROUTINE ADD_BUFCOUNTERS(ielec,nint,nrealk)
implicit none
integer,intent(in) :: ielec
integer(kind=long),intent(in) :: nint,nrealk 
if(mem_InsideOMPsection)THEN
   nODint_tp(ielec) = nODint_tp(ielec) + nint
   nODrealk_tp(ielec) = nODrealk_tp(ielec) + nrealk
else
   nODint(ielec) = nODint(ielec) + nint
   nODrealk(ielec) = nODrealk(ielec) + nrealk
   IF(MOD(ielec,2).NE.0)then!1,3,5 = LHS,PASSQ,PASSF always private
      nODint_tp(ielec) = nODint_tp(ielec) + nint
      nODrealk_tp(ielec) = nODrealk_tp(ielec) + nrealk
   ENDIF
endif
END SUBROUTINE ADD_BUFCOUNTERS

SUBROUTINE MAX_BUFCOUNTERS(ielec,nint,nrealk)
implicit none
integer,intent(in) :: ielec
integer(kind=long),intent(in) :: nint,nrealk 
if(mem_InsideOMPsection)THEN
   nODint_tp(ielec) = MAX(nODint_tp(ielec),nint)
   nODrealk_tp(ielec) = MAX(nODrealk_tp(ielec),nrealk)
else
   nODint(ielec) = MAX(nODint(ielec),nint)
   nODrealk(ielec) = MAX(nODrealk(ielec),nrealk)
   IF(MOD(ielec,2).NE.0)then!1,3,5 = LHS,PASSQ,PASSF always private
      nODint_tp(ielec) = MAX(nODint_tp(ielec),nint)
      nODrealk_tp(ielec) = MAX(nODrealk_tp(ielec),nrealk)
   ENDIF
endif
END SUBROUTINE MAX_BUFCOUNTERS

SUBROUTINE ALLOC_ODLHS_BUFFERS
implicit none
integer(kind=long) :: nsize
nODint_tp(1) = nODint_tp(1)+1
nODrealk_tp(1) = nODrealk_tp(1)+1
CALL MEM_ALLOC(ODLHSintbuffer,nODint_tp(1))
CALL MEM_ALLOC(ODLHSrealkbuffer,nODrealk_tp(1))
LODLHSintbuffer = 1
LODLHSrealkbuffer = 1
nsize = SIZE(ODLHSintbuffer)*mem_intsize + SIZE(ODLHSrealkbuffer)*mem_realsize
call mem_allocated_mem_overlap(nsize)
END SUBROUTINE ALLOC_ODLHS_BUFFERS

SUBROUTINE DEALLOC_ODLHS_BUFFERS
implicit none
integer(kind=long) :: nsize
nsize = SIZE(ODLHSintbuffer)*mem_intsize + SIZE(ODLHSrealkbuffer)*mem_realsize
call mem_deallocated_mem_overlap(nsize)
CALL MEM_DEALLOC(ODLHSintbuffer)
CALL MEM_DEALLOC(ODLHSrealkbuffer)
nODint_tp(1)=0
nODrealk_tp(1)=0
END SUBROUTINE DEALLOC_ODLHS_BUFFERS

SUBROUTINE ALLOC_ODRHS_BUFFERS
implicit none
integer(kind=long) :: nsize
if(mem_InsideOMPsection)THEN
   nODint_tp(2) = nODint_tp(2)+1
   nODrealk_tp(2) = nODrealk_tp(2)+1
   CALL MEM_ALLOC(ODRHSintbuffer_tp,nODint_tp(2))
   CALL MEM_ALLOC(ODRHSrealkbuffer_tp,nODrealk_tp(2))
   LODRHSintbuffer_tp = 1
   LODRHSrealkbuffer_tp = 1
   nsize = SIZE(ODRHSintbuffer_tp)*mem_intsize + SIZE(ODRHSrealkbuffer_tp)*mem_realsize
   call mem_allocated_mem_overlap(nsize)
else
   nODint(2) = nODint(2)+1
   nODrealk(2) = nODrealk(2)+1
   CALL MEM_ALLOC(ODRHSintbuffer,nODint(2))
   CALL MEM_ALLOC(ODRHSrealkbuffer,nODrealk(2))
   LODRHSintbuffer = 1
   LODRHSrealkbuffer = 1
   nsize = SIZE(ODRHSintbuffer)*mem_intsize + SIZE(ODRHSrealkbuffer)*mem_realsize
   call mem_allocated_mem_overlap(nsize)
endif
END SUBROUTINE ALLOC_ODRHS_BUFFERS

SUBROUTINE DEALLOC_ODRHS_BUFFERS
implicit none
integer(kind=long) :: nsize
if(mem_InsideOMPsection)THEN
   nsize = SIZE(ODRHSintbuffer_tp)*mem_intsize + SIZE(ODRHSrealkbuffer_tp)*mem_realsize
   call mem_deallocated_mem_overlap(nsize)
   CALL MEM_DEALLOC(ODRHSintbuffer_tp)
   CALL MEM_DEALLOC(ODRHSrealkbuffer_tp)
   nODint_tp(2)=0
   nODrealk_tp(2)=0
else
   nsize = SIZE(ODRHSintbuffer)*mem_intsize + SIZE(ODRHSrealkbuffer)*mem_realsize
   call mem_deallocated_mem_overlap(nsize)
   CALL MEM_DEALLOC(ODRHSintbuffer)
   CALL MEM_DEALLOC(ODRHSrealkbuffer)
   nODint(2)=0
   nODrealk(2)=0
endif
END SUBROUTINE DEALLOC_ODRHS_BUFFERS

SUBROUTINE REINIT_OD_LHS
implicit none
LODLHSintbuffer = 1
LODLHSrealkbuffer = 1
END SUBROUTINE REINIT_OD_LHS

SUBROUTINE REINIT_OD_PASSF
implicit none
LODPASFintbuffer = 1
LODPASFrealkbuffer = 1
END SUBROUTINE REINIT_OD_PASSF

SUBROUTINE ALLOC_ODPASS_BUFFERS
implicit none
integer(kind=long) :: nsize
nODint_tp(3) = nODint_tp(3)+1
nODrealk_tp(3) = nODrealk_tp(3)+1
CALL MEM_ALLOC(ODPASSintbuffer,nODint_tp(3))
CALL MEM_ALLOC(ODPASSrealkbuffer,nODrealk_tp(3))
LODPASSintbuffer = 1
LODPASSrealkbuffer = 1
nsize = SIZE(ODPASSintbuffer)*mem_intsize + SIZE(ODPASSrealkbuffer)*mem_realsize
call mem_allocated_mem_overlap(nsize)
END SUBROUTINE ALLOC_ODPASS_BUFFERS

SUBROUTINE DEALLOC_ODPASS_BUFFERS
implicit none
integer(kind=long) :: nsize
nsize = SIZE(ODPASSintbuffer)*mem_intsize + SIZE(ODPASSrealkbuffer)*mem_realsize
call mem_deallocated_mem_overlap(nsize)
CALL MEM_DEALLOC(ODPASSintbuffer)
CALL MEM_DEALLOC(ODPASSrealkbuffer)
nODint_tp(3)=0
nODrealk_tp(3)=0
END SUBROUTINE DEALLOC_ODPASS_BUFFERS

SUBROUTINE ALLOC_ODFTUV_BUFFERS
implicit none
integer(kind=long) :: nsize
nODint(4) = nODint(4)+1
nODrealk(4) = nODrealk(4)+1
CALL MEM_ALLOC(ODFTUVintbuffer,nODint(4))
CALL MEM_ALLOC(ODFTUVrealkbuffer,nODrealk(4))
LODFTUVintbuffer = 1
LODFTUVrealkbuffer = 1
nsize = SIZE(ODFTUVintbuffer)*mem_intsize + SIZE(ODFTUVrealkbuffer)*mem_realsize
call mem_allocated_mem_overlap(nsize)
END SUBROUTINE ALLOC_ODFTUV_BUFFERS

SUBROUTINE DEALLOC_ODFTUV_BUFFERS
implicit none
integer(kind=long) :: nsize
nsize = SIZE(ODFTUVintbuffer)*mem_intsize + SIZE(ODFTUVrealkbuffer)*mem_realsize
call mem_deallocated_mem_overlap(nsize)
CALL MEM_DEALLOC(ODFTUVintbuffer)
CALL MEM_DEALLOC(ODFTUVrealkbuffer)
nODint(4)=0
nODrealk(4)=0
END SUBROUTINE DEALLOC_ODFTUV_BUFFERS

SUBROUTINE ALLOC_ODPASSF_BUFFERS
implicit none
integer(kind=long) :: nsize
nODint_tp(5) = nODint_tp(5)+1
nODrealk_tp(5) = nODrealk_tp(5)+1
CALL MEM_ALLOC(ODPASFintbuffer,nODint_tp(5))
CALL MEM_ALLOC(ODPASFrealkbuffer,nODrealk_tp(5))
LODPASFintbuffer = 1
LODPASFrealkbuffer = 1
nsize = SIZE(ODPASFintbuffer)*mem_intsize + SIZE(ODPASFrealkbuffer)*mem_realsize
call mem_allocated_mem_overlap(nsize)
END SUBROUTINE ALLOC_ODPASSF_BUFFERS

SUBROUTINE DEALLOC_ODPASSF_BUFFERS
implicit none
integer(kind=long) :: nsize
nsize = SIZE(ODPASFintbuffer)*mem_intsize + SIZE(ODPASFrealkbuffer)*mem_realsize
call mem_deallocated_mem_overlap(nsize)
CALL MEM_DEALLOC(ODPASFintbuffer)
CALL MEM_DEALLOC(ODPASFrealkbuffer)
nODint_tp(5)=0
nODrealk_tp(5)=0
END SUBROUTINE DEALLOC_ODPASSF_BUFFERS

subroutine real_setODpointer_1dim(A,n,ielec)
implicit none
integer,intent(in)  :: n,ielec
REAL(REALK),pointer :: A(:)
integer(kind=long) :: nsize
nsize = n
nullify(A)
IF(ielec.EQ.2)THEN   
   if(mem_InsideOMPsection)THEN
      IF(LODRHSrealkbuffer_tp+nsize.GT.nODrealk_tp(2))THEN
         print*,'size of RHS real buffer:',nODrealk_tp(2)
         print*,'current first free index:',LODRHSrealkbuffer_tp
         print*,'requested size:',nsize
         call lsquit('OD real buffer full in tp RHS 1dim',-1)
      ENDIF
      A => ODRHSrealkbuffer_tp(LODRHSrealkbuffer_tp:LODRHSrealkbuffer_tp+n-1)
      LODRHSrealkbuffer_tp=LODRHSrealkbuffer_tp+nsize 
   else
      IF(LODRHSrealkbuffer+nsize.GT.nODrealk(2))THEN
         print*,'size of RHS real buffer:',nODrealk(2)
         print*,'current first free index:',LODRHSrealkbuffer
         print*,'requested size:',nsize
         call lsquit('OD real buffer full in RHS 1dim',-1)
      ENDIF
      A => ODRHSrealkbuffer(LODRHSrealkbuffer:LODRHSrealkbuffer+n-1)
      LODRHSrealkbuffer=LODRHSrealkbuffer+nsize 
   endif
ELSEIF(ielec.EQ.1)THEN
   IF(LODLHSrealkbuffer+nsize.GT.nODrealk_tp(1))THEN
      print*,'size of LHS real buffer:',nODrealk_tp(1)
      print*,'current first free index:',LODLHSrealkbuffer
      print*,'requested size:',nsize
      call lsquit('OD real buffer full in LHS 1dim',-1)
   ENDIF
   A => ODLHSrealkbuffer(LODLHSrealkbuffer:LODLHSrealkbuffer+n-1)
   LODLHSrealkbuffer=LODLHSrealkbuffer+nsize
ELSEIF(ielec.EQ.3)THEN
   IF(LODPASSrealkbuffer+nsize.GT.nODrealk_tp(3))THEN
      print*,'size of PASS real buffer:',nODrealk_tp(3)
      print*,'current first free index:',LODPASSrealkbuffer
      print*,'requested size:',nsize
      call lsquit('OD real buffer full in PASS 1dim',-1)
   ENDIF
   A => ODPASSrealkbuffer(LODPASSrealkbuffer:LODPASSrealkbuffer+n-1)
   LODPASSrealkbuffer=LODPASSrealkbuffer+nsize
ELSEIF(ielec.EQ.4)THEN
   IF(LODFTUVrealkbuffer+nsize.GT.nODrealk(4))THEN
      print*,'size of FTUV real buffer:',nODrealk(4)
      print*,'current first free index:',LODFTUVrealkbuffer
      print*,'requested size:',nsize
      call lsquit('OD real buffer full in FTUV 1dim',-1)
   ENDIF
   A => ODFTUVrealkbuffer(LODFTUVrealkbuffer:LODFTUVrealkbuffer+n-1)
   LODFTUVrealkbuffer=LODFTUVrealkbuffer+nsize
ELSEIF(ielec.EQ.5)THEN
   IF(LODPASFrealkbuffer+nsize.GT.nODrealk_tp(5))THEN
      print*,'size of PASS real buffer:',nODrealk_tp(5)
      print*,'current first free index:',LODPASFrealkbuffer
      print*,'requested size:',nsize
      call lsquit('OD real buffer full in PASSF 1dim',-1)
   ENDIF
   A => ODPASFrealkbuffer(LODPASFrealkbuffer:LODPASFrealkbuffer+n-1)
   LODPASFrealkbuffer=LODPASFrealkbuffer+nsize
ELSE
   call lsquit('not recognized option in real_setODpointer_1dim',-1)
ENDIF
end subroutine real_setODpointer_1dim

subroutine int_setODpointer_1dim(A,n,ielec)
implicit none
integer,intent(in)  :: n,ielec
integer,pointer :: A(:)
integer(kind=long) :: nsize
integer :: s,e
nsize = n
nullify(A)
IF(ielec.EQ.2)THEN   
   if(mem_InsideOMPsection)THEN
      IF(LODRHSintbuffer_tp+nsize.GT.nODint_tp(2))THEN
         call lsquit('OD int buffer full in tp RHS 1dim',-1)
      ENDIF
      A => ODRHSintbuffer_tp(LODRHSintbuffer_tp:LODRHSintbuffer_tp+nsize-1)
      LODRHSintbuffer_tp=LODRHSintbuffer_tp+nsize
   else
      IF(LODRHSintbuffer+nsize.GT.nODint(2))THEN
         call lsquit('OD int buffer full in RHS 1dim',-1)
      ENDIF
      A => ODRHSintbuffer(LODRHSintbuffer:LODRHSintbuffer+nsize-1)
      LODRHSintbuffer=LODRHSintbuffer+nsize
   endif
ELSEIF(ielec.EQ.1)THEN
   IF(LODLHSintbuffer+nsize.GT.nODint_tp(1))THEN
      print*,'size(ODLHSintbuffer)',size(ODLHSintbuffer)
      print*,'LODLHSintbuffer',LODLHSintbuffer
      print*,'nsize',nsize
      print*,'nODint_tp(1)',nODint_tp(1)
      call lsquit('OD int buffer full in LHS 1dim',-1)
   ENDIF
   A => ODLHSintbuffer(LODLHSintbuffer:LODLHSintbuffer+nsize-1)
   LODLHSintbuffer=LODLHSintbuffer+nsize
ELSEIF(ielec.EQ.3)THEN
   IF(LODPASSintbuffer+nsize.GT.nODint_tp(3))THEN
      print*,'LODPASSintbuffer',LODPASSintbuffer
      print*,'nsize',nsize
      print*,'nODint_tp(3)',nODint_tp(3)
      print*,'size(ODPASSintbuffer)',size(ODPASSintbuffer)
      call lsquit('OD int buffer full in PASS 1dim',-1)
   ENDIF
   A => ODPASSintbuffer(LODPASSintbuffer:LODPASSintbuffer+nsize-1)
   LODPASSintbuffer=LODPASSintbuffer+nsize
ELSEIF(ielec.EQ.4)THEN
   IF(LODFTUVintbuffer+nsize.GT.nODint(4))THEN
      print*,'size(ODFTUVintbuffer)',size(ODFTUVintbuffer)
      print*,'LODFTUVintbuffer',LODFTUVintbuffer
      print*,'nsize',nsize
      call lsquit('OD int buffer full in FTUV 1dim',-1)
   ENDIF
   A => ODFTUVintbuffer(LODFTUVintbuffer:LODFTUVintbuffer+nsize-1)
   LODFTUVintbuffer=LODFTUVintbuffer+nsize
ELSEIF(ielec.EQ.5)THEN
   IF(LODPASFintbuffer+nsize.GT.nODint_tp(5))THEN
      print*,'size() of PASSF int buffer:',size(ODPASFintbuffer)
      print*,'size of PASSF int buffer:',nODint_tp(5)
      print*,'current first free index:',LODPASFintbuffer
      print*,'requested size:',nsize
      call lsquit('OD int buffer full in PASSF 1dim',-1)
   ENDIF
   A => ODPASFintbuffer(LODPASFintbuffer:LODPASFintbuffer+nsize-1)
   LODPASFintbuffer=LODPASFintbuffer+nsize
ELSE
   call lsquit('not recognized option in int_setODpointer_1dim',-1)
ENDIF
end subroutine int_setODpointer_1dim

subroutine mem_workpointer_alloc(A,n)
implicit none
integer,intent(in)  :: n
REAL(REALK),pointer :: A(:)
nullify(A)
IF(LINTWORK+n.GT.INTWORKLENGTH)THEN
   print*,'n=',n
   print*,'LINTWORK=',LINTWORK   
   print*,'INTWORKLENGTH=',INTWORKLENGTH   
   print*,'LINTWORKALLOC=',LINTWORKALLOC
   call lsquit('Work buffer full in Thermite',-1)
ENDIF
A => INTWORK(LINTWORK:LINTWORK+n-1)
LINTWORK=LINTWORK+n
LINTWORKALLOC=LINTWORKALLOC+n
LMAXWORKALLOC = MAX(LMAXWORKALLOC,LINTWORK)
end subroutine mem_workpointer_alloc

subroutine mem_workpointer_dealloc(A)
implicit none
REAL(REALK),pointer :: A(:)
integer(kind=long) :: n
n = size(A)
LINTWORKALLOC=LINTWORKALLOC-n
IF(LINTWORKALLOC.EQ.0)LINTWORK=1
nullify(A)
end subroutine mem_workpointer_dealloc

END MODULE ThermiteMem_module
