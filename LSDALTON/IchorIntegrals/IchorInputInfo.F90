MODULE IchorInputInfoMod
  use IchorprecisionMod
  use IchorCommonMod
!  use IchorMemory
  public :: InputM,InputMspec,InitIchorInputInfoModule,&
       & FreeIchorInputInfoModule,IchorInputInfoM1,&
       & IchorInputInfoM2,IchorInputInfoM3,IchorInputInfoM4,&
       & IchorInputInfoSpec
  private
  
  TYPE InputMatrixItem
     real(realk),pointer :: M(:,:) 
     integer :: n1,n2
  END type InputMatrixItem
  integer,save :: nInputMatrixItem
  type(InputMatrixItem),pointer :: InputM(:)
  integer,save :: InputMspec(4)
  integer(KIND=long),parameter :: memrealsz=8
  integer(KIND=long),save :: mem_allocated_inputichor
CONTAINS
subroutine InitIchorInputInfoModule()
implicit none
nInputMatrixItem = 0
mem_allocated_inputichor = 0 
end subroutine InitIchorInputInfoModule

subroutine FreeIchorInputInfoModule()
implicit none
integer :: I
do I=1,nInputMatrixItem
!   call mem_ichor_dealloc(InputM(I)%M)
   mem_allocated_inputichor = mem_allocated_inputichor - &
        & size(InputM(I)%M,KIND=long)*memrealsz
   deallocate(InputM(I)%M)
   nullify(InputM(I)%M)
enddo
deallocate(InputM)
nullify(InputM)
nInputMatrixItem = 0
IF(mem_allocated_inputichor.NE.0)THEN
   CALL ICHORQUIT('FreeIchorInputInfoModule: still mem allocated ',-1)
ENDIF
end subroutine FreeIchorInputInfoModule

subroutine IchorInputInfoM1(M,n1,n2)
implicit none
integer,intent(in) :: n1,n2
real(realk),intent(in) :: M(n1,n2)
nInputMatrixItem = 1
allocate(InputM(1))

allocate(InputM(1)%M(n1,n2))
!call mem_ichor_alloc(InputM(1)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(1)%M,KIND=long)*memrealsz
call copyM(InputM(1)%M,M)
InputM(1)%n1 = n1
InputM(1)%n2 = n2
end subroutine IchorInputInfoM1

subroutine IchorInputInfoM2(M1,n11,n12,M2,n21,n22)
implicit none
integer,intent(in) :: n11,n12,n21,n22
real(realk),intent(in) :: M1(n11,n12),M2(n21,n22)
nInputMatrixItem = 2
allocate(InputM(2))

allocate(InputM(1)%M(n11,n12))
!call mem_ichor_alloc(InputM(1)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(1)%M,KIND=long)*memrealsz
call copyM(InputM(1)%M,M1)
InputM(1)%n1 = n11
InputM(1)%n2 = n12

allocate(InputM(2)%M(n21,n22))
!call mem_ichor_alloc(InputM(2)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(2)%M,KIND=long)*memrealsz
call copyM(InputM(2)%M,M2)
InputM(2)%n1 = n21
InputM(2)%n2 = n22

end subroutine IchorInputInfoM2

subroutine IchorInputInfoM3(M1,n11,n12,M2,n21,n22,M3,n31,n32)
implicit none
integer,intent(in) :: n11,n12,n21,n22,n31,n32
real(realk),intent(in) :: M1(n11,n12),M2(n21,n22),M3(n31,n32)
nInputMatrixItem = 3
allocate(InputM(3))

allocate(InputM(1)%M(n11,n12))
!call mem_ichor_alloc(InputM(1)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(1)%M,KIND=long)*memrealsz
call copyM(InputM(1)%M,M1)
InputM(1)%n1 = n11
InputM(1)%n2 = n12

allocate(InputM(2)%M(n21,n22))
!call mem_ichor_alloc(InputM(2)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(2)%M,KIND=long)*memrealsz
call copyM(InputM(2)%M,M2)
InputM(2)%n1 = n21
InputM(2)%n2 = n22

allocate(InputM(3)%M(n31,n32))
!call mem_ichor_alloc(InputM(3)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(3)%M,KIND=long)*memrealsz
call copyM(InputM(3)%M,M3)
InputM(3)%n1 = n31
InputM(3)%n2 = n32

end subroutine IchorInputInfoM3

subroutine IchorInputInfoM4(M1,n11,n12,M2,n21,n22,M3,n31,n32,M4,n41,n42)
implicit none
integer,intent(in) :: n11,n12,n21,n22,n31,n32,n41,n42
real(realk),intent(in) :: M1(n11,n12),M2(n21,n22),M3(n31,n32),M4(n41,n42)
nInputMatrixItem = 4
allocate(InputM(4))

allocate(InputM(1)%M(n11,n12))
!call mem_ichor_alloc(InputM(1)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(1)%M,KIND=long)*memrealsz
call copyM(InputM(1)%M,M1)
InputM(1)%n1 = n11
InputM(1)%n2 = n12

allocate(InputM(2)%M(n21,n22))
!call mem_ichor_alloc(InputM(2)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(2)%M,KIND=long)*memrealsz
call copyM(InputM(2)%M,M2)
InputM(2)%n1 = n21
InputM(2)%n2 = n22

allocate(InputM(3)%M(n31,n32))
!call mem_ichor_alloc(InputM(3)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(3)%M,KIND=long)*memrealsz
call copyM(InputM(3)%M,M3)
InputM(3)%n1 = n31
InputM(3)%n2 = n32

allocate(InputM(4)%M(n41,n42))
!call mem_ichor_alloc(InputM(4)%M)
mem_allocated_inputichor = mem_allocated_inputichor + &
     & size(InputM(4)%M,KIND=long)*memrealsz
call copyM(InputM(4)%M,M4)
InputM(4)%n1 = n41
InputM(4)%n2 = n42

end subroutine IchorInputInfoM4

subroutine copyM(InputM,M)
implicit none
real(realk),intent(in) :: M(:,:)
real(realk),intent(inout) :: InputM(:,:)
integer :: n1,n2,I,J
n1 = size(M,1)
n2 = size(M,2)
do J=1,n2
   do I=1,n1
      InputM(I,J) = M(I,J)
   enddo
enddo
end subroutine copyM

subroutine IchorInputInfoSpec(CenterA,CenterB,CenterC,CenterD)
implicit none
integer :: CenterA,CenterB,CenterC,CenterD

InputMspec(1) = CenterA
IF(CenterA.GT.nInputMatrixItem.OR.CenterA.LT.1)THEN
   call IchorQuit('IchorInputInfoSpecA Error',-1)
ENDIF
InputMspec(2) = CenterB
IF(CenterB.GT.nInputMatrixItem.OR.CenterB.LT.1)THEN
   call IchorQuit('IchorInputInfoSpecA Error',-1)
ENDIF
InputMspec(3) = CenterC
IF(CenterC.GT.nInputMatrixItem.OR.CenterC.LT.1)THEN
   call IchorQuit('IchorInputInfoSpecA Error',-1)
ENDIF
InputMspec(4) = CenterD
IF(CenterD.GT.nInputMatrixItem.OR.CenterD.LT.1)THEN
   call IchorQuit('IchorInputInfoSpecA Error',-1)
ENDIF

end subroutine IchorInputInfoSpec

END MODULE ICHORINPUTINFOMOD
