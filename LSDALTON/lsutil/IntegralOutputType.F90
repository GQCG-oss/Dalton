!> @file 
!> Contains OBJECT CONTAINING INFORMATION ABOUT THE INTEGRAL OUTPUT
MODULE integraloutput_typetype
  use precision
  use lstensor_typetype

type intbatch
integer,pointer :: BATCH(:)
end type intbatch

TYPE INTEGRALOUTPUT
#ifdef VAR_PTR_RESHAPE
REAL(REALK),contiguous,pointer    :: ResultMat(:,:,:,:,:) 
REAL(REALK),contiguous,pointer    :: Result3D(:,:,:) 
#else
REAL(REALK),pointer    :: ResultMat(:,:,:,:,:) 
REAL(REALK),pointer    :: Result3D(:,:,:) 
#endif
type(lstensor),pointer :: screenTensor
type(lstensor),pointer :: resultTensor
type(lstensor)         :: RHScont
Integer,pointer        :: postprocess(:) !see parameters in ls_parameters.f90
LOGICAL                :: memdistResultTensor
LOGICAL                :: doGRAD
Integer                :: ndim(5)
Integer                :: ndim3D(4)
! buffer 
LOGICAL                :: USEBUFMM    ! flag for using/not using a buffer 
                                      !to write multipole moments to file
INTEGER,pointer        :: IBUF(:,:)   ! integer buffer
REAL(REALK),pointer    :: RBUF(:,:)   ! real buffer
REAL(REALK),pointer    :: NBUF(:,:)   ! real buffer for nuclear position information
INTEGER                :: MMBUFLEN    ! length of the integer buffer is : MMBUFLEN*MAXBUFI;  
                                      !        of the real buffer MMBUFLEN*MAXBUFR;
                                      !        of the nuclear position buffer  MMBUFLEN*MAXBUFN
INTEGER                :: MAXBUFI,MAXBUFR,MAXBUFN
INTEGER                :: IBUFI       ! counter for the real and integer buffer
INTEGER                :: IBUFN       ! counter for the nuclear buffer
INTEGER                :: LUITNM, LUITNMR ! logical units for the integer and the real buffer files
! end buffer
logical                :: decpacked
logical                :: decpacked2
logical                :: decpackedK
logical                :: FullAlphaCD
logical                :: RealGabMatrix
real(realk)            :: exchangeFactor
END TYPE INTEGRALOUTPUT

public :: intbatch,INTEGRALOUTPUT

private 

contains

!Added to avoid "has no symbols" linking warning
subroutine integraloutput_typetype_void()
end subroutine integraloutput_typetype_void

end MODULE integraloutput_typetype
