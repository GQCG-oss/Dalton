!> @file 
!> Contains OBJECT CONTAINING INFORMATION ABOUT THE INTEGRAL OUTPUT
MODULE integraloutput_type
  use integraloutput_typetype
  use precision
  use lstensor_typetype
  use memory_handling
CONTAINS
SUBROUTINE nullifyIntegralOutput(IntOut)
implicit none
TYPE(INTEGRALOUTPUT) :: IntOut
NULLIFY(IntOut%resultTensor)
NULLIFY(IntOut%ResultMat)
NULLIFY(IntOut%Result3D)
NULLIFY(IntOut%screenTensor)
NULLIFY(IntOut%postprocess)
NULLIFY(IntOut%IBUF)
NULLIFY(IntOut%RBUF)
NULLIFY(IntOut%NBUF)
IntOut%decpacked = .FALSE.
IntOut%decpacked2 = .FALSE.
IntOut%decpackedK = .FALSE.
IntOut%FullAlphaCD = .FALSE.
IntOut%exchangefactor = 0E0_realk
IntOut%ndim(1) = 0
IntOut%ndim(2) = 0
IntOut%ndim(3) = 0
IntOut%ndim(4) = 0
IntOut%ndim(5) = 0
IntOut%ndim3D(1) = 0
IntOut%ndim3D(2) = 0
IntOut%ndim3D(3) = 0
IntOut%memdistResultTensor = .FALSE.
IntOut%doGRAD = .FALSE.
IntOut%RealGabMatrix = .FALSE.
END SUBROUTINE nullifyIntegralOutput

!> \brief set the dimensions of the integral output structure
!> \author T. Kjaergaard
!> \date 2010
!> \param IntOut the integraloutput structure to be initialised
!> \param dim1 size og dimension 1
!> \param dim2 size og dimension 2
!> \param dim3 size og dimension 3
!> \param dim4 size og dimension 4
!> \param dim5 size og dimension 5
SUBROUTINE initIntegralOutputDims(IntOut,dim1,dim2,dim3,dim4,dim5)
implicit none
TYPE(INTEGRALOUTPUT) :: IntOut
INTEGER              :: dim1,dim2,dim3,dim4,dim5

NULLIFY(IntOut%resultTensor)
IntOut%decpacked = .FALSE.
IntOut%decpacked2 = .FALSE.
IntOut%decpackedK = .FALSE.
IntOut%FullAlphaCD = .FALSE.
call mem_alloc(IntOut%postprocess,dim5)
IntOut%postprocess = 0
call initIntegralOutputDims1(IntOut,dim1,dim2,dim3,dim4,dim5)
IntOut%RealGabMatrix = .FALSE.

END SUBROUTINE initIntegralOutputDims

!> \brief set the dimensions of the integral output structure
!> \author T. Kjaergaard
!> \date 2010
!> \param IntOut the integraloutput structure to be initialised
!> \param dim1 size og dimension 1
!> \param dim2 size og dimension 2
!> \param dim3 size og dimension 3
!> \param dim4 size og dimension 4
!> \param dim5 size og dimension 5
SUBROUTINE initIntegralOutputDims1(IntOut,dim1,dim2,dim3,dim4,dim5)
implicit none
TYPE(INTEGRALOUTPUT) :: IntOut
INTEGER              :: dim1,dim2,dim3,dim4,dim5

IntOut%decpacked = .FALSE.
IntOut%decpacked2 = .FALSE.
IntOut%decpackedK = .FALSE.
IntOut%FullAlphaCD = .FALSE.
IntOut%exchangefactor = 0E0_realk
IntOut%ndim(1) = dim1
IntOut%ndim(2) = dim2
IntOut%ndim(3) = dim3
IntOut%ndim(4) = dim4
IntOut%ndim(5) = dim5
IntOut%memdistResultTensor = .FALSE.
END SUBROUTINE initIntegralOutputDims1

end MODULE integraloutput_type
