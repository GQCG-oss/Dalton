!> @file
!> Contains common routine for the Ichor Code
!> \brief Contains common routine for the Ichor Code
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorBatchToolsMod
use IchorprecisionMod
public:: DetermineBatchSize, ConstructBatchIndexOfType, ObtainMaxGabForType
private
CONTAINS
subroutine DetermineBatchSize(nTypesA,nAtomsOfTypeA,nBatchA)
implicit none
!> nTypesA is the number of different types of shells, each type is defined by 
!> an angular momentum, a number of primitives(nPrim), a number of contracted functions
!> (nCont) a set of exponents and a set of contraction coefficients. [For Center A]
integer,intent(in) :: nTypesA
!> nAtomsOfTypeA is the number of atoms that have the given type. [For Center A]
Integer,intent(in) :: nAtomsOfTypeA(ntypesA)
!> nBatchA the number of integralbatches
Integer,intent(inout) :: nBatchA
! Local variables
integer :: ItypeA
nBatchA = 0
DO ItypeA=1,nTypesA
   nBatchA = nBatchA + nAtomsOfTypeA(ItypeA)
ENDDO
END subroutine DETERMINEBATCHSIZE

subroutine ObtainMaxGabForType(MaxGabForTypeAB,nTypesA,nTypesB,nAtomsOfTypeA,nAtomsOfTypeB,&
     & BatchIndexOfTypeA,BatchIndexOfTypeB,BATCHGAB,nBatchA,nBatchB)
implicit none
integer,intent(in) :: nTypesB,nTypesA,nBatchA,nBatchB
Integer,intent(in) :: nAtomsOfTypeB(ntypesB),nAtomsOfTypeA(ntypesA)
Integer,intent(in) :: BatchIndexOfTypeB(nTypesB)
Integer,intent(in) :: BatchIndexOfTypeA(nTypesA)
real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB)
real(realk),intent(inout) :: MaxGabForTypeAB(nTypesA,nTypesB)
! Local variables
integer :: ItypeA,ItypeB,iBatchA,iBatchB,iAtomA,iAtomB
real(realk) :: MaxGabElm
DO ItypeB=1,nTypesB
   DO ItypeA=1,nTypesA
      MaxGabElm = 0.0E0_realk
      iBatchB = BatchIndexOfTypeB(ItypeB)
      DO IatomB = 1,nAtomsOfTypeB(ItypeB)
         iBatchB = iBatchB + 1
         iBatchA = BatchIndexOfTypeA(ItypeA)
         DO IatomA = 1,nAtomsOfTypeA(ItypeA)
            MaxGabElm = MAX(MaxGabElm,BATCHGAB(iBatchA + IAtomA,iBatchB))
         ENDDO
      ENDDO
      MaxGabForTypeAB(ITypeA,ITypeB) = MaxGabElm
   ENDDO
ENDDO
END subroutine OBTAINMAXGABFORTYPE

subroutine ConstructBatchIndexOfType(BatchIndexOfType,nTypes,nAtomsOfType,nBatch)
implicit none 
integer,intent(in) :: nTypes
integer,intent(in) :: nAtomsOfType(nTypes)
integer,intent(inout) :: nBatch
integer,intent(inout) :: BatchIndexOfType(nTypes)
!local variables 
integer :: Itype
nBatch = 0
DO Itype=1,nTypes
   BatchIndexOfType(Itype) = nBatch
   nBatch = nBatch + nAtomsOfType(Itype)
ENDDO
END subroutine CONSTRUCTBATCHINDEXOFTYPE

END MODULE IchorBatchToolsMod
