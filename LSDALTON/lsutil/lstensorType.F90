!> @file 
!> Contains lstensor structure and associated subroutines
MODULE LSTENSOR_TYPETYPE
  use precision
  TYPE SLSAOTENSOR
     INTEGER :: nelms
     INTEGER(kind=short),pointer :: selms(:)
     !DIM: nAOBATCH(A:B),1:2
     INTEGER,pointer :: nOrb(:)           ! maxBat,2
     INTEGER,pointer :: startLocalOrb(:)  ! maxBat,2
     INTEGER :: maxBat
     INTEGER :: nLocal(2) 
     INTEGER :: ATOM(2)
     INTEGER :: AOBATCH(2)
  END TYPE SLSAOTENSOR

  TYPE LSAOTENSOR
     INTEGER                 :: nelms
     REAL(REALK),pointer     :: elms(:)  
     !DIM: nAngmom(A:D),nAOBATCH(A:D),1:4
     INTEGER,pointer :: nOrb(:)          ! maxAng,maxBat,4
     INTEGER,pointer :: startLocalOrb(:) ! maxAng,maxBat,4
     INTEGER,pointer :: startGlobalOrb(:)! maxAng,maxBat,4
     !DIM: nAOBATCH(A:D),1:4
     INTEGER,pointer :: nAngmom(:)       ! maxBat,4    
     INTEGER :: maxBat 
     INTEGER :: maxAng
     INTEGER :: nLocal(4)
     INTEGER :: ATOM(4)
     INTEGER :: AOBATCH(4)
#ifdef VAR_MPI
     INTEGER :: GLSAOindex
#endif
  END TYPE LSAOTENSOR

  TYPE GLOBALLSAOTENSOR !so far this only works for matrices (2D)
     INTEGER :: G_startGlobalOrb(2)
     INTEGER :: G_nLocal(2)
!     INTEGER :: G_ATOM(2) 
!     INTEGER :: G_AOBATCH(2)
     integer :: mynum
  END TYPE GLOBALLSAOTENSOR

  TYPE LSTENSOR
#ifdef VAR_MPI
     TYPE(GLOBALLSAOTENSOR),pointer :: G_LSAO(:)
#endif
     TYPE(LSAOTENSOR),pointer :: LSAO(:)
     TYPE(SLSAOTENSOR),pointer :: SLSAO(:)
     INTEGER :: nAtom(4)
     INTEGER :: nbast(4)
     INTEGER :: ndim5 !nmat or nderivative*nmat or ...
     INTEGER :: nbatches(4)
     INTEGER :: nLSAO
     INTEGER :: nSLSAO
     INTEGER,pointer :: INDEX(:,:,:,:) !GIVES THE INDEX IN LSAO (or zero)
     INTEGER,pointer :: nAOBATCH(:,:) !maxNatom,4
     logical :: MagGradienttensor
     logical :: Gradienttensor
     logical :: pChargetensor
     !for instance maxgab
     logical :: Screentensor
     logical :: Screenoutput
     logical :: Econtrib ! means (1,1,1,1,ndim5) 
     logical :: LowerDiagZero
     logical :: PermuteResultTensor
     integer(kind=short),pointer :: maxgab(:,:)     !nbatches(1),nbatches(1) 
     integer(kind=short),pointer :: maxprimgab(:,:) !nbatches(1),nbatches(1) 
     integer(kind=short) :: maxgabelm
     integer(kind=short) :: maxprimgabelm
     real(realk),pointer :: MBIE(:,:,:) !nMBIE,nbatches(1),nbatches(1)
     integer :: nMBIE
     integer :: lstmem_index
#ifdef VAR_MPI
     INTEGER :: G_maxnLocal
     INTEGER :: G_nLSAO
     INTEGER :: G_nAtom(2)
     INTEGER :: G_nbast(2)
     INTEGER,pointer :: G_INDEX(:,:) !GIVES THE INDEX IN LSAO (or zero)
     INTEGER,pointer :: G_nAOBATCH(:,:) !maxNatom,2
     INTEGER,pointer :: G_fullatoms1(:) !size nAtom(1) atom in full tensor
     INTEGER,pointer :: G_fullatoms2(:) !size nAtom(2) atom in full tensor
#endif
  END TYPE LSTENSOR

  TYPE LSTENSORP
    TYPE(LSTENSOR),pointer :: p
  END TYPE LSTENSORP

contains

!!Added to avoid "has no symbols" linking warning
!subroutine LSTENSOR_TYPETYPE_void()
!end subroutine LSTENSOR_TYPETYPE_void

subroutine nullifyTENSORLSAO(LSAO)
implicit none
TYPE(LSAOTENSOR) :: LSAO(:)
!
integer :: nsize,I
nsize = SIZE(LSAO)
do I=1,nsize
   LSAO(I)%nelms = 0 
   nullify(LSAO(I)%elms)
   nullify(LSAO(I)%nOrb)
   nullify(LSAO(I)%startLocalOrb)
   nullify(LSAO(I)%startGlobalOrb)
   nullify(LSAO(I)%nAngmom)
   LSAO(I)%maxBat = 0 
   LSAO(I)%maxAng  = 0 
   LSAO(I)%nLocal = 0 
   LSAO(I)%ATOM = 0 
   LSAO(I)%AOBATCH = 0 
#ifdef VAR_MPI
   LSAO(I)%GLSAOindex = 0 
#endif
enddo
end subroutine nullifyTENSORLSAO

END MODULE LSTENSOR_TYPETYPE
