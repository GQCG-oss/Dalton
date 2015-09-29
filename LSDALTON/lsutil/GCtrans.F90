!> @file 
!> contains many structure and associated subroutine
MODULE GCtransMod
 use files
 use precision
 use matrix_module
 use typedeftype
 use matrix_operations
 use memory_handling
 use basis_typetype
private
!STILL NEED TO GET UNRES WORKING AND TEST THE SCALAPACK THING (test without scalapack)
public :: init_AO2GCAO_GCAO2AO, free_AO2GCAO_GCAO2AO, &
     & AO2GCAO_transform_matrixF, AO2GCAO_half_transform_matrix, &
     & GCAO2AO_half_transform_matrix, &!GCAO2AO_transform_matrixD, &
     & GCAO2AO_transform_matrixD2, AO2GCAO_transform_matrixD, &
     & AO2GCAO_transform_fullF, GCAO2AO_transform_fullD, &
!     & GCAO2AO_transform_fullF, AO2GCAO_transform_fullD, &
     & write_GCtransformationmatrix, GCAO2AO_half_transform_matrixfull
type(matrix),save :: GCAOtrans
type(matrix),save :: GCAOtrans_inv
logical :: GCbuild

Contains

subroutine init_AO2GCAO_GCAO2AO()
GCbuild = .FALSE.
call mat_nullify(GCAOtrans)
call mat_nullify(GCAOtrans_inv)
end subroutine INIT_AO2GCAO_GCAO2AO

subroutine free_AO2GCAO_GCAO2AO()
  IF(GCbuild)THEN
     call mat_free(GCAOtrans)
     call mat_free(GCAOtrans_inv)
  ENDIF
end subroutine free_AO2GCAO_GCAO2AO

!> \brief Transform AO Fock matrix to GCAO matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param MAT Matrix to be transformed
!> \param setting the setting structure
!> \param lupri the logical unit number for output
subroutine AO2GCAO_transform_matrixF(F,setting,lupri)
implicit none
integer,intent(in) :: lupri
type(matrix) :: F
TYPE(LSSETTING) :: setting
!
integer :: nrow,ncol
type(Matrix) :: wrk,CC

nrow = F%nrow
ncol = F%ncol
IF(ncol.NE.nrow)THEN
   call lsquit('AO2GCAO_transform_matrixF requires a square matrix',lupri)
ENDIF
IF(GCbuild)THEN
   call MAT_INIT(wrk,nrow,ncol)
   call mat_mul(GCAOtrans,F,'t','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,GCAOtrans,'n','n',1E0_realk,0E0_realk,F)
   call MAT_free(wrk)
ELSE
   call mat_init(CC,ncol,nrow)
   call read_GCtransformationmatrix(CC,nrow,setting,lupri)
   call MAT_INIT(wrk,nrow,ncol)
   call mat_mul(CC,F,'t','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,CC,'n','n',1E0_realk,0E0_realk,F)
   call MAT_free(wrk)
   call MAT_free(CC)
ENDIF
end subroutine AO2GCAO_transform_matrixF

!> \brief Half-transform AO matrix to GCAO matrix
!> \author S. Reine
!> \date Dec 6th 2012
!> \param MAT Matrix to be transformed
!> \param setting the setting structure
!> \param lupri the logical unit number for output
!> \param side index indicating if first or second AO should be transformed (1 or 2)
subroutine AO2GCAO_half_transform_matrix(F,setting,lupri,side)
implicit none
integer,intent(in) :: lupri
type(matrix)       :: F
TYPE(LSSETTING)    :: setting
Integer            :: side
!
integer :: nrow,ncol,ngcao
type(Matrix) :: wrk,CC

nrow = F%nrow
ncol = F%ncol
IF (side.EQ.1) THEN
  ngcao = nrow
ELSEIF (side.EQ.2) THEN
  ngcao = ncol
ELSE
  CALL lsquit('Error in AO2GCAO_half_transform_matrix. Incorrect side.',lupri)
ENDIF

call MAT_init(CC,ngcao,ngcao)
call read_GCtransformationmatrix(CC,ngcao,setting,lupri)
call MAT_INIT(wrk,nrow,ncol)
IF (side.EQ.1) THEN
  call mat_mul(CC,F,'t','n',1E0_realk,0E0_realk,wrk)
ELSE
  call mat_mul(F,CC,'n','n',1E0_realk,0E0_realk,wrk)
ENDIF
call mat_copy(1E0_realk,wrk,F)
call MAT_free(wrk)
call MAT_free(CC)

end subroutine AO2GCAO_half_transform_matrix

!> \brief Half-transform AO matrix to GCAO matrix
!> \author S. Reine
!> \date Dec 6th 2012
!> \param MAT Matrix to be transformed
!> \param setting the setting structure
!> \param lupri the logical unit number for output
!> \param side index indicating if first or second AO should be transformed (1 or 2)
subroutine GCAO2AO_half_transform_matrixFull(fullMat,n1,n2,setting,lupri,side)
implicit none
integer,intent(in) :: lupri,n1,n2
real(realk)       :: fullMat(n1,n2)
TYPE(LSSETTING)    :: setting
Integer            :: side
!
integer :: nrow,ncol,ngcao
real(realk),pointer :: wrk(:,:),CCfull(:,:)
nrow = n1
ncol = n2
IF (side.EQ.1) THEN
  ngcao = nrow
ELSEIF (side.EQ.2) THEN
  ngcao = ncol
ELSE
  CALL lsquit('Error in GCAO2AO_half_transform_matrixfull. Incorrect side.',lupri)
ENDIF

call mem_alloc(CCfull,ngcao,ngcao)
call read_GCtransformationmatrixfull(CCfull,ngcao,setting,lupri)
call mem_alloc(wrk,nrow,ncol)
IF (side.EQ.1) THEN
   !FAO(nao,ncol) = CCfull(ngcao,ngcao)*Ffull(ngcao,ncol)  nao = ngcao
   call DGEMM('n','n',ngcao,ncol,ngcao,1E0_realk,&
        &CCfull,ngcao,fullMAT,ngcao,0E0_realk,WRK,ngcao)
ELSE
   !FAO(nrow,nao) = Ffull(nrow,ngcao)*CCfull(ngcao,ngcao)
   call DGEMM('n','t',nrow,ngcao,ngcao,1E0_realk,&
        &fullMat,nrow,CCfull,ngcao,0E0_realk,WRK,nrow)
ENDIF
call mem_dealloc(CCfull)
CALL DCOPY(n1*n2,wrk,1,fullMat,1)
call mem_dealloc(wrk)

end subroutine GCAO2AO_half_transform_matrixFull

!> \brief Half-transform GCAO matrix to AO matrix
!> \author S. Reine
!> \date May 22nd 2013
!> \param F Matrix to be transformed
!> \param setting the setting structure
!> \param lupri the logical unit number for output
!> \param side index indicating if first or second AO should be transformed (1 or 2)
subroutine GCAO2AO_half_transform_matrix(F,setting,lupri,side)
implicit none
integer,intent(in) :: lupri
type(matrix)       :: F
TYPE(LSSETTING)    :: setting
Integer            :: side
!
integer :: nrow,ncol,ngcao
type(Matrix) :: wrk,CC

nrow = F%nrow
ncol = F%ncol
IF (side.EQ.1) THEN
  ngcao = nrow
ELSEIF (side.EQ.2) THEN
  ngcao = ncol
ELSE
  CALL lsquit('Error in GCAO2AO_half_transform_matrix. Incorrect side.',lupri)
ENDIF

call mat_init(CC,ngcao,ngcao)
call read_GCtransformationmatrix(CC,ngcao,setting,lupri)
call MAT_INIT(wrk,nrow,ncol)
IF (side.EQ.1) THEN
  call mat_mul(CC,F,'n','n',1E0_realk,0E0_realk,wrk)
ELSE
  call mat_mul(F,CC,'n','t',1E0_realk,0E0_realk,wrk)
ENDIF
call mat_copy(1E0_realk,wrk,F)
call MAT_free(wrk)
call MAT_free(CC)

end subroutine GCAO2AO_half_transform_matrix

!!$!> \brief Transform GCAO Density matrix to AO Density matrix
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param MAT Matrix to be transformed
!!$!> \param setting the setting structure
!!$!> \param lupri the logical unit number for output
!!$subroutine GCAO2AO_transform_matrixD(DMAT,setting,lupri)
!!$implicit none
!!$integer,intent(in) :: lupri
!!$type(matrix) :: DMAT
!!$TYPE(LSSETTING) :: setting
!!$!
!!$type(matrix) :: wrk,CC
!!$integer :: ndmat,nbast
!!$nbast=DMAT%nrow
!!$
!!$IF(GCbuild)THEN
!!$   call MAT_INIT(wrk,GCAOtrans%nrow,DMAT%ncol)
!!$   call mat_mul(GCAOtrans,DMAT,'n','n',1E0_realk,0E0_realk,wrk)
!!$   call mat_mul(wrk,GCAOtrans,'n','t',1E0_realk,0E0_realk,DMAT)
!!$   call MAT_free(wrk)
!!$ELSE
!!$   call MAT_INIT(CC,DMAT%ncol,DMAT%nrow)
!!$   call read_GCtransformationmatrix(CC,nbast,setting,lupri)
!!$   call MAT_INIT(wrk,CC%nrow,DMAT%ncol)
!!$   call mat_mul(CC,DMAT,'n','n',1E0_realk,0E0_realk,wrk)
!!$   call mat_mul(wrk,CC,'n','t',1E0_realk,0E0_realk,DMAT)
!!$   call MAT_free(CC)
!!$   call MAT_free(wrk)
!!$ENDIF
!!$end subroutine GCAO2AO_transform_matrixD

!> \brief Transform GCAO Density matrix to AO Density matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param MAT Matrix to be transformed
!> \param setting the setting structure
!> \param lupri the logical unit number for output
subroutine GCAO2AO_transform_matrixD2(DMATGCAO,DMATAO,setting,lupri)
implicit none
integer,intent(in) :: lupri
type(matrix) :: DMATGCAO,DMATAO
TYPE(LSSETTING) :: setting
!
type(matrix) :: wrk,CC
integer :: ndmat,nbast
nbast=DMATAO%nrow

IF(GCbuild)THEN
   call MAT_INIT(wrk,GCAOtrans%nrow,DMATGCAO%ncol)
   call mat_mul(GCAOtrans,DMATGCAO,'n','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,GCAOtrans,'n','t',1E0_realk,0E0_realk,DMATAO)
   call MAT_free(wrk)
ELSE
   call MAT_INIT(CC,DMATGCAO%ncol,DMATGCAO%nrow)
   call read_GCtransformationmatrix(CC,nbast,setting,lupri)
   call MAT_INIT(wrk,CC%nrow,DMATGCAO%ncol)
   call mat_mul(CC,DMATGCAO,'n','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,CC,'n','t',1E0_realk,0E0_realk,DMATAO)
   call MAT_free(wrk)
   call MAT_free(CC)
ENDIF

end subroutine GCAO2AO_transform_matrixD2

!> \brief Transform AO Density matrix to GCAO Density matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param MAT Matrix to be transformed
!> \param setting the setting structure
!> \param lupri the logical unit number for output
subroutine AO2GCAO_transform_matrixD(DMAT,setting,lupri)
implicit none
integer,intent(in) :: lupri
type(matrix) :: DMAT
TYPE(LSSETTING) :: setting
!
integer :: ndmat,nbast
type(matrix) :: wrk,cc_inv
nbast=DMAT%nrow

IF(GCbuild)THEN
   call MAT_INIT(wrk,GCAOtrans_inv%nrow,DMAT%ncol)
   call mat_mul(GCAOtrans_inv,DMAT,'n','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,GCAOtrans_inv,'n','t',1E0_realk,0E0_realk,DMAT)
   call MAT_free(wrk)
ELSE
   call MAT_INIT(CC_inv,DMAT%ncol,DMAT%nrow)
   call read_GCtransInvformationmatrix(CC_inv,nbast,setting,lupri)
   call MAT_INIT(wrk,CC_inv%nrow,DMAT%ncol)
   call mat_mul(CC_inv,DMAT,'n','n',1E0_realk,0E0_realk,wrk)
   call mat_mul(wrk,CC_inv,'n','t',1E0_realk,0E0_realk,DMAT)
   call MAT_free(wrk)
   call MAT_free(CC_inv)
ENDIF

end subroutine AO2GCAO_transform_matrixD

!> \brief Transform 2 dim fortran array from AO to GCAO
!> \author T. Kjaergaard
!> \date 2010
!> \param fullMAT array to be transformed
!> \param nbast number of basis functions
!> \param setting the setting structure
!> \param lupri the logical unit number for output
!> This one works for S,F,H1 and other matrices that transform in the same way
!> For D se next routine
subroutine AO2GCAO_transform_fullF(fullMAT,nbast,setting,lupri)
implicit none
integer,intent(in) :: lupri,nbast
real(realk) :: fullMAT(nbast,nbast)
TYPE(LSSETTING) :: setting
!
real(realk),pointer :: CCfull(:,:)
real(realk),pointer :: WRK(:,:)

call mem_alloc(CCfull,nbast,nbast)
call read_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
call mem_alloc(wrk,nbast,nbast)
call DGEMM('t','n',nbast,nbast,nbast,1E0_realk,&
     &CCfull,nbast,fullMAT,nbast,0E0_realk,WRK,nbast)
call DGEMM('n','n',nbast,nbast,nbast,1E0_realk,&
     &WRK,nbast,CCfull,nbast,0E0_realk,fullMAT,nbast)
call mem_dealloc(wrk)
call mem_dealloc(CCfull)

end subroutine AO2GCAO_transform_fullF

!!$!> \brief Transform Density from AO to GCAO
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!> \param fullMAT array to be transformed
!!$!> \param nbast number of basis functions
!!$!> \param setting the setting structure
!!$!> \param lupri the logical unit number for output
!!$!> This one works for the density matrix and matrices that transform in the same way
!!$!> For S,F,H1  se previous routine
!!$subroutine AO2GCAO_transform_fullD(DmatAO,DmatGCAO,nbast,ndmat,setting,lupri)
!!$implicit none
!!$integer,intent(in) :: lupri,nbast,ndmat
!!$real(realk) :: DmatAO(nbast,nbast,ndmat)
!!$real(realk) :: DmatGCAO(nbast,nbast,ndmat)
!!$TYPE(LSSETTING) :: setting
!!$!
!!$integer :: i
!!$real(realk),pointer :: CCfull_inv(:,:),WRK(:,:)
!!$
!!$call mem_alloc(CCfull_inv,nbast,nbast)
!!$call read_GCtransformationmatrixfull(CCfull_inv,nbast,setting,lupri)
!!$call inv_fullmat(CCfull_inv,nbast)
!!$call mem_alloc(WRK,nbast,nbast)
!!$DO I=1,ndmat
!!$   call DGEMM('n','n',nbast,nbast,nbast,1E0_realk,&
!!$        &CCfull_inv,nbast,DmatAO(:,:,I),nbast,0E0_realk,WRK,nbast)
!!$   call DGEMM('n','t',nbast,nbast,nbast,1E0_realk,&
!!$        &WRK,nbast,CCfull_inv,nbast,0E0_realk,DmatGCAO(:,:,I),nbast)
!!$ENDDO
!!$call mem_dealloc(WRK)
!!$call mem_dealloc(CCfull_inv)
!!$
!!$end subroutine AO2GCAO_transform_fullD

!!$!> \brief Transform GCAO 3dim fortran array to AO
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$!>
!!$!> This one works for S,F,H1 and other matrices that transform in the same way
!!$!> For D se next routine
!!$!>
!!$!> \param MATGCAO GCAO array to be transformed
!!$!> \param MATAO AO the output array
!!$!> \param nbast number of basis functions
!!$!> \param setting the setting structure
!!$!> \param lupri the logical unit number for output
!!$subroutine GCAO2AO_transform_fullF(MATGCAO,MATAO,nbast,ndmat,setting,lupri)
!!$implicit none
!!$integer,intent(in) :: lupri,ndmat,nbast
!!$real(realk) :: MATGCAO(nbast,nbast,ndmat)
!!$real(realk) :: MATAO(nbast,nbast,ndmat)
!!$TYPE(LSSETTING) :: setting
!!$!
!!$integer :: I
!!$real(realk),pointer :: CCfull_inv(:,:),WRK(:,:)
!!$
!!$call mem_alloc(CCfull_inv,nbast,nbast)
!!$call read_GCtransformationmatrixfull(CCfull_inv,nbast,setting,lupri)
!!$call inv_fullmat(CCfull_inv,nbast)
!!$call mem_alloc(WRK,nbast,nbast)
!!$DO I=1,ndmat
!!$   call DGEMM('t','n',nbast,nbast,nbast,1E0_realk,&
!!$        &CCfull_inv,nbast,matgcao(:,:,I),nbast,0E0_realk,WRK,nbast)
!!$   call DGEMM('n','n',nbast,nbast,nbast,1E0_realk,&
!!$        &WRK,nbast,CCfull_inv,nbast,0E0_realk,matao(:,:,I),nbast)
!!$ENDDO
!!$call mem_dealloc(WRK)
!!$call mem_dealloc(CCfull_inv)
!!$
!!$end subroutine GCAO2AO_transform_fullF

!> \brief Transform GCAO 3dim fortran array to AO
!> \author T. Kjaergaard
!> \date 2010
!>
!> This one works for D and other matrices that transform in the same way
!> For F,S,H1 se previous routine
!>
!> \param MATGCAO GCAO array to be transformed
!> \param MATAO AO the output array
!> \param nbast number of basis functions
!> \param setting the setting structure
!> \param lupri the logical unit number for output
subroutine GCAO2AO_transform_fullD(MATGCAO,MATAO,nbast,ndmat,setting,lupri)
implicit none
integer :: lupri,ndmat,nbast
real(realk) :: MATGCAO(nbast,nbast,ndmat)
real(realk) :: MATAO(nbast,nbast,ndmat)
TYPE(LSSETTING) :: setting
!
integer :: I
integer,pointer :: IPVT(:)
real(realk) :: dummy(2),RCOND
real(realk),pointer :: CCfull(:,:),WRK(:,:),WORK1(:)

call mem_alloc(CCfull,nbast,nbast)
call read_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
call mem_alloc(WRK,nbast,nbast)
DO I=1,ndmat
   call DGEMM('n','n',nbast,nbast,nbast,1E0_realk,&
        &CCfull,nbast,matgcao(:,:,I),nbast,0E0_realk,WRK,nbast)
   call DGEMM('n','t',nbast,nbast,nbast,1E0_realk,&
        &WRK,nbast,CCfull,nbast,0E0_realk,matao(:,:,I),nbast)
ENDDO
call mem_dealloc(WRK)
call mem_dealloc(CCfull)

end subroutine GCAO2AO_transform_fullD

!> \brief Build AO,GCAO transformation matrix CC(AO,GCAO)
!> \author T. Kjaergaard
!> \date 2010
!> \param CCfull the transformation matrix to be generated
!> \param nbast the dimension, the number of basis functions
!> \param setting the setting structure
!> \param lupri the logical unit number for output
SUBROUTINE build_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
implicit none
integer,intent(in)   :: lupri
integer,intent(in)   :: nbast
TYPE(LSSETTING) :: setting
!TYPE(MATRIX)    :: CC
real(realk),intent(inout) :: CCfull(nbast,nbast)
!
integer :: nbast1,R1,I,ICHARGE,type,iang,norb,II,JJ,isp
integer :: nOrbComp,cc1,cc2,ccJJ,ccII
real(realk),parameter :: D0=0E0_realk
real(realk),pointer :: TMP(:,:)
real(realk) :: TMPs
logical :: matrix_exsist

!this could be built directly into a coordinate form and then transformed into CSR !
!so introduced special case for CSR 
DO JJ=1,nbast
   DO II=1,nbast
      CCfull(II,JJ) = D0
   ENDDO
ENDDO

nbast1 = 0
R1 = setting%BASIS(1)%p%BINFO(GCTBasParam)%Labelindex
DO I=1,setting%MOLECULE(1)%p%natoms   
   IF(setting%MOLECULE(1)%p%ATOM(I)%pointcharge)cycle
   IF(R1.EQ.0)THEN
      ICHARGE = INT(setting%MOLECULE(1)%p%ATOM(I)%CHARGE)      
      type = setting%BASIS(1)%p%BINFO(GCTBasParam)%CHARGEINDEX(ICHARGE)
   ELSE
      type = setting%MOLECULE(1)%p%ATOM(I)%IDtype(R1)
   ENDIF
   do iang = 1,setting%BASIS(1)%p%BINFO(GCTBasParam)%ATOMTYPE(type)%nAngmom
      norb = setting%BASIS(1)%p%BINFO(GCTBasParam)%ATOMTYPE(type)%SHELL(iang)%segment(1)%ncol
      nOrbComp = (2*(iang-1))+1   !not true if cartesian
      call mem_alloc(TMP,norb,norb)
      DO JJ=1,norb
         DO II=1,norb
            TMP(II,JJ)=&
                 &setting%BASIS(1)%p%BINFO(GCTBasParam)%ATOMTYPE(type)%SHELL(iang)%segment(1)%elms(II+(JJ-1)*norb)
         ENDDO
      ENDDO
      nOrbComp = (2*(iang-1))+1   !not true if cartesian
      DO JJ=1,norb
         ccJJ = nbast1+(JJ-1)*nOrbComp
         DO II=1,norb
            ccII = nbast1+(II-1)*nOrbComp
            TMPs = TMP(II,JJ)
            do isp = 1,nOrbComp
               cc1 = isp+ccII
               cc2 = isp+ccJJ
               CCfull(cc1,cc2)=TMPs
            ENDDO
         ENDDO !II
      ENDDO !JJ
      call mem_dealloc(TMP)
      nbast1 = nbast1 + norb*nOrbComp
   ENDDO
enddo
IF(nbast1.NE.nbast)&
     & call lsquit('dim mismatch in build_GCtransformationmatrix',-1)
end SUBROUTINE build_GCtransformationmatrixfull

!> \brief Build AO,GCAO transformation matrix CC(AO,GCAO)
!> \author T. Kjaergaard
!> \date 2010
!> \param CCfull the transformation matrix to be generated
!> \param nbast the dimension, the number of basis functions
!> \param setting the setting structure
!> \param lupri the logical unit number for output
SUBROUTINE write_GCtransformationmatrixfull(nbast,setting,lupri)
implicit none
integer,intent(in)   :: lupri
integer   :: nbast
TYPE(LSSETTING) :: setting
!
real(realk),pointer :: CCfull(:,:)
logical :: matrix_exsist
integer :: GCAOtrans_lun

INQUIRE(file='GCAOtrans',EXIST=matrix_exsist) 
IF(matrix_exsist)then
   GCAOtrans_lun=-1
   call lsopen(GCAOtrans_lun,'GCAOtrans','OLD','UNFORMATTED')
   call lsclose(GCAOtrans_lun,'DELETE')
ENDIF
call mem_alloc(CCfull,nbast,nbast)
call build_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
GCAOtrans_lun = -1  !initialization
call lsopen(GCAOtrans_lun,'GCAOtrans','UNKNOWN','UNFORMATTED')
rewind GCAOtrans_lun
WRITE(GCAOtrans_lun) nbast
WRITE(GCAOtrans_lun) CCfull

!print*,'write_GCtransformationmatrix nbast',nbast
!call ls_output(CCfull,1,nbast,1,nbast,nbast,nbast,1,6)

call lsclose(GCAOtrans_lun,'KEEP')
call mem_dealloc(CCfull)
end SUBROUTINE write_GCtransformationmatrixfull

!> \brief Build AO,GCAO transformation matrix CC(AO,GCAO)
!> \author T. Kjaergaard
!> \date 2010
!> \param CCfull the transformation matrix to be generated
!> \param nbast the dimension, the number of basis functions
!> \param setting the setting structure
!> \param lupri the logical unit number for output
SUBROUTINE write_GCtransformationmatrix(nbast,setting,lupri)
implicit none
integer,intent(in)   :: lupri
integer   :: nbast
TYPE(LSSETTING) :: setting
real(realk),pointer :: CCfull(:,:)
IF(matrix_type .EQ. mtype_scalapack)THEN
   IF(GCbuild)THEN
      !matrix already buildt
      call mat_free(GCAOtrans)
      call mat_free(GCAOtrans_inv)
      GCbuild = .FALSE.
   ENDIF
   call mem_alloc(CCfull,nbast,nbast)
   call read_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
   
   call mat_init(GCAOtrans,nbast,nbast)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(nbast*nbast,CCfull,1,GCAOtrans%elms,1)
      CALL DCOPY(nbast*nbast,CCfull,1,GCAOtrans%elmsb,1)
   ELSE
      call mat_set_from_full(CCfull,1E0_realk,GCAOtrans)
   ENDIF
   call mat_init(GCAOtrans_inv,nbast,nbast)
   call inv_fullmat(CCfull,nbast)      
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(nbast*nbast,CCfull,1,GCAOtrans_inv%elms,1)
      CALL DCOPY(nbast*nbast,CCfull,1,GCAOtrans_inv%elmsb,1)
   ELSE
      call mat_set_from_full(CCfull,1E0_realk,GCAOtrans_inv)
   ENDIF
   call mem_dealloc(CCfull)
   GCbuild = .TRUE.
!   WRITE(lupri,*)'THE GCAO'
!   call mat_print(GCAOtrans,1,GCAOtrans%nrow,1,GCAOtrans%ncol,lupri)
!   WRITE(lupri,*)'THE GCAO_inv'
!   call mat_print(GCAO,1,GCAO%nrow,1,GCAO%ncol,lupri)
ELSE
   call write_GCtransformationmatrixfull(nbast,setting,lupri)
ENDIF
end SUBROUTINE write_GCtransformationmatrix

!> \brief Build AO,GCAO transformation matrix CC(AO,GCAO)
!> \author T. Kjaergaard
!> \date 2010
!> \param CCfull the transformation matrix to be generated
!> \param nbast the dimension, the number of basis functions
!> \param setting the setting structure
!> \param lupri the logical unit number for output
SUBROUTINE read_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
implicit none
integer,intent(in)   :: lupri
integer,intent(in)   :: nbast
integer   :: nbast2
TYPE(LSSETTING) :: setting
real(realk),intent(inout) :: CCfull(nbast,nbast)
!
logical :: matrix_exsist
integer :: GCAOtrans_lun

INQUIRE(file='GCAOtrans',EXIST=matrix_exsist) 
IF(matrix_exsist)then
   GCAOtrans_lun = -1  !initialization
   call lsopen(GCAOtrans_lun,'GCAOtrans','OLD','UNFORMATTED')
   rewind GCAOtrans_lun
   READ(GCAOtrans_lun) nbast2
   IF(nbast2.NE.nbast) THEN
     write(lupri,'(A,2I5)') 'Error in read_GCtransformationmatrix - dim mismatch:',nbast,nbast2
     call lsquit('dim mismatch read_GCtransformationmatrix',-1)
   ENDIF
   READ(GCAOtrans_lun) CCfull
   call lsclose(GCAOtrans_lun,'KEEP')
ELSE
   call build_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
ENDIF
!print*,'read_GCtransformationmatrixfull nbast',nbast
!call ls_output(CCfull,1,nbast,1,nbast,nbast,nbast,1,6)
end SUBROUTINE read_GCtransformationmatrixfull

SUBROUTINE read_GCtransformationmatrix(CCmatrix,nbast,setting,lupri)
implicit none
integer,intent(in)   :: lupri
integer,intent(in)   :: nbast
integer   :: nbast2
TYPE(LSSETTING) :: setting
type(matrix),intent(inout) :: CCmatrix
real(realk),pointer :: CCfull(:,:)
IF(GCbuild)THEN
   call mat_assign(CCmatrix,GCAOtrans)
ELSE
   call mem_alloc(CCfull,nbast,nbast)
   call read_GCtransformationmatrixfull(CCfull,nbast,setting,lupri)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(CCmatrix%nrow*CCmatrix%ncol,CCfull,1,CCmatrix%elms,1)
      CALL DCOPY(CCmatrix%nrow*CCmatrix%ncol,CCfull,1,CCmatrix%elmsb,1)
   ELSE
      call mat_set_from_full(CCfull,1E0_realk,CCmatrix)
   ENDIF
   call mem_dealloc(CCfull)
ENDIF
end SUBROUTINE read_GCtransformationmatrix

SUBROUTINE read_GCtransInvformationmatrix(CCmatrix_inv,nbast,setting,lupri)
implicit none
integer,intent(in)   :: lupri
integer,intent(in)   :: nbast
integer   :: nbast2
TYPE(LSSETTING) :: setting
type(matrix),intent(inout) :: CCmatrix_inv
real(realk),pointer :: CCfull_inv(:,:)
IF(GCbuild)THEN
   call mat_assign(CCmatrix_inv,GCAOtrans_inv)
ELSE
   call mem_alloc(CCfull_inv,nbast,nbast)
   call read_GCtransformationmatrixfull(CCfull_inv,nbast,setting,lupri)
   call inv_fullmat(CCfull_inv,nbast)
   IF(matrix_type .EQ. mtype_unres_dense)THEN
      CALL DCOPY(nbast*nbast,CCfull_inv,1,CCmatrix_inv%elms,1)
      CALL DCOPY(nbast*nbast,CCfull_inv,1,CCmatrix_inv%elmsb,1)
   ELSE
      call mat_set_from_full(CCfull_inv,1E0_realk,CCmatrix_inv)
   ENDIF
   call mem_dealloc(CCfull_inv)
ENDIF
end SUBROUTINE read_GCtransInvformationmatrix

subroutine inv_fullmat(CCfull_inv,nbast)
implicit none
integer :: nbast
real(realk),intent(inout) :: CCfull_inv(nbast,nbast)
integer :: i
integer,pointer :: IPVT(:)
real(realk) :: dummy(2),RCOND
real(realk),pointer :: WRK(:,:),WORK1(:)
call mem_alloc(IPVT,nbast)
call mem_alloc(WORK1,nbast)
IPVT = 0; RCOND = 0E0_realk; dummy = 0E0_realk
call DGECO(CCfull_inv,nbast,nbast,IPVT,RCOND,work1)
call DGEDI(CCfull_inv,nbast,nbast,IPVT,dummy,work1,01)!01=inverse only
call mem_dealloc(IPVT)
call mem_dealloc(WORK1)
end subroutine inv_fullmat

END MODULE GCTRANSMod



